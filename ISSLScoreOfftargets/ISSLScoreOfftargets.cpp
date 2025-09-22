#include "ISSLScoreOfftargets.hpp"
#include "manprof.hpp"
#ifdef MANPROF_ENABLE
std::FILE* g_profileLogFile = nullptr;
std::chrono::steady_clock::time_point g_progStart{};
#endif

// ---- Manual profiling configuration ----
#ifndef MANPROF_SAMPLE_PERIOD
#define MANPROF_SAMPLE_PERIOD 200  // sample 1 in 200 guides for per-loop timings
#endif

#ifndef MANPROF_ENABLE_LOOP_SAMPLING
#define MANPROF_ENABLE_LOOP_SAMPLING 1  // set to 0 to disable loop timings entirely
#endif
// ----------------------------------------

#ifdef _WIN32
#include <windows.h>
#endif

#include <atomic>
#include <random>
#include <iomanip>

#ifndef USE_OPENCL_SEQ2SIG
#define USE_OPENCL_SEQ2SIG 0
#endif
#if USE_OPENCL_SEQ2SIG
#include "cl_seqsig.hpp"
#endif

struct alignas(64) ProfStats {
    // Times in nanoseconds to keep addition cheap & precise
    uint64_t seq_to_sig_ns = 0;
    uint64_t scoring_total_ns = 0;

    // Work counters
    uint64_t guides_scored = 0;
    uint64_t signatures_seen = 0;    // total entries traversed across all slices
    uint64_t offtargets_scored = 0;  // unique offtargets actually scored (dist<=4 gate)
    uint64_t slices_traversed = 0;

    // Optional sampled loop timings (ns)
    uint64_t sampled_outer_iter_ns = 0;
    uint64_t sampled_inner_loop_ns = 0;
    uint64_t sampled_slice_loop_ns = 0;
    uint64_t samples_outer = 0;
    uint64_t samples_inner = 0;
    uint64_t samples_slice = 0;
};

// Global aggregates (reduced after the parallel region)
static ProfStats g_prof{};

using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::pair;
using std::unordered_map;

// Char to binary encoding
const vector<uint8_t> nucleotideIndex{ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,3 };
// Binary to char encoding
const vector<char> signatureIndex{ 'A', 'C', 'G', 'T' };

uint64_t sequenceToSignature(const std::string& seq, uint64_t seqLen)
{
    uint64_t signature = 0;
    for (uint64_t j = 0; j < seqLen; j++) {
        signature |= static_cast<uint64_t>(nucleotideIndex[seq[j]]) << (j * 2);
    }
    return signature;
}

string signatureToSequence(uint64_t sig, uint64_t seqLen)
{
    string sequence = string(seqLen, ' ');
    for (uint64_t j = 0; j < seqLen; j++) {
        sequence[j] = signatureIndex[(sig >> (j * 2)) & 0x3];
    }
    return sequence;
}

int main(int argc, char** argv)
{
    g_progStart = std::chrono::steady_clock::now();

    // Default log file name
    const char* traceFileName = "profile_trace.txt";
    if (argc >= 7) {
        traceFileName = argv[6];
    }

    g_profileLogFile = std::fopen(traceFileName, "w");
    if (!g_profileLogFile) {
        fprintf(stderr, "Error: could not open trace log file '%s' for writing\n", traceFileName);
        exit(1);
    }

    std::fprintf(g_profileLogFile,
    "metric,value,notes\n"
    "TOTAL_EXECUTION_TIME_S,0,placeholder\n"
    "SEQ_TO_SIG_TOTAL_S,0,\n"
    "SCORING_TOTAL_S,0,\n"
    "GUIDES_SCORED,0,\n"
    "SIGNATURES_SEEN,0,\n"
    "OFFTARGETS_SCORED,0,\n"
    "SLICES_TRAVERSED,0,\n"
    "SAMPLED_OUTER_ITER_AVG_MS,0,\n"
    "SAMPLED_INNER_LOOP_AVG_MS,0,\n"
    "SAMPLED_SLICE_LOOP_AVG_MS,0,\n"
);

    static char s_traceBuf[1 << 20]; // 1 MB
    setvbuf(g_profileLogFile, s_traceBuf, _IOFBF, sizeof(s_traceBuf));

    TRACE_EVT("PROGRAM", "START");
    std::cout << "Trace output will be written to: " << traceFileName << std::endl;

    auto startLoading = std::chrono::high_resolution_clock::now();

    if (argc < 6) {
        fprintf(stderr,
            "Usage: %s [issltable] [query file] [max distance] [score-threshold] [score-method] [optional-tracefile]\n",
            argv[0]);
        exit(1);
    }

    if (g_profileLogFile) {
        auto now   = std::chrono::system_clock::now();
        std::time_t now_c = std::chrono::system_clock::to_time_t(now);

        // Thread-safe localtime: localtime_s on Windows, localtime_r elsewhere
        std::tm tm_buf{};
        #ifdef _WIN32
            localtime_s(&tm_buf, &now_c);
        #else
            localtime_r(&now_c, &tm_buf);
        #endif

        // Format "YYYY-MM-DD HH:MM:SS"
        char ts_base[32];
        std::strftime(ts_base, sizeof(ts_base), "%F %T", &tm_buf);

        // Milliseconds
        auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                    now.time_since_epoch()) % 1000;

        // Final "YYYY-MM-DD HH:MM:SS.mmm"
        char ts_full[40];
        std::snprintf(ts_full, sizeof(ts_full), "%s.%03lld",
                    ts_base, static_cast<long long>(ms.count()));

        // Write to file and flush
        std::fprintf(g_profileLogFile, "Beginning Program Execution, %s\n", ts_full);
        std::fflush(g_profileLogFile);

        // Echo to terminal
        std::cout << "Beginning program execution " << ts_full << std::endl;
    }

    /** The maximum number of mismatches */
    int maxDist = atoi(argv[3]);

    /** The threshold used to exit scoring early */
    double threshold = atof(argv[4]);

    /** Scoring methods. To exit early:
     *      - only CFD must drop below `threshold`
     *      - only MIT must drop below `threshold`
     *      - both CFD and MIT must drop below `threshold`
     *      - CFD or MIT must drop below `threshold`
     *      - the average of CFD and MIT must below `threshold`
     */
    string argScoreMethod = argv[5];
    otScoreMethod scoreMethod;
    bool calcCfd = false;
    bool calcMit = false;
    if (!argScoreMethod.compare("and")) {
        scoreMethod = otScoreMethod::mitAndCfd;
        calcCfd = true;
        calcMit = true;
    }
    else if (!argScoreMethod.compare("or")) {
        scoreMethod = otScoreMethod::mitOrCfd;
        calcCfd = true;
        calcMit = true;
    }
    else if (!argScoreMethod.compare("avg")) {
        scoreMethod = otScoreMethod::avgMitCfd;
        calcCfd = true;
        calcMit = true;
    }
    else if (!argScoreMethod.compare("mit")) {
        scoreMethod = otScoreMethod::mit;
        calcMit = true;
    }
    else if (!argScoreMethod.compare("cfd")) {
        scoreMethod = otScoreMethod::cfd;
        calcCfd = true;
    }
    else
    {
        fprintf(stderr, "Invalid scoring method. Acceptable options are: 'and', 'or', 'avg', 'mit', 'cfd'");
        exit(1);
    }

    /** Begin reading the binary encoded ISSL, structured as:
     *  - The header (3 items)
     *  - All binary-encoded off-target sites
     *  - Slice masks
     *  - Size of slice 1 lists
     *  - Contents of slice 1 lists
     *  ...
     *  - Size of slice N lists (N being the number of slices)
     *  - Contents of slice N lists
     */
    FILE* isslFp = fopen(argv[1], "rb");

    if (isslFp == NULL) {
        throw std::runtime_error("Error reading index: could not open file\n");
    }

    /** The index contains a fixed-sized header
     *      - the number of unique off-targets in the index
     *      - the length of an off-target
     *      - the number of slices
     */
    vector<size_t> slicelistHeader(3);

    if (fread(slicelistHeader.data(), sizeof(size_t), slicelistHeader.size(), isslFp) == 0) {
        throw std::runtime_error("Error reading index: header invalid\n");
    }

    size_t offtargetsCount = slicelistHeader[0];
    size_t seqLength = slicelistHeader[1];
    size_t sliceCount = slicelistHeader[2];

    /** Load in all of the off-target sites */
    vector<uint64_t> offtargets(offtargetsCount);
    if (fread(offtargets.data(), sizeof(uint64_t), offtargetsCount, isslFp) == 0) {
        throw std::runtime_error("Error reading index: loading off-target sequences failed\n");
    }

    /** Read the slice masks and generate 2 bit masks */
    vector<vector<uint64_t>> sliceMasks;
    for (size_t i = 0; i < sliceCount; i++)
    {
        uint64_t maskBinary;
        fread(&maskBinary, sizeof(uint64_t), 1, isslFp);

        vector<uint64_t> mask;
        for (uint64_t j = 0; j < seqLength; j++)
        {
            if (maskBinary & (1ULL << j))
            {
                mask.push_back(j);
            }
        }
        sliceMasks.push_back(mask);
    }

    /** Calculate the total number of lists based on slice count and width */
    size_t sliceListCount = 0;
    for (size_t i = 0; i < sliceCount; i++)
    {
        sliceListCount += 1ULL << (sliceMasks[i].size() * 2);
    }

    /** The number of signatures embedded per slice. Store continguously */
    vector<size_t> allSlicelistSizes(sliceListCount);

    /** The contents of the slices. Stored contiguously
     *  Each signature (64-bit) is structured as:
     *      <occurrences 32-bit><off-target-id 32-bit>
     */
    vector<uint64_t> allSignatures(offtargetsCount * sliceCount);

    /** The number of signatures embedded per slice. Store continguously */
    sliceListCount = 0;
    for (size_t i = 0; i < sliceCount; i++)
    {
        size_t sliceListSize = 1ULL << (sliceMasks[i].size() * 2);
        if (fread(allSlicelistSizes.data() + sliceListCount, sizeof(size_t), sliceListSize, isslFp) == 0) {
            throw std::runtime_error("Error reading index: reading slice list sizes failed\n");
        }

        if (fread(allSignatures.data() + (offtargetsCount * i), sizeof(uint64_t), offtargetsCount, isslFp) == 0) {
            throw std::runtime_error("Error reading index: reading slice contents failed\n");
        }

        sliceListCount += 1ULL << (sliceMasks[i].size() * 2);
    }

    /** End reading the index */
    fclose(isslFp);

    /** Prevent assessing an off-target site for multiple slices
     *
     *      Create enough 1-bit "seen" flags for the off-targets
     *      We only want to score a candidate guide against an off-target once.
     *      The least-significant bit represents the first off-target
     *      0 0 0 1   0 1 0 0   would indicate that the 3rd and 5th off-target have been seen.
     *      The CHAR_BIT macro tells us how many bits are in a byte (C++ >= 8 bits per byte)
     */
    uint64_t numOfftargetToggles = (offtargetsCount / ((size_t)sizeof(uint64_t) * (size_t)CHAR_BIT)) + 1;

    /** Start constructing index in memory
     *
     *      To begin, reverse the contiguous storage of the slices,
     *         into the following:
     *
     *         + Slice 0 :
     *         |---- AAAA : <slice contents>
     *         |---- AAAC : <slice contents>
     *         |----  ...
     *         |
     *         + Slice 1 :
     *         |---- AAAA : <slice contents>
     *         |---- AAAC : <slice contents>
     *         |---- ...
     *         | ...
     */

    vector<vector<uint64_t*>> sliceLists(sliceCount);
    // Assign sliceLists size based on each slice length
    for (size_t i = 0; i < sliceCount; i++)
    {
        sliceLists[i] = vector<uint64_t*>(1ULL << (sliceMasks[i].size() * 2));
    }

    uint64_t* offset = allSignatures.data();
    size_t sliceLimitOffset = 0;
    for (size_t i = 0; i < sliceCount; i++) {
        size_t sliceLimit = 1ULL << (sliceMasks[i].size() * 2);
        for (size_t j = 0; j < sliceLimit; j++) {
            size_t idx = sliceLimitOffset + j;
            sliceLists[i][j] = offset;
            offset += allSlicelistSizes[idx];
        }
        sliceLimitOffset += sliceLimit;
    }

    auto endLoading = std::chrono::high_resolution_clock::now();
    auto startProcessing = std::chrono::high_resolution_clock::now();

    //TODO: rewrite
    /** Load query file (candidate guides)
     *      and prepare memory for calculated global scores
     */
    size_t seqLineLength = seqLength + 1;
    std::filesystem::path queryFile(argv[2]);
    size_t fileSize = std::filesystem::file_size(queryFile);
    if (fileSize % seqLineLength != 0) {
        fprintf(stderr, "Error: query file is not a multiple of the expected line length (%zu)\n", seqLineLength);
        fprintf(stderr, "The sequence length may be incorrect; alternatively, the line endings\n");
        fprintf(stderr, "may be something other than LF, or there may be junk at the end of the file.\n");
        exit(1);
    }
    size_t queryCount = fileSize / seqLineLength;
    FILE* fp = fopen(argv[2], "rb");
    vector<char> queryDataSet(fileSize);
    vector<uint64_t> querySignatures(queryCount);
    vector<double> querySignatureMitScores(queryCount);
    vector<double> querySignatureCfdScores(queryCount);

    if (fread(queryDataSet.data(), fileSize, 1, fp) < 1) {
        fprintf(stderr, "Failed to read in query file.\n");
        exit(1);
    }
    fclose(fp);

    /** Binary encode query sequences */
    auto t_seq_start = std::chrono::steady_clock::now();
    #if USE_OPENCL_SEQ2SIG
    {
        SeqSigEncoderCL enc;                      // one-time init (context, queue, program)
        enc.encode(queryDataSet.data(),           // base pointer to the raw bytes
                seqLineLength,                 // stride (including newline)
                (uint32_t)seqLength,           // guide length (20)
                (uint32_t)queryCount,
                querySignatures);              // output vector<uint64_t>
    }
    #else
    #pragma omp parallel
        {
    #pragma omp for
            for (int i = 0; i < queryCount; i++) {
                char* ptr = &queryDataSet[i * seqLineLength];
                uint64_t signature = sequenceToSignature(ptr, 20);
                querySignatures[i] = signature;
            }
        }
    #endif
    auto t_seq_end = std::chrono::steady_clock::now();

    // Reduce into global at thread exit
    g_prof.seq_to_sig_ns = static_cast<uint64_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(t_seq_end - t_seq_start).count());

    if (g_profileLogFile) {
        auto now   = std::chrono::system_clock::now();
        std::time_t now_c = std::chrono::system_clock::to_time_t(now);

        // Thread-safe localtime: localtime_s on Windows, localtime_r elsewhere
        std::tm tm_buf{};
        #ifdef _WIN32
            localtime_s(&tm_buf, &now_c);
        #else
            localtime_r(&now_c, &tm_buf);
        #endif

        // Format "YYYY-MM-DD HH:MM:SS"
        char ts_base[32];
        std::strftime(ts_base, sizeof(ts_base), "%F %T", &tm_buf);

        // Milliseconds
        auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                    now.time_since_epoch()) % 1000;

        // Final "YYYY-MM-DD HH:MM:SS.mmm"
        char ts_full[40];
        std::snprintf(ts_full, sizeof(ts_full), "%s.%03lld",
                    ts_base, static_cast<long long>(ms.count()));

        // Write to file and flush
        std::fprintf(g_profileLogFile, "ENTERING_SCORING_SECTION, %s\n", ts_full);
        std::fflush(g_profileLogFile);

        // Echo to terminal
        std::cout << "ENTERING_SCORING_SECTION at " << ts_full << std::endl;
    }

    /** Begin scoring */
    auto t_score_start = std::chrono::steady_clock::now();
#pragma omp parallel
    {
        ProfStats local{};
        std::minstd_rand rng((unsigned)omp_get_thread_num() + 1234567u);

        vector<uint64_t> offtargetToggles(numOfftargetToggles);
        uint64_t* offtargetTogglesTail = offtargetToggles.data() + numOfftargetToggles - 1;

        auto t_score_thread_start = std::chrono::steady_clock::now();

        #pragma omp for
        for (int searchIdx = 0; searchIdx < querySignatures.size(); searchIdx++) {
            // -------- sampling decision for per-loop timings ----------
            bool do_sample =
            #if MANPROF_ENABLE_LOOP_SAMPLING
                        (searchIdx % MANPROF_SAMPLE_PERIOD == 0);
            #else
                        false;
            #endif
        // ----------------------------------------------------------
            auto t_outer_start = do_sample ? std::chrono::steady_clock::now() : std::chrono::steady_clock::time_point{};

            auto searchSignature = querySignatures[searchIdx];

            /** Global scores */
            double totScoreMit = 0.0;
            double totScoreCfd = 0.0;

            double maximum_sum = (10000.0 - threshold * 100) / threshold;
            bool checkNextSlice = true;

            size_t sliceLimitOffset = 0;

            // ---------- per-guide work counters ----------
            uint64_t guide_signatures_seen = 0;
            uint64_t guide_offtargets_scored = 0;
            uint64_t guide_slices_traversed = 0;
        // ---------------------------------------------

            /** For each ISSL slice */
            for (size_t i = 0; i < sliceCount; i++) {
                auto t_slice_start = do_sample ? std::chrono::steady_clock::now() : std::chrono::steady_clock::time_point{};
                vector<uint64_t>& sliceMask = sliceMasks[i];
                auto& sliceList = sliceLists[i];

                uint64_t searchSlice = 0ULL;
                for (int j = 0; j < sliceMask.size(); j++)
                {
                    searchSlice |= ((searchSignature >> (sliceMask[j] * 2)) & 3ULL) << (j * 2);
                }

                size_t idx = sliceLimitOffset + searchSlice;

                size_t signaturesInSlice = allSlicelistSizes[idx];
                uint64_t* sliceOffset = sliceList[searchSlice];
                guide_slices_traversed++;

                /** For each off-target signature in slice */
                auto t_inner_start = do_sample ? std::chrono::steady_clock::now() : std::chrono::steady_clock::time_point{};

                for (size_t j = 0; j < signaturesInSlice; j++) {
                    guide_signatures_seen++;
                    auto signatureWithOccurrencesAndId = sliceOffset[j];
                    auto signatureId = signatureWithOccurrencesAndId & 0xFFFFFFFFULL;
                    uint32_t occurrences = (signatureWithOccurrencesAndId >> (32));

                    /** Prevent assessing the same off-target for multiple slices */
                    uint64_t seenOfftargetAlready = 0;
                    uint64_t* ptrOfftargetFlag = (offtargetTogglesTail - (signatureId / 64));
                    seenOfftargetAlready = (*ptrOfftargetFlag >> (signatureId % 64)) & 1ULL;

                    if (!seenOfftargetAlready) {
                        *ptrOfftargetFlag |= (1ULL << (signatureId % 64));

                        /** Find the positions of mismatches
                            *
                            *  Search signature (SS):    A  A  T  T    G  C  A  T
                            *                           00 00 11 11   10 01 00 11
                            *
                            *        Off-target (OT):    A  T  A  T    C  G  A  T
                            *                           00 11 00 11   01 10 00 11
                            *
                            *                SS ^ OT:   00 00 11 11   10 01 00 11
                            *                         ^ 00 11 00 11   01 10 00 11
                            *                  (XORd) = 00 11 11 00   11 11 00 00
                            *
                            *        XORd & evenBits:   00 11 11 00   11 11 00 00
                            *                         & 10 10 10 10   10 10 10 10
                            *                   (eX)  = 00 10 10 00   10 10 00 00
                            *
                            *         XORd & oddBits:   00 11 11 00   11 11 00 00
                            *                         & 01 01 01 01   01 01 01 01
                            *                   (oX)  = 00 01 01 00   01 01 00 00
                            *
                            *         (eX >> 1) | oX:   00 01 01 00   01 01 00 00 (>>1)
                            *                         | 00 01 01 00   01 01 00 00
                            *            mismatches   = 00 01 01 00   01 01 00 00
                            *
                            *   popcount(mismatches):   4
                            */
                        uint64_t xoredSignatures = searchSignature ^ offtargets[signatureId];
                        uint64_t evenBits = xoredSignatures & 0xAAAAAAAAAAAAAAAAULL;
                        uint64_t oddBits = xoredSignatures & 0x5555555555555555ULL;
                        uint64_t mismatches = (evenBits >> 1) | oddBits;
                        uint64_t dist = popcount64(mismatches);

                        if (dist >= 0 && dist <= 4) {
                            guide_offtargets_scored++; 
                            // Begin calculating MIT score
                            if (calcMit) {
                                if (dist > 0) {
                                    totScoreMit += precalculatedMITScores.at(mismatches) * (double)occurrences;
                                }
                            }

                            // Begin calculating CFD score
                            if (calcCfd) {
                                /** "In other words, for the CFD score, a value of 0
                                    *      indicates no predicted off-target activity whereas
                                    *      a value of 1 indicates a perfect match"
                                    *      John Doench, 2016.
                                    *      https://www.nature.com/articles/nbt.3437
                                */
                                double cfdScore = 0;
                                if (dist == 0) {
                                    cfdScore = 1;
                                }
                                else {
                                    cfdScore = cfdPamPenalties[0b1010]; // PAM: NGG, TODO: do not hard-code the PAM

                                    for (size_t pos = 0; pos < 20; pos++) {
                                        size_t mask = pos << 4;

                                        /** Create the mask to look up the position-identity score
                                            *      In Python... c2b is char to bit
                                            *       mask = pos << 4
                                            *       mask |= c2b[sgRNA[pos]] << 2
                                            *       mask |= c2b[revcom(offTaret[pos])]
                                            *
                                            *      Find identity at `pos` for search signature
                                            *      example: find identity in pos=2
                                            *       Recall ISSL is inverted, hence:
                                            *                   3'-  T  G  C  C  G  A -5'
                                            *       start           11 10 01 01 10 00
                                            *       3UL << pos*2    00 00 00 11 00 00
                                            *       and             00 00 00 01 00 00
                                            *       shift           00 00 00 00 01 00
                                            */
                                        uint64_t searchSigIdentityPos = searchSignature;
                                        searchSigIdentityPos &= (3ULL << (pos * 2));
                                        searchSigIdentityPos = searchSigIdentityPos >> (pos * 2);
                                        searchSigIdentityPos = searchSigIdentityPos << 2;

                                        /** Find identity at `pos` for offtarget
                                            *      Example: find identity in pos=2
                                            *      Recall ISSL is inverted, hence:
                                            *                  3'-  T  G  C  C  G  A -5'
                                            *      start           11 10 01 01 10 00
                                            *      3UL<<pos*2      00 00 00 11 00 00
                                            *      and             00 00 00 01 00 00
                                            *      shift           00 00 00 00 00 01
                                            *      rev comp 3UL    00 00 00 00 00 10 (done below)
                                            */
                                        uint64_t offtargetIdentityPos = offtargets[signatureId];
                                        offtargetIdentityPos &= (3ULL << (pos * 2));
                                        offtargetIdentityPos = offtargetIdentityPos >> (pos * 2);

                                        /** Complete the mask
                                            *      reverse complement (^3UL) `offtargetIdentityPos` here
                                            */
                                        mask = (mask | searchSigIdentityPos | (offtargetIdentityPos ^ 3UL));

                                        if (searchSigIdentityPos >> 2 != offtargetIdentityPos) {
                                            cfdScore *= cfdPosPenalties[mask];
                                        }
                                    }
                                }
                                totScoreCfd += cfdScore * (double)occurrences;
                            }

                            /** Stop calculating global score early if possible */
                            if (scoreMethod == otScoreMethod::mitAndCfd) {
                                if (totScoreMit > maximum_sum && totScoreCfd > maximum_sum) {
                                    checkNextSlice = false;
                                    break;
                                }
                            }
                            if (scoreMethod == otScoreMethod::mitOrCfd) {
                                if (totScoreMit > maximum_sum || totScoreCfd > maximum_sum) {
                                    checkNextSlice = false;
                                    break;
                                }
                            }
                            if (scoreMethod == otScoreMethod::avgMitCfd) {
                                if (((totScoreMit + totScoreCfd) / 2.0) > maximum_sum) {
                                    checkNextSlice = false;
                                    break;
                                }
                            }
                            if (scoreMethod == otScoreMethod::mit) {
                                if (totScoreMit > maximum_sum) {
                                    checkNextSlice = false;
                                    break;
                                }
                            }
                            if (scoreMethod == otScoreMethod::cfd) {
                                if (totScoreCfd > maximum_sum) {
                                    checkNextSlice = false;
                                    break;
                                }
                            }
                        }
                    }
                    
                }
                if (do_sample) {
                    auto t_inner_end = std::chrono::steady_clock::now();
                    local.sampled_inner_loop_ns += (uint64_t)std::chrono::duration_cast<std::chrono::nanoseconds>(t_inner_end - t_inner_start).count();
                    local.samples_inner++;
                }

                if (do_sample) {
                    auto t_slice_end = std::chrono::steady_clock::now();
                    local.sampled_slice_loop_ns += (uint64_t)std::chrono::duration_cast<std::chrono::nanoseconds>(t_slice_end - t_slice_start).count();
                    local.samples_slice++;
                }

                if (!checkNextSlice)
                    break;
                sliceLimitOffset += 1ULL << (sliceMasks[i].size() * 2);
            }
            querySignatureMitScores[searchIdx] = 10000.0 / (100.0 + totScoreMit);
            querySignatureCfdScores[searchIdx] = 10000.0 / (100.0 + totScoreCfd);

            memset(offtargetToggles.data(), 0, sizeof(uint64_t) * offtargetToggles.size());

            if (do_sample) {
                auto t_outer_end = std::chrono::steady_clock::now();
                local.sampled_outer_iter_ns += (uint64_t)std::chrono::duration_cast<std::chrono::nanoseconds>(t_outer_end - t_outer_start).count();
                local.samples_outer++;
            }

            local.guides_scored++;
            local.signatures_seen += guide_signatures_seen;
            local.offtargets_scored += guide_offtargets_scored;
            local.slices_traversed += guide_slices_traversed;
        }

        #pragma omp critical
        {
            g_prof.seq_to_sig_ns      += 0; // already added earlier
            g_prof.scoring_total_ns   += 0;
            g_prof.guides_scored      += local.guides_scored;
            g_prof.signatures_seen    += local.signatures_seen;
            g_prof.offtargets_scored  += local.offtargets_scored;
            g_prof.slices_traversed   += local.slices_traversed;

            g_prof.sampled_outer_iter_ns += local.sampled_outer_iter_ns;
            g_prof.sampled_inner_loop_ns += local.sampled_inner_loop_ns;
            g_prof.sampled_slice_loop_ns += local.sampled_slice_loop_ns;
            g_prof.samples_outer += local.samples_outer;
            g_prof.samples_inner += local.samples_inner;
            g_prof.samples_slice += local.samples_slice;
        }
    }
    auto t_score_end = std::chrono::steady_clock::now();
    g_prof.scoring_total_ns = static_cast<uint64_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(t_score_end - t_score_start).count());

    auto progEnd = std::chrono::steady_clock::now();
    auto elapsedUs = std::chrono::duration_cast<std::chrono::microseconds>(progEnd - g_progStart).count();
    double elapsedSec = (double)elapsedUs / 1e6;

    auto ns_to_s = [](uint64_t ns){ return (double)ns / 1e9; };
    auto avg_ms  = [](uint64_t total_ns, uint64_t n){
        return (n == 0) ? 0.0 : ( (double)total_ns / (double)n ) / 1e6;
    };

    if (g_profileLogFile) {
        auto now   = std::chrono::system_clock::now();
        std::time_t now_c = std::chrono::system_clock::to_time_t(now);

        // Thread-safe localtime: localtime_s on Windows, localtime_r elsewhere
        std::tm tm_buf{};
        #ifdef _WIN32
            localtime_s(&tm_buf, &now_c);
        #else
            localtime_r(&now_c, &tm_buf);
        #endif

        // Format "YYYY-MM-DD HH:MM:SS"
        char ts_base[32];
        std::strftime(ts_base, sizeof(ts_base), "%F %T", &tm_buf);

        // Milliseconds
        auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                    now.time_since_epoch()) % 1000;

        // Final "YYYY-MM-DD HH:MM:SS.mmm"
        char ts_full[40];
        std::snprintf(ts_full, sizeof(ts_full), "%s.%03lld",
                    ts_base, static_cast<long long>(ms.count()));

        // Write to file and flush
        std::fprintf(g_profileLogFile, "FINISHED SCORING, %s\n", ts_full);
        std::fflush(g_profileLogFile);

        // Echo to terminal
        std::cout << "FINISHED SCORING at " << ts_full << std::endl;
    }

    if (g_profileLogFile) {
        // rewrite the header lines with real values
        // (for simplicity, append rows again with real data)
        std::fprintf(g_profileLogFile, "TOTAL_EXECUTION_TIME_S, %.6f,\n", elapsedSec);
        std::fprintf(g_profileLogFile, "SEQ_TO_SIG_TOTAL_S, %.6f,\n", ns_to_s(g_prof.seq_to_sig_ns));
        std::fprintf(g_profileLogFile, "SCORING_TOTAL_S, %.6f,\n", ns_to_s(g_prof.scoring_total_ns));
        std::fprintf(g_profileLogFile, "GUIDES_SCORED,% " PRIu64 ",\n", g_prof.guides_scored);
        std::fprintf(g_profileLogFile, "SIGNATURES_SEEN,% " PRIu64 ",counts across all slices\n", g_prof.signatures_seen);
        std::fprintf(g_profileLogFile, "OFFTARGETS_SCORED,% " PRIu64 ",distinct dist<=4 actually scored\n", g_prof.offtargets_scored);
        std::fprintf(g_profileLogFile, "SLICES_TRAVERSED,% " PRIu64 ",\n", g_prof.slices_traversed);

    #if MANPROF_ENABLE_LOOP_SAMPLING
        std::fprintf(g_profileLogFile, "SAMPLED_OUTER_ITER_AVG_MS, %.3f, %" PRIu64 " samples\n",
                    avg_ms(g_prof.sampled_outer_iter_ns, g_prof.samples_outer), g_prof.samples_outer);
        std::fprintf(g_profileLogFile, "SAMPLED_INNER_LOOP_AVG_MS, %.3f, %" PRIu64 " samples\n",
                    avg_ms(g_prof.sampled_inner_loop_ns, g_prof.samples_inner), g_prof.samples_inner);
        std::fprintf(g_profileLogFile, "SAMPLED_SLICE_LOOP_AVG_MS, %.3f, %" PRIu64 " samples\n",
                    avg_ms(g_prof.sampled_slice_loop_ns, g_prof.samples_slice), g_prof.samples_slice);
    #endif

        std::fflush(g_profileLogFile);
        std::fclose(g_profileLogFile);
        g_profileLogFile = nullptr;
    }
    // Also print to console
    std::cout << "Total execution time: " << elapsedSec << " seconds" << std::endl;

    /** Print global scores to stdout */
    for (size_t searchIdx = 0; searchIdx < querySignatures.size(); searchIdx++) {
        auto querySequence = signatureToSequence(querySignatures[searchIdx], 20);
        printf("%s\t", querySequence.c_str());
        if (calcMit)
            printf("%f\t", querySignatureMitScores[searchIdx]);
        else
            printf("-1\t");

        if (calcCfd)
            printf("%f\n", querySignatureCfdScores[searchIdx]);
        else
            printf("-1\n");

    }
    return 0;
}