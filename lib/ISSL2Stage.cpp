#include "../include/ISSL2Stage.hpp"
// TODO: remove
#include <map>
#include <chrono>
#include <fstream>
#include <map>
using std::chrono::time_point;
using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::nanoseconds;
using std::chrono::seconds;

#if defined(_WIN64)
#pragma push_macro("close")
#undef close
#endif

using std::string;
using std::vector;
using std::pair;
using std::unordered_set;
using std::unordered_map;

/** Char to binary encoding */
const vector<uint8_t> nucleotideIndex{ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,3 };
const vector<char> signatureIndex{ 'A', 'C', 'G', 'T' };
enum ScoreMethod { unknown = 0, mit = 1, cfd = 2, mitAndCfd = 3, mitOrCfd = 4, avgMitCfd = 5 };
static std::map<uint64_t, uint64_t> truePerGuideCountTotal;
static std::map<uint64_t, uint64_t> truePerGuideCountUnique;
static std::map<uint64_t, uint64_t> falsePerGuideCountTotal;
static std::map<uint64_t, uint64_t> falsePerGuideCountUnique;

void countTrueTotal(std::filesystem::path& outputPath, vector<vector<uint64_t>>& perGuideCountTotal)
{
    std::ofstream outputFile;

    for (const uint64_t& count : perGuideCountTotal[true])
    {
        truePerGuideCountTotal[count]++;
    }

    outputFile.open(outputPath / "_output" / "truePerGuideCountTotal.txt", std::ios::out | std::ios::binary);
    for (auto const& x : truePerGuideCountTotal)
    {
        outputFile << fmt::format("{}:{}\n", x.first, x.second);
    }
    outputFile.close();
}

void countTrueUnique(std::filesystem::path& outputPath, vector<vector<uint64_t>>& perGuideCountUnique)
{
    std::ofstream outputFile;

    for (const int64_t& count : perGuideCountUnique[true])
    {
        truePerGuideCountUnique[count]++;
    }

    outputFile.open(outputPath / "_output" / "truePerGuideCountUnique.txt", std::ios::out | std::ios::binary);
    for (auto const& x : truePerGuideCountUnique)
    {
        outputFile << fmt::format("{}:{}\n", x.first, x.second);
    }
    outputFile.close();
}

void countFalseTotal(std::filesystem::path& outputPath, vector<vector<uint64_t>>& perGuideCountTotal)
{
    std::ofstream outputFile;

    for (int64_t count : perGuideCountTotal[false])
    {
        falsePerGuideCountTotal[count]++;
    }

    outputFile.open(outputPath / "_output" / "falsePerGuideCountTotal.txt", std::ios::out | std::ios::binary);
    for (auto const& x : falsePerGuideCountTotal)
    {
        outputFile << fmt::format("{}:{}\n", x.first, x.second);
    }
    outputFile.close();
}

void countFalseUnique(std::filesystem::path& outputPath, vector<vector<uint64_t>>& perGuideCountUnique)
{
    std::ofstream outputFile;

    for (int64_t count : perGuideCountUnique[false])
    {
        falsePerGuideCountUnique[count]++;
    }

    outputFile.open(outputPath / "_output" / "falsePerGuideCountUnique.txt", std::ios::out | std::ios::binary);
    for (auto const& x : falsePerGuideCountUnique)
    {
        outputFile << fmt::format("{}:{}\n", x.first, x.second);
    }
    outputFile.close();
}

ISSL2Stage::ISSL2Stage(ConfigManager& cm) :
    toolIsSelected(cm.getBool("offtargetscore", "enabled")),
    optimsationLevel(cm.getString("general", "optimisation")),
    toolCount(cm.getConsensusToolCount()),
    consensusN(cm.getInt("consensus", "n")),
    threadCount(cm.getInt("offtargetscore", "threads")),
    ISSLIndex(cm.getString("input", "offtarget-sites")),
    maxDist(cm.getInt("offtargetscore", "max-distance")),
    scoreMethod(cm.getString("offtargetscore", "method")),
    scoreThreshold(cm.getFloat("offtargetscore", "score-threshold")),
    pageLength(cm.getInt("offtargetscore", "page-length"))
{}

void ISSL2Stage::run(unordered_map<string, unordered_map<string, string>>& candidateGuides)
{
    if (!toolIsSelected)
    {
        printer("Off-target scoring has been configured not to run. Skipping Off-target scoring");
        return;
    }

    printer("Loading ISSL Index.");

    /** Scoring methods. To exit early:
     *      - only CFD must drop below `threshold`
     *      - only MIT must drop below `threshold`
     *      - both CFD and MIT must drop below `threshold`
     *      - CFD or MIT must drop below `threshold`
     *      - the average of CFD and MIT must below `threshold`
     */
    ScoreMethod sm = ScoreMethod::unknown;
    bool calcCfd = false;
    bool calcMit = false;
    if (scoreMethod == "and") {
        sm = ScoreMethod::mitAndCfd;
        calcCfd = true;
        calcMit = true;
    }
    else if (scoreMethod == "or") {
        sm = ScoreMethod::mitOrCfd;
        calcCfd = true;
        calcMit = true;
    }
    else if (scoreMethod == "avg") {
        sm = ScoreMethod::avgMitCfd;
        calcCfd = true;
        calcMit = true;
    }
    else if (scoreMethod == "mit") {
        sm = ScoreMethod::mit;
        calcMit = true;
    }
    else if (scoreMethod == "cfd") {
        sm = ScoreMethod::cfd;
        calcCfd = true;
    }

    /** Begin reading the binary encoded ISSL, structured as:
     *      - a header (5 items)
     *      - precalcuated local MIT scores
     *      - length of all the slices
     *      - the positions within a slice
     *      - all binary-encoded off-target sites
     */
    FILE* isslFp = fopen(fmt::format("{}-1", ISSLIndex).c_str(), "rb");



    /** The index contains a fixed-sized header
     *      - the number of unique off-targets in the index
     *      - the length of an off-target
     *      - the total number of off-targets
     *      - the number of slices
     *      - the number of precalculated MIT scores
     */
    vector<size_t> slicelistHeader(50);

    if (fread(slicelistHeader.data(), sizeof(size_t), slicelistHeader.size(), isslFp) == 0) {
        throw std::runtime_error("Error reading index: header invalid\n");
    }

    size_t offtargetsCount = slicelistHeader[0];
    size_t seqLength = slicelistHeader[1];
    size_t seqCount = slicelistHeader[2];
    size_t sliceCount = slicelistHeader[3];
    size_t scoresCount = slicelistHeader[4];

    /** Read in the precalculated MIT scores
     *      - `mask` is a 2-bit encoding of mismatch positions
     *          For example,
     *              00 01 01 00 01  indicates mismatches in positions 1, 3 and 4
     *
     *      - `score` is the local MIT score for this mismatch combination
     */
    phmap::flat_hash_map<uint64_t, double> precalculatedScores;
    for (int i = 0; i < scoresCount; i++) {
        uint64_t mask = 0;
        double score = 0.0;
        fread(&mask, sizeof(uint64_t), 1, isslFp);
        fread(&score, sizeof(double), 1, isslFp);

        precalculatedScores.insert(pair<uint64_t, double>(mask, score));
    }


    /**
    * Read the slice lengths from header
    */
    vector<size_t> sliceLens;
    for (int i = 0; i < sliceCount; i++)
    {
        size_t sliceLen;
        fread(&sliceLen, sizeof(size_t), 1, isslFp);
        sliceLens.push_back(sliceLen);
    }

    vector<vector<int>> sliceMasks;
    for (int i = 0; i < sliceCount; i++) {
        vector<int> mask;
        for (int j = 0; j < sliceLens[i]; j++)
        {
            int pos;
            fread(&pos, sizeof(int), 1, isslFp);
            mask.push_back(pos);
        }
        sliceMasks.push_back(mask);
    }

    /** Load in all of the off-target sites */
    vector<uint64_t> offtargets(offtargetsCount);
    if (fread(offtargets.data(), sizeof(uint64_t), offtargetsCount, isslFp) == 0) {
        throw std::runtime_error("Error reading index: loading off-target sequences failed\n");
    }

    /** The number of signatures embedded per slice
     *
     *      These counts are stored contiguously
     *
     */
    size_t sliceListCount = 0;
    for (int i = 0; i < sliceCount; i++)
    {
        sliceListCount += 1ULL << (sliceLens[i] * 2);
    }
    vector<size_t> allSlicelistSizes(sliceListCount);

    if (fread(allSlicelistSizes.data(), sizeof(size_t), allSlicelistSizes.size(), isslFp) == 0) {
        throw std::runtime_error("Error reading index: reading slice list sizes failed\n");
    }

    /** The contents of the slices
     *
     *      Stored contiguously
     *
     *      Each signature (64-bit) is structured as:
     *          <occurrences 32-bit><off-target-id 32-bit>
     */
    vector<uint64_t> allSignatures(offtargetsCount * sliceCount);

    if (fread(allSignatures.data(), sizeof(uint64_t), allSignatures.size(), isslFp) == 0) {
        throw std::runtime_error("Error reading index: reading slice contents failed\n");
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
    for (int i = 0; i < sliceCount; i++)
    {
        sliceLists[i] = vector<uint64_t*>(1ULL << (sliceLens[i] * 2));
    }

    uint64_t* offset = allSignatures.data();
    size_t sliceLimitOffset = 0;
    for (size_t i = 0; i < sliceCount; i++) {
        size_t sliceLimit = 1ULL << (sliceLens[i] * 2);
        for (size_t j = 0; j < sliceLimit; j++) {
            size_t idx = sliceLimitOffset + j;
            sliceLists[i][j] = offset;
            offset += allSlicelistSizes[idx];
        }
        sliceLimitOffset += sliceLimit;
    }

    printer("Beginning Off-target scoring.");
    int testedCount = 0;
    int failedCount = 0;
    int pgIdx = 1;
    int guidesInPage = 0;
    auto paginatorIterator = candidateGuides.begin();
    auto pageStart = candidateGuides.begin();
    auto pageEnd = candidateGuides.begin();

    // Outer loop deals with changing iterator start and end points (Pagination)
    while (pageEnd != candidateGuides.end())
    {
        if (pageLength > 0)
        {
            // Advance the pageEnd pointer
            std::advance(pageEnd, (std::min)((int)std::distance(pageEnd, candidateGuides.end()), pageLength));
            // Record page start
            pageStart = paginatorIterator;
            // Print page information
            printer(fmt::format("\tProcessing page {} ({} per page).", commaify(pgIdx), commaify(pageLength)));
        }
        else {
            // Process all guides at once
            pageEnd = candidateGuides.end();
        }


        printer("\tConstructing the Off-target scoring input.");


        // New intergrated method (No need for external file, simply write to char vector)
        vector<char> queryDataSet;
        queryDataSet.reserve(pageLength * seqLength);
        vector<char> pamDataSet;
        pamDataSet.reserve(pageLength * 3);

        guidesInPage = 0;
        while (paginatorIterator != pageEnd)
        {
            string target23 = paginatorIterator->first;
            // Run time filtering
            if (!filterCandidateGuides(paginatorIterator->second, MODULE_SPECIFICITY, optimsationLevel, consensusN, toolCount)) {
                // Advance page end for each filtered out guide
                if (pageEnd != candidateGuides.end())
                {
                    pageEnd++;
                }
                paginatorIterator++;
                continue;
            }

            for (char c : target23.substr(0, 20))
            {
                queryDataSet.push_back(c);
            }

            for (char c : target23.substr(20))
            {
                pamDataSet.push_back(c);
            }
            guidesInPage++;
            paginatorIterator++;
        }



        printer(fmt::format("\t\t{} guides in this page.", commaify(guidesInPage)));

        // Call scoring method
        size_t seqLineLength = seqLength;

        size_t queryCount = guidesInPage;
        vector<uint64_t> querySignatures(queryCount);
        vector<double> querySignatureMitScores(queryCount);
        vector<double> querySignatureCfdScores(queryCount);

        /** Binary encode query sequences */
        omp_set_num_threads(threadCount);
#pragma omp parallel
        {
#pragma omp for
            for (int i = 0; i < queryCount; i++) {
                char* ptr = &queryDataSet[i * seqLineLength];
                uint64_t signature = sequenceToSignature(ptr, seqLength);
                querySignatures[i] = signature;
            }
        }

        // TODO: remove
        // Mutex for thread safety
        std::mutex countsMutex;
        // Timing
        std::chrono::duration<double> usefulTime;
        std::chrono::duration<double> wastedTime;
        // Neighbourhood counts
        vector<uint64_t> neighbourhoodCountTotal(sliceCount);
        vector<uint64_t> neighbourhoodCountUnique(sliceCount);
        // Mismatch counts 
        vector<vector<uint64_t>> offTargetCountTotal(sliceCount, vector<uint64_t>(21));
        vector<vector<uint64_t>> offTargetCountUnique(sliceCount, vector<uint64_t>(21));
        // Per guide counts 
        vector<vector<uint64_t>> perGuideCountTotal(2, vector<uint64_t>(querySignatures.size()));
        vector<vector<uint64_t>> perGuideCountUnique(2, vector<uint64_t>(querySignatures.size()));


        /** Begin scoring */
        omp_set_num_threads(threadCount);
        #pragma omp parallel
        {
            unordered_map<uint64_t, unordered_set<uint64_t>> searchResults;
            vector<uint64_t> offtargetToggles(numOfftargetToggles);

            uint64_t* offtargetTogglesTail = offtargetToggles.data() + numOfftargetToggles - 1;

            // TODO: remove
            // Timing
            std::chrono::duration<double> usefulTimeLocal;
            std::chrono::duration<double> wastedTimeLocal;
            // Neighbourhood counts
            vector<uint64_t> neighbourhoodCountTotalLocal(sliceCount);
            vector<uint64_t> neighbourhoodCountUniqueLocal(sliceCount);
            // Mismatch counts 
            vector<vector<uint64_t>> offTargetCountTotalLocal(sliceCount, vector<uint64_t>(21));
            vector<vector<uint64_t>> offTargetCountUniqueLocal(sliceCount, vector<uint64_t>(21));
            // Per guide counts 
            vector<vector<uint64_t>> perGuideCountTotalLocal(2, vector<uint64_t>(querySignatures.size()));
            vector<vector<uint64_t>> perGuideCountUniqueLocal(2, vector<uint64_t>(querySignatures.size()));

            /** For each candidate guide */
#pragma omp for
            for (int searchIdx = 0; searchIdx < querySignatures.size(); searchIdx++) {

                auto searchSignature = querySignatures[searchIdx];

                /** Global scores */
                double totScoreMit = 0.0;
                double totScoreCfd = 0.0;

                int numOffTargetSitesScored = 0;
                double maximum_sum = (10000.0 - scoreThreshold * 100) / scoreThreshold;
                // TODO: uncommment
                //bool checkNextSlice = true;

                size_t sliceLimitOffset = 0;
                /** For each ISSL slice */
                for (size_t i = 0; i < sliceCount; i++) {
                    vector<int>& sliceMask = sliceMasks[i];
                    auto& sliceList = sliceLists[i];

                    uint64_t searchSlice = 0ULL;
                    for (int j = 0; j < sliceMask.size(); j++)
                    {
                        searchSlice |= ((searchSignature >> (sliceMask[j] * 2)) & 3ULL) << (j * 2);
                    }

                    size_t idx = sliceLimitOffset + searchSlice;

                    size_t signaturesInSlice = allSlicelistSizes[idx];
                    uint64_t* sliceOffset = sliceList[searchSlice];

                    /** For each off-target signature in slice */
                    for (size_t j = 0; j < signaturesInSlice; j++) {
                        // TODO: remove
                        time_point start = high_resolution_clock::now();

                        auto signatureWithOccurrencesAndId = sliceOffset[j];
                        auto signatureId = signatureWithOccurrencesAndId & 0xFFFFFFFFULL;
                        uint32_t occurrences = (signatureWithOccurrencesAndId >> (32));

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
                        /*int dist = popcnt(&mismatches, sizeof(uint64_t));*/
                        int dist = popcount64(mismatches);

                        if (dist >= 0 && dist <= maxDist) {

                            /** Prevent assessing the same off-target for multiple slices */
                            uint64_t seenOfftargetAlready = 0;
                            uint64_t* ptrOfftargetFlag = (offtargetTogglesTail - (signatureId / 64));
                            seenOfftargetAlready = (*ptrOfftargetFlag >> (signatureId % 64)) & 1ULL;


                            if (!seenOfftargetAlready) {
                                // Begin calculating MIT score
                                if (calcMit) {
                                    if (dist > 0) {
                                        totScoreMit += precalculatedScores[mismatches] * (double)occurrences;
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
                                    else if (dist > 0 && dist <= maxDist) {
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
                                            searchSigIdentityPos &= (3UL << (pos * 2));
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
                                            offtargetIdentityPos &= (3UL << (pos * 2));
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

                                *ptrOfftargetFlag |= (1ULL << (signatureId % 64));
                                numOffTargetSitesScored += occurrences;

                                // TODO: uncomment
                                /** Stop calculating global score early if possible */
                                //if (sm == ScoreMethod::mitAndCfd) {
                                //    if (totScoreMit > maximum_sum && totScoreCfd > maximum_sum) {
                                //        checkNextSlice = false;
                                //        break;
                                //    }
                                //}
                                //if (sm == ScoreMethod::mitOrCfd) {
                                //    if (totScoreMit > maximum_sum || totScoreCfd > maximum_sum) {
                                //        checkNextSlice = false;
                                //        break;
                                //    }
                                //}
                                //if (sm == ScoreMethod::avgMitCfd) {
                                //    if (((totScoreMit + totScoreCfd) / 2.0) > maximum_sum) {
                                //        checkNextSlice = false;
                                //        break;
                                //    }
                                //}
                                //if (sm == ScoreMethod::mit) {
                                //    if (totScoreMit > maximum_sum) {
                                //        checkNextSlice = false;
                                //        break;
                                //    }
                                //}
                                //if (sm == ScoreMethod::cfd) {
                                //    if (totScoreCfd > maximum_sum) {
                                //        checkNextSlice = false;
                                //        break;
                                //    }
                                //}
                            }
                            // TODO: remove
                            // Updated counts and timings
                            std::chrono::duration<double> computeTime = high_resolution_clock::now() - start;
                            // Timing
                            usefulTimeLocal += computeTime;
                            // Neighbourhood counts
                            neighbourhoodCountTotalLocal[i] += occurrences;
                            neighbourhoodCountUniqueLocal[i]++;
                            // Mismatch counts 
                            offTargetCountTotalLocal[i][dist] += occurrences;
                            offTargetCountUniqueLocal[i][dist]++;
                            // Per guide counts
                            perGuideCountTotalLocal[true][searchIdx] += occurrences;
                            perGuideCountUniqueLocal[true][searchIdx]++;
                        }
                        else
                        {
                            // TODO: remove
                            // Updated counts and timings
                            std::chrono::duration<double> computeTime = high_resolution_clock::now() - start;
                            // Timing
                            wastedTimeLocal += computeTime;
                            // Neighbourhood counts
                            neighbourhoodCountTotalLocal[i] += occurrences;
                            neighbourhoodCountUniqueLocal[i]++;
                            // Mismatch counts 
                            offTargetCountTotalLocal[i][dist] += occurrences;
                            offTargetCountUniqueLocal[i][dist]++;
                            // Per guide counts
                            perGuideCountTotalLocal[false][searchIdx] += occurrences;
                            perGuideCountUniqueLocal[false][searchIdx]++;
                        }
                    }

                    // TODO: uncomment
                    //if (!checkNextSlice)
                    //    break;
                    sliceLimitOffset += 1ULL << (sliceLens[i] * 2);
                }
                querySignatureMitScores[searchIdx] = 10000.0 / (100.0 + totScoreMit);
                querySignatureCfdScores[searchIdx] = 10000.0 / (100.0 + totScoreCfd);

                memset(offtargetToggles.data(), 0, sizeof(uint64_t) * offtargetToggles.size());
            }
            countsMutex.lock();
            usefulTime += usefulTimeLocal;
            wastedTime += wastedTimeLocal;
            for (size_t i = 0; i < neighbourhoodCountTotalLocal.size(); i++)
            {
                neighbourhoodCountTotal[i] += neighbourhoodCountTotalLocal[i];
                neighbourhoodCountUnique[i] += neighbourhoodCountUniqueLocal[i];
            }

            for (size_t i = 0; i < offTargetCountTotalLocal.size(); i++)
            {
                for (size_t j = 0; j < offTargetCountTotalLocal[i].size(); j++)
                {
                    offTargetCountTotal[i][j] += offTargetCountTotalLocal[i][j];
                    offTargetCountUnique[i][j] += offTargetCountUniqueLocal[i][j];

                }
            }

            for (size_t i = 0; i < perGuideCountTotalLocal.size(); i++)
            {
                for (size_t j = 0; j < perGuideCountTotalLocal[i].size(); j++)
                {
                    perGuideCountTotal[i][j] += perGuideCountTotalLocal[i][j];
                    perGuideCountUnique[i][j] += perGuideCountUniqueLocal[i][j];
                }

            }
            countsMutex.unlock();
        }

        printer("\tStarting to process the Off-target scoring results.");

        // TODO: remove
        printer(fmt::format("Wasted time: {}\tUseful time: {}", wastedTime.count(), usefulTime.count()));

        for (int i = 0; i < neighbourhoodCountTotal.size(); i++)
        {
            printer(fmt::format("Slice {}\tTotal {}\tUnique {}", i + 1, commaify(neighbourhoodCountTotal[i]), commaify(neighbourhoodCountUnique[i])));
            for (int j = 0; j < 21; j++)
            {
                printer(fmt::format("\tMismatch {}\tTotal {}\tUnique {}", j, commaify(offTargetCountTotal[i][j]), commaify(offTargetCountUnique[i][j])));
            }
        }

        std::filesystem::path outputPath(std::filesystem::path(ISSLIndex).parent_path());

        std::thread countTrueTotalThread(countTrueTotal, outputPath, perGuideCountTotal);
        std::thread countTrueUniqueThread(countTrueUnique, outputPath, perGuideCountUnique);
        std::thread countFalseTotalThread(countFalseTotal, outputPath, perGuideCountTotal);
        std::thread countFalseUniqueThread(countFalseUnique, outputPath, perGuideCountUnique);

        countTrueTotalThread.join();
        countTrueUniqueThread.join();
        countFalseTotalThread.join();
        countFalseUniqueThread.join();

        for (size_t searchIdx = 0; searchIdx < querySignatures.size(); searchIdx++) {
            string target20 = signatureToSequence(querySignatures[searchIdx], seqLength);
            char* pamPtr = &pamDataSet[searchIdx * 3];
            string pam(pamPtr, pamPtr + 3);
            string target23 = target20 + pam;
            candidateGuides[target23]["mitOfftargetscore"] = calcMit ? std::to_string(querySignatureMitScores[searchIdx]) : "-1";
            candidateGuides[target23]["cfdOfftargetscore"] = calcCfd ? std::to_string(querySignatureCfdScores[searchIdx]) : "-1";

            if (sm == ScoreMethod::mit)
            {
                if (std::stof(candidateGuides[target23]["mitOfftargetscore"]) < scoreThreshold)
                {
                    candidateGuides[target23]["passedOffTargetScore"] = CODE_REJECTED;
                    failedCount++;
                }
                else { candidateGuides[target23]["passedOffTargetScore"] = CODE_ACCEPTED; }
            }

            else if (sm == ScoreMethod::cfd)
            {
                if (std::stof(candidateGuides[target23]["cfdOfftargetscore"]) < scoreThreshold)
                {
                    candidateGuides[target23]["passedOffTargetScore"] = CODE_REJECTED;
                    failedCount++;
                }
                else { candidateGuides[target23]["passedOffTargetScore"] = CODE_ACCEPTED; }
            }

            else if (sm == ScoreMethod::mitAndCfd)
            {
                if ((std::stof(candidateGuides[target23]["mitOfftargetscore"]) < scoreThreshold) && (std::stof(candidateGuides[target23]["cfdOfftargetscore"]) < scoreThreshold))
                {
                    candidateGuides[target23]["passedOffTargetScore"] = CODE_REJECTED;
                    failedCount++;
                }
                else { candidateGuides[target23]["passedOffTargetScore"] = CODE_ACCEPTED; }
            }

            else if (sm == ScoreMethod::mitOrCfd)
            {
                if ((std::stof(candidateGuides[target23]["mitOfftargetscore"]) < scoreThreshold) || (std::stof(candidateGuides[target23]["cfdOfftargetscore"]) < scoreThreshold))
                {
                    candidateGuides[target23]["passedOffTargetScore"] = CODE_REJECTED;
                    failedCount++;
                }
                else { candidateGuides[target23]["passedOffTargetScore"] = CODE_ACCEPTED; }
            }

            else if (sm == ScoreMethod::avgMitCfd)
            {
                if (((std::stof(candidateGuides[target23]["mitOfftargetscore"]) + std::stof(candidateGuides[target23]["cfdOfftargetscore"])) / 2) < scoreThreshold)
                {
                    candidateGuides[target23]["passedOffTargetScore"] = CODE_REJECTED;
                    failedCount++;
                }
                else { candidateGuides[target23]["passedOffTargetScore"] = CODE_ACCEPTED; }
            }
            testedCount++;
        }

        printer(fmt::format("\t{} of {} failed here.", commaify(failedCount), commaify(testedCount)));
        pageStart = pageEnd;
    }

    printer("Loading ISSL Index.");

    /** Begin reading the binary encoded ISSL, structured as:
     *      - a header (5 items)
     *      - precalcuated local MIT scores
     *      - length of all the slices
     *      - the positions within a slice
     *      - all binary-encoded off-target sites
     */

    printer(fmt::format("{}-2", ISSLIndex).c_str());
    isslFp = fopen(fmt::format("{}-2", ISSLIndex).c_str(), "rb");



    /** The index contains a fixed-sized header
     *      - the number of unique off-targets in the index
     *      - the length of an off-target
     *      - the total number of off-targets
     *      - the number of slices
     *      - the number of precalculated MIT scores
     */
    if (fread(slicelistHeader.data(), sizeof(size_t), slicelistHeader.size(), isslFp) == 0) {
        throw std::runtime_error("Error reading index: header invalid\n");
    }

    offtargetsCount = slicelistHeader[0];
    seqLength = slicelistHeader[1];
    seqCount = slicelistHeader[2];
    sliceCount = slicelistHeader[3];
    scoresCount = slicelistHeader[4];

    /** Read in the precalculated MIT scores
     *      - `mask` is a 2-bit encoding of mismatch positions
     *          For example,
     *              00 01 01 00 01  indicates mismatches in positions 1, 3 and 4
     *
     *      - `score` is the local MIT score for this mismatch combination
     */
    for (int i = 0; i < scoresCount; i++) {
        uint64_t mask = 0;
        double score = 0.0;
        fread(&mask, sizeof(uint64_t), 1, isslFp);
        fread(&score, sizeof(double), 1, isslFp);

        precalculatedScores.insert(pair<uint64_t, double>(mask, score));
    }


    /**
    * Read the slice lengths from header
    */
    sliceLens.clear();
    for (int i = 0; i < sliceCount; i++)
    {
        size_t sliceLen;
        fread(&sliceLen, sizeof(size_t), 1, isslFp);
        sliceLens.push_back(sliceLen);
    }

    sliceMasks.clear();
    for (int i = 0; i < sliceCount; i++) {
        vector<int> mask;
        for (int j = 0; j < sliceLens[i]; j++)
        {
            int pos;
            fread(&pos, sizeof(int), 1, isslFp);
            mask.push_back(pos);
        }
        sliceMasks.push_back(mask);
    }

    /** Load in all of the off-target sites */
    if (fread(offtargets.data(), sizeof(uint64_t), offtargetsCount, isslFp) == 0) {
        throw std::runtime_error("Error reading index: loading off-target sequences failed\n");
    }

    /** The number of signatures embedded per slice
     *
     *      These counts are stored contiguously
     *
     */
    sliceListCount = 0;
    for (int i = 0; i < sliceCount; i++)
    {
        sliceListCount += 1ULL << (sliceLens[i] * 2);
    }
     allSlicelistSizes = vector<size_t>(sliceListCount);

    if (fread(allSlicelistSizes.data(), sizeof(size_t), allSlicelistSizes.size(), isslFp) == 0) {
        throw std::runtime_error("Error reading index: reading slice list sizes failed\n");
    }

    /** The contents of the slices
     *
     *      Stored contiguously
     *
     *      Each signature (64-bit) is structured as:
     *          <occurrences 32-bit><off-target-id 32-bit>
     */
    allSignatures = vector<uint64_t>(offtargetsCount * sliceCount);

    if (fread(allSignatures.data(), sizeof(uint64_t), allSignatures.size(), isslFp) == 0) {
        throw std::runtime_error("Error reading index: reading slice contents failed\n");
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
    numOfftargetToggles = (offtargetsCount / ((size_t)sizeof(uint64_t) * (size_t)CHAR_BIT)) + 1;


    sliceLists = vector<vector<uint64_t*>>(sliceCount);
    // Assign sliceLists size based on each slice length
    for (int i = 0; i < sliceCount; i++)
    {
        sliceLists[i] = vector<uint64_t*>(1ULL << (sliceLens[i] * 2));
    }

    offset = allSignatures.data();
    sliceLimitOffset = 0;
    for (size_t i = 0; i < sliceCount; i++) {
        size_t sliceLimit = 1ULL << (sliceLens[i] * 2);
        for (size_t j = 0; j < sliceLimit; j++) {
            size_t idx = sliceLimitOffset + j;
            sliceLists[i][j] = offset;
            offset += allSlicelistSizes[idx];
        }
        sliceLimitOffset += sliceLimit;
    }

    printer("Beginning Off-target scoring. pass 2");
    testedCount = 0;
    failedCount = 0;
    pgIdx = 1;
    guidesInPage = 0;
    paginatorIterator = candidateGuides.begin();
    pageStart = candidateGuides.begin();
    pageEnd = candidateGuides.begin();

    // Outer loop deals with changing iterator start and end points (Pagination)
    while (pageEnd != candidateGuides.end())
    {
        if (pageLength > 0)
        {
            // Advance the pageEnd pointer
            std::advance(pageEnd, (std::min)((int)std::distance(pageEnd, candidateGuides.end()), pageLength));
            // Record page start
            pageStart = paginatorIterator;
            // Print page information
            printer(fmt::format("\tProcessing page {} ({} per page).", commaify(pgIdx), commaify(pageLength)));
        }
        else {
            // Process all guides at once
            pageEnd = candidateGuides.end();
        }


        printer("\tConstructing the Off-target scoring input.");


        // New intergrated method (No need for external file, simply write to char vector)
        vector<char> queryDataSet;
        queryDataSet.reserve(pageLength * seqLength);
        vector<char> pamDataSet;
        pamDataSet.reserve(pageLength * 3);

        guidesInPage = 0;
        while (paginatorIterator != pageEnd)
        {
            string target23 = paginatorIterator->first;
            // Run time filtering
            if (!filterCandidateGuides(paginatorIterator->second, MODULE_SPECIFICITY, optimsationLevel, consensusN, toolCount)) {
                // Advance page end for each filtered out guide
                if (pageEnd != candidateGuides.end())
                {
                    pageEnd++;
                }
                paginatorIterator++;
                continue;
            }

            for (char c : target23.substr(0, 20))
            {
                queryDataSet.push_back(c);
            }

            for (char c : target23.substr(20))
            {
                pamDataSet.push_back(c);
            }
            guidesInPage++;
            paginatorIterator++;
        }



        printer(fmt::format("\t\t{} guides in this page.", commaify(guidesInPage)));

        // Call scoring method
        size_t seqLineLength = seqLength;

        size_t queryCount = guidesInPage;
        vector<uint64_t> querySignatures(queryCount);
        vector<double> querySignatureMitScores(queryCount);
        vector<double> querySignatureCfdScores(queryCount);

        /** Binary encode query sequences */
        omp_set_num_threads(threadCount);
#pragma omp parallel
        {
#pragma omp for
            for (int i = 0; i < queryCount; i++) {
                char* ptr = &queryDataSet[i * seqLineLength];
                uint64_t signature = sequenceToSignature(ptr, seqLength);
                querySignatures[i] = signature;
            }
        }

        // TODO: remove
        // Mutex for thread safety
        std::mutex countsMutex;
        // Timing
        std::chrono::duration<double> usefulTime;
        std::chrono::duration<double> wastedTime;
        // Neighbourhood counts
        vector<uint64_t> neighbourhoodCountTotal(sliceCount);
        vector<uint64_t> neighbourhoodCountUnique(sliceCount);
        // Mismatch counts 
        vector<vector<uint64_t>> offTargetCountTotal(sliceCount, vector<uint64_t>(21));
        vector<vector<uint64_t>> offTargetCountUnique(sliceCount, vector<uint64_t>(21));
        // Per guide counts 
        vector<vector<uint64_t>> perGuideCountTotal(2, vector<uint64_t>(querySignatures.size()));
        vector<vector<uint64_t>> perGuideCountUnique(2, vector<uint64_t>(querySignatures.size()));


        /** Begin scoring */
        omp_set_num_threads(threadCount);
#pragma omp parallel
        {
            unordered_map<uint64_t, unordered_set<uint64_t>> searchResults;
            vector<uint64_t> offtargetToggles(numOfftargetToggles);

            uint64_t* offtargetTogglesTail = offtargetToggles.data() + numOfftargetToggles - 1;

            // TODO: remove
            // Timing
            std::chrono::duration<double> usefulTimeLocal;
            std::chrono::duration<double> wastedTimeLocal;
            // Neighbourhood counts
            vector<uint64_t> neighbourhoodCountTotalLocal(sliceCount);
            vector<uint64_t> neighbourhoodCountUniqueLocal(sliceCount);
            // Mismatch counts 
            vector<vector<uint64_t>> offTargetCountTotalLocal(sliceCount, vector<uint64_t>(21));
            vector<vector<uint64_t>> offTargetCountUniqueLocal(sliceCount, vector<uint64_t>(21));
            // Per guide counts 
            vector<vector<uint64_t>> perGuideCountTotalLocal(2, vector<uint64_t>(querySignatures.size()));
            vector<vector<uint64_t>> perGuideCountUniqueLocal(2, vector<uint64_t>(querySignatures.size()));

            /** For each candidate guide */
#pragma omp for
            for (int searchIdx = 0; searchIdx < querySignatures.size(); searchIdx++) {

                auto searchSignature = querySignatures[searchIdx];

                /** Global scores */
                double totScoreMit = 0.0;
                double totScoreCfd = 0.0;

                int numOffTargetSitesScored = 0;
                double maximum_sum = (10000.0 - scoreThreshold * 100) / scoreThreshold;
                // TODO: uncommment
                //bool checkNextSlice = true;

                size_t sliceLimitOffset = 0;
                /** For each ISSL slice */
                for (size_t i = 0; i < sliceCount; i++) {
                    vector<int>& sliceMask = sliceMasks[i];
                    auto& sliceList = sliceLists[i];

                    uint64_t searchSlice = 0ULL;
                    for (int j = 0; j < sliceMask.size(); j++)
                    {
                        searchSlice |= ((searchSignature >> (sliceMask[j] * 2)) & 3ULL) << (j * 2);
                    }

                    size_t idx = sliceLimitOffset + searchSlice;

                    size_t signaturesInSlice = allSlicelistSizes[idx];
                    uint64_t* sliceOffset = sliceList[searchSlice];

                    /** For each off-target signature in slice */
                    for (size_t j = 0; j < signaturesInSlice; j++) {
                        // TODO: remove
                        time_point start = high_resolution_clock::now();

                        auto signatureWithOccurrencesAndId = sliceOffset[j];
                        auto signatureId = signatureWithOccurrencesAndId & 0xFFFFFFFFULL;
                        uint32_t occurrences = (signatureWithOccurrencesAndId >> (32));

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
                        /*int dist = popcnt(&mismatches, sizeof(uint64_t));*/
                        int dist = popcount64(mismatches);

                        if (dist >= 0 && dist <= maxDist) {

                            /** Prevent assessing the same off-target for multiple slices */
                            uint64_t seenOfftargetAlready = 0;
                            uint64_t* ptrOfftargetFlag = (offtargetTogglesTail - (signatureId / 64));
                            seenOfftargetAlready = (*ptrOfftargetFlag >> (signatureId % 64)) & 1ULL;


                            if (!seenOfftargetAlready) {
                                // Begin calculating MIT score
                                if (calcMit) {
                                    if (dist > 0) {
                                        totScoreMit += precalculatedScores[mismatches] * (double)occurrences;
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
                                    else if (dist > 0 && dist <= maxDist) {
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
                                            searchSigIdentityPos &= (3UL << (pos * 2));
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
                                            offtargetIdentityPos &= (3UL << (pos * 2));
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

                                *ptrOfftargetFlag |= (1ULL << (signatureId % 64));
                                numOffTargetSitesScored += occurrences;

                                // TODO: uncomment
                                /** Stop calculating global score early if possible */
                                //if (sm == ScoreMethod::mitAndCfd) {
                                //    if (totScoreMit > maximum_sum && totScoreCfd > maximum_sum) {
                                //        checkNextSlice = false;
                                //        break;
                                //    }
                                //}
                                //if (sm == ScoreMethod::mitOrCfd) {
                                //    if (totScoreMit > maximum_sum || totScoreCfd > maximum_sum) {
                                //        checkNextSlice = false;
                                //        break;
                                //    }
                                //}
                                //if (sm == ScoreMethod::avgMitCfd) {
                                //    if (((totScoreMit + totScoreCfd) / 2.0) > maximum_sum) {
                                //        checkNextSlice = false;
                                //        break;
                                //    }
                                //}
                                //if (sm == ScoreMethod::mit) {
                                //    if (totScoreMit > maximum_sum) {
                                //        checkNextSlice = false;
                                //        break;
                                //    }
                                //}
                                //if (sm == ScoreMethod::cfd) {
                                //    if (totScoreCfd > maximum_sum) {
                                //        checkNextSlice = false;
                                //        break;
                                //    }
                                //}
                            }
                            // TODO: remove
                            // Updated counts and timings
                            std::chrono::duration<double> computeTime = high_resolution_clock::now() - start;
                            // Timing
                            usefulTimeLocal += computeTime;
                            // Neighbourhood counts
                            neighbourhoodCountTotalLocal[i] += occurrences;
                            neighbourhoodCountUniqueLocal[i]++;
                            // Mismatch counts 
                            offTargetCountTotalLocal[i][dist] += occurrences;
                            offTargetCountUniqueLocal[i][dist]++;
                            // Per guide counts
                            perGuideCountTotalLocal[true][searchIdx] += occurrences;
                            perGuideCountUniqueLocal[true][searchIdx]++;
                        }
                        else
                        {
                            // TODO: remove
                            // Updated counts and timings
                            std::chrono::duration<double> computeTime = high_resolution_clock::now() - start;
                            // Timing
                            wastedTimeLocal += computeTime;
                            // Neighbourhood counts
                            neighbourhoodCountTotalLocal[i] += occurrences;
                            neighbourhoodCountUniqueLocal[i]++;
                            // Mismatch counts 
                            offTargetCountTotalLocal[i][dist] += occurrences;
                            offTargetCountUniqueLocal[i][dist]++;
                            // Per guide counts
                            perGuideCountTotalLocal[false][searchIdx] += occurrences;
                            perGuideCountUniqueLocal[false][searchIdx]++;
                        }
                    }

                    // TODO: uncomment
                    //if (!checkNextSlice)
                    //    break;
                    sliceLimitOffset += 1ULL << (sliceLens[i] * 2);
                }
                querySignatureMitScores[searchIdx] = 10000.0 / (100.0 + totScoreMit);
                querySignatureCfdScores[searchIdx] = 10000.0 / (100.0 + totScoreCfd);

                memset(offtargetToggles.data(), 0, sizeof(uint64_t) * offtargetToggles.size());
            }
            countsMutex.lock();
            usefulTime += usefulTimeLocal;
            wastedTime += wastedTimeLocal;
            for (size_t i = 0; i < neighbourhoodCountTotalLocal.size(); i++)
            {
                neighbourhoodCountTotal[i] += neighbourhoodCountTotalLocal[i];
                neighbourhoodCountUnique[i] += neighbourhoodCountUniqueLocal[i];
            }

            for (size_t i = 0; i < offTargetCountTotalLocal.size(); i++)
            {
                for (size_t j = 0; j < offTargetCountTotalLocal[i].size(); j++)
                {
                    offTargetCountTotal[i][j] += offTargetCountTotalLocal[i][j];
                    offTargetCountUnique[i][j] += offTargetCountUniqueLocal[i][j];

                }
            }

            for (size_t i = 0; i < perGuideCountTotalLocal.size(); i++)
            {
                for (size_t j = 0; j < perGuideCountTotalLocal[i].size(); j++)
                {
                    perGuideCountTotal[i][j] += perGuideCountTotalLocal[i][j];
                    perGuideCountUnique[i][j] += perGuideCountUniqueLocal[i][j];
                }

            }
            countsMutex.unlock();
        }

        printer("\tStarting to process the Off-target scoring results.");

        //printer(fmt::format("Offtargets tested {}", commaify(neighbourhoodCount)));

        //    for (int i = 0; i < 21; i++)
        //    {
        //        printer(fmt::format("\tMismatch {}\tTotal {}", i, commaify(offTargetCount[i])));
        //    }

        //std::filesystem::path outputPath(std::filesystem::path(ISSLIndex).parent_path());
        //std::ofstream outputFile;

        //std::map<long long, long long> truePerGuideCount;
        //for (long long trueOTCount : perGuideCount[true])
        //{
        //    if (truePerGuideCount.count(trueOTCount))
        //    {
        //        truePerGuideCount[trueOTCount]++;
        //    }
        //    else
        //    {
        //        truePerGuideCount[trueOTCount] = 1;
        //    }
        //}

        //outputFile.open(outputPath / "_output" / "truePerGuideCountTotal.txt");
        //for (auto const& x : truePerGuideCount)
        //{
        //    outputFile << fmt::format("{}:{}\n",x.first,x.second);
        //}
        //outputFile.close();

        //std::map<long long, long long> falsePerGuideCount;
        //for (long long falseOTCount : perGuideCount[false])
        //{
        //    if (falsePerGuideCount.count(falseOTCount))
        //    {
        //        falsePerGuideCount[falseOTCount]++;
        //    }
        //    else
        //    {
        //        falsePerGuideCount[falseOTCount] = 1;
        //    }
        //}

        //outputFile.open(outputPath / "_output" / "falsePerGuideCountTotal.txt");
        //for (auto const& x : falsePerGuideCount)
        //{
        //    outputFile << fmt::format("{}:{}\n", x.first, x.second);
        //}
        //outputFile.close();

        // TODO: remove
        printer(fmt::format("Wasted time: {}\tUseful time: {}", wastedTime.count(), usefulTime.count()));

        for (int i = 0; i < neighbourhoodCountTotal.size(); i++)
        {
            printer(fmt::format("Slice {}\tTotal {}\tUnique {}", i + 1, commaify(neighbourhoodCountTotal[i]), commaify(neighbourhoodCountUnique[i])));
            for (int j = 0; j < 21; j++)
            {
                printer(fmt::format("\tMismatch {}\tTotal {}\tUnique {}", j, commaify(offTargetCountTotal[i][j]), commaify(offTargetCountUnique[i][j])));
            }
        }

        std::filesystem::path outputPath(std::filesystem::path(ISSLIndex).parent_path());

        std::thread countTrueTotalThread(countTrueTotal, outputPath, perGuideCountTotal);
        std::thread countTrueUniqueThread(countTrueUnique, outputPath, perGuideCountUnique);
        std::thread countFalseTotalThread(countFalseTotal, outputPath, perGuideCountTotal);
        std::thread countFalseUniqueThread(countFalseUnique, outputPath, perGuideCountUnique);

        countTrueTotalThread.join();
        countTrueUniqueThread.join();
        countFalseTotalThread.join();
        countFalseUniqueThread.join();

        for (size_t searchIdx = 0; searchIdx < querySignatures.size(); searchIdx++) {
            string target20 = signatureToSequence(querySignatures[searchIdx], seqLength);
            char* pamPtr = &pamDataSet[searchIdx * 3];
            string pam(pamPtr, pamPtr + 3);
            string target23 = target20 + pam;
            candidateGuides[target23]["mitOfftargetscore"] = calcMit ? std::to_string(querySignatureMitScores[searchIdx]) : "-1";
            candidateGuides[target23]["cfdOfftargetscore"] = calcCfd ? std::to_string(querySignatureCfdScores[searchIdx]) : "-1";

            if (sm == ScoreMethod::mit)
            {
                if (std::stof(candidateGuides[target23]["mitOfftargetscore"]) < scoreThreshold)
                {
                    candidateGuides[target23]["passedOffTargetScore"] = CODE_REJECTED;
                    failedCount++;
                }
                else { candidateGuides[target23]["passedOffTargetScore"] = CODE_ACCEPTED; }
            }

            else if (sm == ScoreMethod::cfd)
            {
                if (std::stof(candidateGuides[target23]["cfdOfftargetscore"]) < scoreThreshold)
                {
                    candidateGuides[target23]["passedOffTargetScore"] = CODE_REJECTED;
                    failedCount++;
                }
                else { candidateGuides[target23]["passedOffTargetScore"] = CODE_ACCEPTED; }
            }

            else if (sm == ScoreMethod::mitAndCfd)
            {
                if ((std::stof(candidateGuides[target23]["mitOfftargetscore"]) < scoreThreshold) && (std::stof(candidateGuides[target23]["cfdOfftargetscore"]) < scoreThreshold))
                {
                    candidateGuides[target23]["passedOffTargetScore"] = CODE_REJECTED;
                    failedCount++;
                }
                else { candidateGuides[target23]["passedOffTargetScore"] = CODE_ACCEPTED; }
            }

            else if (sm == ScoreMethod::mitOrCfd)
            {
                if ((std::stof(candidateGuides[target23]["mitOfftargetscore"]) < scoreThreshold) || (std::stof(candidateGuides[target23]["cfdOfftargetscore"]) < scoreThreshold))
                {
                    candidateGuides[target23]["passedOffTargetScore"] = CODE_REJECTED;
                    failedCount++;
                }
                else { candidateGuides[target23]["passedOffTargetScore"] = CODE_ACCEPTED; }
            }

            else if (sm == ScoreMethod::avgMitCfd)
            {
                if (((std::stof(candidateGuides[target23]["mitOfftargetscore"]) + std::stof(candidateGuides[target23]["cfdOfftargetscore"])) / 2) < scoreThreshold)
                {
                    candidateGuides[target23]["passedOffTargetScore"] = CODE_REJECTED;
                    failedCount++;
                }
                else { candidateGuides[target23]["passedOffTargetScore"] = CODE_ACCEPTED; }
            }
            testedCount++;
        }

        printer(fmt::format("\t{} of {} failed here.", commaify(failedCount), commaify(testedCount)));
        pageStart = pageEnd;
    }

}

/// Returns the size (bytes) of the file at `path`
size_t ISSL2Stage::getFileSize(const char* path)
{
    struct p_stat64 statBuf;
    p_stat64(path, &statBuf);
    return statBuf.st_size;
}

/**
 * Binary encode genetic string `ptr`
 *
 * For example,
 *   ATCG becomes
 *   00 11 01 10  (buffer with leading zeroes to encode as 64-bit unsigned int)
 *
 * @param[in] ptr the string containing ATCG to binary encode
 */
uint64_t ISSL2Stage::sequenceToSignature(const char* ptr, const size_t& seqLength)
{
    uint64_t signature = 0;
    for (size_t j = 0; j < seqLength; j++) {
        signature |= (uint64_t)(nucleotideIndex[*ptr]) << (j * 2);
        ptr++;
    }
    return signature;
}

/**
 * Binary encode genetic string `ptr`
 *
 * For example,
 *   00 11 01 10 becomes (as 64-bit unsigned int)
 *    A  T  C  G  (without spaces)
 *
 * @param[in] signature the binary encoded genetic string
 */
string ISSL2Stage::signatureToSequence(uint64_t signature, const size_t& seqLength)
{
    string sequence = string(seqLength, ' ');
    for (size_t j = 0; j < seqLength; j++) {
        sequence[j] = signatureIndex[(signature >> (j * 2)) & 0x3];
    }
    return sequence;
}

#if defined(_WIN64)
#pragma pop_macro("close")
#endif