#include "../include/ISSLOffTargetScoring.hpp"
// TODO: remove
#include <map>
#include <chrono>
#include <fstream>
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

ISSLOffTargetScoring::ISSLOffTargetScoring(ConfigManager& cm) :
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

void ISSLOffTargetScoring::run(unordered_map<string, unordered_map<string, string>>& candidateGuides)
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

    ///** Begin reading the binary encoded ISSL, structured as:
    // *      - a header (5 items)
    // *      - list of slice ranges
    // *      - precalcuated local MIT scores
    // *      - all binary-encoded off-target sites
    // *      - slice list sizes
    // *      - slice contents
    // */
    //FILE* fp = fopen(ISSLIndex.c_str(), "rb");



    ///** The index contains a fixed-sized header
    // *      - the number of unique off-targets in the index
    // *      - the length of an off-target
    // *      - the total number of off-targets
    // *      - the length of the slice range array
    // *      - the number of precalculated MIT scores
    // */
    //vector<size_t> slicelistHeader(50);

    //if (fread(slicelistHeader.data(), sizeof(size_t), slicelistHeader.size(), fp) == 0) {
    //    throw std::runtime_error("Error reading index: header invalid\n");
    //}

    //size_t offtargetsCount = slicelistHeader[0];
    //size_t seqLength = slicelistHeader[1];
    //size_t seqCount = slicelistHeader[2];
    //size_t sliceRangeCount = slicelistHeader[3];
    //size_t scoresCount = slicelistHeader[4];

    ///**
    //* Read the slice ranges from header
    //*/
    //vector<size_t> sliceRanges;
    //for (int i = 0; i < sliceRangeCount; i++)
    //{
    //    size_t sliceRange;
    //    fread(&sliceRange, sizeof(size_t), 1, fp);
    //    sliceRanges.push_back(sliceRange);
    //}

    //vector<size_t> sliceLens;
    //for (int i = 0; i < sliceRanges.size(); i = i + 2) {
    //    sliceLens.push_back((sliceRanges[i + 1] + 2) - (sliceRanges[i]));
    //}

    //size_t sliceCount = sliceLens.size();

    ///** Read in the precalculated MIT scores
    // *      - `mask` is a 2-bit encoding of mismatch positions
    // *          For example,
    // *              00 01 01 00 01  indicates mismatches in positions 1, 3 and 4
    // *
    // *      - `score` is the local MIT score for this mismatch combination
    // */
    //phmap::flat_hash_map<uint64_t, double> precalculatedScores;

    //for (int i = 0; i < scoresCount; i++) {
    //    uint64_t mask = 0;
    //    double score = 0.0;
    //    fread(&mask, sizeof(uint64_t), 1, fp);
    //    fread(&score, sizeof(double), 1, fp);

    //    precalculatedScores.insert(pair<uint64_t, double>(mask, score));
    //}

    ///** Load in all of the off-target sites */
    //vector<uint64_t> offtargets(offtargetsCount);
    //if (fread(offtargets.data(), sizeof(uint64_t), offtargetsCount, fp) == 0) {
    //    throw std::runtime_error("Error reading index: loading off-target sequences failed\n");
    //}

    ///** The number of signatures embedded per slice
    // *
    // *      These counts are stored contiguously
    // *
    // */
    //size_t sliceListCount = 0;
    //for (int i = 0; i < sliceCount; i++)
    //{
    //    sliceListCount += 1ULL << sliceLens[i];
    //}
    //vector<size_t> allSlicelistSizes(sliceListCount);

    //if (fread(allSlicelistSizes.data(), sizeof(size_t), allSlicelistSizes.size(), fp) == 0) {
    //    throw std::runtime_error("Error reading index: reading slice list sizes failed\n");
    //}

    ///** The contents of the slices
    // *
    // *      Stored contiguously
    // *
    // *      Each signature (64-bit) is structured as:
    // *          <occurrences 32-bit><off-target-id 32-bit>
    // */
    //vector<uint64_t> allSignatures(offtargetsCount * sliceCount);

    //if (fread(allSignatures.data(), sizeof(uint64_t), allSignatures.size(), fp) == 0) {
    //    throw std::runtime_error("Error reading index: reading slice contents failed\n");
    //}

    ///** End reading the index */
    //fclose(fp);

    ///** Prevent assessing an off-target site for multiple slices
    // *
    // *      Create enough 1-bit "seen" flags for the off-targets
    // *      We only want to score a candidate guide against an off-target once.
    // *      The least-significant bit represents the first off-target
    // *      0 0 0 1   0 1 0 0   would indicate that the 3rd and 5th off-target have been seen.
    // *      The CHAR_BIT macro tells us how many bits are in a byte (C++ >= 8 bits per byte)
    // */
    //uint64_t numOfftargetToggles = (offtargetsCount / ((size_t)sizeof(uint64_t) * (size_t)CHAR_BIT)) + 1;


    ///** Start constructing index in memory
    // *
    // *      To begin, reverse the contiguous storage of the slices,
    // *         into the following:
    // *
    // *         + Slice 0 :
    // *         |---- AAAA : <slice contents>
    // *         |---- AAAC : <slice contents>
    // *         |----  ...
    // *         |
    // *         + Slice 1 :
    // *         |---- AAAA : <slice contents>
    // *         |---- AAAC : <slice contents>
    // *         |---- ...
    // *         | ...
    // */

    //vector<vector<uint64_t*>> sliceLists(sliceLens.size());
    //// Assign sliceLists size based on each slice length
    //for (int i = 0; i < sliceLens.size(); i++)
    //{
    //    sliceLists[i] = vector<uint64_t*>(1ULL << sliceLens[i]);
    //}

    //uint64_t* offset = allSignatures.data();
    //size_t sliceLimitOffset = 0;
    //for (size_t i = 0; i < sliceCount; i++) {
    //    size_t sliceLimit = 1 << sliceLens[i];
    //    for (size_t j = 0; j < sliceLimit; j++) {
    //        size_t idx = sliceLimitOffset + j;
    //        sliceLists[i][j] = offset;
    //        offset += allSlicelistSizes[idx];
    //    }
    //    sliceLimitOffset += sliceLimit;
    //}

    // TODO: Move this to external header file like CFD penalities
    // Precalculate all the scores
    unordered_map<uint64_t, double> precalculatedScores;

    int maxDist = 4;
    size_t scoresCount = 0;

    for (int i = 1; i <= maxDist; i++) {
        vector<uint64_t> tempMasks;
        tempMasks = computeMasksTwoBit(20, i);
        for (auto mask : tempMasks) {
            double score = predictMITLocalScore(mask);
            precalculatedScores.insert(pair<uint64_t, double>(mask, score));
            scoresCount++;
        }
    }

    // Load cluster results

    int seqLength = 20;

    unordered_map <string, vector<uint64_t>> kClusterLists;
    unordered_map <uint64_t, string> sigToCluster;
    std::ifstream inFile;
    inFile.open(ISSLIndex, std::ios::binary | std::ios_base::in);
    string cluster;

    for (std::string line; std::getline(inFile, line);)
    {
        line = trim(line);
        if (line[0] == '>')
        {
            cluster = line.substr(1);
        }
        else
        {
            uint64_t signature = sequenceToSignature(line.c_str(), line.length());
            kClusterLists[cluster].push_back(signature);
            sigToCluster[signature] = cluster;
        }
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
        // Neighbourhood counts
        long long neighbourhoodCount = 0;
        // Mismatch counts
        vector<long long> offTargetCount(21, 0);
        // Guide counts
        vector<vector<long long>> perGuideCount(2, vector<long long>(querySignatures.size(), 0));

        /** OT scoring loop */
        omp_set_num_threads(threadCount);
        //#pragma omp parallel
        {
            //#pragma omp for
            for (int i = 0; i < queryCount; i++) {
                uint64_t searchSignature = querySignatures[i];

                /** Global scores */
                double totScoreMit = 0.0;
                double totScoreCfd = 0.0;

                int numOffTargetSitesScored = 0;
                double maximum_sum = (10000.0 - scoreThreshold * 100) / scoreThreshold;
                bool checkNextOT = true;

                // TODO: remove
                countsMutex.lock();
                neighbourhoodCount += kClusterLists[sigToCluster[searchSignature]].size();
                countsMutex.unlock();

                // Get Cluster list
                for (auto const& offTarget : kClusterLists[sigToCluster[searchSignature]])
                {

                    uint64_t xoredSignatures = searchSignature ^ offTarget;
                    uint64_t evenBits = xoredSignatures & 0xAAAAAAAAAAAAAAAAULL;
                    uint64_t oddBits = xoredSignatures & 0x5555555555555555ULL;
                    uint64_t mismatches = (evenBits >> 1) | oddBits;
                    int dist = popcount64(mismatches);

                    if (dist >= 0 && dist <= maxDist) {


                        // Begin calculating MIT score
                        if (calcMit) {
                            if (dist > 0) {
                                totScoreMit += precalculatedScores[mismatches];
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
                                    uint64_t offtargetIdentityPos = offTarget;
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
                            totScoreCfd += cfdScore;
                        }


                        countsMutex.lock();
                        perGuideCount[true][i]++;
                        offTargetCount[dist]++;
                        countsMutex.unlock();
                    }
                    else
                    {
                        countsMutex.lock();
                        perGuideCount[false][i]++;
                        offTargetCount[dist]++;
                        countsMutex.unlock();
                    }

                    // TODO: uncomment
                    /** Stop calculating global score early if possible */
                    //if (sm == ScoreMethod::mitAndCfd) {
                    //    if (totScoreMit > maximum_sum && totScoreCfd > maximum_sum) {
                    //        break;
                    //    }
                    //}
                    //if (sm == ScoreMethod::mitOrCfd) {
                    //    if (totScoreMit > maximum_sum || totScoreCfd > maximum_sum) {
                    //        break;
                    //    }
                    //}
                    //if (sm == ScoreMethod::avgMitCfd) {
                    //    if (((totScoreMit + totScoreCfd) / 2.0) > maximum_sum) {
                    //        break;
                    //    }
                    //}
                    //if (sm == ScoreMethod::mit) {
                    //    if (totScoreMit > maximum_sum) {
                    //        break;
                    //    }
                    //}
                    //if (sm == ScoreMethod::cfd) {
                    //    if (totScoreCfd > maximum_sum) {
                    //        break;
                    //    }
                    //}
                }
                querySignatureMitScores[i] = 10000.0 / (100.0 + totScoreMit);
                querySignatureCfdScores[i] = 10000.0 / (100.0 + totScoreCfd);
            }
        }

        ///** Begin scoring */
        //omp_set_num_threads(threadCount);
        //#pragma omp parallel
        //{
        //    unordered_map<uint64_t, unordered_set<uint64_t>> searchResults;
        //    vector<uint64_t> offtargetToggles(numOfftargetToggles);

        //    uint64_t* offtargetTogglesTail = offtargetToggles.data() + numOfftargetToggles - 1;

        //    /** For each candidate guide */
        //    #pragma omp for
        //    for (int searchIdx = 0; searchIdx < querySignatures.size(); searchIdx++) {

        //        auto searchSignature = querySignatures[searchIdx];

        //        /** Global scores */
        //        double totScoreMit = 0.0;
        //        double totScoreCfd = 0.0;

        //        int numOffTargetSitesScored = 0;
        //        double maximum_sum = (10000.0 - scoreThreshold * 100) / scoreThreshold;
        //        // TODO: uncommment
        //        //bool checkNextSlice = true;

        //        size_t sliceLimitOffset = 0;
        //        /** For each ISSL slice */
        //        for (size_t i = 0; i < sliceCount; i++) {
        //            size_t sliceWidth = sliceLens[i];
        //            size_t sliceLimit = 1 << sliceWidth;
        //            uint64_t sliceMask = sliceLimit - 1;
        //            int sliceShift = sliceRanges[i*2];
        //            sliceMask = sliceMask << sliceShift;
        //            auto& sliceList = sliceLists[i];

        //            uint64_t searchSlice = (searchSignature & sliceMask) >> sliceShift;

        //            size_t idx = sliceLimitOffset + searchSlice;

        //            size_t signaturesInSlice = allSlicelistSizes[idx];
        //            uint64_t* sliceOffset = sliceList[searchSlice];

        //            /** For each off-target signature in slice */
        //            for (size_t j = 0; j < signaturesInSlice; j++) {
        //                // TODO: remove
        //                time_point start = high_resolution_clock::now();

        //                auto signatureWithOccurrencesAndId = sliceOffset[j];
        //                auto signatureId = signatureWithOccurrencesAndId & 0xFFFFFFFFULL;
        //                uint32_t occurrences = (signatureWithOccurrencesAndId >> (32));

        //                /** Find the positions of mismatches
        //                 *
        //                 *  Search signature (SS):    A  A  T  T    G  C  A  T
        //                 *                           00 00 11 11   10 01 00 11
        //                 *
        //                 *        Off-target (OT):    A  T  A  T    C  G  A  T
        //                 *                           00 11 00 11   01 10 00 11
        //                 *
        //                 *                SS ^ OT:   00 00 11 11   10 01 00 11
        //                 *                         ^ 00 11 00 11   01 10 00 11
        //                 *                  (XORd) = 00 11 11 00   11 11 00 00
        //                 *
        //                 *        XORd & evenBits:   00 11 11 00   11 11 00 00
        //                 *                         & 10 10 10 10   10 10 10 10
        //                 *                   (eX)  = 00 10 10 00   10 10 00 00
        //                 *
        //                 *         XORd & oddBits:   00 11 11 00   11 11 00 00
        //                 *                         & 01 01 01 01   01 01 01 01
        //                 *                   (oX)  = 00 01 01 00   01 01 00 00
        //                 *
        //                 *         (eX >> 1) | oX:   00 01 01 00   01 01 00 00 (>>1)
        //                 *                         | 00 01 01 00   01 01 00 00
        //                 *            mismatches   = 00 01 01 00   01 01 00 00
        //                 *
        //                 *   popcount(mismatches):   4
        //                 */
        //                uint64_t xoredSignatures = searchSignature ^ offtargets[signatureId];
        //                uint64_t evenBits = xoredSignatures & 0xAAAAAAAAAAAAAAAAULL;
        //                uint64_t oddBits = xoredSignatures & 0x5555555555555555ULL;
        //                uint64_t mismatches = (evenBits >> 1) | oddBits;
        //                /*int dist = popcnt(&mismatches, sizeof(uint64_t));*/
        //                int dist = popcount64(mismatches);

        //                if (dist >= 0 && dist <= maxDist) {

        //                    /** Prevent assessing the same off-target for multiple slices */
        //                    uint64_t seenOfftargetAlready = 0;
        //                    uint64_t* ptrOfftargetFlag = (offtargetTogglesTail - (signatureId / 64));
        //                    seenOfftargetAlready = (*ptrOfftargetFlag >> (signatureId % 64)) & 1ULL;


        //                    if (!seenOfftargetAlready) {
        //                        // Begin calculating MIT score
        //                        if (calcMit) {
        //                            if (dist > 0) {
        //                                totScoreMit += precalculatedScores[mismatches] * (double)occurrences;
        //                            }
        //                        }

        //                        // Begin calculating CFD score
        //                        if (calcCfd) {
        //                            /** "In other words, for the CFD score, a value of 0
        //                             *      indicates no predicted off-target activity whereas
        //                             *      a value of 1 indicates a perfect match"
        //                             *      John Doench, 2016.
        //                             *      https://www.nature.com/articles/nbt.3437
        //                            */
        //                            double cfdScore = 0;
        //                            if (dist == 0) {
        //                                cfdScore = 1;
        //                            }
        //                            else if (dist > 0 && dist <= maxDist) {
        //                                cfdScore = cfdPamPenalties[0b1010]; // PAM: NGG, TODO: do not hard-code the PAM

        //                                for (size_t pos = 0; pos < 20; pos++) {
        //                                    size_t mask = pos << 4;

        //                                    /** Create the mask to look up the position-identity score
        //                                     *      In Python... c2b is char to bit
        //                                     *       mask = pos << 4
        //                                     *       mask |= c2b[sgRNA[pos]] << 2
        //                                     *       mask |= c2b[revcom(offTaret[pos])]
        //                                     *
        //                                     *      Find identity at `pos` for search signature
        //                                     *      example: find identity in pos=2
        //                                     *       Recall ISSL is inverted, hence:
        //                                     *                   3'-  T  G  C  C  G  A -5'
        //                                     *       start           11 10 01 01 10 00
        //                                     *       3UL << pos*2    00 00 00 11 00 00
        //                                     *       and             00 00 00 01 00 00
        //                                     *       shift           00 00 00 00 01 00
        //                                     */
        //                                    uint64_t searchSigIdentityPos = searchSignature;
        //                                    searchSigIdentityPos &= (3UL << (pos * 2));
        //                                    searchSigIdentityPos = searchSigIdentityPos >> (pos * 2);
        //                                    searchSigIdentityPos = searchSigIdentityPos << 2;

        //                                    /** Find identity at `pos` for offtarget
        //                                     *      Example: find identity in pos=2
        //                                     *      Recall ISSL is inverted, hence:
        //                                     *                  3'-  T  G  C  C  G  A -5'
        //                                     *      start           11 10 01 01 10 00
        //                                     *      3UL<<pos*2      00 00 00 11 00 00
        //                                     *      and             00 00 00 01 00 00
        //                                     *      shift           00 00 00 00 00 01
        //                                     *      rev comp 3UL    00 00 00 00 00 10 (done below)
        //                                     */
        //                                    uint64_t offtargetIdentityPos = offtargets[signatureId];
        //                                    offtargetIdentityPos &= (3UL << (pos * 2));
        //                                    offtargetIdentityPos = offtargetIdentityPos >> (pos * 2);

        //                                    /** Complete the mask
        //                                     *      reverse complement (^3UL) `offtargetIdentityPos` here
        //                                     */
        //                                    mask = (mask | searchSigIdentityPos | (offtargetIdentityPos ^ 3UL));

        //                                    if (searchSigIdentityPos >> 2 != offtargetIdentityPos) {
        //                                        cfdScore *= cfdPosPenalties[mask];
        //                                    }
        //                                }
        //                            }
        //                            totScoreCfd += cfdScore * (double)occurrences;
        //                        }

        //                        *ptrOfftargetFlag |= (1ULL << (signatureId % 64));
        //                        numOffTargetSitesScored += occurrences;

        //                        // TODO: uncomment
        //                        /** Stop calculating global score early if possible */
        //                        /*if (sm == ScoreMethod::mitAndCfd) {
        //                            if (totScoreMit > maximum_sum && totScoreCfd > maximum_sum) {
        //                                checkNextSlice = false;
        //                                break;
        //                            }
        //                        }
        //                        if (sm == ScoreMethod::mitOrCfd) {
        //                            if (totScoreMit > maximum_sum || totScoreCfd > maximum_sum) {
        //                                checkNextSlice = false;
        //                                break;
        //                            }
        //                        }
        //                        if (sm == ScoreMethod::avgMitCfd) {
        //                            if (((totScoreMit + totScoreCfd) / 2.0) > maximum_sum) {
        //                                checkNextSlice = false;
        //                                break;
        //                            }
        //                        }
        //                        if (sm == ScoreMethod::mit) {
        //                            if (totScoreMit > maximum_sum) {
        //                                checkNextSlice = false;
        //                                break;
        //                            }
        //                        }
        //                        if (sm == ScoreMethod::cfd) {
        //                            if (totScoreCfd > maximum_sum) {
        //                                checkNextSlice = false;
        //                                break;
        //                            }
        //                        }*/
        //                    }
        //                    // TODO: remove
        //                    // Updated counts and timings
        //                    std::chrono::duration<double> computeTime = high_resolution_clock::now() - start;
        //                    countsMutex.lock();
        //                    // Timing
        //                    usefulTime += computeTime;
        //                    // Neighbourhood counts
        //                    neighbourhoodCountTotal[i] += occurrences;
        //                    neighbourhoodCountUnique[i]++;
        //                    // Mismatch counts 
        //                    offTargetCountTotal[i][dist] += occurrences;
        //                    offTargetCountUnique[i][dist]++;
        //                    // Per guide counts
        //                    perGuideCountTotal[true][searchIdx] += occurrences;
        //                    perGuideCountUnique[true][searchIdx]++;
        //                    countsMutex.unlock();
        //                }
        //                else
        //                {
        //                    // TODO: remove
        //                    // Updated counts and timings
        //                    std::chrono::duration<double> computeTime = high_resolution_clock::now() - start;
        //                    countsMutex.lock();
        //                    // Timing
        //                    wastedTime += computeTime;
        //                    // Neighbourhood counts
        //                    neighbourhoodCountTotal[i] += occurrences;
        //                    neighbourhoodCountUnique[i]++;
        //                    // Mismatch counts 
        //                    offTargetCountTotal[i][dist] += occurrences;
        //                    offTargetCountUnique[i][dist]++;
        //                    // Per guide counts
        //                    perGuideCountTotal[false][searchIdx] += occurrences;
        //                    perGuideCountUnique[false][searchIdx]++;
        //                    countsMutex.unlock();
        //                    
        //                }
        //            }

        //            // TODO: uncomment
        //            //if (!checkNextSlice)
        //            //    break;
        //            sliceLimitOffset += sliceLimit;
        //        }
        //        querySignatureMitScores[searchIdx] = 10000.0 / (100.0 + totScoreMit);
        //        querySignatureCfdScores[searchIdx] = 10000.0 / (100.0 + totScoreCfd);

        //        memset(offtargetToggles.data(), 0, sizeof(uint64_t) * offtargetToggles.size());
        //    }

        //}

        printer("\tStarting to process the Off-target scoring results.");
        
        printer(fmt::format("Offtargets tested {}", commaify(neighbourhoodCount)));

            for (int i = 0; i < 21; i++)
            {
                printer(fmt::format("\tMismatch {}\tTotal {}", i, commaify(offTargetCount[i])));
            }

        std::filesystem::path outputPath(std::filesystem::path(ISSLIndex).parent_path());
        std::ofstream outputFile;

        std::map<long long, long long> truePerGuideCount;
        for (long long trueOTCount : perGuideCount[true])
        {
            if (truePerGuideCount.count(trueOTCount))
            {
                truePerGuideCount[trueOTCount]++;
            }
            else
            {
                truePerGuideCount[trueOTCount] = 1;
            }
        }

        outputFile.open(outputPath / "_output" / "truePerGuideCountTotal.txt");
        for (auto const& x : truePerGuideCount)
        {
            outputFile << fmt::format("{}:{}\n",x.first,x.second);
        }
        outputFile.close();

        std::map<long long, long long> falsePerGuideCount;
        for (long long falseOTCount : perGuideCount[true])
        {
            if (falsePerGuideCount.count(falseOTCount))
            {
                falsePerGuideCount[falseOTCount]++;
            }
            else
            {
                falsePerGuideCount[falseOTCount] = 1;
            }
        }

        outputFile.open(outputPath / "_output" / "falsePerGuideCountTotal.txt");
        for (auto const& x : falsePerGuideCount)
        {
            outputFile << fmt::format("{}:{}\n", x.first, x.second);
        }
        outputFile.close();

        //// TODO: remove
        //printer(fmt::format("Wasted time: {}\tUseful time: {}", wastedTime.count(), usefulTime.count()));

        //for (int i = 0; i < neighbourhoodCountTotal.size(); i++)
        //{
        //    printer(fmt::format("Slice {}\tTotal {}\tUnique {}", i+1, commaify(neighbourhoodCountTotal[i]), commaify(neighbourhoodCountUnique[i])));
        //    for (int j = 0; j < 21; j++)
        //    {
        //        printer(fmt::format("\tMismatch {}\tTotal {}\tUnique {}", j, commaify(offTargetCountTotal[i][j]), commaify(offTargetCountUnique[i][j])));
        //    }
        //}

        //std::filesystem::path outputPath(std::filesystem::path(ISSLIndex).parent_path());

        //std::ofstream outputFile;
        //std::map<long long, long long> truePerGuideCountTotal;
        //for (long long trueOTCount : perGuideCountTotal[true])
        //{
        //    if (truePerGuideCountTotal.count(trueOTCount))
        //    {
        //        truePerGuideCountTotal[trueOTCount]++;
        //    }
        //    else
        //    {
        //        truePerGuideCountTotal[trueOTCount] = 1;
        //    }
        //}

        //outputFile.open(outputPath / "_output" / "truePerGuideCountTotal.txt");
        //for (auto const& x : truePerGuideCountTotal)
        //{
        //    outputFile << fmt::format("{}:{}\n",x.first,x.second);
        //}
        //outputFile.close();

        //std::map<long long, long long> truePerGuideCountUnique;
        //for (long long trueOTCount : perGuideCountUnique[true])
        //{
        //    if (truePerGuideCountUnique.count(trueOTCount))
        //    {
        //        truePerGuideCountUnique[trueOTCount]++;
        //    }
        //    else
        //    {
        //        truePerGuideCountUnique[trueOTCount] = 1;
        //    }
        //}

        //outputFile.open(outputPath / "_output" / "truePerGuideCountUnique.txt");
        //for (auto const& x : truePerGuideCountUnique)
        //{
        //    outputFile << fmt::format("{}:{}\n", x.first, x.second);
        //}
        //outputFile.close();

        //std::map<long long, long long> falsePerGuideCountTotal;
        //for (long long trueOTCount : perGuideCountTotal[false])
        //{
        //    if (falsePerGuideCountTotal.count(trueOTCount))
        //    {
        //        falsePerGuideCountTotal[trueOTCount]++;
        //    }
        //    else
        //    {
        //        falsePerGuideCountTotal[trueOTCount] = 1;
        //    }
        //}

        //outputFile.open(outputPath / "_output" / "falsePerGuideCountTotal.txt");
        //for (auto const& x : falsePerGuideCountTotal)
        //{
        //    outputFile << fmt::format("{}:{}\n", x.first, x.second);
        //}
        //outputFile.close();

        //std::map<long long, long long> falsePerGuideCountUnique;
        //for (long long trueOTCount : perGuideCountUnique[false])
        //{
        //    if (falsePerGuideCountUnique.count(trueOTCount))
        //    {
        //        falsePerGuideCountUnique[trueOTCount]++;
        //    }
        //    else
        //    {
        //        falsePerGuideCountUnique[trueOTCount] = 1;
        //    }
        //}

        //outputFile.open(outputPath / "_output" / "falsePerGuideCountUnique.txt");
        //for (auto const& x : falsePerGuideCountUnique)
        //{
        //    outputFile << fmt::format("{}:{}\n", x.first, x.second);
        //}
        //outputFile.close();


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
    }

}

/// Returns the size (bytes) of the file at `path`
size_t ISSLOffTargetScoring::getFileSize(const char* path)
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
uint64_t ISSLOffTargetScoring::sequenceToSignature(const char* ptr, const size_t& seqLength)
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
string ISSLOffTargetScoring::signatureToSequence(uint64_t signature, const size_t& seqLength)
{
    string sequence = string(seqLength, ' ');
    for (size_t j = 0; j < seqLength; j++) {
        sequence[j] = signatureIndex[(signature >> (j * 2)) & 0x3];
    }
    return sequence;
}


/**
 * computeMasksTwoBit
 *
 * This function generates all of the two bit combinations of mismatches.
 *
 * @param seqLength, The length of the seqeunce
 * @param mismatches, The number of mismatches allowed
 * @return A vector of two bit combinations
 */
vector<uint64_t> ISSLOffTargetScoring::computeMasksTwoBit(int seqLength, int mismatches) {
    vector<uint64_t> masks;

    // There will only be one valid combination when mismatches == seqLength.
    if (mismatches != seqLength) {

        // There are more mismatchs to assign.
        if (mismatches > 0) {
            /**
            * This loop assigns a mismatch to the end of the seqeunce length.
            * It will then recursively call itself, reducing the seqeunce length and mismatch count.
            * (1ULL << (seqLength - 1) * 2) is responsible for assigning a mismatch at the end of the current seq length.
            * By adding the term above to the recursively returned, it will combine them into a unique combination.
            * E.g
            *   Starting with empty seq, seqeunce length 20, mismatch 2:
            *   00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00
            *
            *   Assign (1ULL << (seqLength - 1) * 2):
            *   10 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00
            *
            *   Recursively call function on remaining seq length, seqeunce length 19, mismatch 1:
            *   10 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00
            *
            *   Input to next call, length 19, mismatch 1:
            *   00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00
            *
            *   Eventually returned values will be:
            *   00 10 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00
            *   00 00 10 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00
            *      ....................................................................................................................
            *                                                                                                                     10 00
            *                                                                                                                        10
            *   The loop will add all the recursively genereated values with the original value will give unique combinations:
            *   10 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00
            * +
            *   00 10 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00
            *
            *   10 10 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00
            *
            *   TL:DR This loop reduces both mismatches and total sequence length to generate all combinations from original sequence.
            */
            for (auto mask : computeMasksTwoBit(seqLength - 1, mismatches - 1)) {

                masks.push_back((1ULL << (seqLength - 1) * 2) + mask);
            }

            /**
            * This loop is run after all combinations are generated with the mismatch base on the first sequence length.
            * By shrinking the sequence length it will shift the last possible mismatch and generate a new sequence.
            * This recursive call will result in the process being repeated on this new sequence.
            * E.g
            *   Considered the first mask genereated by (1ULL << (seqLength - 1) * 2):
            *   10 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00
            *
            *   By reducing the sequence size the result of (1ULL << (seqLength - 1) * 2) will be:
            *   00 10 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00
            *
            *   This sequence will then be used as the starting point for the first loop and the process will repeat:
            *   (First loop in recursive call)
            *   Recursively call function on remaining seq length, seqeunce length 19, mismatch 1:
            *   00 10 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00
            *
            *   Input to next call, length 18, mismatch 1:
            *   00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00
            *
            *   Eventually returned values will be:
            *   00 00 10 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00
            *   00 00 00 10 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00
            *      ....................................................................................................................
            *                                                                                                                     10 00
            *                                                                                                                        10
            *   The loop will add all the recursively genereated values with the original value will give unique combinations:
            *   00 00 10 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00
            *
            *   00 00 00 10 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00
            *
            *   00 00 10 10 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00
            *
            * TL:DR This loop only reduces the sequence length, in order to repeat the first loop process with a new sequence.
            */
            for (auto mask : computeMasksTwoBit(seqLength - 1, mismatches)) {
                masks.push_back(mask);
            }
        }
        else {
            // No more mismatches allowed, return a mask of 0 (Exit condition).
            masks.push_back(0ULL);
        }
    }
    /**
    * The only valid combination when mismatches == seqLength will be to assign remaining seq as mismatches.
    * E.g. mismatches=4 and seqLength=4, we want: 10 10 10 10
    * Do not recursivley call self as there are no more valid combinations from this point (Exit condition).
    */
    else {
        uint64_t tempMask = 0;
        for (int i = 0; i < seqLength; i++) {
            tempMask |= (1ULL << i * 2);
        }
        masks.push_back(tempMask);
    }
    return masks;
}

/**
 * calcMITLocalScore
 *
 * This function calculates the local MIT score based on the positions of mismatches.
 * https://dx.doi.org/10.1038/nbt.2647
 *
 * @param mismatch_array, An int array that contains the position of mismatches
 * @param length, The length of the param `mismatch_array`
 * @return The local MIT score
 */
double ISSLOffTargetScoring::calcMITLocalScore(int* mismatch_array, int length) {
    int i;
    double T1 = 1.0, T2, T3, d = 0.0, score;
    /* Mismatch penalty array */
    double M[] = { 0.0, 0.0, 0.014, 0.0, 0.0, 0.395, 0.317, 0.0, 0.389, 0.079, 0.445, 0.508, 0.613, 0.851, 0.732, 0.828, 0.615, 0.804, 0.685, 0.583 };

    /* 1st term */
    for (i = 0; i < length; ++i)
        T1 = T1 * (1.0 - M[mismatch_array[i]]);

    /* 2nd term */
    if (length == 1)
        d = 19.0;
    else {
        for (i = 0; i < length - 1; ++i)
            d += mismatch_array[i + 1] - mismatch_array[i];
        d = d / (length - 1);
    }
    T2 = 1.0 / ((19.0 - d) / 19.0 * 4.0 + 1);

    /* 3rd term */
    T3 = 1.0 / (length * length);

    /* Local score */
    score = T1 * T2 * T3 * 100;
    return score;
}

/**
 * predictMITLocalScore
 *
 * Calculates mistmatch postions before calling scoring equation.
 *
 * @param xoredSignatures, A mask that represents the mismatches between a potential target and off-target. (The result of the XOR operation in practice)
 * @return The local MIT score
 */
double ISSLOffTargetScoring::predictMITLocalScore(uint64_t xoredSignatures)
{
    int mismatch_array[20], m = 0;
    for (size_t j = 0; j < 20; j++) {
        if ((xoredSignatures >> (j * 2)) & 0x3) {
            mismatch_array[m++] = j;
        }
    }
    if (m == 0) return 0.0;
    return calcMITLocalScore(mismatch_array, m);
}


#if defined(_WIN64)
    #pragma pop_macro("close")
#endif