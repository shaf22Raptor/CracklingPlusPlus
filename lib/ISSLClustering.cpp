#include "../include/ISSLClustering.hpp"
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

    for (int64_t i = 0; i < perGuideCountTotal.size(); i++)
    {
        truePerGuideCountTotal[perGuideCountTotal[true][i]]++;
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

    for (int64_t i = 0; i < perGuideCountUnique.size(); i++)
    {
        truePerGuideCountUnique[perGuideCountUnique[true][i]]++;
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

    for (int64_t i = 0; i < perGuideCountTotal.size(); i++)
    {
        uint64_t count = perGuideCountTotal[false][i];
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

    for (int64_t i = 0; i < perGuideCountUnique.size(); i++)
    {
        uint64_t count = perGuideCountUnique[false][i];
        falsePerGuideCountUnique[count]++;
    }

    outputFile.open(outputPath / "_output" / "falsePerGuideCountUnique.txt", std::ios::out | std::ios::binary);
    for (auto const& x : falsePerGuideCountUnique)
    {
        outputFile << fmt::format("{}:{}\n", x.first, x.second);
    }
    outputFile.close();
}

ISSLClustering::ISSLClustering(ConfigManager& cm) :
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

void ISSLClustering::run(unordered_map<string, unordered_map<string, string>>& candidateGuides)
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

        // Mutex
        std::mutex countsMutex;
        // Neighbourhood counts
        uint64_t neighbourhoodCount = 0;
        // Mismatch counts
        vector<uint64_t> offTargetCount(21, 0);
        // Guide counts
        vector<vector<uint64_t>> perGuideCount(2, vector<uint64_t>(querySignatures.size(), 0));

        /** OT scoring loop */
        omp_set_num_threads(threadCount);
        #pragma omp parallel
        {
            #pragma omp for
            for (int i = 0; i < queryCount; i++) {
                uint64_t searchSignature = querySignatures[i];

                /** Global scores */
                double totScoreMit = 0.0;
                double totScoreCfd = 0.0;

                int numOffTargetSitesScored = 0;
                double maximum_sum = (10000.0 - scoreThreshold * 100) / scoreThreshold;

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
                                totScoreMit += precalculatedMITScores.at(mismatches);
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

                        // TODO: remove
                        countsMutex.lock();
                        perGuideCount[true][i]++;
                        offTargetCount[dist]++;
                        countsMutex.unlock();
                    }
                    else
                    {
                        // TODO: remove
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
        for (long long falseOTCount : perGuideCount[false])
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
size_t ISSLClustering::getFileSize(const char* path)
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
uint64_t ISSLClustering::sequenceToSignature(const char* ptr, const size_t& seqLength)
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
string ISSLClustering::signatureToSequence(uint64_t signature, const size_t& seqLength)
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