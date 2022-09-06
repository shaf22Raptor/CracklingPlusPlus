#include "ISSLCreateIndex.hpp"

using std::vector;
using std::map;
using std::pair;
using std::string;
using std::string_view;
using std::regex;
using std::regex_iterator;

const regex extractNumbers("[1234567890]+");

size_t seqLength;
vector<uint8_t> nucleotideIndex(256);
vector<char> signatureIndex(4);

size_t getFileSize(const char* path)
{
    struct p_stat64 statBuf;
    p_stat64(path, &statBuf);
    return statBuf.st_size;
}

uint64_t sequenceToSignature(const char* ptr)
{
    uint64_t signature = 0;
    for (size_t j = 0; j < seqLength; j++) {
        signature |= (uint64_t)(nucleotideIndex[*ptr]) << (j * 2);
        ptr++;
    }
    return signature;
}

string signatureToSequence(uint64_t signature)
{
    string sequence = string(seqLength, ' ');
    for (size_t j = 0; j < seqLength; j++) {
        sequence[j] = signatureIndex[(signature >> (j * 2)) & 0x3];
    }
    return sequence;
}

// combinations
vector<uint64_t> computeMasksTwoBit(int seqLength, int mismatches) {
    vector<uint64_t> masks;

    // there are more positions than mismatches
    if (mismatches < seqLength) {

        if (mismatches > 0) {
            // some mismatches across a long sequence
            for (auto mask : computeMasksTwoBit(seqLength - 1, mismatches - 1)) {
                // orig: masks.push_back((1 << (seqLength - 1)) + mask);
                masks.push_back((1LLU << (seqLength - 1) * 2) + mask);
            }
            for (auto mask : computeMasksTwoBit(seqLength - 1, mismatches)) {
                masks.push_back(mask);
            }
        }
        else {
            // no mismatches at all
            masks.push_back(0LLU);
        }

        // mismatches >= seqLength
        // every position is going to be a mismatch
        // eg: mismatches=4 and seqLength=4, we want: 10 10 10 10
    }
    else {
        // orig: masks.push_back((1LLU << seqLength) - 1);
        uint64_t tempMask = 0;
        for (int i = 0; i < seqLength; i++) {
            tempMask |= (1LLU << i * 2);
        }
        masks.push_back(tempMask);
    }
    return masks;
}

double single_score(int* mismatch_array, int length) {
    int i;
    double T1 = 1.0, T2, T3, d = 0.0, score;
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

    /* Total score */
    score = T1 * T2 * T3 * 100;
    return score;
}

double sscore(uint64_t xoredSignatures)
{
    int mismatch_array[20], m = 0;
    for (size_t j = 0; j < seqLength; j++) {
        if ((xoredSignatures >> (j * 2)) & 0x3) {
            mismatch_array[m++] = j;
        }
    }
    if (m == 0) return 0.0;
    return single_score(mismatch_array, m);
}

int main(int argc, char** argv)
{
    if (argc < 5) {
        fprintf(stderr, "Usage: %s [offtargetSites.txt] [sequence length] [slice width (bits)] [sissltable]\n", argv[0]);
        exit(1);
    }
    size_t fileSize = getFileSize(argv[1]);

    FILE* fp = fopen(argv[1], "rb");
    seqLength = atoi(argv[2]);
    if (seqLength > 32) {
        fprintf(stderr, "Sequence length is greater than 32, which is the maximum supported currently\n");
        exit(1);
    }
    size_t seqLineLength = seqLength + 1; // '\n'
    if (fileSize % seqLineLength != 0) {
        fprintf(stderr, "fileSize: %zu\n", fileSize);
        fprintf(stderr, "Error: file does is not a multiple of the expected line length (%zu)\n", seqLineLength);
        fprintf(stderr, "The sequence length may be incorrect; alternatively, the line endings\n");
        fprintf(stderr, "may be something other than LF, or there may be junk at the end of the file.\n");
        exit(1);
    }
    const string sliceArg = argv[3];
    vector<size_t> sliceRanges;
    for (auto i  = std::sregex_iterator(sliceArg.begin(), sliceArg.end(), extractNumbers);
        i != std::sregex_iterator();
        i++) 
    {
        sliceRanges.push_back(stoi(i->str())*2);
    }
    if (sliceRanges.size() % 2) {
        fprintf(stderr, "Error: Uneven number of slice start and end points provided\n");
        fprintf(stderr, "Please format slice length list like [(start:end)...(start:end)]\n");
        exit(1);
    }
    if (sliceRanges.size() < 2) {
        fprintf(stderr, "Error: Please specific more than 2 slice lengths\n");
        fprintf(stderr, "Please format slice length list like [(start:end)...(start:end)]\n");
        exit(1);
    }
    vector<size_t> sliceLens;
    for (int i = 0; i < sliceRanges.size(); i = i + 2) {
        sliceLens.push_back(sliceRanges[i + 1] - (sliceRanges[i] - 2));
    }

    size_t seqCount = fileSize / seqLineLength;
    fprintf(stderr, "Number of sequences: %zu\n", seqCount);


    nucleotideIndex['A'] = 0;
    nucleotideIndex['C'] = 1;
    nucleotideIndex['G'] = 2;
    nucleotideIndex['T'] = 3;
    signatureIndex[0] = 'A';
    signatureIndex[1] = 'C';
    signatureIndex[2] = 'G';
    signatureIndex[3] = 'T';

    size_t globalCount = 0;

    vector<uint64_t> seqSignatures;
    vector<uint32_t> seqSignaturesOccurrences;

    size_t distinctSites = 0;
    {
        vector<char> entireDataSet(fileSize);

        if (fread(entireDataSet.data(), fileSize, 1, fp) < 1) {
            fprintf(stderr, "Failed to read in file.\n");
            exit(1);
        }
        fclose(fp);

        size_t progressCount = 0;
        size_t offtargetId = 0;
        while (progressCount < seqCount) {
            char* ptr = &entireDataSet[progressCount * seqLineLength];

            uint64_t signature = sequenceToSignature(ptr);

            // check how many times the off-target appears
            // (assumed the list is sorted)
            uint32_t occurrences = 1;
            while (memcmp(ptr, ptr + (seqLineLength * occurrences), seqLength) == 0) {
                occurrences++;
                if ((seqCount - progressCount - occurrences) < 100)
                    fprintf(stderr, "%zu/%zu : %zu\n", (progressCount + occurrences), seqCount, distinctSites);

            }

            seqSignatures.push_back(signature);
            seqSignaturesOccurrences.push_back(occurrences);

            distinctSites++;
            if (progressCount % 10000 == 0)
                fprintf(stderr, "%zu/%zu : %zu\n", progressCount, seqCount, distinctSites);

            progressCount += occurrences;
        }

    }
    printf("Finished counting occurrences, now constructing index...\n");
    size_t offtargetsCount = distinctSites;

    vector<vector<vector<uint64_t>>> sliceLists(sliceLens.size());
    for (int i = 0; i < sliceLens.size(); i++)
    {
        sliceLists[i] = vector<vector<uint64_t>>(1 << sliceLens[i]);
    }
    /*#pragma omp parallel for*/
    for (int i = 0; i < sliceLens.size(); i++) {
        uint64_t sliceMask = (1 << sliceLens[i]) - 1;
        int sliceShift = sliceRanges[i * 2] - 2;
        sliceMask = sliceMask << sliceShift;
        auto& sliceList = sliceLists[i];

        uint32_t signatureId = 0;
        for (uint64_t signature : seqSignatures) {
            uint32_t occurrences = seqSignaturesOccurrences[signatureId];
            uint64_t test = signature & sliceMask;
            uint8_t sliceVal = (signature & sliceMask) >> sliceShift;

            uint64_t seqSigIdVal = (((uint64_t)occurrences) << 32) | (uint64_t)signatureId;
            sliceList[sliceVal].push_back(seqSigIdVal);
            signatureId++;
        }
    }

    printf("Finished constructing index, now precalculating scores...\n");

    // Precalculate all the scores
    map<uint64_t, double> precalculatedScores;

    int maxDist = 4;
    size_t scoresCount = 0;

    for (int i = 1; i <= maxDist; i++) {
        vector<uint64_t> tempMasks;
        tempMasks = computeMasksTwoBit(20, i);
        for (auto mask : tempMasks) {
            double score = sscore(mask);
            precalculatedScores.insert(pair<uint64_t, double>(mask, score));
            scoresCount++;
        }
    }

    printf("Finished calculating scores, now preparing to write to disk...\n");

    fp = fopen(argv[4], "wb");
    vector<size_t> slicelistHeader;
    slicelistHeader.push_back(offtargetsCount);
    slicelistHeader.push_back(seqLength);
    slicelistHeader.push_back(seqCount);
    slicelistHeader.push_back(sliceLens.size());
    slicelistHeader.push_back(scoresCount);
   
    // write the header
    fwrite(slicelistHeader.data(), sizeof(size_t), 50, fp);

    // write slice ranges
    for (const size_t pos : sliceRanges)
    {
        fwrite(&pos, sizeof(size_t), 1, fp);
    }

    // write the precalculated scores
    for (auto const& x : precalculatedScores) {
        fwrite(&x.first, sizeof(uint64_t), 1, fp);
        fwrite(&x.second, sizeof(double), 1, fp);
    }

    // write the offtargets
    fwrite(seqSignatures.data(), sizeof(uint64_t), seqSignatures.size(), fp);

    for (size_t i = 0; i < sliceLens.size(); i++) { // Number of slices
        for (size_t j = 0; j < (1 << sliceLens[i]); j++) { // Slice limit give slice width
            size_t sz = sliceLists[i][j].size();
            fwrite(&sz, sizeof(size_t), 1, fp);
        }
    }

    for (size_t i = 0; i < sliceLens.size(); i++) { // Number of slices
        for (size_t j = 0; j < (1 << sliceLens[i]); j++) { // Slice limit give slice width
            fwrite(sliceLists[i][j].data(), sizeof(uint64_t), sliceLists[i][j].size(), fp); // vector
        }
    }

    printf("Writing to disk...\n");

    fclose(fp);
    printf("Done.\n");
    return 0;
}
