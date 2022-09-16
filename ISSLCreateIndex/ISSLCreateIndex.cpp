#include "ISSLCreateIndex.hpp"

using std::vector;
using std::map;
using std::pair;
using std::string;
using std::string_view;
using std::regex;
using std::regex_iterator;

const regex extractNumbers("[1234567890]+");
const vector<uint8_t> nucleotideIndex{ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,3 };
const vector<char> signatureIndex{ 'A', 'C', 'G', 'T' };
size_t seqLength;

/**
 * getFileSize
 *
 * @param path, The path to the file
 * @return The size of the file (in Bytes) at `path`
 */
size_t getFileSize(const char* path)
{
    struct p_stat64 statBuf;
    p_stat64(path, &statBuf);
    return statBuf.st_size;
}

/**
 * sequenceToSignature
 * 
 * Converts from sequence of letters to two bit encoded signature
 * 
 * @param ptr, A pointer to the beginning of a char seqeuence
 * @return The two bit encoded signature of the sequence at `ptr`
 */
uint64_t sequenceToSignature(const char* ptr)
{
    uint64_t signature = 0;
    for (size_t j = 0; j < seqLength; j++) {
        signature |= (uint64_t)(nucleotideIndex[*ptr]) << (j * 2);
        ptr++;
    }
    return signature;
}

/**
 * sequenceToSignature
 * 
 * Converts from two bit encoded signature to sequence of letters
 * 
 * @param signature, The two bit encoded signature
 * @return The letter sequence of `signature`
 */
string signatureToSequence(uint64_t signature)
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
vector<uint64_t> computeMasksTwoBit(int seqLength, int mismatches) {
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
double calcMITLocalScore(int* mismatch_array, int length) {
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
double predictMITLocalScore(uint64_t xoredSignatures)
{
    int mismatch_array[20], m = 0;
    for (size_t j = 0; j < seqLength; j++) {
        if ((xoredSignatures >> (j * 2)) & 0x3) {
            mismatch_array[m++] = j;
        }
    }
    if (m == 0) return 0.0;
    return calcMITLocalScore(mismatch_array, m);
}

int main(int argc, char** argv)
{
    // Check input args
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
        sliceRanges.push_back(stoi(i->str()) * 2);

    }
    if (sliceRanges.size() % 2) {
        fprintf(stderr, "Error: Uneven number of slice start and end points provided\n");
        fprintf(stderr, "Please format slice length list like [(sliceStart:sliceEnd)...(sliceStart:sliceEnd)]\n");
        exit(1);
    }
    if (sliceRanges.size() < 2) {
        fprintf(stderr, "Error: Please specific more than 2 slice lengths\n");
        fprintf(stderr, "Please format slice length list like [(sliceStart:sliceEnd)...(sliceStart:sliceEnd)]\n");
        exit(1); 
    }
    vector<size_t> sliceLens;
    for (int i = 0; i < sliceRanges.size(); i = i + 2) {
        sliceLens.push_back((sliceRanges[i + 1] + 2) - (sliceRanges[i]));
    }

    size_t seqCount = fileSize / seqLineLength;
    fprintf(stderr, "Number of sequences: %zu\n", seqCount);

    size_t globalCount = 0;

    vector<uint64_t> seqSignatures;
    vector<uint32_t> seqSignaturesOccurrences;

    size_t offtargetsCount = 0;
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
                    fprintf(stderr, "%zu/%zu : %zu\n", (progressCount + occurrences), seqCount, offtargetsCount);

            }

            seqSignatures.push_back(signature);
            seqSignaturesOccurrences.push_back(occurrences);

            offtargetsCount++;
            if (progressCount % 10000 == 0)
                fprintf(stderr, "%zu/%zu : %zu\n", progressCount, seqCount, offtargetsCount);

            progressCount += occurrences;
        }

    }
    printf("Finished counting occurrences, now constructing index...\n");

    vector<vector<vector<uint64_t>>> sliceLists(sliceLens.size());
    // Assign sliceLists size based on each slice length
    for (int i = 0; i < sliceLens.size(); i++)
    {
        sliceLists[i] = vector<vector<uint64_t>>(1ULL << sliceLens[i]);
    }

    // Generate ISSL Index
    #pragma omp parallel for
    for (int i = 0; i < sliceLens.size(); i++) {
        // sliceMask is a mask of 1's the length of the slice, used to extract target slice
        uint64_t sliceMask = (1 << sliceLens[i]) - 1;
        // sliceShift is the amout the target signature needs to be shifted to align with the target slice position 
        int sliceShift = sliceRanges[i*2];
        auto& sliceList = sliceLists[i];
        
        uint32_t signatureId = 0;
        for (uint64_t signature : seqSignatures) {
            uint32_t occurrences = seqSignaturesOccurrences[signatureId];
            // Shift the target signature to align with target slice position
            uint32_t sliceVal = (signature >> sliceShift) & sliceMask;
            // seqSigIdVal represnets the sequence signature ID and number of occurrences of the associated sequence.
            // (((uint64_t)occurrences) << 32), the most significant 32 bits is the count of the occurrences.
            // (uint64_t)signatureId, the index of the sequence in `seqSignatures`
            uint64_t seqSigIdVal = ((uint64_t)occurrences << 32) | (uint64_t)signatureId;
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
            double score = predictMITLocalScore(mask);
            precalculatedScores.insert(pair<uint64_t, double>(mask, score));
            scoresCount++;
        }
    }

    printf("Finished calculating scores, now preparing to write to disk...\n");

    fp = fopen(argv[4], "wb");
    // Generate Header list
    vector<size_t> slicelistHeader;
    slicelistHeader.push_back(offtargetsCount);
    slicelistHeader.push_back(seqLength);
    slicelistHeader.push_back(seqCount);
    slicelistHeader.push_back(sliceRanges.size());
    slicelistHeader.push_back(scoresCount);
   
    // write the header, reserve the first 50 for header information
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

    // write slice list lengths
    for (size_t i = 0; i < sliceLens.size(); i++) { // Number of slices
        for (size_t j = 0; j < (1ULL << sliceLens[i]); j++) { // Slice limit given slice width
            size_t sz = sliceLists[i][j].size();
            fwrite(&sz, sizeof(size_t), 1, fp);
        }
    }

    // write slice list data
    for (size_t i = 0; i < sliceLens.size(); i++) { // Number of slices
        for (size_t j = 0; j < (1ULL << sliceLens[i]); j++) { // Slice limit given slice width
            fwrite(sliceLists[i][j].data(), sizeof(uint64_t), sliceLists[i][j].size(), fp); // vector
        }
    }

    printf("Writing to disk...\n");

    fclose(fp);
    printf("Done.\n");
    return 0;
}
