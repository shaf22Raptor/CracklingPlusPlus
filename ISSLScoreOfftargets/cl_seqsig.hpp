#pragma once
#include <vector>
#include <cstdint>
#include <string>

class SeqSigEncoderCL {
public:
    SeqSigEncoderCL();  // builds context, queue, program, kernel, LUT buffer
    ~SeqSigEncoderCL();

    // Encodes queryCount sequences from queryBytes into outSigs.
    // queryStride = seqLen + 1 if you have '\n' at the end of each line.
    void encode(const char* queryBytes,
                uint32_t    queryStride,
                uint32_t    seqLen,
                uint32_t    queryCount,
                std::vector<uint64_t>& outSigs);

private:
    struct Impl;
    Impl* p_;
};
