// cl_scoring.hpp
#pragma once
#include <cstddef>
#include <cstdint>

// Bring in otScoreMethod so the signatures match your codebase
#include "ISSLScoreOfftargets.hpp"

// Host-side metadata we’ll pass to the scoring backend.
// For now, only the fields you already have on hand.
struct ScoringIndexMeta {
    size_t offtargetsCount = 0;
    size_t seqLength       = 0;
    size_t sliceCount      = 0;

    const uint64_t* pOfftargets           = nullptr; // offtargets.data()
    const uint64_t* pAllSignatures        = nullptr; // allSignatures.data()
    const size_t*   pAllSliceListSizes    = nullptr; // allSlicelistSizes.data()

    // You’ll fill these later when you flatten masks / prefix sums:
    const uint64_t* pSliceMaskPositions   = nullptr;
    const uint32_t* pSliceMaskOffsets     = nullptr;
    const uint32_t* pSliceMaskLengths     = nullptr;
    const uint32_t* pSliceKeyCounts       = nullptr;
    const uint64_t* pSliceKeyBaseOffsets  = nullptr;
};

// Init once after you’ve loaded the index; stub does nothing (returns true).
bool init_scoring_cl(const ScoringIndexMeta&);

// Try to score a whole batch on GPU. For now, stub returns false so CPU path runs.
// Later you’ll use it to fill outMit/outCfd for [0..guideCount).
bool score_batch_cl(const uint64_t* querySigs,
                    size_t guideCount,
                    double* outMit,
                    double* outCfd,
                    otScoreMethod method,
                    double threshold,
                    size_t seqLength);

// Cleanup at program end (stub does nothing).
void shutdown_scoring_cl();
