#pragma once
#include <cstddef>
#include <cstdint>
#include <vector>
#include <string>
#include <cstdio>          // for std::FILE in the extern below
#include "ISSLScoreOfftargets.hpp"

// forward decls
struct AppOptions;               // <-- added
enum otScoreMethod : int;
enum class ClScorePrecision { Float32, Float64 };

#ifdef MANPROF_ENABLE
extern std::FILE* g_profileLogFile;
#endif

// Declarations only
ClScorePrecision cl_get_precision() noexcept;
std::size_t      cl_get_total_global_mem_bytes() noexcept;
std::size_t      cl_get_static_bytes() noexcept;

// GPU-scoring static view of the host-side data
struct ScoringIndexMeta {
    std::size_t offtargetsCount;
    std::size_t seqLength;
    std::size_t sliceCount;

    const std::uint64_t* pOfftargets;        // [offtargetsCount]
    const std::uint64_t* pAllSignatures;     // [sliceCount * offtargetsCount]

    const std::uint32_t* pSliceKeyCounts;       // [sliceKeyTableLen]
    const std::uint64_t* pSliceKeyBaseOffsets;  // [sliceKeyTableLen]
    std::size_t          sliceKeyTableLen;

    const std::uint32_t* pSliceMaskPositions;   // [sliceMaskPositionsLen]
    std::size_t          sliceMaskPositionsLen;
    const std::uint32_t* pSliceMaskOffsets;     // [sliceCount]
    const std::uint32_t* pSliceMaskLengths;     // [sliceCount]
};

// Step-2 API
void  cl_set_precision(ClScorePrecision p);
void  cl_set_mit_lut(const double* data, std::size_t len);
void  cl_set_cfd_pos_penalties(const double* data, std::size_t len);
void  cl_set_cfd_pam_penalties(const double* data, std::size_t len);

// Init/teardown + batch
bool init_scoring_cl(const ScoringIndexMeta& meta, const AppOptions& opts);
bool score_batch_cl(const std::uint64_t* querySigs,
                    std::size_t          guideCount,
                    double*              outMit,
                    double*              outCfd,
                    otScoreMethod        method,
                    double               threshold,
                    std::size_t          seqLength);

void shutdown_scoring_cl();
