// cl_scoring.hpp
#pragma once
#include <cstddef>
#include <cstdint>
#include <vector>
#include <string>

// Bring in otScoreMethod so the signatures match your codebase
#include "ISSLScoreOfftargets.hpp"

enum otScoreMethod : int;
enum class ClScorePrecision { Float32, Float64 };

#ifdef MANPROF_ENABLE
extern std::FILE* g_profileLogFile;
#endif

#if USE_OPENCL_SCORING
ClScorePrecision cl_get_precision() noexcept;
std::size_t      cl_get_total_global_mem_bytes() noexcept;
std::size_t      cl_get_static_bytes() noexcept;
#else
inline ClScorePrecision cl_get_precision() noexcept { return ClScorePrecision::Float32; }
inline std::size_t      cl_get_total_global_mem_bytes() noexcept { return 0; }
inline std::size_t      cl_get_static_bytes() noexcept { return 0; }
#endif


// GPU-scoring static view of the host-side data
struct ScoringIndexMeta {
    size_t offtargetsCount;
    size_t seqLength;
    size_t sliceCount;

    const uint64_t* pOfftargets;         // [offtargetsCount]
    const uint64_t* pAllSignatures;      // [sliceCount * offtargetsCount]

    // flattened tables
    const uint32_t* pSliceKeyCounts;       // [sliceKeyTableLen]
    const uint64_t* pSliceKeyBaseOffsets;  // [sliceKeyTableLen]
    size_t          sliceKeyTableLen;

    const uint32_t* pSliceMaskPositions; // [sliceMaskPositionsLen]
    size_t          sliceMaskPositionsLen;
    const uint32_t* pSliceMaskOffsets;   // [sliceCount]
    const uint32_t* pSliceMaskLengths;   // [sliceCount]
};

// Step-2 API: choose precision and push LUTs
void  cl_set_precision(ClScorePrecision p);     // call once before init/build
void  cl_set_mit_lut(const double* data, std::size_t len); // MIT table lives in host as double; kernel casts if needed
void  cl_set_cfd_pos_penalties(const double* data, std::size_t len);
void  cl_set_cfd_pam_penalties(const double* data, std::size_t len);

// Init/teardown + batch entry points
void init_scoring_cl(const ScoringIndexMeta& meta);
bool score_batch_cl(const uint64_t* querySigs,
                    std::size_t     guideCount,
                    double*         outMit,     // nullable
                    double*         outCfd,     // nullable
                    otScoreMethod   method,     // your enum type is fine too
                    double          threshold,
                    std::size_t     seqLength);

void shutdown_scoring_cl();