// scoring_kernels.cl
// Portable OpenCL 1.2 skeleton for Crackling++ scoring
// - No vendor extensions except optional fp64
// - No program-scope pointers (LUTs passed as kernel args)
// - Safe bounds checks; outputs are zeroed for now

// -------- Optional double precision --------
// Define USE_FP64 at build time if the device supports fp64
#ifdef USE_FP64
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
typedef double ScoreT;
#else
typedef float  ScoreT;
#endif

// 64-bit integers (commonly available; include pragma for safety on old stacks)
#pragma OPENCL EXTENSION cl_khr_int64 : enable

// ---------- Bit helpers (match CPU logic exactly) ----------
inline ulong mismatch_mask(ulong ss, ulong ot)
{
    ulong x = ss ^ ot;
    const ulong EVEN_BITS = (ulong)0xAAAAAAAAAAAAAAAAUL; // ...10 10 pattern
    const ulong ODD_BITS  = (ulong)0x5555555555555555UL; // ...01 01 pattern
    // CPU parity: ((x & even) >> 1) | (x & odd)
    return ((x & EVEN_BITS) >> 1) | (x & ODD_BITS);
}

inline uint dist_from_mask(ulong mm)
{
    // popcount is built-in in OpenCL C
    return (uint)popcount(mm);
}

// Build the per-slice key from a packed guide signature and the mask positions.
//   searchSlice = concat( guide[maskPos[j]] as 2-bit symbols, in-order )
inline ulong build_search_slice(ulong guideSig,
                                __constant const uint* maskPositions,
                                uint maskLen)
{
    ulong acc = (ulong)0;
    // Note: maskLen is typically <= 20; this loop is tiny
    for (uint j = 0; j < maskLen; ++j) {
        uint pos2 = maskPositions[j] * 2u;
        ulong twoBits = (guideSig >> pos2) & (ulong)3u;
        acc |= (twoBits << (2u * j));
    }
    return acc;
}

// ============= Kernel arguments =============
//
// This kernel is a no-op scoring stub for now:
//  • It binds all expected buffers (index + LUTs).
//  • It computes guide index = get_global_id(0) and writes zeros.
//  • Later we’ll implement the full slice/key traversal & scoring.
//
// Global size: launch with at least guideCount work-items.
//
__kernel void score_kernel(
    // Guides (query signatures, 2-bit packed, 64-bit)
    __global const ulong* querySigs,
    ulong                 guideCount,

    // ---- Flattened index/meta ----
    // Off-target sequences (64-bit signatures)
    __global const ulong* offtargets,
    ulong                 offtargetsCount,

    // Concatenated signatures for all slices (layout matches host)
    __global const ulong* allSignatures,

    // Per (slice,key): counts and base offsets (flattened)
    __constant const uint*  sliceKeyCounts,       // length = sum( 4^(maskLen_i) )
    __constant const ulong* sliceKeyBaseOffsets,  // same length

    // Slice mask tables
    __constant const uint*  sliceMaskPositions,   // concatenated positions
    __constant const uint*  sliceMaskOffsets,     // per-slice offset into positions
    __constant const uint*  sliceMaskLengths,     // per-slice length
    uint                    sliceCount,

    // Scoring/LUTs (passed as buffers; constant-qualified)
    __constant const ScoreT* mitLut,
    __constant const ScoreT* cfdPos,
    __constant const ScoreT* cfdPam,

    // Scoring knobs
    int    method,        // enum otScoreMethod (passed as int)
    ScoreT threshold,     // early-exit threshold (same formula as CPU)
    uint   seqLength,     // usually 20

    // Outputs (one per guide)
    __global ScoreT* outMit,
    __global ScoreT* outCfd

    // (Future) per-guide counters/debug could go here
)
{
    ulong gid = (ulong)get_global_id(0);
    if (gid >= guideCount) return;

    // Read this guide’s packed signature
    ulong guide = querySigs[gid];

    // ---- Placeholder “no-op” behavior ----
    // For now, write zeros so the host path can validate plumbing.
    // Later we’ll:
    //   - loop over slices
    //   - build searchSlice via build_search_slice(...)
    //   - locate run via sliceKeyCounts/BaseOffsets
    //   - iterate signatures, dedup offtarget ids with a per-guide bitset
    //   - compute dist and accumulate MIT/CFD using LUTs
    //   - apply early-exit logic by 'method' and 'threshold'
    if (outMit) outMit[gid] = (ScoreT)0;
    if (outCfd) outCfd[gid] = (ScoreT)0;

    // (Optional sanity: touch the LUTs to ensure args wired correctly)
    // This keeps the compiler from optimizing args away in some drivers.
    // Remove once real scoring uses them.
    if ((gid == 0ul) && mitLut && cfdPos && cfdPam) {
        // volatile prevents DCE; write to outMit[0] harmlessly.
        volatile ScoreT sink = mitLut[0];
        sink += cfdPos[0];
        sink += cfdPam[0];
        if (outMit) outMit[0] = (ScoreT)0 + (ScoreT)0 * sink; // still zero
    }
}
