// // scoring_kernels.cl
// // Portable OpenCL 1.2 skeleton for Crackling++ scoring
// // - No vendor extensions except optional fp64
// // - No program-scope pointers (LUTs passed as kernel args)
// // - Safe bounds checks; outputs are zeroed for now

// // -------- Optional double precision --------
// Define USE_FP64 at build time if the device supports fp64
#ifdef USE_FP64
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
typedef double ScoreT;
#else
typedef float  ScoreT;
#endif

// // 64-bit integers (commonly available; include pragma for safety on old stacks)
// #pragma OPENCL EXTENSION cl_khr_int64 : enable

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

inline ulong pack_mismatch_bits(ulong mm, uint seqLength) {
    ulong packed = 0ul;
    for (uint pos = 0; pos < seqLength; ++pos)
        packed |= ((mm >> (2u*pos)) & 1ul) << pos;
    return packed;
}

__kernel void score_kernel(
    __global const ulong* querySigs,
    ulong                 guideCount,

    __global const ulong* offtargets,
    ulong                 offtargetsCount,

    __global const ulong* allSignatures,

    __constant const uint*  sliceKeyCounts,
    __constant const ulong* sliceKeyBaseOffsets,

    __constant const uint*  sliceMaskPositions,
    __constant const uint*  sliceMaskOffsets,
    __constant const uint*  sliceMaskLengths,
    uint                    sliceCount,

    __constant const ScoreT* mitLut,
    __constant const ScoreT* cfdPos,
    __constant const ScoreT* cfdPam,

    int    method,      // enum otScoreMethod as int
    ScoreT threshold,
    uint   seqLength,

    __global ulong*      seenBits,      // guideCount * wordsPerGuide
    ulong                wordsPerGuide, // 64-bit words per guide

    __global ScoreT* outMit,
    __global ScoreT* outCfd
)
{
    ulong gid = (ulong)get_global_id(0);
    if (gid >= guideCount) return;

    ulong guideSig = querySigs[gid];

    // Per-guide bitset window
    __global ulong* seen = seenBits + gid * wordsPerGuide;

    // Accumulators
    ScoreT totMit = (ScoreT)0;
    ScoreT totCfd = (ScoreT)0;

    // Flags for which scores to compute (derive from whether out* was provided)
    const bool haveMitLut = (mitLut != 0);
    const bool calcMit = (outMit != 0) && haveMitLut;
    const bool calcCfd = (outCfd != 0) && (cfdPos != 0) && (cfdPam != 0);

    // Early-exit ceiling (matches CPU)
    const ScoreT maxSum = (((ScoreT)10000) - threshold * (ScoreT)100) / fmax((ScoreT)1e-9, threshold);

    // Walk the concatenated (slice,key) tables with a running index
    uint tbl = 0ul; // start of this slice’s key block in sliceKey{Counts,BaseOffsets}

    for (uint si = 0; si < sliceCount; ++si)
    {
        const uint  maskLen  = sliceMaskLengths[si];
        const uint  posOff   = sliceMaskOffsets[si];

        const ulong keys_i   = (1ul << (2u * maskLen));     // 64-bit!
        const ulong countsStart = tbl;
        const ulong countsEnd   = tbl + keys_i;

        const ulong searchSlice = build_search_slice(guideSig,
                                                    sliceMaskPositions + posOff,
                                                    maskLen);
        const ulong keyIdx   = countsStart + searchSlice;

        const uint  runCount = sliceKeyCounts[keyIdx];       // arrays unchanged
        const ulong base     = sliceKeyBaseOffsets[keyIdx];

        // Iterate signatures for this key
        for (uint j = 0; j < runCount; ++j) {
            const ulong packed = allSignatures[base + (ulong)j];
            const uint  id  = (uint)(packed & 0xFFFFFFFFul);
            const uint  occ = (uint)(packed >> 32);

            // Dedup per guide
            const ulong word = ((ulong)id) >> 6;
            const ulong bit  = (ulong)1 << (id & 63u);
            const ulong oldv = seen[word];
            if ((oldv & bit) != 0ul) continue;     // already seen
            seen[word] = oldv | bit;               // mark seen

            // Hamming distance via your parity trick
            const ulong mm   = mismatch_mask(guideSig, offtargets[id]);
            const uint  dist = dist_from_mask(mm);

            if (dist <= 4u) {
                // MIT
                if (calcMit && dist > 0u) {
                    // NOTE: mm can be up to 2*seqLength bits; mitLut was provisioned to 2^seqLength on host.
                    // Guard array access by masking if needed (optional if host LUT is dense enough).
                    const ulong packed = pack_mismatch_bits(mm, seqLength);   // <-- fix
                    totMit += (ScoreT)mitLut[packed] * (ScoreT)occ;
                }

                // CFD
                if (calcCfd) {
                    ScoreT cfd = (dist == 0u) ? (ScoreT)1 : (ScoreT)cfdPam[0b1010];
                    if (dist != 0u) {
                        const ulong ot = offtargets[id];
                        for (uint pos = 0; pos < seqLength; ++pos) {
                            const uint g2 = (uint)((guideSig >> (2u*pos)) & 3ul);
                            const uint o2 = (uint)((ot       >> (2u*pos)) & 3ul);
                            if (g2 != o2) {
                                const uint m = (pos << 4) | (g2 << 2) | (o2 ^ 3u);
                                // optional debug bound: if (m >= 16u*seqLength) continue;
                                cfd *= cfdPos[m];
                            }
                        }
                    }
                    totCfd += cfd * (ScoreT)occ;
                }

                // Early exit checks (match CPU semantics)
                if (threshold > (ScoreT)0) {
                    if (method == /* mitAndCfd */ 0 && (totMit > maxSum) && (totCfd > maxSum)) break;
                    if (method == /* mitOrCfd  */ 1 && ((totMit > maxSum) || (totCfd > maxSum))) break;
                    if (method == /* avgMitCfd */ 2 && (((totMit + totCfd) * (ScoreT)0.5) > maxSum)) break;
                    if (method == /* mit       */ 3 && (totMit > maxSum)) break;
                    if (method == /* cfd       */ 4 && (totCfd > maxSum)) break;
                }
            }
        }

        // advance to next slice block in the tables
        tbl = countsEnd;

        // if an early-exit was triggered in the inner loop, stop slices
        if (threshold > (ScoreT)0) {
            bool stop = false;
            if (method == 0 && (totMit > maxSum) && (totCfd > maxSum)) stop = true;
            if (method == 1 && ((totMit > maxSum) || (totCfd > maxSum))) stop = true;
            if (method == 2 && (((totMit + totCfd) * (ScoreT)0.5) > maxSum)) stop = true;
            if (method == 3 && (totMit > maxSum)) stop = true;
            if (method == 4 && (totCfd > maxSum)) stop = true;
            if (stop) break;
        }
    }

    // Final transform (same as CPU)
    if (outMit) outMit[gid] = calcMit ? ((ScoreT)10000 / ((ScoreT)100 + totMit)) : (ScoreT)(-1);
    if (outCfd) outCfd[gid] = calcCfd ? ((ScoreT)10000 / ((ScoreT)100 + totCfd)) : (ScoreT)(-1);
}

