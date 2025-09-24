// cl_scoring.cpp
#include "cl_scoring.hpp"

#if USE_OPENCL_SCORING

#include <CL/cl.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <vector>
#include <mutex>
#include <string>

#include <cstdint>
static_assert(sizeof(cl_ulong)   == 8, "Expected cl_ulong to be 64-bit");
static_assert(sizeof(uint64_t)   == 8, "Expected uint64_t to be 64-bit");
static_assert(sizeof(cl_double)  == 8, "Expected cl_double to be 64-bit");
static_assert(sizeof(cl_float)   == 4, "Expected cl_float to be 32-bit");
static_assert(sizeof(cl_uint)    == 4, "Expected cl_uint to be 32-bit");

// ====== lightweight error helpers ======
static inline void cl_die_if(cl_int err, const char* where) {
    if (err != CL_SUCCESS) {
        std::fprintf(stderr, "OpenCL error %d at %s\n", err, where);
        std::fflush(stderr);
        std::exit(1);
    }
}

// ====== global state (scoped to this TU) ======
static std::mutex g_mu;

static cl_platform_id   g_platform = nullptr;
static cl_device_id     g_device   = nullptr;
static cl_context       g_ctx      = nullptr;
static cl_command_queue g_q        = nullptr;
static cl_program       g_prog     = nullptr;
static cl_kernel        g_kernel   = nullptr;

// static (uploaded once) device buffers
static cl_mem d_offtargets           = nullptr;
static cl_mem d_allSignatures        = nullptr;
static cl_mem d_sliceKeyCounts       = nullptr;
static cl_mem d_sliceKeyBaseOffsets  = nullptr;
static cl_mem d_sliceMaskPositions   = nullptr;
static cl_mem d_sliceMaskOffsets     = nullptr;
static cl_mem d_sliceMaskLengths     = nullptr;
static cl_mem d_mitLut               = nullptr;
static cl_mem d_cfdPos               = nullptr;
static cl_mem d_cfdPam               = nullptr;

// host-side precision/LUTs captured via setters (Step 2)
static std::vector<double> g_mitLutHost;
static std::vector<double> g_cfdPosHost;
static std::vector<double> g_cfdPamHost;
static ClScorePrecision g_prec = ClScorePrecision::Float32;


static size_t g_offtargetsCount_cached = 0;
static size_t g_sliceCount_cached = 0;
static size_t g_seqLength_cached = 0;

static size_t g_total_global_mem_bytes = 0; // device property
static size_t g_static_bytes           = 0; // what we've uploaded (sum)

// ====== setters you already declared in the header ======
void cl_set_precision(ClScorePrecision p) {
    std::lock_guard<std::mutex> lk(g_mu);
    g_prec = p;
}
void cl_set_mit_lut(const double* data, std::size_t len) {
    std::lock_guard<std::mutex> lk(g_mu);
    g_mitLutHost.assign(data, data + len);
}
void cl_set_cfd_pos_penalties(const double* data, std::size_t len) {
    std::lock_guard<std::mutex> lk(g_mu);
    g_cfdPosHost.assign(data, data + len);
}
void cl_set_cfd_pam_penalties(const double* data, std::size_t len) {
    std::lock_guard<std::mutex> lk(g_mu);
    g_cfdPamHost.assign(data, data + len);
}

// small util to read kernel source from file
static std::string read_text_file(const char* path) {
    std::ifstream ifs(path, std::ios::binary);
    if (!ifs) {
        std::fprintf(stderr, "Fatal: could not open kernel file: %s\n", path);
        std::exit(1);
    }
    std::ostringstream oss;
    oss << ifs.rdbuf();
    return oss.str();
}

// convert a double-vector to float-vector, no checks
static std::vector<float> to_float(const std::vector<double>& v) {
    std::vector<float> out(v.size());
    for (size_t i = 0; i < v.size(); ++i) out[i] = static_cast<float>(v[i]);
    return out;
}

// ====== init: choose device, build program, upload static data ======
void init_scoring_cl(const ScoringIndexMeta& meta)
{
    std::lock_guard<std::mutex> lk(g_mu);
    if (g_ctx) return;

    g_offtargetsCount_cached = meta.offtargetsCount;
    g_sliceCount_cached      = meta.sliceCount;
    g_seqLength_cached       = meta.seqLength;

    cl_int err = CL_SUCCESS;

    // 1) pick first platform/device (portable & simple)
    cl_uint nplat = 0;
    cl_die_if(clGetPlatformIDs(0, nullptr, &nplat), "clGetPlatformIDs(count)");
    std::vector<cl_platform_id> plats(nplat);
    cl_die_if(clGetPlatformIDs(nplat, plats.data(), nullptr), "clGetPlatformIDs(list)");
    g_platform = plats[0];

    cl_uint ndev = 0;
    cl_die_if(clGetDeviceIDs(g_platform, CL_DEVICE_TYPE_GPU, 0, nullptr, &ndev), "clGetDeviceIDs(count)");
    std::vector<cl_device_id> devs(ndev);
    cl_die_if(clGetDeviceIDs(g_platform, CL_DEVICE_TYPE_GPU, ndev, devs.data(), nullptr), "clGetDeviceIDs(list)");
    g_device = devs[0];

    cl_ulong memBytes = 0;
    cl_die_if(clGetDeviceInfo(g_device, CL_DEVICE_GLOBAL_MEM_SIZE,
                            sizeof(memBytes), &memBytes, nullptr),
                            "clGetDeviceInfo(GLOBAL_MEM_SIZE)");
    g_total_global_mem_bytes = static_cast<size_t>(memBytes);

    // 2) context + queue
#if defined(CL_VERSION_2_0)
    g_ctx = clCreateContext(nullptr, 1, &g_device, nullptr, nullptr, &err);
    cl_die_if(err, "clCreateContext");
    // create a CL 2.0 queue if available, else fall back below
    cl_queue_properties props[] = { 0 };
    g_q = clCreateCommandQueueWithProperties(g_ctx, g_device, props, &err);
    if (err != CL_SUCCESS) {
        g_q = clCreateCommandQueue(g_ctx, g_device, 0, &err);
        cl_die_if(err, "clCreateCommandQueue(fallback)");
    }
#else
    g_ctx = clCreateContext(nullptr, 1, &g_device, nullptr, nullptr, &err);
    cl_die_if(err, "clCreateContext");
    g_q = clCreateCommandQueue(g_ctx, g_device, 0, &err);
    cl_die_if(err, "clCreateCommandQueue");
#endif

    // 3) build program
    const char* kpath = "scoring_kernels.cl"; // ensure it’s next to the exe (see CMake copy step below)
    std::string src = read_text_file(kpath);
    const char* srcPtr = src.c_str();
    size_t      srcLen = src.size();
    g_prog = clCreateProgramWithSource(g_ctx, 1, &srcPtr, &srcLen, &err);
    cl_die_if(err, "clCreateProgramWithSource");

    // Build options: enable fp64 macro if requested
    std::string opts;
    if (g_prec == ClScorePrecision::Float64) {
        // We only *request* fp64 in the kernel via -D USE_FP64.
        // If device lacks cl_khr_fp64, build may fail; you can catch & fallback.
        opts += "-D USE_FP64=1";
    }
    err = clBuildProgram(g_prog, 1, &g_device, opts.c_str(), nullptr, nullptr);
    if (err != CL_SUCCESS) {
        // dump build log to help debugging
        size_t logSz = 0;
        clGetProgramBuildInfo(g_prog, g_device, CL_PROGRAM_BUILD_LOG, 0, nullptr, &logSz);
        std::vector<char> log(logSz + 1, 0);
        clGetProgramBuildInfo(g_prog, g_device, CL_PROGRAM_BUILD_LOG, logSz, log.data(), nullptr);
        std::fprintf(stderr, "Kernel build failed:\n%s\n", log.data());
        std::exit(1);
    }

    g_kernel = clCreateKernel(g_prog, "score_kernel", &err);
    cl_die_if(err, "clCreateKernel(score_kernel)");

    // 4) upload static buffers (index + masks + LUTs)
    //    All sizes are in elements; multiply by sizeof(T) for bytes.

    // offtargets (ulong)
    d_offtargets = clCreateBuffer(g_ctx, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                  meta.offtargetsCount * sizeof(cl_ulong),
                                  const_cast<cl_ulong*>(reinterpret_cast<const cl_ulong*>(meta.pOfftargets)), &err);
    cl_die_if(err, "clCreateBuffer(d_offtargets)");

    // allSignatures (ulong)
    // size is sliceCount * offtargetsCount
    d_allSignatures = clCreateBuffer(g_ctx, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                     (meta.sliceCount * meta.offtargetsCount) * sizeof(cl_ulong),
                                     const_cast<cl_ulong*>(reinterpret_cast<const cl_ulong*>(meta.pAllSignatures)), &err);
    cl_die_if(err, "clCreateBuffer(d_allSignatures)");

    // per (slice,key) tables
    // counts (uint), bases (ulong)
    // NOTE: host must provide the concatenated arrays via meta pointers.
    d_sliceKeyCounts = clCreateBuffer(g_ctx, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                      meta.sliceKeyTableLen * sizeof(cl_uint),
                                      const_cast<cl_uint*>(reinterpret_cast<const cl_uint*>(meta.pSliceKeyCounts)), &err);
    cl_die_if(err, "clCreateBuffer(d_sliceKeyCounts)");

    d_sliceKeyBaseOffsets = clCreateBuffer(g_ctx, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                           meta.sliceKeyTableLen * sizeof(cl_ulong),
                                           const_cast<cl_ulong*>(reinterpret_cast<const cl_ulong*>(meta.pSliceKeyBaseOffsets)), &err);
    cl_die_if(err, "clCreateBuffer(d_sliceKeyBaseOffsets)");

    // mask tables: positions (concat), offsets, lengths
    d_sliceMaskPositions = clCreateBuffer(g_ctx, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                          meta.sliceMaskPositionsLen * sizeof(cl_uint),
                                          const_cast<cl_uint*>(reinterpret_cast<const cl_uint*>(meta.pSliceMaskPositions)), &err);
    cl_die_if(err, "clCreateBuffer(d_sliceMaskPositions)");

    d_sliceMaskOffsets = clCreateBuffer(g_ctx, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                        meta.sliceCount * sizeof(cl_uint),
                                        const_cast<cl_uint*>(reinterpret_cast<const cl_uint*>(meta.pSliceMaskOffsets)), &err);
    cl_die_if(err, "clCreateBuffer(d_sliceMaskOffsets)");

    d_sliceMaskLengths = clCreateBuffer(g_ctx, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                        meta.sliceCount * sizeof(cl_uint),
                                        const_cast<cl_uint*>(reinterpret_cast<const cl_uint*>(meta.pSliceMaskLengths)), &err);
    cl_die_if(err, "clCreateBuffer(d_sliceMaskLengths)");

    // LUTs: convert to float if needed
    if (g_mitLutHost.empty() || g_cfdPosHost.empty() || g_cfdPamHost.empty()) {
        std::fprintf(stderr, "Fatal: LUTs not set before init_scoring_cl\n");
        std::exit(1);
    }

    if (g_prec == ClScorePrecision::Float32) {
        auto mitF = to_float(g_mitLutHost);
        auto posF = to_float(g_cfdPosHost);
        auto pamF = to_float(g_cfdPamHost);

        d_mitLut = clCreateBuffer(g_ctx, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                  mitF.size() * sizeof(cl_float), mitF.data(), &err);
        cl_die_if(err, "clCreateBuffer(d_mitLut)");

        d_cfdPos = clCreateBuffer(g_ctx, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                  posF.size() * sizeof(cl_float), posF.data(), &err);
        cl_die_if(err, "clCreateBuffer(d_cfdPos)");

        d_cfdPam = clCreateBuffer(g_ctx, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                  pamF.size() * sizeof(cl_float), pamF.data(), &err);
        cl_die_if(err, "clCreateBuffer(d_cfdPam)");
    } else {
        d_mitLut = clCreateBuffer(g_ctx, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                  g_mitLutHost.size() * sizeof(cl_double),
                                  const_cast<double*>(g_mitLutHost.data()), &err);
        cl_die_if(err, "clCreateBuffer(d_mitLut)");

        d_cfdPos = clCreateBuffer(g_ctx, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                  g_cfdPosHost.size() * sizeof(cl_double),
                                  const_cast<double*>(g_cfdPosHost.data()), &err);
        cl_die_if(err, "clCreateBuffer(d_cfdPos)");

        d_cfdPam = clCreateBuffer(g_ctx, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                  g_cfdPamHost.size() * sizeof(cl_double),
                                  const_cast<double*>(g_cfdPamHost.data()), &err);
        cl_die_if(err, "clCreateBuffer(d_cfdPam)");
    }
    // <<< put accounting block HERE >>>
    g_static_bytes = 0;
    auto add_bytes = [](size_t& acc, size_t bytes){ acc += bytes; };

    add_bytes(g_static_bytes, meta.offtargetsCount * sizeof(cl_ulong));
    add_bytes(g_static_bytes, meta.sliceCount * meta.offtargetsCount * sizeof(cl_ulong));
    add_bytes(g_static_bytes, meta.sliceKeyTableLen * sizeof(cl_uint));
    add_bytes(g_static_bytes, meta.sliceKeyTableLen * sizeof(cl_ulong));
    add_bytes(g_static_bytes, meta.sliceMaskPositionsLen * sizeof(cl_uint));
    add_bytes(g_static_bytes, meta.sliceCount * sizeof(cl_uint)); // offsets
    add_bytes(g_static_bytes, meta.sliceCount * sizeof(cl_uint)); // lengths

    const bool useF32 = (g_prec == ClScorePrecision::Float32);
    if (useF32) {
        add_bytes(g_static_bytes, g_mitLutHost.size() * sizeof(cl_float));
        add_bytes(g_static_bytes, g_cfdPosHost.size() * sizeof(cl_float));
        add_bytes(g_static_bytes, g_cfdPamHost.size() * sizeof(cl_float));
    } else {
        add_bytes(g_static_bytes, g_mitLutHost.size() * sizeof(cl_double));
        add_bytes(g_static_bytes, g_cfdPosHost.size() * sizeof(cl_double));
        add_bytes(g_static_bytes, g_cfdPamHost.size() * sizeof(cl_double));
    }

    std::fprintf(stderr, "[CL] static_upload_bytes=%zu (of %zu VRAM)\n",
                g_static_bytes, g_total_global_mem_bytes);
}

// ====== run one batch (stub kernel writes zeros) ======
bool score_batch_cl(const uint64_t* querySigs,
                    std::size_t     guideCount,
                    double*         outMit,
                    double*         outCfd,
                    otScoreMethod   method,
                    double          threshold,
                    std::size_t     seqLength)
                    
{
    if (!g_offtargetsCount_cached || !g_sliceCount_cached) return false;
    std::lock_guard<std::mutex> lk(g_mu);
    if (!g_ctx || !g_kernel) return false; // not initialized → caller will use CPU

    cl_int err = CL_SUCCESS;

    // Per-call device buffers (simple & safe for now)
    cl_mem d_query = clCreateBuffer(g_ctx, CL_MEM_READ_ONLY  | CL_MEM_COPY_HOST_PTR,
                                    guideCount * sizeof(cl_ulong),
                                    const_cast<cl_ulong*>(reinterpret_cast<const cl_ulong*>(querySigs)), &err);
    cl_die_if(err, "clCreateBuffer(d_query)");

    // outputs match kernel ScoreT; convert back to double on host if needed
    const bool useF32 = (g_prec == ClScorePrecision::Float32);
    cl_mem d_outMit = nullptr, d_outCfd = nullptr;

    if (useF32) {
        d_outMit = clCreateBuffer(g_ctx, CL_MEM_WRITE_ONLY, guideCount * sizeof(cl_float), nullptr, &err);
        cl_die_if(err, "clCreateBuffer(d_outMit)");
        d_outCfd = clCreateBuffer(g_ctx, CL_MEM_WRITE_ONLY, guideCount * sizeof(cl_float), nullptr, &err);
        cl_die_if(err, "clCreateBuffer(d_outCfd)");
    } else {
        d_outMit = clCreateBuffer(g_ctx, CL_MEM_WRITE_ONLY, guideCount * sizeof(cl_double), nullptr, &err);
        cl_die_if(err, "clCreateBuffer(d_outMit)");
        d_outCfd = clCreateBuffer(g_ctx, CL_MEM_WRITE_ONLY, guideCount * sizeof(cl_double), nullptr, &err);
        cl_die_if(err, "clCreateBuffer(d_outCfd)");
    }

    // ====== set kernel args (order must match the .cl) ======
    cl_uint arg = 0;
    cl_die_if(clSetKernelArg(g_kernel, arg++, sizeof(cl_mem), &d_query), "set arg querySigs");
    cl_ulong guideCountUL = static_cast<cl_ulong>(guideCount);
    cl_die_if(clSetKernelArg(g_kernel, arg++, sizeof(cl_ulong), &guideCountUL), "set arg guideCount");

    cl_die_if(clSetKernelArg(g_kernel, arg++, sizeof(cl_mem), &d_offtargets), "set arg offtargets");
    // we need offtargetsCount and allSignatures too; we cached them in init via ScoringIndexMeta lengths
    // easiest: store them in static variables when you init (add to globals as needed)
    // For now, pass via the counts we can infer: we don't have them cached, so add them to the header if needed.
    // Assuming you added: g_offtargetsCount, g_sliceCount (store at init time)
 
    cl_ulong offtCountUL = (cl_ulong)g_offtargetsCount_cached;

    cl_die_if(clSetKernelArg(g_kernel, arg++, sizeof(cl_ulong), &offtCountUL), "set arg offtargetsCount");

    cl_die_if(clSetKernelArg(g_kernel, arg++, sizeof(cl_mem), &d_allSignatures), "set arg allSignatures");

    cl_die_if(clSetKernelArg(g_kernel, arg++, sizeof(cl_mem), &d_sliceKeyCounts), "set arg sliceKeyCounts");
    cl_die_if(clSetKernelArg(g_kernel, arg++, sizeof(cl_mem), &d_sliceKeyBaseOffsets), "set arg sliceKeyBaseOffsets");

    cl_die_if(clSetKernelArg(g_kernel, arg++, sizeof(cl_mem), &d_sliceMaskPositions), "set arg sliceMaskPositions");
    cl_die_if(clSetKernelArg(g_kernel, arg++, sizeof(cl_mem), &d_sliceMaskOffsets),  "set arg sliceMaskOffsets");
    cl_die_if(clSetKernelArg(g_kernel, arg++, sizeof(cl_mem), &d_sliceMaskLengths),  "set arg sliceMaskLengths");

    cl_uint sliceCountU = (cl_uint)g_sliceCount_cached;
    cl_die_if(clSetKernelArg(g_kernel, arg++, sizeof(cl_uint), &sliceCountU), "set arg sliceCount");

    cl_die_if(clSetKernelArg(g_kernel, arg++, sizeof(cl_mem), &d_mitLut), "set arg mitLut");
    cl_die_if(clSetKernelArg(g_kernel, arg++, sizeof(cl_mem), &d_cfdPos),  "set arg cfdPos");
    cl_die_if(clSetKernelArg(g_kernel, arg++, sizeof(cl_mem), &d_cfdPam),  "set arg cfdPam");

    cl_int methodI = static_cast<cl_int>(method);
    cl_die_if(clSetKernelArg(g_kernel, arg++, sizeof(cl_int), &methodI), "set arg method");

    if (useF32) {
        cl_float thr = (cl_float)threshold;
        cl_die_if(clSetKernelArg(g_kernel, arg++, sizeof(cl_float), &thr), "set arg threshold(f32)");
    } else {
        cl_double thr = (cl_double)threshold;
        cl_die_if(clSetKernelArg(g_kernel, arg++, sizeof(cl_double), &thr), "set arg threshold(f64)");
    }

    cl_uint seqLenU = (cl_uint)seqLength;
    cl_die_if(clSetKernelArg(g_kernel, arg++, sizeof(cl_uint), &seqLenU), "set arg seqLength");

    cl_die_if(clSetKernelArg(g_kernel, arg++, sizeof(cl_mem), &d_outMit), "set arg outMit");
    cl_die_if(clSetKernelArg(g_kernel, arg++, sizeof(cl_mem), &d_outCfd), "set arg outCfd");

    // ====== enqueue ======
    const size_t gsz = guideCount ? guideCount : 1;
    err = clEnqueueNDRangeKernel(g_q, g_kernel, 1, nullptr, &gsz, nullptr, 0, nullptr, nullptr);
    cl_die_if(err, "clEnqueueNDRangeKernel");
    clFinish(g_q);

    std::fprintf(stderr, "[CL] score_kernel launched for %zu guides (stub should write zeros)\n", guideCount);

    if (g_profileLogFile) {
        std::fprintf(g_profileLogFile, "[CL] score_kernel launched for %zu guides\n", guideCount);
        std::fflush(g_profileLogFile);
    }


    // ====== read back ======
    if (outMit) {
        if (useF32) {
            std::vector<float> tmp(guideCount);
            cl_die_if(clEnqueueReadBuffer(g_q, d_outMit, CL_TRUE, 0, tmp.size() * sizeof(float), tmp.data(), 0, nullptr, nullptr), "read outMit f32");
            for (size_t i = 0; i < guideCount; ++i) outMit[i] = (double)tmp[i];
        } else {
            cl_die_if(clEnqueueReadBuffer(g_q, d_outMit, CL_TRUE, 0, guideCount * sizeof(double), outMit, 0, nullptr, nullptr), "read outMit f64");
        }
    }
    if (outCfd) {
        if (useF32) {
            std::vector<float> tmp(guideCount);
            cl_die_if(clEnqueueReadBuffer(g_q, d_outCfd, CL_TRUE, 0, tmp.size() * sizeof(float), tmp.data(), 0, nullptr, nullptr), "read outCfd f32");
            for (size_t i = 0; i < guideCount; ++i) outCfd[i] = (double)tmp[i];
        } else {
            cl_die_if(clEnqueueReadBuffer(g_q, d_outCfd, CL_TRUE, 0, guideCount * sizeof(double), outCfd, 0, nullptr, nullptr), "read outCfd f64");
        }
    }

    // cleanup per-call buffers
    clReleaseMemObject(d_query);
    clReleaseMemObject(d_outMit);
    clReleaseMemObject(d_outCfd);

    // GPU path ran
    return true;
}

// ====== shutdown ======
void shutdown_scoring_cl() {
    std::lock_guard<std::mutex> lk(g_mu);

    auto rel = [](cl_mem& m){ if (m) { clReleaseMemObject(m); m = nullptr; } };

    rel(d_offtargets);
    rel(d_allSignatures);
    rel(d_sliceKeyCounts);
    rel(d_sliceKeyBaseOffsets);
    rel(d_sliceMaskPositions);
    rel(d_sliceMaskOffsets);
    rel(d_sliceMaskLengths);
    rel(d_mitLut);
    rel(d_cfdPos);
    rel(d_cfdPam);

    if (g_kernel)   { clReleaseKernel(g_kernel); g_kernel = nullptr; }
    if (g_prog)     { clReleaseProgram(g_prog);  g_prog   = nullptr; }
    if (g_q)        { clReleaseCommandQueue(g_q); g_q     = nullptr; }
    if (g_ctx)      { clReleaseContext(g_ctx);   g_ctx    = nullptr; }

    g_device = nullptr;
    g_platform = nullptr;

    g_total_global_mem_bytes = 0;
    g_static_bytes = 0;
    g_offtargetsCount_cached = 0;
    g_sliceCount_cached = 0;
    g_seqLength_cached = 0;
}

ClScorePrecision cl_get_precision() noexcept {
    std::lock_guard<std::mutex> lk(g_mu);
    return g_prec;
}
std::size_t cl_get_total_global_mem_bytes() noexcept {
    std::lock_guard<std::mutex> lk(g_mu);
    return g_total_global_mem_bytes;
}
std::size_t cl_get_static_bytes() noexcept {
    std::lock_guard<std::mutex> lk(g_mu);
    return g_static_bytes;
}

#else
// if USE_OPENCL_SCORING == 0, keep empty stubs
void cl_set_precision(ClScorePrecision) {}
void cl_set_mit_lut(const double*, std::size_t) {}
void cl_set_cfd_pos_penalties(const double*, std::size_t) {}
void cl_set_cfd_pam_penalties(const double*, std::size_t) {}
void init_scoring_cl(const ScoringIndexMeta&) {}
bool score_batch_cl(const uint64_t*, std::size_t, double*, double*, otScoreMethod, double, std::size_t) { return false; }
void shutdown_scoring_cl() {}

ClScorePrecision cl_get_precision() { return ClScorePrecision::Float32; }
std::size_t cl_get_total_global_mem_bytes() { return 0; }
std::size_t cl_get_static_bytes() { return 0; }
#endif
