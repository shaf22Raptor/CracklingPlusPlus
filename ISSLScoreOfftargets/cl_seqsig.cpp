#include "cl_seqsig.hpp"
#include <CL/cl.h>
#include <stdexcept>
#include <cstring>
#include <algorithm>
#include "app_options.hpp" 
#include <vector>            // std::vector in pickDevice
#include <string>            // std::string
#include <cctype>            // std::tolower

static const char* kKernelSrc = R"CLC(
__kernel void seq_to_sig(__global const uchar* seqs,
                         uint                  seq_stride,
                         uint                  seq_len,
                         __constant uchar*     lut,
                         __global ulong*       out_sigs,
                         uint                  total_count)   // <-- NEW
{
    size_t gid = get_global_id(0);
    if (gid >= total_count) return;  // <-- guard for padded global size

    const __global uchar* p = seqs + gid * seq_stride;
    ulong sig = 0ul;
    #pragma unroll 32
    for (uint j = 0; j < seq_len; ++j) {
        uint c = (uint)lut[p[j]] & 3u;
        sig |= ((ulong)c) << (j * 2);
    }
    out_sigs[gid] = sig;
}
)CLC";

static std::string to_lower(std::string s){ for(char& c: s) c=(char)std::tolower((unsigned char)c); return s; }
static bool contains_ci(const std::string& hay, const std::string& needle){
    if (needle.empty()) return false;
    auto H = to_lower(hay); auto N = to_lower(needle);
    return H.find(N) != std::string::npos;
}

struct SeqSigEncoderCL::Impl {
    cl_platform_id  platform{};
    cl_device_id    device{};
    cl_context      ctx{};
    cl_command_queue q{};
    cl_program      prog{};
    cl_kernel       kern{};
    cl_mem          lutBuf{};

    // 256-byte ASCII → 2-bit LUT
    unsigned char lut[256]{};

    void initLUT() {
        std::memset(lut, 0, sizeof(lut));
        lut[(unsigned)'A'] = 0; lut[(unsigned)'a'] = 0;
        lut[(unsigned)'C'] = 1; lut[(unsigned)'c'] = 1;
        lut[(unsigned)'G'] = 2; lut[(unsigned)'g'] = 2;
        lut[(unsigned)'T'] = 3; lut[(unsigned)'t'] = 3;
        // everything else stays 0
    }

    // pick device according to: indices > gfxId > nameContains > default(0, gpu then cpu)
    static void pickDevice(const AppOptions& opts, cl_platform_id& outPlat, cl_device_id& outDev)
    {
        cl_int err=0; cl_uint np=0;
        if (clGetPlatformIDs(0,nullptr,&np)!=CL_SUCCESS || !np) throw std::runtime_error("No OpenCL platform");
        std::vector<cl_platform_id> plats(np);
        clGetPlatformIDs(np, plats.data(), nullptr);

        auto query_name = [&](cl_device_id d)->std::string{
            size_t sz=0; clGetDeviceInfo(d, CL_DEVICE_NAME, 0, nullptr, &sz);
            std::string s(sz,'\0');
            clGetDeviceInfo(d, CL_DEVICE_NAME, sz, s.data(), nullptr);
            if (!s.empty() && s.back()=='\0') s.pop_back();
            return s;
        };

        auto pick_by_indices = [&]()->bool{
            if (opts.clPlatformIndex < 0 || opts.clDeviceIndex < 0) return false;
            if ((size_t)opts.clPlatformIndex >= plats.size()) return false;
            cl_platform_id P = plats[(size_t)opts.clPlatformIndex];
            cl_uint nd=0;
            if (clGetDeviceIDs(P, CL_DEVICE_TYPE_ALL, 0, nullptr, &nd)!=CL_SUCCESS || !nd) return false;
            std::vector<cl_device_id> devs(nd);
            clGetDeviceIDs(P, CL_DEVICE_TYPE_ALL, nd, devs.data(), nullptr);
            if ((size_t)opts.clDeviceIndex >= devs.size()) return false;
            outPlat = P; outDev = devs[(size_t)opts.clDeviceIndex];
            return true;
        };

        auto pick_by_gfx = [&]()->bool{
            if (opts.clDeviceGfxId.empty()) return false;
            for (cl_platform_id P : plats) {
                cl_uint nd=0;
                if (clGetDeviceIDs(P, CL_DEVICE_TYPE_ALL, 0, nullptr, &nd)!=CL_SUCCESS || !nd) continue;
                std::vector<cl_device_id> devs(nd);
                clGetDeviceIDs(P, CL_DEVICE_TYPE_ALL, nd, devs.data(), nullptr);
                for (auto d: devs) {
                    std::string name = query_name(d);
                    if (contains_ci(name, opts.clDeviceGfxId)) {
                        outPlat = P; outDev = d; return true;
                    }
                }
            }
            return false;
        };

        auto pick_by_name = [&]()->bool{
            if (opts.clDeviceNameContains.empty()) return false;
            for (cl_platform_id P : plats) {
                cl_uint nd=0;
                if (clGetDeviceIDs(P, CL_DEVICE_TYPE_ALL, 0, nullptr, &nd)!=CL_SUCCESS || !nd) continue;
                std::vector<cl_device_id> devs(nd);
                clGetDeviceIDs(P, CL_DEVICE_TYPE_ALL, nd, devs.data(), nullptr);
                for (auto d: devs) {
                    std::string name = query_name(d);
                    if (contains_ci(name, opts.clDeviceNameContains)) {
                        outPlat = P; outDev = d; return true;
                    }
                }
            }
            return false;
        };

        auto pick_default = [&]()->bool{
            for (cl_platform_id P : plats) {
                cl_uint nd=0;
                if (clGetDeviceIDs(P, CL_DEVICE_TYPE_GPU, 0, nullptr, &nd)==CL_SUCCESS && nd){
                    std::vector<cl_device_id> devs(nd);
                    clGetDeviceIDs(P, CL_DEVICE_TYPE_GPU, nd, devs.data(), nullptr);
                    outPlat = P; outDev = devs[0]; return true;
                }
            }
            // fallback to CPU if no GPU
            for (cl_platform_id P : plats) {
                cl_uint nd=0;
                if (clGetDeviceIDs(P, CL_DEVICE_TYPE_CPU, 0, nullptr, &nd)==CL_SUCCESS && nd){
                    std::vector<cl_device_id> devs(nd);
                    clGetDeviceIDs(P, CL_DEVICE_TYPE_CPU, nd, devs.data(), nullptr);
                    outPlat = P; outDev = devs[0]; return true;
                }
            }
            return false;
        };

        if (pick_by_indices()) return;
        if (pick_by_gfx())     return;
        if (pick_by_name())    return;
        if (pick_default())    return;
        throw std::runtime_error("Unable to select any OpenCL device");
    }

    Impl(const AppOptions& opts){
        cl_int err=0;

        // --- select device per opts ---
        pickDevice(opts, platform, device);

        // Context + queue
        ctx = clCreateContext(nullptr, 1, &device, nullptr, nullptr, &err);
        if (!ctx || err) throw std::runtime_error("clCreateContext failed");
        q = clCreateCommandQueue(ctx, device, 0, &err);
        if (!q || err) throw std::runtime_error("clCreateCommandQueue failed");

        // Program
        const char* src = kKernelSrc; size_t len = std::strlen(kKernelSrc);
        prog = clCreateProgramWithSource(ctx, 1, &src, &len, &err);
        if (!prog || err) throw std::runtime_error("clCreateProgramWithSource failed");
        err = clBuildProgram(prog, 1, &device, "", nullptr, nullptr);
        if (err != CL_SUCCESS) {
            size_t logSize=0;
            clGetProgramBuildInfo(prog, device, CL_PROGRAM_BUILD_LOG, 0, nullptr, &logSize);
            std::string log(logSize, '\0');
            clGetProgramBuildInfo(prog, device, CL_PROGRAM_BUILD_LOG, logSize, log.data(), nullptr);
            throw std::runtime_error("Build failed:\n" + log);
        }

        kern = clCreateKernel(prog, "seq_to_sig", &err);
        if (!kern || err) throw std::runtime_error("clCreateKernel failed");

        // LUT
        initLUT();
        lutBuf = clCreateBuffer(ctx, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                sizeof(lut), lut, &err);
        if (!lutBuf || err) throw std::runtime_error("LUT buffer creation failed");
    }

    ~Impl(){
        if (lutBuf) clReleaseMemObject(lutBuf);
        if (kern)   clReleaseKernel(kern);
        if (prog)   clReleaseProgram(prog);
        if (q)      clReleaseCommandQueue(q);
        if (ctx)    clReleaseContext(ctx);
    }
};

// ctor/dtor
SeqSigEncoderCL::SeqSigEncoderCL(const AppOptions& opts) : p_(new Impl(opts)) {}
SeqSigEncoderCL::~SeqSigEncoderCL() { delete p_; }

void SeqSigEncoderCL::encode(const char* queryBytes,
                             uint32_t    queryStride,
                             uint32_t    seqLen,
                             uint32_t    queryCount,
                             std::vector<uint64_t>& outSigs)
{
    if (seqLen > 32) throw std::runtime_error("seqLen>32 not supported in 64-bit signature");
    if (queryCount == 0) { outSigs.clear(); return; }
    if (seqLen > queryStride) throw std::runtime_error("seqLen must be <= queryStride");

    outSigs.resize(queryCount);

    cl_int err = 0;

    // Heuristic batch size: ~32 MB input → good PCIe amortization
    const size_t bytesPerSeq = (size_t)queryStride;
    const size_t targetBytes = 32ull * 1024 * 1024;
    uint32_t batch = (bytesPerSeq ? (uint32_t)(targetBytes / bytesPerSeq) : 0u);
    if (batch == 0) batch = 1;
    if (batch > queryCount) batch = queryCount;

    cl_mem inBuf  = clCreateBuffer(p_->ctx, CL_MEM_READ_ONLY,  (size_t)batch * bytesPerSeq, nullptr, &err);
    if (!inBuf || err) throw std::runtime_error("inBuf alloc failed");
    cl_mem outBuf = clCreateBuffer(p_->ctx, CL_MEM_WRITE_ONLY, (size_t)batch * sizeof(cl_ulong), nullptr, &err);
    if (!outBuf || err) {
        clReleaseMemObject(inBuf);
        throw std::runtime_error("outBuf alloc failed");
    }

    // Kernel args that do not change across batches
    if (clSetKernelArg(p_->kern, 1, sizeof(cl_uint), &queryStride) != CL_SUCCESS)
        throw std::runtime_error("SetKernelArg #1 failed");
    if (clSetKernelArg(p_->kern, 2, sizeof(cl_uint), &seqLen) != CL_SUCCESS)
        throw std::runtime_error("SetKernelArg #2 failed");
    if (clSetKernelArg(p_->kern, 3, sizeof(cl_mem),  &p_->lutBuf) != CL_SUCCESS)
        throw std::runtime_error("SetKernelArg #3 failed");

    size_t processed = 0;
    while (processed < queryCount) {
        uint32_t thisBatch = std::min<uint32_t>(batch, queryCount - processed);
        const char* src    = queryBytes + (size_t)processed * queryStride;

        // Upload this batch
        err = clEnqueueWriteBuffer(p_->q, inBuf, CL_TRUE, 0,
                                   (size_t)thisBatch * bytesPerSeq, src, 0, nullptr, nullptr);
        if (err != CL_SUCCESS) throw std::runtime_error("WriteBuffer failed");

        // Per-batch args that change
        if (clSetKernelArg(p_->kern, 0, sizeof(cl_mem), &inBuf) != CL_SUCCESS)
            throw std::runtime_error("SetKernelArg #0 failed");
        if (clSetKernelArg(p_->kern, 4, sizeof(cl_mem), &outBuf) != CL_SUCCESS)
            throw std::runtime_error("SetKernelArg #4 failed");
        cl_uint total = thisBatch;
        if (clSetKernelArg(p_->kern, 5, sizeof(cl_uint), &total) != CL_SUCCESS)  // <-- NEW
            throw std::runtime_error("SetKernelArg #5 failed");

        // Launch with padded global size
        size_t gsz = thisBatch;
        const size_t lsz = 256;
        if (gsz % lsz) gsz = ((gsz + lsz - 1) / lsz) * lsz;

        err = clEnqueueNDRangeKernel(p_->q, p_->kern, 1, nullptr, &gsz, &lsz, 0, nullptr, nullptr);
        if (err != CL_SUCCESS) throw std::runtime_error("Enqueue kernel failed");

        // Read back only the valid results
        err = clEnqueueReadBuffer(p_->q, outBuf, CL_TRUE, 0,
                                  (size_t)thisBatch * sizeof(cl_ulong),
                                  outSigs.data() + processed, 0, nullptr, nullptr);
        if (err != CL_SUCCESS) throw std::runtime_error("ReadBuffer failed");

        processed += thisBatch;
    }

    clReleaseMemObject(outBuf);
    clReleaseMemObject(inBuf);
}
