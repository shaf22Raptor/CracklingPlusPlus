#include "cl_seqsig.hpp"
#include <CL/cl.h>
#include <stdexcept>
#include <cstring>
#include <algorithm>

static const char* kKernelSrc = R"CLC(
__kernel void seq_to_sig(__global const uchar* seqs,
                         uint                  seq_stride,
                         uint                  seq_len,
                         __constant uchar*     lut,
                         __global ulong*       out_sigs)
{
    size_t gid = get_global_id(0);
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

    Impl() {
        cl_int err = 0;

        // Platform + device (pick default GPU; fall back to CPU if needed)
        cl_uint np=0; err=clGetPlatformIDs(0,nullptr,&np);
        if (err || !np) throw std::runtime_error("No OpenCL platform");
        std::vector<cl_platform_id> plats(np);
        clGetPlatformIDs(np, plats.data(), nullptr);
        platform = plats[0];

        cl_uint nd=0;
        err = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 1, &device, &nd);
        if (err != CL_SUCCESS) {
            // fallback to CPU device if no GPU
            err = clGetDeviceIDs(platform, CL_DEVICE_TYPE_CPU, 1, &device, &nd);
            if (err != CL_SUCCESS) throw std::runtime_error("No OpenCL device");
        }

        // Context + queue
        ctx = clCreateContext(nullptr, 1, &device, nullptr, nullptr, &err);
        if (!ctx || err) throw std::runtime_error("clCreateContext failed");
        q = clCreateCommandQueue(ctx, device, 0, &err);
        if (!q || err) throw std::runtime_error("clCreateCommandQueue failed");

        // Program
        const char* src = kKernelSrc;
        size_t len = std::strlen(kKernelSrc);
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

        // LUT buffer in constant memory (implemented as read-only buffer)
        initLUT();
        lutBuf = clCreateBuffer(ctx, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                sizeof(lut), lut, &err);
        if (!lutBuf || err) throw std::runtime_error("LUT buffer creation failed");
    }

    ~Impl() {
        if (lutBuf) clReleaseMemObject(lutBuf);
        if (kern)   clReleaseKernel(kern);
        if (prog)   clReleaseProgram(prog);
        if (q)      clReleaseCommandQueue(q);
        if (ctx)    clReleaseContext(ctx);
    }
};

SeqSigEncoderCL::SeqSigEncoderCL() : p_(new Impl) {}
SeqSigEncoderCL::~SeqSigEncoderCL() { delete p_; }

void SeqSigEncoderCL::encode(const char* queryBytes,
                             uint32_t    queryStride,
                             uint32_t    seqLen,
                             uint32_t    queryCount,
                             std::vector<uint64_t>& outSigs)
{
    if (seqLen > 32) throw std::runtime_error("seqLen>32 not supported in 64-bit signature");
    outSigs.resize(queryCount);

    cl_int err = 0;

    // Heuristic batch size: ~32 MB input → good PCIe amortization
    const size_t bytesPerSeq = queryStride;
    const size_t targetBytes = 32ull * 1024 * 1024;
    uint32_t batch = std::max<uint32_t>(1u, (uint32_t)(targetBytes / bytesPerSeq));
    batch = std::min<uint32_t>(batch, queryCount);

    cl_mem inBuf  = clCreateBuffer(p_->ctx, CL_MEM_READ_ONLY,  (size_t)batch * bytesPerSeq, nullptr, &err);
    if (!inBuf || err) throw std::runtime_error("inBuf alloc failed");
    cl_mem outBuf = clCreateBuffer(p_->ctx, CL_MEM_WRITE_ONLY, (size_t)batch * sizeof(cl_ulong), nullptr, &err);
    if (!outBuf || err) { clReleaseMemObject(inBuf); throw std::runtime_error("outBuf alloc failed"); }

    // Kernel args that do not change across batches
    clSetKernelArg(p_->kern, 1, sizeof(cl_uint), &queryStride);
    clSetKernelArg(p_->kern, 2, sizeof(cl_uint), &seqLen);
    clSetKernelArg(p_->kern, 3, sizeof(cl_mem),  &p_->lutBuf);

    size_t processed = 0;
    while (processed < queryCount) {
        uint32_t thisBatch = std::min<uint32_t>(batch, queryCount - processed);
        const char* src = queryBytes + (size_t)processed * queryStride;

        // Upload this batch
        err = clEnqueueWriteBuffer(p_->q, inBuf, CL_TRUE, 0, (size_t)thisBatch * bytesPerSeq, src, 0, nullptr, nullptr);
        if (err != CL_SUCCESS) throw std::runtime_error("WriteBuffer failed");

        // Per-batch args that change
        clSetKernelArg(p_->kern, 0, sizeof(cl_mem), &inBuf);
        clSetKernelArg(p_->kern, 4, sizeof(cl_mem), &outBuf);

        // Launch
        size_t gsz = thisBatch;
        // Round up to a multiple of 256 for better occupancy (optional)
        const size_t lsz = 256;
        if (gsz % lsz) gsz = ((gsz + lsz - 1) / lsz) * lsz;

        err = clEnqueueNDRangeKernel(p_->q, p_->kern, 1, nullptr, &gsz, &lsz, 0, nullptr, nullptr);
        if (err != CL_SUCCESS) throw std::runtime_error("Enqueue kernel failed");

        // Read back
        err = clEnqueueReadBuffer(p_->q, outBuf, CL_TRUE, 0, (size_t)thisBatch * sizeof(cl_ulong),
                                  outSigs.data() + processed, 0, nullptr, nullptr);
        if (err != CL_SUCCESS) throw std::runtime_error("ReadBuffer failed");

        processed += thisBatch;
    }

    clReleaseMemObject(outBuf);
    clReleaseMemObject(inBuf);
}
