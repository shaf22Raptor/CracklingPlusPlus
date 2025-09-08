#pragma once
#ifdef MANPROF_ENABLE
#include <cstdio>
#include <chrono>
#include <cinttypes>
#include <omp.h>

extern std::FILE* g_profileLogFile;
extern std::chrono::steady_clock::time_point g_progStart;

inline uint64_t us_since_start() {
    using namespace std::chrono;
    return (uint64_t)duration_cast<std::chrono::microseconds>(
        std::chrono::steady_clock::now() - g_progStart).count();
}

inline void TRACE_EVT(const char* section, const char* evt, int guideIdx=-1, int sliceIdx=-1) {
    int tid = omp_in_parallel() ? omp_get_thread_num() : -1;
    uint64_t us = us_since_start();
#pragma omp critical(ProfileTraceIo)
    {
        std::fprintf(stderr, "Thread %d | us=%" PRIu64 " | %s | [%s] | guide=%d | slice=%d\n",
                     tid, us, evt, section, guideIdx, sliceIdx);
        if (g_profileLogFile) {
            std::fprintf(g_profileLogFile, "Thread %d | us=%" PRIu64 " | %s | [%s] | guide=%d | slice=%d\n",
                         tid, us, evt, section, guideIdx, sliceIdx);
            std::fflush(g_profileLogFile);
        }
        std::fflush(stderr);
    }
}
#endif