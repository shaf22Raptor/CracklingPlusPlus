// cl_scoring.cpp
#include "cl_scoring.hpp"
#include <string>
#include <mutex>
#include <vector> 

static std::mutex g_mu;
static ClScorePrecision g_prec = ClScorePrecision::Float32;

// Host copies of LUTs (as double for fidelity; we’ll cast at upload)
static std::vector<double> g_mitLut;
static std::vector<double> g_cfdPos;
static std::vector<double> g_cfdPam;

void cl_set_precision(ClScorePrecision p) {
    std::lock_guard<std::mutex> lk(g_mu);
    g_prec = p;
}

void cl_set_mit_lut(const double* data, std::size_t len) {
    std::lock_guard<std::mutex> lk(g_mu);
    g_mitLut.assign(data, data + len);
}

void cl_set_cfd_pos_penalties(const double* data, std::size_t len) {
    std::lock_guard<std::mutex> lk(g_mu);
    g_cfdPos.assign(data, data + len);
}

void cl_set_cfd_pam_penalties(const double* data, std::size_t len) {
    std::lock_guard<std::mutex> lk(g_mu);
    g_cfdPam.assign(data, data + len);
}

void init_scoring_cl(const ScoringIndexMeta& meta)
{
    // current body; remove any `return true/false;`
}

bool score_batch_cl(const uint64_t* querySigs,
                    std::size_t     guideCount,
                    double*         outMit,
                    double*         outCfd,
                    otScoreMethod   method,
                    double          threshold,
                    std::size_t     seqLength)
{
    // stub body for now
    return false;
}
void shutdown_scoring_cl() {
    // stub
}
