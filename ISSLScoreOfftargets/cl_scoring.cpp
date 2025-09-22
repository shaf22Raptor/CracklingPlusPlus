// cl_scoring.cpp
#include "cl_scoring.hpp"

bool init_scoring_cl(const ScoringIndexMeta&) {
    return true; // stub ok
}

bool score_batch_cl(const uint64_t*,
                    size_t,
                    double*,
                    double*,
                    otScoreMethod,
                    double,
                    size_t) {
    return false; // tell caller to use CPU path for now
}

void shutdown_scoring_cl() {
    // stub
}
