#ifndef otScorePenalties
#define otScorePenalties
#include<unordered_map>

extern const std::unordered_map<uint64_t, double> precalculatedMITScores;
extern const double cfdPosPenalties[320];
extern const double cfdPamPenalties[16];

#endif // !otScorePenalties
