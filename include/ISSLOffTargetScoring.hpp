
#pragma once
#include <cfdPenalties.h>
#include <unordered_map>
#include <unordered_set>
#include <phmap/phmap.h>
#include <omp.h>
#include <ConfigManager.hpp>
#include <Constants.hpp>
#include <Helpers.hpp>

#ifdef _WIN32
#define portable_stat64 _stat64
#define portable_popcount __popcnt64
#elif defined __linux__
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#define portable_popcount __builtin_popcountll
#define portable_stat64 stat64
#endif

class ISSLOffTargetScoring
{
public:
	explicit ISSLOffTargetScoring(ConfigManager& cm);

	void run(std::unordered_map<std::string, std::unordered_map<std::string, std::string>>& candidateGuides);

private:
	bool toolIsSelected;
	std::string optimsationLevel;
	int toolCount;
	int consensusN;
	std::string offTargetScoreOutFile;
	std::string offTargetScoreInFile;
	std::string offTargetScoreBin;
	std::string offTargetScoreIndex;
	std::string offTargetScoreMaxDist;
	std::string offTargetScoreMethod;
	float offTagertScoreThreshold;
	int offTargetScorePageLength;

	size_t getFileSize(const char* path);

	uint64_t sequenceToSignature(const char* ptr, const size_t& seqLength);

	std::string signatureToSequence(uint64_t signature, const size_t& seqLength);
};