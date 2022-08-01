// ISSLOffTargetScoring.hpp
#pragma once
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <climits>
#include <stdio.h>
#include <stdint.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <omp.h>
#include <string>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <map>
#include <chrono>
#include <bitset>
#include <iostream>
#include "../include/ConfigManager.hpp"
#include "../include/Constants.hpp"
#include "../include/Helpers.hpp"
#include "../include/phmap/phmap.h"
#include "../include/libpopcnt.h"
#include "../include/cfdPenalties.h"

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
#include "../include/sys/time.h"
#include "../include/unistd.h"
#else
#include <sys/time.h>
#include <unistd.h>
#endif

#ifndef portableStat64
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
#define p_stat64 _stat64
#else
#define p_stat64 stat64
#endif
#endif // !portableStat64


class ISSLOffTargetScoring
{
public:
	explicit ISSLOffTargetScoring(ConfigManager& cm);

	void run(std::map<std::string, std::map<std::string, std::string, std::less<>>, std::less<>>& candidateGuides);

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