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
#include <chrono>
#include <bitset>
#include <iostream>
#include "../include/ConfigManager.hpp"
#include "../include/Constants.hpp"
#include "../include/Helpers.hpp"
#include "../include/phmap/phmap.h"
#include "../include/libpopcnt.h"
#include "../include/cfdPenalties.h"

#ifndef portableStat64
#if (_BSD_SOURCE || _XOPEN_SOURCE >= 500 || _XOPEN_SOURCE && _XOPEN_SOURCE_EXTENDED || /* Since glibc 2.10: */ _POSIX_C_SOURCE >= 200112L)
#include <sys/time.h>
#include <unistd.h>
#define p_stat64 stat64
#elif defined(_WIN64) 
#include "../include/sys/time.h"
#include "../include/unistd.h"
#define p_stat64 _stat64
#endif
#endif // !portableStat64


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
	int threadCount;
	std::string ISSLIndex;
	int maxDist;
	std::string scoreMethod;
	float scoreThreshold;
	int pageLength;

	size_t getFileSize(const char* path);

	uint64_t sequenceToSignature(const char* ptr, const size_t& seqLength);

	std::string signatureToSequence(uint64_t signature, const size_t& seqLength);
};