// ISSLClustering.hpp
#pragma once
#define _POSIX_C_SOURCE 200112L
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
#include <thread>
#include "../include/ConfigManager.hpp"
#include "../include/Constants.hpp"
#include "../include/Helpers.hpp"
#include "../include/phmap/phmap.h"
#include "../include/libpopcnt.h"
#include "../include/otScorePenalties.hpp"

#ifndef portableStat64
#define portableStat64
#if defined(_WIN64)
#include "../include/unistd.h"
#define p_stat64 _stat64
#elif defined(unix) || defined(__unix__) || defined(__unix)
#include <unistd.h>
#define p_stat64 stat64
#else
# error "Error, no stat function"
#endif
#endif // !portableStat64


class ISSLClustering
{
public:
	explicit ISSLClustering(ConfigManager& cm);

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