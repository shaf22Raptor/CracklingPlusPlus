// offTargetScoring.hpp
#pragma once
#include <string>
#include <map>
#include <vector>
#include "../include/ConfigManager.hpp"
#include "../include/Constants.hpp"
#include "../include/Helpers.hpp"

class offTargetScoring
{
public:
	explicit offTargetScoring(ConfigManager& cm);

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

};