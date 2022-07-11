// offTargetScoring.hpp
#pragma once
#include <string>
#include <map>
#include <vector>
#include <ConfigManager.hpp>
#include <Constants.hpp>
#include <Helpers.hpp>

class offTargetScoring
{
public:
	offTargetScoring(ConfigManager cm);

	void run(std::map<std::string, std::map<std::string, std::string, std::less<>>>& candidateGuides);

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