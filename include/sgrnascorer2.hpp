// sgrnascorer2.hpp
#pragma once
#include <string>
#include <map>
#include <ConfigManager.hpp>
#include <Constants.hpp>
#include <Helpers.hpp>

class sgrnascorer2
{
public:
	sgrnascorer2(ConfigManager cm);

	void run(std::map<std::string, std::map<std::string, std::string>>& candidateGuides);

private:
	bool toolIsSelected;
	std::string optimsationLevel;
	int toolCount;
	int consensusN;
	float scoreThreshold;
	char printingBuffer[1024];

};