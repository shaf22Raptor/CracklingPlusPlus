// CHOPCHOP.hpp
#pragma once
#include <string>
#include <map>
#include <ConfigManager.hpp>
#include <Constants.hpp>
#include <Helpers.hpp>

class CHOPCHOP
{
public:
	CHOPCHOP(ConfigManager cm);

	void run(std::map<std::string, std::map<std::string, std::string, std::less<>>>& candidateGuides);
	
	bool static G20(std::string candidateGuide);

private:
	int testedCount;
	int failedCount;
	bool toolIsSelected;
	std::string optimsationLevel;
	int toolCount;
	int consensusN;
};