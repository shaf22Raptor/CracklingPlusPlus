// mm10db.hpp
#pragma once
#include <string>
#include <map>
#include <array>
#include <iterator>
#include <fstream>
#include <list>
#include <regex>
#include <sstream>
#include "../include/ConfigManager.hpp"
#include "../include/Constants.hpp"
#include "../include/Helpers.hpp"

class mm10db
{
public:
	explicit mm10db(ConfigManager& cm);

	void run(std::map<std::string, std::map<std::string, std::string, std::less<>>, std::less<>>& candidateGuides);

	bool static leadingT(std::string_view candidateGuide);

	float static AT_percentage(std::string_view candidateGuide);

	bool static polyT(std::string_view candidateGuide);

	std::string static transToDNA(std::string RNA);

private:
	bool toolIsSelected;
	std::string optimsationLevel;
	int toolCount;
	int consensusN;
	int threadCount;
	std::string RNAFoldOutFile;
	std::string RNAFoldInFile;
	std::string RNAFoldBin;
	float lowEnergyThreshold;
	float highEnergyThreshold;
	int RNAFoldPageLength;
};