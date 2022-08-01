// bowtie2.hpp
#pragma once
#include <string>
#include <map>
#include <vector>
#include "../include/ConfigManager.hpp"
#include "../include/Constants.hpp"
#include "../include/Helpers.hpp"

class bowtie2
{
public:
	explicit bowtie2(ConfigManager& cm);

	void run(std::map<std::string, std::map<std::string, std::string, std::less<>>, std::less<>>& candidateGuides);

private:
	bool toolIsSelected;
	std::string optimsationLevel;
	int toolCount;
	int consensusN;
	int threadCount;
	std::string bowtie2OutFile;
	std::string bowtie2InFile;
	std::string bowtie2Bin;
	std::string bowtie2Index;
	int bowtie2PageLength;

};