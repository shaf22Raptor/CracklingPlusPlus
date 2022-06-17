// bowtie2.hpp
#pragma once
#include <string>
#include <map>
#include <vector>
#include <ConfigManager.hpp>
#include <Constants.hpp>
#include <Helpers.hpp>

class bowtie2
{
public:
	bowtie2(ConfigManager cm);

	void run(std::map<std::string, std::map<std::string, std::string>>& candidateGuides);

private:
	bool toolIsSelected;
	std::string optimsationLevel;
	int toolCount;
	int consensusN;
	std::string bowtie2OutFile;
	std::string bowtie2InFile;
	std::string bowtie2Bin;
	std::string bowtie2Index;
	int bowtie2Threads;
	int bowtie2PageLength;

};