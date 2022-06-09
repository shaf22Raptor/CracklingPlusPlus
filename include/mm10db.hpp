// mm10db.hpp
#pragma once
#include <string>
#include <map>
#include <ConfigManager.hpp>
#include <Constants.hpp>
#include <Helpers.hpp>

class mm10db
{
public:
	mm10db(ConfigManager cm);

	void run(std::map<std::string, std::map<std::string, std::string>> candidateGuides);

	bool leadingT(std::string candidateGuide);

	float AT_percentage(std::string candidateGuide);

	bool polyT(std::string candidateGuide);

	std::string transToDNA(std::string RNA);

private:
	int testedcount;
	int failedcount;
	bool toolIsSelected;
	char printingBuffer[1024];

};