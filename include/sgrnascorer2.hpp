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

	void run(std::map<std::string, std::map<std::string, std::string>> candidateGuides);

private:
	int testedcount;
	int failedcount;
	bool toolIsSelected;
	char printingBuffer[1024];

};