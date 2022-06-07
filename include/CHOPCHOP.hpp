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

	void run(std::map<std::string, std::map<std::string, std::string>> candidateGuides);

private:
	int testedcount;
	int failedcount;
	bool toolIsSelected;

	bool G20(std::string candidateGuide);
};