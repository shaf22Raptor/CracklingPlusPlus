// CHOPCHOP.hpp
#pragma once
#include <string>
#include <ConfigManager.hpp>
#include <Constants.hpp>

class CHOPCHOP
{
public:
	CHOPCHOP(ConfigManager cm);

	void run();

private:
	int testedcount;
	int failedcount;
	bool toolIsSelected;

	bool G20(std::string candidateGuide);
};