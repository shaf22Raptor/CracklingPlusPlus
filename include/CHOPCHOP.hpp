// CHOPCHOP.hpp
#pragma once
#include <string>
#include <unordered_map>
#include <ConfigManager.hpp>
#include <Constants.hpp>
#include <Helpers.hpp>

class CHOPCHOP
{
public:
	explicit CHOPCHOP(ConfigManager& cm);

	void run(std::unordered_map<std::string, std::unordered_map<std::string, std::string>>& candidateGuides);
	
	bool static G20(std::string_view candidateGuide);

private:
	int testedCount;
	int failedCount;
	bool toolIsSelected;
	std::string optimsationLevel;
	int toolCount;
	int consensusN;
};

class G20Input : public std::logic_error
{
public:
	G20Input() : std::logic_error("CHOPCHOP G20: Input lenght must be >= 20!") { };
};