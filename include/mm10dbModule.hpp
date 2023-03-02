#ifndef mm10dbModuleInclude
#define mm10dbModuleInclude
#include <vector>
#include <string>
#include <fstream>
#include <filesystem>
#include <boost/algorithm/string.hpp>
#include <boost/regex.hpp>
#include "../include/consensusModule.hpp"
#include "../include/util.hpp"

class mm10dbModule : private consensusModule
{
public:
	mm10dbModule(cracklingConfig config);
	void run(std::vector<guideResults>& candidateGuides) final;
private:
	rnafoldConfig config;
	bool leadingT(std::string_view guide);
	double AT_percent(std::string_view guide);
	bool polyT(std::string_view guide);
	std::string transToDNA(std::string RNA);
	bool processGuide(const guideResults& guide) final;
	
};

#endif // !mm1odbModuleInclude
