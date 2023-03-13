#ifndef outputModuleInclude
#define outputModuleInclude
#include <vector>
#include <filesystem>
#include <fstream>
#include "../include/util.hpp"
#include "../include/pipelineModule.hpp"

class outputModule : private pipelineModule
{
public:
	outputModule(const cracklingConfig& config);
	void run(std::vector<guideResults>& candidateGuides);
private:
	bool firstWrite;
	std::filesystem::path outputFile;
	void run() final;
};

#endif // !outputModuleInclude
