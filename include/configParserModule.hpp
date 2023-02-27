#ifndef configParserModuleInclude
#define configParserModuleInclude
#include <vector>
#include <string>
#include <fstream>
#include <filesystem>
#include <unordered_map>
#include <boost/algorithm/string.hpp>
#include "../include/util.hpp"
#include "../include/pipelineModule.hpp"

#if defined(_WIN64)
const std::string nullDir("nul");
#elif defined (__unix__) || defined (__unix)
const std::string nullDir("/dev/null");
#endif

class configParserModule : private pipelineModule
{
public:
	cracklingConfig run(const std::string& configFile);
private:
	void run() final;
};

#endif // !configParserModule
