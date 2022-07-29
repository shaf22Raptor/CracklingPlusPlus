// ConfigManager.hpp
#pragma once
#include <string>
#include <list>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <regex>
#include <filesystem>
#include <Helpers.hpp>

#if defined(_POSIX_C_SOURCE) && (_POSIX_C_SOURCE >= 2)
# define portablePopen popen
# define portablePclose pclose
#elif defined(_MSC_VER)
# define portablePopen _popen
# define portablePclose _pclose
#else
# error "Error, no popen"
#endif

class ConfigManager
{
public:
	ConfigManager(std::string configFilePath);
	
	int getConsensusToolCount();

	void set(std::string section, std::string key, std::string value);

	int getInt(std::string section, std::string key);

	float getFloat(std::string section, std::string key);

	double getDouble(std::string section, std::string key);

	std::string getString(std::string section, std::string key);

	const char* getCString(std::string section, std::string key);

	bool getBool(std::string section, std::string key);

	std::filesystem::path getPath(std::string section, std::string key);

	std::list<std::string> getFilesToProcess();

private:
	std::list<std::string> filesToProcess;
	std::unordered_map<std::string, std::unordered_map<std::string, std::string>> configMap;
};