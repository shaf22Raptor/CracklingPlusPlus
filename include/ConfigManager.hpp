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
#include "../include/Helpers.hpp"

#if (_POSIX_C_SOURCE >= 2 || _XOPEN_SOURCE || _BSD_SOURCE || _SVID_SOURCE)
# define p_popen popen
# define p_pclose pclose
#elif defined(_MSC_VER)
# define p_popen _popen
# define p_pclose _pclose
#else
# error "Error, no popen"
#endif

class ConfigManager
{
public:
	explicit ConfigManager(const std::string& configFilePath);

	int getConsensusToolCount();

	void set(const std::string& section, const std::string& key, std::string_view value);

	int getInt(const std::string& section, const std::string& key);

	float getFloat(const std::string& section, const std::string& key);

	double getDouble(const std::string& section, const std::string& key);

	std::string getString(const std::string& section, const std::string& key);

	const char* getCString(const std::string& section, const std::string& key);

	bool getBool(const std::string& section, const std::string& key);

	std::filesystem::path getPath(const std::string& section, const std::string& key);

	std::list<std::string> getFilesToProcess() const;

private:
	std::list<std::string> filesToProcess;
	std::unordered_map<std::string, std::unordered_map<std::string, std::string>> configMap;
};

class InvalidConfiguration : public std::runtime_error
{
	using std::runtime_error::runtime_error;
};