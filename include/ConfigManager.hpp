// ConfigManager.hpp
#pragma once
#include <string>
#include <list>
#include <map>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <regex>
#include <filesystem>

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
	std::map<std::string, std::map<std::string, std::string>> configMap;
	


};