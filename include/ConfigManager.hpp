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


using std::string;
using std::list;

class ConfigManager
{
private:
	list<string> filesToProcess;
	std::map<string, std::map<string, string>> configMap;
	

public:
	
	ConfigManager(string configFilePath);
	
	int getConsensusToolCount();

	void set(string section, string key, string value);

	int getInt(string section, string key);

	float getFloat(string section, string key);

	double getDouble(string section, string key);

	string getString(string section, string key);

	const char* getCString(string section, string key);

	bool getBool(string section, string key);

	std::filesystem::path getPath(string section, string key);

	list<string> getFilesToProcess();
};