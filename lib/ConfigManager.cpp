// ConfigManager Class
#include <ConfigManager.hpp>

using std::string;
using std::list;
using std::filesystem::path;

ConfigManager::ConfigManager(string configFilePath) 
{
	/*
		Check config file exists
	*/ 
	// Create path object
	path configPathObject = configFilePath;
	// Normalise path delimiters to the systems preferred delimiters
	configPathObject = configPathObject.make_preferred();
	// Check file exists
	if (!std::filesystem::exists(configPathObject))
	{
		// File doesn't exist, Throw error
		throw std::runtime_error("Could not find the config file specified");
	}


	/*
		Read config file into map<string, map<string,string>> format
	*/
	// Assign strings to represent appropriate values
	string section, key, value;
	// Regex matching variables
	std::smatch match;
	std::regex expression("\\[.*\\]");
	// Open config file
	std::ifstream configFile(configFilePath);
	// Readlines untill EOF
	for (string currentLine; std::getline(configFile, currentLine);)
	{
		// Remove white spaces and resize string
		currentLine = trim(currentLine);
		// Comment line or empty line, ignore
		if (currentLine[0] == ';' || currentLine == "") { continue; }
		// Section header line, update section variable
		else if (currentLine[0] == '[')
		{
			// Extract section header only 
			std::regex_search(currentLine, match, expression);
			section = match.str().substr(1, match.str().length() - 2);
			std::transform(section.begin(), section.end(), section.begin(), [](unsigned char c) { return std::tolower(c); });
		}
		// Entry line, Key and Value seperated by '='
		else
		{
			// Before '='
			key = trim(currentLine.substr(0, currentLine.find("=")));
			std::transform(key.begin(), key.end(), key.begin(), [](unsigned char c) { return std::tolower(c); });
			// After '='
			value = trim(currentLine.substr(currentLine.find("=")+1, currentLine.length()-1));
			std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) { return std::tolower(c); });
			// Update configMap
			set(section, key, value);
		}
	}


	/*
		Normalise all paths
	*/
	// Ensure the path delimiters are set to the system preferred
	set("input", "exon-sequences", getPath("input", "exon-sequences").make_preferred().string());
	set("input", "offtarget-sites", getPath("input", "offtarget-sites").make_preferred().string());
	set("input", "gff-annotation", getPath("input", "gff-annotation").make_preferred().string());
	set("input", "bowtie2-index", getPath("input", "bowtie2-index").make_preferred().string());
	set("output", "dir", getPath("output", "dir").make_preferred().string());
	set("offtargetscore", "binary", getPath("offtargetscore", "binary").make_preferred().string());
	set("sgrnascorer2", "model", getPath("sgrnascorer2", "model").make_preferred().string());
	set("bowtie2", "binary", getPath("bowtie2", "binary").make_preferred().string());
	set("rnafold", "binary", getPath("rnafold", "binary").make_preferred().string());


	/*
		Check that the config file is valid
	*/
	// Run Validate function
	char printingBuffer[1024];	// Buffer used to format strings
	char stdoutBuffer[1024];	// Buffer used to collet stdout
	int returnCode;				// Return code used to check if the binary was run successfuly
	string resultOutput;		// String that collects the stdout
	FILE* stdoutStream;			// Stream to collect stdout from popen calls

	// Check that binarys are callable

	// Check ISSL
	snprintf(printingBuffer, 1024, "%s 2>&1", getCString("offtargetscore", "binary"));
	stdoutStream = portablePopen(printingBuffer, "r");
	while (fgets(stdoutBuffer, 1024, stdoutStream) != NULL)
	{
		resultOutput.append(stdoutBuffer);
	}
	portablePclose(stdoutStream);

	if (resultOutput != "Usage: isslscoreofftargets.exe [issltable] [query file] [max distance] [score-threshold] [score-method]\n") 
	{ 
		throw std::runtime_error("Could not find Off-target scoring binary");
	}

	// Check bowtie2
	snprintf(printingBuffer, 1024, "%s --version", getCString("bowtie2", "binary"));
	stdoutBuffer[0] = '\0';
	stdoutStream = portablePopen(printingBuffer, "r");
	returnCode = portablePclose(stdoutStream);

	if (returnCode != 0)
	{
		throw std::runtime_error("Could not find Bowtie2 binary");
	}

	// Check rnafold
	snprintf(printingBuffer, 1024, "%s --version", getCString("rnafold", "binary"));
	stdoutBuffer[0] = '\0';
	stdoutStream = portablePopen(printingBuffer, "r");
	returnCode = portablePclose(stdoutStream);

	if (returnCode != 0)
	{
		throw std::runtime_error("Could not find RNAFold binary");
	}


	// Check that the 'n' value for the consensus is valid
	int toolCount = getConsensusToolCount();
	int n = getInt("consensus","n");
	if (n > toolCount)
	{
		throw std::runtime_error("The consensus approach is incorrectly set. You have specified %d tools to be run but the n-value is %d. Change n to be <= %d.");
	}

	// Check that output file doesn't already exist
	// Generate outputfile name
	snprintf(printingBuffer, 1024, "%s-%s", getCString("general", "name"), getCString("output", "filename"));
	// Get output dir path object
	path outputDirPathObject = getPath("output", "dir");
	// Append output file name to output dir
	set("output", "file", (outputDirPathObject / printingBuffer).string());
	if (std::filesystem::exists(getPath("output","file")))
	{
		snprintf(printingBuffer, 1024, "The output file already exists: %s.\nTo avoid loosing data, please rename your output file.", getCString("output", "file"));
		throw std::runtime_error(printingBuffer);
	};

	// Check all the fields have values
	if (
		(getString("general", "name") == "") || 
		(getString("general", "optimisation") == "") ||
		(getString("consensus", "n") == "") ||
		(getString("consensus", "mm10db") == "") ||
		(getString("consensus", "sgrnascorer2") == "") ||
		(getString("consensus", "chopchop") == "") ||
		(getString("input", "exon-sequences") == "") ||
		(getString("input", "offtarget-sites") == "") ||
		(getString("input", "gff-annotation") == "") ||
		(getString("input", "bowtie2-index") == "") ||
		(getString("input", "batch-size") == "") ||
		(getString("output", "dir") == "") ||
		(getString("output", "filename") == "") ||
		(getString("output", "delimiter") == "") ||
		(getString("offtargetscore", "enabled") == "") ||
		(getString("offtargetscore", "binary") == "") ||
		(getString("offtargetscore", "method") == "") ||
		(getString("offtargetscore", "threads") == "") ||
		(getString("offtargetscore", "page-length") == "") ||
		(getString("offtargetscore", "score-threshold") == "") ||
		(getString("offtargetscore", "max-distance") == "") ||
		(getString("sgrnascorer2", "model") == "") ||
		(getString("sgrnascorer2", "score-threshold") == "") ||
		(getString("bowtie2", "binary") == "") ||
		(getString("bowtie2", "threads") == "") ||
		(getString("bowtie2", "page-length") == "") ||
		(getString("rnafold", "binary") == "") ||
		(getString("rnafold", "threads") == "") ||
		(getString("rnafold", "page-length") == "") ||
		(getString("rnafold", "low_energy_threshold") == "") ||
		(getString("rnafold", "high_energy_threshold") == "")
	)
	{
		throw std::runtime_error("Configuration file is missing some fields!");
	}


	/*
		Generate files to process
	*/
	// Create path object for input
	path inputPathObject = getPath("input","exon-sequences");
	// Check for directory or file
	if (std::filesystem::is_directory(inputPathObject))
	{
		// Iterate through directory to find all files to process
		for (const std::filesystem::directory_entry& dir_entry :
			std::filesystem::directory_iterator{ inputPathObject })
		{
			// Normalise path and add to list
			path fileToProcess = dir_entry.path();
			filesToProcess.push_back(fileToProcess.make_preferred().string());
		}
	}
	else if (std::filesystem::is_regular_file(inputPathObject))
	{
		// Only one file to process
		filesToProcess.push_back(inputPathObject.string());
	}

	/*
		Generate temp output file names
	*/
	snprintf(printingBuffer, 1024, "%s-rnafold-input.txt", getCString("general", "name"));
	set("rnafold","input", (outputDirPathObject / printingBuffer).string());

	snprintf(printingBuffer, 1024, "%s-rnafold-output.txt", getCString("general", "name"));
	set("rnafold", "output", (outputDirPathObject / printingBuffer).string());

	snprintf(printingBuffer, 1024, "%s-offtargetscore-input.txt", getCString("general", "name"));
	set("offtargetscore", "input", (outputDirPathObject / printingBuffer).string());

	snprintf(printingBuffer, 1024, "%s-offtargetscore-output.txt", getCString("general", "name"));
	set("offtargetscore", "output", (outputDirPathObject / printingBuffer).string());

	snprintf(printingBuffer, 1024, "%s-bowtie-input.txt", getCString("general", "name"));
	set("bowtie2", "input", (outputDirPathObject / printingBuffer).string());

	snprintf(printingBuffer, 1024, "%s-bowtie-output.txt", getCString("general", "name"));
	set("bowtie2", "output", (outputDirPathObject / printingBuffer).string());

	snprintf(printingBuffer, 1024, "%s-%s.log", getCString("general", "name"), configPathObject.stem().string().c_str());
	set("output", "log", (outputDirPathObject / printingBuffer).string());

	snprintf(printingBuffer, 1024, "%s-%s.errlog", getCString("general", "name"), configPathObject.stem().string().c_str());
	set("output", "error", (outputDirPathObject / printingBuffer).string());
}

int ConfigManager::getConsensusToolCount()
{
	bool mm10db = getBool("consensus", "mm10db");
	bool sgrnascorer2 = getBool("consensus", "sgrnascorer2");
	bool chopchop = getBool("consensus", "chopchop");
	return mm10db+sgrnascorer2+chopchop;
}

void ConfigManager::set(string section, string key, string value)
{
	configMap[section][key] = value;
}

int ConfigManager::getInt(string section, string key)
{
	return std::stoi(this->configMap[section][key]);
}

float ConfigManager::getFloat(string section, string key)
{
	return std::stof(this->configMap[section][key]);
}

double ConfigManager::getDouble(string section, string key)
{
	return std::stod(this->configMap[section][key]);
}

string ConfigManager::getString(string section, string key)
{
	return this->configMap[section][key];
}

const char* ConfigManager::getCString(string section, string key) {
	return this->configMap[section][key].c_str();
}

bool ConfigManager::getBool(string section, string key)
{
	string boolValue = this->configMap[section][key];
	std::transform(boolValue.begin(), boolValue.end(), boolValue.begin(), [](unsigned char c) { return std::tolower(c); });
	if (boolValue == "true") { return true; }
	else if (boolValue == "false") { return false; }
	else 
	{ 
		throw std::invalid_argument("The value selected is not of the type bool!");
	}
}

path ConfigManager::getPath(string section, string key)
{
	return path(this->configMap[section][key]);
}

list<string> ConfigManager::getFilesToProcess() {
	return this->filesToProcess;
}
