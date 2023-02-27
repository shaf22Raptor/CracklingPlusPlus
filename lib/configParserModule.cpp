#include "../include/configParserModule.hpp"


using std::cout;
using std::cerr;
using std::endl;
using std::getline;
using std::vector;
using std::string;
using std::unordered_map;
using std::ifstream;
using std::filesystem::path;
using std::filesystem::directory_entry;
using std::filesystem::directory_iterator;
using std::filesystem::exists;
using std::filesystem::is_directory;
using std::filesystem::create_directory;
using boost::algorithm::trim;
using boost::algorithm::to_lower;
using boost::algorithm::split;



class InvalidConfiguration : public std::runtime_error
{
public:
	InvalidConfiguration(const string& error) : std::runtime_error(error) { };
};


cracklingConfig configParserModule::run(const string& configFile)
{
	cracklingConfig config;

	// Check file exists
	path configPath = configFile;
	if (!exists(configPath))
	{
		throw InvalidConfiguration("Could not find the config file specified");
	}

	// Extract values from config files
	unordered_map<string, unordered_map<string, string>> configMap;
	string section;
	string key;
	string value;

	ifstream inFile(configPath, std::ios::binary, std::ios::in);
	for (string currentLine; getline(inFile, currentLine);)
	{
		trim(currentLine);
		if (currentLine[0] == ';' || currentLine == "") { continue; }
		else if (currentLine[0] == '[')
		{
			section = currentLine.substr(currentLine.find("[") + 1, currentLine.find("]"));
			trim(section);
			to_lower(section);
		}
		else
		{
			vector<string> splitLine;
			split(splitLine, currentLine, boost::is_any_of(" "));
			trim(splitLine[0]);
			to_lower(splitLine[0]);
			trim(splitLine[1]);
			configMap[section][splitLine[0]] = splitLine[1];
		}
	}
	inFile.close();

	// Assign config items
	try
	{
		config.general.name = configMap.at("general").at("name");
		config.general.optimisation = optimisationMap.at(configMap.at("general").at("optimisation"));
		config.consensus.n = stoul(configMap.at("consensus").at("n"));
		config.consensus.mm10db = configMap.at("consensus").at("mm10db") == "True";
		config.consensus.sgrnascorer2 = configMap.at("consensus").at("sgrnascorer2") == "True";
		config.consensus.chopchop = configMap.at("consensus").at("chopchop") == "True";
		config.input.exonSequences = configMap.at("input").at("exon-sequences");
		config.input.gffAnnotation = configMap.at("input").at("gff-annotation");
		config.input.bowtie2Index = configMap.at("input").at("bowtie2-index");
		config.input.batchLen = stoull(configMap.at("input").at("batch-size"));
		config.output.filename = configMap.at("output").at("filename");
		config.output.delimiter = configMap.at("output").at("delimter")[0];
		config.offTarget.enabled = configMap.at("offtargetscore").at("enabled") == "True";
		config.offTarget.method = otScoreMethodMap.at(configMap.at("offtargetscore").at("method"));
		config.offTarget.threads = stoul(configMap.at("offtargetscore").at("threads"));
		config.offTarget.pageLen = stoull(configMap.at("offtargetscore").at("page-length"));
		config.offTarget.scoreThreshold = stod(configMap.at("offtargetscore").at("score-threshold"));
		config.offTarget.maxDist = stoul(configMap.at("offtargetscore").at("max-distance"));
		config.sgrnascorer2.model = configMap.at("sgrnascorer2").at("model");
		config.sgrnascorer2.scoreThreshold = stoi(configMap.at("sgrnascorer2").at("score-threshold"));
		config.bowtie2.binary = configMap.at("bowtie2").at("binary");
		config.bowtie2.threads = stoul(configMap.at("bowtie2").at("threads"));
		config.bowtie2.pageLen = stoull(configMap.at("bowtie2").at("page-length"));
		config.rnafold.binary = configMap.at("rnafold").at("binary");
		config.rnafold.threads = stoul(configMap.at("rnafold").at("threads"));
		config.rnafold.pageLen = stoull(configMap.at("rnafold").at("binary"));
		config.rnafold.lowEngeryThreshold = stoi(configMap.at("rnafold").at("binary"));
		config.rnafold.highEngeryThreshold = stoi(configMap.at("rnafold").at("binary"));
	}
	catch (const std::exception& e)
	{
		cerr << fmt::format("Config file is missing a value.\n\t{}", e.what()) << endl;
		throw InvalidConfiguration("Config file is missing some values!");
	}

	// Check that binarys are callable
	int returnCode;

	returnCode = system(fmt::format("{} --version >{} 2>{}", config.bowtie2.binary, nullDir, nullDir).c_str());
	if (returnCode != 0)
	{
		throw InvalidConfiguration("Could not find Bowtie2 binary");
	}

	returnCode = system(fmt::format("{} --version >{} 2>{}", config.rnafold.binary, nullDir, nullDir).c_str());
	if (returnCode != 0)
	{
		throw InvalidConfiguration("Could not find RNAFold binary");
	}

	// Check that the 'n' value for the consensus is valid
	config.consensus.toolCount = (int)config.consensus.mm10db + (int)config.consensus.sgrnascorer2 + (int)config.consensus.chopchop;
	if (config.consensus.n > config.consensus.toolCount)
	{
		throw InvalidConfiguration(fmt::format("The consensus approach is incorrectly set. You have specified {0} tools to be run but the n-value is {1}. Change n to be <= {0}.", config.consensus.toolCount, config.consensus.n));
	}

	// Check that output file doesn't already exist
	path outputDir(config.output.dir);
	if (!exists(outputDir)) { create_directory(outputDir); }
	config.output.filename = (outputDir / fmt::format("{}-{}", config.general.name, config.output.filename)).string();
	if (path outputFile = path(config.output.filename); exists(outputFile))
	{
		throw InvalidConfiguration(fmt::format("The output file already exists: {}.\nTo avoid loosing data, please rename your output file.", config.output.filename));
	}
	config.output.log = (outputDir / fmt::format("{}-{}.log", config.general.name, configPath.stem())).string();
	config.output.errLog = (outputDir / fmt::format("{}-{}.errlog", config.general.name, configPath.stem())).string();

	// Check that all input files exist
	if (!exists(config.input.exonSequences)) { throw InvalidConfiguration(fmt::format("Could not find input file: {}", config.input.exonSequences)); }
	if (!exists(config.input.offtargetSites)) { throw InvalidConfiguration(fmt::format("Could not find input file: {}", config.input.exonSequences)); }
	if (!exists(path(config.input.bowtie2Index.string() + "1.bt2"))) { throw InvalidConfiguration(fmt::format("Could not find input file: {}", config.input.exonSequences)); }

	// Find files to process
	if (is_directory(config.input.exonSequences))
	{
		for (const directory_entry& dir_entry : directory_iterator{config.input.exonSequences})
		{
			config.input.filesToProcess.push_back(dir_entry.path());
		}
	}
	else
	{
		config.input.filesToProcess.push_back(config.input.exonSequences);
	}

	// Generate tempfile names
	config.rnafold.inFile = (outputDir / fmt::format("{}-rnafold-input.txt", config.general.name));
	config.rnafold.outFile = (outputDir / fmt::format("{}-rnafold-output.txt", config.general.name));
	config.bowtie2.inFile = (outputDir / fmt::format("{}-bowtie2-input.txt", config.general.name));
	config.bowtie2.outFile = (outputDir / fmt::format("{}-bowtie2-output.txt", config.general.name));

	return config;
}

void configParserModule::run() {}
