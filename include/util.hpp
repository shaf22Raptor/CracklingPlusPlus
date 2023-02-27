#ifndef utilInclude
#define utilInclude
#include <vector>
#include <string>
#include <filesystem>
#include <unordered_map>
#include <iostream>
#define FMT_HEADER_ONLY
#include "../include/fmt/format.h"

template <class Char>
class commaFormat : public std::numpunct<Char> {
public:
	std::string do_grouping() const { return "\3"; }
	Char do_thousands_sep() const { return ','; }
};

enum otScoreMethod { mit = 0, cfd = 1, mitAndCfd = 2, mitOrCfd = 3, avgMitCfd = 4 };
enum optimisationLevel { ultralow = 0, low = 1, medium = 2, high = 3};

const static std::unordered_map<std::string, optimisationLevel> optimisationMap = {
	{"ultralow",optimisationLevel::ultralow}, 
	{"low",optimisationLevel::low}, 
	{"medium",optimisationLevel::medium},
	{"high",optimisationLevel::high}
};

const static std::unordered_map<std::string, otScoreMethod> otScoreMethodMap = {
	{"mit",otScoreMethod::mit},
	{"cfd",otScoreMethod::cfd},
	{"mitAndCfd",otScoreMethod::mitAndCfd},
	{"mitOrCfd",otScoreMethod::mitOrCfd},
	{"avgMitCfd",otScoreMethod::avgMitCfd}
};

const char CODE_ACCEPTED = '1';
const char CODE_REJECTED = '0';
const char CODE_UNTESTED = '?';
const char CODE_AMBIGUOUS = '-';
const char CODE_ERROR = '!';

struct generalConfig
{
	std::string name;
	optimisationLevel optimisation;
};

struct consensusConfig
{
	uint8_t n;
	uint8_t toolCount;
	bool mm10db;
	bool sgrnascorer2;
	bool chopchop;
};

struct inputConfig
{
	std::filesystem::path exonSequences;
	std::filesystem::path offtargetSites;
	std::filesystem::path gffAnnotation;
	std::filesystem::path bowtie2Index;
	std::vector<std::filesystem::path> filesToProcess;
	uint64_t batchLen;
};

struct outputConfig
{
	std::filesystem::path dir;
	std::filesystem::path filename;
	std::filesystem::path log;
	std::filesystem::path errLog;
	char delimiter;
};

struct offTargetConfig
{
	bool enabled;
	otScoreMethod method;
	uint8_t threads;
	uint64_t pageLen;
	double scoreThreshold;
	uint8_t maxDist;
};

struct sgrnascorer2Config
{
	std::filesystem::path model;
	int8_t scoreThreshold;
};

struct bowtie2Config
{
	std::filesystem::path binary;
	std::filesystem::path inFile;
	std::filesystem::path outFile;
	uint8_t threads;
	uint64_t pageLen;
};

struct rnafoldConfig
{
	std::filesystem::path binary;
	std::filesystem::path inFile;
	std::filesystem::path outFile;
	uint8_t threads;
	uint64_t pageLen;
	int lowEngeryThreshold;
	int highEngeryThreshold;
};

struct cracklingConfig
{
	generalConfig general;
	consensusConfig consensus;
	inputConfig input;
	outputConfig output;
	offTargetConfig offTarget;
	sgrnascorer2Config sgrnascorer2;
	bowtie2Config bowtie2;
	rnafoldConfig rnafold;
};

struct guideResults
{
	std::string seq;
	std::string header = NULL;
	uint64_t start = NULL;
	uint64_t end = NULL;
	char strand = CODE_AMBIGUOUS;
	bool isUnique = false;
	char passedG20 = CODE_UNTESTED;
	char passedAvoidLeadingT = CODE_UNTESTED;
	char passedATPercent = CODE_UNTESTED;
	char passedTTTT = CODE_UNTESTED;
	char passedSecondaryStructure = CODE_UNTESTED;
	char acceptedByMm10db = CODE_UNTESTED;
	char acceptedBySgRnaScorer = CODE_UNTESTED;
	uint8_t consensusCount = NULL;
	char passedBowtie2 = CODE_UNTESTED;
	char passedOffTargetScore = CODE_UNTESTED;
	double AT = NULL;
	std::string ssL1 = &CODE_UNTESTED;
	std::string ssStructure = &CODE_UNTESTED;
	double ssEnergy = NULL;
	double sgrnascorer2score = NULL;
	std::string bowtie2Chr = &CODE_UNTESTED;
	uint64_t bowtie2Start = NULL;
	uint64_t bowtie2End = NULL;
	double mitOfftargetscore = NULL;
	double cfdOfftargetscore = NULL;
};

class ReturnCode : public std::logic_error
{
public:
	ReturnCode() : std::logic_error("The externally called program returned a non-zero value") { };
};

class tempFileSystemError : public std::runtime_error
{
public:
	tempFileSystemError() : std::runtime_error("Unable to create temp working dir") { };
};

class InvalidConfiguration : public std::runtime_error
{
public:
	InvalidConfiguration(const std::string& error) : std::runtime_error(error) { };
};

void runner(const char* args)
{
	std::cout << fmt::format("| Calling: {}", args) << std::endl;
	try
	{
		int returnCode = system(args);
		if (returnCode != 0)
		{
			throw ReturnCode();
		}
	}
	catch (const ReturnCode& e)
	{
		std::cerr << e.what() << std::endl;
		return;
	}
	std::cout << "| Finished" << std::endl;
	return;
}

#endif // !utilInclude