#ifndef utilInclude
#define utilInclude
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

enum otScoreMethod { unknown = 0, mit = 1, cfd = 2, mitAndCfd = 3, mitOrCfd = 4, avgMitCfd = 5 };
enum optimisationLevel { invalid = 0, ultralow = 1, low = 2, medium = 3, high = 4};
const static std::unordered_map<std::string, optimisationLevel> const table = { 
	{"ultralow",optimisationLevel::ultralow}, 
	{"low",optimisationLevel::low}, 
	{"medium",optimisationLevel::medium},
	{"high",optimisationLevel::high}
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
	uint64_t batchLen;
};

struct outputConfig
{
	std::filesystem::path dir;
	std::string filename;
	char delimiter;
};

struct offtargetConfig
{
	otScoreMethod method;
	uint8_t threads;
	uint64_t pageLen;
	uint8_t scoreThreshold;
	uint8_t maxDist;
};

struct sgrnascorer2Config
{
	std::filesystem::path model;
	int scoreThreshold;
};

struct bowtie2Config
{
	std::filesystem::path binary;
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
	offtargetConfig offtarget;
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
	char passedBowtie = CODE_UNTESTED;
	char passedOffTargetScore = CODE_UNTESTED;
	double AT = NULL;
	std::string ssL1 = &CODE_UNTESTED;
	std::string ssStructure = &CODE_UNTESTED;
	double ssEnergy = NULL;
	double sgrnascorer2score = NULL;
	std::string bowtieChr = &CODE_UNTESTED;
	uint64_t bowtieStart = NULL;
	uint64_t bowtieEnd = NULL;
	double mitOfftargetscore = NULL;
	double cfdOfftargetscore = NULL;
};

class ReturnCode : public std::logic_error
{
public:
	ReturnCode() : std::logic_error("The externally called program returned a non-zero value") { };
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