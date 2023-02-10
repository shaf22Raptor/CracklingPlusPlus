#ifndef utilInclude
#define utilInclude
#include <string>
#include <filesystem>

enum otScoreMethod { unknown = 0, mit = 1, cfd = 2, mitAndCfd = 3, mitOrCfd = 4, avgMitCfd = 5 };
enum optimisationLevel { unknown = 0, ultralow = 1, low = 2, medium = 3, high = 4};

const char CODE_ACCEPTED = '1';
const char CODE_REJECTED = '0';
const char CODE_UNTESTED = '?';
const char CODE_AMBIGUOUS = '-';
const char CODE_ERROR = '!';

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
	uint8_t threads;
	uint64_t pageLen;
	int lowEngeryThreshold;
	int highEngeryThreshold;
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
	std::string ssL1 = CODE_UNTESTED;
	std::string ssStructure = CODE_UNTESTED;
	double ssEnergy = NULL;
	double sgrnascorer2score = NULL;
	std::string bowtieChr = CODE_UNTESTED;
	uint64_t bowtieStart = NULL;
	uint64_t bowtieEnd = NULL;
	double mitOfftargetscore = NULL;
	double cfdOfftargetscore = NULL;
};

#endif // !utilInclude