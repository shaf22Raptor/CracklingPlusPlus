#ifndef inputModuleInclude
#define inputModuleInclude
#include <string>
#include <vector>
#include <fstream>
#include <filesystem>
#include <boost/algorithm/string.hpp>
#include <boost/regex.hpp>
#include "../include/util.hpp"
#include "../include/pipelineModule.hpp"
#include "../include/phmap/phmap.h"

class inputModule : private pipelineModule
{
public:
	inputModule(const cracklingConfig& config);
	void run();
	std::vector<guideResults>* next();
protected:
	std::vector<guideResults> guideBatch;
	std::filesystem::path tempWorkingDir;
	boost::regex fwdExp;
	boost::regex bwdExp;
	uint64_t batchSize;
	uint64_t guidesInBatch;
	uint64_t numDuplicateGuides;
	uint64_t numIdentifiedGuides;
	phmap::flat_hash_set<std::string> candidateGuides;
	phmap::flat_hash_set<std::string> recordedSequences;
	phmap::flat_hash_set<std::string> duplicateGuides;
	std::fstream currentBatchFile;
	std::vector<std::filesystem::path> batchFiles;
	std::vector<std::filesystem::path> filesToProcess;
	void processSeqeunce(const std::string& seqeunce, const std::string& header);
};


#endif // !inputModuleInclude
