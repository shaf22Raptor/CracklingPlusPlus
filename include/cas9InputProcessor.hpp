// cas9InputProcessor.hpp
#pragma once
#include <list>
#include <string>
#include <fstream>
#include <set>
#include <regex>
#include <filesystem>
#include <inputProcessor.hpp>
#include <Helpers.hpp>

class cas9InputProcessor : public inputProcessor
{
public:
	void processInput(std::list<std::string> filesToProcess, int batchSize);

	std::list<std::string> getBatchFiles();

	bool isDuplicateGuide(std::string guide);

	void cleanUp();

private:
	std::set<std::string> duplicateGuides;
	std::list<std::string> batchFiles;

	void processSeqeunce(
		const std::string& seqeunce,
		const std::string& seqHeader,
		std::ofstream& outFile,
		std::filesystem::path& tempWorkingDir,
		int& numIdentifiedGuides,
		int& numDuplicateGuides,
		std::set<std::string>& candidateGuides,
		std::set<std::string>& recordedSequences,
		int& guidesInBatch,
		const int& batchSize
	);
};