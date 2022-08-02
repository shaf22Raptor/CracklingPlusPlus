// cas9InputProcessor.hpp
#pragma once
#include <list>
#include <string>
#include <fstream>
#include <unordered_set>
#include <regex>
#include <filesystem>
#include <numeric>
#include "../include/inputProcessor.hpp"
#include "../include/Helpers.hpp"

class cas9InputProcessor : public inputProcessor
{
public:
	void process(std::list<std::string> const & filesToProcess, int const & batchSize) final;

	const std::list<std::string>& getBatchFiles() const;

	bool isDuplicateGuide(const std::string& guide) const ;

	void cleanUp() final;

	~cas9InputProcessor() final = default;

private:
	std::unordered_set<std::string> duplicateGuides;
	std::list<std::string> batchFiles;
	int numDuplicateGuides;
	int numIdentifiedGuides;
	int guidesInBatch;

	void processSeqeunce(
		std::string_view seqeunce,
		std::string_view seqHeader,
		std::ofstream& outFile,
		std::filesystem::path const& tempWorkingDir,
		std::unordered_set<std::string>& candidateGuides,
		const int& batchSize
	);
};

class FileSystemError : public std::runtime_error
{
public:
	FileSystemError() : std::runtime_error("Unable to create temp working dir") { };
};
