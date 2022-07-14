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

#include <format>

class cas9InputProcessor : public inputProcessor
{
public:
	void process(std::list<std::string> const & filesToProcess, int const & batchSize) final;

	const std::list<std::string>& getBatchFiles() const;

	bool isDuplicateGuide(std::string_view guide) const ;

	void cleanUp() final;

	~cas9InputProcessor() final = default;

private:
	std::set<std::string, std::less<>> duplicateGuides;
	std::list<std::string> batchFiles;
	int numDuplicateGuides;
	int numIdentifiedGuides;
	int guidesInBatch;

	void processSeqeunce(
		std::string_view seqeunce,
		std::string_view seqHeader,
		std::ofstream& outFile,
		std::filesystem::path const& tempWorkingDir,
		std::set<std::string, std::less<>>& candidateGuides,
		const int& batchSize
	);
};

class FileSystemError : public std::runtime_error
{
public:
	FileSystemError() : std::runtime_error("Unable to create temp working dir") { };
};
