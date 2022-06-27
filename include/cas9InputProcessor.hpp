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
	std::set<std::string> duplicateGuides;

	std::list<std::string> processInput(std::list<std::string> filesToProcess, int batchSize);
};