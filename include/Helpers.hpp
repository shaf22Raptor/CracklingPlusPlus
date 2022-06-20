// Helpers.hpp
#pragma once
#include <iostream>
#include <string>
#include <algorithm>
#include <iterator>
#include <array>
#include <ctime>
#include <map>
#include <vector>
#include <Constants.hpp>

std::vector<std::string> split(std::string s, std::string delimiter);

std::string rtrim(const std::string& s);

std::string ltrim(const std::string& s);

std::string trim(const std::string& s);

std::string rc(std::string DNA);

bool filterCandidateGuides(std::map<std::string, std::string> candidateGuideResultMap, std::string selectedModule, std::string optimisation, int consensusN, int toolCount);

void printer(std::string formattedString);

void errPrinter(std::string formattedString);

void runner(char* args);