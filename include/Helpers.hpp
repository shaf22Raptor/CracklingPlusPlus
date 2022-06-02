// Helpers.hpp
#pragma once
#include <iostream>
#include <string>
#include <algorithm>
#include <iterator>
#include <array>
#include <ctime>

std::string rc(std::string DNA);

std::string transToDNA(std::string RNA);

float atPercentage(std::string seq);

void printer(std::string formattedString);

void runner(char* args);