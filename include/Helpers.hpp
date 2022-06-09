// Helpers.hpp
#pragma once
#include <iostream>
#include <string>
#include <algorithm>
#include <iterator>
#include <array>
#include <ctime>

std::string rc(std::string DNA);

float atPercentage(std::string seq);

void printer(std::string formattedString);

void errPrinter(std::string formattedString);

void runner(char* args);