// Helpers.h
#pragma once
#include <iostream>
#include <string>
#include <algorithm>
#include <iterator>
#include <array>
#include <ctime>

using std::string;
using std::array;

string rc(string DNA);

string transToDNA(string RNA);

float atPercentage(string seq);

void printer(string formattedString);