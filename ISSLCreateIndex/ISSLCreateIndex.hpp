// ISSLCreateIndex.hpp
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <regex>
#include <filesystem>
#include <fstream>
#define FMT_HEADER_ONLY
#include "../include/fmt/format.h"
#include "../include/libpopcnt/libpopcnt.h"

const std::regex extractNumbers("[1234567890]+");
const std::vector<uint8_t> nucleotideIndex{ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,3 };
const std::vector<char> signatureIndex{ 'A', 'C', 'G', 'T' };
// 3bit
// const vector<uint8_t> nucleotideIndex{ 1,0,2,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,7 };
// const vector<char> signatureIndex{ '0','A','C','3','G','5','6','T' };
uint64_t seqLength;