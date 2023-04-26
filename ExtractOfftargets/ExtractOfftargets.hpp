#include <iostream>
#include <fstream>
#include <filesystem>
#include <string>
#include <atomic>
#include <omp.h>
#define FMT_HEADER_ONLY
#include "../include/fmt/format.h"
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/regex.hpp>
#include "../include/util.hpp"
#include <chrono>

const boost::regex fwdExp = boost::regex("(?=([ACG][ACGT]{19}[ACGT][AG]G))");
const boost::regex bwdExp = boost::regex("(?=(C[CT][ACGT][ACGT]{19}[TGC]))");
const std::string END = "ZZZZZZZZZZZZZZZZZZZZ";