// ISSLCreateIndex.hpp
/*

Faster and better CRISPR guide RNA design with the Crackling method.
Jacob Bradford, Timothy Chappell, Dimitri Perrin
bioRxiv 2020.02.14.950261; doi: https://doi.org/10.1101/2020.02.14.950261


To compile:

g++ -o isslCreateIndex isslCreateIndex.cpp -O3 -std=c++11 -fopenmp -mpopcnt

*/

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <string>
#include <string_view>
#include <sys/types.h>
#include <sys/stat.h>
#include <map>
#include <regex>

#ifndef portableStat64
#define portableStat64
#if defined(_WIN64)
#include "../include/unistd.h"
#define p_stat64 _stat64
#elif defined(unix) || defined(__unix__) || defined(__unix)
#include <unistd.h>
#define p_stat64 stat64
#else
# error "Error, no stat function"
#endif
#endif // !portableStat64