#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <string>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <map>

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