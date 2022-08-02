// Logger.hpp
#pragma once
#include <iostream>
#include <streambuf>
#include <fstream>
#include <cstring>
#include <ctime>

#if (_POSIX_C_SOURCE >= 1 || _XOPEN_SOURCE || _BSD_SOURCE || _SVID_SOURCE || _POSIX_SOURCE)
# define p_localtime(time_t_ptr, tm_ptr) localtime_r(&tm_ptr, &time_t_ptr)
#elif defined(_MSC_VER)
# define p_localtime(time_t_ptr, tm_ptr) localtime_s(&time_t_ptr, &tm_ptr)
#else
# error "Error, no localtime"
#endif

class logBuffer : public std::streambuf
{
public:
    logBuffer();

    logBuffer(const std::ostream& stdOut, const std::ostream& destFile);

    int sync() final;

    virtual int underflow(int c);

    int overflow(int c) final;

private:
    bool isAtStartOfLine = true;
    std::streambuf* fileBuf = nullptr;
    std::streambuf* stdBuf = nullptr;
};

class Logger : public std::streambuf
{
public:
    Logger(std::ostream& logSource, const std::string& outFile);

    void close();

private:
    std::ostream* origOutputStream;
    std::streambuf* origSrcBuffer;
    std::ofstream outputFileStream;
    logBuffer customBuffer;
};


