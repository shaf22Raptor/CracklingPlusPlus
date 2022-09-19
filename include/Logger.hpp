// Logger.hpp
#pragma once
#define _POSIX_C_SOURCE 1
#include <iostream>
#include <streambuf>
#include <fstream>
#include <cstring>
#include <ctime>


#if defined(_WIN64)
    # define p_localtime(time_t_ptr, tm_ptr) localtime_s(&time_t_ptr, &tm_ptr)
#elif defined(unix) || defined(__unix__) || defined(__unix)
    # define p_localtime(time_t_ptr, tm_ptr) localtime_r(&tm_ptr, &time_t_ptr)
#else
    # error "Error, no localtime function"
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


