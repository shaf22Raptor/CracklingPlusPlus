// Logger.hpp
#pragma once
#include <iostream>
#include <streambuf>
#include <fstream>
#include <cstring>

class logBuffer : public std::streambuf
{
public:
    logBuffer();

    logBuffer(const std::ostream& stdOut ,const std::ostream& destFile);

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


