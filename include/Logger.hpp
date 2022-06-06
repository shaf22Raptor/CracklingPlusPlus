// Logger.hpp
#pragma once
#include <iostream>
#include <streambuf>
#include <fstream>

class logBuffer : public std::streambuf
{
public:
    logBuffer() = default;

    logBuffer(std::ostream& stdOut ,std::ostream& destFile);

    virtual int sync();

    virtual int underflow(int c);

    virtual int overflow(int c);

private:
    bool isAtStartOfLine;
    std::streambuf* fileBuf;
    std::streambuf* stdBuf;
};

class Logger : public std::streambuf
{
public:
    Logger(std::ostream& logSource, std::string outFile);

    void close();

private:
    std::ostream* origOutputStream;
    std::streambuf* origSrcBuffer;
    std::ofstream outputFileStream;
    logBuffer customBuffer;
    std::ostream* loggedOutputStream;
};


