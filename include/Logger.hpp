// Logger.hpp
#pragma once
#include <streambuf>
#include <fstream>

class Logger : public std::streambuf
{
public:

	Logger(std::string fileStream);

	virtual int overflow(int c);

private:
	std::ofstream outFileBuffer;

};