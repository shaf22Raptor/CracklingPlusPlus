#include <Logger.hpp>

Logger::Logger(std::string fileStream)
{
	outFileBuffer.open(fileStream);
};

int Logger::overflow(int c) {
	putchar(c);
	outFileBuffer.put(char(c));
	outFileBuffer.flush();
	return c;
}