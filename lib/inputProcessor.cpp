// inputProcessor class, base class for processing input files for various Cas proteins
#include "../include/inputProcessor.hpp"
using std::string;
using std::list;
using std::logic_error;

void inputProcessor::process(const list<string>& filesToProcess, const int& batchSize)
{
	throw NotImplemented();
}

void inputProcessor::cleanUp()
{
	throw NotImplemented();
}