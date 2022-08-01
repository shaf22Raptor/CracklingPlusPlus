// inputProcessor class, base class for processing input files for various Cas proteins
#include "../include/inputProcessor.hpp"
using std::string;
using std::list;
using std::logic_error;

void inputProcessor::process(list<string> const & filesToProcess, int const & batchSize)
{
	throw NotImplemented();
}

void inputProcessor::cleanUp()
{
	throw NotImplemented();
}	