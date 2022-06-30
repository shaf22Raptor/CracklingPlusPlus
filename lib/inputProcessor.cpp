// inputProcessor class, base class for processing input files for various Cas proteins
#include <inputProcessor.hpp>
using std::string;
using std::list;
using std::logic_error;

list<string> inputProcessor::processInput(list<string>& filesToProcess, int batchSize)
{
	throw logic_error("Function not yet implemented");
}

void inputProcessor::cleanUp()
{
	throw logic_error("Function not yet implemented");
}