// inputProcessor class, base class for processing input files for various Cas proteins
#include <inputProcessor.hpp>
using std::string;
using std::list;
using std::logic_error;

class NotImplemented : public std::logic_error
{
public:
	NotImplemented() : std::logic_error("Function not yet implemented") { };
};

void inputProcessor::process(list<string> const & filesToProcess, int const & batchSize)
{
	throw NotImplemented();
}

void inputProcessor::cleanUp()
{
	throw NotImplemented();
}	