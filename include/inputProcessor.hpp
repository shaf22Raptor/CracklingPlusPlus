// inputProcessor.hpp
#pragma once
#include <list>
#include <string>
#include <stdexcept>

class inputProcessor 
{
public:
	virtual void process(std::list<std::string> const & filesToProcess, int const & batchSize);

	virtual void cleanUp();

	virtual ~inputProcessor() = default;
};