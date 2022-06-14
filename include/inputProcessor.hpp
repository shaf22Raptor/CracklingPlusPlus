// inputProcessor.hpp
#pragma once
#include <list>
#include <string>
#include <stdexcept>

class inputProcessor 
{
public:
	virtual std::list<std::string> processInput(std::list<std::string>& filesToProcess, int batchSize);

};