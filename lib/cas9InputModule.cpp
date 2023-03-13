#include "../include/cas9InputModule.hpp"

using boost::regex;

cas9InputModule::cas9InputModule(const cracklingConfig& config) : inputModule(config)
{
	this->fwdExp = regex("(?=([ATCG]{21}GG))");
	this->bwdExp = regex("(?=(CC[ACGT]{21}))");
}

void cas9InputModule::run()
{
	inputModule::run();
}

std::vector<guideResults>* cas9InputModule::next()
{
	return inputModule::next();
}