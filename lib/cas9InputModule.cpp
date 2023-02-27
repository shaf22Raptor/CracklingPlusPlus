#include "../include/cas9InputModule.hpp"

using boost::regex;

cas9InputModule::cas9InputModule(cracklingConfig config) : inputModule(config)
{
	this->fwdExp = regex("(?=([ATCG]{21}GG))");
	this->bwdExp = regex("(?=(CC[ACGT]{21}))");
}