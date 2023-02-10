#include "../include/consensusModule.hpp"

consensusModule::consensusModule(cracklingConfig config)
{
	this->optimsationLevel = config.general.optimisation;
	this->toolCount = config.consensus.toolCount;
	this->consensusN = config.consensus.n;
}

void consensusModule::run() {}