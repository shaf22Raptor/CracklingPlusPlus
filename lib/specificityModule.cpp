#include "../include/specificityModule.hpp"

specificityModule::specificityModule(cracklingConfig config)
{
	this->toolIsSelected = config.offTarget.enabled;
	this->optimsationLevel = config.general.optimisation;
	this->consensusN = config.consensus.n;
}

void specificityModule::run() {}
