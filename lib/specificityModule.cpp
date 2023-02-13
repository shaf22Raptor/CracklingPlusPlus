#include "../include/specificityModule.hpp"

specificityModule::specificityModule(cracklingConfig config)
{
	this->toolIsSelected = config.offtarget.enabled;
	this->optimsationLevel = config.general.optimisation;
}

void specificityModule::run() {}
