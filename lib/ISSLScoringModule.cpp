#include "ISSLScoringModule.hpp"

ISSLScoringModule::ISSLScoringModule(cracklingConfig config) : specificityModule(config)
{

}

void ISSLScoringModule::run(std::vector<guideResults>& candidateGuides)
{
}

bool ISSLScoringModule::processGuide(const guideResults& guide)
{
	// Process all guides at this level
	if (optimsationLevel == optimisationLevel::ultralow) 
	{ 
		// Reject if off target scoring has failed
		if (guide.pass candidateGuideResultMap["passedOffTargetScore"] == CODE_REJECTED) { return false; }
		else { return true; }
	}

	// For all levels above `ultralow`
	if (!guide.isUnique)
	{
		// Reject all guides that have been seen more than once
		return false;
	}

	// For optimisation levels `medium` and `high`
	if (optimsationLevel == optimisationLevel::medium || optimsationLevel == optimisationLevel::high)
	{
		// Reject if the consensus was not passed
		if (stoi(candidateGuideResultMap["consensusCount"]) < consensusN) { return false; }
		// Reject if bowtie2 was not passed
		if (candidateGuideResultMap["passedBowtie"] == CODE_REJECTED) { return false; }
	}

	// None of the failure conditions have been meet, return true
	return true;
}
