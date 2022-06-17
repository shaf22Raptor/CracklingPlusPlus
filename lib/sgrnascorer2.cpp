#include <sgrnascorer2.hpp>


using std::string;
using std::map;

sgrnascorer2::sgrnascorer2(ConfigManager cm) :
	toolIsSelected(false),
	optimsationLevel(""),
	toolCount(0),
	consensusN(0),
	scoreThreshold(0.0f),
	printingBuffer{ "\0" }
{
	toolIsSelected = cm.getBool("consensus", "mm10db");
	optimsationLevel = cm.getString("general", "optimisation");
	toolCount = cm.getConsensusToolCount();
	consensusN = cm.getInt("consensus", "n");
	scoreThreshold = cm.getFloat("sgrnascorer2", "score-threshold");
}

void sgrnascorer2::run(map<string, map<string, string>>& candidateGuides)
{

	if (!toolIsSelected)
	{
		printer("sgRNAScorer2 has been configured not to run. Skipping sgRNAScorer2");
		return;
	}

	printer("sgRNAScorer2 - score using model.");
	int failedCount = 0;
	int testedCount = 0;
	for (auto const& [target23, resultsMap] : candidateGuides)
	{

		// Run time filtering
		if (!filterCandidateGuides(resultsMap, MODULE_SGRNASCORER2, optimsationLevel, consensusN, toolCount)) { continue; }


	}

}