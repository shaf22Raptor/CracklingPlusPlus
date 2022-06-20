// CHOPCHOP class
#include <CHOPCHOP.hpp>

using std::string;
using std::map;

CHOPCHOP::CHOPCHOP(ConfigManager cm) : 
    failedCount(0),
    testedCount(0),
    toolIsSelected(false),
    optimsationLevel("ultralow"),
    toolCount(0),
    consensusN(0),
    printingBuffer{"\0"}
{
    toolIsSelected = cm.getBool("consensus", "chopchop");
    optimsationLevel = cm.getString("general", "optimisation");
    toolCount = cm.getConsensusToolCount();
    consensusN = cm.getInt("consensus", "n");
}

void CHOPCHOP::run(std::map<std::string, std::map<std::string, std::string>>& candidateGuides)
{
    if (!toolIsSelected)
    {
        printer("CHOPCHOP has been configured not to run. Skipping CHOPCHOP");
        return;
    }

    printer("CHOPCHOP - remove those without G in position 20.");
    failedCount = 0;
    testedCount = 0;
    for (auto const& [target23, resultsMap] : candidateGuides)
    {
        // Run time filtering
        if (!filterCandidateGuides(resultsMap, MODULE_CHOPCHOP, optimsationLevel, consensusN, toolCount)) { continue; }

        if (G20(target23))
        {
            candidateGuides[target23]["passedG20"] = CODE_ACCEPTED;
        }
        else
        {
            candidateGuides[target23]["passedG20"] = CODE_REJECTED;
            failedCount++;
        }
        testedCount++;
    }

    snprintf(printingBuffer, 1024, "\t%d of %d failed here.", failedCount, testedCount);
    printer(printingBuffer);
    return;
}

bool CHOPCHOP::G20(std::string candidateGuide)
{
	return candidateGuide[19] == 'G' || candidateGuide[19] == 'g';
}