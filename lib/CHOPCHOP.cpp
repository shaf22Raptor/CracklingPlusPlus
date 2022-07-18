// CHOPCHOP class
#include <CHOPCHOP.hpp>

using std::string;
using std::map;

CHOPCHOP::CHOPCHOP(ConfigManager& cm) :
    testedCount(),
    failedCount(),
    toolIsSelected(cm.getBool("consensus", "chopchop")),
    optimsationLevel(cm.getString("general", "optimisation")),
    toolCount(cm.getConsensusToolCount()),
    consensusN(cm.getInt("consensus", "n"))
{}

void CHOPCHOP::run(map<string, map<string, string, std::less<>>, std::less<>>& candidateGuides)
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

    printer(fmt::format("\t{} of {} failed here.", commaify(failedCount), commaify(testedCount)));
    return;
}

bool CHOPCHOP::G20(std::string_view candidateGuide)
{
    if (candidateGuide.length() < 20)
    {
        throw G20Input();
    }
	return candidateGuide[19] == 'G';
}