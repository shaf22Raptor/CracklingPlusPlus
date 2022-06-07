// CHOPCHOP class
#include <CHOPCHOP.hpp>

using std::string;
using std::map;

CHOPCHOP::CHOPCHOP(ConfigManager cm) : failedcount(0), testedcount(0), toolIsSelected(false)
{
    toolIsSelected = cm.getBool("consensus", "chopchop");
}

void CHOPCHOP::run(std::map<std::string, std::map<std::string, std::string>> candidateGuides)
{
    if (!toolIsSelected)
    {
        return;
    }

    failedcount = 0;
    testedcount = 0;
    for (auto const& [target23, resultsMap] : candidateGuides)
    {
        if (G20(target23))
        {
            candidateGuides[target23]["passedG20"] = CODE_ACCEPTED;
        }
        else
        {
            candidateGuides[target23]["passedG20"] = CODE_REJECTED;
            failedcount++;
        }
        testedcount++;
    }


    snprintf(printingBuffer, 1024, "%d of %d failed here.", failedcount, testedcount);
    printer(printingBuffer);
    return;
}

bool CHOPCHOP::G20(std::string candidateGuide)
{
	return candidateGuide[19] == 'G' || candidateGuide[19] == 'g';
}