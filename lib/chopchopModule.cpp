#include "../include/chopchopModule.hpp"

using std::string_view;
using std::cout;
using std::endl;

chopchopModule::chopchopModule(cracklingConfig config) : consensusModule(config)
{
	this->toolIsSelected = config.consensus.chopchop;
}

void chopchopModule::run(std::vector<guideResults>& candidateGuides)
{
    if (!toolIsSelected)
    {
        cout << "CHOPCHOP has been configured not to run. Skipping CHOPCHOP" << endl;
        return;
    }


    cout << "CHOPCHOP - remove those without G in position 20." << endl;
    uint64_t failedCount = 0;
    uint64_t testedCount = 0;
    for (guideResults& candidate : candidateGuides)
    {
        // Run time filtering
        if (!processGuide(candidate)) { continue; }

        if (G20(candidate.seq))
        {
            candidate.passedG20 = CODE_ACCEPTED;
        }
        else
        {
            candidate.passedG20 = CODE_REJECTED;
            failedCount++;
        }
        testedCount++;
    }

    cout << fmt::format("\t{} of {} failed here.", failedCount, testedCount) << endl;
    return;
}


bool chopchopModule::processGuide(const guideResults& guide)
{
    // Process all guides at this level
    if (optimsationLevel == optimisationLevel::ultralow) { return true; }

    // For all levels above `ultralow`
    if (!guide.isUnique)
    {
        // Reject all guides that have been seen more than once
        return false;
    }

    // None of the failure conditions have been meet, return true
    return true;
}

bool chopchopModule::G20(string_view guide)
{
	return guide[19] == 'G';
}