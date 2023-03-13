#include "../include/chopchopModule.hpp"

using std::string_view;
using std::cout;
using std::endl;

chopchopModule::chopchopModule(const cracklingConfig& config) : consensusModule(config)
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

bool chopchopModule::G20(string_view guide)
{
	return guide[19] == 'G';
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

    // For optimisation level `high`
    if (optimsationLevel = optimisationLevel::high)
    {
        int countAlreadyAccepted =
            (int)guide.passedG20 == CODE_ACCEPTED +
            (int)guide.acceptedByMm10db == CODE_ACCEPTED +
            (int)guide.acceptedBySgRnaScorer2 == CODE_ACCEPTED;

        int countAlreadyAssessed =
            (int)guide.passedG20 != CODE_UNTESTED +
            (int)guide.acceptedByMm10db != CODE_UNTESTED +
            (int)guide.acceptedBySgRnaScorer2 != CODE_UNTESTED;

        // Reject if the consensus has already been passed
        if (countAlreadyAccepted >= consensusN) { return false; }

        // Reject if there is not enough tests remaining to pass consensus
        if (toolCount - countAlreadyAssessed < consensusN - countAlreadyAccepted) { return false; }
    }

    // None of the failure conditions have been meet, return true
    return true;
}