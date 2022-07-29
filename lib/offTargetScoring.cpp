#include <offTargetScoring.hpp>

using std::string;
using std::map;
using std::vector;

offTargetScoring::offTargetScoring(ConfigManager& cm) :
	toolIsSelected(cm.getBool("offtargetscore", "enabled")),
	optimsationLevel(cm.getString("general", "optimisation")),
	toolCount(cm.getConsensusToolCount()),
	consensusN(cm.getInt("consensus", "n")),
	offTargetScoreOutFile(cm.getString("offtargetscore", "output")),
	offTargetScoreInFile(cm.getString("offtargetscore", "input")),
	offTargetScoreBin(cm.getString("offtargetscore", "binary")),
	offTargetScoreIndex(cm.getString("input", "offtarget-sites")),
	offTargetScoreMaxDist(cm.getString("offtargetscore", "max-distance")),
	offTargetScoreMethod(cm.getString("offtargetscore", "method")),
	offTagertScoreThreshold(cm.getFloat("offtargetscore", "score-threshold")),
	offTargetScorePageLength(cm.getInt("offtargetscore", "page-length"))
{}

void offTargetScoring::run(map<string, map<string, string, std::less<>>, std::less<>>& candidateGuides)
{
	if (!toolIsSelected)
	{
		printer("Off-target scoring has been configured not to run. Skipping Off-target scoring");
		return;
	}

	printer("Beginning Off-target scoring.");
	int testedCount = 0;
	int failedCount = 0;
	int pgIdx = 1;
	int guidesInPage = 0;
	auto paginatorIterator = candidateGuides.begin();
	auto pageStart = candidateGuides.begin();
	auto pageEnd = candidateGuides.begin();

	// Outer loop deals with changing iterator start and end points (Pagination)
	while (pageEnd != candidateGuides.end())
	{
		if (offTargetScorePageLength > 0)
		{
			// Advance the pageEnd pointer
			std::advance(pageEnd, std::min((int)std::distance(pageEnd, candidateGuides.end()), offTargetScorePageLength));
			// Record page start
			pageStart = paginatorIterator;
			// Print page information
			printer(fmt::format("\tProcessing page {} ({} per page).", commaify(pgIdx), commaify(offTargetScorePageLength)));
		}
		else {
			// Process all guides at once
			pageEnd = candidateGuides.end();
		}


		printer("\tConstructing the Off-target scoring input file.");
		// Open input file (Opening in binary to avoid CR+LF on windows)
		std::ofstream inFile;
		inFile.open(offTargetScoreInFile, std::ios_base::binary);

		guidesInPage = 0;
		while (paginatorIterator != pageEnd)
		{
			string target23 = paginatorIterator->first;
			// Run time filtering
			if (!filterCandidateGuides(paginatorIterator->second, MODULE_SPECIFICITY, optimsationLevel, consensusN, toolCount)) {
				// Advance page end for each filtered out guide
				if (pageEnd != candidateGuides.end())
				{
					pageEnd++;
				}
				paginatorIterator++;
				continue;
			}

			inFile << target23.substr(0, 20) << "\n";

			guidesInPage++;
			paginatorIterator++;
		}

		inFile.close();

		printer(fmt::format("\t\t{} guides in this page.", commaify(guidesInPage)));

		// Call scoring method
		runner(fmt::format("{} {} {} {} {} {} > {}",
			offTargetScoreBin,
			offTargetScoreIndex,
			offTargetScoreInFile,
			offTargetScoreMaxDist,
			offTagertScoreThreshold,
			offTargetScoreMethod,
			offTargetScoreOutFile
		).c_str());

		printer("\tStarting to process the Off-target scoring results.");

		map<string, map<string, string, std::less<>>, std::less<>> targetsScored;

		// Open output file 
		std::ifstream outFile;
		outFile.open(offTargetScoreOutFile, std::ios::binary);


		for (string line; std::getline(outFile, line); )
		{
			vector<string> splitLine = split(line, "\t");
			targetsScored[splitLine[0]]["MIT"] = trim(splitLine[1]);
			targetsScored[splitLine[0]]["CFD"] = trim(splitLine[2]);
		}

		// Reset paginatorIterator to page start
		paginatorIterator = pageStart;

		while (paginatorIterator != pageEnd)
		{
			string target23 = paginatorIterator->first;

			if (string target = target23.substr(0, 20); targetsScored.find(target) != targetsScored.end())
			{
				candidateGuides[target23]["mitOfftargetscore"] = targetsScored[target]["MIT"];
				candidateGuides[target23]["cfdOfftargetscore"] = targetsScored[target]["CFD"];

				if (offTargetScoreMethod == "mit")
				{
					if (std::stof(targetsScored[target]["MIT"]) < offTagertScoreThreshold)
					{
						candidateGuides[target23]["passedOffTargetScore"] = CODE_REJECTED;
						failedCount++;
					}
					else { candidateGuides[target23]["passedOffTargetScore"] = CODE_ACCEPTED; }
				}

				else if (offTargetScoreMethod == "cfd")
				{
					if (std::stof(targetsScored[target]["CFD"]) < offTagertScoreThreshold)
					{
						candidateGuides[target23]["passedOffTargetScore"] = CODE_REJECTED;
						failedCount++;
					}
					else { candidateGuides[target23]["passedOffTargetScore"] = CODE_ACCEPTED; }
				}

				else if (offTargetScoreMethod == "and")
				{
					if ((std::stof(targetsScored[target]["MIT"]) < offTagertScoreThreshold) && (std::stof(targetsScored[target]["CFD"]) < offTagertScoreThreshold))
					{
						candidateGuides[target23]["passedOffTargetScore"] = CODE_REJECTED;
						failedCount++;
					}
					else { candidateGuides[target23]["passedOffTargetScore"] = CODE_ACCEPTED; }
				}

				else if (offTargetScoreMethod == "or")
				{
					if ((std::stof(targetsScored[target]["MIT"]) < offTagertScoreThreshold) || (std::stof(targetsScored[target]["CFD"]) < offTagertScoreThreshold))
					{
						candidateGuides[target23]["passedOffTargetScore"] = CODE_REJECTED;
						failedCount++;
					}
					else { candidateGuides[target23]["passedOffTargetScore"] = CODE_ACCEPTED; }
				}

				else if (offTargetScoreMethod == "avg")
				{
					if (((std::stof(targetsScored[target]["MIT"]) + std::stof(targetsScored[target]["MIT"])) / 2) < offTagertScoreThreshold)
					{
						candidateGuides[target23]["passedOffTargetScore"] = CODE_REJECTED;
						failedCount++;
					}
					else { candidateGuides[target23]["passedOffTargetScore"] = CODE_ACCEPTED; }
				}
				testedCount++;
			}
			paginatorIterator++;
		}
		printer(fmt::format("\t{} of {} failed here.", commaify(failedCount), commaify(testedCount)));
	}

}