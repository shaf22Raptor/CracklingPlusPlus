#include <offTargetScoring.hpp>

using std::string;
using std::map;
using std::vector;

offTargetScoring::offTargetScoring(ConfigManager cm) :
	toolIsSelected(false),
	optimsationLevel(""),
	toolCount(0),
	consensusN(0),
	offTargetScoreOutFile(""),
	offTargetScoreInFile(""),
	offTargetScoreBin(""),
	offTargetScoreIndex(""),
	offTargetScoreMaxDist(""),
	offTargetScoreMethod(""),
	offTagertScoreThreshold(0.0f),
	offTargetScorePageLength(0)
{
	toolIsSelected = cm.getBool("offtargetscore", "enabled");
	optimsationLevel = cm.getString("general", "optimisation");
	toolCount = cm.getConsensusToolCount();
	consensusN = cm.getInt("consensus", "n");
	offTargetScoreOutFile = cm.getString("offtargetscore", "output");
	offTargetScoreInFile = cm.getString("offtargetscore", "input");
	offTargetScoreBin = cm.getString("offtargetscore", "binary");
	offTargetScoreIndex = cm.getString("input", "offtarget-sites");
	offTargetScoreMaxDist = cm.getString("offtargetscore", "max-distance");
	offTargetScoreMethod = cm.getString("offtargetscore", "method");
	offTagertScoreThreshold = cm.getFloat("offtargetscore", "score-threshold");
	offTargetScorePageLength = cm.getInt("offtargetscore", "page-length");
}

void offTargetScoring::run(map<string, map<string, string>>& candidateGuides)
{
	if (!toolIsSelected)
	{
		printer("Off-target scoring has been configured not to run. Skipping Off-target scoring");
		return;
	}

	printer("Beginning Off-target scoring.");
	char printingBuffer[1024];
	int testedCount = 0;
	int failedCount = 0;
	int pgIdx = 1;
	int guidesInPage = 0;
	map<string, map<string, string>>::iterator paginatorIterator = candidateGuides.begin();
	map<string, map<string, string>>::iterator pageStart = candidateGuides.begin();
	map<string, map<string, string>>::iterator pageEnd = candidateGuides.begin();

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
			snprintf(printingBuffer, 1024, "\tProcessing page %d (%d per page).", pgIdx, offTargetScorePageLength);
			printer(printingBuffer);
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
		for (paginatorIterator; paginatorIterator != pageEnd; paginatorIterator++)
		{
			string target23 = paginatorIterator->first;
			map<string, string> resultsMap = paginatorIterator->second;
			// Run time filtering
			if (!filterCandidateGuides(resultsMap, MODULE_SPECIFICITY, optimsationLevel, consensusN, toolCount)) {
				// Advance page end for each filtered out guide
				if (pageEnd != candidateGuides.end())
				{
					pageEnd++;
				}
				continue;
			}

			inFile << target23.substr(0, 20) << "\n";

			guidesInPage++;
		}

		inFile.close();

		snprintf(printingBuffer, 1024, "\t\t%d guides in this page.", guidesInPage);
		printer(printingBuffer);

		// Call scoring method
		snprintf(printingBuffer, 1024, "%s %s %s %s %f %s > %s",
			offTargetScoreBin.c_str(),
			offTargetScoreIndex.c_str(),
			offTargetScoreInFile.c_str(),
			offTargetScoreMaxDist.c_str(),
			offTagertScoreThreshold,
			offTargetScoreMethod.c_str(),
			offTargetScoreOutFile.c_str()
		);
		runner(printingBuffer);

		printer("\tStarting to process the Off-target scoring results.");

		map<string, map<string, string>> targetsScored;

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

		for (paginatorIterator; paginatorIterator != pageEnd; paginatorIterator++)
		{
			string target23 = paginatorIterator->first;
			string target = target23.substr(0, 20);

			if (targetsScored.find(target) != targetsScored.end())
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
					if  ( ((std::stof(targetsScored[target]["MIT"]) + std::stof(targetsScored[target]["MIT"]))/2) < offTagertScoreThreshold ) 
					{
						candidateGuides[target23]["passedOffTargetScore"] = CODE_REJECTED;
						failedCount++;
					}
					else { candidateGuides[target23]["passedOffTargetScore"] = CODE_ACCEPTED; }
				}
				testedCount++;
			}
		}
		snprintf(printingBuffer, 1024, "\t%d of %d failed here.", failedCount, testedCount);
		printer(printingBuffer);
	}

}