// mm10db class
#include <mm10db.hpp>

using std::string;
using std::map;
using std::array;
using std::list;

array<char, 4> atArray = { 'a', 't', 'A', 'T' };

mm10db::mm10db(ConfigManager& cm) : 
	toolIsSelected(false), 
	optimsationLevel(""), 
	toolCount(0), 
	consensusN(0), 
	RNAFoldOutFile(""), 
	RNAFoldInFile(""), 
	RNAFoldBin(""),
	RNAFoldPageLength(0),
	lowEnergyThreshold(0.0f),
	highEnergyThreshold(0.0f),
	printingBuffer{"\0"}
{
	toolIsSelected = cm.getBool("consensus", "mm10db");
	optimsationLevel = cm.getString("general", "optimisation");
	toolCount = cm.getConsensusToolCount();
	consensusN = cm.getInt("consensus", "n");
	RNAFoldOutFile = cm.getString("rnafold", "output");
	RNAFoldInFile = cm.getString("rnafold", "input");
	RNAFoldBin = cm.getString("rnafold", "binary");
	RNAFoldPageLength = cm.getInt("rnafold", "page-length");
	lowEnergyThreshold = cm.getFloat("rnafold", "low_energy_threshold");
	highEnergyThreshold = cm.getFloat("rnafold", "high_energy_threshold");
}

void mm10db::run(std::map<std::string, std::map<std::string, std::string>>& candidateGuides)
{
	if (!toolIsSelected)
	{
		printer("mm10db has been configured not to run. Skipping mm10db");
		return;
	}

	printer("mm10db - remove all targets with a leading T(+) or trailing A(-).");
	int failedCount = 0;
	int testedCount = 0;
	for (auto const& [target23, resultsMap] : candidateGuides)
	{
		// Run time filtering
		if (!filterCandidateGuides(resultsMap, MODULE_MM10DB, optimsationLevel, consensusN, toolCount)) { continue; }

		if (leadingT(target23))
		{
			candidateGuides[target23]["passedAvoidLeadingT"] = CODE_REJECTED;
			failedCount++;
		}
		else
		{
			candidateGuides[target23]["passedAvoidLeadingT"] = CODE_ACCEPTED;
			
		}
		testedCount++;
	}
	snprintf(printingBuffer, 1024, "%d of %d failed here.", failedCount, testedCount);
	printer(printingBuffer);

	printer("mm10db - remove based on AT percent.");
	failedCount = 0;
	testedCount = 0;
	for (auto const& [target23, resultsMap] : candidateGuides)
	{
		// Run time filtering
		if (!filterCandidateGuides(resultsMap, MODULE_MM10DB, optimsationLevel, consensusN, toolCount)) { continue; }

		float AT = AT_percentage(target23);
		if (AT < 20.0f || AT > 65.0f)
		{
			candidateGuides[target23]["passedATPercent"] = CODE_REJECTED;
			failedCount++;
		}
		else
		{
			candidateGuides[target23]["passedATPercent"] = CODE_ACCEPTED;
		}
		candidateGuides[target23]["AT"] = std::to_string(AT);
		testedCount++;
	}
	snprintf(printingBuffer, 1024, "%d of %d failed here.", failedCount, testedCount);
	printer(printingBuffer);

	printer("mm10db - remove all targets that contain TTTT.");
	failedCount = 0;
	testedCount = 0;
	for (auto const& [target23, resultsMap] : candidateGuides)
	{
		// Run time filtering
		if (!filterCandidateGuides(resultsMap, MODULE_MM10DB, optimsationLevel, consensusN, toolCount)) { continue; }

		if (polyT(target23))
		{
			candidateGuides[target23]["passedTTTT"] = CODE_REJECTED;
			failedCount++;
		}
		else
		{
			candidateGuides[target23]["passedTTTT"] = CODE_ACCEPTED;
		}
		testedCount++;
	}
	snprintf(printingBuffer, 1024, "%d of %d failed here.", failedCount, testedCount);
	printer(printingBuffer);

	printer("mm10db - check secondary structure.");
	string guide = "GUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGCUUUU";
	std::regex pattern_RNAstructure(".{28}\\({4}\\.{4}\\){4}\\.{3}\\){4}.{21}\\({4}\\.{4}\\){4}\\({7}\\.{3}\\){7}\\.{3}\\s\\((.+)\\)");
	std::regex pattern_RNAenergy("\\s\\((.+)\\)");
	testedCount = 0;
	failedCount = 0;
	int errorCount = 0;
	int notFoundCount = 0;
	int guidesInPage = 0;
	int pgIdx = 1;
	map<string, map<string, string>>::iterator paginatorIterator = candidateGuides.begin();
	map<string, map<string, string>>::iterator pageStart = candidateGuides.begin();
	map<string, map<string, string>>::iterator pageEnd = candidateGuides.begin();
	

	// Outer loop deals with changing iterator start and end points (Pagination)
	while (pageEnd != candidateGuides.end())
	{
		if (RNAFoldPageLength > 0)
		{
			// Advance the pageEnd pointer
			std::advance(pageEnd, std::min( (int)std::distance(pageEnd, candidateGuides.end()), RNAFoldPageLength));
			// Record page start
			pageStart = paginatorIterator;
			// Print page information
			snprintf(printingBuffer, 1024, "\tProcessing page %d (%d per page).", pgIdx, RNAFoldPageLength);
			printer(printingBuffer);
		}
		else {
			// Process all guides at once
			pageEnd = candidateGuides.end();
		}
		printer("\t\tConstructing the RNAfold input file.");

		// Open input file 
		std::ofstream out;
		out.open(RNAFoldInFile);

		guidesInPage = 0;
		for (paginatorIterator; paginatorIterator != pageEnd; paginatorIterator++)
		{
			string target23 = paginatorIterator->first;
			map<string, string> resultsMap = paginatorIterator->second;
			// Run time filtering
			if (!filterCandidateGuides(resultsMap, MODULE_MM10DB, optimsationLevel, consensusN, toolCount)) { continue; }
			
			out << "G" << target23.substr(1, 19) << guide << "\n";
			guidesInPage++;
		}

		out.close();

		// Call RNAFold
		snprintf(printingBuffer, 1024, "%s --noPS -j%d -i %s > %s", RNAFoldBin.c_str(), 1, RNAFoldInFile.c_str(), RNAFoldOutFile.c_str());
		runner(printingBuffer);

		printer("\t\tStarting to process the RNAfold results.");
		// Open output file
		std::ifstream in;
		in.open(RNAFoldOutFile);

		map<string, list<string>> RNAstructures;
		int i = 0;
		string L1, L2, target;
		for (string line; std::getline(in, line); i++)
		{
			if (i % 2 == 0)
			{
				// 0th, 2nd, 4th, etc.
				L1 = rtrim(line);
				target = line.substr(0, 20);
			}
			else
			{
				// 1st, 3rd, 5th, etc.
				L2 = rtrim(line);
				RNAstructures[transToDNA(target.substr(1, 19))] = { L1, L2, target };
			}
		}

		// Reset paginatorIterator to page start
		paginatorIterator = pageStart;

		for (paginatorIterator; paginatorIterator != pageEnd; paginatorIterator++)
		{

			string target23 = paginatorIterator->first;
			string key = target23.substr(1, 19);
			map<string, string> resultsMap = paginatorIterator->second;

			// Run time filtering
			if (!filterCandidateGuides(resultsMap, MODULE_MM10DB, optimsationLevel, consensusN, toolCount)) { continue; }

			if (RNAstructures.find(key) == RNAstructures.end())
			{
				snprintf(printingBuffer, 1024, "Could not find: %s", key.c_str());
				printer(printingBuffer);
				notFoundCount++;
				continue;
			}
			list<string>::iterator listIt = RNAstructures[key].begin();
			L1 = *listIt;
			listIt++;
			L2 = *listIt;
			listIt++;
			target = *listIt;
			string structure = L2.substr(0, L2.find(" "));
			string energy = L2.substr(L2.find(" ")+2);
			candidateGuides[target23]["ssL1"] = L1;
			candidateGuides[target23]["ssStructure"] = structure;
			candidateGuides[target23]["ssEnergy"] = energy;

			if ( (transToDNA(target) != target23.substr(0, 20)) && 
				 ((transToDNA("C" + target.substr(1))) != target23.substr(0, 20)) && 
				 ((transToDNA("A" + target.substr(1))) != target23.substr(0, 20)) )
			{
				candidateGuides[target23]["passedSecondaryStructure"] = CODE_ERROR;
				errorCount++;
				continue;
			}

			std::smatch match_structure;
			std::smatch match_energy;
			if (std::regex_search(L2, match_structure, pattern_RNAstructure))
			{ 
				float energy = std::stof(match_structure[1].str());
				if (energy < lowEnergyThreshold)
				{
					candidateGuides[target23]["passedSecondaryStructure"] = CODE_REJECTED;
					failedCount++;
				}
				else
				{
					candidateGuides[target23]["passedSecondaryStructure"] = CODE_ACCEPTED;
				}
			}
			else if (std::regex_search(L2, match_energy, pattern_RNAenergy))
			{
				float energy = std::stof(match_energy[1].str());
				if (energy <= highEnergyThreshold)
				{
					candidateGuides[target23]["passedSecondaryStructure"] = CODE_REJECTED;
					failedCount++;
				}
				else
				{
					candidateGuides[target23]["passedSecondaryStructure"] = CODE_ACCEPTED;
				}
			}
			testedCount++;
		}

		// Advance paginatorIterator to page end for next loop
		paginatorIterator = pageEnd;
		pgIdx++;
	}
	snprintf(printingBuffer, 1024, "%d of %d failed here.", failedCount, testedCount);
	printer(printingBuffer);
	if (errorCount > 0) 
	{
		snprintf(printingBuffer, 1024, "%d of %d errored here.", errorCount, testedCount);
		printer(printingBuffer);
	}
	if (notFoundCount > 0)
	{
		snprintf(printingBuffer, 1024, "%d of %d not found in RNAfold output.", notFoundCount, testedCount);
		printer(printingBuffer);
	}

	printer("Calculating mm10db final result.");
	int acceptedCount = 0;
	failedCount = 0;
	for (auto const& [target23, resultsMap] : candidateGuides)
	{
		if ((candidateGuides[target23]["passedAvoidLeadingT"] == CODE_ACCEPTED) &&
			(candidateGuides[target23]["passedATPercent"] == CODE_ACCEPTED) &&
			(candidateGuides[target23]["passedTTTT"] == CODE_ACCEPTED) &&
			(candidateGuides[target23]["passedSecondaryStructure"] == CODE_ACCEPTED))
		{
			candidateGuides[target23]["acceptedByMm10db"] = CODE_ACCEPTED;
			acceptedCount++;
		}
		else
		{
			candidateGuides[target23]["acceptedByMm10db"] = CODE_REJECTED;
			failedCount++;
		}
	}

	snprintf(printingBuffer, 1024, "\t%d accepted.", acceptedCount);
	printer(printingBuffer);
	snprintf(printingBuffer, 1024, "\t%d rejected.", failedCount);
	printer(printingBuffer);

}

bool mm10db::leadingT(std::string candidateGuide)
{
	return	( (candidateGuide[0] == 'T' || candidateGuide[0] == 't') && (candidateGuide.substr(candidateGuide.length() - 2, 2) == "GG" || candidateGuide.substr(candidateGuide.length() - 2, 2) == "gg") ) ||
			( (candidateGuide[candidateGuide.length() -1] == 'A' || candidateGuide[candidateGuide.length() - 1] == 'a') && (candidateGuide.substr(0, 2) == "CC" || candidateGuide.substr(0, 2) == "cc"));
}

float mm10db::AT_percentage(std::string candidateGuide)
{
	float total = 0.0f;
	float length = candidateGuide.size();
	array<char, 4>::iterator p;
	for (int i = 0; i < candidateGuide.size(); i++)
	{
		// Check if the char at the current pos is present in the 'AT' array
		p = std::find(atArray.begin(), atArray.end(), candidateGuide[i]);
		if (p != atArray.end())
		{
			total++;
		}
	}

	return (100.0f * total / length);
}

bool mm10db::polyT(std::string candidateGuide)
{
	for (int i = 0; i < candidateGuide.size() - 4; i++)
	{
		if (candidateGuide.substr(i, 4) == "TTTT" || candidateGuide.substr(i, 4) == "tttt")
		{
			return true;
		};
	}
	return false;
}

std::string mm10db::transToDNA(std::string RNA)
{
	// Swap U with T
	std::replace(RNA.begin(), RNA.end(), 'u', 't');
	std::replace(RNA.begin(), RNA.end(), 'U', 'T');
	return RNA;
}