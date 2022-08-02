// mm10db class
#include "../include/mm10db.hpp"

using std::string;
using std::unordered_map;
using std::array;
using std::list;

const array<char, 2> atArray = { 'A', 'T' };

mm10db::mm10db(ConfigManager& cm) :
	toolIsSelected(cm.getBool("consensus", "mm10db")),
	optimsationLevel(cm.getString("general", "optimisation")),
	toolCount(cm.getConsensusToolCount()),
	consensusN(cm.getInt("consensus", "n")),
	threadCount(cm.getInt("rnafold", "threads")),
	RNAFoldOutFile(cm.getString("rnafold", "output")),
	RNAFoldInFile(cm.getString("rnafold", "input")),
	RNAFoldBin(cm.getString("rnafold", "binary")),
	lowEnergyThreshold(cm.getFloat("rnafold", "low_energy_threshold")),
	highEnergyThreshold(cm.getFloat("rnafold", "high_energy_threshold")),
	RNAFoldPageLength(cm.getInt("rnafold", "page-length"))
{}

void mm10db::run(unordered_map<string, unordered_map<string, string>>& candidateGuides)
{

	if (!toolIsSelected)
	{
		printer("mm10db has been configured not to run. Skipping mm10db");
		return;
	}

	printer("mm10db - remove all targets with a leading T(+) or trailing A(-).");
	int failedCount = 0;
	int testedCount = 0;
	for (const auto& [target23, resultsMap] : candidateGuides)
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
	printer(fmt::format("\t{} of {} failed here.", commaify(failedCount), commaify(testedCount)));

	printer("mm10db - remove based on AT percent.");
	failedCount = 0;
	testedCount = 0;
	for (const auto& [target23, resultsMap] : candidateGuides)
	{
		// Run time filtering
		if (!filterCandidateGuides(resultsMap, MODULE_MM10DB, optimsationLevel, consensusN, toolCount)) { continue; }

		float AT = AT_percentage(std::string_view(target23).substr(0, 20));
		if (AT < 20.0f || AT > 65.0f)
		{
			candidateGuides[target23]["passedATPercent"] = CODE_REJECTED;
			failedCount++;
		}
		else
		{
			candidateGuides[target23]["passedATPercent"] = CODE_ACCEPTED;
		}
		std::stringstream converstionStream;
		converstionStream << std::fixed << std::setprecision(1) << AT;
		candidateGuides[target23]["AT"] = converstionStream.str();
		testedCount++;
	}
	printer(fmt::format("\t{} of {} failed here.", commaify(failedCount), commaify(testedCount)));

	printer("mm10db - remove all targets that contain TTTT.");
	failedCount = 0;
	testedCount = 0;
	for (const auto& [target23, resultsMap] : candidateGuides)
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
	printer(fmt::format("\t{} of {} failed here.", commaify(failedCount), commaify(testedCount)));

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
	auto paginatorIterator = candidateGuides.begin();
	auto pageStart = candidateGuides.begin();
	auto pageEnd = candidateGuides.begin();


	// Outer loop deals with changing iterator start and end points (Pagination)
	while (pageEnd != candidateGuides.end())
	{
		if (RNAFoldPageLength > 0)
		{
			// Advance the pageEnd pointer
			std::advance(pageEnd, std::min((int)std::distance(pageEnd, candidateGuides.end()), RNAFoldPageLength));
			// Record page start
			pageStart = paginatorIterator;
			// Print page information
			printer(fmt::format("\tProcessing page {} ({} per page).", commaify(pgIdx), commaify(RNAFoldPageLength)));
		}
		else {
			// Process all guides at once
			pageEnd = candidateGuides.end();
		}
		printer("\t\tConstructing the RNAfold input file.");

		// Open input file 
		std::ofstream out;
		out.open(RNAFoldInFile, std::ios::binary);

		guidesInPage = 0;
		while (paginatorIterator != pageEnd)
		{
			string target23 = paginatorIterator->first;
			// Run time filtering
			if (!filterCandidateGuides(paginatorIterator->second, MODULE_MM10DB, optimsationLevel, consensusN, toolCount)) {
				// Advance page end
				if (pageEnd != candidateGuides.end())
				{
					pageEnd++;
				}
				paginatorIterator++;
				continue;
			}
			out << "G" << target23.substr(1, 19) << guide << "\n";
			guidesInPage++;
			paginatorIterator++;
		}

		out.close();

		printer(fmt::format("\t\t{} guides in this page.", commaify(guidesInPage)));

		// Call RNAFold
		runner(fmt::format("{} --noPS -j{} -i {} > {}", RNAFoldBin, threadCount, RNAFoldInFile, RNAFoldOutFile).c_str());

		printer("\t\tStarting to process the RNAfold results.");
		// Open output file
		std::ifstream in;
		in.open(RNAFoldOutFile, std::ios::binary);

		unordered_map<string, list<string>> RNAstructures;
		int i = 0;
		string L1;
		string L2;
		string target;
		for (string line; std::getline(in, line);)
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
			i++;
		}

		// Reset paginatorIterator to page start
		paginatorIterator = pageStart;

		while (paginatorIterator != pageEnd)
		{

			string target23 = paginatorIterator->first;
			string key = target23.substr(1, 19);

			// Run time filtering
			if (!filterCandidateGuides(paginatorIterator->second, MODULE_MM10DB, optimsationLevel, consensusN, toolCount)) {
				paginatorIterator++;
				continue;
			}

			if (RNAstructures.find(key) == RNAstructures.end())
			{
				printer(fmt::format("Could not find: {}", key));
				notFoundCount++;
				continue;
			}

			auto listIt = RNAstructures[key].begin();
			L1 = *listIt;
			listIt++;
			L2 = *listIt;
			listIt++;
			target = *listIt;
			string ssStructure = L2.substr(0, L2.find(" "));
			string ssEnergy = L2.substr(L2.find(" ") + 2, L2.length() - L2.find(" ") - 3);
			candidateGuides[target23]["ssL1"] = L1;
			candidateGuides[target23]["ssStructure"] = ssStructure;
			candidateGuides[target23]["ssEnergy"] = ssEnergy;

			if ((transToDNA(target) != target23.substr(0, 20)) &&
				((transToDNA("C" + target.substr(1))) != target23.substr(0, 20)) &&
				((transToDNA("A" + target.substr(1))) != target23.substr(0, 20)))
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
			paginatorIterator++;
		}

		// Advance paginatorIterator to page end for next loop
		paginatorIterator = pageEnd;
		pgIdx++;
	}
	printer(fmt::format("\t{} of {} failed here.", commaify(failedCount), commaify(testedCount)));
	if (errorCount > 0)
	{
		printer(fmt::format("\t{} of {} errored here.", commaify(errorCount), commaify(testedCount)));
	}
	if (notFoundCount > 0)
	{
		printer(fmt::format("\t{} of {} not found in RNAfold output.", commaify(notFoundCount), commaify(testedCount)));
	}

	printer("Calculating mm10db final result.");
	int acceptedCount = 0;
	failedCount = 0;
	for (const auto& [target23, resultsMap] : candidateGuides)
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
	printer(fmt::format("\t{} accepted.\n\t{} rejected", commaify(acceptedCount), commaify(failedCount)));
}

bool mm10db::leadingT(std::string_view candidateGuide)
{
	return	(startsWith(candidateGuide, "T") && endsWith(candidateGuide, "GG")) ||
		(startsWith(candidateGuide, "CC") && endsWith(candidateGuide, "A"));
}

float mm10db::AT_percentage(std::string_view candidateGuide)
{
	float total = 0.0f;
	float length = (float)candidateGuide.size();

	for (char c : candidateGuide) {
		// Check if the char is present in the 'AT' array
		if (std::find(atArray.begin(), atArray.end(), c) != atArray.end())
		{
			total += 1.0f;
		}
	}

	return (100.0f * total / length);
}

bool mm10db::polyT(std::string_view candidateGuide)
{
	for (int i = 0; i < candidateGuide.size() - 4; i++)
	{
		if (candidateGuide.substr(i, 4) == "TTTT")
		{
			return true;
		}
	}
	return false;
}

string mm10db::transToDNA(string RNA)
{
	// Swap U with T
	std::replace(RNA.begin(), RNA.end(), 'U', 'T');
	return RNA;
}