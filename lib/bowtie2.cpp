// bowtie2 Class
#include "../include/bowtie2.hpp"

using std::string;
using std::unordered_map;
using std::vector;


bowtie2::bowtie2(ConfigManager& cm) :
	toolIsSelected(cm.getBool("offtargetscore", "enabled")),
	optimsationLevel(cm.getString("general", "optimisation")),
	toolCount(cm.getConsensusToolCount()),
	consensusN(cm.getInt("consensus", "n")),
	threadCount(cm.getInt("bowtie2", "threads")),
	bowtie2OutFile(cm.getString("bowtie2", "output")),
	bowtie2InFile(cm.getString("bowtie2", "input")),
	bowtie2Bin(cm.getString("bowtie2", "binary")),
	bowtie2Index(cm.getString("input", "bowtie2-index")),
	bowtie2PageLength(cm.getInt("bowtie2", "page-length"))
{}

void bowtie2::run(unordered_map<string, unordered_map<string, string>>& candidateGuides)
{
	if (!toolIsSelected)
	{
		printer("bowtie2 has been configured not to run. Skipping bowtie2");
		return;
	}

	printer("Bowtie analysis.");
	int failedCount = 0;
	int testedCount = 0;
	int pgIdx = 1;
	int guidesInPage = 0;
	auto paginatorIterator = candidateGuides.begin();
	auto pageStart = candidateGuides.begin();
	auto pageEnd = candidateGuides.begin();


	// Outer loop deals with changing iterator start and end points (Pagination)
	while (pageEnd != candidateGuides.end())
	{
		if (bowtie2PageLength > 0)
		{
			// Advance the pageEnd pointer
			std::advance(pageEnd, std::min((int)std::distance(pageEnd, candidateGuides.end()), bowtie2PageLength));
			// Record page start
			pageStart = paginatorIterator;
			// Print page information
			printer(fmt::format("\tProcessing page {} ({} per page).", commaify(pgIdx), commaify(bowtie2PageLength)));
		}
		else {
			// Process all guides at once
			pageEnd = candidateGuides.end();
		}
		printer("\t\tConstructing the Bowtie input file.");
		unordered_map<string, string> tempTargetDict_offset;
		// Open input file 
		std::ofstream inFile;
		inFile.open(bowtie2InFile, std::ios::binary);
		guidesInPage = 0;
		vector<char> queryDataSet;
		queryDataSet.reserve(bowtie2PageLength * 23);
		while (paginatorIterator != pageEnd)
		{
			string target23 = paginatorIterator->first;
			// Run time filtering
			if (!filterCandidateGuides(paginatorIterator->second, MODULE_SPECIFICITY, optimsationLevel, consensusN, toolCount)) {
				// Advance page end
				if (pageEnd != candidateGuides.end())
				{
					pageEnd++;
				}
				paginatorIterator++;
				continue;
			}
			// Store original seqeunce for use later when parsing results.
			for (char c : target23) { queryDataSet.push_back(c); }
			
			vector<string> similarTargets = {
				target23.substr(0, 20) + "AGG",
				target23.substr(0, 20) + "CGG",
				target23.substr(0, 20) + "GGG",
				target23.substr(0, 20) + "TGG",
				target23.substr(0, 20) + "AAG",
				target23.substr(0, 20) + "CAG",
				target23.substr(0, 20) + "GAG",
				target23.substr(0, 20) + "TAG"
			};

			for (string bowtieTarget : similarTargets) { inFile << bowtieTarget << "\n"; }

			guidesInPage++;
			paginatorIterator++;
		}

		inFile.close();

		printer(fmt::format("\t\t{} guides in this page.", commaify(guidesInPage)));

		// Call bowtie2
		runner(fmt::format("{} -x {} -p {} --reorder --no-hd -t -r -U {} -S {}", bowtie2Bin, bowtie2Index, threadCount, bowtie2InFile, bowtie2OutFile).c_str());

		printer("\tStarting to process the Bowtie results.");

		// Open output file 
		std::ifstream outFile;
		outFile.open(bowtie2OutFile, std::ios::binary);

		vector<string> bowtie2Results;

		for (string line; std::getline(outFile, line); )
		{
			bowtie2Results.push_back(line);
		}

		for (int i = 0; i < bowtie2Results.size(); i += 8)
		{
			int nb_occurences = 0;

			vector<string> bowtie2Output = split(bowtie2Results[i], "\t");
			string chr = bowtie2Output[2];
			int pos = stoi(bowtie2Output[3]);
			string read = bowtie2Output[9];
			string seq(23, ' ');
			// Bowtie2 results are written in the same order as input, retrieve target seq from `queryDataSet`
			for (int j = 0; j < 23; j++) {
				seq[j] = queryDataSet[((i/8) * 23) + j];
			}
			candidateGuides[seq]["bowtieChr"] = chr;
			candidateGuides[seq]["bowtieStart"] = std::to_string(pos);
			candidateGuides[seq]["bowtieEnd"] = std::to_string(pos + 22);

			// We count how many of the eight reads for this target have a perfect alignment
			for (int j = i; j < i + 8; j++)
			{
				/**
				* http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#sam-output
				* XM : i : <N>    The number of mismatches in the alignment.Only present if SAM record is for an aligned read.
				* XS : i : <N>    Alignment score for the best - scoring alignment found other than the alignment reported.
				*/
				if (bowtie2Results[j].find("XM:i:0") != string::npos)
				{
					nb_occurences++;
					// We also check whether this perfect alignment also happens elsewhere
					if (bowtie2Results[j].find("XS:i:0") != string::npos)
					{
						nb_occurences++;
					}
				}
			}
			// If that number is at least two, the target is removed
			if (nb_occurences > 1)
			{
				candidateGuides[seq]["passedBowtie"] = CODE_REJECTED;
				failedCount++;
			}
			else
			{
				candidateGuides[seq]["passedBowtie"] = CODE_ACCEPTED;
			}
			testedCount++;
		}

		// Point paginatorIterator to page end for next loop
		paginatorIterator = pageEnd;
		pgIdx++;
	}
	printer(fmt::format("\t{} of {} failed here.", commaify(failedCount), commaify(testedCount)));
	return;
}
