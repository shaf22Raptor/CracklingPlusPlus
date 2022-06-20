// bowtie2 Class
#include <bowtie2.hpp>

using std::string;
using std::map;
using std::vector;


bowtie2::bowtie2(ConfigManager cm) :
	toolIsSelected(false),
	optimsationLevel(""),
	toolCount(0),
	consensusN(0),
	bowtie2OutFile(""),
	bowtie2InFile(""),
	bowtie2Bin(""),
	bowtie2PageLength(0),
	bowtie2Index(""),
	bowtie2Threads(0)
{
	toolIsSelected = cm.getBool("offtargetscore", "enabled");
	optimsationLevel = cm.getString("general", "optimisation");
	toolCount = cm.getConsensusToolCount();
	consensusN = cm.getInt("consensus", "n");
	bowtie2OutFile = cm.getString("bowtie2", "output");
	bowtie2InFile = cm.getString("bowtie2", "input");
	bowtie2Bin = cm.getString("bowtie2", "binary");
	bowtie2PageLength = cm.getInt("bowtie2", "page-length");
	bowtie2Index = cm.getString("input", "bowtie2-index");
	bowtie2Threads = cm.getInt("bowtie2", "threads");
}

void bowtie2::run(map<string, map<string, string>>& candidateGuides)
{
	if (!toolIsSelected)
	{
		printer("bowtie2 has been configured not to run. Skipping bowtie2");
		return;
	}	

	printer("Bowtie analysis.");
	char printingBuffer[1024];
	int failedCount = 0;
	int testedCount = 0;
	int pgIdx = 1;
	int guidesInPage = 0;
	map<string, map<string, string>>::iterator paginatorIterator = candidateGuides.begin();
	map<string, map<string, string>>::iterator pageStart = candidateGuides.begin();
	map<string, map<string, string>>::iterator pageEnd = candidateGuides.begin();


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
			snprintf(printingBuffer, 1024, "\tProcessing page %d (%d per page).", pgIdx, bowtie2PageLength);
			printer(printingBuffer);
		}
		else {
			// Process all guides at once
			pageEnd = candidateGuides.end();
		}
		printer("\t\tConstructing the Bowtie input file.");
		map<string, string> tempTargetDict_offset;
		// Open input file 
		std::ofstream inFile;
		inFile.open(bowtie2InFile);

		guidesInPage = 0;
		for (paginatorIterator; paginatorIterator != pageEnd; paginatorIterator++)
		{
			string target23 = paginatorIterator->first;
			map<string, string> resultsMap = paginatorIterator->second;
			// Run time filtering
			if (!filterCandidateGuides(resultsMap, MODULE_SPECIFICITY, optimsationLevel, consensusN, toolCount)) { 
				// Advance page end
				if (pageEnd != candidateGuides.end())
				{
					pageEnd++;
				}
				continue; 
			}

			vector<string> similarTargets ({
				target23.substr(0, 20) + "AGG",
				target23.substr(0, 20) + "CGG",
				target23.substr(0, 20) + "GGG",
				target23.substr(0, 20) + "TGG",
				target23.substr(0, 20) + "AAG",
				target23.substr(0, 20) + "CAG",
				target23.substr(0, 20) + "GAG",
				target23.substr(0, 20) + "TAG"
			});

			for (string bowtieTarget : similarTargets)
			{
				inFile << bowtieTarget << "\n";
				tempTargetDict_offset[bowtieTarget] = target23;
			}

			guidesInPage++;
		}

		inFile.close();

		snprintf(printingBuffer, 1024, "\t\t%d guides in this page.", guidesInPage);
		printer(printingBuffer);

		// Call bowtie2
		snprintf(printingBuffer, 1024, "%s -x %s -p %d --reorder --no-hd -t -r -U %s -S %s", bowtie2Bin.c_str(), bowtie2Index.c_str(), bowtie2Threads, bowtie2InFile.c_str(), bowtie2OutFile.c_str());
		runner(printingBuffer);

		printer("\tStarting to process the Bowtie results.");
		
		// Open output file 
		std::ifstream outFile;
		outFile.open(bowtie2OutFile);

		vector<string> bowtie2Results;

		for (string line; std::getline(outFile, line); )
		{
			bowtie2Results.push_back(line);
		}

		for (int i = 0; i < bowtie2Results.size(); i += 8)
		{
			int nb_occurences = 0;
			bowtie2Results[i] = rtrim(bowtie2Results[i]);
			
			vector<string> line;

			size_t tabPos = 0;
			while ((tabPos = bowtie2Results[i].find('\t')) != std::string::npos) {
				line.push_back(bowtie2Results[i].substr(0, tabPos));
				bowtie2Results[i].erase(0, tabPos + 1);
			}
			line.push_back(bowtie2Results[i].substr(0, tabPos));
			string chr = line[2];
			int pos = stoi(line[3]);
			string read = line[9];
			string seq = "";
			if (tempTargetDict_offset.find(read) != tempTargetDict_offset.end())
			{
				seq = tempTargetDict_offset[read];
			}
			else if (tempTargetDict_offset.find(rc(read)) != tempTargetDict_offset.end())
			{
				seq = tempTargetDict_offset[rc(read)];
			}
			else
			{
				std::cout << "Problem? " << read << std::endl;
			}
			if (seq.substr(seq.length()-2) == "GG")
			{
				candidateGuides[seq]["bowtieChr"] = chr;
				candidateGuides[seq]["bowtieStart"] = pos;
				candidateGuides[seq]["bowtieEnd"] = pos + 22;
			}
			else if (rc(seq).substr(seq.length() - 2) == "CC")
			{
				candidateGuides[seq]["bowtieChr"] = chr;
				candidateGuides[seq]["bowtieStart"] = pos;
				candidateGuides[seq]["bowtieEnd"] = pos + 22;
			}
			else
			{
				std::cout << "Error? " << seq << std::endl;
				exit(-1);
			}

			for (int j = i; j < i + 8; j++)
			{
				if (bowtie2Results[j].find("XM:i:0") != string::npos)
				{
					nb_occurences++;
					if (bowtie2Results[j].find("XM:i:0") != string::npos) 
					{
						nb_occurences++;
					}
				}
			}

			//TODO: Check if this logic makes sense. Later result could undo recorded failure
			if (nb_occurences > 1)
			{
				if (candidateGuides[seq]["passedBowtie"] != CODE_REJECTED)
				{
					failedCount++;
				}
				candidateGuides[seq]["passedBowtie"] = CODE_REJECTED;
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
	snprintf(printingBuffer, 1024, "\t%d of %d failed here.", failedCount, testedCount);
	printer(printingBuffer);

}
