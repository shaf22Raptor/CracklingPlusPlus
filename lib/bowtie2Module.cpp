#include "../include/bowtie2Module.hpp"

using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::ofstream;
using std::unordered_map;
using std::filesystem::remove;
using boost::algorithm::split;

bowtie2Module::bowtie2Module(const cracklingConfig& config) : specificityModule(config)
{
	this->config = config.bowtie2;
	this->indexFile = config.input.bowtie2Index;
}

void bowtie2Module::run(std::vector<guideResults>& candidateGuides)
{

	if (!toolIsSelected)
	{
		cout << "Bowtie2 has been configured not to run. Skipping bowtie2" << endl;
		return;
	}

	cout << "Bowtie analysis." << endl;
	uint64_t failedCount = 0;
	uint64_t testedCount = 0;
	uint64_t pgIdx = 1;
	uint64_t guidesInPage = 0;
	auto paginatorIterator = candidateGuides.begin();
	auto pageStart = candidateGuides.begin();
	auto pageEnd = candidateGuides.begin();

	// Outer loop deals with changing iterator start and end points (Pagination)
	while (pageEnd != candidateGuides.end())
	{
		if (config.pageLen > 0)
		{
			// Advance the pageEnd pointer
			std::advance(pageEnd, std::min(static_cast<uint64_t>(std::distance(pageEnd, candidateGuides.end())), config.pageLen));
			// Record page start
			pageStart = paginatorIterator;
			// Print page information
			cout << fmt::format(comma_locale, "\tProcessing page {:L} ({:L} per page).", pgIdx, config.pageLen) << endl;
		}
		else {
			// Process all guides at once
			pageEnd = candidateGuides.end();
		}

		// Open input file 
		cout << "\t\tConstructing the Bowtie input file." << endl;
		ofstream inFile;
		inFile.open(config.inFile, std::ios::binary | std::ios::out);
		guidesInPage = 0;
		while (paginatorIterator != pageEnd)
		{
			// Run time filtering
			if (!processGuide(*paginatorIterator)) {
				// Advance page end
				if (pageEnd != candidateGuides.end())
				{
					pageEnd++;
				}
				paginatorIterator++;
				continue;
			}

			vector<string> similarTargets = {
				paginatorIterator->seq.substr(0, 20) + "AGG",
				paginatorIterator->seq.substr(0, 20) + "CGG",
				paginatorIterator->seq.substr(0, 20) + "GGG",
				paginatorIterator->seq.substr(0, 20) + "TGG",
				paginatorIterator->seq.substr(0, 20) + "AAG",
				paginatorIterator->seq.substr(0, 20) + "CAG",
				paginatorIterator->seq.substr(0, 20) + "GAG",
				paginatorIterator->seq.substr(0, 20) + "TAG"
			};

			for (string bowtieTarget : similarTargets) { inFile << bowtieTarget << "\n"; }

			guidesInPage++;
			paginatorIterator++;
		}
		inFile.close();

		cout << fmt::format(comma_locale, "\t\t{:L} guides in this page.", guidesInPage) << endl;

		// Call bowtie2
		runner(fmt::format("{} -x {} -p {} --reorder --no-hd -t -r -U {} -S {}", config.binary.string(), indexFile.string(), config.threads, config.inFile.string(), config.outFile.string()).c_str());

		cout << "\tStarting to process the Bowtie results." << endl;

		// Reset page pointer
		paginatorIterator = pageStart;

		std::ifstream outFile;
		outFile.open(config.outFile, std::ios::binary | std::ios::in);
		while (paginatorIterator != pageEnd)
		{
			// Run time filtering
			if (!processGuide(*paginatorIterator)) {
				paginatorIterator++;
				continue;
			}

			uint64_t nb_occurences = 0;

			vector<string> bowtie2Results;
			for (uint8_t j = 0; j < 8; j++)
			{
				string line;
				std::getline(outFile, line);
				bowtie2Results.push_back(line);
			}

			vector<string> bowtie2Output;
			split(bowtie2Output, bowtie2Results[0], boost::is_any_of("\t"));
			string chr = bowtie2Output[2];
			int pos = stoi(bowtie2Output[3]);
			string read = bowtie2Output[9];
			string seq(23, ' ');

			paginatorIterator->bowtie2Chr = bowtie2Output[2];
			paginatorIterator->bowtie2Start = stoull(bowtie2Output[3]);
			paginatorIterator->bowtie2End = paginatorIterator->bowtie2Start + 22;

			// We count how many of the eight reads for this target have a perfect alignment
			for (const string& outputLine : bowtie2Results)
			{
				/**
				* http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#sam-output
				* XM : i : <N>    The number of mismatches in the alignment.Only present if SAM record is for an aligned read.
				* XS : i : <N>    Alignment score for the best - scoring alignment found other than the alignment reported.
				*/
				if (outputLine.find("XM:i:0") != string::npos)
				{
					nb_occurences++;
					// We also check whether this perfect alignment also happens elsewhere
					if (outputLine.find("XS:i:0") != string::npos)
					{
						nb_occurences++;
					}
				}
			}
			// If that number is at least two, the target is removed
			if (nb_occurences > 1)
			{
				paginatorIterator->passedBowtie2 = CODE_REJECTED;
				failedCount++;
			}
			else
			{
				paginatorIterator->passedBowtie2 = CODE_ACCEPTED;
			}
			testedCount++;
			paginatorIterator++;
		}
		outFile.close();
		// Clean up intermediate files
		remove(config.inFile);
		remove(config.outFile);
		// Point paginatorIterator to page end for next loop
		paginatorIterator = pageEnd;
		pgIdx++;
	}
	cout << fmt::format(comma_locale, "\t{:L} of {:L} failed here.", failedCount, testedCount) << endl;
}

bool bowtie2Module::processGuide(const guideResults& guide)
{
	// Process all guides at this level
	if (optimsationLevel == optimisationLevel::ultralow) { return true; }

	// For all levels above `ultralow`
	if (!guide.isUnique)
	{
		// Reject all guides that have been seen more than once
		return false;
	}

	// For optimisation levels `medium` and `high`
	if (optimsationLevel == optimisationLevel::medium || optimsationLevel == optimisationLevel::high)
	{
		// Reject if the consensus was not passed
		if (guide.consensusCount < consensusN) { return false; }
	}

	// None of the failure conditions have been meet, return true
	return true;
}
