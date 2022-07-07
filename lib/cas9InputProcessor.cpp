// cas9InputProcessor class, implementation of inputProcessor for Cas9.
#include <cas9InputProcessor.hpp>

using std::string;
using std::list;
using std::set;
using std::regex;
using std::regex_iterator;
using std::string_view;
using std::match_results;
using std::ifstream;
using std::ofstream;
using std::filesystem::path;
using std::runtime_error;


// Guide searching
const regex patternForward("(?=([ATCG]{21}GG))");
const regex patternReverse("(?=(CC[ACGT]{21}))");

void cas9InputProcessor::process(list<string> const & filesToProcess, int const & batchSize)
{
	printer("Analysing files...");

	// Progress reporting
	int totalSizeBytes = 0;
	int completedSizeBytes = 0;
	double completedPercent = 0.0;
	for (const string& file : filesToProcess)
	{
		totalSizeBytes += std::filesystem::file_size(path(file));
	}

	// Formatting prints
	string printingBuffer;

	// In file
	ifstream inFile;
	string inputLine;

	// Batched files
	ofstream outFile;
	path outFileName;
	guidesInBatch = 0;

	// File processing
	string seqHeader;
	list<string> seq;

	// Duplicate tracking
	set<string, std::less<>> candidateGuides;
	set<string, std::less<>> recordedSequences;
	numDuplicateGuides = 0;
	numIdentifiedGuides = 0;

	// Setup temp working dir
	path systemTempDir = std::filesystem::temp_directory_path();
	if (std::filesystem::is_directory(systemTempDir / "Crackling")) { std::filesystem::remove_all(systemTempDir / "Crackling"); }
	if (!std::filesystem::create_directory(systemTempDir / "Crackling")) { throw runtime_error("Unable to create temp working dir!"); }
	path tempWorkingDir(systemTempDir / "Crackling");

	
	printer(std::format("Storing batch files in: {0}", tempWorkingDir.string()));


	// Create first batch file
	outFileName = tempWorkingDir / (std::to_string(batchFiles.size()) + "_batchFile.txt");
	batchFiles.push_back(outFileName.string());
	outFile.open(outFileName, std::ios::binary);

	// Begin processing files
	for (const string& file : filesToProcess)
	{
		// Identify file format
		printer(std::format("Identifying possible target sites in : {0}", file));


		inFile.open(file, std::ios::binary);
		std::getline(inFile, inputLine);

		// file is Fasta formatted
		if (inputLine[0] == '>')
		{
			seqHeader = inputLine.substr(1);
			seq = {};
			while (std::getline(inFile, inputLine)) 
			{
				inputLine = trim(inputLine);
				if (inputLine[0] == '>')
				{
					if (!recordedSequences.contains(seqHeader))
					{
						recordedSequences.insert(seqHeader);
						string concatanatedSeq;
						for (string seqFragment : seq) { concatanatedSeq += makeUpper(seqFragment); }

						processSeqeunce(
							concatanatedSeq,
							seqHeader,
							outFile,
							tempWorkingDir,
							candidateGuides,
							batchSize
						);
					}
					seqHeader = inputLine.substr(1);
					seq = {};
				}
				else
				{
					seq.push_back(inputLine);
				}
			}
			// EOF process last seq
			if (!recordedSequences.contains(seqHeader))
			{
				recordedSequences.insert(seqHeader);
				string concatanatedSeq;
				for (string seqFragment : seq) { concatanatedSeq += makeUpper(seqFragment); }

				processSeqeunce(
					concatanatedSeq,
					seqHeader,
					outFile,
					tempWorkingDir,
					candidateGuides,
					batchSize
				);
			}
			outFile.close();
		}
		// file is plain text, assume one sequence per line
		else
		{
			do
			{
				// Clean up input line
				inputLine = makeUpper(trim(inputLine));

				processSeqeunce(
					inputLine,
					file,
					outFile,
					tempWorkingDir,
					candidateGuides,
					batchSize
				);

			} while (std::getline(inFile, inputLine));

		}

		// Report overall progress before processing next file
		completedSizeBytes += std::filesystem::file_size(path(file));
		completedPercent = completedSizeBytes / totalSizeBytes * 100.0;

		std::cout.precision(2);
		printer(std::format("\tProcessed {}% of input.", completedPercent));

	}

	// Finished processing files, report results
	float duplicatePercent = ((float)numDuplicateGuides / (float)numIdentifiedGuides) * 100.0f;

	printer(std::format(
		"\tIdentified {} possible target sites.\n"
		"\tOf these, {} are not unique. These sites occur a total of {} times.\n"
		"\t{} of {} ({}%) of guides will be ignored for optimisation levels over ultralow.\n"
		"\t{} distinct guides were identified.",
		numIdentifiedGuides, (int)duplicateGuides.size(), numDuplicateGuides, numDuplicateGuides, numIdentifiedGuides, duplicatePercent, (int)candidateGuides.size())
	);

	return;
}

void cas9InputProcessor::processSeqeunce(
	string_view seqeunce, 
	string_view seqHeader,
	ofstream& outFile,
	path const& tempWorkingDir,
	set<string, std::less<>>& candidateGuides,
	const int& batchSize
	)
{

	for (auto regexItr = regex_iterator<string_view::const_iterator>(seqeunce.begin(), seqeunce.end(), patternForward);
		regexItr != regex_iterator<string_view::const_iterator>();
		regexItr++)
	{
		numIdentifiedGuides++;
		match_results<string_view::const_iterator> m = *regexItr;
		string guide = m[1].str();
		long long matchPos = m.position();
		if (!candidateGuides.contains(guide))
		{
			candidateGuides.insert(guide);
			if (++guidesInBatch > batchSize)
			{
				outFile.close();
				path outFileName = tempWorkingDir / (std::to_string(batchFiles.size()) + "_batchFile.txt");
				batchFiles.push_back(outFileName.string());
				outFile.open(outFileName, std::ios::binary);
				guidesInBatch = 1;
			}
			outFile << guide << "," << seqHeader << "," << matchPos << "," << (matchPos + 23) << "," << "+" << "\n";
		}
		else
		{
			numIdentifiedGuides++;
			duplicateGuides.insert(guide);
		}
	}

	for (auto regexItr = regex_iterator<string_view::const_iterator>(seqeunce.begin(), seqeunce.end(), patternReverse);
		regexItr != regex_iterator<string_view::const_iterator>();
		regexItr++)
	{
		numIdentifiedGuides++;
		match_results<string_view::const_iterator> m = *regexItr;
		string guide = rc(m[1].str());
		long long matchPos = m.position();
		if (!candidateGuides.contains(guide))
		{
			candidateGuides.insert(guide);
			if (++guidesInBatch > batchSize)
			{
				outFile.close();
				path outFileName = tempWorkingDir / (std::to_string(batchFiles.size()) + "_batchFile.txt");
				batchFiles.push_back(outFileName.string());
				outFile.open(outFileName, std::ios::binary);
				guidesInBatch = 1;
			}
			outFile << guide << "," << seqHeader << "," << matchPos << "," << (matchPos + 23) << "," << "+" << "\n";
		}
		else
		{
			numDuplicateGuides++;
			duplicateGuides.insert(guide);
		}
	}
	return;
}

list<string> cas9InputProcessor::getBatchFiles()
{
	return batchFiles;
}

bool cas9InputProcessor::isDuplicateGuide(string guide)
{
	return duplicateGuides.contains(guide);
}

void cas9InputProcessor::cleanUp()
{
	std::filesystem::remove_all(path(std::filesystem::temp_directory_path() / "Crackling"));
	return;
}