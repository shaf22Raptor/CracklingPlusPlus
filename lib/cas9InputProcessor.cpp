// cas9InputProcessor class, implementation of inputProcessor for Cas9.
#include "../include/cas9InputProcessor.hpp"

using std::string;
using std::list;
using std::unordered_set;
using std::regex;
using std::regex_iterator;
using std::string_view;
using std::match_results;
using std::ifstream;
using std::ofstream;
using std::filesystem::path;
using std::runtime_error;


// Guide searching
const regex patternForward("(?=([ATCG]{21}GG))", std::regex::optimize);
const regex patternReverse("(?=(CC[ACGT]{21}))", std::regex::optimize);

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
	std::vector<string> seq;

	// Duplicate tracking
	unordered_set<string> candidateGuides;
	unordered_set<string> recordedSequences;
	numDuplicateGuides = 0;
	numIdentifiedGuides = 0;

	// Setup temp working dir
	path systemTempDir = std::filesystem::temp_directory_path();
	if (std::filesystem::is_directory(systemTempDir / "Crackling")) { std::filesystem::remove_all(systemTempDir / "Crackling"); }
	if (!std::filesystem::create_directory(systemTempDir / "Crackling")) { throw FileSystemError(); }
	path tempWorkingDir(systemTempDir / "Crackling");

	
	printer(fmt::format("Storing batch files in: {}", tempWorkingDir.string()));


	// Create first batch file
	outFileName = tempWorkingDir / (std::to_string(batchFiles.size()) + "_batchFile.txt");
	batchFiles.push_back(outFileName.string());
	outFile.open(outFileName, std::ios::binary);

	int seqLength = 0;

	// Begin processing files
	for (const string& file : filesToProcess)
	{
		// Identify file format
		printer(fmt::format("Identifying possible target sites in : {}", file));


		inFile.open(file, std::ios::binary);
		std::getline(inFile, inputLine);

		// file is Fasta formatted
		if (inputLine[0] == '>')
		{
			seqHeader = trim(inputLine).substr(1);
			seq = {};
			while (std::getline(inFile, inputLine)) 
			{
				inputLine = trim(inputLine);
				if (inputLine[0] == '>')
				{
					if (recordedSequences.find(seqHeader) == recordedSequences.end())
					{
						recordedSequences.insert(seqHeader);
						string concatanatedSeq = makeUpper(std::accumulate(seq.begin(), seq.end(), std::string{}));

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
					seqLength += inputLine.length();
				}
			}
			// EOF process last seq
			if (recordedSequences.find(seqHeader) == recordedSequences.end())
			{
				recordedSequences.insert(seqHeader);
				string concatanatedSeq = makeUpper(std::accumulate(seq.begin(), seq.end(), std::string{}));

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

		printer(fmt::format("\tProcessed {:.2f}% of input.", completedPercent));

	}

	// Finished processing files, report results
	float duplicatePercent = ((float)numDuplicateGuides / (float)numIdentifiedGuides) * 100.0f;

	printer(fmt::format(
		"\tIdentified {} possible target sites.\n"
		"\tOf these, {} are not unique. These sites occur a total of {} times.\n"
		"\t{} of {} ({:.2f}%) of guides will be ignored for optimisation levels over ultralow.\n"
		"\t{} distinct guides were identified.",
		commaify(numIdentifiedGuides),
		commaify((int)duplicateGuides.size()),
		commaify(numDuplicateGuides),
		commaify(numDuplicateGuides),
		commaify(numIdentifiedGuides),
		duplicatePercent,
		commaify((int)candidateGuides.size()))
	);

	return;
}

void cas9InputProcessor::processSeqeunce(
	string_view seqeunce,
	string_view seqHeader,
	ofstream& outFile,
	const path& tempWorkingDir,
	unordered_set<string>& candidateGuides,
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
		if (candidateGuides.find(guide) == candidateGuides.end())
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

	for (auto regexItr = regex_iterator<string_view::const_iterator>(seqeunce.begin(), seqeunce.end(), patternReverse);
		regexItr != regex_iterator<string_view::const_iterator>();
		regexItr++)
	{
		numIdentifiedGuides++;
		match_results<string_view::const_iterator> m = *regexItr;
		string guide = rc(m[1].str());
		long long matchPos = m.position();
		if (candidateGuides.find(guide) == candidateGuides.end())
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
			outFile << guide << "," << seqHeader << "," << matchPos << "," << (matchPos + 23) << "," << "-" << "\n";
		}
		else
		{
			numDuplicateGuides++;
			duplicateGuides.insert(guide);
		}
	}
	return;
}

const list<string>& cas9InputProcessor::getBatchFiles() const
{
	return batchFiles;
}

bool cas9InputProcessor::isDuplicateGuide(const string& guide) const
{
	return duplicateGuides.find(guide) != duplicateGuides.end();
}

void cas9InputProcessor::cleanUp()
{
	std::filesystem::remove_all(path(std::filesystem::temp_directory_path() / "Crackling"));
	return;
}