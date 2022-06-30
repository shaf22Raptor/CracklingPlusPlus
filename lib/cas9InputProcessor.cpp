// cas9InputProcessor class, implementation of inputProcessor for Cas9.
#include <cas9InputProcessor.hpp>

using std::string;
using std::list;
using std::set;
using std::regex;
using std::ifstream;
using std::ofstream;
using std::filesystem::path;

// Guide searching
regex patternForward("(?=([ATCG]{21}GG))");
regex patternReverse("(?=(CC[ACGT]{21}))");

void cas9InputProcessor::processInput(list<string> filesToProcess, int batchSize)
{
	printer("Analysing files...");

	// Progress reporting
	int totalSizeBytes = 0;
	int completedSizeBytes = 0;
	double completedPercent = 0.0;
	for (string file : filesToProcess)
	{
		totalSizeBytes += std::filesystem::file_size(path(file));
	}

	// Formatting prints
	char printingBuffer[1024];

	// In file
	ifstream inFile;
	string inputLine;

	// Batched files
	ofstream outFile;
	path outFileName;
	int guidesInBatch = 0;

	// File processing
	string seqHeader;
	list<string> seq;

	// Duplicate tracking
	set<string> candidateGuides;
	set<string> recordedSequences;
	int numDuplicateGuides = 0;
	int numIdentifiedGuides = 0;

	// Setup temp working dir
	path systemTempDir = std::filesystem::temp_directory_path();
	if (std::filesystem::is_directory(systemTempDir / "Crackling")) { std::filesystem::remove_all(systemTempDir / "Crackling"); }
	if (!std::filesystem::create_directory(systemTempDir / "Crackling")) { throw std::runtime_error("Unable to create temp working dir!"); }
	path tempWorkingDir(systemTempDir / "Crackling");

	snprintf(printingBuffer, 1024, "Storing batch files in: %s", tempWorkingDir.string().c_str());
	printer(printingBuffer);

	// Create first batch file
	outFileName = tempWorkingDir / (std::to_string(batchFiles.size()) + "_batchFile.txt");
	batchFiles.push_back(outFileName.string());
	outFile.open(outFileName, std::ios::binary);

	for (const string& file : filesToProcess)
	{
		snprintf(printingBuffer, 1024, "Identifying possible target sites in : %s", file.c_str());
		printer(printingBuffer);

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
					if (recordedSequences.find(seqHeader) == recordedSequences.end())
					{
						recordedSequences.insert(seqHeader);
						string concatanatedSeq;
						for (string seqFragment : seq) { concatanatedSeq += makeUpper(seqFragment); }

						processSeqeunce(
							concatanatedSeq,
							seqHeader,
							outFile,
							tempWorkingDir,
							numIdentifiedGuides,
							numDuplicateGuides,
							candidateGuides,
							recordedSequences,
							guidesInBatch,
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
			if (recordedSequences.find(seqHeader) == recordedSequences.end())
			{
				recordedSequences.insert(seqHeader);
				string concatanatedSeq;
				for (string seqFragment : seq) { concatanatedSeq += makeUpper(seqFragment); }

				processSeqeunce(
					concatanatedSeq,
					seqHeader,
					outFile,
					tempWorkingDir,
					numIdentifiedGuides,
					numDuplicateGuides,
					candidateGuides,
					recordedSequences,
					guidesInBatch,
					batchSize
				);
			}
			outFile.close();
		}
		// TODO: Other file formats here

		// Report overall progress before processing next file
		completedSizeBytes += std::filesystem::file_size(path(file));
		completedPercent = completedSizeBytes / totalSizeBytes * 100.0;

		snprintf(printingBuffer, 1024, "\tProcessed %.2f%% of input.", completedPercent);
		printer(printingBuffer);

	}

	// Finished processing files, report results
	float duplicatePercent = ((float)numDuplicateGuides / (float)numIdentifiedGuides) * 100.0f;
	snprintf(printingBuffer, 1024,	"\tIdentified %d possible target sites.\n"
									"\tOf these, %d are not unique. These sites occur a total of %d times.\n"
									"\t%d of %d (%.2f%%) of guides will be ignored for optimisation levels over ultralow.\n"
									"\t%d distinct guides were identified.",
									numIdentifiedGuides, (int)duplicateGuides.size(), numDuplicateGuides, numDuplicateGuides, numIdentifiedGuides, duplicatePercent, (int)candidateGuides.size());
	printer(printingBuffer);

	return;
}

void cas9InputProcessor::processSeqeunce(
	const string& seqeunce, 
	const string& seqHeader, 
	ofstream& outFile,
	path& tempWorkingDir,
	int&numIdentifiedGuides,
	int& numDuplicateGuides,
	set<string>& candidateGuides,
	set<string>& recordedSequences,
	int& guidesInBatch,
	const int& batchSize
	)
{

	for (std::sregex_iterator regexItr = std::sregex_iterator(seqeunce.begin(), seqeunce.end(), patternForward);
		regexItr != std::sregex_iterator();
		regexItr++)
	{
		numIdentifiedGuides++;
		std::smatch m = *regexItr;
		string guide = m[1].str();
		int matchPos = m.position();
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
			numIdentifiedGuides++;
			duplicateGuides.insert(guide);
		}
	}

	for (std::sregex_iterator regexItr = std::sregex_iterator(seqeunce.begin(), seqeunce.end(), patternReverse);
		regexItr != std::sregex_iterator();
		regexItr++)
	{
		numIdentifiedGuides++;
		std::smatch m = *regexItr;
		string guide = rc(m[1].str());
		int matchPos = m.position();
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
}

list<string> cas9InputProcessor::getBatchFiles()
{
	return batchFiles;
}

bool cas9InputProcessor::isDuplicateGuide(std::string guide)
{
	return duplicateGuides.find(guide) != duplicateGuides.end();
}

void cas9InputProcessor::cleanUp()
{
	std::filesystem::remove_all(path(std::filesystem::temp_directory_path() / "Crackling"));
}