#include <cas9InputProcessor.hpp>

using std::string;
using std::list;
using std::set;

void cas9InputProcessor::processInput(list<string> filesToProcess, int batchSize)
{
	printer("Analysing files...");

	// Guide searching
	std::regex patternForward("(?=([ATCG]{21}GG))");
	std::regex patternReverse("(?=(CC[ACGT]{21}))");

	// Formatting prints
	char printingBuffer[1024];

	// In file
	std::ifstream inFile;
	string inputLine;

	// Batched files
	std::ofstream outFile;
	std::filesystem::path outFileName;
	int guidesInBatch = 0;

	// File processing
	string seqHeader;
	list<string> seq;

	// Duplicate tracking
	set<string> candidateGuides;
	set<string> recordedSequences;
	int numDuplicateGuides = 0;
	int	numIdentifiedGuides = 0;

	// Setup temp working dir
	std::filesystem::path systemTempDir = std::filesystem::temp_directory_path();
	if (std::filesystem::is_directory(systemTempDir / "Crackling")) { std::filesystem::remove_all(systemTempDir / "Crackling"); }
	if (!std::filesystem::create_directory(systemTempDir / "Crackling")) { throw std::runtime_error("Unable to create temp working dir!"); }
	std::filesystem::path tempWorkingDir(systemTempDir / "Crackling");


	snprintf(printingBuffer, 1024, "Storing batch files in: %s", tempWorkingDir.string().c_str());
	printer(printingBuffer);

	// Create first batch file
	outFileName = tempWorkingDir / (std::to_string(batchFiles.size()) + "_batchFile.txt");
	batchFiles.push_back(outFileName.string());
	outFile.open(outFileName);

	for (const string& file : filesToProcess)
	{
		snprintf(printingBuffer, 1024, "Identifying possible target sites in : %s", file.c_str());
		printer(printingBuffer);

		inFile.open(file);
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

						for (std::sregex_iterator regexItr = std::sregex_iterator(concatanatedSeq.begin(), concatanatedSeq.end(), patternForward);
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
									outFileName = tempWorkingDir / (std::to_string(batchFiles.size()) + "_batchFile.txt");
									batchFiles.push_back(outFileName.string());
									outFile.open(outFileName);
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

						for (std::sregex_iterator regexItr = std::sregex_iterator(concatanatedSeq.begin(), concatanatedSeq.end(), patternReverse);
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
									outFileName = tempWorkingDir / (std::to_string(batchFiles.size()) + "_batchFile.txt");
									batchFiles.push_back(outFileName.string());
									outFile.open(outFileName);
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

				for (std::sregex_iterator regexItr = std::sregex_iterator(concatanatedSeq.begin(), concatanatedSeq.end(), patternForward);
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
							outFileName = tempWorkingDir / (std::to_string(batchFiles.size()) + "_batchFile.txt");
							batchFiles.push_back(outFileName.string());
							outFile.open(outFileName);
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

				for (std::sregex_iterator regexItr = std::sregex_iterator(concatanatedSeq.begin(), concatanatedSeq.end(), patternReverse);
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
							outFileName = tempWorkingDir / (std::to_string(batchFiles.size()) + "_batchFile.txt");
							batchFiles.push_back(outFileName.string());
							outFile.open(outFileName);
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
			outFile.close();
		}
		// TODO: Other file formats here

		snprintf(printingBuffer, 1024, "\tExtracted from %d%% of input.", numDuplicateGuides);
		printer(printingBuffer);

	}

	float duplicatePercent = ((float)numDuplicateGuides / (float)numIdentifiedGuides) * 100.0f;
	snprintf(printingBuffer, 1024, "\tIdentified %d possible target sites.", numIdentifiedGuides);
	printer(printingBuffer);
	snprintf(printingBuffer, 1024, "\tOf these, %d are not unique. These sites occur a total of %d times.", (int)duplicateGuides.size(), numDuplicateGuides);
	printer(printingBuffer);
	snprintf(printingBuffer, 1024, "\t%d of %d (%.2f%%) of guides will be ignored for optimisation levels over ultralow",numDuplicateGuides, numIdentifiedGuides, duplicatePercent);
	printer(printingBuffer);
	snprintf(printingBuffer, 1024, "\t%d distinct guides were identified.", (int)candidateGuides.size());
	printer(printingBuffer);

	return;
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
	std::filesystem::remove_all(std::filesystem::path(std::filesystem::temp_directory_path() / "Crackling"));
}