#include <inputProcessor.hpp>

using std::string;
using std::list;
using std::set;

list<string> processInput(list<string>& filesToProcess)
{
	// Guide searching
	std::regex patternForward("[ATCG]{21}GG");
	std::regex patternReverse("CC[ACGT]{21}");

	// Formatting prints
	char printingBuffer[1024];

	// In file
	std::ifstream inFile;
	string inputLine;

	// Batched files
	std::ofstream outFile;
	string outFileName;
	list<string> batchFiles;


	// File processing
	string seqHeader;
	list<string> seq;
	
	// Duplicate tracking
	set<string> candidateGuides;
	set<string> duplicateGuides;
	set<string> recordedSequences;

	// TODO: move these assingments to constructor
	int guidesInBatch;
	int batchSize;

	

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
			while (std::getline(inFile, inputLine)) {
				inputLine = inputLine.substr(0, inputLine.find("\n"));
				if (inputLine[0] == '>')
				{
					if (recordedSequences.find(seqHeader) == recordedSequences.end())
					{
						string concatanatedSeq;
						for (string seqFragment : seq) { concatanatedSeq += seqFragment; }

						std::regex_iterator<string::iterator> regexItr(concatanatedSeq.begin(), concatanatedSeq.end(), patternForward);
						std::regex_iterator<string::iterator> regexItrEnd;
						while (regexItr != regexItrEnd) {
							outFile << regexItr->str() << "," << seqHeader << "," << regexItr->position() << "," << (regexItr->position() + 23) << "," << "+" << "\n";
							regexItr++;
						}

						regexItr = std::regex_iterator<string::iterator>(concatanatedSeq.begin(), concatanatedSeq.end(), patternReverse);
						while (regexItr != regexItrEnd) {
							outFile << regexItr->str() << "," << seqHeader << "," << regexItr->position() << "," << (regexItr->position() + 23) << "," << "-" << "\n";
							regexItr++;
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
				string concatanatedSeq;
				for (string seqFragment : seq)	{ concatanatedSeq += seqFragment; }

				std::regex_iterator<string::iterator> regexItr(concatanatedSeq.begin(), concatanatedSeq.end(), patternForward);
				std::regex_iterator<string::iterator> regexItrEnd;
				while (regexItr != regexItrEnd) {
					outFile << regexItr->str() << "," << seqHeader << "," << regexItr->position() << "," << (regexItr->position() + 23) << "," << "+" << "\n";
					regexItr++;
				}

				regexItr = std::regex_iterator<string::iterator>(concatanatedSeq.begin(), concatanatedSeq.end(), patternReverse);
				while (regexItr != regexItrEnd) {
					outFile << regexItr->str() << "," << seqHeader << "," << regexItr->position() << "," << (regexItr->position() + 23) << "," << "-" << "\n";
					regexItr++;
				}
			}
		}
		// TODO: Other file formats here

	}
	return {};
}