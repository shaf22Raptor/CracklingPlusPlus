// CracklingPlusPlus.cpp : Defines the entry point for the application.

#include <Helpers.hpp>
#include <Constants.hpp>
#include <Logger.hpp>
#include <ConfigManager.hpp>
#include <cas9InputProcessor.hpp>
#include <CHOPCHOP.hpp>
#include <mm10db.hpp>
#include <sgrnascorer2.hpp>
#include <bowtie2.hpp>
#include <offTargetScoring.hpp>


int main(int argc, char** argv)
{
	try
	{
		// Check input arguments
		if (argc != 2) {
			fprintf(stderr, "Usage: %s [Crackling Config File]\n", argv[0]);
			exit(1);
		}

		// Load config
		ConfigManager cm(argv[1]);

		// Create logger objects
		Logger coutLogger(std::cout, cm.getString("output", "log"));
		Logger cerrLogger(std::cerr, cm.getString("output", "error"));

		// Process input
		cas9InputProcessor ip;
		std::list<std::string> batchFiles = ip.processInput(cm.getFilesToProcess(), cm.getInt("input", "batch-size"));

		// Create pipeline objects
		CHOPCHOP CHOPCHOPModule(cm);
		mm10db mm10dbModule(cm);
		sgrnascorer2 sgRNAScorer2Module(cm);
		bowtie2 bowtie2Module(cm);
		offTargetScoring otsModule(cm);

		// Add header line to output file
		std::ofstream outFile(cm.getString("output", "file"), std::ios_base::binary);
		std::string headerLine;
		for (std::string guideProperty : DEFAULT_GUIDE_PROPERTIES_ORDER)
		{
			headerLine += guideProperty + ",";
		}
		headerLine = headerLine.substr(0, headerLine.length() - 1) + "\n";

		outFile << headerLine;

		outFile.close();

		// Start of pipeline
		for (std::string fileName : batchFiles)
		{
			std::map <std::string, std::map<std::string, std::string>> candidateGuides;
			std::ifstream inFile;
			inFile.open(fileName);

			for (std::string line; std::getline(inFile, line);)
				//while (std::getline(inFile, inputLine))
			{
				std::string guideInfo[5];
				for (int i = 0; i < 5; i++)
				{
					size_t pos = line.find(',');
					guideInfo[i] = line.substr(0, pos);
					line.erase(0, pos + 1);
				}

				candidateGuides[guideInfo[0]] = DEFAULT_GUIDE_PROPERTIES;
				candidateGuides[guideInfo[0]]["seq"] = guideInfo[0];
				if (ip.duplicateGuides.find(guideInfo[0]) != ip.duplicateGuides.end())
				{
					candidateGuides[guideInfo[0]]["header"] = CODE_AMBIGUOUS;
					candidateGuides[guideInfo[0]]["start"] = CODE_AMBIGUOUS;
					candidateGuides[guideInfo[0]]["end"] = CODE_AMBIGUOUS;
					candidateGuides[guideInfo[0]]["strand"] = CODE_AMBIGUOUS;
					candidateGuides[guideInfo[0]]["isUnique"] = CODE_REJECTED;
				}
				else
				{
					candidateGuides[guideInfo[0]]["header"] = guideInfo[1];
					candidateGuides[guideInfo[0]]["start"] = guideInfo[2];
					candidateGuides[guideInfo[0]]["end"] = guideInfo[3];
					candidateGuides[guideInfo[0]]["strand"] = guideInfo[4];
				}
			}

			CHOPCHOPModule.run(candidateGuides);


			mm10dbModule.run(candidateGuides);

			
			sgRNAScorer2Module.run(candidateGuides);

			printer("Evaluating efficiency via consensus approach.");
			int failedCount = 0;
			int testedCount = 0;
			for (auto const& [target23, resultsMap] : candidateGuides)
			{
				candidateGuides[target23]["consensusCount"] = std::to_string((candidateGuides[target23]["acceptedByMm10db"] == CODE_ACCEPTED) +
					(candidateGuides[target23]["acceptedBySgRnaScorer"] == CODE_ACCEPTED) +
					(candidateGuides[target23]["passedG20"] == CODE_ACCEPTED));
				if (std::stoi(candidateGuides[target23]["consensusCount"]) < cm.getInt("consensus", "n")) { failedCount++; }
				testedCount++;
			}
			char printingBuffer[1024];
			snprintf(printingBuffer, 1024, "\t%d of %d failed here.", failedCount, testedCount);
			printer(printingBuffer);

			bowtie2Module.run(candidateGuides);



			otsModule.run(candidateGuides);



			printer("Writing results to file.");


			std::ofstream outFile(cm.getString("output", "file"), std::ios_base::app | std::ios_base::binary);
			std::string headerLine;

			for (std::map<std::string, std::map<std::string, std::string>>::iterator iter = candidateGuides.begin(); iter != candidateGuides.end(); iter++)
			{
				std::string target23 = iter->first;
				std::string line;
				for (std::string guideProperty : DEFAULT_GUIDE_PROPERTIES_ORDER)
				{
					line += candidateGuides[target23][guideProperty] + ",";
				}
				line = line.substr(0, line.length() - 1) + "\n";
				outFile << line;
			}

			outFile.close();

			printer("Cleaning auxiliary files.");

			std::filesystem::remove(cm.getString("rnafold", "input"));
			std::filesystem::remove(cm.getString("rnafold", "output"));
			std::filesystem::remove(cm.getString("offtargetscore", "input"));
			std::filesystem::remove(cm.getString("offtargetscore", "output"));
			std::filesystem::remove(cm.getString("bowtie2", "input"));
			std::filesystem::remove(cm.getString("bowtie2", "output"));

			printer("Done.");

			snprintf(printingBuffer, 1024, "%d guides evaluated.", (int)candidateGuides.size());
			printer(printingBuffer);

		}


		// Clean up
		coutLogger.close();
		cerrLogger.close();
		
		return 0;
	}
	catch (const std::exception& error)
	{
		errPrinter(error.what());
		// Clean up after error and close gracefully
		return -1;
	}
}
