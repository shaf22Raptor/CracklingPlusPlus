// CracklingPlusPlus.cpp : Defines the entry point for the application.

#include "../include/Helpers.hpp"
#include "../include/Constants.hpp"
#include "../include/Logger.hpp"
#include "../include/ConfigManager.hpp"
#include "../include/cas9InputProcessor.hpp"
#include "../include/CHOPCHOP.hpp"
#include "../include/mm10db.hpp"
#include "../include/sgrnascorer2.hpp"
#include "../include/bowtie2.hpp"
#include "../include/ISSLOffTargetScoring.hpp"

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
#pragma push_macro(close)
#undef close
#endif

int main(int argc, char** argv)
{
	try
	{
		// Check input arguments
		if (argc != 2) {
			std::cout << fmt::format("Usage: {} [Crackling Config File]\n", argv[0]);
			exit(1);
		}

		// Load config
		ConfigManager cm(argv[1]);

		// Create logger objects
		Logger coutLogger(std::cout, cm.getString("output", "log"));
		Logger cerrLogger(std::cerr, cm.getString("output", "error"));

		// Record start time
		auto start = std::chrono::high_resolution_clock::now();

		// Process input
		cas9InputProcessor ip;
		ip.process(cm.getFilesToProcess(), cm.getInt("input", "batch-size"));

		// Create pipeline objects
		CHOPCHOP CHOPCHOPModule(cm);
		mm10db mm10dbModule(cm);
		sgrnascorer2 sgRNAScorer2Module(cm);
		bowtie2 bowtie2Module(cm);
		ISSLOffTargetScoring OTSModule(cm);

		// Add header line to output file
		std::ofstream outFile(cm.getString("output", "file"), std::ios_base::binary | std::ios_base::out);
		std::string headerLine;
		for (const std::string& guideProperty : DEFAULT_GUIDE_PROPERTIES_ORDER)
		{
			headerLine += guideProperty + ",";
		}

		outFile << headerLine.substr(0, headerLine.length() - 1) + "\n";;

		outFile.close();

		// Start of pipeline
		for (const std::string& fileName : ip.getBatchFiles())
		{
			// Record batch start time
			auto batchStart = std::chrono::high_resolution_clock::now();

			std::unordered_map <std::string, std::unordered_map<std::string, std::string>> candidateGuides;
			std::ifstream inFile;
			inFile.open(fileName, std::ios::binary | std::ios_base::in);

			for (std::string line; std::getline(inFile, line);)
			{
				std::array<std::string, 5> guideInfo;
				for (int i = 0; i < 5; i++)
				{
					size_t pos = line.find(',');
					guideInfo[i] = line.substr(0, pos);
					line.erase(0, pos + 1);
				}

				candidateGuides[guideInfo[0]] = DEFAULT_GUIDE_PROPERTIES;
				candidateGuides[guideInfo[0]]["seq"] = guideInfo[0];
				if (ip.isDuplicateGuide(guideInfo[0]))
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
					candidateGuides[guideInfo[0]]["isUnique"] = CODE_ACCEPTED;
				}
			}

			CHOPCHOPModule.run(candidateGuides);

			mm10dbModule.run(candidateGuides);

			sgRNAScorer2Module.run(candidateGuides);

			printer("Evaluating efficiency via consensus approach.");
			int failedCount = 0;
			int testedCount = 0;
			for (const auto& [target23, resultsMap] : candidateGuides)
			{
				candidateGuides[target23]["consensusCount"] = std::to_string((int)(candidateGuides[target23]["acceptedByMm10db"] == CODE_ACCEPTED) +
					(int)(candidateGuides[target23]["acceptedBySgRnaScorer"] == CODE_ACCEPTED) +
					(int)(candidateGuides[target23]["passedG20"] == CODE_ACCEPTED));
				if (std::stoi(candidateGuides[target23]["consensusCount"]) < cm.getInt("consensus", "n")) { failedCount++; }
				testedCount++;
			}

			printer(fmt::format("\t{} of {} failed here.", commaify(failedCount), commaify(testedCount)));

			bowtie2Module.run(candidateGuides);

			OTSModule.run(candidateGuides);

			printer("Writing results to file.");

			std::ofstream resultsFile(cm.getString("output", "file"), std::ios_base::app | std::ios_base::binary);

			for (const auto& [target23, resultsMap] : candidateGuides) {
				std::string line;
				for (std::string guideProperty : DEFAULT_GUIDE_PROPERTIES_ORDER)
				{
					line += candidateGuides[target23][guideProperty] + ",";
				}
				line = line.substr(0, line.length() - 1) + "\n";
				resultsFile << line;
			}

			resultsFile.close();

			printer("Cleaning auxiliary files.");

			std::filesystem::remove(cm.getString("rnafold", "input"));
			std::filesystem::remove(cm.getString("rnafold", "output"));
			std::filesystem::remove(cm.getString("offtargetscore", "input"));
			std::filesystem::remove(cm.getString("offtargetscore", "output"));
			std::filesystem::remove(cm.getString("bowtie2", "input"));
			std::filesystem::remove(cm.getString("bowtie2", "output"));

			printer("Done.");

			printer(fmt::format("{} guides evaluated.", commaify((int)candidateGuides.size())));

			auto totalSeconds = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - batchStart);

			int days = (int)totalSeconds.count() / 86400;
			int hours = (totalSeconds.count() % 86400) / 3600;
			int minutes = ((totalSeconds.count() % 86400) % 3600) / 60;
			int seconds = ((totalSeconds.count() % 86400) % 3600) % 60;

			printer(fmt::format("This batch ran in {:02} {:02}:{:02}:{:02} (dd hh:mm:ss) or {} seconds", days, hours, minutes, seconds, (int)totalSeconds.count()));

		}

		auto totalSeconds = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - start);

		int days = (int)totalSeconds.count() / 86400;
		int hours = (totalSeconds.count() % 86400) / 3600;
		int minutes = ((totalSeconds.count() % 86400) % 3600) / 60;
		int seconds = ((totalSeconds.count() % 86400) % 3600) % 60;


		printer(fmt::format("Total run time {:02} {:02}:{:02}:{:02} (dd hh:mm:ss) or {} seconds", days, hours, minutes, seconds, (int)totalSeconds.count()));

		// Clean up
		coutLogger.close();
		cerrLogger.close();

		ip.cleanUp();

		return 0;
	}
	catch (const std::exception& error)
	{
		errPrinter(error.what());
		// Clean up after error and close gracefully
		return -1;
	}
}

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
#pragma pop_macro(close)
#endif