// CracklingPlusPlus.cpp : Defines the entry point for the application.
//#include "../include/Helpers.hpp"
//#include "../include/Constants.hpp"
#include "../include/Logger.hpp"
//#include "../include/ConfigManager.hpp"
//#include "../include/cas9InputProcessor.hpp"
//#include "../include/CHOPCHOP.hpp"
//#include "../include/mm10db.hpp"
//#include "../include/sgrnascorer2.hpp"
//#include "../include/bowtie2.hpp"
//#include "../include/ISSLOffTargetScoring.hpp"
//#include "../include/ISSLClustering.hpp"
//#include "../include/ISSL2Stage.hpp"
#include "../include/configParserModule.hpp"
#include "../include/cas9InputModule.hpp"
#include "../include/chopchopModule.hpp"
#include "../include/mm10dbModule.hpp"
#include "../include/sgrnascorer2Module.hpp"
#include "../include/bowtie2Module.hpp"
#include "../include/ISSLScoringModule.hpp"
#include "../include/ISSLScoringModuleMMF.hpp"
#include "../include/outputModule.hpp"

int main(int argc, char** argv)
{
	// Check input arguments
	if (argc != 2) {
		std::cout << fmt::format("Usage: {} [Crackling Config File]\n", argv[0]);
		exit(1);
	}

	// Load config
	configParserModule cmNew; 
	cracklingConfig config = cmNew.run(argv[1]);

	// Create logger objects and apply number formatting
	Logger coutLogger(std::cout, config.output.log.string());
	Logger cerrLogger(std::cerr, config.output.errLog.string());
	std::cout.imbue(comma_locale);
	std::cerr.imbue(comma_locale);

	// Create pipeline modules
	cas9InputModule		cas9IM(config);
	chopchopModule		chopchop(config);
	mm10dbModule		mm10db(config);
	sgrnascorer2Module	sgrnascorer2(config);
	bowtie2Module		bowtie2(config);
	ISSLScoringModule	ISSLScoring(config);
	ISSLScoringModuleMMF ISSLScoringMMF(config);
	outputModule		output(config);

	// Record start time
	auto startTime = std::chrono::steady_clock::now();

	cas9IM.run();

	uint64_t i = 0;
	// Process guides in batches
	for (std::vector<guideResults>* currentBatch; currentBatch = cas9IM.next();)
	{
		// Record batch start time
		auto batchStartTime = std::chrono::steady_clock::now();

		std::cout << "Processing batch " << ++i << std::endl;

		// Consensus scoring
		//chopchop.run(*currentBatch);
		//mm10db.run(*currentBatch);
		//sgrnascorer2.run(*currentBatch);

		// Complete consensus evaluation
		//std::cout << "Evaluating efficiency via consensus approach." << std::endl;
		//uint64_t failedCount = 0;
		//uint64_t testedCount = 0;
		//for (guideResults& candidate : *currentBatch)
		//{
		//	candidate.consensusCount = (candidate.passedG20 == CODE_ACCEPTED) + (candidate.acceptedByMm10db == CODE_ACCEPTED) + (candidate.acceptedBySgRnaScorer2 == CODE_ACCEPTED);
		//	if (candidate.consensusCount < config.consensus.n) { failedCount++; }
		//	testedCount++;
		//}
		//std::cout << fmt::format("\t{:L} of {:L} failed here.", failedCount, testedCount) << std::endl;

		// Specificity scoring
		//bowtie2.run(*currentBatch);
		//ISSLScoring.run(*currentBatch);
		ISSLScoringMMF.run(*currentBatch);

		// Print output
		output.run(*currentBatch);

		std::chrono::nanoseconds nanoSec = std::chrono::steady_clock::now() - batchStartTime;
		std::chrono::duration<uint64_t, std::ratio<86400>> days = std::chrono::duration_cast<std::chrono::duration<uint64_t, std::ratio<86400>>>(nanoSec);
		std::chrono::hours hours = std::chrono::duration_cast<std::chrono::hours>(nanoSec - days);
		std::chrono::minutes minutes = std::chrono::duration_cast<std::chrono::minutes>(nanoSec - days - hours);
		std::chrono::seconds sec = std::chrono::duration_cast<std::chrono::seconds>(nanoSec - days - hours - minutes);
		std::chrono::milliseconds milliSec = std::chrono::duration_cast<std::chrono::milliseconds>(nanoSec - days - hours - minutes - sec);
		std::chrono::microseconds microSec = std::chrono::duration_cast<std::chrono::microseconds>(nanoSec - days - hours - minutes - sec - milliSec);
		std::cout << fmt::format("This batch ran in {:02} {:02}:{:02}:{:02} (dd hh:mm:ss) or {} seconds", days.count(), hours.count(), minutes.count(), sec.count(), std::chrono::duration_cast<std::chrono::seconds>(nanoSec).count()) <<  std::endl;

	}
	std::chrono::nanoseconds nanoSec = std::chrono::steady_clock::now() - startTime;
	std::chrono::duration<uint64_t, std::ratio<86400>> days = std::chrono::duration_cast<std::chrono::duration<uint64_t, std::ratio<86400>>>(nanoSec);
	std::chrono::hours hours = std::chrono::duration_cast<std::chrono::hours>(nanoSec - days);
	std::chrono::minutes minutes = std::chrono::duration_cast<std::chrono::minutes>(nanoSec - days - hours);
	std::chrono::seconds sec = std::chrono::duration_cast<std::chrono::seconds>(nanoSec - days - hours - minutes);
	std::chrono::milliseconds milliSec = std::chrono::duration_cast<std::chrono::milliseconds>(nanoSec - days - hours - minutes - sec);
	std::chrono::microseconds microSec = std::chrono::duration_cast<std::chrono::microseconds>(nanoSec - days - hours - minutes - sec - milliSec);
	std::cout << fmt::format("Total run time {:02} {:02}:{:02}:{:02} (dd hh:mm:ss) or {} seconds", days.count(), hours.count(), minutes.count(), sec.count(), std::chrono::duration_cast<std::chrono::seconds>(nanoSec).count()) << std::endl;

	coutLogger.close();
	cerrLogger.close();
	// cas9IM.cleanup();

	//try
	//{
	//	// Check input arguments
	//	if (argc != 2) {
	//		std::cout << fmt::format("Usage: {} [Crackling Config File]\n", argv[0]);
	//		exit(1);
	//	}

	//	// Load config
	//	ConfigManager cm(argv[1]);

	//	configParserModule cmNew; 
	//	cracklingConfig config = cmNew.run(argv[1]);

	//	// Create logger objects
	//	Logger coutLogger(std::cout, cm.getString("output", "log"));
	//	Logger cerrLogger(std::cerr, cm.getString("output", "error"));

	//	// TODO: remove line
	//	cm.set("general", "optimisation", "ultralow");
	//	cm.set("offtargetscore", "method", "mit");

	//	// Record start time
	//	auto start = std::chrono::high_resolution_clock::now();

	//	// Process input
	//	cas9InputProcessor ip;
	//	ip.process(cm.getFilesToProcess(), cm.getInt("input", "batch-size"));

	//	// Create pipeline objects
	//	CHOPCHOP CHOPCHOPModule(cm);
	//	mm10db mm10dbModule(cm);
	//	sgrnascorer2 sgRNAScorer2Module(cm);
	//	bowtie2 bowtie2Module(cm);
	//	ISSLOffTargetScoring OTSModule(cm);
	//	//ISSL2Stage OTSModule(cm);
	//	//ISSLClustering OTSModule(cm);

	//	// Add header line to output file
	//	std::ofstream outFile(cm.getString("output", "file"), std::ios_base::binary | std::ios_base::out);
	//	std::string headerLine;
	//	for (const std::string& guideProperty : DEFAULT_GUIDE_PROPERTIES_ORDER)
	//	{
	//		headerLine += guideProperty + ",";
	//	}

	//	outFile << headerLine.substr(0, headerLine.length() - 1) + "\n";;

	//	outFile.close();

	//	// Start of pipeline
	//	for (const std::string& fileName : ip.getBatchFiles())
	//	{
	//		// Record batch start time
	//		auto batchStart = std::chrono::high_resolution_clock::now();

	//		// Load candidate guides from batch file
	//		std::unordered_map <std::string, std::unordered_map<std::string, std::string>> candidateGuides;
	//		std::ifstream inFile;
	//		inFile.open(fileName, std::ios::binary | std::ios_base::in);

	//		for (std::string line; std::getline(inFile, line);)
	//		{
	//			std::array<std::string, 5> guideInfo;
	//			for (int i = 0; i < 5; i++)
	//			{
	//				size_t pos = line.find(',');
	//				guideInfo[i] = line.substr(0, pos);
	//				line.erase(0, pos + 1);
	//			}

	//			candidateGuides[guideInfo[0]] = DEFAULT_GUIDE_PROPERTIES;
	//			candidateGuides[guideInfo[0]]["seq"] = guideInfo[0];
	//			if (ip.isDuplicateGuide(guideInfo[0]))
	//			{
	//				candidateGuides[guideInfo[0]]["header"] = CODE_AMBIGUOUS;
	//				candidateGuides[guideInfo[0]]["start"] = CODE_AMBIGUOUS;
	//				candidateGuides[guideInfo[0]]["end"] = CODE_AMBIGUOUS;
	//				candidateGuides[guideInfo[0]]["strand"] = CODE_AMBIGUOUS;
	//				candidateGuides[guideInfo[0]]["isUnique"] = CODE_REJECTED;
	//			}
	//			else
	//			{
	//				candidateGuides[guideInfo[0]]["header"] = guideInfo[1];
	//				candidateGuides[guideInfo[0]]["start"] = guideInfo[2];
	//				candidateGuides[guideInfo[0]]["end"] = guideInfo[3];
	//				candidateGuides[guideInfo[0]]["strand"] = guideInfo[4];
	//				candidateGuides[guideInfo[0]]["isUnique"] = CODE_ACCEPTED;
	//			}
	//		}

	//		// TODO: uncomment
	//		// Run scoring modules
	//		//CHOPCHOPModule.run(candidateGuides);

	//		//mm10dbModule.run(candidateGuides);

	//		//sgRNAScorer2Module.run(candidateGuides);

	//		// Complete consensus evaluation
	//		//printer("Evaluating efficiency via consensus approach.");
	//		//int failedCount = 0;
	//		//int testedCount = 0;
	//		//for (const auto& [target23, resultsMap] : candidateGuides)
	//		//{
	//		//	candidateGuides[target23]["consensusCount"] = std::to_string((int)(candidateGuides[target23]["acceptedByMm10db"] == CODE_ACCEPTED) +
	//		//		(int)(candidateGuides[target23]["acceptedBySgRnaScorer"] == CODE_ACCEPTED) +
	//		//		(int)(candidateGuides[target23]["passedG20"] == CODE_ACCEPTED));
	//		//	if (std::stoi(candidateGuides[target23]["consensusCount"]) < cm.getInt("consensus", "n")) { failedCount++; }
	//		//	testedCount++;
	//		//}

	//		//printer(fmt::format("\t{} of {} failed here.", commaify(failedCount), commaify(testedCount)));

	//		// TODO: uncomment
	//		// Run Specifity (Off target) modules 
	//		//bowtie2Module.run(candidateGuides);

	//		OTSModule.run(candidateGuides);

	//		// Write results to output file
	//		printer("Writing results to file.");

	//		std::ofstream resultsFile(cm.getString("output", "file"), std::ios_base::app | std::ios_base::binary);

	//		for (const auto& [target23, resultsMap] : candidateGuides) {
	//			std::string line;
	//			for (std::string guideProperty : DEFAULT_GUIDE_PROPERTIES_ORDER)
	//			{
	//				line += candidateGuides[target23][guideProperty] + ",";
	//			}
	//			line = line.substr(0, line.length() - 1) + "\n";
	//			resultsFile << line;
	//		}

	//		resultsFile.close();

	//		// Clean up
	//		printer("Cleaning auxiliary files.");

	//		std::filesystem::remove(cm.getString("rnafold", "input"));
	//		std::filesystem::remove(cm.getString("rnafold", "output"));
	//		std::filesystem::remove(cm.getString("offtargetscore", "input"));
	//		std::filesystem::remove(cm.getString("offtargetscore", "output"));
	//		std::filesystem::remove(cm.getString("bowtie2", "input"));
	//		std::filesystem::remove(cm.getString("bowtie2", "output"));

	//		printer("Done.");

	//		printer(fmt::format("{} guides evaluated.", commaify((int)candidateGuides.size())));

	//		// Update timing
	//		auto totalSeconds = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - batchStart);

	//		int days = (int)totalSeconds.count() / 86400;
	//		int hours = (totalSeconds.count() % 86400) / 3600;
	//		int minutes = ((totalSeconds.count() % 86400) % 3600) / 60;
	//		int seconds = ((totalSeconds.count() % 86400) % 3600) % 60;

	//		printer(fmt::format("This batch ran in {:02} {:02}:{:02}:{:02} (dd hh:mm:ss) or {} seconds", days, hours, minutes, seconds, (int)totalSeconds.count()));

	//	}

	//	auto totalSeconds = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - start);

	//	int days = (int)totalSeconds.count() / 86400;
	//	int hours = (totalSeconds.count() % 86400) / 3600;
	//	int minutes = ((totalSeconds.count() % 86400) % 3600) / 60;
	//	int seconds = ((totalSeconds.count() % 86400) % 3600) % 60;


	//	printer(fmt::format("Total run time {:02} {:02}:{:02}:{:02} (dd hh:mm:ss) or {} seconds", days, hours, minutes, seconds, (int)totalSeconds.count()));

	//	// Clean up
	//	coutLogger.close();
	//	cerrLogger.close();

	//	ip.cleanUp();

	//	return 0;
	//}
	//catch (const std::exception& error)
	//{
	//	errPrinter(error.what());
	//	// Clean up after error and close gracefully
	//	return -1;
	//}
}

#if defined(_WIN64)
	#pragma pop_macro("close")
#endif