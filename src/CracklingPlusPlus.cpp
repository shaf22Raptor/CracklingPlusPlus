// CracklingPlusPlus.cpp : Defines the entry point for the application.

#include <Helpers.hpp>
#include <ConfigManager.hpp>
#include <Constants.hpp>
#include <Logger.hpp>
#include <CHOPCHOP.hpp>
#include <mm10db.hpp>
#include <cas9InputProcessor.hpp>

int main(int argc, char** argv)
{
	// Check input arguments

	// Load config
	
	ConfigManager cm("data/test_config.ini");

	// Create logger objects
	Logger coutLogger( std::cout , cm.getString("output", "log"));
	Logger cerrLogger( std::cerr , cm.getString("output", "error"));

	// Process input
	cas9InputProcessor ip;
	std::list<std::string> batchFiles =  ip.processInput(cm.getFilesToProcess(), cm.getInt("input", "batch-size"));

	// Create pipeline objects
	CHOPCHOP CHOPCHOPModule(cm);
	mm10db mm10dbModule(cm);

	// Start of pipeline
	for (std::string fileName : batchFiles)
	{
		std::map <std::string, std::map<std::string, std::string>> candidateGuides;
		std::ifstream inFile;
		std::string inputLine;
		inFile.open(fileName);
		while (std::getline(inFile, inputLine))
		{
			std::string guideInfo[5];
			for (int i = 0; i < 5; i++)
			{
				size_t pos = inputLine.find(',');
				guideInfo[i] = inputLine.substr(0, pos);
				inputLine.erase(0, pos + 1);
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

		std::cout << "Finished CHOPCHOP" << std::endl;

		mm10dbModule.run(candidateGuides);

		std::cout << "Finished mm10db" << std::endl;

	}


	// Clean up
	coutLogger.close();
	cerrLogger.close();

	return 0;
}
