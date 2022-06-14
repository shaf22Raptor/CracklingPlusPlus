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

	// std::list<std::string> batchFiles = processInput(cm.getFilesToProcess(), cm.getInt("input", "batch-size"));
	cas9InputProcessor ip;
	std::list<std::string> batchFiles =  ip.processInput(cm.getFilesToProcess(), cm.getInt("input", "batch-size"));


	// Clean up
	coutLogger.close();
	cerrLogger.close();

	return 0;
}
