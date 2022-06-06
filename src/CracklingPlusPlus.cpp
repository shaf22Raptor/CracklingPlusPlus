// CracklingPlusPlus.cpp : Defines the entry point for the application.

#include <Helpers.hpp>
#include <ConfigManager.hpp>
#include <Constants.hpp>
#include <Logger.hpp>

// "C:\\Users\\n9725059\\Documents\\Crackling-fork\\Crackling\\log.log"



int main(int argc, char** argv)
{
	// Check input arguments

	// Load config
	
	ConfigManager cm("test_config.ini");

	// Create logger objects
	Logger coutLogger( std::cout , cm.getString("output", "log"));
	Logger cerrLogger( std::cerr , cm.getString("output", "error"));


	// Run crackling
	printer("Starting crackling");
	printer(CODE_ACCEPTED);
	printer(CODE_REJECTED);
	printer(CODE_UNTESTED);
	printer(CODE_AMBIGUOUS);
	printer("Finished");

	// Clean up
	coutLogger.close();

	return 0;
}
