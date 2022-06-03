// CracklingPlusPlus.cpp : Defines the entry point for the application.

#include <Helpers.hpp>
#include <ConfigManager.hpp>
#include <Constants.hpp>
#include <Logger.hpp>
#include <Windows.h>

// "C:\\Users\\n9725059\\Documents\\Crackling-fork\\Crackling\\log.log"



int main(int argc, char** argv)
{
	// Check input arguments

	// Load config
	
	ConfigManager cm("test_config.ini");

	// Redirect cout and cerr
	
	// Store old buffers
	std::streambuf* oldCout = std::cout.rdbuf();
	std::streambuf* oldCerr = std::cerr.rdbuf();

	// Open file streams
	//std::ofstream stdLog;
	//stdLog.open(cm.getString("output", "log"));
	//std::ofstream errLog;
	//errLog.open(cm.getString("output", "error"));

	// Create logger objects
	Logger coutLogger(cm.getString("output", "log"));
	std::ostream newCout(&coutLogger);
	Logger cerrLogger(cm.getString("output", "error"));
	std::ostream newCerr(&cerrLogger);

	// Set cout and cerr to custom logger
	std::cout.rdbuf(newCout.rdbuf());
	std::cerr.rdbuf(newCerr.rdbuf());

	// Run crackling
	printer("Starting crackling");

	printer(CODE_ACCEPTED);
	printer(CODE_REJECTED);
	printer(CODE_UNTESTED);
	printer(CODE_AMBIGUOUS);
	Sleep(5000);
	printer("Finished");

	// Reset cout and cerr
	std::cout.rdbuf(oldCout); 
	std::cout.rdbuf(oldCerr);

	return 0;
}
