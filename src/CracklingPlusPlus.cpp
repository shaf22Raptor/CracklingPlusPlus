// CracklingPlusPlus.cpp : Defines the entry point for the application.

#include <Helpers.hpp>
#include <ConfigManager.hpp>
#include <Constants.hpp>


int main(int argc, char** argv)
{
	// Check input arguments

	// Load config
	
	ConfigManager cm("test-test_config.ini");

	// Run crackling
	printer("Starting crackling");
	printer("Finished");
	printer(CODE_ACCEPTED);
	printer(CODE_REJECTED);
	printer(CODE_UNTESTED);
	printer(CODE_AMBIGUOUS);
}
