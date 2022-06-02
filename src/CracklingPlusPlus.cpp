// CracklingPlusPlus.cpp : Defines the entry point for the application.

#include <Helpers.hpp>
#include <ConfigManager.hpp>


int main(int argc, char** argv)
{
	// Check input arguments

	// Load config
	
	ConfigManager cm("test-test_config.ini");

	// Run crackling
	printer("Starting crackling");
	printer("Finished");
}
