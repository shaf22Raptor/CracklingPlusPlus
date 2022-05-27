// CracklingPlusPlus.cpp : Defines the entry point for the application.

#include <Helpers.hpp>


int main(int argc, char** argv)
{
	printer("Starting OUTPUT");
	printer(rc("ATGC"));
	printer(transToDNA("AAAATTTTGGGGUUUU"));
	printer(std::to_string(atPercentage("AAAATTTTGGGGCCC")));
	char argsBuffer[1024];
	snprintf(argsBuffer, 1024, "%s", "dir");
	runner(argsBuffer);
	runner(argsBuffer);
}
