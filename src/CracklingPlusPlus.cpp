// CracklingPlusPlus.cpp : Defines the entry point for the application.

#include <Helpers.h>


int main()
{
	printer("Starting OUTPUT");
	printer(rc("ATGC"));
	printer(transToDNA("AAAATTTTGGGGUUUU"));
	printer(std::to_string(atPercentage("AAAATTTTGGGGCCC")));
}
