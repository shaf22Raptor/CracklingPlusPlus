// CracklingPlusPlus.cpp : Defines the entry point for the application.

#include <Helpers.h>


int main()
{
	printer("Starting OUTPUT");
	std::cout << rc("ATGC") << "\n";
	std::cout << transToDNA("AAAATTTTGGGGUUUU") << "\n";
	std::cout << atPercentage("AAAATTTTGGGGCCC") << "\n";
	return 0;
}
