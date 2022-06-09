// Helpers Library
#include <Helpers.hpp>

using std::string;
using std::array;

array<char, 24> nulceotideArray = { 'a', 'c', 'g', 't', 'r', 'y', 'm', 'k', 'b', 'd', 'h', 'v', 'A', 'C', 'G', 'T', 'R', 'Y', 'M', 'K', 'B', 'D', 'H', 'V' };
array<char, 24> complementArray = { 't', 'g', 'c', 'a', 'y', 'r', 'k', 'm', 'v', 'h', 'd', 'b', 'T', 'G', 'C', 'A', 'Y', 'R', 'K', 'M', 'V', 'H', 'D', 'B' };

string rc(string DNA)
{
	// Ensure input lenght is greater than 0
	if (DNA.length() < 1)
	{
		throw std::length_error("Type Error, Seqeunce length must be greater than 0!");
	}
	// Ensure input lenght is not greater than 1024
	else if (DNA.length() > 1024)
	{
		throw std::length_error("Type Error, Seqeunce length must be less than 1024!");
	}
	// Reverse the input seqeuence 
	std::reverse(DNA.begin(), DNA.end());
	// Convert each character to the complement
	for (int i = 0; i < DNA.size(); i++)
	{
		array<char, 24>::iterator nulceotidePos = std::find(nulceotideArray.begin(), nulceotideArray.end(), DNA[i]);
		int complementPos = std::distance(nulceotideArray.begin(), nulceotidePos);
		DNA[i] = complementArray[complementPos];
	}
	// Return reverse compliment
	return DNA;
}

void printer(string formattedString)
{
	std::cout << formattedString << std::endl;
}

void errPrinter(string formattedString)
{
	std::cerr << formattedString << std::endl;
}

void runner(char* args)
{
	char buffer[1024];
	snprintf(buffer, 1024, "| Calling: %s", args);
	printer(buffer);
	try 
	{
		int returnCode = system(args);
		if (returnCode != 0)
		{
			throw std::runtime_error("Runtime Error, returned a normal 0 value!");
		}
	}
	catch (string error)
	{
		errPrinter(error);
		return;
	}
    printer("| Finished");
	return;
}
