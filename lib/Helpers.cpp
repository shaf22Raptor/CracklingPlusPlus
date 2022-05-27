// Helpers Library
#include <Helpers.hpp>

array<char, 24> nulceotideArray = { 'a', 'c', 'g', 't', 'r', 'y', 'm', 'k', 'b', 'd', 'h', 'v', 'A', 'C', 'G', 'T', 'R', 'Y', 'M', 'K', 'B', 'D', 'H', 'V' };
array<char, 24> complementArray = { 't', 'g', 'c', 'a', 'y', 'r', 'k', 'm', 'v', 'h', 'd', 'b', 'T', 'G', 'C', 'A', 'Y', 'R', 'K', 'M', 'V', 'H', 'D', 'B' };
array<char, 4> atArray = { 'a', 't', 'A', 'T' };

string rc(string DNA)
{
	// Reverse the input seqeuence 
	std::reverse(DNA.begin(), DNA.end());
	// Convert each character to the complement
	for (int i = 0; i < DNA.size(); i++)
	{
		array<char, 24>::iterator nulceotidePos = std::find(nulceotideArray.begin(), nulceotideArray.end(), DNA[i]);
		int complementPos = std::distance(nulceotideArray.begin(), nulceotidePos);
		DNA[i] = complementArray[complementPos];
	}
	return DNA;
}

string transToDNA(string RNA)
{
	// Swap U with T
	std::replace(RNA.begin(), RNA.end(), 'u', 't');
	std::replace(RNA.begin(), RNA.end(), 'U', 'T');
	return RNA;
}

float atPercentage(string seq)
{
	float total = 0.0f;
	float length = seq.size();
	array<char, 4>::iterator p;
	for (int i = 0; i < seq.size(); i++)
	{
		// Check if the char at the current pos is present in the 'AT' array
		p = std::find(atArray.begin(), atArray.end(), seq[i]);
		if (p != atArray.end())
		{
			total++;
		}
	}

	return (100.0f * total / length);
}

void printer(string formattedString)
{
	time_t rawtime = time(0);
	struct tm* timeinfo = localtime(&rawtime);
	char timestampBuffer[32];
	strftime(timestampBuffer, 32, ">>> %Y-%m-%d %H:%M:%S:\t", timeinfo);
	
	std::cout << timestampBuffer << formattedString << std::endl;
	return;
}

void runner(char* args)
{
	int returnCode;

	returnCode = system(args);

	printer(std::to_string(returnCode));
}
