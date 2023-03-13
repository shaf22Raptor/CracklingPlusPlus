#include "../include/util.hpp"

const std::locale comma_locale(std::locale(), new commaFormat());


const std::unordered_map<std::string, optimisationLevel> optimisationMap = {
	{"ultralow",optimisationLevel::ultralow},
	{"low",optimisationLevel::low},
	{"medium",optimisationLevel::medium},
	{"high",optimisationLevel::high}
};

const std::unordered_map<std::string, otScoreMethod> otScoreMethodMap = {
	{"mit",otScoreMethod::mit},
	{"cfd",otScoreMethod::cfd},
	{"mitAndCfd",otScoreMethod::mitAndCfd},
	{"mitOrCfd",otScoreMethod::mitOrCfd},
	{"avgMitCfd",otScoreMethod::avgMitCfd}
};

const char CODE_ACCEPTED = '1';
const char CODE_REJECTED = '0';
const char CODE_UNTESTED = '?';
const char CODE_AMBIGUOUS = '-';
const char CODE_ERROR = '!';

void runner(const char* args)
{
	std::cout << fmt::format("| calling: {}", args) << std::endl;	
	try
	{
		int returncode = system(args);
		if (returncode != 0)
		{
			throw ReturnCode();
		}
	}
	catch (const std::exception& e)
	{
		std::cerr << e.what() << std::endl;
		return;
	}
	std::cout << "| finished" << std::endl;
	return;
}

const std::vector<char> nucleotideArray = { 'a', 'c', 'g', 't', 'r', 'y', 'm', 'k', 'b', 'd', 'h', 'v', 'A', 'C', 'G', 'T', 'R', 'Y', 'M', 'K', 'B', 'D', 'H', 'V' };
const std::vector<char> complementArray = { 't', 'g', 'c', 'a', 'y', 'r', 'k', 'm', 'v', 'h', 'd', 'b', 'T', 'G', 'C', 'A', 'Y', 'R', 'K', 'M', 'V', 'H', 'D', 'B' };

std::string rc(std::string DNA)
{
	// Reverse the input seqeuence 
	std::reverse(DNA.begin(), DNA.end());
	// Convert each character to the complement
	std::for_each(DNA.begin(), DNA.end(), [](char& c) {
		auto nulceotidePos = std::find(nucleotideArray.begin(), nucleotideArray.end(), c);
		long long complementPos = std::distance(nucleotideArray.begin(), nulceotidePos);
		c = complementArray[complementPos];
		}
	);

	// Return reverse compliment
	return DNA;
}
