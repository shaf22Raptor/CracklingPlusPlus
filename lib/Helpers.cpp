// Helpers Library
#include <Helpers.hpp>

using std::string;
using std::array;
using std::map;
using std::stoi;
using std::vector;
using std::string_view;

const array<char, 24> nulceotideArray = { 'a', 'c', 'g', 't', 'r', 'y', 'm', 'k', 'b', 'd', 'h', 'v', 'A', 'C', 'G', 'T', 'R', 'Y', 'M', 'K', 'B', 'D', 'H', 'V' };
const array<char, 24> complementArray = { 't', 'g', 'c', 'a', 'y', 'r', 'k', 'm', 'v', 'h', 'd', 'b', 'T', 'G', 'C', 'A', 'Y', 'R', 'K', 'M', 'V', 'H', 'D', 'B' };
const string WHITESPACE = " \n\r\t\f\v";


string makeUpper(const string& s)
{
	string upper = s;
	std::transform(upper.begin(), upper.end(), upper.begin(), [](unsigned char c) { return std::toupper(c); });
	return upper;
}

string makeLower(const string& s)
{
	string lower = s;
	std::transform(lower.begin(), lower.end(), lower.begin(), [](unsigned char c) { return std::tolower(c); });
	return lower;
}

vector<string> split(string& s, string_view delimiter)
{
	vector<string> splitLine;

	size_t delimiterPos = 0;
	while ((delimiterPos = s.find(delimiter)) != std::string::npos) {
		splitLine.push_back(s.substr(0, delimiterPos));
		s.erase(0, delimiterPos + 1);
	}
	splitLine.push_back(s.substr(0, delimiterPos));

	return splitLine;
}


string rtrim(string_view s)
{
	size_t end = s.find_last_not_of(WHITESPACE);
	if (end == string::npos) { return ""; }
	return { s.begin(), s.begin() + end + 1 };
}

string ltrim(string_view s)
{
	size_t start = s.find_first_not_of(WHITESPACE);
	if (start == string::npos) { return ""; }
	return { s.begin() + start, s.end()};
}

string trim(const string& s)
{
	return ltrim(rtrim(s));
}

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
	std::for_each(DNA.begin(), DNA.end(), [](char& c) { 
		auto nulceotidePos = std::find(nulceotideArray.begin(), nulceotideArray.end(), c);
		long long complementPos = std::distance(nulceotideArray.begin(), nulceotidePos);
		c = complementArray[complementPos];
		}
	);

	// Return reverse compliment
	return DNA;
}

bool filterCandidateGuides(map<string, string, std::less<>> candidateGuideResultMap, string_view selectedModule, string_view optimisation, int consensusN, int toolCount)
{
	// Ultralow optimisation, process all guides
	if (optimisation == "ultralow") { return true; }

	// For all modules
	if (optimisation == "low" || optimisation == "medium" || optimisation == "high")
	{	
		// Reject all guides that have been seen more than once
		if (candidateGuideResultMap["isUnique"] == CODE_REJECTED) { return false; }
	}

	// mm10db filtering
	if ( (selectedModule == MODULE_MM10DB) && (optimisation == "medium" || optimisation == "high") )
	{
		// Reject if any mm10db test has failed
		if (candidateGuideResultMap["passedAvoidLeadingT"] == CODE_REJECTED ||
			candidateGuideResultMap["passedATPercent"] == CODE_REJECTED ||
			candidateGuideResultMap["passedTTTT"] == CODE_REJECTED ||
			candidateGuideResultMap["passedSecondaryStructure"] == CODE_REJECTED ||
			candidateGuideResultMap["acceptedByMm10db"] == CODE_REJECTED)
		{
			return false;
		}
	}

	// For all consensus scoring tools
	if ( (selectedModule == MODULE_CHOPCHOP || selectedModule == MODULE_MM10DB || selectedModule == MODULE_SGRNASCORER2) && 
		(optimisation == "high") )
	{
		int countAlreadyAccepted = 
			(candidateGuideResultMap["acceptedByMm10db"] == CODE_ACCEPTED) + 
			(candidateGuideResultMap["passedG20"] == CODE_ACCEPTED) +
			(candidateGuideResultMap["acceptedBySgRnaScorer"] == CODE_ACCEPTED);

		int countAlreadyAssessed = 
			(candidateGuideResultMap["acceptedByMm10db"] != CODE_UNTESTED) +
			(candidateGuideResultMap["passedG20"] != CODE_UNTESTED) +
			(candidateGuideResultMap["acceptedBySgRnaScorer"] != CODE_UNTESTED);

		// Reject if the consensus has already been passed
		if (countAlreadyAccepted >= consensusN) { return false; }

		// Reject if there is not enough tests remaining to pass consensus
		if (toolCount - countAlreadyAssessed < consensusN - countAlreadyAccepted) { return false; }
	}

	// For specificty modules
	if ((selectedModule == MODULE_SPECIFICITY) && (optimisation == "medium" || optimisation == "high"))
	{
		// Reject if the consensus was not passed
		if (stoi(candidateGuideResultMap["consensusCount"]) < consensusN) { return false; }
		// Reject if bowtie2 was not passed
		if (candidateGuideResultMap["passedBowtie"] == CODE_REJECTED) { return false; }
	}

	// Given none of the failure conditions have been meet, return true
	return true;
}

void printer(string_view formattedString)
{
	std::cout << formattedString << std::endl;
}

void errPrinter(string_view formattedString)
{
	std::cerr << formattedString << std::endl;
}

void runner(char* args)
{

	printer((std::format("| Calling: {}", args)));
	try 
	{
		int returnCode = system(args);
		if (returnCode != 0)
		{
			throw std::runtime_error("Runtime Error, returned a normal 0 value!");
		}
	}
	catch (std::exception e)
	{
		errPrinter(string_view(e.what()));
		return;
	}
    printer(string_view("| Finished"));
	return;
}
