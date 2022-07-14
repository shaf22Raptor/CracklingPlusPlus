// Helpers.hpp
#pragma once
#include <iostream>
#include <string>
#include <algorithm>
#include <iterator>
#include <array>
#include <ctime>
#include <map>
#include <vector>
#include <format>
#include <locale>
#include <iomanip>
#include <sstream>
 
#include <Constants.hpp>

std::string makeUpper(const std::string& s);

std::string makeLower(const std::string& s);

std::vector<std::string> split(std::string& s, std::string_view delimiter);

std::string rtrim(std::string_view s);

std::string ltrim(std::string_view s);

std::string trim(const std::string& s);

std::string rc(std::string DNA);

bool filterCandidateGuides(std::map<std::string, std::string, std::less<>> candidateGuideResultMap, std::string_view selectedModule, std::string_view optimisation, const int& consensusN, const int& toolCount);

void printer(std::string_view formattedString);

void errPrinter(std::string_view formattedString);

void runner(const char* args);

class ReturnCode : public std::logic_error
{
public:
	ReturnCode() : std::logic_error("The externally called program returned a non-zero value") { };
};

std::string commaify(int);

class comma_numpunct : public std::numpunct<char>
{
protected:
    char do_thousands_sep() const final
    {
        return ',';
    }

    std::string do_grouping() const final
    {
        return "\03";
    }
};
