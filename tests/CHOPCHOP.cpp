#include "../include/CHOPCHOP.hpp"
#include "../include/doctest.h"

TEST_CASE("G20" * doctest::description("Ensure that the G20 function is working correctly") * doctest::timeout(5))
{
    SUBCASE("Test G20 on 'AT' seq")
    {
        CHECK(CHOPCHOP::G20("ATATATATATATATATATATA") == false);
    }
    SUBCASE("Test G20 on 'A' seq")
    {
        CHECK(CHOPCHOP::G20("AAAAAAAAAAAAAAAAAAAAA") == false);
    }
    SUBCASE("Test G20 on 'T' seq")
    {
        CHECK(CHOPCHOP::G20("TTTTTTTTTTTTTTTTTTTTT") == false);
    }
    SUBCASE("Test G20 on 'G' seq")
    {
        CHECK(CHOPCHOP::G20("GGGGGGGGGGGGGGGGGGGGG") == true);
    }
    SUBCASE("Test G20 on 'C' seq")
    {
        CHECK(CHOPCHOP::G20("CCCCCCCCCCCCCCCCCCCCC") == false);
    }
    SUBCASE("Test G20 on repeating 'ATGC' seq")
    {
        CHECK(CHOPCHOP::G20("ATGCATGCATGCATGCATGCA") == false);
    }
    SUBCASE("Test G20 on random seq")
    {
        CHECK(CHOPCHOP::G20("TTTGTGTCATATTCTTCCTGT") == true);
    }
    SUBCASE("Test G20 on mixed uppercase and lowercase seq with no G @ pos 20")
    {
        CHECK(CHOPCHOP::G20("AtGgtcatGAactgcaAGAtc") == false);
    }
    SUBCASE("Test G20 on mixed uppercase and lowercase seq with G @ pos 20")
    {
        CHECK(CHOPCHOP::G20("AtGgtGAactcgcaAGatAgc") == false);
    }
    SUBCASE("Test G20 on short seq (length < 20)")
    {
        CHECK_THROWS_WITH_AS(CHOPCHOP::G20("TTTGTGTCAT"), "CHOPCHOP G20: Input lenght must be >= 20!", G20Input);
    }
    SUBCASE("Test G20 on long seq (length > 20) with G @ pos 20")
    {
        CHECK(CHOPCHOP::G20("TTTGTGTCATATTCTTCCTGTTTTGTGTCATATTCTTCCTGTTTTGTGTCATATTCTTCCTGT") == true);
    }
    SUBCASE("Test G20 on long seq (length > 20) with NO G @ pos 20")
    {
        CHECK(CHOPCHOP::G20("ATGCATGCATGCATGCATGCAATGCATGCATGCATGCATGCAATGCATGCATGCATGCATGCA") == false);
    }
}

TEST_CASE("CHOPCHOP Module" * doctest::description("Ensure CHOPCHOP module is working correctly") * doctest::timeout(5))
{
    ConfigManager cm("data/test_config.ini");
    cm.set("general", "optimisation", "ultralow");
    CHOPCHOP testModule(cm);
    
    std::unordered_map<std::string, std::unordered_map<std::string, std::string>> result = {
        {"ACTCCTCATGCTGGACATTCTGG", {{"passedG20" , CODE_UNTESTED}} },
        {"ATTCTGGTTCCTAGTATATCTGG", {{"passedG20" , CODE_UNTESTED}} },
        {"GTATATCTGGAGAGTTAAGATGG", {{"passedG20" , CODE_UNTESTED}} }
    };

    std::unordered_map<std::string, std::unordered_map<std::string, std::string>> expected = {
    {"ACTCCTCATGCTGGACATTCTGG", {{"passedG20" , CODE_REJECTED}} },
    {"ATTCTGGTTCCTAGTATATCTGG", {{"passedG20" , CODE_REJECTED}} },
    {"GTATATCTGGAGAGTTAAGATGG", {{"passedG20" , CODE_REJECTED}} }
    };

    testModule.run(result);
    CHECK(result == expected);
}