#include <CHOPCHOP.hpp>
#include <doctest.h>

TEST_CASE("G20" * doctest::description("Ensure that the G20 function is working correctly") * doctest::timeout(5))
{
    SUBCASE("Test G20 'AT' seq")
    {
        CHECK(CHOPCHOP::G20("ATATATATATATATATATATA") == false);
    }
    SUBCASE("Test G20 'A' seq")
    {
        CHECK(CHOPCHOP::G20("AAAAAAAAAAAAAAAAAAAAA") == false);
    }
    SUBCASE("Test G20 'T' seq")
    {
        CHECK(CHOPCHOP::G20("TTTTTTTTTTTTTTTTTTTTT") == false);
    }
    SUBCASE("Test G20 'G' seq")
    {
        CHECK(CHOPCHOP::G20("GGGGGGGGGGGGGGGGGGGGG") == true);
    }
    SUBCASE("Test G20 'C' seq")
    {
        CHECK(CHOPCHOP::G20("CCCCCCCCCCCCCCCCCCCCC") == false);
    }
    SUBCASE("Test G20 repeating 'ATGC' seq")
    {
        CHECK(CHOPCHOP::G20("ATGCATGCATGCATGCATGCA") == false);
    }
    SUBCASE("Test G20 random seq")
    {
        CHECK(CHOPCHOP::G20("TTTGTGTCATATTCTTCCTGT") == true);
    }
    SUBCASE("Test G20 mixed uppercase and lowercase seq with no G @ pos 20")
    {
        CHECK(CHOPCHOP::G20("AtGgtcatGAactgcaAGAtc") == false);
    }
    SUBCASE("Test G20 mixed uppercase and lowercase seq with G @ pos 20")
    {
        CHECK(CHOPCHOP::G20("AtGgtGAactcgcaAGatAgc") == true);
    }
}

//TEST_CASE("CHOPCHOP Module" * doctest::description("Ensure CHOPCHOP module is working correctly") * doctest::timeout(5))
//{
//    ConfigManager cm("data/test_config.ini");
//    cm.set("general", "optimisation", "ultralow");
//    CHOPCHOP testModule(cm);
//    
//    std::map<std::string, std::map<std::string, std::string>> result = {
//        {"ACTCCTCATGCTGGACATTCTGG", {{"passedG20" , CODE_UNTESTED}} },
//        {"ATTCTGGTTCCTAGTATATCTGG", {{"passedG20" , CODE_UNTESTED}} },
//        {"GTATATCTGGAGAGTTAAGATGG", {{"passedG20" , CODE_UNTESTED}} }
//    };
//
//    std::map<std::string, std::map<std::string, std::string>> expected = {
//    {"ACTCCTCATGCTGGACATTCTGG", {{"passedG20" , CODE_REJECTED}} },
//    {"ATTCTGGTTCCTAGTATATCTGG", {{"passedG20" , CODE_REJECTED}} },
//    {"GTATATCTGGAGAGTTAAGATGG", {{"passedG20" , CODE_REJECTED}} }
//    };
//
//    testModule.run(result);
//    CHECK(result == expected);
//}