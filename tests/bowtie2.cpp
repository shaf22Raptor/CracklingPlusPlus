#include "../include/bowtie2.hpp"
#include "../include/doctest.h"


TEST_CASE("bowtie2 Module" * doctest::description("Ensure bowtie2 module is working correctly") * doctest::timeout(5))
{
    SUBCASE("Run bowtie2 and check results")
    {
        ConfigManager cm("data/test_config.ini");
        cm.set("general", "optimisation", "ultralow");
        bowtie2 testModule(cm);

        std::unordered_map<std::string, std::unordered_map<std::string, std::string>> result = {
            {"ACTCCTCATGCTGGACATTCTGG", {
                {"bowtieChr"    , CODE_UNTESTED},
                {"bowtieStart"  , CODE_UNTESTED},
                { "bowtieEnd"   , CODE_UNTESTED},
                {"passedBowtie" , CODE_UNTESTED}
            } },
            {"ATTCTGGTTCCTAGTATATCTGG", {
                {"bowtieChr"    , CODE_UNTESTED},
                {"bowtieStart"  , CODE_UNTESTED},
                { "bowtieEnd"   , CODE_UNTESTED},
                {"passedBowtie" , CODE_UNTESTED}
            } },
            {"GTATATCTGGAGAGTTAAGATGG", {
                {"bowtieChr"    , CODE_UNTESTED},
                {"bowtieStart"  , CODE_UNTESTED},
                { "bowtieEnd"   , CODE_UNTESTED},
                {"passedBowtie" , CODE_UNTESTED}
            }}
        };

        std::unordered_map<std::string, std::unordered_map<std::string, std::string>> expected = {
            {"ACTCCTCATGCTGGACATTCTGG", {
                {"bowtieChr"    , "*"},
                {"bowtieStart"  , "0"},
                { "bowtieEnd"   , "22"},
                {"passedBowtie" , CODE_ACCEPTED}
            } },
            {"ATTCTGGTTCCTAGTATATCTGG", {
                {"bowtieChr"    , "test-genome"},
                {"bowtieStart"  , "28"},
                { "bowtieEnd"   , "50"},
                {"passedBowtie" , CODE_REJECTED}
            } },
            {"GTATATCTGGAGAGTTAAGATGG", {
                {"bowtieChr"    , "*"},
                {"bowtieStart"  , "0"},
                { "bowtieEnd"   , "22"},
                {"passedBowtie" , CODE_ACCEPTED}
            }}
        };

        testModule.run(result);
        CHECK(result == expected);
    }
    SUBCASE("bowtie2 test clean up")
    {
        std::filesystem::remove_all("data/output");
    }
}
