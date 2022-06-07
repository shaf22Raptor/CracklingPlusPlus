#include <CHOPCHOP.hpp>
#include <doctest.h>

TEST_CASE("G20" * doctest::description("Ensure that the G20 function is working correctly") * doctest::timeout(5))
{
    ConfigManager cm("data/test_config.ini");
    CHOPCHOP testModule(cm);

    SUBCASE("Test G20 'AT' seq")
    {
        CHECK(testModule.G20("ATATATATATATATATATATA") == false);
    }
    SUBCASE("Test G20 'A' seq")
    {
        CHECK(testModule.G20("AAAAAAAAAAAAAAAAAAAAA") == false);
    }
    SUBCASE("Test G20 'T' seq")
    {
        CHECK(testModule.G20("TTTTTTTTTTTTTTTTTTTTT") == false);
    }
    SUBCASE("Test G20 'G' seq")
    {
        CHECK(testModule.G20("GGGGGGGGGGGGGGGGGGGGG") == true);
    }
    SUBCASE("Test G20 'C' seq")
    {
        CHECK(testModule.G20("CCCCCCCCCCCCCCCCCCCCC") == false);
    }
    SUBCASE("Test G20 repeating 'ATGC' seq")
    {
        CHECK(testModule.G20("ATGCATGCATGCATGCATGCA") == false);
    }
    SUBCASE("Test G20 random seq")
    {
        CHECK(testModule.G20("TTTGTGTCATATTCTTCCTGT") == true);
    }
    SUBCASE("Test G20 mixed uppercase and lowercase seq with no G @ pos 20")
    {
        CHECK(testModule.G20("AtGgtcatGAactgcaAGAtc") == false);
    }
    SUBCASE("Test G20 mixed uppercase and lowercase seq with G @ pos 20")
    {
        CHECK(testModule.G20("AtGgtGAactcgcaAGatAgc") == true);
    }
}