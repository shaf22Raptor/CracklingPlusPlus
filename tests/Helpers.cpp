#include "../include/Helpers.hpp"
#include "../include/doctest.h"

TEST_CASE("rc" * doctest::description("Ensure that the Reverse Compliment function is working correctly") * doctest::timeout(5))
{
    SUBCASE("Test rc on palindromic seq")
    {
        CHECK(rc("ACGTACGTACGTACGTACGT") == "ACGTACGTACGTACGTACGT");
    }
    SUBCASE("Test rc on 'AT' seq")
    {
        CHECK(rc("ATATATATATATATATATATAGG") == "CCTATATATATATATATATATAT");
    }
    SUBCASE("Test rc on 'A' seq")
    {
        CHECK(rc("AAAAAAAAAAAAAAAAAAAAAGG") == "CCTTTTTTTTTTTTTTTTTTTTT");
    }
    SUBCASE("Test rc on 'T' seq")
    {
        CHECK(rc("TTTTTTTTTTTTTTTTTTTTTGG") == "CCAAAAAAAAAAAAAAAAAAAAA");
    }
    SUBCASE("Test rc on 'G' seq")
    {
        CHECK(rc("GGGGGGGGGGGGGGGGGGGGGGG") == "CCCCCCCCCCCCCCCCCCCCCCC");
    }
    SUBCASE("Test rc on 'C' seq")
    {
        CHECK(rc("CCCCCCCCCCCCCCCCCCCCCGG") == "CCGGGGGGGGGGGGGGGGGGGGG");
    }
    SUBCASE("Test rc on short seq ( < 23)")
    {
        CHECK(rc("GACTGG") == "CCAGTC");
    }
    SUBCASE("Test rc on long seq ( > 23)")
    {
        CHECK(rc("GACTGGTTTGTGTCATATTCTTCCTGTGG") == "CCACAGGAAGAATATGACACAAACCAGTC");
    }
    SUBCASE("Test rc on empty seq")
    {
        CHECK_THROWS_WITH_AS(rc(""), "Type Error, Seqeunce length must be greater than 0!", std::length_error);
    }
    SUBCASE("Test rc on extremely large seq")
    {
        std::string extremelyLargeInput = "AAATCTAGCTGCTACGATCGATCGATCGCTAGCATGATCGATCGATCGATCGATCGACTGACTAGCTAGCATGCTAGCATCGATCGATGCTAGCATGCATGCTAGCATGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTGATCGATGACTGCTAGTACGTCGTACGTAGTCGTAGTAGCTGATCGATGCTACGTAGCATGCTAGCATGCATAGCTAGCTAGCTAGCTAGCTGCATAGCTAGCTAGCTAGCTGATCGATCGATCGATCGTCGCTAGAATCTAGCTGCTACGATCGATCGATCGCTAGCATGATCGATCGATCGATCGATCGACTGACTAGCTAGCATGCTAGCATCGATCGATGCTAGCATGCATGCTAGCATGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTGATCGATGACTGCTAGTACGTCGTACGTAGTCGTAGTAGCTGATCGATGCTACGTAGCATGCTAGCATGCATAGCTAGCTAGCTAGCTAGCTGCATAGCTAGCTAGCTAGCTGATCGATCGATCGATCGTCGCTAGAATCTAGCTGCTACGATCGATCGATCGCTAGCATGATCGATCGATCGATCGATCGACTGACTAGCTAGCATGCTAGCATCGATCGATGCTAGCATGCATGCTAGCATGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTGATCGATGACTGCTAGTACGTCGTACGTAGTCGTAGTAGCTGATCGATGCTACGTAGCATGCTAGCATGCATAGCTAGCTAGCTAGCTAGCTGCATAGCTAGCTAGCTAGCTGATCGATCGATCGATCGTCGCTAGAATCTAGCTGCTACGATCGATCGATCGCTAGCATGATCGATCGATCGATCGATCGACTGACTAGCTAGCATGCTAGCATCGATCGATGCTAGCATGCATGCTAGCATGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTGATCGATGACTGCTAGTACGTCGTACGTAGTCGTAGTAGCTGATCGAT";
        CHECK_THROWS_WITH_AS(rc(extremelyLargeInput), "Type Error, Seqeunce length must be less than 1024!", std::length_error);
    }
}


TEST_CASE("commaify" * doctest::description("Ensure that the commaify function is working correctly") * doctest::timeout(5))
{
    SUBCASE("Test commaify on 0")
    {
        CHECK(commaify(0) == "0");
    }

    SUBCASE("Test commaify on -1")
    {
        CHECK(commaify(-1) == "-1");
    }

    SUBCASE("Test commaify on number < 1,000")
    {
        CHECK(commaify(999) == "999");
    }

    SUBCASE("Test commaify on number >= 1,000")
    {
        CHECK(commaify(1000) == "1,000");
    }

    SUBCASE("Test commaify on number > -1,000")
    {
        CHECK(commaify(-999) == "-999");
    }

    SUBCASE("Test commaify on number <= -1,000")
    {
        CHECK(commaify(-1000) == "-1,000");
    }

    SUBCASE("Test commaify on largest signed long long")
    {
        CHECK(commaify(9223372036854775807) == "9,223,372,036,854,775,807");
    }

    SUBCASE("Test commaify on smallest signed long long")
    {
        CHECK(commaify(-9223372036854775807) == "-9,223,372,036,854,775,807");
    }

    SUBCASE("Test commaify on largest signed int")
    {
        CHECK(commaify(2147483647) == "2,147,483,647");
    }

    SUBCASE("Test commaify on smallest signed int")
    {
        CHECK(commaify(-2147483647) == "-2,147,483,647");
    }
}