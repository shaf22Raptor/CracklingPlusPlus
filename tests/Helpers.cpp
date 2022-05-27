#include <Helpers.hpp>
#include <doctest.h>

TEST_CASE("Reverse Compliments are computed correctly")
{
    SUBCASE("Test rc on palindromic seq"){
        CHECK( rc("ACGTACGTACGTACGTACGT") == "ACGTACGTACGTACGTACGT");
    }
    SUBCASE("Test rc on 'AT' seq"){
        CHECK( rc("ATATATATATATATATATATAGG") == "CCTATATATATATATATATATAT");
    }
    SUBCASE("Test rc on 'A' seq"){
        CHECK( rc("AAAAAAAAAAAAAAAAAAAAAGG") == "CCTTTTTTTTTTTTTTTTTTTTT");
    }
    SUBCASE("Test rc on 'T' seq"){
        CHECK( rc("TTTTTTTTTTTTTTTTTTTTTGG") == "CCAAAAAAAAAAAAAAAAAAAAA");
    }
    SUBCASE("Test rc on 'G' seq"){
        CHECK( rc("GGGGGGGGGGGGGGGGGGGGGGG") == "CCCCCCCCCCCCCCCCCCCCCCC");
    }
    SUBCASE("Test rc on 'C' seq"){
        CHECK( rc("CCCCCCCCCCCCCCCCCCCCCGG") == "CCGGGGGGGGGGGGGGGGGGGGG");
    }
    SUBCASE("Test rc on empty seq"){
        CHECK_THROWS_AS(rc(""), const std::exception&);
    }
    SUBCASE("Test rc on short seq ( < 23)"){
        CHECK( rc("GACTGG") == "CCAGTC");
    }
    SUBCASE("Test rc on long seq ( > 23)"){
        CHECK( rc("GACTGGTTTGTGTCATATTCTTCCTGTGG") == "CCACAGGAAGAATATGACACAAACCAGTC");
    }
}
    