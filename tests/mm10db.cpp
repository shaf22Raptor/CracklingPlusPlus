#include "../include/mm10db.hpp"
#include "../include/doctest.h"
#include <fstream>
#include <vector>

TEST_CASE("LeadingT" * doctest::description("Ensure that the LeadingT function is working correctly") * doctest::timeout(5))
{
    SUBCASE("Test LeadingT on 'AT' seq")
    {
        CHECK(mm10db::leadingT("ATATATATATATATATATATAGG") == false);
    }
    SUBCASE("Test LeadingT on 'A' seq")
    {
        CHECK(mm10db::leadingT("AAAAAAAAAAAAAAAAAAAAAGG") == false);
    }
    SUBCASE("Test LeadingT on 'T' seq")
    {
        CHECK(mm10db::leadingT("TTTTTTTTTTTTTTTTTTTTTGG") == true);
    }
    SUBCASE("Test LeadingT on 'G' seq")
    {
        CHECK(mm10db::leadingT("GGGGGGGGGGGGGGGGGGGGGGG") == false);
    }
    SUBCASE("Test LeadingT on 'C' seq")
    {
        CHECK(mm10db::leadingT("CCCCCCCCCCCCCCCCCCCCCGG") == false);
    }
    SUBCASE("Test LeadingT on repeating 'ATGC' seq")
    {
        CHECK(mm10db::leadingT("ATGCATGCATGCATGCATGCAGG") == false);
    }
    SUBCASE("Test LeadingT on random seq")
    {
        CHECK(mm10db::leadingT("TTTGTGTCATATTCTTCCTGTGG") == true);
    }
    SUBCASE("Test LeadingT on mixed uppercase and lowercase seq with NO leading T")
    {
        CHECK(mm10db::leadingT("AtGgtcatGAactgcaAGAtcGG") == false);
    }
    SUBCASE("Test LeadingT on mixed uppercase and lowercase seq with leading T")
    {
        CHECK(mm10db::leadingT("ttGgtGAactcgcaAGatAgcgg") == false);
    }
}

TEST_CASE("AT_percentage" * doctest::description("Ensure that the AT_percentage function is working correctly") * doctest::timeout(5))
{
    SUBCASE("Test AT_percentage on 'AT' seq")
    {
        CHECK(mm10db::AT_percentage("ATATATATATATATATATATA") == 100.00f);
    }
    SUBCASE("Test AT_percentage on 'A' seq")
    {
        CHECK(mm10db::AT_percentage("AAAAAAAAAAAAAAAAAAAAA") == 100.00f);
    }
    SUBCASE("Test AT_percentage on 'T' seq")
    {
        CHECK(mm10db::AT_percentage("TTTTTTTTTTTTTTTTTTTTT") == 100.00f);
    }
    SUBCASE("Test AT_percentage on 'G' seq")
    {
        CHECK(mm10db::AT_percentage("GGGGGGGGGGGGGGGGGGGGG") == 0.00f);
    }
    SUBCASE("Test AT_percentage on 'C' seq")
    {
        CHECK(mm10db::AT_percentage("CCCCCCCCCCCCCCCCCCCCC") == 0.00f);
    }
    SUBCASE("Test AT_percentage on repeating 'ATGC' seq")
    {
        CHECK(mm10db::AT_percentage("ATGCATGCATGCATGCATGC") == 50.00f);
    }
    SUBCASE("Test AT_percentage on random seq")
    {
        CHECK(mm10db::AT_percentage("TTTGTGTCATATTCTTCCTA") == 70.00f);
    }
    SUBCASE("Test AT_percentage on mixed uppercase and lowercase seq")
    {
        CHECK(mm10db::AT_percentage("AtgcaTGcatGcatGcatGC") == 10.00f);
    }
    SUBCASE("Test AT_percentage on mixed uppercase and lowercase seq with leading T")
    {
        CHECK(mm10db::AT_percentage("tttGtGtcAtAttCttCctA") == 15.00f);
    }
}

TEST_CASE("polyT" * doctest::description("Ensure that the polyT function is working correctly") * doctest::timeout(5))
{
    SUBCASE("Test polyT on 'AT' seq")
    {
        CHECK(mm10db::polyT("ATATATATATATATATATATA") == false);
    }
    SUBCASE("Test polyT on 'A' seq")
    {
        CHECK(mm10db::polyT("AAAAAAAAAAAAAAAAAAAAA") == false);
    }
    SUBCASE("Test polyT on 'T' seq")
    {
        CHECK(mm10db::polyT("TTTTTTTTTTTTTTTTTTTTT") == true);
    }
    SUBCASE("Test polyT on 'G' seq")
    {
        CHECK(mm10db::polyT("GGGGGGGGGGGGGGGGGGGGG") == false);
    }
    SUBCASE("Test polyT on 'C' seq")
    {
        CHECK(mm10db::polyT("CCCCCCCCCCCCCCCCCCCCC") == false);
    }
    SUBCASE("Test polyT on repeating 'ATGC' seq")
    {
        CHECK(mm10db::polyT("ATGCATGCATGCATGCATGC") == false);
    }
    SUBCASE("Test polyT on random seq")
    {
        CHECK(mm10db::polyT("TTTGTGTCATATTCTTCCTA") == false);
    }
    SUBCASE("Test polyT on mixed uppercase and lowercase seq")
    {
        CHECK(mm10db::polyT("AtgcaTGcatGcatGcatGC") == false);
    }
    SUBCASE("Test polyT on mixed uppercase and lowercase seq with leading T")
    {
        CHECK(mm10db::polyT("tttGtGtcAtAttCttCctA") == false);
    }
}

TEST_CASE("transToDNA" * doctest::description("Ensure that the transToDNA function is working correctly") * doctest::timeout(5))
{
    SUBCASE("Test transToDNA on 'AU' seq")
    {
        CHECK(mm10db::transToDNA("AUAUAUAUAUAUAUAUAUAUA") == "ATATATATATATATATATATA");
    }
    SUBCASE("Test transToDNA on 'A' seq")
    {
        CHECK(mm10db::transToDNA("AAAAAAAAAAAAAAAAAAAAA") == "AAAAAAAAAAAAAAAAAAAAA");
    }
    SUBCASE("Test transToDNA on 'T' seq")
    {
        CHECK(mm10db::transToDNA("UUUUUUUUUUUUUUUUUUUUU") == "TTTTTTTTTTTTTTTTTTTTT");
    }
    SUBCASE("Test transToDNA on 'G' seq")
    {
        CHECK(mm10db::transToDNA("GGGGGGGGGGGGGGGGGGGGG") == "GGGGGGGGGGGGGGGGGGGGG");
    }
    SUBCASE("Test transToDNA on 'C' seq")
    {
        CHECK(mm10db::transToDNA("CCCCCCCCCCCCCCCCCCCCC") == "CCCCCCCCCCCCCCCCCCCCC");
    }
    SUBCASE("Test transToDNA on repeating 'ATGC' seq")
    {
        CHECK(mm10db::transToDNA("AUGCAUGCAUGCAUGCAUGCA") == "ATGCATGCATGCATGCATGCA");
    }
    SUBCASE("Test transToDNA on random seq")
    {
        CHECK(mm10db::transToDNA("UUUGUGUCAUAUUCUUCCUAU") == "TTTGTGTCATATTCTTCCTAT");
    }
    SUBCASE("Test transToDNA on mixed uppercase and lowercase seq")
    {
        CHECK(mm10db::transToDNA("AugcaUGcauGcauGcauGC") == "AugcaTGcauGcauGcauGC");
    }
}

//TEST_CASE("RNAFold" * doctest::description("Ensure that RNAFold is working correctly") * doctest::timeout(5))
//{
//    std::vector<std::string> expected{ {
//        "GCUCCUCAUGCUGGACAUUCGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGCUUUU",
//        "((.......))(((((..(((((((((.((((....))))...)))))))..))...)))))......((((....))))(((((((...)))))))... (-22.40)",
//        "GUUCUGGUUCCUAGUAUAUCGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGCUUUU",
//        ".(((((((((..((........))..)))))))))....((((((...((((.(......).)))).)))))).......(((((((...)))))))... (-25.20)",
//        "GUAUAUCUGGAGAGUUAAGAGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGCUUUU",
//        ".......((((.((((....(((((((.((((....))))...)))))))..))))..))))......((((....))))(((((((...)))))))... (-22.50)"
//    } };
//
//    std::vector<std::string> result;
//
//    std::ofstream outFile("input.txt", std::ios::binary);
//
//    outFile << "GCTCCTCATGCTGGACATTCGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGCUUUU\n"
//            << "GTTCTGGTTCCTAGTATATCGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGCUUUU\n"
//            << "GTATATCTGGAGAGTTAAGAGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGCUUUU\n";
//
//    outFile.close();
//
//    system("RNAfold --noPS -j16 -i input.txt > output.txt");
//
//    std::ifstream inFile("output.txt", std::ios::binary);
//
//    for (std::string line; std::getline(inFile, line);)
//    {
//        result.push_back(line);
//    }
//
//    inFile.close();
//
//    remove("input.txt");
//    remove("output.txt");
//
//    CHECK(result == expected);
//}
//
//TEST_CASE("mm10db" * doctest::description("Ensure that mm10db module is working correctly") * doctest::timeout(5))
//{
//    ConfigManager cm("data/test_config.ini");
//    cm.set("general", "optimisation", "ultralow");
//    mm10db testModule(cm);
//    
//    std::map<std::string, std::map<std::string, std::string>> result = {
//        {"ACTCCTCATGCTGGACATTCTGG", {
//            {"passedAvoidLeadingT" , CODE_UNTESTED},
//            {"passedATPercent" , CODE_UNTESTED},
//            {"AT" , CODE_UNTESTED},
//            {"passedTTTT" , CODE_UNTESTED},
//            {"ssL1" , CODE_UNTESTED},
//            {"ssStructure" , CODE_UNTESTED},
//            {"ssEnergy" , CODE_UNTESTED},
//            {"passedSecondaryStructure" , CODE_UNTESTED},
//            {"acceptedByMm10db" , CODE_UNTESTED}
//        }},
//        {"ATTCTGGTTCCTAGTATATCTGG", {
//            {"passedAvoidLeadingT" , CODE_UNTESTED},
//            {"passedATPercent" , CODE_UNTESTED},
//            {"AT" , CODE_UNTESTED},
//            {"passedTTTT" , CODE_UNTESTED},
//            {"ssL1" , CODE_UNTESTED},
//            {"ssStructure" , CODE_UNTESTED},
//            {"ssEnergy" , CODE_UNTESTED},
//            {"passedSecondaryStructure" , CODE_UNTESTED},
//            {"acceptedByMm10db" , CODE_UNTESTED}
//        }},
//        {"GTATATCTGGAGAGTTAAGATGG", {
//            {"passedAvoidLeadingT" , CODE_UNTESTED},
//            {"passedATPercent" , CODE_UNTESTED},
//            {"AT" , CODE_UNTESTED},
//            {"passedTTTT" , CODE_UNTESTED},
//            {"ssL1" , CODE_UNTESTED},
//            {"ssStructure" , CODE_UNTESTED},
//            {"ssEnergy" , CODE_UNTESTED},
//            {"passedSecondaryStructure" , CODE_UNTESTED},
//            {"acceptedByMm10db" , CODE_UNTESTED}
//        }}
//    };
//
//    std::map<std::string, std::map<std::string, std::string>> expected = {
//        {"ACTCCTCATGCTGGACATTCTGG", {
//            {"passedAvoidLeadingT" , CODE_ACCEPTED},
//            {"passedATPercent" , CODE_ACCEPTED},
//            {"AT" , "50.0"},
//            {"passedTTTT" , CODE_ACCEPTED},
//            {"ssL1" , "GCUCCUCAUGCUGGACAUUCGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGCUUUU"},
//            {"ssStructure" , "((.......))(((((..(((((((((.((((....))))...)))))))..))...)))))......((((....))))(((((((...)))))))..."},
//            {"ssEnergy" , "-22.40"},
//            {"passedSecondaryStructure" , CODE_ACCEPTED},
//            {"acceptedByMm10db" , CODE_ACCEPTED}
//        }},
//        {"ATTCTGGTTCCTAGTATATCTGG", {
//            {"passedAvoidLeadingT" , CODE_ACCEPTED},
//            {"passedATPercent" , CODE_ACCEPTED},
//            {"AT" , "65.0"},
//            {"passedTTTT" , CODE_ACCEPTED},
//            {"ssL1" , "GUUCUGGUUCCUAGUAUAUCGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGCUUUU"},
//            {"ssStructure" , ".(((((((((..((........))..)))))))))....((((((...((((.(......).)))).)))))).......(((((((...)))))))..."},
//            {"ssEnergy" , "-25.20"},
//            {"passedSecondaryStructure" , CODE_REJECTED},
//            {"acceptedByMm10db" , CODE_REJECTED}
//        }},
//        {"GTATATCTGGAGAGTTAAGATGG", {
//            {"passedAvoidLeadingT" , CODE_ACCEPTED},
//            {"passedATPercent" , CODE_ACCEPTED},
//            {"AT" , "65.0"},
//            {"passedTTTT" , CODE_ACCEPTED},
//            {"ssL1" , "GUAUAUCUGGAGAGUUAAGAGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGCUUUU"},
//            {"ssStructure" , ".......((((.((((....(((((((.((((....))))...)))))))..))))..))))......((((....))))(((((((...)))))))..."},
//            {"ssEnergy" , "-22.50"},
//            {"passedSecondaryStructure" , CODE_ACCEPTED},
//            {"acceptedByMm10db" , CODE_ACCEPTED}
//        }}
//    };
//
//    testModule.run(result);
//    CHECK(result == expected);
//}