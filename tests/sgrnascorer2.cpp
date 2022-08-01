#include "../include/sgrnascorer2.hpp"
#include "../include/doctest.h"

//TEST_CASE("sgRNAScorer2" * doctest::description("Ensure that the sgRNAScorer2 Module is working correctly") * doctest::timeout(5))
//{
//    ConfigManager cm("data/test_config.ini");
//    cm.set("general", "optimisation", "ultralow");
//    sgrnascorer2 testModule(cm);
//
//    std::map<std::string, std::map<std::string, std::string>> result = {
//        {"ACTCCTCATGCTGGACATTCTGG", {
//            {"acceptedBySgRnaScorer"    , CODE_UNTESTED},
//            {"sgrnascorer2score"        , CODE_UNTESTED}
//        }},
//        {"ATTCTGGTTCCTAGTATATCTGG", {
//            {"acceptedBySgRnaScorer"    , CODE_UNTESTED},
//            {"sgrnascorer2score"        , CODE_UNTESTED}
//        }},
//        {"GTATATCTGGAGAGTTAAGATGG", {
//            {"acceptedBySgRnaScorer"    , CODE_UNTESTED},
//            {"sgrnascorer2score"        , CODE_UNTESTED}
//        }}
//    };
//
//    std::map<std::string, std::map<std::string, std::string>> expected = {
//        {"ACTCCTCATGCTGGACATTCTGG", {
//            {"acceptedBySgRnaScorer"    , CODE_REJECTED},
//            {"sgrnascorer2score"        , "-1.8394958342445178"}
//        }},
//        {"ATTCTGGTTCCTAGTATATCTGG", {
//            {"acceptedBySgRnaScorer"    , CODE_REJECTED},
//            {"sgrnascorer2score"        , "-3.3132386532897384"}
//        }},
//        {"GTATATCTGGAGAGTTAAGATGG", {
//            {"acceptedBySgRnaScorer"    , CODE_REJECTED},
//            {"sgrnascorer2score"        , "-0.328729309372379"}
//        }}
//    };
//
//    testModule.run(result);
//
//    CHECK(result == expected);
//}
