#include <Constants.hpp>
#include <doctest.h>

TEST_CASE("Constants" * doctest::description("Ensure that the Constant values have not changed") * doctest::timeout(5))
{
    SUBCASE("Test CODE_ACCEPTED")
    {
        CHECK(CODE_ACCEPTED == "1");
    }
    SUBCASE("Test CODE_REJECTED")
    {
        CHECK(CODE_REJECTED == "0");
    }
    SUBCASE("Test CODE_AMBIGUOUS")
    {
        CHECK(CODE_AMBIGUOUS == "-");
    }
    SUBCASE("Test CODE_UNTESTED")
    {
        CHECK(CODE_UNTESTED == "?");
    }
    SUBCASE("Test MODULE_CHOPCHOP")
    {
        CHECK(MODULE_CHOPCHOP == "chopchop");
    }
    SUBCASE("Test MODULE_MM10DB")
    {
        CHECK(MODULE_MM10DB == "mm10db");
    }
    SUBCASE("Test MODULE_SGRNASCORER2")
    {
        CHECK(MODULE_SGRNASCORER2 == "sgrnascorer2");
    }
    SUBCASE("Test MODULE_CONSENSUS")
    {
        CHECK(MODULE_CONSENSUS == "consensus");
    }
    SUBCASE("Test MODULE_SPECIFICITY")
    {
        CHECK(MODULE_SPECIFICITY == "specificity");
    }
    SUBCASE("Test DEFAULT_GUIDE_PROPERTIES")
    {
        std::map < std::string, std::string> result = {
            {"seq"                      , ""},
            {"header"                   , ""},
            {"isUnique"                 , CODE_ACCEPTED},
            {"start"                    , CODE_UNTESTED},
            {"end"                      , CODE_UNTESTED},
            {"strand"                   , CODE_UNTESTED},
            {"passedTTTT"               , CODE_UNTESTED},
            {"passedATPercent"          , CODE_UNTESTED},
            {"passedG20"                , CODE_UNTESTED},
            {"passedSecondaryStructure" , CODE_UNTESTED},
            {"ssL1"                     , CODE_UNTESTED},
            {"ssStructure"              , CODE_UNTESTED},
            {"ssEnergy"                 , CODE_UNTESTED},
            {"acceptedByMm10db"         , CODE_UNTESTED},
            {"acceptedBySgRnaScorer"    , CODE_UNTESTED},
            {"consensusCount"           , CODE_UNTESTED},
            {"passedBowtie"             , CODE_UNTESTED},
            {"passedOffTargetScore"     , CODE_UNTESTED},
            {"sgrnascorer2score"        , CODE_UNTESTED},
            {"AT"                       , CODE_UNTESTED},
            {"bowtieChr"                , CODE_UNTESTED},
            {"bowtieStart"              , CODE_UNTESTED},
            {"bowtieEnd"                , CODE_UNTESTED},
            {"mitOfftargetscore"        , CODE_UNTESTED},
            {"cfdOfftargetscore"        , CODE_UNTESTED},
            {"passedAvoidLeadingT"      , CODE_UNTESTED},
        };
        CHECK(DEFAULT_GUIDE_PROPERTIES == result);
    }
    SUBCASE("Test DEFAULT_GUIDE_PROPERTIES_ORDER")
    {
        std::list<std::string> result({
            "seq",
            "sgrnascorer2score",
            "header",
            "start",
            "end",
            "strand",
            "isUnique",
            "passedG20",
            "passedTTTT",
            "passedATPercent",
            "passedSecondaryStructure",
            "ssL1",
            "ssStructure",
            "ssEnergy",
            "acceptedByMm10db",
            "acceptedBySgRnaScorer",
            "consensusCount",
            "passedBowtie",
            "passedOffTargetScore",
            "AT",
            "bowtieChr",
            "bowtieStart",
            "bowtieEnd",
            "mitOfftargetscore",
            "cfdOfftargetscore",
            "passedAvoidLeadingT",
        });
        CHECK(DEFAULT_GUIDE_PROPERTIES_ORDER == result);
    }

}