// Constants.hpp
#pragma once
#include <string>
#include <list>
#include <unordered_map>

const std::string CODE_ACCEPTED = "1";
const std::string CODE_REJECTED = "0";
const std::string CODE_UNTESTED = "?";
const std::string CODE_AMBIGUOUS = "-";
const std::string CODE_ERROR = "!";

const std::string MODULE_MM10DB = "mm10db";
const std::string MODULE_SGRNASCORER2 = "sgrnascorer2";
const std::string MODULE_CHOPCHOP = "chopchop";
const std::string MODULE_CONSENSUS = "consensus";
const std::string MODULE_SPECIFICITY = "specificity";

const std::unordered_map < std::string, std::string> DEFAULT_GUIDE_PROPERTIES = {
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

const std::list<std::string> DEFAULT_GUIDE_PROPERTIES_ORDER({
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
    //"passedReversePrimer"
    });