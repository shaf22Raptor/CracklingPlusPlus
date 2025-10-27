#pragma once
#include <string>

struct AppOptions {
    bool        useCLSeq2Sig = false;
    bool        useCLScoring = false;

    int         clPlatformIndex = -1;      // prefer indices >
    int         clDeviceIndex   = -1;      // >
    std::string clDeviceGfxId;             // > gfx#### >
    std::string clDeviceNameContains;      // > name substring

    std::string traceFile = "profile_trace.txt";
};