#ifndef chopchopModuleInclude
#define chopchopModuleInclude

#include <iostream>
#include <string>
#include <string_view>
#include "../include/util.hpp"
#include "../include/consensusModule.hpp"
#define FMT_HEADER_ONLY
#include "../include/fmt/format.h"

class chopchopModule : private consensusModule
{
public:
	chopchopModule(cracklingConfig config);
	void run(std::vector<guideResults>& candidateGuides) final;
private:
	bool G20(std::string_view guide);
	bool processGuide(const guideResults& guide) final;
};

#endif // !chopchopModuleInclude
