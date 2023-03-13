#ifndef chopchopModuleInclude
#define chopchopModuleInclude

#include <string>
#include "../include/util.hpp"
#include "../include/consensusModule.hpp"

class chopchopModule : private consensusModule
{
public:
	chopchopModule(const cracklingConfig& config);
	void run(std::vector<guideResults>& candidateGuides) final;
private:
	bool G20(std::string_view guide);
	bool processGuide(const guideResults& guide) final;
};

#endif // !chopchopModuleInclude
