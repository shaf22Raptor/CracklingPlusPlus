#ifndef specificityModuleInclude
#define specificityModuleInclude
#include <vector>
#include "../include/util.hpp"
#include "../include/pipelineModule.hpp"

class specificityModule : private pipelineModule
{
protected:
	bool toolIsSelected;
	optimisationLevel optimsationLevel;
	uint8_t consensusN;
	specificityModule(const cracklingConfig& config);
	virtual void run(std::vector<guideResults>& candidateGuides) = 0;
	virtual bool processGuide(const guideResults& guide) = 0;
private:
	void run() final;
};

#endif // !specificityModuleInclude
