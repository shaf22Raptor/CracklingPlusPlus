#ifndef consensusModuleInclude
#define consensusModuleInclude
#include "../include/util.hpp"
#include <vector>

class consensusModule
{
protected:
	consensusModule(cracklingConfig config);
	run(const std::vector<guideResults>& candidateGuides) = 0;
private:
	run() {};
};

consensusModule::consensusModule(cracklingConfig config)
{

}

#endif // !consensusModuleInclude
