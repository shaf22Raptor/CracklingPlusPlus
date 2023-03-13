#ifndef sgrnascorer2ModuleInclude
#define sgrnascorer2ModuleInclude

#include <string>
#include <array>
#include <unordered_map>
#include "../include/libsvm/svm.h"
#include "../include/util.hpp"
#include "../include/consensusModule.hpp"

class sgrnascorer2Module : private consensusModule
{
public:
	sgrnascorer2Module(const cracklingConfig& config);
	void run(std::vector<guideResults>& candidateGuides) final;
private:
	sgrnascorer2Config sgrnascorer2Config;
	bool processGuide(const guideResults& guide) final;
};

#endif // !sgrnascorer2ModuleInclude
