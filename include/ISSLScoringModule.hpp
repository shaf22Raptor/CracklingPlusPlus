#ifndef ISSLScoringModuleInclude
#define ISSLScoringModuleInclude

#include "../include/specificityModule.hpp"
#include "../include/util.hpp"

class ISSLScoringModule : private specificityModule
{
public:
	ISSLScoringModule(cracklingConfig config);
	void run(std::vector<guideResults>& candidateGuides) final;
private:
	bool processGuide(const guideResults& guide) final;
};


#endif // !ISSLScoringModuleInclude
