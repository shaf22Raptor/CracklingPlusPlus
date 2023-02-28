#ifndef cas9InputModuleInclude
#define cas9InputMoudleinclude
#include <boost/regex.hpp>
#include "../include/inputModule.hpp"

class cas9InputModule : private inputModule
{
public:
	cas9InputModule(cracklingConfig config);
	void run();
	std::vector<guideResults>* next();
};


#endif // !cas9InputModuleInclude
