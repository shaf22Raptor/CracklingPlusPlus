#ifndef cas9InputModuleInclude
#define cas9InputMoudleinclude
#include <boost/regex.hpp>
#include "../include/inputModule.hpp"

class cas9InputModule : private inputModule
{
public:
	cas9InputModule(cracklingConfig config);
};


#endif // !cas9InputModuleInclude
