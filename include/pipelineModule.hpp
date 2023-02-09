#ifndef pipelineModuleInclude
#define pipelineModuleInclude

/*
The simplest form of a pipeline module in Crackling.
Simply define by it's ability to be `run`.
*/

class pipeLineModule
{
protected:
	run() = 0;
};


#endif // !pipelineModuleInclude