#ifndef bowtie2ModuleInclude
#define bowtie2ModuleInclude
#include <string>
#include <vector>
#include <fstream>
#include <filesystem>
#include <unordered_map>
#include <boost/algorithm/string.hpp>
#include "../include/specificityModule.hpp"
#include "../include/util.hpp"

class bowtie2Module : private specificityModule
{
public:
	bowtie2Module(cracklingConfig config);
	void run(std::vector<guideResults>& candidateGuides) final;
private:
	bowtie2Config config;
	std::filesystem::path indexFile;
	bool processGuide(const guideResults& guide) final;
};

#endif // !bowtie2ModuleInclude
