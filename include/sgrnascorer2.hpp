// sgrnascorer2.hpp
#include <string>
#include <map>
#include <array>
#include "../include/svm.h"
#include "../include/ConfigManager.hpp"
#include "../include/Constants.hpp"
#include "../include/Helpers.hpp"

// Variables for model training

class sgrnascorer2
{
public:
	//TODO: add rule of 5(6) to round out this class

	explicit sgrnascorer2(ConfigManager& cm);

	//TODO: Custom deconstructor to properly handle pointer destruction

	void run(std::map<std::string, std::map<std::string, std::string, std::less<>>, std::less<>>& candidateGuides) const;

private:
	bool toolIsSelected;
	std::string optimsationLevel;
	int toolCount;
	int consensusN;
	float scoreThreshold;
	struct svm_model* sgRNAScorer2Model;

};
