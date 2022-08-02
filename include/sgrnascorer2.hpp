// sgrnascorer2.hpp
#include <string>
#include <unordered_map>
#include <array>
#include "../include/svm.h"
#include "../include/ConfigManager.hpp"
#include "../include/Constants.hpp"
#include "../include/Helpers.hpp"

// Variables for model training

class sgrnascorer2
{
public:

	explicit sgrnascorer2(ConfigManager& cm);

	void run(std::unordered_map<std::string, std::unordered_map<std::string, std::string>>& candidateGuides) const;

private:
	bool toolIsSelected;
	std::string optimsationLevel;
	int toolCount;
	int consensusN;
	float scoreThreshold;
	struct svm_model* sgRNAScorer2Model;

};
