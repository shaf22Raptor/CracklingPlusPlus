#include <sgrnascorer2.hpp>

using std::string;
using std::map;
using std::array;

map <char, string> encoding = {
	{'A' , "0001"},
	{'C' , "0010"},
	{'T' , "0100"},
	{'G' , "1000"},
	{'K' , "1100"},
	{'M' , "0011"},
	{'R' , "1001"},
	{'Y' , "0110"},
	{'S' , "1010"},
	{'W' , "0101"},
	{'B' , "1110"},
	{'V' , "1011"},
	{'H' , "0111"},
	{'D' , "1101"},
	{'N' , "1111"},
};

sgrnascorer2::sgrnascorer2(ConfigManager cm) :
	toolIsSelected(false),
	optimsationLevel(""),
	toolCount(0),
	consensusN(0),
	scoreThreshold(0.0f)
{
	toolIsSelected = cm.getBool("consensus", "mm10db");
	optimsationLevel = cm.getString("general", "optimisation");
	toolCount = cm.getConsensusToolCount();
	consensusN = cm.getInt("consensus", "n");
	scoreThreshold = cm.getFloat("sgrnascorer2", "score-threshold");
	sgRNAScorer2Model = svm_load_model(cm.getCString("sgrnascorer2", "model"));
}

void sgrnascorer2::run(map<string, map<string, string>>& candidateGuides)
{

	if (!toolIsSelected)
	{
		printer("sgRNAScorer2 has been configured not to run. Skipping sgRNAScorer2");
		return;
	}

	printer("sgRNAScorer2 - score using model.");
	int failedCount = 0;
	int testedCount = 0;
	for (auto const& [target23, resultsMap] : candidateGuides)
	{

		// Run time filtering
		if (!filterCandidateGuides(resultsMap, MODULE_SGRNASCORER2, optimsationLevel, consensusN, toolCount)) { continue; }

		string seqUpper = makeUpper(target23);

		// TODO: Helper function to encode. All data should eventually be bit encoded.
		array<double, 80> encodedSeq;
		for (int i = 0; i < 20; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				encodedSeq[(i*4)+j] = (double)encoding[seqUpper[i]][j] - 48;
			}
		}

		struct svm_node* nodeToTest = (struct svm_node*)malloc((encodedSeq.size() + 1) * sizeof(struct svm_node));

		int j = 0;
		for (int k = 0; k < encodedSeq.size(); k++, j++)
		{
			nodeToTest[j].index = k + 1; //index of value
			nodeToTest[j].value = encodedSeq[k]; //value

		}
		nodeToTest[j].index = -1;//state the end of data vector

		double predictedValue;

		svm_predict_values(sgRNAScorer2Model, nodeToTest, &predictedValue);

		free(nodeToTest);

		candidateGuides[target23]["sgrnascorer2score"] = std::to_string(predictedValue);

		if (predictedValue < scoreThreshold)
		{
			candidateGuides[target23]["acceptedBySgRnaScorer"] = CODE_REJECTED;
			failedCount++;
		}
		else
		{
			candidateGuides[target23]["acceptedBySgRnaScorer"] = CODE_ACCEPTED;
		}
		testedCount++;
	}
	char printingBuffer[1024];
	snprintf(printingBuffer, 1024, "\t%d of %d failed here.", failedCount, testedCount);
	printer(printingBuffer);
}

