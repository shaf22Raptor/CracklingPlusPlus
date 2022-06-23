#include <sgrnascorer2.hpp>
#include <svmData.hpp>

using std::string;
using std::map;
using std::array;

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

	struct svm_parameter svmParameters;
	struct svm_problem svmProblem;
	struct svm_node* svmTrainingNodes;

	svmProblem.l = (int)data.size();
	//svmProblem.y = (double*)malloc((svmProblem.l) * sizeof(double));
	svmProblem.y = &labels[0];
	svmProblem.x = (struct svm_node**)malloc((svmProblem.l) * sizeof(struct svm_node*));
	svmTrainingNodes = (struct svm_node*)malloc(((data[0].size() + 1) * svmProblem.l) * sizeof(struct svm_node));

	int j = 0; 
	for (int i = 0; i < svmProblem.l; i++)
	{
		svmProblem.x[i] = &svmTrainingNodes[j];
		for (int k = 0; k < data[i].size(); k++, j++)
		{
			svmTrainingNodes[j].index = k + 1;
			svmTrainingNodes[j].value = data[i][k];

		}
		svmTrainingNodes[j].index = -1;
		j++;
	}

	svmParameters.svm_type = C_SVC;
	svmParameters.kernel_type = LINEAR;
	svmParameters.degree = 3;
	svmParameters.gamma = 0.1;
	svmParameters.coef0 = 0.0;
	svmParameters.cache_size = 200.;
	svmParameters.eps = 1e-9;
	svmParameters.C = 1;
	svmParameters.nr_weight = 0;
	svmParameters.weight_label = NULL;
	svmParameters.weight = NULL;
	svmParameters.nu = 0.0;
	svmParameters.p = 0.1;
	svmParameters.shrinking = 1;
	svmParameters.probability = 0;

	sgRNAScorer2Model = svm_train(&svmProblem, &svmParameters);

	for (int i = 0; i < svmProblem.l; i++)
	{
		free(svmProblem.x[i]);
	}
	free(svmTrainingNodes);
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

		string seqUpper = target23;
		std::transform(seqUpper.begin(), seqUpper.end(), seqUpper.begin(), [](unsigned char c) { return std::toupper(c); });

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

