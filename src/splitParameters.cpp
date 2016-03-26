#include "splitParameters.h"

SplitParams::SplitParams()
{

}

SplitParams::~SplitParams()
{

}


void SplitParams::UpdateMissingNode(double predIncrement, double trainWIncrement, long numIncrement)
{
	// Move data point from right node to missing
	MissingWeightResiduals += predIncrement;
	MissingTotalWeight += trainWIncrement;
	MissingNumObs += numIncrement;

	RightWeightResiduals -= predIncrement;
	RightTotalWeight -= trainWIncrement;
	RightNumObs -= numIncrement;
}

void SplitParams::UpdateLeftNode(double predIncrement, double trainWIncrement, long numIncrement)
{
	// Move data point from right node to left node
	LeftWeightResiduals += predIncrement;
	LeftTotalWeight += trainWIncrement;
	LeftNumObs += numIncrement;

	RightWeightResiduals -= predIncrement;
	RightTotalWeight -= trainWIncrement;
	RightNumObs -= numIncrement;

}

void SplitParams::NodeGradResiduals()
{
	// Returns weighted

	double dTemp = 0.0;
	double dResult = 0.0;

	// Only need to look at left and right
	if(MissingNumObs == 0.0)
	{
		dTemp = LeftWeightResiduals/LeftTotalWeight - RightWeightResiduals/RightTotalWeight;
		dResult = LeftTotalWeight*RightTotalWeight*dTemp*dTemp/(LeftTotalWeight+RightTotalWeight);
	}
	else
	{
		// Grad - left/right
		dTemp = LeftWeightResiduals/LeftTotalWeight - RightWeightResiduals/RightTotalWeight;
		dResult += LeftTotalWeight*RightTotalWeight*dTemp*dTemp;

		// Grad - left/missing
		dTemp = LeftWeightResiduals/LeftTotalWeight - MissingWeightResiduals/MissingTotalWeight;
		dResult += LeftTotalWeight*MissingTotalWeight*dTemp*dTemp;

		// Grad - right/missing
		dTemp = RightWeightResiduals/RightTotalWeight - MissingWeightResiduals/MissingTotalWeight;
		dResult += RightTotalWeight*MissingTotalWeight*dTemp*dTemp;
		dResult /= (LeftTotalWeight + RightTotalWeight + MissingTotalWeight);
	}

	// Update current residuals
	ImprovedResiduals = dResult;
}

bool SplitParams::SplitIsCorrMonotonic(long specifyMonotone)
{
	double weightedGrad = RightWeightResiduals * LeftTotalWeight- LeftWeightResiduals * RightTotalWeight;
	return (specifyMonotone == 0 || specifyMonotone * weightedGrad > 0);
}

bool SplitParams::HasMinNumOfObs(long minObsInNode)
{
	return ((LeftNumObs >= minObsInNode) &&
				(RightNumObs >= minObsInNode));
}
