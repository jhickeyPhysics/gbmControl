//-----------------------------------
//
// File: varsplitter.cpp
//
// Description: class that implements the splitting of a node on a variable.
//
//-----------------------------------

//-----------------------------------
// Includes
//-----------------------------------
#include "varsplitter.h"


VarSplitter::VarSplitter(unsigned long minNumObs)
{

	hasBestSplit = false;
	setForNode = false;

	InitTotalWeight = 0.0;
	InitWeightResiduals = 0.0;
	InitNumObs = 0;

	dLastXValue = -HUGE_VAL;
	minObsInNode = minNumObs;

	// Hard Code Max
	adGroupSumZ.resize(1024);
	adGroupW.resize(1024);
	acGroupN.resize(1024);
	groupdMeanAndCategory.resize(1024);

}

VarSplitter::~VarSplitter()
{
}

void VarSplitter::IncorporateObs
(
    double dX,
    double dZ,
    double dW,
    long lMonotone
)
{
	if(hasBestSplit) return;
	if(ISNA(dX))
	{
		proposedSplit.UpdateMissingNode(dW*dZ, dW);
	}
	else if(proposedSplit.SplitClass == 0)
	{
		if(dLastXValue > dX)
		{
			throw GBM::failure("Observations are not in order. gbm() was unable to build an index for the design matrix. Could be a bug in gbm or an unusual data type in data.");
		}

		// Evaluate the current split
		// the newest observation is still in the right child
		proposedSplit.SplitValue = 0.5*(dLastXValue + dX);

		if((dLastXValue != dX) &&
			proposedSplit.HasMinNumOfObs(minObsInNode) &&
			proposedSplit.SplitIsCorrMonotonic(lMonotone))
		{
			proposedSplit.NodeGradResiduals();

			if(proposedSplit.HasMinNumOfObs(minObsInNode) &&
					(proposedSplit.ImprovedResiduals > bestSplit.ImprovedResiduals))
			{
				bestSplit = proposedSplit;
				WrapUpSplit();

			}

		}

		// now move the new observation to the left
		// if another observation arrives we will evaluate this
		proposedSplit.UpdateLeftNode(dW*dZ, dW);
		dLastXValue = dX;
	}
	else // variable is categorical, evaluates later
	{
		adGroupSumZ[(unsigned long)dX] += dW*dZ;
		adGroupW[(unsigned long)dX] += dW;
		acGroupN[(unsigned long)dX] ++;
	}
}


void VarSplitter::EvaluateCategoricalSplit()
{
  long i=0;
  unsigned long cFiniteMeans = 0;

  // Check if already evaluated
  if(hasBestSplit) return;

  if(proposedSplit.SplitClass == 0)
	{
	  throw GBM::invalid_argument();
	}

  cFiniteMeans = 0;
  for(i=0; i < proposedSplit.SplitClass; i++)
    {
      groupdMeanAndCategory[i].second = i;

      if(adGroupW[i] != 0.0)
      {
    	  groupdMeanAndCategory[i].first = adGroupSumZ[i]/adGroupW[i];
    	  cFiniteMeans++;
      }
      else
      {
    	  groupdMeanAndCategory[i].first = HUGE_VAL;
      }
    }

  std::sort(groupdMeanAndCategory.begin(), groupdMeanAndCategory.begin() + proposedSplit.SplitClass);


  // if only one group has a finite mean it will not consider
  // might be all are missing so no categories enter here
  for(i=0; (cFiniteMeans>1) && ((ULONG)i<cFiniteMeans-1); i++)
    {
      proposedSplit.SplitValue = (double)i;
      proposedSplit.UpdateLeftNode(adGroupSumZ[groupdMeanAndCategory[i].second], adGroupW[groupdMeanAndCategory[i].second],
    		  	  	  	  	  	  acGroupN[groupdMeanAndCategory[i].second]);
      proposedSplit.NodeGradResiduals();
      proposedSplit.setBestCategory(groupdMeanAndCategory);

      if(proposedSplit.HasMinNumOfObs(minObsInNode)
    		  && (proposedSplit.ImprovedResiduals > bestSplit.ImprovedResiduals))
      {
    	  bestSplit = proposedSplit;
		  WrapUpSplit();

      }

    }
}

void VarSplitter::SetForNode(CNode& nodeToSplit)
{
	// If not set for this node then
	if(!setForNode)
	{
		InitWeightResiduals = nodeToSplit.dPrediction * nodeToSplit.dTrainW;
		InitTotalWeight = nodeToSplit.dTrainW;
		InitNumObs = nodeToSplit.cN;
	}
	hasBestSplit = !(nodeToSplit.splitType == none);
}

void VarSplitter::SetForVariable(unsigned long iWhichVar, long cVarClasses)
{

	if(hasBestSplit) return;
	if (int(cVarClasses) > adGroupSumZ.size())
	{
	throw GBM::failure("too many variable classes");
	}

	std::fill(adGroupSumZ.begin(), adGroupSumZ.begin() + cVarClasses, 0);
	std::fill(adGroupW.begin(), adGroupW.begin() + cVarClasses, 0);
	std::fill(acGroupN.begin(), acGroupN.begin() + cVarClasses, 0);

	bestSplit.ResetSplitProperties(InitWeightResiduals, InitTotalWeight, InitNumObs);
	proposedSplit.ResetSplitProperties(InitWeightResiduals, InitTotalWeight, InitNumObs,
		  proposedSplit.SplitValue,	cVarClasses, iWhichVar);

	dLastXValue = -HUGE_VAL;

}

void VarSplitter::WrapUpSplit()
{
	if(proposedSplit.MissingNumObs <= 0)
	{
		bestSplit.MissingWeightResiduals   = InitWeightResiduals;
		bestSplit.MissingTotalWeight = InitTotalWeight;
		bestSplit.MissingNumObs      = 0;
	}
	else
	{
		bestSplit.MissingWeightResiduals = proposedSplit.MissingWeightResiduals;
		bestSplit.MissingTotalWeight = proposedSplit.MissingTotalWeight;
		bestSplit.MissingNumObs = proposedSplit.MissingNumObs;

	}
}
