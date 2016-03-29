//------------------------------------------------------------------------------
//  GBM by Greg Ridgeway  Copyright (C) 2003
//
//  File:       node_search.cpp
//
//------------------------------------------------------------------------------
#include "node_search.h"

CNodeSearch::CNodeSearch()
{
    bestSplit.SplitVar = 0;

    bestSplit.SplitValue = 0.0;
    fIsSplit = false;

    bestSplit.MissingTotalWeight = 0.0;
    proposedSplit.MissingTotalWeight = 0.0;
    bestSplit.MissingWeightResiduals = 0.0;
    proposedSplit.MissingWeightResiduals = 0.0;

    adGroupSumZ.resize(1024);
    adGroupW.resize(1024);
    acGroupN.resize(1024);
    adGroupMean.resize(1024);
    aiCurrentCategory.resize(1024);
    aiBestCategory.resize(1024);
}


CNodeSearch::~CNodeSearch()
{
}


void CNodeSearch::Initialize
(
    unsigned long cMinObsInNode
)
{
    this->cMinObsInNode = cMinObsInNode;
}


void CNodeSearch::IncorporateObs
(
    double dX,
    double dZ,
    double dW,
    long lMonotone
)
{
    static double dWZ = 0.0;

    if(fIsSplit) return;

    dWZ = dW*dZ;

    if(ISNA(dX))
    {
        proposedSplit.UpdateMissingNode(dWZ, dW);
    }
    else if(proposedSplit.SplitClass == 0)   // variable is continuous
    {
        if(dLastXValue > dX)
        {
        	throw GBM::failure("Observations are not in order. gbm() was unable to build an index for the design matrix. Could be a bug in gbm or an unusual data type in data.");
        }

        // Evaluate the current split
        // the newest observation is still in the right child
        proposedSplit.SplitValue = 0.5*(dLastXValue + dX);
        if((dLastXValue != dX) &&
            proposedSplit.HasMinNumOfObs(cMinObsInNode) &&
            proposedSplit.SplitIsCorrMonotonic(lMonotone))
        {
        	proposedSplit.NodeGradResiduals();
            if(proposedSplit.ImprovedResiduals > bestSplit.ImprovedResiduals )
            {
            	bestSplit = proposedSplit;
            }
        }

        // now move the new observation to the left
        // if another observation arrives we will evaluate this
        proposedSplit.UpdateLeftNode(dWZ, dW);
        dLastXValue = dX;
    }
    else // variable is categorical, evaluates later
    {
        adGroupSumZ[(unsigned long)dX] += dWZ;
        adGroupW[(unsigned long)dX] += dW;
        acGroupN[(unsigned long)dX] ++;
    }
}



void CNodeSearch::Set(CNode nodeToSplit)
{
    dInitSumZ = nodeToSplit.dPrediction * nodeToSplit.dTrainW;
    dInitTotalW = nodeToSplit.dTrainW;
    cInitN = nodeToSplit.cN;

    bestSplit.ResetSplitProperties(dInitSumZ, dInitTotalW, cInitN);
    proposedSplit.ResetSplitProperties(0.0, dInitTotalW, cInitN);
    fIsSplit = !nodeToSplit.isTerminal;

}


void CNodeSearch::ResetForNewVar
(
    unsigned long iWhichVar,
    long cCurrentVarClasses
)
{
  if(fIsSplit) return;

  if (int(cCurrentVarClasses) > adGroupSumZ.size()) {
    throw GBM::failure("too many variable classes");
  }

  std::fill(adGroupSumZ.begin(), adGroupSumZ.begin() + cCurrentVarClasses, 0);
  std::fill(adGroupW.begin(), adGroupW.begin() + cCurrentVarClasses, 0);
  std::fill(acGroupN.begin(), acGroupN.begin() + cCurrentVarClasses, 0);
  
  proposedSplit.ResetSplitProperties(dInitSumZ, dInitTotalW, cInitN,
		  	  	  	  	  	  	  	  proposedSplit.SplitValue, cCurrentVarClasses, iWhichVar);
  dLastXValue = -HUGE_VAL;
}



void CNodeSearch::WrapUpCurrentVariable()
{
  if(proposedSplit.SplitVar == bestSplit.SplitVar)
    {
      if(proposedSplit.MissingNumObs > 0)
        {
		  bestSplit.MissingWeightResiduals   = proposedSplit.MissingWeightResiduals;
		  bestSplit.MissingTotalWeight = proposedSplit.MissingTotalWeight;
		  bestSplit.MissingNumObs      = proposedSplit.MissingNumObs;
        }
      else // DEBUG: consider a weighted average with parent node?
        {
		  bestSplit.MissingWeightResiduals   = dInitSumZ;
		  bestSplit.MissingTotalWeight = dInitTotalW;
		  bestSplit.MissingNumObs      = 0;
        }
    }
}



void CNodeSearch::EvaluateCategoricalSplit()
{
  long i=0;
  unsigned long cFiniteMeans = 0;
  
  if(fIsSplit) return;
  
  if(proposedSplit.SplitClass == 0)
    {
      throw GBM::invalid_argument();
    }

  cFiniteMeans = 0;
  for(i=0; i < proposedSplit.SplitClass; i++)
    {
      aiCurrentCategory[i] = i;
      if(adGroupW[i] != 0.0)
      {
    	  adGroupMean[i] = adGroupSumZ[i]/adGroupW[i];
    	  cFiniteMeans++;
      }
      else
        {
	  adGroupMean[i] = HUGE_VAL;
        }
    }
  
  rsort_with_index(&adGroupMean[0],&aiCurrentCategory[0],proposedSplit.SplitClass);
    
  // if only one group has a finite mean it will not consider
  // might be all are missing so no categories enter here
  for(i=0; (cFiniteMeans>1) && ((ULONG)i<cFiniteMeans-1); i++)
    {
      proposedSplit.SplitValue = (double)i;
      proposedSplit.UpdateLeftNode(adGroupSumZ[aiCurrentCategory[i]], adGroupW[aiCurrentCategory[i]],
    		  	  	  	  	  	  acGroupN[aiCurrentCategory[i]]);

      proposedSplit.NodeGradResiduals();

      if(proposedSplit.HasMinNumOfObs(cMinObsInNode) &&
	 (proposedSplit.ImprovedResiduals > bestSplit.ImprovedResiduals ))
        {
	  bestSplit.SplitValue = proposedSplit.SplitValue;
	  if(bestSplit.SplitVar!= proposedSplit.SplitVar)
	  {
		  bestSplit.SplitVar = proposedSplit.SplitVar;
		  bestSplit.SplitClass = proposedSplit.SplitClass;

	      std::copy(aiCurrentCategory.begin(),
			aiCurrentCategory.end(),
			aiBestCategory.begin());
	  }

	  bestSplit.LeftWeightResiduals      = proposedSplit.LeftWeightResiduals;
	  bestSplit.LeftTotalWeight    = proposedSplit.LeftTotalWeight;
	  bestSplit.LeftNumObs         = proposedSplit.LeftNumObs;
	  bestSplit.RightWeightResiduals     = proposedSplit.RightWeightResiduals;
	  bestSplit.RightTotalWeight   = proposedSplit.RightTotalWeight;
	  bestSplit.RightNumObs        = proposedSplit.RightNumObs;
	  bestSplit.ImprovedResiduals    = proposedSplit.ImprovedResiduals;
        }
    }
}




void CNodeSearch::SetupNewNodes(CNode& nodeToSplit)
{
	nodeToSplit.SplitNode();

}


void CNodeSearch::AssignToNode(CNode& terminalNode)
{
	terminalNode.childrenParams =  bestSplit;
	terminalNode.splitCategory = aiBestCategory;
}
