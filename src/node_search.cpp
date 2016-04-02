//------------------------------------------------------------------------------
//  GBM by Greg Ridgeway  Copyright (C) 2003
//
//  File:       node_search.cpp
//
//------------------------------------------------------------------------------
#include "node_search.h"

CNodeSearch::CNodeSearch()
{
    fIsSplit = false;
    cTerminalNodes = 1;

    adGroupSumZ.resize(1024);
    adGroupW.resize(1024);
    acGroupN.resize(1024);
    groupdMeanAndCategory.resize(1024);
    proposedSplits.resize(1);
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

void CNodeSearch::Reset(const CDataset& data)
{
	cTerminalNodes = 1;
	proposedSplits.resize(data.get_numFeatures());
}

template<>
void CNodeSearch::IncorporateObs<false>
(
    double dX,
    double dZ,
    double dW,
    long lMonotone,
    SplitParams& proposedSplit
)
{

    if(fIsSplit) return;
    if(ISNA(dX))
    {
        proposedSplit.UpdateMissingNode(dW*dZ, dW);
    }
    else
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
			std::cout << "L Prop:" << proposedSplit.LeftWeightResiduals  << endl;
			std::cout << "Improv Prop: " << proposedSplit.ImprovedResiduals << endl;
			std::cout << endl;
        }

        // now move the new observation to the left
        // if another observation arrives we will evaluate this
        proposedSplit.UpdateLeftNode(dW*dZ, dW);
        dLastXValue = dX;
    }

}

template<>
void CNodeSearch::IncorporateObs<true>
(
    double dX,
    double dZ,
    double dW,
    long lMonotone,
    SplitParams& proposedSplit
)
{
	if(fIsSplit) return;
	if(ISNA(dX))
	{
		proposedSplit.UpdateMissingNode(dW*dZ, dW);
	}
	else // variable is categorical, evaluates later
	{
		adGroupSumZ[(unsigned long)dX] += dW*dZ;
		adGroupW[(unsigned long)dX] += dW;
		acGroupN[(unsigned long)dX] ++;
	}
}

void CNodeSearch::GenerateAllSplits
(
		vector<CNode*>& vecpTermNodes,
		const CDataset& data,
		double* adZ,
		vector<unsigned long>& aiNodeAssign
)
{
	unsigned long iWhichObs = 0;
	const CDataset::index_vector colNumbers(data.random_order());
	const CDataset::index_vector::const_iterator final = colNumbers.begin() + data.get_numFeatures();

	// Loop over terminal nodes
	for(long iNode = 0; iNode < cTerminalNodes; iNode++)
	{
	  Set(*vecpTermNodes[iNode]);

	  // Loop over variables - Generate splits
	  for(CDataset::index_vector::const_iterator it=colNumbers.begin();
			  it != final;
			  it++)
	  {
		  ResetForNewVar(*vecpTermNodes[iNode], *it, data.varclass(*it), proposedSplits[*it]);
		  std::cout << "VAR:" << *it << endl;
		 /* std::cout << "Left:" << proposedSplits[*it].LeftWeightResiduals << endl;
		  std::cout << "Right:" << proposedSplits[*it].RightWeightResiduals << endl;
		  std::cout << "Missing:" << proposedSplits[*it].MissingWeightResiduals << endl;
		  std::cout << "SplitVal: " << proposedSplits[*it].SplitValue << endl;
		  std::cout << "Residuals: " << proposedSplits[*it].ImprovedResiduals << endl;
		  std::cout << endl;*/

		  bool varIsCategorical = (bool) data.varclass(*it);
		  for(long iOrderObs=0; iOrderObs < data.get_trainSize(); iOrderObs++)
		  {
			  //Get Observation and add to split if needed
			  iWhichObs = data.order_ptr()[(*it)*data.get_trainSize() + iOrderObs];
			  if((aiNodeAssign[iWhichObs] == iNode) && data.GetBag()[iWhichObs])
			  {
				  const double dX = data.x_value(iWhichObs, *it);
				  if(data.varclass(*it) !=0)
				  {
					  IncorporateObs<true>(dX,
										adZ[iWhichObs],
										data.weight_ptr()[iWhichObs],
										data.monotone(*it), proposedSplits[*it]);
				  }
				  else
				  {
					  IncorporateObs<false>(dX,
										adZ[iWhichObs],
										data.weight_ptr()[iWhichObs],
										data.monotone(*it), proposedSplits[*it]);
				  }

			  }

		  }
		  if(data.varclass(*it) != 0) // evaluate if categorical split
		  {
			  EvaluateCategoricalSplit(proposedSplits[*it]);
		  }
	  }
	  // Assign best split to node
	  throw  GBM::failure("EXIT");
	  AssignToNode(*vecpTermNodes[iNode]);
	}

}


double CNodeSearch::SplitAndCalcImprovement
(
		vector<CNode*>& vecpTermNodes,
		const CDataset& data,
		vector<unsigned long>& aiNodeAssign
)
{
	// search for the best split
	long iBestNode = 0;
	double dBestNodeImprovement = 0.0;
	for(long iNode=0; iNode < cTerminalNodes; iNode++)
	{
		if(vecpTermNodes[iNode]->SplitImprovement() > dBestNodeImprovement)
		{
			iBestNode = iNode;
			dBestNodeImprovement = vecpTermNodes[iNode]->SplitImprovement();
		}
	}

	// Split Node if improvement is non-zero
	if(dBestNodeImprovement != 0.0)
	{
		//Split Node
		vecpTermNodes[iBestNode]->SplitNode();
		cTerminalNodes += 2;

		// Move data to children nodes
		ReAssignData(iBestNode, vecpTermNodes, data, aiNodeAssign);

		// Add children to terminal node list
		vecpTermNodes[cTerminalNodes-2] = vecpTermNodes[iBestNode]->pRightNode;
		vecpTermNodes[cTerminalNodes-1] = vecpTermNodes[iBestNode]->pMissingNode;
		vecpTermNodes[iBestNode] = vecpTermNodes[iBestNode]->pLeftNode;
	}


	return dBestNodeImprovement;
}

void CNodeSearch::ReAssignData
(
		long splittedNodeIndex,
		vector<CNode*>& vecpTermNodes,
		const CDataset& data,
		vector<unsigned long>& aiNodeAssign
)
{
	// assign observations to the correct node
	for(long iObs=0; iObs < data.get_trainSize(); iObs++)
	{
		if(aiNodeAssign[iObs]==splittedNodeIndex)
		{
		  signed char schWhichNode = vecpTermNodes[splittedNodeIndex]->WhichNode(data,iObs);
		  if(schWhichNode == 1) // goes right
		  {
			  aiNodeAssign[iObs] = cTerminalNodes-2;
		  }
		  else if(schWhichNode == 0) // is missing
		  {
			  aiNodeAssign[iObs] = cTerminalNodes-1;
		  }
		  // those to the left stay with the same node assignment
		  }
	}
}




void CNodeSearch::Set(CNode nodeToSplit)
{
    dInitSumZ = nodeToSplit.dPrediction * nodeToSplit.dTrainW;
    dInitTotalW = nodeToSplit.dTrainW;
    cInitN = nodeToSplit.cN;

    bestSplit.ResetSplitProperties(dInitSumZ, dInitTotalW, cInitN);
    fIsSplit = !(nodeToSplit.splitType == none);

}


void CNodeSearch::ResetForNewVar
(
	CNode nodeToSplit,
    unsigned long iWhichVar,
    long cCurrentVarClasses,
    SplitParams& proposedSplit
)
{
  if(fIsSplit) return;

  if (int(cCurrentVarClasses) > adGroupSumZ.size()) {
    throw GBM::failure("too many variable classes");
  }

  std::fill(adGroupSumZ.begin(), adGroupSumZ.begin() + cCurrentVarClasses, 0);
  std::fill(adGroupW.begin(), adGroupW.begin() + cCurrentVarClasses, 0);
  std::fill(acGroupN.begin(), acGroupN.begin() + cCurrentVarClasses, 0);
  proposedSplit.ResetSplitProperties(dInitSumZ, dInitTotalW, cInitN, proposedSplit.SplitValue,
										cCurrentVarClasses, iWhichVar);

  dLastXValue = -HUGE_VAL;
}



void CNodeSearch::WrapUpProposedSplit(SplitParams& proposedSplit)
{

	if(proposedSplit.MissingNumObs <= 0)
	{
		proposedSplit.MissingWeightResiduals   = dInitSumZ;
		proposedSplit.MissingTotalWeight = dInitTotalW;
		proposedSplit.MissingNumObs      = 0;
	}

}



void CNodeSearch::EvaluateCategoricalSplit(SplitParams& proposedSplit)
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
    }
}


void CNodeSearch::GenerateAllSplitsForVar(SplitParams& storeProposedSplit)
{

}

void CNodeSearch::AssignToNode(CNode& terminalNode)
{
	// Find the best split
	long bestSplitInd = 0;
	double bestErrImprovement = 0.0;
	double currErrImprovement = 0.0;
	for(long it = 0; it < proposedSplits.size(); it++)
	{
		currErrImprovement = proposedSplits[it].GetImprovement();
		if(currErrImprovement > bestErrImprovement)
		{
			bestErrImprovement = currErrImprovement;
			bestSplitInd = it;
		}
	}
	// Wrap up variable
	WrapUpProposedSplit(proposedSplits[bestSplitInd]);
	std::cout << "Left Residuals : " << proposedSplits[1].LeftWeightResiduals << endl;
	std::cout << "Right Residuals : " << proposedSplits[1].RightWeightResiduals << endl;
	std::cout << "Missing Residuals :" << proposedSplits[1].MissingWeightResiduals << endl;
	std::cout << bestSplitInd << endl;
	terminalNode.childrenParams =  proposedSplits[bestSplitInd];
}
