//  GBM by Greg Ridgeway  Copyright (C) 2003

#include "node.h"

CNode::CNode(double nodePrediction,
		double trainingWeight, long numObs, bool terminalFlag):aiLeftCategory()
{
    dPrediction = nodePrediction;
    dTrainW = trainingWeight;
    isTerminal = terminalFlag;
    isContinuous = true;
    cN = numObs;

    dSplitValue = 0.0;
    iSplitVar = 0;
    dImprovement = 0.0;

    // Set children to NULL
	pLeftNode = NULL;
	pRightNode = NULL;
	pMissingNode = NULL;
}


CNode::~CNode()
{
	// Each node is responsible for deleting its
	// children
    delete pLeftNode;
    delete pRightNode;
    delete pMissingNode;
}

void CNode::Adjust
(
    unsigned long cMinObsInNode
)
{
	// Only adjust if node is not terminal
	if(!isTerminal)
	{
		pLeftNode->Adjust(cMinObsInNode);
		pRightNode->Adjust(cMinObsInNode);

		if(pMissingNode->isTerminal && (pMissingNode->cN < cMinObsInNode))
		{
			dPrediction = ((pLeftNode->dTrainW)*(pLeftNode->dPrediction) +
				 (pRightNode->dTrainW)*(pRightNode->dPrediction))/
			(pLeftNode->dTrainW + pRightNode->dTrainW);
			pMissingNode->dPrediction = dPrediction;
		}
		else
		{
			pMissingNode->Adjust(cMinObsInNode);
			dPrediction =
			((pLeftNode->dTrainW)*   (pLeftNode->dPrediction) +
			(pRightNode->dTrainW)*  (pRightNode->dPrediction) +
			(pMissingNode->dTrainW)*(pMissingNode->dPrediction))/
			(pLeftNode->dTrainW + pRightNode->dTrainW + pMissingNode->dTrainW);
		}
	}

}



void CNode::Predict
(
    const CDataset &data,
    unsigned long iRow,
    double &dFadj
)
{
	// If node is terminal set the function adjustment to the
	// prediction.  Else move down the tree.
	if(isTerminal)
	{
		dFadj = dPrediction;
	}
	else
	{
		signed char schWhichNode = WhichNode(data,iRow);
		if(schWhichNode == -1)
		{
		  pLeftNode->Predict(data, iRow, dFadj);
		}
		else if(schWhichNode == 1)
		{
		  pRightNode->Predict(data, iRow, dFadj);
		}
		else
		{
		  pMissingNode->Predict(data, iRow, dFadj);
		}
	}

}


void CNode::Predict
(
    double *adX,
    unsigned long cRow,
    unsigned long cCol,
    unsigned long iRow,
    double &dFadj
)
{
	// If node is terminal set the function adjustment to the
	// prediction.  Else move down the tree.
	if(isTerminal)
	{
		dFadj = dPrediction;
	}
	else
	{
		signed char schWhichNode = WhichNode(adX,cRow,cCol,iRow);
		if(schWhichNode == -1)
		{
		  pLeftNode->Predict(adX,cRow,cCol,iRow,dFadj);
		}
		else if(schWhichNode == 1)
		{
		  pRightNode->Predict(adX,cRow,cCol,iRow,dFadj);
		}
		else
		{
		  pMissingNode->Predict(adX,cRow,cCol,iRow,dFadj);
		}
	}

}


void CNode::GetVarRelativeInfluence
(
    double *adRelInf
)
{
	// Relative influence of split variable only updated in non-terminal nodes
	if(!isTerminal)
	{
		adRelInf[iSplitVar] += dImprovement;
		pLeftNode->GetVarRelativeInfluence(adRelInf);
		pRightNode->GetVarRelativeInfluence(adRelInf);
	}

}

void CNode::ApplyShrinkage(double dLambda)
{
	dPrediction *= dLambda;
}


void CNode::PrintSubtree
(
 unsigned long cIndent
)
{
  unsigned long i = 0;

  if(this->isTerminal)
  {
	  for(i=0; i< cIndent; i++) Rprintf("  ");
	  Rprintf("N=%f, Prediction=%f *\n",
		  dTrainW,
		  dPrediction);
  }
  else
  {
	  if(this->isContinuous)
	    {
	  	  for(i=0; i< cIndent; i++) Rprintf("  ");
	  	    Rprintf("N=%f, Improvement=%f, Prediction=%f, NA pred=%f\n",
	  	  	  dTrainW,
	  	  	  dImprovement,
	  	  	  dPrediction,
	  	  	  (pMissingNode == NULL ? 0.0 : pMissingNode->dPrediction));

	  	    for(i=0; i< cIndent; i++) Rprintf("  ");
	  	    Rprintf("V%d < %f\n",
	  	  	  iSplitVar,
	  	  	  dSplitValue);
	  	    pLeftNode->PrintSubtree(cIndent+1);

	  	    for(i=0; i< cIndent; i++) Rprintf("  ");
	  	    Rprintf("V%d > %f\n",
	  	  	  iSplitVar,
	  	  	  dSplitValue);
	  	    pRightNode->PrintSubtree(cIndent+1);

	  	    for(i=0; i< cIndent; i++) Rprintf("  ");
	  	    Rprintf("missing\n");
	  	    pMissingNode->PrintSubtree(cIndent+1);
	    }
	    else
	    {
	  	  const std::size_t cLeftCategory = aiLeftCategory.size();

	  	    for(i=0; i< cIndent; i++) Rprintf("  ");
	  	    Rprintf("N=%f, Improvement=%f, Prediction=%f, NA pred=%f\n",
	  	  	  dTrainW,
	  	  	  dImprovement,
	  	  	  dPrediction,
	  	  	  (pMissingNode == NULL ? 0.0 : pMissingNode->dPrediction));

	  	    for(i=0; i< cIndent; i++) Rprintf("  ");
	  	    Rprintf("V%d in ",iSplitVar);
	  	    for(i=0; i<cLeftCategory; i++)
	  	      {
	  	        Rprintf("%d",aiLeftCategory[i]);
	  	        if(i<cLeftCategory-1) Rprintf(",");
	  	      }
	  	    Rprintf("\n");
	  	    pLeftNode->PrintSubtree(cIndent+1);

	  	    for(i=0; i< cIndent; i++) Rprintf("  ");
	  	    Rprintf("V%d not in ",iSplitVar);
	  	    for(i=0; i<cLeftCategory; i++)
	  	      {
	  	        Rprintf("%d",aiLeftCategory[i]);
	  	        if(i<cLeftCategory-1) Rprintf(",");
	  	      }
	  	    Rprintf("\n");
	  	    pRightNode->PrintSubtree(cIndent+1);

	  	    for(i=0; i< cIndent; i++) Rprintf("  ");
	  	    Rprintf("missing\n");
	  	    pMissingNode->PrintSubtree(cIndent+1);
	    }
  }


}

void CNode::SplitNode(SplitParams bestSplit, std::vector<long>& bestCategory)
{
	// set up a continuous split
	if(bestSplit.SplitClass==0)
	{
		isTerminal = false;
		dSplitValue = bestSplit.SplitValue;
		iSplitVar = bestSplit.SplitVar;
	}
	else
	{
		isTerminal = false;
		isContinuous = false;
		iSplitVar = bestSplit.SplitVar;
		aiLeftCategory.resize(1 + (ULONG)bestSplit.SplitValue);
							  std::copy(bestCategory.begin(),
									  	  bestCategory.begin() + aiLeftCategory.size(),
										 aiLeftCategory.begin());
	}

	dImprovement = bestSplit.ImprovedResiduals ;
	pLeftNode    = new CNode(bestSplit.LeftWeightResiduals/bestSplit.LeftTotalWeight, bestSplit.LeftTotalWeight,
									bestSplit.LeftNumObs, true);
	pRightNode   = new CNode(bestSplit.RightWeightResiduals/bestSplit.RightTotalWeight,
							bestSplit.RightTotalWeight, bestSplit.RightNumObs, true);

	pMissingNode = new CNode(bestSplit.MissingWeightResiduals/bestSplit.MissingTotalWeight,
							bestSplit.MissingTotalWeight, bestSplit.MissingNumObs, true);

}

signed char CNode::WhichNode
(
    const CDataset &data,
    unsigned long iObs
)
{
    signed char ReturnValue = 0;
    double dX = data.x_value(iObs, iSplitVar);

    if(isContinuous)
    {
    	 if(!ISNA(dX))
    	    {
    	        if(dX < dSplitValue)
    	        {
    	            ReturnValue = -1;
    	        }
    	        else
    	        {
    	            ReturnValue = 1;
    	        }
    	    }
    	    // if missing value returns 0

    	    return ReturnValue;
    }
    else
    {
    	if(!ISNA(dX))
    	    {
    	      if(std::find(aiLeftCategory.begin(),
    			   aiLeftCategory.end(),
    			   (ULONG)dX) != aiLeftCategory.end())
    	        {
    	            ReturnValue = -1;
    	        }
    	        else
    	        {
    	            ReturnValue = 1;
    	        }
    	    }
    	    // if missing value returns 0

    	    return ReturnValue;
    }

}


signed char CNode::WhichNode
(
    double *adX,
    unsigned long cRow,
    unsigned long cCol,
    unsigned long iRow
)
{
    signed char ReturnValue = 0;
    double dX = adX[iSplitVar*cRow + iRow];

    if(isContinuous)
    {
    	if(!ISNA(dX))
    	    {
    	        if(dX < dSplitValue)
    	        {
    	            ReturnValue = -1;
    	        }
    	        else
    	        {
    	            ReturnValue = 1;
    	        }
    	    }
    	    // if missing value returns 0

    	    return ReturnValue;
    }
    else
    {
    	 if(!ISNA(dX))
    	    {
    	      if(std::find(aiLeftCategory.begin(),
    			   aiLeftCategory.end(),
    			   (ULONG)dX) != aiLeftCategory.end())
    	        {
    	            ReturnValue = -1;
    	        }
    	        else
    	        {
    	            ReturnValue = 1;
    	        }
    	    }
    	    // if missing value returns 0

    	    return ReturnValue;
    }

}

void CNode::TransferTreeToRList
(
    int &iNodeID,
    const CDataset &data,
    int *aiSplitVar,
    double *adSplitPoint,
    int *aiLeftNode,
    int *aiRightNode,
    int *aiMissingNode,
    double *adErrorReduction,
    double *adWeight,
    double *adPred,
    VEC_VEC_CATEGORIES &vecSplitCodes,
    int cCatSplitsOld,
    double dShrinkage
)
{
	if(this->isTerminal)
	{
		aiSplitVar[iNodeID] = -1;
		adSplitPoint[iNodeID] = dShrinkage*dPrediction;
		aiLeftNode[iNodeID] = -1;
		aiRightNode[iNodeID] = -1;
		aiMissingNode[iNodeID] = -1;
		adErrorReduction[iNodeID] = 0.0;
		adWeight[iNodeID] = dTrainW;
		adPred[iNodeID] = dShrinkage*dPrediction;

		iNodeID++;
	}
	else if(this->isContinuous)
	{
		int iThisNodeID = iNodeID;
		aiSplitVar[iThisNodeID] = iSplitVar;
		adSplitPoint[iThisNodeID] = dSplitValue;
		adErrorReduction[iThisNodeID] = dImprovement;
		adWeight[iThisNodeID] = dTrainW;
		adPred[iThisNodeID] = dShrinkage*dPrediction;

		iNodeID++;
		aiLeftNode[iThisNodeID] = iNodeID;
		pLeftNode->TransferTreeToRList(iNodeID,
					 data,
					 aiSplitVar,
					 adSplitPoint,
					 aiLeftNode,
					 aiRightNode,
					 aiMissingNode,
					 adErrorReduction,
					 adWeight,
					 adPred,
					 vecSplitCodes,
					 cCatSplitsOld,
					 dShrinkage);

		aiRightNode[iThisNodeID] = iNodeID;
		pRightNode->TransferTreeToRList(iNodeID,
					  data,
					  aiSplitVar,
					  adSplitPoint,
					  aiLeftNode,
					  aiRightNode,
					  aiMissingNode,
					  adErrorReduction,
					  adWeight,
					  adPred,
					  vecSplitCodes,
					  cCatSplitsOld,
					  dShrinkage);

		aiMissingNode[iThisNodeID] = iNodeID;
		pMissingNode->TransferTreeToRList(iNodeID,
						data,
						aiSplitVar,
						adSplitPoint,
						aiLeftNode,
						aiRightNode,
						aiMissingNode,
						adErrorReduction,
						adWeight,
						adPred,
						vecSplitCodes,
						cCatSplitsOld,
						dShrinkage);
	}
	else
	{
		int iThisNodeID = iNodeID;
			    unsigned long cCatSplits = vecSplitCodes.size();
			    unsigned long i = 0;
			    int cLevels = data.varclass(iSplitVar);
			    const std::size_t cLeftCategory = aiLeftCategory.size();

			    aiSplitVar[iThisNodeID] = iSplitVar;
			    adSplitPoint[iThisNodeID] = cCatSplits+cCatSplitsOld; // 0 based
			    adErrorReduction[iThisNodeID] = dImprovement;
			    adWeight[iThisNodeID] = dTrainW;
			    adPred[iThisNodeID] = dShrinkage*dPrediction;

			    vecSplitCodes.push_back(VEC_CATEGORIES());

			    vecSplitCodes[cCatSplits].resize(cLevels,1);
			    for(i=0; i<cLeftCategory; i++)
			      {
			        vecSplitCodes[cCatSplits][aiLeftCategory[i]] = -1;
			      }

			    iNodeID++;
			    aiLeftNode[iThisNodeID] = iNodeID;
			    pLeftNode->TransferTreeToRList(iNodeID,
			  				 data,
			  				 aiSplitVar,
			  				 adSplitPoint,
			  				 aiLeftNode,
			  				 aiRightNode,
			  				 aiMissingNode,
			  				 adErrorReduction,
			  				 adWeight,
			  				 adPred,
			  				 vecSplitCodes,
			  				 cCatSplitsOld,
			  				 dShrinkage);
			    aiRightNode[iThisNodeID] = iNodeID;
			    pRightNode->TransferTreeToRList(iNodeID,
			  				  data,
			  				  aiSplitVar,
			  				  adSplitPoint,
			  				  aiLeftNode,
			  				  aiRightNode,
			  				  aiMissingNode,
			  				  adErrorReduction,
			  				  adWeight,
			  				  adPred,
			  				  vecSplitCodes,
			  				  cCatSplitsOld,
			  				  dShrinkage);

			    aiMissingNode[iThisNodeID] = iNodeID;
			    pMissingNode->TransferTreeToRList(iNodeID,
			  				    data,
			  				    aiSplitVar,
			  				    adSplitPoint,
			  				    aiLeftNode,
			  				    aiRightNode,
			  				    aiMissingNode,
			  				    adErrorReduction,
			  				    adWeight,
			  				    adPred,
			  				    vecSplitCodes,
			  				    cCatSplitsOld,
			  				    dShrinkage);
	}

}






