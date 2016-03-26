//  GBM by Greg Ridgeway  Copyright (C) 2003
#include <algorithm>

#include "tree.h"

CCARTTree::CCARTTree()
{
    pRootNode = NULL;
    dShrink = 1.0;
}


CCARTTree::~CCARTTree()
{
	delete pRootNode;
}


void CCARTTree::Initialize
(
)
{
}


void CCARTTree::Reset() {
  delete pRootNode;
  
  iBestNode = 0;
  dBestNodeImprovement = 0.0;

  schWhichNode = 0;
  
}



//------------------------------------------------------------------------------
// Grows a regression tree
//------------------------------------------------------------------------------
void CCARTTree::grow
(
 double *adZ,
 const CDataset& data,
 const double *adW,
 const double *adF,
 unsigned long nBagged,
 double dLambda,
 unsigned long cMinObsInNode,
 std::vector<unsigned long>& aiNodeAssign,
 CNodeSearch *aNodeSearch
)
{
#ifdef NOISY_DEBUG
  Rprintf("Growing tree\n");
#endif

	if((adZ==NULL) || (adW==NULL) || (adF==NULL) ||
	 (depthOfTree < 1))
	{
	  throw GBM::invalid_argument();
	}

  double dSumZ = 0.0;
  double dSumZ2 = 0.0;
  double dTotalW = 0.0;
  
#ifdef NOISY_DEBUG
  Rprintf("initial tree calcs\n");
#endif

	for(long iObs=0; iObs<data.get_trainSize(); iObs++)
	{
		// aiNodeAssign tracks to which node each training obs belongs
		aiNodeAssign[iObs] = 0;

		if(data.GetBag()[iObs])
		{
			// get the initial sums and sum of squares and total weight
			dSumZ += adW[iObs]*adZ[iObs];
			dSumZ2 += adW[iObs]*adZ[iObs]*adZ[iObs];
			dTotalW += adW[iObs];
		}
	}

  dError = dSumZ2-dSumZ*dSumZ/dTotalW;
  pRootNode = new CNode(dSumZ/dTotalW, dTotalW, true);


  vecpTermNodes.resize(2*depthOfTree + 1,NULL); // accounts for missing nodes
  vecpTermNodes[0] = pRootNode;
  aNodeSearch[0].Set(dSumZ, dTotalW, nBagged);
  
  // build the tree structure
#ifdef NOISY_DEBUG
  Rprintf("Building tree 1 ");
#endif
  cTotalNodeCount = 1;
  cTerminalNodes = 1;
  for(long cDepth=0; cDepth<depthOfTree; cDepth++)
    {
#ifdef NOISY_DEBUG
      Rprintf("%d ",cDepth);
#endif
      GetBestSplit(data,
		   aNodeSearch,
		   aiNodeAssign,
		   adZ,
		   adW);
      
      if(dBestNodeImprovement == 0.0)
        {
	  break;
        }
      
      // setup the new nodes and add them to the tree
      aNodeSearch[iBestNode].SetupNewNodes(*vecpTermNodes[iBestNode]);
      cTotalNodeCount += 3;
      cTerminalNodes += 2;

        // assign observations to the correct node
      for(long iObs=0; iObs < data.get_trainSize(); iObs++)
        {
	  signed char iWhichNode = aiNodeAssign[iObs];
	  if(iWhichNode==iBestNode)
            {
	      schWhichNode = vecpTermNodes[iBestNode]->WhichNode(data,iObs);
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
      
      // set up the node search for the new right node
      aNodeSearch[cTerminalNodes-2].Set(aNodeSearch[iBestNode].bestSplit.RightWeightResiduals,
					aNodeSearch[iBestNode].bestSplit.RightTotalWeight,
					aNodeSearch[iBestNode].bestSplit.RightNumObs);
      // set up the node search for the new missing node
      aNodeSearch[cTerminalNodes-1].Set(aNodeSearch[iBestNode].bestSplit.MissingWeightResiduals,
					aNodeSearch[iBestNode].bestSplit.MissingTotalWeight,
					aNodeSearch[iBestNode].bestSplit.MissingNumObs);
      // set up the node search for the new left node
      // must be done second since we need info for right node first
      aNodeSearch[iBestNode].Set(aNodeSearch[iBestNode].bestSplit.LeftWeightResiduals,
				 aNodeSearch[iBestNode].bestSplit.LeftTotalWeight,
				 aNodeSearch[iBestNode].bestSplit.LeftNumObs);

      vecpTermNodes[cTerminalNodes-2] = vecpTermNodes[iBestNode]->pRightNode;
            vecpTermNodes[cTerminalNodes-1] = vecpTermNodes[iBestNode]->pMissingNode;
            vecpTermNodes[iBestNode] = vecpTermNodes[iBestNode]->pLeftNode;

    } // end tree growing

    // DEBUG
    // Print();
}


void CCARTTree::GetBestSplit
(
 const CDataset &data,
 CNodeSearch *aNodeSearch,
 std::vector<unsigned long>& aiNodeAssign,
 double *adZ,
 const double *adW
 )
{
  
  unsigned long iNode = 0;
  unsigned long iOrderObs = 0;
  unsigned long iWhichObs = 0;
  
  const CDataset::index_vector colNumbers(data.random_order());
  const CDataset::index_vector::const_iterator final = colNumbers.begin() + data.get_numFeatures();
  
  for(CDataset::index_vector::const_iterator it=colNumbers.begin();
      it != final;
      it++)
    {
      const int iVar = *it;
      const int cVarClasses = data.varclass(iVar);
      
      for(iNode=0; iNode < cTerminalNodes; iNode++)
        {
	  aNodeSearch[iNode].ResetForNewVar(iVar, cVarClasses);
        }

      // distribute the observations in order to the correct node search
      for(iOrderObs=0; iOrderObs < data.get_trainSize(); iOrderObs++)
        {
	  iWhichObs = data.order_ptr()[iVar*data.get_trainSize() + iOrderObs];
	  if(data.GetBag()[iWhichObs])
            {
	      const int iNode = aiNodeAssign[iWhichObs];
	      const double dX = data.x_value(iWhichObs, iVar);
	      aNodeSearch[iNode].IncorporateObs(dX,
						adZ[iWhichObs],
						adW[iWhichObs],
						data.monotone(iVar));
            }
        }
        for(iNode=0; iNode<cTerminalNodes; iNode++)
        {
            if(cVarClasses != 0) // evaluate if categorical split
            {
	      aNodeSearch[iNode].EvaluateCategoricalSplit();
            }
            aNodeSearch[iNode].WrapUpCurrentVariable();
        }
    }

    // search for the best split
    iBestNode = 0;
    dBestNodeImprovement = 0.0;
    for(iNode=0; iNode<cTerminalNodes; iNode++)
    {
        aNodeSearch[iNode].SetToSplit();
        if(aNodeSearch[iNode].BestImprovement() > dBestNodeImprovement)
        {
            iBestNode = iNode;
            dBestNodeImprovement = aNodeSearch[iNode].BestImprovement();
        }
    }

}


void CCARTTree::GetNodeCount
(
    int &cNodes
)
{
    cNodes = cTotalNodeCount;
}



void CCARTTree::PredictValid
(
 const CDataset &data,
 unsigned long nValid,
 double *adFadj
 )
{
  int i=0;
  
  for(i=data.nrow() - nValid; i<data.nrow(); i++)
    {
      pRootNode->Predict(data, i, adFadj[i]);
      adFadj[i] *= dShrink;
    }
}



void CCARTTree::Predict
(
    double *adX,
    unsigned long cRow,
    unsigned long cCol,
    unsigned long iRow,
    double &dFadj
)
{

    if(pRootNode)
      {
        pRootNode->Predict(adX,cRow,cCol,iRow,dFadj);
        dFadj *= dShrink;
      }
    else
      {
        dFadj = 0.0;
      }
}



void CCARTTree::Adjust
(
 const std::vector<unsigned long>& aiNodeAssign,
 double *adFadj,
 unsigned long cMinObsInNode
)
{
  unsigned long iObs = 0;
  
  pRootNode->Adjust(cMinObsInNode);
  
  // predict for the training observations
  for(iObs=0; iObs<aiNodeAssign.size(); iObs++)
    {
      adFadj[iObs] = vecpTermNodes[aiNodeAssign[iObs]]->dPrediction;
    }
}


void CCARTTree::Print()
{
    if(pRootNode)
    {
      pRootNode->PrintSubtree(0);
      Rprintf("shrinkage: %f\n",dShrink);
      Rprintf("initial error: %f\n\n",dError);
    }
}



void CCARTTree::GetVarRelativeInfluence
(
    double *adRelInf
)
{
  if(pRootNode)
    {
      pRootNode->GetVarRelativeInfluence(adRelInf);
    }
}



void CCARTTree::TransferTreeToRList
(
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

    int iNodeID = 0;

    if(pRootNode)
    {
        pRootNode->TransferTreeToRList(iNodeID,
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
      throw GBM::failure();
    }
}


