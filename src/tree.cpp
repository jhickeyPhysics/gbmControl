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


void CCARTTree::Initialize()
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
 const double *adF,
 unsigned long cMinObsInNode,
 std::vector<unsigned long>& aiNodeAssign,
 CNodeSearch *aNodeSearch
)
{
#ifdef NOISY_DEBUG
  Rprintf("Growing tree\n");
#endif

	if((adZ==NULL) || (data.weight_ptr()==NULL) || (adF==NULL) ||
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
			dSumZ += data.weight_ptr()[iObs]*adZ[iObs];
			dSumZ2 += data.weight_ptr()[iObs]*adZ[iObs]*adZ[iObs];
			dTotalW += data.weight_ptr()[iObs];
		}
	}

  dError = dSumZ2-dSumZ*dSumZ/dTotalW;
  pRootNode = new CNode(dSumZ/dTotalW, dTotalW, data.GetTotalInBag(), true);


  vecpTermNodes.resize(2*depthOfTree + 1,NULL); // accounts for missing nodes
  vecpTermNodes[0] = pRootNode;
  aNodeSearch[0].Set(*pRootNode);
  
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
      
      unsigned long iNode = 0;
      unsigned long iOrderObs = 0;
      unsigned long iWhichObs = 0;

      const CDataset::index_vector colNumbers(data.random_order());
      const CDataset::index_vector::const_iterator final = colNumbers.begin() + data.get_numFeatures();

      // Loop over nodes
      for(iNode = 0; iNode < cTerminalNodes; iNode++)
      {
    	  // Loop over variables
    	  for(CDataset::index_vector::const_iterator it=colNumbers.begin();
    			  it != final;
    			  it++)
    	  {
    		  aNodeSearch[iNode].ResetForNewVar(*it, data.varclass(*it));


    		  // Loop over observations
    		  for(iOrderObs=0; iOrderObs < data.get_trainSize(); iOrderObs++)
    		  {
    			  //Get Observation
    			  iWhichObs = data.order_ptr()[(*it)*data.get_trainSize() + iOrderObs];
				  if(aiNodeAssign[iWhichObs] == iNode && data.GetBag()[iWhichObs])
				  {
					  const double dX = data.x_value(iWhichObs, *it);
					  aNodeSearch[iNode].IncorporateObs(dX,
									adZ[iWhichObs],
									data.weight_ptr()[iWhichObs],
									data.monotone(*it));
				  }

    		  }

			  if(data.varclass(*it) != 0) // evaluate if categorical split
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
      aNodeSearch[cTerminalNodes-2].Set(*(vecpTermNodes[iBestNode]->pRightNode));
      // set up the node search for the new missing node
      aNodeSearch[cTerminalNodes-1].Set(*(vecpTermNodes[iBestNode]->pMissingNode));
      // set up the node search for the new left node
      // must be done second since we need info for right node first
      aNodeSearch[iBestNode].Set(*(vecpTermNodes[iBestNode]->pLeftNode));

      vecpTermNodes[cTerminalNodes-2] = vecpTermNodes[iBestNode]->pRightNode;
      vecpTermNodes[cTerminalNodes-1] = vecpTermNodes[iBestNode]->pMissingNode;
      vecpTermNodes[iBestNode] = vecpTermNodes[iBestNode]->pLeftNode;

    } // end tree growing

    // DEBUG
    // Print();
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


