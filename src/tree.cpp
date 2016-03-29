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

void CCARTTree::Reset()
{
  delete pRootNode;
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
  long cTerminalNodes = 1;
  for(long cDepth=0; cDepth < depthOfTree; cDepth++)
  {
#ifdef NOISY_DEBUG
      Rprintf("%d ",cDepth);
#endif
      
      unsigned long iWhichObs = 0;
      const CDataset::index_vector colNumbers(data.random_order());
      const CDataset::index_vector::const_iterator final = colNumbers.begin() + data.get_numFeatures();

      // Loop over terminal nodes
      for(long iNode = 0; iNode < cTerminalNodes; iNode++)
      {
    	  aNodeSearch[0].Set(*vecpTermNodes[iNode]);
    	  // Loop over variables
    	  for(CDataset::index_vector::const_iterator it=colNumbers.begin();
    			  it != final;
    			  it++)
    	  {

    		  aNodeSearch[0].ResetForNewVar(*it, data.varclass(*it));
    		  // Loop over observations
    		  for(long iOrderObs=0; iOrderObs < data.get_trainSize(); iOrderObs++)
    		  {
    			  //Get Observation
    			  iWhichObs = data.order_ptr()[(*it)*data.get_trainSize() + iOrderObs];
				  if(aiNodeAssign[iWhichObs] == iNode && data.GetBag()[iWhichObs])
				  {
					  const double dX = data.x_value(iWhichObs, *it);
					  aNodeSearch[0].IncorporateObs(dX,
									adZ[iWhichObs],
									data.weight_ptr()[iWhichObs],
									data.monotone(*it));
				  }

    		  }

			  if(data.varclass(*it) != 0) // evaluate if categorical split
			  {
				  aNodeSearch[0].EvaluateCategoricalSplit();
			  }
			  aNodeSearch[0].WrapUpCurrentVariable();
		  }

    	  // Assign best split to node
    	  aNodeSearch[0].AssignToNode(*vecpTermNodes[iNode]);
	  }

	// search for the best split
	long iBestNode = 0;
	double dBestNodeImprovement = 0.0;
	for(long iNode=0; iNode < cTerminalNodes; iNode++)
	{
		aNodeSearch[0].SetToSplit();
		if(vecpTermNodes[iNode]->SplitImprovement() > dBestNodeImprovement)
		{
			iBestNode = iNode;
			dBestNodeImprovement = vecpTermNodes[iNode]->SplitImprovement();
		}
	}

	if(dBestNodeImprovement == 0.0)
	{
	  break;
	}
      
      // setup the new nodes and add them to the tree
      vecpTermNodes[iBestNode]->SplitNode();
      cTotalNodeCount += 3;
      cTerminalNodes += 2;


        // assign observations to the correct node
      for(long iObs=0; iObs < data.get_trainSize(); iObs++)
      {
    	  signed char iWhichNode = aiNodeAssign[iObs];

	  if(iWhichNode==iBestNode)
	  {
	      signed char schWhichNode = vecpTermNodes[iBestNode]->WhichNode(data,iObs);
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

      vecpTermNodes[cTerminalNodes-2] = vecpTermNodes[iBestNode]->pRightNode;
	  vecpTermNodes[cTerminalNodes-1] = vecpTermNodes[iBestNode]->pMissingNode;
	  vecpTermNodes[iBestNode] = vecpTermNodes[iBestNode]->pLeftNode;

    } // end tree growing

    // DEBUG
    // Print();
}

void CCARTTree::GetNodeCount(int &cNodes)
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


CNode* CCARTTree::GetRootNode()
{
	return pRootNode;
}
const CNode* CCARTTree::GetRootNode() const
{
	return pRootNode;
}


