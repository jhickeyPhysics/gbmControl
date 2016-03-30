//-----------------------------------
//
// File: gbmTreeComps.cpp
//
// Description: class that implements the public methods of the
//    gbm engine tree components.
//
//-----------------------------------

//------------------------------
// Includes
//------------------------------
#include "gbmTreeComps.h"

//----------------------------------------
// Function Members - Public
//----------------------------------------
//-----------------------------------
// Function: CTreeComps
//
// Returns: none
//
// Description: Default constructor for the gbm engine tree components.
//
// Parameters: double - the shrinkage parameter
//    unsigned long - number of training examples
//    unsigned long - number of features
//    double - fraction of training examples in bag
//    unsigned long - depth of the tree to grow
//    unsigned long - minimum number of observations in node
//    int - number of groups in data
//
//-----------------------------------
CTreeComps::CTreeComps(double dLambda,
	    unsigned long cDepth,
	    unsigned long cMinObsInNode)
{
	this-> dLambda = dLambda;
	this-> cMinObsInNode = cMinObsInNode;
	aNodeSearch.Initialize(cMinObsInNode);
	ptreeTemp.reset(new CCARTTree);
	ptreeTemp->SetDepth(cDepth);

}

//-----------------------------------
// Function: ~CTreeComps()
//
// Returns: none
//
// Description: default destructor for the gbmTreeComps
//
// Parameters: none
//
//-----------------------------------
CTreeComps::~CTreeComps()
{
}

//-----------------------------------
// Function: Initialize
//
// Returns: none
//
// Description: initializes the tree components
//
// Parameters: const CDataset* - ptr to the data object in GBM
//
//-----------------------------------
void CTreeComps::TreeInitialize(const CDataset* pData)
{
	// aiNodeAssign tracks to which node each training obs belongs
	aiNodeAssign.resize(pData->get_trainSize());
}


//-----------------------------------
// Function: GrowTrees
//
// Returns: none
//
// Description: grows the tree
//
// Parameters: const CDataset ptr - pointer to the gbm data
//    int& - reference to  number of nodes in tree
//
//-----------------------------------
void CTreeComps::GrowTrees(const CDataset* pData, int& cNodes, double* adZ, const double* adFadj)
{
	#ifdef NOISY_DEBUG
	  Rprintf("Reset tree\n");
	#endif

	  //Reset tree and searcher
	  ptreeTemp->Reset();
	  aNodeSearch.Reset();

	#ifdef NOISY_DEBUG
	  Rprintf("grow tree\n");
	#endif

	ptreeTemp->grow(&(adZ[0]),
	                *(pData),
	                &(adFadj[0]),
	                cMinObsInNode,
	                aiNodeAssign,
	                 aNodeSearch);

	#ifdef NOISY_DEBUG
	  tempTree->Print();
	#endif

	  ptreeTemp->GetNodeCount(cNodes);
	#ifdef NOISY_DEBUG
	  Rprintf("get node count=%d\n",cNodes);
	#endif
}

//-----------------------------------
// Function: AdjustAndShrink
//
// Returns: none
//
// Description: adjusts the tree and shrinks.
//
// Parameters: none
//
//-----------------------------------
void CTreeComps::AdjustAndShrink(double * adFadj)
{
	ptreeTemp->Adjust(aiNodeAssign,
	                  &(adFadj[0]),
	                  cMinObsInNode);
	ptreeTemp->SetShrinkage(dLambda);
	#ifdef NOISY_DEBUG
	  ptreeTemp->Print();
	#endif
}

//-----------------------------------
// Function: TransferTreeToRList
//
// Returns: none
//
// Description:
//
// Parameters:
//
//-----------------------------------
void CTreeComps::TransferTreeToRList(const CDataset &pData,
	     int *aiSplitVar,
	     double *adSplitPoint,
	     int *aiLeftNode,
	     int *aiRightNode,
	     int *aiMissingNode,
	     double *adErrorReduction,
	     double *adWeight,
	     double *adPred,
	     VEC_VEC_CATEGORIES &vecSplitCodes,
	     int cCatSplitsOld)
{
	int iNodeID = 0;

	if(ptreeTemp->GetRootNode())
	{
		ptreeTemp->GetRootNode()->TransferTreeToRList(iNodeID,
													   pData,
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
													   dLambda);
	}
	else
	{
	  throw GBM::failure("Can't transfer to list - RootNode does not exist.");
	}
}

//-----------------------------------
// Function: PredictValid
//
// Returns: none
//
// Description: makes predictions on validation set
//
// Parameters: const CDataset ptr - ptr to gbm data
//
void CTreeComps::PredictValid(const CDataset* pData, double* adFadj)
{
	ptreeTemp->PredictValid(*(pData), pData->GetValidSize(), &(adFadj[0]));
}

//-----------------------------------
// Function: GetNodeAssign
//
// Returns: vector<unsigned long>
//
// Description: getter node assignments
//
// Parameters: none
//
//-----------------------------------
std::vector<unsigned long> CTreeComps::GetNodeAssign()
{
	return aiNodeAssign;
}

//-----------------------------------
// Function: GetTermNodes
//
// Returns: VEC_P_NODETERMINAL
//
// Description: getter for terminal nodes
//
// Parameters: none
//
//-----------------------------------
vector<CNode*> CTreeComps::GetTermNodes()
{
	return ptreeTemp->GetTermNodes();
}

//-----------------------------------
// Function: GetLambda
//
// Returns: double
//
// Description: get shrinkage
//
// Parameters: none
//
//-----------------------------------
double CTreeComps::GetLambda()
{
	return dLambda;
}

const double CTreeComps::GetLambda() const
{
	return dLambda;
}

//-----------------------------------
// Function: GetMinNodeObs
//
// Returns: unsigned long
//
// Description: get min no of observation in node
//
// Parameters: none
//
//-----------------------------------
unsigned long CTreeComps::GetMinNodeObs()
{
	return cMinObsInNode;
}




