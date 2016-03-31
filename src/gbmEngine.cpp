//  GBM by Greg Ridgeway  Copyright (C) 2003
//#define NOISY_DEBUG
#include <algorithm>
#include "gbmEngine.h"

CGBM::CGBM()
{
	// Set up checks for initialization
    fInitialized = false;
    hasDataAndDist = false;
    hasTreeContainer = false;

    // Initialize distributions




}


CGBM::~CGBM()
{
	delete pDataCont;
	delete pTreeComp;
}

void CGBM::SetDataAndDistribution(SEXP radY, SEXP radOffset, SEXP radX, SEXP raiXOrder,
        SEXP radWeight, SEXP racVarClasses,
        SEXP ralMonotoneVar, SEXP radMisc, const std::string& family,
		const int cTrain, const int cFeatures, int& cGroups, double bagFraction)
{

	pDataCont=new CGBMDataContainer(radY, radOffset, radX, raiXOrder,
	        radWeight, racVarClasses, ralMonotoneVar, radMisc, family,
	        cTrain, cFeatures, cGroups, bagFraction);

	// Set up residuals container
	adZ.assign(pDataCont->getData()->nrow(), 0);
	hasDataAndDist = true;
}

void CGBM::SetTreeContainer(double dLambda,
   	    unsigned long cDepth,
   	    unsigned long cMinObsInNode)
{
	pTreeComp = new CTreeComps(dLambda, cDepth, cMinObsInNode);
	hasTreeContainer = true;
}

void CGBM::Initialize()
{
	// Throw error if initialization called incorrectly
	if(!(hasDataAndDist && hasTreeContainer))
	{
		throw GBM::failure("GBM object could not be built - missing: "
				"data or distribution or weak learner");
	}
	pDataCont-> Initialize();
	pTreeComp -> TreeInitialize(pDataCont->getData());
	fInitialized = true;
}

void CGBM::FitLearner
(
  double *adF,
  double &dTrainError,
  double &dValidError,
  double &dOOBagImprove,
  int &cNodes
)
{
  if(!fInitialized)
  {
    throw GBM::failure("GBM not initialized");
  }

  dTrainError = 0.0;
  dValidError = 0.0;
  dOOBagImprove = 0.0;

  // Initialize adjustments to function estimate
  std::vector<double> adFadj(pDataCont->getData()->nrow(), 0);

  // Bag data
  pDataCont->BagData();

#ifdef NOISY_DEBUG
  Rprintf("Compute working response\n");
#endif

  // Compute Residuals and fit tree
  pDataCont->ComputeResiduals(&adF[0], &adZ[0]);
  pTreeComp->GrowTrees(pDataCont->getData(), cNodes, &adZ[0], &adFadj[0]);

  // Now I have adF, adZ, and vecpTermNodes (new node assignments)
  // Fit the best constant within each terminal node
#ifdef NOISY_DEBUG
  Rprintf("fit best constant\n");
#endif

  // Adjust terminal node predictions and shrink
  pDataCont->ComputeBestTermNodePreds(&adF[0], &adZ[0], pTreeComp, cNodes);
  pTreeComp->AdjustAndShrink(&adFadj[0]);

  // Compute the error improvement within bag
  dOOBagImprove = pDataCont->ComputeBagImprovement(&adF[0], pTreeComp->GetLambda(), &adFadj[0]);

  // Update the function estimate
  unsigned long i = 0;
  for(i=0; i < pDataCont->getData()->get_trainSize(); i++)
  {
    adF[i] += pTreeComp->GetLambda() * adFadj[i];

  }

  // Make validation predictions
  dTrainError = pDataCont->ComputeDeviance(&adF[0], false);
  pTreeComp->PredictValid(pDataCont->getData(), &adFadj[0]);

  for(i=pDataCont->getData()->get_trainSize();
	  i < pDataCont->getData()->get_trainSize()+pDataCont->getData()->GetValidSize();
	  i++)
  {
    adF[i] += adFadj[i];
  }
  dValidError = pDataCont->ComputeDeviance(&adF[0], true);

}


void CGBM::GBMTransferTreeToRList
(
 int *aiSplitVar,
 double *adSplitPoint,
 int *aiLeftNode,
 int *aiRightNode,
 int *aiMissingNode,
 double *adErrorReduction,
 double *adWeight,
 double *adPred,
 VEC_VEC_CATEGORIES &vecSplitCodes,
 int cCatSplitsOld
 )
{
	pTreeComp->TransferTreeToRList(*(pDataCont->getData()),
				 aiSplitVar,
				 adSplitPoint,
				 aiLeftNode,
				 aiRightNode,
				 aiMissingNode,
				 adErrorReduction,
				 adWeight,
				 adPred,
				 vecSplitCodes,
				 cCatSplitsOld);
}

double CGBM::InitF()
{
	return pDataCont->InitialFunctionEstimate();
}
