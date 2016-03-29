//-----------------------------------
//
// File: tdist.cpp
//
// Description: t distribution implementation for GBM.
//
//-----------------------------------

//-----------------------------------
// Includes
//-----------------------------------
#include "locationm.h"
#include "tdist.h"
#include <vector>

//----------------------------------------
// Function Members - Private
//----------------------------------------
CTDist::CTDist(SEXP radMisc):
  CDistribution(radMisc), mpLocM("tdist", CDistribution::misc_ptr(true)[0])
{
	mdNu = CDistribution::misc_ptr(true)[0];
}


//----------------------------------------
// Function Members - Public
//----------------------------------------
CDistribution* CTDist::Create(SEXP radMisc,
									const char* szIRMeasure,
									int& cGroups, int& cTrain)
{
	return new CTDist(radMisc);
}

CTDist::~CTDist()
{

}

void CTDist::ComputeWorkingResponse
(
	const CDataset* pData,
    const double *adF,
    double *adZ
)
{
    unsigned long i = 0;
    double dU = 0.0;

    if(pData->offset_ptr(false) == NULL)
    {
        for(i=0; i<pData->get_trainSize(); i++)
        {
	  dU = pData->y_ptr()[i] - adF[i];
	  adZ[i] = (2 * dU) / (mdNu + (dU * dU));
        }
    }
    else
    {
        for(i=0; i<pData->get_trainSize(); i++)
        {
 		    dU = pData->y_ptr()[i] - pData->offset_ptr(false)[i] - adF[i];
			adZ[i] = (2 * dU) / (mdNu + (dU * dU));
        }
    }
}


void CTDist::InitF
(
	const CDataset* pData,
    double &dInitF,
    unsigned long cLength
)
{

	// Local variables
	int ii;

	// Get objects to pass into the LocM function
	int iN = int(cLength);
	std::vector<double> adArr(cLength);

	for (ii = 0; ii < iN; ii++)
	{
		double dOffset = (pData->offset_ptr(false)==NULL) ? 0.0 : pData->offset_ptr(false)[ii];
		adArr[ii] = pData->y_ptr()[ii] - dOffset;
	}

	dInitF = mpLocM.LocationM(iN, &adArr[0], pData->weight_ptr(), 0.5);
}

double CTDist::Deviance
(
	const CDataset* pData,
    const double *adF,
    bool isValidationSet
)
{
    unsigned long i=0;
    double dL = 0.0;
    double dW = 0.0;
	double dU = 0.0;

	// Switch to validation set if necessary
	long cLength = pData->get_trainSize();
	if(isValidationSet)
	{
	   pData->shift_to_validation();
	   cLength = pData->GetValidSize();
	}

    if(pData->offset_ptr(false) == NULL)
    {
        for(i=0; i<cLength; i++)
        {
			dU = pData->y_ptr()[i] - adF[i];
			dL += pData->weight_ptr()[i] * std::log(mdNu + (dU * dU));
            dW += pData->weight_ptr()[i];
        }
    }
    else
    {
        for(i=0; i<cLength; i++)
        {
			dU = pData->y_ptr()[i] - pData->offset_ptr(false)[i] - adF[i];
		    dL += pData->weight_ptr()[i] * std::log(mdNu + (dU * dU));
            dW += pData->weight_ptr()[i];
        }
    }

    // Switch back to training set if necessary
    if(isValidationSet)
    {
 	   pData->shift_to_train();
    }

    return dL/dW;
}


void CTDist::FitBestConstant
(
	const CDataset* pData,
    const double *adF,
    unsigned long cTermNodes,
    CTreeComps* pTreeComps
)
{
   	// Local variables
    unsigned long iNode = 0;
    unsigned long iObs = 0;


    std::vector<double> adArr, adW;
	// Call LocM for the array of values on each node
    for(iNode=0; iNode<cTermNodes; iNode++)
    {
      if(pTreeComps->GetTermNodes()[iNode]->cN >= pTreeComps->GetMinNodeObs())
        {
	  adArr.clear();
	  adW.clear();

	  for (iObs = 0; iObs < pData->get_trainSize(); iObs++)
	    {
	      if(pData->GetBagElem(iObs) && (pTreeComps->GetNodeAssign()[iObs] == iNode))
                {
		  const double dOffset = (pData->offset_ptr(false)==NULL) ? 0.0 : pData->offset_ptr(false)[iObs];
		  adArr.push_back(pData->y_ptr()[iObs] - dOffset - adF[iObs]);
		  adW.push_back(pData->weight_ptr()[iObs]);
                }
	    }

	  pTreeComps->GetTermNodes()[iNode]->dPrediction = mpLocM.LocationM(adArr.size(), &adArr[0],
							       &adW[0], 0.5);

        }
    }
}

double CTDist::BagImprovement
(
	const CDataset& data,
    const double *adF,
    const bag& afInBag,
    const CTreeComps* pTreeComps
)
{
    double dReturnValue = 0.0;
    unsigned long i = 0;
    double dW = 0.0;

    for(i=0; i<data.get_trainSize(); i++)
    {
        if(!data.GetBagElem(i))
        {
            const double dF = adF[i] + ((data.offset_ptr(false)==NULL) ? 0.0 : data.offset_ptr(false)[i]);
	    const double dU = (data.y_ptr()[i] - dF);
	    const double dV = (data.y_ptr()[i] - dF - pTreeComps->GetLambda() * pTreeComps->GetRespAdj()[i]) ;

            dReturnValue += data.weight_ptr()[i] * (std::log(mdNu + (dU * dU)) - log(mdNu + (dV * dV)));
            dW += data.weight_ptr()[i];
        }
    }

    return dReturnValue/dW;
}




