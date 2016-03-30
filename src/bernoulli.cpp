//-----------------------------------
//
// File: bernoulli.cpp
//
// Description: bernoulli distribution.
//
//-----------------------------------

//-----------------------------------
// Includes
//-----------------------------------

#include "bernoulli.h"
#include <memory>

//----------------------------------------
// Function Members - Private
//----------------------------------------
CBernoulli::CBernoulli(SEXP radMisc): CDistribution(radMisc)
{
  // Used to issue warnings to user that at least one terminal node capped
  fCappedPred = false;
}

//----------------------------------------
// Function Members - Public
//----------------------------------------
CDistribution* CBernoulli::Create(SEXP radMisc,
											const char* szIRMeasure,
											int& cGroups, int& cTrain)
{
	return new CBernoulli(radMisc);
}

CBernoulli::~CBernoulli()
{
}

void CBernoulli::ComputeWorkingResponse
(
	const CDataset* pData,
    const double *adF,
    double *adZ
)
{
  unsigned long i = 0;
  double dProb = 0.0;
  double dF = 0.0;

  for(i=0; i<pData->get_trainSize(); i++)
  {
    dF = adF[i] + ((pData->offset_ptr(false)==NULL) ? 0.0 : pData->offset_ptr(false)[i]);
    dProb = 1.0/(1.0+std::exp(-dF));

    adZ[i] = pData->y_ptr()[i] - dProb;
#ifdef NOISY_DEBUG
//  Rprintf("dF=%f, dProb=%f, adZ=%f, pData->y_ptr()=%f\n", dF, dProb, adZ[i], pData->y_ptr()[i]);
    if(dProb<  0.0001) Rprintf("Small prob(i=%d)=%f Z=%f\n",i,dProb,adZ[i]);
    if(dProb>1-0.0001) Rprintf("Large prob(i=%d)=%f Z=%f\n",i,dProb,adZ[i]);
#endif
  }
}


void CBernoulli::InitF
(
	const CDataset* pData,
    double &dInitF,
    unsigned long cLength
)
{
    unsigned long i=0;
    double dTemp=0.0;

    if(pData->offset_ptr(false)==NULL)
    {
        double dSum=0.0;
        for(i=0; i<cLength; i++)
        {
            dSum += pData->weight_ptr()[i]*pData->y_ptr()[i];
            dTemp += pData->weight_ptr()[i];
        }
        dInitF = std::log(dSum/(dTemp-dSum));
    }
    else
    {
        // Newton method for solving for F
        // should take about 3-6 iterations.
        double dNum=0.0;         // numerator
        double dDen=0.0;         // denominator
        double dNewtonStep=1.0;  // change
        dInitF = 0.0;
        while(dNewtonStep > 0.0001)
        {
            dNum=0.0;
            dDen=0.0;
            for(i=0; i<cLength; i++)
            {
                dTemp = 1.0/(1.0+std::exp(-(pData->offset_ptr(false)[i] + dInitF)));
                dNum += pData->weight_ptr()[i]*(pData->y_ptr()[i]-dTemp);
                dDen += pData->weight_ptr()[i]*dTemp*(1.0-dTemp);
            }
            dNewtonStep = dNum/dDen;
            dInitF += dNewtonStep;
        }
    }
}

double CBernoulli::Deviance
(
	const CDataset* pData,
    const double *adF,
    bool isValidationSet
)
{
   unsigned long i=0;
   double dL = 0.0;
   double dF = 0.0;
   double dW = 0.0;

   // Switch to validation set if necessary
   long cLength = pData->get_trainSize();
   if(isValidationSet)
   {
	   pData->shift_to_validation();
	   cLength = pData->GetValidSize();
   }

   if(pData->offset_ptr(false)==NULL)
   {
      for(i=0; i!=cLength; i++)
      {
         dL += pData->weight_ptr()[i]*(pData->y_ptr()[i]*adF[i] - std::log(1.0+std::exp(adF[i])));
         dW += pData->weight_ptr()[i];
      }
   }
   else
   {
      for(i=0; i!=cLength; i++)
      {
         dF = adF[i] + pData->offset_ptr(false)[i];
         dL += pData->weight_ptr()[i]*(pData->y_ptr()[i]*dF - std::log(1.0+std::exp(dF)));
         dW += pData->weight_ptr()[i];
      }
   }
   // Switch back to trainig set if necessary
   if(isValidationSet)
   {
	   pData->shift_to_train();
   }

   return -2*dL/dW;
}


void CBernoulli::FitBestConstant
(
  const CDataset* pData,
  const double *adF,
  unsigned long cTermNodes,
  double* adZ,
  CTreeComps* pTreeComps
)
{
  unsigned long iObs = 0;
  unsigned long iNode = 0;
  double dTemp = 0.0;
  
  vecdNum.resize(cTermNodes);
  vecdNum.assign(vecdNum.size(),0.0);
  vecdDen.resize(cTermNodes);
  vecdDen.assign(vecdDen.size(),0.0);

  for(iObs=0; iObs<pData->get_trainSize(); iObs++)
  {
    if(pData->GetBagElem(iObs))
    {
      vecdNum[pTreeComps->GetNodeAssign()[iObs]] += pData->weight_ptr()[iObs]*adZ[iObs];
      vecdDen[pTreeComps->GetNodeAssign()[iObs]] +=
          pData->weight_ptr()[iObs]*(pData->y_ptr()[iObs]-adZ[iObs])*(1-pData->y_ptr()[iObs]+adZ[iObs]);
#ifdef NOISY_DEBUG
/*
      Rprintf("iNode=%d, dNum(%d)=%f, dDen(%d)=%f\n",
              aiNodeAssign[iObs],
              iObs,vecdNum[aiNodeAssign[iObs]],
              iObs,vecdDen[aiNodeAssign[iObs]]);
*/
#endif
    }
  }

  for(iNode=0; iNode<cTermNodes; iNode++)
  {
    if(pTreeComps->GetTermNodes()[iNode]!=NULL)
    {
      if(vecdDen[iNode] == 0)
      {
          pTreeComps->GetTermNodes()[iNode]->dPrediction = 0.0;
      }
      else
      {
        dTemp = vecdNum[iNode]/vecdDen[iNode];
        // avoid large changes in predictions on log odds scale
        if(std::abs(dTemp) > 1.0)
        {
          if(!fCappedPred)
          {
            // set fCappedPred=true so that warning only issued once
            fCappedPred = true;  
            Rcpp::warning("Some terminal node predictions were excessively large for Bernoulli and have been capped at 1.0. Likely due to a feature that separates the 0/1 outcomes. Consider reducing shrinkage parameter.");
          }
          if(dTemp>1.0) dTemp = 1.0;
          else if(dTemp<-1.0) dTemp = -1.0;
        }
        pTreeComps->GetTermNodes()[iNode]->dPrediction = dTemp;
      }
    }
  }
}


double CBernoulli::BagImprovement
(
	const CDataset& data,
    const double *adF,
    const bag& afInBag,
    const double shrinkage,
    const double* adFadj
)
{
    double dReturnValue = 0.0;
    double dF = 0.0;
    double dW = 0.0;
    unsigned long i = 0;

    for(i=0; i<data.get_trainSize(); i++)
    {
        if(!data.GetBagElem(i))
        {
            dF = adF[i] + ((data.offset_ptr(false)==NULL) ? 0.0 : data.offset_ptr(false)[i]);

            if(data.y_ptr()[i]==1.0)
            {
                dReturnValue += data.weight_ptr()[i]*shrinkage*adFadj[i];
            }
            dReturnValue += data.weight_ptr()[i]*
                            (std::log(1.0+std::exp(dF)) -
                             std::log(1.0+std::exp(dF+shrinkage*adFadj[i])));
            dW += data.weight_ptr()[i];
        }
    }

    return dReturnValue/dW;
}
