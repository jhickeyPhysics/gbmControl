//------------------------------------------------------------------------------
//
//  File:       gbmDataContainer.h
//
//  Description:   Header file for a class that stores the data and dist for gbm.
//
//------------------------------------------------------------------------------


#ifndef __gbmDataContainer_h__
#define __gbmDataContainer_h__

//------------------------------
// Includes
//------------------------------
#include "buildinfo.h"
#include "dataset.h"
#include "distribution.h"
#include "distributionFactory.h"
#include "gbmTreeComps.h"
#include <Rcpp.h>
#include <vector>
#include <memory>

//------------------------------
// Class definition
//------------------------------
class CGBMDataContainer
{
public:
	//----------------------
	// Public Constructors
	//----------------------
    CGBMDataContainer(SEXP radY, SEXP radOffset, SEXP radX, SEXP raiXOrder,
            SEXP radWeight, SEXP racVarClasses,
            SEXP ralMonotoneVar, SEXP radMisc, const std::string& family,
            int cTrain, int cFeatures, int& cGroups, double bagFraction);

	//---------------------
	// Public destructor
	//---------------------
    ~CGBMDataContainer();

    //---------------------
	// Public Functions
	//---------------------
    void Initialize();
    void InitializeFunctionEstimate(double &dInitF, unsigned long cLength);
    void ComputeResiduals(const double* adF, double* adZ);
    void ComputeBestTermNodePreds(const double* adF, double* adZ, CTreeComps* pTreeComp, int& cNodes);
    double ComputeDeviance(const double *adF, bool isValidationSet=false);
    double ComputeBagImprovement(const double* adF, const double shrinkage, const double* adFadj);
    void BagData();
    CDistribution* getDist();
    CDataset* getData();


private:
	//-------------------
	// Private Variables
	//-------------------
    CDataset data;
    CDistribution* pDist;
    DistributionFactory* DistFactory; // currently a singleton - does not need to be now - TODO:remove property.

};

#endif //  __gbmDataContainer_h__
