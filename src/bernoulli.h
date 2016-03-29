//------------------------------------------------------------------------------
//
//  File:       bernoulli.h
//
//  Description:   bernoulli distribution class used in GBM
//
//  History:    3/26/2001   gregr created
//              2/14/2003   gregr: adapted for R implementation
//
//------------------------------------------------------------------------------

#ifndef __bernoulli_h__
#define __bernoulli_h__

//------------------------------
// Includes
//------------------------------

#include "distribution.h"
#include "buildinfo.h"
#include <memory>

//------------------------------
// Class definition
//------------------------------
class CBernoulli : public CDistribution
{

public:
	//---------------------
	// Factory Function
	//---------------------
	static CDistribution* Create(SEXP radMisc,
											const char* szIRMeasure,
											int& cGroups, int& cTrain);

	//---------------------
	// Public destructor
	//---------------------
    virtual ~CBernoulli();

    //---------------------
    // Public Functions
    //---------------------
    void ComputeWorkingResponse(const CDataset* pData,
    		const double *adF,
				double *adZ);

    double Deviance(const CDataset* pData,
    				const double *adF,
                    bool isValidationSet=false);

    void InitF(const CDataset* pData,
    		double &dInitF,
	       unsigned long cLength);

    void FitBestConstant(const CDataset* pData,
    		const double *adF,
			 unsigned long cTermNodes,
				CTreeComps* pTreeComps);
    
    double BagImprovement(const CDataset& data,
    					  const double *adF,
    					  const bag& afInBag,
                          const CTreeComps* pTreeComps);

private:
    //----------------------
    // Private Constructors
    //----------------------
    CBernoulli(SEXP radMisc);

    //-------------------
    // Private Variables
    //-------------------
    vector<double> vecdNum;
    vector<double> vecdDen;
    bool fCappedPred;
};

#endif // BERNOULLI_H



