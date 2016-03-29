//------------------------------------------------------------------------------
//
//  File:       gbmTreeComps.h
//
//  Description:   Header file for a class containing the tree components used by
//    the gbm engine.
//
//------------------------------------------------------------------------------

#ifndef __gbmTreeComps_h__
#define __gbmTreeComps_h__
//------------------------------
// Includes
//------------------------------
#include "buildinfo.h"
#include "tree.h"
#include "dataset.h"
#include <vector>
#include <memory>

//------------------------------
// Class definition
//------------------------------
class CTreeComps
{
public:
	//----------------------
	// Public Constructors
	//----------------------
    CTreeComps(double dLambda,
    	    unsigned long cDepth,
    	    unsigned long cMinObsInNode);


	//---------------------
	// Public destructor
	//---------------------
    ~CTreeComps();

    //---------------------
	// Public Functions
	//---------------------
    void TreeInitialize(const CDataset* pData);
    void GrowTrees(const CDataset* pData, int& cNodes);
    void AdjustAndShrink();
    void PredictValid(const CDataset* pData);
    void TransferTreeToRList(const CDataset &pData,
		     int *aiSplitVar,
		     double *adSplitPoint,
		     int *aiLeftNode,
		     int *aiRightNode,
		     int *aiMissingNode,
		     double *adErrorReduction,
		     double *adWeight,
		     double *adPred,
		     VEC_VEC_CATEGORIES &vecSplitCodes,
		     int cCatSplitsOld);

    // getters
	std::vector<unsigned long> GetNodeAssign();
	vector<CNode*> GetTermNodes();

	double* GetGrad();
	double* GetRespAdj();
	const double* GetRespAdj() const;
	const double  RespAdjElem(int ind);

	double GetLambda();
	const double GetLambda() const;
	unsigned long GetMinNodeObs();

private:
	//-------------------
	// Private Variables
	//-------------------

    // these objects are for the tree growing
    // allocate them once here for all trees to use
    std::vector<unsigned long> aiNodeAssign;
    std::vector<CNodeSearch> aNodeSearch;
    std::auto_ptr<CCARTTree> ptreeTemp;

    std::vector<double> adZ;
    std::vector<double> adFadj;

    double dLambda;
    unsigned long cMinObsInNode;
};

#endif //  __gbmTreeComps_h__
