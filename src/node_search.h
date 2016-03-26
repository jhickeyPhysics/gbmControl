//------------------------------------------------------------------------------
//  GBM by Greg Ridgeway  Copyright (C) 2003
//
//  File:       node_search.h
//
//  License:    GNU GPL (version 2 or later)
//
//  Contents:   does the searching for where to split a node
//
//  Owner:      gregr@rand.org
//
//  History:    3/26/2001   gregr created
//              2/14/2003   gregr: adapted for R implementation
//
//------------------------------------------------------------------------------

#ifndef NODESEARCH_H
#define NODESEARCH_H

#include <vector>
#include "dataset.h"
#include "node.h"
#include "splitParameters.h"

using namespace std;

class CNodeSearch
{
public:

    CNodeSearch();
    ~CNodeSearch();
    void Initialize(unsigned long cMinObsInNode);

    void IncorporateObs(double dX,
			double dZ,
			double dW,
			long lMonotone);

    void Set(double dSumZ,
	     double dTotalW,
	     unsigned long cTotalN);
    void ResetForNewVar(unsigned long iWhichVar,
			long cVarClasses);
    
    static double Improvement
        (
            double dLeftW,
            double dRightW,
            double dMissingW,
            double dLeftSum,
            double dRightSum,
            double dMissingSum
        )
        {
            double dTemp = 0.0;
            double dResult = 0.0;

            if(dMissingW == 0.0)
            {
                dTemp = dLeftSum/dLeftW - dRightSum/dRightW;
                dResult = dLeftW*dRightW*dTemp*dTemp/(dLeftW+dRightW);
            }
            else
            {
                dTemp = dLeftSum/dLeftW - dRightSum/dRightW;
                dResult += dLeftW*dRightW*dTemp*dTemp;
                dTemp = dLeftSum/dLeftW - dMissingSum/dMissingW;
                dResult += dLeftW*dMissingW*dTemp*dTemp;
                dTemp = dRightSum/dRightW - dMissingSum/dMissingW;
                dResult += dRightW*dMissingW*dTemp*dTemp;
                dResult /= (dLeftW + dRightW + dMissingW);
            }

            return dResult;
        }

    double BestImprovement() { return bestSplit.ImprovedResiduals ; }
    void SetToSplit()
    {
        fIsSplit = true;
    };
    void SetupNewNodes(CNode& nodeToSplit);

    void EvaluateCategoricalSplit();
    void WrapUpCurrentVariable();

    double dInitTotalW;
    double dInitSumZ;
    unsigned long cInitN;
    SplitParams bestSplit;

private:
    // Split Parameters
    SplitParams proposedSplit;

    bool fIsSplit;

    unsigned long cMinObsInNode;
    double dLastXValue;

    std::vector<double> adGroupSumZ;
    std::vector<double> adGroupW;
    std::vector<unsigned long> acGroupN;
    std::vector<double> adGroupMean;

    // Has to be int for r_sort_index
    std::vector<int> aiCurrentCategory;
    std::vector<long> aiBestCategory;

};

#endif // NODESEARCH_H
