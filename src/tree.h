//------------------------------------------------------------------------------
//  GBM by Greg Ridgeway  Copyright (C) 2003
//
//  File:       tree.h
//
//  License:    GNU GPL (version 2 or later)
//
//  Contents:   regression tree
//
//  Owner:      gregr@rand.org
//
//  History:    3/26/2001   gregr created
//              2/14/2003   gregr: adapted for R implementation
//
//------------------------------------------------------------------------------

#ifndef TREGBM_H
#define TREGBM_H

#include <cstdio>
#include <cfloat>
#include <algorithm>
#include <vector>
#include "dataset.h"
#include "node_search.h"
#include <ctime>


class CCARTTree
{
public:

    CCARTTree();
    ~CCARTTree();

    void Initialize();
    void grow(double *adZ,
	      const CDataset& data,
	      const double *adAlgW,
	      const double *adF,
	      unsigned long nBagged,
	      double dLambda,
	      unsigned long cMinObsInNode,
	      std::vector<unsigned long>& aiNodeAssign,
	      CNodeSearch *aNodeSearch);
    void Reset();

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
			     int cCatSplitsOld,
			     double dShrinkage);

    void PredictValid(const CDataset &pData,
		      unsigned long nValid,
		      double *adFadj);
    
    void Predict(double *adX,
		 unsigned long cRow,
		 unsigned long cCol,
		 unsigned long iRow,
		 double &dFadj);
    void Adjust(const std::vector<unsigned long>& aiNodeAssign,
		double *adFadj,
		unsigned long cMinObsInNode);
    
    void GetNodeCount(int &cNodes);
    void SetShrinkage(double dShrink)
    {
        this->dShrink = dShrink;
    }
    double GetShrinkage() {return dShrink;}
    long GetDepth(){return depthOfTree;}
    vector<CNode*> GetTermNodes(){return vecpTermNodes;}
    void SetDepth(long depth){ depthOfTree = depth;}
    void Print();
    void GetVarRelativeInfluence(double *adRelInf);


    double dError; // total squared error before carrying out the splits
private:
    void GetBestSplit(const CDataset &pData,
		      CNodeSearch *aNodeSearch,
		      std::vector<unsigned long>& aiNodeAssign,
		      double *adZ,
		      const double *adW);
    
    // Definition of a tree
    CNode* pRootNode;
    vector<CNode*> vecpTermNodes;
    long depthOfTree;
    double dShrink;


    // objects used repeatedly
    unsigned long cTerminalNodes;
    unsigned long cTotalNodeCount;

    unsigned long iBestNode;
    double dBestNodeImprovement;

    signed char schWhichNode;

};

#endif // TREGBM_H



