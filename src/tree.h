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

    void grow(double *adZ,
	      const CDataset& data,
	      const double *adF,
	      unsigned long cMinObsInNode,
	      std::vector<unsigned long>& aiNodeAssign,
	      CNodeSearch& aNodeSearch);
    void Reset();

    CNode* GetRootNode();
    const CNode* GetRootNode() const;

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
    long GetDepth(){return depthOfTree;}
    vector<CNode*> GetTermNodes(){return vecpTermNodes;}
    void SetDepth(long depth){ depthOfTree = depth;}
    void Print();

private:
    
    // Definition of a tree
    CNode* pRootNode;
    vector<CNode*> vecpTermNodes;
    long depthOfTree;
    double dShrink;
    double dError; // total squared error before carrying out the splits
    unsigned long cTotalNodeCount;

};

#endif // TREGBM_H



