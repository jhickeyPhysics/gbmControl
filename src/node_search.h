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
    void GenerateAllSplits(vector<CNode*>& vecpTermNodes, const CDataset& data,
    						double* residuals, vector<unsigned long>& aiNodeAssign);
    double SplitAndCalcImprovement(vector<CNode*>& vecpTermNodes, const CDataset& data,
    					vector<unsigned long>& aiNodeAssign);

    void IncorporateObs(double dX,
			double dZ,
			double dW,
			long lMonotone);

    void Set(CNode nodeToSplit);
    void ResetForNewVar(unsigned long iWhichVar,
			long cVarClasses);
    
    void Reset();

    void SetupNewNodes(CNode& nodeToSplit);
    void AssignToNode(CNode& terminalNode);
    void EvaluateCategoricalSplit();
    void WrapUpCurrentVariable();

    double dInitTotalW;
    double dInitSumZ;
    unsigned long cInitN;
    SplitParams bestSplit;

private:
    //Private methods
    void ReAssignData(long splittedNodeIndex, vector<CNode*>& vecpTermNodes,
    					const CDataset& data, vector<unsigned long>& aiNodeAssign);

    // Split Parameters
    SplitParams proposedSplit;

    bool fIsSplit;
    long cTerminalNodes;
    unsigned long cMinObsInNode;
    double dLastXValue;

    // Move into proposed split
    std::vector<double> adGroupSumZ;
    std::vector<double> adGroupW;
    std::vector<unsigned long> acGroupN;
    std::vector<double> adGroupMean;

    // Has to be int for r_sort_index
    // Make both int? - Can we remove both of these?!
    std::vector<int> aiCurrentCategory;
    std::vector<int> aiBestCategory;

};

#endif // NODESEARCH_H
