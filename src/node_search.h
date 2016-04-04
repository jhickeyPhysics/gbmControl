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
    double SplitAndCalcImprovement(vector<CNode*>& vecpTermNodes,
    					const CDataset& data,
    					vector<unsigned long>& aiNodeAssign);

    template <volatile bool xIsFactor>
    void IncorporateObs(double dX,
			double dZ,
			double dW,
			long lMonotone, SplitParams& proposedSplit, SplitParams& bestSplit);
    
    void Reset(const CDataset& data);


    // Remove InitTotalW if possible
    double dInitTotalW;
    double dInitSumZ;
    unsigned long cInitN;

private:
    //Private methods
    void ReAssignData(long splittedNodeIndex, vector<CNode*>& vecpTermNodes,
    					const CDataset& data, vector<unsigned long>& aiNodeAssign);
    void AssignToNode(CNode& terminalNode);
	void EvaluateCategoricalSplit(SplitParams& proposedSplit, SplitParams& bestSplit);
	void GenerateAllSplitsForVar(SplitParams& storeProposedSplit);
	void WrapUpProposedSplit(SplitParams& proposedSplit, SplitParams& bestSplit);
	void Set(CNode nodeToSplit);
	void ResetForNewVar(CNode nodeToSplit, unsigned long iWhichVar,
				long cVarClasses, SplitParams& proposedSplit, SplitParams& bestSplit);

    // Split Parameters -
    std::vector<SplitParams> proposedSplits;
    std::vector<SplitParams> bestSplits;
    SplitParams bestSplit;


    // Clean up if possible
    bool fIsSplit;
    long cTerminalNodes;
    unsigned long cMinObsInNode;
    double dLastXValue;

    // Move into proposed split
    std::vector<double> adGroupSumZ;
    std::vector<double> adGroupW;
    std::vector<unsigned long> acGroupN;
    std::vector<pair<double, int> > groupdMeanAndCategory;





};

#endif // NODESEARCH_H
