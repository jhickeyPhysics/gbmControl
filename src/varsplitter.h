//------------------------------------------------------------------------------
//
//  File:       varsplitter.h
//
//  Description: header for class that splits a node on a particular variable.
//
//------------------------------------------------------------------------------

#ifndef __varsplitter_h__
#define __varsplitter_h__

//------------------------------
// Includes
//------------------------------
#include "node.h"
#include "nodeParameters.h"
#include <Rcpp.h>

//------------------------------
// Class Definition
//------------------------------
class VarSplitter
{
public:
	VarSplitter(unsigned long minNumObs);
	~VarSplitter();

	void SetForNode(CNode& nodeToSet);
	void SetForVariable(unsigned long iWhichVar, long cVarClasses);

	double GetBestImprovement() { return bestSplit.GetImprovement(); };
	void IncorporateObs(double dX,
			double dZ,
			double dW,
			long lMonotone);
	void EvaluateCategoricalSplit();
	NodeParams GetBestSplit() { return bestSplit;};


private:
	void WrapUpSplit();

	bool hasBestSplit, setForNode;
	unsigned long minObsInNode;
	double InitTotalWeight, InitWeightResiduals, dLastXValue;
	unsigned long InitNumObs;
	NodeParams bestSplit, proposedSplit;



};
#endif // __varplitter_h__
