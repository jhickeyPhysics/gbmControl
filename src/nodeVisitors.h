//------------------------------------------------------------------------------
//
//  File:       nodeVisitors.h
//
//  Description:  header  contains the parameters used to split a node
//
//	Author: 	James Hickey
//------------------------------------------------------------------------------

#ifndef __nodeVisitors_h__
#define __nodeVisitors_h__

//------------------------------
// Includes
//------------------------------
#include "dataset.h"
#include "node.h"
#include <Rcpp.h>

//------------------------------
// Generic Visitor Definition
//------------------------------
template<typename DerivedNodeVisitor>
class GenericNodeVisitor
{
public:

	//---------------------
	// Public Functions
	//---------------------
	void printSubTree(CNode* node, unsigned long cIndent)
	{
		dispatchNodeToPrint(node, cIndent);
	}


private:

	DerivedNodeVisitor& derivedVisitor()
	{
		return *static_cast<DerivedNodeVisitor*>(this);
	}


	void dispatchNodeToPrint(CNode* node, unsigned long cIndent)
	{
		switch (node->splitType)
		{
			case CNode::categorical:
				derivedVisitor().printCategorical(node, cIndent);
				break;

			case CNode::continuous:
				derivedVisitor().printContinuous(node, cIndent);
				break;

			case CNode::none:
				derivedVisitor().printTerminal(node, cIndent);
		}
	}

};

//------------------------------
// Categorical Split Visitor Def.
//------------------------------
class CategoricalSplitVisitor: public GenericNodeVisitor<CategoricalSplitVisitor>
{
public:
	void printCategorical(CNode* node, unsigned long cIndent)
	{

		const std::size_t cLeftCategory = node->aiLeftCategory.size();

		for(long i=0; i< cIndent; i++) Rprintf("  ");
		Rprintf("N=%f, Improvement=%f, Prediction=%f, NA pred=%f\n",
		  node->dTrainW,
		  node->dImprovement,
		  node->dPrediction,
		  (node->pMissingNode == NULL ? 0.0 : node->pMissingNode->dPrediction));

		for(long i=0; i< cIndent; i++) Rprintf("  ");
		Rprintf("V%d in ",node->iSplitVar);
		for(long i=0; i<cLeftCategory; i++)
		  {
			Rprintf("%d", node->aiLeftCategory[i]);
			if(i<cLeftCategory-1) Rprintf(",");
		  }
		Rprintf("\n");
		node->pLeftNode->PrintSubtree((cIndent+1));

		for(long i=0; i< cIndent; i++) Rprintf("  ");
		Rprintf("V%d not in ", node->iSplitVar);
		for(long i=0; i<cLeftCategory; i++)
		  {
			Rprintf("%d", node->aiLeftCategory[i]);
			if(i<cLeftCategory-1) Rprintf(",");
		  }
		Rprintf("\n");
		node->pRightNode->PrintSubtree(cIndent+1);

		for(long i=0; i< cIndent; i++) Rprintf("  ");
		Rprintf("missing\n");
		node->pMissingNode->PrintSubtree(cIndent+1);
	}
};

//------------------------------
// Continuous Split Visitor Def.
//------------------------------
class ContinuousSplitVisitor: public GenericNodeVisitor<ContinuousSplitVisitor>
{
public:
	void printContinuous(CNode* node, unsigned long cIndent)
	{
		for(long i=0; i< cIndent; i++) Rprintf("  ");

		Rprintf("N=%f, Improvement=%f, Prediction=%f, NA pred=%f\n",
		  node->dTrainW,
		  node->dImprovement,
		  node->dPrediction,
		  (node->pMissingNode == NULL ? 0.0 : node->pMissingNode->dPrediction));

		for(long i=0; i< cIndent; i++) Rprintf("  ");

		Rprintf("V%d < %f\n",
		  node->iSplitVar,
		  node->dSplitValue);
		node->pLeftNode->PrintSubtree(cIndent+1);

		for(long i=0; i< cIndent; i++) Rprintf("  ");
		Rprintf("V%d > %f\n",
		  node->iSplitVar,
		  node->dSplitValue);
		node->pRightNode->PrintSubtree(cIndent+1);

		for(long i=0; i< cIndent; i++) Rprintf("  ");
		Rprintf("missing\n");
		node->pMissingNode->PrintSubtree(cIndent+1);
	}

};

//------------------------------
// NoSplit Visitor Def.
//------------------------------
class NoSplitVisitor: public GenericNodeVisitor<NoSplitVisitor>
{
public:
	void printTerminal(CNode* node, unsigned long cIndent)
	{
		  for(long i=0; i< cIndent; i++) Rprintf("  ");
		  Rprintf("N=%f, Prediction=%f *\n",
			  node->dTrainW,
			  node->dPrediction);
	}
};
