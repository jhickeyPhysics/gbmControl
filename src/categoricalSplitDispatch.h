//------------------------------------------------------------------------------
//
//  File:       categoricalSplitDispatch.h
//
//  Description: dispatch for categorical splits.
//
//	Author: 	James Hickey
//------------------------------------------------------------------------------

#ifndef __categoricalSplitDispatch_h__
#define __categoricalSplitDispatch_h__

//------------------------------
// Includes
//------------------------------
#include "dataset.h"
#include "node.h"
#include <Rcpp.h>


//------------------------------
// Class Definition
//------------------------------
class CategoricalSplitDispatch
{
public:
	//----------------------
	// Public Constructors
	//----------------------
	CategoricalSplitDispatch(){};

	//---------------------
	// Public destructor
	//---------------------
	~CategoricalSplitDispatch(){};

	//---------------------
	// Public Functions
	//---------------------
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
#endif //__categoricalSplitDispatch_h__
