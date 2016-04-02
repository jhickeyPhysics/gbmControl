//------------------------------------------------------------------------------
//
//  File:       continuousSplitDispatch.h
//
//  Description: dispatch for no split terminal nodes.
//
//	Author: 	James Hickey
//------------------------------------------------------------------------------

#ifndef __noSplitDispatch_h__
#define __noSplitDispatch_h__

//------------------------------
// Includes
//------------------------------
#include "dataset.h"
#include "node.h"
#include <Rcpp.h>


//------------------------------
// Class Definition
//------------------------------
class NoSplitDispatch
{
public:
	//----------------------
	// Public Constructors
	//----------------------
	NoSplitDispatch(){};

	//---------------------
	// Public destructor
	//---------------------
	~NoSplitDispatch(){};

	//---------------------
	// Public Functions
	//---------------------
	void printTerminal(CNode* node, unsigned long cIndent)
	{
		  for(long i=0; i< cIndent; i++) Rprintf("  ");
		  Rprintf("N=%f, Prediction=%f *\n",
			  node->dTrainW,
			  node->dPrediction);
	}
};
#endif //__noSplitDispatch_h__
