//------------------------------------------------------------------------------
//
//  File:       genericNodeDispatch.h
//
//  Description: generic dispatcher for all split types - defines the interface.
//
//	Author: 	James Hickey
//------------------------------------------------------------------------------

#ifndef __genericNodeDispatch_h__
#define __genericNodeDispatch_h__

//------------------------------
// Includes
//------------------------------
#include "categoricalSplitDispatch.h"
#include "continuousSplitDispatch.h"
#include "dataset.h"
#include "node.h"
#include "noSplitDispatch.h"
#include <Rcpp.h>


//------------------------------
// Generic Dispatch Definition
//------------------------------
class GenericNodeDispatch
{
public:
	//----------------------
	// Public Constructors
	//----------------------
	GenericNodeDispatch(){};

	//---------------------
	// Public destructor
	//---------------------
	virtual ~GenericNodeDispatch(){};

	//---------------------
	// Public Functions
	//---------------------
	void printSubTree(CNode* node, unsigned long cIndent)
	{
		dispatchNodeToPrint(node, cIndent);
	}


private:


	//---------------------
	// Private Functions
	//---------------------
	// Dispatch to correct delegate class

	void dispatchNodeToPrint(CNode* node, unsigned long cIndent)
	{
		switch (node->splitType)
		{
			case categorical:
				catDispatch.printCategorical(node, cIndent);
				break;

			case continuous:
				contDispatch.printContinuous(node, cIndent);
				break;

			case none:
				termDispatch.printTerminal(node, cIndent);
				break;

			default:
				throw GBM::failure("Can't dispatch node of unknown type - print fnc.");
		}
	}

	//---------------------
	// Private Variables
	//---------------------
	CategoricalSplitDispatch catDispatch;
	ContinuousSplitDispatch contDispatch;
	NoSplitDispatch termDispatch;
};


#endif // __genericNodeDispatch_h__
