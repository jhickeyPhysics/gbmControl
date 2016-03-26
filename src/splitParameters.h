//------------------------------------------------------------------------------
//
//  File:       splitParameters.h
//
//  Description:   struct that contains the parameters used to split a node
//
//------------------------------------------------------------------------------

#ifndef __splitParameters_h__
#define __splitParameters_h__



typedef std::vector<int> bag;

//------------------------------
// Struct Definition
//------------------------------
struct SplitParams
{

	// Left Node Definition
	double LeftWeightResiduals;
	double LeftTotalWeight;
	long LeftNumObs;

	// Right Node Definition
	double RightWeightResiduals;
	double RightTotalWeight;
	long RightNumObs;

	// Missing Node Definition
	double MissingWeightResiduals;
	double MissingTotalWeight;
	long MissingNumObs;

	// Splitting values
	double SplitValue; // Continuous Split Value
	long SplitVar; // Which feature to split on
	long SplitClass; // Categorical Split Value
	double ImprovedResiduals;

};

#endif // __splitParameters_h__
