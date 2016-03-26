//------------------------------------------------------------------------------
//
//  File:       splitParameters.h
//
//  Description:  header  contains the parameters used to split a node
//
//------------------------------------------------------------------------------

#ifndef __splitParameters_h__
#define __splitParameters_h__

//------------------------------
// Class Definition
//------------------------------
class SplitParams
{
public:
	//----------------------
	// Public Constructors
	//----------------------
	SplitParams();

	//---------------------
	// Public destructor
	//---------------------
    ~SplitParams();

	//---------------------
	// Public Functions
	//---------------------
	void UpdateMissingNode(double predIncrement, double trainWIncrement, long numIncrement = 1);
	void UpdateLeftNode(double predIncrement, double trainWIncrement, long numIncrement = 1);
	bool SplitIsCorrMonotonic(long specifyMonotone);
	void NodeGradResiduals();
	bool HasMinNumOfObs(long minObsInNode);

	SplitParams& operator=(const SplitParams rhs)
	{
		RightWeightResiduals = rhs.RightWeightResiduals;
		RightTotalWeight = rhs.RightTotalWeight;
		RightNumObs = rhs.RightNumObs;

		LeftWeightResiduals = rhs.LeftWeightResiduals;
		LeftTotalWeight = rhs.LeftTotalWeight;
		LeftNumObs = rhs.LeftNumObs;

		MissingWeightResiduals = rhs.MissingWeightResiduals;
		MissingTotalWeight = rhs.MissingTotalWeight;
		MissingNumObs = rhs.MissingNumObs;

		SplitValue = rhs.SplitValue;
		SplitVar = rhs.SplitVar;
		SplitClass = rhs.SplitClass;
		ImprovedResiduals = rhs.ImprovedResiduals;


	}
	//---------------------
	// Public Variables
	//---------------------
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
