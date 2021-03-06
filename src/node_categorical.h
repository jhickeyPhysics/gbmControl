//------------------------------------------------------------------------------
//  GBM by Greg Ridgeway  Copyright (C) 2003
//
//  File:       node_categorical.h
//
//  License:    GNU GPL (version 2 or later)
//
//  Contents:   a node with a categorical split
//
//  Owner:      gregr@rand.org
//
//  History:    3/26/2001   gregr created
//              2/14/2003   gregr: adapted for R implementation
//
//------------------------------------------------------------------------------

#ifndef NODECATEGORICAL_H
#define NODECATEGORICAL_H

#include <float.h>
#include <algorithm>
#include "node_nonterminal.h"

class CNodeCategorical : public CNodeNonterminal
{
public:

 CNodeCategorical() : aiLeftCategory() {};
  ~CNodeCategorical();

  void PrintSubtree(unsigned long cIndent);
  void TransferTreeToRList(int &iNodeID,
			   const CDataset &data,
			   int *aiSplitVar,
			   double *adSplitPoint,
			   int *aiLeftNode,
			   int *aiRightNode,
			   int *aiMissingNode,
			   double *adErrorReduction,
			   double *adWeight,
			   double *adPred,
			   VEC_VEC_CATEGORIES &vecSplitCodes,
			   int cCatSplitsOld,
			   double dShrinkage);

  signed char WhichNode(const CDataset &data,
			unsigned long iObs);
  signed char WhichNode(double *adX,
			unsigned long cRow,
			unsigned long cCol,
			unsigned long iRow);

  void RecycleSelf(CNodeFactory *pNodeFactory);

  void reset() {
    CNodeNonterminal::reset();
    aiLeftCategory.resize(0);
  }

  std::vector<unsigned long> aiLeftCategory;
};

typedef CNodeCategorical *PCNodeCategorical;

#endif // NODECATEGORICAL_H



