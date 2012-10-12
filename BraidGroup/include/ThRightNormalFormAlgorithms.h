// Copyright (C) 2005 Alexander Ushakov
// Contents: Definition of class ThRightNormalForm
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//


#ifndef _ThRightNormalFormAlgorithms_h_
#define _ThRightNormalFormAlgorithms_h_



#include "set"
using namespace std;


#include "ThRightNormalForm.h"



//! Compute summit set representatives for a tuple of braids
pair< vector< ThRightNormalForm > , ThRightNormalForm > getSummitSetRepresentative( int rank , const vector< ThRightNormalForm >& elts );


//! Compute a simple conjugator for "tuple" starting from "start". 
/*!
  Algorithm from Gonzalez-Menese, "Improving an algorithm to solve Multiple Simultaneous Conjugacy problem for braid groups".
  Conjugating by a simple element does not decrease the infimum of the tuple.
*/
Permutation getSimpleConjugator( int rank , const vector< ThRightNormalForm >& tuple , const Permutation& start );


set< Permutation > getSimpleConjugators( int rank , const vector< ThRightNormalForm >& tuple );


#endif
