// Copyright (C) 2007 Alexander Ushakov
// Contents: Definition of class FreeMetabelianGroup
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//


#ifndef _FreeMetabelianGroup_h_
#define _FreeMetabelianGroup_h_



#include "map"
#include "list"
#include "set"
using namespace std;

#include "tuples.h"
#include "Word.h"


//---------------------------------------------------------------------------//
//-------------------------- FreeMetabelianGroup ----------------------------//
//---------------------------------------------------------------------------//


//! Static class FreeMetabelianGroupAlgorithms encapsulates algorithms for FreeMetabelianGroup groups
/*!
  The class FreeMetabelianGroupAlgorithms is static, i.e., all member functions are static and there is
  no constructor defined. 
 */


class FreeMetabelianGroupAlgorithms
{

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Constructors:                                      //
  //                                                     //
  /////////////////////////////////////////////////////////

public:

  //! Default constructor is not instantiated to protect from creating the obects of this class
  FreeMetabelianGroupAlgorithms( );


  /////////////////////////////////////////////////////////
  //                                                     //
  //  Accessors:                                         //
  //                                                     //
  /////////////////////////////////////////////////////////

public:

  //! Decide if a word represents the identity of a free metabelian group
  /*!
    \f$N\f$ is the rank of the corresponding metabelian group,
    \f$w\f$ is the given word.
  */
  static bool trivial( int N , const Word& w );


  //! Decide if words represent conjugate elements in a free metabelian group
  /*!
    \f$N\f$ is the rank of the corresponding metabelian group,
    \f$w1, w2\f$ are the given word.
    The function returns a pair \f$(A,B)\f$ where \f$A\f$ is true if and only if
    words are conjugate and \f$B\f$ is a conjugator, i.e., \f$B^{-1} w1 B = w2\f$ 
    holds.
  */
  static pair< bool , Word > conjugate( int N , Word w1 , Word w2 );


  /////////////////////////////////////////////////////////
  //                                                     //
  //  Internal functions:                                //
  //                                                     //
  /////////////////////////////////////////////////////////

 private:
 public:

  //! Compute a word (not a shortest) defining the same edge map as EM
  /*!
    It is assumed here that EM defines an element of \f$M_r'\f$.
  */
  static Word getWordFromEdgeMap( int N , const map< vector< int > , int >& EM );


  /////////////////////////////////////////////////////////
  //                                                     //
  //  Data:                                              //
  //                                                     //
  /////////////////////////////////////////////////////////

private:


};


#endif
