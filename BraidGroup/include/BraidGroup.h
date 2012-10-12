// Copyright (C) 2005 Alexander Ushakov
// Contents: Definition of class BraidGroup
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//

#ifndef _BraidGroup_H_
#define _BraidGroup_H_

#include "Word.h"


//---------------------------------------------------------------------------//
//-------------------------------- BraidGroup -------------------------------//
//---------------------------------------------------------------------------//


//! Class BraidGroup (defines a representation of a Braid Group)//
/*!
  A BraidGroup is uniquely defined by its index.
*/


class BraidGroup
{
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Constructors                                       //
  //                                                     //
  /////////////////////////////////////////////////////////
 public:

  //! Constructor (creates the braid group on n strands).
  BraidGroup( int n ) : theRank(n) { }

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Accessors:                                         //
  //                                                     //
  /////////////////////////////////////////////////////////

 public:

  //! get the rank of the braid group (number of strands)
  int getRank( ) const { return theRank; }

  //! get the braid word called a half-twist (its square generates the center of the braid group)
  Word twist( const Word& w ) const;
  
  //! Get a list (non-symmetrized) of braid relators.
  static list< Word > getBraidRelators( int N );

  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Data members                                       //
  //                                                     //
  /////////////////////////////////////////////////////////
private:

  int theRank;
  // specifies the number of strands (number of generators + 1)

};

#endif
