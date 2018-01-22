// Copyright (C) 2005 Alexander Ushakov
// Contents: Implementation of class BraidGroup
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//

#include "BraidGroup.h"


//---------------------------------------------------------------------------//
//-------------------------------- BraidGroup -------------------------------//
//---------------------------------------------------------------------------//


Word BraidGroup::twist( const Word& w ) const
{
  list< int > result;

	Word::const_iterator w_it = w.begin( );
  for( ; w_it!=w.end( ) ; ++w_it) {
    int g = *w_it;
    result.push_back( g<0 ? -g-theRank : theRank-g );
  }

  return Word(std::move(result));
}
