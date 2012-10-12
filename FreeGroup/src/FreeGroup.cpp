// Copyright (C) 2005 Alexander Ushakov
// Contents: Implementation of class FreeGroup
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//

#include "FreeGroup.h"

//---------------------------------------------------------------------------//
//-------------------------------- FreeGroup --------------------------------//
//---------------------------------------------------------------------------//


FreeGroup::FreeGroup( int rank ) :
  theRank( rank ),
  theAlphabet(  rank )
{
  
}

FreeGroup::FreeGroup(  const FiniteAlphabet& a ) :
  theRank( a.size() ),
  theAlphabet( a )
{
  
}

bool FreeGroup::isPrimitive( const Word& w ) const
{
  return false;
}


bool FreeGroup::isAlmostPrimitive( const Word& w ) const
{
  return false;
}


bool FreeGroup::doesContain( const SubgroupFG& sbgp , const Word& w ) const
{
  return false;
}
