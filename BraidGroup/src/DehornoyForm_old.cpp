// Copyright (C) 2005 Alexander Ushakov
// Contents: Implementation of class DehornoyForm
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//

#include "DehornoyForm.h"


DehornoyForm::DehornoyForm( const int N , const Word& w ) :
  theLinkedStructure( N ),
  theIndex( N )
{
  theDehornoyForm = computeDehornoyForm( w );
}


Word DehornoyForm::computeDehornoyForm( const Word& w )
{
  theLinkedStructure.push_back( w.begin( ) , w.end( ) );
  theLinkedStructure.transformToDehornoyForm( );
  return theLinkedStructure.getWord( );
}
