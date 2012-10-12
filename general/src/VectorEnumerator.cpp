// Copyright (C) 2000 Alexander Ushakov
// Contents: Implementation of class VectorEnumerator
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//

#include "VectorEnumerator.h"

#include <iostream>

//---------------------------------------------------------------------------//
//--------------------------- VectorEnumerator ------------------------------//
//---------------------------------------------------------------------------//


vector< int > VectorEnumerator::operator* ( )
{
  if( seqComplete( ) )
    return vector< int >( curVector.begin( ) , curVector.begin( )+curLength );
  
  if( end( ) ) 
    return vector< int >();
  
  ++(*this);
  return curVector;
}


//=========================================================


VectorEnumerator& VectorEnumerator::operator++( )
{
  int cur;
  
  if( seqLimit( ) ) {
    curLength--;
    if( curLength>=0 ) {
      stepBack( curVector[curLength] );
      cur = next( curVector[curLength] );
    }
  }
  else
    cur = start( );
  
  while( curLength>=0 ) {
    
    if( !finish( cur ) ) {
      
      if( curLength>=curVector.capacity( ) )
	curVector.resize( curLength+5 );
      
      curVector[curLength] = cur;
      ++curLength;
      
      if( seqOK( ) ) {
	stepTo( );
	if( seqComplete( ) ) {
	  return *this;
	}
	else {
	  if( seqLimit( ) ) {
	    curLength--;
	    stepBack( curVector[curLength] );
	    cur = next( curVector[curLength] );
	  }
	  else {
	    cur = start( );
	  }
	}
      }
      else {
	cur = next( curVector[--curLength] );
      }
    } // if( !finish( cur ) ) {
    else{
      --curLength;
      if( curLength>=0 ) {
	stepBack( curVector[curLength] );
	cur = curVector[curLength];
      }
      cur = next( cur );
    } // if( !finish( cur ) ) { } else {
    
  } // while( curLength>=0 ) {
  
  return *this;
}


//=========================================================


ostream& VectorEnumerator::printSeq( ostream& os ) const
{
  os << "{ ";
  for( int t=0 ; t<curLength ; ++t ) {
    if( t ) 
      os << " , ";
    os << curVector[t];
    
  }
  os << " }";
  
  return os;
}


//=========================================================

