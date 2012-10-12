// Copyright (C) 2005 Alexander Ushakov
// Contents: Implementation of class PowerWord
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//

#include "stdlib.h"
#include "PowerWord.h"
#include <strstream>
using namespace std;


PowerWord& PowerWord::freelyReduce( const_iterator beg , const_iterator end )
{
  list< PII >::const_iterator b = beg.theIterator;
  list< PII >::const_iterator e = end.theIterator;
  
  for( const_iterator cur_it=beg ; cur_it!=end ; ++cur_it ) {
    
    const_iterator cur_it2 = cur_it;
    if( ++cur_it2==this->end( ) )
      break;
  }
  
  return *this;
}


//=========================================================


PowerWord PowerWord::power( int t ) const
{
  if( t==0 ) return PowerWord( );
  
  PowerWord result = ( t<0 ? this->inverse() : *this );
  for( int i=1 ; i<abs( t ) ; ++i )
    result *= ( t<0 ? this->inverse() : *this );
  return result;
}


//=========================================================


PowerWord PowerWord::randomPowerWord( int gens , int wLen )
{
  if( wLen==0 )	return PowerWord( );
  
  int old = 0;
  list< int > result;
  for( int i=0 ; i<wLen ; ++i ) {
    int div = i==0 ? 2*gens : 2*gens-1;
    int g = ::rand()%div-gens;
    g = g>=0 ? g+1 : g;
    if( g+old==0 )
      g = gens;
    result.push_back( old = g );
  }
  
  return result;
}


//=========================================================
