// Copyright (C) 2005 Alexander Ushakov
// Contents: Implementation of short braid form
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//

#include "ShortBraidForm.h"
#include "braid_group.h"
#include "DehornoyForm.h"
#include "ThRightNormalForm.h"


Word shortenBraid( int N , const Word& w )
{
  BraidGroup B( N );

  int n = -1;
  DehornoyForm DF( N , w );
  Word w3 = DF.getDehornoyForm( );
  Word w1 = ( w3.length( ) < w.length( ) ? w3 : w );
  Word result = w1;
  
  for( int i=0 ; i<3 ; ++i ) {
    DehornoyForm DF( N , B.twist( w1 ) );
    Word w2 = DF.getDehornoyForm( );
    if( w2.length( ) < result.length( ) ) {
      result = w2;
      n = i;
    }
    w1 = w2;
  }

  if( n==-1 )
    return w;

  if( n%2==0 )
    result = B.twist( result );

  return result;
}


Word shortBraidForm( int N , const Word& w )
{
  BraidGroup B( N );
  ThRightNormalForm NF( B , w );
  return shortenBraid( N , NF.getShortWord( ) );
}


vector< Word > shortBraidSbgpForm( int N , const vector< Word >& sbgp )
{
  vector< Word > result( sbgp.size() );
  for( int i=0 ; i<sbgp.size() ; ++i )
    result[i] = shortBraidForm( N , sbgp[i] );

  return result;
}
