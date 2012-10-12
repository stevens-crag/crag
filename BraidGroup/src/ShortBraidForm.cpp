
#include "ShortBraidForm.h"
#include "LinkedBraidStructure.h"
#include "BraidGroup.h"
#include "DehornoyForm.h"
#include "ThRightNormalForm.h"


LinkedBraidStructure shortenLBS( LinkedBraidStructure& lbs )
{
  LinkedBraidStructure result = lbs;
  for( int i=0 ; i<4 ; ++i ) {
    if( i%2==0 )
      lbs.removeRightHandles( );
    else
      lbs.removeLeftHandles( );
    if( result.size( )>lbs.size( ) )
      result = lbs;
  }
  
  return result;
}


Word shortenBraid( int N , const Word& w )
{
  LinkedBraidStructure df( N-1 , w );
  LinkedBraidStructure result = df;
  for( int i=0 ; i<4 ; ++i ) {
    if( i%2==0 )
      df.removeRightHandles( );
    else
      df.removeLeftHandles( );
    if( result.size( )>df.size( ) )
      result = df;
  }  

  return result.translateIntoWord( );
}

Word shortBraidForm( int N , const Word& w )
{
  BraidGroup B( N );
  ThRightNormalForm NF( B , w );
  Word w1 = NF.getShortWord( );
  //  cout << "int len = " << w1.length() << endl;
  return shortenBraid( N , w1 );
}

vector< Word > shortBraidSbgpForm( int N , const vector< Word >& sbgp )
{
  vector< Word > result( sbgp.size() );
  for( int i=0 ; i<sbgp.size() ; ++i )
    result[i] = shortBraidForm( N , sbgp[i] );

  return result;
}
