
#include "DehornoyForm.h"
#include "LinkedBraidStructure.h"


DehornoyForm::DehornoyForm( const int N , const Word& w ) :
  theIndex( N )
{
  theDehornoyForm = computeDehornoyForm( w );
}


Word DehornoyForm::computeDehornoyForm( const Word& w )
{
  // cout << "|w| = " << w.length( ) << endl;
  LinkedBraidStructure LBS( theIndex-1 );
  for( auto w_it=w.begin( ) ; w_it!=w.end( ) ; ++w_it )
    LBS.push_back( *w_it );
  LBS.removeLeftHandles( );
  // LBS.removeRightHandles( );
  return LBS.translateIntoWord( );
}


