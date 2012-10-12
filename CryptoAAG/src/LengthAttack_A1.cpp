// Copyright (C) 2005 Alexander Ushakov
// Contents: Implementation of length attack
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//

// A1

#include <set>
#include "LengthAttack.h"
#include "ShortBraidForm.h"
#include <time.h>


int LengthAttack_A1::sbgpGeneratorsWeight( const vector< Word >& A )
{
  int result = 0;
  for( int i=0 ; i<A.size() ; ++i )
    result += A[i].length( );
  return result;
}

void LengthAttack_A1::addNewElt( const vector< Word >& A , set< ELT >& checkedElements , set< ELT >& uncheckedElements )
{
  int weight = sbgpGeneratorsWeight( A );
  ELT new_elt( weight , A );

  if(   checkedElements.find( new_elt )!=  checkedElements.end( ) ) return;
  if( uncheckedElements.find( new_elt )!=uncheckedElements.end( ) ) return;

  uncheckedElements.insert( new_elt );
}


void LengthAttack_A1::tryElt( int N , const ELT& cur , const vector< Word >& B , set< ELT >& checkedElements , set< ELT >& uncheckedElements )
{
  // we better vary this value, depending on parameters of A and B
  int MAX_DELTA = 80;
  int maxDecrease = 0;
  vector< Word > maxTuple;

  for( int i=0 ; i<B.size( ) ; ++i ) {
    
    Word b = B[i];
    for( int d=0 ; d<2 ; ++d ) {
    
      //cout << "   #" << i+1 << "," << d+1;
      if( d==1 ) b = -b;
      bool candidate = true;
      int delta = 0;
      vector< Word > A = cur.second;
      for( int t=0 ; t<A.size( ) && candidate ; ++t ) {
	delta -= A[t].length( );
	A[t] = shortenBraid( N , -b*A[t]*b );
	delta += A[t].length( );

	if( delta>MAX_DELTA && t>=4 ) {
	  candidate = false;
	  //cout << " stopped @ " << t << endl;
	}
      }
      if ( -delta > maxDecrease ){
	maxDecrease = -delta;
	maxTuple = A;
      }
      /*      
	      if( candidate ) {
	      addNewElt( A , checkedElements , uncheckedElements );
	      //cout << " accepted with " << delta << endl;
	      }
      */
    }
  }
  
  if ( maxDecrease > 0 ){
    addNewElt( maxTuple , checkedElements , uncheckedElements );
  }
  
}


bool LengthAttack_A1::check_ifVectorsEqual( int N , const vector< Word >& A1 , const vector< Word >& A2 )
{
  int len = A1.size( );
  for( int i=0 ; i<len ; ++i ) {
    Word s = shortenBraid( N , A1[i]*-A2[i]  );
    if( s.length( )>0 )
      return false;
  }

  return true;
}


findKey_LengthBasedResult LengthAttack_A1::findKey_LengthBased( int N , const vector< Word >& A1 , const vector< Word >& A2 , const vector< Word >& B , int sec, ostream& out )
{
  set< ELT > checkedElements;
  set< ELT > uncheckedElements;

  int init_time = time( 0 );

  int init_weight1 = sbgpGeneratorsWeight( A1 );
  int init_weight2 = sbgpGeneratorsWeight( A2 );
  ELT init( init_weight2 , A2 );
  uncheckedElements.insert( init );
  out << "Initial weights: " << init_weight1 << ", " << init_weight2 << endl;
  int best_result = 999999;

  for( int c=0 ; uncheckedElements.size( ) && c<100000 ; ++c ) {

    out << "Elts to try: " << uncheckedElements.size() << endl;

    ELT cur = *uncheckedElements.begin( );
    uncheckedElements.erase( uncheckedElements.begin( ) );
    checkedElements.insert( cur );

    int cur_time = time( 0 );
    if( best_result>cur.first ) best_result = cur.first;
    out << "Current (best) weight: " << cur.first << " (" << best_result << ")" << ", tm = " << cur_time << endl;

    if( cur_time-init_time > sec )
      return TIME_EXPIRED;

    if( check_ifVectorsEqual( N , A1 , cur.second ) )
      return SUCCESSFULL;
    tryElt( N , cur , B , checkedElements , uncheckedElements );
  }
  
  return FAILED;
}
