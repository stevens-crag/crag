// Copyright (C) 2005 Alexander Ushakov
// Contents: Implementation of length attack
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//

// A2

#include <set>
#include "LengthAttack.h"
#include "ShortBraidForm.h"
#include <time.h>

int LengthAttack_A2::sbgpGeneratorsWeight( const vector<Word>& A )
{
  int result = 0;
  for( int i=0 ; i<A.size() ; ++i )
    result += A[i].length( );
  return result;
}

void LengthAttack_A2::addNewElt( const vector<Word>& A , set< ELT >& checkedElements , set< ELT >& uncheckedElements )
{
  int weight = sbgpGeneratorsWeight( A );
  ELT new_elt( weight , A );

  if(   checkedElements.find( new_elt )!=  checkedElements.end( ) ) return;
  if( uncheckedElements.find( new_elt )!=uncheckedElements.end( ) ) return;

  uncheckedElements.insert( new_elt );
}


void LengthAttack_A2::tryElt( int N , const ELT& cur , const vector< Word >& B , set< ELT >& checkedElements , set< ELT >& uncheckedElements, ostream& out )
{
  // we better vary this value, depending on parameters of A and B
  int MAX_DELTA = 0;

  int max_delta_observed = 0;
  int max_neg_criteria = 0;
  for( int i=0 ; i<B.size( ) ; ++i ) {

    Word b = B[i];
    for( int d=0 ; d<2 ; ++d ) {
    
      //cout << "   #" << i+1 << "," << d+1;
      if( d==1 ) b = -b;
      bool candidate = true;
      int delta = 0;
      int neg_sum = 0;
      
      vector< Word > A = cur.second;
      for( int t=0 ; t<A.size( ) && candidate ; ++t ) {
	int old_len = A[t].length( );
	delta -= A[t].length( );
	A[t] = shortenBraid( N , -b*A[t]*b );
	delta += A[t].length( );
	
	if (old_len <=  A[t].length( )) neg_sum++;

	if( delta>MAX_DELTA && t>=4 ) {
	  candidate = false;
	  //  out << " stopped @";
	}
      }
      
      if (delta < 0 ){
	out << endl << "Try " << "B_" << i+1 << " : ";
	out << " d = " << delta << " ";
      }

      if ( max_delta_observed > delta ){
	max_neg_criteria = neg_sum;
	max_delta_observed = delta;
      }
      if( candidate ) {
	addNewElt( A , checkedElements , uncheckedElements );
	//cout << " accepted with " << delta << endl;
      }
    }
  }

  out << " " << max_delta_observed << " " << max_neg_criteria << " ";
  // if not enough, do conjugations
  // || max_neg_criteria > 1
  if ((max_delta_observed > -100  ) && B.size() < cur.second.size()*cur.second.size()){
    if (max_delta_observed > -100 )
      out << "Maximal decrease is less then 200. ";
    if (max_neg_criteria > 1)
      out << "Too many non-positive in the maximal decrease. ";
  
    out << endl << "Try conjugations ..." << endl;
    // extend B to conjugates
    vector<Word> ext_B;
    for (int i=0;i<B.size();i++)
      for (int j=0;j<B.size();j++)
	if (i!=j){
	  ext_B.push_back(B[i]*B[j]*-B[i]);
	  ext_B.push_back(-B[i]*B[j]*B[i]);
	}

    tryElt( N , cur , ext_B , checkedElements , uncheckedElements, out );

    
  }
  
  out << endl;
  
}


bool LengthAttack_A2::check_ifVectorsEqual( int N , const vector< Word >& A1 , const vector< Word >& A2 )
{
  int len = A1.size( );
  for( int i=0 ; i<len ; ++i ) {
    Word s = shortenBraid( N , A1[i]*-A2[i]  );
    if( s.length( )>0 )
      return false;
  }

  return true;
}


findKey_LengthBasedResult LengthAttack_A2::findKey_LengthBased( int N , const vector< Word >& A1 , const vector< Word >& A2 , const vector< Word >& B , int sec, ostream& out )
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
    
    ELT cur = *uncheckedElements.begin( );
    
    int cur_time = time( 0 );
    out << "Current (best) weight: " << cur.first << " (" << best_result << ")" << ", tm = " << cur_time << endl;
    for( vector< Word >::const_iterator c_it=cur.second.begin( ) ; c_it!=cur.second.end( ) ; ++c_it )
      out << (*c_it).length() << "  ";
    out << endl;
    
    
    uncheckedElements.erase( uncheckedElements.begin( ) );
    checkedElements.insert( cur );

    
    if( best_result>cur.first ) 
      best_result = cur.first;
    
    // Check time
    if( cur_time-init_time > sec )
      return TIME_EXPIRED;

    if( check_ifVectorsEqual( N , A1 , cur.second ) )
      return SUCCESSFULL;
    
    tryElt( N , cur , B , checkedElements , uncheckedElements, out );
  }
  
  return FAILED;
}
