// Copyright (C) 2007 Alexander Ushakov
// Contents: Implementation of "Random FSA Generation" methods.
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//


#include "RandomFSA.h"
#include "Permutation.h"
#include <time.h>

#include "gmpxx.h"
typedef mpz_class LongInteger;


//---------------------------------------------------------------------------//
//--------------------------- mixPartialInjection ---------------------------//
//---------------------------------------------------------------------------//


// Take an output of randomPartialInjection, mix it, and output in a more convenient form
vector< int > mixPartialInjection( int N , const list< pair< int , bool > >& PI , gmp_randclass& R )
{
  vector< int > result( N , -1 );
  
  // Generate random permutation 
  int cur_pos = 0;
  Permutation P = Permutation::random( N );
  cout << "N = " << N << endl;
  cout << "P = " << P << endl;
  for( list< pair< int , bool > >::const_iterator PI_it=PI.begin( ) ; PI_it!=PI.end( ) ; ++PI_it ) {
    int l = (*PI_it).first;
    cout << "l = " << l << endl;
    for( int i=0 ; i<l-1 ; ++i ) {
      result[P[cur_pos]] = P[cur_pos+1];
      cur_pos++;
    }
    if( (*PI_it).second )
      result[P[cur_pos]] = P[cur_pos-l+1];
    cur_pos++;
  }
  
  cout << "( ";
  for( vector< int >::const_iterator r_it=result.begin( ) ; r_it!=result.end( ) ; ++r_it )
    cout << " " << *r_it << " ";
  cout << " )" << endl;


  return result;
}


//---------------------------------------------------------------------------//
//------------------------- randomPartialInjection --------------------------//
//---------------------------------------------------------------------------//


// Generate a random partial injection on N elements
// The result is a partition of N where each part is labeled 
// either false - sequence
// or     true  - cycle
vector< int > randomPartialInjection( const int N , vector< LongInteger > I , vector< LongInteger > K , gmp_randclass& R )
{
  list< pair< int , bool > > result;

  
  // cout << "=========================" << endl;
  for( int n=N ; n>0 ; ) {
    
    LongInteger max = I[n]/K[n-1];
    LongInteger dice = R.get_z_range( max );
    // cout << "I = " << I[n] << endl;
    // cout << "K = " << K[n] << endl;
    // cout << "M = " << max << endl;
    // cout << "D = " << dice << endl;
    
    int k = 1;
    LongInteger S = (I[n-1]*2)/K[n-1];
    // cout << "S = " << S << endl;
    while( dice>=S ) {
      k++;
      S += (I[n-k]*(k+1))/K[n-k];
    }
    // cout << "k = " << k << endl;
    // cout << endl;
    
    result.push_back( pair< int , bool >( k , R.get_z_range( k )==0  ) );
    n -= k;
  }
  
  const list< pair< int , bool > >& PI = result;
  for( list< pair< int , bool > >::const_iterator PI_it=PI.begin( ) ; PI_it!=PI.end( ) ; ++PI_it )
    cout << "(" << (*PI_it).first << "," << (*PI_it).second << ") ";
  cout << endl;

  return mixPartialInjection( N , result , R );
}


//---------------------------------------------------------------------------//
//------------------------------ RandomFSA ----------------------------------//
//---------------------------------------------------------------------------//


FSA randomFSA( int N , int L )
{
  // Precompute the list of I-values
  vector< LongInteger > I;
  I.push_back( 1 );
  I.push_back( 2 );
  for( int i=2 ; i<=N ; ++i ) {
    LongInteger A = *----I.end( );
    LongInteger B = *--I.end( );
    I.push_back( 2*i*B - (i-1)*(i-1)*A );
  }
  
  // Precompute the list of K-values
  vector< LongInteger > K;
  K.push_back( 1 );
  for( int i=1 ; i<=N ; ++i )
    K.push_back( *--K.end( ) * i );
  
  // Prepare the set of vertices
  FSA fsa;
  for( int i=0 ; i<N ; ++i ) {
    fsa.newState( );
  }

  gmp_randclass R( gmp_randinit_default );
  R.seed( time( 0 ) );

  // Generate the set of edges
  for( int l=0 ; l<L ; ++l ) {
    vector< int > PI = randomPartialInjection( N , I , K , R );
    int i=0;
    for( vector< int >::const_iterator pi_it=PI.begin( ) ; pi_it!=PI.end( ) ; ++pi_it, ++i ) {
      if( *pi_it!=-1 ) {
	fsa.newEdge( i , *pi_it ,  l+1 );
	fsa.newEdge( *pi_it , i , -l-1 );
      }
    }
  }

  return fsa;
}
