// Copyright (C) 2006 Aleksey Myasnikov
// Contents: Example for class Permutation
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//


#include "Permutation.h"

#include "iostream"
using namespace std;


int main( )
{
  //---------------------------------------------------------------------------//
  //------------------------- Examples: Word ----------------------------------//
  //---------------------------------------------------------------------------//




  //& Permutation; How do I create a trivial permutation on N symbols?
  int N = 4;
  Permutation p1(N);
  cout << "Trivial permutation on " << N << " symbols: " << p1 << endl;


  

  //& Permutation; How do I create a permutation on N symbols from array?
  int p_array[] = { 3 , 0 , 2 , 1 };
  Permutation p2( 4 , p_array );




  //& Permutation; How do I create a permutation on N symbols from a vector?
  vector<int> v(4,0);                         // vector of 4 elements
  v[0] = 1; v[1] = 0; v[2] = 4; v[3] = 3;     // initialize elements
  Permutation p3( v );



  
  //& Permutation; How do I create a permutation on N symbols by a pair of iterators?
  Permutation p4( v.begin() , v.end() );
  
  


  //& Permutation; How do I get the 1st element of a permutation?
  int first_elt = p1[0];      // everything starts from 0!!!
  



  //& Permutation; How do I compare 2 permutations?
  if( p1==p2 )
    cout << "Permutations " << p1 << " and " << p2 << " are equal" << endl;
  if( p1!=p2 )
    cout << "Permutations " << p1 << " and " << p2 << " are not equal" << endl;
  if( p1<p2 )
    cout << "Permutation " << p1 << " is less than " << p2 << endl;



  
  //& Permutation; How do I perform algebraic operations with permutations?
  // 1) How do I multiply permutations?
  Permutation p23 = p2 * p3;
  // 2) How do I invert a permutation?
  Permutation p2_inv = -p2;
  // 3) How do I raise a permutation into a power?
  Permutation p2_pwr = p2.power( -2 );
  // 4) How do I extend a permutation to more symbols?
  Permutation p2_ext = p2.increaseSize( 6 );
  


  
  //& Permutation; How do I get a geodesic word for a permutation?
  vector< int > geodesic = p2.geodesic( );
  
  
  
  
  //& Permutation; How do I perform lattice operations with permutations?
  Permutation r_gcd = p2.RightGCD( p3 );
  Permutation r_lcm = p2.RightLCM( p3 );
  Permutation l_gcd = p2. LeftGCD( p3 );
  Permutation l_lcm = p2. LeftLCM( p3 );
  
}
