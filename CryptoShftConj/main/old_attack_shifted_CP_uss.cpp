// Copyright (C) 2006 Aleksey Myasnikov
// Contents: Attack on Dehornoy shifted conjugacy based authentication protocol
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//


#include "braid_group.h"
#include "ShortBraidForm.h"
#include "ThLeftNormalForm.h"
#include "ShftConjKeyGeneration.h"


//---------------------------------------------------------------------------//
//-------------------------------------- main -------------------------------//
//---------------------------------------------------------------------------//


typedef ThLeftNormalForm NF;
set< Word > partofCentralizer( int N , const Word& p1 );
Word getSmallDelta( int rank );


//---------------------------------------------------------------------------//
//-------------------------------------- main -------------------------------//
//---------------------------------------------------------------------------//


Word centr_attack( int N , const Word& p1 , const Word& p2 , const Word& candidate )
{
  Word delta = getSmallDelta( N+1 );


  // Elements from the centralizer
  Word c1 = NF( N+1 ,  2 , list< Permutation >( ) ).getWord( );
  Word c2 =  generatorShift( p1 ) * Word(1) * -delta;
  Word c3 = Word( N )*Word( N );
  for( int i=N-1 ; i>=1; --i ) {
    c3.push_front( i );
    c3.push_back ( i );
  }
  
  typedef triple< int , int , int > TIII;
  set< pair< int , TIII > > new_elts; 
  set< pair< int , TIII > > old_elts; 
  set< TIII > explored;
  
  int L = shortBraidForm( N+2 , candidate * generatorShift(p1) * Word(1) * generatorShift(-candidate) * -p2 ).length( );
  new_elts.insert( pair< int , TIII >( L , TIII(0,0,0) ) );
  explored.insert( TIII(0,0,0) );
  
  while( new_elts.size( ) ) {
    
    pair< int , TIII > cur_elt = *new_elts.begin( );
    new_elts.erase( *new_elts.begin( ) );
    old_elts.insert( cur_elt ); 

    int w = cur_elt.first;
    int i = cur_elt.second.first;
    int j = cur_elt.second.second;
    int k = cur_elt.second.third;
    cout << "  (" << i << "," << j << "," << k << ") -> " << w << endl;
    
    if( w==0 ) {
      cout << "Key is found" << endl;
      return candidate * c1.power(i) * c2.power(j) * c3.power(k);
    }
    
    for( int d=-1 ; d<=1 ; d+=2 ) {
      if( explored.find( TIII( i+d , j , k ) )==explored.end( ) ) {
	Word new_candidate = candidate * c1.power(i+d) * c2.power(j) * c3.power(k);
	Word c = shortBraidForm( N+2, new_candidate * generatorShift(p1) * Word(1) * generatorShift(-new_candidate) * -p2 );
	new_elts.insert( pair<int,TIII>( c.length( ) , TIII( i+d , j , k ) ) );
	explored.insert( TIII( i+d , j , k ) );
      }
      if( explored.find( TIII( i , j+d , k ) )==explored.end( ) ) {
	Word new_candidate = candidate * c1.power(i) * c2.power(j+d) * c3.power(k);
	Word c = shortBraidForm( N+2, new_candidate * generatorShift(p1) * Word(1) * generatorShift(-new_candidate) * -p2 );
	new_elts.insert( pair<int,TIII>( c.length( ) ,TIII( i , j+d , k ) ) );
	explored.insert( TIII( i , j+d , k ) );
      }
      if( explored.find( TIII( i , j , k+d ) )==explored.end( ) ) {
	Word new_candidate = candidate * c1.power(i) * c2.power(j) * c3.power(k+d);
	Word c = shortBraidForm( N+2, new_candidate * generatorShift(p1) * Word(1) * generatorShift(-new_candidate) * -p2 );
	new_elts.insert( pair<int,TIII>( c.length( ) , TIII( i , j , k+d ) ) );
	explored.insert( TIII( i , j , k+d ) );
      }
    }
  }
}


//---------------------------------------------------------------------------//
//-------------------------------------- main -------------------------------//
//---------------------------------------------------------------------------//


Word getSmallDelta( int rank )
{
  Word result;
  for( int i=rank-1 ; i>=1 ; --i )
    result.push_back( i );
  return result;
}


//---------------------------------------------------------------------------//
//-------------------------------------- main -------------------------------//
//---------------------------------------------------------------------------//


int main( )
{
  int N = 80;           //  The rank of the braid group
  int baseLenth = 1000;   //  The length of the base word 
  int keyLength = 1000;   //  The length of the key
  
  // Generate  instance
  
  int exp = 1000;

  for( int e=0 ; e<exp ; ++e ) {

    cout << "++++++++++++++++++++++++++++++" << endl;
    cout << "e = " << e << endl;
    ShftConjKeyInstance SCK = ShftConjKeyInstance::random( N , baseLenth , keyLength );
    
    pair< Word , Word > public_key = SCK.getPublicKey( );
    Word p1 = public_key.first;
    Word p2 = public_key.second;
    
    Word delta = getSmallDelta( N+1 );
    NF nf  = NF( N+1 , generatorShift( p1 ) * Word(1) * -delta );
    NF nf2 = NF( N+1 , p2 * -delta );
    
    pair< bool , NF > res = nf.areConjugate_uss( nf2 );
    Word candidate = shortenBraid( N+1 , res.second.getWord( ) );
    candidate = centr_attack( N , p1 , p2 , candidate );
    
    if( NF( N+2 , candidate * generatorShift(p1) * Word(1) * generatorShift(-candidate) * -p2 ).isTrivial( ) ) {
      cout << "The key is correct" << endl;
    } else {
      cout << "The key is not correct" << endl;
    }
  }
  
  return 0;
}
