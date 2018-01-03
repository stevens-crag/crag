// Copyright (C) 2006 Aleksey Myasnikov
// Contents: Attack on Dehornoy shifted conjugacy based authentication protocol
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//


#include "fstream"
#include "strstream"
#include "BraidGroup.h"
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


void test_inc_rank( )
{
  int N = 20;

  for( int i=0 ; i<100 ; ++i ) {
    Word w = Word::randomWord( N-1 , 200 ); 
    NF nf( N , w );
    NF nf2 = nf.increaseRank( N+1 );
    
    Word u = -w * nf2.getWord( );
    NF nf_check( N+1 , -w * nf2.getWord( ) );
    if( nf_check.isTrivial( ) ) {
      cout << "The result is correct" << endl;
    } else {
      cout << "The result is not correct" << endl;
      exit(1);
    }
  }
}


//---------------------------------------------------------------------------//
//-------------------------------------- main -------------------------------//
//---------------------------------------------------------------------------//


int getLength( const NF& nf )
{
  return abs(nf.getPower())+nf.getDecomposition( ).size( );
}


//---------------------------------------------------------------------------//
//-------------------------------------- main -------------------------------//
//---------------------------------------------------------------------------//


NF centr_attack( int N , const NF& p1_2 , const NF& p2_2 , const NF& candidate2 , const NF& delta2 )
{
  NF delta = NF( N+2 , getSmallDelta( N+1 ) );
  NF s2( N+2 , Word(1) );
  
  // Elements from the centralizer
  NF c1( N+1 ,  2 , list< Permutation >( ) );
  c1 = c1.increaseRank( N+2 );
  NF c2 = (-delta2*p1_2*delta2) * s2 * -delta;
  Word wc3 = Word( N )*Word( N );
  for( int i=N-1 ; i>=1; --i ) {
    wc3.push_front( i );
    wc3.push_back ( i );
  }
  NF c3( N+2 , wc3 );

  typedef triple< int , int , int > TIII;
  set< triple< int , TIII , NF > > new_elts; 
  set< pair< int , TIII > > old_elts; 
  set< TIII > explored;

  int L = getLength( candidate2 * (-delta2*p1_2*delta2) * s2 * (-delta2*-candidate2*delta2) * -p2_2 );
  new_elts.insert( triple< int , TIII , NF >( L , TIII(0,0,0) , candidate2 ) );
  explored.insert( TIII(0,0,0) );
  
  while( new_elts.size( ) ) {
    
    triple< int , TIII , NF > cur_elt = *new_elts.begin( );
    new_elts.erase( *new_elts.begin( ) );
    old_elts.insert( pair< int , TIII >( cur_elt.first , cur_elt.second ) ); 

    int w = cur_elt.first;
    int i = cur_elt.second.first;
    int j = cur_elt.second.second;
    int k = cur_elt.second.third;
    NF conj = cur_elt.third;
    cout << "  (" << i << "," << j << "," << k << ") -> " << w << endl;
    
    if( w==0 ) {
      cout << "Key is found" << endl;
      return conj;
    }
    
    for( int d=-1 ; d<=1 ; d+=2 ) {
      if( explored.find( TIII( i+d , j , k ) )==explored.end( ) ) {
	NF new_conj = conj * ( d==-1 ? -c1 : c1 );
	int L = getLength( new_conj * (-delta2*p1_2*delta2) * s2 * (-delta2*-new_conj*delta2) * -p2_2 );
	new_elts.insert( triple<int,TIII,NF>( L , TIII( i+d , j , k ) , new_conj ) );
	explored.insert( TIII( i+d , j , k ) );
      }
      if( explored.find( TIII( i , j+d , k ) )==explored.end( ) ) {
	NF new_conj = conj * ( d==-1 ? -c2 : c2 );
	int L = getLength( new_conj * (-delta2*p1_2*delta2) * s2 * (-delta2*-new_conj*delta2) * -p2_2 );
	new_elts.insert( triple<int,TIII,NF>( L , TIII( i , j+d , k ) , new_conj ) );
	explored.insert( TIII( i , j+d , k ) );
      }
      if( explored.find( TIII( i , j , k+d ) )==explored.end( ) ) {
	NF new_conj = conj * ( d==-1 ? -c3 : c3 );
	int L = getLength( new_conj * (-delta2*p1_2*delta2) * s2 * (-delta2*-new_conj*delta2) * -p2_2 );
	new_elts.insert( triple<int,TIII,NF>( L , TIII( i , j , k+d ) , new_conj ) );
	explored.insert( TIII( i , j , k+d ) );
      }
    }
  }
}


//---------------------------------------------------------------------------//
//-------------------------------------- main -------------------------------//
//---------------------------------------------------------------------------//


Word getSmallDelta( int N )
{
  Word result;
  for( int i=N-1 ; i>=1 ; --i )
    result.push_back( i );
  return result;
}


//---------------------------------------------------------------------------//
//-------------------------------------- main -------------------------------//
//---------------------------------------------------------------------------//


int main( )
{
  int Numbers[] = {10,40,80};      //  The rank of the braid group
  int exp = 100;                   //  The number of experiments to perform for each parameter value
  int Lengths[] = {100,400,800};   //  The length of the base word 
  
  // Generate  instance
  
  for( int n=0 ; n<sizeof(Numbers)/sizeof(int) ; ++n ) {
    for( int l=0 ; l<sizeof(Lengths)/sizeof(int) ; ++l ) {

      int N = Numbers[n];
      int L = Lengths[l];
  
      int success_orig = 0;
      int success_new = 0;
      for( int e=0 ; e<exp ; ++e ) {
	
	cout << "++++++++++++++++++++++++++++++" << endl;
	cout << "e = " << e << endl;
	
	ShftConjKeyInstance SCK = ShftConjKeyInstance::random( N , L , L );
	
	pair< Word , Word > public_key = SCK.getPublicKey( );
	Word p1 = public_key.first;
	Word p2 = public_key.second;
	
	Word delta = getSmallDelta( N+1 );
	NF nf  = NF( N+1 , generatorShift( p1 ) * Word(1) * -delta );
	NF nf2 = NF( N+1 , p2 * -delta );
	
	int time_sec_bound = 60*60  ; // time bound is 60 minutes
	pair< bool , NF > res = nf.areConjugate_uss( nf2 , time_sec_bound );
	if( !res.first ) {
	  cout << "USS failure" << endl;
	  continue;
	}
	NF candidate = res.second;
	
	NF candidate2 = candidate.increaseRank( N+2 );
	NF s2( N+2 , Word(1) );
	NF delta2( N+2 , getSmallDelta( N+2 ) );
	NF p1_2( N+2 , p1 );
	NF p2_2( N+2 , p2 );
	
	candidate2 = centr_attack( N , p1_2 , p2_2 , candidate2 , delta2);
	
	if( ( candidate2 * (-delta2*p1_2*delta2) * s2 * (-delta2*-candidate2*delta2) * -p2_2 ).isTrivial( ) ) {
	  cout << "The key is correct" << endl;

	  // Now check with the original private key
	  NF priv = NF( N+2 , SCK.getPrivateKey( ) );
	  if( priv==candidate2 ) {
	    cout << "The original key obtained" << endl;
	    success_orig++;
	  } else {
	    cout << "The obtained key is new" << endl;
	    success_new++;
	  }
	} else {
	  cout << "The key is not correct" << endl;
	}
	
      }


      const int filename_sz = 100;
      char filename[filename_sz];
      ostrstream ostr( filename , filename_sz );
      ostr << "results_N" << N << "_L" << L << ".txt" << ends;
      
      ofstream of( filename );
      of << "N = " << n << endl;
      of << "baseLenth = " << l << endl;
      of << "keyLength = " << l << endl;
      of << "Total experiments: " << exp << endl;
      of << "Successful experiments: " << success_new+success_orig << endl;
      of << "Percentage_orig: " << 100*success_orig/exp << endl;
      of << "Percentage_new: " << 100*success_new/exp << endl;
      // test_inc_rank( );
      
    }
  }
  
  return 0;
}
