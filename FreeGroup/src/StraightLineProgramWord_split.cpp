// Copyright (C) 2007 Alexander Ushakov
// Contents: Implementation of class StraightLineProgramWord
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//

#define DEBUG_WP_AUT

#include "tuples.h"
#include "StraightLineProgramWord.h"

#include "set"
using namespace std;


//---------------------------------------------------------------------------//
//------------------------------- splitAssertion ----------------------------//
//---------------------------------------------------------------------------//


triple< bool , bool , LongInteger > needRule( LongInteger L1 , LongInteger L2 , LongInteger L3 , LongInteger L4 )
{
  if( L2-L1==0 || L4-L3==0 || L1>=L4 || L3>=L2 )
    return triple< bool , bool , LongInteger >( false , true , 0 );
  return triple< bool , bool , LongInteger >( true , L1<=L3 , abs(L1-L3) );
}


set< StraightLineProgramWord::Assertion > 
StraightLineProgramWord::splitAssertion( const Assertion& A , bool firstTerm , const StraightLineProgramWord& SLP ) const
{
  const StraightLineProgramWord& SLP1 = *this;
  const StraightLineProgramWord& SLP2 = SLP;
  int V1    = A.theVertex1;
  int V2    = A.theVertex2;
  

  // Get the productions for the first rule
  map< bool , vector< int > > V;
  if( abs(V1)<=theTerminals || !firstTerm ) {
    V[false].push_back( V1 );
    V[false].push_back( 0 );
  } else {
    const Production& pr1 = (*SLP1.theRules.find(abs(V1))).second;
    V[false].push_back( pr1.theTerm1 );
    V[false].push_back( pr1.theTerm2 );
    if( V1<0 ) 
      invertProductionPair( V[false][0] , V[false][1] );
  }


  // Get the productions for the second rule
  if( abs(V2)<=theTerminals  || firstTerm ) {
    V[true].push_back( V2 );
    V[true].push_back( 0 );
  } else {
    const Production& pr2 = (*SLP2.theRules.find(abs(V2))).second;
    V[true].push_back( pr2.theTerm1 );
    V[true].push_back( pr2.theTerm2 );
    if( V2<0 ) 
      invertProductionPair( V[true][0] , V[true][1] );
  }
#ifdef DEBUG_WP_AUT
  cout << V[false][0] << " / " << V[false][1] << endl;
  cout << V[true][0] << " / " << V[true][1] << endl;
#endif
  
  // Get the lengths of the production parts
  LongInteger LA = SLP1.length_rule( V[false][0] );
  LongInteger LB = SLP1.length_rule( V[false][1] );
  LongInteger LC = SLP2.length_rule( V[true][0] );
  LongInteger LD = SLP2.length_rule( V[true][1] );
  
#ifdef DEBUG_WP_AUT
  cout << LA << " * " << LB << endl;
  cout << LC << " * " << LD << endl;
#endif 
  
  // Get the ralitive positions of production parts
  map< bool , vector< LongInteger > > L;
  if( A.theBase1 ) {
    L[false].push_back( 0 );
    L[false].push_back( LA );
    L[false].push_back( LA+LB );
    L[true] .push_back( A.theLength );
    L[true] .push_back( A.theLength+LC );
    L[true] .push_back( A.theLength+LC+LD );
  } else {
    L[false].push_back( A.theLength );
    L[false].push_back( A.theLength+LA );
    L[false].push_back( A.theLength+LA+LB );
    L[true] .push_back( 0 );
    L[true] .push_back( LC );
    L[true] .push_back( LC+LD );
  }

  // Find the result
  set< Assertion > result;
  for( int i=0 ; i<=1 ; ++i ) {
    for( int j=0 ; j<=1 ; ++j ) {
      triple< bool , bool , LongInteger > nr = needRule( L[false][i], L[false][i+1], L[true][j], L[true][j+1] );
#ifdef DEBUG_WP_AUT
      cout << "   Consider intervals: (" 
	   << L[false][i] << "," << L[false][i+1] << ") + (" 
	   << L[true][j] << "," << L[true][j+1] << ") -> " 
	   << nr.first << ", " << nr.second << ", " << nr.third << endl;
#endif
      if( nr.first )
	result.insert( Assertion( nr.second , V[false][i] , V[true][j] , nr.third ) );
    }
  }

  return result;
}
