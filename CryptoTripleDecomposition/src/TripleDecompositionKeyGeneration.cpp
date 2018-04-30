// Copyright (C) 2005 Alexander Ushakov
// Contents: Implementation of class TripleDecompositionProtocolInstance
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//


#include "TripleDecompositionKeyGeneration.h"

#include "RanlibCPP.h"
#include "DehornoyForm.h"
#include "braid_group.h"


//---------------------------------------------------------------------------//
//-------------------- TripleDecompositionProtocolInstance ------------------//
//---------------------------------------------------------------------------//


TripleDecompositionProtocolInstance::TripleDecompositionProtocolInstance
( int braid_rank , 
  quadruple< Word , Word , Word , Word > conjugators ,
  quintuple< Word , Word , Word , Word , Word > privateKeyA ,
  quintuple< Word , Word , Word , Word , Word > privateKeyB ) :
  theRank( braid_rank ),
  theConjugators( conjugators ),
  thePrivateKeyA( privateKeyA ),
  thePrivateKeyB( privateKeyB )
{
  thePublicKeyA.first  = DehornoyForm( braid_rank ,  privateKeyA.first  * privateKeyA.fourth ).getDehornoyForm( );
  thePublicKeyA.second = DehornoyForm( braid_rank , -privateKeyA.fourth * privateKeyA.second * privateKeyA.fifth ).getDehornoyForm( );
  thePublicKeyA.third  = DehornoyForm( braid_rank , -privateKeyA.fifth  * privateKeyA.third ).getDehornoyForm( );
  
  thePublicKeyB.first  = DehornoyForm( braid_rank ,  privateKeyB.first  * privateKeyB.fourth ).getDehornoyForm( );
  thePublicKeyB.second = DehornoyForm( braid_rank , -privateKeyB.fourth * privateKeyB.second * privateKeyB.fifth ).getDehornoyForm( );
  thePublicKeyB.third  = DehornoyForm( braid_rank , -privateKeyB.fifth  * privateKeyB.third ).getDehornoyForm( );
  
  theSharedKey = ThRightNormalForm( braid_rank , privateKeyA.first * privateKeyB.first * privateKeyA.second * privateKeyB.second * privateKeyA.third * privateKeyB.third );
  
  
}


//---------------------------------------------------------------------------//
//-------------------------------- random -----------------------------------//
//---------------------------------------------------------------------------//


TripleDecompositionProtocolInstance TripleDecompositionProtocolInstance::random
( int braid_rank , int conjugator_length , int element_length )
{
  int lowerIndex = braid_rank/3;
  int upperIndex = 2*lowerIndex;
  
  // A. Choose randomly the set of conjugators
  quadruple< Word , Word , Word , Word > C;
  C.first  = Word::randomWord( braid_rank-1 , conjugator_length );
  C.second = Word::randomWord( braid_rank-1 , conjugator_length );
  C.third  = Word::randomWord( braid_rank-1 , conjugator_length );
  C.fourth = Word::randomWord( braid_rank-1 , conjugator_length );
  
  // B. Choose randomly the private keys for both parties
  quintuple< Word , Word , Word , Word , Word > thePrivateKeyA;
  thePrivateKeyA.first  = randomWord( 1 , braid_rank-1 , element_length );                        // a1
  thePrivateKeyA.fourth = -C.first  * randomWord( 1 , lowerIndex-1 , element_length ) * C.first;  // x1
  thePrivateKeyA.second = -C.second * randomWord( 1 , lowerIndex-1 , element_length ) * C.second; // a2
  thePrivateKeyA.fifth  = -C.third  * randomWord( 1 , upperIndex-1 , element_length ) * C.third;  // x2
  thePrivateKeyA.third  = -C.fourth * randomWord( 1 , upperIndex-1 , element_length ) * C.fourth; // a3

  quintuple< Word , Word , Word , Word , Word > thePrivateKeyB;
  thePrivateKeyB.first  = -C.first  * randomWord( lowerIndex+1 , braid_rank-1 , element_length ) * C.first;  // b1
  thePrivateKeyB.fourth = -C.second * randomWord( lowerIndex+1 , braid_rank-1 , element_length ) * C.second; // y1
  thePrivateKeyB.second = -C.third  * randomWord( upperIndex+1 , braid_rank-1 , element_length ) * C.third;  // b2
  thePrivateKeyB.fifth  = -C.fourth * randomWord( upperIndex+1 , braid_rank-1 , element_length ) * C.fourth; // y2
  thePrivateKeyB.third  = randomWord( 1            , braid_rank-1 , element_length );                        // b3

  // return the result
  return TripleDecompositionProtocolInstance( braid_rank , C , thePrivateKeyA , thePrivateKeyB );
}


//---------------------------------------------------------------------------//
//------------------------------ randomWord ---------------------------------//
//---------------------------------------------------------------------------//


Word TripleDecompositionProtocolInstance::randomWord( int lowerIndex , int upperIndex , int len )
{
  Word result;

  int r = upperIndex - lowerIndex + 1;

  int old = 0;
  for( int i=0 ; i<len ; ++i ) {
    
    int div = i==0 ? 2*r : 2*r-1;
    int g = RandLib::ur.irand( 0 , div-1 )-r;
    g = g>=0 ? g+1 : g;
    if( g+old==0 )
      g = r;
    old = g;
    result.push_back( g<0 ? g-lowerIndex+1 : g+lowerIndex-1 );
  }
  
  return result;
}
