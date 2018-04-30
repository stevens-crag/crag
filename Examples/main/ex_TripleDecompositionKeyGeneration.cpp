// Copyright (C) 2005 Alexander Ushakov
// Contents: Example for class TripleDecompositionProtocolInstance
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//

#include "TripleDecompositionKeyGeneration.h"
#include "braid_group.h"
#include "ThRightNormalForm.h"


int main( )
{
  // Define parameters
  int N = 15;
  int baseLenth = 10;
  int keyLength = 10;
  
  
  //& Triple Decomposition Key Exchange protocol ; How do I generate an instance of Triple Decomposition protocol with specific parameters?
  TripleDecompositionProtocolInstance TD = TripleDecompositionProtocolInstance::random( N , baseLenth , keyLength );
  
  
  //& Triple Decomposition Key Exchange protocol ; How do I retrieve Alice's and Bob's public information?
  triple< Word , Word , Word > publicKeyA = TD.getPublicKeyA( );
  triple< Word , Word , Word > publicKeyB = TD.getPublicKeyB( );

  
  
  //& Triple Decomposition Key Exchange protocol ; How do I retrieve Alice's and Bob's private information?
  quintuple< Word , Word , Word , Word , Word > privateKeyA = TD.getPrivateKeyA( );
  quintuple< Word , Word , Word , Word , Word > privateKeyB = TD.getPrivateKeyB( );
  
  
  
  // Check the shared key correctness:
  ThRightNormalForm c1( N , privateKeyA.first * publicKeyB.first * privateKeyA.second * publicKeyB.second * privateKeyA.third * publicKeyB.third );
  
  ThRightNormalForm c2( N , publicKeyA.first * privateKeyB.first * publicKeyA.second * privateKeyB.second * publicKeyA.third * privateKeyB.third );
  
  if( c1==c2 && c1==TD.getSharedKey( ) ) {
    cout << "Shared key: ok" << endl;
  } else {
    cout << "Shared key: error" << endl;
  }
  
  return 0;
}
