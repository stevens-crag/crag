// Copyright (C) 2005 Alexander Ushakov
// Contents: Example for class KLProtocolInstance
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//

#include "KLKeyGeneration.h"


int main( )
{
  // Parameters:
  // N - the rank of tHe braid group
  int N = 10;
  // - the lengh of the base word
  int baseLenth = 40;
  // - the length of the private keys
  int keyLenth = 20;
  
  //& Ko-Lee protocol ; How do I generate an instance of Ko-Lee protocol with specific parameters?
  KLProtocolInstance KL = KLProtocolInstance::random( N , baseLenth , keyLenth );
  
  // Get the keys
  int  theRank = KL.getBraidRank  ( );
  Word theKeyA = KL.getPrivateKeyA( );
  Word theKeyB = KL.getPrivateKeyB( );
  Word pubKeyA = KL.getPublicKeyA ( );
  Word pubKeyB = KL.getPublicKeyB ( );
  
  
  
  return 0;
}
