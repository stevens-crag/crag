// Copyright (C) 2006 Aleksey Myasnikov
// Contents: Example for class ShftConjKeyInstance
//
// Principal Authors: Aleksey Myasnikov
//
// Revision History:
//

#include "ShftConjKeyGeneration.h"


int main( )
{

  // Define parameters

  int N = 80;           //  The rank of the braid group
  int baseLenth = 20;   //  The length of the base word 
  int keyLength = 20;   //  The length of the key
 


  //& Shifted conjugacy authentication ; How do I generate a random instance of Shifted Conjugacy  protocol with specific parameters?


  // Generate  instance

  ShftConjKeyInstance SCK = ShftConjKeyInstance::random( N,baseLenth,keyLength );

  pair< Word , Word > keyA = SCK.getPublicKey( );

}
