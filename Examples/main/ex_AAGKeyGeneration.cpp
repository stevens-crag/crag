// Copyright (C) 2006 Aleksey Myasnikov
// Contents: Example for class AAGProtocolInstance
//
// Principal Authors: Aleksey Myasnikov
//
// Revision History:
//

#include "AAGKeyGeneration.h"


int main( )
{

  // Define parameters

  int N = 80;                             //  The rank of the braid group
  int num_gens = 20;                      //  Number of subgroup generators 
  int min_len = 20;                       //  Minimal length of subgroup generators
  int max_len = 23;                       //  Maximal length of subgroup generators
  int AliceDecompositionLength = 100;     //  
  int BobDecompositionLength = 100;



  //& Arithmetica Key Exchange protocol ; How do I generate an instance of AAG protocol with specific parameters?


  // Generate AAG instance
  AAGProtocolInstance AAG = AAGProtocolInstance::random( N , num_gens , min_len , max_len, 
							 AliceDecompositionLength , BobDecompositionLength );
  
  vector< Word > Sbgp_A  = AAG.getAlicePublicSbgp( );
  vector< Word > Sbgp_A2 = AAG.getAliceConjSbgp( );
  vector< Word > Sbgp_B  = AAG.getBobPublicSbgp( );
  vector< Word > Sbgp_B2 = AAG.getBobConjSbgp( );


  // Get keys
  Word keyA = AAG.getAliceKey ( );
  Word keyB = AAG.getBobKey( );
  Word key  = AAG.getSharedKey( );

  cout << "The shared key : " << key << endl;

}
