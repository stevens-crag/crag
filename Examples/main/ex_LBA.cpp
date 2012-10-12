// Copyright (C) 2006 Aleksey Myasnikov
// Contents: Example for class AAGProtocolInstance
//
// Principal Authors: Aleksey Myasnikov
//
// Revision History:
//


#include <iostream>
using namespace std;

#include "LengthAttack.h"
#include "AAGKeyGeneration.h"


//---------------------------------------------------------------------------//
//---------------------------------- main -----------------------------------//
//---------------------------------------------------------------------------//


int main( int argc, char** argv )
{


  // Define parameters
  int N = 80;           // Rank of Braid group
  int num_gens = 20;    // Number of SG generators
  int min_len = 20;     // Minimal generator length
  int max_len = 23;     // Maximal generator length
  int AliceDecompositionLength = 50; // Length of the Alice's secret
  int BobDecompositionLength = 50;   // Length of the Bob's secret
  
  // Set time limit
  int sec_in_day = 86400;
  int TIME_LIMIT_SEC = 1*sec_in_day;
  

  //& Arithmetica Key Exchange protocol ; How do I execute the Length-Based attack on an instance of AAG protocol with specific parameters?
 

  // Generate AAG instance
  AAGProtocolInstance AAG = AAGProtocolInstance::random( N , num_gens , min_len, max_len, 
							 AliceDecompositionLength , BobDecompositionLength );
  
  vector< Word > Sbgp_A = AAG.getAlicePublicSbgp( );
  vector< Word > Sbgp_A2 = AAG.getAliceConjSbgp();
  vector< Word > Sbgp_B = AAG.getBobPublicSbgp( );
  
  // Execute the  attack 
  LengthAttack_A2 A;
  cout << "Executing algorithm : A" << A.type() << endl;
  
  // try to find Bob' key
  cout << "Start attack ... " << endl;
  switch( A.findKey_LengthBased( N , Sbgp_A , Sbgp_A2 , Sbgp_B , TIME_LIMIT_SEC, cout ) ) {
  case SUCCESSFULL:
    cout << "Success" << endl;
    break;
  case FAILED:
    cout << "Failed" << endl;
    break;
  case TIME_EXPIRED:
    cout << "Time expired" << endl;
    break;
  }
  
  return 0;
}
