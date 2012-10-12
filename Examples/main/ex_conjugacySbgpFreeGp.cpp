// Copyright (C) 2006 Alexander Ushakov
// Contents: Example for SubgroupFG::areConjugate( )
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//

#include <fstream>
#include "SubgroupFG.h"


int main( )
{
  //& Subgroups of free groups ; How do I check if two subgroup are conjugate?
  int N = 2;
  vector< Word > gens;
  SubgroupFG Sbgp( N );
  for( int i=0 ; i<2 ; ++i )
    Sbgp += Word::randomWord( N , 5 );
  
  Word c1 = Word::randomWord( N , 5 );
  SubgroupFG Sbgp1 = Sbgp^c1;
  
  Word c2 = Word::randomWord( N , 5 );
  SubgroupFG Sbgp2 = Sbgp^c2;

  // Determine if subgroups Sbgp1 and Sbgp2 are conjugate
  pair< bool , Word > C = Sbgp1.areConjugate( Sbgp2 );

  if( C.first ) {                                   // If subgroups are conjugate
    cout << "Subgroups are conjugate" << endl;

    // Check the correctness of the conjugator
    if( Sbgp2==(Sbgp1^C.second) ) {
      cout << "conjugator is correct" << endl;
    } else {
      cout << "conjugator is incorrect" << endl;
    }
    
  } else {                                          // If subgroups are not conjugate
    cout << "Subgroups are not conjugate" << endl;
  }
  
  return 0;
}

