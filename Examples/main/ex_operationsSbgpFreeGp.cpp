// Copyright (C) 2006 Alexander Ushakov
// Contents: Example for SubgroupFG::centralizer( )
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//

#include <fstream>
#include "SubgroupFG.h"


int main( )
{
  //& Subgroups of free groups ; How do I create trivial subgroup of a free group?
  int N = 2;
  SubgroupFG Sbgp( N );
  

  //& Subgroups of free groups ; How do I extend the subgroup basis?
  // Extend the basis with 3 random words of length 7
  for( int i=0 ; i<3 ; ++i )
    Sbgp += Word::randomWord( N , 7 );
  

  //& Subgroups of free groups ; How do I conjugate a subgroup?
  Word conjugator = Word::randomWord( N , 7 );
  SubgroupFG Sbgp2 = Sbgp^conjugator;
  
  
  //& Subgroups of free groups ; How do I trim the subgroup graph (erase the tail from the subgroup graph)?
  pair< SubgroupFG , Word > pr = Sbgp2.trim( );
  // Check the returned conjugator
  if( pr.first==(Sbgp2^pr.second) ) {
    cout << "Correct trimming" << endl;
  } else {
    cout << "Problem with trimming" << endl;
  }
  
  
  return 0;
}
