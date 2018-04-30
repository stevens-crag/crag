// Copyright (C) 2005 Alexander Ushakov
// Contents: Example for class BraidGroup
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//


#include "braid_group.h"


int main( )
{
  // Create the braid group on N strands object
  BraidGroup B( 10 );
  
  // Get the rank of the braid group
  int N = B.getRank( );
  
  // Generate a random braid word of length L ( the first parameter is the number of generators)
  int L = 20;
  Word bw = Word::randomWord( N-1 , L );
  
  // Twist the braid word (conjugate by a half-twist permutation) in B
  Word h_twist = B.twist( bw );
  
  
  return 0;
}
