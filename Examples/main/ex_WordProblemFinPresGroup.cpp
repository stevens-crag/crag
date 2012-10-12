// Copyright (C) 2005 Alexander Ushakov
// Contents: Example for Advanced Dehn algorithm
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//


#include "AdvDehnAlgorithm.h"


int main( )
{

  //& Finitely Presented Group ;  How do I check whether a word represents a trivial element of a finitely presented group


  // Generate a random freely reduced word over alphabet on 2 symbols of length 7 and push it into a vector relators
  vector< Word > relators;
  relators.push_back( Word::randomWord( 2 , 7 ) );
  // Create a group presentation G with 2 generators and the relator set - relators
  FPGroup G( 2 , relators );
  // Generate a random freely reduced word over alphabet on 2 symbols of length 20
  Word w = Word::randomWord( 2 , 20 );

  // Output the group G and the word w
  cout << G << endl;
  cout << w << endl;

  // Create the object incapsulating the Word Problem algorithm for G
  AdvDehnAlgorithm ADA( G , w );
  
  int depth;
  int max_depth = 3;
  bool loop = ADA.isLoop( w );
  bool coset_limit_reached;
  for( depth=0 ; !loop && depth<max_depth ; ++depth ) {
    cout << "iter = " << depth << endl;
    coset_limit_reached = !ADA.builtup( 0 , 2000000 );
    loop = ADA.isLoop( w );
  }
  if( coset_limit_reached )
    cout << "The limit on the coset set was reached" << endl;
  if( !loop ) {
    cout << "Word is not trivial or has depth more than " << max_depth << "." << endl;
    exit( 1 );
  } else {
    cout << "Depth = " << depth << endl;
  }


  
  return 0;
}
