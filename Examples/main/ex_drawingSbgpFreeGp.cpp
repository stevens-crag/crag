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
  //& Subgroups of free groups ; How can I visualize the subgroup graph?
  // Create a random subgroup of a free group
  int N = 2;
  SubgroupFG Sbgp( N );
  for( int i=0 ; i<3 ; ++i )
    Sbgp += Word::randomWord( N , 7 );
  ofstream OF( "fsa.txt" );
  OF << Sbgp.graphviz_format( ) << endl;;
  
  // Now, you can open file "fsa.txt" in Graphviz
  
  return 0;
}
