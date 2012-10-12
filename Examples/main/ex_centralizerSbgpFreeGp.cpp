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
  // Create a subgroup of a free group (it will be normal)
  vector< Word > gens;

  Word x(1);
  Word y(2);
  gens.push_back( x^2 );
  gens.push_back( y^2 );
  gens.push_back( -x*-y*x*y );
  gens.push_back( x^2^y );
  gens.push_back( y^2^x );
  SubgroupFG Sbgp( 2 , gens );
  
  
  //& Subgroups of free groups ; How do I compute a centralizer of a subgroup of a free group?
  SubgroupFG Z = Sbgp.centralizer( );
  vector< Word > zs = Z.getGenerators( );
  cout << "Generators of the centralizer:" << endl;
  for( vector< Word >::const_iterator g_it = zs.begin( ) ; g_it!=zs.end( ) ; ++g_it )
    cout << *g_it << endl;
  
  return 0;
}
