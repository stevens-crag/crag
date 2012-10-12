// Copyright (C) 2006 Alexander Ushakov
// Contents: Example for Subgroups of Free Groups.
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//


#include "Word.h"
#include "SubgroupFG.h"


int main( )
{
  int N = 5;
  int K = 6;
  int L = 7;

  // Create a random vector of K freely reduced words each of length L over the alphabet {x_1,...,x_N}
  vector< Word > generators;
  for( int i=0 ; i<K ; ++i )
    generators.push_back( Word::randomWord( N , L ) );


  //& Subgroups of free groups ; How do I create a subgroup of a free group?
  SubgroupFG sbgp( N , generators );


  //& Subgroups of free groups ; How do I oupput a subgroup of a free group?
  cout << "SBGP = " << endl << sbgp << endl;
  

  //& Subgroups of free groups ; How do I check whether two subgroups coincide?
  // Create another random subgroup with the same parameters
  vector< Word > generators2;
  for( int i=0 ; i<K ; ++i )
    generators2.push_back( Word::randomWord( N , L ) );
  SubgroupFG sbgp2( N , generators2 );
  if( sbgp==sbgp2 ) {
    cout << "Soubgroups sbgp and sbgp2 coincide" << endl;
  } else {
    cout << "Soubgroups sbgp and sbgp2 are different" << endl;
  }

  
  //& Subgroups of free groups ; How do I find intersection of two subgroups?
  SubgroupFG sbgp3 = sbgp*sbgp2;
  

  //& Subgroups of free groups ; How do I find a set of Nielsen Generators for a subgroup?
  vector< Word > niels_gens = sbgp.getNielsenGenerators( );

  
  //& Subgroups of free groups ; How do I get a FSA corresponing to a subgroup?
  IntLabeledGraph FSA = sbgp.getFSA( );

  
  //& Subgroups of free groups ; How do I check whether a word belongs to a subgroup?
  Word w = Word::randomWord( N , 15 );        // generate random word
  if( sbgp.doesBelong( w ) ) {                // check if it belongs to a subgroup
    cout << "A word belongs to a subgroup" << endl;
  } else {
    cout << "A word does not belong to a subgroup" << endl;
  }

  //& Subgroups of free groups ; If a word belongs to a subgroup how do I find actual product of initial generators giving the word?
  if( sbgp.doesBelong( w ) ) {
    Word decomposition = sbgp.express( w );
    // output the decomposition
    Word::const_iterator w_it=decomposition.begin( );
    for( int c=0 ; w_it!=decomposition.end( ) ; ++w_it, ++c )
      cout << ( c>0 ? "  *  " : "" ) << ( *w_it>0 ? generators[*w_it-1] : -generators[-*w_it-1] );
    cout << endl;
  }


  //& Subgroups of free groups ; How do I compute the index of a subgroup in a free group?
  int index = sbgp.getIndex( );


  //& Subgroups of free groups ; How do I compute the rank of a subgroup?
  int rank = sbgp.getRank( );




  return 0;
}
