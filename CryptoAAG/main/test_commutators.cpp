
#include "time.h"
#include "stdlib.h"
#include "braid_group.h"
#include "ThRightNormalForm.h"
#include "ShortBraidForm.h"

#include <set>
#include <map>
#include <iterator>
#include <iostream>
#include <fstream>
using namespace std;

#include "AAGKeyGeneration.h"


bool addElement( int N , const Word& w , set< Word >& oldWords , set< Word >& newWords , ofstream& OF_details )
{
  Word c = shortBraidForm( N , w );
  if( c.length( ) && c.length( ) < 40 ) {
    if( oldWords.find( c )==oldWords.end( ) && newWords.find( c )==newWords.end( ) ) {
      OF_details << c.length( ) << "  :  " << c << endl;
      if( c.length()==2 ) {
	newWords.clear( );
	newWords.insert( c );
	cout << c << endl;
	return true;
      }
      newWords.insert( c );
    }
  }
  return false;
}


int main( )
{
  ofstream OF_details( "details.txt" );
  
  int N = 80;
  int num_gens = 20;
  int min_len = 30;
  // choose_parameters( res_filename , N  , num_gens , min_len );
  int max_len = min_len+2;
  int AliceDecompositionLength = 10;
  int BobDecompositionLength = 10;
  
  
  OF_details << "------------------------------------------------" << endl;
  OF_details << "Rank of the group = " << N << endl;
  OF_details << "Minimium generator length = " << min_len << endl;
  OF_details << "Number of subgroup generators = " << num_gens << endl;
  OF_details << "------------------------------------------------" << endl;
  
  // Generate AAG instance
  AAGProtocolInstance AAG = AAGProtocolInstance::random( N , num_gens , min_len , max_len , 
							 AliceDecompositionLength , BobDecompositionLength );

  set< Word > oldWords;    
  set< Word > newWords;

  OF_details << "Words:" << endl;
  vector< Word > sbgp  = AAG.getAlicePublicSbgp( );
  for( int i=0 ; i<num_gens ; ++i ) {
    OF_details << sbgp[i].length( ) << "  ";
    newWords.insert( sbgp[i] );
  }
  OF_details << endl;
  OF_details << "------------------------------------------------" << endl;

  while( !newWords.empty( ) ) {
    
    Word a = *newWords.begin( );
    for( set< Word >::const_iterator w_it = oldWords.begin( ) ; w_it!=oldWords.end( ) ; ++w_it ) {
      Word b = *w_it;
      
      addElement( N ,  a * b * -a * -b , oldWords , newWords , OF_details );
      addElement( N , -a * b *  a      , oldWords , newWords , OF_details );
    }
    newWords.erase( a );
    oldWords.insert( a );
  }

  return 0;
}
