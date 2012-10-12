// Copyright (C) 2006 Alexander Ushakov
// Contents: Example for class Word
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//


#include "TheGrigorchukGroupAlgorithms.h"


//---------------------------------------------------------------------------//
//------------------------- Examples: Word ----------------------------------//
//---------------------------------------------------------------------------//



int main( )
{
  
  //& The Grigorchuk Group ; How do I randomly generate a random reduced word?
  int L = 10;
  Word w = TheGrigorchukGroupAlgorithms::randomWord( L );
  


  //& The Grigorchuk Group ; How do I check that a word over the Grigorchuk group alphabet is trivial?
  if( TheGrigorchukGroupAlgorithms::trivial( w ) )
    cout << w << "  is trivial" << endl;
  else
    cout << w << "  is not trivial" << endl;


  
  //& The Grigorchuk Group ; How do I check if two words over the Grigorchuk group alphabet are conjugate?
  Word w2 = TheGrigorchukGroupAlgorithms::randomWord( L );
  if( TheGrigorchukGroupAlgorithms::conjugate( w , w2 ) ) {
    cout << "Words are conjugate" << endl;
  } else {
    cout << "Words are not conjugate" << endl;
  }
  
  
  //& The Grigorchuk Group ; How do I check if two words over the Grigorchuk group alphabet are conjugate and if so find a conjuator?
  set< Word > conjugators = TheGrigorchukGroupAlgorithms::findConjugator_Kcosets( w , w2 );
  if( conjugators.size( )!=0 ) {
    cout << "Words are conjugate" << endl;
    cout << "A few conjugators are:" << endl;
    for( set< Word >::const_iterator c_it=conjugators.begin( ) ; c_it!=conjugators.end( ) ; ++c_it )
      cout << "  " << *c_it << endl;
  } else {
    cout << "Words are not conjugate" << endl;
  }
  
  
  return 0;
}
