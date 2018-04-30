// Copyright (C) 2005 Alexander Ushakov
// Contents: Example for Thurston's left normal form of braids.
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//

#include "Word.h"
#include "LinkedBraidStructure.h"
#include "braid_group.h"
#include "ThLeftNormalForm.h"
#include "ShortBraidForm.h"

#include "iostream"
using namespace std;


//---------------------------------------------------------------------------//
//---------------------- Examples: ThLeftNormalForm ------------------------//
//---------------------------------------------------------------------------//


int main( )
{
  // Fix the number of strands and the length of a braid word
  int N = 5;
  int L = 20;
  
  // Generate two random braid words from B_N of length L
  Word w1 = Word::randomWord( N-1 , L );
  Word w2 = Word::randomWord( N-1 , L );
  
  
  //& Thurston's left normal form : How do I create Thurston's left normal form of a braid word?
  BraidGroup B( N );                    // a) create an object "braid group" representing B_N
  ThLeftNormalForm nf1( B , w1 );      // b) compute the normal form of w1
  ThLeftNormalForm nf2( B , w2 );      // b) compute the normal form of w2
  
  
  //& Thurston's left normal form : How do I output Thurston's left normal form?
  cout << "NF1 = " << nf1 << endl;
  cout << "NF2 = " << nf2 << endl;


  //& Thurston's left normal form : How do I check that Thurston's left normal form is trivial?
  //& Braid word : How do I check that a braid word represents trivial braid (using Thurston's left normal forms)?
  if( nf1.isTrivial( ) )
    cout << "Braid word is trivial" << endl;
  else
    cout << "Braid word is not trivial" << endl;
  
  
  
  //& Thurston's left normal form : How do I check equality of two Thurston's left normal forms?
  if( nf1==nf2 )
    cout << "Words w1 and w2 represent the same element of B" << endl;
  else 
    cout << "Words w1 and w2 represent different elements of B" << endl;
  

  //& Thurston's left normal form : How can I get a word (some word) represented by a Thurston's left normal form?
  Word nfw  = nf1.getWord( );
 
  
  //& Thurston's left normal form : How do I multiply two Thurston's left normal forms?
  ThLeftNormalForm nf3 = nf1 * nf2;


  //& Thurston's left normal form : How do I invert a Thurston's left normal form?
  ThLeftNormalForm nf4 = -nf1;
  

  //& Thurston's left normal form : How do I compute a representative of the super summit set for a Thurston's left normal form?
  pair< ThLeftNormalForm , ThLeftNormalForm > sss_res = nf1.findSSSRepresentative( );
  ThLeftNormalForm representative = sss_res.first;
  ThLeftNormalForm conjugator = sss_res.second;
  // Now representative = a representative of the super summit set of nf1
  // Also, representative = -conjugator * nf1 * conjugator.
  
  
  
  //& Thurston's left normal form : How do I check if two Thurston's left normal forms are conjugate using Super Summit Set construction?
  pair<bool,ThLeftNormalForm> conj = nf1.areConjugate( nf2 );
  // if nf1 and nf2 are not conjugate then conj.first = false
  // if nf1 and nf2 are conjugate then conj.first = true and conj.second is an actual conjugator
  


  //& Thurston's left normal form : How do I check if two Thurston's left normal forms are conjugate using Ultra Summit Set construction?
  // IN GENERAL USS-CONSTRUCTION IS MUCH MORE EFFICIENT THAN SSS-COSNTRUCTION
  try{
    pair<bool,ThLeftNormalForm> conj_uss = nf1.areConjugate_uss( nf2 );
    // if nf1 and nf2 are not conjugate then conj_uss.first = false
    // if nf1 and nf2 are conjugate then conj_uss.first = true and conj_uss.second is an actual conjugator
    if( conj_uss.first ) {
      cout << "Braids are conjugate" << endl;
    } else {
      cout << "Braids are not conjugate" << endl;
    }
  } 
  catch( const exception& e ) {
    // Getting into this block means that there is not enough memory available to perform the conjugacy test
    cout << "Got an exception: " << e.what( ) << endl;
  }
  
  
  return 0;
}
