// Copyright (C) 2005 Alexander Ushakov
// Contents: Example for Thurston's right normal form of braids.
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//

#include "Word.h"
#include "LinkedBraidStructure.h"
#include "braid_group.h"
#include "ThRightNormalForm.h"
#include "ShortBraidForm.h"

#include "iostream"
using namespace std;


//---------------------------------------------------------------------------//
//---------------------- Examples: ThRightNormalForm ------------------------//
//---------------------------------------------------------------------------//


int main( )
{
  // Fix the number of strands and the length of a braid word
  int N = 5;
  int L = 20;
  
  // Generate two random braid words from B_N of length L
  Word w1 = Word::randomWord( N-1 , L );
  Word w2 = Word::randomWord( N-1 , L );
  
  
  //& Thurston's right normal form : How do I create Thurston's right normal form of a braid word?
  BraidGroup B( N );                    // a) create an object "braid group" representing B_N
  ThRightNormalForm nf1( B , w1 );      // b) compute the normal form of w1
  ThRightNormalForm nf2( B , w2 );      // b) compute the normal form of w2
  
  
  //& Thurston's right normal form : How do I output Thurston's right normal form?
  cout << "NF1 = " << nf1 << endl;
  cout << "NF2 = " << nf2 << endl;


  //& Thurston's right normal form : How do I check that Thurston's right normal form is trivial?
  //& Braid word : How do I check that a braid word represents trivial braid (using Thurston's right normal forms)?
  if( nf1.isTrivial( ) )
    cout << "Braid word is trivial" << endl;
  else
    cout << "Braid word is not trivial" << endl;
  
  
  
  //& Thurston's right normal form : How do I check equality of two Thurston's right normal forms?
  if( nf1==nf2 )
    cout << "Words w1 and w2 represent the same element of B" << endl;
  else 
    cout << "Words w1 and w2 represent different elements of B" << endl;
  

  //& Thurston's right normal form : How can I get a word (some word) represented by a Thurston's right normal form?
  // Basically there are two ways to get a word represented by a Thurston's right normal form
  // The second way give shorter words
  Word nfw  = nf1.getWord( );
  Word nfw2 = nf1.getShortWord( );
  
  
  //& Thurston's right normal form : How do I multiply two Thurston's right normal forms?
  ThRightNormalForm nf3 = nf1 * nf2;


  //& Thurston's right normal form : How do I invert a Thurston's right normal form?
  ThRightNormalForm nf4 = -nf1;
  

  //& Thurston's right normal form : How do I compute a representative of the super summit set for a Thurston's right normal form?
  pair< ThRightNormalForm , ThRightNormalForm > sss_res = nf1.findSSSRepresentative( );
  ThRightNormalForm representative = sss_res.first;
  ThRightNormalForm conjugator = sss_res.second;
  // Now representative = a representative of the super summit set of nf1
  // Also, representative = -conjugator * nf1 * conjugator.
  
  
  //& Thurston's right normal form : How do I compute a normalizer of a Thurston's right normal form?
  set<ThRightNormalForm> normalizer = nf1.computeCentralizer( );
  // After execution normalizer contains generators of the normalizer of nf1.
  
  
  //& Thurston's right normal form : How do I check if two Thurston's right normal forms are conjugate?
  pair<bool,ThRightNormalForm> conj = nf1.areConjugate( nf2 );
  // if nf1 and nf2 are not conjugate then conj.first = false
  // if nf1 and nf2 are conjugate then conj.first = true and conj.second is an actual conjugator
  
  
  return 0;
}
