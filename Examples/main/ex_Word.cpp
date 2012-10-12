// Copyright (C) 2006 Aleksey Myasnikov
// Contents: Example for class Word
//
// Principal Authors: Aleksey Myasnikov
//
// Revision History:
//


#include "Word.h"
#include <iostream>
#include <algorithm>
#include <iterator>
#include <vector>
#include <list>

using namespace std;



int main()
{
  //---------------------------------------------------------------------------//
  //------------------------- Examples: Word ----------------------------------//
  //---------------------------------------------------------------------------//
  
  /*
    Reduced words are  represented  as a list of non-trivial integers list< int >.
    Each generator  x is represented by the unique number  n_x.
    For each x  it is assumed that n_{x^{-1}} = -n_x.
  */
  
  
  //& Word; How do I create  a word?
  
  // Create an empty word (identity)
  Word e;
  cout << "Create an empty word : " << e << endl;





  // Create a word from a vector of generators  <x1, ..., xk>, 
  // where xi is a positive or negative  integer s.t. 
  // abs(xi) = 1, ..., n and n is the rank of a group
  vector<int> v(4,0); // 4 generators in the vector
  
  // Enter <1,1,-2,-2> which corresponds to a word x1^2 x2^-2 in <x1,x2>
  v[0] = 1; v[1] = 1; v[2] = -2; v[3] = -2; 
  





  // create a word from v
  Word wv(v);
  cout << "Create a word from a vector of generators : " << wv << endl;
  
  // Create a word from a list of generators  <g_1, ..., g_k>, 
  list<int> l; // 4 generators in the vector
  
  // Create a list which corresponds to a word (x1 x2)^10 in <x1,x2>
  for (int i=0;i<10;i++){
    l.push_back(1); // append the first generator  (see below)
    l.push_back(2); // append the second generator
  }
  
  // create a word from l
  Word wl(l);
  cout << "Create a word from a list of generators : " << wl << endl;
  
  




  //& Word : How do I modify words; How do I add letters  and subwords to a word
  

  // Multiply the word by a one-letter word defined by gen on the right. The result is being reduced.
  Word w_p = e;
  e.push_back(1);
  cout << "Concatenate " << w_p << " and x1 = " << e << endl;



  // Multiply the word by a one-letter word defined by gen on the left. The result is being reduced.
  w_p = e;
  e.push_front(2);
  cout << "Concatenate x2 and " << w_p << " = " << e << endl;



  // Multiply the word by a word on the right. The result is being reduced.
  w_p = e;
  e.push_back(wv);
  cout << "Concatenate " <<  w_p << " and  " << wv << " = " << e << endl;



  // Multiply the word by a word on the left. The result is being reduced.
  w_p = e;
  e.push_front(wv);
  cout << "Concatenate " <<  wv << " and  " << w_p << " = " << e << endl;
  




  
  //& Word : How do I modify words; How do I insert letters or subwords into a word


  // Insert a sequence of generators [B,E) into a word at th position 5
  Word w_tmp = wl;
  //  w_tmp.insert<ConstWordIterator>( 5 , wv.begin(),wv.end() );
  //  cout << "Insert :" << wl << " -> " << w_tmp << endl;
  



  // Insert a generator x5 into a word after the 5th letter
  w_tmp = wl;
  w_tmp.insert( 5,5 );
  cout <<  "Insert :" << wl << " -> " << w_tmp << endl;
  



  // Insert a sequence of generators [B,E) into a word before the second position
  w_tmp = wl;
  WordIterator wI = wl.begin();
  //  w_tmp.insert<list<int>::const_iterator>( ++wI,wv.getList().begin(),wv.getList().end() );
  //cout <<  "Insert :" << wl << " -> " << w_tmp << endl;
  



  // Insert a generator x5 into a word before the second position
  w_tmp = wl;
  wI = w_tmp.begin();
  w_tmp.insert( ++wI,5 );
  cout <<  "Insert :" << wl << " -> " << w_tmp << endl;
  
  


  //& Word : How do I modify words; How do I replacing letters and subwords
  

  //Replace a generator at the second position by x10
  w_tmp = wl;
  wI = w_tmp.begin();
  w_tmp.replace( ++wI,5 );
  cout <<  "Replace :" << wl << " -> " << w_tmp << endl;


  
  // Replace a subword of a word starting at the second position by a word [B,E). 
  // The length of the word does not increase if [B,E) is longer than the terminal
  // segment of the word [it,end()). 
  //In that case terminal symbols of [B,E) are ignored.
  w_tmp = wl;
  wI = w_tmp.begin();
  //  w_tmp.replace<list<int>::const_iterator>( ++wI,wv.getList().begin(), wv.getList().end() );
  // cout << "Replace :" << wl << " -> " << w_tmp << endl;
  
  


  //& Word : How do I modify words; How do I multiply words
  

  // Multiply two words. The result is reduced.
  cout << "Multiply :" << wv << " * " << wl << " = " << wv*wl << endl;
  
  // Multiply a word on the right by another word. The result is reduced.
  Word wmul(wv);
  wmul *= wl;
  cout << "Multiply : " << wmul << endl;;
  
  


  //& Word; How to invert a word
  cout << "Inverse of " << wv << " = " << -wv << endl;
  cout << "Inverse of " << wv << " = " << wv.inverse() << endl;
  
  

  //& Word : How do I modify words; How to reduce words
  
  // Freely reduce a word.
  wv.freelyReduce( );
  
  // Freely reduce a segment of a word defined by [B,E).
  wv.freelyReduce( wv.begin(),wv.end() );
  
  // Cyclically reduce a word
  Word cw = wv*wl*-wv;
  cw.cyclicallyReduceWord();
  cout << "Cyclically reduce " << wv*wl*-wv << " = " << cw << endl;
  
  // Cyclically reduces a word and returns the corresponding conjugator.
  Word conjugator;
  cw = wv*wl*-wv;
  cw.cyclicallyReduceWord(conjugator);
  cout << "Cyclically reduce " << wv*wl*-wv << " = " << cw << endl;
  cout << "Conjugator " << conjugator << endl;
  
  // Returns the cyclically reduced word. The original word is unchanged
  cout << "Cyclically reduce " << wv*wl*-wv << " = " << (wv*wl*-wv).cyclicallyReduce(  ) << endl;
  
  // Returns the cyclically reduced word and the corresponding conjugator. The original word is unchanged
  cout << "Cyclically reduce " << wv*wl*-wv << " = " << (wv*wl*-wv).cyclicallyReduce( conjugator ) << endl;
  cout << "Conjugator " << conjugator << endl;
  
  


  //& Word;  How do I iterate through the letters of a word
  
  // Iterate and print generators  of the word. 
  cout << "Generators of " << wv << " : " << endl;
  for (ConstWordIterator I=wv.begin();I!=wv.end();I++)
    cout << *I << " ";
  cout << endl;
  


  // Get a constant representation of a word as a list of integers corresponding to the generators.
  list< int > l1 = wv.getList( );
  cout << "List representation of " << wv << " is ";
  copy(l1.begin(),l1.end(), ostream_iterator<int>( cout," "));
  cout << endl;
  

  // Get a representation of a word as a list of integers corresponding to the generators.
  // This allows direct manipulation with the representation (which requires caution).
  cout << "Change second letter of " << wv << " -> ";
  list< int >& l2 = wv.getList( );
  // change second letter to x10
  list<int>::iterator lI = l2.begin();
  *(++lI) = 10;
  cout << wv << endl;
  
  
  
  //& Word; How do I get words' length
  cout << "Length of " << wv << " is " << wv.length() << endl;
  


  //& Word; How to extract subwords
  
  // Get an initial segment of the word of length 2
  cout << "Initial segment of " << wv << " is " << wv.initialSegment( 2 ) << endl;
  
  // Get a terminal segment of the word of length 2.
  cout << "Terminal segment of " << wv << " is " << wv.terminalSegment( wv.length() - 2 ) << endl;
  
  // Get a segment of the word defined 
  cout << "Segment of " << wv << " is " << wv.segment( 1,3 ) << endl;
  
  
  //& Word; How do I compare words
  
  // Words are compared using lexicographical order 
  if (wv < wl )
    cout << wv << " is less than " << wl << endl;
  
  if (wv > wl )
    cout << wv << " is greater than " << wl << endl;
  
  if (wv == wl )
    cout << wv << " is equal to " << wl << endl;
  
  if (wv != wl )
    cout << wv << " is not equal to " << wl << endl;
  



  //& Word: Operations with cyclic  words; How do I cyclically shift a word

  // Shifts the word one position to the left
  Word wvls = wv;
  wvls.cyclicLeftShift();
  cout << "Cyclic left shift : " << wv << " -> " << wvls  << endl;


  // Shifts the word one position to the right
  Word wvrs = wv;
  wvrs.cyclicRightShift();  
  cout  << "Cyclic right shift : " << wv << " -> " <<  wvrs  << endl;
  
  



  //& Word: Operations with cyclic  words; How do I cyclically permute a word


  //Cyclically left-shift permute the word by 5 positions.
  cout << "Cyclically left-shift permute : "  << wv << " -> " << wv.cyclicallyPermute( 5 ) << endl;



  
  //Cyclically righ-shift permute the word by 5 positions.
  cout << "Cyclically right-shift permute : "  << wv << " -> " << wv.cyclicallyPermute( -5 ) << endl;
  

    
  //& Word; How to generate a pseudo-random word
  
  // Generate a pseudo randomly reduced word of the length 10 in 2 generators
  cout << "Generate pseudo-random word of length 10 " <<  Word::randomWord( 2 , 10 ) << endl;
  
  // Generates a pseudo randomly reduced word of a length in [10,15] and  2 generators
  cout << "Generate pseudo-random word of length in [10,15] " <<  Word::randomWord( 2,10,15 ) << endl;
  


  //& Word; How do I print a word

  // Print a word into the standard output stream. 
  // Default operator will print a word in alphabet x1 ... xN
  // See Alphabet class description too change the output alphabet 
  cout << "Print word : " << wv << endl;
  

  //& Word; How do I  input  a word
  
  // Read a word from the standard input stream. 
  // Default operator will accept a word in alphabet x1 ... xN
  // See Alphabet class description too read words in arbitrary alphabet 
  cout << "Enter a word in generators x1, x2 ... " << endl
       << "Example: x1^2 [x2,x3^-2] (x5 x10)^-1 ; " << flush; 
  cin >> wv;
  cout << "You entered " << wv << endl;
  
  

  //& Word; How to determine the power of the word as an element of a free monoid
  // Determines the power of the word (as an element of a free monoid, not as an element of a free group).
  Word base;
  cout << "Power of " << wv << " is " << wv.getPower( base ) << endl;
  cout << "Base : " << base << endl;
    
  
  //& Word; How to check if a word contains a generator x2
  Word wr = Word::randomWord(5,5);
  cout << wr << ( wr.doesContain( 2 ) ? " contains " : " doesnot contain " ) << " x2 " << endl;
  
  
  //& Word; How do I compute  the  exponent sum of a generator
  cout << "Exponent sum of x2 in " << wr << " is " << wr.exponentSum( 2 ) << endl;
  
  //& Word; How do I compute a power of a word 
  // Compute wr^2
  cout << wr << " squared = " << wr.power(2) << endl;
  
  //& Word;  How to get the number of times a generator occurs in a word  
  // Print how many times x2 occur in wr 
  cout << "x2 occurs " << wr.isIn( 2 )  << " times in " << wr << endl;
  
  
  //! Compute the minimal equivalent word. 
  /*! permutableGenerators must be positive.
   */
  //Word minimalEquivalentForm( const set< int >& permutableGenerators , bool inverses , bool cyclicPermutations ) const;
  
  
}
