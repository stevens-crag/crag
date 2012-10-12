// Copyright (C) 2005 Alexander Ushakov
// Contents: Example for class LinkedBraidStructure and DehornoyForm
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//

#include "Word.h"

#include "iostream"
using namespace std;


//---------------------------------------------------------------------------//
//------------------------- Examples: Word ----------------------------------//
//---------------------------------------------------------------------------//

int main( )
{
  // Fix the number of strands and the length of a braid word
  int N = 5;
  int L = 20;
  
  
  //& Braid word : How do I generate a random freely reduced braid word from B_N of length L?
  Word w1 = Word::randomWord( N-1 , L );
  Word w2 = Word::randomWord( N-1 , L );
  
  
  //& Braid word : How do I invert a braid word?
  Word u1 = -w1;
  

  //& Braid word : How do I multiply braid words?
  Word p = w1 * w2 * u1 * -w2;    // here p = [w1,w2]; the result of the product is reduced whenever each factor is
  
  
  return 0;
}

