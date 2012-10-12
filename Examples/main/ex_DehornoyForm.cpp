// Copyright (C) 2005 Alexander Ushakov
// Contents: Example for class LinkedBraidStructure and DehornoyForm
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//

#include "Word.h"
#include "LinkedBraidStructure.h"
#include "ShortBraidForm.h"
#include "DehornoyForm.h"

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
  
  
  // Generate random freely reduced braid word from B_N of length L.
  Word w1 = Word::randomWord( N-1 , L );
  Word w2 = Word::randomWord( N-1 , L );


  //& Dehornoy Form : How do I compute Dehornoy Form of a braid word?
  DehornoyForm DF1( N , w1 );          // a) compute dehornoy form
  Word df1 = DF1.getDehornoyForm( );   // b) get corresponding word
  cout << "Dehornoy Form of w1 is " << df1 << endl;
  
  //& Braid word : How do I check that a braid word represents trivial braid (using Dehornoy forms)?
  DehornoyForm DF2( N , w2 );        // a) compute dehornoy form of w2
  Word df2 = DF2.getDehornoyForm( ); // b) get corresponding word
  if( df2.length()==0 )              // c) check if it is trivial
    cout << "Braid word w2 is trivial" << endl;
  else
    cout << "Braid word w2 is not trivial" << endl;
  
  
  return 0;
}
