// Copyright (C) 2005 Alexander Ushakov
// Contents: Example for class FreeGroup
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//

#include "FreeGroup.h"
#include "Alphabet.h"

int main( )
{

  int R = 10;

  //& Free Group ; How do I create a  free group of rank R
  FreeGroup F1( R );

  //& Free Group ; How do I create a free group from an alphabet
  FiniteAlphabet a( R );
  
  FreeGroup F2( a );

  //& Free Group; How do I print a free group
  cout << F2 << endl;

  //& Free Group; How do I input a free group from a  stream
  cout << "Enter a presentation of a free group (Example: <a,b,c,d>) : "; cin >> F1;
  cout << F1 << endl;
  
  return 0;
}
