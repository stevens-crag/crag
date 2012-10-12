// Copyright (C) 2005 Alexander Ushakov
// Contents: Example for class FreeGroup
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//

#include "FPGroup.h"
#include "Alphabet.h"

int main( )
{

  FiniteAlphabet a;
  cout << "Enter an alphabet (Example: {a,b,c})" << endl;
  cin >> a;

  cout << "Enter relators (Example: { [a,b] , (c b)^2 })" << endl;
  vector<Word> r = a.readVector(cin);


  //&  Finitely Presented Group  ; How do I create a Finitely Presented Group 
  FPGroup G( a,r );

  //&  Finitely Presented Group ; How do I output a Finitely Presented Group into a stream
  cout << G << endl;

  //&  Finitely Presented Group ; How do I read a Finitely Presented Group from a stream
  cout << "Enter a presentation (Example: <a,b,c | [a,b], a^2 >) " << endl; 
  cin >> G;
  cout << G << endl;
  
  return 0;
}
