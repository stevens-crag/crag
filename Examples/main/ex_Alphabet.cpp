// Copyright (C) 2005 Aleksey Myasnikov
// Contents: Example for class FiniteAlphabet
//
// Principal Authors: Aleksey Myasnikov
//
// Revision History:
//

#include "Alphabet.h"
#include "Word.h"
#include <vector>

int main( )
{
  int R = 3;

  //& Alphabet ; How do I create a finite alphabet of size R
  FiniteAlphabet A1( R );

  //& Alphabet ; How do I create a free group from a list of letters
  vector<string> letters(R);
  cout << " Enter " << R << " letters : " << flush;  
  for (int i=0;i<R;i++)
    cin >> letters[i];
  
  FiniteAlphabet A2( letters );
  
  //& Alphabet; How do I print an alphabet
  cout << A2 << endl;

  //& Alphabet; How do I input an alphabet from a stream
  cout << "Enter a finite alphabet (Example: {a,b,c,d}) : "; cin >> A1;
  cout << A1 << endl;
  

  //& Alphabet; How do I input a word in the alphabet
  cout << "Input a word in " << A1 << " ending with ';' " <<  flush;
  Word w1 = A1.readWord( cin );

  //& Alphabet; How do I print word in letters of  the alphabet
  A1.printWord( cout, w1 ); cout << endl;


  
  //& Alphabet; How do I input a vector of words in the alphabet
  cout << "Input a vector of words in " << A1 << " Example : { a^2 b , a b } " <<  flush;
  vector<Word> v1 = A1.readVector( cin );

  //& Alphabet; How do I print a vector of words in letters of  the alphabet
  A1.printVector( cout, v1 ); cout << endl;

  
  return 0;
}
