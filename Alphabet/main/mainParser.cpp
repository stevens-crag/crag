#include <iostream>
#include "Alphabet.h"
#include <iterator>
#include "Word.h"

int main()
{
  
  //readFPPresentation( cin );
  
  FiniteAlphabet a;
  cin  >> a;
    
  while (1) {
        cout << endl << "Enter a vector of words in the alphabet " << a << ": ";
	//    vector<Word> w = a.readVector( cin );
	Word w = a.readWord( cin );
        cout << endl;
        //a.printVector( cout,w );
	a.printWord( cout, w);

	//cout << endl << "Enter a word in the default alphabet " 
	//	 << InfiniteAlphabet::defaultAlphabet << ": ";
	//    vector<Word> v = InfiniteAlphabet::defaultAlphabet.readVector( cin );
	//    InfiniteAlphabet::defaultAlphabet.printVector(cout,v); cout << endl;
    
  }

  //  copy( p.getWord().begin(), p.getWord().end(),ostream_iterator<int>(cout, " ")); 

  return 0;
}
