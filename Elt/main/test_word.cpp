
#include "Word.h"


int main( )
{
	Word w( 2 );
	cout << "w = " << w << endl;

	Word u = Word::randomWord( 3 , 10 );
	cout << "u = " << u << endl;

	for( Word::const_iterator u_it = u.begin( ) ; u_it!=u.end( ) ; ++u_it ) {
		cout << *u_it << endl;
	}

	return 0;
}

