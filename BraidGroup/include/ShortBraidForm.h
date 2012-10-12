
#ifndef _ShortBraidForm_h_
#define _ShortBraidForm_h_

#include <vector>
using namespace std;

class Word;
class LinkedBraidStructure;

Word shortenBraid( int N , const Word& w );

Word shortBraidForm( int N , const Word& w );

vector< Word > shortBraidSbgpForm( int N , const vector< Word >& w );

LinkedBraidStructure shortenLBS( LinkedBraidStructure& lbs );


#endif
