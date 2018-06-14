
#ifndef _ShortBraidForm_h_
#define _ShortBraidForm_h_

#include <vector>
using namespace std;

class Word;
class LinkedBraidStructure;

//! Attempt to reduce |w| using Dehornoy handle free form
Word shortenBraid(int N, const Word &w);
//! Attempt to reduce |w| using Dehornoy handle free form (for longer words)
Word shortenBraid2(int n, const Word &w);
//! Compute normal form for w and then shorten the result
Word shortBraidForm( int N , const Word& w );

Word dehornoy(int N, const Word& w);

Word garsideDehornoy(int N, const Word& w);

vector<Word> shortBraidSbgpForm(int N, const vector<Word> &w);

LinkedBraidStructure shortenLBS(LinkedBraidStructure &lbs);

#endif
