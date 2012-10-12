// Copyright (C) 2005 Alexander Ushakov
// Contents: Definition of short braid form
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//

#ifndef _ShortBraidForm_h_
#define _ShortBraidForm_h_

#include <vector>
using namespace std;

class Word;

Word shortenBraid( int N , const Word& w );

Word shortBraidForm( int N , const Word& w );

vector< Word > shortBraidSbgpForm( int N , const vector< Word >& w );

#endif
