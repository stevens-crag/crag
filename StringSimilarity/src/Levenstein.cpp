// Copyright (C) 2003 Dmitry Bormotov
// 
// Contents: Implementation of class Levenstein
//
// Principal Author: Dmitry Bormotov
//
// Status: in progress
//
// Revision History:
//
//


#include "Levenstein.h"
#include <string.h>


//#define DEBUG_LEVENSTEIN


// ------------------------ Levenstein -------------------------- //


int    addition = 1;
int    change   = 1;
int    deletion = 1;

#define COMP_LEN     10000
#define ARR_SIZE     COMP_LEN + 1
#define SMALLEST_OF(x,y,z) ( (x<y) ? min(x,z) : min(y,z) )
#define ZERO_IF_EQUAL(x,y) (requested[x-1] == found[y-1] ? 0 : change)

int    distnce[ARR_SIZE][ARR_SIZE];

int ldistance( char* requested, char* found)
{
   int i,j;
   int r_len, f_len;

#if SAFETY > 0
   if( strlen(requested)>COMP_LEN || strlen(found)>COMP_LEN )
     cout << "*** The length is longer than expected! ***" << endl;
#endif  

   r_len = (strlen(requested)>COMP_LEN ? COMP_LEN : strlen(requested));
   f_len = (strlen(found)>COMP_LEN ? COMP_LEN : strlen(found));
   
   distnce[0][0] = 0;
   for (j = 1; j <= ARR_SIZE; j++)
      distnce[0][j] = distnce[0][j-1] + addition;
   for (j = 1; j <= ARR_SIZE; j++)
      distnce[j][0] = distnce[j-1][0] + deletion;

   for (i = 1; i <= r_len; i++)
      for (j = 1; j <= f_len; j++)
         distnce[i][j] = SMALLEST_OF(
	     (distnce[i-1][j-1] + ZERO_IF_EQUAL(i,j)),
	     (distnce[i][j-1]   + addition),
	     (distnce[i-1][j]   + deletion) );

    return( distnce[r_len][f_len] );
}

int Levenstein::compute( const Word& u, const Word& v ) {
  
  char* s = wordToString(u);
  char* t = wordToString(v);
  int d = ldistance(s,t);
  delete [] s;
  delete [] t;
  return d;
}

char* Levenstein::wordToString( const Word& w ) {

  int wLen = w.length();
  char* s = new char[wLen+1];

  {
    int i=0;
    for( ConstWordIterator I = w.begin(); I != w.end();I++,i++ ) {
      int o = *I;
      if( o > 0 )
	s[i] = 2*o;
      else
	s[i] = -2*o+1;
    }
  }
  s[wLen] = 0;

#ifdef DEBUG_LEVENSTEIN
  for( int i = 0; i < wLen; ++i )
    cout << " " << int(s[i]);
  cout << endl;
#endif

  return s;
}
