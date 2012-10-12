// Contents: 
//
// Principal Author: Alexei Miasnikov (2004)
//
// Status:
//
// Revision History:
//

#include "Word.h"
#include "StringScramblers.h"
#include "RanlibCPP.h"

////////////////////////////////////////////////////////////////////////////
//
// UNIFORM SCRAMBLER
//
///////////////////////////////////////////////////////////////////////////


Word UniformScrambler::scramble(const Word& w)const
{
  double changed = 0.0;
  Word returnWord;
  
  for (ConstWordIterator I = w.begin(); I!=w.end();I++)
    if (RandLib::ur.rand() <= theChangeProb){
      Generator gen = 0;
      do {
	gen = RandLib::ur.irand(-numGens, numGens);
      } while ( gen == 0 || gen == *I );
      returnWord.push_back(gen);
      changed += 1.0;
    } else {
      returnWord.push_back(*I);
    }
  ((UniformScrambler*)this)->theFracChanged = changed / double(w.length());
  
  return returnWord.freelyReduce();
}

////////////////////////////////////////////////////////////////////////////
//
// SUBWORD SCRAMBLER
//
///////////////////////////////////////////////////////////////////////////

// Function requires revision!!! Need more efficient one.
Word SubwordScrambler::scramble(const Word& w)const
{
  double changed = 0.0;
  
  int sw_length = int(double(w.length()) * theChangeFrac);
  int init_pos = RandLib::ur.irand(0, w.length() - sw_length);
  Word subw =  Word::randomWord( numGens , sw_length );

  Word returnWord = w;
    
  WordIterator I = returnWord.begin();
  for (int i=0;i<init_pos;i++) I++;
  returnWord.replace( I , subw.begin( ) , subw.end( ) );

  return returnWord.freelyReduce();

}

////////////////////////////////////////////////////////////////////////////
//
// WORD MULTIPLY SCRAMBLER
//
///////////////////////////////////////////////////////////////////////////


Word WordMultiplyScrambler::scramble(const Word& w)const
{
  int sw_length = int(double(w.length()) * theChangeFrac);
  Word mulw = Word::randomWord( numGens, sw_length );
  Word returnWord = w * mulw;
  return returnWord.freelyReduce();
}
