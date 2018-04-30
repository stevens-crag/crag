
#include "braid_group.h"
#include "ThRightNormalForm.h"

#include <iterator>
#include <iostream>
#include <fstream>
#include "ShortBraidForm.h"

#ifndef AAG_Challenge
#define AAG_Challenge

using namespace std;

#include "LengthAttack.h"
#include "AAGKeyGeneration.h"

namespace AAGChallenge
{
  
  //
  // RETURNS GENERATORS of Dp obtained from generators of two free factors in genComp
  //                    and relators in Prelts. 
  //
  vector<Word> getDpSubgroup( int N, int Pgens, 
			      const vector<Word>& Prelts, 
			      const pair< vector<Word>, vector<Word> >& genComps, int conjLen );
  
  //
  //
  // RETURNS i'th delta obtained from tubes of size c
  //
  //
  Word getDelta( int i,int c );
  
  //
  //
  //  CREATES GENERATORS of the free factors
  //
  //
  pair< vector<Word>,vector<Word> >  getSgGenComponentsRandom( int Pgens, int c, int k, int len );
  pair< vector<Word>,vector<Word> >  getSgGenComponentsSquares( int Pgens, int c, int k );

  
  //
  //
  //  GENERATES a random word from the subgroup sg in 
  //            in generators of Bn
  //
  //
  Word randomSubgroupWord( int N,const vector<Word>& sg );
  
  //
  //
  //  GENERATES word Wn
  //
  //
  Word specialSubgroupWord( const Word& a, const Word& t, int n );
  Word specialSubgroupWordRandom( const Word& a, const Word& t, int len );

  
  //---------------------------------------------------------------------------//
  //---------------------------------- generateSubgroup -----------------------//
  //---------------------------------------------------------------------------//
  
  vector<Word> generateSubgroup( int c, int k, int comp_len  );
  
  Word generateKeyDecomp(  int n );

}

#endif
