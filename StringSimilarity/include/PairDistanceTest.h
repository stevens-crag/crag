
// Contents: 
//           
// Principal Author: Copyright (2005) Alexei Miasnikov 
//
// Status:
//
// Revision History:
//
// 

#include "Word.h"
#include "RanlibCPP.h"


#ifndef PAIR_DIST_TEST_H
#define PAIR_DIST_TEST_H

#include "errormsgs.h"
#include "SimilarityMeasures.h"

#include <strstream>
#include <vector>
#include <algorithm>
#include <sstream>

//////////////////////////////////////////////////////////////////////////
//
//
//  PairGenerators
//
//
/////////////////////////////////////////////////////////////////////////

//! Abstract interface for classes defining sword  pairs similarities.
class PairGenerator
{
public:
  //! Returns a pair of words  similar under specified criteria.
  virtual pair<Word,Word> getTruePair() = 0;
  //! Returns a pair of words that are not similar under specified criteria.
  virtual pair<Word,Word> getFalsePair() = 0;
};



//! Implements a class for generating pseudo-random pairs of words.
class RandomPairGenerator : public PairGenerator
{
public:
  //! Constructor.
  /*!
    \param N - number of generators
    \param len - the length of each word
   */
  RandomPairGenerator(int N, int len ):
    nGens( N ),
    Len( len ) {}
  
  //! Returns a pair of pseudo-random, independently generated words.
  pair<Word,Word> getTruePair(){
    pair <Word,Word> p;
    p.first = Word::randomWord( nGens,Len );
    p.second = Word::randomWord( nGens,Len ); 

    return p;
  }

  //! Returns a pair of correlated  words 
  /*!
    There are many ways to create correlated words. This is just for testing.
    First word \c w1 is generated randomly and the second \c w2 = w1. Thereofre
    they are equal.
   */
  pair<Word,Word> getFalsePair(){
    pair <Word,Word> p;
    p.first = Word::randomWord( nGens,Len );
    p.second = p.first; //-p.first*Word::randomWord( nGens,Len )*p.first; 

    return p;

  }
private:
  int nGens;
  int Len;
};

//////////////////////////////////////////////////////////////////////////
//
//
//  TEST ON THE SIMILARITY OF PAIRS
//
//
/////////////////////////////////////////////////////////////////////////

//! Implements a hypothesis testing that two words a similar according to some given criteria.
class PairDistanceSimilarityTest
{
public:
  //! Constructor. 
  PairDistanceSimilarityTest( int n,PairGenerator* g, StringSimilarityMeasure* sm) : 
    sampleSize( n ),
    measureDistr( n ),
    theGenerator( g ),
    theSimilarity( sm ){}
  
  void   estimateTrueDistribution();
  double testTruePair();
  double testFalsePair();
  
  double testPair(const pair<Word,Word>& p);
  
private:
  // Vector of measure values
  vector<double> measureDistr;
  StringSimilarityMeasure* theSimilarity;

  int sampleSize;
  PairGenerator* theGenerator;
};

#endif
