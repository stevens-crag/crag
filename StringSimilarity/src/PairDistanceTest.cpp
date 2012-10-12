
// Contents: 
//           
// Principal Author: Alexei Miasnikov (2004)
//
// Status:
//
// Revision History:
//
// 

#include "Word.h"
#include "RanlibCPP.h"
#include "ConfigFile.h"
#include "errormsgs.h"
#include "SimilarityMeasures.h"
#include "PairDistanceTest.h"
#include "FormatOutput.h"


#include <iostream>
#include <vector>
#include <algorithm>
#include <sstream>
#include <iterator>



void PairDistanceSimilarityTest::estimateTrueDistribution()
{
  // generate pairs and compute corresponding measure values
  for ( int lCount=0;lCount<sampleSize;lCount++ ){
    
    pair<Word,Word> p = theGenerator->getTruePair( );
    measureDistr[lCount] = theSimilarity->measure( p.first,p.second );
    
    cout << PBar( double(lCount+1)/double(sampleSize) ) << " "  << measureDistr[lCount];
  }
  cout << endl;
  
  sort(measureDistr.begin(),measureDistr.end());
  copy( measureDistr.begin(), measureDistr.end(),ostream_iterator<double>(cout," "));
  
}

double PairDistanceSimilarityTest::testTruePair()
{
  pair<Word,Word> p = theGenerator->getTruePair( );
  return testPair( p );
}

double PairDistanceSimilarityTest::testFalsePair()
{
  pair<Word,Word> p = theGenerator->getFalsePair( );
  return testPair( p );
}
  
double PairDistanceSimilarityTest::testPair(const pair<Word,Word>& p)
{
  if (measureDistr.size() == 0)
    msgs::error("PairDistanceSimilarityTest::testPair(...): Need to estiamte true distribution first.");
  return  getLeftTaleConfidenceValue(measureDistr,theSimilarity->measure( p.first,p.second ));
}

