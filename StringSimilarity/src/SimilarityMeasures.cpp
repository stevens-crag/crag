// Contents: 
//
// Principal Author: Alexei Miasnikov (2004)
//
// Status:
//
// Revision History:
//

#include "Word.h"
#include "SimilarityMeasures.h"
#include "Levenstein.h"
#include <limits.h>
#include <algorithm>

const double ALMOST_ZERO = 0.00000001;

double getLeftTaleConfidenceValue(const vector<double>& distr, double value)
{
  int p_value = 0;
  for (int i=0;i<distr.size();i++)
    if (value > distr[i]){
      p_value = i;
    }
  
  return double(p_value) / double( distr.size());
  
}


vector<double> WordPairComparison::distrEstimate(int minLen, int maxLen, int nSamples)const
{
  // Vector of measure values
  vector<double> measureDistr(nSamples);
  
  // generate "number" of pairs
  for ( int lCount=0;lCount<nSamples;lCount++ ){
    
    // Randomly generate a pair
    Word w1 = Word::randomWord( theRank,  minLen,maxLen );
    Word w2 = Word::randomWord( theRank,  minLen,maxLen );
    
    measureDistr[lCount] = pSSM->measure( w1,w2 );
    
  }
  
  sort(measureDistr.begin(),measureDistr.end());

  return measureDistr;

}


double WordPairComparison::comparePair( const Word& w1, const Word& w2, const vector<double>& measureDistr)const
{
  cout << "Check if it is  proper!!!!!" << endl;
  return  getLeftTaleConfidenceValue(measureDistr,pSSM->measure( w1,w2 ));
  //  return  1.0 - getLeftTaleConfidenceValue(measureDistr,pSSM->measure( w1,w2 ));
}

double WordPairComparison::comparePair( const Word& w1, const Word& w2)const
{
   return  pSSM->measure( w1,w2 );
}


/////////////////////////////////////////////////////
//
//
//  MEASURES
//
//
/////////////////////////////////////////////////////



double HammingDistance::measure(const Word& w1, const Word& w2) const
{
  int minmalLength = min(w1.length(),w2.length());
  int distance = w1.length()-minmalLength+w2.length()-minmalLength;
  for (ConstWordIterator I1 = w1.begin(), I2 = w2.begin(); I1!=w1.end() && I2!=w2.end();I1++,I2++)
    if (*I1 != *I2)
      distance++;

  return double(distance) / double(max(w1.length(),w2.length()));
}

double HammingDistanceCyclic::measure(const Word& w1, const Word& w2) const
{

  HammingDistance hd;
  double min_distance = INT_MAX;
  Word w1cr = w1.cyclicallyReduce();
  Word w2cr = w2.cyclicallyReduce();
  
  Word minWord;
  Word maxWord;
  if ( w1cr.length() < w2cr.length() ){
    minWord = w1cr;
    maxWord  = w2cr;
  } else {
    minWord = w2cr;
    maxWord  = w1cr;   
  }
  

  for (int i=0;i<maxWord.length();i++){
    maxWord.cyclicLeftShift();
    double d = hd.measure(maxWord, minWord);
    //cout << d << endl;
    if ( d < min_distance )
      min_distance = d;
  }     
  //  return (min_distance - abs( w1.length() - w2.length()))/ double(min(w1.length(),w2.length()));
  return min_distance;
}

double SubwordHammingDistanceCyclic::measure(const Word& w1, const Word& w2) const
{
  HammingDistance hd;
  Word w1cr = w1.cyclicallyReduce();
  Word w2cr = w2.cyclicallyReduce();
  
  Word minWord;
  Word maxWord;
  if ( w1cr.length() < w2cr.length() ){
    minWord = w1cr;
    maxWord  = w2cr;
  } else {
    minWord = w2cr;
    maxWord  = w1cr;   
  }
  
  double min_distance = INT_MAX;
  for (int i=0;i<maxWord.length();i++){
    maxWord.cyclicLeftShift();
    double d = hd.measure(maxWord.initialSegment(min(maxWord.length()-1,minWord.length()-1)), minWord);
    if ( d < min_distance )
      min_distance = d;
  }     
  return min_distance;
};


double EditingDistance::measure(const Word& w1, const Word& w2) const
{
  Levenstein ld;
  //return double(ld.compute(w1,w2));
  return double(ld.compute(w1,w2)) / double(max(w1.length(),w2.length()));
};


double SubwordEditingDistanceCyclic::measure(const Word& w1, const Word& w2) const
{

  EditingDistance ed;
  Word w1cr = w1.cyclicallyReduce();
  Word w2cr = w2.cyclicallyReduce();

  Word minWord;
  Word maxWord;
  if ( w1cr.length() < w2cr.length() ){
    minWord = w1cr;
    maxWord  = w2cr;
  } else {
    minWord = w2cr;
    maxWord  = w1cr;   
  }
  

  double min_distance = INT_MAX;
  for (int i=0;i<maxWord.length();i++){
    maxWord.cyclicLeftShift();
    double d = ed.measure(maxWord.initialSegment(min(maxWord.length()-1,minWord.length()-1)), minWord);
    if ( d < min_distance )
      min_distance = d;
  }     
  return min_distance;
}


/*
double WhiteheadGraphSimilarity::measure(const Word& w1, const Word& w2) const
{
  WhiteheadSimpleGraph wg1( w1,theGroup );
  WhiteheadSimpleGraph wg2( w2,theGroup );
  
  vector<double> wFeatures1 = wg1.getWeightVector();
  vector<double> wFeatures2 = wg2.getWeightVector();

  // NORMALIZE??????
  for (int i=0;i<wFeatures1.size();i++){
    wFeatures1[i] /= double(wFeatures1.size());
    wFeatures2[i] /= double(wFeatures2.size());
  }
  
*/
  /* Computes Euclidean distance
  double fsum = 0;
  for (int i=0;i<wFeatures1.size();i++){
    fsum += (wFeatures1[i] - wFeatures2[i])*(wFeatures1[i] - wFeatures2[i]);
  }
  return sqrt(fsum);
  */
  /*
  // Compute Kullback-Leibler distance
  vector<double>&  minWord = wFeatures1;
  vector<double>&  maxWord = wFeatures2;
  if (w1.length() > w2.length()){
    minWord = wFeatures2;
    maxWord = wFeatures1;
  }
  
  double kl_sum = 0;
  for (int i=0;i<maxWord.size();i++){
    double dMinWordi = (minWord[i] == 0) ? ALMOST_ZERO : double(minWord[i]);
    double dMaxWordi = (maxWord[i] == 0) ? ALMOST_ZERO : double(maxWord[i]);
   
    kl_sum += double(dMinWordi) * log(double(dMinWordi) / double(dMaxWordi));
    //    cout << kl_sum << " " << double(dMinWordi) << " " << double(dMaxWordi) << endl;
  }
  
  return kl_sum;
  
}


double MotiveSimilarity::measure(const Word& w1, const Word& w2) const
{
  MotivePatternWrapper mpw( theGroup );
  
  vector<Word> seq(2);
  seq[0] = w1;
  seq[1] = w2;
  mpw.initializeSequence( seq );
  mpw.compute();
  vector<MotivePattern> ps = mpw.getPatterns();

  int maxPatternLength = 0;

  for ( int i=0;i<ps.size();i++)
    if (ps[i].sequences >= 2 && ps[i].thePattern.length() > maxPatternLength )
      maxPatternLength = ps[i].thePattern.length();

  //  cout << "Motive pattern size : " << maxPatternLength << endl;
  return maxPatternLength;
}
*/
