// Contents: Classes implementing similarity measures for strings

//
// Principal Author:Copyright (2005)  Alexei Miasnikov
//
// Status:
//
// Revision History:
//

#ifndef STRING_SIMILARITY_H
#define STRING_SIMILARITY_H


#include "Word.h"
#include <vector>
//#include "MotivePatternWrapper.h"


//! Returns the left tale confidence value.
/*!
  Given an estimate (as a list of sorted samples) of a probability
  distributions, returns the confidence value for a given \c p-value.
  \param dist - distribution estimate.
  \param p - the \c p -value.
  \return the confidence value (probability).
 */
double getLeftTaleConfidenceValue(const vector<double>& dist, double p);

class StringSimilarityMeasure;


//! Implements a probabilistic measure for comparing two words
/*!
  Basic Idea: ...
 */
class WordPairComparison
{
 public:
  //! Constructor. 
  /*!
    \param rank - the rank of a free group
    \param sm   - pointer to the corresponding similarity measure.
  */
  WordPairComparison( int rank,  const StringSimilarityMeasure* sm): theRank( rank ), pSSM( sm ) {}
  
  //! Estimates the distribution of the distances (similarities) between two randomly generated  words.
  /*!
    \param minLen - the minimal length of a randomly generated word
    \param maxLen - the maximal length of a randomly generated word
    \param nSamples - the number of pair samples to be generated (1000 is the default value).
    \return vector of the sorted distances (similarity measures)  
    
   */
  vector<double> distrEstimate(int minLen, int maxLen, int nSamples = 1000)const;
  
  //! Compare the pair of words using probabilistic measure.
  /*!
    \param w1 - the first word
    \param w2 - the second word.
    \param measureDistr - distribution of the distance between two random words
    \return an estimate of the probability  of words \c w1 and \c w2 been generated independently
   */
  double comparePair( const Word& w1, const Word& w2, const vector<double>& measureDistr)const;

  //! Compare the pair of words using similarity measure.
  /*!
    \param w1 - the first word
    \param w2 - the second word.
    \return the distance (value of the similarity measure) 
   */
  double comparePair( const Word& w1, const Word& w2)const;
 private:
  WordPairComparison( const WordPairComparison& );
  WordPairComparison& operator = (const WordPairComparison&);
  const StringSimilarityMeasure* pSSM;
  int theRank;
};

////////////////////////////////////////////////////////////////////////////
//
// MEASURE ABSTRACT CLASS
//
///////////////////////////////////////////////////////////////////////////

//! Abstract interface for a class implementing string (Word)  similarity measure
class StringSimilarityMeasure
{
public:
  //! Returns a distance (measure) between two words.
  /*!
    \param w1 - the first word
    \param w2 - the second word
    \return the distance. It is assumed that greater values indicate less similarities.
   */
  virtual double measure(const Word& w1, const Word& w2)const  = 0;
};


////////////////////////////////////////////////////////////////////////////
//
// HAMMING DISTANCE
//
///////////////////////////////////////////////////////////////////////////


//! Implements the Hamming distance between two words
class HammingDistance : public StringSimilarityMeasure
{
public:
  HammingDistance() {}
  
  //! Returns the hamming distance scaled by the length of the longest word.
  /*!
    \param w1 - the first word
    \param w2 - the second word
    \return the scaled Hamming distance.
   */
  double measure(const Word& w1, const Word& w2) const;
};

//! Implements the Hamming distance between two cyclic words
class HammingDistanceCyclic : public StringSimilarityMeasure
{
public:
  HammingDistanceCyclic() {}
  //! Returns the hamming distance computed for all cyclic permutations of given words and scaled by the length of the longest word.
  /*!
    \param w1 - the first word
    \param w2 - the second word
    \return the scaled Hamming distance between cyclic words.
   */
  double measure(const Word& w1, const Word& w2) const;
};

//! Implements the Hamming distance between two cyclic words. 
/*! This is similar to \link HammingDistanceCyclic  HammingDistanceCyclic \endlink, 
except no penalty is given if words have different length. Let \f$ l_m = \min\{|w_1|,|w_2|\}\f$
be the length of shortest word, then distance is computed between initial segments of length \f$\_,\f$.
*/
class SubwordHammingDistanceCyclic : public StringSimilarityMeasure
{
 public:
  SubwordHammingDistanceCyclic() {}
  //! Returns the hamming distance computed for all cyclic permutations of initial segments of  given words.
  /*!
    \param w1 - the first word
    \param w2 - the second word
    \return the scaled Hamming distance between initial segments of cyclic words.
   */
  double measure(const Word& w1, const Word& w2) const;
};


////////////////////////////////////////////////////////////////////////////
//
// EDITING DISTANCE
//
///////////////////////////////////////////////////////////////////////////

//! Implements the Editing (Levenstein)  distance between two words
class EditingDistance : public StringSimilarityMeasure
{
public:
  EditingDistance() {}
  
  //! Returns the editing  distance.
  /*!
    \param w1 - the first word
    \param w2 - the second word
    \return the scaled Editing distance.
   */
  double measure(const Word& w1, const Word& w2) const;
};

//! Implements the Editing distance between two cyclic words. 
/*! No penalty is given if words have different length. Let \f$ l_m = \min\{|w_1|,|w_2|\}\f$
  be the length of shortest word, then distance is computed between initial segments of length \f$\_,\f$.
*/
class SubwordEditingDistanceCyclic : public StringSimilarityMeasure
{
 public:
  SubwordEditingDistanceCyclic() {}
  //! Returns the editing distance computed for all cyclic permutations of initial segments of  given words and scaled by the length.
  /*!
    \param w1 - the first word
    \param w2 - the second word
    \return the scaled Editing distance.
  */
 double measure(const Word& w1, const Word& w2) const;
};




////////////////////////////////////////////////////////////////////////////
//
// WHITEHEAD GRAPH
//
///////////////////////////////////////////////////////////////////////////
				     /*
class WhiteheadGraphSimilarity : public StringSimilarityMeasure
{
public:
  WhiteheadGraphSimilarity( const FreeGroup& F): theGroup( F ) {}
  double measure(const Word&, const Word&) const;
private:
  FreeGroup theGroup;
};
				     */
////////////////////////////////////////////////////////////////////////////
//
// MOTIVES
//
///////////////////////////////////////////////////////////////////////////
				 /*
class MotiveSimilarity : public StringSimilarityMeasure
{
public:
  MotiveSimilarity( const FreeGroup& F): theGroup( F ) {}
  double measure(const Word&, const Word&) const;
private:
  FreeGroup theGroup;
};
				 */


#endif
