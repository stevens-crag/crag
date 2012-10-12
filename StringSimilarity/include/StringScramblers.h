// Contents: Classes implementing scrambling operators acting on group words (or strings)

//
// Principal Author: Alexei Miasnikov (2004)
//
// Status:
//
// Revision History:
//

#ifndef STRING_SCRAMBLER_H
#define STRING_SCRAMBLER_H


#include "Word.h"
#include <vector>


////////////////////////////////////////////////////////////////////////////
//
// SCRAMBLER ABSTRACT CLASS
//
///////////////////////////////////////////////////////////////////////////


class StringScrambler
{
public:
  virtual Word scramble( const Word& )const  = 0;
  virtual double fracChanged()const = 0;
};


////////////////////////////////////////////////////////////////////////////
//
// UNIFORM SCRAMBLER
//
///////////////////////////////////////////////////////////////////////////

class UniformScrambler : public StringScrambler
{
public:
  UniformScrambler( int gs, double prob ): numGens( gs ), theChangeProb( prob ),theFracChanged(0.0)  {}
  Word scramble( const Word& )const;
  double fracChanged()const { return theFracChanged; }

 private:
  double theChangeProb;
  double theFracChanged;
  int    numGens;
};

////////////////////////////////////////////////////////////////////////////
//
// SUBWORD SCRAMBLER
//
///////////////////////////////////////////////////////////////////////////

class SubwordScrambler : public StringScrambler
{
public:
  SubwordScrambler( int gs, double f ): numGens( gs ), theChangeFrac( f )  {}
  Word scramble( const Word& )const;
  double fracChanged()const { return theChangeFrac; }

 private:
  double theChangeFrac;
  int    numGens;
};

////////////////////////////////////////////////////////////////////////////
//
// MULTIPLY WORD SCRAMBLER
//
///////////////////////////////////////////////////////////////////////////

class WordMultiplyScrambler : public StringScrambler
{
public:
  WordMultiplyScrambler( int gs, double f ): numGens( gs ), theChangeFrac( f )  {}
  Word scramble( const Word& )const;
  double fracChanged()const { return theChangeFrac; }
  
 private:
  double theChangeFrac;
  int    numGens;
};


#endif
