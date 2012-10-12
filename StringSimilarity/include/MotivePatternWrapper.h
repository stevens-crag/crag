/*
 *   $Id: MotivePatternWrapper.h,v 1.1 2005/11/28 23:25:56 amiasnik Exp $
 */
 
// Contents: Classes implementing a wrapper for IBM implementation of the Motive Pattern Extractor

//
// Principal Author: Alexei Miasnikov (2004)
//
// Status:
//
// Revision History:
//

#ifndef MOTIVE_PATTERN_WRAPPER_H
#define MOTIVE_PATTERN_WRAPPER_H


#include "global.h"
#include "FreeGroup.h"
#include "Word.h"
#include <vector>
#include <stdlib.h>

////////////////////////////////////////////////////////////////////////////  
//
// MOTIVE PATTERN 
//
///////////////////////////////////////////////////////////////////////////                                                     

struct MotivePattern
{
  int occurrences;
  int sequences;

  Word thePattern;
};

////////////////////////////////////////////////////////////////////////////
//
// MOTIVE PATTERN WRAPPER
//
///////////////////////////////////////////////////////////////////////////

class MotivePatternWrapper
{
public:
  MotivePatternWrapper( const FreeGroup& f): theGroup( f ) {}
  
  void initializeSequence( const vector<Word>&  );
  bool compute();
  vector<MotivePattern> getPatterns();
 private:
  string input_file;
  string output_file;


  FreeGroup theGroup;
};


#endif
