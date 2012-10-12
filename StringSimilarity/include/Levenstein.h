// Copyright (C) 2003 Dmitry Bormotov

// Contents: Definition of class Levenstein
//
// Principal Authors: Dmitry Bormotov
//
// Status: in progress
//
// Description: 
//
//   A class for computing the Levenstein distance between two
//   words. It converts words to strings (char*) and then calls a
//   known string algorithm.
//
// Revision History:
//

#ifndef Levenstein_h_
#define Levenstein_h_

//#include "CFGWord.h"
#include "Word.h"


//!   A class for computing the Levenstein distance between two words. 
/*!
  Converts words to strings (char*) and then calls a
  known string algorithm.
*/
class Levenstein {
  
public:

  ///////////////////////////////////////////////////////
  //                                                   //
  //  Constructors                                     //
  //                                                   //
  ///////////////////////////////////////////////////////

  Levenstein( ) { }


  ///////////////////////////////////////////////////////
  //                                                   //
  //  Public functions                                 //
  //                                                   //
  ///////////////////////////////////////////////////////

  //! Computes the Levenstein distance between two words
  /*!
    \param w1 - first word.
    \param w2 - second word
    \return non-scaled distance, i.e. the actual minimal  number of operations
    required to rewrite one word into another. 
   */
  int compute( const Word& w1, const Word& w2);

 
private:

  ///////////////////////////////////////////////////////
  //                                                   //
  //  Private functions                                //
  //                                                   //
  ///////////////////////////////////////////////////////

  char* wordToString( const Word& );


  ///////////////////////////////////////////////////////
  //                                                   //
  //  Data Members                                     //
  //                                                   //
  ///////////////////////////////////////////////////////

};

#endif
