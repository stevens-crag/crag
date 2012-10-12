// Copyright (C) 2005 Alexander Ushakov
// Contents: Definition of class FPGroup
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//

#ifndef _FPGROUP_H_
#define _FPGROUP_H_


#include <iostream>
#include <vector>
#include <string>
#include <strstream>


#include "Word.h"
#include "Alphabet.h"

using namespace std;


//---------------------------------------------------------------------------//
//--------------------------------- FPGroup ---------------------------------//
//---------------------------------------------------------------------------//


//! Class FPGroup (finitely presented group)

class FPGroup
{
  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Constructors                                       //
  //                                                     //
  /////////////////////////////////////////////////////////

 public:
  
  FPGroup( int num=0 );
  //  FPGroup( const vector< string >& genNames );
  //  FPGroup( const vector< string >& genNames , const vector< Word >& relators );
  FPGroup( int numOfGen , const vector< Word >& relators );
  FPGroup( const FiniteAlphabet& a );
  FPGroup( const FiniteAlphabet& a, const vector< Word >& relators );  
  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Accessors:                                         //
  //                                                     //
  /////////////////////////////////////////////////////////

 public:

  //! Get the number of generators
  inline int numberOfGenerators( ) const { return numOfGenerators; }
  //! Get the names of the generators
  inline const vector< string >& getGeneratorsNames( ) const { return theAlphabet.getLetters(); }
  //! Get the alphabet
  inline const FiniteAlphabet& getAlphabet( ) const { return theAlphabet; }  
  // Get the relators
  inline const vector< Word >& relators( ) const { return theRelators; }

  
  //! Generate random equivalent word
  /*! Function inserts into random positions in the given word \f$w\f$ 
    relators (and their cyclic permutations and inverses) of the group 
    conjugated by randomly chosen words. The lengths of conjugators
    is chosen using geometric distribution with parameter = conj_param.
    When the required length is reached the word is being reduced and output.
    So, the result often is shorter than the given parameter length.
   */
  Word randomEqWord_Baltimore( const Word& w , int length , float conj_param ) const;

  
  
  //! Generate random trivial word
  /*! Function starts with a pair of trivial words \f$w_1 = \varepsilon\f$, \f$w_2 = \varepsilon\f$.
    On each iteration randomly takes a relators \f$r\f$, takes a random cyclic permutation \f$r'\f$ 
    of it or its inverse, then randomly cuts in two pieces  \f$r' = r_1 \circ r_2\f$ (one can be trivial) 
    and multiplies \f$w_1\f$ on the right by \f$r_1\f$ and \f$w_2\f$ on the right by \f$r_2^{-1}\f$.
    When \f$|w_1|+|w_2|\f$ reaches the length output freely reduced \f$w_1 w_2^{-1}\f$.
   */
  Word randomIdentity_Stack( int length ) const;

  //! Generate random trivial word
  /*! Function starts with a trivial word \f$w = \varepsilon\f$. On each iteration
    it multiplies \f$w\f$ on the right by a relator (and their cyclic permutations 
    and inverses) of the group 
    conjugated by randomly chosen words. Lengths of conjugators
    are chosen using geometric distribution with parameter = conj_param.
    When the required length is reached the word is being reduced and output.
    So, the result often is shorter than the given parameter length.
   */
  Word randomIdentity_Classic( int length , float conj_param ) const;
  
  //! Generate random trivial word (using randomEqWord_Baltimore)
  Word randomIdentity_Baltimore( int length , float conj_param ) const;

  //! Triangulate the set of relations. 
  /*! 
    \return The result is a finite group presentation \f$\langle Y;S \rangle\f$
    with all relators of length at most 3. Notice that the names of the original 
    generators change to \f$x_1,\ldots,x_n\f$ to avoid possible name collisions. 
    New generator names are \f$x_{n+1},\ldots\f$.
  */
  FPGroup triangulatePresentation( ) const;

 
  /////////////////////////////////////////////////////////
  //                                                     //
  //  I/O                                                //
  //                                                     //
  /////////////////////////////////////////////////////////

 public:
  

  /////////////////////////////////////////////////////////
  //                                                     //
  //  SWIG functions:                                    //
  //                                                     //
  /////////////////////////////////////////////////////////

   /*
#ifdef SWIG
%extend {

  char* __str__( ) const {
    strstream ostr;
    ostr << *self << ends;
    return ostr.str( );
  }

  Word read( const char* str ) {
    Word result;
    self->readWord( str , result );
    return result;
  }

  char* write( const Word& word ) const {
    strstream ostr;
    self->writeWord( ostr , word ) << ends;
    return ostr.str( );
  }

  char* write( const Mapping& m ) const {
    strstream ostr;
    self->writeMapping( ostr , m ) << ends;
    return ostr.str( );
  }

}
#endif
   */

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Internal functions                                 //
  //                                                     //
  /////////////////////////////////////////////////////////
  
 public:
  
  static vector< string > initializeGenNames( int num );
  
 private:
  
  friend ostream& operator << ( ostream& os , const FPGroup& group );
  friend istream& operator >> ( istream& is , FPGroup& group );
  
  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Data members                                       //
  //                                                     //
  /////////////////////////////////////////////////////////
  
 protected:
  
  int numOfGenerators;
  // does it make sense to have this variable? You can change it by generatorsNames.size( )
  //  vector< string > generatorsNames;
  vector< Word > theRelators;
  
  FiniteAlphabet theAlphabet;
  bool useDefaultAlphabet;  
};


#endif

