// Copyright (C) 2005 Alexander Ushakov
// Contents: Definition of class BKLLeftNormalForm
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//

#ifndef _BKLLeftNormalForm_h_
#define _BKLLeftNormalForm_h_

#include "tuples.h"
#include "Permutation.h"

// #include "vector"
#include "list"
using namespace std ;

class BraidGroup;
class Word;


//---------------------------------------------------------------------------//
//--------------------------- BKLLeftNormalForm -----------------------------//
//---------------------------------------------------------------------------//


//! Birman-Ko-Lee left normal form. 
/*! 
  The standard presentation for BKL is a pair \f$(p,(\xi_1,\ldots,\xi_m))\f$, where \f$p\in \mathbb{Z}\f$ and
  each \f$\xi_i \in S_n\f$ is a permutation.
*/

class BKLLeftNormalForm
{
  
 public:

  //! Presentation of a normal form.
  /*!
    The first component specifies the rank of a braid group.
    The second component specifies  the power of the twist.
    The third component specifies the list of braid permutations.
  */
  typedef triple< int , int , list< Permutation > > NF;
  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Constructors:                                      //
  //                                                     //
  /////////////////////////////////////////////////////////


  //! Create trivial normal form  
  BKLLeftNormalForm( ) :
    theOmegaPower( 0 ) { }
    
    
  //! Copy constructor
  BKLLeftNormalForm( const BKLLeftNormalForm& bkl ) : 
    theRank( bkl.theRank ) ,
    theOmegaPower( bkl.theOmegaPower ) ,
    theDecomposition( bkl.theDecomposition ) { }
    
  //! Create normal form by its presentation (no check that given presentation is correct, use adjustDecomposition if the presentation is incorrect).
  BKLLeftNormalForm( int rank , int p , const list< Permutation >& d ) : 
    theRank( rank ) ,
    theOmegaPower( p ) ,
    theDecomposition( d ) { }
    
  //! Create normal form by its presentation (no check that given presentation is correct, use adjustDecomposition if the presentation is incorrect).
  BKLLeftNormalForm( const NF& bkl ) : 
    theRank( bkl.first ) ,
    theOmegaPower( bkl.second ) ,
    theDecomposition( bkl.third ) { }


  //! Create a normal form of a braid word.
  BKLLeftNormalForm( const BraidGroup& G , const Word& w );


  //! Assignment operator.
  BKLLeftNormalForm& operator = ( const BKLLeftNormalForm& bkl ) {
    theRank = bkl.theRank;
    theOmegaPower = bkl.theOmegaPower;
    theDecomposition = bkl.theDecomposition;
  }
  
  
  //! Assignment operator.
  BKLLeftNormalForm& operator = ( const NF& bkl ) {
    theRank = bkl.first;
    theOmegaPower = bkl.second;
    theDecomposition = bkl.third;
  }
  
  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Operators:                                         //
  //                                                     //
  /////////////////////////////////////////////////////////

 public:

  BKLLeftNormalForm& operator *= ( const BKLLeftNormalForm& bkl );
  BKLLeftNormalForm  operator *  ( const BKLLeftNormalForm& bkl ) const;
  BKLLeftNormalForm  operator -  ( ) const;

  bool operator ==( const BKLLeftNormalForm& bkl ) const;
  bool operator !=( const BKLLeftNormalForm& bkl ) const;
  bool operator < ( const BKLLeftNormalForm& bkl ) const;


  /////////////////////////////////////////////////////////
  //                                                     //
  //  Accessor functions:                                //
  //                                                     //
  /////////////////////////////////////////////////////////

 public:

  //! Get the power of omega.  
  int getPower( ) const { return theOmegaPower; }
  //! Get the decomposition.
  const list< Permutation >& getDecomposition( ) const { return theDecomposition; }

  //! Get the presentation of the normal form
  operator NF( ) const { return NF( theRank , theOmegaPower , theDecomposition ); }

  //! Check if the normal form os trivial.
  bool isTrivial( ) const { return theOmegaPower==0 && theDecomposition.size()==0; }

  //! Returns a cyclic permutation \f$\delta = (1,2,3,\ldots,n-1,0)\f$.
  static Permutation getTinyTwistPermutation( int theIndex ) {
    Permutation result( theIndex );
    for( int i=0 ; i<theIndex ; ++i )
      result[i] = (i+theIndex-1)%theIndex;
    return result;
  }


  //! Return a word represented by the normal form.
  Word getWord( ) const;
  
  
  //! Check if two normal forms are conjugate.
  pair< bool , Word > areConjugate( const BKLLeftNormalForm& bkl ) const;
  
  //! Cycle operation.
  pair< BKLLeftNormalForm , BKLLeftNormalForm > cycle( ) const;
  //! Decycle operation.
  pair< BKLLeftNormalForm , BKLLeftNormalForm > decycle( ) const;
  //! Find a representative of the supper summit set.
  pair< BKLLeftNormalForm , BKLLeftNormalForm > findSSSRepresentative( ) const;


  /////////////////////////////////////////////////////////
  //                                                     //
  //  Modifiers:                                         //
  //                                                     //
  /////////////////////////////////////////////////////////

 public:

  //! Set the power of \f$\Delta\f$. The result is always a correct normal form.
  inline void setPower( int p )
    { theOmegaPower = p; }

  //! Set the decomposition \f$(\xi_1,\ldots,\xi_m)\f$. The result might be an incorrect form (not satisfying some gready conditions of this BKL form).
  inline void setDecomposition( const list< Permutation >& d )
    { theDecomposition = d; }
  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Internal functions:                                //
  //                                                     //
  /////////////////////////////////////////////////////////

 private:
  
  
  //! Invert the normal form
  NF inverse( ) const;
  
  
  //! Multiply the normal form by another normal form on the right.
  NF multiply( const BKLLeftNormalForm& bkl ) const;


  //! The result of the transformation.
  enum transformationResult { TWO_MULTIPLIERS , ONE_MULTIPLIER , NO_CHANGE };
  //! Main function for computing the normal form of the word.
  static transformationResult transform ( int theIndex , Permutation& p1 , Permutation& p2 );


  //! Function returns a set of permutation braids (well, permutations) such that ... Must be checked
  static set<Permutation> getSimpleConjugators( const NF& bkl );


  //! Function returns a set of permutation braids (well, permutations) such that ... Must be checked
  static set<Permutation> getSimpleSummitConjugators( const NF& bkl );

  
  //! Function adjusting the decomposition in a normal form.
  static void adjustDecomposition( int rank ,
				   int& power ,
				   list<Permutation>& decomp );
  // this function transforms any (reasonable) decomposition
  // to a decomposition of a normal form 

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Data members:                                      //
  //                                                     //
  /////////////////////////////////////////////////////////

 private:

  //! The rank of the braid group (number of strands).
  int theRank;

  //! Power of omega.
  int theOmegaPower;

  //! Sequence of permutations.
  list< Permutation > theDecomposition;
  
};


ostream& operator << ( ostream& os, const BKLLeftNormalForm& bkl );

#endif
