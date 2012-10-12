// Copyright (C) 2005 Alexander Ushakov
// Contents: Definition of class BKLRightNormalForm
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//

#ifndef _BKLRightNormalForm_h_
#define _BKLRightNormalForm_h_

#include "tuples.h"
#include "Permutation.h"

#include "vector"
#include "list"
using namespace std ;

class BraidGroup;
class Word;


#define ERROR_EXISTS

//---------------------------------------------------------------------------//
//--------------------------- BKLRightNormalForm ----------------------------//
//---------------------------------------------------------------------------//


//! (Contains errors!!! Not to be used yet) Birman-Ko-Lee right normal form. 
/*! 
  The standard presentation for BKL is a pair \f$(p,(\xi_1,\ldots,\xi_m))\f$, where \f$p\in \mathbb{Z}\f$ and
  each \f$\xi_i \in S_n\f$ is a permutation.
*/


class BKLRightNormalForm
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

 public:
  
  //! Create trivial normal form  
  BKLRightNormalForm( ) : theRank( 0 ) { 
#ifdef ERROR_EXISTS
    exit(1);
#endif
  }
    
    
  //! Copy constructor
  BKLRightNormalForm( const BKLRightNormalForm& bkl ) : 
    theOmegaPower( bkl.theOmegaPower ) ,
    theDecomposition( bkl.theDecomposition ) { 
#ifdef ERROR_EXISTS
    exit(1);
#endif
}

    
  //! Create normal form by its presentation (no check that given presentation is correct, use adjustDecomposition if the presentation is incorrect).
  BKLRightNormalForm( int rank , int power , const list< Permutation >& decomp ) : 
    theRank( rank ) ,
    theOmegaPower( power ) ,
    theDecomposition( decomp ) { 
#ifdef ERROR_EXISTS
    exit(1);
#endif
}
  
    
  //! Create a normal form of a braid word.
  BKLRightNormalForm( const BraidGroup& G , const Word& w );

  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Operatos:                                          //
  //                                                     //
  /////////////////////////////////////////////////////////

 public:

  //! Assignment operator.
  BKLRightNormalForm& operator = ( const BKLRightNormalForm& bkl ) {
    theRank = bkl.theRank;
    theOmegaPower = bkl.theOmegaPower;
    theDecomposition = bkl.theDecomposition;
  }

  //! Assignment operator.
  BKLRightNormalForm& operator = ( const NF& bkl ) {
    theRank = bkl.first;
    theOmegaPower = bkl.second;
    theDecomposition = bkl.third;
  }

  //! Invert a normal form.
  BKLRightNormalForm operator- ( ) const { return inverse( ); }
  
  
  //! Multiply two normal forms.
  BKLRightNormalForm operator* ( const BKLRightNormalForm& bkl ) const { return multiply( bkl ); }


  bool operator ==( const BKLRightNormalForm& bkl ) const;
  bool operator !=( const BKLRightNormalForm& bkl ) const;
  bool operator < ( const BKLRightNormalForm& bkl ) const;


  /////////////////////////////////////////////////////////
  //                                                     //
  //  Accessors:                                         //
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

  
  //! Function adjusting the decomposition in a normal form.
  static void adjustDecomposition( int rank ,
				   int& power ,
				   list<Permutation>& decomp );
  // this function transforms any (reasonable) decomposition
  // to a decomposition of a normal form 


  /////////////////////////////////////////////////////////
  //                                                     //
  //  Modifiers:                                         //
  //                                                     //
  /////////////////////////////////////////////////////////

 public:
  
  //! Set the power of \f$\Delta\f$. The result is always a correct normal form.
  inline void setPower( int p ) { theOmegaPower = p; }

  //! Set the decomposition \f$(\xi_1,\ldots,\xi_m)\f$. The result might be an incorrect form (not satisfying some gready conditions of this BKL form).
  inline void setDecomposition( const list< Permutation >& d ) { theDecomposition = d; }
  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Internal functions:                                //
  //                                                     //
  /////////////////////////////////////////////////////////

 private:

  //! The result of the transformation.
  enum   transformationResult { TWO_MULTIPLIERS , ONE_MULTIPLIER , NO_CHANGE };
  //! Main function for computing the normal form of the word.
  static transformationResult transform ( int theIndex , Permutation& p1 , Permutation& p2 );

  //! Invert the normal form
  BKLRightNormalForm inverse( ) const;
  
  
  //! Multiply the normal form by another normal form on the right.
  BKLRightNormalForm multiply( const BKLRightNormalForm& bkl ) const;


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


ostream& operator << ( ostream& os, const BKLRightNormalForm& rep );


#endif
