
#ifndef _PermutationEnumerator_H_
#define _PermutationEnumerator_H_

#include "VectorEnumerator.h"
#include "Permutation.h"


//---------------------------------------------------------------------------//
//------------------------ PermutationEnumerator ----------------------------//
//---------------------------------------------------------------------------//


class PermutationEnumerator : public VectorEnumerator
{

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Constructors                                       //
  //                                                     //
  /////////////////////////////////////////////////////////
  
 public:

  //! Constructor
  PermutationEnumerator( int L ) : theLength(L), theUsedElements( L , false ) { }


  //! No assignment operator
  PermutationEnumerator operator= ( const PermutationEnumerator& );
  
  
  //! No copy constructor
  PermutationEnumerator( const PermutationEnumerator& );

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Accessors                                          //
  //                                                     //
  /////////////////////////////////////////////////////////
  
 public:

  Permutation getPermutation( ) const { return Permutation( getSeq( ).begin( ) , getSeq( ).begin( )+theLength ); }

  virtual bool seqOK            ( ) const { 
    int last_value = getSeq( )[getLength( )-1];
    return !theUsedElements[last_value]  &&  getLength( )<=theLength;
  }

  virtual bool seqLimit         ( ) const {
    return getLength( ) > theLength;
  }
  virtual bool seqComplete      ( ) const {
    return theLength==getLength( ); 
  }

  virtual int  start (         ) const {
    return 0;
  }
  virtual int  next  ( int cur ) const {
    return cur+1;
  }
  virtual bool finish( int cur ) const {
    return cur>=theLength;
  }

  virtual void stepTo( ) { 
    int last_value = getSeq( )[getLength( )-1];
    theUsedElements[last_value] = true;
  }
  virtual void stepBack( int cur ) { 
    int last_value = getSeq( )[getLength( )];
    theUsedElements[last_value] = false;
  }


  /////////////////////////////////////////////////////////
  //                                                     //
  //  Data members                                       //
  //                                                     //
  /////////////////////////////////////////////////////////
  
 private:
  
  //! The length of a permutation
  int theLength;
  
  //! Array of flags - elements used in already constructed part of a permutation
  vector< bool > theUsedElements;
};

#endif
