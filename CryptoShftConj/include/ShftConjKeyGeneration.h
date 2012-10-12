// Copyright (C) 2005 Alexander Ushakov
// Contents: Definition of class ShftConjKeyInstance
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//

#ifndef _ShftConjKeyGeneration_h_
#define _ShftConjKeyGeneration_h_

#include "Word.h"

#include <vector>
using namespace std;


//---------------------------------------------------------------------------//
//-------------------------- ShftConjKeyInstance ----------------------------//
//---------------------------------------------------------------------------//


class ShftConjKeyInstance
{
  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Constructors:                                      //
  //                                                     //
  /////////////////////////////////////////////////////////
  
 public:

  //! Create an instance of the protocol.
  ShftConjKeyInstance( int braid_rank , Word publicKeyA , Word privateKey );
  
  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Accessors:                                         //
  //                                                     //
  /////////////////////////////////////////////////////////

 public:
  
  //! Generate a random instance of the protocol
  /*!
    \param braid_rank - rank of the braid group;
    \param baseLenth - the length of thePublicKey.first
    \param keyLength - the length of thePrivateKey
   */
  static ShftConjKeyInstance random( int braid_rank , int baseLenth , int keyLength );

  //! (accessor function) get the rank of the braid group
  int  getBraidRank  ( ) const { return theRank       ; }
  //! (accessor function) get the private key
  Word getPrivateKey( ) const { return thePrivateKey; }
  //! (accessor function) get the public key
  pair< Word , Word > getPublicKey ( ) const { return thePublicKey ; }
  
  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Data members:                                      //
  //                                                     //
  /////////////////////////////////////////////////////////

 private:

  //! the rank of the braid group
  int theRank;

  //! the public key of the first party \f$P_B = D(sh(A)^{-1} \sigma_1 sh(w) A)\f$
  pair< Word , Word > thePublicKey;

  //! the private key of the second party (denoted by \f$B\f$)
  Word thePrivateKey;
};


//---------------------------------------------------------------------------//
//------------------------------- Algorithms --------------------------------//
//---------------------------------------------------------------------------//


Word generatorShift( const Word& w );


#endif
