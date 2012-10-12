// Copyright (C) 2005 Alexander Ushakov
// Contents: Definition of class KLProtocolInstance
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//

#ifndef _KLKeyGeneration_h_
#define _KLKeyGeneration_h_

#include "Word.h"

#include <vector>
using namespace std;


//---------------------------------------------------------------------------//
//--------------------------- KLProtocolInstance ----------------------------//
//---------------------------------------------------------------------------//


class KLProtocolInstance
{
  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Constructors:                                      //
  //                                                     //
  /////////////////////////////////////////////////////////
  
 public:

  //! Create an instance of the protocol.
  KLProtocolInstance( int braid_rank , Word base_word , Word keyA , Word KeyB );
  
  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Accessors:                                         //
  //                                                     //
  /////////////////////////////////////////////////////////

 public:
  
  //! Generate a random instance of the protocol
  /*!
    \param braid_rank - rank of the braid group;
    \param baseLenth - the length of the base word
    \param keyLength - the length of the key
   */
  static KLProtocolInstance random( int braid_rank , int baseLenth , int keyLength );

  //! (accessor function) get the rank of the braid group
  int  getBraidRank  ( ) const { return theRank       ; }
  //! (accessor function) get the private key of the first  party
  Word getPrivateKeyA( ) const { return thePrivateKeyA; }
  //! (accessor function) get the private key of the second party
  Word getPrivateKeyB( ) const { return thePrivateKeyB; }
  //! (accessor function) get the public  key of the first  party
  Word getPublicKeyA ( ) const { return thePublicKeyA ; }
  //! (accessor function) get the public  key of the second party
  Word getPublicKeyB ( ) const { return thePublicKeyB ; }
  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Data members:                                      //
  //                                                     //
  /////////////////////////////////////////////////////////

 private:

  //! the rank of the braid group
  int theRank;

  //! the base element (denoted by \f$w\f$)
  Word theBase;

  //! the private key of the first party (denoted by \f$A\f$)
  Word thePrivateKeyA;

  //! the public key of the first party \f$P_B = D(A^{-1} w A)\f$
  Word thePublicKeyA;

  //! the private key of the second party (denoted by \f$B\f$)
  Word thePrivateKeyB;

  //! the public key of the second party \f$P_B = D(B^{-1} w B)\f$
  Word thePublicKeyB;

};


#endif
