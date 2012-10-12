// Copyright (C) 2005 Alexander Ushakov
// Contents: Definition of class ShftConjKeyInstance. 
// This a version for Dehornoy's authentification protocol which works with Garside normal forms.
// Personally, I believe that it produces more interesting instances.
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//

#ifndef _ShftConjKeyGeneration_h_
#define _ShftConjKeyGeneration_h_

#include "Word.h"
#include "ThRightNormalForm.h"

#include <vector>
using namespace std;


//---------------------------------------------------------------------------//
//-------------------------- ShftConjKeyInstance ----------------------------//
//---------------------------------------------------------------------------//


//! Container for public and private information in Dehornoy's authentification protocol.
/*!
  Use this class to generate and keep an instance of Dehornoy's authentification protocol.
 */

class ShftConjKeyInstanceGarside
{
  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Constructors:                                      //
  //                                                     //
  /////////////////////////////////////////////////////////
  
 public:

  //! Create an instance of the protocol.
  /*!
    
  */
  ShftConjKeyInstanceGarside( int braid_rank , ThRightNormalForm publicKeyA , ThRightNormalForm privateKey );
  
  
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
  static ShftConjKeyInstanceGarside random( int braid_rank , int baseLenth , int keyLength );

  //! (accessor function) get the rank of the braid group
  int  getBraidRank  ( ) const { return theRank; }
  //! (accessor function) get the private key
  ThRightNormalForm getPrivateKey( ) const { return thePrivateKey; }
  //! (accessor function) get the public key
  pair< ThRightNormalForm , ThRightNormalForm > getPublicKey ( ) const { return thePublicKey; }

  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Data members:                                      //
  //                                                     //
  /////////////////////////////////////////////////////////

 private:

  //! the rank of the braid group
  int theRank;

  //! the private key of the first party (denoted by \f$A\f$)
  ThRightNormalForm thePrivateKey;

  //! the public key of the first party \f$P_B = D(sh(A)^{-1} \sigma_1 sh(w) A)\f$
  pair< ThRightNormalForm , ThRightNormalForm > thePublicKey;

};


//---------------------------------------------------------------------------//
//------------------------------- Algorithms --------------------------------//
//---------------------------------------------------------------------------//


ThRightNormalForm shiftedConjugation( const ThRightNormalForm& w , const ThRightNormalForm& c );


#endif
