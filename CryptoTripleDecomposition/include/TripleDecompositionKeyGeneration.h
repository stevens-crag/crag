// Copyright (C) 2005 Alexander Ushakov
// Contents: Definition of class TripleDecompositionProtocolInstance
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//

#ifndef _TripleDecompositionProtocolInstance_h_
#define _TripleDecompositionProtocolInstance_h_

#include "tuples.h"
#include "Word.h"
#include "ThRightNormalForm.h"


#include <vector>
using namespace std;


//---------------------------------------------------------------------------//
//-------------------- TripleDecompositionProtocolInstance ------------------//
//---------------------------------------------------------------------------//


//! Definition of the class TripleDecompositionProtocolInstance
/*!
  Objects of this class contain public and private information of two parties 
  as used in Kurt' key exchange protocol. Also it contains a routine for random
  generation of keys as proposed in 
  Y. Kurt, J. Koh, "A New Key Exchange Primitive Based on the Triple Decomposition Problem".  
*/

class TripleDecompositionProtocolInstance
{
  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Constructors:                                      //
  //                                                     //
  /////////////////////////////////////////////////////////
  
 public:

  //! Create an instance of the protocol for a particular choice of certain public/private information.
  TripleDecompositionProtocolInstance( int braid_rank ,
				       quadruple< Word , Word , Word , Word > conjugators ,
				       quintuple< Word , Word , Word , Word , Word > privateKeyA ,
				       quintuple< Word , Word , Word , Word , Word > privateKeyB );
  
  
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
  static TripleDecompositionProtocolInstance random( int braid_rank , int baseLenth , int keyLength );

  //! (accessor function) Get the rank of the braid group.
  int  getBraidRank  ( ) const { return theRank; }


  //! (accessor function) Get the private key of the first party (Alice).
  quintuple< Word , Word , Word , Word , Word > getPrivateKeyA( ) const { return thePrivateKeyA; }
  //! (accessor function) Get the private key of the second party (Bob).
  quintuple< Word , Word , Word , Word , Word > getPrivateKeyB( ) const { return thePrivateKeyB; }


  //! (accessor function) Get the public key of the first party (Alice).
  triple< Word , Word , Word > getPublicKeyA( ) const { return thePublicKeyA; }
  //! (accessor function) Get the public key of the second party (Bob).
  triple< Word , Word , Word > getPublicKeyB( ) const { return thePublicKeyB; }


  //! (accessor function) Get the shared key of two parties.
  ThRightNormalForm getSharedKey( ) const { return theSharedKey; }

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Internal functions:                                //
  //                                                     //
  /////////////////////////////////////////////////////////

 private:
  
  //! Generate a random reduced word over the alphabet \f$x_l,\ldots,x_u\f$ of length len.
  static Word randomWord( int lowerIndex , int upperIndex , int len );
  
  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Data members:                                      //
  //                                                     //
  /////////////////////////////////////////////////////////

 private:

  //! the rank of the braid group
  int theRank;
  
  
  //! The private key of the first party (Alice).
  /*!
    The private key \f$Private_A = (a_1,a_2,a_3,x_1,x_2)\f$. 
    See notation of Y. Kurt, J. Koh, "A New Key Exchange Primitive Based on the Triple Decomposition Problem".
  */
  quintuple< Word , Word , Word , Word , Word > thePrivateKeyA;


  //! The private key of the second party (Bob) 
  /*!
    The private key \f$Private_B = (b_1,b_2,b_3,y_1,y_2)\f$.
    See notation of Y. Kurt, J. Koh, "A New Key Exchange Primitive Based on the Triple Decomposition Problem".
  */
  quintuple< Word , Word , Word , Word , Word > thePrivateKeyB;
  
  
  //! The public key of the first party (Alice).
  /*!
    The public key \f$Public_A = (u,v,w)\f$.
    See notation of Y. Kurt, J. Koh, "A New Key Exchange Primitive Based on the Triple Decomposition Problem".
  */
  triple< Word , Word , Word > thePublicKeyA;

  //! The public key of the second party (Bob).
  /*!
    The public key \f$Public_B = (p,q,r)\f$.
    See notation of Y. Kurt, J. Koh, "A New Key Exchange Primitive Based on the Triple Decomposition Problem".
  */
  triple< Word , Word , Word > thePublicKeyB;
  
  
  //! The conjugators for the elementary commuting subgroups (public information).
  /*!
    The conjugators \f$s_1,s_2,s_3,s_4\f$.
    See notation of Y. Kurt, J. Koh, "A New Key Exchange Primitive Based on the Triple Decomposition Problem".
  */
  quadruple< Word , Word , Word , Word > theConjugators;
  

  //! The shared key \f$a_1 b_1 a_2 b_2 a_3 b_3\f$.
  ThRightNormalForm theSharedKey;
};


#endif
