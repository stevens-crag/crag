// Copyright (C) 2005 Alexander Ushakov
// Contents: Definition of class AAGProtocolInstance
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//

#ifndef _AAGKeyGeneration_h_
#define _AAGKeyGeneration_h_

#include "Word.h"

#include <vector>
using namespace std;


//---------------------------------------------------------------------------//
//-------------------------- AAGProtocolInstance ----------------------------//
//---------------------------------------------------------------------------//


class AAGProtocolInstance
{

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Constructors:                                      //
  //                                                     //
  /////////////////////////////////////////////////////////

  public:

  AAGProtocolInstance( int N ,
		       const vector< Word >& aSbgp,
		       const vector< Word >& bSbgp, 
		       const Word& aDecomp,
		       const Word& bDecomp,
		       bool useForm = true );
  // Constructor for AAG-Instance
  // parameters must be provided by some generating procedure


  /////////////////////////////////////////////////////////
  //                                                     //
  //  Accessors:                                         //
  //                                                     //
  /////////////////////////////////////////////////////////

  public:

  vector< Word > getAlicePublicSbgp( ) const { return AlicePublicSbgp; }
  vector< Word > getBobPublicSbgp  ( ) const { return BobPublicSbgp; }
  vector< Word > getAliceConjSbgp  ( ) const { return AliceConjugatedSbgp; }
  vector< Word > getBobConjSbgp    ( ) const { return BobConjugatedSbgp; }

  Word getAliceKey ( ) const { return AliceKey; }
  Word getBobKey   ( ) const { return BobKey; }
  Word getSharedKey( ) const { return theSharedKey; }


  void printStats( ostream& out );

  static AAGProtocolInstance random( int N , int num_gens , int min_len , int max_len ,
				     int AliceDecompositionLength , int BobDecompositionLength );
  
  static AAGProtocolInstance challenge( int N ,  int sg_conj_len, int c, int k, int key_len );


  /////////////////////////////////////////////////////////
  //                                                     //
  //  Internal functions:                                //
  //                                                     //
  /////////////////////////////////////////////////////////

  private:

  //! Function generates a hard key.
  static Word generateHardProductOfGenerators( int num_gens , int product_length );
  

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Data members:                                      //
  //                                                     //
  /////////////////////////////////////////////////////////

  private:

  vector< Word > AlicePublicSbgp;
  vector< Word >   BobPublicSbgp;

  Word AliceKeyDecomposition;
  Word   BobKeyDecomposition;

  Word AliceKey;
  Word   BobKey;
  
  vector< Word > AliceConjugatedSbgp;
  vector< Word >   BobConjugatedSbgp;

  Word theSharedKey;
};


//---------------------------------------------------------------------------//
//------------------------------- Algorithms --------------------------------//
//---------------------------------------------------------------------------//

vector< Word > conjugateSubgroup_PGBF( int N , const vector< Word >& sbgp , Word w , bool useForm = true );


#endif
