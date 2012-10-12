// Copyright (C) 2005 Alexander Ushakov
// Contents: Implementation of class ShftConjKeyGeneration
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//


#include "DehornoyForm.h"
#include "ShftConjKeyGeneration.h"


//---------------------------------------------------------------------------//
//--------------------------- ShftConjKeyInstance ---------------------------//
//---------------------------------------------------------------------------//


ShftConjKeyInstance::ShftConjKeyInstance( int braid_rank , Word publicKeyA , Word privateKey ) :
  theRank( braid_rank ),
  thePrivateKey( privateKey )
{
  DehornoyForm DF1( braid_rank+1 , thePrivateKey * generatorShift(publicKeyA) * Word(1) * generatorShift(-thePrivateKey) );
  thePublicKey.first = publicKeyA;
  thePublicKey.second = DF1.getDehornoyForm( );
}


ShftConjKeyInstance ShftConjKeyInstance::random( int braid_rank , int baseLenth , int keyLength )
{
  return ShftConjKeyInstance( braid_rank , 
			      Word::randomWord( braid_rank-1 , baseLenth ) ,
			      Word::randomWord( braid_rank-1 , keyLength ) );
}


Word generatorShift( const Word& w )
{
  Word result;
  Word::const_iterator w_it = w.begin( );
  for( ; w_it!=w.end( ) ; ++w_it ) {
    int g = *w_it;
    result.push_back( g>0 ? g+1 : g-1 );
  }
  return result;
}
