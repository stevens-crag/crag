// Copyright (C) 2005 Alexander Ushakov
// Contents: Implementation of class ShftConjKeyInstanceGarside
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//


#include "RanlibCPP.h"
#include "braid_group.h"
#include "ShftConjKeyGenerationGarside.h"


//---------------------------------------------------------------------------//
//------------------------ ShftConjKeyInstanceGarside -----------------------//
//---------------------------------------------------------------------------//


ShftConjKeyInstanceGarside::ShftConjKeyInstanceGarside( int braid_rank , ThRightNormalForm publicKeyA , ThRightNormalForm privateKey ) :
  theRank( braid_rank ),
  thePrivateKey( privateKey )
{
  thePublicKey.first = publicKeyA;
  thePublicKey.second = shiftedConjugation( publicKeyA , privateKey );
}


//---------------------------------------------------------------------------//
//--------------------------------- random ----------------------------------//
//---------------------------------------------------------------------------//


ShftConjKeyInstanceGarside ShftConjKeyInstanceGarside::random( int braid_rank , int baseLength , int keyLength )
{
  ThRightNormalForm publicKeyA = ThRightNormalForm::randomPositive( braid_rank , baseLength );
  ThRightNormalForm privateKey = ThRightNormalForm::randomPositive( braid_rank , keyLength );
  if( baseLength>0 )
    publicKeyA.setPower( -RandLib::ur.irand( 0 , (baseLength+1)/2 ) );
  if( keyLength>0 )
    privateKey.setPower( -RandLib::ur.irand( 0 , (keyLength+1)/2 ) );
  return ShftConjKeyInstanceGarside( braid_rank , publicKeyA , privateKey );
}


//---------------------------------------------------------------------------//
//------------------------------- Algorithms --------------------------------//
//---------------------------------------------------------------------------//


ThRightNormalForm shiftedConjugation( const ThRightNormalForm& w , const ThRightNormalForm& c )
{
  int r1 = w.getRank( );
  int r2 = c.getRank( );
  
  int rank = 1 + ( r1>r2 ? r1:r2 );
  
  ThRightNormalForm iw = w.increaseRank( rank );
  ThRightNormalForm ic = c.increaseRank( rank );
  ThRightNormalForm d = Permutation::getCyclePermutation( rank );
  ThRightNormalForm s1( rank , Word(1) );
  
  return ic * (-d* iw *d) * s1 * (-d* -ic *d);
}
