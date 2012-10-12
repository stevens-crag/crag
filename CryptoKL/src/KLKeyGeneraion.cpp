// Copyright (C) 2005 Alexander Ushakov
// Contents: Implementation of class KLProtocolInstance
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//


#include "DehornoyForm.h"
#include "KLKeyGeneration.h"


//---------------------------------------------------------------------------//
//--------------------------- KLProtocolInstance ----------------------------//
//---------------------------------------------------------------------------//


KLProtocolInstance::KLProtocolInstance( int braid_rank , Word base_word , Word keyA , Word KeyB ) :
  theRank( braid_rank ),
  theBase( base_word ),
  thePrivateKeyA( keyA ),
  thePrivateKeyB( KeyB )
{
  DehornoyForm DF1( braid_rank-1 , -thePrivateKeyA * theBase * thePrivateKeyA );
  thePublicKeyA = DF1.getDehornoyForm( );
  DehornoyForm DF2( braid_rank-1 , -thePrivateKeyB * theBase * thePrivateKeyB );
  thePublicKeyB = DF2.getDehornoyForm( );
}


KLProtocolInstance KLProtocolInstance::random( int braid_rank , int baseLenth , int keyLength )
{
  return KLProtocolInstance( braid_rank , 
			     Word::randomWord( braid_rank , baseLenth ) ,
			     Word::randomWord( braid_rank , keyLength ) ,
			     Word::randomWord( braid_rank , keyLength ) );
}
