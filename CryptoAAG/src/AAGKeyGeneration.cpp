// Copyright (C) 2005 Alexander Ushakov
// Contents: Implementation of class AAGProtocolInstance
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//

#include "RanlibCPP.h"
#include "AAGKeyGeneration.h"
#include "ShortBraidForm.h"
#include "braid_group.h"
// #include "FreeGroup.h"
#include "tuples.h"
#include <set>
#include "AAGChallengeGeneration.h"
#include "MajorDump.h"

//---------------------------------------------------------------------------//
//------------------------------- randomKey ---------------------------------//
//---------------------------------------------------------------------------//


Word constructKey( int N , const vector< Word >& sbgp , const Word& decomposition, bool useForm = true )
{
  Word result;
  int len = decomposition.length( );
  auto d_it = decomposition.begin( );
  for( ; d_it!=decomposition.end( ) ; ++d_it ) {
    int ngen = *d_it;
    if( ngen<0 )
      result *= sbgp[-ngen-1].inverse( );
    else
      result *= sbgp[ngen-1];
  }

  if ( useForm )
    return shortenBraid( N , result.freelyReduce( ) );
    // return shortBraidForm( N , result.freelyReduce( ) );
  else
    return result.freelyReduce( );
}


//---------------------------------------------------------------------------//
//---------------------------- generateAAGInstance --------------------------//
//---------------------------------------------------------------------------//


vector< Word > conjugateSubgroup_PGBF( int N , const vector< Word >& sbgp , Word w , bool useForm )
{
  vector< Word > result( sbgp.size( ) );
  for( int i=0 ; i<sbgp.size( ) ; ++i )
    if ( useForm )
      result[i] = shortenBraid( N , -w * sbgp[i] * w );
  //result[i] = shortBraidForm( N , w.inverse( )*sbgp[i]*w );
    else
      result[i] = -w * sbgp[i] * w ;

  return result;
}


//---------------------------------------------------------------------------//
//----------------------------- generateSubgroup ----------------------------//
//---------------------------------------------------------------------------//


vector< Word > generateSubgroup( int N , int num_gens , int min_len , int max_len )
{
  // choose random lengths of a generating set
  int total_length = 0;
  //  cout << "@am WARNING !!! length hack" << endl;
  vector< int > lenghts( num_gens , 0 );
  for( int i=0 ; i<num_gens ; ++i ) {
    // lenghts[i] = min_len + ::rand( )%(max_len-min_len+1);
    lenghts[i] = min_len + RandLib::ur.irand( 0 , max_len-min_len );
    total_length += lenghts[i];
    // cout << "  len" << i << " -> " << lenghts[i] << endl;
  }

  // generate generating set until all generators are involved
  bool progress = true;
  vector< Word > result( num_gens );
  while( progress ) {

    set< triple< int , int , int > > gens;
    for( int i=0 ; i<N-1 ; ++i )
      gens.insert( triple< int , int , int >( RandLib::ur.irand(0,3*N-1) , i , i ) );
    for( int i=N-1 ; i<total_length ; ++i )
      gens.insert( triple< int , int , int >( RandLib::ur.irand(0,3*N-1) , rand( )%(N-1) , i ) );

    result = vector< Word >( num_gens );
    set< triple< int , int , int > >::iterator g_it = gens.begin( );
    for( int i=0 ; i<num_gens ; ++i ) {
      for( int k=0 ; k<lenghts[i] ; ++k , ++g_it )
        result[i] *= Word( RandLib::ur.irand(0,1)==0 ? 1+(*g_it).second : -1-(*g_it).second );
    }
    result = shortBraidSbgpForm( N , result );

    // check if all generators are involved
    vector< bool > involved( N , false );
    for( int i=0 ; i<num_gens ; ++i ) {
      Word& w = result[i];
      auto w_it = w.begin( );
      for( ; w_it!=w.end( ) ; ++w_it )
        involved[abs( *w_it )] = true;
    }

    // if some generators are not involved then start all over again
    progress = false;
    for( int i=1 ; i<N ; ++i )
      if( !involved[i] )
        progress = true;
    //if( progress )
    //  cout << "Weak generating set" << endl;
  }

  return result;
}


//---------------------------------------------------------------------------//
//---------------------------- generateAAGInstance --------------------------//
//---------------------------------------------------------------------------//


AAGProtocolInstance AAGProtocolInstance::random( int N , int num_gens , int min_len , int max_len ,
						 int AliceDecompositionLength , int BobDecompositionLength )
{
  return AAGProtocolInstance( N , 
			      generateSubgroup( N , num_gens , min_len , max_len ),
			      generateSubgroup( N , num_gens , min_len , max_len ),
			      Word::randomWord( num_gens , AliceDecompositionLength ),
			      //generateHardProductOfGenerators( num_gens , AliceDecompositionLength ),
			      Word::randomWord( num_gens , BobDecompositionLength )
			      //generateHardProductOfGenerators( num_gens , BobDecompositionLength )
                            );
}

//---------------------------------------------------------------------------//
//---------------------------- generateAAGIChallenge ------------------------//
//---------------------------------------------------------------------------//


AAGProtocolInstance AAGProtocolInstance::challenge( int N ,  int sg_conj_len, int c, int k, int key_len )
{
  return AAGProtocolInstance( N , 
			      AAGChallenge::generateSubgroup(  c,k, sg_conj_len ),
			      AAGChallenge::generateSubgroup(  c,k, sg_conj_len ),
			      AAGChallenge::generateKeyDecomp(  key_len  ),
			      AAGChallenge::generateKeyDecomp(  key_len  )
			      //, false
                            );
}


//---------------------------------------------------------------------------//
//---------------------------- AAGProtocolInstance --------------------------//
//---------------------------------------------------------------------------//


AAGProtocolInstance::AAGProtocolInstance( int N ,
					  const vector< Word >& aSbgp,
					  const vector< Word >& bSbgp, 
					  const Word& aDecomp,
					  const Word& bDecomp,
					  bool useForm ) :
  AlicePublicSbgp( aSbgp ),
  BobPublicSbgp( bSbgp ),
  AliceKeyDecomposition( aDecomp ),
  BobKeyDecomposition( bDecomp )
{

  AliceKey = constructKey( N , AlicePublicSbgp , aDecomp , useForm );
  //  *Dump::dump_out << "A" << endl;
  BobKey   = constructKey( N , BobPublicSbgp   , bDecomp , useForm );
  //  *Dump::dump_out << "B" << endl;
  BobConjugatedSbgp = conjugateSubgroup_PGBF( N , BobPublicSbgp   , AliceKey , useForm );
  //   *Dump::dump_out << "C" << endl;
  AliceConjugatedSbgp   = conjugateSubgroup_PGBF( N , AlicePublicSbgp , BobKey , useForm   );
  //   *Dump::dump_out << "D" << endl;

  if (useForm )
    theSharedKey = shortenBraid( N , 
				 AliceKey * 
				 BobKey * 
				 AliceKey.inverse( ) * 
				 BobKey.inverse( ) );
  else
    theSharedKey = 		 AliceKey * 
				 BobKey * 
				 AliceKey.inverse( ) * 
                                 BobKey.inverse( );

}



//---------------------------------------------------------------------------//
//---------------------------- AAGProtocolInstance --------------------------//
//---------------------------------------------------------------------------//


void AAGProtocolInstance::printStats( ostream& out ) 
{
  out << "# ASG gens          : " << getAlicePublicSbgp().size() << endl
      << "# BSG gens          : " << getBobPublicSbgp().size() << endl;

  int ASgLength = 0;
  int BSgLength = 0;
  int AConSgLength = 0;
  int BConSgLength = 0;

  for (int i=0;i<AlicePublicSbgp.size();i++){
    ASgLength    += AlicePublicSbgp[i].length();
    AConSgLength += AliceConjugatedSbgp[i].length();
  }

  for (int i=0;i<BobPublicSbgp.size();i++){
    BSgLength    += BobPublicSbgp[i].length();
    BConSgLength += BobConjugatedSbgp[i].length();
  }

  out << "ASG / AConSG length : " <<      ASgLength << " / " <<  AConSgLength << endl
      << "BSG / BConSG length : " <<      BSgLength << " / " <<  BConSgLength << endl;
  
  
  out << "AKey length         : " << getAliceKey().length() << endl
      << "BKey length         : " << getBobKey().length() << endl;
  
  out << "Shared key length   : " << getSharedKey().length() << endl;
    
}





//---------------------------------------------------------------------------//
//---------------------------- AAGProtocolInstance --------------------------//
//---------------------------------------------------------------------------//


Word AAGProtocolInstance::generateHardProductOfGenerators( int num_gens , int product_length )
{
  if( product_length<4 )
    return Word::randomWord( num_gens , product_length );
  
  auto result = Word::randomWord( num_gens , 2 ).toList( );
  result.push_back( -(*  result.begin( ) ) );
  result.push_back( -(*++result.begin( ) ) );
  
  while( product_length - result.size( ) > 1 ) {
    
    int g = RandLib::ur.irand( 1 , num_gens );
    g = RandLib::ur.irand( 0 , 1 )==0 ? -g : g;
    
    if( RandLib::ur.irand( 0 , 1 )==0 ) {
      result.insert( --result.end( ) , g );
      result.push_back( -g );
    } else {
      result.insert( --(--result.end( ) ) , g );
      result.push_back( -g );
    }
  }
  
  Word product = Word(std::move(result));
  while( product.length( )<product_length ) {
    int g = RandLib::ur.irand( 1 , num_gens );
    g = RandLib::ur.irand( 0 , 1 )==0 ? -g : g;
    product.push_front( g );
  }
  
  return product;
}
