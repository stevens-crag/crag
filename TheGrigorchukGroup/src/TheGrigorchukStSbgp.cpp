// Copyright (C) 2006 Alexander Ushakov
// Contents: Part of implementation of class TheGrigorchukGroupAlgorithms
// Algorithms for a normal subgroup St(1).
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//


#include "Word.h"
#include "TheGrigorchukGroupAlgorithms.h"


//---------------------------------------------------------------------------//
//----------------------- TheGrigorchukGroupAlgorithms ----------------------//
//---------------------------------------------------------------------------//


pair< Word , Word > TheGrigorchukGroupAlgorithms::liftToSTone( const Word& w1 , const Word& w2 )
{
  pair< Word , Word > result;
  Word r1 = reduce(w1);
  Word r2 = reduce(w2);
  
  // A. Clean w2
  for( Word::const_iterator w_it=r2.begin( ) ; w_it!=r2.end( ) ; ++w_it ) {

    int g1 = r1.length( )==0 ? 0 : *r1.begin( );
    int g2 = *w_it;

    switch( g2 ) {
    case 1:
      if( g1==3 ) {
	push_back( result.first , 1 );
	push_back( result.first , 2 );
	push_back( result.first , 1 );
	r1.pop_front( );
	break;
      }
      push_back( result.first , 1 );
      push_back( result.first , 3 );
      push_back( result.first , 1 );
      push_front( r1 , 4 );
      break;
    case 2:
      push_back( result.first , 4 );
      break;
    case 3:
      push_back( result.first , 2 );
      push_front( r1 , 1 );
      break;
    case 4:
      push_back( result.first , 3 );
      push_front( r1 , 1 );
      break;
    }
  }
  
  // B. find the decomposition w1 as an element of <b>
  pair< Word , list< Word > > D = TheGrigorchukGroupAlgorithms::decompositionBSbgp( r1 );
  for( list< Word >::const_iterator l_it=D.second.begin( ) ; l_it!=D.second.end() ; ++l_it ) {
    
    // current conjugator for b
    Word l = *l_it;
    Word c;
    for( Word::const_iterator w_it=l.begin() ; w_it!=l.end( ) ; ++w_it ) {
      
      switch( *w_it ) {
      case 1:
	push_back( c , 2 );
	break;
      case 4:
	push_back( c , 1 );
	push_back( c , 3 );
	push_back( c , 1 );
	break;
      }
    }
    result.first *= reduce( -c * Word(1) * Word(4) * Word(1) * c );
  }
  
  result.first = reduce( result.first );
  result.second = D.first;
  
  return result;
}
