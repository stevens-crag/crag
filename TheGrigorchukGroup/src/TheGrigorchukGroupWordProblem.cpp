// Copyright (C) 2006 Alexander Ushakov
// Contents: Part of implementation of class TheGrigorchukGroupAlgorithms
// Word Problem algorithm
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//


#include "Word.h"
#include "RanlibCPP.h"
#include "TheGrigorchukGroupAlgorithms.h"


//---------------------------------------------------------------------------//
//----------------------- TheGrigorchukGroupAlgorithms ----------------------//
//---------------------------------------------------------------------------//


void TheGrigorchukGroupAlgorithms::push_front( Word& w , int g )
{
  g = g>0 ? g : -g;
  if( w.length( )==0 ) {
    w.push_front( g );
    return;
  }

  int first_symbol = *w.begin( );
  
  // cases when a new symbol stack onto w
  if( g==1 ^ first_symbol==1 ) {
    w.push_front( g );
    return;
  }

  // cases when 2 symbols cancel out
  if( first_symbol==g ) {
    w.pop_front( );
    return;
  }
  
  // cases when 2 symbols merge
  w.pop_front( );
  w.push_front( 9-first_symbol-g );
}


//---------------------------------------------------------------------------//
//----------------------- TheGrigorchukGroupAlgorithms ----------------------//
//---------------------------------------------------------------------------//


void TheGrigorchukGroupAlgorithms::push_back( Word& w , int g )
{
  g = g>0 ? g : -g;
  if( w.length( )==0 ) {
    w.push_back( g );
    return;
  }

  int last_symbol = *--w.end( );

  // cases when a new symbol stack onto w
  if( last_symbol==1 ^ g==1 ) {
    w.push_back( g );
    return;
  }


  // cases when 2 symbols cancel out
  if( last_symbol==g ) {
    w.pop_back( );
    return;
  }
  
  // cases when 2 symbols merge
  w.pop_back( );
  w.push_back( 9-last_symbol-g );
}


//---------------------------------------------------------------------------//
//----------------------- TheGrigorchukGroupAlgorithms ----------------------//
//---------------------------------------------------------------------------//


pair< Word , Word > TheGrigorchukGroupAlgorithms::split( const Word& w )
{
  pair< Word , Word > result;
  Word& w0 = result.first;
  Word& w1 = result.second;
  for( Word::const_iterator w_it=w.begin( ) ; w_it!=w.end( ) ; ++w_it ) {
    
    int g0 = *w_it;
    if( g0==1 ) {
      int g1 = *++w_it;
      ++w_it;
      if( g1!=4 ) {
	push_back( w1 , 1 );
	push_back( w0 , g1+1 );
      } else {
	push_back( w0 , 2 );
      }
    } else {
      if( g0!=4 ) {
	push_back( w0 , 1 );
	push_back( w1 , g0+1 );
      } else {
	push_back( w1 , 2 );
      }
    }
  }
  
  return result;
}


//---------------------------------------------------------------------------//
//----------------------- TheGrigorchukGroupAlgorithms ----------------------//
//---------------------------------------------------------------------------//


bool TheGrigorchukGroupAlgorithms::trivial( const Word& w )
{
  Word r = reduce( w );
  return trivial_reduced( r );
}


//---------------------------------------------------------------------------//
//----------------------- TheGrigorchukGroupAlgorithms ----------------------//
//---------------------------------------------------------------------------//


bool TheGrigorchukGroupAlgorithms::trivial_reduced( const Word& w )
{
  if( w.length( )==0 )
    return true;

  // check that w belongs to ST(1)
  triple< int , int , int > a = abelianImage( w);
  if( a.first || a.second || a.third )
    return false;
  
  pair< Word , Word > sp = split( w );
  return trivial(sp.first) && trivial(sp.second);
}


//---------------------------------------------------------------------------//
//----------------------- TheGrigorchukGroupAlgorithms ----------------------//
//---------------------------------------------------------------------------//


Word TheGrigorchukGroupAlgorithms::reduce( const Word& w )
{
  Word result;
  for( Word::const_iterator w_it=w.begin( ) ; w_it!=w.end( ) ; ++w_it )
    push_back( result , *w_it );
  
  return result;
}


//---------------------------------------------------------------------------//
//----------------------- TheGrigorchukGroupAlgorithms ----------------------//
//---------------------------------------------------------------------------//


Word TheGrigorchukGroupAlgorithms::randomWord( int len )
{
  Word result;
  
  int a_symbol = RandLib::ur.irand( 0 , 1 );
  for( int i=0 ; i<len ; ++i ) {
    if( a_symbol==1 )
      result.push_back(1);
    else
      result.push_back( RandLib::ur.irand( 2 , 4 ) );
    a_symbol = 1-a_symbol;
  }
  
  return result;
}


//---------------------------------------------------------------------------//
//----------------------- TheGrigorchukGroupAlgorithms ----------------------//
//---------------------------------------------------------------------------//


int TheGrigorchukGroupAlgorithms::findOrder( const Word& w )
{
  if( w.length( )<=1 )
    return w.length( )+1;
  
  Word r = reduce( w );
  triple< int , int , int > A = abelianImage( w );
  if( A.first==1 )
    r = reduce( r*r );
  pair< Word , Word > S = split( r );
  int O1 = findOrder( S.first  );
  int O2 = findOrder( S.second );
  int O = O1>O2 ? O1 : O2;
  O = A.first==1 ? 2*O : O;
  
  return O;
}
