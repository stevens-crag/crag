// Copyright (C) 2005 Alexander Ushakov
// Contents: Implementation of class Equation
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//

#include "Equation.h"


//---------------------------------------------------------------------------//
//--------------------------------- Equation --------------------------------//
//---------------------------------------------------------------------------//


Equation::Equation( int nGen , int nVar , const Word& eq ) :
  theNumberOfGenerators(nGen),
  theNumberOfVariables(nVar),
  theEquation(eq)
{
  
}


//---------------------------------------------------------------------------//
//-------------------------------- isQuadratic ------------------------------//
//---------------------------------------------------------------------------//


bool Equation::isQuadratic( ) const
{
  map< int , int > var_count;
  
  Word::const_iterator w_it = theEquation.begin( );
  for( ; w_it!=theEquation.end( ) ; ++w_it ) {
    int ag = abs( *w_it );
    if( ag>theNumberOfGenerators ) {
      var_count[ag]++;
      if( var_count[ag]==3 )
	return false;
    }
  }
  
  map< int , int >::const_iterator c_it = var_count.begin( );
  for( ; c_it!=var_count.end( ) ; ++c_it )
    if( (*c_it).second<2 )
      return false;

  return true;
}


//---------------------------------------------------------------------------//
//-------------------------------- isGenerator ------------------------------//
//---------------------------------------------------------------------------//


bool Equation::isGenerator( int g ) const
{
  return abs( g )<=theNumberOfGenerators;
}


//---------------------------------------------------------------------------//
//-------------------------------- isVariable -------------------------------//
//---------------------------------------------------------------------------//


bool Equation::isVariable( int g ) const
{
  return abs( g )>theNumberOfGenerators;
}


//---------------------------------------------------------------------------//
//----------------------------- trivialSolution -----------------------------//
//---------------------------------------------------------------------------//


bool Equation::trivialSolution( ) const
{
  Word image;
  Word::const_iterator w_it = theEquation.begin( );
  for( ; w_it!=theEquation.end( ) ; ++w_it )
    if( isGenerator(*w_it) )
      image.push_back( *w_it );
  
  return image.length()==0;
}


//---------------------------------------------------------------------------//
//------------------------------ operator << --------------------------------//
//---------------------------------------------------------------------------//


ostream& operator << ( ostream& os , const Equation& eq )
{
  os << "{" 
     << eq.getTheNumberOfGenerators( ) << " , " 
     << eq.getTheNumberOfVariables( ) << " , " 
     << eq.getTheEquation( ) << "}";
  
  return os;
}


//---------------------------------------------------------------------------//
//-------------------------- randomQuadraticEquation ------------------------//
//---------------------------------------------------------------------------//



Equation Equation::randomQuadraticEquation( int nGen , int nVar , int len )
{
  int l = len+2*nVar;
  vector< int > result( l , 0 );
  
  // assign positions to variables
  for( int v=0 ; v<nVar ; ++v ) {    // for all variables
    for( int p=0 ; p<2 ; ++p ) {     // for each variable twice

      // loop until the position and power for a variable is found
      bool posFound = false;
      while( !posFound ) {
	int pos = RandLib::ur.irand( 0 , l-1 );
	if( result[pos]!=0 )
	  continue;
	
	int l_pos = pos==0 ? l-1 : pos-1;
	int r_pos = pos==l-1 ? 0 : pos+1;
	
	int var = RandLib::ur.irand( 0 , 1 ) ? v+nGen+1 : -v-nGen-1;
	if( result[l_pos]+var==0 || result[r_pos]+var==0 )
	  continue;
	result[pos] = var;
	posFound = true;
      }
    }
  }
  
  
  for( int pos=0 ; pos<l ; ++pos ) {
    
    if( result[pos]!=0 )
      continue;
    
    int l_pos = pos==0 ? l-1 : pos-1;
    int r_pos = pos==l-1 ? 0 : pos+1;
    
    while( 1 ) {
      int gen = RandLib::ur.irand( 0 , nGen-1 );
      gen = RandLib::ur.irand( 0 , 1 ) ? gen+1 : -gen-1;
      if( result[l_pos]+gen==0 || result[r_pos]+gen==0 )
	continue;
      result[pos] = gen;
      break;
    }
    
  }
  
  return Equation( nGen , nVar , Word(result) );
}
