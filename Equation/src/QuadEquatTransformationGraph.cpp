// Copyright (C) 2005 Alexander Ushakov
// Contents: Implementation of class QuadEquationTranformationGraph
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//

#include "QuadEquatTransformationGraph.h"


//---------------------------------------------------------------------------//
//---------------------- QuadEquationTranformationGraph ---------------------//
//---------------------------------------------------------------------------//


QuadEquationTranformationGraph::QuadEquationTranformationGraph( const Equation& eq ) :
  theEquation( eq )
{
  // initialize the graph
  equationsInProcess[ theEquation.getTheEquation( ) ] = theGraph.newVertex( );
  hasSolution = isTrivialSolution( theEquation.getTheEquation( ) );
}


//---------------------------------------------------------------------------//
//----------------------------------- isDone --------------------------------//
//---------------------------------------------------------------------------//


bool QuadEquationTranformationGraph::isDone( ) const
{
  return equationsInProcess.size( )==0;
}


//---------------------------------------------------------------------------//
//--------------------------------- hasSolution -----------------------------//
//---------------------------------------------------------------------------//


const IntLabeledGraph& QuadEquationTranformationGraph::getGraph( ) const
{
  return theGraph;
}


//---------------------------------------------------------------------------//
//--------------------------------- hasSolution -----------------------------//
//---------------------------------------------------------------------------//


bool QuadEquationTranformationGraph::solutionFound( ) const
{
  return hasSolution;
}


//---------------------------------------------------------------------------//
//------------------------------ isTrivialSolution --------------------------//
//---------------------------------------------------------------------------//


bool QuadEquationTranformationGraph::isTrivialSolution( const Word& e ) const
{
  Word r;
  for( Word::const_iterator e_it=e.begin() ; e_it!=e.end( ) ; ++e_it ) {
    int g = *e_it;
    if( theEquation.isGenerator( g ) )
      r.push_back( g );
  }
  return r.length( )==0;
}


//---------------------------------------------------------------------------//
//----------------------------------- extend --------------------------------//
//---------------------------------------------------------------------------//


void QuadEquationTranformationGraph::extend( )
{
  typedef IntLabeledGraph::edge_type edge_type;
  typedef IntLabeledGraph::vertex_type vertex_type;

  // cout << "Extend" << endl;
  
  // A. choose any non-processed equation
  if( isDone( ) )
    return;
  
  
  Word eq     = (*equationsInProcess.begin( )).first;
  int  eq_num = (*equationsInProcess.begin( )).second;
  equationsInProcess.erase( equationsInProcess.begin( ) );
  processedEquations[eq] = eq_num;
  
  
  // B. Find the neighbours of the current equation and add them into the graph
  set< Word > n = getNeighbours( eq );

  for( set< Word >::const_iterator n_it=n.begin( ) ; n_it!=n.end( ) ; ++n_it ) {

    Word new_eq = *n_it;
    
    
    int target;
    map< Word , int >::iterator e_it;
    if( ( e_it = processedEquations.find(new_eq) ) !=processedEquations.end( ) ) {
      target = (*e_it).second;
    } else if( ( e_it = equationsInProcess.find(new_eq) ) !=equationsInProcess.end( ) ) {
      target = (*e_it).second;
    } else {
      target = equationsInProcess[ new_eq ] = theGraph.newVertex( );
      if( isTrivialSolution(new_eq) )
	hasSolution = true;
    }
    
    theGraph.newEdge( eq_num , edge_type( target , 1 ) );
    
  }
}


//---------------------------------------------------------------------------//
//-------------------------- applyAdjointTransformation ---------------------//
//---------------------------------------------------------------------------//


Word QuadEquationTranformationGraph::applyAdjointTransformation( const Word& eq , int x , const Word& im )
{
  Word image;
  Word::const_iterator eq_it = eq.begin( );
  for( ; eq_it!=eq.end( ) ; ++eq_it ) {
    int g = *eq_it;
    
    if( g==x ) {
      image.push_back( im );
    } else if( g==-x) {
      image.push_back( -im );
    } else {
      image.push_back( g );
    }
  }

  // image = image.minimalEquivalentForm( permutableGenerators , true , true );
  return image;
}


//---------------------------------------------------------------------------//
//-------------------------------- getNeighbours ----------------------------//
//---------------------------------------------------------------------------//


set< Word > QuadEquationTranformationGraph::getNeighbours( const Word& eq )
{
  set< Word > result;
  
  // construct all of its descendants
  Word::const_iterator eq_it = eq.begin( );
  for( ; eq_it!=eq.end( ) ; ++eq_it ) {
    
    // skip generators
    if( theEquation.isGenerator( *eq_it ) )
      continue;
    
    // look to the left
    Word::const_iterator eq_it2 = eq_it;
    if( eq_it2!=eq.begin( ) ) {
      int lg = *(--eq_it2);
      // cout << *eq_it << " , " << lg << endl;
      if( *eq_it-lg!=0 )
	result.insert( applyAdjointTransformation( eq , *eq_it , Word( -lg ) * Word( *eq_it ) ) );
    }
    
    // look to the right
    eq_it2 = eq_it;
    if( ++eq_it2!=eq.end( ) ) {
      int lg = *eq_it2;
      // cout << *eq_it << " , " << lg << endl;
      if( *eq_it-lg!=0 )
	result.insert( applyAdjointTransformation( eq , *eq_it , Word( *eq_it ) * Word( -lg ) ) );
    }
  }
  
  return result;
}
