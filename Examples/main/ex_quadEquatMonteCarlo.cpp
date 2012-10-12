// Copyright (C) 2005 Alexander Ushakov
// Contents: Example for Monte Carlo method for solving 
// the Diophantine problem for Quadratic Equations over Free Groups
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//


#include "Equation.h"
#include "QuadEquatTransformationGraph.h"


//---------------------------------------------------------------------------//
//----------------------------------- main ----------------------------------//
//---------------------------------------------------------------------------//


int main( )
{
  //& Quadratic equation : How do I generate a random quadratic equation?
  int nGen = 2;
  int nVar = 3;
  int len = 8;
  Equation eq = Equation::randomQuadraticEquation( nGen , nVar , len );
  cout << "Eq = " << eq << endl;
  
  
  //& Quadratic equation : How do I solve a quadratic equation using Monte-Carlo method?
  QuadEquationTranformationGraph MC( eq );
  int max_iter = 10000;
  for( int iter=0 ; iter<max_iter ; ++iter ) {
    MC.extend( );
    if( MC.solutionFound( ) ) {
      cout << "The trivial solution is found" << endl;
      cout << "Constructed part contains " << MC.getGraph( ).getVertices( ).size( ) << " vertices." << endl;
      break;
    }
    if( MC.isDone( ) ) {
      cout << "The graph is fully constructed" << endl;
      break;
    }
    if( iter==max_iter-1 )
      cout << "Solution is not found within " << max_iter << " iterations." << endl;
  }
  cout << "Final size = " << MC.getGraph( ).getVertices( ).size( ) << " vertices." << endl;
  
  return 0;
}
