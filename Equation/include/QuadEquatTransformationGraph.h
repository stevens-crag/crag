// Copyright (C) 2005 Alexander Ushakov
// Contents: Definition of class QuadEquationTranformationGraph
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//


#ifndef _QuadEquationTranformationGraph_h_
#define _QuadEquationTranformationGraph_h_


#include "Equation.h"


//---------------------------------------------------------------------------//
//---------------------- QuadEquationTranformationGraph ---------------------//
//---------------------------------------------------------------------------//


class QuadEquationTranformationGraph
{
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Constructors:                                      //
  //                                                     //
  /////////////////////////////////////////////////////////
  
 public:
  
  QuadEquationTranformationGraph( const Equation& eq );
  
  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Accessors:                                         //
  //                                                     //
  /////////////////////////////////////////////////////////

 public:

  void extend( );
  

  //! Get the current graph
  const IntLabeledGraph& getGraph( ) const;

  
  //! Determine if the graph is completely constructed
  bool isDone( ) const;


  //! Determine if the equation has solutions
  bool solutionFound( ) const;
  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  I/O:                                               //
  //                                                     //
  /////////////////////////////////////////////////////////

 private:

  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Internal functions:                                //
  //                                                     //
  /////////////////////////////////////////////////////////

 private:
  
  //! Function returns all neighbours of the vertex in the graph
  set< Word > getNeighbours( const Word& eq );
  
  
  //! Function is used in getNeighbours( )
  Word applyAdjointTransformation( const Word& eq , int x , const Word& im );

  
  //! Check if an equation e has trivial solution
  bool isTrivialSolution( const Word& e ) const;
  
  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Data members:                                      //
  //                                                     //
  /////////////////////////////////////////////////////////

 private:
  
  //! Equation under consideration
  const Equation theEquation;
  
  
  //! The transformation graph
  /*!
    A graph with edges labelled by integers. 
    Labels of edges are all ones at the moment. In the future they will encode
    corresponding transformations.
   */
  IntLabeledGraph theGraph;

  
  //! Already processed equations and their numbers in the graph 
  map< Word , int > processedEquations;
  
  
  //! New equations and their numbers in the graph
  map< Word , int > equationsInProcess;


  //! Flag. True if the equation with the trivial solution is found
  bool hasSolution;
  
};


#endif
