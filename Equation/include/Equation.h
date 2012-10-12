// Copyright (C) 2005 Alexander Ushakov
// Contents: Definition of class Equation
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//


#ifndef _Equation_h__
#define _Equation_h__


#include "GraphType.h"
#include "GraphConcept.h"
#include "GraphConceptAlgorithms.h"
using namespace Graphs;

#include "Word.h"


//---------------------------------------------------------------------------//
//--------------------------------- Equation --------------------------------//
//---------------------------------------------------------------------------//


class Equation
{

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Constructors:                                      //
  //                                                     //
  /////////////////////////////////////////////////////////
  
 public:

  Equation( int nGen , int nVar , const Word& eq );


  /////////////////////////////////////////////////////////
  //                                                     //
  //  Accessors:                                         //
  //                                                     //
  /////////////////////////////////////////////////////////

 public:

  //! Determine if g is a letter in this equation.
  bool isGenerator( int g ) const;


  //! Determine if g is a variable in this equation
  bool isVariable ( int g ) const;


  //! Get the number of the generators of the group
  int getTheNumberOfGenerators( ) const { return theNumberOfGenerators; }


  //! Get the number of variables in the formula
  int getTheNumberOfVariables( ) const { return theNumberOfVariables; }


  //! Get the word presentation of the equation
  const Word& getTheEquation( ) const { return theEquation; }

  
  //! Determine if the equation is quadratic
  bool isQuadratic( ) const;

  
  //! Determine if the equation has trivial solution
  bool trivialSolution( ) const;
  
  
  //! Generate random (strictly) quadratic equation of length \f$len + 2 \cdot nVars\f$ with \f$nGen\f$ generators, \f$nVar\f$ variables.
  /*!
    Routine  "arranges" len generators and \f$2 \cdot nVars\f$ variables (each variable twice) into a reduced equation.
    The distribution is not uniform among equations of this type (even though for large values of len
    I think it will be close to uniform). For uniform distribution one has
    to construct a FSA accepting all equations of this type, then assign weights to edges using dynamic
    programming, and finally choose words from that FSA accoring to the weights. 
   */
  static Equation randomQuadraticEquation( int nGen , int nVar , int len );


  /////////////////////////////////////////////////////////
  //                                                     //
  //  I/O:                                               //
  //                                                     //
  /////////////////////////////////////////////////////////

 private:
  
  friend ostream& operator << ( ostream& os , const Equation& eq );
  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Internal functions:                                //
  //                                                     //
  /////////////////////////////////////////////////////////

 private:
  
  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Data members:                                      //
  //                                                     //
  /////////////////////////////////////////////////////////

 private:

  //! the number of the generators of the group
  int theNumberOfGenerators;
  //! the number of variables in the formula (the actual number of variables in the equation can be lesser)
  int theNumberOfVariables;
  //! the presentation of the equation
  /*!
    A word theEquation is a sequence of generators. Each generator \f$g\f$ is interpreted the following way:
    - if \f$|g| \le theNumberOfGenerators\f$ then \f$g\f$ is the corresponding generator of the group
    - if \f$|g| > theNumberOfGenerators\f$ then \f$g\f$ is a variable with index \f$|g| - theNumberOfGenerators\f$
    raised in the power \f$\pm 1\f$ depending on the sign of \f$g\f$.
   */
  Word theEquation;

};


#endif

