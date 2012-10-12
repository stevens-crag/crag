// Copyright (C) 2005 Alexander Ushakov
// Contents: Definition of class Graph
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//

#ifndef _Graph_H_
#define _Graph_H_


#include "ObjectOf.h"
#include "GraphRep.h"
#include "vector"
using namespace std;



//---------------------------------------------------------------------------//
//---------------------------------- Graph ----------------------------------//
//---------------------------------------------------------------------------//


class Graph : public ObjectOf< GraphRep >
{

public:
  typedef GraphRep::state_type state_type;
  typedef GraphRep::edge_type edge_type;
  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Constructors                                       //
  //                                                     //
  /////////////////////////////////////////////////////////
public:

  Graph( ) : ObjectOf< GraphRep >( new GraphRep() ) { }

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Accessors                                          //
  //                                                     //
  /////////////////////////////////////////////////////////
public:

  const map< int, state_type >& getStates( ) const { return look( )->getStates( ); }
  int  newState( ) { return change( ) -> newState( ); }
  void newEdge( int v1 , int v2 ) { change( ) -> newEdge( v1 , v2 ); }
  void clear( ) { change( ) -> clear( ); }


  /////////////////////////////////////////////////////////
  //                                                     //
  //  I/O                                                //
  //                                                     //
  /////////////////////////////////////////////////////////
private:
  
  // ostream& operator << ( ostream& os , const Graph& g );
  // Function is defined outside of the scope of GraphRep

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Operators                                          //
  //                                                     //
  /////////////////////////////////////////////////////////

};


ostream& operator << ( ostream& os , const Graph& g );

//! Function generates a random (non-directed) graph on N vertices. Each edge has equal probability edge_param to appear in the result.
Graph randomGraph( int N , float edge_param );

//! Compute a table of all lengths in the directed graph G. 
vector< vector< int > > lengthTable( const Graph& G );

//! Compute a table of inner (Gromov's) products in the directed graph G. 
vector< vector< int > > innerProductTable( const Graph& G , int origin );

//! For a finite directed graph G compute a constant of hyperbolisity.
float getHyperbolicityConst( const Graph& G );

#endif

