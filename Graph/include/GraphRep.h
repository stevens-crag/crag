// Copyright (C) 2005 Alexander Ushakov
// Contents: Definition of class GraphRep
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//

#ifndef _GraphRep_H_
#define _GraphRep_H_

#include "RefCounter.h"

#include "set"
#include "map"
using namespace std;




//---------------------------------------------------------------------------//
//-------------------------------- FSAEdge ----------------------------------//
//---------------------------------------------------------------------------//


struct GraphEdge 
{

  GraphEdge( ) : target(-1) { }
  // dummy constructor, not to be used, required for STL::map

  GraphEdge( int t ) : target(t) { }

  bool operator < ( const GraphEdge& e ) const { return target< e.target; }
  bool operator ==( const GraphEdge& e ) const { return target==e.target; }

  int target;
};


//---------------------------------------------------------------------------//
//------------------------------- FSAState ----------------------------------//
//---------------------------------------------------------------------------//


struct GraphState 
{
  typedef GraphEdge edge_type;
  
  GraphState( int s ) : theState(s) { }
  GraphState( ) : theState(0) { }
  // doomy constructor for class map
  
  bool operator== ( const GraphState& s ) const { 
    return theState==s.theState; 
  }
  bool operator< ( const GraphState& s ) const { 
    return theState<s.theState; 
  }

  set< edge_type > in;
  set< edge_type > out;
  int theState;
};


//---------------------------------------------------------------------------//
//-------------------------------- GraphRep ---------------------------------//
//---------------------------------------------------------------------------//


class GraphRep : public RefCounter
{

public:
	typedef GraphEdge  edge_type;
	typedef GraphState state_type;

	friend class Graph;


  /////////////////////////////////////////////////////////
  //                                                     //
  //  Constructors                                       //
  //                                                     //
  /////////////////////////////////////////////////////////
private:
	
	GraphRep( ) : maxState( 0 ) { }

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Operators                                          //
  //                                                     //
  /////////////////////////////////////////////////////////
private:

	

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Accessors                                          //
  //                                                     //
  /////////////////////////////////////////////////////////

public:
  GraphRep* clone( ) const { return new GraphRep( *this ); }

private:

	int  newState( );
	void newEdge( int v1 , int v2 );
	void clear( );

	const map< int, state_type >& getStates( ) const { return theStates; }

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Internal functions                                 //
  //                                                     //
  /////////////////////////////////////////////////////////
	private:

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Data members                                       //
  //                                                     //
  /////////////////////////////////////////////////////////
  private:

	map< int, state_type > theStates;
	// the set of vertices in the graph

	int maxState;
	// the maximal unused number of a vertex (used in addNewVertex)
};



#endif
