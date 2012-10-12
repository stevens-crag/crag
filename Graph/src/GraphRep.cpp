// Copyright (C) 2005 Alexander Ushakov
// Contents: Implementation of class GraphRep
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//

#include "GraphRep.h"


//---------------------------------------------------------------------------//
//-------------------------------- GraphRep ---------------------------------//
//---------------------------------------------------------------------------//

int GraphRep::newState( )
{
	theStates[maxState] = state_type( maxState );
	return maxState++;
}


void GraphRep::newEdge( int v1 , int v2 )
{
	state_type& vert1 = theStates[v1];
	vert1.out.insert( GraphEdge( v2 ) );

	state_type& vert2 = theStates[v2];
	vert2.in.insert( GraphEdge( v1 ) );
}


void GraphRep::clear( )
{
	theStates.clear( );
	maxState = 0;
}
