// Copyright (C) 2005 Alexander Ushakov
// Contents: Implementation of class Graph
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//

#include "Graph.h"
#include "GraphAlgorithms.h"
#include "RanlibCPP.h"


//---------------------------------------------------------------------------//
//---------------------------------- Graph ----------------------------------//
//---------------------------------------------------------------------------//


ostream& operator << ( ostream& os , const Graph& g )
{
  typedef Graph::state_type state_type;
  typedef Graph::edge_type edge_type;
  const map< int, state_type >& theStates = g.getStates( );
  
  os << "{" << endl;
  map< int, state_type >::const_iterator v_it = theStates.begin( );
  for( ; v_it!=theStates.end( ) ; ++v_it ) {
    const set< edge_type >& out = (*v_it).second.out;
    os << "  " << (*v_it).second.theState << " -> ";
    for( set< edge_type >::const_iterator e_it=out.begin( ) ; e_it!=out.end( ) ; ++e_it ) {
      if( e_it!=out.begin( ) )
	os << ", ";
      os << (*e_it).target;
    }
    os << ";" << endl;
  }
  os << "}" << endl;
  
  return os;
}


//---------------------------------------------------------------------------//
//---------------------------------- Graph ----------------------------------//
//---------------------------------------------------------------------------//


Graph randomGraph( int N , float edge_param )
{
  Graph result;

  for( int i=0 ; i<N ; ++i )
    result.newState( );

  for( int i=0 ; i<N ; ++i ) {
    for( int j=i+1 ; j<N ; ++j ) {
      if( RandLib::ur.rand()<edge_param ) {
	result.newEdge( i , j );
	result.newEdge( j , i );
      }
    }
  }
  
  return result;
}


//---------------------------------------------------------------------------//
//---------------------------------- Graph ----------------------------------//
//---------------------------------------------------------------------------//


vector< vector< int > > lengthTable( const Graph& G )
{
  typedef Graph::state_type state_type;

  // Get the adjacency structure
  const map< int, state_type >& theStates = G.getStates( );

  // Get the number of states
  int N = theStates.size( );

  // Prepare enumeration of vertices
  map< int , int > stateNumber;
  map< int, state_type >::const_iterator s_it = theStates.begin( );
  for( int i=0 ; s_it!=theStates.end( ) ; ++i , ++s_it )
    stateNumber[(*s_it).first] = i;
  
  // Compute the result
  vector< vector< int > > result( N , vector< int >( ) );
  s_it = theStates.begin( );
  for( int i=0 ; s_it!=theStates.end( ) ; ++i , ++s_it ) {
    result[i] = vector< int >( N , -1 );
    map< int , int > dist = getDistances_out( G , (*s_it).first );
    map< int , int >::iterator d_it = dist.begin( );
    for( ; d_it!=dist.end( ) ; ++d_it )
      result[i][stateNumber[(*d_it).first]] = (*d_it).second;
  }
  
  return result;
}


//---------------------------------------------------------------------------//
//---------------------------------- Graph ----------------------------------//
//---------------------------------------------------------------------------//


vector< vector< int > > innerProductTable( const Graph& G , int origin )
{
  int N = G.getStates( ).size( );
  vector< vector< int > > dist = lengthTable( G );
  vector< vector< int > > result( N , vector< int >( ) );
  for( int i=0 ; i<N ; ++i ) {
    result[i] = vector< int >( N , -1 );
    for( int j=0 ; j<N ; ++j ) {
      if( dist[origin][i]==-1 || dist[origin][j]==-1 || dist[i][j]==-1 )
	continue;
      result[i][j] = dist[origin][i] + dist[origin][j] - dist[i][j];
    }
  }
  
  return result;
}


//---------------------------------------------------------------------------//
//---------------------------------- Graph ----------------------------------//
//---------------------------------------------------------------------------//


float getHyperbolicityConst( const Graph& G )
{
  // Get the adjacency structure
  typedef Graph::state_type state_type;
  const map< int, state_type >& theStates = G.getStates( );
  if( theStates.empty( ) )
    return 0;
  
  float result = 0;
  int N = theStates.size( );
  int origin = (*theStates.begin()).first;
  vector< vector< int > > IP = innerProductTable( G , origin );
  for( int a=0 ; a<N ; ++a ) {
    for( int b=0 ; b<N ; ++b ) {
      for( int c=0 ; c<N ; ++c ) {
	if( IP[a][b]==-1 || IP[a][c]==-1 || IP[b][c]==-1 )
	  continue;
	int val = IP[a][c] < IP[b][c] ? IP[a][c] : IP[b][c];
	val -= IP[a][b];
	if( result<val )
	  result = val;
      }
    }
  }
  
  return result/2;
}
