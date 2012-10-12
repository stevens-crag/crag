// Copyright (C) 2005 Alexander Ushakov
// Contents: Example for class Graph
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//


#include "GraphType.h"
#include "GraphConcept.h"


using namespace Graphs;


int main( )
{
  //& Graph ; How do I create an empty graph?
  Graph G;
  
  //& Graph ; How do I add new vertices to a graph?
  int v1 = G.newVertex( );
  int v2 = G.newVertex( );
  int v3 = G.newVertex( );
  int v4 = G.newVertex( );
  
  //& Graph ; How do I add new edges to a graph?
  G.newEdge( v1 , GraphEdge( v2 ) );      // add an edge v1 -> v2
  G.newEdge( v2 , GraphEdge( v3 ) );      // add an edge v2 -> v3
  G.newEdge( v1 , GraphEdge( v4 ) );      // add an edge v2 -> v3
  
  //& Graph ; How do I compute a geodesic tree with a fixed root?
  map< int , pres_type::edge_type > tree_out = getGeodesicTree_out( G , v1 );
  
  
  //& Graph ; How do I find geodesic path between two vertices?
  // using spanning tree already computed above to find a geodesic path between v1 and v3
  pair< bool , list< pres_type::edge_type > > path = readoffGeodesicTree( tree_out , v3 );
  if( !path.first ) {
    cout << "There is no path from " << v1 << " to " << v3 << endl;
  } else {
    cout << "Distance from " << v1 << " to " << v3 << " is " << path.second.size( ) << endl;
    cout << "The actual path passes through the vertices : ";
    for( list< pres_type::edge_type >::const_iterator p_it=path.second.end() ; p_it!=path.second.begin() ; )
      cout << *--p_it << " ";
    cout << endl;
  }
  
  //& Graph ; How do I find distances from a fixed vertex to all other vertices?
  map< int , int > D = getDistances_in( G , 0 );
  map< int , int >::const_iterator D_it = D.begin( );
  for( ; D_it!=D.end() ; ++D_it )
    cout << (*D_it).first << " -> " << (*D_it).second << endl;
  
  
  return 0;
}
