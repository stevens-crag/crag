
#include "fstream"
#include "GraphType.h"
#include "GraphConcept.h"

using namespace std;
using namespace Graphs;

void graph_example()
{
  typedef GraphVertex< GraphEdge >                vertex_type;
  typedef GraphEdge                               edge_type;
  typedef GraphConcept< vertex_type , edge_type > pres_type;
  
  GraphConcept< vertex_type , edge_type > G;
  
  // generate a random graph
  int N = 16;
  for( int i=0 ; i<16 ; ++i ) G.newVertex( vertex_type( ) );
  float edge_param = 2.0/N;
  for( int i=0 ; i<N ; ++i )
    for( int j=0 ; j<N ; ++j )
      if( RandLib::ur.rand()<edge_param )
	G.newEdge( i , edge_type( j ) );
  
  ofstream OF( "g.txt" );
  OF << graphviz_format( G ) << endl;
}


void labeled_graph_example( )
{
  typedef IntLabeledGraph::edge_type edge_type;
  typedef IntLabeledGraph::vertex_type vertex_type;

  IntLabeledGraph G;

  int N = 16;
  for( int i=0 ; i<N ; ++i ) G.newVertex( vertex_type( ) );
  float edge_param = 2.0/N;
  for( int i=0 ; i<N ; ++i )
    for( int j=0 ; j<N ; ++j )
      if( RandLib::ur.rand()<edge_param )
	G.newEdge( i , edge_type( j , (2*RandLib::ur.irand(0,1)-1) * RandLib::ur.irand( 1 , 2 ) ) );

  ofstream OF( "f.txt" );
  OF << graphviz_format( G ) << endl;
}


int main( )
{
  // graph_example( );
  labeled_graph_example( );
  
  
  return 0;
}
