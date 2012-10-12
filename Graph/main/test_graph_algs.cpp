
#include "Word.h"

#include "GraphType.h"
#include "GraphConcept.h"
#include "GraphConceptAlgorithms.h"
using namespace Graphs;

/*
void test_express( )
{
	int N=2;
	vector< Word > gens;
	for( int i=0 ; i<5; ++i )
		gens.push_back( Word::randomWord( N , 3 ) );
	SubgroupFG sbgp( N , gens );
	if( sbgp.getIndex( )!=1 )
		return;

	Word decomp = sbgp.express( Word(1) );
}
*/


void test_unfold( )
{
  typedef IntLabeledGraph Graph;
  typedef Graph::edge_type edge_type;
  typedef Graph::vertex_type vertex_type;
  Graph F;

  int o = F.newVertex( );
  Word x1( 1 );
  Word x2( 2 );
  for( int i=0 ; i<5; ++i ) {
    Word w = Word::randomWord( 2 , 3 );
    addLoop( F , o , w.begin( ) , w.end( ) );
  }
  
  set< int > candidates;
  candidates.insert( o );
  list< FoldDetails< vertex_type , edge_type > > details;
  fold( F , candidates , &details );
  if( F.getVertices( ).size( )>1 )
    return;
  
  pair< bool , list< edge_type > > path = trace_path( F , o , x2.begin( ) , x2.end( ) );
  liftup( F , o , path.second , details.begin( ) , details.end( ) );
}


int main( )
{
  srand( time(0) );
  test_unfold( );
  return 0;
  
  
  Graph G;
  int v1 = G.newVertex( );
  int v2 = G.newVertex( );
  int v3 = G.newVertex( );
  int v4 = G.newVertex( );
  G.newEdge( v1 , v2 );
  G.newEdge( v1 , v3 );
  G.newEdge( v3 , v4 );
  
  cout << "G = " << G << endl;
  
  {
    map< int , Graph::edge_type > spanning_tree = getGeodesicTree_in( G , v1 );
    map< int , Graph::edge_type >::iterator t_it = spanning_tree.begin( );
    for( ; t_it!=spanning_tree.end( ) ; ++t_it ) {
      cout << (*t_it).first << " -> " << (*t_it).second.theTarget << endl;
    }
  }
  {
    map< int , Graph::edge_type > spanning_tree = getGeodesicTree_out( G , v1 );
    map< int , Graph::edge_type >::iterator t_it = spanning_tree.begin( );
    for( ; t_it!=spanning_tree.end( ) ; ++t_it ) {
      cout << (*t_it).first << " -> " << (*t_it).second.theTarget << endl;
    }
    // finds a path from v4 to v1
    pair< bool , list< Graph::edge_type > > path = readoffGeodesicTree( spanning_tree , v4 );
    for( list< Graph::edge_type >::iterator p_it = path.second.begin( ) ; p_it!=path.second.end( ) ; ++p_it )
      cout << (*p_it).theTarget << ", ";
    cout << endl;
  }
  
  
  
  return 0;
}
