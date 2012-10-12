
#include "GraphType.h"
#include "GraphConcept.h"
#include "GraphConceptAlgorithms.h"
#include "Word.h"

#include "RanlibCPP.h"


using namespace Graphs;


//---------------------------------------------------------------------------//
//----------------------------- graph_example -------------------------------//
//---------------------------------------------------------------------------//


void graph_example( )
{
  typedef GraphVertex< GraphEdge >                vertex_type;
  typedef GraphEdge                               edge_type;
  typedef GraphConcept< vertex_type , edge_type > pres_type;
  
  GraphConcept< vertex_type , edge_type > G;
  
  // generate a random graph
  int N = 16;
  for( int i=0 ; i<16 ; ++i ) G.newVertex( vertex_type( ) );
  float edge_param = 3.0/N;
  for( int i=0 ; i<N ; ++i )
    for( int j=0 ; j<N ; ++j )
      if( RandLib::ur.rand()<edge_param )
	G.newEdge( i , edge_type( j ) );

  // Output G
  cout << "==========================================" << endl;
  cout << "Graph G" << endl;
  cout << G << endl;

  cout << "==========================================" << endl;
  cout << "In tree" << endl;
  // Compute geodesic tree "inward" to the vertex 0 and output it
  map< int , pres_type::edge_type > tree_in = getGeodesicTree_in( G , 0 );
  map< int , pres_type::edge_type >::const_iterator t_it = tree_in.begin( );
  for( ; t_it!=tree_in.end( ) ; ++t_it ) {
    cout << "  " << (*t_it).first << " -> " << (*t_it).second << endl;
  }
  
  cout << "==========================================" << endl;
  cout << "Out tree" << endl;
  // Compute geodesic leaving the vertex n1 and output it
  map< int , pres_type::edge_type > tree_out = getGeodesicTree_out( G , 0 );
  t_it = tree_out.begin( );
  for( ; t_it!=tree_out.end( ) ; ++t_it ) {
    cout << "  " << (*t_it).first << " -> " << (*t_it).second << endl;
  }

  cout << "==========================================" << endl;
  cout << "Distances from the vertex" << endl;
  map< int , int > D = getDistances_out( G , 0 );
  map< int , int >::const_iterator D_it = D.begin( );
  for( ; D_it!=D.end() ; ++D_it )
    cout << (*D_it).first << " -> " << (*D_it).second << endl;

  cout << "==========================================" << endl;
  cout << "Distances to the vertex" << endl;
  D = getDistances_in( G , 0 );
  D_it = D.begin( );
  for( ; D_it!=D.end() ; ++D_it )
    cout << (*D_it).first << " -> " << (*D_it).second << endl;

  cout << "==========================================" << endl;
  cout << "Geodesic paths from the vertex" << endl;
  for( int i=0 ; i<N ; ++i ) {
    pair< bool , list< pres_type::edge_type > > path = readoffGeodesicTree( tree_out , i );
    if( !path.first ) {
      cout << "Distance from " << 0 << " to " << i << " is -1" << endl;
    } else {
      cout << "Distance from " << 0 << " to " << i << " is " << path.second.size( ) << ": ";
      for( list< pres_type::edge_type >::const_iterator p_it=path.second.end() ; p_it!=path.second.begin() ; )
	cout << *--p_it << " ";
      cout << endl;
    }
  }

  cout << "==========================================" << endl;
  cout << "Geodesic paths to the vertex" << endl;
  for( int i=0 ; i<N ; ++i ) {
    pair< bool , list< pres_type::edge_type > > path = readoffGeodesicTree( tree_in , i );
    if( !path.first ) {
      cout << "Distance from " << i << " to " << 0 << " is -1" << endl;
    } else {
      cout << "Distance from " << i << " to " << 0 << " is " << path.second.size( ) << ": ";
      for( list< pres_type::edge_type >::const_iterator p_it=path.second.begin() ; p_it!=path.second.end() ; ++p_it )
	cout << *p_it << " ";
      cout << endl;
    }
  }
  
  
  
}


//---------------------------------------------------------------------------//
//------------------------- labeled_graph_example ---------------------------//
//---------------------------------------------------------------------------//


void labeled_graph_example( )
{
  // typedef                             IntLabeledEdge  edge_type;
  // typedef              GraphVertex< IntLabeledEdge >  vertex_type;
  // typedef GraphConceptRep< vertex_type , edge_type >  pres_type;
  // GraphConcept< vertex_type , edge_type > G;

  typedef IntLabeledGraph Graph;
  typedef Graph::edge_type edge_type;
  typedef Graph::vertex_type vertex_type;

  Graph G;

  int N = 16;
  for( int i=0 ; i<N ; ++i ) G.newVertex( vertex_type( ) );
  float edge_param = 3.0/N;
  for( int i=0 ; i<N ; ++i )
    for( int j=0 ; j<N ; ++j )
      if( RandLib::ur.rand()<edge_param )
	G.newEdge( i , edge_type( j , (2*RandLib::ur.irand(0,1)-1) * RandLib::ur.irand( 1 , 2 ) ) );
  
  cout << G << endl;
  cout << "==========================================" << endl;

  // Compute geodesic tree "inward" to the vertex 0 and output it
  map< int , edge_type > tree = getGeodesicTree_in( G , 0 );
  map< int , edge_type >::const_iterator t_it = tree.begin( );
  for( ; t_it!=tree.end( ) ; ++t_it ) {
    cout << "  " << (*t_it).first << " -> " << (*t_it).second << endl;
  }
  cout << "==========================================" << endl;

  // Compute geodesic leaving the vertex n1 and output it
  map< int , edge_type > tree_out = getGeodesicTree_out( G , 0 );
  t_it = tree_out.begin( );
  for( ; t_it!=tree_out.end( ) ; ++t_it ) {
    cout << "  " << (*t_it).first << " -> " << (*t_it).second << endl;
  }

  cout << "==========================================" << endl;
  cout << "       Perform Folding:" << endl;
  set< int > candidates;
  for( int i=0 ; i<N ; ++i ) candidates.insert( i );
  list< FoldDetails< vertex_type , edge_type > > details;
  fold( G , candidates , &details );
  cout << G << endl;
  cout << "==========================================" << endl;
  cout << "       Perform Unfolding:" << endl;
  unfold( G , details.begin( ) , details.end( ) );
  cout << G << endl;
  cout << "==========================================" << endl;
  
  
}


void labeled_graph_folding_example( )
{
  typedef IntLabeledGraph Graph;
  typedef Graph::edge_type edge_type;
  typedef Graph::vertex_type vertex_type;

  Graph G;
  G.newVertex( vertex_type( ) );
  int alphabet = 2;
  
  Word u = Word::randomWord( alphabet , 5 );
  cout << u << endl;
  addRay( G , 0 , u.begin() , u.end( ) );
  for( int i=0 ; i<15 ; ++i ) {
    Word w = Word::randomWord( alphabet , 5 );
    cout << w << endl;
    addLoop( G , 0 , w.begin() , w.end( ) );
  }
  cout << G << endl;
  
  set< int > candidates;
  list< FoldDetails< vertex_type , edge_type > > details;
  candidates.insert( 0 );
  fold( G , candidates , &details );
  cout << G << endl;
  
  cout << "  Size of details is " << details.size( ) << endl;

  pair< bool , list< edge_type > > path = trace_path( G , 0 , u.begin( ) , u.end( ) );
  if( !path.first ) {
    cout << "Error! No path!" << endl;
  } else {
    list< edge_type >::const_iterator e_it = path.second.end( );
    --e_it;
    if( (*e_it).theTarget!=0 ) {
      cout << "Error! Path is not a loop" << endl;
    } else {
      cout << "Loop found" << endl;
      liftup( G , 0 , path.second , details.begin( ) , details.end( ) );
      cout << "Path length = " << path.second.size( ) << endl;
      list< edge_type >::const_iterator p_it = path.second.begin( );
      for( ; p_it!=path.second.end( ) ; ++p_it ) 
	cout << *p_it << "  ";
      cout << endl;
    }
  }
}


//---------------------------------------------------------------------------//
//--------------------------------- main ------------------------------------//
//---------------------------------------------------------------------------//


int main( )
{
  // graph_example( );
  // labeled_graph_example( );
  labeled_graph_folding_example( );

  return 0;
}
