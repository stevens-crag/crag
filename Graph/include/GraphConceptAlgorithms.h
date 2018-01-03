// Copyright (C) 2005 Alexander Ushakov
// Contents: Definition of the graph concept algorithms
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//

#include <iostream>
#include <map>
#include <set>
#include <list>
#include <stdlib.h>

using namespace std;


#ifndef _GraphConceptAlgorithms_H_
#define _GraphConceptAlgorithms_H_


namespace Graphs
{

  //---------------------------------------------------------------------------//
  //-------------------------- getGeodesicTree_in -----------------------------//
  //---------------------------------------------------------------------------//
  //! Compute directed geodesic tree "inward" to the vertex init_v.
  /*!
    Returns a list of edges that form a a geodesic subtree of a graph directed to the vertex init_vert. 
    If init_vert is reachable from V if and only if result[V] is defined.
   */
  template < class Graph > map< int , typename Graph::edge_type > getGeodesicTree_in( const Graph& graph , int init_v )
  {
    typedef typename Graph::vertex_type vertex_type;
    typedef typename Graph::edge_type edge_type;
    map< int , edge_type > result;
    
    const map< int , vertex_type >& theVertices = graph.getVertices( );
    
    if( theVertices.find( init_v )==theVertices.end( ) ) {
      cout << "Nonexisting initial state. (9482)" << endl;
      exit( 1 );
    }
    
    set< int > checked_verts;
    set< int > cur_verts;
    set< int > new_verts;
    new_verts.insert( init_v );
    result[init_v] = edge_type( );
    
    while( new_verts.size( ) ) {
      
      cur_verts = new_verts;
      new_verts.clear( );
      while( cur_verts.size( ) ) {
	
	int cur_v = *cur_verts.begin( );
	cur_verts.erase( cur_verts.begin( ) );
	checked_verts.insert( cur_v );
	
	typename map< int , vertex_type >::const_iterator st_it = theVertices.find( cur_v );
	if( st_it==theVertices.end( ) ) {
	  cout << "Unexpected situation in getGeodesicTree. (0276)" << endl;
	  exit( 1 );
	}
	const vertex_type& cur_V = (*st_it).second;
	const set< edge_type >& in = cur_V.in;
	typename set< edge_type >::const_iterator in_it = in.begin( );
	for( ; in_it!=in.end( ) ; ++in_it ) {
	  
	  int new_v = (*in_it).theTarget;
	  if( checked_verts.find(new_v)==checked_verts.end() && cur_verts.find(new_v)==cur_verts.end() && new_verts.find(new_v)==new_verts.end() ) {
	    new_verts.insert( new_v );
	    edge_type edge = *in_it;
	    edge.theTarget = cur_v;
	    result[new_v] = edge;
	  }
	}
      }
    }
    
    return result;
  }
  
  //---------------------------------------------------------------------------//
  //-------------------------- getGeodesicTree_out ----------------------------//
  //---------------------------------------------------------------------------//
  //! Compute geodesic directed tree starting the vertex init_v.
  /*!
    Returns a list of edges that form a a geodesic subtree of a graph directed to the vertex init_vert. 
    If init_vert is reachable from V if and only if result[V] is defined.
   */
  template < class Graph >
    map< int , typename Graph::edge_type > getGeodesicTree_out( const Graph& graph , int init_v )
  {
    typedef typename Graph::vertex_type vertex_type;
    typedef typename Graph::edge_type   edge_type;
    map< int , edge_type > result;
    
    const map< int , vertex_type >& theVertices = graph.getVertices( );
    
    if( theVertices.find( init_v )==theVertices.end( ) ) {
      cout << "Nonexisting initial state. (9482)" << endl;
      exit( 1 );
    }
    
    set< int > checked_verts;
    set< int > cur_verts;
    set< int > new_verts;
    new_verts.insert( init_v );
    result[init_v] = edge_type( );
    
    while( new_verts.size( ) ) {
      
      cur_verts = new_verts;
      new_verts.clear( );
      while( cur_verts.size( ) ) {
	
	int cur_v = *cur_verts.begin( );
	cur_verts.erase( cur_verts.begin( ) );
	checked_verts.insert( cur_v );
	
	typename map< int , vertex_type >::const_iterator st_it = theVertices.find( cur_v );
	if( st_it==theVertices.end( ) ) {
	  cout << "Unexpected situation in getGeodesicTree. (0276)" << endl;
	  exit( 1 );
	}
	const vertex_type& cur_V = (*st_it).second;
	const set< edge_type >& out = cur_V.out;
	typename set< edge_type >::const_iterator out_it = out.begin( );
	for( ; out_it!=out.end( ) ; ++out_it ) {

	  int new_v = (*out_it).theTarget;
	  if( checked_verts.find(new_v)==checked_verts.end() && cur_verts.find(new_v)==cur_verts.end() && new_verts.find(new_v)==new_verts.end() ) {
	    new_verts.insert( new_v );
	    edge_type edge = *out_it;
	    edge.theTarget = cur_v;
	    result[new_v] = edge;
	  }
	}
      }
    }
    
    return result;
  }
  

  //---------------------------------------------------------------------------//
  //---------------------------- getDistances_out -----------------------------//
  //---------------------------------------------------------------------------//
  
  //! Compute distances from the vertex init_v.
  template < class Graph >
    map< int , int > getDistances_out( const Graph& graph , int init_v )
    {
      typedef typename Graph::vertex_type vertex_type;
      typedef typename Graph::edge_type edge_type;
      map< int , int > result;
      
      const map< int , vertex_type >& theVertices = graph.getVertices( );
      
      if( theVertices.find( init_v )==theVertices.end( ) ) {
	cout << "Nonexisting initial state. (9432)" << endl;
	exit( 1 );
      }
  
      set< int > checked_verts;
      set< int > cur_verts;
      set< int > new_verts;
      new_verts.insert( init_v );
      result[init_v] = 0;
  
      for( int dist=1 ; new_verts.size( ) ; ++dist ) {
	cur_verts = new_verts;
	new_verts.clear( );
	while( cur_verts.size( ) ) {
	  
	  int cur_v = *cur_verts.begin( );
	  cur_verts.erase( cur_verts.begin( ) );
	  checked_verts.insert( cur_v );
	  
	  typename map< int , vertex_type >::const_iterator st_it = theVertices.find( cur_v );
	  if( st_it==theVertices.end( ) ) {
	    cout << "Unexpected situation in getGeodesicTree. (89766)" << endl;
	    exit( 1 );
	  }
	  const vertex_type& cur_V = (*st_it).second;
	  const set< edge_type >& out = cur_V.out;
	  typename set< edge_type >::const_iterator out_it = out.begin( );
	  for( ; out_it!=out.end( ) ; ++out_it ) {
	    
	    int new_v = (*out_it).theTarget;
	    if( checked_verts.find(new_v)==checked_verts.end() && cur_verts.find(new_v)==cur_verts.end() && new_verts.find(new_v)==new_verts.end() ) {
	      new_verts.insert( new_v );
	      result[new_v] = dist;
	    }
	  }
	}
      }
      
      return result;
    }
  

  //---------------------------------------------------------------------------//
  //----------------------------- getDistances_in -----------------------------//
  //---------------------------------------------------------------------------//
  

  //! Compute distances to the vertex init_v.
  template < class Graph >
    map< int , int > getDistances_in( const Graph& graph , int init_v )
    {
      typedef typename Graph::vertex_type vertex_type;
      typedef typename Graph::edge_type edge_type;
      map< int , int > result;
      
      const map< int , vertex_type >& theVertices = graph.getVertices( );
      
      if( theVertices.find( init_v )==theVertices.end( ) ) {
	cout << "Nonexisting initial state. (9432)" << endl;
	exit( 1 );
      }
  
      set< int > checked_verts;
      set< int > cur_verts;
      set< int > new_verts;
      new_verts.insert( init_v );
      result[init_v] = 0;
  
      for( int dist=1 ; new_verts.size( ) ; ++dist ) {
	cur_verts = new_verts;
	new_verts.clear( );
	while( cur_verts.size( ) ) {
	  
	  int cur_v = *cur_verts.begin( );
	  cur_verts.erase( cur_verts.begin( ) );
	  checked_verts.insert( cur_v );
	  
	  typename map< int , vertex_type >::const_iterator st_it = theVertices.find( cur_v );
	  if( st_it==theVertices.end( ) ) {
	    cout << "Unexpected situation in getGeodesicTree. (89766)" << endl;
	    exit( 1 );
	  }
	  const vertex_type& cur_V = (*st_it).second;
	  const set< edge_type >& in = cur_V.in;
	  typename set< edge_type >::const_iterator in_it = in.begin( );
	  for( ; in_it!=in.end( ) ; ++in_it ) {
	    
	    int new_v = (*in_it).theTarget;
	    if( checked_verts.find(new_v)==checked_verts.end() && cur_verts.find(new_v)==cur_verts.end() && new_verts.find(new_v)==new_verts.end() ) {
	      new_verts.insert( new_v );
	      result[new_v] = dist;
	    }
	  }
	}
      }
      
      return result;
    }


  //---------------------------------------------------------------------------//
  //-------------------------- readoffGeodesicTree ----------------------------//
  //---------------------------------------------------------------------------//
  
  
  //! Function finds a path in the tree starting from the state st_num to the root. 
  template < class edge_type >
    pair< bool , list< edge_type > > readoffGeodesicTree( const map< int , edge_type >& tree , int v ) 
    {
      list< edge_type > result;
      
      while( 1 ) {
	typename map< int , edge_type >::const_iterator t_it = tree.find( v );
	if( t_it==tree.end( ) )
	  return pair< bool , list< edge_type > >( false , list< edge_type >( ) );
	if( (*t_it).second.theTarget==-1 )
	  break;
	result.push_back( (*t_it).second );
	v = (*t_it).second.theTarget;
      }
      
      return pair< bool , list< edge_type > >( true , result );
    }
  
  //---------------------------------------------------------------------------//
  //---------------------------------- trace ----------------------------------//
  //---------------------------------------------------------------------------//
  
  
  template < class LabelledGraph , class ConstIterator >
    int trace( const LabelledGraph& LG , int init_v , ConstIterator B , ConstIterator E )
  {
    typedef typename LabelledGraph::vertex_type vertex_type;
    typedef typename LabelledGraph::edge_type edge_type;
    const map< int , vertex_type >& theVertices = LG.getVertices( );
    
    int cur_v = init_v;
    for( ; B!=E ; ++B ) {
      
      typename map< int , vertex_type >::const_iterator v_it = theVertices.find( cur_v );
      if( v_it==theVertices.end( ) )
	return -1;
      
      bool success = false;
      int label = *B;
      const vertex_type& cur_V = (*v_it).second;
      const set< edge_type >& out = cur_V.out;
      typename set< edge_type >::const_iterator out_it = out.begin( );
      for( ; out_it!=out.end( ) ; ++out_it ) {
	if( (*out_it).theLabel == label ) {
	  cur_v = (*out_it).theTarget;
	  success = true;
	  break;
	}
      }
      if( !success )
	return -1;
    }
    
    return cur_v;
  }


  //---------------------------------------------------------------------------//
  //-------------------------------- trace_path -------------------------------//
  //---------------------------------------------------------------------------//
  
  
  template < class LabelledGraph , class ConstIterator >
    pair< bool , list< typename LabelledGraph::edge_type > > trace_path( const LabelledGraph& LG , int init_v , ConstIterator B , ConstIterator E )
    {
      typedef typename LabelledGraph::vertex_type vertex_type;
      typedef typename LabelledGraph::edge_type edge_type;
      const map< int , vertex_type >& theVertices = LG.getVertices( );
      
      int cur_v = init_v;
      pair< bool , list< edge_type > > result;
      result.first = true;
      for( ; B!=E ; ++B ) {
	
	typename map< int , vertex_type >::const_iterator v_it = theVertices.find( cur_v );
	if( v_it==theVertices.end( ) )
	  return pair< bool , list< edge_type > >( false , list< edge_type >( ) );
	
	bool success = false;
	int label = *B;
	const vertex_type& cur_V = (*v_it).second;
	const set< edge_type >& out = cur_V.out;
	typename set< edge_type >::const_iterator out_it = out.begin( );
	for( ; out_it!=out.end( ) ; ++out_it ) {
	  if( (*out_it).theLabel == label ) {
	    cur_v = (*out_it).theTarget;
	    success = true;
	    result.second.push_back( *out_it );
	    break;
	  }
	}
	if( !success )
	  return pair< bool , list< edge_type > >( false , list< edge_type >( ) );
      }
      
      return result;
    }
  
  
  //---------------------------------------------------------------------------//
  //-------------------------------- FoldDetails ------------------------------//
  //---------------------------------------------------------------------------//

  /*! 
    It is assumed that theEdge1.theTarget < theEdge2.theTarget.
   */
  template < class VertexType , class EdgeType >
    struct FoldDetails {
      
      FoldDetails( int o , const EdgeType& e1 , const EdgeType& e2 , const VertexType& V1 , const VertexType& V2 ) : 
	theOrigin(o), theEdge1(e1), theEdge2(e2), theVertex1(V1), theVertex2(V2) { }
      
      // int newVertex;  // = min{ theEdge1.theTarget , theEdge2.theTarget }
      // int oldVertex1; // = theEdge1.theTarget;
      // int oldVertex2; // = theEdge2.theTarget;
      
      int theOrigin;
      EdgeType theEdge1; // edge: theOrigin -> theVertex1 
      EdgeType theEdge2; // edge: theOrigin -> theVertex2
      
      VertexType theVertex1;
      VertexType theVertex2;
    };
  
  //---------------------------------------------------------------------------//
  //------------------------------------ fold ---------------------------------//
  //---------------------------------------------------------------------------//
  
  //! Fold a labelled graph. To aplly this function the graph edges must have theLabel defined. 
  /*!
    The graph is not folded if there is a vertex and two different edges leaving it equally labelled.
    In that event "fold" pinches end-points of these edges.
  */
  template < class LabelledGraph >
    void fold( LabelledGraph& G , set< int > candidates , 
	       list< FoldDetails< typename LabelledGraph::vertex_type , typename LabelledGraph::edge_type > >* details = NULL )
    {
      typedef typename LabelledGraph::vertex_type vertex_type;
      typedef typename LabelledGraph::edge_type edge_type;
      typedef FoldDetails< vertex_type , edge_type > FD;
      const map< int , vertex_type >& theVertices = G.getVertices( );


      while( !candidates.empty( ) ) {
	int v = *candidates.begin( );
	candidates.erase( v );
	if( theVertices.find(v)==theVertices.end( ) ) continue;
	const vertex_type& V = (*theVertices.find( v )).second;

	const set< edge_type >& out = V.out;
	typename set< edge_type >::const_iterator e_it = out.begin( );
	for( ; e_it!=out.end( ) ; ++e_it ) {
	  typename set< edge_type >::const_iterator e_it2 = e_it;
	  ++e_it2;
	  if( e_it2==out.end( ) )
	    break;
	  if( (*e_it).theLabel==(*e_it2).theLabel ) {
	    int t1 = (*e_it).theTarget;
	    int t2 = (*e_it2).theTarget;
	    if( details ) {
	      const vertex_type& V1 = (*theVertices.find( t1 )).second;
	      const vertex_type& V2 = (*theVertices.find( t2 )).second;
	      if( t1<t2 )
		details->push_back( FD( v , *e_it , *e_it2 , V1 , V2 ) );
	      else 
		details->push_back( FD( v , *e_it2 , *e_it , V2 , V1 ) );
	    }
	    G.pinch( t1 , t2 );
	    candidates.insert( t1<t2 ? t1:t2 );
	    candidates.insert( v );
	    break;
	  }
	}
      }
      
    }

  template < class LabelledGraph >
    void fold( const LabelledGraph& G , int candidate , 
	       list< FoldDetails< typename LabelledGraph::vertex_type , typename LabelledGraph::edge_type > >* details = NULL )
    {
      set< int > candidates;
      candidates.insert( candidate );
      fold( G , candidates , details );
    }

  template < class LabelledGraph >
    void fold( const LabelledGraph& G , 
	       list< FoldDetails< typename LabelledGraph::vertex_type , typename LabelledGraph::edge_type > >* details = NULL )
    {
      typedef typename LabelledGraph::vertex_type vertex_type;
      typedef typename LabelledGraph::edge_type edge_type;
      
      set< int > candidates;
      const map< int , vertex_type >& theVertices = G.getVertices( );
      typename map< int , vertex_type >::const_iterator v_it = theVertices.begin( );
      for( ; v_it!=theVertices.end( ) ; ++v_it )
	candidates.insert( *v_it );

      fold( G , candidates , details );
    }
  
  
  //---------------------------------------------------------------------------//
  //----------------------------------- unfold --------------------------------//
  //---------------------------------------------------------------------------//
  
  
  //! Revert the folding.
  template < class LabelledGraph , class FoldDetailsConstIterator >
    void unfold( LabelledGraph& G , FoldDetailsConstIterator B , FoldDetailsConstIterator E )
  {
    for( ; B!=E ; ) {
      --E;
      int oldVertex1 = (*E).theEdge1.theTarget;
      int oldVertex2 = (*E).theEdge2.theTarget;
      G.eraseVertex( oldVertex1 );
      G.replaceVertex( oldVertex1 , (*E).theVertex1 );
      G.replaceVertex( oldVertex2 , (*E).theVertex2 );
    }
  }  

  //---------------------------------------------------------------------------//
  //----------------------------------- liftup --------------------------------//
  //---------------------------------------------------------------------------//
  
  template < class edge_type > void reduce_path( int init_vertex , list< edge_type >& path );
  //! Revert the folding.
  template < class LabelledGraph , class FoldDetailsConstIterator >
    void liftup( const LabelledGraph& graph , int init_vertex , 
		 list< typename LabelledGraph::edge_type >& path ,
		 FoldDetailsConstIterator B , FoldDetailsConstIterator E )
  {
    LabelledGraph G = graph;
    typedef typename LabelledGraph::vertex_type vertex_type;
    typedef typename LabelledGraph::edge_type edge_type;
    
    for( ; B!=E ; ) {

      // unfold
      --E;
      int theOrigin = (*E).theOrigin;
      int oldVertex1 = (*E).theEdge1.theTarget;
      int oldVertex2 = (*E).theEdge2.theTarget;
      G.eraseVertex( oldVertex1 );
      G.replaceVertex( oldVertex1 , (*E).theVertex1 );
      G.replaceVertex( oldVertex2 , (*E).theVertex2 );

      const map< int , vertex_type >& theVertices = G.getVertices( );

      // liftup the current path edge by edge
      int cur_vertex = init_vertex;
      typename list< edge_type >::iterator path_it = path.begin( );
      for( ; path_it!=path.end( ) ; ++path_it ) {

	bool edge_found = false;
	edge_type edge = *path_it;
	int new_vertex = (*path_it).theTarget;
	for( int i=0 ; i<(cur_vertex==oldVertex1?2:1) && !edge_found ; ++i ) { 
	  // if the current vertex (origin of the current edge) is the split vertex then we make 2 steps
	  int v1 = (i==0 ? cur_vertex : oldVertex2);
	  
	  for( int j=0 ; j<(new_vertex==oldVertex1?2:1) && !edge_found ; ++j ) { 
	    // if the new vertex (target of the current edge) is the split vertex then we make 2 steps
	    int v2 = (j==0 ? new_vertex : oldVertex2);
	    
	    edge.theTarget = v2;
	    const vertex_type& V1 = (*theVertices.find( v1 )).second;
	    if( V1.out.find( edge )!=V1.out.end( ) ) { // the edge is found
	      
	      (*path_it).theTarget = v2;
	      if( i==1 ) { // initial vertex has jumped
		// insert a 2-edge path before the current edge connecting cur_vertex=oldVertex1 and oldVertex2
		path.insert( path_it , (*E).theEdge1.inverse( theOrigin ) );
		path.insert( path_it , (*E).theEdge2 );
	      }
	      
	      if( j==1 ) { // terminal vertex has jumped
		// insert a 2-edge path after the current edge connecting cur_vertex and removedVertex
		++path_it;
		path.insert( path_it , (*E).theEdge2.inverse( theOrigin ) );
		path.insert( path_it , (*E).theEdge1 );
		--path_it;
	      }
	      edge_found = true;
	    }
	  }
	}
	
	if( !edge_found ) {
	  cout << "Big shit" << endl;
	  exit(1);
	}
	cur_vertex = new_vertex;
      }
      reduce_path( init_vertex , path );
    }
  }

  //---------------------------------------------------------------------------//
  //----------------------------------- liftup --------------------------------//
  //---------------------------------------------------------------------------//
  
  //! Reduce a path (remove all pairs of consequent inverse edges).
  template < class edge_type >
    void reduce_path( int init_vertex , list< edge_type >& path )
  {
    list< int > vertices;
    vertices.push_back( init_vertex );
    typename list< edge_type >::iterator path_it = path.begin( );
    for( ; path_it!=path.end( ) ; ) {
      
      vertices.push_back( (*path_it).theTarget );
      if( path_it==path.begin( ) ) {
	++path_it;
	continue;
      }
      int old_vertex = *-- -- --vertices.end( );

      typename list< edge_type >::iterator path_it2 = path_it;
      --path_it2;
      if( *path_it!=(*path_it2).inverse( old_vertex ) ) {
	++path_it;
	continue;
      }

      vertices.pop_back( );
      vertices.pop_back( );
      path.erase( path_it2 );
      path_it = path.erase( path_it );
    }
  }
  
}


#endif
