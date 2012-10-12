// Copyright (C) 2005 Alexander Ushakov
// Contents: Definitions of Graph algorithms
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//

#ifndef _GraphAlgorithms_h_
#define _GraphAlgorithms_h_


#include "iostream"
#include "map"
#include "set"
#include "list"
#include <stdlib.h>

using namespace std;


//---------------------------------------------------------------------------//
//-------------------------- getGeodesicTree_in -----------------------------//
//---------------------------------------------------------------------------//


template < class Graph >
map< int , typename Graph::edge_type > getGeodesicTree_in( const Graph& graph , int init_st )
{
  typedef typename Graph::state_type state_type;
  typedef typename Graph::edge_type edge_type;
  map< int , edge_type > result;

  const map< int , state_type >& theStates = graph.getStates( );

  if( theStates.find( init_st )==theStates.end( ) ) {
   cout << "Nonexisting initial state. (9482)" << endl;
   exit( 1 );
  }

  set< int > checked_pts;
  set< int > cur_pts;
  set< int > new_pts;
  new_pts.insert( init_st );
  result[init_st] = edge_type( );

  while( new_pts.size( ) ) {

    cur_pts = new_pts;
    new_pts.clear( );
    while( cur_pts.size( ) ) {

      int cur_st = *cur_pts.begin( );
      cur_pts.erase( cur_pts.begin( ) );
      checked_pts.insert( cur_st );
	
      typename map< int , state_type >::const_iterator st_it = theStates.find( cur_st );
      if( st_it==theStates.end( ) ) {
        cout << "Unexpected situation in getGeodesicTree. (0276)" << endl;
        exit( 1 );
      }
      const state_type& cur_state = (*st_it).second;
      const set< edge_type >& in = cur_state.in;
      typename set< edge_type >::const_iterator in_it = in.begin( );
      for( ; in_it!=in.end( ) ; ++in_it ) {

        int new_st = (*in_it).target;
        if( checked_pts.find(new_st)==checked_pts.end() && cur_pts.find(new_st)==cur_pts.end() && new_pts.find(new_st)==new_pts.end() ) {
          new_pts.insert( new_st );
          edge_type edge = *in_it;
          edge.target = cur_st;
          result[new_st] = edge;
        }
      }
    }
  }

  return result;
}


//---------------------------------------------------------------------------//
//-------------------------- getGeodesicTree_out ----------------------------//
//---------------------------------------------------------------------------//


template < class Graph >
map< int , typename Graph::edge_type > getGeodesicTree_out( const Graph& graph , int init_st )
{
  typedef typename Graph::state_type state_type;
  typedef typename Graph::edge_type edge_type;
  map< int , edge_type > result;

  const map< int , state_type >& theStates = graph.getStates( );

  if( theStates.find( init_st )==theStates.end( ) ) {
   cout << "Nonexisting initial state. (9482)" << endl;
   exit( 1 );
  }

  set< int > checked_pts;
  set< int > cur_pts;
  set< int > new_pts;
  new_pts.insert( init_st );
  result[init_st] = edge_type( );

  while( new_pts.size( ) ) {

    cur_pts = new_pts;
    new_pts.clear( );
    while( cur_pts.size( ) ) {

      int cur_st = *cur_pts.begin( );
      cur_pts.erase( cur_pts.begin( ) );
      checked_pts.insert( cur_st );

      typename map< int , state_type >::const_iterator st_it = theStates.find( cur_st );
      if( st_it==theStates.end( ) ) {
        cout << "Unexpected situation in getGeodesicTree. (0276)" << endl;
        exit( 1 );
      }
      const state_type& cur_state = (*st_it).second;
      const set< edge_type >& out = cur_state.out;
      typename set< edge_type >::const_iterator out_it = out.begin( );
      for( ; out_it!=out.end( ) ; ++out_it ) {

        int new_st = (*out_it).target;
        if( checked_pts.find(new_st)==checked_pts.end() && cur_pts.find(new_st)==cur_pts.end() && new_pts.find(new_st)==new_pts.end() ) {
          new_pts.insert( new_st );
          edge_type edge = *out_it;
          edge.target = cur_st;
          result[new_st] = edge;
        }
      }
    }
  }

  return result;
}


//---------------------------------------------------------------------------//
//-------------------------- getGeodesicTree_out ----------------------------//
//---------------------------------------------------------------------------//


template < class Graph >
map< int , int > getDistances_out( const Graph& graph , int init_st )
{
  typedef typename Graph::state_type state_type;
  typedef typename Graph::edge_type edge_type;
  map< int , int > result;
  
  const map< int , state_type >& theStates = graph.getStates( );
  
  if( theStates.find( init_st )==theStates.end( ) ) {
   cout << "Nonexisting initial state. (9432)" << endl;
   exit( 1 );
  }
  
  set< int > checked_pts;
  set< int > cur_pts;
  set< int > new_pts;
  new_pts.insert( init_st );
  result[init_st] = 0;
  
  for( int dist=1 ; new_pts.size( ) ; ++dist ) {
    cur_pts = new_pts;
    new_pts.clear( );
    while( cur_pts.size( ) ) {

      int cur_st = *cur_pts.begin( );
      cur_pts.erase( cur_pts.begin( ) );
      checked_pts.insert( cur_st );

      typename map< int , state_type >::const_iterator st_it = theStates.find( cur_st );
      if( st_it==theStates.end( ) ) {
        cout << "Unexpected situation in getGeodesicTree. (89766)" << endl;
        exit( 1 );
      }
      const state_type& cur_state = (*st_it).second;
      const set< edge_type >& out = cur_state.out;
      typename set< edge_type >::const_iterator out_it = out.begin( );
      for( ; out_it!=out.end( ) ; ++out_it ) {

        int new_st = (*out_it).target;
        if( checked_pts.find(new_st)==checked_pts.end() && cur_pts.find(new_st)==cur_pts.end() && new_pts.find(new_st)==new_pts.end() ) {
          new_pts.insert( new_st );
          result[new_st] = dist;
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
list< edge_type > readoffGeodesicTree( const map< int , edge_type >& tree , int st_num )
{
  list< edge_type > result;

  while( 1 ) {
    // cout << "st_num = " << st_num << endl;
    typename map< int , edge_type >::const_iterator t_it = tree.find( st_num );
    if( t_it==tree.end( ) ) {
      cout << "Error. Unpredicted situation in readoffGeodesicTree( ... )" << endl;
      exit( 1 );
    }
    if( (*t_it).second.target==-1 )
      break;
    result.push_back( (*t_it).second );
    st_num = (*t_it).second.target;
  }

  return result;
}


//---------------------------------------------------------------------------//
//---------------------------------- trace ----------------------------------//
//---------------------------------------------------------------------------//


template < class LabelledGraph , class ConstIterator >
int trace( const LabelledGraph& LG , int init , ConstIterator B , ConstIterator E )
{
  typedef typename LabelledGraph::state_type state_type;
  typedef typename LabelledGraph::edge_type edge_type;
  const map< int , state_type >& theStates = LG.getStates( );

	int cur_st = init;
	for( ; B!=E ; ++B ) {

    typename map< int , state_type >::const_iterator st_it = theStates.find( cur_st );
		if( st_it==theStates.end( ) )
			return -1;

		bool success = false;
		int label = *B;
    const state_type& cur_state = (*st_it).second;
    const set< edge_type >& out = cur_state.out;
		typename set< edge_type >::const_iterator out_it = out.begin( );
    for( ; out_it!=out.end( ) ; ++out_it ) {
			if( (*out_it).label == label ) {
				cur_st = (*out_it).target;
				success = true;
				break;
			}
		}
		if( !success )
			return -1;
	}

	return cur_st;
}


//---------------------------------------------------------------------------//
//---------------------------------- trace ----------------------------------//
//---------------------------------------------------------------------------//


template < class LabelledGraph , class ConstIterator >
pair< bool , list< typename LabelledGraph::edge_type > > trace_path( const LabelledGraph& LG , int init , ConstIterator B , ConstIterator E )
{
  typedef typename LabelledGraph::state_type state_type;
  typedef typename LabelledGraph::edge_type edge_type;
  const map< int , state_type >& theStates = LG.getStates( );

	int cur_st = init;
	pair< bool , list< edge_type > > result;
	result.first = true;
	for( ; B!=E ; ++B ) {

    typename map< int , state_type >::const_iterator st_it = theStates.find( cur_st );
		if( st_it==theStates.end( ) )
			return pair< bool , list< edge_type > >( false , list< edge_type >( ) );

		bool success = false;
		int label = *B;
    const state_type& cur_state = (*st_it).second;
    const set< edge_type >& out = cur_state.out;
		typename set< edge_type >::const_iterator out_it = out.begin( );
    for( ; out_it!=out.end( ) ; ++out_it ) {
			if( (*out_it).label == label ) {
				cur_st = (*out_it).target;
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


#endif

