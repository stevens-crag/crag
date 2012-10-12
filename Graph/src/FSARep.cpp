// Copyright (C) 2005 Alexander Ushakov
// Contents: Implementation of class FSARep
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//

#include "Word.h"
#include "FSARep.h"
#include "iostream"


//---------------------------------------------------------------------------//
//-------------------------------- FSARep -----------------------------------//
//---------------------------------------------------------------------------//


FSARep::FSARep( ) :
  maxState( 0 )
{
  
}


int FSARep::newState( )
{
  theStates[maxState] = maxState;
  return maxState++;
}

void FSARep::newEdge( int state1 , int state2 , int label )
{
  theStates[state1].out.insert( FSAEdge(state2, label) );
  theStates[state2]. in.insert( FSAEdge(state1, label) );
}

void FSARep::eraseEdge( int state1 , int state2 , int label )
{
  theStates[state1].out.erase( FSAEdge(state2, label) );
  theStates[state2]. in.erase( FSAEdge(state1, label) );
}

void FSARep::eraseState( int state )
{
  set< FSAEdge > in  = theStates[state]. in;
  set< FSAEdge > out = theStates[state].out;

  set< FSAEdge >::iterator it2;

  // state <- (*it2).first by (*it2).second
  for( it2 = in.begin( ) ; it2!=in.end( ) ; ++it2 ) {
    FSAEdge p( state , (*it2).label );
    theStates[(*it2).target].out.erase( p );
  }

  // state2 -> (*it2).first by (*it2).second
  for( it2 = out.begin( ) ; it2!=out.end( ) ; ++it2 ) {
    FSAEdge p( state , (*it2).label );
    theStates[(*it2).target].in.erase( p );
  }  
  theStates.erase( state );
}


//---------------------------------------------------------------------------//
//--------------------------------- fold ------------------------------------//
//---------------------------------------------------------------------------//


void FSARep::fold( const set< int >* candidates , list< FoldDetails >* details )
{
  set< int > states_to_check;
  if( candidates ) {
    states_to_check = *candidates;
  } else {
    for( map< int , FSAState >::iterator s_it = theStates.begin( ) ; s_it!=theStates.end( ) ; ++s_it )
      states_to_check.insert( (*s_it).first );
  }

  while( !states_to_check.empty( ) ) {
    
    // take a first state in the set and remove it
    int state = *states_to_check.begin( );
    states_to_check.erase( states_to_check.begin( ) );
    map< int , FSAState >::iterator s_it = theStates.find( state );
    if(  s_it==theStates.end( ) )
      continue;
    
    // consider all edges leaving the state
    set< FSAEdge >& out = (*s_it).second.out;
    set< FSAEdge >::iterator o_it1 = out.begin( );
    set< FSAEdge >::iterator o_it2 = ++(out.begin( ));
    while( o_it2!=out.end( ) ) {
      
      if( (*o_it1).label==(*o_it2).label ) {
	if( details )
	  details->push_back( FoldDetails( (*s_it).first  ,  (*o_it1).label  ,  (*o_it1).target , theStates[(*o_it1).target] , theStates[(*o_it2).target] ) );
	pinch( (*o_it1).target , (*o_it2).target );
	states_to_check.insert( (*o_it1).target );
	states_to_check.insert( state );
	break;
      }
      ++o_it1;
      ++o_it2;
    }
    
  }
}


//---------------------------------------------------------------------------//
//-------------------------------- pinch ------------------------------------//
//---------------------------------------------------------------------------//


void FSARep::pinch( int state1 , int state2 )
{
  if( state1==state2 ) return;
  if( state1 >state2 )
    swap( state1 , state2 );

  set< FSAEdge > in  = theStates[state2]. in;
  set< FSAEdge > out = theStates[state2].out;

  set< FSAEdge >::iterator it2;

  // state2 <- (*it2).first by (*it2).second
  for( it2 = in.begin( ) ; it2!=in.end( ) ; ++it2 ) {
    if( (*it2).target!=state2 ) {
      theStates[state1].in.insert( (*it2) );
      FSAEdge p( state1 , (*it2).label );
      theStates[(*it2).target].out.insert( p );
    } else {
      FSAEdge p( state1 , (*it2).label );
      theStates[state1].in.insert( p );
      theStates[state1].out.insert( p );
    }
  }

  // state2 -> (*it2).first by (*it2).second
  for( it2 = out.begin( ) ; it2!=out.end( ) ; ++it2 ) {
    if( (*it2).target!=state2 ) {
      theStates[state1].out.insert( (*it2) );
      FSAEdge p( state1 , (*it2).label );
      theStates[(*it2).target].in.insert( p );
    } else {
      FSAEdge p( state1 , (*it2).label );
      theStates[state1].out.insert( p );
      theStates[state1].in.insert( p );
    }
  }

  eraseState( state2 );
}


//---------------------------------------------------------------------------//
//-------------------------------- unfold -----------------------------------//
//---------------------------------------------------------------------------//


void FSARep::unfold( const list< FoldDetails >& details )
{
	list< FoldDetails >::const_iterator det_it = details.end( );
	for( ; det_it!=details.begin( ) ; ) {

		--det_it;
		const FoldDetails& FD = *det_it;

		// erase folded state
		eraseState( FD.state_num );

		// insert two old states
	  theStates[FD.state1.theState] = FD.state1.theState;
	  theStates[FD.state2.theState] = FD.state2.theState;

		// insert old edges
		set< FSAEdge >::iterator edge_it = FD.state1.out.begin( );
		for( ; edge_it!=FD.state1.out.end( ) ; ++edge_it )
			newEdge( FD.state1.theState , (*edge_it).target , (*edge_it).label );
		edge_it = FD.state1.in.begin( );
		for( ; edge_it!=FD.state1.in.end( ) ; ++edge_it )
			newEdge( (*edge_it).target , FD.state1.theState , (*edge_it).label );

		edge_it = FD.state2.out.begin( );
		for( ; edge_it!=FD.state2.out.end( ) ; ++edge_it )
			newEdge( FD.state2.theState , (*edge_it).target , (*edge_it).label );
		edge_it = FD.state2.in.begin( );
		for( ; edge_it!=FD.state2.in.end( ) ; ++edge_it )
			newEdge( (*edge_it).target , FD.state2.theState , (*edge_it).label );
	}
}


//---------------------------------------------------------------------------//
//-------------------------------- liftup -----------------------------------//
//---------------------------------------------------------------------------//


void FSARep::liftup( const list< FoldDetails >& details , list< FSAEdge >& path , int init_state )
{
	int term_state = (*--path.end()).target;

	list< FoldDetails >::const_iterator det_it = details.end( );
	for( ; det_it!=details.begin( ) ; ) {

		--det_it;
		FoldDetails FD = *det_it;

		// cout << " +++ " << FD.state1.theState << ", " << FD.state2.theState << " +++ " << endl;
		// erase folded state
		eraseState( FD.state_num );

		// insert two old states
	  theStates[FD.state1.theState] = FD.state1.theState;
	  theStates[FD.state2.theState] = FD.state2.theState;

		// insert old edges
		set< FSAEdge >::iterator edge_it = FD.state1.out.begin( );
		for( ; edge_it!=FD.state1.out.end( ) ; ++edge_it )
			newEdge( FD.state1.theState , (*edge_it).target , (*edge_it).label );
		edge_it = FD.state1.in.begin( );
		for( ; edge_it!=FD.state1.in.end( ) ; ++edge_it )
			newEdge( (*edge_it).target , FD.state1.theState , (*edge_it).label );

		edge_it = FD.state2.out.begin( );
		for( ; edge_it!=FD.state2.out.end( ) ; ++edge_it )
			newEdge( FD.state2.theState , (*edge_it).target , (*edge_it).label );
		edge_it = FD.state2.in.begin( );
		for( ; edge_it!=FD.state2.in.end( ) ; ++edge_it )
			newEdge( (*edge_it).target , FD.state2.theState , (*edge_it).label );

		// liftup the path
		list< FSAEdge > path2;
		int cur_state = init_state;
		list< FSAEdge >::iterator path_it = path.begin( );
		// cout << "  *" << endl;
		for( ; path_it!=path.end( ) ; ) {

			int label = (*path_it).label;
			int target = (*path_it).target;
			// cout << "  " << cur_state << " -> " << label << " -> " << target << endl;
			const FSAState& S = theStates[cur_state];
			const set< FSAEdge >& out = S.out;

			// 1. check if we can lift up an edge "directly"
			set< FSAEdge >::const_iterator out_it = out.find( *path_it );
			if( out_it!=out.end( ) ) {
				path2.push_back( *out_it );
				cur_state = target;
				++path_it;
				continue;
			}

			// 2. check if we can switch a target
			if( target==FD.state1.theState || target==FD.state2.theState ) {
				FSAEdge mirror_edge( target==FD.state1.theState ? FD.state2.theState : FD.state1.theState , label );
				out_it = out.find( mirror_edge );
				if( out_it!=out.end( ) ) {
					path2.push_back( mirror_edge );
					cur_state = mirror_edge.target;
					++path_it;
					continue;
				}
			}
			
			// 3. check if we can switch a base (current) state
			if( cur_state==FD.state1.theState ) {
				path2.push_back( FSAEdge( FD.base , -FD.label ) );
				path2.push_back( FSAEdge( FD.state2.theState , FD.label ) );
				cur_state = FD.state2.theState;
			} else if( cur_state==FD.state2.theState ) {
				path2.push_back( FSAEdge( FD.base , -FD.label ) );
				path2.push_back( FSAEdge( FD.state1.theState , FD.label ) );
				cur_state = FD.state1.theState;
			}

		}

		path = path2;
		if( term_state!=(*--path.end()).target ) {
			if( term_state==FD.state1.theState ) {
				path.push_back( FSAEdge( FD.base , -FD.label ) );
				path.push_back( FSAEdge( FD.state1.theState , FD.label ) );
			} else {
				path.push_back( FSAEdge( FD.base , -FD.label ) );
				path.push_back( FSAEdge( FD.state2.theState , FD.label ) );
			}
		}
	}
}

/*
void reducePath( list< FSAEdge >& path , int init_state )
{
	if( path.size()<2 )
		return;

	int cur_state = init_state;
	list< int > states;
	list< FSAEdge >::iterator path_it1 = path.begin( );
	for( ; path_it1!=path.end( ) ; ) {

		// find the edge after the current one
		list< FSAEdge >::iterator path_it2 = path_it1;
		++path_it2;
		if( path_it2==path.end( ) )
			break;

		if( cur_state==(*path_it2).target && (*path_it1).label==-(*path_it2).label ) {

			// cout << "     ~~~R~~~" << endl;
			path.erase( path_it1 );
			path_it1 = path.erase( path_it2 );
			if( path_it1!=path.begin( ) ) {
				--path_it1;
				states.pop_back( );
			}
		} else {
			++path_it1;
			states.push_back( cur_state );
			cur_state = (*path_it1).target;
		}
	}

}
*/
