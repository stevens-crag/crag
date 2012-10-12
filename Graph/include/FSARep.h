// Copyright (C) 2005 Alexander Ushakov
// Contents: Definition of class FSARep
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//

#ifndef _FSARep_h_
#define _FSARep_h_


#include "RefCounter.h"

#include "map"
#include "set"
#include "list"
#include "vector"
#include <stdlib.h>

using namespace std;


//---------------------------------------------------------------------------//
//-------------------------------- FSAEdge ----------------------------------//
//---------------------------------------------------------------------------//


struct FSAEdge 
{

  FSAEdge( ) : target(-1), label(0) { }
  // dummy constructor, not to be used, required for STL::map

  FSAEdge( int t , int l ) : target(t), label(l) { }

  bool operator ==( const FSAEdge& e ) const {
    return target==e.target && label==e.label;
  }

  bool operator < ( const FSAEdge& e ) const {
    if( label<e.label )
      return true;
    if( label>e.label )
      return false;

    if( target<e.target )
      return true;
    return false;
  }

  int target;
  int label;
};


//---------------------------------------------------------------------------//
//------------------------------- FSAState ----------------------------------//
//---------------------------------------------------------------------------//


struct FSAState 
{
  typedef FSAEdge edge_type;
  
  FSAState( int s ) : theState(s) { }
  FSAState( ) : theState(0) { }
  // doomy constructor for class map
  
  bool operator== ( const FSAState& s ) const { 
    return theState==s.theState; 
  }
  bool operator< ( const FSAState& s ) const { 
    return theState<s.theState; 
  }

  set< edge_type > in;
  set< edge_type > out;
  int theState;
};


//---------------------------------------------------------------------------//
//------------------------------- FSAState ----------------------------------//
//---------------------------------------------------------------------------//


struct FoldDetails
{
	FoldDetails( int b , int l , int n , const FSAState& s1 , const FSAState& s2 ) :
		base(b), label(l), state_num( n ), state1(s1), state2(s2) { }

	int base;
	int label;
	int state_num;
	FSAState state1;
	FSAState state2;
};


//---------------------------------------------------------------------------//
//-------------------------------- FSARep -----------------------------------//
//---------------------------------------------------------------------------//


class FSARep : public RefCounter
{
  friend class FSA;

  typedef FSAState::edge_type edge_type;

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Constructors                                       //
  //                                                     //
  /////////////////////////////////////////////////////////

 private:

  FSARep( );


  /////////////////////////////////////////////////////////
  //                                                     //
  //  Accessors:                                         //
  //                                                     //
  /////////////////////////////////////////////////////////
public:

  
  FSARep* clone( ) const { return new FSARep( *this ); }

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Subgroup graph functions:                          //
  //                                                     //
  /////////////////////////////////////////////////////////
public:

  //! Fold a FSA. Outputs details of the fold, that will allow to unfold latter. Argument candidates is a list of vertices to check for fold (all others won't be checked). If candidates=NULL then all vertices will be checked.
  void fold( const set< int >* candidates = NULL ,  list< FoldDetails >* details = NULL );
  
  //! Pinch two vertices
  void pinch( int state1 , int state2 );
  
  //! Unfold 2 vertices using details of the previous fold.
  void unfold( const list< FoldDetails >& details );
  
  //! Function allows to lift a path in the folded FSA up to unfolded one.
  void liftup( const list< FoldDetails >& details , list< FSAEdge >& path , int init_state );
  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Adding and removing elements                       //
  //                                                     //
  /////////////////////////////////////////////////////////

public:
  //! Add a new state to a FSA.
  int    newState( );
  //! Erase a state from a FSA.
  void eraseState( int state );
  //! Add a new edge to a FSA.
  void   newEdge( int state1 , int state2 , int label );
  //! Erase an edge from a FSA.
  void eraseEdge( int state1 , int state2 , int label );
  
  //! Add a loop to a FSA at a vertex vert labelled with [F,L).
  template< class ConstIntIterator > void addLoop( int vert  , ConstIntIterator F , ConstIntIterator L ) {
    if( F==L ) return;
    int cur_vert = vert;
    ConstIntIterator L2 = L;
    --L2;
    for( ConstIntIterator IT=F ; IT!=L2 ; ++IT ) {
      int new_vert = newState( );
      newEdge( cur_vert , new_vert ,  *IT );
      newEdge( new_vert , cur_vert , -*IT );
      cur_vert = new_vert;
    }
    
    newEdge( cur_vert , vert ,  *L2 );
    newEdge( vert , cur_vert , -*L2 );
  }

  //! Add a ray to a FSA at a vertex vert labelled with [F,L).
  template< class ConstIntIterator > void addRay ( int vert , ConstIntIterator F , ConstIntIterator L ) {
    int cur_vert = vert;
    for( ConstIntIterator IT=F ; IT!=L ; ++IT ) {
      int new_vert = newState( );
      newEdge( cur_vert , new_vert ,  *IT );
      newEdge( new_vert , cur_vert , -*IT );
      cur_vert = new_vert;
    }
  }
  
  //! not implemented yet
  void addFSA ( int vert1 , int vert2 , const FSARep& fsa ) { exit( 1 ); }

  const map< int , FSAState >& getStates( ) const { return theStates; }
        map< int , FSAState >& getStates( )       { return theStates; }

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Initial and terminal states:                       //
  //                                                     //
  /////////////////////////////////////////////////////////
public:
  void makeInitial    ( int s ) { initStates.insert( s ); }
  void makeTerminal   ( int s ) { termStates.insert( s ); }
  void makeNonInitial ( int s ) { initStates.erase ( s ); }
  void makeNonTerminal( int s ) { termStates.erase ( s ); }
  const set< int >& getInitStates( ) const { return initStates; }
  const set< int >& getTermStates( ) const { return termStates; }

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Data members                                       //
  //                                                     //
  /////////////////////////////////////////////////////////
  
private:

  // friend class AdvCPAlgorithm;
  
  map< int , FSAState > theStates;
  // here argument (int) coincide with FSAState.theState

  int maxState;
  
  set< int > initStates;
  // the set of initial states
  set< int > termStates;
  // the set of terminal states
};

//! to be checked!!!
void reducePath( list< FSAEdge >& path , int init_state );


#endif

