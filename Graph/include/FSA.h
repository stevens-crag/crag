// Copyright (C) 2005 Alexander Ushakov
// Contents: Definition of class FSA
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//

#ifndef _FSA_h_
#define _FSA_h_

#include "FSARep.h"
#include "ObjectOf.h"


//---------------------------------------------------------------------------//
//--------------------------------- FSA -------------------------------------//
//---------------------------------------------------------------------------//


class FSA : public ObjectOf< FSARep >
{
public:

  typedef FSAState state_type;
  typedef FSAState::edge_type edge_type;


  /////////////////////////////////////////////////////////
  //                                                     //
  //  Constructors                                       //
  //                                                     //
  /////////////////////////////////////////////////////////
public:

  FSA( ) : ObjectOf< FSARep >( new FSARep( ) ) { }
  // copy constructor supplied by compiler
  // destructor supplied by compiler

private:
  
  FSA( const FSARep& rep ) : ObjectOf< FSARep >( new FSARep( rep ) ) { }



  /////////////////////////////////////////////////////////
  //                                                     //
  //  Operators:                                         //
  //                                                     //
  /////////////////////////////////////////////////////////
public:
  
  FSA operator * ( const FSA& F ) const;
  bool operator== ( const FSA& F ) const; // "complete" equality

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Subgroup graph functions:                          //
  //                                                     //
  /////////////////////////////////////////////////////////
public:


  void fold( const set< int >* candidates = NULL , list< FoldDetails >* details = NULL ) { 
    change( ) -> fold( candidates , details ); 
  }
  void pinch( int state1 , int state2 ) { change( ) -> pinch( state1 , state2 ); }
  
  void unfold( const list< FoldDetails >& details ) { change( ) -> unfold( details ); }
  void liftup( const list< FoldDetails >& details , list< FSAEdge >& path , int init_state ) { change( ) -> liftup( details , path , init_state ); }

  bool isDeterministic( ) const;
  FSA    deterministic( ) const;

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Adding and removing elements                       //
  //                                                     //
  /////////////////////////////////////////////////////////

public:
  int    newState( ) { return change( ) -> newState( ); }
  void eraseState( int state ) { change( ) -> eraseState( state ); }
  void   newEdge( int state1 , int state2 , int label ) { change( ) -> newEdge( state1 , state2 , label ); }
  void eraseEdge( int state1 , int state2 , int label ) { change( ) -> eraseEdge( state1 , state2 , label ); }

  template< class ConstIntIterator > void addLoop( int vert , ConstIntIterator F , ConstIntIterator L ) 
    { change( )->addLoop( vert , F , L ); }
  template< class ConstIntIterator > void addRay ( int vert , ConstIntIterator F , ConstIntIterator L ) 
    { change( )->addRay ( vert , F , L ); }
  void addFSA ( int vert1 , int vert2 , const FSA& fsa ) { change( )->addFSA( vert1 , vert2 , *fsa.look( ) ); }

  const map< int , FSAState >& getStates( ) const { return look()->getStates( ); }
        map< int , FSAState >& getStates( )       { return change()->getStates( ); }

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Initial and terminal states:                       //
  //                                                     //
  /////////////////////////////////////////////////////////
public:


  void makeInitial    ( int s ) { change( ) -> makeInitial( s ); }
  void makeTerminal   ( int s ) { change( ) -> makeTerminal( s ); }
  void makeNonInitial ( int s ) { change( ) -> makeNonInitial( s ); }
  void makeNonTerminal( int s ) { change( ) -> makeNonTerminal( s ); }
  const set< int >& getInitStates( ) const { return look( ) -> getInitStates( ); }
  const set< int >& getTermStates( ) const { return look( ) -> getTermStates( ); }


};


ostream& operator << ( ostream& os , const FSA& g );


#endif
