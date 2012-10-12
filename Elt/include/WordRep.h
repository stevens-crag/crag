// Copyright (C) 2005 Alexander Ushakov
// Contents: Definition of class WordRep
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//

#ifndef _WordRep_H_
#define _WordRep_H_



#include "RefCounter.h"

#include "vector"
#include "list"
#include "iostream"
using namespace std;

class WordIterator;

typedef int Generator;

//---------------------------------------------------------------------------//
//-------------------------------- WordRep ----------------------------------//
//---------------------------------------------------------------------------//


class WordRep : public RefCounter
{
  friend class Word;
  typedef pair< int , int > PII;

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Constructors                                       //
  //                                                     //
  /////////////////////////////////////////////////////////

private:

  WordRep( );
  WordRep( const WordRep& wr );
  WordRep( const list< int >& gens );
  WordRep( const vector< int >& gens );
	// special constructors, here we assume that -x is the inverse of x

  template< class IntIterator > WordRep( const IntIterator& B , const IntIterator& E ) {
    for( IntIterator C=B ; C!=E ; ++C )
      push_back( *C );
  }


  WordRep( int g );

  WordRep operator=( const WordRep wr );

public:


  /////////////////////////////////////////////////////////
  //                                                     //
  //  Operators                                          //
  //                                                     //
  /////////////////////////////////////////////////////////
private:

  WordRep& operator *= ( const WordRep& w );
  bool     operator <  ( const WordRep& wr ) const;
  bool     operator >  ( const WordRep& wr ) const;
  bool     operator == ( const WordRep& wr ) const;
  
  
  //! Conjugate a word by another word
  WordRep& operator ^= ( const WordRep& conjugator );
  WordRep  operator ^  ( const WordRep& conjugator ) const {
    WordRep result( *this );
    result ^= conjugator;
    return result;
  }


  //! Raise a word into a power
  WordRep& operator ^= ( int power );
  WordRep  operator ^  ( int power ) const {
    WordRep result( *this );
    result ^= power;
    return result;
  }
  
  
  inline WordRep  operator * ( const WordRep& w ) const {
    WordRep result( *this );
    result *= w;
    return result;
  }
  
  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Accessors:                                         //
  //                                                     //
  /////////////////////////////////////////////////////////

public:

  WordRep* clone( ) const { return new WordRep( *this ); }

  const list< int >& getList( ) const { return theElements; }
  list< int >& getList( ) { return theElements; }

	
private:

  bool doesContain( int gen ) const;

  inline int length( ) const { return theElements.size( ); }

  int exponentSum( int gen ) const;
  int isIn( int gen ) const;

  int getPower( WordRep& base ) const;


  /////////////////////////////////////////////////////////
  //                                                     //
  //  Manipulators:                                      //
  //                                                     //
  /////////////////////////////////////////////////////////
	
private:

  WordRep inverse( ) const;
  
  
  //! Make a word trivial
  void clear( ) { theElements.clear( ); }
  
  
  void freelyReduce( WordIterator beg , WordIterator end );
  void cyclicallyReduce( );
  void cyclicallyReduce( WordRep& conjugator );

  void cyclicLeftShift( );
  // the leftmost symbol goes to the rightmost position
  void cyclicRightShift( );
  // the rightmost symbol goes to the leftmost position

  void cyclicallyPermute( int n );
  // n>0 => left-shift permute
  // n<0 => rigth-shift permute

  void segment( int from , int to );
  // cut [from, to) part of the word
  // whole word is [0,length) part of itself
  void  initialSegment( int to   );
  void terminalSegment( int from );

  template< class ConstIntIterator > 
    void insert( int pos , ConstIntIterator B , ConstIntIterator E );
  template< class ConstIntIterator > 
    void insert( WordIterator it , ConstIntIterator B , ConstIntIterator E );

  void insert( int pos , int g );
  void insert( WordIterator it , int g );
  
  void replace( WordIterator it , const Generator& g );
  void replace( int pos , const Generator& g );

  template< class ConstIntIterator > 
    void replace( WordIterator it , ConstIntIterator B , ConstIntIterator E );

  /////////////////////////////////////////////////////////
  //                                                     //
  //  I/O:                                               //
  //                                                     //
  /////////////////////////////////////////////////////////

public:

  // friend ostream& operator << ( ostream& os , const WordRep& wr ) {
  ostream& printOn( ostream& os ) const;


  /////////////////////////////////////////////////////////
  //                                                     //
  //  Internal functions:                                //
  //                                                     //
  /////////////////////////////////////////////////////////

private:

  void push_back( int g );
  void push_front( int g );
  void pop_back( ) { theElements.pop_back( ); }
  void pop_front( ) { theElements.pop_front( ); }

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Data members                                       //
  //                                                     //
  /////////////////////////////////////////////////////////

private:

  list< int > theElements;
  // list of generators, negative integers represent inverses of positive integers
};


#endif
