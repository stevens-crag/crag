// Copyright (C) 2005 Alexander Ushakov
// Contents: Definition of class PowerWordRep
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//

#ifndef _PowerWordRep_H_
#define _PowerWordRep_H_



#include "RefCounter.h"

#include "vector"
#include "list"
#include "iostream"
using namespace std;



//---------------------------------------------------------------------------//
//------------------------------ PowerWordRep -------------------------------//
//---------------------------------------------------------------------------//


class PowerWordRep : public RefCounter
{
  friend class PowerWord;
  typedef pair< int , int > PII;

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Constructors                                       //
  //                                                     //
  /////////////////////////////////////////////////////////

private:

  PowerWordRep( );
  PowerWordRep( const list< PII >& gens );
  PowerWordRep( const vector< PII >& gens );
  PowerWordRep( const PowerWordRep& wr );

  PowerWordRep( const list< int >& gens );
  PowerWordRep( const vector< int >& gens );
  // special constructors, here we assume that -x is the inverse of x

  PowerWordRep( int g , int p = 1 );

  PowerWordRep operator=( const PowerWordRep wr );

public:


  /////////////////////////////////////////////////////////
  //                                                     //
  //  Operators                                          //
  //                                                     //
  /////////////////////////////////////////////////////////
private:

  PowerWordRep& operator *= ( const PowerWordRep& w );
  bool     operator <  ( const PowerWordRep& wr ) const;
  bool     operator >  ( const PowerWordRep& wr ) const;
  bool     operator == ( const PowerWordRep& wr ) const;

  inline PowerWordRep  operator * ( const PowerWordRep& w ) const {
    PowerWordRep result( *this );
    result *= w;
    return result;
  }
	
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Accessors:                                         //
  //                                                     //
  /////////////////////////////////////////////////////////

public:

  PowerWordRep* clone( ) const { return new PowerWordRep( *this ); }

  const list< PII >& getList( ) const { return theElements; }
  list< PII >& getList( ) { return theElements; }

	
private:

  bool doesContain( int gen ) const;

  inline int length( ) const { return theLength; }

  int exponentSum( int gen ) const;
  int isIn( int gen ) const;

  int getPower( PowerWordRep& base ) const;


  ////////////////////////////////////////////
  //                                        //
  //  Manipulators                          //
  //                                        //
  ////////////////////////////////////////////
	
private:

  PowerWordRep inverse( ) const;

  void cyclicallyReduce( );
  void cyclicallyReduce( PowerWordRep& conjugator );

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

  void insert( const PowerWordRep& wr , int pos );


  /////////////////////////////////////////////////////////
  //                                                     //
  //  I/O:                                               //
  //                                                     //
  /////////////////////////////////////////////////////////

public:

  // friend ostream& operator << ( ostream& os , const PowerWordRep& wr ) {
  ostream& printOn( ostream& os ) const {
    list< PII >::const_iterator w_it = theElements.begin( );
    for( ; w_it!=theElements.end( ) ; ++w_it ) {
			int g = (*w_it).first;
			int p = (*w_it).second;
      if( w_it!=theElements.begin( ) ) os << " ";
      os << 'x' << g;
			if( p!=1 ) os << '^' << p;
    }
    return os;
  }


  /////////////////////////////////////////////////////////
  //                                                     //
  //  Internal functions:                                //
  //                                                     //
  /////////////////////////////////////////////////////////

private:

	void pushGeneratorBack( int g , int p );
	void pushGeneratorFront( int g , int p );

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Data members                                       //
  //                                                     //
  /////////////////////////////////////////////////////////

private:

  list< PII > theElements;
	// list of pairs (gen#,power)

	int theLength;
	// the total length of the word (the total sum of power)
};


#endif
