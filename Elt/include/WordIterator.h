// Copyright (C) 2005 Alexander Ushakov
// Contents: Definition of class WordIterator
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//

#ifndef _WordIterator_H_
#define _WordIterator_H_


#include "list"
using namespace std;

class Word;

//---------------------------------------------------------------------------//
//------------------------------- WordIterator ------------------------------//
//---------------------------------------------------------------------------//


class WordIterator
{
	friend class ConstWordIterator;
	friend class WordRep;
	friend class Word;

public:
	typedef pair< int , int > PII;

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Constructors                                       //
  //                                                     //
  /////////////////////////////////////////////////////////

public:
	WordIterator( );
	WordIterator( Word& w , bool begin=true );

private:
	WordIterator( list<int>& lst , list< int >::iterator it );

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Operators                                          //
  //                                                     //
  /////////////////////////////////////////////////////////
public:
	bool operator!= ( const WordIterator& WI ) const;
	bool operator== ( const WordIterator& WI ) const;

	const WordIterator& operator++ ( );
	      WordIterator  operator++ ( int doomy );
	const WordIterator& operator-- ( );
	      WordIterator  operator-- ( int doomy );
	int   operator*  ( ) const;

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Accessors:                                         //
  //                                                     //
  /////////////////////////////////////////////////////////

public:

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Data members                                       //
  //                                                     //
  /////////////////////////////////////////////////////////
private:

	list< int >* theList;
  list< int >::iterator theIterator;
};


//---------------------------------------------------------------------------//
//----------------------------- ConstWordIterator ---------------------------//
//---------------------------------------------------------------------------//


class ConstWordIterator
{
	friend class WordRep;

public:
	typedef pair< int , int > PII;

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Constructors                                       //
  //                                                     //
  /////////////////////////////////////////////////////////
public:
	ConstWordIterator( );
	ConstWordIterator( const Word& w , bool begin=true );
	ConstWordIterator( const WordIterator& WI );

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Operators                                          //
  //                                                     //
  /////////////////////////////////////////////////////////
public:
	ConstWordIterator& operator= ( const WordIterator& WI );

	bool operator!= ( const ConstWordIterator& WI ) const;
	bool operator== ( const ConstWordIterator& WI ) const;

	ConstWordIterator& operator++ ( );
	      ConstWordIterator  operator++ ( int doomy );
	ConstWordIterator& operator-- ( );
	      ConstWordIterator  operator-- ( int doomy );
	int   operator*  ( ) const;

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Accessors:                                         //
  //                                                     //
  /////////////////////////////////////////////////////////

public:


  /////////////////////////////////////////////////////////
  //                                                     //
  //  Data members                                       //
  //                                                     //
  /////////////////////////////////////////////////////////
private:

	const list< int >* theList;
  list< int >::const_iterator theIterator;
};


#endif


