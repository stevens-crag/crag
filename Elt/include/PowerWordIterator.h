// Copyright (C) 2005 Alexander Ushakov
// Contents: Definition of class PowerWordIterator
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
	PII  operator*  ( ) const;

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

	Word* theWord;
  list< PII >::iterator theIterator;
	int theOffset;

};


//---------------------------------------------------------------------------//
//----------------------------- ConstWordIterator ---------------------------//
//---------------------------------------------------------------------------//


class ConstWordIterator
{
	frient class Word;

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

	const ConstWordIterator& operator++ ( );
	      ConstWordIterator  operator++ ( int doomy );
	const ConstWordIterator& operator-- ( );
	      ConstWordIterator  operator-- ( int doomy );
	PII  operator*  ( ) const;

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

	const Word* theWord;
  list< PII >::const_iterator theIterator;
	int theOffset;

};


#endif


