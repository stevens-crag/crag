// Copyright (C) 2007 Alexander Ushakov
// Contents: Definition of class ThompsonGroupFNormalForm
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//


#ifndef _ThompsonGroupFNormalForm_H_
#define _ThompsonGroupFNormalForm_H_


#include "Word.h"


//---------------------------------------------------------------------------//
//------------------------- ThompsonGroupFNormalForm ------------------------//
//---------------------------------------------------------------------------//


//! Class ThompsonGroupFNormalForm (defines a representation of a normal form of an element of the Thompson's group F (infinitely generated))//
/*!
  A normal form of an element of the Thompson's group F is a word \f$w = n\cdot p\f$
  To be finished.
 */


class ThompsonGroupFNormalForm : protected Word
{
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Constructors                                       //
  //                                                     //
  /////////////////////////////////////////////////////////

 public:

  ThompsonGroupFNormalForm( );

  ThompsonGroupFNormalForm( const Word& w );

  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Operators                                          //
  //                                                     //
  /////////////////////////////////////////////////////////

public:
  
  bool operator== ( const ThompsonGroupFNormalForm& nf ) const {
    return (Word)*this==(Word)nf;
  }
  
  
  //! Multiply the word on the right by another word (the result is reduced)
  inline ThompsonGroupFNormalForm& operator *= ( const ThompsonGroupFNormalForm& nf ) { 
    list< int > unit1 = getList( );
    list< int > unit2 = nf.getList( );
    list< Generator > new_unit = mergeUnits( unit1 , unit2 );
    *this = removeBadPairs( new_unit );
    return *this;
  }
  
  
  //! Multiply two normal forms (the result is a normal form)
  inline ThompsonGroupFNormalForm operator * ( const ThompsonGroupFNormalForm& nf ) const { 
    list< int > unit1 = getList( );
    list< int > unit2 = nf.getList( );
    list< Generator > new_unit = mergeUnits( unit1 , unit2 );
    return removeBadPairs( new_unit );
  }
  
  //! Multiply two normal forms (the result is a normal form)
  inline ThompsonGroupFNormalForm operator - ( ) const { 
    return ThompsonGroupFNormalForm( -(Word)*this );
  }
  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Accessors                                          //
  //                                                     //
  /////////////////////////////////////////////////////////
	
public:
  

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Internal functions                                 //
  //                                                     //
  /////////////////////////////////////////////////////////

 private:

  //! Function returns a "sorted peak-like word" which can contain bad pairs of generators.
  static Word semiNormalFormFor( const Word& w );
  
  //! Function merges 2 seminormal forms.
  static list< int > mergeUnits( list< int >& unit1 , list< int >& unit2 );
  
  //! Auxiliary function used in "mergeUnits".
  static int nextOperation( int n1 , int n2 , int n3 , int n4 );

  //! Remove pairs of inverse generators in a seminormal form
  static Word removeBadPairs( const Word& w );

	 
  /////////////////////////////////////////////////////////
  //                                                     //
  //  I/O                                                //
  //                                                     //
  /////////////////////////////////////////////////////////

   friend ostream& operator << ( ostream& os , const ThompsonGroupFNormalForm& nf );
   
  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Data members                                       //
  //                                                     //
  /////////////////////////////////////////////////////////

  private:
		

  
};


#endif
