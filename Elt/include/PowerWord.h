// Copyright (C) 2005 Alexander Ushakov
// Contents: Definition of class PowerWord
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//

#ifndef _Word_H_
#define _Word_H_


#include "ObjectOf.h"
#include "PowerWord.h"
#include "WordRep.h"
#include "PowerWordIterator.h"
#include <string>
#include <strstream>
using namespace std;


//---------------------------------------------------------------------------//
//-------------------------------- PowerWord --------------------------------//
//---------------------------------------------------------------------------//


class PowerWord : public ObjectOf< PowerWordRep >
{
  // friend class WordIterator;
  
 public:
  typedef ConstPowerWordIterator const_iterator;
  typedef PowerWordIterator iterator;
  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Constructors                                       //
  //                                                     //
  /////////////////////////////////////////////////////////

public:

  typedef pair< int , int > PII;

  PowerWord( ) : ObjectOf< PowerWordRep >( new WordRep() ) { }
  PowerWord( const vector< int >& gens ) : ObjectOf< PowerWordRep >( new PowerWordRep( gens ) ) { }
  PowerWord( const list< int >& gens ) : ObjectOf< PowerWordRep >( new PowerWordRep( gens ) ) { }
  PowerWord( const vector< PII >& gens ) : ObjectOf< PowerWordRep >( new PowerWordRep( gens ) ) { }
  PowerWord( const list< PII >& gens ) : ObjectOf< PowerWordRep >( new PowerWordRep( gens ) ) { }
 
  PowerWord( int g , int p=1 ) : ObjectOf< PowerWordRep >( new PowerWordRep( g , p ) ) { }
  // Cast constructor.
  
  // copy constructor supplied by compiler
  // destructor supplied by compiler
  
 
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Operators                                          //
  //                                                     //
  /////////////////////////////////////////////////////////
public:

  bool operator < ( const PowerWord& wr ) const {
    return (*look() < *wr.look( ) );
  }
  bool operator > ( const PowerWord& wr ) const {
    return (*look() > *wr.look( ) );
  }
  
  bool operator == ( const PowerWord& wr ) const {
    return (*look() == *wr.look( ) );
  }
  bool operator != ( const PowerWord& wr ) const {
    return !(*look() == *wr.look( ) );
  }
  
  inline PowerWord& operator *= ( const PowerWord& w ) { 
    *change( ) *= *w.look(); 
    return *this;
  }
  
  inline PowerWord operator * ( const PowerWord& w ) const { 
    PowerWord result( *this );
    result *= w;
    return result;
  }

  inline PowerWord operator - ( ) const { 
    PowerWord result;
    (*result.change( )) = look( )->inverse( );
    return result;
  }


  /////////////////////////////////////////////////////////
  //                                                     //
  //  Accessors                                          //
  //                                                     //
  /////////////////////////////////////////////////////////
	
 public:
  
  const_iterator begin( ) const { return const_iterator(*this); }
  const_iterator end  ( ) const { return const_iterator(*this,false); }
  iterator       begin( )       { return iterator(*this); }
  iterator       end  ( )       { return iterator(*this,false); }
  
  PowerWord& freelyReduce( ) { return freelyReduce( begin() , end( ) ); }
  PowerWord& freelyReduce( const_iterator beg , const_iterator end );
  // freely reduces a segment of a word specified by [beg,end]. 
  // Most of the operations produce freely reduced words,
  // except operations insert, remove. The function takes linear time to performs
  // which is not efficient. So, either try to avoid uses of functions insert and remove
  // or try to bound the changes in variables beg and end
  
  static PowerWord randomWord( int gens , int wLen );
  // generates a pseudo random reduced word 
  
  inline const list< PII >& getList( ) const { return look( )->getList( ); }
  // get a list representing the Word
  inline list< PII >& getList( ) { return change( )->getList( ); }
  // get a list representing the Word
  
  inline int length( ) const { return look( )->length( ); }
  // the length of the Word

  inline void push_back ( int gen , int power ) { change()->pushGeneratorBack ( gen , power ); }
  inline void push_front( int gen , int power ) { change()->pushGeneratorFront( gen , power ); }
  inline void push_back ( const pair< int , int >& g ) { change()->pushGeneratorBack ( g.first , g.second ); }
  inline void push_front( const pair< int , int >& g ) { change()->pushGeneratorFront( g.first , g.second ); }

  int getPower( PowerWord& base ) const { return look( )->getPower( *base.change( ) ); }
  // determines the power of the Word 
  // (not as an element of the corresponding Free Group)

  bool doesContain( const int& gen ) const {
    return look( )->doesContain( gen );
  }
  // Checks if the Word contains a generator

  inline void cyclicLeftShift( )  { change( )->cyclicLeftShift(  ); }
  // Shifts the word one position to the left
  inline void cyclicRightShift( ) { change( )->cyclicRightShift( ); }
  // Shifts the word one position to the right

  inline PowerWord cyclicallyReduce( ) const {
    Word result( *this );
    result.change( ) -> cyclicallyReduce( );
    return result;
  }
  // Returns the cyclically reduced Word
  
  inline void cyclicallyReduceWord( ) { change( ) -> cyclicallyReduce( ); }
  // Cyclically reduces the Word

  inline PowerWord cyclicallyReduce( PowerWord& conjugator ) const {
    PowerWord result( *this );
    result.change( ) -> cyclicallyReduce( *conjugator.change( ) );
    return result;
  }
  inline void cyclicallyReduceWord( PowerWord& conjugator ) { 
    change( ) -> cyclicallyReduce( *conjugator.change( ) ); 
  }
  
  inline PowerWord inverse( ) const {
    PowerWord result;
    (*result.change( )) = look( )->inverse( );
    return result;
  }
  
  inline PowerWord cyclicallyPermute( int n ) const {
    Word result( *this );
    result.change( ) -> cyclicallyPermute( n );
    return result;
  }
  // n>0 => left-shift permute
  // n<0 => rigth-shift permute
  
  inline PowerWord initialSegment( int len ) const {
    PowerWord result( *this );
    result.change( ) -> initialSegment( len );
    return result;		
  }
  
  inline PowerWord terminalSegment( int len ) const {
    PowerWord result( *this );
    result.change( ) -> terminalSegment( len );
    return result;			
  }
  
  inline PowerWord segment( int from , int to ) const {
    Word result( *this );
    result.change( ) -> segment( from , to );
    return result;	
  }
  
  inline int exponentSum( const int& gen ) const {
    return look( )->exponentSum( gen );
  }
  
  int isIn( const int& gen ) const {
    return look( )->isIn( gen );
  }
  
  Word power( int t ) const;

  void insert( const PowerWord& wr , int pos ) {
    change( ) -> insert( *wr.look( ) , pos );
  }

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Internal functions                                 //
  //                                                     //
  /////////////////////////////////////////////////////////

 private:

	 
  /////////////////////////////////////////////////////////
  //                                                     //
  //  I/O                                                //
  //                                                     //
  /////////////////////////////////////////////////////////

  
  friend ostream& operator << ( ostream& os , const PowerWord& w ) {
    return w.look( )->printOn( os );
  }
  
  ostream& printOn ( ostream& os ) const {
    return look( )->printOn( os );
  }
  
  /*
    friend ostream& operator << ( ostream& os , const Word& w ) {
    os << *w.look( );
    return os;
    }
  */
  


  /////////////////////////////////////////////////////////
  //                                                     //
  //  Data members                                       //
  //                                                     //
  /////////////////////////////////////////////////////////

  private:
		
};


#endif
