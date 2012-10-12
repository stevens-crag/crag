// Copyright (C) 2005 Alexander Ushakov
// Contents: Definition of class Word
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//

#ifndef _Word_H_
#define _Word_H_



#include "ObjectOf.h"
#include "Word.h"
#include "WordRep.h"
#include "WordIterator.h"
#include <string>
#include "Alphabet.h"
// #include <strstream>

#include <set>
using namespace std;


//---------------------------------------------------------------------------//
//---------------------------------- Word -----------------------------------//
//---------------------------------------------------------------------------//

//! Class Word (defines a representation of a Word over a group alphabet)//
/*!
  A reduced word is the one which does not involve subwords of the type \f$x x^{-1}\f$ or \f$x^{-1} x\f$.
  We represent a reduced word over a group alphabet \f$X\f$ as a list of non-trivial integers list< int >.
  Each generator \f$x \in X\f$ is represented by the unique number  \f$n_x\f$.
  For each \f$x \in X\f$ it is assumed that \f$n_{x^{-1}} = -n_x\f$.
*/

class Word : public ObjectOf< WordRep >
{

 public:
  typedef ConstWordIterator const_iterator;
  typedef WordIterator iterator;
  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Constructors                                       //
  //                                                     //
  /////////////////////////////////////////////////////////

 public:

  typedef pair< int , int > PII;

  //! Default constructor (creates the empty word \f$w = \varepsilon\f$).
  Word( ) : ObjectOf< WordRep >( new WordRep() ) { }
  //! Cast constructor. Constructs a word by its presentation (if the word defined in gens is not reduced then it reduces it).
  Word( const vector< int >& gens ) : ObjectOf< WordRep >( new WordRep( gens ) ) { }
  //! Cast constructor. Constructs a word by its presentation (if the word defined in gens is not reduced then it reduces it).
  Word( const list< int >& gens ) : ObjectOf< WordRep >( new WordRep( gens ) ) { }
  //! Cast constructor. Constructs a one letter word.
  Word( int g ) : ObjectOf< WordRep >( new WordRep( g ) ) { }
    
  template< class IntIterator > 
    Word( const IntIterator& B , const IntIterator& E ) : ObjectOf< WordRep >( new WordRep( B , E ) ) { }

  // copy constructor supplied by compiler
  // destructor supplied by compiler
  
 
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Operators                                          //
  //                                                     //
  /////////////////////////////////////////////////////////
public:

  //! Comparison operator
  /*! The result of comparison of \f$u = u_1 \ldots u_k\f$ and \f$v = v_1 \ldots v_m\f$ is defined by
    \f$ u<v ~\Leftrightarrow~ k<m ~\vee~ ( k=m ~\wedge~ u_i<v_i \mbox{ where $i$ is the smallest index such that $u_i\ne v_i$}) \f$
  */
  bool operator < ( const Word& wr ) const { return (*look() < *wr.look( ) ); }
  //! Comparison operator
  bool operator > ( const Word& wr ) const { return (*look() > *wr.look( ) ); }
  //! Comparison operator
  bool operator == ( const Word& wr ) const { return (*look() == *wr.look( ) ); }
  //! Comparison operator
  bool operator != ( const Word& wr ) const { return !(*look() == *wr.look( ) ); }
  
  //! Multiply the word on the right by another word (the result is reduced)
  inline Word& operator *= ( const Word& w ) { 
    *change( ) *= *w.look(); 
    return *this;
  }
  //! Multiply two words (the result is reduced)
  inline Word operator * ( const Word& w ) const { 
    Word result( *this );
    result *= w;
    return result;
  }
  
  
  //! Conjugate a word by another word (the result is reduced)
  Word& operator ^= ( const Word& conjugator ) {
    *change( ) ^= *conjugator.look(); 
    return *this;
  }
  Word  operator ^  ( const Word& conjugator ) const {
    Word result( *this );
    result ^= conjugator;
    return result;
  }


  //! Conjugate a word by another word (the result is reduced)
  Word& operator ^= ( int power ) {
    *change( ) ^= power; 
    return *this;
  }
  Word  operator ^  ( int power ) const {
    Word result( *this );
    result ^= power;
    return result;
  }
  

  //! Invert a word (works the same as inverse).
  inline Word operator - ( ) const { 
    Word result;
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
  
  //! Freely reduce a word.
  Word freelyReduce( ) { return freelyReduce( begin() , end( ) ); }
  //! Freely reduce a segment of a word defined by [B,E).
  Word freelyReduce( iterator B , iterator E );
  // freely reduces a segment of a word specified by [beg,end]. 
  // Most of the operations produce freely reduced words,
  // except operations insert, remove. The function takes linear time to performs
  // which is not efficient. So, either try to avoid uses of functions insert and remove
  // or try to bound the changes in variables beg and end
  
  //! Generates a pseudo randomly reduced word of the length wLen over the alphabet \f$X = \{x_1,\ldots,x_n\}\f$.
  static Word randomWord( int n , int wLen );
  //! Generates a pseudo randomly reduced word of a length in [wLenMin,wLenMax] and  over the alphabet \f$X = \{x_1,\ldots,x_n\}\f$.
  static Word randomWord( int n , int wLenMin, int wLenMax );

  //! Get a constant representation the word.
  inline const list< int >& getList( ) const { return look( )->getList( ); }
  //! Get a representation the word. This allows direct manipulation with the representation (which requires caution).
  inline list< int >& getList( ) { return change( )->getList( ); }
  
  //! Get the length of the word.
  inline int length( ) const { return look( )->length( ); }

  //! Multiply the word by a one-letter word defined by gen on the right. The result is being reduced.
  inline Word& push_back ( int gen ) { change()->push_back ( gen ); return *this; }
  //! Multiply the word by a one-letter word defined by gen on the left. The result is being reduced.
  inline Word& push_front( int gen ) { change()->push_front( gen ); return *this; }
  //! Multiply the word by a word on the right. The result is being reduced.
  Word& push_back ( const Word& w );
  //! Multiply the word by a word on the left. The result is being reduced.
  Word& push_front( const Word& w );

  //! Remove the last symbol
  void pop_back ( ) { change()->pop_back ( ); }
  //! Remove the first symbol
  void pop_front( ) { change()->pop_front( ); }


  //! Determines the power of the word (as an element of a free monoid, not as an element of a free group).
  int getPower( Word& base ) const { return look( )->getPower( *base.change( ) ); }
  
  //! Checks if the word contains a generator given by gen.
  bool doesContain( const int& gen ) const { return look( )->doesContain( gen ); }
  
  //! Shifts the word one position to the left, i.e., \f$cyclicLeftShift( x_1 x_2 \ldots x_{k-1} x_k ) = x_2 \ldots x_{k-1} x_k x_1 \f$.
  inline void cyclicLeftShift( )  { change( )->cyclicLeftShift(  ); }
  //! Shifts the word one position to the right, i.e., \f$cyclicRightShift( x_1 x_2 \ldots x_{k-1} x_k ) = x_k x_1 x_2 \ldots x_{k-1} \f$.
  inline void cyclicRightShift( ) { change( )->cyclicRightShift( ); }
  
  //! Returns the cyclically reduced word.
  inline Word cyclicallyReduce( ) const {
    Word result( *this );
    result.change( ) -> cyclicallyReduce( );
    return result;
  }
  
  //! Cyclically reduces the Word
  inline void cyclicallyReduceWord( ) { change( ) -> cyclicallyReduce( ); }

  //! Returns the cyclically reduced word and the corresponding conjugator.
  inline Word cyclicallyReduce( Word& conjugator ) const {
    Word result( *this );
    result.change( ) -> cyclicallyReduce( *conjugator.change( ) );
    return result;
  }
  //! Cyclically reduces the Word and returns the corresponding conjugator.
  inline void cyclicallyReduceWord( Word& conjugator ) { 
    change( ) -> cyclicallyReduce( *conjugator.change( ) ); 
  }

  //! Invert the word.
  inline Word inverse( ) const {
    Word result;
    (*result.change( )) = look( )->inverse( );
    return result;
  }

  //! Cyclically permute the word and return the result.
  /*! n>0 => left-shift, n<0 => rigth-shift permute. */
  inline Word cyclicallyPermute( int n ) const {
    Word result( *this );
    result.change( ) -> cyclicallyPermute( n );
    return result;
  }
  //! Cyclically permute the word.
  inline Word& _cyclicallyPermute( int n ) {
    change( ) -> cyclicallyPermute( n );
    return *this;
  }

  //! Get an initial segment of the word of length len.
  inline Word initialSegment( int len ) const {
    Word result( *this );
    result.change( ) -> initialSegment( len );
    return result;		
  }

  //! Get a terminal segment of the word of length len.
  inline Word terminalSegment( int len ) const {
    Word result( *this );
    result.change( ) -> terminalSegment( len );
    return result;			
  }

  // Get a segment of the word defined by ???
  inline Word segment( int from , int to ) const {
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
  

  //! Insert a sequence of generators [B,E) into a word at a position pos
  template< class ConstIntIterator > 
  void insert( int pos , ConstIntIterator B , ConstIntIterator E );

  //! Insert a generator g into a word at a position pos
  void insert( int pos , int g );

  //! Insert a sequence of generators [B,E) into a word before a position it
  template< class ConstIntIterator > 
  void insert( WordIterator it , ConstIntIterator B , ConstIntIterator E );
  //! Insert a generator g into a word before a position it
  void insert( WordIterator it , int g );

  //! Replace a generator at a position it by g
  void replace( WordIterator it , const Generator& g );

  //! Replace a generator at a position pos by g
  void replace( int pos , const Generator& g );

  //! Replace a subword of a word starting at a position it by a word [B,E). The length of the word does not increase if [B,E) is longer than the terminal segment of the word [it,end()). In that case terminal symbols of [B,E) are ignored.
  template< class ConstIntIterator > 
    void replace( WordIterator it , ConstIntIterator B , ConstIntIterator E );

  //! Returns a word in which generators are replaced by words
  /*! Replace generators with the words (images) contained 
    in the vector \param images
    */  
  Word replaceGenerators( const vector<Word>& images )const;


  //! Compute the minimal equivalent word. 
  /*! permutableGenerators must be positive.
   */
  Word minimalEquivalentForm( const set< int >& permutableGenerators , bool inverses , bool cyclicPermutations ) const;


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

  
  friend ostream& operator << ( ostream& os , const Word& w ) {
    InfiniteAlphabet::defaultAlphabet.printWord( os, w );
    return os;
    //return w.look( )->printOn( os );
  }
  
  friend istream& operator >> ( istream& is ,  Word& w ) {
    w = InfiniteAlphabet::defaultAlphabet.readWord( is );
    return is;
    //return w.look( )->printOn( os );
  }

  ostream& printOn ( ostream& os ) const {
    InfiniteAlphabet::defaultAlphabet.printWord( os, *this );
    return os;
    //return look( )->printOn( os );
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
