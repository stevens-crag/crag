// Copyright (C) 2007 Alexander Ushakov
// Contents: Definition of class StraightLineProgramWord
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//


#ifndef _StraightLineProgramWord_H_
#define _StraightLineProgramWord_H_
#include <gmpxx.h>

#include "map"
#include "list"
using namespace std;


typedef mpz_class LongInteger;
#include "Map.h"


//---------------------------------------------------------------------------//
//------------------------- StraightLineProgramWord -------------------------//
//---------------------------------------------------------------------------//


//! class StraightLineProgramWord.
/*!
  A straight line program is basically a program computing a single word.
  We used "Polynomial-time Word problems" by Saul Schleimer when implementing
  this class.
  
  In this implementation the class StraightLineProgramWord serves 1 purpose -
  it is a special container class for words. In other words it is just
  another representation for words, used mainly to efficiently solve
  the Word problem for the automorphism group of a free group.

  Notice that we do not reduce SLP when constructing, i.e., initially
  a word represented by SLP is not reduced and length() is not the same
  as the length of getWord(). You need to use reduce() for that.

  Another thing to remember... I do not assume that for all rules \f$x_i \rightarrow x_j x_k\f$ 
  \f$i > j,k\f$. Though it is assumed that \f$x_j \ne 1\f$. In a situation when \f$x_j=1,~ x_k=1\f$
  \f$x_i\f$ must be removed from the system.

  After all is done need to check that all those assumptions are really satisfied.
*/


class StraightLineProgramWord
{
private:

  //! A production for a rule of a composition system.
  /*! 
    Each rule of a composition system is of the form 
    \f$X_i \rightarrow X_j^{\pm 1} \cdot X_k^{\pm 1}\f$ or \f$X_i \rightarrow X_j^{\pm 1}\f$. 
    If \f$j\f$ or \f$k\f$ is zero then that part does not exist.
  */
  struct Production {
    
    Production( int t1=0, int t2=0 , LongInteger l=0 , int h=0 ) :
      theTerm1(t1), 
      theTerm2(t2),
      theLength(l),
      theHeight(h),
      reduced(false)
    { }
    
    int theTerm1;
    int theTerm2;
    LongInteger theLength;
    int theHeight;
    bool reduced;
  };


  //! Structure used in function equal() only
  /*!
    This structure is used for comparison of 2 straight line program words (as elements of a semigroup).
    theVertex1 is a number of a vertex in the first graph, theVertex2 is a number of a vertex in the second graph,
    theBase1 is true if theVertex1 is a base vertex, theLength is the distance between the vertices.
   */
  struct Assertion {
    
    Assertion( bool b , int v1 , int v2 , LongInteger l ) : theBase1(b), theVertex1(v1), theVertex2(v2), theLength(l) { }
    
    bool operator < ( const Assertion& A ) const {
      if( theBase1<A.theBase1)
	return true;
      if( theBase1>A.theBase1)
	return false;
      if( theVertex1<A.theVertex1 )
	return true;
      if( theVertex1>A.theVertex1 )
	return false;
      if( theVertex2<A.theVertex2 )
	return true;
      if( theVertex2>A.theVertex2 )
	return false;
      if( theLength<A.theLength )
	return true;
      return false;
    }
    
    bool similar( const Assertion& A ) const {
      return (theBase1==A.theBase1 && theVertex1==A.theVertex1 && theVertex2==A.theVertex2);
    }
    
    friend ostream& operator << ( ostream& os , const Assertion& A ) {
      os << "(" << A.theBase1 << "," << A.theVertex1 << "," << A.theVertex2 << "," << A.theLength << ")";
      return os;
    }
    
    //! true if and only if the vertex 1 comes from the first SLP
    bool theBase1;
    
    //! the number of the vertex 1
    int theVertex1;
    
    //! the number of the vertex 2
    int theVertex2;
    
    //! shift in theVertex 1
    LongInteger theLength;
  };
  
  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Constructors                                       //
  //                                                     //
  /////////////////////////////////////////////////////////
  
public:

  //! Default constructor
  StraightLineProgramWord( ) : theRoot(0) , theTerminals(0) { }


  //! Create a composition system producing a 1-generator word 'init'.
  StraightLineProgramWord( int init , int gens );
  
  
  //! Construct a composition system for a generator init acted on by a sequence of mappings.
  /*!
    Here if \f$M_1,\ldots,M_n\f$ is a sequence of mappings then
    the action of a sequence is a sequence of actions where \f$M_1\f$ is
    the first mapping to apply, and so on up to \f$M_n\f$. 
  */
  template< class ConstMapIterator > 
  StraightLineProgramWord( int init , ConstMapIterator B , ConstMapIterator E ) {
    int ai = abs( init );
    theRules[ai+1] = Production( init , 0 , 1 , 1 );
    theRoot = ai+1;
    theTerminals = ai;
    for( ; B!=E ; ++B )
      *this = applyMapping( *B );
  }
  

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Operators                                          //
  //                                                     //
  /////////////////////////////////////////////////////////

public:
  
  //! Assignement operator.
  StraightLineProgramWord& operator= ( const StraightLineProgramWord& SLP ) {
    theRoot      = SLP.theRoot;
    theRules     = SLP.theRules;
    theTerminals = SLP.theTerminals;
  }

  
  //! Index operator. Obtain a generator at the ith position in the word.
  int operator[] ( LongInteger i ) const {
    return getGenerator( theRoot , i );
  }
  
  
  //! Compute a composition system which produces the inverse of the current word
  StraightLineProgramWord operator- ( ) const;
  
  
  //! Compute an SLP which produces a product of words produced by *this and SLP (as semigroup elements)
  StraightLineProgramWord operator * ( const StraightLineProgramWord& SLP ) const;

  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Accessors                                          //
  //                                                     //
  /////////////////////////////////////////////////////////
	
public:

  //! Check if words produced by \f$n_1\f$ in *this and \f$n_2\f$ in SLP are the same (as elements of semigroup)
  bool equal( int n1 , const StraightLineProgramWord& SLP , int n2 ) const;


  //! Add a vertex \f$n'\f$ into SLP (if required) which represents the word \f$w\f$ which the terminal segment of the word represented by \f$n\f$ of required length.
  int terminalSegment( int n , LongInteger length ) {
    return truncateVertexLeft( n , length_rule(n)-length );
  }


  //! Add a vertex \f$n'\f$ into SLP (if required) which represents the word \f$w\f$ which the initial segment of the word represented by \f$n\f$ of required length.
  int initialSegment( int n , LongInteger length ) {
    return truncateVertexRight( n , length_rule(n)-length );
  }


  //! Add a vertex \f$n'\f$ into SLP (if required) which represents the word \f$w_n\f$ truncated on the left.
  /*!
    If \f$w_n = u \circ v\f$, \f$length = |u|\f$, and \f$t\f$ is the output of this function
    then  \f$w_t = v\f$. The vertex \f$n\f$ can be negative. Foolproof.
   */
  int truncateVertexLeft( int n , LongInteger length );
  
  
  //! Add a vertex \f$n'\f$ into SLP (if required) which represents the word \f$w_n\f$ truncated on the left.
  /*!
    If \f$w_n = u \circ v\f$, \f$length = |v|\f$, and \f$t\f$ is the output of this function
    then  \f$w_t = u\f$. The vertex \f$n\f$ can be negative. Foolproof.
   */
  int truncateVertexRight( int n , LongInteger length );
  
  
  //! Compute the length of the longest common initial segment of 2 words given by SLPs (as elements of semigroup)
  LongInteger leftGCDLength( const StraightLineProgramWord& CS ) const {
    return leftGCDLength( theRoot , CS , CS.theRoot );
  }
  
  
  //! Get a word represented (produced) by the SLP.
  /*!
    Be careful!!! Relatively small SLPs can represent huge words.
  */
  Word getWord( ) const { return getWord( theRoot ); }
  
  
  //! Get a word represented (produced) by (non-)terminal symbol.
  /*!
    The current implementation is the most straightforward approach (not efficient).
    Notice that \f$N\f$ can be negative.
   */
  Word getWord( int N ) const;
  
  
  //! Compute the length of the word produced by the SLP.
  LongInteger length( ) const { return length_rule( theRoot ); }
  
  
  //! Reduce the SLP.
  void reduce( ) { if( theRoot!=0 ) reduce_rule( theRoot ); }
  
  
  //! Simplify the structure of SLP
  /*!
    Update vertices of out_degree 1 and remove all vertices that can not be reached
    from the root.
   */
  void simplify( );
  
  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Internal functions                                 //
  //                                                     //
  /////////////////////////////////////////////////////////

private:

  //! Order the descendants of \f$n\f$.
  /*!
    Operation returns 2 structures: order - contains the ordered set of descendants
    and closure which contains the same vertices ordered as integers.
    It is assumed that order and closure are initially empty.
  */
  void order_vertices( int n , list< int >& order , set< int >& closure ) const;
  void order_vertices( int n , list< int >& order ) const {
    set< int > closure;
    order_vertices( n , order , closure );
  }

  
  //! Compute the length of the word corresponding to a rule.
  LongInteger length_rule( int n ) const { 
    if( n==0 )
      return 0;
    if( abs( n )<=theTerminals )
      return 1;
    return (*theRules.find(abs(n))).second.theLength;
  }
  
  
  //! Compute the height of the rule.
  int height_rule( int n ) const { 
    if( n==0 )
      return 0;
    if( abs( n )<=theTerminals )
      return 0;
    return (*theRules.find(abs(n))).second.theHeight;
  }


  //! Function inverts a production pair (A,B)
  static void invertProductionPair( int& A , int& B );


  //! Determine the type of the assertion
  /*!
    The function determines the type of the assertion.
    The output==true  if A is an overlap assertion.
    The output==false if A is an subword assertion.
  */
  bool assertionType( const Assertion& A , const StraightLineProgramWord& SLP ) const;

  
  //! Extend the set of non-terminals. The word produced by the result is the same as the original.
  void extendTerminals( int nt );


  //! Reduce a subgraph of the SLP reachable from the vertex n.
  void reduce_rule( int n );
  

  //! Compute the length of the longest common initial segment of 2 words given by composition systems
  LongInteger leftGCDLength( int n1 , const StraightLineProgramWord& CS , int n2 ) const;

  
  // Split an assertion.
  set< Assertion > splitAssertion( const Assertion& A , bool firstTerm , const StraightLineProgramWord& SLP ) const;
  
  
  //! Update the lengths of all rules in the system.
  void update_lengths( );
  //! Recursively update the length of the rule.
  void update_length( int n );
  
  //! Update the heights of all rules in the system.
  void update_heights( );
  //! Recursively update the height of the rule.
  void update_height( int n );

  
  //! Get a generator at ith position of the word \f$w_n\f$.
  int getGenerator( int n , LongInteger pos ) const;
  
  
  //! Find the amount of cancellation in a product \f$w_{n_1} w_{n_2}\f$.
  LongInteger cancellationLength( int n1 , int n2 ) const;
  
  
  //! Construct the composition system rules for the map M.
  /*!
    Function prepares a block of rules for the given map.
    The range terminals have the lowest numbers \f$1,\ldots,t\f$.
    The domain non-terminals are given in the vector.
   */
  static pair< map< int , Production > , vector< int > > map_rules( const Map& M );


  //! For a composition system \f$C\f$ and a mapping \f$\alpha\f$ compute a composition system representing \f$w_C^\alpha\f$
  /*!
    The result is computed in a very straightforward way. 
    There is no check that \f$\alpha\f$ is applicable. 
    In case when the domain alphabet of \f$\alpha\f$ is strictly smaller than 
    the number of terminals the output is an empty composition sequence.
    In case when the domain alphabet of \f$\alpha\f$ is strictly greater than 
    the number of terminals the result is a correct composition sequence.
   */
  StraightLineProgramWord applyMapping( const Map& M ) const;
  
  
  //! Get the maximal number of the rule in SLP
  int max_rule_number( ) const;
  
  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  I/O                                                //
  //                                                     //
  /////////////////////////////////////////////////////////

public:  

  friend ostream& operator << ( ostream& os , const StraightLineProgramWord& CS );

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Data members                                       //
  //                                                     //
  /////////////////////////////////////////////////////////

  private:

  //! The number of the root non-terminal. Theroot==0 if the system is empty (default constructor).
  /*!
    The root is always non-negative!!!
   */
  int theRoot;
  
  
  //! The number of terminals in the system.
  int theTerminals;
  
  
  //! The production rules of the system.
  map< int , Production > theRules;
  
};


#endif
