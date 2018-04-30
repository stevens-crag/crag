// Copyright (C) 2005 Alexander Ushakov
// Contents: Definition of class ThRightNormalForm
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//

#ifndef _ThRightNormalForm_h_
#define _ThRightNormalForm_h_

#include <vector>
#include <list>
#include <set>

#include "tuples.h"
#include "Permutation.h"
#include "braid_group.h"

using namespace std;

class ThLeftNormalForm;


//---------------------------------------------------------------------------//
//--------------------------- ThRightNormalForm -----------------------------//
//---------------------------------------------------------------------------//

//! Defines a representation of a right Garside normal form
/*!
  Right Garside normal form has the following presentation \f$ \xi_1 \ldots \xi_{t-1} \xi_t \Delta^\omega\f$
  where \f$\Delta\f$ is a half-twist permutation and \f$\xi_1 \ldots \xi_{t-1} \xi_t\f$ is a right-weighted 
  sequence of permutation braids. A sequence of permutation braids is called right-weighted if for any pair 
  \f$\xi_i\xi_{i+1}\f$ \f$F(\xi_i) \subseteq S(\xi_{i+1})\f$. 
  Here \f$S(\xi) = \{i \mid \xi = x_i \xi' \mbox{ for some permutation } \xi' \}\f$ and
  \f$F(\xi) = \{i \mid \xi = \xi' x_i \mbox{ for some permutation } \xi' \}\f$ are starting and finishing sets.
  
  Permutations \f$\xi\f$ are translated to braids "from right to left". So, you construct a minimal positive
  braid where points on the right will go to positions defined by \f$\xi\f$ on the left.
 */

class ThRightNormalForm
{
public:
  //! Presentation of a normal form. 
  /*!
    The first  component specifies the rank  of a braid group.
    The second component specifies the power of the half twist.
    The third  component specifies the list  of braid permutations.
  */
  typedef triple< int , int , list< Permutation > > NF;
  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Constructors:                                      //
  //                                                     //
  /////////////////////////////////////////////////////////

public:
  
  //! Default constructor (rank specifies the rank of a braid group the form belongs to)
  /*!
    Default rank value is for map<...> only
    make sure to provide the argument for this constructor
  */
  ThRightNormalForm( int rank=0 ) : 
    theRank( rank ) ,
    theOmegaPower( 0 ) ,
    theDecomposition( ) { }
    
    
  //! Constructor (rank specifies the rank of a braid group the form belongs to, p - a power of a half twist, and d - a list of permutations. So the pair (p,d) is a presentation)
  /*!
    There is no check that the pair (p,d) defines a correct normal form representation.
    If you are not sure if (p,d) is correct apply static function adjustDecomposition first.
  */
  ThRightNormalForm( int rank , int p , const list< Permutation >& d ) : 
    theRank( rank ),
    theOmegaPower( p ) ,
    theDecomposition( d ) { }


  //! Constructor (tr specifies a presentation of the normal form)
  ThRightNormalForm( const NF& tr ) : 
    theRank( tr.first ),
    theOmegaPower( tr.second ) ,
    theDecomposition( tr.third ) { }
    
    
  //! Constructs the normal form of a braid word w.
  ThRightNormalForm( const crag::braidgroup::BraidGroup& G , const Word& w );
  
  
  //! Construct a positive braid from a permutation
  ThRightNormalForm( const Permutation& p ) {
    theRank = p.size( );
    Permutation omega = Permutation::getHalfTwistPermutation( theRank );
    if( p.isTrivial( ) ) {
      theOmegaPower = 0;
    } else if( p==omega ) {
      theOmegaPower = 1;
    } else {
      theOmegaPower = 0;
      theDecomposition.push_back( p );
    }
  }
  
  
  //! Assignment operator
  ThRightNormalForm& operator = ( const NF& tr ) {
    theRank = tr.first;
    theOmegaPower = tr.second;
    theDecomposition = tr.third;
    return *this;
  }
  

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Operators                                          //
  //                                                     //
  /////////////////////////////////////////////////////////
public:


  //! Cast operator (returns a representation of a normal form)
  operator NF( ) const {
    return NF( theRank , theOmegaPower , theDecomposition );
  }
  
  
  //! Compare (check if equal)
  bool operator == ( const ThRightNormalForm& rep ) const {
    return 
      theRank==rep.theRank && 
      theOmegaPower==rep.theOmegaPower && 
      theDecomposition==rep.theDecomposition;
  }
  
  
  //! Compare (check if not equal)
  bool operator != ( const ThRightNormalForm& rep ) const {
    return 
      theRank!=rep.theRank ||
      theOmegaPower!=rep.theOmegaPower || 
      theDecomposition!=rep.theDecomposition;
  }
  
  //! Compare (check if less, performed componentwise)
  bool operator < ( const ThRightNormalForm& rep ) const {
    if( theRank<rep.theRank )
      return true;
    if( theRank>rep.theRank )
      return false;
    if( theOmegaPower<rep.theOmegaPower )
      return true;
    if( theOmegaPower>rep.theOmegaPower )
      return false;
    return theDecomposition<rep.theDecomposition;
  }


  // Invert a normal form.  
  ThRightNormalForm operator - ( ) const {
    return inverse( );
  }
  
  //! Multiply two normal forms
  ThRightNormalForm operator * ( const ThRightNormalForm& rep ) const {
    return multiply( rep );
  }
  
  //! Multiply a normal form by the other on the right
  ThRightNormalForm& operator *= ( const ThRightNormalForm& rep ) {
    *this = multiply( rep );
    return *this;
  }

  //! Cast operator
  operator ThLeftNormalForm( ) const;

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Accessors:                                         //
  //                                                     //
  /////////////////////////////////////////////////////////

public:
  
  //! Get the rank of a corresponding braid group
  inline int getRank ( ) const { return theRank; }
  //! Get half-twist power.
  inline int getPower( ) const { return theOmegaPower; }
  //! Get a list of permutations.
  inline const list< Permutation >& getDecomposition( ) const { return theDecomposition; }
  

  //! Check if a normal form is trivial
  bool isTrivial( ) const { 
    return theOmegaPower==0 && theDecomposition.size( )==0;
  }


  //! Computes a word represented by a normal form.
  /*!
    For each permutation-braid computes a braid-word it represents
    and, finally, concatenate all the results.
   */
  Word getWord( ) const;


  //! Computes a "shorter" word represented by a normal form.
  /*!
    Try to use the fact that the square of a half twist generates a center.
    If the power of half-twist is negative then it cancel them with the other
    permutation braids. After that, for newly obtained permutations computes
    braid words and concatenates them.
   */
  Word getShortWord( ) const;


  //! Generate random positive normal form 
  /*!
    Function returns a right normal form with canonical length decomp_length (it can be smaller due to possible cancellations).
   */
  static ThRightNormalForm randomPositive( int rank , int decomp_length );


  //! Increase the rank of the braid.
  /*!
    The current braid does not change as an element of \f$F_\infty\f$.
    We simply increase the rank and recompute the normal form presentation.
    If N<theRank then output *this.
   */
  ThRightNormalForm increaseRank( int N ) const;
  

  //! Check whether two braids are conjugate
  /*!
    Very slow function. Checks if super summit sets of two braids contain the same element (i.e., coincide).
    ... and the super summit sets are typically huge even for decent parameters.
    @return - a pair (R,C), where R is a boolean value (true if braids are conjugate, false otherwise), and C is a conjugator (if braids conjugate)
   */
  pair< bool , ThRightNormalForm > areConjugate( const ThRightNormalForm& rep ) const;
  

  //! Adjust a normal form
  /*!
    Use this function if the triple of arguments does not satisfy "greedy conditions".
   */
  static void adjustDecomposition( int rank , int& power , list<Permutation>& decomp );

  //! Adjust a normal form
  /*!
    Use this function if there is a possiblitiy that the sequence of permutations does not 
    satisfy the "greedy conditions". Right now this might happen only after use of a 
    constructor ThRightNormalForm (const NF &tr) or, more likely, in function setDecomposition(...). 
    This function uses adjustDecomposition( ... ).
  */
  void adjust( ) {
    adjustDecomposition( theRank , theOmegaPower , theDecomposition );
  }

  //! Compute generators of the centralizer of an element
  /*! Very slow function. Function computes the graph with the super summit set 
    of the element as the vertex set. Even for a small rank of a braid group (say 10)
    and a short braid word (say 100) the super summit set is typically huge.  
   */
  set<ThRightNormalForm> computeCentralizer( ) const;


  
  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  SSS computation:                                   //
  //                                                     //
  /////////////////////////////////////////////////////////
public:
  
  //! (Low level function) Compute super summit set representative.
  /*
    Function uses cyclings and decyclings to obtain a representative of SSS.
    @return a pair (R,C), where R is the result and C is a conjugator (if A is an input then \f$R = C^{-1} A C\f$)
  */
  pair< ThRightNormalForm , ThRightNormalForm > findSSSRepresentative( ) const;
  
  
  //! Cycle a normal form.
  /*!
    Operation: \f$ \xi_1 \ldots \xi_{t-1} \xi_t \Delta^s ~\mapsto~ \xi_t \Delta^s \xi_1 \ldots \xi_{t-1} \Delta^s\f$.
    Main purpose of the function is to increase the infimum (half twist power) of a normal form.
    If \f$n\f$ cyclings do not increase infimum then the infimum is already maximal.
    @return a pair (R,C), where R is the result of a cycling and C is a conjugator (if A is an input then \f$R = C^{-1} A C\f$)
   */
  pair< ThRightNormalForm , ThRightNormalForm > cycle( ) const;
  

  //! Decycle a normal form.
  /*!
    Operation: \f$ \xi_1 \xi_2 \ldots \xi_t \Delta^s ~\mapsto~ \Delta^s \xi_2 \ldots \xi_t \Delta^s \xi_1\f$.
    Main purpose of the function is to decrease the supremum (size of the decomposition) of a normal form.
    If \f$n\f$ decyclings do not decrease supremum then the supremum is already minimal.
    @return a pair (R,C), where R is the result of a cycling and C is a conjugator (if A is an input then \f$R = C^{-1} A C\f$)
   */
  pair<ThRightNormalForm,ThRightNormalForm> decycle( ) const;
  
  
  //! Compute the minimal simple conjugator (permutation) starting from "start". 
  /*!
    Algorithm 1 from Franco, Gonzalez-Menese, "Conjugacy problem for braid and Garside groups".
    Conjugating by a simple element does not decrease the infimum of the element.
   */
  Permutation getSimpleConjugator( const Permutation& start ) const;
  
  
  //! Find a list of simple conjugators (permutations) conjugation by which does not decrease the infimum of the element.
  /*!
    The current version is very simple. For each permutation cycle \f$(i,i+1)\f$ it computes the minimal
    permutation \f$p\f$ which starts with \f$(i,i+1)\f$ such that conjugation by \f$p\f$ does not decrease the infimum 
    of the element and outputs all of such permutations. The output is not guaranteed to be minimal, 
    i.e., some of the permutations can be prefixes of the other.
  */
  set< Permutation > getSimpleConjugators( ) const;
  
  
  //! Compute a minimal simple summit conjugator for starting from "start". 
  /*!
    Algorithm 4 from Franco, Gonzalez-Menese, "Conjugacy problem for braid and Garside groups".
    Conjugating by a simple summit element does not decrease the infimum of the element and does 
    not increase the canonical length.
   */
  Permutation getSimpleSummitConjugator( const Permutation& start ) const;
  


  //! Find a list of simple conjugators (permutations) conjugation by which does not change infimum and supremum of the given element
  /*!
    The current version is very simple. For each permutation cycle \f$(i,i+1)\f$ it computes the minimal
    permutation \f$p\f$ which starts with \f$(i,i+1)\f$ such that conjugation by \f$p\f$ does not decrease the infimum 
    of the element and does not increase the canonical length,
    and outputs all of such permutations. The output is not guaranteed to be minimal, i.e., some of the permutations
    can be prefixes of the other.
  */
  set< Permutation > getSimpleSummitConjugators( ) const;




  /////////////////////////////////////////////////////////
  //                                                     //
  //  USS computation:                                   //
  //                                                     //
  /////////////////////////////////////////////////////////
public:


  //! (Low level function) Compute ultra summit set representative.
  /*!
    @return a triple (R,C,P), where R is the result, C is a conjugator (if A is an input then \f$R = C^{-1} A C\f$),
    and P is the period of the trajectory.
   */
  triple< ThRightNormalForm , ThRightNormalForm , int > findUSSRepresentative( ) const;

  
  //! Compute a simple ultra conjugator for starting from "startSymbol". 
  /*!
    Algorithm 4.9 from Gebhardt, "A New Approach to the Conjugacy Problem in Garside Groups".
    Conjugating by a simple ultra element does not decrease the infimum of the element and does 
    not increase the canonical length. Also, the result of the conjugaction belongs to the USS.
    @return a triple (R,M,P), where R is the resulting permutation, M is true if R is minimal, and P is a period of a conjugated braid.
   */
  triple< Permutation , bool , int > getSimpleUltraConjugator( int period , const Permutation& start );
  
  
  //! Compute all minimal simple ultra conjugators for "rep" 
  set< pair< Permutation , int > > getSimpleUltraConjugators( int period );
  

  //! (Aux, for USS simple) Compute a set of transports for a pair *this and (*this)^u where p is a permutation. 
  /*!
    F-set, see Definition 4.2 in Gebhardt, "A New Approach to the Conjugacy Problem in Garside Groups".
    It is important that braids (*this) and (*this)^u belong to their SSS.
    Also, it is important that (*this) is not a power of \f$\Delta\f$ (USS is "trivial" for such elements).
   */
  set< Permutation > getTransports( const Permutation& u , int period ) const;
  
  
  //! (Aux, for USS simple) Compute a transport for *this and B = (*this)^u where p is a permutation. 
  /*!
    See Proposition 2.1 and Definition 2.4 of Gebhardt "A New Approach to the Conjugacy Problem in Garside Groups".
    It is important that braids (*this) and B = (*this)^u belong to their SSS.
    Also, it is important that (*this) is not a power of \f$\Delta\f$ (USS is "trivial" for such elements).
  */
  Permutation getTransport( const ThRightNormalForm& B , const Permutation& u ) const;
  
  
  //! (Aux, for USS simple) Compute a pullback for *this and B = (*this)^s where s is a permutation. 
  /*!
    See Definition 4.6 in Gebhardt, "A New Approach to the Conjugacy Problem in Garside Groups".
    It is important that braids (*this) and (*this)^u belong to their SSS.
    Also, it is important that (*this) is not a power of \f$\Delta\f$ (USS is "trivial" for such elements).
   */
  Permutation getPullback( const Permutation& s ) const;


  //! (Aux, for USS simple) Compute the main pullback for *this and a permutation u
  /*!
   */
  Permutation getPullback( const Permutation& u , int period ) const;



  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Modifiers:                                         //
  //                                                     //
  /////////////////////////////////////////////////////////
public:

  //! Set a power of a half-twist.
  inline void setPower( int p ) { theOmegaPower = p; }
  //! Set a decomposition of a normal form (without any check of "greedy conditions"). 
  inline void setDecomposition( const list< Permutation >& d ) { theDecomposition = d; }
  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Internal functions:                                //
  //                                                     //
  /////////////////////////////////////////////////////////
private:

  
  //! (Low level function) Invert a normal form.
  NF inverse( ) const;
  
  
  //! (Low level function) Multiply normal forms.
  NF multiply( const ThRightNormalForm& rep ) const;
  
  
  //! The main function which "canonify" a product of two permutations.
  enum transformationResult { TWO_MULTIPLIERS , ONE_MULTIPLIER , NO_CHANGE };
  static transformationResult transform ( int theIndex , Permutation& p1 , Permutation& p2 );
  

  //! I think function processes a vertex in a Super Summit Set graph.
  pair< bool , ThRightNormalForm > sssConstructionIteration( map< ThRightNormalForm , ThRightNormalForm >& sss_new1 ,
							     map< ThRightNormalForm , ThRightNormalForm >& sss_checked1 ,
							     const map< ThRightNormalForm , ThRightNormalForm >& sss_new2 ,
							     const map< ThRightNormalForm , ThRightNormalForm >& sss_checked2 ) const;
  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Data members:                                      //
  //                                                     //
  /////////////////////////////////////////////////////////

private:

  //! The rank of the braid group (number of strands).
  int theRank;

  //! Power of omega.
  int theOmegaPower;

  //! Sequence of permutations.
  list< Permutation > theDecomposition;
  
};


ostream& operator << ( ostream& os, const ThRightNormalForm& nf );


#endif
