// Copyright (C) 2006 Alexander Ushakov
// Contents: Definition of class TheGrigorchukGroupAlgorithms
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//


#ifndef _TheGrigorchukGroupAlgorithms_h_
#define _TheGrigorchukGroupAlgorithms_h_



#include "map"
#include "list"
#include "set"
using namespace std;

#include "tuples.h"
#include "Word.h"


//---------------------------------------------------------------------------//
//----------------------- TheGrigorchukGroupAlgorithms ----------------------//
//---------------------------------------------------------------------------//


//! Static class TheGrigorchukGroupAlgorithms encapsulates algorithms for the original Grigorchuk group
/*!
  The class TheGrigorchukGroupAlgorithms is static, i.e., all member functions are static and there is
  no constructor defined. For introduction to the original Grigorchuk group see P. de la Harpe
  "Topics in Geometric Group Theory". For more on algorithmic properties see survey by R. Grigorchuk
  "Solved and Unsolved Problems Around One Group".
 */
class TheGrigorchukGroupAlgorithms
{
  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Constructors:                                      //
  //                                                     //
  /////////////////////////////////////////////////////////

public:
  
  //! Default constructor is not instantiated to protect from creating the obects of this class
  TheGrigorchukGroupAlgorithms( );
  
  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Accessors:                                         //
  //                                                     //
  /////////////////////////////////////////////////////////

public:


  //! Solve the Identity Problem for a non-reduced word
  /*!
    Check if a word over the original alphabet represents the identity.
    See VIII.E(47) of de la Harpe.
  */
  static bool trivial( const Word& w );


  //! Solve the Identity Problem for a reduced word
  /*!
    Check if a word over the original alphabet represents the identity.
    See VIII.E(47) of de la Harpe.
  */
  static bool trivial_reduced( const Word& w );


  //! Compute the abelian image of an element. 
  /*
    Result belongs to \f$\mathbb{Z}_2 \times \mathbb{Z}_2 \times \mathbb{Z}_2\f$.
    \f$(1,0,0)\f$ corresponds to \f$a\f$, \f$(0,1,0)\f$ corresponds to \f$b\f$, 
    \f$(0,0,1)\f$ corresponds to \f$c\f$.
  */
  static triple< int , int , int > abelianImage( const Word& w );


  //! Reduce a word
  /*!
    See VIII.B(12) of de la Harpe.
   */
  static Word reduce( const Word& w );
  
  
  //! Find the order of an element
  /*!
    Compute the smallest positive number \f$n\f$ such that \f$w^n\f$ is trivial in the Grigorchuk group.
   */
  static int findOrder( const Word& w );

  
  //! Determine if 2 words represent conjugate elements of the Grigorchuk group.
  /*
    See Section 5 of "Solved and Unsolved Problems Around One Group".
  */
  static bool conjugate( const Word& w1 , const Word& w2 ) {
    return conjugate_Kcosets( w1 , w2 ).size()!=0;
  }


  //! Compute the \f$Q\f$-set for a pair of words.
  /*
    For a pair of words \f$w_1\f$, \f$w_2\f$, \f$Q(w_1,w_2) = \{ xK \mid w_1^x = w_2\}\f$, 
    where \f$K = ncl(abab)\f$.
    See Section 5 of "Solved and Unsolved Problems Around One Group".
  */
  static set< int > conjugate_Kcosets( const Word& w1 , const Word& w2 );
  
  
  //! Determine if 2 words represent conjugate elements of the Grigorchuk group and find the actual conjugators, one for each K-coset.
  /*
    See Section 5 of "Solved and Unsolved Problems Around One Group".
  */
  static set< Word > findConjugator_Kcosets( const Word& w1 , const Word& w2 );
  

  //! Split an element of \f$St(1)\f$.
  /*!
    Make sure \f$w \in St(1)\f$! The output is a pair of words \f$(w_0,w_1)\f$, \f$w_0\f$ acts on the left subtree,
    and \f$w_1\f$ acts on the right subtree.
  */
  static pair< Word , Word > split( const Word& w );


  /////////////////////////////////////////////////////////
  //                                                     //
  //  Element routines:                                  //
  //                                                     //
  /////////////////////////////////////////////////////////

public:

  //! Multiply an element by a generator on the right and perform a reduction if possible
  static void push_back( Word& w , int g );


  //! Multiply an element by a generator on the left and perform a reduction if possible
  static void push_front( Word& w , int g );
  
  
  //! Generate a random reduced word.
  static Word randomWord( int len );

  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Certain subgroups:                                 //
  //                                                     //
  /////////////////////////////////////////////////////////

public:

  //! Compute the decomposition of a word as an element of a subgroup \f$ncl_G(b)\f$.
  /*!
    The function returns a pair \f$(c,L)\f$ where \f$c\f$ is a right coset for \f$w\f$
    and \f$L\f$ is a sequence of conjugators for \f$b\f$ comprising the decomposition for
    \f$w\f$. The equality \f$w = (\prod_{i=1}^n L_i^{-1} b L_i ) c\f$ holds.
  */
  static pair< Word , list< Word > > decompositionBSbgp( const Word& w );
    


  //! Compute the preimage of \f$(w_1,w_2)\f$ along \f$\psi : St_{\Gamma}(1) \rightarrow \Gamma \times \Gamma\f$.
  /*!
    \f$\Gamma \times \Gamma \ne \psi(St_{\Gamma}(1))\f$. The function returns a pair \f$(L,C)\f$,
    where \f$C\f$ is such that \f$(w_1 C^{-1},w_2) \in \psi(St_{\Gamma}(1))\f$ 
    and \f$L\f$ is a preimage of \f$(w_1 C^{-1},w_2)\f$.
  */
  static pair< Word , Word > liftToSTone( const Word& w1 , const Word& w2 );
  
  
  //! Find a number of a right K-coset \f$Kw\f$ where \f$K = ncl_\Gamma(abab)\f$.
  /*
    The index of \f$K\f$ in \f$\Gamma\f$ is \f$8\f$.
    Cosets are numbered by 0,...,7 - ground level, 8,...,15 above level (involving b).
    On the level elements are numbered from the trivial element, d-element next, and so on.
  */
  static int cosetRepresentativeKSbgp( const Word& w );


  //! Get a word representative of a right K-coset with number c
  /*
    The index of \f$K\f$ in \f$\Gamma\f$ is \f$8\f$.
    Cosets are numbered by 0,...,7 - ground level, 8,...,15 above level (involving b).
    On the level elements are numbered from the trivial element, d-element next, and so on.
  */
  static Word cosetRepresentativeKSbgp( int c );
  
  
  //! Lift the pair of right K-cosets \f$(Kx,Ky)\f$ to \f$Kz\f$ if possible.
  /*!
    If \f$z = (x,y)\f$ then knowing cosets \f$xK\f$ and \f$yK\f$ it is possible
    to find a coset \f$xK\f$. This operation is defined not for all pairs \f$x,y\f$.
   */
  static int liftPairKcosetsUP( int x , int y );
  
  
  //! For two sets of K-cosets find all K-lifts.
  /*!
    Function uses liftPairKcosetsUP() for each pair of cosets.
   */
  static set< int > liftPairsKcosetsUP( const set< int >& K1 , const set< int >& K2 );


  /////////////////////////////////////////////////////////
  //                                                     //
  //  Internal functions:                                //
  //                                                     //
  /////////////////////////////////////////////////////////

 private:
 public:

  //! For each pair of K-cosets that can be lifted up to \f$ST_\Gamma(1)\f$ it associates such a lift.
  static map< pair< int , int > , int > KCosetLiftTable( );
  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Data:                                              //
  //                                                     //
  /////////////////////////////////////////////////////////

private:

};


#endif
