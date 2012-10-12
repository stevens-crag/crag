// Copyright (C) 2005 Alexander Ushakov
// Contents: Definition of class Permutation
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//


#ifndef _Permutation_h_
#define _Permutation_h_

#include <set>
#include <map>
#include <vector>
#include <iostream>

using namespace std;


//---------------------------------------------------------------------------//
//----------------------------- Permutation ---------------------------------//
//---------------------------------------------------------------------------//


//! Permutation.
/*!
  We use the standard representation of a presentation on \f$n\f$ symbols - a sequence of numbers
  \f$( x_0, \ldots, x_{n-1} )\f$. Notice that indexes start at 0.
*/


class Permutation
{

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Constructors                                       //
  //                                                     //
  /////////////////////////////////////////////////////////
  
public:
  
  //! Default constructor (size specifies the size of the permutation)
  Permutation( int size=0 );
  
  //! Construct a permutation of the specified size given by the array p
  Permutation( int size , const int* p );

  //! Construct a permutation by a vector of numbers
  Permutation( const vector< int >& p );
  

  //! Construct a permutation by a sequence of numbers
  template< class IntIterator > Permutation( const IntIterator& B , const IntIterator& E ) : theValue( B , E ) { }

  
  //! Copy constructor provided by compiler
  // Permutation( const Permutation& p );
  
  
  //! What is it?
  static Permutation CYCLE( int N , const vector< int >& cycle );
  
  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Operators:                                         //
  //                                                     //
  /////////////////////////////////////////////////////////
  
public:

  
  //! Get the ith element of the permutation
  int  operator[] ( int ind ) const { return theValue[ind]; }
  
  
  //! Get the reference to the ith element of the permutation
  int& operator[] ( int ind ) { return theValue[ind]; }

  //! Check if 2 permutations are equal
  bool operator == ( const Permutation& p ) const;


  //! Check if 2 permutations are not equal
  bool operator != ( const Permutation& p ) const;
  
  
  //! Check if *this is strictly less than p (compared lexicographically)
  bool operator < ( const Permutation& p ) const;


  //! Multiple 2 permutations
  Permutation  operator *  ( const Permutation& p ) const;


  //! Multiple 1 permutation by another on the right
  Permutation& operator *= ( const Permutation& p );


  //! Invert the permutation
  Permutation operator - ( ) const { return inverse( ); }


  //! Invert the permutation
  Permutation inverse( ) const;


  //! A function used for computation of BKL normal forms of braids
  Permutation& left_mult_by_cycle( const vector< int >& cycle );


  //! A function used for computation of BKL normal forms of braids
  static void lr_multiply_by_cycles( Permutation& P , Permutation& I , const vector< int >& M1 , const vector< int >& M2 );
  // a function for fast multiplication of a permutation by cycles on the left and on the right
  // P = -I !!! a must


  /////////////////////////////////////////////////////////
  //                                                     //
  //  Accessors:                                         //
  //                                                     //
  /////////////////////////////////////////////////////////
  
public:
  
  
  //! Swap the values at uth and jth position (multiply this by a cycle (i,j) on the left)
  void change( int i , int j ) { swap( theValue[i] , theValue[j] ); }

  
  //! Get the vector of numbers which represents the permutation.
  const vector< int >& getVector( ) const { return theValue; }


  //! Raise the permutation into the power p.
  Permutation power( int p ) const;
  

  //! Get the size of the permutation
  int size( ) const { return theValue.size( ); }
  
  //! Increase the size of a permutation.
  Permutation increaseSize( int N ) const;
  
  
  //! Compute the length of a permutation (length of a geodesic). Current version is slow, need to update.
  int length( ) const;


  //! Compute the number of positions with different elements
  int difference( const Permutation& p ) const;

  
  //! Compute the conjugacy class representative (2 elements conjugate iff they have the same result of this function)
  Permutation computeConjugacyClassRepresentative( Permutation& conj ) const;


  //! Compute a conjugator for a couple of permutations (if permutations are not conjugate then the result makes no sense)
  Permutation computeConjugator( const Permutation& p ) const;
  // bool computeConjugator( const Permutation& p , Permutation& res ) const;
  
  
  //! Flip the permutation (conjugate by a half-twist permutation \f$\Delta\f$)
  Permutation flip( ) const;
  
  
  //! Perform a tiny flip (conjugate by \f$\delta\f$)
  Permutation tinyFlip( int sh ) const;
  
  
  //! Generate a random permutation of specified size
  static Permutation random( int size );
  
  
  //! Get half twist permutation \f$\Delta = (n-1,n-2,\ldots,2,1,0)\f$
  static Permutation getHalfTwistPermutation( int size );
  //! Get half twist permutation \f$\Delta = (1,2,\ldots,n-2,n-1,0)\f$
  static Permutation getCyclePermutation( int size );
  
  
  //! Check if 2 permutations are mixable (it is not used anywhere in the system, god knows that means)
  static bool mixable( const Permutation& p1 , const Permutation& p2 );
  

  //! Check if the permutation is trivial
  bool isTrivial( ) const;


  //! Find a geodesic word representing the permutation (indices start with 0)
  vector<int> geodesic( ) const;

  //! Find a geodesic word representing the permutation (indices start with 1)
  vector<int> geodesicWord( ) const;

  
  //! Compute the word presentation for a permutation which is "parallel descending cycles" (for other permutations it makes no sense)
  vector< int > getWordPresentation( ) const;
  
  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Lattice operations:                                //
  //                                                     //
  /////////////////////////////////////////////////////////
  
public:
  
  //! Compute RightGCD of 2 permutations
  /*!
    Let \f$p_1\f$ and \f$p_2\f$ be two permutations. The maximal permutation \f$P\f$
    which ends \f$p_1\f$ and \f$p_2\f$ is called the
    right greatest common divisor of \f$p_1\f$ and \f$p_2\f$, i.e., \f$P\f$ is maximal such that
    \f$p_1 = d_1 \circ P \f$ and \f$p_2 = d_2 \circ P \f$ for some permutations \f$d_1\f$ and \f$d_2\f$.
  */
  Permutation RightGCD( const Permutation& p ) const;
  
  
  //! Compute RightLCM of 2 permutations
  /*!
    Let \f$p_1\f$ and \f$p_2\f$ be two permutations. The minimal permutation \f$P\f$
    which ends with \f$p_1\f$ and \f$p_2\f$ is called the
    right least common multiple, i.e., \f$P\f$ is minimal such that
    \f$P = d_1 \circ p_1 \f$ and \f$P = d_2 \circ p_2 \f$ .
  */
  Permutation RightLCM( const Permutation& p ) const;
  
  
  //! Compute LeftGCD of 2 permutations
  /*!
    Let \f$p_1\f$ and \f$p_2\f$ be two permutations. The maximal permutation \f$P\f$
    which ends \f$p_1\f$ and \f$p_2\f$ is called the
    left greatest common divisor of \f$p_1\f$ and \f$p_2\f$, i.e., \f$P\f$ is maximal such that
    \f$p_1 = P \circ d_1 \f$ and \f$p_2 = P \circ d_1 \f$ for some permutations \f$d_1\f$ and \f$d_2\f$.
  */
  Permutation LeftGCD( const Permutation& p ) const;
  
  
  //! Compute LeftLCM of 2 permutations
  /*!
    Let \f$p_1\f$ and \f$p_2\f$ be two permutations. The minimal permutation \f$P\f$
    which starts with \f$p_1\f$ and \f$p_2\f$ is called the
    left least common multiple, i.e., \f$P\f$ is minimal such that
    \f$P = p_1 \circ d_1 \f$ and \f$P = p_2 \circ d_2 \f$ .
  */
  Permutation LeftLCM( const Permutation& p ) const;
  
  
  //! Compute GCD for 2 permutations that are "parallel descending cycles" (used for BKL Normal Forms of braids)
  Permutation meet2( const Permutation& p ) const;
  
  
  //! Compute LCM for 2 permutations that are "parallel descending cycles" (used for BKL Normal Forms of braids)
  Permutation join2( const Permutation& p ) const;
  

  /////////////////////////////////////////////////////////
  //                                                     //
  //  I/O:                                               //
  //                                                     //
  /////////////////////////////////////////////////////////

public:

  friend ostream& operator << ( ostream& os , const Permutation& p );

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Internal functions:                                //
  //                                                     //
  /////////////////////////////////////////////////////////
  
private:
  
  //! (Aux) The main operation to compute RightGCD and all other lattice functions
  void _sub_meet( const Permutation& p , 
		  const Permutation& ip1 ,
		  const Permutation& ip2 ,
		  Permutation& cur , 
		  int* left_indeces_a ,
		  int* left_indeces_b ,
		  int* right_indeces_a ,
		  int* right_indeces_b ,
		  int beg , int end ) const;


  struct triple {
    triple( int p1=0 , int p2=0 , int p3=0 ) : c1(p1), c2(p2), c3(p3) { }
    int c1, c2, c3;
  };
  
  void prepare_pairs( int N, 
		      Permutation& P, 
		      vector< pair<int,int> >& pairs ) const;
 
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Data members                                       //
  //                                                     //
  /////////////////////////////////////////////////////////

  friend class BraidGroup;
  
private:

  //! A vector representing the permutation.
  vector< int > theValue;

};




#endif
