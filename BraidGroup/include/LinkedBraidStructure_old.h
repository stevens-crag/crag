// Copyright (C) 2005 Alexander Ushakov
// Contents: Definition of class BraidNode
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//

#ifndef _LinkedBraidStructure_h_
#define _LinkedBraidStructure_h_


#include "list"
#include "vector"
using namespace std;


//---------------------------------------------------------------------------//
//------------------------------ BraidNode ----------------------------------//
//---------------------------------------------------------------------------//

//! Defines a crossing in a linked braid structure//
/*!
  Represented by the orientation (type) of the crossing
  and 6 pointers to other crossings
*/


struct BraidNode
{
  
  // BraidNode( ) : left( 0 ), ahead( 0 ), right( 0 ),
  // back_left( 0 ), back( 0 ), back_right( 0 ) , marked( 0 ) { }
  
  //! Default constructor (creates a disconnected node - all pointers are NULL)
  BraidNode( bool m=false ) : left( 0 ), ahead( 0 ), right( 0 ),
			back_left( 0 ), back( 0 ), back_right( 0 ) , marked( m ) { }
  
  //! pointer to the ahead-left node
/*!
  if \f$w = w_1 \circ x_i \circ w_2\f$ and the node corresponds to \f$x_i\f$
  then left points to the rightmost \f$x_{i-1}\f$ in \f$w_1\f$. Except one
  case! If between \f$x_{i-1}\f$ and \f$x_{i}\f$ there is another \f$x_{i}\f$
  then left = NULL.
*/
  BraidNode* left;
  
  //! pointer to the ahead node
/*!
  if \f$w = w_1 \circ x_i \circ w_2\f$ and the node corresponds to \f$x_i\f$
  then ahead points to the rightmost \f$x_{i}\f$ in \f$w_1\f$.
*/
  BraidNode* ahead;

  //! pointer to the ahead-right node
/*!
  if \f$w = w_1 \circ x_i \circ w_2\f$ and the node corresponds to \f$x_i\f$
  then right points to the rightmost \f$x_{i+1}\f$ in \f$w_1\f$. Except one
  case! If between \f$x_{i+1}\f$ and \f$x_{i}\f$ there is another \f$x_{i}\f$
  then right = NULL.
*/
  BraidNode* right;

  //! pointer to the back-left node
/*!
  if \f$w = w_1 \circ x_i \circ w_2\f$ and the node corresponds to \f$x_i\f$
  then back_left points to the leftmost \f$x_{i-1}\f$ in \f$w_2\f$. Except one
  case! If between \f$x_{i}\f$ and \f$x_{i-1}\f$ there is another \f$x_{i}\f$
  then back_left = NULL.
*/
  BraidNode* back_left;

  //! pointer to the back node
/*!
  if \f$w = w_1 \circ x_i \circ w_2\f$ and the node corresponds to \f$x_i\f$
  then ahead points to the leftmost \f$x_{i}\f$ in \f$w_2\f$.
*/
  BraidNode* back;

  //! pointer to the back-right node
/*!
  if \f$w = w_1 \circ x_i \circ w_2\f$ and the node corresponds to \f$x_i\f$
  then back_right points to the leftmost \f$x_{i+1}\f$ in \f$w_2\f$. Except one
  case! If between \f$x_{i}\f$ and \f$x_{i+1}\f$ there is another \f$x_{i}\f$
  then back_right = NULL.
*/
  BraidNode* back_right;

  // bool ahead_left_order;
  // bool ahead_right_order;
  // shows the order of links ahead (0 - center first, 1 - sides first)
  // bool back_left_order;
  // bool back_right_order;
  // shows the order of links back (0 - center first, 1 - sides first)

  // static void addNode( vector< BraidNode* >& frontNodes , int g );
  static bool isHandle( vector< BraidNode* >& backNodes , int g );
  static void removeHandle( vector< BraidNode* >& frontNodes , vector< BraidNode* >& backNodes , int g , list< int >& processedHandle );
  static void removeNodeTerminalClosure( vector< BraidNode* >& frontNodes , vector< BraidNode* >& backNodes , int g , BraidNode* node , list< int >& removedClosure );
  static void removeTerminalNode( vector< BraidNode* >& frontNodes , vector< BraidNode* >& backNodes , int g , BraidNode* node , list< int >& removedClosure );

  // static list< int > translateIntoWord( vector< BraidNode* >& frontNodes );

  //! crossing orientation (true=positive, false=negative)
  bool type;
  
  //! Auxiliary flag. Used in LinkedBraidStructure::getWord( ) to mark used crossings.
  mutable bool marked;
};


//---------------------------------------------------------------------------//
//--------------------------- LinkedBraidStructure --------------------------//
//---------------------------------------------------------------------------//


class LinkedBraidStructure
{
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Constructors                                       //
  //                                                     //
  /////////////////////////////////////////////////////////
 public:

  LinkedBraidStructure( int N );
  ~LinkedBraidStructure( );
  
  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Accessors:                                         //
  //                                                     //
  /////////////////////////////////////////////////////////

 public:
  
  void push_back( int i );
  template< class ConstIntIterator >
    void push_back( ConstIntIterator B , ConstIntIterator E ) {
    for( ; B!=E ; ++B )
      push_back( *B );
  }
  
  void clear( );
  void push_front( int i );

  void transformToDehornoyForm( );
  list< int > getWord( ) const;

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Internal functions:                                //
  //                                                     //
  /////////////////////////////////////////////////////////
private:

  //! Function used in getWord( )
  void processNode( list< int >& result , int index , const BraidNode* bn ) const;

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Data members                                       //
  //                                                     //
  /////////////////////////////////////////////////////////
private:

  int theRank;
  // specifies the number of strands (number of generators + 1)
  
  vector< BraidNode* > frontNodes;
  vector< BraidNode* > backNodes;
  
  // Auxiliary variable. Used in getWord to specify the type of the last mark.
  mutable bool theMark;
};


#endif

