// Copyright (C) 2005 Alexander Ushakov
// Contents: Definition of class BalancedTree
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//


#ifndef _BalancedTree_H_
#define _BalancedTree_H_

#include "list"
#include "iostream"
using namespace std;

//---------------------------------------------------------------------------//
//------------------------------ BalancedTree -------------------------------//
//---------------------------------------------------------------------------//


//! Class BalancedTree (container class intended to keep an ordered sequence of objects, supports fast \f$n \log n \f$ ninsertions)
/*!
  Typical usage for this class is generation of trivial elements of groups, where a lot
  of insertions of relators is performed, i.e., start from a trivial word and cosequently
  insert relators (with their inverses and cyclic permutations). The result will be a word
  representing the identity of the group. The main feature of this class - fast \f$n \log n \f$ 
  asymptotic insertions which cannot be achieved using standard vector<...> and list<...>.
*/



template< class Obj >
class BalancedTree
{
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Internal structures:                               //
  //                                                     //
  /////////////////////////////////////////////////////////

 private:
  
  struct BTNode {
    BTNode(         ) : subtree1(0), subtree2(0), height1(0), height2(0), weight1(0), weight2(0) { }
    BTNode( Obj obj ) : subtree1(0), subtree2(0), height1(0), height2(0), weight1(0), weight2(0), theObject(obj) { }
    ~BTNode( ) {
      if( subtree1 ) delete subtree1;
      if( subtree2 ) delete subtree2;
    }

    Obj     theObject;
    BTNode* subtree1;   int height1;  int weight1;
    BTNode* subtree2;   int height2;  int weight2;
  };

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Constructors                                       //
  //                                                     //
  /////////////////////////////////////////////////////////

 public:

  //! Default constructor.
  BalancedTree( ) : theRoot(0) { }
  //! Constructor. Creates a tree for a sequence in range [B,E)
  template< class ConstObjIter > BalancedTree( const ConstObjIter& B , const ConstObjIter& E ) : theRoot(0) {
    insert( 0 , B , E );
  }
  //! Destructor recursively frees memory
  ~BalancedTree( ) {
    if( theRoot ) delete theRoot;
  }
  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Operators                                          //
  //                                                     //
  /////////////////////////////////////////////////////////
 public:

  int size( ) const {
    if( !theRoot ) return 0;
    return theRoot->weight1 + theRoot->weight2 + 1;
  }

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Accessors                                          //
  //                                                     //
  /////////////////////////////////////////////////////////
	
public:

  //! Insert a sequence of objects bounded by [B,E) into a position pos.
  template< class ConstObjIter > 
  void insert( int pos , const ConstObjIter& B , const ConstObjIter& E ) {
    int i=pos;
    for( ConstObjIter it=B ; it!=E ; ++it, ++i )
      insert( i , *it );
  }
  
  //! Insert an object obj into position pos.
  void insert( int pos , const Obj& obj ) {
    if( theRoot==0 )
      theRoot = new BTNode( obj );
    else
      insert( NULL , theRoot, pos, obj );
  }

  //! Get an ordered list of objects stored in a balanced tree.
  list< Obj > getList( ) const {
    list< Obj > result;
    getList( theRoot , result );
    return result;
  }

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Internal functions                                 //
  //                                                     //
  /////////////////////////////////////////////////////////

 private:

  //! Get an ordered list of objects in a subtree with a root node.
  void getList( BTNode* node , list< Obj >& result ) const {
    if( !node ) return;
    getList( node->subtree1 , result );
    result.push_back( node->theObject );
    getList( node->subtree2 , result );
  }

  //! Insert an onject obj into a subtree with a root node at the offset position offset
  void insert( BTNode* parent , BTNode* node , int offset , const Obj& obj ) {

    // going down to insert a new element
    if( offset<=node->weight1 ) {
      if( node->weight1==0 )
	node->subtree1 = new BTNode( obj );
      else
	insert( node , node->subtree1 , offset , obj );
      node->weight1++;
      node->height1 = 1 +
	( node->subtree1->height1 > node->subtree1->height2 ? 
	  node->subtree1->height1 : node->subtree1->height2 ); 
    } else {
      if( !node->subtree2 )
	node->subtree2 = new BTNode( obj );
      else
	insert( node , node->subtree2 , offset-node->weight1-1 , obj );
      node->weight2++;
      node->height2 = 1 +
	( node->subtree2->height1 > node->subtree2->height2 ? 
	  node->subtree2->height1 : node->subtree2->height2 ); 
    }

    // going up and balance the branches if required
    if( node->height1<=node->height2-2 ) {
      rotateLeft( parent , node );
    } else if( node->height1-2>=node->height2 ) {
      rotateRight( parent , node );
    }
  }

  //! Left rotation.
  void rotateLeft( BTNode* parent , BTNode* child ) {
    BTNode* right = child->subtree2;
    if( right->height1<=right->height2 ) { // 1st type of rotation
      if( parent ) {
	if( child==parent->subtree1 )
	  parent->subtree1 = right;
	else
	  parent->subtree2 = right;
      } else
	theRoot = right;
      child->subtree2 = right->subtree1;
      right->subtree1 = child;

      child->height2 = right->height1;
      child->weight2 = right->weight1;
      right->height1 = ( child->height1 > child->height2 ? child->height1 : child->height2 ) + 1;
      right->weight1 = child->weight1 + child->weight2 + 1;
    } else { // 2nd type of rotation
      BTNode* right_left = right->subtree1;

      if( parent ) {
	if( child==parent->subtree1 )
	  parent->subtree1 = right_left;
	else
	  parent->subtree2 = right_left;
      } else
	theRoot = right_left;
      
      // child
      child->subtree2 = right_left->subtree1;
      child->weight2 = right_left->weight1;
      child->height2 = right_left->height1;
      
      // right
      right->subtree1 = right_left->subtree2;
      right->weight1 = right_left->weight2;
      right->height1 = right_left->height2;
      
      // right_left
      right_left->subtree1 = child;
      right_left->height1 = ( child->height1 > child->height2 ? child->height1 : child->height2 ) + 1;
      right_left->weight1 = child->weight1 + child->weight2 + 1;
      right_left->subtree2 = right;
      right_left->height2 = ( right->height1 > right->height2 ? right->height1 : right->height2 ) + 1;
      right_left->weight2 = right->weight1 + right->weight2 + 1;
    }
  }

  //! Right rotation.
  void rotateRight( BTNode* parent , BTNode* child ) {
    BTNode* left = child->subtree1;
    if( left->height2<=left->height1 ) { // 1st type of rotation
      
      if( parent ) {
	if( child==parent->subtree1 )
	  parent->subtree1 = left;
	else
	  parent->subtree2 = left;
      } else
      theRoot = left;
      child->subtree1 = left->subtree2;
      left->subtree2 = child;
      
      child->height1 = left->height2;
      child->weight1 = left->weight2;
      left->height2 = ( child->height1 > child->height2 ? child->height1 : child->height2 ) + 1;
      left->weight2 = child->weight1 + child->weight2 + 1;

    } else {
      BTNode* left_right = left->subtree2;
      if( parent ) {
	if( child==parent->subtree1 )
	  parent->subtree1 = left_right;
	else
	  parent->subtree2 = left_right;
      } else
	theRoot = left_right;
      
      // child
      child->subtree1 = left_right->subtree2;
      child->weight1  = left_right->weight2;
      child->height1  = left_right->height2;
      
      // left
      left->subtree2 = left_right->subtree1;
      left->weight2  = left_right->weight1;
      left->height2  = left_right->height1;
      
      // right_left
      left_right->subtree1 = left;
      left_right->height1 = ( left->height1 > left->height2 ? left->height1 : left->height2 ) + 1;
      left_right->weight1 = left->weight1 + left->weight2 + 1;
      left_right->subtree2 = child;
      left_right->height2 = ( child->height1 > child->height2 ? child->height1 : child->height2 ) + 1;
      left_right->weight2 = child->weight1 + child->weight2 + 1;
    }
  }

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Data members                                       //
  //                                                     //
  /////////////////////////////////////////////////////////

private:
  
  BTNode* theRoot;

};


#endif
