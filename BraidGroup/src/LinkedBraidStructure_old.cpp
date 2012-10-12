// Copyright (C) 2005 Alexander Ushakov
// Contents: Implementation of class BraidNode
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//

#include "Word.h"
#include "iostream"
#include "LinkedBraidStructure.h"


bool BraidNode::isHandle( vector< BraidNode* >& backNodes , int l )
{
  int al = abs( l );

  BraidNode* node = backNodes[al-1];
  if( !node || node->left )
    return false;

  BraidNode* ahead = node->ahead;
  if( !ahead || ahead->back_left )
    return false;

  return node->type!=ahead->type;
}


void BraidNode::removeHandle( vector< BraidNode* >& frontNodes , vector< BraidNode* >& backNodes , int l , list< int >& processedHandle )
{
  list< int > handle;

  int al = abs( l );
  BraidNode* node = backNodes[al-1];
  BraidNode* ahead = node->ahead;
  
  removeNodeTerminalClosure( frontNodes , backNodes , al , node , handle );
  removeNodeTerminalClosure( frontNodes , backNodes , al , ahead , handle );
  // cout << "  handle = " << Word( handle ) << endl;
  handle.pop_front( );
  handle.pop_back( );

  // cout << "  handle2 = " << Word( handle ) << endl;

  // process the handle
  list< int >::iterator handle_it = handle.begin( );
  for( ; handle_it!=handle.end( ) ; ++handle_it ) {
    int l_1 = *handle_it;
    int al_1 = abs( l_1 );
    if( al_1==al+1 ) {
      handle.insert( handle_it , l>0 ? l+1 : l-1 );
      *handle_it = l_1>0 ? l_1-1 : l_1+1;
      handle.insert( ++handle_it , l>0 ? -l-1 : -l+1 );
      --handle_it;
    }
  }

  // cout << "  handle3 = " << Word( handle ) << endl;
  for( handle_it = handle.end( ) ; handle_it!=handle.begin( ) ; )
    processedHandle.push_front( *--handle_it );
}


void BraidNode::removeNodeTerminalClosure( vector< BraidNode* >& frontNodes , vector< BraidNode* >& backNodes , int al , BraidNode* node , list< int >& removedClosure )
{
  if( node->back_left )
    removeNodeTerminalClosure( frontNodes , backNodes , al-1 , node->back_left , removedClosure );
  if( node->back_right ) 
    removeNodeTerminalClosure( frontNodes , backNodes , al+1 , node->back_right , removedClosure );
  if( node->back )
    removeNodeTerminalClosure( frontNodes , backNodes , al , node->back , removedClosure );

  removeTerminalNode( frontNodes , backNodes , al , node , removedClosure );
}


void BraidNode::removeTerminalNode( vector< BraidNode* >& frontNodes , vector< BraidNode* >& backNodes , int al , BraidNode* node , list< int >& removedClosure )
{
  if( node->left )
    node->left->back_right = 0;
  if( node->right )
    node->right->back_left = 0;
  if( node->ahead )
    node->ahead->back = 0;

  backNodes[al-1] = node->ahead;
  if( node==frontNodes[al-1] ) frontNodes[al-1] = NULL;
  removedClosure.push_front( node->type ? al : -al );
  delete node;
}

/*
list< int > BraidNode::translateIntoWord( vector< BraidNode* >& frontNodes )
{
  list< int > result;

  for( int i=0 ; i<frontNodes.size( ) ; ++i )
    while( frontNodes[i] )
      removeNodeTerminalClosure( frontNodes , i+1 , frontNodes[i] , result );

  return result;
}
*/

//---------------------------------------------------------------------------//
//--------------------------- LinkedBraidStructure --------------------------//
//---------------------------------------------------------------------------//


LinkedBraidStructure::LinkedBraidStructure( int N ) :
  theRank( N ) ,
  theMark( 0 ) ,
  frontNodes( N-1 , NULL ) ,
  backNodes ( N-1 , NULL )
{

}


LinkedBraidStructure::~LinkedBraidStructure( )
{
  clear( );
}

//---------------------------------------------------------------------------//
//--------------------------------- push_back -------------------------------//
//---------------------------------------------------------------------------//


void LinkedBraidStructure::push_back( int l )
{
  int al = abs(l);
  BraidNode* left  = al>1 ? backNodes[al-2] : 0;
  BraidNode* ahead = backNodes[al-1];
  BraidNode* right = al<theRank-1 ? backNodes[al] : 0;

  BraidNode* newNode = new BraidNode( theMark );
  
  newNode->ahead = ahead;
  if( ahead ) 
    ahead->back = newNode;

  if( !ahead || ahead && ahead->back_left ) {
    newNode->left = left;
    if( left )
      left->back_right = newNode;
  }
  if( !ahead || ahead && ahead->back_right ) {
    newNode->right = right;
    if( right )
      right->back_left = newNode;
  }
  newNode->type = l>0;
  backNodes[al-1] = newNode;
  
  if( !frontNodes[al-1] )
    frontNodes[al-1] = newNode;
}


//---------------------------------------------------------------------------//
//--------------------------------- push_front ------------------------------//
//---------------------------------------------------------------------------//


void LinkedBraidStructure::push_front( int l )
{
  int al = abs(l);
  BraidNode* back_left  = al>1 ? frontNodes[al-2] : 0;
  BraidNode* back = frontNodes[al-1];
  BraidNode* back_right = al<theRank-1 ? frontNodes[al] : 0;
  
  BraidNode* newNode = new BraidNode( theMark );
  
  newNode->back = back;
  if( back ) 
    back->ahead = newNode;

  if( !back || back && back->left ) {
    newNode->back_left = back_left;
    if( back_left )
      back_left->right = newNode;
  }
  if( !back || back && back->right ) {
    newNode->back_right = back_right;
    if( back_right )
      back_right->left = newNode;
  }
  newNode->type = l>0;
  frontNodes[al-1] = newNode;
  
  if( !backNodes[al-1] )
    backNodes[al-1] = newNode;
}


//---------------------------------------------------------------------------//
//----------------------------------- getWord -------------------------------//
//---------------------------------------------------------------------------//


list< int > LinkedBraidStructure::getWord( ) const
{
  list< int > result;
  
  theMark = !theMark;
  // we consider each node starting from the back
  // for( int i=theRank-2 ; i>=0 ; --i )
  for( int i=0 ; i<theRank-1 ; ++i )
    for( BraidNode* cur_node = backNodes[i] ; cur_node ; cur_node=cur_node->ahead )
      processNode( result , i+1 , cur_node );
  
  return result;
}


//---------------------------------------------------------------------------//
//-------------------------------- processNode ------------------------------//
//---------------------------------------------------------------------------//


void LinkedBraidStructure::processNode( list< int >& result , int index , const BraidNode* bn ) const
{
  if( bn->marked==theMark )
    return;
  
  if( bn->back_left )
    processNode( result , index-1 , bn->back_left );
  if( bn->back_right )
    processNode( result , index+1 , bn->back_right );
  if( !bn->back_left && !bn->back_right && bn->back )
    processNode( result , index , bn->back );
  result.push_front( bn->type ? index : -index );
  bn->marked = theMark;
}

void LinkedBraidStructure::clear( )
{
  for( int i=0 ; i<theRank-1 ; ++i ) {
    for( BraidNode* cur_node = backNodes[i] ; cur_node ; ) {
      BraidNode* next_node = cur_node->ahead;
      delete cur_node;
      cur_node = next_node;
    }
  }
  
  frontNodes = vector< BraidNode* >( theRank-1 , NULL );
  backNodes  = vector< BraidNode* >( theRank-1 , NULL );
}


void LinkedBraidStructure::transformToDehornoyForm( )
{
  Word w = getWord( );
  clear( );
  
  list< int > lw = w.getList( );
  while( lw.begin()!=lw.end( ) ) {
    int l = *lw.begin( );
    // cout << "  " << l << endl;
    lw.erase( lw.begin( ) );
    push_back( l );
    if( BraidNode::isHandle( backNodes , l ) ) {
      // cout << "    Handle" << endl;
      BraidNode::removeHandle( frontNodes , backNodes , l , lw );
    }
  }
  
}
