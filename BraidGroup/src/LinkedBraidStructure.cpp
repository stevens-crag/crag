
#include "Word.h"
#include "iostream"
#include "LinkedBraidStructure.h"



ostream& operator << ( ostream& os , const BraidNode& bn )
{
  os << &bn << "{ " 
     << bn.type << "| "
     << bn.left << ", "
     << bn.ahead << ", "
     << bn.right << ", "
     << bn.back_left << ", "
     << bn.back << ", "
     << bn.back_right << " }";

  return os;
}


//---------------------------------------------------------------------------//
//---------------------------- makeTransform --------------------------------//
//---------------------------------------------------------------------------//


LinkedBraidStructureTransform LinkedBraidStructure::make_EraseTransform( BraidNode* bn , int pos ) const 
{ 
  return LinkedBraidStructureTransform( bn->theNumber , 
					pos , 
					LinkedBraidStructureTransform::ERASED , 
					bn->type ,
					bn->left       ? bn->left      ->theNumber : -1 ,
					bn->ahead      ? bn->ahead     ->theNumber : -1 ,
					bn->right      ? bn->right     ->theNumber : -1 ,
					bn->back_left  ? bn->back_left ->theNumber : -1 ,
					bn->back       ? bn->back      ->theNumber : -1 ,
					bn->back_right ? bn->back_right->theNumber : -1 );
  
}

LinkedBraidStructureTransform LinkedBraidStructure::make_AddTransform( BraidNode* bn , int pos ) const 
{ 
  return LinkedBraidStructureTransform( bn->theNumber , 
					pos , 
					LinkedBraidStructureTransform::ADDED , 
					bn->type );
}


LinkedBraidStructureTransform LinkedBraidStructure::make_ChangeType( BraidNode* bn , int pos ) const
{
  return LinkedBraidStructureTransform( bn->theNumber , 
					pos , 
					LinkedBraidStructureTransform::CHANGE_TYPE , 
					bn->type );
}


//---------------------------------------------------------------------------//
//----------------------------- operator <  ---------------------------------//
//---------------------------------------------------------------------------//


bool LinkedBraidStructure::operator < ( const LinkedBraidStructure& lbs ) const
{
  if( size()<lbs.size() )
    return true;
  if( size()>lbs.size() )
    return false;

  if( theIndex<lbs.theIndex )
    return true;
  if( theIndex>lbs.theIndex )
    return false;

  return translateIntoWord( )<lbs.translateIntoWord( );
}


//---------------------------------------------------------------------------//
//------------------------- LinkedBraidStructure ----------------------------//
//---------------------------------------------------------------------------//


LinkedBraidStructure::LinkedBraidStructure( int N ) :
  maxNodeNumber( 0 ),
  theIndex( N ),
  backNodes ( vector< BraidNode* >( N , (BraidNode*)NULL ) ),
  frontNodes( vector< BraidNode* >( N , (BraidNode*)NULL ) )
{
  
}


//---------------------------------------------------------------------------//
//------------------------- LinkedBraidStructure ----------------------------//
//---------------------------------------------------------------------------//


LinkedBraidStructure::LinkedBraidStructure( int N , const Word& w ) :
  maxNodeNumber( 0 ),
  theIndex( N ),
  backNodes ( vector< BraidNode* >( N , (BraidNode*)NULL ) ),
  frontNodes( vector< BraidNode* >( N , (BraidNode*)NULL ) )
{
  for( ConstWordIterator w_it=w.begin( ) ; w_it!=w.end( ) ; ++w_it )
    push_back( *w_it );
}


//---------------------------------------------------------------------------//
//------------------------- LinkedBraidStructure ----------------------------//
//---------------------------------------------------------------------------//


LinkedBraidStructure& LinkedBraidStructure::operator= ( const LinkedBraidStructure& LBS )
{
  // cout << "   #" << endl;

  theIndex = LBS.theIndex;
  theNodes = LBS.theNodes;
  maxNodeNumber = LBS.maxNodeNumber;
  backNodes  = LBS.backNodes;
  frontNodes = LBS.frontNodes;
  
  map< int , BraidNode >::const_iterator n_it1 = LBS.theNodes.begin( );
  map< int , BraidNode >::iterator       n_it2 =     theNodes.begin( );
  for( ; n_it1!=LBS.theNodes.end( ) ; ++n_it1, ++n_it2 )
    (*n_it1).second.link = &(*n_it2).second;
  
  n_it2 = theNodes.begin( );
  for( ; n_it2!=theNodes.end( ) ; ++n_it2 ) {
    BraidNode& cur = (*n_it2).second;
    if( cur.left       ) cur.left       = cur.left->link;
    if( cur.ahead      ) cur.ahead      = cur.ahead->link;
    if( cur.right      ) cur.right      = cur.right->link;
    if( cur.back_left  ) cur.back_left  = cur.back_left->link;
    if( cur.back       ) cur.back       = cur.back->link;
    if( cur.back_right ) cur.back_right = cur.back_right->link;
  }

  for( int i=0 ; i<theIndex ; ++i ) {
    if( backNodes [i] ) backNodes [i] = backNodes [i]->link;
    if( frontNodes[i] ) frontNodes[i] = frontNodes[i]->link;
  }
 
  LBS.clearLinks( );

  return *this;
}


//---------------------------------------------------------------------------//
//------------------------- LinkedBraidStructure ----------------------------//
//---------------------------------------------------------------------------//


LinkedBraidStructure::LinkedBraidStructure( const LinkedBraidStructure& LBS ) :
  theIndex( LBS.theIndex ),
  theNodes( LBS.theNodes ),
  maxNodeNumber( LBS.maxNodeNumber ),
  backNodes( LBS.backNodes ),
  frontNodes( LBS.frontNodes )
{
  // cout << "   &" << endl;
  map< int , BraidNode >::const_iterator n_it1 = LBS.theNodes.begin( );
  map< int , BraidNode >::iterator       n_it2 =     theNodes.begin( );
  for( ; n_it1!=LBS.theNodes.end( ) ; ++n_it1, ++n_it2 )
    (*n_it1).second.link = &(*n_it2).second;
  
  n_it2 = theNodes.begin( );
  for( ; n_it2!=theNodes.end( ) ; ++n_it2 ) {
    BraidNode& cur = (*n_it2).second;
    if( cur.left       ) cur.left       = cur.left->link;
    if( cur.ahead      ) cur.ahead      = cur.ahead->link;
    if( cur.right      ) cur.right      = cur.right->link;
    if( cur.back_left  ) cur.back_left  = cur.back_left->link;
    if( cur.back       ) cur.back       = cur.back->link;
    if( cur.back_right ) cur.back_right = cur.back_right->link;
  }

  for( int i=0 ; i<theIndex ; ++i ) {
    if( backNodes [i] ) backNodes [i] = backNodes [i]->link;
    if( frontNodes[i] ) frontNodes[i] = frontNodes[i]->link;
  }

  LBS.clearLinks( );
}


//---------------------------------------------------------------------------//
//--------------------------------- push_back -------------------------------//
//---------------------------------------------------------------------------//


LinkedBraidStructureTransform LinkedBraidStructure::push_back( int g )
{
  int ag = abs(g);
  BraidNode* left  = ag>1 ? backNodes[ag-2] : NULL;
  BraidNode* ahead = backNodes[ag-1];
  BraidNode* right = ag<theIndex ? backNodes[ag] : NULL;

  BraidNode& newNode = theNodes[maxNodeNumber] = BraidNode( maxNodeNumber , g>0 );
  maxNodeNumber++;
  
  newNode.ahead = ahead;
  if( ahead ) 
    ahead->back = &newNode;

  if( !ahead || ahead && ahead->back_left ) {
    newNode.left = left;
    if( left )
      left->back_right = &newNode;
  }
  if( !ahead || ahead && ahead->back_right ) {
    newNode.right = right;
    if( right )
      right->back_left = &newNode;
  }
  if( !frontNodes[ag-1] )
    frontNodes[ag-1] = &newNode;
  backNodes[ag-1] = &newNode;

  return make_AddTransform( &newNode , ag-1 );
}


//---------------------------------------------------------------------------//
//-------------------------------- push_front -------------------------------//
//---------------------------------------------------------------------------//


LinkedBraidStructureTransform LinkedBraidStructure::push_front( int g )
{
  int ag = abs(g);
  BraidNode* left  = ag>1 ? frontNodes[ag-2] : 0;
  BraidNode* ahead = frontNodes[ag-1];
  BraidNode* right = ag<theIndex ? frontNodes[ag] : 0;

  BraidNode& newNode = theNodes[maxNodeNumber] = BraidNode( maxNodeNumber , g>0 );
  maxNodeNumber++;

  newNode.back = ahead;
  if( ahead ) 
    ahead->ahead = &newNode;

  if( !ahead || ahead && ahead->left ) {
    newNode.back_left = left;
    if( left )
      left->right = &newNode;
  }
  if( !ahead || ahead && ahead->right ) {
    newNode.back_right = right;
    if( right )
      right->left = &newNode;
  }
  
  if( !backNodes[ag-1] )
    backNodes[ag-1] = &newNode;
  frontNodes[ag-1] = &newNode;

  return make_AddTransform( &newNode , ag-1 );
}


//---------------------------------------------------------------------------//
//---------------------------- removeLeftHandles ----------------------------//
//---------------------------------------------------------------------------//


void LinkedBraidStructure::removeLeftHandles( list< LinkedBraidStructureTransform >* result )
{
  // set of nodes to check (1st value = number of uppercrossings, 2nd value = -position in the Linked Structure)
  typedef triple< int , int , BraidNode* > NODE;
  set< NODE > to_check;

  // form the initial set of nodes to check (all nodes from the Linked Structure)
  for( int i=0 ; i<theIndex ; ++i ) {
    for( BraidNode* c = frontNodes[i] ; c!=backNodes[i] ; c = c->back ) {
      int weight = checkIfStartsLeftHandle( i , c );
      if( weight!=-1 ) {
	c->weight=weight;
	to_check.insert( NODE( weight , i , c ) );
      }
    }
  }
  
  // check the nodes one by one starting from the rightmost
  for( ; !to_check.empty( ) ; ) {
    
    // a. choose the node and remove it from to_check
    NODE cur_node = *to_check.begin( );
    to_check.erase( to_check.begin( ) );
    
    // b. remove the handle
    int weight = checkIfStartsLeftHandle( cur_node.second , cur_node.third );
    if( weight!=-1 )
      removeLeftHandle( cur_node , to_check , result );
  }
}


//---------------------------------------------------------------------------//
//------------------------ checkIfStartsLeftHandle --------------------------//
//---------------------------------------------------------------------------//


int LinkedBraidStructure::checkIfStartsLeftHandle( int pos , BraidNode* bn )
{
  BraidNode* back = bn->back;
  // if there is no more x_i
  if( !back )
    return -1;
  
  // if there is an obstacle for handle
  if( bn->back_right )
    return -1;

  // if crossings have the same orientation
  if( bn->type==back->type )
    return -1;

  // compute the number of upper crossings
  BraidNode* l1 = bn->back_left;
  BraidNode* l2 = back->left;
  int counter=0;
  bool crossingType;
  for( ; l1 ; l1 = l1->back ) {
    if( counter==0 ) {
      crossingType = l1->type;
    } else {
      if( crossingType!=l1->type )
	return -1;
    }
    ++counter;
    if( l1==l2 )
      break;
  }
  
  return counter;
}


//---------------------------------------------------------------------------//
//----------------------------- removeLeftHandle ----------------------------//
//---------------------------------------------------------------------------//


void LinkedBraidStructure::removeLeftHandle( triple< int , int , BraidNode* > node , set< triple< int , int , BraidNode* > >& to_check , list< LinkedBraidStructureTransform >* lst )
{
  typedef triple< int , int , BraidNode* > NODE;
  
  int pos = node.second;
  BraidNode* n1 = node.third;
  BraidNode* n2 = n1->back;
  bool type = n1->type;

  BraidNode* l1 = n1->back_left;
  BraidNode* l2 = n2->left;

  // Removal of a handle can introduce new handles. Here we store some nodes to check. A few will be added later.
  set< pair< int , BraidNode* > > to_check2;
  if( n1->ahead ) to_check2.insert( pair< int , BraidNode* >( pos   , n1->ahead ) );
  if( n1->left  ) to_check2.insert( pair< int , BraidNode* >( pos-1 , n1->left  ) );
  for( BraidNode* cn=n1 ; cn ; cn=cn->ahead ) {
    if( !cn->right )
      continue;
    to_check2.insert( pair< int , BraidNode* >( pos+1 , cn->right ) );
    break;
  }
  
  // A. process left nodes
  int counter=0;
  for( ; l1 ; ) {
    // ++counter;
    BraidNode* l3=l1->back;
    BraidNode* new_node  = insertBackRight ( l1 , pos , l1->type );
    BraidNode* new_node2 = insertBackLeft( new_node , pos-1 , type );

    if( lst ) {
      lst->push_back( make_AddTransform( new_node , pos ) );
      lst->push_back( make_AddTransform( new_node2 , pos-1 ) );
      lst->push_back( make_ChangeType( l1 , pos ) );
    }
    
    l1->type = !type;
    to_check2.insert( pair< int , BraidNode* >( pos-1 , new_node2 ) );
    if( l1==l2 ) {
      to_check2.insert( pair< int , BraidNode* >( pos   , new_node  ) );
      break;
    }
    l1 = l3;
  }

  // B. Remove "boundary" crossings that formed a handle.
  to_check.erase( NODE( n1->weight , pos , n1 ) );
  to_check.erase( NODE( n2->weight , pos , n2 ) );
  if( lst ) {
    lst->push_back( removeNode( n1 , pos ) ); 
    lst->push_back( removeNode( n2 , pos ) );
  } else {
    removeNode( n1 , pos );
    removeNode( n2 , pos );
  }

  // Check if new handles were introduced.
  for( set< pair< int , BraidNode* > >::iterator c_it = to_check2.begin( ) ; c_it!=to_check2.end( ) ; ++c_it ) {
    int pos = (*c_it).first;
    BraidNode* bn = (*c_it).second;
    int weight = checkIfStartsLeftHandle( pos , bn );
    if( weight!=-1 ) {
      // remove old version
      to_check.erase( NODE( bn->weight , pos , bn ) );
      // add a new
      bn->weight=weight;
      to_check.insert( NODE( weight , pos , bn ) );
    }
  }
}



//---------------------------------------------------------------------------//
//---------------------------- removeRightHandles ---------------------------//
//---------------------------------------------------------------------------//


void LinkedBraidStructure::removeRightHandles( list< LinkedBraidStructureTransform >* result )
{
  // set of nodes to check (1st value = number of uppercrossings, 2nd value = position in the Linked Structure) 
  typedef triple< int , int , BraidNode* > NODE;
  set< NODE > to_check;
  
  // form the initial set of nodes to check (all nodes from the Linked Structure)
  for( int i=0 ; i<theIndex ; ++i ) {
    for( BraidNode* c = frontNodes[i] ; c!=backNodes[i] ; c = c->back ) {
      int weight = checkIfStartsRightHandle( i , c );
      if( weight!=-1 ) {
	c->weight=weight;
	to_check.insert( NODE( weight , i , c ) );
      }
    }
  }

  // check the nodes one by one starting from the rightmost
  for( ; !to_check.empty( ) ; ) {

    // choose the node and remove it from to_check
    NODE cur_node = *to_check.begin( );
    to_check.erase( to_check.begin( ) );
    
    // b. remove the handle
    int weight = checkIfStartsRightHandle( cur_node.second , cur_node.third );
    if( weight!=-1 )
      removeRightHandle( cur_node , to_check , result );
  }
}


//---------------------------------------------------------------------------//
//------------------------ checkIfStartsRightHandle -------------------------//
//---------------------------------------------------------------------------//


int LinkedBraidStructure::checkIfStartsRightHandle( int pos , BraidNode* bn )
{
  BraidNode* back = bn->back;
  // if there is no more x_i
  if( !back )
    return -1;
  
  // if there is an obstacle for handle
  if( bn->back_left )
    return -1;

  // if crossings have the same orientation
  if( bn->type==back->type )
    return -1;

  // compute the number of upper crossings
  BraidNode* l1 = bn->back_right;
  BraidNode* l2 = back->right;
  int counter=0;
  bool crossingType;
  for( ; l1 ; l1 = l1->back ) {
    if( counter==0 ) {
      crossingType = l1->type;
    } else {
      if( crossingType!=l1->type )
	return -1;
    }
    ++counter;
    if( l1==l2 )
      break;
  }
  
  return counter;
}


//---------------------------------------------------------------------------//
//----------------------------- removeRightHandle ---------------------------//
//---------------------------------------------------------------------------//


void LinkedBraidStructure::removeRightHandle( triple< int , int , BraidNode* > node , set< triple< int , int , BraidNode* > >& to_check , list< LinkedBraidStructureTransform >* lst )
{
  typedef triple< int , int , BraidNode* > NODE;
  
  int pos = node.second;
  BraidNode* n1 = node.third;
  BraidNode* n2 = n1->back;
  bool type = n1->type;
  
  BraidNode* r1 = n1->back_right;
  BraidNode* r2 = n2->right;

  // Removal of a handle can introduce new handles. Here we store some nodes to check. A few will be added later.
  set< pair< int , BraidNode* > > to_check2;
  if( n1->ahead ) to_check2.insert( pair< int , BraidNode* >( pos   , n1->ahead ) );
  if( n1->right ) to_check2.insert( pair< int , BraidNode* >( pos+1 , n1->right  ) );
  for( BraidNode* cn=n1 ; cn ; cn=cn->ahead ) {
    if( !cn->left )
      continue;
    to_check2.insert( pair< int , BraidNode* >( pos-1 , cn->left ) );
    break;
  }
  
  // B. process right nodes
  int counter=0;
  for( ; r1 ; ) {
    // ++counter;

    // cout << "Process right node" << endl;
    BraidNode* r3=r1->back;
    BraidNode* new_node  = insertBackLeft ( r1 , pos , r1->type );
    BraidNode* new_node2 = insertBackRight( new_node , pos+1 , type );

    if( lst ) {
      lst->push_back( make_AddTransform( new_node , pos ) );
      lst->push_back( make_AddTransform( new_node2 , pos+1 ) );
      lst->push_back( make_ChangeType( r1 , pos ) );
    }

    r1->type = !type;

    to_check2.insert( pair< int , BraidNode* >( pos+1 , new_node2 ) );

    if( r1==r2 ) {
      to_check2.insert( pair< int , BraidNode* >( pos   , new_node  ) );
      break;
    }
    r1 = r3;
  }

  // B. Remove "boundary" crossings that formed a handle.
  to_check.erase( NODE( n1->weight , pos , n1 ) );
  to_check.erase( NODE( n2->weight , pos , n2 ) );
  if( lst ) {
    lst->push_back( removeNode( n1 , pos ) );
    lst->push_back( removeNode( n2 , pos ) );
  } else {
    removeNode( n1 , pos );
    removeNode( n2 , pos );
  }

  // Check if new handles were introduced.
  for( set< pair< int , BraidNode* > >::iterator c_it = to_check2.begin( ) ; c_it!=to_check2.end( ) ; ++c_it ) {
    int pos = (*c_it).first;
    BraidNode* bn = (*c_it).second;
    int weight = checkIfStartsRightHandle( pos , bn );
    if( weight!=-1 ) {
      // remove old version
      to_check.erase( NODE( bn->weight , pos , bn ) );
      // add a new
      bn->weight=weight;
      to_check.insert( NODE( weight , pos , bn ) );
    }
  }
}


//---------------------------------------------------------------------------//
//-------------------------------- removeNode -------------------------------//
//---------------------------------------------------------------------------//


LinkedBraidStructureTransform LinkedBraidStructure::removeNode( BraidNode* bn , int pos )
{
  LinkedBraidStructureTransform result = make_EraseTransform( bn , pos );
  
  // A. process left
  if( bn->left )
    bn->left->back_right = bn->back_left ? NULL : bn->back;

  // B. process ahead
  if( bn->ahead ) {
    bn->ahead->back = bn->back;
    if( !bn->ahead->back_left )
      bn->ahead->back_left = bn->back_left;
    if( !bn->ahead->back_right )
      bn->ahead->back_right = bn->back_right;
  }

  // C. process right
  if( bn->right )
    bn->right->back_left = bn->back_right ? NULL : bn->back;
  
  // D. process back_left
  if( bn->back_left )
    bn->back_left->right = bn->left ? NULL : bn->ahead;

  // E. process back
  if( bn->back ) {
    bn->back->ahead = bn->ahead;
    if( !bn->back->left )
      bn->back->left = bn->left;
    if( !bn->back->right )
      bn->back->right = bn->right;
  }
  
  // F. process back_right
  if( bn->back_right )
    bn->back_right->left = bn->right ? NULL : bn->ahead;
  
  // G. process frontNodes and backNodes
  if( frontNodes[pos]==bn )
    frontNodes[pos] = bn->back;
  if( backNodes[pos]==bn )
    backNodes[pos] = bn->ahead;

  // H. finally delet the memory allocated for the node
  theNodes.erase( bn->theNumber );

  return result;
}


//---------------------------------------------------------------------------//
//------------------------------ insertBackLeft -----------------------------//
//---------------------------------------------------------------------------//


BraidNode* LinkedBraidStructure::insertBackLeft( BraidNode* bn , int pos , bool type )
{
  // A. determine the "main" nodes around the new node
  BraidNode* back = bn->back;

  BraidNode* left = NULL;
  for( BraidNode* c=bn ; c!=NULL && !left ; c=c->ahead )
    if( c->left )
      left = c->left;

  BraidNode* back_left = NULL;
  for( BraidNode* c=bn ; c!=NULL && !back_left ; c=c->back )
    if( c->back_left )
      back_left = c->back_left;

  BraidNode* left_back_left;
  if( left )
    left_back_left = left->back_left;
  else if( back_left ) {
    if( !back_left->left )
      left_back_left = NULL;
    else
      left_back_left = pos>0 ? frontNodes[pos-1] : NULL;
  } else {
    left_back_left = pos>0 ? frontNodes[pos-1] : NULL;
  }

  BraidNode* back_left_right;
  if( back_left )
    back_left_right = back_left->right!=bn ? bn->back : NULL;
  else
    back_left_right = back;

  // B. create a new node and link it to other nodes
  BraidNode& newNode = theNodes[maxNodeNumber] = BraidNode( maxNodeNumber , type , NULL , left , bn , left_back_left , back_left , back_left_right );
  maxNodeNumber++;
  
  bn->back_left = &newNode;

  if( left ) {
    left->back = &newNode;
    left->back_left = NULL;
  }
  if( back_left ) {
    back_left->ahead = &newNode;
    back_left->right = back_left->right==bn ? NULL : back_left->right;
  }
  if( left_back_left )
    left_back_left->right = &newNode;
  if( back_left_right )
    back_left_right->left = &newNode;

  if( !newNode.back )
    backNodes[pos] = &newNode;
  if( !newNode.ahead )
    frontNodes[pos] = &newNode;

  return &newNode;
}


//---------------------------------------------------------------------------//
//------------------------------ insertBackRight ----------------------------//
//---------------------------------------------------------------------------//


BraidNode* LinkedBraidStructure::insertBackRight( BraidNode* bn , int pos , bool type )
{
  // A. determine the "main" nodes around the new node
  BraidNode* back  = bn->back;
  
  BraidNode* right = NULL;
  for( BraidNode* c=bn ; c!=NULL && !right; c=c->ahead )
    if( c->right )
      right = c->right;

  BraidNode* back_right = NULL;
  for( BraidNode* c=bn ; c!=NULL && !back_right ; c=c->back )
    if( c->back_right )
      back_right = c->back_right;

  BraidNode* right_back_right;
  if( right )
    right_back_right = right->back_right;
  else if( back_right ) {
    if( !back_right->right )
      right_back_right = NULL;
    else
      right_back_right = pos<theIndex-1 ? frontNodes[pos+1] : NULL;
  } else {
    right_back_right = pos<theIndex-1 ? frontNodes[pos+1] : NULL;
  }
  
  BraidNode* back_right_left;
  if( back_right )
    back_right_left = back_right->left!=bn ? bn->back : NULL;
  else
    back_right_left = back;

  // B. create a new node and link it to other nodes
  BraidNode& newNode = theNodes[maxNodeNumber] = BraidNode( maxNodeNumber , type , bn , right , NULL , back_right_left , back_right , right_back_right );
  maxNodeNumber++;

  bn->back_right = &newNode;
  
  if( right ) {
    right->back = &newNode;
    right->back_right = NULL;
  }
  if( back_right ) {
    back_right->ahead = &newNode;
    back_right->left = back_right->left==bn ? NULL : back_right->left;
  }
  
  if( back_right_left )
    back_right_left->right = &newNode;
  if( right_back_right )
    right_back_right->left = &newNode;

  // update backNodes and frontNodes
  if( !newNode.back )
    backNodes[pos] = &newNode;
  if( !newNode.ahead )
    frontNodes[pos] = &newNode;

  return &newNode;
}


//---------------------------------------------------------------------------//
//---------------------------------- insert ---------------------------------//
//---------------------------------------------------------------------------//


BraidNode* LinkedBraidStructure::insert( const LinkedBraidStructureTransform& lbst )
{
  BraidNode& newNode = theNodes[lbst.theNumber] = BraidNode( lbst.theNumber , 
							     lbst.type , 
							     lbst.left==-1       ? NULL : &theNodes[lbst.left] ,
							     lbst.ahead==-1      ? NULL : &theNodes[lbst.ahead] ,
							     lbst.right==-1      ? NULL : &theNodes[lbst.right] ,
							     lbst.back_left==-1  ? NULL : &theNodes[lbst.back_left] ,
							     lbst.back==-1       ? NULL : &theNodes[lbst.back] ,
							     lbst.back_right==-1 ? NULL : &theNodes[lbst.back_right] );
  
  // A. determine the "main" nodes around the new node
  if( newNode.left ) newNode.left->back_right = &newNode;
  if( newNode.ahead ) {
    newNode.ahead->back = &newNode;
    if( !newNode.left  ) newNode.ahead->back_left  = NULL;
    if( !newNode.right ) newNode.ahead->back_right = NULL;
  }
  if( newNode.right ) newNode.right->back_left = &newNode;
  if( newNode.back_left ) newNode.back_left->right = &newNode;
  if( newNode.back ) {
    newNode.back->ahead = &newNode;
    if( !newNode.back_left  ) newNode.back->left  = NULL;
    if( !newNode.back_right ) newNode.back->right = NULL;
  }
  if( newNode.back_right ) newNode.back_right->left = &newNode;
  
  // B. update backNodes and frontNodes
  if( !newNode.back )
    backNodes[lbst.thePosition] = &newNode;
  if( !newNode.ahead )
    frontNodes[lbst.thePosition] = &newNode;

  return &newNode;
}


//---------------------------------------------------------------------------//
//-------------------------------- push_front -------------------------------//
//---------------------------------------------------------------------------//


void LinkedBraidStructure::processTree( int al , BraidNode* node , Word& result ) const
{
  if( node->back_left  && !node->back_left ->link ) processTree( al-1 , node->back_left  , result );
  if( node->back_right && !node->back_right->link ) processTree( al+1 , node->back_right , result );
  if( node->back       && !node->back->link       ) processTree( al   , node->back       , result );
  result.push_front( node->type ? al : -al );
  node->link = node;
}


Word LinkedBraidStructure::translateIntoWord( ) const
{
  Word result;

  for( int i=0 ; i<theIndex ; ++i )
    if( frontNodes[i] && !frontNodes[i]->link )
      processTree( i+1 , frontNodes[i] , result );

  clearLinks( );

  return result;
}


void LinkedBraidStructure::clearLinks( ) const
{
  map< int , BraidNode >::const_iterator n_it = theNodes.begin( );
  for( ;  n_it!=theNodes.end( ) ; ++n_it )
    (*n_it).second.link = NULL;
}


//---------------------------------------------------------------------------//
//---------------------------------- clear ----------------------------------//
//---------------------------------------------------------------------------//


void LinkedBraidStructure::clear( )
{
  theNodes.clear( );
  frontNodes = vector< BraidNode* >( theIndex , (BraidNode*)NULL );
  backNodes  = vector< BraidNode* >( theIndex , (BraidNode*)NULL );
}


//---------------------------------------------------------------------------//
//---------------------------------- undo -----------------------------------//
//---------------------------------------------------------------------------//


void LinkedBraidStructure::undo( const list< LinkedBraidStructureTransform >& lbst_seq )
{
  list< LinkedBraidStructureTransform >::const_iterator s_it=lbst_seq.end( );
  for( ; s_it!=lbst_seq.begin( ) ; )
    undo( *--s_it );
}


void LinkedBraidStructure::undo( const LinkedBraidStructureTransform& lbst )
{
  if( lbst.theTransform==LinkedBraidStructureTransform::ADDED ) {
    removeNode( &theNodes[lbst.theNumber] , lbst.thePosition );
  } else if( lbst.theTransform==LinkedBraidStructureTransform::ERASED ) {
    insert( lbst );
  } else if( lbst.theTransform==LinkedBraidStructureTransform::CHANGE_TYPE ) {
    theNodes[lbst.theNumber].type = lbst.type;
  }
}
