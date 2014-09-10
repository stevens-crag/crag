
#ifndef _LinkedBraidStructure_h_
#define _LinkedBraidStructure_h_


#include "set"
#include "map"
#include "list"
#include "vector"
using namespace std;

#include "tuples.h"


//---------------------------------------------------------------------------//
//------------------------------ BraidNode ----------------------------------//
//---------------------------------------------------------------------------//


struct BraidNode
{

  ///////////////////////////////////////////////////////// 
  //                                                     //
  //  Constructors:                                      //
  //                                                     //
  /////////////////////////////////////////////////////////
    
  public:

  BraidNode( int num=0, bool tp=true, BraidNode* l=NULL, BraidNode* a=NULL, BraidNode* r=NULL, BraidNode* bl=NULL, BraidNode* b=NULL, BraidNode* br=NULL) : 
    left( l ), ahead( a ), right( r ),
    back_left( bl ), back( b ), back_right( br ) ,
    type( tp ) ,
    link( NULL ),
    weight(0),
    theNumber( num )
  { }

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Data members:                                      //
  //                                                     //
  /////////////////////////////////////////////////////////
  
  public:
  
  BraidNode* left;
  BraidNode* ahead;
  BraidNode* right;

  BraidNode* back_left;
  BraidNode* back;
  BraidNode* back_right;
  
  bool type;
  // shows whether a crossing positive or negative
  
  mutable BraidNode* link;
  // auxiliary member, used in LinkedBraidStructure copy constructor and in translateIntoWord
  
  mutable int weight;
  // auxiliary member, used in remove Handle functions
  
  int theNumber;
  // a unique number associated with the node in a LBS
};

ostream& operator << ( ostream& os , const BraidNode& bn );



//---------------------------------------------------------------------------//
//---------------------- LinkedBraidStructureTransform ----------------------//
//---------------------------------------------------------------------------//


struct LinkedBraidStructureTransform
{
  enum TRANSFORM{ ERASED , ADDED , CHANGE_TYPE };

  LinkedBraidStructureTransform( int n , int p , TRANSFORM tr , bool t=false , int l=0 , int a=0 , int r=0 , int bl=0 , int b=0 , int br=0 ) :
    theTransform(tr), left(l), ahead(a), right(r), back_left(bl), back(b), back_right(br), type(t), theNumber( n ), thePosition(p)
  { }

  TRANSFORM theTransform;

  int left;
  int ahead;
  int right;

  int back_left;
  int back;
  int back_right;
  
  bool type;
  int theNumber;
  int thePosition;
};


//---------------------------------------------------------------------------//
//------------------------- LinkedBraidStructure ----------------------------//
//---------------------------------------------------------------------------//


class LinkedBraidStructure
{

  ///////////////////////////////////////////////////////// 
  //                                                     //
  //  Constructors:                                      //
  //                                                     //
  /////////////////////////////////////////////////////////
    
 public:

  LinkedBraidStructure( int N );
  LinkedBraidStructure( int N , const Word& w );
  LinkedBraidStructure( const LinkedBraidStructure& LBS );

  LinkedBraidStructure& operator= ( const LinkedBraidStructure& LBS );

  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Accessors:                                         //
  //                                                     //
  /////////////////////////////////////////////////////////

 public:

  bool operator < ( const LinkedBraidStructure& lbs ) const;
  // function check if the braid-word represented by *this is shortlex smaller than the one given by lbs
  // usually the computation of the words is not required, so the function is often fast

  int size( ) const { return theNodes.size( ); }
  void clear( );

  LinkedBraidStructureTransform push_back ( int g );
  LinkedBraidStructureTransform push_front( int g );
  
  void removeLeftHandles ( list< LinkedBraidStructureTransform >* result = NULL );
  void removeRightHandles( list< LinkedBraidStructureTransform >* result = NULL );
  
  Word translateIntoWord( ) const;

  void undo( const list< LinkedBraidStructureTransform >& lbst_seq );
  void undo( const LinkedBraidStructureTransform& lbst );

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Internal functions:                                //
  //                                                     //
  /////////////////////////////////////////////////////////

 private:

  LinkedBraidStructureTransform make_EraseTransform( BraidNode* bn , int pos ) const;
  LinkedBraidStructureTransform make_AddTransform( BraidNode* bn , int pos ) const;
  LinkedBraidStructureTransform make_ChangeType( BraidNode* bn , int pos ) const;

  int checkIfStartsLeftHandle ( int pos , BraidNode* bn );
  int checkIfStartsRightHandle( int pos , BraidNode* bn );
  void removeLeftHandle( triple< int , int , BraidNode* > node , set< triple< int , int , BraidNode* > >& to_check , list< LinkedBraidStructureTransform >* lst );
  void removeRightHandle(triple< int , int , BraidNode* > node , set< triple< int , int , BraidNode* > >& to_check , list< LinkedBraidStructureTransform >* lst );

  LinkedBraidStructureTransform removeNode( BraidNode* bn , int pos );
  BraidNode* insertBackRight( BraidNode* bn , int pos , bool type );
  BraidNode* insertBackLeft( BraidNode* bn , int pos , bool type );
  BraidNode* insert( const LinkedBraidStructureTransform& lbst );
  
  void processTree( int al , BraidNode* node , Word& result ) const;

  void clearLinks( ) const;

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Data members:                                      //
  //                                                     //
  /////////////////////////////////////////////////////////

 private:
  
  //! Number of generators!!!
  int theIndex;
  vector< BraidNode* > frontNodes;
  vector< BraidNode* >  backNodes;
  map< int , BraidNode > theNodes;
  
  int maxNodeNumber;
};


#endif

