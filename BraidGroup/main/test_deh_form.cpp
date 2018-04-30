
#include "braid_group.h"
#include "LinkedBraidStructure.h"
#include "DehornoyForm.h"
#include "ShortBraidForm.h"
#include "ThRightNormalForm.h"

#include <time.h>
#include <iostream>
using namespace std;


//---------------------------------------------------------------------------//
//---------------------------------- main -----------------------------------//
//---------------------------------------------------------------------------//


void test_LBS_history( )
{
  int N = 4;

  LinkedBraidStructure lbs( N );
  list< LinkedBraidStructureTransform > history;
  for( int i=0 ; i<10 ; ++i ) {
    if( rand()%2==0 ) {
      history.push_back( lbs.push_back( 1+rand()%(N-1) ) );
    } else {
      history.push_back( lbs.push_front( 1+rand()%(N-1) ) );
    }
    Word w = lbs.translateIntoWord( );
    cout << w << endl;
  }

  list< LinkedBraidStructureTransform >::iterator h_it = history.end( );
  for( ; h_it!=history.begin( ) ; ) {
    --h_it;
    lbs.undo( *h_it );
    Word w1 = lbs.translateIntoWord( );
    cout << w1 << endl;
  }
}

//---------------------------------------------------------------------------//
//---------------------------------- main -----------------------------------//
//---------------------------------------------------------------------------//


void test_LBS_history2( )
{
  int N = 60;

  for( int i=0 ; i<100 ; ++i ) {

    cout << "i=" << i << endl;
    // Word w = Word(-1) * Word(-2) * Word(1);
    Word w = Word::randomWord( N-1 , 30000 );
    LinkedBraidStructure lbs( N , w );
    w = lbs.translateIntoWord( );
    // cout << "==============================" << endl;
    // cout << "w  = " << w << endl;

    int t1 = time(0);
    list< LinkedBraidStructureTransform > history;
    lbs.removeRightHandles( &history );
    Word w1 = lbs.translateIntoWord( );
    // cout << "w1 = " << w1 << endl;
    
    int t2 = time(0);
    lbs.undo( history );
    Word w2 = lbs.translateIntoWord( );
    // cout << "w2 = " << w2 << endl;
    int t3 = time(0);

    if( w!=w2 ) {
      cout << "Incorrect" << endl;
      exit(1);
    } else {
      cout << "  ok, " << t2-t1 << ", " << t3-t2 << endl;
    }
  }
}


//---------------------------------------------------------------------------//
//---------------------------------- main -----------------------------------//
//---------------------------------------------------------------------------//


void test_LBS_history3( )
{
  int N = 10;
  Word w = Word::randomWord( N-1 , 1000 );
  LinkedBraidStructure lbs( N , w );
  lbs.removeRightHandles( );
  cout << lbs.size( ) << endl;

  for( int i=0 ; i<100 ; ++i ) {
    Word c = Word::randomWord( N-1 , 25 );
    list< LinkedBraidStructureTransform > history;
    for( Word::const_iterator c_it=c.begin( ) ; c_it!=c.end( ) ; ++c_it ) {
      // multiply lbs on the right and remember the transformation
      history.push_back( lbs.push_back(   *c_it ) );
      // multiply lbs on the left and remember  the transformation
      history.push_back( lbs.push_front( -*c_it ) );
    }
    list< LinkedBraidStructureTransform > history2;
    lbs.removeRightHandles( &history2 );
    cout << lbs.size( ) << endl;
    lbs.undo( history2 );
    lbs.undo( history );
  }
}


//---------------------------------------------------------------------------//
//---------------------------------- main -----------------------------------//
//---------------------------------------------------------------------------//


void test_short_forms( )
{
  int N = 100;
  int L = 8000;
  
  for( int i=0 ; i<100 ; ++i ) {
    
    cout << i << ": " << endl;
    Word w = Word::randomWord( N-1 , L );
    int t = time(0);
    Word w1 = shortBraidForm( N , w );
    // Word w1 = shortenBraid( N , w );
    cout << "     " << w.length( ) << " -> " << w1.length( ) << endl;
    cout << "          time = " << time(0)-t << endl;
  }
}

void test_short_forms_correct( )
{
  int N = 10;
  int L = 100;
  
  for( int i=0 ; i<100 ; ++i ) {
    
    cout << i << ": " << endl;
    Word w = Word::randomWord( N-1 , L );
    Word w1 = shortBraidForm( N , w );

    crag::braidgroup::BraidGroup B( N );
    ThRightNormalForm nf( B , w*-w1 );
    if( nf.isTrivial( ) )
      cout << "Correct: " << w.length( ) << ", " << w1.length( ) << endl;
    else {
      cout << "Wrong" << endl;
      // cout << w << endl;
      exit(1);
    }
    
  }
}


//---------------------------------------------------------------------------//
//---------------------------------- main -----------------------------------//
//---------------------------------------------------------------------------//


void test_deh_form_efficiency( )
{
  int N = 100;

  int exp = 20;
  int total_length = 0;
  for( int i=0 ; i<exp ; ++i ) {
    cout << i << ": " << endl;
    Word w = Word::randomWord( N-1 , 40000 );
    cout << "b" << endl;
    DehornoyForm DF( N , w );
    cout << "c" << endl;
    Word w1 = DF.getDehornoyForm( );
    cout << "d" << endl;
    cout << w1.length( ) << endl;
    total_length += w1.length( );
  }
  cout << "avg = " << total_length/exp << endl;
}


//---------------------------------------------------------------------------//
//---------------------------------- main -----------------------------------//
//---------------------------------------------------------------------------//


void test_nf_shorten_efficiency( )
{
  int N = 100;
  int L = 2000;
  
  int exp = 20;
  int avg = 0;
  int total_length = 0;
  for( int i=0 ; i<exp ; ++i ) {
    cout << i << ": " << endl;
    cout << "a" << endl;
    Word w = Word::randomWord( N-1 , L );
    cout << "b" << endl;
    crag::braidgroup::BraidGroup B( N );
    ThRightNormalForm nf( B , w );
    cout << "c" << endl;
    DehornoyForm DF( N , nf.getShortWord( ) );
    cout << "d" << endl;
    Word w1 = DF.getDehornoyForm( );
    cout << "e" << endl;
    cout << L << " -> " << w1.length() << endl;
    avg += w1.length( );
  }
  cout << "avg = " << total_length/exp << endl;
}


//---------------------------------------------------------------------------//
//---------------------------------- main -----------------------------------//
//---------------------------------------------------------------------------//


void test_deh_form_correct( )
{
  int N = 10;

  for( int i=0 ; i<100 ; ++i ) {

    cout << i << ": " << endl;
    Word w = Word::randomWord( N-1 , 1000 );
    DehornoyForm DF( N , w );
    Word w1 = DF.getDehornoyForm( );

    crag::braidgroup::BraidGroup B( N );
    ThRightNormalForm nf( B , w*-w1 );
    if( nf.isTrivial( ) )
      cout << "Correct: " << w.length( ) << ", " << w1.length( ) << endl;
    else {
      cout << "Wrong" << endl;
      cout << w << endl;
      exit(1);
    }
  }
}


//---------------------------------------------------------------------------//
//---------------------------------- main -----------------------------------//
//---------------------------------------------------------------------------//

void test_deh_form_correct2( )
{
  int N = 5;

  for( int i=0 ; i<100 ; ++i ) {

    cout << i << ": " << endl;
    Word w = Word::randomWord( N-1 , 100 );
    DehornoyForm DF( N , w );
    DehornoyForm DF2( DF );

    Word w1 = DF.getDehornoyForm( );
    Word w2 = DF2.getDehornoyForm( );
    // cout << w1 << endl;
    // cout << w2 << endl;

    if( w1==w2 )
      cout << "Correct: " << w.length( ) << ", " << w1.length( ) << endl;
    else {
      cout << "Wrong" << endl;
      cout << w << endl;
      exit(1);
    }
  }
}


//---------------------------------------------------------------------------//
//---------------------------------- main -----------------------------------//
//---------------------------------------------------------------------------//

void test_deh_form_correct3( )
{
  cout << "test_deh_form_correct3" << endl;
  int N = 20;
  
  for( int i=0 ; i<100 ; ++i ) {

  // cout << i << ": " << endl;
    Word w = Word::randomWord( N-1 , 1000 );
    // Word x = Word(1), y = Word(2);
    // Word w = y*-x*-x*-y*-x*y;
    // cout << w << endl;
    DehornoyForm DF( N , w );
    Word w1 = DF.getDehornoyForm( );
    // cout << "===" << endl;
    DehornoyForm DF2( N , w1 );
    Word w2 = DF2.getDehornoyForm( );
    // cout << w1 << endl;
    // cout << w2 << endl;

    if( w1==w2 )
      cout << "Correct: " << w.length( ) << ", " << w1.length( ) << endl;
    else {
      cout << "Wrong: " << w1.length( ) << ", " << w2.length( ) << endl;
      cout << w << endl;
      cout << w1 << endl;
      cout << w2 << endl;
      exit(1);
    }
  }
}


//---------------------------------------------------------------------------//
//---------------------------------- main -----------------------------------//
//---------------------------------------------------------------------------//


void PB3_relations_check( )
{
  Word a(1);
  Word b(2);

  Word a12 = a^2;
  Word a23 = b^2;
  Word a13 = b*(a^2)*-b;
  
  Word t = a12*a13*a23;

  DehornoyForm DF1( 3 , a23*t*-a23*-t );
  DehornoyForm DF2( 3 , a13*t*-a13*-t );
  cout << "L1 = " << DF1.getDehornoyForm( ).length() << endl;
  cout << "L2 = " << DF2.getDehornoyForm( ).length() << endl;
}


//---------------------------------------------------------------------------//
//---------------------------------- main -----------------------------------//
//---------------------------------------------------------------------------//


int main( )
{
  srand( time(0) );

  PB3_relations_check( );
  // test_short_forms_correct( );
  // test_LBS_history3( );
  // test_deh_form_correct( );
  // test_short_forms( );
  // test_deh_form_efficiency( );
  // test_nf_shorten_efficiency( );
  return 0;

  int N = 10;

  for( int i=0 ; i<100 ; ++i ) {

    cout << i << ": " << endl;
    Word w = Word::randomWord( N-1 , 100 );
    // Word w = Word(-1) * Word(2) * Word(2);

    
    DehornoyForm DF( N , w );
    Word w1 = DF.getDehornoyForm( );
    // cout << "------------------------" << endl;
    // cout << w << endl;
    // cout << w1 << endl;
    cout << "     " << w.length( ) << " -> " << w1.length( ) << endl;

    // cout << w << endl;
    // cout << w1 << endl;
    // cout << w2 << endl;


    crag::braidgroup::BraidGroup B( N );
    ThRightNormalForm nf( B , w*-w1 );
    if( nf.isTrivial( ) )
      cout << "Correct: " << w.length( ) << ", " << w1.length( ) << endl;
    else {
      cout << "Wrong" << endl;
      cout << w << endl;
      exit(1);
    }
    /*
    DehornoyForm DF2( N , w1 );
    Word w2 = DF2.getDehornoyForm( );
    if( w1!=w2 ) {
      cout << "   Shit" << endl;
      exit(2);
    } else
      cout << "   ok" << endl;
    */

  }

  return 0;
}
