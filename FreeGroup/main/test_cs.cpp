// Copyright (C) 2007 Alexander Ushakov
// Contents: Test for class StraightLineProgramWord
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//


#include "StraightLineProgramWord.h"
#include "WhiteheadAutoSet.h"
#include <time.h>

#include "fstream"
using namespace std;


//---------------------------------------------------------------------------//
//------------------------------------ main ---------------------------------//
//---------------------------------------------------------------------------//

void test3( )
{
  ofstream of( "out.txt" );
  NielsenAutoSet NAS( 2 );
  const SetOfMaps& NielsSet = NAS.getSet();
  vector< Map > maps;
  of << "==================" << endl;
  for( SetOfMaps::const_iterator ns_it=NielsSet.begin( ) ; ns_it!=NielsSet.end( ) ; ++ns_it ) {
    of << maps.size() << ": " << *ns_it << endl;
    maps.push_back( *ns_it );
  }
  of << "==================" << endl;
  
  int aut = 20;
  int exp = 10000;

  for( int e=0 ; e<exp ; ++e ) {

    of << "Exp N " << e << endl;
    vector< Map > T1;
    for( int i=0 ; i<aut ; ++i ) {
      int m = rand()%maps.size( );
      of << m << ", ";
      T1.push_back( maps[m] );
    }
    of << endl;

    /*
    int M[] = { 5,0,4 };
    for( int i=0 ; i<sizeof(M)/sizeof(int) ; ++i ) {
      int m = M[i];
      of << m << ", ";
      T1.push_back( maps[m] );
    }
    of << endl;
    */

    StraightLineProgramWord CS1( 1 , T1.begin( ) , T1.end( ) );
    // CS1.simplify( );
    cout << CS1 << endl;

    
    CS1.reduce( );
    if( CS1.length( )!=CS1.getWord( ).length( ) ) {
      cout << "w = " << CS1.getWord( ) << endl;
      cout << "L = " << CS1.length( ) << endl;
      cout << "|w| = " << CS1.getWord( ).length( ) << endl;
      cout << "Failure" << endl;

      cout << CS1 << endl;

      exit(1);
    } else {
      cout << "Success" << endl;
    }
  }
}


//---------------------------------------------------------------------------//
//------------------------------------ main ---------------------------------//
//---------------------------------------------------------------------------//


void test2( )
{
  NielsenAutoSet NAS( 2 );
  const SetOfMaps& NielsSet = NAS.getSet();
  vector< Map > maps;
  cout << "==================" << endl;
  for( SetOfMaps::const_iterator ns_it=NielsSet.begin( ) ; ns_it!=NielsSet.end( ) ; ++ns_it ) {
    cout << maps.size() << ": " << *ns_it << endl;
    maps.push_back( *ns_it );
  }
  cout << "==================" << endl;
  
  int aut = 15;
  int exp = 10;

  for( int e=0 ; e<exp ; ++e ) {

    vector< Map > T1;
    for( int i=0 ; i<aut ; ++i )
      T1.push_back( maps[rand()%maps.size( )] );
    StraightLineProgramWord CS1( 1 , T1.begin( ) , T1.end( ) );
    
    vector< Map > T2;
    for( int i=0 ; i<aut ; ++i )
      T2.push_back( maps[rand()%maps.size( )] );
    StraightLineProgramWord CS2( 1 , T2.begin( ) , T2.end( ) );
    
    StraightLineProgramWord CS = CS1*CS2;
    if( CS.getWord()!=CS1.getWord()*CS2.getWord() ) {
      cout << "Failure" << endl;
      cout << "=========================" << endl;
      cout << CS1 << endl;
      cout << "=========================" << endl;
      cout << CS2 << endl;
      cout << "=========================" << endl;
      cout << CS << endl;
      exit(1);
    }
  }
}


//---------------------------------------------------------------------------//
//------------------------------------ main ---------------------------------//
//---------------------------------------------------------------------------//


int main( )
{
  srand(time(0));
  // test2( );
  test3( );

  return 0;
}
