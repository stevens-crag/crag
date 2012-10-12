// Copyright (C) 2007 Alexander Ushakov
// Contents: Test of class ThompsonGroupFNormalForm
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//


#include "ThompsonGroupFNormalForm.h"


void test_wp( )
{
  int E=100;
  for( int e=0 ; e<E ; ++e ) {
    Word w1 = Word::randomWord(2,100);
    Word w2 = Word::randomWord(2,100);
    
    ThompsonGroupFNormalForm nf1(w1);
    ThompsonGroupFNormalForm nf2(w2);
    ThompsonGroupFNormalForm nf (w1*w2);
    ThompsonGroupFNormalForm nf3 = nf1*nf2;
    if( nf==nf3 ) {
      cout << "#" << e << " success" << endl;
    } else {
      cout << "#" << e << " failure" << endl;
      break;
    }
  }
}


//---------------------------------------------------------------------------//
//--------------------------------- main ------------------------------------//
//---------------------------------------------------------------------------//


void test_wp2( )
{
  int E=100;
  for( int e=0 ; e<E ; ++e ) {
    
    Word w = Word::randomWord(2,100);
    ThompsonGroupFNormalForm nf( w);
    ThompsonGroupFNormalForm nf2(-w);
    ThompsonGroupFNormalForm nf3 = -nf;
    if( nf2==nf3 ) {
      cout << "#" << e << " success" << endl;
    } else {
      cout << "#" << e << " failure" << endl;
      break;
    }
  }
}


//---------------------------------------------------------------------------//
//--------------------------------- main ------------------------------------//
//---------------------------------------------------------------------------//


void exp_commuting_pairs( )
{
  int S = 0;
  int L = 50;
  int E = 100000;
  for( int e=0 ; e<E ; ++e ) {
    
    ThompsonGroupFNormalForm nf1( Word::randomWord(2,L) );
    ThompsonGroupFNormalForm nf2( Word::randomWord(2,L) );
    if( nf1*nf2==nf2*nf1 ) {
      cout << "#" << e << " yes" << endl;
      ++S;
    } else {
      cout << "#" << e << " no" << endl;
    }
  }
  cout << "S = " << S << endl;
}


int main( )
{
  // test_wp( );
  // test_wp2( );

  exp_commuting_pairs( );

  return 0;
}
