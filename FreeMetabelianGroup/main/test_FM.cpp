// Copyright (C) 2007 Alexander Ushakov
// Contents: Test for class FreeMetabelianGroupAlgorithms
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//


#include "Word.h"
#include "FreeMetabelianGroupAlgorithms.h"

#include "fstream"
#include "strstream"
#include "string"
using namespace std;


vector< int > getTail( int N , const Word& w );
map< vector< int > , int > getEdgeMap( int N , const Word w );
map< vector< int > , int > dropTrivialEdgesInMap( const map< vector< int > , int >& EM );


//---------------------------------------------------------------------------//
//------------------------------------ main ---------------------------------//
//---------------------------------------------------------------------------//


string arrow_ps_code( int x , int y , int BX , int BY , int dir , int arity , int ar1 , int ar2 )
{
  char str[1024];
  ostrstream OS( str, 1024 );
  
  int px = x*BX;
  int py = y*BY;

  int ad = abs(dir);
  int x2 = dir==1 ? x+1 : x;
  int y2 = dir==2 ? y+1 : y;
  int px2 = BX*x2;
  int py2 = BY*y2;
  // cout << x << " , " << y << " -> " << ad << endl;
  
  OS << px << " " << py << " moveto" << endl;
  OS << px2 << " " << py2 << " lineto" << endl << endl;
  if( dir==1 ) {
    if( arity<0 ) {
      OS << px << " " << py << " moveto" << endl;
      OS << px+ar1 << " " << py-ar2 << " lineto" << endl << endl;
      OS << px << " " << py << " moveto" << endl;
      OS << px+ar1 << " " << py+ar2 << " lineto" << endl << endl;
    } else {
      OS << px2 << " " << py2 << " moveto" << endl;
      OS << px2-ar1 << " " << py2-ar2 << " lineto" << endl << endl;
      OS << px2 << " " << py2 << " moveto" << endl;
      OS << px2-ar1 << " " << py2+ar2 << " lineto" << endl << endl;
    }
  } else {
    if( arity<0 ) {
      OS << px << " " << py << " moveto" << endl;
      OS << px-ar2 << " " << py+ar2 << " lineto" << endl << endl;
      OS << px << " " << py << " moveto" << endl;
      OS << px+ar2 << " " << py+ar2 << " lineto" << endl << endl;
    } else {
      OS << px2 << " " << py2 << " moveto" << endl;
      OS << px2-ar2 << " " << py2-ar1 << " lineto" << endl << endl;
      OS << px2 << " " << py2 << " moveto" << endl;
      OS << px2+ar2 << " " << py2-ar1 << " lineto" << endl << endl;
    }
  }
  
  OS << ends;

  return str;
}


//---------------------------------------------------------------------------//
//------------------------------------ main ---------------------------------//
//---------------------------------------------------------------------------//


void visualize( const string& FN , const Word& w )
{
  // N = 2;
  
  ofstream OF( FN.c_str( ) );
  OF << "%!PS" << endl << endl;
  OF << ".5 setlinewidth" << endl << endl;

  // move the origin to the center
  int LX = 590;
  int LY = 840;
  int CX = LX/2;
  int CY = LY/2;
  int BX = 5;
  int BY = 5;
  
  OF << CX << " " << CY << " translate" << endl;
  int ax_x = 200;
  int ax_y = 200;
  int ar1 = 10;
  int ar2 = 4;
  OF << "newpath" << endl;
  OF << -ax_x << " 0 moveto" << endl;
  OF <<  ax_x << " 0 lineto" << endl << endl;
  OF <<  ax_x << " 0 moveto" << endl;
  OF <<  ax_x-ar1 << " " << -ar2 << " lineto" << endl << endl;
  OF <<  ax_x << " 0 moveto" << endl;
  OF <<  ax_x-ar1 << " " << ar2 << " lineto" << endl << endl;

  OF << "0 " << -ax_y << " moveto" << endl;
  OF << "0 " <<  ax_y << " lineto" << endl << endl;
  OF << "0 " <<  ax_y << " moveto" << endl;
  OF <<  -ar2 << " " << ax_x-ar1 << " lineto" << endl << endl;
  OF << "0 " <<  ax_y << " moveto" << endl;
  OF <<   ar2 << " " << ax_x-ar1 << " lineto" << endl << endl;

  for( int i=-6 ; i<=6 ; ++i ) {
    OF << -ax_x << " " << 5*i*BX << " moveto" << endl;
    OF <<  ax_x << " " << 5*i*BX << " lineto" << endl;
    OF << 5*i*BX << " " << -ax_y << " moveto" << endl;
    OF << 5*i*BX << " " <<  ax_y << " lineto" << endl;
  }

  OF << "/Times-Roman findfont" << endl;
  OF << "12 scalefont" << endl;
  OF << "setfont" << endl;
  OF <<  ax_x+5 << " 0 moveto" << endl;
  OF << "(x) show" << endl;
  OF <<  5 << " " << ax_y << " moveto" << endl;
  OF << "(y) show" << endl;

  OF << "0.5 setgray" << endl;
  OF << "stroke" << endl;
  
  // the origin
  OF << "newpath" << endl;
  OF << "0 0 2 0 360 arc" << endl;
  OF << "closepath fill" << endl;

  // the end point
  vector< int > T = getTail( 2 , w );
  OF << "newpath" << endl;
  OF << BX*T[0] << " " << BY*T[1] << " 2 0 360 arc" << endl;
  OF << "closepath fill" << endl;
  
  map< vector< int > , int > EM = dropTrivialEdgesInMap( getEdgeMap( 2 , w ) );
  for( map< vector< int > , int >::const_iterator E_it=EM.begin( ) ; E_it!=EM.end( ) ; ++E_it ) {

    pair< vector< int > , int > C = *E_it;
    int x = C.first[0];
    int y = C.first[1];
    OF << arrow_ps_code( x , y , BX , BY , C.first[2] , C.second , 1 , 1 );
  }

  OF << "0 setgray" << endl;
  OF << "stroke" << endl;
  OF << "showpage" << endl;
}


//---------------------------------------------------------------------------//
//------------------------------------ main ---------------------------------//
//---------------------------------------------------------------------------//


int main( )
{
  Word w1 = Word::randomWord( 2 , 1000 );
  // Word w2 = Word::randomWord( 2 , 1000 );
  // Word w2 = Word(1);
  Word w2 = Word(-1)*Word(-2)*Word(1)*Word(2);
  // visualize( "w.ps" , w1 );
  visualize( "w1.ps" , w1 );
  visualize( "w2.ps" , w2 );
  visualize( "w.ps" , w1*w2*-w1*-w2 );
  return 0;


  
  int N = 3;
  
  for( int i=0 ; i<1000 ; ++i ) {
    
    cout << endl << "Exp N " << i << endl;
    Word w1 = Word::randomWord( N , 159 );
    // Word w1 = Word( 1 ) * Word( 2 ) * Word( -1 ) * Word( -2 );
    Word c = Word::randomWord( N , 257 );
    Word w2 = -c* w1 *c;
    
    pair< bool , Word > T = FreeMetabelianGroupAlgorithms::conjugate( N , w1 , w2 );
    if( T.first==0 )
      exit(0);

    if( FreeMetabelianGroupAlgorithms::trivial( N , -T.second*w1*T.second * -w2 ) )
      cout << "Test ok" << endl;
    else {
      cout << "Test failure" << endl;
      exit(0);
    }
    
  }

  return 0;
}
