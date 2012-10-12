

#include "Word.h"
#include "TheGrigorchukGroupAlgorithms.h"

#include "map"
#include "vector"
#include "fstream"
using namespace std;


#include "GraphType.h"
#include "GraphConcept.h"

using namespace Graphs;


//---------------------------------------------------------------------------//
//----------------------- TheGrigorchukGroupAlgorithms ----------------------//
//---------------------------------------------------------------------------//

void act( int sz , vector< int >::iterator p_it , int x )
{
  int half = sz >> 1;
  
  if( x==1 )
    for( int i=0 ; i<half ; ++i )
      swap( *(p_it+i) , *(p_it+half+i) );
  
  if( sz<=2 ) return;
  
  if( x==2 ) {
    act( half , p_it , 1 );
    act( half , p_it+half , 3 );
  }
  if( x==3 ) {
    act( half , p_it , 1 );
    act( half , p_it+half , 4 );
  }
  if( x==4 )
    act( half , p_it+half , 2 );
}


void constructLevelGraph( )
{
  int L = 3;
  int B = 1 << L;
  cout << "B = " << B << endl;
  vector< int > P(B);
  for( int i=0 ; i<B ; ++i )
    P[i] = i;
  
  map< vector< int > , int > new_pt;
  map< vector< int > , int > old_pt;

  Graph G;
  new_pt[P] = G.newVertex( );

  while( !new_pt.empty( ) ) {
    
    pair< vector< int > , int > pr = *new_pt.begin( );
    new_pt.erase( new_pt.begin( ) );
    
    vector< int > P = pr.first;
    int num = pr.second;
    old_pt[P] = num;

    
    for( int action=1 ; action<=4 ; ++action ) {

      vector< int > P2 = P;
      act( B , P2.begin( ) , action );

     
      map< vector< int > , int >::iterator f = new_pt.find( P2 );
      if( f!=new_pt.end( ) ) {
	int num2 = (*f).second;
	G.newEdge( num , GraphEdge( num2 ) );
	continue;
      }

      f = old_pt.find( P2 );
      if( f!=old_pt.end( ) ) {
	int num2 = (*f).second;
	G.newEdge( num , GraphEdge( num2 ) );
	continue;
      }
      
      new_pt[P2] = G.newVertex( );
      cout << new_pt[P2] << "   ";
      for( int i=0 ; i<B ; ++i ) cout << P2[i] << " "; cout << endl;
    }
  }
  ofstream OF( "g.txt" );
  OF << graphviz_format( G ) << endl;
}


//---------------------------------------------------------------------------//
//----------------------- TheGrigorchukGroupAlgorithms ----------------------//
//---------------------------------------------------------------------------//


void test_lifts( )
{
  // Generators of K:
  Word G[3];
  G[0] = (Word(1)*Word(2)).power(2);
  G[1] = (Word(2)*Word(1)*Word(4)*Word(1)).power(2);
  G[2] = (Word(1)*Word(2)*Word(1)*Word(4)).power(2);

  typedef pair< int , int > PII;
  map< PII , int > preimages;
  
  
  // A. Consider 8 of 16 representatives belonging to ST(1)
  for( int b=0 ; b<=1 ; ++b ) {
    for( int a=0 ; a<=3 ; ++a ) {
      for( int d=0 ; d<=1 ; ++d ) {

	Word representative = TheGrigorchukGroupAlgorithms::reduce( Word(2).power(b) * (Word(4)*Word(1)).power(a) * Word(4).power(d) );
	if( a%2==0 ) {
	  cout << "==========================================" << endl;
	  pair< Word , Word > S = TheGrigorchukGroupAlgorithms::split( representative );
	  int n  = TheGrigorchukGroupAlgorithms::cosetRepresentativeKSbgp( representative );
	  int n1 = TheGrigorchukGroupAlgorithms::cosetRepresentativeKSbgp( S.first );
	  int n2 = TheGrigorchukGroupAlgorithms::cosetRepresentativeKSbgp( S.second );
	  cout << n << "  -> (" << n1 << " , " << n2 << ")" << endl;
	  Word c1 = TheGrigorchukGroupAlgorithms::cosetRepresentativeKSbgp( n1 );
	  Word c2 = TheGrigorchukGroupAlgorithms::cosetRepresentativeKSbgp( n2 );
	  cout << "      -> (" << c1 << " , " << c2 << ")" << endl;
	  
	  // now starting from (n1,n2) compute the image of the coset C r -> (C r1, C r2)
	  map< PII , Word > I;
	  map< PII , Word > I2;
	  preimages[PII(n1,n2)] = n;
	  I2[PII(n1,n2)] = representative;
	  while( !I2.empty( ) ) {
	    
	    pair< PII , Word > P = *I2.begin( );
	    I2.erase( I2.begin( ) );
	    I[P.first] = P.second;
	    // cout << "   +++  Current pair ( " << P.first.first << "," << P.first.second << ")" << endl;
	    
	    for( int i=0 ; i<3 ; ++i ) {
	      Word nw = P.second * G[i];
	      // cout << nw << endl;
	      pair< Word , Word > S = TheGrigorchukGroupAlgorithms::split( nw );
	      int nw1 = TheGrigorchukGroupAlgorithms::cosetRepresentativeKSbgp( S.first );
	      int nw2 = TheGrigorchukGroupAlgorithms::cosetRepresentativeKSbgp( S.second );

	      if( I.find(PII(nw1,nw2))==I.end( ) && I2.find(PII(nw1,nw2))==I2.end( ) ) {
		cout << "   -> (" << nw1 << " , " << nw2 << ")" << endl;
		Word c1 = TheGrigorchukGroupAlgorithms::cosetRepresentativeKSbgp( nw1 );
		Word c2 = TheGrigorchukGroupAlgorithms::cosetRepresentativeKSbgp( nw2 );
		Word lift = TheGrigorchukGroupAlgorithms::liftToSTone( c1 , c2 ).first;
		cout << "      -> (" << c1 << " , " << c2 << ")" << endl;
		cout << "          " << lift << endl;
		
		I2[PII(nw1,nw2)] = nw;
		preimages[PII(nw1,nw2)] = n;
	      }
	    }
	  }
	}
      }
    }
  }
  
  
  


}


//---------------------------------------------------------------------------//
//-------------------------------------- main -------------------------------//
//---------------------------------------------------------------------------//


void test_st_lift( )
{
  int exp = 100;
  for( int e=0 ; e<exp ; ++e ) {
    
    Word w1 = TheGrigorchukGroupAlgorithms::randomWord( 100 );
    Word w2 = TheGrigorchukGroupAlgorithms::randomWord( 100 );
    
    pair< Word , Word > L = TheGrigorchukGroupAlgorithms::liftToSTone( w1 , w2 );
    pair< Word , Word > S = TheGrigorchukGroupAlgorithms::split( L.first );
    if( TheGrigorchukGroupAlgorithms::trivial(-w2*S.second) && TheGrigorchukGroupAlgorithms::trivial(L.second*-w1*S.first) ) {
      cout << "ST lift - ok" << endl;
      cout << "   " << L.first.length() << endl;
    } else {
      cout << "ST lift - error" << endl;
      exit( 1 );
    }
    
  }
}


//---------------------------------------------------------------------------//
//-------------------------------------- main -------------------------------//
//---------------------------------------------------------------------------//


void test_Bsbgp( )
{
  int exp = 100;
  for( int e=0 ; e<exp ; ++e ) {
    
    Word w = TheGrigorchukGroupAlgorithms::randomWord( 100 );
    Word r = TheGrigorchukGroupAlgorithms::reduce( w );
    pair< Word , list< Word > > D = TheGrigorchukGroupAlgorithms::decompositionBSbgp( w );
    
    Word to_check;
    for( list< Word >::const_iterator l_it=D.second.begin( ) ; l_it!=D.second.end() ; ++l_it ) {
      Word factor = -*l_it * Word(2) * *l_it;
      to_check *= factor;
    }
    to_check *= D.first;
    to_check *= -w;
    
    if( TheGrigorchukGroupAlgorithms::trivial( to_check ) ) {
      cout << "Decomposition - ok" << endl;
    } else {
      cout << "Decomposition - error" << endl;
      exit(1);
    }
  }
}


//---------------------------------------------------------------------------//
//-------------------------------------- main -------------------------------//
//---------------------------------------------------------------------------//


void test_CP( )
{
  int exp = 10000;
  for( int e=0 ; e<exp ; ++e ) {
    cout << "----------------------------------------------------" << endl;
    Word w = TheGrigorchukGroupAlgorithms::randomWord( 200 );
    Word c = TheGrigorchukGroupAlgorithms::randomWord( 400 );
    
    /*
    set< Word > solutions = TheGrigorchukGroupAlgorithms::findConjugator_Kcosets( w , -c * w * c );
    if( !solutions.size( ) ) {
      cout << "Conjugacy: error" << endl;
      exit(1);
    } else {
      cout << e << " - Conjugacy: ok" << endl;
      for( set< Word >::const_iterator s_it=solutions.begin( ) ; s_it!=solutions.end( ) ; ++s_it ) {
	Word conj = *s_it;
	// cout << TheGrigorchukGroupAlgorithms::reduce( -conj*w*conj * -c*-w*c ) << endl;
	// cout << conj << endl;
	// cout << c << endl;
	if( TheGrigorchukGroupAlgorithms::trivial( -conj*w*conj * -c*-w*c ) ) {
	  // cout << "      " << c << " - ok" << endl;
	} else {
	  cout << "      " << c << " - error" << endl;
	  exit(1);
	}
      }
    }
    */
    
    if( !TheGrigorchukGroupAlgorithms::conjugate( w , -c * w * c ) ) {
      cout << "Conjugacy: error" << endl;
      exit(1);
    } else {
      cout << e << " - Conjugacy: ok" << endl;
    }

  }
}


//---------------------------------------------------------------------------//
//-------------------------------------- main -------------------------------//
//---------------------------------------------------------------------------//


void test_WP( )
{
  int exp = 100;
  for( int e=0 ; e<exp ; ++e ) {
    Word w = Word::randomWord( 4 , 1000 );
    Word r = TheGrigorchukGroupAlgorithms::reduce( w );
    if( TheGrigorchukGroupAlgorithms::trivial( -r*w ) ) {
      cout << "Trivial - ok" << endl;
    } else {
      cout << "Trivial - error" << endl;
      exit(1);
    }
  }
}

//---------------------------------------------------------------------------//
//-------------------------------------- main -------------------------------//
//---------------------------------------------------------------------------//


void test_reduction( )
{
  int exp = 1;
  for( int e=0 ; e<exp ; ++e ) {
    Word w = Word::randomWord( 4 , 10 );
    cout << "w = " << w << endl;

    Word r = TheGrigorchukGroupAlgorithms::reduce( w );
    cout << "r = " << r << endl;
  }
}


//---------------------------------------------------------------------------//
//-------------------------------------- main -------------------------------//
//---------------------------------------------------------------------------//


void test_order( )
{
  int L = 1000;
  int exp = 100;
  for( int e=0 ; e<exp ; ++e ) {
    Word w = TheGrigorchukGroupAlgorithms::randomWord( L );
    // cout << "w = " << w << endl;
    
    int O = TheGrigorchukGroupAlgorithms::findOrder( w );
    cout << L << "  " << O << endl;
  }
}


//---------------------------------------------------------------------------//
//-------------------------------------- main -------------------------------//
//---------------------------------------------------------------------------//


int main( )
{
  // test_reduction( );
  // test_WP( );
  // test_Bsbgp( );
  // test_st_lift( );
  test_lifts( );

  // test_CP();
  // test_order( );

  // constructLevelGraph( );

  return 0;
}
