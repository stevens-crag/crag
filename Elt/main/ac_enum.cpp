
#include <set>
#include <list>
#include <fstream>
#include "Word.h"


// typedef unsigned int UI;
typedef unsigned long long int UI;
int MAX_LENGTH = 35;

typedef pair< Word , Word > PWW;
typedef pair< UI , UI > PII;


//---------------------------------------------------------------------------//
//------------------------- smallestEquivalentWord --------------------------//
//---------------------------------------------------------------------------//


void smallestEquivalentWord( Word& w )
{
  w.cyclicallyReduceWord( );
  Word result = w;
  for( int d=0 ; d<2 ; ++d ) {
    if( d ) w = -w;
    for( int s=0 ; s<w.length( ) ; ++s ) {
      w.cyclicLeftShift( );
      if( result<w ) result = w;
    }
  }
  w = result;
}


//---------------------------------------------------------------------------//
//------------------------------ number_2_word ------------------------------//
//---------------------------------------------------------------------------//

Word number_2_word( UI n )
{
  Word result;
  
  list< int > decomp;
  while( n>7 ) {
    decomp.push_front( n%3 );
    n /= 3;
  }

  n = n%4;
  int prev_gen = n<2 ? n-2 : n-1;
  result.push_back( prev_gen );

  while( !decomp.empty( ) ) {
    n = *decomp.begin( );
    decomp.pop_front( );
    int gen = n<2 ? n-2 : n-1;
    if( gen+prev_gen==0 )
      gen = 2;
    prev_gen = gen;
    result.push_back( prev_gen );
  }

  return result;
}


//---------------------------------------------------------------------------//
//------------------------------ word_2_number ------------------------------//
//---------------------------------------------------------------------------//


UI word_2_number( const Word& w )
{
  // w is assumed to be reduced and non-trivial
  UI result = 4;
  auto w_it = w.begin( );
  int prev_gen = *(w_it++);
  result += prev_gen<0 ? prev_gen+2 : prev_gen+1;

  for( ; w_it!=w.end( ) ; ++w_it ) {

    int gen = *w_it;
    int num = ( gen<0 ? gen+2 : gen+1 );
    if( num==3 ) 
      num = ( prev_gen>0 ? -prev_gen+2 : -prev_gen+1 );
    result *= 3;
    result += num;
    prev_gen = gen;
  }

  return result;
}


//---------------------------------------------------------------------------//
//------------------------------ add_new_pair -------------------------------//
//---------------------------------------------------------------------------//


void add_new_pair( PWW p , set< PII >& checked_pairs , set< PII >& unchecked_pairs )
{
  smallestEquivalentWord( p.first  );
  smallestEquivalentWord( p.second );
  if( p.first.length( )>MAX_LENGTH || p.second.length( )>MAX_LENGTH )
    return;

  UI n1 = word_2_number( p.first  );
  UI n2 = word_2_number( p.second );
  if( n1<n2 ) swap( n1 , n2 );
  PII pii( n1 , n2 );
  if( checked_pairs.find( pii )!=checked_pairs.end( ) )
    return;

  unchecked_pairs.insert( pii );
}


//---------------------------------------------------------------------------//
//------------------------------ process_pair -------------------------------//
//---------------------------------------------------------------------------//


void process_pair( PII& p , set< PII >& checked_pairs , set< PII >& unchecked_pairs )
{
  Word w1 = number_2_word( p.first );
  Word w2 = number_2_word( p.second );
  static ofstream OF( "checked.txt" );
  OF << w1 << " , " << w2 << endl;

  for( int t1=0 ; t1<w1.length() ; ++t1 ) {
    w1.cyclicLeftShift( );
    for( int t2=0 ; t2<w2.length() ; ++t2 ) {
      w2.cyclicLeftShift( );
      add_new_pair( PWW( w1     , w2* w1 ) , checked_pairs , unchecked_pairs );
      add_new_pair( PWW( w1     , w2*-w1 ) , checked_pairs , unchecked_pairs );
      add_new_pair( PWW( w1* w2 , w2     ) , checked_pairs , unchecked_pairs );
      add_new_pair( PWW( w1*-w2 , w2     ) , checked_pairs , unchecked_pairs );
    }
  }
}


//---------------------------------------------------------------------------//
//----------------------------- enumerate_pairs -----------------------------//
//---------------------------------------------------------------------------//


void enumerate_pairs( PWW p )
{
  set< PII > checked_pairs;
  set< PII > unchecked_pairs;

  smallestEquivalentWord( p.first  );
  smallestEquivalentWord( p.second );
  int n1 = word_2_number( p.first  );
  int n2 = word_2_number( p.second );
  if( n1<n2 ) swap( n1 , n2 );
  unchecked_pairs.insert( PII(n1,n2) );
  
  for( int i=0 ; !unchecked_pairs.empty( ) && i<100000 ; ++i ) {
    PII p = *unchecked_pairs.begin( );
    unchecked_pairs.erase( unchecked_pairs.begin( ) );
    checked_pairs.insert( p );
    process_pair( p , checked_pairs , unchecked_pairs );
  }
}


//---------------------------------------------------------------------------//
//----------------------------------- main ----------------------------------//
//---------------------------------------------------------------------------//


int main( )
{
  Word x( 1 );
  Word y( 2 );
  enumerate_pairs( PWW(x,y) );

  return 0;
}
