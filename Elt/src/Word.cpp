// Copyright (C) 2005 Alexander Ushakov
// Contents: Implementation of class Word
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//

#include <set>
#include <map>
using namespace std;


#include "RanlibCPP.h"
#include "stdlib.h"
#include "Word.h"
// #include <strstream>


//---------------------------------------------------------------------------//
//---------------------------------- Word -----------------------------------//
//---------------------------------------------------------------------------//


Word& Word::push_back ( const Word& w )
{
  Word::const_iterator w_it = w.begin( );
  for( ; w_it!=w.end( ) ; ++w_it )
    push_back( *w_it );
  
  return *this;
}


//=========================================================


Word& Word::push_front( const Word& w )
{
  Word::const_iterator w_it = w.end( );
  for( ; w_it!=w.begin( ) ; )
    push_front( *--w_it );
  
  return *this;
}


//=========================================================


Word Word::power( int t ) const
{
  if( t==0 ) return Word( );
  
  Word result = ( t<0 ? this->inverse() : *this );
  for( int i=1 ; i<abs( t ) ; ++i )
    result *= ( t<0 ? this->inverse() : *this );
  return result;
}


//=========================================================


Word Word::randomWord( int gens , int wLen )
{
  if( wLen==0 )	return Word( );
  
  int old = 0;
  list< int > result;
  for( int i=0 ; i<wLen ; ++i ) {
    int div = i==0 ? 2*gens : 2*gens-1;
    // int g = ::rand()%div-gens;
    int g = RandLib::ur.irand( 0 , div-1 )-gens;
    g = g>=0 ? g+1 : g;
    if( g+old==0 )
      g = gens;
    result.push_back( old = g );
  }

  return result;
}

//=========================================================

Word Word::randomWord( int gens , int wLenMin, int wLenMax )
{

  int wLen = RandLib::ur.irand( wLenMin , wLenMax );
  return randomWord(  gens ,  wLen );
}

//=========================================================

Word Word::freelyReduce( iterator beg , iterator end )
{
  Word result( *this );
  result.change( ) -> freelyReduce( beg , end );
  return result;
}

//=========================================================

template< class ConstIntIterator > 
void Word::insert( int pos , ConstIntIterator B , ConstIntIterator E )
{
  change( )->insert( pos , B , E );
}

void Word::insert( int pos , int g )
{
  change( )->insert( pos , g );
}

template< class ConstIntIterator > 
void Word::insert( WordIterator it , ConstIntIterator B , ConstIntIterator E )
{
  change( )->insert( it , B , E );
}


void Word::insert( WordIterator it , int g )
{
  change( )->insert( it , g );
}

//=========================================================

void Word::replace( WordIterator it , const Generator& g )
{
  change( )->replace( it , g );
}


void Word::replace( int pos , const Generator& g )
{
  change( )->replace( pos , g );
}

//=========================================================

template< class ConstIntIterator > 
void Word::replace( WordIterator it , ConstIntIterator B , ConstIntIterator E )
{
  change( )->replace( it , B , E );
}


//=========================================================

Word Word::replaceGenerators( const vector<Word>& images ) const
{
	Word newWord;
  
	for ( ConstWordIterator I = this->begin(); I!=this->end(); I++){
    	Generator g = *I;
//    	cout << abs(g)-1  << " " << images.size() << endl;
		if (static_cast<unsigned int>(abs(g)) >= images.size() + 1){
	  cout <<"Word::replaceGenerators() : image vector index overflow." <<endl; 
	  exit(1);
	}
    	newWord *= ( g > 0 ? images[abs(g)-1] :  images[abs(g)-1].inverse() );
  	}

	return newWord; 
}


//=========================================================


Word Word::minimalEquivalentForm( const set< int >& permutableGenerators , bool inverses , bool cyclicPermutations ) const
{
  Word w;
  int P = 1;
  if( cyclicPermutations ) {
    w = this->cyclicallyReduce( );
    P = w.length( );
  } else {
    w = *this;
  }
  Word result = w;
  int I = inverses ? 2:1;

  // cout << "   >>>   " << cur << endl;
  
  for( int i=0;  i<I ; ++i ) {
    Word cur = i==0 ? w : -w;
    for( int p=0;  p<P ; ++p , cur.cyclicLeftShift( ) ) {

      // prepare a permutation
      map< int , int > permutation;
      for( set< int >::const_iterator g_it = permutableGenerators.begin( ) ; g_it!=permutableGenerators.end( ) ; ++g_it )
	permutation[*g_it] = permutation[-*g_it] = 0;
      
      // slide along the word and replace permutable generators to guarantee minimality
      Word rep;
      set< int >::const_iterator counter_it = permutableGenerators.end( );
      for( Word::iterator w_it=cur.begin( ) ; w_it!=cur.end() ; ++w_it ) {
	int g = *w_it;
	// check if g is permutable
	map< int , int >::iterator p_it = permutation.find(g);
	if( p_it==permutation.end( ) ) { // non-permutable
	  rep.push_back(g);
	} else { // permutable
	  if( (*p_it).second==0 ) {
	    (*p_it).second  = -*(--counter_it);
	    permutation[-g] =  *counter_it;
	  }
	  rep.push_back( (*p_it).second );
	}
      }
      // cout << rep << endl;
      if( rep<result )
	result = rep;

      //for( set< int >::const_iterator g_it = permutableGenerators.begin( ) ; g_it!=permutableGenerators.end( ) ; ++g_it )
      // cout << "  " << permutation[*g_it] << " , " << permutation[-*g_it] << endl;

    }
  }
  
  // cout << "   <<<   " << result << endl;

  return result;
}
