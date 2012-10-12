// Copyright (C) 2005 Alexander Ushakov
// Contents: Implementation of class WordRep
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//


#include <math.h>
#include "WordRep.h"
#include "WordIterator.h"
#include <stdlib.h>

using namespace std;

//---------------------------------------------------------------------------//
//-------------------------------- WordRep ----------------------------------//
//---------------------------------------------------------------------------//


ostream& WordRep::printOn( ostream& os ) const
{
  int p_base = 0;
  int deg = 0;

  if( theElements.size()==0 ) {
    os << "1";
    return os;
  }

  int pb;
	for( list< int >::const_iterator e_it = theElements.begin( ) ; e_it!=theElements.end() ; ++e_it ) {
    if( *e_it==p_base )
      ++deg;
    else {
      pb = p_base;
      if( pb!=0 ) {
        os << 'x' << abs( p_base );
        if( deg!=1 || pb<0 )
          os << "^" << ( pb<0 ? -deg : deg );
        os << " ";
      }
      pb = p_base = *e_it;
      deg = 1;
    }
  }
	
	os << 'x' << abs( p_base );
  if( deg!=1 || pb<0 )
    os << "^" << ( pb<0 ? -deg : deg );

	return os;
}


//=========================================================


void WordRep::push_back( int g )
{
  if( !theElements.size( ) )
    theElements.push_back( g );
  else {
    int e = *(--theElements.end());
    if( e+g==0 )
      theElements.pop_back( );
    else
      theElements.push_back( g );
  }
}


//=========================================================


void WordRep::push_front( int g )
{
  if( !theElements.size( ) )
    theElements.push_front( g );
  else {
    int e = *(theElements.begin());
    if( e+g==0 )
      theElements.pop_front( );
    else
      theElements.push_front( g );
  }
}


//=========================================================


WordRep::WordRep( )
{ 
  
}

//=========================================================

WordRep::WordRep( const list< int >& gens )
{
  for( list< int >::const_iterator g_it = gens.begin( ) ; g_it!=gens.end( ) ; ++g_it )
    push_back( *g_it );
}

//=========================================================

WordRep::WordRep( const vector< int >& gens )
{
  for( vector< int >::const_iterator g_it = gens.begin( ) ; g_it!=gens.end( ) ; ++g_it )
    push_back( *g_it );
}

//=========================================================


WordRep::WordRep( const WordRep& wr ) :
  theElements( wr.theElements )
{
  
}


//=========================================================


WordRep::WordRep( int g )
{ 
  push_back( g );
}


//=========================================================


WordRep WordRep::operator=( const WordRep wr ) 
{
  theElements = wr.theElements;
  return *this;
}

//=========================================================

WordRep& WordRep::operator ^= ( int power )
{
  WordRep base = ( power>0 ? *this : this->inverse( ) );
  this->clear( );
  int ap = abs( power );
  for( int i=0 ; i<ap ; ++i )
    *this *= base;
  return *this;
}

//=========================================================

WordRep& WordRep::operator ^= ( const WordRep& w )
{
  *this *= w;
  for( list< int >::const_iterator m_it = w.theElements.begin( ) ; m_it!=w.theElements.end( ) ; ++m_it )
    push_front( -*m_it );
  return *this;
}

//=========================================================


WordRep& WordRep::operator *= ( const WordRep& w ) 
{
  for( list< int >::const_iterator m_it = w.theElements.begin( ) ; m_it!=w.theElements.end( ) ; ++m_it )
    push_back( *m_it );
  return *this;
}

//=========================================================

bool WordRep::doesContain( int gen ) const
{
  for( list< int >::const_iterator g_it = theElements.begin( ) ; g_it!=theElements.end( ) ; ++g_it )
    if( *g_it==gen || *g_it==-gen )
      return true;

  return false;
}


//=========================================================


void WordRep::cyclicLeftShift( )
{
  if( theElements.size( )<2 ) 
    return;
  push_back( *theElements.begin( ) );
  theElements.pop_front( );
}


//=========================================================


void WordRep::cyclicRightShift( )
{
  if( theElements.size( )<2 ) 
    return;
  push_front( *--theElements.end( ) );
  theElements.pop_back( );
}


//=========================================================


int WordRep::isIn( int gen ) const
{
  int result = 0;
  list< int >::const_iterator g_it = theElements.begin( );
  for( ; g_it!=theElements.end( ) ; ++g_it )
    if( *g_it==gen || *g_it==-gen )
      result += 1;
  return result;
}


//=========================================================


int WordRep::exponentSum( int gen ) const
{
  int result = 0;
  list< int >::const_iterator g_it = theElements.begin( );
  for( ; g_it!=theElements.end( ) ; ++g_it )
    if( *g_it==gen || *g_it==-gen )
      result += 1;
  return result;
}


//=========================================================


WordRep WordRep::inverse( ) const
{
  WordRep result;
  list< int >::const_iterator g_it = theElements.begin( );
  for( ; g_it!=theElements.end( ) ; ++g_it )
    result.push_front( -*g_it );
  return result;
}


//=========================================================


bool WordRep::operator == ( const WordRep& wr ) const
{
  if( theElements.size( )!=wr.theElements.size( ) )
    return false;

  list< int >::const_iterator g_it1 = theElements.begin( );
  list< int >::const_iterator g_it2 = wr.theElements.begin( );
  while( g_it1!=theElements.end( ) )
    if( *(g_it1++) != *(g_it2++) )
      return false;

  return true;
}


//=========================================================


bool WordRep::operator < ( const WordRep& wr ) const
{
  int len1 = theElements.size( );
  int len2 = wr.theElements.size( );
  if( len1<len2 ) return true;
  if( len1>len2 ) return false;

  list< int >::const_iterator g_it1 = theElements.begin( );
  list< int >::const_iterator g_it2 = wr.theElements.begin( );
  for( int t=0 ; t<len1 ; ++t ) {
    if( *g_it1 < *g_it2 )
      return true;
    if( *(g_it1++) > *(g_it2++) )
      return false;
  }
  
  return false;
}


//=========================================================


bool WordRep::operator > ( const WordRep& wr ) const
{
  int len1 = theElements.size( );
  int len2 = wr.theElements.size( );
  if( len1<len2 ) return true;
  if( len1>len2 ) return false;
  
  list< int >::const_iterator g_it1 = theElements.begin( );
  list< int >::const_iterator g_it2 = wr.theElements.begin( );
  for( int t=0 ; t<len1 ; ++t ) {
    if( *g_it1 > *g_it2 )
      return true;
    if( *(g_it1++) < *(g_it2++) )
      return false;
  }

  return false;
}


//=========================================================


void WordRep::insert( WordIterator it , int g )
{
  theElements.insert( it.theIterator , g );
}

template< class ConstIntIterator > 
void WordRep::insert( WordIterator it , ConstIntIterator B , ConstIntIterator E )
{
  theElements.insert( it.theIterator , B , E );
}

void WordRep::insert( int pos , int g )
{
  if( pos<0 || theElements.size( )<pos ) return;
  list< int >::iterator l_it = theElements.begin( );
  for( ; pos>0 ; --pos, ++l_it );
  theElements.insert( l_it , g );
}

template< class ConstIntIterator > 
void WordRep::insert( int pos , ConstIntIterator B , ConstIntIterator E )
{
  if( pos<0 || theElements.size( )<pos ) return;
  list< int >::iterator l_it = theElements.begin( );
  for( ; pos>0 ; --pos, ++l_it );
  theElements.insert( l_it , B , E );
}


//=========================================================

void WordRep::replace( WordIterator it , const Generator& g )
{
  *(it.theIterator) = g;
}

void WordRep::replace( int pos , const Generator& g )
{
  if( pos<0 || theElements.size( )<pos ) return;
  list< int >::iterator l_it = theElements.begin( );
  for( ; pos>0 ; --pos, ++l_it );
  *l_it = g;
}

//=========================================================

template< class ConstIntIterator > 
void WordRep::replace( WordIterator it , ConstIntIterator B , ConstIntIterator E )
{
  // list< int >::iterator l_it = 
  list< int >::iterator c_it = it.theIterator;
  for( ; c_it!=theElements.end( ) && B!=E ; ++c_it, ++B )
    *c_it = *B;
}

//=========================================================


void WordRep::cyclicallyPermute( int n )
{
  int len = theElements.size( );
  if( n==0 || len==0 ) return;

  // make it positive
  n %= len;
  n = n<0 ? n+len : n;
  
  WordRep segm1 = *this;
  segm1.initialSegment( n );
  WordRep segm2 = *this;
  segm2.terminalSegment( n );
  *this = segm2*segm1;
}


//=========================================================


void WordRep::segment( int from , int to )
{
  initialSegment( to );
  terminalSegment( from );
}

//=========================================================

void WordRep::initialSegment( int len )
{
  if( len<=0 ) {
    theElements.clear( );
    return;
  }
  
  for( int left=theElements.size( ) - len ; left>0 ; --left )
    theElements.pop_back( );
}


//=========================================================


void WordRep::terminalSegment( int len )
{
  int l = len < theElements.size( ) ? len : theElements.size( );
  for( int t=0 ; t<l ; ++t )
    theElements.pop_front( );
}

//=========================================================


void WordRep::cyclicallyReduce( )
{
  while( theElements.size( )>1 ) {
    int b = *theElements.begin( );
    int e = *--theElements.end( );
    if( e+b )
      break;
    theElements.pop_back( );
    theElements.pop_front( );
  }
}


//=========================================================


void WordRep::cyclicallyReduce( WordRep& conjugator )
{
  conjugator.theElements.clear( );
  while( theElements.size( )>1 ) {
    int b = *theElements.begin( );
    int e = *--theElements.end( );
    if( e+b )
      break;
    theElements.pop_back( );
    theElements.pop_front( );
    conjugator.theElements.push_front( e );
  }
}


//=========================================================


int WordRep::getPower( WordRep& base ) const
{
  int len = theElements.size( );
  if( len<=1 ) return 1;
  
  int up = sqrt((double)len)+1;
  
  // Decompose len into a product of primes 
  list< pair< int , int > > primes;
  list< pair< int , int > > primes_to_use;
  
  for( int i=2 ; i<up ; ++i ) {
    bool prime = true;
    list< pair< int , int > >::iterator p_it = primes.begin( );
    for( ; p_it!=primes.end( ) ; ++p_it ) {
      int a = (*p_it).first;
      int b = (*p_it).second;
      (*p_it).second = b+1==a ? (prime=false) : b+1;
    }
    if( prime ) {
      primes.push_back( pair< int , int >(i,0) );
      int count = 0;
      while( len%i==0 ) {
	++count;
	len /= i;
      }
      if( count ) {
	primes_to_use.push_back( pair<int,int>(i,count) );
	up = sqrt((double)len)+1;
      }
    }
  }
  if( len>1 )
    primes_to_use.push_back( pair<int,int>(len,1) );

  /*{
	list< pair< int , int > >::iterator p_it = primes_to_use.begin( );
	for( ; p_it!=primes_to_use.end( ) ; ++p_it ) {
		int a = (*p_it).first;
		int b = (*p_it).second;
		cout << " " << a << "," << b << endl;
	}
	}*/

  len = theElements.size( );
  
  // Construct a vector of generators (to have indexed access)
  vector< int > wrd_vct(len);
  list< int >::const_iterator g_it = theElements.begin( );
  for( int j=0 ; g_it!=theElements.end() ; ++j, ++g_it )
    wrd_vct[j] = *g_it;
  
  // Find the power
  list< pair< int , int > >::iterator p_it = primes_to_use.begin( );
  for( ; p_it!=primes_to_use.end( ) ; ++p_it ) {
    int parts = (*p_it).first;
    int count = (*p_it).second;
    int shrt_len = len/parts;
    
    bool progress = true;
    for( int t=0 ; t<count && progress ; ++t ) {
      for( int offset=0 ; offset<shrt_len  && progress ; ++offset ) {
	int e = wrd_vct[offset];
	for( int j=offset ; j<len && progress ; j+=shrt_len )
	  if( wrd_vct[j]!=e )
	    progress = false;
      }
      if( progress )
	shrt_len = (len = shrt_len) / parts;
    }
  }
  
  base = *this;
  base.initialSegment( len );
  return theElements.size()/len;
}


//=========================================================


void WordRep::freelyReduce( WordIterator beg , WordIterator end )
{
  list< int >::iterator b = beg.theIterator;
  list< int >::iterator e = end.theIterator;
  
  for( list< int >::iterator cur_it=b ; cur_it!=e && cur_it!=theElements.end( ) ; ) {
    
    list< int >::iterator cur_it2 = cur_it;
    if( ++cur_it2==theElements.end( ) )
      break;
    if( *cur_it+*cur_it2 ) {
      ++cur_it;
      continue;
    }
    
    if( cur_it2==e )
      ++e;
    theElements.erase( cur_it );
    cur_it = theElements.erase( cur_it2 );
    if( cur_it!=theElements.begin( ) )
      --cur_it;
  }
  
  
}
