// Copyright (C) 2005 Alexander Ushakov
// Contents: Implementation of class PowerWordRep
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//


#include "math.h"
#include "PowerWordRep.h"


//---------------------------------------------------------------------------//
//-------------------------------- PowerWordRep ----------------------------------//
//---------------------------------------------------------------------------//


void PowerWordRep::pushGeneratorBack( int g , int p )
{
  if( !theLength ) {
    theElements.push_back( PII(g,p) );
    theLength += abs(p);
  } else {
    PII& e = *(--theElements.end());
    int g1 = e.first;
    int p1 = e.second;
    if( g1!=g ) {
      theElements.push_back( PII(g,p) );
      theLength += abs(p);
    } else {
      theLength -= abs( p1 );
      if( p1+p ) {
	e.second = p1+p;
	theLength += abs( p1+p );
      } else {
	theElements.pop_back( );
      }
    }
  }
}


//=========================================================


void PowerWordRep::pushGeneratorFront( int g , int p )
{
  if( !theLength ) {
    theElements.push_front( PII(g,p) );
    theLength += abs(p);
  } else {
    PII& e = *theElements.begin();
    int g1 = e.first;
    int p1 = e.second;
    if( g1!=g ) {
      theElements.push_front( PII(g,p) );
      theLength += abs(p);
    } else {
      theLength -= abs( p1 );
      if( p1+p ) {
	e.second = p1+p;
	theLength += abs( p1+p );
      } else
	theElements.pop_front( );
    }
  }
}


//=========================================================


PowerWordRep::PowerWordRep( ) :
  theLength( 0 )
{ 

}

//=========================================================

PowerWordRep::PowerWordRep( const list< int >& gens ) :
  theLength(0)
{
  for( list< int >::const_iterator g_it = gens.begin( ) ; g_it!=gens.end( ) ; ++g_it )
    if( *g_it>0 )
      pushGeneratorBack( *g_it , 1 );
    else
      pushGeneratorBack( -*g_it , -1 );
}

PowerWordRep::PowerWordRep( const list< pair< int , int > >& gens ) :
  theLength(0)
{
  for( list< pair< int , int > >::const_iterator g_it = gens.begin( ) ; g_it!=gens.end( ) ; ++g_it )
    pushGeneratorBack( (*g_it).first , (*g_it).second );
}

//=========================================================

PowerWordRep::PowerWordRep( const vector< int >& gens ) :
  theLength(0)
{
  for( vector< int >::const_iterator g_it = gens.begin( ) ; g_it!=gens.end( ) ; ++g_it )
    if( *g_it>0 )
      pushGeneratorBack( *g_it , 1 );
    else
      pushGeneratorBack( -*g_it , -1 );
  
}

PowerWordRep::PowerWordRep( const vector< pair< int , int > >& gens ) :
  theLength(0)
{
  for( vector< pair< int , int > >::const_iterator g_it = gens.begin( ) ; g_it!=gens.end( ) ; ++g_it )
    pushGeneratorBack( (*g_it).first , (*g_it).second );
}

//=========================================================


PowerWordRep::PowerWordRep( const PowerWordRep& wr ) :
  theElements( wr.theElements ),
  theLength( wr.theLength )
{

}


//=========================================================


PowerWordRep::PowerWordRep( int g , int p )
{ 
  pushGeneratorBack( g , p );
  theLength = abs(p);
}


//=========================================================


PowerWordRep PowerWordRep::operator=( const PowerWordRep wr ) 
{
  theElements = wr.theElements;
  theLength = wr.theLength;
  return *this;
}


//=========================================================


PowerWordRep& PowerWordRep::operator *= ( const PowerWordRep& w ) 
{
  for( list< PII >::const_iterator m_it = w.theElements.begin( ) ; m_it!=w.theElements.end( ) ; ++m_it )
    pushGeneratorBack( (*m_it).first , (*m_it).second );
  
  return *this;
}

//=========================================================

bool PowerWordRep::doesContain( int gen ) const
{
  for( list< PII >::const_iterator g_it = theElements.begin( ) ; g_it!=theElements.end( ) ; ++g_it )
    if( (*g_it).first==gen )
      return true;

  return false;
}


//=========================================================


void PowerWordRep::cyclicLeftShift( )
{
  if( theLength<=1 ) return;
  int& g = (*theElements.begin( )).first;
  int& p = (*theElements.begin( )).second;
  int s = p<0 ? -1 : 1;
  
  if( abs(p)==1 )
    theElements.pop_front( );
  else
    p -= s;
  
  theLength -= 1;
  pushGeneratorBack( g , s );
}


//=========================================================


void PowerWordRep::cyclicRightShift( )
{
  if( theLength<=1 ) return;
  int& g = (*--theElements.end( )).first;
  int& p = (*--theElements.end( )).second;
  int s = p<0 ? -1 : 1;
  
  if( abs(p)==1 )
    theElements.pop_back( );
  else
    p -= s;
  
  theLength -= 1;
  pushGeneratorFront( g , s );
}


//=========================================================


int PowerWordRep::isIn( int gen ) const
{
  int result = 0;
  list< PII >::const_iterator g_it = theElements.begin( );
  for( ; g_it!=theElements.end( ) ; ++g_it )
    if( (*g_it).first==gen )
      result += abs((*g_it).second);
  return result;
}


//=========================================================


int PowerWordRep::exponentSum( int gen ) const
{
  int result = 0;
  list< PII >::const_iterator g_it = theElements.begin( );
  for( ; g_it!=theElements.end( ) ; ++g_it )
    if( (*g_it).first==gen )
      result += (*g_it).second;
  return result;
}


//=========================================================


PowerWordRep PowerWordRep::inverse( ) const
{
	PowerWordRep result;
  list< PII >::const_iterator g_it = theElements.begin( );
  for( ; g_it!=theElements.end( ) ; ++g_it )
    result.pushGeneratorFront( (*g_it).first , -(*g_it).second );

  return result;
}


//=========================================================


bool PowerWordRep::operator == ( const PowerWordRep& wr ) const
{
  if( theElements.size( )!=wr.theElements.size( ) || theLength!=wr.theLength )
    return false;

  list< PII >::const_iterator g_it1 = theElements.begin( );
  list< PII >::const_iterator g_it2 = wr.theElements.begin( );
  for( int t=0 ; t<theElements.size( ) ; ++t )
    if( *(g_it1++) != *(g_it2++) )
      return false;

  return true;
}


//=========================================================


bool PowerWordRep::operator < ( const PowerWordRep& wr ) const
{
  if( theLength<wr.theLength )
    return true;
  if( theLength>wr.theLength )
    return false;
  int len1 = theElements.size( );
  int len2 = wr.theElements.size( );
  if( len1<len2 ) return true;
  if( len1>len2 ) return false;

  list< PII >::const_iterator g_it1 = theElements.begin( );
  list< PII >::const_iterator g_it2 = wr.theElements.begin( );
  for( int t=0 ; t<len1 ; ++t ) {
    if( *g_it1 < *g_it2 )
      return true;
    if( *(g_it1++) > *(g_it2++) )
      return false;
  }
  
  return false;
}


//=========================================================


bool PowerWordRep::operator > ( const PowerWordRep& wr ) const
{
  if( theLength>wr.theLength )
    return true;
  if( theLength<wr.theLength )
    return false;
  int len1 = theElements.size( );
  int len2 = wr.theElements.size( );
  if( len1<len2 ) return true;
  if( len1>len2 ) return false;


  list< PII >::const_iterator g_it1 = theElements.begin( );
  list< PII >::const_iterator g_it2 = wr.theElements.begin( );
  for( int t=0 ; t<len1 ; ++t ) {
    if( *g_it1 > *g_it2 )
      return true;
    if( *(g_it1++) < *(g_it2++) )
      return false;
  }

  return false;
}


//=========================================================

void PowerWordRep::insert( const PowerWordRep& wr , int pos )
{
  /*
    if( pos<0 || theElements.size( )<pos ) return;
    list< int >::iterator l_it = theElements.begin( );
    for( int i=0 ; i<pos+1 ; ++i )
    ++l_it;
    theElements.insert( l_it , wr.theElements.begin( ) , wr.theElements.end( ) );
  */
}


//=========================================================



void PowerWordRep::cyclicallyPermute( int n )
{
  int len = theElements.size( );
  if( n==0 || len==0 ) return;

  // make it positive
  n %= len;
  n = n<0 ? n+len : n;
  
  PowerWordRep segm1 = *this;
  segm1.initialSegment( n );
  PowerWordRep segm2 = *this;
  segm2.terminalSegment( n );
  *this = segm2*segm1;
}


//=========================================================


void PowerWordRep::segment( int from , int to )
{
  initialSegment( to );
  terminalSegment( from );
}

//=========================================================

void PowerWordRep::initialSegment( int len )
{
  if( len<=0 ) {
    theElements.clear( );
    theLength = 0;
    return;
  }

  // remove a few segments
  while( 1 ) {
    PII& e = *(--theElements.end( ));
    if( theLength-abs(e.second)<len )
      break;
    theElements.pop_back( );
    theLength -= abs(e.second);
  }
  
  // remove a part of the last segment
  PII& e = *(--theElements.end( ));
  e.second -= e.second>0 ? theLength-len : len-theLength;
  theLength = len;
}


//=========================================================


void PowerWordRep::terminalSegment( int len )
{
  int l = len < theElements.size( ) ? len : theElements.size( );
  for( int t=0 ; t<l ; ++t )
    theElements.pop_front( );
}

//=========================================================


void PowerWordRep::cyclicallyReduce( )
{
  while( theElements.size( )>1 ) {

    list< PII >::iterator it1 = theElements.begin( );
    list< PII >::iterator it2 = theElements.end( );
    it2--;

    if( (*it1).first==(*it2).first ) {
      if( (*it1).second+(*it2).second==0 ) {
	theLength -= 2*abs((*it1).second);
	theElements.pop_front( );
	theElements.pop_back( );
      } else {
	if( ((*it1).second<0) == ((*it2).second<0) )
	  break;
	if( abs((*it1).second)>abs((*it2).second) ) {
	  theLength -= 2*abs((*it2).second);
	  (*it1).second += (*it2).second;
	  theElements.pop_back( );
	} else {
	  theLength -= 2*abs((*it1).second);
	  (*it2).second += (*it1).second;
	  theElements.pop_front( );
	}
      }
    } else 
      break;
  }
}


//=========================================================


void PowerWordRep::cyclicallyReduce( PowerWordRep& conjugator )
{
  conjugator.theElements.clear( );
  
  while( theElements.size( )>1 ) {

    list< PII >::iterator it1 = theElements.begin( );
    list< PII >::iterator it2 = theElements.end( );
    it2--;

    if( (*it1).first==(*it2).first ) {
      if( (*it1).second+(*it2).second==0 ) {
	theLength -= 2*abs((*it1).second);
	conjugator.pushGeneratorFront( (*it2).first , (*it2).second );
	theElements.pop_front( );
	theElements.pop_back( );
      } else {
	if( ((*it1).second<0) == ((*it2).second<0) )
	  break;
	if( abs((*it1).second)>abs((*it2).second) ) {
	  theLength -= 2*abs((*it2).second);
	  conjugator.pushGeneratorFront( (*it2).first , (*it2).second );
	  (*it1).second += (*it2).second;
	  theElements.pop_back( );
	} else {
	  theLength -= 2*abs((*it1).second);
	  conjugator.pushGeneratorFront( (*it2).first , -(*it1).second );
	  (*it2).second += (*it1).second;
	  theElements.pop_front( );
	}
      }
    } else 
      break;
  }
}


//=========================================================


int PowerWordRep::getPower( PowerWordRep& base ) const
{
  int len = theElements.size( );
  if( len<=1 ) return 1;
  
  int up = sqrt(len)+1;
  
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
	up = sqrt(len)+1;
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
  vector< PII > wrd_vct(len);
  list< PII >::const_iterator g_it = theElements.begin( );
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
	PII e = wrd_vct[offset];
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

