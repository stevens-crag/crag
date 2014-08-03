// Copyright (C) 2005 Alexander Ushakov
// Contents: Implementation of class Permutation
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//


#include "RanlibCPP.h"
#include "Permutation.h"
#include <list>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <cassert>

using namespace std;


//---------------------------------------------------------------------------//
//----------------------------- Permutation ---------------------------------//
//---------------------------------------------------------------------------//



Permutation::Permutation( int l ) : theValue( l )
{
  for( int i=0 ; i<l ; ++i )
    theValue[i] = i;
}

Permutation::Permutation( int size , const int* p ) :
  theValue( size )
{
  for( int i=0 ; i<size ; ++i )
    theValue[i] = p[i];
}

Permutation::Permutation( const vector< int >& p ) :
  theValue( p )
{ 
  
}


Permutation Permutation::power( int p ) const
{
	Permutation base = p<0 ? -*this : *this;
	int sz = theValue.size( );
	Permutation result( sz );

	int ap = abs(p);
	while( ap ) {
		if( ap%2==1 ) 
			result *= base;
		ap = ap>>1;
		base = base*base;
	}
	return result;
}


Permutation Permutation::CYCLE( int N , const vector< int >& cycle )
{
	Permutation result( N );
	for( int i=1 ; i<cycle.size() ; ++i )
		swap( result.theValue[cycle[i-1]] , result.theValue[cycle[i]] );

	return result;
}

Permutation& Permutation::left_mult_by_cycle( const vector< int >& cycle )
{
	for( int i=1 ; i<cycle.size() ; ++i )
		swap( theValue[cycle[i-1]] , theValue[cycle[i]] );

	return *this;
}

void Permutation::lr_multiply_by_cycles( Permutation& P , Permutation& I , const vector< int >& M1 , const vector< int >& M2 )
{
	for( int i=1 ; i<M1.size() ; ++i ) {
		int a = M1[i-1];
		int b = M1[i];
		int c = P.theValue[a];
		int d = P.theValue[b];
		swap( P.theValue[a] , P.theValue[b] );
		swap( I.theValue[c] , I.theValue[d] );
	}

	for( int i=M2.size()-1 ; i>0 ; --i ) {
		int a = M2[i-1];
		int b = M2[i];
		int c = I.theValue[a];
		int d = I.theValue[b];
		swap( I.theValue[a] , I.theValue[b] );
		swap( P.theValue[c] , P.theValue[d] );
	}
}


//---------------------------------------------------------------------------//
//------------------- computeConjugacyClassRepresentative -------------------//
//---------------------------------------------------------------------------//


Permutation
Permutation::computeConjugacyClassRepresentative( Permutation& conj ) const
{
  int i;
  int size = theValue.size( );
  Permutation p1( *this );
  conj = Permutation( size );

  // 1. arrange cycles
  set< pair< int,int > > cycles;
  for( i=0 ; i<size ; ++i ) {
    
    if( p1.theValue[i]!=i ) {
      
      int r = i;
      do {
	
	if( p1.theValue[r]!=r+1 ) {

	  Permutation c( size );
	  c.change( r+1 , p1.theValue[r] );
	  
	  p1 = c * p1 * c;
	  conj *= c;
	}
	r = p1.theValue[r];
      } while( p1.theValue[r]!=i );
      cycles.insert( pair< int,int >( r+1-i , i ) );
      i = r;
    } else 
      cycles.insert( pair< int,int >( 1 , i ) );
  }
  
  set< pair< int,int > >::iterator c_it = cycles.begin( );
  
  int pos = 0;
  while( cycles.size( ) ) {

    c_it = cycles.begin( );
    int len = (*c_it).first;
    int pt = (*c_it).second;
    cycles.erase( c_it );
    
    if( pt!=pos ) {
      
      // 2.1. prepare shift conjugators
      Permutation c1(size);
      Permutation c2(size);
      Permutation c3(size);

      for( i=pos ; i<pt+len ; ++i )
	c1.theValue[i] = pt+len-1 - (i-pos);

      for( i=pos ; i<pos+len ; ++i )
	c2.theValue[i] = pos+len-1 - (i-pos);
      
      for( i=pos+len ; i<pt+len ; ++i )
	c3.theValue[i] = pt+len-1 - (i-pos-len);
      
      // 2.2. shift cycle to a position pos
      p1 = c1 * p1 * c1.inverse( );
      p1 = c2 * p1 * c2.inverse( );
      p1 = c3 * p1 * c3.inverse( );
      conj *= c1.inverse( );
      conj *= c2.inverse( );
      conj *= c3.inverse( );
      
      set< pair< int,int > > _cycles;
      c_it = cycles.begin( );
      for( ; c_it!=cycles.end( ) ; ++c_it ) {
	if( (*c_it).second<pt )
	  _cycles.insert( pair< int,int >( (*c_it).first , (*c_it).second+len ) );
	else
	  _cycles.insert( pair< int,int >( (*c_it).first , (*c_it).second ) );
      }
      cycles = _cycles;
    }
    pos += len;
  }
  
  return p1;
}


//---------------------------------------------------------------------------//
//-------------------------- computeConjugator ------------------------------//
//---------------------------------------------------------------------------//




Permutation
Permutation::computeConjugator( const Permutation& p ) const
{
  int size = theValue.size( );
  assert(size == p.theValue.size( ));
//  if( size!=p.theValue.size( ) ) Do something with this!!!
//    return false;
  
  Permutation result( size );
  
  Permutation conj1( size );
  Permutation p1 = computeConjugacyClassRepresentative( conj1 );

  Permutation conj2( size );
  Permutation p2 = p.computeConjugacyClassRepresentative( conj2 );
  
  // cout << "_p1 = " << p1 << endl;
  // cout << "_p2 = " << p2 << endl;
  
  result = conj1 * conj2.inverse( );
  return result;
}

/*
Permutation
Permutation::computeConjugator( const Permutation& p ) const
{
  int size = theValue.size( );
  if( size!=p.theValue.size( ) )
    return false;
  
  Permutation result( size );
  
  vector< bool > h( size , false );
  for( int i=size-1 ; i>=0 ; --i ) {
    int r = i;
    int s = i;
    while( !h[r] ) {
      h[r] = true;
      r = theValue[r];
      s = p.theValue[s];
      if( r!=s )
	result.theValue[r] = s;
    }
  }
  
  // cout << " ---- " << result << endl;
  
  vector< bool > k( size , false );
  for( int i=size-1 ; i>=0 ; --i ) {
    if( !k[i] ) {
      k[i] = true;
      int r = i;
      while( !k[result.theValue[r]] && result.theValue[r]!=r ) {
	r = result.theValue[r];
	k[r] = true;
      }
      
      result.theValue[r] = i;
    }
  }
  
  // cout << " ---- " << result << endl;
  
  return result;
}
*/

/*
bool 
Permutation::computeConjugator( const Permutation& p , Permutation& res ) const
{
  Permutation _p( *this );

  int size = _p.theValue.size( );
  if( size!=p.theValue.size( ) )
    return false;
  
  res = Permutation( size );

  cout << "------------------" << endl;
  cout << _p << endl;
  cout << p << endl;
  
  for( int i=0 ; i<size ; ++i ) {
    
    int y = p.theValue[i];
    Permutation inv = _p.inverse( );
    int k = inv[y];
    if( k==i )
      continue;
    
    cout << "i = " << i << endl;
    cout << "k = " << k << endl;
    Permutation c( size );
    c.change( i , k );
    
    _p = c * _p * c;
    
    if( _p.theValue[i]!=y )
      return false;
    cout << "------------------" << endl;
    cout << _p << endl;
    cout << p << endl;
    
    res *= c;
  }
  
  return true;
}
*/

vector<int> Permutation::getWordPresentation( ) const
{
  int j;
  vector<int> result;

  /*
  // first implementation of this function
  // much worse than the second one
  vector<int> val( theValue );
  for( int i=0 ; i<val.size() ; ++i ) {
    if( val[i]!=i ) {
      int t = val[i];
      int s = i;
      
      for( int j=t-1 ; j>s ; --j )
	result.push_back( -j-1 );
      for( int j=s ; j<t ; ++j )
	result.push_back( j+1 );
      
      swap( val[i] , val[val[i]] );
      i--;
    }
  }
  */


  // second implementation  
  vector<int> val( theValue );
  Permutation inv = inverse( );
  vector<int> inv_val( inv.theValue );
  for( int i=0 ; i<theValue.size( ) ; ++i ) {
    if( val[i]!=i ) {
      
      int r = val[i];
      int t = inv_val[i];
      int s = i;
      
      for( j=t-1 ; j>s ; --j )
	result.push_back( -j-1 );
      for( j=s ; j<t ; ++j )
	result.push_back( j+1 );
      
      swap( val[s] , val[t] );
      swap( inv_val[s] , inv_val[r] );
      
    }
  }
  
  return result;
}


vector<int> Permutation::geodesic( ) const
{
  int i;
  vector<int> result;
  
  Permutation cur( theValue.size( ) );
  Permutation inv( theValue.size( ) );
  
  for( i=0 ; i<theValue.size( ) ; ++i ) {
    
    int pos = inv.theValue[theValue[i]];
    for( int j=pos-1 ; j>=i ; --j ) {
      result.push_back( j );
      inv.change( cur.theValue[j] , cur.theValue[j+1] );
      cur.change( j , j+1 );
    }
  }

  for( i=0 ; i<result.size( )/2 ; ++i )
    swap( result[i] , result[result.size( )-i-1] );
  
  return result;
}

vector< int > Permutation::geodesicWord( ) const
{
  vector<int> V = geodesic( );
  for( int i=0 ; i<V.size( ) ; ++i ) ++V[i];
  return V;
}

/*
void Permutation::change( int i , int j )
{
  int len = i<j ? j : i;
  
  if( len<=theValue.size() )
    for( int t=theValue.size() ; t<=len ; ++t ) 
      theValue.push_back( t );

  swap( theValue[i] , theValue[j] );
}
*/


bool Permutation::isTrivial( ) const
{
  for( int i=0 ; i<theValue.size( ) ; ++i )
    if( theValue[i]!=i )
      return false;
  return true;
}

bool Permutation::operator == ( const Permutation& p ) const 
{
  return theValue==p.theValue;
}

bool Permutation::operator != ( const Permutation& p ) const 
{
  return theValue!=p.theValue;
}

bool Permutation::operator < ( const Permutation& p ) const 
{
  if( theValue.size( )<p.theValue.size( ) )
    return true;
  if( theValue.size( )>p.theValue.size( ) )
    return false;

  for( int t=0 ; t<theValue.size( ) ; ++t ) {
    if( theValue[t]<p.theValue[t] )
      return true;
    if( theValue[t]>p.theValue[t] )
      return false;
  }
  return false;
}


Permutation Permutation::operator * ( const Permutation& p ) const 
{
  return (Permutation(*this) *= p);
}


Permutation& Permutation::operator *= ( const Permutation& p ) 
{
  int l1 = theValue.size( );
  int l2 = p.theValue.size( );
  
  if( theValue.size()<p.theValue.size() ) {
    for( int i=l1 ; i<l2 ; ++i ) 
      theValue.push_back( i );
  }
  
  for( int t=0 ; t < theValue.size( ) ; ++t ) {
    if (theValue[t] < p.theValue.size()) {
      theValue[t] = p.theValue[theValue[t]];
    }
  }
  return *this;
}


Permutation Permutation::inverse( ) const 
{
  int len = theValue.size( );
  Permutation res( len );
  
  for( int t=0 ; t<len ; ++t )
    res.theValue[theValue[t]] = t;
  return res;
}

int Permutation::difference( const Permutation& p ) const
{
	if( theValue.size( )!=p.theValue.size( ) )
		return 999999;

	int result = 0;
	for( int i=0 ; i<theValue.size( ) ; ++i )
		if( theValue[i]!=p.theValue[i] )
			++result;

	return result;
}


Permutation Permutation::random( int l ) 
{
  Permutation res(l);
  
  for( int t=0 ; t<l-1 ; ++t ) {
    // int pos = rand( )%(l-t);
		int pos = RandLib::ur.irand( 0 , l-t-1 );
    swap( res.theValue[t] , res.theValue[t+pos] );
  }
  
  return res;
}

ostream& operator << ( ostream& os , const Permutation& p ) 
{
  int len = p.theValue.size( );
  
  os << "{";
  for( int t=0 ; t<len ; ++t ) {
    if( t>0 )
      os << ",";
    os << (int)p.theValue[t];
  }
  os << "}";
  return os;
}


Permutation Permutation::join2( const Permutation& p ) const
{
  int i;
  int l1 = theValue.size( );
  int l2 = p.theValue.size( );
  if( l1!=l2 ) {
    cerr << "Check dimensions in join2 operation" << endl;
    exit( 1 );
  }
  

  /*
  Permutation omega( l1 );
  for( int i=0 ; i<l1 ; ++i )
    omega.theValue[i] = (i+l1-1)%l1;
  
  Permutation p1 = inverse( )*omega;
  Permutation p2 = p.inverse( )*omega;
  Permutation p3 = p1.meet( p2 );

  return p3.inverse( )*omega;
  */

  vector< int > cycleN( l1 , -1 );
  vector< int > foundPts( l1, -1 );
  for( i=l1-1 ; i>=0 ; --i ) {
    
    if( cycleN[i]!=-1 ) 
      continue;
    
    list< int > toCheck;
    toCheck.push_back( i );
    foundPts[i] = i;
    while( toCheck.begin( )!=toCheck.end( ) ) {
      int cur = *toCheck.begin( );
      toCheck.pop_front( );
      cycleN[cur] = i;
      
      int next1 = theValue[cur];
      int next2 = p.theValue[cur];
      
      if( foundPts[next1]!=i ) {
	foundPts[next1] = i;
	toCheck.push_back( next1 );
      }
      if( foundPts[next2]!=i ) {
	foundPts[next2] = i;
	toCheck.push_back( next2 );
      }
    }
  }

  vector< int > perm( l1 , -1 );

  vector< int > firstPt( l1 , -1 );
  vector< int > prevPt ( l1 , -1 );
  for( i=l1-1 ; i>=0 ; --i ) {
    int c = cycleN[i];
    perm[i] = i;
    if( firstPt[c]==-1 )
      firstPt[c] = i;
    else
      swap( perm[i] , perm[prevPt[c]] );
    prevPt[c] = i;
  }  
  
  return perm;
}


Permutation Permutation::meet2( const Permutation& p ) const
{
  int i;
  int l1 = theValue.size( );
  int l2 = p.theValue.size( );
  if( l1!=l2 ) {
    cerr << "Check dimensions in meet2 operation" << endl;
    exit( 1 );
  }

  // 1. extract cycles from the second permutation
  vector< int > cycle2N( l1 , 0 );
  for( i=0 ; i<l1 ; ++i ) {
    if( p.theValue[i]!=i && cycle2N[i]==0 ) {
      cycle2N[i] = i+1;
      for( int t=p.theValue[i] ; t!=i ; t=p.theValue[t] )
	cycle2N[t] = i+1;
    }
    // cout << cycle2N[i] << ",";
  }
  // cout << endl;
  
  // 2. extract cycles from the first permutation
  //    and compute the result
  Permutation result( l1 );
  vector< int > cycle1N( l1 , 0 );
  
  for( i=0 ; i<l1 ; ++i ) {
    if( theValue[i]!=i && cycle1N[i]==0 ) {
      
      cycle1N[i] = i+1;
      vector< pair<int,int> > pairs;
      if( cycle2N[i] )
	pairs.push_back( pair<int,int>( cycle2N[i] , i ) );
      for( int t=theValue[i] ; t!=i ; t=theValue[t] ) {
	cycle1N[t] = i+1;
	if( cycle2N[t] )
	  pairs.push_back( pair<int,int>( cycle2N[t] , t ) );
      }

      if( pairs.size( )<1 )
	continue;

      // cout << "> " << result << endl;
      prepare_pairs( l1+1 , result , pairs );
      // cout << "< " << result << endl;
    }
    // cout << cycle1N[i] << ",";
  }
  // cout << endl;
  
  return result;
}


void
Permutation::prepare_pairs
( int N, 
  Permutation& P,
  vector< pair<int,int> >& pairs 
  ) const
{
  
  set< pair<int,int> > pairs1;
  for( int i=0 ; i<pairs.size( ) ; ++i )
    pairs1.insert( pairs[i] );
  
  set< pair<int,int> >::reverse_iterator it = pairs1.rbegin( );
  int len = 1;
  //int beg = (*it).second;
  pair< int , int > prev_pair = (*it);
  it++;
  for( ; it!=pairs1.rend( ) ; ++it ) {
    // cout << "(" << (*it).first << "," << (*it).second << ")" << endl;
    
    if( prev_pair.first==(*it).first ) {
      P.change( prev_pair.second , (*it).second );
      // cout << prev_pair.second << "," << (*it).second << endl;
    } else {
      len = 0;
      //beg = (*it).second;
    }
    len++;
    prev_pair = (*it);
  }
  
  /*
  vector< int > counts( N , 0 );
  for( int i=0 ; i<pairs.size( ) ; ++i )
    counts[pairs[i].first]++;

  for( int i=1 ; i<N ; ++i )
    counts[i] += counts[i-1];

  vector< pair<int,int> > pairs1( pairs.size( ) );
  vector< int > cur_num( N , 0 );
  for( int i=0 ; i<pairs.size( ) ; ++i ) {
    pairs1[--counts[pairs[i].first]] = pairs[i];
    // counts[triples[i].c1]--;
  }
  
  for( int i=1 ; i<N ; ++i ) {
    if( count[i]-count[i-1]>1 ) {
      vector< int > counts( N , 0 );
      
    }
  }
  pairs = pairs1;
  */
}


Permutation Permutation::getHalfTwistPermutation( int size )
{
  Permutation result( size );
  for( int i=0 ; i<size ; ++i )
    result[i] = size-i-1;
  return result;
}


//---------------------------------------------------------------------------//
//------------------------------ RightGCD -----------------------------------//
//---------------------------------------------------------------------------//


Permutation Permutation::RightGCD( const Permutation& p ) const
{
  int l1 = theValue.size( );
  int l2 = p.theValue.size( );
  if( l1!=l2 ) {
    cerr << "Check dimensions in meet operation" << endl;
    exit( 1 );
  }

  int* l_ind_a = new int[l1];
  int* l_ind_b = new int[l1];
  int* r_ind_a = new int[l1];
  int* r_ind_b = new int[l1];
  
  Permutation result( l1 );
  _sub_meet( p, inverse(), p.inverse(), result, l_ind_a, l_ind_b, r_ind_a , r_ind_b , 0 , l1 );
  
  delete[] l_ind_a;
  delete[] l_ind_b;
  delete[] r_ind_a;
  delete[] r_ind_b;
  
  return result;
}


//---------------------------------------------------------------------------//
//------------------------------ RightLCM -----------------------------------//
//---------------------------------------------------------------------------//


Permutation Permutation::RightLCM( const Permutation& p ) const
{
  Permutation Delta = getHalfTwistPermutation( size( ) );
  return (*this * Delta).RightGCD( p*Delta ) * Delta;
}


//---------------------------------------------------------------------------//
//------------------------------- LeftGCD -----------------------------------//
//---------------------------------------------------------------------------//


Permutation Permutation::LeftGCD( const Permutation& p ) const
{
  return -((-*this).RightGCD( -p ));
}


//---------------------------------------------------------------------------//
//------------------------------- LeftLCM -----------------------------------//
//---------------------------------------------------------------------------//


Permutation Permutation::LeftLCM( const Permutation& p ) const
{
  Permutation Delta = getHalfTwistPermutation( size( ) );
  Permutation P3 = Delta * inverse( );
  Permutation P4 = Delta * p.inverse( );
  return (Delta*P3.RightGCD( P4 )).inverse( );
}


//---------------------------------------------------------------------------//
//------------------------------ _sub_meet ----------------------------------//
//---------------------------------------------------------------------------//


void Permutation::_sub_meet ( const Permutation& p , 
			      const Permutation& ip1 ,
			      const Permutation& ip2 ,
			      Permutation& cur , 
			      int* l_ind_a ,
			      int* l_ind_b ,
			      int* r_ind_a ,
			      int* r_ind_b ,
			      int beg , int end ) const
{
  int i;
  if( end-beg==1 )
    return;

  const int A = beg;
  const int B = beg+(end-beg)/2;
  const int C = end;
  
  // I. reorder left and right parts of permutation according to meet operation
  _sub_meet( p, ip1, ip2, cur, l_ind_a, l_ind_b, r_ind_a , r_ind_b , A , B );
  _sub_meet( p, ip1, ip2, cur, l_ind_a, l_ind_b, r_ind_a , r_ind_b , B , C );
  
  
  // II. merge left and right parts
  
  // 1. find left indeces
  for( i=B-1 ; i>=A; --i ) {
    int a = ip1.theValue[cur.theValue[i]];
    int b = ip2.theValue[cur.theValue[i]];
    if( i==B-1 ) {
      l_ind_a[i] = a;
      l_ind_b[i] = b;
    } else {
      l_ind_a[i] = a<l_ind_a[i+1] ? a : l_ind_a[i+1];
      l_ind_b[i] = b<l_ind_b[i+1] ? b : l_ind_b[i+1];
    }
  }
  
  // 2. find right indeces
  for( i=B ; i<C; ++i ) {
    int a = ip1.theValue[cur.theValue[i]];
    int b = ip2.theValue[cur.theValue[i]];
    if( i==B ) {
      r_ind_a[i] = a;
      r_ind_b[i] = b;
    } else {
      r_ind_a[i] = a<r_ind_a[i-1] ? r_ind_a[i-1] : a;
      r_ind_b[i] = b<r_ind_b[i-1] ? r_ind_b[i-1] : b;
    }
  }

  // 3. merge lists
  int i1 = A;
  int i2 = B;
  int* new_sublist = new int[C-A];
  for( i=0 ; i<C-A ; ++i ) {
    if( i1==B ) {
      new_sublist[i] = cur.theValue[i2++];
      continue;
    }
    if( i2==C ) {
      new_sublist[i] = cur.theValue[i1++];
      continue;
    }
    
    if( l_ind_a[i1]>r_ind_a[i2] && 
	l_ind_b[i1]>r_ind_b[i2] )
      new_sublist[i] = cur.theValue[i2++];
    else
      new_sublist[i] = cur.theValue[i1++];
  }
  for( i=0 ; i<C-A ; ++i )
    cur.theValue[A+i] = new_sublist[i];

  delete[] new_sublist;
}


Permutation Permutation::tinyFlip( int sh ) const
{
  int len = theValue.size( );
  Permutation result( len );
  
  if( sh<0 )
    sh = sh-(sh/len-1)*len;
  
  for( int t=0 ; t<len ; ++t )
    result.theValue[(t+sh)%len] = (theValue[t]+sh)%len;
  
  return result;
}


//---------------------------------------------------------------------------//
//------------------------------- mixable -----------------------------------//
//---------------------------------------------------------------------------//


bool Permutation::mixable( const Permutation& p1 , const Permutation& p2 )
{
  int l1 = p1.theValue.size( );
  int l2 = p2.theValue.size( );
  if( l1!=l2 ) {
    cerr << "Check dimensions in mixable operation" << endl;
    exit( 1 );
  }
  
  Permutation ip1 = p1.inverse( );
  for( int i=1 ; i<l1 ; ++i )
    if( ip1.theValue[i-1]>ip1.theValue[i] && p2.theValue[i-1]<p2.theValue[i] )
      return false;
  
  return true;
}


//---------------------------------------------------------------------------//
//-------------------------------- flip -------------------------------------//
//---------------------------------------------------------------------------//


Permutation Permutation::flip( ) const
{
  int sz = size( );
  Permutation result(sz);
  
  for( int i=0 ; i<sz ; ++i )
    result.theValue[sz-i-1] = sz-theValue[i]-1;
  
  return result;
}


//---------------------------------------------------------------------------//
//-------------------------------- length -----------------------------------//
//---------------------------------------------------------------------------//


int Permutation::length( ) const
{
  return geodesic( ).size( );
}


//---------------------------------------------------------------------------//
//------------------------------ increaseSize -------------------------------//
//---------------------------------------------------------------------------//


Permutation Permutation::increaseSize( int N ) const
{
  int sz = theValue.size( );
  if( N<=sz )
    return *this;
  
  Permutation result = *this;
  result.theValue.resize( N );
  for( int i=sz ; i<N ; ++i )
    result.theValue[i] = i;
  return result;
}


//---------------------------------------------------------------------------//
//--------------------------- getCyclePermutation ---------------------------//
//---------------------------------------------------------------------------//


Permutation Permutation::getCyclePermutation( int size )
{
  Permutation result;
  
  for( int i=0 ; i<size-1 ; ++i )
    result.theValue[i] = i+1;
  result.theValue[size-1] = 0;
  
  return result;
}
