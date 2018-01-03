// Copyright (C) 2007 Alexander Ushakov
// Contents: Implementation of class StraightLineProgramWord
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//

#define DEBUG_WP_AUT


#include "StraightLineProgramWord.h"

#include <algorithm>
#include "set"
using namespace std;


//---------------------------------------------------------------------------//
//---------------------------- truncateVertexLeft ---------------------------//
//---------------------------------------------------------------------------//


int StraightLineProgramWord::truncateVertexLeft( int n , LongInteger length )
{
  // cout << "   n = " << n << ", " << "l = " << length << endl;
  
  // If the word is empty or there is nothing to truncate
  if( n==0 || length<=0 )
    return n;
  
  // If the word is a generator (here \f$length\ne 0\f$)
  if( abs(n)<=theTerminals )
    return 0;
  
  // Get the current rule
  Production pr = theRules[abs(n)];
  LongInteger L = pr.theLength;
  
  // if length is at least the length of the word
  if( length>=L )
    return 0;
  
  // If the word = another word
  if( pr.theTerm1!=0 && pr.theTerm2==0 )
    return n>0 ? truncateVertexLeft( pr.theTerm1 , length ) : truncateVertexLeft( -pr.theTerm1 , length );
  
  // If the word is a product of two words
  int T1 = pr.theTerm1;
  int T2 = pr.theTerm2;
  LongInteger L1 = length_rule(T1);
  LongInteger L2 = length_rule(T2);
  if( n<0 ) {
    if( length<L2 ) {
      int t = truncateVertexLeft( -T2 , length );
      int new_rl = max_rule_number( )+1;
      theRules[new_rl] = Production( t , -T1 , L1+L2-length , 1+max( height_rule(t) ,height_rule(T1) ) );
      return new_rl;
    } else if( length==L2 ) {
      return -T1;
    } else {
      return truncateVertexLeft( -T1 , length-L2 );
    }
  } else {
    if( length<L1 ) {
      int t = truncateVertexLeft( T1 , length );
      int new_rl = max_rule_number( )+1;
      theRules[new_rl] = Production( t , T2 , L1+L2-length , 1+max( height_rule(t) , height_rule(T2) ) );
      return new_rl;
    } else if( length==L1 ) {
      return T2;
    } else {
      return truncateVertexLeft( T2 , length-L1 );
    }
  }
  
  return 0;
}


//---------------------------------------------------------------------------//
//---------------------------- truncateVertexRight --------------------------//
//---------------------------------------------------------------------------//


int StraightLineProgramWord::truncateVertexRight( int n , LongInteger length )
{
  // cout << "   n = " << n << ", " << "l = " << length << endl;
  
  // If the word is empty or there is nothing to truncate
  if( n==0 || length<=0 )
    return n;
  
  // If the word is a generator (here \f$length> 0\f$)
  if( abs(n)<=theTerminals )
    return 0;
  
  // Get the current rule
  Production pr = theRules[abs(n)];
  LongInteger L = pr.theLength;
  
  // if length is at least the length of the word
  if( length>=L )
    return 0;
  
  // If the word = another word
  if( pr.theTerm1!=0 && pr.theTerm2==0 )
    return n>0 ? truncateVertexRight( pr.theTerm1 , length ) : truncateVertexRight( -pr.theTerm1 , length );
  
  // If the word is a product of two words
  int T1 = pr.theTerm1;
  int T2 = pr.theTerm2;
  LongInteger L1 = length_rule(T1);
  LongInteger L2 = length_rule(T2);
  if( n<0 ) {
    if( length<L1 ) {
      int t = truncateVertexRight( -T1 , length );
      int new_rl = max_rule_number( )+1;
      theRules[new_rl] = Production( -T2 , t , L1+L2-length , 1+max( height_rule(t) , height_rule(T2) ) );
      return new_rl;
    } else if( length==L1 ) {
      return -T2;
    } else {
      return truncateVertexRight( -T2 , length-L1 );
    }
  } else {
    if( length<L2 ) {
      int t = truncateVertexRight( T2 , length );
      int new_rl = max_rule_number( )+1;
      theRules[new_rl] = Production( T1 , t , L1+L2-length , 1+max( height_rule(t) , height_rule(T1) ) );
      return new_rl;
    } else if( length==L2 ) {
      return T1;
    } else {
      return truncateVertexRight( T1 , length-L2 );
    }
  }
  
  return 0;
}


//---------------------------------------------------------------------------//
//------------------------------ leftGCDLength ------------------------------//
//---------------------------------------------------------------------------//


LongInteger StraightLineProgramWord::leftGCDLength( int n1 , const StraightLineProgramWord& CS , int n2 ) const
{
  LongInteger L1 = length_rule( n1 );
  LongInteger L2 = CS.length_rule( n2 );
#ifdef DEBUG_WP_AUT
  cout << "L1 = " << L1 << endl;
  cout << "L2 = " << L2 << endl;
  cout << "n1 = " << n1 << endl;
  cout << "n2 = " << n2 << endl;
#endif

  LongInteger a = 0;
  LongInteger b = L1<L2 ? L1 : L2;
  
  StraightLineProgramWord SLP1 = *this;
  StraightLineProgramWord SLP2 = CS;
  while( a!=b ) {
    
#ifdef DEBUG_WP_AUT
    cout << "Test (" << a << "," << b << ")" << endl;
#endif
    if( b-a==1 ) {
      int v1 = SLP1.truncateVertexRight( n1 , L1-b );
      int v2 = SLP2.truncateVertexRight( n2 , L2-b );
#ifdef DEBUG_WP_AUT
      cout << "U1 = " << SLP1.getWord( v1 ) << endl;
      cout << "U2 = " << SLP2.getWord( v2 ) << endl;
#endif
      int e = SLP1.equal( v1 , SLP2 , v2 );
#ifdef DEBUG_WP_AUT
      cout << "e = " << e << endl;
#endif
      return e ? b : a;
    }
    
    LongInteger c = (a+b)/2;
#ifdef DEBUG_WP_AUT
    cout << "c = " << c << endl;
#endif
    int v1 = SLP1.truncateVertexRight( n1 , L1-c );
    int v2 = SLP2.truncateVertexRight( n2 , L2-c );
#ifdef DEBUG_WP_AUT
    cout << "   v1 = " << SLP1.getWord( v1 ) << endl;
    cout << "   v2 = " << SLP2.getWord( v2 ) << endl;
    cout << SLP1 << endl;
    cout << SLP2 << endl;
#endif
    bool e = SLP1.equal( v1 , SLP2 , v2 );
    if( e ) a = c;
    else    b = c;
    // if( e ) cout << "SUCCESS" << endl;
    // else    cout << "FAILURE" << endl;
  }
  return b;
}


//---------------------------------------------------------------------------//
//--------------------------- cancellationLength ----------------------------//
//---------------------------------------------------------------------------//


LongInteger StraightLineProgramWord::cancellationLength( int n1 , int n2 ) const
{
  return leftGCDLength( -n1 , *this , n2 );
}

