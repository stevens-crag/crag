// Copyright (C) 2007 Alexander Ushakov
// Contents: Implementation of class FreeMetabelianGroupAlgorithms
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//

#include "Map.h"
#include "FreeMetabelianGroupAlgorithms.h"




//---------------------------------------------------------------------------//
//------------------------------- printEdgeMap ------------------------------//
//---------------------------------------------------------------------------//


//! Get a word given by a vector
Word getTailWord( int N , const vector< int >& T )
{
  Word C;
  for( int i=0 ; i<N ; ++i )
    C *= Word(i+1).power( T[i] );
  return C;
}


//---------------------------------------------------------------------------//
//------------------------------- printEdgeMap ------------------------------//
//---------------------------------------------------------------------------//


void printEdgeMap( const map< vector< int > , int >& EM )
{
  for( map< vector< int > , int >::const_iterator E_it=EM.begin( ) ; E_it!=EM.end( ) ; ++E_it ) {
    pair< vector< int > , int > C = *E_it;
    for( int i=0 ; i<C.first.size( ) ; ++i )
      cout << C.first[i] << " ";
    cout << " -> " << C.second << endl;
  }
}


//---------------------------------------------------------------------------//
//--------------------------------- getTail ---------------------------------//
//---------------------------------------------------------------------------//


vector< int > getTail( int N , const Word& w )
{
  vector< int > V( N , 0 );
  for( Word::const_iterator w_it=w.begin( ) ; w_it!=w.end( ) ; ++w_it ) {
    int  x = *w_it;
    int ax = abs(x);
    V[ax-1] += x>0 ? 1:-1;
  }
  return V;
}


//---------------------------------------------------------------------------//
//-------------------------------- getEdgeMap -------------------------------//
//---------------------------------------------------------------------------//


map< vector< int > , int > getEdgeMap( int N , const Word w )
{
  vector< int > V( N , 0 );
  map< vector< int > , int > E;
  for( Word::const_iterator w_it=w.begin( ) ; w_it!=w.end( ) ; ++w_it ) {
    int  x = *w_it;
    int ax = abs(x);
    if( x>0 ) {
      vector< int > W = V;
      W.push_back( x );
      E[W]++;
    }
    V[ax-1] += x>0 ? 1:-1;
    if( x<0 ) {
      vector< int > W = V;
      W.push_back( -x );
      E[W]--;
    }
  }
  
  return E;
}


//---------------------------------------------------------------------------//
//--------------------------- dropTrivialEdgesInMap -------------------------//
//---------------------------------------------------------------------------//


map< vector< int > , int > dropTrivialEdgesInMap( const map< vector< int > , int >& EM )
{
  map< vector< int > , int > result;

  for( map< vector< int > , int >::const_iterator E_it=EM.begin( ) ; E_it!=EM.end( ) ; ++E_it )
    if( (*E_it).second!=0 )
      result.insert( *E_it );
  
  return result;
}


//---------------------------------------------------------------------------//
//----------------------------- getModConjugator ----------------------------//
//---------------------------------------------------------------------------//


//! Compute a conjugator for a given edge which shifts it into the strip specified by \f$S\f$.
/*!
  Here \f$S\f$ is a tail vector for \f$w\f$.
*/
void getModEdgeConjugator( int N, const vector< int >& S, vector< int > edge , int arity , map< vector< int > , int >& EM , bool positive )
{
  int direction = edge[N];
  if( direction==1 )
    return;

  arity = positive ? arity:-arity;

  while( edge[0]<0 ) {

    // left edge
    EM[edge] += arity;
    
    // lower path
    edge[N] = 1;
    int x1 = edge[0];
    for( int i=x1 ; edge[0]<0 ; ++edge[0] )
      EM[edge] -= arity;
    edge[N] = direction;

    // right edge
    EM[edge] -= arity;
    
    // upper path
    edge[N] = 1;
    edge[direction-1]++;
    for( int i=0 ; edge[0]>x1 ; ) {
      --edge[0];
      EM[edge] += arity;
    }
    edge[direction-1]--;
    edge[N] = direction;
    
    for( int i=0 ; i<N ; ++i )
      edge[i] += S[i];
  }

  int u1 = S[0];
  while( edge[0]>u1 ) {

    for( int i=0 ; i<N ; ++i )
      edge[i] -= S[i];
    
    // right edge
    EM[edge] -= arity;
    
    // lower path
    edge[N] = 1;
    int x1 = edge[0];
    for( ; edge[0]>0 ; ) {
      --edge[0];
      EM[edge] -= arity;
    }
    edge[N] = direction;

    // left edge
    EM[edge] += arity;

    // upper path
    edge[N] = 1;
    edge[direction-1]++;
    for( ; edge[0]<x1 ; ++edge[0] )
      EM[edge] += arity;
    edge[direction-1]--;
    edge[N] = direction;
  }
}


//---------------------------------------------------------------------------//
//----------------------------- getModConjugator ----------------------------//
//---------------------------------------------------------------------------//


//! Compute a conjugator for a given word which "compactifies" the word to its main strip.
/*!
  Here \f$S\f$ is a tail vector for \f$w\f$.
*/
Word getModConjugator( int N , const vector< int >& S , const Word& w )
{
  // Shift all the edges to the main strip
  map< vector< int > , int > EM = dropTrivialEdgesInMap( getEdgeMap( N , w ) );
  map< vector< int > , int > EM_r;
  for( map< vector< int > , int >::const_iterator E_it=EM.begin( ) ; E_it!=EM.end( ) ; ++E_it )
    getModEdgeConjugator( N , S , (*E_it).first , (*E_it).second , EM_r , true );
  return FreeMetabelianGroupAlgorithms::getWordFromEdgeMap( N , dropTrivialEdgesInMap( EM_r ) );
}


//---------------------------------------------------------------------------//
//----------------------------- getLoopConjugator ---------------------------//
//---------------------------------------------------------------------------//


//! Compute a conjugator which belongs to a derived subgroup (assuming such exists)
/*!
  The procedure outputs an element \f$c\f$ such that \f$c^{-1} w1 c = w2 \f$.
*/
Word getLoopConjugator( int N , const vector< int >& S , const Word& w1 , const Word& w2 )
{
  Word result;
  
  Word c1 = getModConjugator( N , S , w1 );
  Word c2 = getModConjugator( N , S , w2 );
  
  map< vector< int > , int > EM1 = dropTrivialEdgesInMap( getEdgeMap( N , -c1*w1*c1 ) );
  map< vector< int > , int > EM2 = dropTrivialEdgesInMap( getEdgeMap( N , -c2*w2*c2 ) );
  map< vector< int > , int > EM;
  for( map< vector< int > , int >::iterator E_it=EM1.begin( ) ; E_it!=EM1.end( ) ; ++E_it ) {
    pair< vector< int > , int > C = *E_it;
    if( C.first[0]==0 && C.first[N]!=1 ) 
      EM[C.first] = C.second;
  }
  for( map< vector< int > , int >::iterator E_it=EM2.begin( ) ; E_it!=EM2.end( ) ; ++E_it ) {
    pair< vector< int > , int > C = *E_it;
    if( C.first[0]==0 && C.first[N]!=1 ) 
      EM[C.first] -= C.second;
  }
  Word c3 = FreeMetabelianGroupAlgorithms::getWordFromEdgeMap( N , dropTrivialEdgesInMap( EM ) );
  return c1 * c3 * -c2;
}


//---------------------------------------------------------------------------//
//------------------------------- shiftEdgeMap ------------------------------//
//---------------------------------------------------------------------------//


map< vector< int > , int > shiftEdgeMap( int N , const vector< int >& S , const map< vector< int > , int >& EM )
{
  map< vector< int > , int > result;
  
  for( map< vector< int > , int >::const_iterator E_it=EM.begin( ) ; E_it!=EM.end( ) ; ++E_it ) {
    pair< vector< int > , int > C = *E_it;
    for( int i=0 ; i<N ; ++i )
      C.first[i] += S[i];
    result.insert( C );
  }

  return result;
}


//---------------------------------------------------------------------------//
//-------------------------------- modEdgeMap -------------------------------//
//---------------------------------------------------------------------------//


//! Take an edge map EM modulo a vector S. This shifts all edges to a strip \f$[0,S[0]) \times \mathbb{Z}^{N-1}\f$.
/*!
  S[0] is assumed to be positive.
*/
map< vector< int > , int > modEdgeMap( int N , const vector< int >& S , const map< vector< int > , int >& EM )
{
  map< vector< int > , int > result;
  int u1 = S[0];
  for( map< vector< int > , int >::const_iterator E_it=EM.begin( ) ; E_it!=EM.end( ) ; ++E_it ) {
    pair< vector< int > , int > C = *E_it;
    int x1 = C.first[0];
    int r = x1>=0 ? x1/u1 : (x1-u1+1)/u1;
    
    for( int i=0 ; i<N ; ++i ) 
      C.first[i] -= r*S[i];
    result[C.first] += C.second;
  }
  
  return result;
}


//---------------------------------------------------------------------------//
//------------------------------ x1shiftEdgeMap -----------------------------//
//---------------------------------------------------------------------------//


//! Take an edge map EM modulo a vector S. This shifts all edges to a strip \f$[0,S[0]) \times \mathbb{Z}^{N-1}\f$.
/*!
  S[0] is assumed to be positive. EM is assumed to be taken mod S.
*/
map< vector< int > , int > x1shiftEdgeMap( int N , const vector< int >& S , const map< vector< int > , int >& EM )
{
  map< vector< int > , int > result;
  int u1 = S[0];
  for( map< vector< int > , int >::const_iterator E_it=EM.begin( ) ; E_it!=EM.end( ) ; ++E_it ) {
    pair< vector< int > , int > C = *E_it;
    if( ++C.first[0]==u1 )
      for( int i=0 ; i<N ; ++i )
	C.first[i] -= S[i];
    result.insert( C );
  }
  
  return result;
}


//---------------------------------------------------------------------------//
//--------------------------- suggestTailConjugator -------------------------//
//---------------------------------------------------------------------------//


//! Suggest a tail-conjugator for two elements
/*!
  Function returns a pair \f$(A,B)\f$ where \f$A\f$ is false 
  if it is clear that EM1 and EM2
  are not conjugate, \f$B\f$ is a vector one should try 
  (only if \f$A\f$ is true).
*/
pair< bool , vector< int > > suggestTailConjugator( int N , const map< vector< int > , int >& EM1 , const map< vector< int > , int >& EM2 )
{
  // 1. Construct the interiors without x1-edges
  map< vector< int > , int > _EM1;
  for( map< vector< int > , int >::const_iterator E_it=EM1.begin( ) ; E_it!=EM1.end( ) ; ++E_it ) {
    pair< vector< int > , int > C = *E_it;
    if( 1==*--C.first.end( ) ) {
      C.first.pop_back( );
      _EM1.insert(C);
    }
  }

  map< vector< int > , int > _EM2;
  for( map< vector< int > , int >::const_iterator E_it=EM2.begin( ) ; E_it!=EM2.end( ) ; ++E_it ) {
    pair< vector< int > , int > C = *E_it;
    if( 1==*--C.first.end( ) ) {
      C.first.pop_back( );
      _EM2.insert(C);
    }
  }

  // 2. Sizes must be the same
  if( _EM1.size()!=_EM2.size() )
    return pair< bool , vector< int > >( false , vector< int >( ) );
  
  // 3. Find the appropriate vector
  pair< vector< int > , int > C1 = *_EM1.begin( );
  pair< vector< int > , int > C2 = *_EM2.begin( );
  if( C1.second!=C2.second )
    return pair< bool , vector< int > >( false , vector< int >( ) );
  vector< int > S( N , 0 );
  for( int i=0 ; i<N ; ++i )
    S[i] = C2.first[i]-C1.first[i];
  
  return pair< bool , vector< int > >( true , S );
}


//---------------------------------------------------------------------------//
//-------------------------------- testVector -------------------------------//
//---------------------------------------------------------------------------//


//! Test the given tail-vector.
/*!
  Test if \f$S(V^{-1} w1 V) = w2\f$, where \f$S\f$ is modEdgeMap( ... ).
*/
bool testVector( int N , const vector< int >& S , const Word& w1 , const Word& w2 , const vector< int >& V )
{
  Word C = getTailWord( N , V );
  map< vector< int > , int > EM1 = dropTrivialEdgesInMap( getEdgeMap( N , C * w1 * -C ) );
  map< vector< int > , int > EM2 = dropTrivialEdgesInMap( getEdgeMap( N , w2 ) );
  EM1 = dropTrivialEdgesInMap( modEdgeMap( N , S , EM1 ) );
  EM2 = dropTrivialEdgesInMap( modEdgeMap( N , S , EM2 ) );
  
  return EM1==EM2;
}


//---------------------------------------------------------------------------//
//----------------------------------- trivial -------------------------------//
//---------------------------------------------------------------------------//


bool FreeMetabelianGroupAlgorithms::trivial( int N , const Word& w )
{
  return dropTrivialEdgesInMap( getEdgeMap( N , w ) ).size()==0;
}


//---------------------------------------------------------------------------//
//------------------------------- conjugate ---------------------------------//
//---------------------------------------------------------------------------//


pair< bool , Word > FreeMetabelianGroupAlgorithms::conjugate( int N , Word w1 , Word w2 )
{
  // A. Compute the tails for w1 and w2
  vector< int > T1 = getTail( N , w1 );
  vector< int > T2 = getTail( N , w2 );
  
  // B. If tails are different then the answer is no
  if( T1!=T2 )
    return pair< bool , Word >( false , Word() );

  //----------------------------------------------------------
  // C. If the tails are trivial
  if( T1==vector< int >( N, 0 ) ) {
    // 1. Compute the edge-maps and check if they are of the same size
    map< vector< int > , int > EM1 = dropTrivialEdgesInMap( getEdgeMap( N , w1 ) );
    map< vector< int > , int > EM2 = dropTrivialEdgesInMap( getEdgeMap( N , w2 ) );
    if( EM1.size( )!=EM2.size( ) )
      return pair< bool , Word >( false , Word() );
    
    // 2. Check if the elements are trivial
    if( EM1.size( )==0 )
      return pair< bool , Word >( true , Word() );
    
    // 3. Determine the required shift 
    pair< vector< int > , int > C1 = *EM1.begin( );
    pair< vector< int > , int > C2 = *EM2.begin( );
    vector< int > S( N , 0 );
    for( int i=0 ; i<N ; ++i )
      S[i] = C2.first[i]-C1.first[i];
    EM1 = shiftEdgeMap( N , S , EM1 );
    if( EM1==EM2)
      return pair< bool , Word >( true , -getTailWord( N , S ) );
    else
      return pair< bool , Word >( false , Word() );
  }

  //----------------------------------------------------------
  // D. Move the space so that x1-component is positive and others are non-negative
  vector< Word > Images(N);
  for( int j=0 ; j<N ; ++j )
    Images[j] = Word( j+1 );
  if( T1[0]==0 ) {

    // 1. Find the first non-trivial component xi and swap(x1,xi)
    for( int i=1 ; i<N ; ++i ) {
      if( T1[i]!=0 ) {
	Images[0] = Word( i+1 );
	Images[i] = Word( 1 );
	swap( T1[0] , T1[i] );
	break;
      }
    }
  }
  Map M1( N , N , Images );
  w1 = M1.imageOf( w1 );
  w2 = M1.imageOf( w2 );

  //----------------------------------------------------------
  // E. Make negative components positive
  vector< Word > Images2( N );
  for( int i=0 ; i<N ; ++i ) {
    if( T1[i]>=0 )
      Images2[i] = Word(+i+1);
    else
      Images2[i] = Word(-i-1);
    T1[i] = abs(T1[i]);
  }
  Map M2( N , N , Images2 );
  w1 = M2.imageOf( w1 );
  w2 = M2.imageOf( w2 );
  
  //----------------------------------------------------------
  // F. Compute the strip-forms for elements
  map< vector< int > , int > EM1 = dropTrivialEdgesInMap( modEdgeMap( N , T1 , getEdgeMap( N , w1 ) ) );
  map< vector< int > , int > EM2 = dropTrivialEdgesInMap( modEdgeMap( N , T1 , getEdgeMap( N , w2 ) ) );
 
  //----------------------------------------------------------
  // G. For all possible x1-transformations:
  int u1 = T1[0];
  for( int i=0 ; i<u1 ; ++i ) {

    // 1. Apply an x1-shift
    EM1 = x1shiftEdgeMap( N , T1 , EM1 );

    // 2. Suggest a tail-conjugator
    pair< bool , vector< int > > r = suggestTailConjugator( N , EM1 , EM2 );
    r.second[0] += i+1;

    // 3. Test a tail conjugator
    if( r.first ) {
      if( testVector( N , T1 , w1 , w2 , r.second ) ) {
	// Find an actual conjugator
	Word C = getTailWord( N , r.second );
	Word D = getLoopConjugator( N , T1 , w1 , -C*w2*C );
	return pair< bool , Word >( true , M1.imageOf( M2.imageOf( D*-C ) ) );
      }
    }
  }
  
  return pair< bool , Word >( false , Word() );
}


//---------------------------------------------------------------------------//
//---------------------------- getWordFromEdgeMap ---------------------------//
//---------------------------------------------------------------------------//


void exhaust( int N , map< vector< int > , set< pair< int , int > > >& adjacencyList , vector< int > I , vector< int > T , list< int >& R , list< int >::iterator R_it )
{
  Word result;

  // pick an edge leaving the current vertex - I
  pair< int , int > P = *adjacencyList[I].begin( );
  int direction = P.first;
  int arity = P.second;
  
  // update the adjacency list
  adjacencyList[I].erase( P );
  if( arity>1 )
    adjacencyList[I].insert( pair< int , int >( direction , arity-1 ) );

  // save the value of I and compute a new value of I
  vector< int > I0 = I;
  I[abs(direction)-1] += direction<0 ? -1 : 1;
  
  // add a new letter to a word
  R_it = R.insert( R_it , direction );
  list< int >::iterator R_it2 = R_it;
  ++R_it2;
  
  // continue the process if we are not at the initial vertex
  if( I!=T )
    exhaust( N , adjacencyList , I , T , R , R_it2 );
  
  ++R_it;
  while( !adjacencyList[I].empty( ) )
    exhaust( N , adjacencyList , I , I , R , R_it );
  adjacencyList.erase( I );
}


Word FreeMetabelianGroupAlgorithms::getWordFromEdgeMap( int N , const map< vector< int > , int >& EM )
{
  Word result;
  
  // adjacencyList = point in Z^n -> (edge direction, arity)
  map< vector< int > , set< pair< int , int > > > adjacencyList;
  for( map< vector< int > , int >::const_iterator E_it=EM.begin( ) ; E_it!=EM.end( ) ; ++E_it ) {
    pair< vector< int > , int > C = *E_it;
    int direction = C.first[N];
    C.first.pop_back( );
    if( C.second<0 ) {
      ++C.first[direction-1];
      adjacencyList[C.first].insert( pair< int , int >( -direction , -C.second ) );
    } else {
      adjacencyList[C.first].insert( pair< int , int >(  direction ,  C.second ) );
    }
  }

  // Cover the components
  while( !adjacencyList.empty( ) ) {
    pair< vector< int > , set< pair< int , int > > > C = *adjacencyList.begin( );
    list< int > R;
    exhaust( N , adjacencyList , C.first , C.first , R , R.begin( ) );
    Word TW = getTailWord( N , C.first );
    result *= TW*Word( R )*-TW;
  }
  
  return result;
}
