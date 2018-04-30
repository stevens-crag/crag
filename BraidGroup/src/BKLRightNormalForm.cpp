
#include "BKLRightNormalForm.h"
#include "braid_group.h"
#include "Word.h"


//---------------------------------------------------------------------------//
//--------------------------- BKLRightNormalForm ----------------------------//
//---------------------------------------------------------------------------//


ostream& operator << ( ostream& os, const BKLRightNormalForm& rep )
{
  const list< Permutation >& decomposition = rep.getDecomposition( );
  list< Permutation >::const_iterator d_it = decomposition.begin( );
  
  for( int i=0 ; d_it!=decomposition.end( ) ; ++d_it )
    os << *d_it << endl;
  
  os << "Power = " << rep.getPower( ) << endl;
  return os;
}

//---------------------------------------------------------------------------//
//--------------------------- BKLRightNormalForm ----------------------------//
//---------------------------------------------------------------------------//


bool BKLRightNormalForm::operator !=( const BKLRightNormalForm& bkl ) const
{
  if( theOmegaPower!=bkl.theOmegaPower )
    return false;
  return theDecomposition!=bkl.theDecomposition;
}


bool BKLRightNormalForm::operator ==( const BKLRightNormalForm& bkl ) const
{
  return theOmegaPower==bkl.theOmegaPower && theDecomposition==bkl.theDecomposition;
}


bool BKLRightNormalForm::operator < ( const BKLRightNormalForm& bkl ) const
{
  if( theOmegaPower<bkl.theOmegaPower )
    return true;
  if( theOmegaPower>bkl.theOmegaPower )
    return false;
  return theDecomposition<bkl.theDecomposition;
}


//---------------------------------------------------------------------------//
//--------------------------- BKLRightNormalForm ----------------------------//
//---------------------------------------------------------------------------//


BKLRightNormalForm::BKLRightNormalForm( const BraidGroup& G , const Word& w ) :
  theRank( G.getRank( ) ),
  theOmegaPower( 0 )
{
#ifdef ERROR_EXISTS
    exit(1);
#endif
  const Permutation omega = getTinyTwistPermutation( theRank );
  cout << "omega = " << omega << endl;
  
  // 1. compute reverse permutation decomposition of a given braid word
  int shift = 0;
  for( Word::const_iterator w_it=w.begin( ) ; w_it!=w.end( ) ; w_it++ ) {
    int g = *w_it;
    if( g>0 ) {
      
      int letter  =  g-1;
      int strand1 = (letter+shift  )%theRank;
      int strand2 = (letter+shift+1)%theRank;
      
      Permutation transvection( theRank );
      transvection.change( strand1 , strand2 );
      
      theDecomposition.push_back( transvection );
    } else {
      int letter = -g-1;
      shift       = (shift+theRank-1) % theRank;
      int strand1 = (letter+shift   ) % theRank;
      int strand2 = (letter+shift+1 ) % theRank;
      
      Permutation transvection( theRank );
      transvection.change( strand1 , strand2 );
      
      theDecomposition.push_back( omega*transvection );
      theOmegaPower -= 1;
    }
  }
  
  
  // for( list< Permutation >::iterator l_it = theDecomposition.begin( ) ; l_it!=theDecomposition.end( ) ; ++l_it )
  cout << "------------------------------" << endl;
  cout << *this << endl;
  cout << "------------------------------" << endl;
  
  adjustDecomposition( theRank , theOmegaPower , theDecomposition );
}


//---------------------------------------------------------------------------//
//------------------------------- compute -----------------------------------//
//---------------------------------------------------------------------------//


//! Function adjusting the decomposition in a normal form.
void BKLRightNormalForm::adjustDecomposition( int rank , int& power , list<Permutation>& decomp )
{
  const Permutation omega = getTinyTwistPermutation( rank );
  
  int shift = 0;
  list< Permutation>::iterator it1 =  decomp.begin( );
  for( ; it1!=decomp.end( ) ; it1++ ) {
    
    *it1 = (*it1).tinyFlip( shift );
    
    list< Permutation >::iterator it2 = it1;
    for( ; it2!=decomp.begin( ) ; it2-- ) {
      list< Permutation >::iterator it3 = it2;
      it3--;
      
      transformationResult tr_res = transform( rank , *it3 , *it2 );
      switch( tr_res ) {
      case ONE_MULTIPLIER:
	cout << "ONE_MULTIPLIER" << endl;
	if( it1==it2 )
	  it1 = it2 = decomp.erase( it3 );
	else
	  it2 = decomp.erase( it3 );
	it2++;
	break;
      case NO_CHANGE:
	cout << "NO_CHANGE" << endl;
	it2 = decomp.begin( );
	it2++;
	break;
      }
    }

    if( *it1==omega ) {
      power++;
      it1 = decomp.erase( it1 );
      it1--;
      shift = (shift+1)%rank;
    }
  }
}


//---------------------------------------------------------------------------//
//------------------------------- compute -----------------------------------//
//---------------------------------------------------------------------------//

/*
BKLRightNormalForm* 
BKLRightNormalForm::compute( const BraidGroup& G , const Word& w )
{
  int theIndex = G.getIndex( );
  
  list< Permutation > result;
  
  const Permutation omega = getTinyTwistPermutation( theIndex );
  int omega_power = 0;

  // 1. compute reverse permutation decomposition of a given braid word
  int shift = 0;
  for( int i=0 ; i<w.length( ) ; i++ ) {
    if( ord(w[i])>0 ) {
      
      int letter = ord(w[i])-1;
      int strand1 = (letter+shift  )%theIndex;
      int strand2 = (letter+shift+1)%theIndex;
      
      Permutation transvection( theIndex );
      transvection.change( strand1 , strand2 );
      
      result.push_back( transvection );
    } else {
      int letter = -ord(w[i])-1;
      shift = (shift+theIndex-1)%theIndex;
      int strand1 = (letter+shift  )%theIndex;
      int strand2 = (letter+shift+1)%theIndex;
      
      Permutation transvection( theIndex );
      transvection.change( strand1 , strand2 );
      
      result.push_back( omega*transvection );
      omega_power -= 1;
    }
  }
  

  // 2. compute normal form of a decomposition
  shift = 0;
  list< Permutation>::iterator it1 =  result.begin( );
  for( ; it1!=result.end( ) ; it1++ ) {
    
    *it1 = (*it1).tinyFlip( shift );
    
    list< Permutation >::iterator it2 = it1;
    for( ; it2!=result.begin( ) ; it2-- ) {
      list< Permutation >::iterator it3 = it2;
      it3--;
      
      transformationResult tr_res = transform( theIndex , *it3 , *it2 );
      switch( tr_res ) {
      case ONE_MULTIPLIER:
	if( it1==it2 )
	  it1 = it2 = result.erase( it3 );
	else
	  it2 = result.erase( it3 );
	it2++;
	break;
      case NO_CHANGE:
	it2 = result.begin( );
	it2++;
	break;
      }
    }

    if( *it1==omega ) {
      omega_power++;
      it1 = result.erase( it1 );
      it1--;
      shift = (shift+1)%theIndex;
    }

  }
  
  vector< Permutation > result1( result.size( ) );
  it1 =  result.begin( );
          for( int i=0 ; it1!=result.end( ) ; ++it1 , i++ )
    result1[i] = *it1;
  
  return new BKLRightNormalForm( omega_power , result1 );
}
*/


//---------------------------------------------------------------------------//
//------------------------------- transform ---------------------------------//
//---------------------------------------------------------------------------//


BKLRightNormalForm::transformationResult 
BKLRightNormalForm::transform( int theIndex , Permutation& p1 , Permutation& p2 )
{
  Permutation p4 = p1*p2;

  transformationResult res = TWO_MULTIPLIERS;
  const Permutation omega = getTinyTwistPermutation( theIndex );
  
  cout << p1 << "  ---  " << p2 << endl;
  
  Permutation p3 = p1.meet2( p2.inverse( )*omega );
  cout << "   p3 = " << p3 << endl;
  if( p3==p1 )
    res = ONE_MULTIPLIER;
  if( p3.isTrivial( ) )
    res = NO_CHANGE;
  p2 = p3*p2;
  p1 *= p3.inverse( );
  
  Permutation p5 = p1*p2;
  cout << "    p1, p2: " << p1 << " , " << p2 << endl;
  cout << "    p4, p5: " << p4 << " , " << p5 << endl;
  
  return res;
}


//---------------------------------------------------------------------------//
//--------------------------------- inverse ---------------------------------//
//---------------------------------------------------------------------------//


BKLRightNormalForm BKLRightNormalForm::inverse( ) const
{
  int power = -theOmegaPower-theDecomposition.size( );
  const Permutation omega = getTinyTwistPermutation( theRank );
  list<Permutation> result;
  
  int flip = power;
  for( list<Permutation>::const_iterator d_it=theDecomposition.end( ) ; d_it!=theDecomposition.end( ) ; ) {
    --d_it;
    Permutation p = (*d_it).inverse( ) * omega;
    result.push_back( p.tinyFlip( -flip++ ) );
  }
  
  return BKLRightNormalForm( theRank , power , result );
}


//---------------------------------------------------------------------------//
//--------------------------------- multiply --------------------------------//
//---------------------------------------------------------------------------//


BKLRightNormalForm BKLRightNormalForm::multiply( const BKLRightNormalForm& rep ) const
{
  const Permutation omega = getTinyTwistPermutation( theRank );
  
  // 1. shift omegas to the left
  int power = theOmegaPower + rep.theOmegaPower;

  // 2. create initial decomposition
  list< Permutation > blocks;
  for( list< Permutation >::const_iterator d_it=theDecomposition.begin( ) ; d_it!=theDecomposition.end( ) ; ++d_it ) {
    blocks.push_back( (*d_it).tinyFlip( -rep.theOmegaPower ) );
    // cout << theDecomposition[t] << " -> " 
    // << theDecomposition[t].tinyFlip( rep.theOmegaPower ) << endl;
  }
  for( list< Permutation >::const_iterator d_it=rep.theDecomposition.begin( ) ; d_it!=rep.theDecomposition.end( ) ; ++d_it )
    blocks.push_back( *d_it );

  /*
  list< Permutation >::iterator it = blocks.begin( );
  cout << "+++++++++++++++++++++++++++++++" << endl;
  for( ; it!=blocks.end( ) ; ++it )
    cout << *it << endl;
  cout << "+++++++++++++++++++++++++++++++" << endl;
  */

  adjustDecomposition( theRank , power , blocks );
  return BKLRightNormalForm( theRank , power , blocks );
}


//---------------------------------------------------------------------------//
//--------------------------------- getWord ---------------------------------//
//---------------------------------------------------------------------------//


Word BKLRightNormalForm::getWord( ) const
{
  Word result;

  vector< int > geodesic;
  for( list< Permutation >::const_iterator d_it = theDecomposition.begin( ) ; d_it!=theDecomposition.end( ) ; ++d_it ) {
    geodesic = (*d_it).geodesic( );
    for( int j=0 ; j<geodesic.size( ) ; ++j ) 
      result *= Generator( geodesic[j]+1 );
  }
  
  if( theOmegaPower ) {
    const Permutation omega = getTinyTwistPermutation( theRank );
    vector< int > geodesic = omega.geodesic( );
    Word omegaWord;
    for( int i=0 ; i<geodesic.size( ) ; ++i ) omegaWord *= Generator( geodesic[i]+1 );
    
    if( theOmegaPower<0 )
      omegaWord = omegaWord.inverse( );
    cout << "         /// " << omegaWord << endl;
    for( int i=0 ; i<abs(theOmegaPower) ; ++i )
      result *= omegaWord;
  }
  
  return result;
}
