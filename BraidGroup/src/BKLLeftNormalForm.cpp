
#include "BKLLeftNormalForm.h"
#include "braid_group.h"
#include "Word.h"


ostream& operator << ( ostream& os, const BKLLeftNormalForm& bkl )
{
  const list< Permutation >& decomposition = bkl.getDecomposition( );
  
  list< Permutation >::const_iterator it = decomposition.begin( );
  for( ; it!=decomposition.end( ) ; ++it )
    os << (*it) << endl;

  os << "Power = " << bkl.getPower( ) << endl;
  return os;
}


bool BKLLeftNormalForm::operator !=( const BKLLeftNormalForm& bkl ) const
{
  if( theOmegaPower!=bkl.theOmegaPower )
    return false;
  return theDecomposition!=bkl.theDecomposition;
}


bool BKLLeftNormalForm::operator ==( const BKLLeftNormalForm& bkl ) const
{
  return theOmegaPower==bkl.theOmegaPower && theDecomposition==bkl.theDecomposition;
}


bool BKLLeftNormalForm::operator < ( const BKLLeftNormalForm& bkl ) const
{
  if( theOmegaPower<bkl.theOmegaPower )
    return true;
  if( theOmegaPower>bkl.theOmegaPower )
    return false;
  return theDecomposition<bkl.theDecomposition;
}


//---------------------------------------------------------------------------//
//------------------------------- compute -----------------------------------//
//---------------------------------------------------------------------------//


BKLLeftNormalForm::BKLLeftNormalForm( const BraidGroup& G , const Word& w ) :
  theRank( G.getRank( ) ),
  theOmegaPower( 0 )
{
  const Permutation omega = getTinyTwistPermutation( theRank );

  // 1. compute reverse permutation decomposition of a given braid word
  int shift = 0;
  for( Word::const_iterator w_it=w.end( ) ; w_it!=w.begin( ) ; ) {
    int g = *(--w_it);
    if( g>0 ) {
      
      int letter = g-1;
      int strand1 = (letter+shift  )%theRank;
      int strand2 = (letter+shift+1)%theRank;
      
      Permutation transvection( theRank );
      transvection.change( strand1 , strand2 );
      
      theDecomposition.push_front( transvection );
    } else {
      int letter = -g-1;
      shift = (shift+1)%theRank;
      int strand1 = (letter+shift  )%theRank;
      int strand2 = (letter+shift+1)%theRank;
      
      Permutation transvection( theRank );
      transvection.change( strand1 , strand2 );
      
      theDecomposition.push_front( transvection*omega );
      theOmegaPower -= 1;
    }
  }
  
  adjustDecomposition( theRank , theOmegaPower , theDecomposition );
}


//---------------------------------------------------------------------------//
//------------------------------- transform ---------------------------------//
//---------------------------------------------------------------------------//


BKLLeftNormalForm::transformationResult 
BKLLeftNormalForm::transform( int rank , Permutation& p1 , Permutation& p2 )
{
  transformationResult res = TWO_MULTIPLIERS;
  const Permutation omega = getTinyTwistPermutation( rank );
  
  Permutation p3 = p2.meet2( p1.inverse( )*omega );
  if( p3==p2 )
    res = ONE_MULTIPLIER;
  if( p3.isTrivial( ) )
    res = NO_CHANGE;
  p2 = p3.inverse( )*p2;
  p1 *= p3;

  return res;
}


//---------------------------------------------------------------------------//
//--------------------------------- inverse ---------------------------------//
//---------------------------------------------------------------------------//



BKLLeftNormalForm BKLLeftNormalForm::operator - ( ) const
{
  return inverse( );
}

//---------------------------------------------------------------------------//
//--------------------------------- inverse ---------------------------------//
//---------------------------------------------------------------------------//


BKLLeftNormalForm::NF BKLLeftNormalForm::inverse( ) const
{
  int dec_size = theDecomposition.size( );
  
  int power = -theOmegaPower-dec_size;
  const Permutation omega = getTinyTwistPermutation( theRank );
  list<Permutation> result;
  
  int flip = theOmegaPower+1;
  list<Permutation>::const_iterator it = theDecomposition.begin( );
  for( ; it!=theDecomposition.end( ) ; ++it ) {
    Permutation p = (*it).inverse( ) * omega;
    result.push_front( p.tinyFlip( flip ) );
    // cout << -flip << "," << (*it) << " : " << p << " => " << p.tinyFlip( -flip ) << endl;
    flip++;
  }
  
  return BKLLeftNormalForm::NF( theRank , power , result );
}


//---------------------------------------------------------------------------//
//--------------------------------- operator* -------------------------------//
//---------------------------------------------------------------------------//


BKLLeftNormalForm& BKLLeftNormalForm::operator *= ( const BKLLeftNormalForm& bkl )
{
  NF res = multiply( bkl );
  theOmegaPower = res.second;
  theDecomposition = res.third;
  return *this;
}


BKLLeftNormalForm  BKLLeftNormalForm::operator *  ( const BKLLeftNormalForm& bkl ) const
{
  return multiply( bkl );
}


//---------------------------------------------------------------------------//
//--------------------------------- multiply --------------------------------//
//---------------------------------------------------------------------------//


BKLLeftNormalForm::NF
BKLLeftNormalForm::multiply( const BKLLeftNormalForm& bkl ) const
{
  int theRank = bkl.theRank;
  
  const Permutation omega = getTinyTwistPermutation( theRank );
  
  // 1. shift omegas to the left
  int power = theOmegaPower + bkl.theOmegaPower;

  // 2. create initial decomposition
  list< Permutation > blocks;
  list< Permutation >::const_iterator it = theDecomposition.begin( );
  for( ; it!=theDecomposition.end( ) ; ++it )
    blocks.push_back( (*it).tinyFlip( -bkl.theOmegaPower ) );
  it = bkl.theDecomposition.begin( );
  for( ; it!=bkl.theDecomposition.end( ) ; ++it )
    blocks.push_back( *it );
  
  adjustDecomposition( theRank , power , blocks );
  
  return BKLLeftNormalForm::NF( theRank , power , blocks );
}


//---------------------------------------------------------------------------//
//--------------------------- adjustDecomposition ---------------------------//
//---------------------------------------------------------------------------//


void 
BKLLeftNormalForm::adjustDecomposition
( int rank , int& power , list<Permutation>& decomp )
{
  const Permutation omega = getTinyTwistPermutation( rank );
  
  int shift = 0;
  list< Permutation >::iterator it1 = decomp.end( );
  while( it1!=decomp.begin( ) ) {
    --it1;

    *it1 = (*it1).tinyFlip( shift );
    list< Permutation >::iterator it2 = it1;
    for( ; it2!=decomp.end( ) ; ++it2 ) {
      
      list< Permutation >::iterator it3 = it2;
      if( ++it3==decomp.end( ) )
	continue;
      
      transformationResult tr_res = transform( rank , *it2 , *it3 );
      switch( tr_res ) {
      case ONE_MULTIPLIER:
	if( it1==it2) {
	  it2 = decomp.erase( it3 );
	  it2--;
	  it1 = it2;
	} else {
	  it2 = decomp.erase( it3 );
	  it2--;
	}
	it2--;
	break;
      case NO_CHANGE:
	// means we can stop the inner loop
	it2 = decomp.end( );
	it2--;
	break;
      }
    }
    
    if( *it1==omega ) {
      power++;
      it1 = decomp.erase( it1 );
      shift = (shift+rank-1)%rank;
    }
  }
}


//---------------------------------------------------------------------------//
//--------------------------------- getWord ---------------------------------//
//---------------------------------------------------------------------------//


Word BKLLeftNormalForm::getWord( ) const
{
  Word result;

  if( theOmegaPower ) {
    Word omegaWord;
    const Permutation omega = getTinyTwistPermutation( theRank );
    vector< int > geodesic = omega.geodesic( );
    for( int i=0 ; i<geodesic.size( ) ; ++i )
      omegaWord *= Generator( geodesic[i]+1 );
    
    if( theOmegaPower<0 )
      omegaWord = omegaWord.inverse( );
    for( int i=0 ; i<abs(theOmegaPower) ; ++i )
      result *= omegaWord;
  }
  
  vector< int > geodesic;
  for( list< Permutation >::const_iterator d_it = theDecomposition.begin( ) ; d_it!=theDecomposition.end( ) ; ++d_it ) {
    geodesic = (*d_it).getWordPresentation( );
    for( int j=0 ; j<geodesic.size( ) ; ++j )
      result *= Generator( geodesic[j] );
  }
  
  return result;
}


//---------------------------------------------------------------------------//
//------------------------------ cycle/decycle ------------------------------//
//---------------------------------------------------------------------------//


pair< BKLLeftNormalForm , BKLLeftNormalForm >
BKLLeftNormalForm::cycle( ) const
{
  typedef pair< BKLLeftNormalForm , BKLLeftNormalForm > PNF;
  
  if( theDecomposition.empty( ) )
    return PNF( *this , BKLLeftNormalForm( ) );
  
  BKLLeftNormalForm result = *this;
  
  // create initial decomposition of a result
  Permutation conj = (*result.theDecomposition.begin( )).tinyFlip( result.theOmegaPower );
  result.theDecomposition.push_back( conj );
  result.theDecomposition.pop_front( );
  
  // transform to a normal form
  adjustDecomposition( result.theRank , result.theOmegaPower , result.theDecomposition );
  return PNF( result , BKLLeftNormalForm( theRank , 0 , list<Permutation>( 1 , conj ) ) );
}


//---------------------------------------------------------------------------//
//------------------------------ cycle/decycle ------------------------------//
//---------------------------------------------------------------------------//


pair< BKLLeftNormalForm , BKLLeftNormalForm >
BKLLeftNormalForm::decycle( ) const
{
  typedef pair< BKLLeftNormalForm , BKLLeftNormalForm > PNF;

  if( theDecomposition.empty( ) )
    return PNF( *this , BKLLeftNormalForm( ) );

  BKLLeftNormalForm result = *this;

  // create initial decomposition of a result
  Permutation conj = *--result.theDecomposition.end( );
  result.theDecomposition.push_front( conj.tinyFlip( -theOmegaPower ) );
  result.theDecomposition.pop_back( );
  
  BKLLeftNormalForm conjugator( theRank , 0 , list<Permutation>( 1 , conj ) );
  conjugator = -conjugator;
  
  // transform to a normal form
  adjustDecomposition( theRank , result.theOmegaPower , result.theDecomposition );
  return PNF( result , conjugator );
}


//---------------------------------------------------------------------------//
//--------------------------- findSSSRepresentative -------------------------//
//---------------------------------------------------------------------------//


pair< BKLLeftNormalForm , BKLLeftNormalForm >
BKLLeftNormalForm::findSSSRepresentative( ) const
{
  typedef pair< BKLLeftNormalForm , BKLLeftNormalForm > PNF;
  BKLLeftNormalForm cur_conj;
  BKLLeftNormalForm conjugator;
  
  BKLLeftNormalForm result = *this;
  BKLLeftNormalForm cur_nf = *this;
  for( bool progress=true ; progress ; ) {
    progress = false;
    for( int i=0; i<theRank ; ++i ) {
      PNF pr = cur_nf.cycle( );
      cur_nf = pr.first;
      cur_conj = cur_conj.multiply( pr.second );
      if( cur_nf.theOmegaPower>result.theOmegaPower ) {
	result = cur_nf;
        conjugator = cur_conj;
	progress = true;
	break;
      }
    }
  }
  
  cur_nf = result;
  cur_conj = conjugator;
  for( bool progress=true ; progress ; ) {
    progress = false;
    for( int i=0; i<theRank ; ++i ) {
      PNF pr = cur_nf.decycle( );
      cur_nf = pr.first;
      cur_conj = cur_conj.multiply( pr.second );
      if( cur_nf.theDecomposition.size( )<result.theDecomposition.size( ) ) {
	result = cur_nf;
        conjugator = cur_conj;
	progress = true;
	break;
      }
    }
  }

  return pair< BKLLeftNormalForm::NF , BKLLeftNormalForm >( result , conjugator );
}


//---------------------------------------------------------------------------//
//--------------------------- getSimpleConjugators --------------------------//
//---------------------------------------------------------------------------//


set<Permutation> BKLLeftNormalForm::getSimpleConjugators( const NF& bkl )
{
  int theRank = bkl.first;
  const list< Permutation >& decomp = bkl.third;
  
  set<Permutation> result;
  
  for( int i=0 ; i<theRank-1 ; ++i ) {
    Permutation c = Permutation( theRank );
    c.change( i , i+1 );

    // cout << "i = " << i << endl;
    
    bool done = false;
    if( bkl.second>0 )
      done = true;
    
    Permutation v = c;
    list< Permutation >::const_iterator it = decomp.begin( );
    while( !done ) {
      for( ; it!=decomp.end( ) ; ++it ) {
	v = (*it).inverse() * v.join2( *it );
	if( v.isTrivial( ) ) {
	  done = true;
	  break;
	}
      }
      if( !done ) {
	v = v.join2( c );
	if( v==c )
	  done = true;
	else
	  c = v;
      }
      // cout << c << endl;
    }
    
    result.insert( c );
  }
  
  return result;
}


//---------------------------------------------------------------------------//
//------------------------- getSimpleSummitConjugators ----------------------//
//---------------------------------------------------------------------------//


// Is it finished?
set<Permutation> BKLLeftNormalForm::getSimpleSummitConjugators( const NF& bkl )
{
  set<Permutation> result;
  
  set<Permutation> conjs = getSimpleConjugators( bkl );
  
  
  return conjs;
}

  
//---------------------------------------------------------------------------//
//------------------------------- areConjugate ------------------------------//
//---------------------------------------------------------------------------//


pair<bool,Word>
BKLLeftNormalForm::areConjugate( const BKLLeftNormalForm& bkl ) const
{
  pair<bool,Word> result;

  // typedef pair< int , vector< Permutation > > NF;

  // 1. find representative of SSS1
  pair< NF , BKLLeftNormalForm > pr1 = findSSSRepresentative( );
  NF nf1 = pr1.first;
  pair< NF , BKLLeftNormalForm > pr2 = bkl.findSSSRepresentative( );
  NF nf2 = pr2.first;
  
  if( nf1.second!=nf2.second || nf1.third.size( )!=nf2.third.size( ) )
    return pair<bool,Word>( false , Word( ) );
  
  if( nf1==nf2 ) {
    cout << "Equal" << endl;
    return pair<bool,Word>( true , Word( ) );
  }
  
  // 2. multiply by the n-th power of omega
  if( nf1.first<0 ) {
    int power = (1-nf1.first/theRank)*theRank;
    nf1.first += power;
    nf2.first += power;
  }
  
  // 3. construct super summit sets
  set< NF > sss_new1;
  set< NF > sss_checked1;
  set< NF > sss_new2;
  set< NF > sss_checked2;
  
  sss_new1.insert( nf1 );
  sss_new2.insert( nf2 );
  
  for( int i=0 ; i<1 && sss_new1.size( ) && sss_new2.size( ) ; ++i ) {
    
    // 4. select unchecked element of SSS1
    NF el1 = *sss_new1.begin( );
    sss_checked1.insert( el1 );
    sss_new1.erase( sss_new1.begin( ) );
    
    BKLLeftNormalForm tmp_nf = el1;
    cout << "====================================" << endl;
    cout << tmp_nf << endl;
    cout << "====================================" << endl;

    // 5. compute the set of simple summit conjugators 
    set<Permutation> conj1 = getSimpleSummitConjugators( el1 );
    
    // 6. conjugate selected element by conjugators
    set<Permutation>:: iterator it = conj1.begin( );
    for( ; it!=conj1.end( ) ; ++it ) {
      BKLLeftNormalForm r_mult( theRank , 0 , list<Permutation>( 1 , *it ) );
      BKLLeftNormalForm l_mult = r_mult.inverse( );
      BKLLeftNormalForm new_el1 = l_mult.multiply( el1 );
      new_el1 = new_el1.multiply( r_mult );
      
      if( sss_new1.find( new_el1 )!=sss_new1.end( ) )
	continue;
      if( sss_checked1.find( new_el1 )!=sss_checked1.end( ) )
	continue;
      
      if( sss_new2.find( new_el1 )!=sss_new2.end( ) )
	return pair<bool,Word>( true , Word( ) );
      if( sss_checked2.find( new_el1 )!=sss_checked2.end( ) )
	return pair<bool,Word>( true , Word( ) );

      cout << "====================================" << endl;
      cout << "N = " << i << endl;
      cout << new_el1 << endl;
      cout << "====================================" << endl;
      
      sss_new1.insert( new_el1 );
    }
  }
  
  cout << "Similar" << endl;
  return pair<bool,Word>( true , Word( ) );
}
