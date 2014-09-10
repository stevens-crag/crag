// Copyright (C) 2005 Alexander Ushakov
// Contents: Implementation of class ThRightNormalForm
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//

#include <cassert>
#include <fstream>

#include "BraidGroup.h"
#include "ShortBraidForm.h"
#include "ThRightNormalForm.h"
#include "ThLeftNormalForm.h"
#include "Word.h"


ostream& printOn( ostream& os , const list<Permutation>& lp  )
{
  os << "(";
  list<Permutation>::const_iterator lp_it = lp.begin( );
  for( ; lp_it!=lp.end( ) ; ++lp_it ) {
    // os << (*lp_it) << " -> " << (*lp_it).geodesic( ).size( ) << endl;
    if( lp_it!=lp.begin( ) )
      os << ",";
    os.width(2);
    os << (*lp_it).geodesic( ).size( );
  }
  os << ")";
  
  return os;
}



ostream& operator << ( ostream& os, const ThRightNormalForm& rep )
{
  const list< Permutation >& decomposition = rep.getDecomposition( );
  
  list< Permutation >::const_iterator it = decomposition.begin( );
  for( ; it!=decomposition.end( ) ; ++it )
    os << *it << endl;
  
  os << "Power = " << rep.getPower( ) << endl;
  os << "Rank = "  << rep.getRank ( ) << endl;
  return os;
}


//---------------------------------------------------------------------------//
//-------------------------- ThRightNormalForm ------------------------------//
//---------------------------------------------------------------------------//


ThRightNormalForm::operator ThLeftNormalForm( ) const
{
  if( theOmegaPower%2==0 ) {
    ThLeftNormalForm nf( theRank , theOmegaPower , theDecomposition );
    nf.adjust( );
    return nf;
  }
  
  list< Permutation > D;
  for( list< Permutation >::const_iterator d_it=theDecomposition.begin( ) ; d_it!=theDecomposition.end( ) ; ++d_it )
    D.push_back( (*d_it).flip( ) );
  
  ThLeftNormalForm nf( theRank , theOmegaPower , D );
  nf.adjust( );
  return nf;
}


//---------------------------------------------------------------------------//
//--------------------------- ThRightNormalForm -----------------------------//
//---------------------------------------------------------------------------//


ThRightNormalForm::ThRightNormalForm( const BraidGroup& G , const Word& w ) :
	theRank( G.getRank( ) ),
	theOmegaPower( 0 )
{
  const Permutation omega = Permutation::getHalfTwistPermutation( theRank );

  // 1. compute permutation decomposition of a given braid word
  Permutation curMult( theRank );
  //bool trivialMult = true;
  for( Word::const_iterator w_it=w.end( ) ; w_it!=w.begin( ) ; ) {

    // find the current block
    Word::const_iterator w_it2 = --w_it;
    for( ; (*w_it2<0)==(*w_it<0) && w_it2!=w.begin( ) ; --w_it2 );
    if( (*w_it2<0)!=(*w_it<0) ) ++w_it2;
    w_it++;
    
    // cout << " ***** " << *w_it2 << endl;
    // check the sign of the block
    if( *w_it2>0 ) {
      
      // process the positive block
      for( Word::const_iterator w_it3=w_it ; w_it3!=w_it2 ; ) {
	
	--w_it3;
	int gen = *w_it3;
	// cout << " +++++ " << gen << endl;
	if( curMult[gen-1]>curMult[gen] ) {
	  theDecomposition.push_front( curMult );
	  curMult = Permutation( theRank );
	}
	swap( curMult[gen-1] , curMult[gen] );
	//trivialMult = false;
      }
      theDecomposition.push_front( curMult );
      curMult = Permutation( theRank );
      //trivialMult = true;
      
    } else {
      
      // process the negative block
      list< Permutation > mult;
      for( Word::const_iterator w_it3 =w_it2 ; w_it3!=w_it ; ++w_it3 ) {
	
	int gen = -*w_it3;
	// cout << " ------ " << gen << endl;
	if( curMult[gen-1]>curMult[gen] ) {
	  mult.push_front( curMult );
	  curMult = Permutation( theRank );
	}
	swap( curMult[gen-1] , curMult[gen] );
	//trivialMult = false;
      }
      mult.push_front( curMult );
      curMult = Permutation( theRank );
      //trivialMult = true;
      
      list< Permutation >::iterator it = mult.begin( );
      for( ; it!=mult.end( ) ; ++it ) {
	theDecomposition.push_front( omega );
	if( !((*it).inverse( )*omega).isTrivial( ) )
	  theDecomposition.push_front( (*it).inverse( )*omega );
	theOmegaPower -= 2;
      }
    }
    
    w_it = w_it2;
  }
	
  adjustDecomposition( theRank , theOmegaPower , theDecomposition );
}


//---------------------------------------------------------------------------//
//--------------------------- adjustDecomposition ---------------------------//
//---------------------------------------------------------------------------//


void 
ThRightNormalForm::adjustDecomposition
( int rank , int& power , list<Permutation>& decomp )
{
  const Permutation omega = Permutation::getHalfTwistPermutation( rank );

  bool flip = false;
  list< Permutation>::iterator it1 =  decomp.begin( );
  for( ; it1!=decomp.end( ) ; it1++ ) {
    
    if( flip )
      *it1 = omega*(*it1)*omega;
		
    list< Permutation >::iterator it2 = it1;
    for( ; it2!=decomp.begin( ) ; it2-- ) {
      list< Permutation >::iterator it3 = it2;
      it3--;
      
      transformationResult tr_res = transform( rank , *it3 , *it2 );
      switch( tr_res ) {
      case TWO_MULTIPLIERS:
          assert(false);
      case ONE_MULTIPLIER:
				if( it1==it2 )
					it1 = it2 = decomp.erase( it3 );
				else
					it2 = decomp.erase( it3 );
				it2++;
				break;
      case NO_CHANGE:
				it2 = decomp.begin( );
				it2++;
				break;
      }
    }
		
    if( *it1==omega ) {
      power++;
      it1 = decomp.erase( it1 );
      it1--;
      flip = !flip;
    }
		
  }
}


//---------------------------------------------------------------------------//
//------------------------------- transform ---------------------------------//
//---------------------------------------------------------------------------//


// finds longest head of p2 
// that can be multiplied by p1 on the left
ThRightNormalForm::transformationResult
ThRightNormalForm::transform( int theRank , Permutation& p1 , Permutation& p2 )
{
  transformationResult res = TWO_MULTIPLIERS;
  const Permutation omega = Permutation::getHalfTwistPermutation( theRank );

  Permutation p3 = p1.RightGCD( omega * -p2 );
  if( p3==p1 )
    res = ONE_MULTIPLIER;
  if( p3.isTrivial( ) )
    res = NO_CHANGE;
  p2 = p3 * p2;
  p1 *= -p3;

  return res;
}


//---------------------------------------------------------------------------//
//------------------------------- inverse -----------------------------------//
//---------------------------------------------------------------------------//


ThRightNormalForm::NF
ThRightNormalForm::inverse( ) const
{
  int dec_size = theDecomposition.size( );
  
  int power = -theOmegaPower-dec_size;
  const Permutation omega = Permutation::getHalfTwistPermutation( theRank );
  list<Permutation> result;
  
  list<Permutation>::const_iterator it = --theDecomposition.end( );
  for( int i=0 ; i<dec_size ; ++i , --it )
    if( ( i%2==0 )==( theOmegaPower%2!=0 ) )
      result.push_back( omega * -(*it) );
    else
      result.push_back( -(*it) * omega );
  
  return NF( theRank , power , result );
}


//---------------------------------------------------------------------------//
//------------------------------- multiply ----------------------------------//
//---------------------------------------------------------------------------//


ThRightNormalForm::NF
ThRightNormalForm::multiply( const ThRightNormalForm& rep ) const
{
  const Permutation omega = Permutation::getHalfTwistPermutation( theRank );
  
  // 1. shift omegas to the right
  int power = theOmegaPower + rep.theOmegaPower;
  
  list< Permutation > blocks( theDecomposition );
  
  // 2. if power of the omega on the left is odd we must flip 
  //    permutations of the second form
  list< Permutation >::const_iterator it = rep.theDecomposition.begin( );
  if( theOmegaPower%2 ) {
    for( size_t t=0 ; t<rep.theDecomposition.size( ) ; ++t , ++it )
      blocks.push_back( omega * (*it) * omega );
  } else
    blocks.insert( blocks.end( ) , 
		   rep.theDecomposition.begin( ) , 
		   rep.theDecomposition.end( ) );
  
  adjustDecomposition( theRank , power , blocks );
  return NF( theRank , power , blocks );
}


//---------------------------------------------------------------------------//
//-------------------------------- getWord ----------------------------------//
//---------------------------------------------------------------------------//


Word ThRightNormalForm::getWord( ) const
{
  Word result;
  vector< int > decomposition;
  
  vector< int > geodesic;
  list<Permutation>::const_iterator it = theDecomposition.begin( );
  for( ; it!=theDecomposition.end( ) ; ++it ) {
    geodesic = (*it).geodesic( );
    for( size_t j=0 ; j<geodesic.size( ) ; ++j )
      result.push_back( geodesic[j]+1 );
  }

  if( theOmegaPower==0 )
    return result;
  
  Word omegaWord;
  const Permutation omega = Permutation::getHalfTwistPermutation( theRank );
  geodesic = omega.geodesic( );
  for( size_t i=0 ; i<geodesic.size( ) ; ++i )
    omegaWord.push_back( geodesic[i]+1 );
  
  if( theOmegaPower<0 )
    omegaWord = -omegaWord;
    
  for( int i=0 ; i<abs(theOmegaPower) ; ++i )
    result *= omegaWord;
  
  return result;
}


Word ThRightNormalForm::getShortWord( ) const
{
  Word result;
  Permutation omega = Permutation::getHalfTwistPermutation( getRank( ) );
  int power = getPower( );

  if( power<0 ) {

    const list< Permutation >& decomp = getDecomposition( );
    list<Permutation>:: const_iterator it = decomp.begin( );
    for( int j=0 ; it!=decomp.end( ) ; ++it, ++j ) {
      int n = j - decomp.size( ) - power;
      if( n<0 ) {
	vector< int > gd = (*it).geodesic();
	for( size_t t=0 ; t<gd.size() ; ++t )
	  result.push_back( gd[t]+1 );
      } else {
	
	Permutation p = ( n%2 == 1 ? (*it).inverse() * omega : omega * (*it).inverse() );
	
	vector<int> gd = p.geodesic();
	for( int t=gd.size( )-1 ; t>=0 ; --t )
	  result.push_back( -gd[t]-1 );
      }
    }
    Word omega_w = omega.geodesicWord( );
    omega_w = -omega_w;
    for( int j=decomp.size( ) ; j<-power ; ++j )
      result = omega_w*result;
  } else
    result = getWord( );
  
  return result;
}


//---------------------------------------------------------------------------//
//---------------------------- computeCentralizer ---------------------------//
//---------------------------------------------------------------------------//


set<ThRightNormalForm> 
ThRightNormalForm::computeCentralizer( ) const
{
  const Permutation omega = Permutation::getHalfTwistPermutation( theRank );
  
  set<ThRightNormalForm> result;
  
  map< ThRightNormalForm , ThRightNormalForm > new_states;
  map< ThRightNormalForm , ThRightNormalForm > unchecked_states;
  map< ThRightNormalForm , ThRightNormalForm > checked_states;
  
  if( theOmegaPower<0 ) {
    ThRightNormalForm rep = *this;
    rep.theOmegaPower = -theOmegaPower%2;
    cout << "Power = " << rep.theOmegaPower << endl;
    new_states[rep] = ThRightNormalForm( theRank );
  } else 
    new_states[*this] = ThRightNormalForm( theRank );
/*
  cout << endl << "  init:" << endl;
  cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
  cout << (*new_states.begin( )).first;
  cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
*/  
  
  while( !new_states.empty( ) ) {
    
    { // was unsupported by MSVC++
      // checked_states.insert( unchecked_states.begin(), unchecked_states.end() );
      map< ThRightNormalForm , ThRightNormalForm >::iterator m_it = unchecked_states.begin( );
      for( ; m_it!=unchecked_states.end( ) ; ++m_it )
        checked_states.insert( *m_it );
    }

    // process all new vertices:
    unchecked_states = new_states;
    new_states.clear( );
    
    size_t counter = 0;
    map< ThRightNormalForm , ThRightNormalForm >::const_iterator 
      it = unchecked_states.begin( );
    for( ; it!=unchecked_states.end( ) ; ++it ) {
      
      ThRightNormalForm cur_el = (*it).first;
      ThRightNormalForm cur_conj = (*it).second;

/*
			cout << endl << "  cur_el:" << endl;
			cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
			cout << cur_el;
			cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
*/
      // cout << "c1" << endl;
      // cout << cur_el << endl;
      // set<Permutation> conjugators = getSimpleSummitConjugators( G , cur_el );
      set<Permutation> conjugators = cur_el.getSimpleConjugators( );
      // cout << "c2" << endl;
      
      for( set<Permutation>::const_iterator conj_it = conjugators.begin( ) ; conj_it!=conjugators.end( ) ; ++conj_it ) {
				
	ThRightNormalForm r_mult( theRank , 0 , list<Permutation>( 1 , *conj_it ) );
	if( *conj_it==omega )
	  r_mult = ThRightNormalForm( theRank , 1 , list<Permutation>( 0 ) );
	ThRightNormalForm new_el = -r_mult * cur_el * r_mult;
	ThRightNormalForm new_conj = cur_conj * r_mult;

/*
				cout << "**********************************************" << endl;
				cout << endl << "  new_el:" << endl;
				cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
				cout << new_el;
				cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;

				cout << endl << "  new_conj:" << endl;
				cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
				cout << new_conj;
				cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
*/
				
	// check whether we made a loop
	map< ThRightNormalForm , ThRightNormalForm >::const_iterator 
	  fit = new_states.find( new_el );
	if( fit!=new_states.end( ) ) {
	  result.insert( new_conj * -(*fit).second );
	  continue;
	}
	
	fit = unchecked_states.find( new_el );
	if( fit!=unchecked_states.end( ) ) {
	  result.insert( new_conj * -(*fit).second );
	  continue;
	}
	
	fit = checked_states.find( new_el );
	if( fit!=checked_states.end( ) ) {
	  result.insert( new_conj * -(*fit).second );
	  continue;
	}
	
	cout << "Size = " << new_states.size( ) << "," << unchecked_states.size( ) << "," << checked_states.size( ) << " (" << result.size( ) << ")" << endl;
	new_states[new_el] = new_conj;
	// cout << cur_el << endl;
	// cout << new_el;
	// cout << "+++++++++++++++++++++++++++++++" << endl;

	if( counter!=result.size( ) ) {
	  ofstream of( "out.txt" );
	  set<ThRightNormalForm>::iterator r_it = result.begin();
	  for( ; r_it!=result.end( ) ; ++r_it ) {
	    Word w = shortBraidForm( theRank , (*r_it).getShortWord() );
	    of << w << endl;
	  }
	  counter = result.size( );
	}
	
	
	
	if( new_el.theOmegaPower<0 ) {
	  cerr << "Fail" << endl;
	  exit(1);
	}
      }
    }
  }


  cout << "+++++++++++++++++++++++++++++++" << endl;
  set<ThRightNormalForm>::iterator r_it = result.begin( );
	for( ; r_it!=result.end( ) ; ++r_it ) {
		cout << *r_it << endl;
		ThRightNormalForm n = -(*r_it) * (*this) * (*r_it);
		if( n!=*this ) {
			cout << "Error!" << endl;
			exit( 23 );
		}
	}


  return result;
}


//---------------------------------------------------------------------------//
//------------------------------ randomPositive -----------------------------//
//---------------------------------------------------------------------------//


ThRightNormalForm ThRightNormalForm::randomPositive( int rank , int decomp_length )
{
  list< Permutation > D;
  for( int i=0 ; i<decomp_length ; ++i )
    D.push_back( Permutation::random( rank ) );
  ThRightNormalForm result( rank , 0 , D );
  result.adjust( );
  
  return result;
}


//---------------------------------------------------------------------------//
//------------------------------- increaseRank ------------------------------//
//---------------------------------------------------------------------------//


ThRightNormalForm ThRightNormalForm::increaseRank( int N ) const
{
  if( N<=theRank )
    return *this;
  
  list< Permutation > D1;
  for( list< Permutation >::const_iterator d_it=theDecomposition.begin( ) ; d_it!=theDecomposition.end( ) ; ++d_it )
    D1.push_back( (*d_it).increaseSize(N) );
  
  list< Permutation > D2;
  for( int i=0 ; i<abs(theOmegaPower) ; ++i )
    D2.push_back( Permutation::getHalfTwistPermutation( theRank ).increaseSize(N) );
  
  ThRightNormalForm r1( N , 0 , D1 );
  ThRightNormalForm r2( N , 0 , D2 );
  r1.adjust( );
  r2.adjust( );
  
  return r1 * (theOmegaPower<0 ? -r2 : r2 );
}
