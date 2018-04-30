// Copyright (C) 2006 Alexander Ushakov
// Contents: Implementation of class ThRightNormalForm / USS algorithms
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//


#include "ThRightNormalForm.h"
#include "braid_group.h"
#include "Word.h"


//---------------------------------------------------------------------------//
//-------------------------------- getTransport -----------------------------//
//---------------------------------------------------------------------------//


Permutation ThRightNormalForm::getTransport( const ThRightNormalForm& B , const Permutation& u ) const
{
  // we assume that the braid is not a power of \Delta
  Permutation A1 = *(--  theDecomposition.end( ));
  Permutation B1 = *(--B.theDecomposition.end( ));

  if( getPower()%2!=0 ) {
    A1 = A1.flip( );
    B1 = B1.flip( );
  }
  
  return A1*u*(-B1);
}


//---------------------------------------------------------------------------//
//------------------------------ getTransports ------------------------------//
//---------------------------------------------------------------------------//


set< Permutation > ThRightNormalForm::getTransports( const Permutation& u , int period ) const
{
  Permutation cur_perm = u;

  // to be cycled below
  ThRightNormalForm A = *this;
  
  
  // conjugate and get B to be cycled above
  ThRightNormalForm c( u );
  ThRightNormalForm B = -c * A * c;
  
  
  // set of transports to be considered
  set< Permutation > F;
  
  
  // if( A==B ) cout << "==" << endl;

  // trajectory contains elements obtained by consequent cycling of the current element
  // first=trajectory element, second.first=transport, second.second=loop number
  map< ThRightNormalForm , pair< Permutation , int > > trajectory;
  trajectory[B] = pair< Permutation , int >( cur_perm , 0 );
  for( int i=1 ; 1 ; ++i ) {
    
    // cycle and compute transports
    for( int t=0 ; t<period ; ++t ) {
      cur_perm = A.getTransport( B , cur_perm );
      B = A.cycle( ).first;
      A = B.cycle( ).first;
    }
    
    // check if got into a loop
    map< ThRightNormalForm , pair< Permutation , int > >::iterator t_it = trajectory.find( B );
    if( t_it!=trajectory.end( ) ) {
      
      int loop_start = (*t_it).second.second;
      // cout << "   : " << loop_start << " , " << i << endl;
      for( t_it=trajectory.begin( ) ; t_it!=trajectory.end( ) ; ++t_it ) {
	
	if( (*t_it).second.second<loop_start )
	  continue;
	// cout << "       Candidate = " << (*t_it).second.first << endl;
	F.insert( (*t_it).second.first );
      }
      break;
    }
    trajectory[B] = pair< Permutation , int >( cur_perm , i );
  }
  
  return F;
}


//---------------------------------------------------------------------------//
//-------------------------------- getPullback ------------------------------//
//---------------------------------------------------------------------------//


Permutation ThRightNormalForm::getPullback( const Permutation& s ) const
{
  const Permutation omega = Permutation::getHalfTwistPermutation( theRank );

  // we assume that the braid is not a power of \Delta
  Permutation A1 = *(--theDecomposition.end( ));
  if( getPower()%2!=0 )
    A1 = A1.flip( );
  
  // Compute the minimal permutation such that u1 will start from s
  Permutation s1 = -A1 * A1.LeftLCM( s );

  Permutation p = s;
  for( list< Permutation >::const_iterator d_it=--theDecomposition.end( ) ; d_it!=theDecomposition.begin( ) ; ) {
    Permutation A = *--d_it;
    transform( theRank , A , p );
    p = A;
  }
  cout << "      s1 = " << s1 << endl;
  cout << "      p  = " << p << endl;
  s1 = s1.LeftLCM( p );

  cout << "   " << s << " + " << A1 << " = " << s1 << endl;

  // Complete s1 to a simple summit conjugator
  Permutation s2 = getSimpleSummitConjugator( s1 );
  
  return s2;
}


//---------------------------------------------------------------------------//
//------------------------------ getPullback --------------------------------//
//---------------------------------------------------------------------------//


Permutation ThRightNormalForm::getPullback( const Permutation& u , int period ) const
{
  // A. Compute the cycling trajectory for *this
  ThRightNormalForm cur = *this;
  list< ThRightNormalForm > trajectory1;
  for( int i=0 ; i<period ; ++i )
    trajectory1.push_back( cur = cur.cycle( ).first );

  // B. Compute pullback sequence (test)
  int j=0;
  Permutation pb = u;
  cout << "0) " << pb << endl;
  for( list< ThRightNormalForm >::iterator t_it=trajectory1.end( ) ; t_it!=trajectory1.begin( ) ; ) {
    ThRightNormalForm base = *--t_it;
    pb = base.getPullback( pb );
    cout << ++j << ") " << pb << endl;
  }
  
  if( pb==u )
    cout << "coincide" << endl;
  else
    cout << "does not coincide" << endl;
  
  ThRightNormalForm a( pb );
  cout << (-a * *this * a) << endl;
  
  return pb;
}




//---------------------------------------------------------------------------//
//------------------------ getSimpleUltraConjugator -------------------------//
//---------------------------------------------------------------------------//


triple< Permutation , bool , int > ThRightNormalForm::getSimpleUltraConjugator( int period , const Permutation& start )
{
  // Compute transports
  set< Permutation > F = getTransports( start , period );
  
  cout << "c = " << start << endl;
  cout << "    F = " << endl;
  for( set< Permutation >::iterator F_it = F.begin( ) ; F_it!=F.end( ) ; ++F_it ) {
    cout << "     " << *F_it << endl;
    
    Permutation p = -start * *F_it;
    if( start.length( )+p.length( )==(*F_it).length( ) ) {
      
      int p=1;
      ThRightNormalForm f = *F_it;
      ThRightNormalForm B = -f * *this * f;
      for( ThRightNormalForm B1 = B.cycle( ).first ; B1!=B ; ++p )
	B1 = B1.cycle( ).first;

      cout << "         success: period = " << p << endl;
      // need to find period before returning
      return triple< Permutation , bool , int >( *F_it , true , p );
    }
  }
  
  if( F.size()==1 && *F.begin( )==Permutation( theRank ) ) {
    cout << "         need pull_back" << endl;
    
    // Compute a pullbacks
    Permutation pb = getPullback( start , period );

    // Find and element in the F set starting from s

  } else {
    cout << "         not-minimal, need pull_back" << endl;
  }
  
  
  // Compute pull-back
  
  
  return triple< Permutation , bool , int >( );
  //return triple< Permutation , bool , int >( result , true , 0 );
}


//---------------------------------------------------------------------------//
//------------------------ getSimpleUltraConjugators ------------------------//
//---------------------------------------------------------------------------//


set< pair< Permutation , int > > ThRightNormalForm::getSimpleUltraConjugators( int period )
{
  set< pair< Permutation , int > > result;

  // Check if rep is a power of Delta
  if( theDecomposition.size( )==0 )
    return result;
  
  // Otherwise we compute the set of simple summit conjugators and try to lift them up to simple ultra conjugators
  set< Permutation > conjugators = getSimpleSummitConjugators( );
  for( set<Permutation>::const_iterator c_it=conjugators.begin( ) ; c_it!=conjugators.end( ) ; ++c_it ) {

    // cout << "c = " << (*c_it) << endl;
    triple< Permutation , bool , int > suc = getSimpleUltraConjugator( period , *c_it );
  }
  
  return result;
}


//---------------------------------------------------------------------------//
//--------------------------- findUSSRepresentative -------------------------//
//---------------------------------------------------------------------------//


triple< ThRightNormalForm , ThRightNormalForm , int > ThRightNormalForm::findUSSRepresentative( ) const
{
  ThRightNormalForm result = *this;
  ThRightNormalForm conjugator( theRank );
  
  // Compute SSS representative
  pair< ThRightNormalForm , ThRightNormalForm > sss = findSSSRepresentative( );
  conjugator = sss.second;
  result = sss.first;
  
  // Cycle until get into a loop
  map< ThRightNormalForm , pair< ThRightNormalForm , int > > trajectory;
  trajectory[result] = pair<ThRightNormalForm,int>(conjugator,0);
  for( int c=0 ; 1 ; ++c ) {
    
    pair< ThRightNormalForm , ThRightNormalForm > pr = result.cycle( );
    result = pr.first;
    conjugator *= pr.second;
    map< ThRightNormalForm , pair< ThRightNormalForm , int > >::iterator t_it = trajectory.find( result );
    if( t_it==trajectory.end( ) )
      trajectory[result] = pair<ThRightNormalForm,int>(conjugator,c+1);
    else
      return triple< ThRightNormalForm , ThRightNormalForm , int >( result , (*t_it).second.first, c+1-(*t_it).second.second );
  }
}


