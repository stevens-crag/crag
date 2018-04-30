// Copyright (C) 2005 Alexander Ushakov
// Contents: Implementation of SSS algorithms for class ThRightNormalForm
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//

#include "fstream"
#include "ShortBraidForm.h"


#include "ThRightNormalForm.h"
#include "braid_group.h"
#include "Word.h"


//---------------------------------------------------------------------------//
//---------------------------------- cycle ----------------------------------//
//---------------------------------------------------------------------------//


pair< ThRightNormalForm , ThRightNormalForm > ThRightNormalForm::cycle( ) const
{
  if( theDecomposition.empty( ) )
    return pair< ThRightNormalForm , ThRightNormalForm > ( *this ,   ThRightNormalForm( theRank ) );
  
  ThRightNormalForm result = *this;
  Permutation conj = *(--theDecomposition.end( ));
  if( theOmegaPower%2!=0 )
    conj = conj.flip( );
  result.theDecomposition.push_front( conj );
  result.theDecomposition.pop_back( );
  result.adjust( );
  
  return pair<NF,ThRightNormalForm> ( result , -ThRightNormalForm( conj ) );
}


//---------------------------------------------------------------------------//
//---------------------------------- decycle --------------------------------//
//---------------------------------------------------------------------------//


pair< ThRightNormalForm , ThRightNormalForm > ThRightNormalForm::decycle( ) const
{
  if( theDecomposition.empty( ) )
    return pair< ThRightNormalForm , ThRightNormalForm >( *this , ThRightNormalForm( theRank ) );
  
  ThRightNormalForm result = *this;
  Permutation conj = *theDecomposition.begin( );
  if( theOmegaPower%2!=0 )
    result.theDecomposition.push_back( conj.flip( ) );
  else
    result.theDecomposition.push_back( conj );
  result.theDecomposition.pop_front( );
  result.adjust( );

  return pair< ThRightNormalForm , ThRightNormalForm > ( result , conj );
}




//---------------------------------------------------------------------------//
//--------------------------- findSSSRepresentative -------------------------//
//---------------------------------------------------------------------------//


pair< ThRightNormalForm , ThRightNormalForm > ThRightNormalForm::findSSSRepresentative( ) const 
{
  bool progress;
  
  ThRightNormalForm cur_nf = *this;
  ThRightNormalForm result = *this;
  ThRightNormalForm cur_conj( theRank );
  ThRightNormalForm conjugator( theRank );
  for( progress=true ; progress ; ) {

    progress = false;
    for( int i=0; i<theRank ; ++i ) {

      pair< ThRightNormalForm , ThRightNormalForm > pr = cur_nf.cycle( );
      cur_nf    = pr.first;
      cur_conj *= pr.second;
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
  for( progress=true ; progress ; ) {

    progress = false;
    for( int i=0; i<theRank ; ++i ) {

      pair< ThRightNormalForm , ThRightNormalForm > pr = cur_nf.decycle( );
      cur_nf = pr.first;
      cur_conj *= pr.second;
      if( cur_nf.theDecomposition.size( )<result.theDecomposition.size( ) ) {

	result = cur_nf;
	conjugator = cur_conj;
	progress = true;
	break;
      }
    }
  }
  
  return pair< ThRightNormalForm , ThRightNormalForm >( result , conjugator );
}


//---------------------------------------------------------------------------//
//--------------------------- getSimpleConjugator ---------------------------//
//---------------------------------------------------------------------------//


Permutation ThRightNormalForm::getSimpleConjugator( const Permutation& start ) const
{
  const Permutation omega = Permutation::getHalfTwistPermutation( theRank );
  
  Permutation  c = start;
  Permutation ic = theOmegaPower%2 ? omega * c * omega : c;
    
  Permutation v = c;
  list< Permutation >::const_iterator it = theDecomposition.begin( );
  while( 1 ) {
      
    // move through permutations to the right
    for( it = theDecomposition.begin( ) ; it!=theDecomposition.end( ) ; ++it )
      if( ( v = (*it).inverse() * v.LeftLCM( *it ) ).isTrivial( ) )
	return c;
    
    // if "cur" does not start with c then we make last check if "cur*c" starts with c 
    if( ( v = v.LeftLCM( ic ) )==ic )
      return c;
    
    ic = v;
    c = v = theOmegaPower%2 ? omega * ic * omega : ic;
  }
}


//---------------------------------------------------------------------------//
//------------------------ getSimpleSummitConjugator ------------------------//
//---------------------------------------------------------------------------//


Permutation ThRightNormalForm::getSimpleSummitConjugator( const Permutation& start ) const
{
  const Permutation omega = Permutation::getHalfTwistPermutation( theRank );
  
  Permutation c = getSimpleConjugator( start );
  
  while( 1 ) {
      
    if( c==omega )
      return c;
    
    // Conjugate *this by the current permutation
    ThRightNormalForm conj( c );
    ThRightNormalForm R = -conj * *this * conj;


    // If the current conjugator does not change the cannonical length then stop
    if( R.theDecomposition.size( )==theDecomposition.size( ) )
      return c;


    // If the cannonical length increased then we cycle R
    if( R.theDecomposition.size( )+R.theOmegaPower > 
	theDecomposition.size( )+theOmegaPower ) {
      c = c * (*R.theDecomposition.begin( ));
      continue;
    }
    
    Permutation v = *--R.theDecomposition.end( );
    if( 1-theOmegaPower%2 )
      v = v.flip( );
    v = v.inverse( ) * omega;
    c = c*v;
  }
}



//---------------------------------------------------------------------------//
//--------------------------- getSimpleConjugators --------------------------//
//---------------------------------------------------------------------------//


set< Permutation > ThRightNormalForm::getSimpleConjugators( ) const
{
  set<Permutation> result;
  for( int i=0 ; i<theRank-1 ; ++i ) {
    Permutation c( theRank );
    c.change( i , i+1 );
    result.insert( getSimpleConjugator( c ) );
  }
  return result;
}



//---------------------------------------------------------------------------//
//------------------------- getSimpleSummitConjugators ----------------------//
//---------------------------------------------------------------------------//
// WARNING!!! DOESN'T WORK FOR NON-SUMMIT ELEMENTS!!! (GETS INTO A LOOP)


set< Permutation > ThRightNormalForm::getSimpleSummitConjugators( ) const
{
  set< Permutation > result;
  for( int i=0 ; i<theRank-1 ; ++i ) {
    Permutation c( theRank );
    c.change( i , i+1 );
    result.insert( getSimpleSummitConjugator( c ) );
  }
  return result;
}


//---------------------------------------------------------------------------//
//------------------------- sssConstructionIteration ------------------------//
//---------------------------------------------------------------------------//


pair< bool , ThRightNormalForm > 
ThRightNormalForm::sssConstructionIteration( map< ThRightNormalForm , ThRightNormalForm >& sss_new1 , 
					     map< ThRightNormalForm , ThRightNormalForm >& sss_checked1 , 
					     const map< ThRightNormalForm , ThRightNormalForm >& sss_new2 , 
					     const map< ThRightNormalForm , ThRightNormalForm >& sss_checked2 ) const
{
  const Permutation omega = Permutation::getHalfTwistPermutation( theRank );

  const pair< ThRightNormalForm , ThRightNormalForm > pr = *sss_new1.begin( );
  sss_checked1[pr.first] = pr.second;
  sss_new1.erase( sss_new1.begin( ) );

  // 1. compute the set of simple summit conjugators 
  set< Permutation > conj = pr.first.getSimpleSummitConjugators( );
  for( set< Permutation >::iterator it=conj.begin( ) ; it!=conj.end( ) ; ++it ) {
    
    // Conjugate the current element
    ThRightNormalForm conj( *it );
    ThRightNormalForm res = -conj * pr.first * conj;
    
    // Check if the element is new
    if( sss_new1.find( res )!=sss_new1.end( ) ||
      sss_checked1.find( res )!=sss_checked1.end( ) )
      continue;
    
    // Compute new conjugator
    ThRightNormalForm new_conjugator = pr.second * conj;
    
    // a. Check if the new element belongs to the other SSS. Stop if true.
    map< ThRightNormalForm , ThRightNormalForm >::const_iterator sss_it = sss_new2.find( res );
    if( sss_it!=sss_new2.end( ) )
      return pair< bool , ThRightNormalForm >( true , (*sss_it).second * -new_conjugator );
    
    // b. Check if the new element belongs to the other SSS. Stop if true.
    sss_it = sss_checked2.find( res );
    if( sss_it!=sss_checked2.end( ) )
      return pair< bool , ThRightNormalForm >( true , (*sss_it).second * -new_conjugator );
    
    sss_new1[res] = new_conjugator;
  }
  
  return pair<bool,ThRightNormalForm>( false , ThRightNormalForm( theRank ) );
}


//---------------------------------------------------------------------------//
//------------------------------ areConjugate -------------------------------//
//---------------------------------------------------------------------------//


pair< bool , ThRightNormalForm > ThRightNormalForm::areConjugate( const ThRightNormalForm& rep ) const
{
  if( theRank!=rep.theRank )
    return pair<bool,ThRightNormalForm>( false , ThRightNormalForm( theRank ) );
  
  
  // 1. find representative of SSS1
  pair< ThRightNormalForm , ThRightNormalForm > pr1 =     findSSSRepresentative( );
  pair< ThRightNormalForm , ThRightNormalForm > pr2 = rep.findSSSRepresentative( );
  
  
  // Infimum and supremum of nf1 and nf2 must be the same
  if( pr1.first.theOmegaPower           !=pr2.first.theOmegaPower || 
      pr1.first.theDecomposition.size( )!=pr2.first.theDecomposition.size( ) )
    return pair< bool , ThRightNormalForm >( false , ThRightNormalForm( theRank ) );

  
  // If nf1==nf2 then we are lucky and can stop
  if( pr1.first==pr2.first )
    return pair< bool , ThRightNormalForm >( true , pr2.second * -pr1.second );

  
  // Multiply nf1 and nf2 by omega^2 to make infimums positive
  if( pr1.first.theOmegaPower<0 ) {
    int power = (1-pr1.first.theOmegaPower/2)*2;
    pr1.first.theOmegaPower += power;
    pr2.first.theOmegaPower += power;
  }
  
  
  // 3. construct super summit sets
  map< ThRightNormalForm , ThRightNormalForm > sss_new1;
  map< ThRightNormalForm , ThRightNormalForm > sss_checked1;
  map< ThRightNormalForm , ThRightNormalForm > sss_new2;
  map< ThRightNormalForm , ThRightNormalForm > sss_checked2;
  
  sss_new1[pr1.first] = pr1.second;
  sss_new2[pr2.first] = pr2.second;
  
  
  for( int i=0 ; i<1000000 && sss_new1.size( ) && sss_new2.size( ) ; ++i ) {
    
    pair< bool , ThRightNormalForm > pr = sssConstructionIteration( sss_new1 , sss_checked1 , sss_new2 , sss_checked2 );
    if( pr.first )
      return pair< bool , ThRightNormalForm >( true , pr.second );

    
    pr = sssConstructionIteration( sss_new2 , sss_checked2 , sss_new1 , sss_checked1 );
    if( pr.first )
      return pair< bool , ThRightNormalForm >( true , pr.second.inverse( ) );
  }
  
  return pair<bool,ThRightNormalForm>( false , ThRightNormalForm( theRank ) );
}


