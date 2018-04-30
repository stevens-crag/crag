// Copyright (C) 2005 Alexander Ushakov
// Contents: Implementation of class ThLeftNormalForm
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//


#include <algorithm>
#include "ThLeftNormalForm.h"
#include "ThRightNormalForm.h"
#include "braid_group.h"
#include "Word.h"
#include <time.h>


ostream& operator << ( ostream& os, const ThLeftNormalForm& rep )
{
  const list< Permutation >& decomposition = rep.getDecomposition( );
  
  os << "Power = " << rep.getPower( ) << endl;
  list< Permutation >::const_iterator it = decomposition.begin( );
  for( ; it!=decomposition.end( ) ; ++it )
    os << (*it) << endl;
  os << "Rank = "  << rep.getRank ( ) << endl;

  return os;
}


//---------------------------------------------------------------------------//
//-------------------------- ThRightNormalForm ------------------------------//
//---------------------------------------------------------------------------//


ThLeftNormalForm::operator ThRightNormalForm( ) const
{
  if( theOmegaPower%2==0 ) {
    ThRightNormalForm nf( theRank , theOmegaPower , theDecomposition );
    nf.adjust( );
    return nf;
  }
  
  list< Permutation > D;
  for( list< Permutation >::const_iterator d_it=theDecomposition.begin( ) ; d_it!=theDecomposition.end( ) ; ++d_it )
    D.push_back( (*d_it).flip( ) );
  
  ThRightNormalForm nf( theRank , theOmegaPower , D );
  nf.adjust( );
  return nf;
}


//---------------------------------------------------------------------------//
//------------------------------- compute -----------------------------------//
//---------------------------------------------------------------------------//


bool ThLeftNormalForm::operator==( const ThLeftNormalForm& rep ) const
{
  return 
    theRank          == rep.theRank &&
    theOmegaPower    == rep.theOmegaPower &&
    theDecomposition == rep.theDecomposition;
}


//---------------------------------------------------------------------------//
//--------------------------- ThLeftNormalForm ------------------------------//
//---------------------------------------------------------------------------//


ThLeftNormalForm::ThLeftNormalForm( const crag::braidgroup::BraidGroup& G , const Word& w ) :
  theRank( G.getRank( ) )
{
  Word reversed_w;
  for( Word::const_iterator w_it=w.begin( ) ; w_it!=w.end( ) ; ++w_it )
    reversed_w.push_front( *w_it );
  ThRightNormalForm F( G , reversed_w );
  
  NF pr ( theRank , F.getPower() , F.getDecomposition( ) );
  reverse( pr );
  
  theOmegaPower = pr.second;
  theDecomposition = pr.third;
}


//---------------------------------------------------------------------------//
//------------------------------- reverse -----------------------------------//
//---------------------------------------------------------------------------//


void ThLeftNormalForm::reverse( ThLeftNormalForm::NF& pr )
{
  list<Permutation> rev_list;
  list<Permutation>::iterator it = pr.third.begin( );
  for( ; it!=pr.third.end( ) ; ++it )
    rev_list.push_front( -*it );
  pr.third = rev_list;
}


//---------------------------------------------------------------------------//
//-------------------------- adjustDecomposition ----------------------------//
//---------------------------------------------------------------------------//


void ThLeftNormalForm::adjustDecomposition( int rank , int& power , list<Permutation>& decomp )
{
  NF pr( rank , power , decomp );
  reverse( pr );
  ThRightNormalForm::adjustDecomposition( rank , pr.second , pr.third );
  reverse( pr );
  power = pr.second;
  decomp = pr.third;
}

//---------------------------------------------------------------------------//
//------------------------------- inverse -----------------------------------//
//---------------------------------------------------------------------------//


ThLeftNormalForm ThLeftNormalForm::inverse( ) const
{
  NF pr( theRank , theOmegaPower , theDecomposition );
  reverse( pr );
  ThRightNormalForm rightNF = pr;
  pr = -rightNF;
  reverse( pr );
  return pr;
}


//---------------------------------------------------------------------------//
//------------------------------- multiply ----------------------------------//
//---------------------------------------------------------------------------//


ThLeftNormalForm ThLeftNormalForm::multiply( const ThLeftNormalForm& rep ) const
{
  NF pr1 = *this;
  reverse( pr1 );
  ThRightNormalForm rightNF1 = pr1;

  NF pr2 = rep;
  reverse( pr2 );
  ThRightNormalForm rightNF2 = pr2;
  
  pr1 = rightNF2 * rightNF1;
  reverse( pr1 );
  
  return pr1;
}


//---------------------------------------------------------------------------//
//-------------------------------- getWord ----------------------------------//
//---------------------------------------------------------------------------//

Word ThLeftNormalForm::getWord() const {
  Word result;
  // A. Take care of Deltas
  if( theOmegaPower!=0 ) {
    const Permutation omega = Permutation::getHalfTwistPermutation( theRank );
    Word omegaWord = Word(omega.geodesicWord());
    if( theOmegaPower<0 )
      omegaWord = -omegaWord;
    result = omegaWord.power(abs(theOmegaPower));
  }
  
  // B. Take care of permutation braids
  for (const auto &d : theDecomposition) {
    result *= Word(d.geodesicWord());
  }
  
  return result;
}

Word ThLeftNormalForm::getReducedWord2() const {
  if (theOmegaPower >= 0 || theDecomposition.size() == 0)
    return getWord();

  const auto p = -theOmegaPower;
  const auto d = theDecomposition.size();
  auto a = p < d ? p : d;

  Word result;

  // 1. Process omega
  const Permutation omega = Permutation::getHalfTwistPermutation(theRank);
  Word omegaWord = Word(omega.geodesicWord());
  omegaWord = -omegaWord;
  result = omegaWord.power(p - a);

  // 2. Sort permutations and find the cut_value
  vector<int> sizes;
  for (const auto &perm : getDecomposition())
    sizes.push_back(perm.geodesic().size());
  sort(sizes.begin(), sizes.end(), [](const int &lhs, const int &rhs) { return lhs > rhs; });
  const int cut_value = sizes[a - 1];

  // 3. Cancel omega^-1 with positive permutations
  for (auto it = theDecomposition.begin(); it != theDecomposition.end(); ++it) {
    auto perm = *it;
    if (a > 0 && perm.geodesic().size() >= cut_value) {
      perm = (-perm) * omega;
      if ((--a) % 2 != 0)
        perm = perm.flip();
      result *= -Word(perm.geodesicWord());
    } else {
      if (a % 2 != 0)
        perm = perm.flip();
      result *= Word(perm.geodesicWord());
    }
  }
  return result;
}

Word ThLeftNormalForm::getReducedWord() const {
  if (theOmegaPower >= 0 || theDecomposition.size() == 0)
    return getWord();

  const auto p = -theOmegaPower;
  const auto d = theDecomposition.size();
  const auto a = p < d ? p : d;

  Word result;

  // 1. Process omega
  const Permutation omega = Permutation::getHalfTwistPermutation(theRank);
  Word omegaWord = Word(omega.geodesicWord());
  omegaWord = -omegaWord;
  result = omegaWord.power(p - a);

  // 2. Cancel omega^-1 with positive permutations
  auto it = theDecomposition.begin();
  for (int i = 0; i < a; ++i, ++it) {
    auto perm = (-(*it)) * omega;
    if ((a - i - 1) % 2 != 0)
      perm = perm.flip();
    result *= -Word(perm.geodesicWord());
  }
  // 3. process the rest of positive permutations
  for (; it != theDecomposition.end(); ++it) {
    result *= Word((*it).geodesicWord());
  }

  return result;
}

//---------------------------------------------------------------------------//
//----------------------------------- cycle ---------------------------------//
//---------------------------------------------------------------------------//


pair< ThLeftNormalForm , ThLeftNormalForm > ThLeftNormalForm::cycle( ) const
{
  if( theDecomposition.empty( ) )
    return pair< ThLeftNormalForm , ThLeftNormalForm >( *this , ThLeftNormalForm( theRank ) );

  ThLeftNormalForm result = *this;
  Permutation first = *(theDecomposition.begin( ));
  if( theOmegaPower%2!=0 )
    first = first.flip( );
  result.theDecomposition.push_back( first );
  result.theDecomposition.pop_front( );
  result.adjust( );
  
  return pair< ThLeftNormalForm , ThLeftNormalForm >( result , ThLeftNormalForm(first) );
}


//---------------------------------------------------------------------------//
//---------------------------------- decycle --------------------------------//
//---------------------------------------------------------------------------//


pair< ThLeftNormalForm , ThLeftNormalForm > ThLeftNormalForm::decycle( ) const
{
  if( theDecomposition.empty( ) )
    return pair< ThLeftNormalForm , ThLeftNormalForm >( *this , ThLeftNormalForm( theRank ) );
  
  ThLeftNormalForm result = *this;
  Permutation last = *(--theDecomposition.end( ));
  if( theOmegaPower%2==0 )
    result.theDecomposition.push_front( last );
  else 
    result.theDecomposition.push_front( last.flip( ) );
  result.theDecomposition.pop_back( );
  result.adjust( );
  
  return pair< ThLeftNormalForm , ThLeftNormalForm >( result , -ThLeftNormalForm(last) );
}


//---------------------------------------------------------------------------//
//--------------------------- findSSSRepresentative -------------------------//
//---------------------------------------------------------------------------//


pair< ThLeftNormalForm , ThLeftNormalForm > ThLeftNormalForm::findSSSRepresentative( ) const
{
  bool progress;
  
  ThLeftNormalForm result = *this;
  ThLeftNormalForm cur_nf = *this;
  ThLeftNormalForm cur_conj( theRank );
  ThLeftNormalForm conjugator( theRank );
  for( progress=true ; progress ; ) {
    
    progress = false;
    for( int i=0; i<theRank ; ++i ) {
      
      pair< ThLeftNormalForm , ThLeftNormalForm > pr = cur_nf.cycle( );
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
      
      pair< ThLeftNormalForm , ThLeftNormalForm > pr = cur_nf.decycle( );
      cur_nf    = pr.first;
      cur_conj *= pr.second;
      if( cur_nf.theDecomposition.size( )<result.theDecomposition.size( ) ) {
	
        result = cur_nf;
	conjugator = cur_conj;
        progress = true;
        break;
      }
    }
  }

  return pair< ThLeftNormalForm , ThLeftNormalForm >( result , conjugator );
}


//---------------------------------------------------------------------------//
//---------------------------- getSimpleConjugator --------------------------//
//---------------------------------------------------------------------------//


Permutation ThLeftNormalForm::getSimpleConjugator( const Permutation& start ) const
{
  ThRightNormalForm rightNF = *this;
  return rightNF.getSimpleConjugator( start );
}


//---------------------------------------------------------------------------//
//--------------------------- getSimpleConjugators --------------------------//
//---------------------------------------------------------------------------//


set< Permutation > ThLeftNormalForm::getSimpleConjugators( ) const
{
  ThRightNormalForm rightNF = *this;
  return rightNF.getSimpleConjugators( );
}


//---------------------------------------------------------------------------//
//------------------------ getSimpleSummitConjugator ------------------------//
//---------------------------------------------------------------------------//


Permutation ThLeftNormalForm::getSimpleSummitConjugator( const Permutation& start ) const
{
  ThRightNormalForm rightNF = *this;
  return rightNF.getSimpleSummitConjugator( start );
}


//---------------------------------------------------------------------------//
//------------------------- getSimpleSummitConjugators ----------------------//
//---------------------------------------------------------------------------//


set< Permutation > ThLeftNormalForm::getSimpleSummitConjugators( ) const
{
  ThRightNormalForm rightNF = *this;
  return rightNF.getSimpleSummitConjugators( );
}


//---------------------------------------------------------------------------//
//------------------------------ areConjugate -------------------------------//
//---------------------------------------------------------------------------//


pair< bool , ThLeftNormalForm > ThLeftNormalForm::areConjugate( const ThLeftNormalForm& rep ) const
{
  NF pr1 = *this;
  reverse( pr1 );
  ThRightNormalForm rightNF1 = pr1;

  NF pr2 = rep;
  reverse( pr2 );
  ThRightNormalForm rightNF2 = pr2;
  
  pair< bool , ThRightNormalForm > res = rightNF1.areConjugate( rightNF2 );
  NF conj = res.second;
  reverse( conj );
  return pair< bool , ThLeftNormalForm >( res.first , -ThLeftNormalForm(conj) );
}


//---------------------------------------------------------------------------//
//------------------------------ getTransport -------------------------------//
//---------------------------------------------------------------------------//


Permutation ThLeftNormalForm::getTransport( const ThLeftNormalForm& B , const Permutation& u ) const
{
  // we assume that the braid is not a power of \Delta
  Permutation A1 = *  theDecomposition.begin( );
  Permutation B1 = *B.theDecomposition.begin( );
  
  if( getPower()%2!=0 ) {
    A1 = A1.flip( );
    B1 = B1.flip( );
  }
  // return A1*u*-B1;
  return -A1*u*B1;
}


//---------------------------------------------------------------------------//
//-------------------------- findUSSRepresentative --------------------------//
//---------------------------------------------------------------------------//


triple< ThLeftNormalForm , ThLeftNormalForm , int > ThLeftNormalForm::findUSSRepresentative( ) const
{
  ThLeftNormalForm result = *this;
  ThLeftNormalForm conjugator( theRank );
  
  // Compute SSS representative
  pair< ThLeftNormalForm , ThLeftNormalForm > sss = findSSSRepresentative( );
  conjugator = sss.second;
  result = sss.first;
  
  // Cycle until get into a loop
  map< ThLeftNormalForm , pair< ThLeftNormalForm , int > > trajectory;
  trajectory[result] = pair< ThLeftNormalForm , int >( conjugator , 0 );
  for( int c=0 ; 1 ; ++c ) {
    
    pair< ThLeftNormalForm , ThLeftNormalForm > pr = result.cycle( );
    result = pr.first;
    conjugator *= pr.second;
    map< ThLeftNormalForm , pair< ThLeftNormalForm , int > >::iterator t_it = trajectory.find( result );
    if( t_it==trajectory.end( ) )
      trajectory[result] = pair< ThLeftNormalForm , int >( conjugator , c+1 );
    else
      return triple< ThLeftNormalForm , ThLeftNormalForm , int >( result , (*t_it).second.first, c+1-(*t_it).second.second );
  }
}


//---------------------------------------------------------------------------//
//------------------------------- getPullback -------------------------------//
//---------------------------------------------------------------------------//


Permutation ThLeftNormalForm::getPullback( const Permutation& s ) const
{
  Permutation result;
  Permutation Delta = Permutation::getHalfTwistPermutation( theRank );

  Permutation B1 = *theDecomposition.begin( );
  Permutation b1 = s;
  if( theOmegaPower%2!=0 ) {
    B1 = B1.flip( );
    b1 = b1.flip( );
  }
  B1 = Delta * -B1;
  Permutation b0 = -B1 * B1.LeftLCM( s.flip( ) );

  for( list< Permutation >::const_iterator d_it=++theDecomposition.begin( ) ; d_it!=theDecomposition.end( ) ; ++d_it ) {
    const Permutation& B = *d_it;
    b1 = -B * B.LeftLCM( b1 );
  }
  
  return getSimpleSummitConjugator( b0.LeftLCM( b1 ) );
}


//---------------------------------------------------------------------------//
//------------------------------ getTransports ------------------------------//
//---------------------------------------------------------------------------//


set< Permutation > ThLeftNormalForm::getTransports( const Permutation& u , int period ) const
{
  Permutation cur_perm = u;
  ThLeftNormalForm A = *this;
  ThLeftNormalForm c( u );
  ThLeftNormalForm B = -c * A * c;
  
  // set of transports to be considered
  set< Permutation > F;
  
  // trajectory contains elements obtained by consequent cycling of the current element
  // first=trajectory element, second.first=transport, second.second=loop number
  map< Permutation , int > trajectory;
  trajectory[cur_perm] = 0;
  for( int i=1 ; 1 ; ++i ) {
    
    // cycle and compute transports
    for( int t=0 ; t<period ; ++t ) {
      cur_perm = A.getTransport( B , cur_perm );
      A = A.cycle( ).first;
      B = B.cycle( ).first;
    }
    
    // check if got into a loop
    map< Permutation , int >::iterator t_it = trajectory.find( cur_perm );
    if( t_it!=trajectory.end( ) ) {
      
      int loop_start = (*t_it).second;
      for( t_it=trajectory.begin( ) ; t_it!=trajectory.end( ) ; ++t_it )
	if( (*t_it).second>=loop_start )
	  F.insert( (*t_it).first );
      break;
    }
    trajectory[cur_perm] = i;
  }
  
  return F;
}


//---------------------------------------------------------------------------//
//------------------------------ getMainPullback ----------------------------//
//---------------------------------------------------------------------------//


Permutation ThLeftNormalForm::getMainPullback( const Permutation& s , int period ) const
{
  // A. Compute the cycling trajectory for *this
  vector< ThLeftNormalForm > cycle_trajectory( period );
  int j = 1;
  cycle_trajectory[0] = *this;
  for( ThLeftNormalForm nf=cycle( ).first ; nf!=*this ; nf=nf.cycle( ).first , ++j )
    cycle_trajectory[j] = nf;
  
  // B. Compute the permutation pullback trajectory
  map< Permutation , int > trajectory;
  Permutation cur_perm = s;
  trajectory[cur_perm] = 0;
  for( j=1 ; ; ++j ) {
    
    for( int i=0 ; i<period ; ++i )
      cur_perm = cycle_trajectory[period-i-1].getPullback( cur_perm );
    
    map< Permutation , int >::iterator t_it=trajectory.find( cur_perm );
    if( t_it!=trajectory.end( ) ) {
      
      int loop_start = (*t_it).second;
      int l = j-loop_start;
      for( t_it=trajectory.begin( ) ; t_it!=trajectory.end( ) ; ++t_it )
	if( (*t_it).second>=loop_start && (*t_it).second%l==0 )
	  return (*t_it).first;
    }
    trajectory[cur_perm] = j;
  }
}


//---------------------------------------------------------------------------//
//------------------------------ computePeriod ------------------------------//
//---------------------------------------------------------------------------//


int ThLeftNormalForm::computePeriod( ) const
{
  int result = 1;
  for( ThLeftNormalForm nf=cycle( ).first ; nf!=*this ; nf=nf.cycle( ).first, ++result );
  return result;
}


//---------------------------------------------------------------------------//
//------------------------- getSimpleUltraConjugator ------------------------//
//---------------------------------------------------------------------------//


pair< Permutation , bool > ThLeftNormalForm::getSimpleUltraConjugator( int period , const Permutation& start ) const
{
  Permutation c = getSimpleSummitConjugator( start );
  set< Permutation > F = getTransports( c , period );
  
  for( set< Permutation >::const_iterator F_it=F.begin( ) ; F_it!=F.end( ) ; ++F_it ) {
    Permutation v = *F_it;
    if( (-start*v).length( )+start.length()==v.length( ) )
      return pair< Permutation , bool >( v , true );
  }
  
  bool minimal = false;
  if( F.size( )==1 && (*F.begin( )).isTrivial( ) )
    minimal = true;
  
  // Compute the main pullback
  Permutation main_pullback = getMainPullback( start , period );
  set< Permutation > G = getTransports( main_pullback , period );
  
  for( set< Permutation >::const_iterator G_it=G.begin( ) ; G_it!=G.end( ) ; ++G_it ) {
    Permutation v = *G_it;
    if( (-start*v).length( )+start.length()==v.length( ) )
      return pair< Permutation , bool >( v , minimal );
  }
}


//---------------------------------------------------------------------------//
//------------------------- getSimpleUltraConjugators -----------------------//
//---------------------------------------------------------------------------//


set< pair< Permutation , bool > > ThLeftNormalForm::getSimpleUltraConjugators( int period ) const
{
  set< pair< Permutation , bool > > result;
  for( int i=0 ; i<theRank-1 ; ++i ) {
    Permutation start( theRank );
    start.change( i , i+1 );
    result.insert( getSimpleUltraConjugator( period , start ) );
  }
  
  return result;
}


//---------------------------------------------------------------------------//
//---------------------------------- ussAdd ---------------------------------//
//---------------------------------------------------------------------------//


pair< bool , ThLeftNormalForm > 
ThLeftNormalForm::ussAddTrajectory( const ThLeftNormalForm& new_elt , const ThLeftNormalForm& new_conj , 
				    map< ThLeftNormalForm , pair< ThLeftNormalForm , int > >& uss_new1 , 
				    map< ThLeftNormalForm , pair< ThLeftNormalForm , int > >& uss_checked1 , 
				    const map< ThLeftNormalForm , pair< ThLeftNormalForm , int > >& uss_new2 , 
				    const map< ThLeftNormalForm , pair< ThLeftNormalForm , int > >& uss_checked2 ) const
{
  pair< bool , ThLeftNormalForm >  result;
  
  // Add the whole trajectory to already constructed part of USS
  ThLeftNormalForm cur_res = new_elt;
  ThLeftNormalForm new_conjugator = new_conj;
  do {
    
    // a. Check if the new element belongs to the other USS. Stop if true.
    map< ThLeftNormalForm , pair< ThLeftNormalForm , int > >::const_iterator uss_it = uss_new2.find( cur_res );
    if( uss_it!=uss_new2.end( ) )
      return pair< bool , ThLeftNormalForm >( true , (*uss_it).second.first * -new_conjugator );
    
    // b. Check if the new element belongs to the other USS. Stop if true.
    uss_it = uss_checked2.find( cur_res );
    if( uss_it!=uss_checked2.end( ) )
      return pair< bool , ThLeftNormalForm >( true , (*uss_it).second.first * -new_conjugator );
    
    uss_new1[cur_res] = pair< ThLeftNormalForm , int >( new_conjugator , cur_res.computePeriod( ) );
    
    pair< ThLeftNormalForm , ThLeftNormalForm > pr_res = cur_res.cycle( );
    cur_res = pr_res.first;
    new_conjugator *= pr_res.second;
  } while( cur_res!=new_elt );
  
  return result;
}

//---------------------------------------------------------------------------//
//------------------------- sssConstructionIteration ------------------------//
//---------------------------------------------------------------------------//


pair< bool , ThLeftNormalForm > 
ThLeftNormalForm::ussConstructionIteration( map< ThLeftNormalForm , pair< ThLeftNormalForm , int > >& uss_new1 , 
					    map< ThLeftNormalForm , pair< ThLeftNormalForm , int > >& uss_checked1 , 
					    const map< ThLeftNormalForm , pair< ThLeftNormalForm , int > >& uss_new2 , 
					    const map< ThLeftNormalForm , pair< ThLeftNormalForm , int > >& uss_checked2 ) const
{
  const Permutation omega = Permutation::getHalfTwistPermutation( theRank );
  
  const pair< ThLeftNormalForm , pair< ThLeftNormalForm , int > > pr = *uss_new1.begin( );
  uss_checked1[pr.first] = pr.second;
  uss_new1.erase( uss_new1.begin( ) );
  
  // 1. compute the set of simple summit conjugators 
  set< pair< Permutation , bool > > conj = pr.first.getSimpleUltraConjugators( pr.second.second );
  for( set< pair< Permutation , bool > > ::iterator it=conj.begin( ) ; it!=conj.end( ) ; ++it ) {
    
    // Conjugate the current element
    ThLeftNormalForm conj( (*it).first );
    ThLeftNormalForm res = -conj * pr.first * conj;
    
    // Check if the element is new
    if( uss_new1.find( res )!=uss_new1.end( ) ||
      uss_checked1.find( res )!=uss_checked1.end( ) )
      continue;
    
    // Compute new conjugator
    ThLeftNormalForm new_conjugator = pr.second.first * conj;
    
    pair< bool , ThLeftNormalForm > traj_add_result = 
      ussAddTrajectory( res , new_conjugator , uss_new1 , uss_checked1 , uss_new2 , uss_checked2 );
    if( traj_add_result.first )
      return traj_add_result;
  }
  
  return pair<bool,ThLeftNormalForm>( false , ThLeftNormalForm( theRank ) );
}


//---------------------------------------------------------------------------//
//----------------------------- areConjugate_uss ----------------------------//
//---------------------------------------------------------------------------//


pair< bool , ThLeftNormalForm > ThLeftNormalForm::areConjugate_uss( const ThLeftNormalForm& rep , int time_sec_bound ) const
{
  if( theRank!=rep.theRank )
    return pair< bool , ThLeftNormalForm >( false , ThLeftNormalForm( theRank ) );
  
  int init_time = time(0);
  
  // 1. find representative of SSS1
  triple< ThLeftNormalForm , ThLeftNormalForm , int > tr1 =     findUSSRepresentative( );
  triple< ThLeftNormalForm , ThLeftNormalForm , int > tr2 = rep.findUSSRepresentative( );
  
  
  // Infimum and supremum of nf1 and nf2 must be the same
  if( tr1.first.theOmegaPower           !=tr2.first.theOmegaPower || 
      tr1.first.theDecomposition.size( )!=tr2.first.theDecomposition.size( ) )
    return pair< bool , ThLeftNormalForm >( false , ThLeftNormalForm( theRank ) );
  
  
  // If nf1==nf2 then we are lucky and can stop
  if( tr1.first==tr2.first )
    return pair< bool , ThLeftNormalForm >( true , tr2.second * -tr1.second );
  
  // 3. construct ultra summit sets
  map< ThLeftNormalForm , pair< ThLeftNormalForm , int > > uss_new1;
  map< ThLeftNormalForm , pair< ThLeftNormalForm , int > > uss_checked1;
  map< ThLeftNormalForm , pair< ThLeftNormalForm , int > > uss_new2;
  map< ThLeftNormalForm , pair< ThLeftNormalForm , int > > uss_checked2;
  
  ussAddTrajectory( tr1.first , tr1.second , uss_new1 , uss_checked1 , uss_new2 , uss_checked2 );
  pair< bool , ThLeftNormalForm > traj_add_result = 
    ussAddTrajectory( tr2.first , tr2.second , uss_new2 , uss_checked2 , uss_new1 , uss_checked1 );
  if( traj_add_result.first )
    return pair< bool , ThLeftNormalForm >( true , -traj_add_result.second );

  // uss_new1[tr1.first] = pair< ThLeftNormalForm , int >( tr1.second , tr1.third );
  // uss_new2[tr2.first] = pair< ThLeftNormalForm , int >( tr2.second , tr2.third );

  for( int i=0 ; i<1000000 && uss_new1.size( ) && uss_new2.size( ) ; ++i ) {

    // Check if we have time to proceed
    int cur_time = time(0);
    if( cur_time-init_time>time_sec_bound )
      return pair< bool , ThLeftNormalForm >( false , ThLeftNormalForm( theRank ) );
    
    pair< bool , ThLeftNormalForm > pr = ussConstructionIteration( uss_new1 , uss_checked1 , uss_new2 , uss_checked2 );
    if( pr.first ) {
      // cout << "Size: " << uss_new1.size()+uss_checked1.size() << "  -  " << uss_new2.size()+uss_checked2.size() << endl;
      return pair< bool , ThLeftNormalForm >( true , pr.second );
    }

    // Check if we have time to proceed
    cur_time = time(0);
    if( cur_time-init_time>time_sec_bound )
      return pair< bool , ThLeftNormalForm >( false , ThLeftNormalForm( theRank ) );

    pr = ussConstructionIteration( uss_new2 , uss_checked2 , uss_new1 , uss_checked1 );
    if( pr.first ) {
      // cout << "Size: " << uss_new1.size()+uss_checked1.size() << "  -  " << uss_new2.size()+uss_checked2.size() << endl;
      return pair< bool , ThLeftNormalForm >( true , pr.second.inverse( ) );
    }
  }
  
  return pair< bool , ThLeftNormalForm >( false , ThLeftNormalForm( theRank ) );
}


//---------------------------------------------------------------------------//
//------------------------------- increaseRank ------------------------------//
//---------------------------------------------------------------------------//


ThLeftNormalForm ThLeftNormalForm::increaseRank( int N ) const
{
  if( N<=theRank )
    return *this;
  
  list< Permutation > D1;
  for( list< Permutation >::const_iterator d_it=theDecomposition.begin( ) ; d_it!=theDecomposition.end( ) ; ++d_it )
    D1.push_back( (*d_it).increaseSize(N) );
  
  list< Permutation > D2;
  for( int i=0 ; i<abs(theOmegaPower) ; ++i )
    D2.push_back( Permutation::getHalfTwistPermutation( theRank ).increaseSize(N) );
  
  ThLeftNormalForm r1( N , 0 , D1 );
  ThLeftNormalForm r2( N , 0 , D2 );
  r1.adjust( );
  r2.adjust( );
  
  return (theOmegaPower<0 ? -r2 : r2 ) * r1;
}
