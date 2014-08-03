// Copyright (C) 2007 Alexander Ushakov
// Contents: Implementation of class StraightLineProgramWord
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//

#define DEBUG_WP_AUT


#include "StraightLineProgramWord.h"

#include "set"
using namespace std;


LongInteger gcd(const LongInteger& a , const LongInteger& b )
{
  LongInteger result;
  mpz_gcd(result.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
  return result;
}


//---------------------------------------------------------------------------//
//------------------------------- simplify ----------------------------------//
//---------------------------------------------------------------------------//


void StraightLineProgramWord::extendTerminals( int nt )
{
  if( nt<=theTerminals )
    return;
  
  int shift = nt-theTerminals;
  map< int , Production > newRules;
  for( map< int , Production >::iterator r_it=theRules.begin( ) ; r_it!=theRules.end( ) ; ++r_it ) {
    
    int rl = (*r_it).first;
    Production prod = (*r_it).second;
    
    if( prod.theTerm1<-theTerminals )
      prod.theTerm1 -= shift;
    else if( prod.theTerm1>theTerminals )
      prod.theTerm1 += shift;
    
    if( prod.theTerm2<-theTerminals )
      prod.theTerm2 -= shift;
    else if( prod.theTerm2>theTerminals )
      prod.theTerm2 += shift;
    
    newRules[rl+shift] = prod;
  }

  if( theRoot<-theTerminals )
    theRoot -= shift;
  if( theRoot>theTerminals )
    theRoot += shift;
  theRules = newRules;
  theTerminals = nt;
}


//---------------------------------------------------------------------------//
//------------------------------- simplify ----------------------------------//
//---------------------------------------------------------------------------//


void StraightLineProgramWord::simplify( )
{
  // update vertices of degree 1
  for( map< int , Production >::iterator r_it=theRules.begin( ) ; r_it!=theRules.end( ) ; ++r_it ) {
    int n = (*r_it).first;
    Production& pr = theRules[abs(n)];
    if( pr.theTerm1!=0 && pr.theTerm2==0 ) {
      if( abs(pr.theTerm1)>theTerminals ) {
	Production pr1 = theRules[abs(pr.theTerm1)];
	if( pr.theTerm1>0 )
	  pr = pr1;
	else {
	  pr = pr1;
	  swap( pr.theTerm1 , pr.theTerm2 );
	  pr.theTerm1 = -pr.theTerm1;
	  pr.theTerm2 = -pr.theTerm2;
	}
	pr.theHeight--;
      }
    }
  }

  // find vertices reachable from the root
  set< int > reach;
  reach.insert( theRoot );
  for( set< int >::const_iterator rch_it=reach.end( ) ; rch_it!=reach.begin( ) ; ) {
    int n = *--rch_it;
    Production pr = theRules[abs(n)];
    if( abs(pr.theTerm1)>theTerminals ) reach.insert( abs(pr.theTerm1) );
    if( abs(pr.theTerm2)>theTerminals ) reach.insert( abs(pr.theTerm2) );
  }
  
  // remove vertices that can not be reached from the root
  set< int > non_reach;
  for( map< int , Production >::iterator r_it=theRules.begin( ) ; r_it!=theRules.end( ) ; ++r_it )
    if( reach.find((*r_it).first)==reach.end( ) )
      non_reach.insert((*r_it).first);
  for( set< int >::iterator nr_it=non_reach.begin( ) ; nr_it!=non_reach.end() ; ++nr_it )
    theRules.erase( *nr_it );
}


//---------------------------------------------------------------------------//
//------------------------------ max_rule_number ----------------------------//
//---------------------------------------------------------------------------//


int StraightLineProgramWord::max_rule_number( ) const
{
  if( theRules.empty( ) )
    return theTerminals;
  return (*--theRules.end( )).first;
}


//---------------------------------------------------------------------------//
//------------------------------- order_vertices ----------------------------//
//---------------------------------------------------------------------------//


void StraightLineProgramWord::order_vertices( int n , list< int >& order , set< int >& closure ) const
{
  n = abs(n);
  if( n<theTerminals )
    return;
  if( closure.find(n)!=closure.end() )
    return;
  
  const Production& pr = (*theRules.find(n)).second;
  int T1 = pr.theTerm1;
  int T2 = pr.theTerm2;
  order_vertices( T1 , order , closure );
  order_vertices( T2 , order , closure );
  
  closure.insert( n );
  order.push_front( n );
}


//---------------------------------------------------------------------------//
//------------------------------- assertionType -----------------------------//
//---------------------------------------------------------------------------//


bool StraightLineProgramWord::assertionType( const Assertion& A , const StraightLineProgramWord& SLP ) const
{
  LongInteger L1 = length_rule( A.theVertex1 );
  LongInteger L2 = SLP.length_rule( A.theVertex2 );
  if( A.theBase1 )
    return L1<=A.theLength+L2;
  else
    return L2<=A.theLength+L1;
}


//---------------------------------------------------------------------------//
//---------------------------------- equal ----------------------------------//
//---------------------------------------------------------------------------//


bool StraightLineProgramWord::equal( int n1 , const StraightLineProgramWord& SLP , int n2 ) const
{
  // A. A few fast checks
  LongInteger L1 = length_rule( n1 );
  LongInteger L2 = SLP.length_rule( n2 );
  if( L1!=L2 )
    return false;
  
  if( L1==0 && L2==0 )
    return true;
  
  // B. Order the set of vertices relative to the length of the corresponding words
  typedef pair< bool , int > RI; // (*this?, non-terminal rule number)
  typedef pair< LongInteger , int > PII; // (length, height)
  typedef pair< PII , RI > ORI;  // (length, rule information)
  set< ORI > orderedRules;
  for( map< int , Production >::const_iterator r_it=theRules.begin() ; r_it!=theRules.end( ) ; ++r_it ) {
    int r = (*r_it).first;
    orderedRules.insert( ORI( PII( length_rule(r) , height_rule(r) ) , RI(true,r) ) );
  }
  for( map< int , Production >::const_iterator r_it=SLP.theRules.begin() ; r_it!=SLP.theRules.end( ) ; ++r_it ) {
    int r = (*r_it).first;
    orderedRules.insert( ORI( PII( SLP.length_rule(r) , SLP.height_rule(r) ) , RI(false,r) ) );
  }
  
  // C. Construct the initial set of assertions
  set< Assertion > assertions;
  assertions.insert( Assertion( true , n1 , n2 , 0 ) );
#ifdef DEBUG_WP_AUT
  cout << "(" << n1 << " , " << n2 << ")" << endl;
#endif

  
  // D. Process assertions
  while( !orderedRules.empty( ) ) {
    
    // D.1. Pick a longest unprocessed rule and remove it from the system
    bool firstTerm = (*--orderedRules.end( )).second.first;
    int  theVertex      = (*--orderedRules.end( )).second.second;
    orderedRules.erase( --orderedRules.end( ));


    bool print_flag = false;

    // D.2. Go through all assertions and split those containing the chosen longest rule
    set< Assertion > new_assertions;
    set< Assertion > copy_assertions = assertions;
    for( set< Assertion >::iterator a_it=copy_assertions.begin( ) ; a_it!=copy_assertions.end( ) ; ++a_it ) {
      
      const Assertion& A = *a_it;
      // D.2.a. the assertion has the rule on the 1st place
      if( (  firstTerm && abs(A.theVertex1)==theVertex ) || 
	  ( !firstTerm && abs(A.theVertex2)==theVertex ) ) { 
#ifdef DEBUG_WP_AUT
	if( !print_flag ) {
	  cout << endl;
	  cout << theVertex << ", " << firstTerm << endl;
	  cout << "==== Assertions =======" << endl;
	  for( set< Assertion >::iterator a_it=assertions.begin( ) ; a_it!=assertions.end( ) ; ++a_it )
	    cout << "   " << *a_it << endl;
	  cout << "=======================" << endl;
	  print_flag = true;
	}
	cout << "Assertion to process: " << A << " | " << firstTerm << endl;
#endif
	set< Assertion > na = splitAssertion( A , firstTerm , SLP );
#ifdef DEBUG_WP_AUT
	for( set< Assertion >::const_iterator na_it=na.begin( ) ; na_it!=na.end( ) ; ++na_it )
	  cout << "   -> " << *na_it << endl;
	if( na.empty( ) )
	  cout << "   empty" << endl;
#endif
	new_assertions.insert( na.begin( ) , na.end( ) );
	assertions.erase( A );
	continue;
      }
    }
    
    // D.3. Check if new assertions give a non compact system, compactify, and add them into the system
    for( set< Assertion >::iterator na_it=new_assertions.begin( ) ; na_it!=new_assertions.end( ) ; ++na_it ) {
      
      // check that a new assertion does not belong to the current system of assertions
      Assertion a = *na_it;
      if( assertions.find(a)!=assertions.end() )
	continue;
      // Subword assertion cannot be compactified, so I add it
      if( !assertionType( a , SLP ) ) {
	assertions.insert( a );
	continue;
      }
      
      // find similar overlap assertions
      vector< Assertion > similar_assertions;
      for( set< Assertion >::const_iterator as_it=assertions.begin( ) ; as_it!=assertions.end( ) ; ++as_it )
	if( a.similar( *as_it ) && assertionType( *as_it , SLP ) )
	  similar_assertions.push_back( *as_it );
      
      // if there are 2 assertions, we need to remove one and update the new assertion before adding into the system
      if( similar_assertions.size( )==2 ) {
	LongInteger p1 = similar_assertions[0].theLength;
	LongInteger p2 = similar_assertions[1].theLength;
	LongInteger p3 = a.theLength;

#ifdef DEBUG_WP_AUT
	  cout << "A1 = " << similar_assertions[0] << endl;
	  cout << "A2 = " << similar_assertions[1] << endl;
	  cout << "a  = " << a << endl;
	  cout << "p1 = " << p1 << endl;
	  cout << "p2 = " << p2 << endl;
	  cout << "p3 = " << p3 << endl;
	  cout << "abs(p3-p2) = " << abs(p3-p2) << endl;
	  cout << "abs(p3-p1) = " << abs(p3-p1) << endl;
#endif
	LongInteger p = gcd( abs(p3-p2), abs(p3-p1) );
	assertions.erase(similar_assertions[1]);
	a.theLength = p;
      }
      assertions.insert( a );
    }
  }

#ifdef DEBUG_WP_AUT
  cout << "==== Final Assertions =======" << endl;
  for( set< Assertion >::iterator a_it=assertions.begin( ) ; a_it!=assertions.end( ) ; ++a_it )
    cout << "   " << *a_it << endl;
  cout << "=======================" << endl;
#endif

  for( set< Assertion >::iterator a_it=assertions.begin( ) ; a_it!=assertions.end( ) ; ++a_it )
    if( (*a_it).theVertex1!=(*a_it).theVertex2 )
      return false;

  return true;
}


//---------------------------------------------------------------------------//
//------------------------------- reduce_rule -------------------------------//
//---------------------------------------------------------------------------//


void StraightLineProgramWord::reduce_rule( int n )
{
  // A. Check if there is no rule of the rule is a terminal
  if( abs(n)<=theTerminals )
    return;
  
  // B. Get the rule and check if it is already reduced
  Production& pr = theRules[n];
  if( pr.reduced )
    return;
  
  // C. Make sure that the children are reduced
  if( pr.theTerm1!=0 )
    reduce_rule( abs(pr.theTerm1) );
  if( pr.theTerm2!=0 )
    reduce_rule( abs(pr.theTerm2) );

#ifdef DEBUG_WP_AUT
  cout << "Process: " << n << endl;
  Word W1 = getWord( pr.theTerm1 );
  Word W2 = getWord( pr.theTerm2 );
#endif
  
  LongInteger c = cancellationLength( pr.theTerm1 , pr.theTerm2 );
  if( c!=0 ) {
#ifdef DEBUG_WP_AUT
    cout << "Cancellation: " << n << " -> " << pr.theTerm1 << ", " << pr.theTerm2 << endl;
    cout << "c = " << c << endl;
    cout << "W1 = " << W1 << endl;
    cout << "W2 = " << W2 << endl;
#endif
    pr.theTerm1 = truncateVertexRight( pr.theTerm1 , c );
    pr.theTerm2 = truncateVertexLeft ( pr.theTerm2 , c );
    if( pr.theTerm1==0 ) {
      pr.theTerm1 = pr.theTerm2;
      pr.theTerm2 = 0;
    }

#ifdef DEBUG_WP_AUT
    Word W1_new = getWord( pr.theTerm1 );
    Word W2_new = getWord( pr.theTerm2 );
    if( W1*W2!=W1_new*W2_new || (W1*W2).length( )!=W1_new.length()+W2_new.length() ) {
      cout << (W1*W2).length( ) << " | " << W1_new.length()+W2_new.length() << endl;
      cout << "W1' = " << W1_new << endl;
      cout << "W2' = " << W2_new << endl;
      cout << "Cancellation failure" << endl;
      cout << endl;
      cout << *this << endl;
      exit(1);
    }
#endif
  }
  pr.theLength = length_rule( pr.theTerm1 ) + length_rule( pr.theTerm2 );
  int h1 = height_rule( pr.theTerm1 );
  int h2 = height_rule( pr.theTerm2 );
  pr.theHeight = 1 + ( h1>h2 ? h1 : h2 );
  pr.reduced = true;
}



//---------------------------------------------------------------------------//
//----------------------------- update_lengths ------------------------------//
//---------------------------------------------------------------------------//


void StraightLineProgramWord::update_heights( )
{
  for( map< int , Production >::iterator r_it=theRules.begin( ) ; r_it!=theRules.end( ) ; ++r_it ) 
    (*r_it).second.theHeight = -1;
  
  for( map< int , Production >::iterator r_it=theRules.begin( ) ; r_it!=theRules.end( ) ; ++r_it ) 
    if( (*r_it).second.theHeight==-1 )
      update_height( (*r_it).first );
}

void StraightLineProgramWord::update_height( int n )
{
  // n is a non-terminal
  Production& pr = theRules[n];
  pr.theHeight = 0;

  //check if heights of children are computed
  int T1 = abs(pr.theTerm1);
  if( T1==0 ) {
  } else if( T1<=theTerminals ) {
    pr.theHeight++;
  } else if( T1>theTerminals ) {
    Production& pr1 = theRules[T1];
    if( pr1.theHeight==-1 )
      update_height( T1 );
    pr.theHeight += pr1.theHeight;
  }
  
  int T2 = abs(pr.theTerm2);
  if( T2==0 ) {
  } else if( T2<=theTerminals ) {
    pr.theHeight++;
  } else if( T2>theTerminals ) {
    Production& pr2 = theRules[T2];
    if( pr2.theHeight==-1 )
      update_height( T2 );
    pr.theHeight += pr2.theHeight;
  }
}

//---------------------------------------------------------------------------//
//----------------------------- update_lengths ------------------------------//
//---------------------------------------------------------------------------//


void StraightLineProgramWord::update_lengths( )
{
  for( map< int , Production >::iterator r_it=theRules.begin( ) ; r_it!=theRules.end( ) ; ++r_it ) 
    (*r_it).second.theLength = -1;
  
  for( map< int , Production >::iterator r_it=theRules.begin( ) ; r_it!=theRules.end( ) ; ++r_it ) 
    if( (*r_it).second.theLength==-1 )
      update_length( (*r_it).first );
}


void StraightLineProgramWord::update_length( int n )
{
  // n is a non-terminal
  Production& pr = theRules[n];
  pr.theLength = 0;

  //check if lengths of children are computed
  int T1 = abs(pr.theTerm1);
  if( T1==0 ) {
  } else if( T1<=theTerminals ) {
    pr.theLength++;
  } else if( T1>theTerminals ) {
    Production& pr1 = theRules[T1];
    if( pr1.theLength==-1 )
      update_length( T1 );
    pr.theLength += pr1.theLength;
  }
  
  int T2 = abs(pr.theTerm2);
  if( T2==0 ) {
  } else if( T2<=theTerminals ) {
    pr.theLength++;
  } else if( T2>theTerminals ) {
    Production& pr2 = theRules[T2];
    if( pr2.theLength==-1 )
      update_length( T2 );
    pr.theLength += pr2.theLength;
  }
}


//---------------------------------------------------------------------------//
//----------------------------- applyMapping  -------------------------------//
//---------------------------------------------------------------------------//


StraightLineProgramWord StraightLineProgramWord::applyMapping( const Map& M ) const
{
  StraightLineProgramWord result;
  
  // A. If the CS produces the empty word
  if( theRoot==0 )
    return result;
  
  // B. If not all images are defined
  if( M.domainSize( )<theTerminals )
    return result;
  
  // C. Construct the result
  result.theTerminals = M.rangeSize( );

  // C.1. add new rules
  pair< map< int , Production > , vector< int > > R = map_rules( M );
  result.theRules = R.first;
  
  // C.2 add old rules
  int shift = (*--result.theRules.end()).first-theTerminals;
  for( map< int , Production >::const_iterator r_it=theRules.begin( ) ; r_it!=theRules.end( ) ; ++r_it ) {
    
    int rl = (*r_it).first;
    Production prod = (*r_it).second;
#ifdef DEBUG_WP_AUT
    cout << "      Old rule: " << rl << " -> " << prod.theTerm1 << " , " << prod.theTerm2 << endl;
#endif
    if( prod.theTerm1<-theTerminals ) {
      prod.theTerm1 -= shift;
    } else if( prod.theTerm1<0 ) {
      prod.theTerm1 = -R.second[-1-prod.theTerm1];
    } else if( prod.theTerm1>theTerminals )
      prod.theTerm1 += shift;
    else if( prod.theTerm1>0 ) {
      prod.theTerm1 = R.second[prod.theTerm1-1];
      // cout << "Here: " << prod.theTerm1 << endl;
    }
    
    if( prod.theTerm2<-theTerminals )
      prod.theTerm2 -= shift;
    else if( prod.theTerm2<0 )
      prod.theTerm2 = -R.second[-1-prod.theTerm2];
    else if( prod.theTerm2>theTerminals )
      prod.theTerm2 += shift;
    else if( prod.theTerm2>0 ) {
      prod.theTerm2 = R.second[prod.theTerm2-1];
    }
    
    result.theRules[rl+shift] = prod;
  }
  
  // update the root
  result.theRoot = theRoot+shift;
  
  // update the lengths 
  result.update_lengths( );
  result.update_heights( );
  
  for( map< int , Production >::const_iterator r_it=result.theRules.begin( ) ; r_it!=result.theRules.end( ) ; ++r_it ) {
    int rl = (*r_it).first;
    Production pr = (*r_it).second;
#ifdef DEBUG_WP_AUT
    cout << "   -  " << rl << " -> " << pr.theTerm1 << " , " << pr.theTerm2 << endl;
#endif
  }
  
  return result;
}


//---------------------------------------------------------------------------//
//---------------------------------- map_rules ------------------------------//
//---------------------------------------------------------------------------//


pair< map< int , StraightLineProgramWord::Production > , vector< int > > StraightLineProgramWord::map_rules( const Map& M )
{
  map< int , Production > rules;
  vector< int > domain_elts;
  
  // the current number of the largest rule in the system
  int last_rule = M.rangeSize( );

  // construct rules
  const vector< Word >& I = M.generatingImages( );
  for( int i=0 ; i<I.size( ) ; ++i ) {
    
    const Word& w = I[i];
    if( w.length( )==0 ) {
      rules[++last_rule] = Production( );
      domain_elts.push_back( last_rule );
    } else if( w.length( )==1 ) {
      rules[++last_rule] = Production( *w.begin( ) , 0 , 1 );
      domain_elts.push_back( last_rule );
    } else {
      int cur_len = 2;
      rules[++last_rule] = Production( *w.begin( ) , *++w.begin( ) , cur_len );
      Word::const_iterator w_it = w.begin( );
      for( ++++w_it ; w_it!=w.end( ) ; ++w_it ) {
	rules[last_rule+1] = Production( last_rule , *w_it , ++cur_len );
	++last_rule;
      }
      domain_elts.push_back( last_rule );
    }
  }
  
  /*
  for( map< int , Production >::const_iterator r_it=rules.begin( ) ; r_it!=rules.end( ) ; ++r_it ) {
    int rl = (*r_it).first;
    Production pr = (*r_it).second;
    cout << rl << " -> " << pr.theTerm1 << " , " << pr.theTerm2 << endl;
  }
  */
  
  return pair< map< int , Production > , vector< int > >( rules , domain_elts );
}
