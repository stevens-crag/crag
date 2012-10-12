// Copyright (C) 2007 Alexander Ushakov
// Contents: Implementation of class StraightLineProgramWord
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//


#include "StraightLineProgramWord.h"

#include "set"
using namespace std;


//---------------------------------------------------------------------------//
//-------------------------- StraightLineProgramWord ------------------------//
//---------------------------------------------------------------------------//


StraightLineProgramWord::StraightLineProgramWord( int init , int gens ) : 
  theRoot(gens+1),
  theTerminals(gens) 
{ 
  theRules[gens+1] = Production( init , 0 , 1 , 1 );
}


//---------------------------------------------------------------------------//
//------------------------------ operator * ---------------------------------//
//---------------------------------------------------------------------------//


StraightLineProgramWord StraightLineProgramWord::operator * ( const StraightLineProgramWord& SLP ) const
{
  StraightLineProgramWord result = *this;
  StraightLineProgramWord SLP2 = SLP;
  if( theTerminals>SLP.theTerminals )
    SLP2.extendTerminals( theTerminals );
  if( theTerminals<SLP.theTerminals )
    result.extendTerminals( SLP.theTerminals );
  
  SLP2.extendTerminals( result.max_rule_number( ) );
  SLP2.theTerminals = result.theTerminals;
  for( map< int , Production >::iterator r_it=SLP2.theRules.begin() ; r_it!=SLP2.theRules.end( ) ; ++r_it )
    result.theRules.insert( *r_it );
  int h1 = height_rule( theRoot );
  int h2 = SLP2.height_rule( SLP2.theRoot );
  int h = 1+(h1>h2 ? h1:h2);
  result.theRules[1+result.max_rule_number( )] = Production( theRoot , 
							     SLP2.theRoot , 
							     length_rule(theRoot)+SLP2.length_rule(SLP2.theRoot) ,
							     h );
  
  result.theRoot = result.max_rule_number( );
  return result;
}


//---------------------------------------------------------------------------//
//-------------------------------- operator- --------------------------------//
//---------------------------------------------------------------------------//


StraightLineProgramWord StraightLineProgramWord::operator- ( ) const
{
  StraightLineProgramWord result = *this;
  // Empty word
  if( theRoot==0 )
    return result;
  
  // Terminal (generator)
  if( abs(theRoot)<=theTerminals ) {
    result.theRoot = -result.theRoot;
    return result;
  }
  
  // Get the rule
  Production& pr = result.theRules[theRoot];
  
  // If the production contains 1 element
  if( pr.theTerm2==0 ) {
    pr.theTerm1 = -pr.theTerm1;
    return result;
  }
  
  // If the production contains 2 elements
  pr.theTerm1 = -pr.theTerm1;
  pr.theTerm2 = -pr.theTerm2;
  swap( pr.theTerm1 , pr.theTerm2 );
  
  return result;
}


//---------------------------------------------------------------------------//
//-------------------------------- operator << ------------------------------//
//---------------------------------------------------------------------------//


ostream& operator << ( ostream& os , const StraightLineProgramWord& CS )
{
  os << "{" << endl
     << "  R = "
     << CS.theRoot << endl
     << "  T = " 
     << CS.theTerminals << endl;
  
  map< int , StraightLineProgramWord::Production >::const_iterator r_it=CS.theRules.begin( );
  for( ; r_it!=CS.theRules.end( ) ; ++r_it ) {
    const StraightLineProgramWord::Production& pr = (*r_it).second;
    os << "  " 
       << (*r_it).first << " -> (" 
       << pr.theTerm1 << ","
       << pr.theTerm2 << ") " << pr.theLength << endl;
  }
  os << "}" << endl;
  
  return os; 
}


//---------------------------------------------------------------------------//
//----------------------------------- getWord -------------------------------//
//---------------------------------------------------------------------------//


Word StraightLineProgramWord::getWord( int N ) const
{
  // A. If the word is empty.
  if( N==0 )
    return Word( );
  
  // B. Get the rule for the current (non-)terminal.
  map< int , Production >::const_iterator r_it = theRules.find( abs(N) );
  
  // C. If N is a non-terminal.
  if( r_it==theRules.end( ) )
    return Word( N );

  // D. get indices of descendants
  int t1 = (*r_it).second.theTerm1;
  int t2 = (*r_it).second.theTerm2;

  // E. Compute the values for descendants
  Word w1 = getWord( t1 );
  Word w2 = getWord( t2 );
  
  return N>0 ? w1*w2 : -w2*-w1;
}


//---------------------------------------------------------------------------//
//----------------------------- getGenerator --------------------------------//
//---------------------------------------------------------------------------//


int StraightLineProgramWord::getGenerator( int n , LongInteger pos ) const
{
  // cout << "n = " << n << ", p = " << pos << endl;
  if( n==0 )
    return 0;
  if( abs(n)<=theTerminals ) {
    if( pos==0 )
      return n;
    else
      return 0;
  }
  
  Production pr = (*theRules.find(abs(n))).second;
  if( pos<0 || pr.theLength<=pos )
    return 0;

  int T1 = pr.theTerm1;
  int T2 = pr.theTerm2;
  if( n<0 ) {
    if( T2!=0 ) {
      swap( T1 , T2 );
      T2 = -T2;
    }
    T1 = -T1;
  }
  LongInteger L1 = length_rule( T1 );
  LongInteger L2 = length_rule( T2 );
  
  if( T2==0 )
    return getGenerator( T1 , pos );

  if( pos<L1 )
    return getGenerator( T1 , pos );
  else
    return getGenerator( T2 , pos-L1 );
}

//---------------------------------------------------------------------------//
//-------------------------- StraightLineProgramWord ------------------------//
//---------------------------------------------------------------------------//


void StraightLineProgramWord::invertProductionPair( int& A , int& B )
{
  if( B!=0 )
    swap( A , B );
  A = -A;
  B = -B;
}
