// Copyright (C) 2005 Alexander Ushakov
// Contents: Implementation of class FPGroup
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//


#include "RanlibCPP.h"
#include "FPGroup.h"
#include "BalancedTree.h"
#include "errormsgs.h"

#include <set>


//---------------------------------------------------------------------------//
//--------------------------------- FPGroup ---------------------------------//
//---------------------------------------------------------------------------//


FPGroup::FPGroup( int numOfGen ) :
  numOfGenerators( numOfGen ),
  theAlphabet( numOfGen ),
  theRelators( )
{
  
}

//=========================================================

	
FPGroup::FPGroup(  const FiniteAlphabet& a  ) :
  numOfGenerators( a.size() ),
  theAlphabet( a ),
  theRelators( )
{
  
}


//=========================================================


FPGroup::FPGroup(  const FiniteAlphabet& a, const vector< Word >& relators ) :
  numOfGenerators( a.size() ),
  theAlphabet( a ),
  theRelators( relators )
{
  
}


//=========================================================


FPGroup::FPGroup( int numOfGen , const vector< Word >& relators ) :
  numOfGenerators( numOfGen ),
  theAlphabet( numOfGen ),
  theRelators( relators )
{
  
}


//=========================================================


vector< string > FPGroup::initializeGenNames( int num )
{
  vector< string > result( num );
  
  char str[10];
  for( int t=0 ; t<num ; ++t ) {
    strstream(str,10) << t << ends;
    // itoa( t , str , 10 );
    result[t] = "x";
    result[t].append( str );
  }
  return result;
}


//=========================================================

FPGroup FPGroup::triangulatePresentation( ) const
{
  int new_gen_num = numOfGenerators;
  char str[10];
  
  // rename old generators to avoid collisions
  vector< string > new_gen_names = initializeGenNames( numOfGenerators );
  vector< Word > new_relators;
  
  int new_gens_num = 0;
  for( vector< Word >::const_iterator r_it=theRelators.begin( ) ; r_it!=theRelators.end( ) ; ++r_it ) {
    
    const Word& rel = *r_it;
    int len = rel.length( );
    if( len>3 ) {

      // initial setup (first 3 letters)
      Word::const_iterator r_it = rel.begin( );
      int g1 = *r_it++;
      int g2 = *r_it++;
      int g3 = *r_it++;
      strstream(str,10) << (++new_gens_num) << ends;
      new_gen_names.push_back( string( "x" ).append( str ) );
      new_relators.push_back( Word( g1 )*Word( g2 )*Word( -new_gens_num ) );
      
      // iterations
      for( int i=0 ; i<len-4 ; ++i, ++r_it ) {
        strstream(str,10) << (++new_gens_num) << ends;
        new_gen_names.push_back( string( "x" ).append( str ) );
        int new_gen = new_gen_names.size( );
        new_relators.push_back( Word( new_gens_num-1 ) * Word(g3) * Generator( -new_gens_num ) );
        g3 = *r_it;
      }
      
      // the last relator
      new_relators.push_back( Word( g3 ) * Generator( *r_it ) * Word( -new_gens_num ) );
    } else
      new_relators.push_back( rel );
  }

  return FPGroup( FiniteAlphabet( new_gen_names) , new_relators );
}


//=========================================================


Word FPGroup::randomEqWord_Baltimore( const Word& w , int length , float conj_param ) const
{
  BalancedTree< int > BT( w.begin( ), w.end( ) );
  
  while( BT.size( )<length ) {
    
    int pos = RandLib::ur.irand( 0 , BT.size( ) );

    int conj_length = 0;
    for( ; RandLib::ur.rand()<conj_param ; ++conj_length );
    Word conjugator = Word::randomWord( numOfGenerators , conj_length );
    
    Word rel( theRelators[RandLib::ur.irand( 0 , theRelators.size( )-1 )] );
    rel = ( RandLib::ur.irand( 0 , 1 )==0 ? rel : rel.inverse( ) );
    rel.cyclicallyPermute( RandLib::ur.irand( 0 , rel.length( )-1 ) );

    Word to_insert = -conjugator * rel * conjugator;
    BT.insert( pos , to_insert.begin() , to_insert.end( ) );
  }
  
  return Word( BT.getList( ) );
}


//=========================================================


Word FPGroup::randomIdentity_Baltimore( int length , float conj_param ) const
{
  return randomEqWord_Baltimore( Word( ) , length , conj_param );
}


//=========================================================


Word FPGroup::randomIdentity_Stack( int length ) const
{
  if( theRelators.size( )==0 )
    return Word();

  // 1. Construct a set of rules
  int ngens = numberOfGenerators( );
  int nrels = theRelators.size( );
  vector< pair< Word , Word > > rules;
  for( int r=0 ; r<nrels ; ++r ) {
    Word rel = theRelators[r];
    for( int l=0 ; l<rel.length() ; ++l ) {
      rel = rel.cyclicallyPermute( 1 );
      for( int p=0 ; p<=rel.length() ; ++p ) {
        Word left = rel.initialSegment( p );
        Word right = rel.terminalSegment( p );
        rules.push_back( pair< Word , Word >(  left ,  right ) );
        rules.push_back( pair< Word , Word >( -left , -right ) );
      }
    }
  }
  for( int i=0 ; i<ngens ; ++i ) {
    rules.push_back( pair< Word , Word >( Generator( i+1 ) , Generator( -i-1 ) ) );
    rules.push_back( pair< Word , Word >( Generator(-i-1 ) , Generator(  i+1 ) ) );
  }

  int rule_num = rules.size( );

  // 2. Construct a word using the rules
  Word result;
  while( result.length() < length ) {
    while( result.length() < length ) {
      int rn = RandLib::ur.irand( 0 , rule_num-1 );
      result = rules[rn].first * result * rules[rn].second;
    }
  }

  return result;
}


//=========================================================


Word FPGroup::randomIdentity_Classic( int length , float conj_param ) const
{
  if( theRelators.size( )==0 )
    return Word();

  int ngens = numberOfGenerators( );
  int nrels = theRelators.size( );
  Word result;

  while( result.length()<length ) {
    while( result.length()<length ) {

      // generate conjugator
      int conj_length = 0;
      for( ; RandLib::ur.rand()<conj_param ; ++conj_length );
      Word conjugator = Word::randomWord( numOfGenerators , conj_length );

      Word rel( theRelators[RandLib::ur.irand( 0 , theRelators.size( )-1 )] );
      rel = ( RandLib::ur.irand( 0 , 1 )==0 ? rel : rel.inverse( ) );
      rel.cyclicallyPermute( RandLib::ur.irand( 0 , rel.length( )-1 ) );
      result *= -conjugator * rel * conjugator;
    }
  }

  return result;
}

//=========================================================


ostream& operator << ( ostream& os , const FPGroup& group )
{
  os << "< ";
  
  int i;
  const vector< string >& names = group.getGeneratorsNames();
  for( i=0 ; i<names.size( ) ; ++i ) {
    if( i )
      os << " , ";
    os << names[i];
  }
  
  const vector< Word >& relators = group.relators();
  if( relators.size( )>0 )
    os << " ; ";
  for( i=0 ; i<relators.size( ) ; ++i ) {
    if( i )
      os << " , ";
    // group.writeWord( os , relators[i] );
    group.getAlphabet().printWord( os,relators[i] );
  }
  
  os << " >";
  return os;
}


//=========================================================

istream& operator >> ( istream& in ,  FPGroup& group )
{
  
  AParser ap( in );
  ap.parse();
  if ( ap.getType() != REL_PRES) msgs::error("readFPPresentation() : syntax error");
  group.theAlphabet = ap.getAlphabet();
  
  
  list<Word> l;
  
  while (1){
    Parser p(in,&group.theAlphabet);
    p.parse();
    Word w(p.getWord());
    l.push_back(w);
    if (p.getWordTerminalSymbol() == '>') break;
    if (p.getWordTerminalSymbol() != ',' ){ //&& p.getWordTerminalSymbol() != ';'){
      msgs::error("readFPPresentation(): syntax error.");
    }
  }
  
  group.theRelators = vector<Word>(l.size());
  copy(l.begin(),l.end(),group.theRelators.begin());
  
  //  cout << " Al : " << a << endl << "Rels : ";
  //  copy(v.begin(),v.end(),ostream_iterator<Word>(cout,", "));
  
  
}

//=========================================================

