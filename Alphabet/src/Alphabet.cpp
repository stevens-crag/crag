
// Contents: Implementation of Alphabet related  classes
//
// Principal Author:   Alexei Miasnikov
// Copyright (C) 2005
//
// Status: in progress
//
// Revision History:
//

#include <iostream>
#include <vector>
#include <string>
#include "Word.h"
#include "Parser.h"
#include <stdlib.h>
#include <sstream>
#include "Alphabet.h"
#include "errormsgs.h"
#include <iterator>

void Alphabet::printWord( ostream& out, const Word& w )const {
  int p_base = 0;
  int deg = 0;
  
  if( w.length()==0 ) {
    out << "1";
    return;
  }
  
  int pb;
  for( ConstWordIterator  e_it = w.begin( ) ; e_it!=w.end() ; ++e_it ) {
    if( *e_it==p_base )
      ++deg;
    else {
      pb = p_base;
      if( pb!=0 ) {
	out << getLetter( abs( p_base ) );
	if( deg!=1 || pb<0 )
	  out << "^" << ( pb<0 ? -deg : deg );
	out << " ";
      }
      pb = p_base = *e_it;
      deg = 1;
    }
  }
  
  out << getLetter( abs( p_base ) );
  if( deg!=1 || pb<0 )
    out << "^" << ( pb<0 ? -deg : deg );
  
}

Word Alphabet::readWord( istream& in )const {
  Parser p( cin,(const Alphabet*)this );
  p.parse();
  
  return Word( p.getWord() );
}  


void Alphabet::printVector( ostream& out, const vector<Word>& v )const
{
  out << "{ ";
  for (size_t i=0;i<v.size()-1;i++){
    printWord(out,v[i]); out << ", ";
  }
  printWord(out,v[v.size()-1]); out << " }" << flush;
}
  

vector<Word> Alphabet::readVector( istream& in )const
{
  list<Word> l;

  // read the first symbol which supposed to be '}'
  char ch = ' ';
  while ( ch == ' ' || ch=='\t' || ch == '\n' || ch == '\r' ){
    in >> ch;
  }
  
  if (ch!='{') {
    msgs::error("Alphabet::readVector(): syntax error.");
  }
  while (1){
    Parser p(in,this);
    p.parse();
    Word w(p.getWord());
    l.push_back(w);
    if (p.getWordTerminalSymbol() == '}') break;
    if (p.getWordTerminalSymbol() != ',' && p.getWordTerminalSymbol() != ';'){
      msgs::error("Alphabet::readVector(): syntax error.");
    }
  }
  
  vector<Word> v(l.size());
  copy(l.begin(),l.end(),v.begin());
  return v;
}



/**********************************************************
 *
 *  FINITE ALPHABET
 *
 ***********************************************************/


int FiniteAlphabet::getNum( const string& letter ) const {
  for (size_t i=0;i<theLetters.size();i++)
    if ( theLetters[i] == letter ) return i+1;
  
  return 0;
}


string FiniteAlphabet::getLetter( int index ) const {
  if (index == 0 || static_cast<unsigned int>(abs(index)) > theLetters.size() ){
    cout << "Index " << index << " is out of bounds" << endl;
    exit(0);
  }
  return theLetters[index-1];
}

const vector<string>& FiniteAlphabet::getLetters( ) const {
  return theLetters;
}


/**********************************************************
 *
 *  INFINITE ALPHABET
 *
 ***********************************************************/

int InfiniteAlphabet::getNum( const string& letter ) const {
  // match the letter prefix first
  if ( thePrefix.size() >= letter.size() ) return 0;
  for (size_t i=0;i<thePrefix.size();i++)
    if ( thePrefix[i] != letter[i] ) return 0;
  
  // extract the index
  int index = atoi((letter.substr(thePrefix.size(),letter.size()-thePrefix.size())).c_str());
  
  return index;
}
  

string InfiniteAlphabet::getLetter( int index ) const {
  if (abs(index) == 0 ){
    cout << "Index " << index << " is out of bounds" << endl;
    exit(0);
  }
  stringstream ss;
  ss << thePrefix << index  << flush;
  return ss.str();
}

InfiniteAlphabet InfiniteAlphabet::defaultAlphabet = InfiniteAlphabet();


