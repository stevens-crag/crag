
#ifndef _PARSER_H_
#define _PARSER_H_

#include <iostream>
#include <list>

class Alphabet;
class FiniteAlphabet;

using namespace std;

enum AInputType { NONE, SET, FREE_PRES, REL_PRES }; 

 
class Parser
{
 public:
  Parser(istream& in, const Alphabet* a);
  ~Parser();
  void parse();

  char getWordTerminalSymbol( ) const;
 
  const list<int>& getWord()const;
 private:
  class WordFlexLexer* localFlexLexer;
  list<int> theWord;
  char theTS;
};

class AParser
{
 public:
  AParser(istream& in);
  ~AParser();
  void parse();
  AInputType getType()const;
  FiniteAlphabet getAlphabet() const;
 private:
  class AlphabetFlexLexer* localFlexLexer;
  AInputType theType;
  FiniteAlphabet* theAlphabet;
};

#endif
