/*
 *   
 *
 * Contents: Bison grammars for alphabet parser
 *         
 *
 *
 *
 *  Author: Alexei Miasnikov (2006)
 *
 * Status: in progress
 *
*/


%{


#include <iostream>
#include <sstream>
#include <list>
#include <string>
#include <algorithm>
#include <iterator>
#include "Alphabet.h"
#include <vector>

#ifndef			yyFlexLexer 
#define 		yyFlexLexer aaFlexLexer 
#include		 <FlexLexer.h> 
#elif 	 	 	yyFlexLexer != aaFlexLexer 
#undef			yyFlexLexer 
#define			yyFlexLexer wwFlexLexer 
#include	 	 <FlexLexer.h> 
#endif 
 

#include "Parser.h"

  
  //  typedef list<int> ParseWord;
  
  class AlphabetFlexLexer{
  public:

    
    AlphabetFlexLexer( istream* in , ostream* out  ) {
      iType = NONE;
      lexer = new aaFlexLexer( in, out );   
      lexer->yyrestart(in);
    }
    
    ~AlphabetFlexLexer() { delete lexer; }
    
    void addGenerator( const string& g) { theAlphabet.addGenerator( g ); }
    const FiniteAlphabet& getAlphabet() { return theAlphabet; }
    
    FlexLexer& getLexer() { return *lexer; }

    AInputType getType()const { return iType; }
    void setType( AInputType t) {  iType = t; }
  private:
    FlexLexer* lexer;  
    FiniteAlphabet theAlphabet;
    AInputType iType;
  };
  
  AlphabetFlexLexer* alphabetFlexLexer;
  
  inline int aalex() { return alphabetFlexLexer->getLexer().yylex(); }
  

   void stop_parser(AInputType t) { alphabetFlexLexer->setType(t); }
  
 //  list<string> theParseWord;

  using namespace std;

  extern int aalex();
  void aaerror(char*);
  //  extern void yyterminate( );

  string a_p_buf_string;


  //#define Info( S ) cout << "[" << S << "]:" << flush;
  #define Info( S );

  #define Flush_STRING_POINTER( S ) { cout << *(S) << flush; delete S; }


%}

/* Declarations */

%union{
  void* lst;
  void* str;
  int   num;
}

%token ALETTER

%type <str> ALETTER_STR
/*%type <num> GENERATOR_STR*/
/*%type <lst> generator*/
%%

/***************************************************************************
*
* general input grammars
*
****************************************************************************/



input:  '{' generator_complete_list { if (alphabetFlexLexer->getType() !=  SET ) aaerror("Letter set syntax error"); }
| '<' generator_complete_list { 
  if (alphabetFlexLexer->getType() != FREE_PRES && alphabetFlexLexer->getType() != REL_PRES   ) 
    aaerror("Presentation syntax error"); 
}
;


generator_complete_list: generator_list ALETTER_STR { alphabetFlexLexer->addGenerator( *(string*)$2 ); delete (string*)$2; }

generator_list: /*empty*/
| generator_list ALETTER_STR ',' { alphabetFlexLexer->addGenerator( *(string*)$2 ); delete (string*)$2; }
;


ALETTER_STR: ALETTER {  $$ = new string( a_p_buf_string ); }
;

%%

void aaerror( char* s ){
  cout << " Error. Parser stoped: { " << s << " } " <<  flush;
  exit(0);
  //  yyerrok;
  //  yyclearin;
}

AParser::AParser(istream& in)
  : theAlphabet( NULL )
{
  localFlexLexer = new AlphabetFlexLexer(&in,&cout);
}
AParser::~AParser()
{
  if (theAlphabet) delete theAlphabet;
  delete localFlexLexer;
}

void AParser::parse() { 
  AlphabetFlexLexer* savFlexLexer = alphabetFlexLexer;
  alphabetFlexLexer = localFlexLexer;  
  aaparse(); 
  theAlphabet = new FiniteAlphabet(alphabetFlexLexer->getAlphabet());
  theType =  alphabetFlexLexer->getType();
  alphabetFlexLexer = savFlexLexer;
}

FiniteAlphabet AParser::getAlphabet()const { return *theAlphabet; }

AInputType AParser::getType()const { return theType;  }


/*
int main()
{
  Parser p( cin );
  //  p.parse();

  return 0;
}
*/
