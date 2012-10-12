/*
 *   
 *
 * Contents: Bison grammars for word parsers
 *         
 *
 *
 *
 *  Author: Alexei Miasnikov (2005)
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

#include <FlexLexer.h> 
 

#include "Parser.h"


#define formulaLatex2HTML( X, A ) ()
#define TOGGLE( A ) ( A = (A + 1 )% 2 )


  //#include "FlexLexer.h"

  typedef list<int> ParseWord;
  
  
  class WordFlexLexer{
  public:
    WordFlexLexer( istream* in , ostream* out, const Alphabet* a  ): theAlphabet( a ) {
      lexer = new yyFlexLexer( in, out );   
      lexer->yyrestart(in);
      theTerminalSymbol = '0';
    }
    
    ~WordFlexLexer() { delete lexer; }
    
    
    const Alphabet* getAlphabet() { return theAlphabet; }

    void setWordTerminalSymbol( char s ) { theTerminalSymbol = s; }
    char getWordTerminalSymbol(  ) const { return theTerminalSymbol; } 

    FlexLexer& getLexer() { return *lexer; }
    const ParseWord& getWord() { return theParsedWord; } 
    void setWord( const ParseWord& w)
      {
	theParsedWord = w;
      }
  private:
    ParseWord theParsedWord;
    FlexLexer* lexer;  
    const Alphabet* theAlphabet;
    char theTerminalSymbol;
  };
  
  WordFlexLexer* cragFlexLexer = NULL;
  
  inline int yylex() { return cragFlexLexer->getLexer().yylex(); }
  

  void word_terminate_symbol( char s ) { cragFlexLexer->setWordTerminalSymbol(s); }


  //  list<string> theParseWord;

  using namespace std;

  extern char* yytext;
  //  extern int yylex();
  void yyerror(char*);
  //  extern void yyterminate( );

  string p_buf_string;


  //#define Info( S ) cout << "[" << S << "]:" << flush;
  #define Info( S );

  #define Flush_STRING_POINTER( S ) { cout << *(S) << flush; delete S; }


  ParseWord inv( const ParseWord& w ){
    ParseWord invw;
    for ( ParseWord::const_reverse_iterator I=w.rbegin();I!=w.rend();I++ )
      invw.push_back(-(*I));

    return invw;
  }


%}

/* Declarations */

%union{
  void* lst;
  void* str;
  int   num;
}

%token GENERATOR POWER

/*%type <str> GENERATOR*/
%type <num> POWER_NUM
%type <num> GENERATOR_STR
%type <lst> word
%type <lst> generator
%%

/***************************************************************************
*
* general input grammars
*
****************************************************************************/

input:       /* empty */ { Info("Empty input"); }
| input word {  
  cragFlexLexer->setWord( *((ParseWord*)$2) );
}
| error   {  yyerrok; }
;

/***************************************************************************
*
* word grammar
*
****************************************************************************/

word: generator {  Info("GEN"); 
 ParseWord ls = *((ParseWord*)$1);
 //cout << "Param: " << ls.size() << endl;
 
 $$ = new ParseWord();
 ((ParseWord*)$$)->splice( ((ParseWord*)$$)->end(), *(ParseWord*)$1);
 
 delete (ParseWord*)$1;
}

|     word word { Info("GEN WORD"); 
 
 $$ = new ParseWord(*(ParseWord*)$1);
 ((ParseWord*)$$)->splice( ((ParseWord*)$$)->end(), *(ParseWord*)$2); 
 
 // cout <<  ((ParseWord*)$$)->size() << endl;

 delete (ParseWord*)$1;
 delete (ParseWord*)$2;

}
|     '(' word ')' '^' POWER_NUM {   Info( "W_POWER"); 
  int p = $5;
  int ap = abs(p);
  ParseWord w;
  if ( p > 0 )
    w = *((ParseWord*)$2);
  else
    w = inv(*((ParseWord*)$2));

  $$ = new ParseWord();
  for (int i=0;i<ap-1;i++)
    ((ParseWord*)$$)->insert( ((ParseWord*)$$)->end(), w.begin(),w.end());
  
  ((ParseWord*)$$)->splice( ((ParseWord*)$$)->end(),w);
  
 
  //cout <<  ((ParseWord*)$$)->size() << endl;
  
  delete (ParseWord*)$2;
    
   
}
|     '[' word ',' word ']' '^' POWER_NUM { Info(" COMM_POWER "); 
	ParseWord comm(*(ParseWord*)$2); 
	comm.insert( comm.end(), ((ParseWord*)$4)->begin(),((ParseWord*)$4)->end());
	ParseWord iw2 =  inv(*(ParseWord*)$2);
	comm.splice( comm.end(),iw2);
	ParseWord iw4 =  inv(*(ParseWord*)$4);
	comm.splice( comm.end(), iw4);	
  	//cout << "End comm_power" << endl;
	

	int p = $7;
	int ap = abs($7);
	if ( p < 0)
	  comm = inv(comm);
	
	$$ = new ParseWord();
	for (int i=0;i<ap-1;i++)
	  ((ParseWord*)$$)->insert( ((ParseWord*)$$)->end(), comm.begin(),comm.end());
	
	((ParseWord*)$$)->splice( ((ParseWord*)$$)->end(),comm);	
	
	delete (ParseWord*)$2;
	delete (ParseWord*)$4;		
}
|     '[' word ',' word ']'  {  Info(" COMM "); // [a,b] = a b A B
	$$ = new ParseWord(*(ParseWord*)$2);
	((ParseWord*)$$)->insert( ((ParseWord*)$$)->end(), ((ParseWord*)$4)->begin(),((ParseWord*)$4)->end());
	ParseWord iw2 =  inv(*(ParseWord*)$2);
	((ParseWord*)$$)->splice( ((ParseWord*)$$)->end(),iw2);
	ParseWord iw4 =  inv(*(ParseWord*)$4);
	((ParseWord*)$$)->splice( ((ParseWord*)$$)->end(), iw4);
	
 	//cout <<  ((ParseWord*)$$)->size() << endl;

	delete (ParseWord*)$2;
	delete (ParseWord*)$4;	
}
;


generator: GENERATOR_STR { 
  $$ = new ParseWord();
  ((ParseWord*)$$)->push_back($1);
  //cout << ((ParseWord*)$$)->size() << endl;
}
|   GENERATOR_STR '^' POWER_NUM {
  int p = $3;
  int ap = abs(p);
   $$ = new ParseWord();
   for (int i=0;i<ap;i++) {
     if ( p > 0)
       ((ParseWord*)$$)->push_back($1);
     else
       ((ParseWord*)$$)->push_back(-($1));
   }

   //cout << ((ParseWord*)$$)->size() << endl;
}
;

GENERATOR_STR: GENERATOR {  
  int il = cragFlexLexer->getAlphabet()->getNum( p_buf_string );
  if (!il) {
    cout << "Cannot recognize the letter " <<  p_buf_string << endl;
    exit(0);
  }
  $$ = il;
}

POWER_NUM: POWER  { $$ = atoi(p_buf_string.c_str()); }
|      '-' POWER  { $$ = -atoi(p_buf_string.c_str()); }
;

%%


void yyerror( char* s ){
  cout << " Error. Parser stoped: { " << s << " } " <<  flush;
  exit(0);
  //  yyerrok;
  //  yyclearin;
}

Parser::Parser(istream& in, const Alphabet* a)
{
  localFlexLexer = new WordFlexLexer(&in,&cout,a);
}
Parser::~Parser()
{
  delete localFlexLexer;
}

const list<int>& Parser::getWord()const { return theWord; }

char  Parser::getWordTerminalSymbol(  ) const { return theTS;  } 


void Parser::parse() { 
  WordFlexLexer* savFlexLexer = cragFlexLexer;
  cragFlexLexer = localFlexLexer;  
  yyparse(); 
  theWord =  cragFlexLexer->getWord();
  theTS =  cragFlexLexer->getWordTerminalSymbol();
  cragFlexLexer = savFlexLexer;
}

/*
int main()
{
  Parser p( cin );
  //  p.parse();

  return 0;
}
*/
