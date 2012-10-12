/*
*   Author: Alexei Miasnikov (2006)
*/
%{


#include  "AlphabetBisonGrammar.h"
#include  "Parser.h"
using namespace std;


#include <iostream>
#include <sstream>
#include <string>




/*
//extern int yyinput();


//strstream*  inputString = new strstream();
//strstream inputString;
//ostream& out = cout;
//ostream& out = *inputString;
//InputBuffer inputString(100000000);
//InputBuffer& out  = inputString;

*/

extern string a_p_buf_string;
extern bool stop_parser( AInputType t );

%}
%option prefix="aa"
%option noyywrap

SP [ ]+
GEN    [a-zA-Z_][a-zA-Z_0-9]*
ANYTHING  .

%%

"{"  return '{';
"}"  {stop_parser( SET ); yyterminate(); }
"<"  return '<';
">"  {stop_parser( FREE_PRES );yyterminate();}
"|"  {stop_parser( REL_PRES );yyterminate();}
";"  {stop_parser( REL_PRES );yyterminate();}
","  return ',';
{GEN}        {  a_p_buf_string = string(yytext); return ALETTER; }


[\r\t\n ]+    { /*cout << yytext << flush ;*/}

%%


/*

int main(int argc,char** argv)
{
	yylex();
	return 0;
}
*/



