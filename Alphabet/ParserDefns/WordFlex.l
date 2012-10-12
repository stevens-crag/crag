/*
*   Author: Alexei Miasnikov (2005)
*/
%{


#include  "WordBisonGrammar.h"
using namespace std;


#include <iostream>
#include <sstream>
#include <string>
#include <stack>



/*
//extern int yyinput();


//strstream*  inputString = new strstream();
//strstream inputString;
//ostream& out = cout;
//ostream& out = *inputString;
//InputBuffer inputString(100000000);
//InputBuffer& out  = inputString;

*/

stack<int> bracketStack;
int powerFlag = 0;

extern string p_buf_string;
extern void word_terminate_symbol( char s );

%}
%option noyywrap

INTEGER [1-9][0-9]*
SP [ ]+
GEN    [a-zA-Z_][a-zA-Z_0-9]*
ANYTHING  .

%%

"["  { bracketStack.push(1); return '[';}
"]"  { bracketStack.pop(); return  ']'; }
"^"  { powerFlag = 1; return '^'; }
"("  return '(';
")"  return ')';
"-"  { if (powerFlag)
         return '-';
       else {
         word_terminate_symbol( '-' );
         yyterminate();
       }
      }
","  { if (bracketStack.empty() ){
        word_terminate_symbol( ',' );
        yyterminate();
     } else
        return ',';
     }
";"  {  word_terminate_symbol( ';' );  yyterminate(); } 
"}"  {  word_terminate_symbol( '}' );  yyterminate(); } 
">"  {  word_terminate_symbol( '>' );  yyterminate(); } 


{INTEGER} {  powerFlag = 0;  p_buf_string = string(yytext); return POWER; }
{GEN}        {   p_buf_string = string(yytext); return GENERATOR; }


[\r\t\n ]+    { /*cout << yytext << flush ;*/}

%%


/*

int main(int argc,char** argv)
{
	yylex();
	return 0;
}
*/



