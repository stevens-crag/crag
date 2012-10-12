
// Contents: Implementation of ConfigDile class
//
// Principal Author: Alexei Miasnikov
//
// Status: in progress
//
// Revision History:
//
//


#include <string>
#include "ConfigFile.h"
#include <iostream>
#include <fstream>

ConfigFile::ConfigFile( )
{
}

ConfigFile::ConfigFile( const string& f_name)
{  
  ifstream in(f_name.c_str());

  if (!in){
    cout << "ConfigFile( const string& f_name) : can't open the configuration file." << endl;
    exit(0);
  }

  while( !in.eof() ) {
    
    char s[500],c;
    in.get(s,499,'\n');
    in.get(c);
    if( strlen(s) == 0 || s[0] == '#' ) continue;

    istringstream in(s);
    char var[500];
    char value[500];
    in.get(var,499,':');
    in.get(c);
    in.get(value,499,'\n');

    trim(var);
    trim(value);

    parameters[string(var)] = Value(value);
  }
}


void ConfigFile::setVariable( const char* varName, const Value& value)
{
  parameters[string(varName)] = value;
}

const Value& ConfigFile::getValue( const char* varName)
{
  ParameterType::iterator ret  = parameters.find(string(varName));
  if (ret == parameters.end()){
    cout << "ConfigFile::getValue( ... ) : requested parameter \"" << varName << "\""
	 << "has not been found." << endl;
    exit(0);
  }
  
  return ret->second;
}


void ConfigFile::printOn( ostream& ostr ) const
{
  for (ParameterType::const_iterator I=parameters.begin();I != parameters.end();I++){
    //ostr << "\"" << I->first << "\" : " << (string)I->second << endl;
    ostr << I->first << " : " << (string)I->second << endl;
 }
}



void ConfigFile::rtrim(char* s)
{
  int l = strlen(s);
  for (int i=l-1;i>=0;i--)
    if ((s[i] != '\0') && (s[i] != '\n') && (s[i] != ' ')){
      s[i+1] = '\0';
      break;
    }
}

void ConfigFile::ltrim(char* s)
{
  int l = strlen(s);
  int space_offset = 0;
  for ( space_offset=0;s[space_offset]==' ' && space_offset < l;space_offset++);
  
  if (space_offset>0)
  for (int i=space_offset;i<l;i++){
    s[i-space_offset] = s[i];
    s[i] = '\0';
  }

}

void ConfigFile::trim(char* s)
{
  ltrim(s);
  rtrim(s);
}
