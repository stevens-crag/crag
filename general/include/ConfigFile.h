// Contents: Definition of ConfigFile class
//
// Principal Author:   Alexei Miasnikov
// Copyright (C) 2005
//
// Status: in progress
//
// Revision History:
//


#ifndef _WSET_CONFIG_H_
#define _WSET_CONFIG_H_

#include <string.h>
#include <sstream>
#include <iostream>
#include <map>
#include <stdlib.h>

using namespace std;

//! Implements a comparison operator on two strings.
struct ltstr
{
  bool operator()(const string& s1, const string& s2) const
  {
    return strcmp(s1.c_str(), s2.c_str()) < 0;
  }
};


//! Implements data types of parameters obtained from a configuration file.  
/*!
  Allows universal handling of parameters with different data types.
  Usually an instance of this class is returned by \link ConfigFile::getValue() ConfigFile::getValue()\endlink function.
*/
class Value
{
public:
  //! Default constructor
  Value() {}
  //! Constructor 
  /*!
    \param v - a string  value.
  */
  Value( const char* v) : theValue(v) { }	
  //! Constructor 
  /*!
    \param n - an integer value.
  */    
  Value( int n ) {
    stringstream s;
    s << n << '\0' << flush;
    theValue = string(s.str());
  }
  
  //! Converts value into an integer.
  /*!
    \return an integer value. 0 if an error.
  */
  operator int( ) const { return atoi(theValue.c_str()); }
  //! Converts value into a double precision number.
  /*!
    \return a double  value. 0 if an error.
  */
  operator double( ) const { return atof(theValue.c_str()); }
  //! Converts value into a string.
  /*!
    \return a string  value.
  */  
  operator string( ) const { return theValue; }
  //! Input operator.
  friend istream& operator >> (istream& in, Value& v){
    in >> v.theValue;
    return in;
  }
  
private:
  string theValue;
};


//! Implements a mapping from parameter's name into its value.
typedef map<string, Value, ltstr> ParameterType;


//------------------------------ ConfigFile ----------------------------------//


//! Implements a mechanism for passing parameters from a configuration file.
/*!
  This class generates a list of parameters and their values from a given
  configuration file. The format of the configuration file is follows:
  
  ParameterName1 : ParameterValue1  <br>
  ParameterName2 : ParameterValue2  <br>
  . . .            <br>
  ParameterNameN : ParameterValueN <br>
  
  Parameter names are identifiers containing letters and numbers. Parameter values
  can be an integer or float point numbers and strings.
  
  It is possible to add comments. Everything after a symbol \# and until the end of a line
  is considered to be a comment.
  
  \bug On some platforms scanning may hang if an empty string is contained in the
  input file.
  
*/
class ConfigFile
{
  
public:
  
  /////////////////////////////////////////////////////////////////////////
  //                                                                     //
  // Constructors:                                                       //
  //                                                                     //
  /////////////////////////////////////////////////////////////////////////
  
  //! Default constructor. 
  ConfigFile();
  //! Constructor. Collects parameters and their values from a configuration file.
  /*!
    Will scan the file whose name is given as a parameter and create a list
    of parameters together with the corresponding values.
    \param f_name - name of the configuration file.
  */
  ConfigFile(const string& f_name);
  
  // copy constructor supplied by compiler.
  // destructor supplied by compiler.
  
  
  /////////////////////////////////////////////////////////////////////////
  //                                                                     //
  // Accessors/Modifiers:                                                          //
  //                                                                     //
  /////////////////////////////////////////////////////////////////////////
  
  //! Set a value of a parameter.
  /*!
    Sets a value of a parameter with the given name. Will create the parameter
    if not in the list.
    
    \param p_name - the name of the parameter.
    \param p_value - the value of the parameter.
  */
  void setVariable( const char* p_name, const Value& p_value );
  // used by readFrom() to set variables
  
  
  //! Get a value of a parameter.
  /*!
    Returns a value of the parameter with the name p_name. <br>
    Note, since class \link Value Value \endlink is equipped with the 
    operators converting the parameter value into standard types,  it is not necessary to create instances
    of the class Value. Typical use of this function: <br>
    \code int v = configFile.getValue("THENAME"); \endcode
    or with explicit conversion:
    \code string v = string(configFile.getValue("THENAME")); \endcode
  */
  const Value& getValue(const char* p_name);
  // returns value of the parameter

  /////////////////////////////////////////////////////////////////////////
  //                                                                     //
  // OI:                                                                 //
  //                                                                     //
  /////////////////////////////////////////////////////////////////////////
  
  // assignment operator supplied by compiler

  //! Output operator
  /*!
    Prints a list of parameters together with their values.
  */
  friend ostream& operator << ( ostream& ostr, const ConfigFile& C )
  {
    C.printOn(ostr);
    return ostr;
  }
    

private:

  /////////////////////////////////////////////////////////////////////////
  //                                                                     //
  // Private functions:                                                  //
  //                                                                     //
  /////////////////////////////////////////////////////////////////////////

  void readFrom( istream& istr );

  void printOn( ostream& ostr ) const;

  void rtrim(char* ch);
  void ltrim(char* ch);
  void trim(char* ch);

  /////////////////////////////////////////////////////////////////////////
  //                                                                     //
  // Data members:                                                       //
  //                                                                     //
  /////////////////////////////////////////////////////////////////////////

  ParameterType parameters;

  
};
#endif
