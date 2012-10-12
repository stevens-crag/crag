
// Contents: Inlined global functions error, warn.
//
// Principal Author: Copyright (2005) Alexei Miasnikov
//
// Status: Useable.

#ifndef _ERROR_H_
#define _ERROR_H_

#include <iostream>

using namespace std;

//! Errors output handling
namespace msgs {

  //! Print a message  to the standard error output.
  /*!
    This function will go  into infinite loop (freeze) so one can get a stack backtrace when debugging.    
    \param msg - the error message.
   */ 
  inline void error(const char* msg) {
    
    cerr << endl << "Fatal error: " << msg << endl;
    
    while ( 1 );
  }
  
  //! Print a message  to the standard error output.
  /*!
    This function will go  into infinite loop (freeze) so one can get a stack backtrace when debugging.    
    \param msg - the error message.
   */ 
  inline void error(const string& msg) {
    
    cerr << endl << "Fatal error: " << msg << endl;
    
    while ( 1 );
  }

  //! Conditional output of  an error  message  to the standard error output.
  /*!
    Prints an error message if \c exp is not \c true. 
    This function will go  into infinite loop (freeze) so one can get a stack backtrace when debugging.    
    \param exp - boolean expression.
    \param msg - the error message.
   */ 
  inline bool error(bool exp, const char* msg) {
    
    if (!exp){
      cerr << endl << "Fatal error: " << msg << endl;
      
      while ( 1 );
      return true;
    } else
    return false;
  }
  
  
  //! Print a message  to the standard error output.
  /*!
    \param msg - the error message.
   */ 
  inline void warn(const char* msg) {
    cerr << endl << "Warning: " << msg << endl;
  }
  
}

#endif
