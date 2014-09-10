// Contents: Wrapper for ranlib functions
//           
//
// 
//
// Principal Author: Alexei Miasnikov (2002)
//
// Status: in progress
//

//
//
//   $Id: RanlibCPP.h,v 1.5 2007/09/27 16:44:03 amiasnik Exp $
//
//

#ifndef _RANLIB_H_
#define _RANLIB_H_


#include "ranlib.h"
#include "cdflib.h"
#include <time.h>
#ifdef _MSC_VER
//for getpid
#include <process.h>
#else
#include <unistd.h>
#endif
#include <stdlib.h>
#include <iostream>

//! A wrapper class for  RANLIB Library
/*! This is a wrapper class for functions implemented in <br>
RANLIB  Library of Fortran Routines for Random Number Generation, Compiled and Written by: <br>
    Barry W. Brown and   James Lovato <br>
Department of Biomathematics, Box 237,
The University of Texas, M.D. Anderson Cancer Center,
1515 Holcombe Boulevard,
Houston, TX      77030
*/
class RandLibURG
{
 public:
  //! Constructor.
  /*! Set the seeds of random generator to \c i1 and \c u2 */
  RandLibURG(long i1,long i2) {
    setall(i1,i2);
  }

  //! Default constructor. The seeds of the random number generator a set using current time values*/
  RandLibURG() {
    setall(time(NULL),time(NULL));
    long low = irand(0,2147483560);
    long high = irand(0,2147483560);
    setall(low,high);
  }

  //! Set the seeds of random generator to \c i1 and \c u2 
  void reset(long i1,long i2) {
    setall(i1,i2);
  }

  //! Set the seeds of the random number generator  using current time values  
  void reset(){
    setall(time(NULL),time(NULL));
    long low = irand(0,2147483560);
    long high = irand(0,2147483560);
    setall(low,high); 
  }

  //! Set the seeds of the random number generator  using current time and process id values   
  void setSeedPID(){
#ifdef _MSC_VER
    long pid = _getpid();
#else
    long pid = getpid();
#endif
    long low  = labs( (pid + (pid % 131) * time(NULL)) % 2147483560); // pid +  ... is to make shure that 
    long high = labs( (pid + (pid % 17)  * time(NULL)) % 2147483560); // that seed will not reduce to 0
    setall(low,high); 
    low = irand(0,2147483560); 
    high = irand(0,2147483560); 
    setall(low,high);  
  } 

  //! Get the seeds of the random number generator
  void getseed( long& s1, long& s2 ){
    long iseed1,iseed2;
    getsd(&iseed1,&iseed2);
    s1 = iseed1; s2 = iseed2;
  }

  //! Returns a float point pseudo random number in [0,1] 
  float rand() {
    return ranf(); 
  }

  
  //! Returns a float point pseudo random number in [low,high] 
  float rand(float low,float high) {
    return genunf( low, high);
  }

  
  //! Returns an integer pseudo random number in [low,high] 
  long irand(long low,long high){
   return  ignuin( low, high);
  }
};


//! Static wrapper class for RANDLIB Library
class RandLib
{
 public:
  //! The static instance of the wrapper class. 
  /*! This is done to make sure that the random number generator is seeded and  only once. */
  static RandLibURG ur;

  
  //! Returns a p-value for the chi distribution. 
  static double chiPValue( int idf, double alpha ){
    int which = 2;
    double p = 1.0 - alpha;
    double q = alpha;
    double x;
    double df = idf;
    int status;
    double bound;
   
    cdfchi(&which,&p,&q,&x,&df,&status,&bound);

    return x;
  }

};
#endif
