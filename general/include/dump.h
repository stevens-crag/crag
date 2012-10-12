// Contents:  Implements helper class to dump container
//            objects into a stream           
//
// 
//
// Principal Author: Alexei Miasnikov (2002)
//
// Status: in progress
//
// Revision History:


//    $Id: dump.h,v 1.1 2005/11/11 15:18:35 amiasnik Exp $

#ifndef _DUMP_H_
#define _DUMP_H_


#include <iostream>
#include <vector>




template<class T> class dump
{
public:
  typedef class T::const_iterator const_iterator;
  dump( const T& s ): theCont( s ){}
  
  void printOn( ostream& o) const{
    o << "{ ";
    for (const_iterator I=theCont.begin();I != theCont.end();I++) 
      o << *I << ", "; 
    o << "}" << flush;
  }
  
  inline friend ostream& operator << ( ostream& o, const dump& d ) {
    d.printOn( o );
    return o;
  }

 private:
  const T& theCont; 
};


#endif
