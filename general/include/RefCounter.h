
#ifndef _REF_COUNTER_H_
#define _REF_COUNTER_H_

#ifdef _MSC_VER
#pragma warning(disable:4786)
#endif

class RefCounter {

public:

  /////////////////////////////////////////////////////////////////////////
  //                                                                     //
  // Exported Type Handles                                               //
  //                                                                     //  
  /////////////////////////////////////////////////////////////////////////

  typedef unsigned int refCounterType;
  
  /////////////////////////////////////////////////////////////////////////
  //                                                                     //
  // Constructors                                                        //
  //                                                                     //  
  /////////////////////////////////////////////////////////////////////////

  RefCounter( ) : xrefs(0) { }
  // a new reference counter is initialised
  // with a ref count of 1 (ie 0 extra refs)

  RefCounter( const RefCounter& rc ) : xrefs(0) { }
  // for derived representation classes whose copy constructor is generated
  // by the compiler: a new (copied) representation also has an initial
  // extra ref count of 0

  // Destructor provided by compiler

  /////////////////////////////////////////////////////////////////////////
  //                                                                     //
  // Reference Count Handles                                             //
  //                                                                     //  
  /////////////////////////////////////////////////////////////////////////

  bool lastRef( ) const { return xrefs == 0; }
  
  bool sharedRef( ) const { return xrefs != 0; }
  
  void addRef( ) const { ++((RefCounter*)this)->xrefs; }
  // addRef is logically const
  
  void delRef( ) const { --((RefCounter*)this)->xrefs; }
  // delRef is logically const

  #ifdef DEBUG
  refCounterType nxrefs( ) const { return xrefs; }
  #endif
  
private:
  
  /////////////////////////////////////////////////////////////////////////
  //                                                                     //
  // Data Members                                                        //
  //                                                                     //  
  /////////////////////////////////////////////////////////////////////////

  refCounterType xrefs; // extra references (ie 0 means one ref)

  /////////////////////////////////////////////////////////////////////////
  //                                                                     //
  // Hidden Members                                                      //
  //                                                                     //  
  /////////////////////////////////////////////////////////////////////////

  RefCounter& operator = ( const RefCounter& );
  // disable assignment: normally representations are never assigned

};

#endif

