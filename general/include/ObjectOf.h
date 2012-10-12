
#ifndef _OBJECT_OF_H_
#define _OBJECT_OF_H_


#include <iostream>
using namespace std;

#include "RefCounter.h"


template<class Rep> class ObjectOf
{
	
public:
	
	///////////////////////////////////////////////////////////
	//                                                       //
	// Constructors:                                         //
	//                                                       //
	///////////////////////////////////////////////////////////
	
	ObjectOf( const ObjectOf& o ) { theRep = o.theRep; theRep->addRef(); }
	
	~ObjectOf() { if (theRep->lastRef()) delete theRep; else theRep->delRef(); }
	
	///////////////////////////////////////////////////////////
	//                                                       //
	// Standard Operators:                                   //
	//                                                       //
	///////////////////////////////////////////////////////////
	
#ifndef SWIG
	ObjectOf& operator = ( const ObjectOf& o )
	{
		o.theRep->addRef();
		if ( theRep->lastRef() ) delete theRep; else theRep->delRef();
		theRep = o.theRep;
		return *this;
	}
#endif

	///////////////////////////////////////////////////////////
	//                                                       //
	// Debugging Stuff:                                      //
	//                                                       //
	///////////////////////////////////////////////////////////
	
#ifdef DEBUG
	void debugPrint( ostream& ostr ) {
		ostr << "this : " << this << "; theRep : " << theRep << "; xrefs : "
			<< look()->nxrefs();
	}
	
	//friend int main( );
#endif
	
	
protected:
	
	///////////////////////////////////////////////////////////
	//                                                       //
	// Representation Access:                                //
	//                                                       //
	///////////////////////////////////////////////////////////
	
	const Rep* look( ) const { return theRep; }
	// For safe read-only access.
	
	Rep* enhance( ) const { return theRep; }
	// DANGEROUS: for altering an object without triggering cloning.
	// Use to change theRep without altering semantics.
	
	Rep* change( ) {
		if ( theRep->lastRef() ) return theRep;
		else { theRep->delRef(); return theRep = (Rep*)theRep->clone(); }
	}
	// For safe read/write access.
	
	void acquireRep( const Rep* rep )
	{
		((Rep*&)rep)->addRef();
		// cast away physical constness to permit logically const
		// incrementation of ref count
		if ( theRep->lastRef() ) delete theRep; else theRep->delRef();
		theRep = ((Rep*&)rep);
		// cast away physical constness of representation for
		// acquisition through new object; semantics of look() and
		// and change() guarantee that logical constness is maintained
	}
	// special `assignment-like' member: to permit an object to acquire
	// another's representation consistently
	// useful in delegation of object-level methods
	
	///////////////////////////////////////////////////////////
	//                                                       //
	// Special Constructors:                                 //
	//                                                       //
	///////////////////////////////////////////////////////////
	
	ObjectOf( Rep* newrep ) { theRep = newrep; }
	// To wrap new representations
	
private:
	
	///////////////////////////////////////////////////////////
	//                                                       //
	// Data Members:                                         //
	//                                                       //
	///////////////////////////////////////////////////////////
	
	Rep* theRep;
	
	///////////////////////////////////////////////////////////
	//                                                       //
	// Safety:                                               //
	//                                                       //
	///////////////////////////////////////////////////////////
	
	void force_derivation( ) { RefCounter* rc = theRep; }
	// With this member RefCounter is forced to be 
	// an accessible base class of Rep
	
};


#endif
