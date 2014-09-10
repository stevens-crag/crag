// Copyright (C) 2005 Alexander Ushakov
// Contents: Implementation of class WordIterator
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//

#include "WordIterator.h"
#include "Word.h"


//---------------------------------------------------------------------------//
//------------------------------- WordIterator ------------------------------//
//---------------------------------------------------------------------------//


ConstWordIterator::ConstWordIterator( ) :
theList( 0 )
{

}

ConstWordIterator::ConstWordIterator( const WordIterator& WI ) :
theList(WI.theList),
theIterator(WI.theIterator)
{

}

ConstWordIterator& ConstWordIterator::operator= ( const WordIterator& WI )
{
	theList = WI.theList;
	theIterator = WI.theIterator;
	return *this;
}


ConstWordIterator::ConstWordIterator( const Word& w , bool begin ) : 
theList( &w.getList() ),
theIterator( begin ? w.getList( ).begin( ) : w.getList( ).end() )
{

}


ConstWordIterator& ConstWordIterator::operator++ ( )
{
	++theIterator;
	return *this;
}


ConstWordIterator ConstWordIterator::operator++ ( int doomy )
{
	ConstWordIterator WI = *this;
	++theIterator;
	return WI;
}


ConstWordIterator& ConstWordIterator::operator-- ( )
{
	--theIterator;
	return *this;
}


ConstWordIterator ConstWordIterator::operator-- ( int doomy )
{
	ConstWordIterator WI = *this;
	--theIterator;
	return WI;
}


int ConstWordIterator::operator*  ( ) const
{
	return *theIterator;
}


bool ConstWordIterator::operator== ( const ConstWordIterator& WI ) const
{
	return theIterator==WI.theIterator;
}


bool ConstWordIterator::operator!= ( const ConstWordIterator& WI ) const
{
	return theIterator!=WI.theIterator;
}


//---------------------------------------------------------------------------//
//------------------------------- WordIterator ------------------------------//
//---------------------------------------------------------------------------//


WordIterator::WordIterator( ) :
theList( 0 )
{

}


WordIterator::WordIterator( Word& w , bool begin ) : 
theList( &w.getList( ) ),
theIterator( begin ? w.getList( ).begin( ) : w.getList( ).end( ) )
{

}

WordIterator::WordIterator( list<int>& lst , list< int >::iterator it ) :
theList( &lst ),
theIterator( it )
{

}


const WordIterator& WordIterator::operator++ ( )
{
	++theIterator;
	return *this;
}


WordIterator WordIterator::operator++ ( int doomy )
{
	WordIterator WI = *this;
	++theIterator;
	return WI;
}


const WordIterator& WordIterator::operator-- ( )
{
	--theIterator;
	return *this;
}


WordIterator WordIterator::operator-- ( int doomy )
{
	WordIterator WI = *this;
	--theIterator;
	return WI;
}


int WordIterator::operator*  ( ) const
{
	return *theIterator;
}


bool WordIterator::operator== ( const WordIterator& WI ) const
{
	return theIterator==WI.theIterator;
}


bool WordIterator::operator!= ( const WordIterator& WI ) const
{
	return theIterator!=WI.theIterator;
}
