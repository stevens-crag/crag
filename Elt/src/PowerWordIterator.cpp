// Copyright (C) 2005 Alexander Ushakov
// Contents: Implementation of class PowerWordIterator
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
	theWord( 0 ),
	theOffset( 0 )
{

}

ConstWordIterator::ConstWordIterator( const WordIterator& WI ) :
	theWord(WI.theWord),
	theIterator(WI.theIterator),
	theOffset(WI.theOffset)
{

}

ConstWordIterator& ConstWordIterator::operator= ( const WordIterator& WI )
{
	theWord = WI.theWord;
	theIterator = WI.theIterator;
	theOffset = WI.theOffset;
}


ConstWordIterator::ConstWordIterator( const Word& w , bool begin ) : 
	theWord( &w ),
	theIterator( begin ? w.getList( ).begin( ) : w.getList( ).end() ),
	theOffset( 0 )
{

}


const ConstWordIterator& ConstWordIterator::operator++ ( )
{
	if( theOffset<abs( (*theIterator).second )-1 )
		++theOffset;
	else {
		++theIterator;
		theOffset = 0;
	}
	return *this;
}


ConstWordIterator ConstWordIterator::operator++ ( int doomy )
{
	ConstWordIterator WI = *this;
	if( theOffset<abs( (*theIterator).second )-1 )
		++theOffset;
	else {
		++theIterator;
		theOffset = 0;
	}
	return WI;
}


const ConstWordIterator& ConstWordIterator::operator-- ( )
{
	if( theOffset==0 )
		theOffset = abs( (*--theIterator).second)-1;
	else
		--theOffset;
	return *this;
}


ConstWordIterator ConstWordIterator::operator-- ( int doomy )
{
	ConstWordIterator WI = *this;
	if( theOffset==0 )
		theOffset = abs( (*--theIterator).second)-1;
	else
		--theOffset;
	return WI;
}


pair< int , int > ConstWordIterator::operator*  ( ) const
{
	int p = (*theIterator).second;
	p = p<0 ? -1 : 1;
	return PII( (*theIterator).first , p );
}


bool ConstWordIterator::operator== ( const ConstWordIterator& WI ) const
{
	if( theIterator!=WI.theIterator )
		return false;
	if( theOffset!=WI.theOffset )
		return false;
	return true;
}


bool ConstWordIterator::operator!= ( const ConstWordIterator& WI ) const
{
	if( operator==( WI ) )
		return false;
	return true;
}


//---------------------------------------------------------------------------//
//------------------------------- WordIterator ------------------------------//
//---------------------------------------------------------------------------//


WordIterator::WordIterator( ) :
	theWord( 0 ),
	theOffset( 0 )
{

}


WordIterator::WordIterator( Word& w , bool begin ) : 
	theWord( &w ),
	theIterator( begin ? w.getList( ).begin( ) : w.getList( ).end( ) ),
	theOffset( 0 )
{

}


const WordIterator& WordIterator::operator++ ( )
{
	if( theOffset<abs( (*theIterator).second )-1 )
		++theOffset;
	else {
		++theIterator;
		theOffset = 0;
	}
	return *this;
}


WordIterator WordIterator::operator++ ( int doomy )
{
	WordIterator WI = *this;
	if( theOffset<abs( (*theIterator).second )-1 )
		++theOffset;
	else {
		++theIterator;
		theOffset = 0;
	}
	return WI;
}


const WordIterator& WordIterator::operator-- ( )
{
	if( theOffset==0 )
		theOffset = abs( (*--theIterator).second)-1;
	else
		--theOffset;
	return *this;
}


WordIterator WordIterator::operator-- ( int doomy )
{
	WordIterator WI = *this;
	if( theOffset==0 )
		theOffset = abs( (*--theIterator).second)-1;
	else
		--theOffset;
	return WI;
}


pair< int , int > WordIterator::operator*  ( ) const
{
	int p = (*theIterator).second;
	p = p<0 ? -1 : 1;
	return PII( (*theIterator).first , p );
}


bool WordIterator::operator== ( const WordIterator& WI ) const
{
	if( theIterator!=WI.theIterator )
		return false;
	if( theOffset!=WI.theOffset )
		return false;
	return true;
}


bool WordIterator::operator!= ( const WordIterator& WI ) const
{
	if( operator==( WI ) )
		return false;
	return true;
}
