// Copyright (C) 2000 Alexander Ushakov
// Contents: Definition of class VectorEnumerator
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//


#ifndef _VECTORENUMERATOR_H_
#define _VECTORENUMERATOR_H_


#include <vector>
#include <iostream>

using namespace std;


//---------------------------------------------------------------------------//
//--------------------------- VectorEnumerator ------------------------------//
//---------------------------------------------------------------------------//


//! Class VectorEnumerator
/*!
  VectorEnumerator provides the ground-level interface for vector-enumerators
  (for instance it can be used for enumeration of permutations, finite sequences over finite alphabet, etc.).
*/

class VectorEnumerator
{

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Constructors                                       //
  //                                                     //
  /////////////////////////////////////////////////////////
  
 public:

  //! Default constructor
  VectorEnumerator( ) : curLength( 0 ) { }
    
    
  // copy constructor, assignment operator and destructor are
  // supplied by compiler
  
  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Accessors                                          //
  //                                                     //
  /////////////////////////////////////////////////////////
  
 public:
  
  //! Compute the next element.
  VectorEnumerator& operator++( );

  
  //! Get the current element
  vector< int > operator* ( );
  

  //! Returns true if there are no more elements to enumerate
  bool end( ) { 

    if( curLength==-1 )
      return true;
    
    // if sequence is not completed we need to complete it
    if( !seqComplete( ) )
      ++(*this);	
    
    return (curLength==-1); 
  }
  
  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Debug:                                             //
  //                                                     //
  /////////////////////////////////////////////////////////
  
  ostream& printSeq( ostream& os ) const;
  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Internal functions:                                //
  //                                                     //
  /////////////////////////////////////////////////////////
  
 protected:
  
  //! Check if the currently constructed sequence is legitimate
  virtual bool seqOK            ( ) const = 0;
  //! Returns true if the length-limit on the sequence is passed
  virtual bool seqLimit         ( ) const = 0;
  //! Returns true if the current sequence is complete
  virtual bool seqComplete      ( ) const = 0;

  //! Returns the initial value for the current element of the sequence (usually 0)
  virtual int  start (         ) const = 0;
  //! Returns the next value for the current element of the sequence (usually cur+1)
  virtual int  next  ( int cur ) const = 0;
  //! 
  virtual bool finish( int cur ) const = 0;
  
  //! Get the length of the current sequence
  int                  getLength( ) const { return curLength; }
  //! Get the current sequence
  const vector< int >& getSeq   ( ) const { return curVector; }

  //! Function is invoked when the current value of the current element of the sequence is accepted and we go for the next element
  /*!
    The value of curElement is already increased, so you are at the next element.
   */
  virtual void stepTo( ) { }
  //! Function is invoked when we go backward from the current element to the previous element 
  /*! 
    The value of curElement is already decreased, so you are at the previous element
   */
  virtual void stepBack( int cur ) { }
  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Data members                                       //
  //                                                     //
  /////////////////////////////////////////////////////////
  
 private:
  
  int           curLength;
  vector< int > curVector;
  
};



#endif
