// Copyright (C) 2005 Alexander Ushakov
// Contents: Definition of class DehornoyForm
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//

#ifndef _DehornoyForm_h_
#define _DehornoyForm_h_


#include "Word.h"
#include "LinkedBraidStructure.h"


//---------------------------------------------------------------------------//
//----------------------------- DehornoyForm --------------------------------//
//---------------------------------------------------------------------------//

//! Dehornoy Form of a braid word (aka/ handle free form)
/*!
  This class uses BraidLinkedStructure to compute the Dehornoy form a braid-word.
  The object of this class is just a container to keep the form.
 */

class DehornoyForm
{

  ///////////////////////////////////////////////////////// 
  //                                                     //
  //  Constructors:                                      //
  //                                                     //
  /////////////////////////////////////////////////////////
    
 public:

  //! Create a Dehornoy form of a braid word
  /*!
    Note that the form will be computed right here. So, if you are not sure that 
    the form will be used later do not create this object.
    \param N - rank of a braid group;
    \param w - braid word
   */
  DehornoyForm( const int N , const Word& w );

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Accessors:                                         //
  //                                                     //
  /////////////////////////////////////////////////////////

 public:

  //! (Accessor function) Get a braid word representing a form.
  Word getDehornoyForm( ) const { return theDehornoyForm; }


  /////////////////////////////////////////////////////////
  //                                                     //
  //  Internal functions:                                //
  //                                                     //
  /////////////////////////////////////////////////////////

 private:

  //! Function which computes the form (using LinkedBraidStructure).
  Word computeDehornoyForm( const Word& w );

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Data members:                                      //
  //                                                     //
  /////////////////////////////////////////////////////////

 private:

  //! The rank of the corresponding braid group
  const int theIndex;

  //! The Dehornoy form.
  Word theDehornoyForm;

  //! Structure in which the form is kept
  LinkedBraidStructure theLinkedStructure;
};

#endif

