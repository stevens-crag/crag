
#ifndef _DehornoyForm_h_
#define _DehornoyForm_h_


#include "Word.h"


//---------------------------------------------------------------------------//
//----------------------------- DehornoyForm --------------------------------//
//---------------------------------------------------------------------------//


class DehornoyForm
{

  ///////////////////////////////////////////////////////// 
  //                                                     //
  //  Constructors:                                      //
  //                                                     //
  /////////////////////////////////////////////////////////
    
 public:

  DehornoyForm( const int N , const Word& w );


  /////////////////////////////////////////////////////////
  //                                                     //
  //  Accessors:                                         //
  //                                                     //
  /////////////////////////////////////////////////////////

 public:

  Word getDehornoyForm( ) const { return theDehornoyForm; }


  /////////////////////////////////////////////////////////
  //                                                     //
  //  Internal functions:                                //
  //                                                     //
  /////////////////////////////////////////////////////////

 private:

  Word computeDehornoyForm( const Word& w );

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Data members:                                      //
  //                                                     //
  /////////////////////////////////////////////////////////

 private:
  
  //! the rank of the braid group = number of strands = number of generators + 1
  const int theIndex;

  Word theDehornoyForm;
};

#endif

