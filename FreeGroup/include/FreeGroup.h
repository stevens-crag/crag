// Copyright (C) 2005 Alexander Ushakov
// Contents: Definition of class FreeGroup
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//


#ifndef _FreeGroup_h_
#define _FreeGroup_h_


#include "SubgroupFG.h"
#include "Alphabet.h"
#include "Word.h"


//---------------------------------------------------------------------------//
//-------------------------------- FreeGroup --------------------------------//
//---------------------------------------------------------------------------//


class FreeGroup
{
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Constructors:                                      //
  //                                                     //
  /////////////////////////////////////////////////////////
  
 public:

  //! (Constructor) Free group of the given rank.
  FreeGroup( int rank );
  FreeGroup( const FiniteAlphabet& a );
  

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Accessors:                                         //
  //                                                     //
  /////////////////////////////////////////////////////////

 public:

  //! Determine if the word w is primitive or not
  bool isPrimitive( const Word& w ) const;

  //! Determine if the word w is almost primitive or not
  bool isAlmostPrimitive( const Word& w ) const;

  //! Determine whether a subgroup sbgp of a free group contains a word w
  bool doesContain( const SubgroupFG& sbgp , const Word& w ) const;

  const FiniteAlphabet& getAlphabet()const { return theAlphabet; }
  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  I/O operators                                      //
  //                                                     //
  /////////////////////////////////////////////////////////

   //! Output the presentation into a stream 
  friend ostream& operator << ( ostream& out, const FreeGroup& g ){
    cout << "< ";
    for (int i=0;i<g.theRank-1;i++)
      out << g.theAlphabet.getLetter(i+1) << ", ";

    cout << g.theAlphabet.getLetter(g.theRank) << " >" << flush;
    return out;
  }
  
  //! Read an alphabet from a string 
  friend istream& operator >> ( istream& in,  FreeGroup& g ){
    FiniteAlphabet a;
    in >> a;
    g.theAlphabet = a;
    g.theRank = a.size();
    return in;
  }

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Data members:                                      //
  //                                                     //
  /////////////////////////////////////////////////////////

 private:

  int theRank;
  FiniteAlphabet theAlphabet;
  bool useDefaultAlphabet;

};

#endif

