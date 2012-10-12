
#ifndef _AdvDehnAlgorithm_h_
#define _AdvDehnAlgorithm_h_

#include "FPGroup.h"

#include "GraphType.h"
#include "GraphConcept.h"
#include "GraphConceptAlgorithms.h"
using namespace Graphs;


//---------------------------------------------------------------------------//
//--------------------------- AdvDehnAlgorithm ------------------------------//
//---------------------------------------------------------------------------//

/*!
  Requires a revision.
  set< int > checkedStates;
  must be replaced with 
  int lastCheckedVertex;
 */

class AdvDehnAlgorithm 
{

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Constructor:                                       //
  //                                                     //
  /////////////////////////////////////////////////////////

 public:
  
  AdvDehnAlgorithm( const FPGroup& G , const Word& w );
  AdvDehnAlgorithm( const FPGroup& G , const set< Word >& gens , const Word& w );


  /////////////////////////////////////////////////////////
  //                                                     //
  //  Accessors:                                         //
  //                                                     //
  /////////////////////////////////////////////////////////

 public:
  
  const IntLabeledGraph& getFSA( ) const { return theFSA; }
  bool builtup( set< Word >* conj=0 , int coset_limit=100000 );
  // we need conjugators to construct a generating set for the automaton
  
  bool isLoop( const Word& w ) const;
  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Manipulators:                                      //
  //                                                     //
  /////////////////////////////////////////////////////////

 private:
  
  void addCycle( const Word& w , int origin );
  // function ported from my dehn.C file
  // I should unify it somehow
  
  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Data members:                                      //
  //                                                     //
  /////////////////////////////////////////////////////////

 private:
  
  const FPGroup theGroup;
  const Word theWord;  

  set< int > checkedStates;
  
  IntLabeledGraph theFSA;
};




#endif
