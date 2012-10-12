// Copyright (C) 2005 Alexander Ushakov
// Contents: Definition of class WhiteheadGraph
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//

#ifndef _WhiteheadGraph_h_
#define _WhiteheadGraph_h_


#include "GraphType.h"
#include "GraphConcept.h"
#include "GraphConceptAlgorithms.h"
using namespace Graphs;

#include "Word.h"


//---------------------------------------------------------------------------//
//---------------------------- WhiteheadGraph -------------------------------//
//---------------------------------------------------------------------------//


class WhiteheadGraph
{
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Constructors:                                      //
  //                                                     //
  /////////////////////////////////////////////////////////

public:
  WhiteheadGraph( int n_gens , const Word& w ) { assign( n_gens , w.begin( ) , w.end( ) ); }

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Accessors:                                         //
  //                                                     //
  /////////////////////////////////////////////////////////

public:
  template< class ConstIterator > void assign( int n_gens , ConstIterator B , ConstIterator E ) {

    if( B==E ) return;
    theGraph.clear( );
    for( int i=0 ; i<n_gens ; ++i ) {
      theGraph.newVertex( );
      theGraph.newVertex( );
    }
    
    int old_l = *--E;
    int old_num = old_l>0 ? 2*old_l-1 : -2*old_l-2;
    for( ++E ; B!=E ; ++B ) {
      int l = *B;
      int num = l>0 ? 2*l-1 : -2*l-2;
      theGraph.newEdge( old_num , num );
      // cout << old_num << ", " << num << endl;
      old_l = l;
      old_num = num;
    }
  }
  
  const Graph& getGraph( ) { return theGraph; }
	

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Data members:                                      //
  //                                                     //
  /////////////////////////////////////////////////////////

private:

  Graph theGraph;
};


#endif
