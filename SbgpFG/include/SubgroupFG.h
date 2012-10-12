// Copyright (C) 2005 Alexander Ushakov
// Contents: Definition of class SubgroupFG
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//

#ifndef _SubgroupFG_h_
#define _SubgroupFG_h_


#include "Word.h"
#include "GraphType.h"
#include "GraphConcept.h"
#include "GraphConceptAlgorithms.h"
using namespace Graphs;


#include "vector"
using namespace std;


//---------------------------------------------------------------------------//
//------------------------------ SubgroupFG ---------------------------------//
//---------------------------------------------------------------------------//


class SubgroupFG
{
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Constructors                                       //
  //                                                     //
  /////////////////////////////////////////////////////////
public:
  SubgroupFG( int n_gens = 0 );
  SubgroupFG( int n_gens , const vector< Word >& gens );
  
  template< class ConstWordIterator > SubgroupFG( int n_gens , ConstWordIterator B , ConstWordIterator E ) : 
    theNumberOfGenerators( n_gens ),
    fsaDone( false ),
    nielsDone( false )
      {
	int sz = 0;
	ConstWordIterator C = B;
	for( ; C!=E ; ++C, ++sz );
	theGenerators = vector< Word >( sz );
	for( int i=0 ; B!=E ; ++B, ++i )
	  theGenerators[i] = *B;
      }

 protected:
  SubgroupFG( int n_gens , const IntLabeledGraph& fsa );


  /////////////////////////////////////////////////////////
  //                                                     //
  //  Operators:                                         //
  //                                                     //
  /////////////////////////////////////////////////////////
public:

  //! Check the equality of two subgroups
  bool operator== ( const SubgroupFG& S ) const;


  //! Compute the intersection of two subgroups
  SubgroupFG operator* ( const SubgroupFG& S ) const;
  
  
  //! Conjugate a subgroup
  SubgroupFG& operator^= ( const Word& conjugator );
  
  
  //! Conjugate a subgroup
  SubgroupFG  operator^  ( const Word& conjugator ) const {
    SubgroupFG result = *this;
    result ^= conjugator;
    return result;
  }
  

  //! Extend subgroup basis
  SubgroupFG& operator+= ( const SubgroupFG& sbgp );


  //! Extend subgroup basis
  SubgroupFG  operator+  ( const SubgroupFG& sbgp ) const {
    SubgroupFG result = *this;
    result += sbgp;
    return result;
  }


  //! Extend subgroup basis
  SubgroupFG& operator+= ( const Word& w );


  //! Extend subgroup basis
  SubgroupFG  operator+  ( const Word& w ) const {
    SubgroupFG result = *this;
    result += w;
    return result;
  }

  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Accessors:                                         //
  //                                                     //
  /////////////////////////////////////////////////////////
public:
  
  //! Check if two subgroup graphs *this and S are isomorphic as automata, restricted that v1 -> v2
  bool checkIsomorphism( const SubgroupFG& S , int vert1 , int vert2 ) const;
  
  
  //! Get the initial generators of the subgroup.
  const vector< Word >& getGenerators( ) const { return theGenerators; }


  //! Get a set of Nielsen generators of the subgroup.
  const vector< Word >& getNielsenGenerators( ) const;
  
  
  //! Get an automaton corresponding to a subgroup of a free froup
  const IntLabeledGraph& getFSA( ) const;
  
  
  //! Compute the index of a subgroup of a free group, -1 means infinite
  int getIndex( ) const;


  //! Compute the rank of a subgroup of a free group
  int getRank( ) const;
  
  
  //! check if a word belongs to a subgroup of a free group
  bool doesBelong( const Word& w ) const;


  //! Express a word in a generators of a subgroup of a free group (word must belong to subgroup)
  Word express( const Word& w ) const;
  

  //! Compute the centralizer of a subgroup
  SubgroupFG centralizer( ) const;
  
  
  //! Compute the normalizer of a subgroup
  SubgroupFG normalizer( ) const;
  

  //! Trim the subgroup graph. (Cut off the "tail".)
  pair< SubgroupFG , Word > trim( ) const;
  
  
  //! Check if two subgroups are conjugate
  pair< bool , Word > areConjugate( const SubgroupFG& sbgp ) const;
  
  
  //! Returns the graphviz graphical description of the subgroup graph (save it to file and use Graphviz to visualize the graph)
  string graphviz_format( ) const;



  /////////////////////////////////////////////////////////
  //                                                     //
  //  Internal functions:                                //
  //                                                     //
  /////////////////////////////////////////////////////////
  
private:


  void computeFSA( ) const;
  void computeNielsenGenerators( ) const;
  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Data members                                       //
  //                                                     //
  /////////////////////////////////////////////////////////
  
private:

  vector< Word > theGenerators;
  int theNumberOfGenerators;
  
  mutable bool fsaDone;
  mutable IntLabeledGraph theFSA;
  mutable list< FoldDetails< IntLabeledGraph::vertex_type , IntLabeledGraph::edge_type > > foldDetails;
  
  mutable bool nielsDone;
  mutable vector< Word > theNielsenGenerators;
};


//---------------------------------------------------------------------------//
//------------------------------ Substitute ---------------------------------//
//---------------------------------------------------------------------------//


template < class ConstIterator >
Word substitute( ConstIterator B , ConstIterator E , const vector< Word >& wrds ) 
{
  Word result;
  for( ; B!=E ; ++B ) {
    int letter = *B;
    if( letter>0 )
      result *=  wrds[ letter-1];
    else
      result *= -wrds[-letter-1];
  }
  
  return result;
}


Word substitute( const Word& w , const vector< Word >& wrds );


//---------------------------------------------------------------------------//
//------------------------------ operator << --------------------------------//
//---------------------------------------------------------------------------//

ostream& operator << ( ostream& os , const SubgroupFG& sbgp );



#endif

