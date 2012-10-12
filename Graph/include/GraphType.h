// Copyright (C) 2005 Alexander Ushakov
// Contents: Definition of the GraphVertex and GraphEdge
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//


#ifndef _GraphType_H_
#define _GraphType_H_


#include "set"
#include "vector"
#include "iostream"
#include "string"
using namespace std;

#include "GraphConcept.h"
#include "GraphConceptAlgorithms.h"
#include "GraphDrawingAttributes.h"
using namespace Graphs;

#include "RanlibCPP.h"


/*!
  The order on edges must be defined so that the edges with the least labels go first for folding to be correct.
 */

//---------------------------------------------------------------------------//
//------------------------------- GraphEdge ---------------------------------//
//---------------------------------------------------------------------------//


//! Defines non-labelled edge of the directed graph.
struct GraphEdge
{
  //! Constructor for an edge. Argument t is a number of the target.
  GraphEdge( int t ) : theTarget(t) { }

  //! Dummy constructor, not to be used, required for STL::map. Also, can be used to denote "failure" or "dead-end" edge.
  GraphEdge( ) : theTarget(-1) { }


 
  //! Invert an edge.
  GraphEdge inverse( int origin ) { return GraphEdge( origin ); }


  //! Check if one edge is less than the other.
  bool operator < ( const GraphEdge& e ) const { return theTarget< e.theTarget; }


  //! Check if two edges equal.
  bool operator ==( const GraphEdge& e ) const { return theTarget==e.theTarget; }


  //! Check if two edges are not equal.
  bool operator!= ( const GraphEdge& e ) const { return theTarget!=e.theTarget; }

  
  //! The target of the edge (for out edges), or origin of the vertex (for in edges).
  int theTarget;
};

ostream& operator << ( ostream& os , const GraphEdge& e );


//---------------------------------------------------------------------------//
//------------------------------- GraphVertex -------------------------------//
//---------------------------------------------------------------------------//


//! Graph vertex.
template< class EdgeType >
struct GraphVertex
{
  //! Type of graph edges.
  typedef EdgeType edge_type;
  
  
  //! Default vertex.
  GraphVertex( ) { }

  
  //! Incoming edges.
  set< edge_type > in;

  
  //! Leaving edges.
  set< edge_type > out;
};


//---------------------------------------------------------------------------//
//----------------------------- IntLabeledEdge ------------------------------//
//---------------------------------------------------------------------------//


//! Defines labelled (by integer) edge of the directed graph.
struct IntLabeledEdge : virtual public GraphEdge
{
  //! Constructor for an edge. Argument t is a number of the target.
  IntLabeledEdge( int target , int label ) : GraphEdge(target), theLabel( label ) { }

  //! Dummy constructor, not to be used, required for STL::map. Also, can be used to denote "failure" or "dead-end" edge.
  IntLabeledEdge( ) : GraphEdge(-1), theLabel( 0 ) { }
  
  
  //! Invert an edge.
  IntLabeledEdge inverse( int origin ) { return IntLabeledEdge( origin , -theLabel ); }
  
  //! Check if one edge is less than the other.
  bool operator < ( const IntLabeledEdge& e ) const { 
    if( theLabel!=e.theLabel )
      return theLabel<e.theLabel;
    return theTarget<e.theTarget; 
  }
  
  
  //! Check if two edges equal.
  bool operator ==( const IntLabeledEdge& e ) const { return this->GraphEdge::operator==( e ) && theLabel==e.theLabel; }


  //! Check if two edges are not equal.
  bool operator!= ( const IntLabeledEdge& e ) const { return this->GraphEdge::operator!=( e ) || theLabel!=e.theLabel; }


  //! The label of the edge.  
  int theLabel;
};


//! Output operator for IntLabeledEdge.
ostream& operator << ( ostream& os , const IntLabeledEdge& e );


//! Function preparing two edges to a fold (for IntLabeledEdge does nothing).
void prepareToFold( const IntLabeledEdge& e1 , const IntLabeledEdge& e2 );


//---------------------------------------------------------------------------//
//----------------------------- PlanarGraphEdge -----------------------------//
//---------------------------------------------------------------------------//


//! Defines non-labelled edge of the directed planar graph.
struct PlanarGraphEdge : virtual public GraphEdge
{
  
  //! Constructor for an edge. Argument t is a number of the target.
  PlanarGraphEdge( int target , int cell1=-1 , int cell2=-1 ) : GraphEdge(target), theCell1( cell1 ), theCell2( cell2 ) { }
  
  
  //! Dummy constructor, not to be used, required for STL::map. Also, can be used to denote "failure" or "dead-end" edge.
  PlanarGraphEdge( ) : GraphEdge(-1), theCell1( -1 ), theCell2( -1 ) { }
  
  
  //! Invert an edge.
  PlanarGraphEdge inverse( int origin ) { return PlanarGraphEdge( origin , theCell2 , theCell1 ); }
  
  
  //! Check if one edge is less than the other.
  bool operator < ( const PlanarGraphEdge& e ) const { 
    if( theCell1!=e.theCell1 )
      return theCell1<e.theCell1;
    if( theCell2!=e.theCell2 )
      return theCell2<e.theCell2;
    return theTarget<e.theTarget;
  }
  
  
  //! Check if two edges equal.
  bool operator ==( const PlanarGraphEdge& e ) const { return theTarget==e.theTarget && theCell1==e.theCell1 && theCell2==e.theCell2; }
  
  
  //! Check if two edges are not equal.
  bool operator!= ( const PlanarGraphEdge& e ) const { return theTarget!=e.theTarget || theCell1!=e.theCell1 || theCell2!=e.theCell2; }
  
  
  //! The left neighbor-cell of the edge. Numbering of cells starts from 0, -1 denotes the outer face.
  int theCell1;
  
  
  //! The right neighbor-cell of the edge. Numbering of cells starts from 0, -1 denotes the outer face.
  int theCell2;
};


//---------------------------------------------------------------------------//
//------------------------ PlanarGraphIntLabelledEdge -----------------------//
//---------------------------------------------------------------------------//


//! Defines non-labelled edge of the directed planar graph.
struct PlanarGraphIntLabelledEdge : public PlanarGraphEdge, public IntLabeledEdge
{
  
  //! Constructor for an edge. Argument t is a number of the target.
  PlanarGraphIntLabelledEdge( int target , int label , int cell1=-1 , int cell2=-1 ) : PlanarGraphEdge(target,cell1,cell2), IntLabeledEdge(target,label) { }
  
  
  //! Dummy constructor, not to be used, required for STL::map. Also, can be used to denote "failure" or "dead-end" edge.
  PlanarGraphIntLabelledEdge( ) { }
  
  
  //! Invert an edge.
  PlanarGraphIntLabelledEdge inverse( int origin ) { return PlanarGraphIntLabelledEdge( origin , -theLabel , theCell2 , theCell1 ); }
  
  
  //! Check if one edge is less than the other.
  bool operator < ( const PlanarGraphIntLabelledEdge& e ) const { 
    if( theLabel!=e.theLabel )
      return theLabel<e.theLabel;
    if( theCell1!=e.theCell1 )
      return theCell1<e.theCell1;
    if( theCell2!=e.theCell2 )
      return theCell2<e.theCell2;
    return theTarget<e.theTarget;
  }
  
  
  //! Check if two edges equal.
  bool operator ==( const PlanarGraphIntLabelledEdge& e ) const { return theLabel==e.theLabel && theTarget==e.theTarget && theCell1==e.theCell1 && theCell2==e.theCell2; }
  
  
  //! Check if two edges are not equal.
  bool operator!= ( const PlanarGraphIntLabelledEdge& e ) const { return theLabel!=e.theLabel || theTarget!=e.theTarget || theCell1!=e.theCell1 || theCell2!=e.theCell2; }
};


//---------------------------------------------------------------------------//
//----------------------------- IntLabeledEdge ------------------------------//
//---------------------------------------------------------------------------//



//! Directed Graph
typedef GraphConcept< GraphVertex< GraphEdge > , GraphEdge > Graph;



//! Directed Graph with integer labelled edges
typedef GraphConcept< GraphVertex< IntLabeledEdge > , IntLabeledEdge > IntLabeledGraph;




//---------------------------------------------------------------------------//
//---------------------------- Specific functions ---------------------------//
//---------------------------------------------------------------------------//


//! Function attaches a sequence of edges to a graph consequently labelled by numbers [B,E)
template< class ConstIntIterator >
int addRay( IntLabeledGraph& G , int v , ConstIntIterator B , ConstIntIterator E , bool direct=false )
{
  int cur_v = v;
  for( ; B!=E ; ) {
    int l = *(B++);
    int new_v = G.newVertex( IntLabeledGraph::vertex_type( ) );
    G.newEdge( cur_v , IntLabeledGraph::edge_type( new_v ,  l ) );
    if( !direct )
      G.newEdge( new_v , IntLabeledGraph::edge_type( cur_v , -l ) );
    cur_v = new_v;
  }
  return cur_v;
}


//! Function attaches a loop of edges to a graph consequently labelled by numbers [B,E)
template< class ConstIntIterator >
void addLoop( IntLabeledGraph& G , int v , ConstIntIterator B , ConstIntIterator E , bool direct=false )
{
  int cur_v = v;
  for( ; B!=E ; ) {
    int l = *(B++);
    if( B!=E ) {
      int new_v = G.newVertex( IntLabeledGraph::vertex_type( ) );
      G.newEdge( cur_v , IntLabeledGraph::edge_type( new_v ,  l ) );
      if( !direct )
	G.newEdge( new_v , IntLabeledGraph::edge_type( cur_v , -l ) );
      cur_v = new_v;
    } else {
      G.newEdge( cur_v , IntLabeledGraph::edge_type( v ,  l ) );
      if( !direct )
	G.newEdge( v , IntLabeledGraph::edge_type( cur_v , -l ) );
    }
  }
}


//! Create random graph on N verteices where each edge has probability to appear = edge_param.
Graph randomGraph( int N , float edge_param );


//! Compute distances between all pairs of vertices.
vector< vector< int > > lengthTable( const Graph& G );


//! Compute Gromov's inner product for all pairs of vertices (the origin is the vertex with the smallest number)
vector< vector< int > > innerProductTable( const Graph& G , int origin );


//! Compute constant of hyperbolicity of a graph (uses Gromov's inner product).
float getHyperbolicityConst( const Graph& G );


//! Get Graphviz graph description.
string graphviz_format( const Graph& G , const GraphDrawingAttributes& GDA = GraphDrawingAttributes( ) );


//! Get Graphviz graph description.
string graphviz_format( const IntLabeledGraph& G , const GraphDrawingAttributes& GDA = GraphDrawingAttributes( ) );


#endif


