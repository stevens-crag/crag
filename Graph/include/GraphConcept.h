// Copyright (C) 2005 Alexander Ushakov
// Contents: Definition of the graph concept
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//


#ifndef _GraphConcept_h_
#define _GraphConcept_h_


#include "RefCounter.h"
#include "ObjectOf.h"


#include "set"
#include "map"
using namespace std;

/*!
  EdgeType must support: 
  e.inverse( v ), where v is the number of the origin
 */


//---------------------------------------------------------------------------//
//------------------------------- GraphConcept ------------------------------//
//---------------------------------------------------------------------------//



namespace Graphs
{
  
  template< class VertexType , class EdgeType > class GraphConcept;
  template< class VertexType , class EdgeType > class GraphConceptRep;
  
  
  //---------------------------------------------------------------------------//
  //-------------------------------- GraphRep ---------------------------------//
  //---------------------------------------------------------------------------//
  
  
  //! Representation class for graph types.
  template< class VertexType , class EdgeType > class GraphConceptRep : public RefCounter 
  {
      /////////////////////////////////////////////////////////
      //                                                      //
      //  Constructors                                       //
      //                                                     //
      /////////////////////////////////////////////////////////
      public:
    
    typedef VertexType vertex_type;
    typedef EdgeType   edge_type;

    friend class GraphConcept< vertex_type , edge_type >;

      /////////////////////////////////////////////////////////
      //                                                     //
      //  Constructors                                       //
      //                                                     //
      /////////////////////////////////////////////////////////
      private:
    
    GraphConceptRep( ) : maxVertex( 0 ) { }


      /////////////////////////////////////////////////////////
      //                                                     //
      //  Accessors                                          //
      //                                                     //
      /////////////////////////////////////////////////////////
      public:

    GraphConceptRep* clone( ) const { return new GraphConceptRep( *this ); }
  
      /////////////////////////////////////////////////////////
      //                                                     //
      //  Internal functions                                 //
      //                                                     //
      /////////////////////////////////////////////////////////
      private:
    
    //! Get a set of all vertices
    const map< int, vertex_type >& getVertices( ) const { return theVertices; }
    
    //! Add a new vertex to a graph disconnected from all others.
    int newVertex( ) {
      theVertices[maxVertex] = vertex_type( );
      return maxVertex++;
    }
    //! Add a new vertex to a graph. All edges in "in" and "out" sets will be incorporated to the graph.
    int newVertex( const vertex_type& V ) {
      _newVertex( maxVertex , V );
      return maxVertex++;
    }
    //! Replace a vertex v with V. If v does not exist then V simply becomes a vertex with the number v.
    void replaceVertex( int v , const vertex_type& V ) {
      eraseVertex( v );
      _newVertex( v , V );
    }

    //! Erase a vertex from a graph.
    void eraseVertex( int v ) {
      const vertex_type& V = theVertices[v];
      set< edge_type >  in = V. in;
      set< edge_type > out = V.out;
      
      typename set< edge_type >::iterator it;
      
      // state <- (*it).first by (*it).second
      for( it=in.begin( ) ; it!=in.end( ) ; ++it ) {
	edge_type edge = *it;
	int origin = edge.theTarget;
	edge.theTarget = v;
	theVertices[origin].out.erase( edge );
      }
      
      // state2 -> (*it).first by (*it).second
      for( it=out.begin( ) ; it!=out.end( ) ; ++it ) {
	edge_type edge = *it;
	int target = edge.theTarget;
	edge.theTarget = v;
	theVertices[target].in.erase( edge );
      }
      theVertices.erase( v );
    }
    
    //! Add a new edge to a graph.
    void newEdge( int v , const edge_type& E ) {
      int t = E.theTarget;
      theVertices[v].out.insert( E );
      edge_type back = E;
      back.theTarget = v;
      theVertices[t]. in.insert( back );
    }

    //! Erase an edge from a graph.
    void eraseVertex( int v , const edge_type& E ) {
      int t = E.theTarget;
      theVertices[v].out.erase( E );
      edge_type back = E;
      back.theTarget = v;
      theVertices[t]. in.erase( back );
    }
    
    
    //! Clear a graph.
    void clear( ) {
      theVertices.clear( );
      maxVertex = 0;
    }
    
    void pinch( int v1 , int v2 ) {
      if( v1==v2 ) return;
      if( v1>v2 ) swap( v1 , v2 );
      const vertex_type& V2 = theVertices[v2];
      set< edge_type >  in = V2. in;
      set< edge_type > out = V2.out;
      
      typename set< edge_type >::iterator it;
      
      // V2 <- *it
      for( it=in.begin( ) ; it!=in.end( ) ; ++it ) {
	edge_type edge = *it;
	int origin = edge.theTarget;
	if( origin!=v2 ) {
	  edge.theTarget = v2;
	  // cout << "    Erase  from out N " << origin << " : " << edge << endl;
	  theVertices[origin].out.erase( edge );
	  edge.theTarget = v1;
	  // cout << "    Insert to   out N " << origin << " : " << edge << endl;
	  newEdge( origin , edge );
	  // theVertices[origin].out.insert( edge );
	} else {
	  edge.theTarget = v1;
	  theVertices[v1].out.insert( edge );
	}
      }
      
      // V2 -> *it
      for( it=out.begin( ) ; it!=out.end( ) ; ++it ) {
	edge_type edge = *it;
	int target = edge.theTarget;
	
	if( target!=v2 ) {
	  // add a edge v1 -> target
	  newEdge( v1 , edge );
	  // remove old incoming edge v2 -> target
	  edge.theTarget = v2;
	  // cout << "    Erase  from in N " << target << " : " << edge << endl;
	  theVertices[target].in.erase( edge );
	  // cout << "    Insert to   in N " << target << " : " << edge << endl;
	} else {
	  edge.theTarget = v1;
	  theVertices[v1].in.insert( edge );
	}
      }
      theVertices.erase( v2 );
    }

      /////////////////////////////////////////////////////////
      //                                                     //
      //  Internal Functions:                                //
      //                                                     //
      /////////////////////////////////////////////////////////
      private:
    
    //! This function assumes that the vertex v does not exist and V has a correct structure to be incorporated to a graph.
    int _newVertex( int v , const vertex_type& V ) {
      
      // add V into the graph
      theVertices[v] = V;
      
      // link new edges incoming to V and leaving V
      set< edge_type >  in = V. in;
      set< edge_type > out = V.out;
      typename set< edge_type >::iterator it;

      // V <- *it
      for( it=in.begin( ) ; it!=in.end( ) ; ++it ) {
	edge_type edge = *it;
	int origin = edge.theTarget;
	if( origin!=v ) {
	  edge.theTarget = v;
	  theVertices[origin].out.insert( edge );
	}
      }

      // V -> *it
      for( it=out.begin( ) ; it!=out.end( ) ; ++it ) {
	edge_type edge = *it;
	int target = edge.theTarget;
	if( target!=v ) {
	  edge.theTarget = v;
	  theVertices[target].in.insert( edge );
	}
      }
    }

    
      /////////////////////////////////////////////////////////
      //                                                     //
      //  Data members                                       //
      //                                                     //
      /////////////////////////////////////////////////////////
      private:

    // the set of vertices in the graph
    map< int, vertex_type > theVertices;
    
    //! The maximal unused number of a vertex (used in newVertex)
    int maxVertex;
    
    
  };


  //---------------------------------------------------------------------------//
  //---------------------------------- Graph ----------------------------------//
  //---------------------------------------------------------------------------//
  

  //! The main class for graph types.
  /* Graph concept provides main functions and data structures for directed graph types.
     Examples of directed graph types are: directed graph, directed graph with labelled vertices and edges
     and different combinations. By definition Graph = (Vertices,Edges). Changing types of template parameters
     VertexType and EdgeType will give different types of graphs.
   */

  template< class VertexType , class EdgeType > class GraphConcept : public ObjectOf< GraphConceptRep< VertexType , EdgeType > >
  {
    
    
  public:
    //! The type of the graph vertex.
    typedef VertexType                               vertex_type;
    //! The type of the graph edge.
    typedef EdgeType                                 edge_type;
    
    
  private:
    typedef GraphConceptRep< VertexType , EdgeType > presentation_type;

      /////////////////////////////////////////////////////////
      //                                                     //
      //  Constructors                                       //
      //                                                     //
      /////////////////////////////////////////////////////////
      public:
    
    //! Create empty graph.
    GraphConcept( ) : ObjectOf< presentation_type >( new presentation_type( ) ) { }

      /////////////////////////////////////////////////////////
      //                                                     //
      //  Accessors                                          //
      //                                                     //
      /////////////////////////////////////////////////////////
      public:


    //! Get the vertex set.
    const map< int, vertex_type >& getVertices( ) const { return this->look( )->getVertices( ); }


    //! Create new vertex (default). Function returns the number of the new vertex.
    int  newVertex( ) { return this->change( ) -> newVertex( ); }
    
    
    //! Add new vertex into a graph. Function returns the number of the new vertex.
    int  newVertex( const vertex_type& V ) { return this->change( ) -> newVertex( V ); }

    
    //! Erase a vertex from the graph
    void eraseVertex( int v ) { return this->change( ) -> eraseVertex( v ); }
    
    
    //! Function locally changes graph by replacing the vertex with number v by the vertex V.
    void replaceVertex( int v , const vertex_type& V ) { return this->change( ) -> replaceVertex( v , V ); }


    //! Create new edge with the origin v1 and target [and, maybe other parameters] provided in E
    void newEdge( int v1 , const EdgeType& E ) { this->change( ) -> newEdge( v1 , E ); }


    //! Make the graph trivial.
    void clear( ) { this->change( ) -> clear( ); }


    //! Pinch two vertices v1 and v2.
    void pinch( int v1 , int v2 ) { this->change( ) -> pinch( v1 , v2 ); }
      
  };


  //! Output the graph.
  /*! To use this function a function 
    ostream& operator << ( ostream& os , const EdgeType& E ) 
    must be defined. Otherwise you will get an error when linking.
   */
  template< class VertexType , class EdgeType > 
    ostream& operator << ( ostream& os , const GraphConcept< VertexType , EdgeType >& G ) {

    typedef VertexType                               vertex_type;
    typedef EdgeType                                 edge_type;
    const map< int, vertex_type >& theVertices = G.getVertices( );
    typename map< int, vertex_type >::const_iterator v_it = theVertices.begin( );
    for( ; v_it!=theVertices.end( ) ; ++v_it ) {

      int v = (*v_it).first;
      const vertex_type& V = (*v_it).second;
      const set< edge_type >& out = V.out;
      os << v << ":";
      typename set< edge_type >::const_iterator out_it = out.begin( );
      for( ; out_it!=out.end( ) ; ++out_it )
	os << *out_it << " ";
      os << endl;
    }

    return os;
  }

}


#endif
