// Copyright (C) 2005 Alexei Miasnikov
// Contents: Definition of classes which implement different versions
//           of Whitehead graphs
//
// Principal Author: Alexei Miasnikov (2005)
//
// Status: Useable.
//



#ifndef _WGRAPH_H_
#define _WGRAPH_H_


#include "Word.h"
#include <vector>
#include "errormsgs.h"

typedef int Generator;

////////////////////////////////////////////////////////////
//
//  Class to compute all cut vertexes of a graph  ( there is an edge (i,j) in the graph iff G[i][j]\f$\neq$\f  0).
//
////////////////////////////////////////////////////////////

//! Implements an algorithm to find all articulation points of a graph//
class CutVertices
{
 public:
   //! Constructor
  /*!
    \param G - the adjacency matrix of a graph ( there is an edge from i to j in the graph iff G[i][j] \f$\neq\f$ 0  ). 
    \param n - the number of vertices.
  */
  CutVertices( const int** G, int n);
  ~CutVertices();

  //! Finds articulation  points.
  /*!
    Finds all articulation points of the  graph. The result can be obtained through \link  getCutVertices() getCutVertices()\endlink function.
    This algorithm  requires linear time in the size of the graph, or \f$O(n^2)\f$. 
   */
  void compute();

  
  //! Finds articulation  points (Brute force algorithm.
  /*!
    Finds all articulation points of the  graph. The result can be obtained through \link  getCutVertices() getCutVertices()\endlink function.
    Complexity: \f$O(n^2)\f$. 
   */  
  void computeBruteForce();

  
  
  //! Computes the number of connected components of the graph.
  /*!
    Computes the number of connected components of the graph.
    \param ignore_single_vertices - if set to \c true, then components containing only one vertex are ignored. Default value is \c false.
    \return The number of connected components
   */  
  int numberOfComponents( bool ignore_single_vertices = false );

  
  //! Get the list of articulation points.
  /*!
    Return the list of articulation points obtained using one of the algorithms. \link compute() compute() \endlink or
    \link computeBruteForce() computeBruteForce() \endlink must be called first.
    \return The list of articulation vertices.
   */ 
  vector<int> getCutVertices()const {return articulation_point;}
  
 private:
  void init();

  
   //! Execute depth first search.
  void DepthFirstSearch();

  //! Execute recursive depth first search.
  void RecursiveDepthFirstSearch(int v);

  //! Compute the low subtree recursively.
  void RDFS_Compute_Low(int v);

  //! Extract articulation points from the labels
  void ArticulationPoints();

  //! Adjacency matrix of the graph
  int** theGraph;

  //! Number of vertices
  int N;

  //! Current time label (internal use only)
  int time;
  //! List of visited points (internal use only)
  vector<int> visit;
  //! List of predecessors (internal use only)
  vector<int> pred;
  //! Labels of discovered vertices (internal use only)
  vector<int>  discover;
  vector<int>  Low;

  
  //! List of articulation points
  vector<int> articulation_point;

};


////////////////////////////////////////////////////////////
//
//  Base (interface)  class 
//
////////////////////////////////////////////////////////////


//! Interface class for Whitehead graphs
class WhiteheadGraph 
{
 public:

  //! Constructor
  /*!
    \param w - the input word.
    \param num_of_gens - the number of generators in the corresponding group (i.e. \f$ w \in F_{num_of_gens}\f$.
  */
  WhiteheadGraph( const Word& w, int num_of_gens ): 
    theWord( w ), 
    nOfGenerators(num_of_gens  ) {}

  //! Get the input word
  /*
    \return Reference to the input word. 
   */
 const Word& getWord()const { return theWord; }
 protected:
 //! The input word 
 Word theWord;
 //! The number of generators in the group
 int nOfGenerators;
};


////////////////////////////////////////////////////////////
//
//  Implementation of Whitehead Simple Graph (i.e. no multiple
//                    edges or loops allowed)  
//
////////////////////////////////////////////////////////////

//! Implements  Whitehead Simple graph (i.e. no multiple edges or loops allowed)  
class WhiteheadSimpleGraph : public WhiteheadGraph
{
 public:

  //! Constructor
  /*!
    \param w - input word.
    \param num_of_gen - the number of generators in the corresponding group
  */
  WhiteheadSimpleGraph( const Word& w, int num_of_gens );
  
  //! Output operator
  friend ostream& operator << ( ostream& out, const WhiteheadSimpleGraph& g){
    g.printOn( out );
    return out;
  }
  
  ~WhiteheadSimpleGraph();
    
  //! Get the label of the edge (i,j). 
  /*! The label is the number of times (counts) the subword \f$ i j^{-1}\f$  occurred in w.
     \param i - first vertex (letter)
     \param j - second vertex (letter)
     \return Label of the edge (i,j)
  */
  int getCount(int i, int j)const {
    if ( i<0 || j<0 || i>=theSize || j>=theSize )
      msgs::error("WhiteheadSimpleGraph::getCount(...): index is out of bounds.");
    return theAdjMatrix[i][j];
  }
  
  //! Get the number of vertices of the graph
  /*!
  \return The number of vertices
  */
  int getSize()const { return theSize; }
  //! Get a vector containing weights of the edges of the  graph
  /*!
    \return Vector of the weights of the edges of the graph. 
    NOTE: the graph is treated as undirected Whitehead graph when weights are computed.
    See \link makeUndirected() makeUndirected() \endlink for more information.
  */
  vector<double> getWeightVector() const;
  //! Get a vector of words, corresponding to edges of the graph
  /*!
    \return The vector of words, corresponding to the edges of the Whitehead graph.
  */
  vector<Word> getWeightNames()const;

  //! Check if undirected graph.
  /*!
    By default the graph is undirected Whitehead graph. Run \link makeUndirected() makeUndirected() \endlink 
    to convert to undirected graph. 
    \return \c true if the graph is undirected Whitehead graph, \c false otherwise.
  */
  bool isUndirected()const { return undirected; }

  //! Convert to undirected Whitehead graph.
  /*!
  Initially a Directed Whitehead graph is constructed. Use this function to convert to
  undirected Whitehead graph. 

  Edges are added, and the new edge weights \f$ \omega^*(i,j) = \omega^*(j,i) =  \omega(i,j) + \omega(j,i)\f$, 
  where \f$\omega(m,n)\f$ are the weights in the directed graph.
  */
  void makeUndirected();
  
  //! Get the number of connected components
  /*!
    \return The number of connected components.
   */
  int numberOfComponents() const;
  //! Get the  number of cut vertices (articulation points)
  /*!
    \return The number of cut vertices.
   */  
  int numberOfCutVertices()const;
  //! Get the list  of cut vertices (articulation points)
  /*!
    \return The list of cut vertices.
  */    
  vector<int> cutVertices()const;  
  
  int numberOfCutVerticesBruteForce()const;
  vector<int> cutVerticesBruteForce()const;

 private:
  //! Outputs the graph into a stream
  void printOn( ostream& out ) const;
  
  //! Transfers a letter into a graph vertex
  int genToIndex( const Generator& g )const{
    int ind = g;
    if ( ind > 0) 
      return ind-1;
    else 
      return -ind + nOfGenerators - 1;
  }
  
  //! Transfers a vertex into the corresponding letter
  int indToGenerator( int i )const{
    if (i < nOfGenerators)
      // if positive
      return i+1;
    else {
      // if negative
      return -( i - nOfGenerators + 1 );
    }
  }

  //! computes the number of connected components
  int nOfComponents() const;

  //! The number of vertices
  int theSize;
  
  //! The adjacensy matrix of the graph
  int** theAdjMatrix;

  //! True if undirected
  bool undirected;
};


////////////////////////////////////////////////////////////
//
//  Implementation of Whitehead Multi-Graph
//
////////////////////////////////////////////////////////////

//! Whitehead Multi-Graph (Not implemented yet)
class WhiteheadMultiGraph : public WhiteheadGraph
{
 public:
  WhiteheadMultiGraph( const Word& w, int n) : WhiteheadGraph( w,n ) {} 
 private:
};

#endif
