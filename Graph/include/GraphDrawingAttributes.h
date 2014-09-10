// Copyright (C) 2005 Alexander Ushakov
// Contents: Definition of the graph concept
//
// Principal Authors: Alexander Ushakov
//


#ifndef _GraphDrawingAttributes_H_
#define _GraphDrawingAttributes_H_


#include "tuples.h"
#include "map"
#include "set"
using namespace std;

//---------------------------------------------------------------------------//
//-------------------------- GraphDrawingAttributes -------------------------//
//---------------------------------------------------------------------------//


class GraphDrawingAttributes
{
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Internal types:                                    //
  //                                                     //
  /////////////////////////////////////////////////////////

 public:

  typedef triple< int , int , int > COLOR;
  enum NODESHAPE { 
    box , polygon  , ellipse , circle , point , egg , triangle , plaintext , diamond , trapezium ,
    parallelogram , house  , pentagon , hexagon , septagon , octagon , doublecircle , doubleoctagon ,
    tripleoctagon , invtriangle , invtrapezium , invhouse , Mdiamond , Msquare , Mcircle , rect ,
    rectangle , none };
  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Constructors                                       //
  //                                                     //
  /////////////////////////////////////////////////////////
 public:
  
  GraphDrawingAttributes( ) : 
    theDefaultNodeColor( 0 , 120 , 120 ),
    theDefaultNodeShape( circle ) { }
  
  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Accessors                                          //
  //                                                     //
  /////////////////////////////////////////////////////////
 public:

  //! Set the color for a node.
  void setNodeColor( int v , const COLOR& c ) { nodeColor[v] = c; }
  //! Get the color of a node.
  COLOR getNodeColor( int v ) const { 
    if( nodeColor.find(v)==nodeColor.end( ) ) 
      return theDefaultNodeColor;
    return nodeColor.at(v);
  }

  //! Set the shape for a node.
  void setNodeShape( int v , const NODESHAPE& s ) { nodeShape[v] = s; }
  //! Get the shape of a node.
  NODESHAPE getNodeShape( int v ) const {
    if( nodeShape.find(v)==nodeShape.end( ) ) 
      return theDefaultNodeShape;
    return nodeShape.at(v);
  }

  void setDefaultNodeShape( const NODESHAPE& s ) { theDefaultNodeShape = s; }
  
  static const map< int , string >& getNodeShapeNames( ) { return nodeShapeNames; }
  
  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Internal functions                                 //
  //                                                     //
  /////////////////////////////////////////////////////////
 private:
  
  static map< int , string > initializeNodeShapeNames( );
  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Data members                                       //
  //                                                     //
  /////////////////////////////////////////////////////////
 private:

  // The default color for nodes (of nodes that are not out of "nodeColor").
  COLOR theDefaultNodeColor;
  // Colors of nodes.
  map< int , COLOR > nodeColor;

  static map< int , string > nodeShapeNames;
  // The default shape for nodes (of nodes that are not out of "nodeShape").
  NODESHAPE theDefaultNodeShape;
  // Shapes of nodes.
  map< int , NODESHAPE > nodeShape;
};


#endif

