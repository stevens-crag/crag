
#include "GraphDrawingAttributes.h"

map< int , string > GraphDrawingAttributes::nodeShapeNames = GraphDrawingAttributes::initializeNodeShapeNames( );


map< int , string > GraphDrawingAttributes::initializeNodeShapeNames( )
{
  map< int , string > result;
  result[box] = "box";
  result[polygon] = "polygon";
  result[ellipse] = "ellipse";
  result[circle] = "circle";
  result[point] = "point";
  result[egg] = "egg";
  result[triangle] = "triangle";
  result[plaintext] = "plaintext";
  result[diamond] = "diamond";
  result[trapezium] = "trapezium";
  result[parallelogram] = "parallelogram";
  result[house] = "house";
  result[pentagon] = "pentagon";
  result[hexagon] = "hexagon";
  result[septagon] = "septagon";
  result[octagon] = "octagon";
  result[doublecircle] = "doublecircle";
  result[doubleoctagon] = "doubleoctagon";
  result[invtriangle] = "invtriangle";
  result[invtrapezium] = "invtrapezium";
  result[invhouse] = "invhouse";
  result[Mdiamond] = "Mdiamond";
  result[Msquare] = "Msquare";
  result[Mcircle] = "Mcircle";
  result[rect] = "rect";
  result[rectangle] = "rectangle";
  result[none] = "none";
  
  return result;
}
