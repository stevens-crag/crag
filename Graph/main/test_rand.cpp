// Copyright (C) 2007 Alexander Ushakov
// Contents: Test of "Random FSA Generation" methods.
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//

#include "string"
#include <iomanip>
#include "strstream"
#include "fstream"
using namespace std;

#include "RandomFSA.h"
#include "GraphDrawingAttributes.h"


string graphviz_format( const FSA& theFSA )
{
  typedef GraphDrawingAttributes::COLOR COLOR;
  typedef GraphDrawingAttributes::NODESHAPE NODESHAPE;
  const map< int , string >& nodeShapeNames = GraphDrawingAttributes::getNodeShapeNames( );
  
  typedef FSA::state_type vertex_type;
  typedef FSA::edge_type edge_type;
  const map< int, vertex_type >& theVertices = theFSA.getStates( );


  GraphDrawingAttributes GDA;
  GDA.setNodeColor( 0 , COLOR(255,0,0) );
  GDA.setNodeShape( 0 , GraphDrawingAttributes::doublecircle );
  
  ostrstream ostr;
  string vertices;
  ostr << "digraph G { ";
  for( map< int, vertex_type >::const_iterator v_it = theVertices.begin( ) ; v_it!=theVertices.end( ) ; ++v_it ) {
    
    COLOR c = GDA.getNodeColor( (*v_it).first );
    NODESHAPE s = GDA.getNodeShape( (*v_it).first );
    string shape = nodeShapeNames.at(s);
    
    // Vertex:
    ostr << setbase(10) << (*v_it).first;
    // Attrubutes of vertex: Shape
    ostr << "[shape=" << shape;

    ostr << ",style=filled,color=\"#";
    ostr << setbase(16);
    // Attrubutes of vertex: Color
    if( c. first<16 ) ostr << "0";
    ostr << c.first;
    if( c.second<16 ) ostr << "0";
    ostr << c.second;
    if( c. third<16 ) ostr << "0";
    ostr << c.third;
    ostr << "\"]; ";
  }

  for( map< int, vertex_type >::const_iterator v_it = theVertices.begin( ) ; v_it!=theVertices.end( ) ; ++v_it ) {
    ostr << setbase(10);
    int v = (*v_it).first;
    const vertex_type& V = (*v_it).second;
    const set< edge_type >& out = V.out;
    set< edge_type >::const_iterator out_it = out.begin( );
    for( ; out_it!=out.end( ) ; ++out_it ) {
      if( (*out_it).label>0 )
	ostr << v << " -> " << (*out_it).target << " [label=\"" << (*out_it).label << "\"];";
    }
  }
  ostr << "}";
  
  return ostr.str( );
}


//---------------------------------------------------------------------------//
//-------------------------------- main -------------------------------------//
//---------------------------------------------------------------------------//


int main( )
{
  int N = 50;
  int L = 2;
  FSA fsa = randomFSA( N , L );
  
  cout << "======================" << endl;
  cout << fsa << endl;
  cout << "======================" << endl;

  ofstream of( "g.txt" );
  of << graphviz_format( fsa ) << endl;


  return 0;
}
