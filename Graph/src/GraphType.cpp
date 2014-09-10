
#include "GraphType.h"
#include <strstream>
#include <iomanip>
using namespace std;


ostream& operator << ( ostream& os , const GraphEdge& e ) 
{
  os << e.theTarget;
  return os;
}

ostream& operator << ( ostream& os , const IntLabeledEdge& e )
{
  os << "(" << e.theTarget << "," << e.theLabel << ")";
  return os;
}

void prepareToFold( const IntLabeledEdge& e1 , const IntLabeledEdge& e2 )
{
  return;
}


//---------------------------------------------------------------------------//
//---------------------------------- Graph ----------------------------------//
//---------------------------------------------------------------------------//


Graph randomGraph( int N , float edge_param ) 
{
  Graph result;

  for( int i=0 ; i<N ; ++i )
    result.newVertex( );

  for( int i=0 ; i<N ; ++i ) {
    for( int j=i+1 ; j<N ; ++j ) {
      if( RandLib::ur.rand()<edge_param ) {
	result.newEdge( i , j );
	result.newEdge( j , i );
      }
    }
  }
  
  return result;
}


//---------------------------------------------------------------------------//
//---------------------------------- Graph ----------------------------------//
//---------------------------------------------------------------------------//


vector< vector< int > > lengthTable( const Graph& G )
{
  typedef Graph::vertex_type vertex_type;

  // Get the adjacency structure
  const map< int, vertex_type >& theVertices = G.getVertices( );

  // Get the number of states
  int N = theVertices.size( );

  // Prepare enumeration of vertices
  map< int , int > vertexNumber;
  map< int, vertex_type >::const_iterator v_it = theVertices.begin( );
  for( int i=0 ; v_it!=theVertices.end( ) ; ++i , ++v_it )
    vertexNumber[(*v_it).first] = i;
  
  // Compute the result
  vector< vector< int > > result( N , vector< int >( ) );
  v_it = theVertices.begin( );
  for( int i=0 ; v_it!=theVertices.end( ) ; ++i , ++v_it ) {
    result[i] = vector< int >( N , -1 );
    map< int , int > dist = getDistances_out( G , (*v_it).first );
    map< int , int >::iterator d_it = dist.begin( );
    for( ; d_it!=dist.end( ) ; ++d_it )
      result[i][vertexNumber[(*d_it).first]] = (*d_it).second;
  }
  
  return result;
}


//---------------------------------------------------------------------------//
//---------------------------------- Graph ----------------------------------//
//---------------------------------------------------------------------------//


vector< vector< int > > innerProductTable( const Graph& G , int origin )
{
  int N = G.getVertices( ).size( );
  vector< vector< int > > dist = lengthTable( G );
  vector< vector< int > > result( N , vector< int >( ) );
  for( int i=0 ; i<N ; ++i ) {
    result[i] = vector< int >( N , -1 );
    for( int j=0 ; j<N ; ++j ) {
      if( dist[origin][i]==-1 || dist[origin][j]==-1 || dist[i][j]==-1 )
	continue;
      result[i][j] = dist[origin][i] + dist[origin][j] - dist[i][j];
    }
  }
  
  return result;
}


float getHyperbolicityConst( const Graph& G )
{
  // Get the adjacency structure
  typedef Graph::vertex_type vertex_type;
  const map< int, vertex_type >& theVertices = G.getVertices( );
  if( theVertices.empty( ) )
    return 0;
  
  float result = 0;
  int N = theVertices.size( );
  int origin = (*theVertices.begin()).first;
  vector< vector< int > > IP = innerProductTable( G , origin );
  for( int a=0 ; a<N ; ++a ) {
    for( int b=0 ; b<N ; ++b ) {
      for( int c=0 ; c<N ; ++c ) {
	if( IP[a][b]==-1 || IP[a][c]==-1 || IP[b][c]==-1 )
	  continue;
	int val = IP[a][c] < IP[b][c] ? IP[a][c] : IP[b][c];
	val -= IP[a][b];
	if( result<val )
	  result = val;
      }
    }
  }
  
  return result/2;
}


//---------------------------------------------------------------------------//
//----------------------------- graphviz_format -----------------------------//
//---------------------------------------------------------------------------//


string graphviz_format( const Graph& G , const GraphDrawingAttributes& GDA )
{
  typedef GraphDrawingAttributes::COLOR COLOR;
  typedef GraphDrawingAttributes::NODESHAPE NODESHAPE;
  const map< int , string >& nodeShapeNames = GraphDrawingAttributes::getNodeShapeNames( );

  typedef Graph::vertex_type vertex_type;
  typedef Graph::edge_type edge_type;
  const map< int, vertex_type >& theVertices = G.getVertices( );

  ostrstream ostr;
  
  string vertices;
  ostr << "digraph G {  node [shape=circle,style=filled,color=\".7 .3 1.0\"]; ";
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
    ostr << setbase(10);
    // ostr << (*v_it).first << ";";
  }

  for( map< int, vertex_type >::const_iterator v_it = theVertices.begin( ) ; v_it!=theVertices.end( ) ; ++v_it ) {
    int v = (*v_it).first;
    const vertex_type& V = (*v_it).second;
    const set< edge_type >& out = V.out;
    set< edge_type >::const_iterator out_it = out.begin( );
    for( ; out_it!=out.end( ) ; ++out_it )
      ostr << v << " -> " << (*out_it).theTarget << ";";
  }

  ostr << "}";
  
  return ostr.str( );
}

//---------------------------------------------------------------------------//
//----------------------------- graphviz_format -----------------------------//
//---------------------------------------------------------------------------//


string graphviz_format( const IntLabeledGraph& G , const GraphDrawingAttributes& GDA )
{
  typedef GraphDrawingAttributes::COLOR COLOR;
  typedef GraphDrawingAttributes::NODESHAPE NODESHAPE;
  const map< int , string >& nodeShapeNames = GraphDrawingAttributes::getNodeShapeNames( );
  
  typedef IntLabeledGraph::vertex_type vertex_type;
  typedef IntLabeledGraph::edge_type edge_type;
  const map< int, vertex_type >& theVertices = G.getVertices( );

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
      if( v<=(*out_it).theTarget )
	ostr << v << " -> " << (*out_it).theTarget << " [label=\"" << (*out_it).theLabel << "\"];";
    }
  }

  ostr << "}";
  
  return ostr.str( );
}
