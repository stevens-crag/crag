// Copyright (C) 2005 Alexander Ushakov
// Contents: Implementation of class SubgroupFG
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//


#include <fstream>
#include <iomanip>
#include "iterator"
#include "strstream"
#include "SubgroupFG.h"
#include "GraphDrawingAttributes.h"


//---------------------------------------------------------------------------//
//------------------------------ SubgroupFG ---------------------------------//
//---------------------------------------------------------------------------//


SubgroupFG::SubgroupFG( int n_gens ) :
  theNumberOfGenerators( n_gens ),
  fsaDone( false ),
  nielsDone( false )
{

}


SubgroupFG::SubgroupFG( int n_gens , const vector< Word >& gens ) :
  theNumberOfGenerators( n_gens ),
  theGenerators( gens ),
  fsaDone( false ),
  nielsDone( false )
{

}


SubgroupFG::SubgroupFG( int n_gens , const IntLabeledGraph& fsa ) :
theNumberOfGenerators( n_gens )
{
  // 1. make fsa deterministic
  theFSA = fsa;
  fsaDone = true;
  
  // 2. find a Nielsen basis
  computeNielsenGenerators( );
  theGenerators = theNielsenGenerators;
}


//---------------------------------------------------------------------------//
//------------------------------ operator << --------------------------------//
//---------------------------------------------------------------------------//


ostream& operator << ( ostream& os , const SubgroupFG& sbgp )
{
  os << "[";
  const vector< Word >& gens = sbgp.getGenerators( );
  for( int i=0 ; i<gens.size( ) ; ++i ) {
    if( i ) os << ", ";
    os << gens[i];
  }
  os << "]";
  
  return os;
}


//---------------------------------------------------------------------------//
//------------------------------ operator^=  --------------------------------//
//---------------------------------------------------------------------------//


SubgroupFG& SubgroupFG::operator^= ( const Word& conjugator )
{
  for( vector< Word >::iterator g_it = theGenerators.begin( ) ; g_it!=theGenerators.end( ) ; ++g_it )
    *g_it ^= conjugator;
  
  fsaDone = false;
  nielsDone = false;
  
  theFSA.clear( );
  theNielsenGenerators.clear( );
  
  return *this;
}


//---------------------------------------------------------------------------//
//------------------------------ operator^=  --------------------------------//
//---------------------------------------------------------------------------//


SubgroupFG& SubgroupFG::operator+= ( const Word& w )
{
  theGenerators.push_back( w );
  
  fsaDone = false;
  nielsDone = false;
  
  theFSA.clear( );
  theNielsenGenerators.clear( );
  
  return *this;
}


SubgroupFG& SubgroupFG::operator+= ( const SubgroupFG& sbgp )
{
  for( vector< Word >::const_iterator g_it = sbgp.theGenerators.begin( ) ; g_it!=sbgp.theGenerators.end( ) ; ++g_it )
    theGenerators.push_back( *g_it );
  
  fsaDone = false;
  nielsDone = false;
  
  theFSA.clear( );
  theNielsenGenerators.clear( );
  
  return *this;
}


//---------------------------------------------------------------------------//
//-------------------------------- getFSA -----------------------------------//
//---------------------------------------------------------------------------//


const IntLabeledGraph& SubgroupFG::getFSA( ) const
{
  if( !fsaDone ) computeFSA( );
  return theFSA;
}


void SubgroupFG::computeFSA( ) const
{
  int is = theFSA.newVertex( );
  for( int i=0 ; i<theGenerators.size( ) ; ++i ) {
    const list< int >& lst = theGenerators[i].getList( );
    addLoop( theFSA , is , lst.begin( ) , lst.end( ) );
  }
  fold( theFSA , 0 , &foldDetails );
  fsaDone = true;
}


//---------------------------------------------------------------------------//
//------------------------------ doesBelong ---------------------------------//
//---------------------------------------------------------------------------//


bool SubgroupFG::doesBelong( const Word& w ) const
{
  if( !fsaDone ) computeFSA( );
  return trace( getFSA( ) , 0 , w.begin( ) , w.end( ) )==0;
}


//---------------------------------------------------------------------------//
//------------------------------- getIndex ----------------------------------//
//---------------------------------------------------------------------------//


int SubgroupFG::getIndex( ) const
{
  typedef IntLabeledGraph::edge_type edge_type;
  typedef IntLabeledGraph::vertex_type vertex_type;

  if( !fsaDone ) computeFSA( );
  const map< int , vertex_type >& theVertices = theFSA.getVertices( );
  
  for( map< int , vertex_type >::const_iterator v_it = theVertices.begin( ) ; v_it!=theVertices.end( ) ; ++v_it ) {
    const set< edge_type >& out = (*v_it).second.out;
    if( out.size( )!=2*theNumberOfGenerators )
      return -1;
    set< edge_type >::const_iterator out_it = out.begin( );
    for( int i=-theNumberOfGenerators ; out_it!=out.end( ) ; ++out_it, ++i ) {
      if( i==0 ) ++i;
      if( (*out_it).theLabel!=i )
	return -1;
    }	
  }
  
  return theVertices.size( );
}


//---------------------------------------------------------------------------//
//----------------------------- checkIsomorphism ----------------------------//
//---------------------------------------------------------------------------//


bool SubgroupFG::checkIsomorphism( const SubgroupFG& S , int vert1 , int vert2 ) const
{
  typedef IntLabeledGraph::edge_type edge_type;
  typedef IntLabeledGraph::vertex_type vertex_type;

  if( !fsaDone ) computeFSA( );
  if( !S.fsaDone ) S.computeFSA( );
  const map< int , vertex_type >& V1 =   theFSA.getVertices( );
  const map< int , vertex_type >& V2 = S.theFSA.getVertices( );
  if( V1.size( )!=V2.size( ) )
    return false;
  
  set< int > images;
  map< int , int > M;
  images.insert( vert2 );
  M[vert1] = vert2;
  
  typedef pair< int , int > PII;
  set< PII > vertices;
  vertices.insert( PII( 0 , vert1 ) );
  
  while( 1 ) {
    
    // cout << "A1" << endl;
    // process the first element
    PII p = *vertices.begin( );
    if( p.first==1 ) break;
    int n1 = p.second;
    vertices.erase( vertices.begin( ) );
    vertices.insert( PII( 1 , n1 ) );
    
    // cout << "A2" << endl;
    // get states and edges
    int n2 = M[n1];
    const vertex_type& v1 = (*V1.find( n1 )).second;
    const vertex_type& v2 = (*V2.find( n2 )).second;
    const set< edge_type >& out1 = v1.out;
    const set< edge_type >& out2 = v2.out;
    if( out1.size( )!=out2.size( ) )
      return false;
    
    // cout << "A3" << endl;
    // process leaving edges
    set< edge_type >::const_iterator out_it1 = out1.begin( );
    set< edge_type >::const_iterator out_it2 = out2.begin( );
    for( ; out_it1!=out1.end( ) ; ++out_it1, ++out_it2 ) {
      
      if( (*out_it1).theLabel!=(*out_it2).theLabel )
	return false;
      int t1 = (*out_it1).theTarget;
      int t2 = (*out_it2).theTarget;
      
      if( M.find(t1)==M.end() ) {
	// without this you'd a morphism
	if( images.find( t2 )!=images.end( ) )
	  return false;
	images.insert( t2 );
	M[t1] = t2;
      } else {
	if( M[t1]!=t2 )
	  return false;
      }
      
      if( vertices.find( PII(0,t1) )==vertices.end( ) && vertices.find( PII(1,t1) )==vertices.end( ) )
	vertices.insert( PII(0,t1) );
    }
    // cout << "A4" << endl;
    
  }
  
  return true;
}


//---------------------------------------------------------------------------//
//-------------------------------- operator== -------------------------------//
//---------------------------------------------------------------------------//


bool SubgroupFG::operator== ( const SubgroupFG& S ) const
{
  return checkIsomorphism( S , 0 , 0 );
}


//---------------------------------------------------------------------------//
//-------------------------------- operator* --------------------------------//
//---------------------------------------------------------------------------//


SubgroupFG SubgroupFG::operator* ( const SubgroupFG& S ) const
{
  typedef IntLabeledGraph::edge_type edge_type;
  typedef IntLabeledGraph::vertex_type vertex_type;
  
  if( !fsaDone ) computeFSA( );
  if( !S.fsaDone ) S.computeFSA( );
  const map< int , vertex_type >& V1 =   theFSA.getVertices( );
  const map< int , vertex_type >& V2 = S.theFSA.getVertices( );
  
  IntLabeledGraph fsa;
  typedef pair< int , int > PII;
  map< PII , int > vertices;
  vertices[PII(0,0)] = fsa.newVertex( );
  list< PII > to_check;
  to_check.push_back( PII(0,0) );
  
  while( !to_check.empty( ) ) {
    
    PII p = *to_check.begin( );
    to_check.erase( to_check.begin( ) );
    
    int n1 = p.first;
    int n2 = p.second;
    int n = vertices[PII(n1,n2)];
    const vertex_type& v1 = (*V1.find( n1 )).second;
    const vertex_type& v2 = (*V2.find( n2 )).second;
    const set< edge_type >& out1 = v1.out;
    const set< edge_type >& out2 = v2.out;
    set< edge_type >::const_iterator out_it1 = out1.begin( );
    set< edge_type >::const_iterator out_it2 = out2.begin( );
    
    while( out_it1!=out1.end( ) && out_it2!=out2.end( ) ) {
      int l1 = (*out_it1).theLabel;
      int l2 = (*out_it2).theLabel;
      if( l1<l2 ) { ++out_it1; continue; }
      if( l2<l1 ) { ++out_it2; continue; }
      
      int t1 = (*out_it1).theTarget;
      int t2 = (*out_it2).theTarget;
      int t;
      map< PII , int >::iterator vertices_it = vertices.find( PII(t1,t2) );
      if( vertices_it==vertices.end( ) ) {
	to_check.push_back( PII(t1,t2) );
	vertices[PII(t1,t2)] = t = fsa.newVertex( );
	// cout << "  " << t << " -> (" << t1 << "," << t2 << ")" << endl;
      } else 
	t = (*vertices_it).second;
      
      fsa.newEdge( n , edge_type( t ,  l1 ) );
      fsa.newEdge( t , edge_type( n , -l1 ) );
      ++out_it1;
      ++out_it2;
    }
    
  }
  
  return SubgroupFG( theNumberOfGenerators>S.theNumberOfGenerators ? theNumberOfGenerators : S.theNumberOfGenerators , fsa );
}


//---------------------------------------------------------------------------//
//-------------------------- computeNielsenGenerators -----------------------//
//---------------------------------------------------------------------------//


void SubgroupFG::computeNielsenGenerators( ) const
{
  typedef IntLabeledGraph::edge_type edge_type;
  typedef IntLabeledGraph::vertex_type vertex_type;
  
  if( !fsaDone ) computeFSA( );
  
  // compute geodesic subtree
  map< int , edge_type > GT = getGeodesicTree_in( theFSA , 0 );

  // go over all states
  const map< int , vertex_type >& theVertices = theFSA.getVertices( );
  map< int , vertex_type >::const_iterator v_it = theVertices.begin( );
  for( ; v_it!=theVertices.end( ) ; ++v_it ) {
    
    // process all vertices incoming into a state
    int cur_v = (*v_it).first;
    const vertex_type& cur_V = (*v_it).second;
    const edge_type tree_edge = GT[cur_v];
    const set< edge_type >& out = cur_V.out;
    set< edge_type >::const_iterator out_it = out.begin( );
    for( ; out_it!=out.end( ) ; ++out_it ) {

      if( tree_edge==*out_it || (*out_it).theLabel<0 )
	continue;
      
      int from_v = (*out_it).theTarget;
      const edge_type tree_edge2 = GT[from_v];
      if( tree_edge2==edge_type( cur_v , -(*out_it).theLabel ) )
	continue;
      
      // cout << "   ! " << from_state << " " << tree_edge2.label << " " << tree_edge2.target << endl;
      // cout << "   - " << cur_state.theState << " " << (*out_it).label << " " << from_state << endl;
      list< edge_type > L1 = readoffGeodesicTree( GT , cur_v ).second;
      list< edge_type > L2 = readoffGeodesicTree( GT , from_v ).second;
      L2.push_front( *out_it );
      
      Word w;
      list< edge_type >::iterator L_it = L1.begin( );
      for( ; L_it!=L1.end( ) ; ++L_it )
	w.push_front( -(*L_it).theLabel );
      L_it = L2.begin( );
      for( ; L_it!=L2.end( ) ; ++L_it )
	w.push_back( (*L_it).theLabel );
      theNielsenGenerators.push_back( w );
    }
  }
  
  nielsDone = true;
}


//---------------------------------------------------------------------------//
//--------------------------- getNielsenGenerators --------------------------//
//---------------------------------------------------------------------------//


const vector< Word >& SubgroupFG::getNielsenGenerators( ) const
{
  if( !nielsDone )
    computeNielsenGenerators( );
  
  return theNielsenGenerators;
}


//---------------------------------------------------------------------------//
//--------------------------------- getRank ---------------------------------//
//---------------------------------------------------------------------------//


int SubgroupFG::getRank( ) const
{
  if( !nielsDone )
    computeNielsenGenerators( );
  
  return theNielsenGenerators.size( );
}


//---------------------------------------------------------------------------//
//--------------------------------- express ---------------------------------//
//---------------------------------------------------------------------------//


Word SubgroupFG::express( const Word& w ) const
{
  typedef IntLabeledGraph::edge_type edge_type;
  typedef IntLabeledGraph::vertex_type vertex_type;
  
  Word result;
  
  IntLabeledGraph fsa = theFSA;
  pair< bool , list< edge_type > > path = trace_path( fsa , 0 , w.begin( ) , w.end( ) );
  if( !path.first || (path.second.size( ) && (*--path.second.end()).theTarget!=0 ) ) {
    cout << "Error" << endl;
    exit( 1 );
  }
  
  liftup( fsa , 0 , path.second , foldDetails.begin( ) , foldDetails.end( ) );
  map< Word , int > numbered_gens;
  vector< Word >::const_iterator g_it = theGenerators.begin( );
  for( int t=1 ; g_it!=theGenerators.end( ) ; ++g_it, ++t ) {
    numbered_gens[*g_it]  =  t;
    numbered_gens[-*g_it] = -t;
  }

  Word cur_gen;
  int cur_state = 0;
  list< edge_type >::iterator path_it = path.second.begin( );
  for( ; path_it!=path.second.end( ) ; ++path_it ) {
    cur_gen.push_back( (*path_it).theLabel );
    if( (*path_it).theTarget==0 ) {
      int n = numbered_gens[cur_gen];
      if( n )	result.push_back( n );
      cur_gen = Word( );
    }
  }
  
  return result;
}


//---------------------------------------------------------------------------//
//------------------------------- centralizer -------------------------------//
//---------------------------------------------------------------------------//


SubgroupFG SubgroupFG::centralizer( ) const
{
  vector< Word > result;
  int rank = getRank( );
  if( rank==0 )
    for( int i=0 ; i<theNumberOfGenerators ; ++i )
      result.push_back( Word(i+1) );
  if( rank!=1 )
    return SubgroupFG( theNumberOfGenerators , result );
  
  if( !nielsDone )
    computeNielsenGenerators( );
  
  Word base;
  theNielsenGenerators[0].getPower( base );
  result.push_back( base );
  return SubgroupFG( theNumberOfGenerators , result );
}


//---------------------------------------------------------------------------//
//------------------------------- centralizer -------------------------------//
//---------------------------------------------------------------------------//


SubgroupFG SubgroupFG::normalizer( ) const
{
  SubgroupFG result;
  
  if( !fsaDone ) computeFSA( );
  
  typedef IntLabeledGraph::edge_type edge_type;
  typedef IntLabeledGraph::vertex_type vertex_type;

  // compute geodesic subtree
  map< int , edge_type > GT = getGeodesicTree_in( theFSA , 0 );
  
  const map< int , vertex_type >& theVertices = theFSA.getVertices( );
  map< int , vertex_type >::const_iterator v_it = theVertices.begin( );
  for( ; v_it!=theVertices.end( ) ; ++v_it ) {

    if( checkIsomorphism( *this , 0 , (*v_it).first ) ) {
      list< edge_type > L1 = readoffGeodesicTree( GT , (*v_it).first ).second;
      Word w;
      list< edge_type >::iterator L_it = L1.begin( );
      for( ; L_it!=L1.end( ) ; ++L_it )
	w.push_front( -(*L_it).theLabel );
      
      result += *this^w;
    }
  }
  
  result.computeNielsenGenerators( );
  result.theGenerators = result.theNielsenGenerators;
  return result;
}


//---------------------------------------------------------------------------//
//--------------------------------- getRank ---------------------------------//
//---------------------------------------------------------------------------//


pair< bool , Word > SubgroupFG::areConjugate( const SubgroupFG& sbgp ) const
{
  typedef IntLabeledGraph::edge_type edge_type;
  typedef IntLabeledGraph::vertex_type vertex_type;
  
  pair< SubgroupFG , Word > pr1 = this->trim( );
  pair< SubgroupFG , Word > pr2 = sbgp. trim( );
  const SubgroupFG& Sbgp1 = pr1.first;
  const SubgroupFG& Sbgp2 = pr2.first;

  if( !Sbgp1.fsaDone ) Sbgp1.computeFSA( );
  if( !Sbgp2.fsaDone ) Sbgp2.computeFSA( );
  const map< int , vertex_type >& V1 = Sbgp1.theFSA.getVertices( );
  const map< int , vertex_type >& V2 = Sbgp2.theFSA.getVertices( );
  if( V1.size( )!=V2.size( ) )
    return pair< bool , Word >( false , Word( ) );
  
  // compute geodesic subtree
  map< int , edge_type > GT = getGeodesicTree_in( Sbgp2.theFSA , 0 );
  
  map< int , vertex_type >::const_iterator v_it = V2.begin( );
  for( ; v_it!=V2.end( ) ; ++v_it ) {
    
    if( !Sbgp1.checkIsomorphism( Sbgp2 , 0 , (*v_it).first ) )
      continue;
    
    list< edge_type > L1 = readoffGeodesicTree( GT , (*v_it).first ).second;
    Word w;
    list< edge_type >::iterator L_it = L1.begin( );
    for( ; L_it!=L1.end( ) ; ++L_it )
      w.push_front( -(*L_it).theLabel );
    return pair< bool , Word >( true , pr1.second * -w * -pr2.second );
  }
  return pair< bool , Word >( false , Word( ) );
}


//---------------------------------------------------------------------------//
//---------------------------------- trim -----------------------------------//
//---------------------------------------------------------------------------//


pair< SubgroupFG , Word > SubgroupFG::trim( ) const
{
  typedef IntLabeledGraph::edge_type edge_type;
  typedef IntLabeledGraph::vertex_type vertex_type;
  
  if( !fsaDone ) computeFSA( );
  
  int prev_v = 0;
  int cur_v = 0;
  Word tail;
  const map< int , vertex_type >& theVertices = theFSA.getVertices( );
  for( int i=0 ; ; ++i ) {
    
    vertex_type V = (*theVertices.find(cur_v)).second;
    const set< edge_type >& out = V.out;
    if( i==0 && out.size( )==1 ) {
      const edge_type& E = *out.begin( );
      prev_v = cur_v;
      cur_v = E.theTarget;
      tail.push_back( E.theLabel );
      continue;
    }
    if( i>0 && out.size( )==2 ) {
      const edge_type& E1 = *out.begin( );
      const edge_type& E2 = *(++out.begin( ));
      if( E1.theTarget!=prev_v ) {
	prev_v = cur_v;
	cur_v = E1.theTarget;
	tail.push_back( E1.theLabel );
      } else {
	prev_v = cur_v;
	cur_v = E2.theTarget;
	tail.push_back( E2.theLabel );
      }
      continue;
    }
    
    break;  
  }
  
  return pair< SubgroupFG , Word >( *this^tail , tail );
}


//---------------------------------------------------------------------------//
//--------------------------------- getRank ---------------------------------//
//---------------------------------------------------------------------------//


string SubgroupFG::graphviz_format( ) const
{
  typedef GraphDrawingAttributes::COLOR COLOR;
  typedef GraphDrawingAttributes::NODESHAPE NODESHAPE;
  const map< int , string >& nodeShapeNames = GraphDrawingAttributes::getNodeShapeNames( );
  
  typedef IntLabeledGraph::vertex_type vertex_type;
  typedef IntLabeledGraph::edge_type edge_type;
  const map< int, vertex_type >& theVertices = theFSA.getVertices( );

  if( !fsaDone ) computeFSA( );

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
      if( (*out_it).theLabel>0 )
	ostr << v << " -> " << (*out_it).theTarget << " [label=\"" << (*out_it).theLabel << "\"];";
    }
  }
  ostr << "}";
  
  return ostr.str( );
}


//---------------------------------------------------------------------------//
//--------------------------------- getRank ---------------------------------//
//---------------------------------------------------------------------------//


Word substitute( const Word& w , const vector< Word >& wrds ) 
{
  return substitute( w.begin() , w.end() , wrds );
}

