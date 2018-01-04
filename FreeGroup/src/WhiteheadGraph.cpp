/*
 *   $Id: WhiteheadGraph.cpp,v 1.1 2005/12/20 15:46:22 amiasnik Exp $
 */

// Contents: Implementation of classes which implement different versions
//           of Whitehead graphs
//
// Principal Author: Alexei Miasnikov (2003)
//
// Status: Useable.
//


#include "WhiteheadGraph.h"

//#include <unistd.h>


/*********************************************************************************/
//
//    CutVertices
//
/*********************************************************************************/
CutVertices::CutVertices(const int** G, int n):
  N(n),
  time(0),
  visit(N),
  pred(N),
  discover(N),
  Low(N),
  articulation_point(N)
{
  theGraph = new int*[n];
  for ( int i=0;i<n;i++){
    theGraph[i] =  new int[n];
    // copy the graph
    for ( int j=0;j<n;j++)
      theGraph[i][j]  =  G[i][j];    
  }
}


CutVertices::~CutVertices()
{
  for (int i=0;i<N;i++)
    delete [] theGraph[i];
  
  delete [] theGraph;
}

void  CutVertices::init()
{
 // initialize
  time = 0;
  
  for (int i=0;i<N;i++) {
    visit[i] = 0;
    pred[i] = -1;
    discover[i] = 0;
    Low[i] = 0;
  }
}

void CutVertices::compute()
{
  // initialize
  init();

  for (int i=0;i<N;i++)
    articulation_point[i] = 0;

  for (int v=0;v<N;v++)
    if (visit[v] == 0)  RDFS_Compute_Low(v);
  
  ArticulationPoints();


}

void CutVertices::computeBruteForce()
{

  init();
  
 for (int i=0;i<N;i++)
    articulation_point[i] = 0;

 int nc = numberOfComponents(  );
  
  for (int V=0;V<N;V++){
    // remove the vertex from the graph
    vector<int> save_column(N);
    vector<int> save_row(N);
    for (int i=0;i<N;i++){
      save_column[i] = theGraph[i][V];
      save_row[i] = theGraph[V][i];
      theGraph[i][V] = theGraph[V][i] = 0;
    }
    
    // if after a removal we got more connected componenets, then cut vertex
    int nnc =  numberOfComponents( );

    //cout << nc << " ? " << nnc << endl;
    if ( nc < (nnc - 1) )
      articulation_point[V]=1;
    
    // restore the graph
    for (int i=0;i<N;i++){
      theGraph[i][V] = save_column[i];
      theGraph[V][i] = save_row[i];
    }
    
  }
  
}

int CutVertices::numberOfComponents( bool ignore_single_vertices )
{
  // initialize
  init();
  
  int nc = 0;
  for (int v=0;v<N;v++)
    if (visit[v] == 0) {
      RecursiveDepthFirstSearch(v);
       nc++;
    }

  if (ignore_single_vertices) {
    for (int v=0;v<N;v++){
     bool has_edge = false;
     for (int w=0;w<N;w++)
	if ( theGraph[v][w] || theGraph[w][v] )
	  has_edge = true;
      
      if (!has_edge) nc--;
    }
  }
  return nc;
  
}

void CutVertices::DepthFirstSearch()
{
  for (int v=0;v<N;v++)
    if (visit[v] == 0) RecursiveDepthFirstSearch(v);
}


void CutVertices::RecursiveDepthFirstSearch(int v) 
{
  visit[v] = 1;
  discover[v] = ++time;
  for (int w=0;w<N;w++)
    if ( theGraph[v][w] ) { // if connected
      if (visit[w] == false) {
	pred[w]=v;
	RecursiveDepthFirstSearch(w);
      }
    }
}

template<typename T>
T min(const T& a, const T& b) {
  return (a < b) ? a : b;
}

void CutVertices::RDFS_Compute_Low(int v) 
{
  visit[v]=1;
  Low[v]=discover[v] = ++time;
  for (int w=0;w<N;w++) {
    if (theGraph[v][w]){
      if (visit[w] == 0) {
	pred[w]=v;
	RDFS_Compute_Low(w);
	// When RDFS_Compute_Low(w) returns,
	// Low[w] stores the
	// lowest value it can climb up
	// for a subtree rooted at w.
	Low[v] = std::min(Low[v], Low[w]);
      } else if (w != pred[v]) {
	// {v, w} is a back edge
	Low[v] = std::min(Low[v], discover[w]);
      }
    }
  }
}



void CutVertices::ArticulationPoints()
{
  for (int v=0;v<N;v++){
    if (pred[v] == -1) { //v is a root

      // compute the number of children in the TREE
      int nChildren = 0;
      for (int w=0;w<N;w++)
	if ( (pred[w] == v) && theGraph[v][w] ) nChildren++;
      
      if (nChildren > 1) // if more then one then we hava a cut vertex
	articulation_point[v]=1;
      
    } else{

      for (int  w=0;w<N;w++) {
	if ( (pred[w] == v) && (theGraph[v][w]) ){
	  if (Low[w] >= discover[v])
	    articulation_point[v]=1;
	}
      }

    }
  }
}


/*********************************************************************************/
//
//    WhiteheadSimpleGraph
//
/*********************************************************************************/

WhiteheadSimpleGraph::WhiteheadSimpleGraph( const Word& w, int num_of_gens ): WhiteheadGraph( w,num_of_gens), undirected( false )
{
  theSize = 2*num_of_gens;

  // allocate space for the adjacency matrix
  theAdjMatrix = new int*[theSize];
  for ( int i=0;i<theSize;i++){
    theAdjMatrix[i] =  new int[theSize];
    // initialize with zeros
    for ( int j=0;j<theSize;j++)
      theAdjMatrix[i][j]  =  0;    
  }
  
  // fill out the edge information
  // trace throught the word and add edge (g1,g2^-1)
  // where g1, g2 are two letters subsequently occuring in
  // the word
  //for (int a=1;a<w.length();a++)
  //    theAdjMatrix[genToIndex(w[a-1])][genToIndex(-(w[a]))] += 1;
  for ( ConstWordIterator I = ++(w.begin()); I!=w.end();I++){
    ConstWordIterator I1 = I;
    Generator g1 = *(--I1);
    Generator g2 = *I;
    theAdjMatrix[genToIndex(g1)][genToIndex(-(g2))] += 1;
  }

  // connect the last with the first
  if ( w.length() > 0 )
    theAdjMatrix[genToIndex(*(--(w.end())))][genToIndex(-(*(w.begin())))] += 1;

}

WhiteheadSimpleGraph::~WhiteheadSimpleGraph()
{
  for (int i=0;i<theSize;i++)
    delete [] theAdjMatrix[i];
  
  delete theAdjMatrix;
}

vector<double> WhiteheadSimpleGraph::getWeightVector( ) const
{
  
  int s = getSize();
  vector<double> theVector( (2*nOfGenerators*(2*nOfGenerators-1))/2 , 0.0 );
  int count = 0;
  for (int i=0;i<s;i++)
    for (int j=i;j<s;j++)
      if (i!=j){
	//cout << i << "," << j << ":" << wg.getCount( i,j ) << "," << wg.getCount( j,i ) << endl;
	if (undirected)
	  theVector[count] = double(getCount( i,j ));
	else 
	  theVector[count] = double(getCount( i,j ) + getCount( j,i ));
	count++;
      }
  return theVector; 
}

vector<Word> WhiteheadSimpleGraph::getWeightNames()const
{
  
  int s = getSize();
  vector<Word> theVector( (2*nOfGenerators*(2*nOfGenerators-1))/2 );
  int count = 0;
  for (int i=0;i<s;i++)
    for (int j=i;j<s;j++)
      if (i!=j){
	//cout << i << "," << j << ":" << wg.getCount( i,j ) << "," << wg.getCount( j,i ) << endl;
	theVector[count] = indToGenerator(i) * indToGenerator(j);
	count++;
      }
  return theVector; 
}

 void WhiteheadSimpleGraph::makeUndirected() {
    if ( !undirected ){
      undirected = true;
      for (int i=0;i<theSize;i++)
	for (int j=i;j<theSize;j++){
	  theAdjMatrix[i][j] += theAdjMatrix[j][i];
	  theAdjMatrix[j][i] =  theAdjMatrix[i][j];
	}
    }
  }
  
  int WhiteheadSimpleGraph::numberOfComponents() const {
    if (!isUndirected()) msgs::error("Cut vertices for undirected WG only.");
    
    CutVertices cv( (const int**)theAdjMatrix, theSize );
    return cv.numberOfComponents();
  }
  
  int WhiteheadSimpleGraph::numberOfCutVertices()const {
    vector<int> vertices = cutVertices();
    
    int nv = 0;
    for (int i=0;i<theSize;i++)
      if (vertices[i] > 0)
	nv++;
    
    return nv;
  }

  vector<int> WhiteheadSimpleGraph::cutVertices()const {
    if (!isUndirected()) msgs::error("Cut vertices for undirected WG only.");
    
    CutVertices cv( (const int**)theAdjMatrix, theSize );
    cv.compute();
    return  cv.getCutVertices();
    
  }
  
  int WhiteheadSimpleGraph::numberOfCutVerticesBruteForce()const {

    vector<int> vertices = cutVerticesBruteForce();
    
    int nv = 0;
    for (int i=0;i<theSize;i++)
      if (vertices[i] > 0)
	nv++;

    return nv;
    
  }

  vector<int> WhiteheadSimpleGraph::cutVerticesBruteForce()const {
    if (!isUndirected()) msgs::error("Cut vertices for undirected WG only.");
    
    CutVertices cv( (const int**)theAdjMatrix, theSize );
    cv.computeBruteForce();
    return cv.getCutVertices();
  }


/*
void WhiteheadSimpleGraph::getEdgeMaps(const Word& theEdge, Map& m1, Map& m2)const
{
  // extract the max weight maps
  VectorOf<Word> image( nOfGenerators );
  for (int i=0;i<nOfGenerators;i++)
    image[i] = Generator( i+1 );
  
  // get the first map
  int firstMapGenIndex = abs(theEdge[0]) - 1;
  image[firstMapGenIndex] = ord( theEdge[0] ) > 0 ?  theEdge : theEdge.inverse();
  m1 = Map(theGroup,theGroup,image);
  
  // get the second map
  image[firstMapGenIndex] = Generator(firstMapGenIndex+1);
  int secondMapGenIndex = abs(theEdge[1]) - 1;
  Word edgeReverse = theEdge[1]*theEdge[0];
  image[secondMapGenIndex] = ord( theEdge[1] ) > 0 ?  edgeReverse : edgeReverse.inverse();
  m2 = Map(theGroup,theGroup,image);
}

*/

void WhiteheadSimpleGraph::printOn( ostream& out ) const
{
  for (int i=0;i<theSize;i++){
    for (int j=0;j<theSize;j++)
      out << theAdjMatrix[i][j] << " ";
    out << endl;
  }
}


