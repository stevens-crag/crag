
#include "AdvDehnAlgorithm.h"


//---------------------------------------------------------------------------//
//--------------------------- AdvDehnAlgorithm ------------------------------//
//---------------------------------------------------------------------------//


AdvDehnAlgorithm::AdvDehnAlgorithm( const FPGroup& G , const Word& w ) :
  theGroup( G ),
  theWord( w )
{
  addRay( theFSA , theFSA.newVertex( ) , w.begin( ) , w.end( ) );
}


AdvDehnAlgorithm::AdvDehnAlgorithm( const FPGroup& G , const set< Word >& gens , const Word& w ) :
  theGroup( G ),
  theWord( w )
{
  theFSA.newVertex( );
  addRay( theFSA , 0 , w.begin( ) , w.end( ) );
  for( set< Word >::const_iterator w_it=gens.begin( ) ; w_it!=gens.end( ) ; ++w_it ) 
    addLoop( theFSA , 0 , (*w_it).begin( ) , (*w_it).end( ) );
  
  // cout << theFSA << "<br><br>" << endl;
  
  set< int > candidates;
  candidates.insert( 0 );
  fold( theFSA , candidates );

  // cout << theFSA << "<br><br>" << endl;
}


//---------------------------------------------------------------------------//
//-------------------------------- builtup ----------------------------------//
//---------------------------------------------------------------------------//


bool AdvDehnAlgorithm::builtup( set< Word >* conj , int coset_limit )
{
  typedef IntLabeledGraph Graph;
  typedef Graph::edge_type edge_type;
  typedef Graph::vertex_type vertex_type;

  const map< int , vertex_type >& theVertices = theFSA.getVertices( );
  map< int , vertex_type >::const_iterator v_it = theVertices.begin( );
  
  // find a set of unchecked states (could be implemented better)
  set< int > unchecked_v;
  for( ; v_it!=theVertices.end( ) ; ++v_it )
    if( checkedStates.find( (*v_it).first )==checkedStates.end( ) )
      unchecked_v.insert( (*v_it).first );

  map< int , edge_type > tree = getGeodesicTree_out( theFSA , 0 );
  
  // add all cycles to all unchecked states
  set< int >::iterator unch_v_it = unchecked_v.begin( );
  for( ; unch_v_it!=unchecked_v.end( ) ; ++unch_v_it ) {
    
    // check if the coset limit is reached
    if( theVertices.size( )>coset_limit )
      return false;
    
    // here we check if the current state was not removed (folded with the vertex of lesser number)
    if( theVertices.find( *unch_v_it )!=theVertices.end( ) ) {
      
      // if we need to construct the set of conjugators:
      if( conj ) {
        // find the shortest path from 0 to *unch_pts_it, and add to conj
	
        list< int > labells;
	list< edge_type > path = readoffGeodesicTree( tree , *unch_v_it ).second;
	for( list< edge_type >::iterator p_it=path.begin() ; p_it!=path.end() ; ++p_it )
	  labells.push_front( (*p_it).theLabel );
	conj->insert( Word( labells ) );
      }
      
      // add new loops at the current vertex
      vector< Word > relators = theGroup.relators( );
      for( int r=0 ; r<relators.size( ) ; ++r ) {
        Word relator = relators[r];
        int len = relator.length( );
        for( int l=0 ; l<len ; ++l , relator.cyclicLeftShift( ) )
          addLoop( theFSA , *unch_v_it , relator.begin( ) , relator.end( ) );
      }

      // fold the obtained graph at the current vertex if necessary
      set< int > candidates;
      candidates.insert( *unch_v_it );
      fold( theFSA , candidates );
      
      checkedStates.insert( *unch_v_it );
    }
  } // for( ; unch_pts_it

  return true;
}


//---------------------------------------------------------------------------//
//--------------------------------- isLoop ----------------------------------//
//---------------------------------------------------------------------------//


bool AdvDehnAlgorithm::isLoop( const Word& w ) const
{
  return trace( theFSA , 0 , w.begin( ) , w.end( ) )==0;
}
