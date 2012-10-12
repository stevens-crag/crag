// Copyright (C) 2005 Alexander Ushakov
// Contents: Implementation of class FSA
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//

#include "iterator"
#include "algorithm"
#include "FSA.h"


//---------------------------------------------------------------------------//
//--------------------------------- FSA -------------------------------------//
//---------------------------------------------------------------------------//


ostream& operator << ( ostream& os , const FSA& g )
{
	typedef FSA::state_type state_type;
	typedef FSA::edge_type edge_type;
	const map< int, state_type >& theStates = g.getStates( );

	os << "{" << endl;
	map< int, state_type >::const_iterator v_it = theStates.begin( );
	for( ; v_it!=theStates.end( ) ; ++v_it ) {
		const set< edge_type >& out = (*v_it).second.out;
		os << "  " << (*v_it).second.theState << " -> ";
		for( set< edge_type >::const_iterator e_it=out.begin( ) ; e_it!=out.end( ) ; ++e_it ) {
			if( e_it!=out.begin( ) )
				os << ", ";
			os << "(" << (*e_it).label << "," << (*e_it).target << ")";
		}
		os << ";" << endl;
	}

	os << "  IS = ";
	const set< int >& IS = g.getInitStates( );
	copy( IS.begin() , IS.end() , ostream_iterator<int>( os , ", " ) );
	os << endl;

	os << "  TS = ";
	const set< int >& TS = g.getTermStates( );
	copy( TS.begin() , TS.end() , ostream_iterator<int>( os , ", " ) );
	os << endl;

	os << "}" << endl;

	return os;
}


//---------------------------------------------------------------------------//
//------------------------------ operator * ---------------------------------//
//---------------------------------------------------------------------------//

FSA FSA::operator * ( const FSA& F ) const
{
	FSA result;

	typedef pair< set< int > , set< int > > NODE;
	map< NODE , int > S;
	list< NODE > toCheck;
	const map< int , FSAState >& S1 =   getStates( );
	const map< int , FSAState >& S2 = F.getStates( );
	const set< int >& T1 =   getTermStates( );
	const set< int >& T2 = F.getTermStates( );

	// the initial vertex
	NODE init( getInitStates( ) , F.getInitStates( ) );
	if( init.first.empty( ) || init.second.empty( ) )
		return result;
	result.makeInitial( S[init] = result.newState( ) );
	toCheck.push_back( init );
	
	while( !toCheck.empty( ) ) {

		NODE curNode = *toCheck.begin( );
		toCheck.erase( toCheck.begin( ) );
		int curNodeNum = S[curNode];
		map< int , NODE > target;

		set< int >::iterator st_it = curNode.first.begin( );
		for( ; st_it!=curNode.first.end( ) ; ++st_it ) {
			const FSAState& st = (*S1.find( *st_it )).second;
			set< FSAEdge >::const_iterator out_it = st.out.begin( );
			for( ; out_it!=st.out.end( ) ; ++out_it ) {
				target[(*out_it).label].first.insert( (*out_it).target );
			}
		}

		st_it = curNode.second.begin( );
		for( ; st_it!=curNode.second.end( ) ; ++st_it ) {
			const FSAState& st = (*S2.find( *st_it )).second;
			set< FSAEdge >::const_iterator out_it = st.out.begin( );
			for( ; out_it!=st.out.end( ) ; ++out_it ) {
				target[(*out_it).label].second.insert( (*out_it).target );	
			}
		}
	
		map< int , NODE >::iterator target_it = target.begin( );
		for( ; target_it!=target.end( ) ; ++target_it ) {

			int label = (*target_it).first;
			NODE& node = (*target_it).second;
			if( (*target_it).second.first.empty( ) || (*target_it).second.second.empty( ) )
				continue;

			map< NODE , int >::iterator node_it = S.find( node );
			if( node_it==S.end( ) ) {
				S[node] = result.newState( );
				node_it = S.find( node );
				toCheck.push_back( node );

				list< int > intersection1;
				set_intersection( T1.begin( ) , T1.end( ) , node.first.begin( ) , node.first.end( ) , intersection1.begin( ) );
				if( !intersection1.empty( ) ) {
					list< int > intersection2;
					set_intersection( T2.begin( ) , T2.end( ) , node.second.begin( ) , node.second.end( ) , intersection2.begin( ) );
					if( !intersection2.empty( ) )
						result.makeTerminal( (*node_it).second );
				}
			}
			result.newEdge( curNodeNum , (*node_it).second , label );
		}
	}
	

	return result;
}


//---------------------------------------------------------------------------//
//---------------------------- isDeterministic ------------------------------//
//---------------------------------------------------------------------------//


bool FSA::isDeterministic( ) const
{
	const set< int >& IS = getInitStates( );
	if( IS.size( )!=1 )
		return false;

	const map< int , FSAState >& theStates = getStates( );
	map< int , FSAState >::const_iterator s_it = theStates.begin( );
	for( ; s_it!=theStates.end() ; ++s_it ) {
		const FSAState& s = (*s_it).second;
		const set< FSAEdge >& out = s.out;
		set< FSAEdge >::const_iterator out_it1 = out.begin( ), out_it2 = out.begin( );
		for( ++out_it2 ; out_it2!=out.end( ) ; ++out_it1, ++out_it2 )
			if( (*out_it1).label==(*out_it2).label )
				return false;
	}

	return true;
}


//---------------------------------------------------------------------------//
//------------------------------ deterministic ------------------------------//
//---------------------------------------------------------------------------//


FSA FSA::deterministic( ) const
{
	typedef set< int > NODE;
	map< NODE , int > S;
	list< NODE > toCheck;
	const map< int , FSAState >& S1 =   getStates( );
	const set< int >& I1 =   getInitStates( );
	const set< int >& T1 =   getTermStates( );

	FSA result;
	if( I1.empty( ) )
		return result;

	// the initial vertex
	NODE init( I1 );
	result.makeInitial( S[init] = result.newState( ) );
	toCheck.push_back( init );
	
	while( !toCheck.empty( ) ) {

		NODE curNode = *toCheck.begin( );
		toCheck.erase( toCheck.begin( ) );
		int curNodeNum = S[curNode];
		map< int , NODE > target;

		set< int >::iterator st_it = curNode.begin( );
		for( ; st_it!=curNode.end( ) ; ++st_it ) {
			const FSAState& st = (*S1.find( *st_it )).second;
			set< FSAEdge >::const_iterator out_it = st.out.begin( );
			for( ; out_it!=st.out.end( ) ; ++out_it ) {
				target[(*out_it).label].insert( (*out_it).target );
			}
		}

		map< int , NODE >::iterator target_it = target.begin( );
		for( ; target_it!=target.end( ) ; ++target_it ) {

			int label = (*target_it).first;
			NODE& node = (*target_it).second;
			if( (*target_it).second.empty( ) )
				continue;

			map< NODE , int >::iterator node_it = S.find( node );
			if( node_it==S.end( ) ) {
				S[node] = result.newState( );
				node_it = S.find( node );
				toCheck.push_back( node );

				// cout << "  " << S[node] << "  ->  ";
				// copy( node.begin() , node.end() , ostream_iterator< int > ( cout , " " ) );
				// cout << endl;

				list< int > intersection1;
				set_intersection( T1.begin( ) , T1.end( ) , node.begin( ) , node.end( ) , intersection1.begin( ) );
				if( !intersection1.empty( ) )
					result.makeTerminal( (*node_it).second );
			}
			result.newEdge( curNodeNum , (*node_it).second ,  label );
		}

	}

	return result;
}


//---------------------------------------------------------------------------//
//-------------------------------- operator== -------------------------------//
//---------------------------------------------------------------------------//


bool FSA::operator== ( const FSA& F ) const
{
	return getStates( )==F.getStates( ) && 
		getInitStates( )==F.getInitStates( ) &&
		getTermStates( )==F.getTermStates( );
}
