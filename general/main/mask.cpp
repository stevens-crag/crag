
#include "Permutation.h"
#include "RanlibCPP.h"


#include "tuples.h"
#include "set"
#include "list"
#include "map"
#include "vector"
#include "iostream"
using namespace std;


typedef pair< int , int > PII;
typedef vector< int > CYCLE;
typedef triple< PII , int , triple< int , int , int > > EDU;
// error-detection-unit
typedef list< EDU > EDC; 
//error-detection-code


ostream& print( ostream& os , list< vector< int > >& partition )
{
	list< vector< int > >::iterator pn_it = partition.begin( );
	for( int i=0 ; pn_it!=partition.end( ) ; ++pn_it, ++i ) {
		vector< int >& part = *pn_it;
		vector< int >::iterator p_it = part.begin( );
		cout << i << ") ";
		for( ; p_it!=part.end() ; ++p_it )
			cout << *p_it << ", ";
		cout << endl;
	}
	return os;
}


//---------------------------------------------------------------------------//
//--------------------------- generatePartition -----------------------------//
//---------------------------------------------------------------------------//


list< vector< int > > generatePartition( int N , int min , int max )
{
	list< vector< int > > partition;

	Permutation p = Permutation::random( N );
	// cout << p << endl;

	vector< int > part( min , 0 );
	for( int i=0 ; i<min ; ++i )
		part[i] = p[i];
	partition.push_back( part );
	for( int i=min ; i<N ; ) {
		int _max = max<N-i ? max : N-i;
		int _min = min<N-i ? min : N-i;
		int sz = RandLib::ur.irand( _min , _max );
		if( sz==1 ) {
			(*partition.begin()).push_back( p[i] );
		} else {
			vector< int > part( sz , 0 );
			for( int j=0 ; j<sz ; ++j )
				part[j] = p[i+j];
			partition.push_back( part );
		}
		i += sz;
	}

	// vector< vector< int > > result;
	// result.insert( result.end() , partition.begin() , partition.end( ) );
	return partition;
}


//---------------------------------------------------------------------------//
//---------------------------- generateSubgroup -----------------------------//
//---------------------------------------------------------------------------//

/*
list< Permutation > generateSubgroup( int N , const list< vector< int > >& partition )
{
	list< Permutation > result;

	for( list< vector< int > >::const_iterator pn_it = partition.begin( ) ; pn_it!=partition.end( ) ; ++pn_it )
		result.push_back( Permutation::CYCLE( N , *pn_it ) );

	return result;
}
*/

//---------------------------------------------------------------------------//
//------------------------------ generateKey --------------------------------//
//---------------------------------------------------------------------------//


Permutation generateKey( int N , const list< vector< int > >& partition )
{
	Permutation result( N );

	list< vector< int > >::const_iterator pn_it = partition.begin( );
	for( ; pn_it!=partition.end( ) ; ++pn_it ) {
		int power = RandLib::ur.irand( 0 , (*pn_it).size( )-1 );
		for( int i=0 ; i<power ; ++i )
			result.left_mult_by_cycle( *pn_it );
	}

	return result;
}


//---------------------------------------------------------------------------//
//------------------------------- perturb -----------------------------------//
//---------------------------------------------------------------------------//


void perturb( int N , list< vector< int > >& partition , int crit )
{
	set< int > to_perturb;
	while( to_perturb.size( )<crit && partition.size( )>0 ) {
		int sz = partition.size( );
		int n = RandLib::ur.irand( 0 , sz-1 );
		list< vector< int > >::iterator pn_it = partition.begin( );
		for( int i=0 ; i<n ; ++i, ++pn_it );
		// cout << "    n = " << n << endl;
		to_perturb.insert( (*pn_it).begin( ) , (*pn_it).end( ) );
		partition.erase( pn_it );
	}

	list< vector< int > > partition2 = generatePartition( to_perturb.size( ) , 2 , 7 );
	vector< int > mapping;
	mapping.insert( mapping.end( ) , to_perturb.begin() , to_perturb.end( ) );
	for( list< vector< int > >::iterator pn_it2 = partition2.begin( ) ; pn_it2!=partition2.end( ) ; ++pn_it2 ) {
		vector< int >& part = *pn_it2;
		for( int i=0 ; i<part.size() ; ++i ) part[i] = mapping[part[i]];
	}
	partition.insert( partition.end( ) , partition2.begin( ) , partition2.end( ) );
}


//---------------------------------------------------------------------------//
//------------------------------ generateKey --------------------------------//
//---------------------------------------------------------------------------//


quadruple< Permutation , set< Permutation > , Permutation , set< Permutation > > 
generateKey( int N )
{
	//quadruple< Permutation , set< Permutation > , Permutation , set< Permutation > > result;

	// generate initial setup
	list< vector< int > > partition = generatePartition( N , 2 , 7 );
	Permutation A = generateKey( N , partition );
	// print( cout , partition ) << endl;
	// cout << "A = " << A << endl << endl;

	// perturb subgroup and partition, and construct perturbed permutation
	list< vector< int > > partition2 = partition;
	perturb( N , partition2 , N/10 );
	Permutation B = generateKey( N , partition2 );
	set< Permutation > _sbgp, _sbgp2;

	// print( cout , partition2 ) << endl;
	// cout << "B = " << B << endl << endl;
	
/*
	_sbgp .insert( sbgp .begin() , sbgp .end( ) );
	_sbgp2.insert( sbgp2.begin() , sbgp2.end( ) );
*/
	return quadruple< Permutation , set< Permutation > , Permutation , set< Permutation > >( A , _sbgp , B , _sbgp2 );
}


//---------------------------------------------------------------------------//
//------------------------------ generateBase -------------------------------//
//---------------------------------------------------------------------------//


Permutation generateBase( int N )
{
	return Permutation::random( N );
}


//---------------------------------------------------------------------------//
//------------------------------- exchange ----------------------------------//
//---------------------------------------------------------------------------//


Permutation perturb( int N , const Permutation& p , int crit )
{
	Permutation result = p;

	Permutation r = Permutation::random( N );
	for( int i=0 ; i<crit-1 ; ++i )
		result.change( r[i] , r[i+1] );

	return result;
}


//---------------------------------------------------------------------------//
//------------------------------- exchange ----------------------------------//
//---------------------------------------------------------------------------//


vector< vector< int > > generate_mixing_cycles( int N , int sz , int it )
{
	list< vector< int > > partition1 = generatePartition( N , sz , sz );
	list< vector< int > > partition2 = generatePartition( N , sz , sz );

	vector< vector< int > > result;
	result.insert( result.begin( ) , partition1.begin( ) , partition1.end( ) );
	result.insert( result.begin( ) , partition2.begin( ) , partition2.end( ) );
	return result;
}


//---------------------------------------------------------------------------//
//------------------------------- exchange ----------------------------------//
//---------------------------------------------------------------------------//


list< triple< pair< int , int > , int , triple< int , int , int > > > 
generate_error_detection_code(int N , int loops,
															const vector< vector< int > >& mixing_cycles ,
															const Permutation& K )
{
	EDC result;

	int sz = mixing_cycles.size( );
	Permutation cur_perm = K;
	Permutation cur_inv  = -K;
	Permutation left_mult( N );
	Permutation inv_mult ( N );
	for( int i=0 ; i<loops ; ++i ) {
		Permutation P = Permutation::random( N );
		const vector< int >& V = P.getVector( );
		for( int j=0 ; j<N ; j+=3 ) {
			int c1 = RandLib::ur.irand( 0 , sz-1 );
			int c2 = RandLib::ur.irand( 0 , sz-1 );
			Permutation::lr_multiply_by_cycles( cur_perm , cur_inv  , mixing_cycles[c1] , mixing_cycles[c2] );
			Permutation::lr_multiply_by_cycles( left_mult, inv_mult , mixing_cycles[c1] , vector< int >( ) );

			int a = inv_mult[V[j]];
			int b = inv_mult[V[j+1]];
			int c = inv_mult[V[j+2]];
			int sum = cur_perm[a]+cur_perm[b]+cur_perm[c];
			result.push_back( EDU( pair< int , int >( c1 , c2 ) , sum , triple< int , int , int >( a , b , c ) ) );
		}
	}

	return result;
}


//---------------------------------------------------------------------------//
//----------------------------- detect_errors -------------------------------//
//---------------------------------------------------------------------------//


set< int > detect_errors( int N , const vector< vector< int > >& mixing_cycles , const EDC& D , const Permutation& K , int max )
{
	Permutation cur_perm = K;
	Permutation cur_inv  = -K;
	Permutation left_mult( N );

	vector< int > weights( N , 0 );

	EDC::const_iterator d_it = D.begin( );
	for( ; d_it!=D.end( ) ; ++d_it ) {
		const EDU& edu = *d_it;
		int c1 = edu.first.first;
		int c2 = edu.first.second;
		int sum = edu.second;
		int a = edu.third.first;
		int b = edu.third.second;
		int c = edu.third.third;
		left_mult.left_mult_by_cycle( mixing_cycles[c1] );
		Permutation::lr_multiply_by_cycles( cur_perm , cur_inv , mixing_cycles[c1] , mixing_cycles[c2] );
		int new_sum = cur_perm[a]+cur_perm[b]+cur_perm[c];
		if( sum==new_sum ) {
			// cout << "x";
		} else {
			// cout << ".";
			weights[left_mult[a]]++;
			weights[left_mult[b]]++;
			weights[left_mult[c]]++;
		}
	}
	// cout << endl;

	set< int > result;
	for( int i=0 ; i<N ; ++i ) 
		if( weights[i]<max-1 )
			result.insert( i );

/*
	set< PII >::iterator w_it = ordered_weight.begin( );
	for( int i=0 ; i<N ; ++i,++w_it )
		cout << "(" << (*w_it).first << "," << (*w_it).second << ") ";
	cout << endl;
*/

	return result;
}


//---------------------------------------------------------------------------//
//------------------------------- exchange ----------------------------------//
//---------------------------------------------------------------------------//


void exchange( )
{
	int N = 12000;

	// Generate the left pair
	quadruple< Permutation , set< Permutation > , Permutation , set< Permutation > > left_elts = generateKey( N );
	quadruple< Permutation , set< Permutation > , Permutation , set< Permutation > > right_elts = generateKey( N );
	Permutation base = generateBase( N );

	Permutation A_l = left_elts.first;
	Permutation B_l = left_elts.third;
	Permutation B_r = right_elts.first;
	Permutation A_r = right_elts.third;

	Permutation base1 = perturb( N , base , N/10 );
	Permutation base2 = perturb( N , base , N/10 );
	cout << "diff1 = " << base.difference( base1 ) << endl;
	cout << "diff2 = " << base.difference( base2 ) << endl;

	Permutation K_A = A_l * B_l * base2 * B_r * A_r;
	Permutation K_B = B_l * A_l * base1 * A_r * B_r;
	cout << "diff = " << K_A.difference( K_B ) << endl;
	// for( int i=0 ; i<N ; ++i ) if( K_A[i]!=K_B[i] ) cout << i << " ";
	// cout << endl;

	int cycles = 15;
	vector< vector< int > > mixing_cycles = generate_mixing_cycles( N , 10 , 2 );
	EDC D = generate_error_detection_code( N , cycles , mixing_cycles , K_A );

	set< int > correct_pos = detect_errors( N , mixing_cycles , D , K_B , cycles );

	// check
	set< int >::iterator cp_it = correct_pos.begin( );
	for( ; cp_it!=correct_pos.end( ) ; ++cp_it ) {
		if( K_A[*cp_it]!=K_B[*cp_it] )
			cout << "Mistake" << endl;
	}
	cout << "   Removed = " << N-correct_pos.size( ) << endl;
}


//---------------------------------------------------------------------------//
//--------------------------------- main ------------------------------------//
//---------------------------------------------------------------------------//


int main( )
{
	exchange( );

	return 0;
}
