
#include <time.h>
#include <stdlib.h>
#include <string.h>

#include "braid_group.h"
#include "ThRightNormalForm.h"
#include "ShortBraidForm.h"

#include <set>
#include <map>
#include <iterator>
#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;


//#define MPI


#ifdef MPI
#include <mpi.h>
#endif


#include "AAGKeyGeneration.h"

// #define WITH_CONJUGATED_SBGP
// when this constant is defined the program works in parallel 
// with the initial set of generators and the conjugated set of
// generators
// This is not important when checking whether S ~ B_n it just
// consumes processor time (a lot of it since conjugated elements are longer)


//---------------------------------------------------------------------------//
//------------------------------- areSeparate -------------------------------//
//---------------------------------------------------------------------------//


bool areSeparate( int N , const Word& w1 , const Word& w2 )
{
  vector< bool > letters( N , false );
  for( auto w_it=w1.begin( ) ; w_it!=w1.end( ) ; ++w_it ) {
    int g = *w_it;
    int ag = abs(g)-1;
    if( ag ) letters[ag-1] = true;
    letters[ag] = true;
    if( ag<N-1 ) letters[ag+1] = true;
  }

  for( auto w_it=w2.begin( ) ; w_it!=w2.end( ) ; ++w_it ) {
    int g = *w_it;
    int ag = abs(g)-1;
    if( letters[ag] )
      return false;
  }

  return true;
}


//---------------------------------------------------------------------------//
//------------------------------- reducePair --------------------------------//
//---------------------------------------------------------------------------//


bool reducePair( int N , set< Word >& new_elts , map< Word , bool >& checked_elts , map< Word , Word >& conj_elts , Word w1 , Word w2 , Word cw2 , ostream& os=cout )
{
  if( areSeparate( N , w1 , w2 ) ) 
    return false;

  Word nw, cnw;
  Word cw1 = conj_elts[w1];

  for( int d1=0 ; d1<2 ; ++d1 ) {

    for( int d2=0 ; d2<2 ; ++d2 ) {

      nw  = shortBraidForm( N , ( d1 ?  w1 :  w1.inverse( ) ) * ( d2 ?  w2 :  w2.inverse( ) ) );

#ifdef WITH_CONJUGATED_SBGP
      cnw = shortBraidForm( N , ( d1 ? cw1 : cw1.inverse( ) ) * ( d2 ? cw2 : cw2.inverse( ) ) );
#endif

      // if w2 or w2^-1 already exists
      if( nw.length( )==0 )
        return true;

      // if w2 is reducible
      if( nw.length( )+w1.length()==w2.length( ) ) {
        new_elts.insert( nw );
        conj_elts[nw] = cnw;
        return true;
      }

      // if w1 is reducible
      if( nw.length( )+w2.length()==w1.length( ) ) {

        checked_elts.erase( w1 );
        new_elts.insert( nw );
        conj_elts[nw] = cnw;

        os << "      " << w1 << " ~> " << nw << endl;
        return false;
      }
    }
  }

  return false;
}


//---------------------------------------------------------------------------//
//------------------------ reduceSbgpPresentation ---------------------------//
//---------------------------------------------------------------------------//


void reduceSbgpPresentation( int N , set< Word >& _new_elts , set< Word >& _checked_elts , map< Word , Word >& conj_elts , const Word& w , ostream& os=cout )
{
  set< Word > new_elts;
  new_elts.insert( w );

  map< Word , bool > checked_elts;
  set< Word >::iterator el_it = _new_elts.begin( );
  for( ; el_it!=_new_elts.end( ) ; ++el_it )
    checked_elts[*el_it] = 0;
  el_it = _checked_elts.begin( );
  for( ; el_it!=_checked_elts.end( ) ; ++el_it )
    checked_elts[*el_it] = 1;

  while( new_elts.size( ) ) {

    Word w2 = *new_elts.begin( );
    new_elts.erase( w2 );
    Word cw2 = conj_elts[w2];

    bool add = true;
    map< Word , bool > checked_elts_copy = checked_elts;
    map< Word , bool >::iterator el_it = checked_elts_copy.begin( );
    for( ; el_it!=checked_elts_copy.end( ) ; ++el_it ) {

      Word w1 = (*el_it).first;
      if( reducePair( N , new_elts , checked_elts , conj_elts , w1 , w2 , cw2 , os ) ) {
        add = false;
        break;
      }
    }

    if( add )
      checked_elts[w2] = 0;
  }

  // clean conj_elts;
  map< Word , Word > conj_elts_copy = conj_elts;
  map< Word , Word >::iterator ce_it = conj_elts_copy.begin( );
  for( ; ce_it!=conj_elts_copy.end( ) ; ++ce_it ) {
    Word w = (*ce_it).first;
    if( checked_elts.find(w)==checked_elts.end( ) )
      conj_elts.erase(w);
  }

  static int print_counter = 0;
  print_counter++;
  const int periodicity = 50;
  bool to_print = (print_counter%periodicity)==0;
  
  // update _checked_elts and _new_elts
  _new_elts.clear( );
  _checked_elts.clear( );
  map< Word , bool >::iterator el_it2 = checked_elts.begin( );
  if( to_print )
    os << "------------------------------------------------" << endl;
  for( int i=0 ; el_it2!=checked_elts.end( ) ; ++el_it2, ++i ) {
    if( (*el_it2).second )
      _checked_elts.insert( (*el_it2).first );
    else
      _new_elts.insert( (*el_it2).first );
    if( to_print ) {
      os.width( 4 );
      os << i << " : (" << (*el_it2).second << ") " << (*el_it2).first << endl;
    }
  }
  if( to_print )
    os << "------------------------------------------------" << endl;
}


//---------------------------------------------------------------------------//
//-------------------------------- add_elt ----------------------------------//
//---------------------------------------------------------------------------//


int check_if_belongs( int N , const set< Word >& new_elts , const set< Word >& checked_elts , const Word& w )
{
  if( new_elts.find( w )!=new_elts.end( ) )
    return 1;
  if( checked_elts.find( w )!=checked_elts.end( ) )
    return 2;
  
  Word w1 = shortBraidForm( N , -w );
  if( new_elts.find( w1 )!=new_elts.end( ) )
    return 3;
  if( checked_elts.find( w1 )!=checked_elts.end( ) )
    return 4;
  
  return 0;
}

//---------------------------------------------------------------------------//
//-------------------------------- add_elt ----------------------------------//
//---------------------------------------------------------------------------//


void add_elt( int N , set< Word >& new_elts , set< Word >& checked_elts , map< Word , Word >& conj_elts , Word w , Word conj_w , int limit , ostream& os=cout ) 
{
  Word nw = shortBraidForm( N , w );

  if( nw.length( )==0 || nw.length( )>6 || nw.length( )>=limit )
    return;

  Word n_conj_w;
#ifdef WITH_CONJUGATED_SBGP
  n_conj_w = shortBraidForm( N , conj_w );
#endif
  
  if( check_if_belongs( N , new_elts , checked_elts , nw )==0 ) {
    os << nw << " -> " << nw.length( ) << endl;
    conj_elts[nw] = n_conj_w;
    reduceSbgpPresentation( N , new_elts , checked_elts , conj_elts , nw , os );
  }
}


//---------------------------------------------------------------------------//
//-------------------------- precomputeShortWords ---------------------------//
//---------------------------------------------------------------------------//


vector< Word > precomputeShortWords( int N , set< Word >& new_elts , set< Word >& checked_elts , map< Word , Word >& conj_elts , ostream& os=cout )
{
  vector< Word > result;

  vector< Word > sbgp;
  set< Word >::iterator el_it = checked_elts.begin( );
  for( ; el_it!=checked_elts.end( ) ; ++el_it )
    sbgp.push_back( *el_it );
  el_it = new_elts.begin( );
  for( ; el_it!=new_elts.end( ) ; ++el_it )
    sbgp.push_back( *el_it );


  vector< set< int > > indices( N );
  for( int i=0 ; i<sbgp.size( ) ; ++i ) {
    const Word& w = sbgp[i];
    for( auto w_it=w.begin( ) ; w_it!=w.end( ) ; ++w_it ) {
      int g = *w_it;
      int ag = abs( g )-1;
      indices[ag].insert( i );
    }
  }

  for( int i=1 ; i<N-1 ; ++i ) {

    Word desired1 = Word(  i ) * Word(-i-1 );
    Word desired2 = Word( -i ) * Word( i+1 );

    set< int >::iterator s_it1 = indices[i-1].begin( );
    for( ; s_it1!=indices[i-1].end( ) ; ++s_it1 ) {

      if( check_if_belongs( N , new_elts , checked_elts , desired1 )!=0 || check_if_belongs( N , new_elts , checked_elts , desired2 )!=0 )
	break;

      Word w1 = sbgp[*s_it1];
      Word cw1 = conj_elts[w1];
      set< int >::iterator s_it2 = indices[i].begin( );
      for( ; s_it2!=indices[i].end( ) ; ++s_it2 ) {

        Word w2 = sbgp[*s_it2];
	if( w1.length( )<=2 && w2.length( )<=2 )
	  continue;
	
        Word cw2 = conj_elts[w2];
        Word nw = shortBraidForm( N , -w1 * -w1 * -w2 * w1 * w2 * w1 );
        if( nw.length()<5 )
          add_elt( N , new_elts , checked_elts , conj_elts , nw , -cw1 * -cw1 * -cw2 * cw1 * cw2 * cw1 , 5 , os );

        nw = shortBraidForm( N , w1 * w1 * w2 * -w1 * -w2 * -w1 );
        if( nw.length()<5 )
          add_elt( N , new_elts , checked_elts , conj_elts , nw , cw1 * cw1 * cw2 * -cw1 * -cw2 * -cw1 , 5 , os );
      }
    }
  }

  return result;
}


//---------------------------------------------------------------------------//
//--------------------------- completeShortWords ----------------------------//
//---------------------------------------------------------------------------//


void locally_perturbate_elts( int N , set< Word >& new_elts , set< Word >& checked_elts , map< Word , Word >& conj_elts , ostream& os=cout )
{
  vector< Word > result;
  
  vector< Word > sbgp;
  set< Word >::iterator el_it = checked_elts.begin( );
  for( ; el_it!=checked_elts.end( ) ; ++el_it )
    sbgp.push_back( *el_it );
  el_it = new_elts.begin( );
  for( ; el_it!=new_elts.end( ) ; ++el_it )
    sbgp.push_back( *el_it );


  vector< set< int > > indices( N );
  for( int i=0 ; i<sbgp.size( ) ; ++i ) {
    const Word& w = sbgp[i];
    for( auto w_it=w.begin( ) ; w_it!=w.end( ) ; ++w_it ) {
      int g = *w_it;
      int ag = abs( g )-1;
      indices[ag].insert( i );
    }
  }

  for( int i=0 ; i<N ; ++i ) {

    set< int >::iterator s_it1 = indices[i].begin( );
    for( ; s_it1!=indices[i].end( ) ; ++s_it1 ) {

      Word w1 = sbgp[*s_it1];
      if( w1.length( )<=2 ) continue;
      
      Word cw1 = conj_elts[w1];
      set< int >::iterator s_it2 = indices[i].begin( );
      for( ; s_it2!=indices[i].end( ) ; ++s_it2 ) {
	
        Word w2 = sbgp[*s_it2];
	if( w2.length( )!=2 ) continue;
	
        Word cw2 = conj_elts[w2];
        Word nw = shortBraidForm( N , w1*w2 );
        if( nw.length( )==w1.length( ) && nw!=w1 ) {
          add_elt( N , new_elts , checked_elts , conj_elts , nw , cw1*cw2 , w1.length( )+1 , os );
	  break;
	}

        nw = shortBraidForm( N , w1*-w2 );
        if( nw.length( )==w1.length( ) && nw!=w1 ) {
          add_elt( N , new_elts , checked_elts , conj_elts , nw , cw1 * -cw2 , w1.length( )+1 , os );
	  break;
	}
	
        nw = shortBraidForm( N , -w1*w2 );
        if( nw.length( )==w1.length( ) && nw!=w1 ) {
          add_elt( N , new_elts , checked_elts , conj_elts , nw , -cw1*cw2 , w1.length( )+1 , os );
	  break;
	}

        nw = shortBraidForm( N , -w1*-w2 );
        if( nw.length( )==w1.length( ) && nw!=w1 ) {
          add_elt( N , new_elts , checked_elts , conj_elts , nw , -cw1*-cw2 , w1.length( )+1 , os );
	  break;
	}
	
      }
    }
  }
}


//---------------------------------------------------------------------------//
//--------------------------- completeShortWords ----------------------------//
//---------------------------------------------------------------------------//


void completeShortWords( int N , set< Word >& new_elts , set< Word >& checked_elts , map< Word , Word >& conj_elts , ostream& os=cout )
{
  typedef pair< int , Word > PIW;
  set< int > to_check;
  map< int , set< PIW > > words;

  set< Word >::iterator w_it = new_elts.begin( );
  for( ; w_it!=new_elts.end( ) ; ++w_it ) {
    Word w = *w_it;
    if( w.length( )!=2 )
      continue;

    // int g1 = ord(w[0]);
    // int g2 = ord(w[1]);
    int g1 = *w.begin( );
    int g2 = *(++w.begin( ));
    Word cw = conj_elts[w];
    words[g1].insert( PIW(g2,cw) );
    to_check.insert(g1);
    words[-g2].insert( PIW(-g1,cw.inverse()) );
    to_check.insert(-g2);
  }
  w_it = checked_elts.begin( );
  for( ; w_it!=checked_elts.end( ) ; ++w_it ) {
    Word w = *w_it;
    if( w.length( )!=2 )
      continue;

    // int g1 = ord(w[0]);
    // int g2 = ord(w[1]);
    int g1 = *w.begin( );
    int g2 = *(++w.begin( ));
    Word cw = conj_elts[w];
    words[g1].insert( PIW(g2,cw) );
    to_check.insert(g1);
    words[-g2].insert( PIW(-g1,cw.inverse()) );
    to_check.insert(-g2);
  }


  while( to_check.size( ) ) {

    int c = *to_check.begin();
    to_check.erase( to_check.begin() );
    set< PIW > s = words[c];

    for( set< PIW >::iterator s_it1 = s.begin( ) ; s_it1!=s.end( ) ; ++s_it1 ) {
      set< PIW >::iterator s_it2 = s_it1;
      for( ++s_it2 ; s_it2!=s.end( ) ; ++s_it2 ) {
        words[-(*s_it1).first].insert( PIW( (*s_it2).first , (*s_it1).second.inverse() * (*s_it2).second ) );
        words[-(*s_it2).first].insert( PIW( (*s_it1).first , (*s_it2).second.inverse() * (*s_it1).second ) );
      }
    }
  }

  for( map< int , set< PIW > >::iterator words_it = words.begin( ) ; words_it!=words.end( ) ; ++words_it ) {

    int g1 = (*words_it).first;
    set< PIW >& ws = (*words_it).second;

    Word w1( g1 );
    if( check_if_belongs( N , new_elts , checked_elts , w1 )>0 )
      continue;
    
    for( set< PIW >::iterator ws_it=ws.begin( ) ; ws_it!=ws.end( ) ; ++ws_it ) {
      
      Word w2( (*ws_it).first );
      if( check_if_belongs( N , new_elts , checked_elts , w2 )>0 )
	continue;
      
      Word w( g1 );
      w *= Word( (*ws_it).first );
      add_elt( N , new_elts , checked_elts , conj_elts , w , (*ws_it).second , 3 , os );
    }
  }

}


//---------------------------------------------------------------------------//
//----------------------------- randomEltOfSet ------------------------------//
//---------------------------------------------------------------------------//


Word randomEltOfSet( const set< Word >& S )
{
  int n = rand( ) % S.size( );
  if( rand( )%3==0 ) 
    n = S.size( )-1;
  set< Word >::const_iterator w_it = S.begin( );
  for( int i=0 ; i<n ; ++i ) ++w_it;
  return *w_it;
}


//---------------------------------------------------------------------------//
//--------------------------- simplify_braid_sbgp ---------------------------//
//---------------------------------------------------------------------------//


vector< Word > simplify_braid_sbgp( int N , vector< Word > sbgp , const vector< Word >& conjSbgp , ostream& os = cout , bool thorough_check=true )
{
  // union of these to guys is a generating set
  set< Word > new_elts;
  set< Word > checked_elts;

  // here we store conjugated to the generating set
  map< Word , Word > conj_elts;
  // initial setup
  for( int i=0 ; i<sbgp.size( ) ; ++i ) {
    new_elts.insert( sbgp[i] );
    conj_elts[sbgp[i]] = conjSbgp[i];
  }

  if( thorough_check ) {
    
    precomputeShortWords( N , new_elts , checked_elts , conj_elts );
    while( new_elts.size( ) ) {
      
      // Word w1 = *new_elts.begin( );
      Word w1 = randomEltOfSet( new_elts );
      Word cw1 = conj_elts[w1];
      int l1 = w1.length( );
      new_elts.erase( w1 );
      checked_elts.insert( w1 );
      
      set< Word > checked_elts_copy = checked_elts;
      set< Word >::iterator ch_it = checked_elts_copy.begin( );
      for( ; ch_it!=checked_elts_copy.end( ) ; ++ch_it ) {
	
	Word w2 = *ch_it;
	Word cw2 = conj_elts[w2];
	if( !areSeparate( N , w1 , w2 ) ) {
	  int l2 = w2.length( );
	  for( int d=0 ; d<2 ; ++d ) {
	    w2 = w2.inverse( );
	    add_elt( N , new_elts , checked_elts , conj_elts ,  w2*w1*-w2 ,  cw2*cw1*-cw2 , l2+l1+l2 , os );
	    add_elt( N , new_elts , checked_elts , conj_elts , -w2*w1* w2 , -cw2*cw1* cw2 , l2+l1+l2 , os );
	    add_elt( N , new_elts , checked_elts , conj_elts ,  w1*w2*-w1 ,  cw1*cw2*-cw1 , l1+l2+l1 , os );
	    add_elt( N , new_elts , checked_elts , conj_elts , -w1*w2* w1 , -cw1*cw2* cw1 , l1+l2+l1 , os );
	    add_elt( N , new_elts , checked_elts , conj_elts ,  w2* w1 ,  cw2* cw1 , l1+l2 , os );
	    add_elt( N , new_elts , checked_elts , conj_elts ,  w2*-w1 ,  cw2*-cw1 , l1+l2 , os );
	    add_elt( N , new_elts , checked_elts , conj_elts ,  w1* w2 ,  cw1* cw2 , l1+l2 , os );
	    add_elt( N , new_elts , checked_elts , conj_elts , -w1* w2 , -cw1* cw2 , l1+l2 , os );
	    
	  }
	}
      }
      if( rand( )%3==0 ) {
	precomputeShortWords( N , new_elts , checked_elts , conj_elts );
	completeShortWords( N , new_elts , checked_elts , conj_elts );
      }
    } 
  } else {
    bool progress = true;
    int perturbations = 0;
    while( progress ) {
      progress = false;
      set< Word > new_elts_old = new_elts;
      set< Word > checked_elts_old = checked_elts;
      os << endl << "Precompute:" << endl << endl;
      precomputeShortWords( N , new_elts , checked_elts , conj_elts , os );
      os << endl << "Complition:" << endl << endl;
      completeShortWords( N , new_elts , checked_elts , conj_elts , os );
      if( new_elts_old!=new_elts || checked_elts_old!=checked_elts )
	progress = true;
      else {
	if( perturbations<1 ) {
	  os << endl << "Perturbations:" << endl << endl;
	  locally_perturbate_elts( N , new_elts , checked_elts , conj_elts , os );
	  ++perturbations;
	  progress = true;
	}
      }
    }
    return vector< Word >( new_elts.begin( ) , new_elts.end( ) );
  }
  
  return vector< Word >( checked_elts.begin( ) , checked_elts.end( ) );
}


//---------------------------------------------------------------------------//
//---------------------------------- main -----------------------------------//
//---------------------------------------------------------------------------//


double average_sbgp_generator_length( const vector< Word >& sbgp )
{
  if( sbgp.size()==0 ) return 0;
  double result = 0;
  for( int i=0 ; i<sbgp.size( ) ; ++i )
    result += sbgp[i].length( );
  
  return result/sbgp.size( );
}


//---------------------------------------------------------------------------//
//---------------------------------- main -----------------------------------//
//---------------------------------------------------------------------------//


void choose_parameters( const char* res_filename , int& N  , int& num_gens , int& min_len )
{
  int RANK      [] = {80};
  int NUM_GENS  [] = {24,26,28,30,32,34,36};
  int MIN_LENGTH[] = {11,13,15,17,19,21,23};
  
  typedef triple< int , int , int > TIII;
  map< TIII , int > done_exp;
  for( int r=0 ; r<sizeof(RANK)/sizeof(int) ; ++r )
    for( int i=0 ; i<sizeof(NUM_GENS)/sizeof(int) ; ++i )
      for( int j=0 ; j<sizeof(MIN_LENGTH)/sizeof(int) ; ++j )
	done_exp[ TIII( RANK[r] , NUM_GENS[i] , MIN_LENGTH[j] ) ] = 0;
  
  const int buffer_size = 1024;
  char buffer[buffer_size];

  ifstream IF( res_filename );
  while( !IF.eof( ) ) {

    IF.getline( buffer , buffer_size );
    if( strlen(buffer)<=1 )
      break;
    // cout << "* " << buffer << endl;
    int n = atoi( buffer );
    int g = atoi( strchr( buffer , ',' )+1 );
    int l = atoi( strchr( strchr( buffer , ',' )+1 , ',' )+1 );
    // cout << "    " << n << " " << g << " " << l << endl;
    done_exp[ TIII( n , g , l ) ]++;
  }
  
  typedef quadruple< int , int , int , int > QIIII;
  set< QIIII > ordered_exp;
  map< TIII , int >::iterator de_it = done_exp.begin( );
  for( ; de_it!=done_exp.end( ) ; ++de_it )
    ordered_exp.insert( QIIII( (*de_it).second , (*de_it).first.first , (*de_it).first.second , (*de_it).first.third ) );

  QIIII to_check = *ordered_exp.begin( );
  
  N = to_check.second;
  num_gens = to_check.third;
  min_len = to_check.fourth;
  
  // cout << "    " << N << " " << num_gens << " " << min_len << endl;
}


//---------------------------------------------------------------------------//
//---------------------------------- main -----------------------------------//
//---------------------------------------------------------------------------//


int main( int argc, char** argv )
{

#ifdef MPI

  /// Initialize MPI

  MPI_Status status;      /* return status for receive */

  /* Initialize MPI */
  MPI_Init(&argc, &argv);

  /* Get my process rank */
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  /* How many processes are being used? */
  int pn;
  MPI_Comm_size(MPI_COMM_WORLD, &pn);

  /// End initialize MPI

  /// Init OUTPUT
  stringstream f_name;
  f_name << "Details" << time(NULL) << ".p"<<my_rank<<'\0'<<flush;
  ofstream  OF_details(f_name.str().c_str());
  OF_details << "Rank - " << my_rank << " N - " << pn <<endl;


  //  MAIN PROCESS 0
  if (my_rank == 0) {
    for (int source = 1; source < pn; source++) {
      int res;
      int tag = 0;
      MPI_Recv(&res, 1, MPI_FLOAT, source, tag,
               MPI_COMM_WORLD, &status);
      OF_details << "Recieved from  " << source << " result : " << res << endl;
    }
  } else {

#else


  ofstream OF_details( "details.txt" );

#endif

  const char* res_filename = "results.txt";

  const int PROGRAM_LIFE_TIME = 7*24*60*60;
  int time1 = time(0);

  while( time(0)-time1<PROGRAM_LIFE_TIME ) {
    
    // Define parameters
    int N;
    int num_gens;
    int min_len;
    choose_parameters( res_filename , N  , num_gens , min_len );
    int max_len = min_len+2;
    int AliceDecompositionLength = 10;
    int BobDecompositionLength = 10;
    
    
    OF_details << "------------------------------------------------" << endl;
    OF_details << "Rank of the group = " << N << endl;
    OF_details << "Minimium generator length = " << min_len << endl;
    OF_details << "Number of subgroup generators = " << num_gens << endl;
    OF_details << "------------------------------------------------" << endl;
    
    // Generate AAG instance
    AAGProtocolInstance AAG = AAGProtocolInstance::random( N , num_gens , min_len , max_len , 
							   AliceDecompositionLength , BobDecompositionLength );
    
    // simplify the subgroup
    vector< Word > Sbgp_A  = AAG.getAlicePublicSbgp( );
    vector< Word > Sbgp_A2 = AAG.getAliceConjSbgp  ( );
    int time1 = time(0);
    vector< Word > new_sbgp_A =  simplify_braid_sbgp( N , Sbgp_A , Sbgp_A2 , OF_details , false );
    int time2 = time(0);
    
    double av_length = average_sbgp_generator_length( new_sbgp_A );
    OF_details << "Init . Average length = " << average_sbgp_generator_length( Sbgp_A ) << endl;
    OF_details << "Final. Average length = " << av_length << endl;
    
    {
      ofstream OF_result( res_filename , ios::app );
      OF_result.width(4); OF_result << N << ",";
      OF_result.width(3); OF_result << num_gens << ",";
      OF_result.width(3); OF_result << min_len << ",";
      OF_result.width(9); OF_result << average_sbgp_generator_length( Sbgp_A ) << ",";
      OF_result.width(9); OF_result << average_sbgp_generator_length( new_sbgp_A ) << ",";
      OF_result.width(6); OF_result << time2-time1 << endl;
    }
  }
  
#ifdef MPI
  /// SEND RESULTS TO PROCESS 0
    int dest = 0;
    int tag = 0;
    int success = 1;
    MPI_Send(&success, 1, MPI_FLOAT, dest, tag, MPI_COMM_WORLD);

  }

  /* Shut down MPI */
  MPI_Finalize();
#endif



  return 0;
}
