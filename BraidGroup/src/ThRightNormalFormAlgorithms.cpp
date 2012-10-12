
#include "Word.h"
#include "ShortBraidForm.h"
#include "ThRightNormalFormAlgorithms.h"


//---------------------------------------------------------------------------//
//----------------------- ThRightNormalFormAlgorithms -----------------------//
//---------------------------------------------------------------------------//


ostream& printon( ostream& os , const vector< ThRightNormalForm >& elts ) 
{
  int sz = elts.size( );
  
  cout << "( ";
  for( int i=0 ; i<sz ; ++i )
    cout << elts[i].getPower( ) << " ";
  cout << ")";
  
  return os;
}


//---------------------------------------------------------------------------//
//----------------------- ThRightNormalFormAlgorithms -----------------------//
//---------------------------------------------------------------------------//


map< vector< ThRightNormalForm > , ThRightNormalForm > processTuple( int rank , pair< vector< ThRightNormalForm > , ThRightNormalForm > pr )
{
  vector< ThRightNormalForm > tuple = pr.first;
  ThRightNormalForm tuple_conjugator = pr.second;
  
  const Permutation omega = ThRightNormalForm::getHalfTwistPermutation( rank );

  
  int sz = tuple.size( );
  map< vector< ThRightNormalForm > , ThRightNormalForm > result;
  

  // compute the set of simple conjugators for the current tuple
  set< Permutation > conjugators = getSimpleConjugators( rank , tuple );

  /*
  cout << "|conjugators| = " << conjugators.size() << endl;
  int lr = 0;
  for( set< Permutation >::iterator c_it=conjugators.begin( ) ; c_it!=conjugators.end( ) ; ++c_it )
    lr += (*c_it).geodesic().size( );
  cout << "     lr = " << lr << endl;
  
  
  for( int i=0 ; i<sz && conjugators.size( )>1 ; ++i ) {
    for( int j=0 ; j<i && conjugators.size( )>1 ; ++j ) {

      
      ThRightNormalForm nf;
      
      nf = -tuple[i] * -tuple[j] * tuple[i] * tuple[j];
      if( !nf.isTrivial( ) ) {
	vector< ThRightNormalForm > tuple2 = tuple;
	tuple2.push_back( nf );
	set< Permutation > conjugators2 = getSimpleConjugators( rank , tuple2 );
	if( conjugators!=conjugators2 ) {
	  conjugators = conjugators2;
	  tuple = tuple2;
	  ++sz;
	  cout << "extension : [" << i << "," << j << "]" << endl;

	  lr = 0;
	  cout << "|conjugators| = " << conjugators.size() << endl;
	  for( set< Permutation >::iterator c_it=conjugators.begin( ) ; c_it!=conjugators.end( ) ; ++c_it )
	    lr += (*c_it).geodesic().size( );
	  cout << "     lr = " << lr << endl;
	}
      }

      nf = -tuple[j] * -tuple[i] * tuple[j] * tuple[i];
      if( !nf.isTrivial( ) ) {
	vector< ThRightNormalForm > tuple2 = tuple;
	tuple2.push_back( nf );
	set< Permutation > conjugators2 = getSimpleConjugators( rank , tuple2 );
	if( conjugators!=conjugators2 ) {
	  conjugators = conjugators2;
	  tuple = tuple2;
	  ++sz;
	  cout << "extension : [" << j << "," << i << "]" << endl;

	  lr = 0;
	  cout << "|conjugators| = " << conjugators.size() << endl;
	  for( set< Permutation >::iterator c_it=conjugators.begin( ) ; c_it!=conjugators.end( ) ; ++c_it )
	    lr += (*c_it).geodesic().size( );
	  cout << "     lr = " << lr << endl;
	}
      }
      
      nf = tuple[i] * tuple[j];
      if( !nf.isTrivial( ) ) {
	vector< ThRightNormalForm > tuple2 = tuple;
	tuple2.push_back( nf );
	set< Permutation > conjugators2 = getSimpleConjugators( rank , tuple2 );
	if( conjugators!=conjugators2 ) {
	  conjugators = conjugators2;
	  tuple = tuple2;
	  ++sz;
	  cout << "extension : " << i << " * " << j << endl;

	  lr = 0;
	  cout << "|conjugators| = " << conjugators.size() << endl;
	  for( set< Permutation >::iterator c_it=conjugators.begin( ) ; c_it!=conjugators.end( ) ; ++c_it )
	    lr += (*c_it).geodesic().size( );
	  cout << "     lr = " << lr << endl;
	}
      }

      nf = tuple[j] * tuple[i];
      if( !nf.isTrivial( ) ) {
	vector< ThRightNormalForm > tuple2 = tuple;
	tuple2.push_back( nf );
	set< Permutation > conjugators2 = getSimpleConjugators( rank , tuple2 );
	if( conjugators!=conjugators2 ) {
	  conjugators = conjugators2;
	  tuple = tuple2;
	  ++sz;
	  cout << "extension : " << j << " * " << i << endl;

	  lr = 0;
	  cout << "|conjugators| = " << conjugators.size() << endl;
	  for( set< Permutation >::iterator c_it=conjugators.begin( ) ; c_it!=conjugators.end( ) ; ++c_it )
	    lr += (*c_it).geodesic().size( );
	  cout << "     lr = " << lr << endl;
	}
      }


    }
  }
  */

  // conjugate
  for( set< Permutation >::const_iterator c_it=conjugators.begin( ) ; c_it!=conjugators.end( ) ; ++c_it ) {
    
    vector< ThRightNormalForm > new_tuple = tuple;
    ThRightNormalForm conjugator( rank , 0 , list<Permutation>(1,*c_it) );
    if( *c_it==omega )
      conjugator = ThRightNormalForm( rank , 1 , list<Permutation>( ) );
    

    for( int i=0 ; i<sz ; ++i )
      new_tuple[i]  = -conjugator*new_tuple[i]*conjugator;
    result[new_tuple] = tuple_conjugator*conjugator;
  }
  
  return result;
}


//---------------------------------------------------------------------------//
//----------------------- ThRightNormalFormAlgorithms -----------------------//
//---------------------------------------------------------------------------//


pair< vector< ThRightNormalForm > , ThRightNormalForm > getSummitSetRepresentative( int rank , const vector< ThRightNormalForm >& elts )
{

  int sz = elts.size( );

  cout << "Given : "; printon( cout , elts ) << endl;
  
  map< vector< ThRightNormalForm > , ThRightNormalForm > C;  // checked
  map< vector< ThRightNormalForm > , ThRightNormalForm > N;  // new
  N[elts] = ThRightNormalForm( rank , 0 , list<Permutation>( ) );
  
  int counter=0;
  for( ; !N.empty( ) && counter<100000 ; ++counter ) {

    cout << "counter = " << counter << endl;
    cout << "  C.size() = " << C.size() << endl;
    cout << "  N.size() = " << N.size() << endl;

    // take a new tuple from the summit set
    pair< vector< ThRightNormalForm > , ThRightNormalForm > cur_pr = *N.begin( );
    const vector< ThRightNormalForm >& cur_tuple = cur_pr.first;
    const ThRightNormalForm& cur_conjugator = cur_pr.second;
    N.erase( N.begin( ) );
    C.insert( cur_pr );

    map< vector< ThRightNormalForm > , ThRightNormalForm > new_elts = processTuple( rank , cur_pr );
    
    // process new tuples
    for( map< vector< ThRightNormalForm > , ThRightNormalForm >::const_iterator ne_it=new_elts.begin( ) ; ne_it!=new_elts.end( ) ; ++ne_it ) {
      
      // check if the new element is really new
      if( C.find((*ne_it).first)!=C.end( ) || N.find((*ne_it).first)!=N.end( ) )
	continue;

      cout << "new_tuple : "; printon( cout , (*ne_it).first ) << endl;

      /*
      // DEBUG
      vector< ThRightNormalForm > elts2 = elts;
      for( int i=0 ; i<sz ; ++i ) {
	elts2[i]  = -(*ne_it).second * elts2[i] * (*ne_it).second;
	cout << shortBraidForm( rank , elts2[i].getShortWord( ) ) << "  ,  ";
      }
      cout << endl;
      
      if( elts2!=(*ne_it).first ) {
	cout << "Problem" << endl;
	exit(1);
      }
      */
      
      // Check if the new element has smaller component than the original
      bool increase = false;
      for( int i=0 ; i<sz ; ++i )
	if( cur_tuple[i].getPower( )<(*ne_it).first[i].getPower( ) )
	  increase = true;
      
      if( increase ) {
	cout << "================================" << endl;
	N.clear( );
	C.clear( );
	N[(*ne_it).first] = (*ne_it).second;
	break;
      }
      N[(*ne_it).first] = (*ne_it).second;
    }
  }
  
  pair< vector< ThRightNormalForm > , ThRightNormalForm > result;
  return result;
}


//---------------------------------------------------------------------------//
//----------------------- ThRightNormalFormAlgorithms -----------------------//
//---------------------------------------------------------------------------//


set< Permutation > getSimpleConjugators( int rank , const vector< ThRightNormalForm >& tuple )
{
  set< Permutation > result;
  for( int i=0 ; i<rank-1 ; ++i ) {
    Permutation c( rank );
    c.change( i , i+1 );
    result.insert( getSimpleConjugator( rank , tuple , c ) );
  }
  return result;
}


Permutation getSimpleConjugator( int rank , const vector< ThRightNormalForm >& tuple , const Permutation& start )
{
  Permutation c = start;
  int sz = tuple.size( );
  
  // cout << "start = " << start << endl;
  for( bool progress=true ; progress ; ) {
    progress = false;
    for( int i=0 ; i<sz ; ++i ) {
      Permutation new_c = ThRightNormalForm::getSimpleConjugator( rank , ThRightNormalForm::NF(tuple[i]) , c );
      if( new_c!=c ) {
	c = new_c;
	// cout << "    c = " << c << endl;
	progress = true;
      }
    }
  }
  
  // cout << "    r = " << c << endl;
  return c;
}
