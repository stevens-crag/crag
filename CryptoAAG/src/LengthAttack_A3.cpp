// Copyright (C) 2005 Alexander Ushakov
// Contents: Implementation of length attack
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//

// A3

#include <set>
#include "LengthAttack.h"
#include "ShortBraidForm.h"
#include "FormatOutput.h"
#include <time.h>


void LengthAttack_A3::addProducts(  const vector<Word>& elem_set, vector<Word>& ext_set, vector<Word>& ext_set_sg_gens, const Word& sel_gen, int sel_gen_sg )
{
  
  const Word& M = sel_gen;
  for (int i=0;i<elem_set.size();i++)
    if ( (i+1)!= abs(sel_gen_sg) ){
      // add M*x
      ext_set_sg_gens.push_back(Word(sel_gen_sg)*Word(i+1));
      ext_set.push_back(M*elem_set[i]);
      // add M*-x
      ext_set_sg_gens.push_back(Word(sel_gen_sg)*Word(-(i+1)));
      ext_set.push_back(M*-elem_set[i]);
      // add x*M
      ext_set_sg_gens.push_back(Word(i+1)*Word(sel_gen_sg));
      ext_set.push_back(elem_set[i]*M);
      // add -x*M
      ext_set_sg_gens.push_back(Word(-(i+1))*Word(sel_gen_sg));
      ext_set.push_back(-elem_set[i]*M);      
      // add -x*M*x
      ext_set_sg_gens.push_back(Word(-(i+1))*Word(sel_gen_sg)*Word(i+1) );
      ext_set.push_back(-elem_set[i]*M*elem_set[i]);      
      // add x*M*-x
      ext_set_sg_gens.push_back(Word(i+1)*Word(sel_gen_sg)*Word(-(i+1)) );
      ext_set.push_back(elem_set[i]*M*-elem_set[i]);          
    } else {
      // add M*M
      ext_set_sg_gens.push_back(Word(sel_gen_sg)*Word(sel_gen_sg));
      ext_set.push_back(M*M);
      // add -M*-M
      ext_set_sg_gens.push_back(Word(-sel_gen_sg)*Word(-sel_gen_sg));
      ext_set.push_back(-M*-M);      
    }
  
}

void LengthAttack_A3::addAllProducts(  const vector<Word>& elem_set, vector<Word>& ext_set, vector<Word>& ext_set_sg_gens )
{
  
  
  for (int j=0;j<elem_set.size();j++){

    Word M = elem_set[j];
    int sel_gen_sg = j+1;
    
    for (int i=0;i<elem_set.size();i++){
      if ( (i+1)!= abs(sel_gen_sg) ){
	// add M*x
	ext_set_sg_gens.push_back(Word(sel_gen_sg)*Word(i+1));
	ext_set.push_back(M*elem_set[i]);
	// add M*-x
	ext_set_sg_gens.push_back(Word(sel_gen_sg)*Word(-(i+1)));
	ext_set.push_back(M*-elem_set[i]);
	// add -M*x
	ext_set_sg_gens.push_back(Word(-sel_gen_sg)*Word(i+1));
	ext_set.push_back(-M*elem_set[i]);
	// add -M*-x
	ext_set_sg_gens.push_back(Word(-sel_gen_sg)*Word(-(i+1)));
	ext_set.push_back(-M*-elem_set[i]);      
	// add -x*M*x
	ext_set_sg_gens.push_back(Word(-(i+1))*Word(sel_gen_sg)*Word(i+1) );
	ext_set.push_back(-elem_set[i]*M*elem_set[i]);      
	// add x*M*-x
	ext_set_sg_gens.push_back(Word(i+1)*Word(sel_gen_sg)*Word(-(i+1)) );
	ext_set.push_back(elem_set[i]*M*-elem_set[i]);          
	// add -x*-M*x
	ext_set_sg_gens.push_back(Word(-(i+1))*Word(-sel_gen_sg)*Word(i+1) );
	ext_set.push_back(-elem_set[i]*-M*elem_set[i]);      
	// add x*-M*-x
	ext_set_sg_gens.push_back(Word(i+1)*Word(-sel_gen_sg)*Word(-(i+1)) );
	ext_set.push_back(elem_set[i]*-M*-elem_set[i]);      
      } else {
	// add M*M
	ext_set_sg_gens.push_back(Word(sel_gen_sg)*Word(sel_gen_sg));
	ext_set.push_back(M*M);
      }
    }
    
  }
}


int LengthAttack_A3::sbgpGeneratorsWeight( const vector<Word>& A )
{
  int result = 0;
  for( int i=0 ; i<A.size() ; ++i )
    result += A[i].length( );
  return result;
}

void LengthAttack_A3::addNewElt( const vector<Word>& A , set< ELT >& checkedElements , set< ELT >& uncheckedElements )
{
  int weight = sbgpGeneratorsWeight( A );
  ELT new_elt( weight , A );
  
  if(   checkedElements.find( new_elt )!=  checkedElements.end( ) ) return;
  if( uncheckedElements.find( new_elt )!=uncheckedElements.end( ) ) return;
  
  uncheckedElements.insert( new_elt );
}


void LengthAttack_A3::tryElt( int N , const ELT& cur , const vector< Word >& B , const vector<Word>& B_sg_gens,
	     set< ELT >& checkedElements , set< ELT >& uncheckedElements, 
	     bool is_B_extended,
	     ostream& out )
{
  // we better vary this value, depending on parameters of A and B
  int MAX_DELTA = 0;

  int max_delta_observed = 0;
  int max_neg_criteria = 0;
  int max_delta_sg_gen = 0;
  Word max_delta_gen;


  int n_of_conj = 2;
  if (is_B_extended) n_of_conj = 1; // skip b*A[i]*-b when applying extended set of transformations

  for( int i=0 ; i<B.size( ) ; ++i ) {
    Word b = B[i];
    for( int d=0 ; d<n_of_conj ; ++d ) {
    
      //cout << "   #" << i+1 << "," << d+1;
      if( d==1 ) b = -b;
      bool candidate = true;
      int delta = 0;
      int neg_sum = 0;
      
      out << endl << "Try " <<  (( d==1 ) ? "-" : "") <<  "B_" << i+1
	  << " ( " << (( d==1 ) ? -B_sg_gens[i] : B_sg_gens[i]) << " ) : [ ";
        

      vector< Word > A = cur.second;
      for( int t=0 ; t<A.size( ) && candidate ; ++t ) {
	int old_len = A[t].length( );
	delta -= A[t].length( );
	A[t] = shortenBraid( N , -b*A[t]*b );
	delta += A[t].length( );
	
	if (old_len <=  A[t].length( )) neg_sum++;

        out << old_len - A[t].length() << " ";

	if( delta>MAX_DELTA && t>=4 ) {
	  candidate = false;
	  // out << " stopped @";
	}
      }
      
      
      //      if ( delta < 0 ) {
      //	out << endl << "Try " <<  (( d==1 ) ? "-" : "") <<  "B_" << i+1 
      //	    << " ( " << (( d==1 ) ? -B_sg_gens[i] : B_sg_gens[i]) << " ) : ";
      out << "] d = " << -delta << " ";
      //      }

      //      cout << i << " " << d << " " << max_delta_observed << " ---> ";
      if ( max_delta_observed > delta ){
	max_neg_criteria = neg_sum;
	max_delta_observed = delta;
	max_delta_gen = b;
	max_delta_sg_gen = ( d==1 ) ? -(i+1) : i+1;
      }

      //      cout<< max_delta_observed <<  endl;

      int weight = delta;
      if( candidate ) {
	addNewElt( A, checkedElements , uncheckedElements );
	//cout << " accepted with " << delta << endl;
      }
    }
  }

  
  out << endl << -max_delta_observed << " " << max_neg_criteria << " ";

  // max_delta_observed == 0 if no decreases of length occured
  // in this case do all extended products and conjugations
  if (max_delta_observed != 0 ) {
    out << max_delta_sg_gen << " ( " 
	<< ((max_delta_sg_gen < 0) ? -B_sg_gens[abs(max_delta_sg_gen)-1] : B_sg_gens[abs(max_delta_sg_gen)-1]) 
	<< " ) ";

    // Do extended set of transformations (products and conjugations with the maximal)
    if ( !is_B_extended && max_delta_observed < 0 ) { // > -200 
      out << endl << "Try extended set ..." << endl;
      vector<Word> B_ext;
      vector<Word> B_ext_sg_gens;
      
      Word& M = max_delta_gen;

      addProducts(  B, B_ext, B_ext_sg_gens, M, max_delta_sg_gen );
      
      tryElt( N , cur , B_ext , B_ext_sg_gens, checkedElements , uncheckedElements, true, out );
    }
    
  } else {
    
    // Do extended set of transformations (products and conjugations with all)
    if ( !is_B_extended ) { // > -200 
      out << endl << "Try all products ..." << endl;
      vector<Word> B_ext;
      vector<Word> B_ext_sg_gens;

      addAllProducts(  B, B_ext, B_ext_sg_gens );
      
      tryElt( N , cur , B_ext , B_ext_sg_gens, checkedElements , uncheckedElements, true, out );
    }

  }
  
  out << " <-----------------------------------------------------" << endl;
  
}


bool LengthAttack_A3::check_ifVectorsEqual( int N , const vector< Word >& A1 , const vector< Word >& A2 )
{
  int len = A1.size( );
  for( int i=0 ; i<len ; ++i ) {
    Word s = shortenBraid( N , A1[i]*-A2[i]  );
    if( s.length( )>0 )
      return false;
  }

  return true;
}


findKey_LengthBasedResult LengthAttack_A3::findKey_LengthBased( int N , const vector< Word >& A1 , const vector< Word >& A2 , const vector< Word >& B , int sec, ostream& out )
{
  set< ELT > checkedElements;
  set< ELT > uncheckedElements;

  int init_time = time( 0 );

  int init_weight1 =  sbgpGeneratorsWeight( A1 );
  int init_weight2 =  sbgpGeneratorsWeight( A2 );
  ELT init( init_weight2 , A2 );
  uncheckedElements.insert( init );
  out << "Initial weights: " << init_weight1 << ", " << init_weight2 << endl;
  int best_result = 999999;

  vector<Word> B_sg_gens(B.size());
  for (int i=0;i<B.size();i++)
    B_sg_gens[i] = Word(i+1);

  for( int c=0 ; uncheckedElements.size( ) && c<100000 ; ++c ) {
    
    ELT cur = *uncheckedElements.begin( );
    
    int cur_time = time( 0 );
    out << "Current (best) weight: " << cur.first << " (" << best_result << ")" << ", tm = " << cur_time << endl;
    for( vector< Word >::const_iterator c_it=cur.second.begin( ) ; c_it!=cur.second.end( ) ; ++c_it )
      out << (*c_it).length() << "  ";
    out << endl;
    
    
    uncheckedElements.erase( uncheckedElements.begin( ) );
    checkedElements.insert( cur );

    
    if( best_result>cur.first ) 
      best_result = cur.first;
    
    // Check time
    if( cur_time-init_time > sec )
      return TIME_EXPIRED;

    if( check_ifVectorsEqual( N , A1 , cur.second ) )
      return SUCCESSFULL;
    
    tryElt( N , cur , B , B_sg_gens, checkedElements , uncheckedElements, false, out );
  }
  
  return FAILED;
}
