// Copyright (C) 2006 Aleksey Myasnikov
// Contents: Tests data to test pick correlation
//
// Principal Authors: Aleksey Myasnikov
//
// Revision History:
//



#include "AAGKeyGeneration.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <iterator>
#include <algorithm>
#include "ShortBraidForm.h"


//#define FIXED_B
#define FIXED_A

int main( )
{

  // Define parameters

  int N = 80;                             //  The rank of the braid group
  int num_gens = 20;                      //  Number of subgroup generators 
  int min_len = 5;                       //  Minimal length of subgroup generators
  int max_len = 8;                       //  Maximal length of subgroup generators
  int AliceDecompositionLength = 50;     //  
  int BobDecompositionLength = 50;


  // initialize default rank to 1 (for no MPI execution)
  int my_rank = 1;
  int length_incr = 1;


	

  cout << "Parameters:" << endl
	<< "     N = " << N << endl
	<< "     num_gens = " << num_gens << endl
	<< "     min_len = " << min_len + (my_rank-1)*length_incr << endl
	<< "     max_len = " << max_len + (my_rank-1)*length_incr << endl 
	<< "     Decomp  = "<< AliceDecompositionLength << endl
	<< endl;
   
 
  vector< Word > Sbgp_B;
  vector< Word > Sbgp_A;

  bool first_time = true;
   
  for (int test_count=0;test_count<100;test_count++) { 
	
	

  // Generate AAG instance
  AAGProtocolInstance AAG = AAGProtocolInstance::random( N , num_gens , min_len , max_len, 
							 							AliceDecompositionLength, 
							 							BobDecompositionLength );

#ifdef FIXED_A  
	if ( first_time )
#endif
	Sbgp_A  = AAG.getAlicePublicSbgp( );

//  vector< Word > Sbgp_A2 = AAG.getAliceConjSbgp( );
  
#ifdef FIXED_B
	if ( first_time )
#endif
	Sbgp_B  = AAG.getBobPublicSbgp( );
  
  first_time = false;
  
  
//  vector< Word > Sbgp_B2 = AAG.getBobConjSbgp( );
  
  
  // Get keys
//  Word keyA = AAG.getAliceKey ( );
//  Word keyB = AAG.getBobKey( );
//  Word key  = AAG.getSharedKey( );




  // construct conjugator
  Word decomposition = Word::randomWord(num_gens,BobDecompositionLength);


	stringstream ss;
	ss << "PICK_CORR_Data_" << time(0) << "_" << test_count+1 << flush;
	cout << "Processing " << ss.str() << endl;
	
   
   ofstream out(ss.str().c_str());  
 
  Word conj;
  int len_A = 0;
  for (int i=0;i<Sbgp_A.size();i++)
  	len_A += Sbgp_A[i].length();
  
  out << conj.length() << " " << len_A << endl;
  
  for (ConstWordIterator d_it = decomposition.begin( ); d_it!=decomposition.end( ) ; ++d_it ) {
    int ngen = *d_it;
    
	Word conj_tmp;
    if( ngen<0 )
      conj_tmp = Sbgp_B[-ngen-1].inverse( );
    else
      conj_tmp = Sbgp_B[ngen-1];

    conj *= shortenBraid( N , conj_tmp );
    
    // construct the subgroup
    int len_A3 = 0;
    for (int i=0;i<Sbgp_A.size();i++){
      Word tmp_A = shortenBraid( N, conj*Sbgp_A[i]*-conj );
      len_A3 += tmp_A.length();
    }
	out << conj.length() << " " << len_A3 << endl;
    
  	}
 	out.close();
  }

   return 0;
       
}
