// Copyright (C) 2006 Aleksey Myasnikov
// Contents: Tests for AAG
//
// Principal Authors: Aleksey Myasnikov
//
// Revision History:
//

//#define MPI


#ifdef MPI
#include <mpi.h>
#endif


#include "AAGKeyGeneration.h"
#include <iostream>
#include <iterator>
#include <algorithm>
#include "ShortBraidForm.h"
#include "LengthAttack.h"

int main( )
{

  // Define parameters

  int N = 80;                             //  The rank of the braid group
  int num_gens = 20;                      //  Number of subgroup generators 
  int min_len = 20;                       //  Minimal length of subgroup generators
  int max_len = 23;                       //  Maximal length of subgroup generators
  int AliceDecompositionLength = 50;     //  
  int BobDecompositionLength = 50;


  // initialize default rank to 1 (for no MPI execution)
  int my_rank = 1;
  int length_incr = 1;

#ifdef MPI
 
 /// Initializa MPI

  MPI_Status status;      /* return status for receive */

  /* Initialize MPI */
  MPI_Init(&argc, &argv);

  /* Get my process rank */
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  /* How many processes are being used? */
  int pn;
  MPI_Comm_size(MPI_COMM_WORLD, &pn);

  /// End initialize MPI


 /// Init OUTPUT
  stringstream f_name;
  f_name << "Out_TEST_AAG_" << time(NULL) << ".p"<<my_rank<<'\0'<<flush;
  ofstream out(f_name.str().c_str());
  out << "Rank - " << my_rank << " N - " << pn <<endl;


  //  MAIN PROCESS 0
  if (my_rank == 0) {
    for (int source = 1; source < pn; source++) {
      int res;
      int tag = 0;
      MPI_Recv(&res, 1, MPI_FLOAT, source, tag,
	       MPI_COMM_WORLD, &status);
      out << "Recieved from  " << source << " result : " << res << endl;
    }
  } else {

#else
    
    // IF not multiprocess, then output to the console
    ostream&  out = cout;
#endif


  out << "Parameters:" << endl
	<< "     N = " << N << endl
	<< "     num_gens = " << num_gens << endl
	<< "     min_len = " << min_len + (my_rank-1)*length_incr << endl
	<< "     max_len = " << max_len + (my_rank-1)*length_incr << endl 
	<< "     Decomp  = "<< AliceDecompositionLength << endl
	<< endl;


  // Generate AAG instance
  AAGProtocolInstance AAG = AAGProtocolInstance::random( N , num_gens , min_len , max_len, 
							 AliceDecompositionLength, 
							 BobDecompositionLength );
  
  vector< Word > Sbgp_A  = AAG.getAlicePublicSbgp( );
  vector< Word > Sbgp_A2 = AAG.getAliceConjSbgp( );
  vector< Word > Sbgp_B  = AAG.getBobPublicSbgp( );
  vector< Word > Sbgp_B2 = AAG.getBobConjSbgp( );
  
  
  // Get keys
  Word keyA = AAG.getAliceKey ( );
  Word keyB = AAG.getBobKey( );
  Word key  = AAG.getSharedKey( );



  // construct conjugator
  Word decomposition = Word::randomWord(num_gens,BobDecompositionLength);
  
  int len = decomposition.length( );
  vector<Word> conjv(len);
  vector<int> gens(len);
  Word conj;
  int i=0;
  for (ConstWordIterator d_it = decomposition.begin( ); d_it!=decomposition.end( ) ; ++d_it,i++ ) {
    int ngen = *d_it;
    
    gens[i] = ngen;

    if( ngen<0 )
      conjv[i]= Sbgp_B[-ngen-1].inverse( );
    else
      conjv[i]= Sbgp_B[ngen-1];
    
    conj *= conjv[i];
    
  }
  conj = shortenBraid( N , conj );


  // construct the subgroup
  vector<Word> Sbgp_A3(Sbgp_A.size());
  for (int i=0;i<Sbgp_A.size();i++){
    Sbgp_A3[i] = shortenBraid( N, conj*Sbgp_A[i]*-conj );
  }

  for (int I=0;I<conjv.size();I++){
    out << " STEP : " << I+1 << " -------------------------------------------- " << endl;
    int old_len = 0;
    int new_len = 0;
    int old4_len = 0;
    int new4_len = 0;    
    stringstream s;


    int maxInc = -9999999;
    int maxInc4 =  0;
    int maxIncGen = 0;
    
  
    vector<Word> Sbgp_A3_sav = Sbgp_A3;

    
    for (int i=0;i<Sbgp_A3.size();i++){
      
      old_len += Sbgp_A3_sav[i].length();
      Sbgp_A3[i]  = shortenBraid( N, -conjv[I]*Sbgp_A3_sav[i]*conjv[I] );
      new_len += Sbgp_A3[i].length();
      s << Sbgp_A3_sav[i].length() - Sbgp_A3[i].length() << " ";
      if (i < 5){
	old4_len += Sbgp_A3_sav[i].length();
	new4_len += Sbgp_A3[i].length();      
      }
      
    }
    out << "TRUE : "<< " ( " << gens[I] << " ) "
	<< old_len - new_len << " " <<  old4_len - new4_len << " [ " << s.str() << " ] " << endl;
  
    //    if (old_len - new_len < 0)
      for (int bi=0;bi<Sbgp_B.size();bi++){
	
	int old_len1 = 0;
	int new_len1 = 0;
	int old_len2 = 0;
	int new_len2 = 0;
	int old4_len1 = 0;
	int new4_len1 = 0;
	int old4_len2 = 0;
	int new4_len2 = 0;
	int bb = 0;
	stringstream s1, s2;

	if (abs(gens[I]) != bi+1)
	for (int d=0;d<2;d++) {
	  
	  for (int i=0;i<Sbgp_A3_sav.size();i++){
	    if (d){
	      
	      old_len1 += Sbgp_A3_sav[i].length();
	      Word bw1  = shortenBraid( N, -Sbgp_B[bi]*Sbgp_A3_sav[i]*Sbgp_B[bi] );
	      new_len1 += bw1.length();
	      s1 << Sbgp_A3_sav[i].length() - bw1.length() << " ";
	      if (i < 5){
		old4_len1 += Sbgp_A3_sav[i].length();
		new4_len1 += bw1.length();      
	      }
	      
	    } else {
	      
	      old_len2 += Sbgp_A3_sav[i].length();
	      Word bw2  = shortenBraid( N, Sbgp_B[bi]*Sbgp_A3_sav[i]*-Sbgp_B[bi] );
	      new_len2 += bw2.length();
	      s2 << Sbgp_A3_sav[i].length() - bw2.length() << " ";
	      if (i < 5){
		old4_len2 += Sbgp_A3_sav[i].length();
		new4_len2 += bw2.length();      
	      }	    
	    }
	  }
	}
	
	if (max(old_len1 - new_len1, old_len2 - new_len2) > 0)
	  if (old_len1 - new_len1 > old_len2 - new_len2){
	    out << "FALSE : "<< " abs( " << bi+1 << " ) "
		<< old_len1 - new_len1 << " " <<  old4_len1 - new4_len1 << " [ " << s1.str() << " ] " << endl;
	    if (abs(gens[I]) != bi+1) {
	      if(maxInc < old_len1 - new_len1){
		maxInc =  old_len1 - new_len1;
		maxInc4 =  old4_len1 - new4_len1;
		maxIncGen = bi+1;
	      }
	    }
	  } else {
	    out << "FALSE : "<< " abs( " << -(bi+1) << " ) "
		<< old_len2 - new_len2 << " " <<  old4_len2 - new4_len2 << " [ " << s2.str() << " ] " << endl;
	    if (abs(gens[I]) != bi+1) {
	      if(maxInc < old_len2 - new_len2){
		maxInc =  old_len2 - new_len2;
		maxInc4 =  old4_len2 - new4_len2;
		maxIncGen = bi+1;
	      }
	    }
	  }
	//	cout << Sbgp_B[bi] << endl;
	
    }
    
    out << "MAX : "<< " abs( " << maxIncGen << " ) " << maxInc << " " << maxInc4 << endl;
  }
  
  copy(Sbgp_A3.begin(),Sbgp_A3.end(),ostream_iterator<Word>(out,"\n"));
  out << "------------------>" << endl;
  copy(Sbgp_A.begin(),Sbgp_A.end(),ostream_iterator<Word>(out,"\n"));
  
  for (int i=0;i<Sbgp_A3.size();i++)
    out << shortenBraid( N, -Sbgp_A[i]*Sbgp_A3[i] ) << endl;
  
  // findKey_LengthBased( N , Sbgp_A , Sbgp_A3 , Sbgp_B , 48000, cout );



#ifdef MPI
    /// SEND RESULTS TO PROCESS 0
    int dest = 0;
    int tag = 0;
    MPI_Send(&success, 1, MPI_FLOAT, dest, tag, MPI_COMM_WORLD);
    
  }

  out.close();

  /* Shut down MPI */
  MPI_Finalize();
#endif


  return 0;
       
}
