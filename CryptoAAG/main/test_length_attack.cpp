
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include "RanlibCPP.h"

#include "braid_group.h"
#include "ThRightNormalForm.h"

#include <sstream>
#include <iterator>
#include <iostream>
#include <fstream>


//
//  Uncomment to use MPI parallel computing interface
//
//#define USE_MPI


#ifdef USE_MPI
#include <mpi.h>
#endif

using namespace std;

#include "LengthAttack.h"
#include "AAGKeyGeneration.h"


//---------------------------------------------------------------------------//
//---------------------------------- main -----------------------------------//
//---------------------------------------------------------------------------//


int main( int argc, char** argv )
{

  srand( time(0) );

  // Define parameters
  int N = 80;           // Rank of Braid group
  int num_gens = 20;    // Number of SG generators
  int min_len = 20;     // Minimal generator length
  int max_len = 23;     // Maximal generator length
  int AliceDecompositionLength = 50; // Length of the Alice's secret
  int BobDecompositionLength = 50;   // Length of the Bob's secret

  // number of samples
  int n = 100;

  // Set time limit
  int sec_in_day = 86400;
  int TIME_LIMIT_SEC = 1*sec_in_day;
 
  // initialize default rank to 1 (for no MPI execution)
  int my_rank = 1;
  int length_incr = 0; // 1;


#ifdef USE_MPI
 
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
  
  
  //  MAIN PROCESS 0
  if (my_rank == 0) {
    
    // send info to the  processors 
    for (int dest = 0; dest < pn; dest++) {
      int tag = 0;
      
      // send random seeds
      long seed1 = RandLib::ur.irand(0,32000);
      long seed2 = RandLib::ur.irand(0,32000);
      
      MPI_Send(&seed1, 1, MPI_INTEGER, dest, tag, MPI_COMM_WORLD);
      MPI_Send(&seed2, 1, MPI_INTEGER, dest, tag, MPI_COMM_WORLD);
      
      
      // send number of samples to compute
      int k = int(n / pn);
      if ( dest < n % pn ) k++;
      
      MPI_Send(&k, 1, MPI_INTEGER, dest, tag, MPI_COMM_WORLD);

      
    }
  } 
  
  
  {
    
    int source = 0;
    long seed1 = 0;
    long seed2 = 0;
    int k;
    
    MPI_Recv(&seed1, 1, MPI_INTEGER, source, 0,  MPI_COMM_WORLD, &status);
    MPI_Recv(&seed2, 1, MPI_INTEGER, source, 0, MPI_COMM_WORLD, &status);
    RandLib::ur.reset(seed1,seed2);


    MPI_Recv(&k, 1, MPI_INTEGER, source, 0, MPI_COMM_WORLD, &status);
    
    for (int sample_count=0;sample_count < k;sample_count++){
      // Init OUTPUT  
      stringstream f_name; 
      f_name << "Out" << time(NULL) << "_P"<<my_rank +1 << "_S" << sample_count+1<<'\0'<<flush;
      ofstream out(f_name.str().c_str());
      out << "Rank - " << my_rank << " N - " << pn << " Sample #: " << sample_count <<endl;
      
      
#else
      
      // IF not multiprocess, then output to the console
      ostream&  out = cout;
#endif

      // Construct the attack 
      LengthAttack_A2 A;
      out << "Executing algorithm : A" << A.type() << endl;

      long s1, s2;
      RandLib::ur.getseed(s1,s2);

      out << "Using random seed : "<< s1 << " " << s2 << endl;

      out << "Parameters:" << endl
	  << "     N = " << N << endl
	  << "     num_gens = " << num_gens << endl
	  << "     min_len = " << min_len + (my_rank-1)*length_incr << endl
	  << "     max_len = " << max_len + (my_rank-1)*length_incr << endl 
	  << "     decomp  = " << AliceDecompositionLength << endl << endl;
      
      out << "Start generation of keys ..." << endl;
      
      // Generate AAG instance
      AAGProtocolInstance AAG = AAGProtocolInstance::random( N , num_gens , min_len + (my_rank-1)*length_incr , 
							     max_len + (my_rank-1)*length_incr, 
							     AliceDecompositionLength , BobDecompositionLength );
      
      vector< Word > Sbgp_A = AAG.getAlicePublicSbgp( );
      vector< Word > Sbgp_A2 = AAG.getAliceConjSbgp();
      vector< Word > Sbgp_B = AAG.getBobPublicSbgp( );
      
      // output
      // cout << "=============================================" << endl;
      // copy( Sbgp_A.begin( ) , Sbgp_A.end( ) , ostream_iterator< Word >( cout , "\n" ) );
      // cout << "=============================================" << endl;
      // copy( Sbgp_A2.begin( ) , Sbgp_A2.end( ) , ostream_iterator< Word >( cout , "\n" ) );
      // cout << "=============================================" << endl;
      // copy( Sbgp_B.begin( ) , Sbgp_B.end( ) , ostream_iterator< Word >( cout , "\n" ) );
      
      // try to find Bob' key
      out << "Start attack ... " << endl;
      int success = 999;
      switch( A.findKey_LengthBased( N , Sbgp_A , Sbgp_A2 , Sbgp_B , TIME_LIMIT_SEC, out ) ) {
      case SUCCESSFULL:
	out << "Success" << endl;
	success = 1;
	break;
      case FAILED:
	out << "Failed" << endl;
	success = 0;
	break;
      case TIME_EXPIRED:
	out << "Time expired" << endl;
	success = -1;
	break;
      }
      
#ifdef USE_MPI
      /*    /// SEND RESULTS TO PROCESS 0
	    int dest = 0;
	    int tag = 0;
	    MPI_Send(&success, 1, MPI_FLOAT, dest, tag, MPI_COMM_WORLD);
      */
      out.close();
    }
  }
  
  
  /* Shut down MPI */
  MPI_Finalize();
#endif
  
  
  return 0;
}
