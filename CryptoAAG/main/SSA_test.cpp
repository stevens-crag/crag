#include "Word.h"
#include "braid_group.h"
#include "ShortBraidForm.h"
#include "AAGKeyGeneration.h"
#include "braid.h"
#include "RanlibCPP.h"
#include <iostream>
#include "MajorDump.h"


#define USE_MPI


#ifdef USE_MPI
#include <mpi.h>
#endif


int main( int argc, char** argv)
{


  srand( time(0) );


  // Define parameters
  // parameters defined to be used in braid.h !!!
  // Very inconvinient should be changed
  extern int  STRANDS;

  int len = 10;
  int Total_N   = 100;
  int c   = 3;
  int k   = 6;
  int sg_conj_len = 0; //100;
  int key_len = 5;
  STRANDS = c*k;

  // initialize default rank to 1 (for no MPI execution)
  int my_rank = 1;


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
      int N_exp_send = int(Total_N / pn);
      if ( dest < Total_N % pn ) N_exp_send++;
      
      MPI_Send(&N_exp_send, 1, MPI_INTEGER, dest, tag, MPI_COMM_WORLD);

      
    }
  } 
  
  
  {
    
    int source = 0;
    long seed1 = 0;
    long seed2 = 0;
    int N_exp;
    
    MPI_Recv(&seed1, 1, MPI_INTEGER, source, 0,  MPI_COMM_WORLD, &status);
    MPI_Recv(&seed2, 1, MPI_INTEGER, source, 0, MPI_COMM_WORLD, &status);
    RandLib::ur.reset(seed1,seed2);



    MPI_Recv(&N_exp, 1, MPI_INTEGER, source, 0, MPI_COMM_WORLD, &status);
    
    // Init OUTPUT  
    stringstream f_name; 
    f_name << "Out" << time(NULL) << "_SSA"<<my_rank +1 <<'\0'<<flush;
    (Dump::dump_out) = new ofstream(f_name.str().c_str());
    

   *Dump::dump_out << "Rank - " << my_rank << " N - " << pn << " N exp - " << N_exp << endl;

    
    for (int exp_count=0;exp_count < N_exp; exp_count++){
      
      
#else
      
      // IF not multiprocess, then output to the console
      Dump::dump_out = &cout;
      
      for (int exp_count=0;exp_count < n;exp_count++){

#endif
	
	*(Dump::dump_out) << exp_count+1 << " =========================================================== " << endl;

	*(Dump::dump_out) << "Parameters  " << endl
	    << "                 c : " << c << endl
	    << "                 k : " << k << endl
	    << "                 N : " << STRANDS << endl
	    << "    Len of SG conj : " << sg_conj_len << endl
	    << "    Random key len : " <<  key_len << endl << endl;
	
	
	//
	// Generate AAG challenge
	//  
	
	*(Dump::dump_out) << "Generate AAG Instance ..." << endl;
	
	AAGProtocolInstance AAG = AAGProtocolInstance::challenge( STRANDS, sg_conj_len, c, k, key_len );
	// AAGProtocolInstance AAG = AAGProtocolInstance::random( STRANDS, 10 , 10 , 13 , 50 , 50 );
	

	vector< Word > Sbgp_B      = AAG.getBobPublicSbgp( );
	vector< Word > Sbgp_BConj  = AAG.getBobConjSbgp();
	Word           A_key       = AAG.getAliceKey( );
	
	*Dump::dump_out << "Finished generating AAG " << endl;
	
	AAG.printStats(*(Dump::dump_out) );
	
	
	//
	// Convert to cbraid format
	//
	
	*(Dump::dump_out) << "Convert to cbraid format ..." << endl;
	
	BraidType cbraid_A_key(STRANDS); craig2cbraid( A_key.getList() , cbraid_A_key ); cbraid_A_key.MakeLCF();
	cbraid_A_key.LeftDelta=0;
	*(Dump::dump_out) << "Key conversion complete." << endl;
	
	vector<BraidType> cbraid_Sbgp_B( Sbgp_B.size(),BraidType(STRANDS) );
	vector<BraidType> cbraid_Sbgp_BConj( Sbgp_BConj.size(),BraidType(STRANDS) );  
	for (int i=0;i<cbraid_Sbgp_B.size();i++){
	  BraidType tmp_B(STRANDS); craig2cbraid( Sbgp_B[i].getList() , tmp_B ); tmp_B.MakeLCF();
	  //    BraidType tmp_BConj(STRANDS); craig2cbraid( Sbgp_BConj[i].getList() , tmp_BConj ); tmp_BConj.MakeLCF();
	  
	  
	  //    BraidType aag_check  = !cbraid_A_key*tmp_B*cbraid_A_key; aag_check.MakeLCF();
	  //    if ( aag_check != tmp_BConj ){
	  //      cout << "Error in translation!!!" << endl;
	  //      exit(0);
	  //    } else {
	  //      cout << "G"<<i+1<< " conversion complete" << endl;
	  //    }
	  
	  cbraid_Sbgp_B[i] = tmp_B;
	  
	  cbraid_Sbgp_BConj[i] = (!cbraid_A_key)*tmp_B*cbraid_A_key;
	  cbraid_Sbgp_BConj[i].MakeLCF();
	  
	  //    cbraid_Sbgp_BConj[i] = tmp_BConj;
	  
	 *(Dump::dump_out) << i+1 << " " << flush;
	  
	}
	*(Dump::dump_out) << endl;
	
	*(Dump::dump_out) << "Start SSA ..." << endl; 
	BraidType mya(STRANDS); 
	BraidType tmp(STRANDS); 
	
	BraidType         a = cbraid_A_key;
	vector<BraidType> b = cbraid_Sbgp_B;
	vector<BraidType> x = cbraid_Sbgp_BConj;
	
	
	//    BraidType         a(STRANDS); // = cbraid_A_key;
	//    vector<BraidType> b; // = cbraid_Sbgp_B;
	//    vector<BraidType> x; // = cbraid_Sbgp_BConj;
	
	bool isgood; 
	
	
	// generate system instance and call attack
	//    cout << "generating AAFG system..." << endl;
	//    initAAFG(&x,&b,&a,10,10,50,false, true);
	*(Dump::dump_out) << "testing generated system..." << endl;
	for (int i=0;i<b.size();i++) {
	  tmp=(!a)*b[i]*a*(!x[i]); tmp.MakeLCF();
	  if (!(tmp.CompareWithIdentity())) {
	   *(Dump::dump_out) << "not correct!" << endl;
	  } else {
	   *(Dump::dump_out) << "G" << i+1 << " is good." << endl;
	  }
	}
	*(Dump::dump_out) << "attacking system..." << endl;
	mya=attackAAFG(&x,&b,&a,Dump::dump_out); // cf. comments above
	// mya=attackAAFG(&x,&b);
	
	
	// check if attack was successful
	isgood=true;
	for (int i=0;(i<b.size()) && isgood;i++) {
	  BraidType myx=(!mya)*b[i]*mya;
	  myx.MakeLCF();
	  isgood=isgood && (myx==x[i]);
	}
	if (isgood) {

	  *(Dump::dump_out) << "result: found good a :)" << endl;
	  *(Dump::dump_out) << "Bad instance information:" << endl
			    << "Key:" << endl
			    << A_key << endl
			    << " Subgroup generators : " << endl;
	  for (int ii=0;ii<Sbgp_B.size(); ii++)
	    *(Dump::dump_out) << Sbgp_B[ii] << endl;
	    
	    
	} else
	  *(Dump::dump_out) << "result: did not succeed :(" << endl;
	//	cout << "`difference' LCF(a*estimated_a^-1):" << endl;
	//    tmp=a*!mya; tmp.MakeLCF(); showBraid(&tmp);
	
      }
      
#ifdef USE_MPI
      delete Dump::dump_out;
      /*    /// SEND RESULTS TO PROCESS 0
	    int dest = 0;
	    int tag = 0;
	    MPI_Send(&success, 1, MPI_FLOAT, dest, tag, MPI_COMM_WORLD);
      */


  }
  
  
  /* Shut down MPI */
  MPI_Finalize();
#endif

  return 0;
}
  
