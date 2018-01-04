#include "Word.h"
#include <iostream>
#include <iterator>
#include <algorithm>
#include "AEProtocol.h"
#include "TTPAttack.h"
#include "RanlibCPP.h"
#include "ProgressBar.h"
#include "ThLeftNormalForm.h"
#include "BraidGroup.h"
#include "ShortBraidForm.h"
#include <time.h>

using namespace std;

int main() {
  RandLib::ur.reset();
  long s1, s2;
  RandLib::ur.getseed( s1, s2 );
  cout << "Seed : " << s1 << " " << s2 << endl;
  
  TTP_Conf ttp_conf;
  

  // AE suggested parameters
  ttp_conf.nBL    =  6; // 5;     // # Generators in BL
  ttp_conf.nBR    =  6; // 5;     // # Generators in BR
  ttp_conf.N      = 14; // 12;    // Group rank
  ttp_conf.nGamma = 10; // Tuple size
  
  ttp_conf.len_z  = 150; // 18;  // Conjugator's length
  ttp_conf.len_w  = 150;  // Word's length
  
  cout << ttp_conf << endl;


  int nExp = 100;
  
  ////////////////////////////////////////////////////////////////////////
  //
  //   TTP attack
  //
  ////////////////////////////////////////////////////////////////////////
  
  
  int suc_count = 0;


  for (int i=0;i<nExp;i++){
    
    BSets bs = BSets::generateEqual( ttp_conf.N );
    //BSets bs = BSets::generateRandom( ttp_conf.N );
    cout << "BS generated : " << bs << endl;
    TTPTuple dw = AEKeyExchange::generateTuples( ttp_conf, bs );

    
    TTPAttack ttp( ttp_conf.N, bs );
    //    if ( ttp.testTuples( dw ) )
    clock_t begin_t = clock();
    if ( ttp.run(dw) ) {
    	cout << "SUCC" << endl;
      	suc_count++;
    } else
    	cout << "FAIL" << endl;
    
  	cout << "TIME : " << double(clock() - begin_t)/ double(CLOCKS_PER_SEC) << endl;   
  }
  
  cout << endl << "Success : " << double(suc_count) / double(nExp) << endl;



  return 0;

}
