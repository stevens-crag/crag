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
#include <fstream>

using namespace std;

void getAnticipatedDeltaPowerIdea() {
  for (int n = 5; n <= 30; n += 5) {
    BraidGroup B(n);
    for (int l = 40; l <= 800; l += 40) {
      int total_delta = 0;
      int total_dec = 0;
      for (int e = 0; e < 10; e++) {
        const auto w = Word::randomWord(n - 1, l);
        ThLeftNormalForm nf(B, w);
        total_delta += nf.getPower();
        total_dec += nf.getDecomposition().size();
      }
      cout << "(" << total_delta << "," << total_dec << ") ";
    }
    cout << endl;
  }
}

static bool wayToSort(int i, int j) { return i > j; }

void getDecompositionDetails() {
  const int n = 16;
  const int l = 1000;
  const BraidGroup B(n);
  const auto w = Word::randomWord(n - 1, l);
  ThLeftNormalForm nf(B, w);
  vector<int> sizes;
  for (const auto &perm : nf.getDecomposition()) {
    cout << perm.geodesic().size() << "-";
    sizes.push_back(perm.geodesic().size());
  }
  cout << endl;
  const auto p = -nf.getPower();
  cout << "Power = " << p << endl;
  sort(sizes.begin(), sizes.end(), wayToSort);
  for (const auto &s : sizes) {
    cout << s << "-";
  }
  cout << endl;
  cout << "Last: " << sizes[p - 1] << endl;
}

void compareGenerateWordFunctions() {
  const int n = 16;
  const int l = 10000;
  const BraidGroup B(n);
  for (int i = 0; i < 100; ++i) {
    const auto w = Word::randomWord(n - 1, l);
    ThLeftNormalForm nf(B, w);
    const auto w1 = nf.getWord();
    const auto w2 = nf.getReducedWord();
    const auto w3 = nf.getReducedWord2();
    cout << w1.length() << " vs " << w2.length() << " vs " << w3.length() << endl;
    ThLeftNormalForm nf3(B, w3);
    if (nf != nf3) {
      exit(1);
    }
  }

}

static int abelinization(const Word& w) {
  int result = 0;
  for (const auto &g : w) {
    result += g > 0 ? 1 : -1;
  }
  return result;
}

void abelinizationTest() {
  const int n = 16;
  const int l = 2000;
  const BraidGroup B(n);
  for (int i = 0; i < 1000; ++i) {
    const auto w = Word::randomWord(n - 1, l);
    cout << abelinization(w) << ", ";
  }
  cout << endl;
  cout << 2 * Permutation::getHalfTwistPermutation(n).geodesic().size() << endl;
}

void test_shortenBraid2() {
  const int n = 16;
  typedef ThLeftNormalForm NF;
  BraidGroup B(n);

  const int l = 200;
  for (int i = 0; i < 100; ++i) {
    const auto w1 = Word::randomWord(n - 1, l);
    const auto w2 = shortenBraid(n, w1);
    const auto w3 = shortenBraid2(n, w1);
    cout << w1.length() << " -> " << w2.length() << " -> " << w3.length() << endl;
    NF nf1(B, w1);
    NF nf3(B, w3);
    if (nf1 != nf3) {
      cout << "Failure!" << endl;
      exit(1);
    }
  }
}

int main() {
  RandLib::ur.reset();
  long s1, s2;
  RandLib::ur.getseed( s1, s2 );
  cout << "Seed : " << s1 << " " << s2 << endl;
  
  // test_shortenBraid2();
  // getDecompositionDetails();
  // compareGenerateWordFunctions();
  // abelinizationTest();
  // return 0;

  TTP_Conf ttp_conf;
  
  // AE suggested parameters
  ttp_conf.nBL    =  7; // 5;     // # Generators in BL
  ttp_conf.nBR    =  7; // 5;     // # Generators in BR
  ttp_conf.N      = ttp_conf.nBL + ttp_conf.nBR + 2; // 12;    // Group rank
  ttp_conf.nGamma = 10; // Tuple size
  
  ttp_conf.len_z  = 600; // 18;  // Conjugator's length
  ttp_conf.len_w  = 600;  // Word's length
  
  cout << ttp_conf << endl;


  const int nExp = 100;
  
  ////////////////////////////////////////////////////////////////////////
  //
  //   TTP attack
  //
  ////////////////////////////////////////////////////////////////////////
  
  
  double suc_count = 0;

  fstream of("experiments.txt", ios::app);
  for (int i = 0; i < nExp; i++) {
    BSets bs = BSets::generateEqual( ttp_conf.N );
    //BSets bs = BSets::generateRandom( ttp_conf.N );
    cout << "==============================================" << endl;
    cout << "Experiment #" << i + 1 << endl;
    cout << "BS generated : " << bs << endl;
    TTPTuple dw = AEKeyExchange::generateTuples(ttp_conf, bs);

    of << "Experiment #" << i + 1 << endl;

    TTPAttack ttp(ttp_conf.N, bs);
    clock_t begin_t = clock();
    if ( ttp.run(dw) ) {
    	cout << "SUCCESS" << endl;
      	suc_count++;
    } else
    	cout << "FAIL" << endl;

    cout << "TIME : " << double(clock() - begin_t) / double(CLOCKS_PER_SEC) << endl;
    cout << endl << " >>> Success: " << suc_count << " out of " << (i + 1) << " <<< " << endl;
    of << "TIME : " << double(clock() - begin_t) / double(CLOCKS_PER_SEC) << endl;
    of << endl << " >>> Success: " << suc_count << " out of " << (i + 1) << " <<< " << endl;
  }

  return 0;

}
