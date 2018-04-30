#include "Word.h"
#include <iostream>
#include <iterator>
#include <algorithm>
#include "AEProtocol.h"
#include "TTPAttack.h"
#include "RanlibCPP.h"
#include "ProgressBar.h"
#include "ThLeftNormalForm.h"
#include "braid_group.h"
#include "ShortBraidForm.h"
#include <time.h>
#include <fstream>
#include <boost/program_options.hpp>

using namespace std;
namespace po = boost::program_options;

void getAnticipatedDeltaPowerIdea() {
  for (int n = 5; n <= 30; n += 5) {
    crag::braidgroup::BraidGroup B(n);
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
  const crag::braidgroup::BraidGroup B(n);
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
  const crag::braidgroup::BraidGroup B(n);
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
  const crag::braidgroup::BraidGroup B(n);
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
  crag::braidgroup::BraidGroup B(n);

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

void prepareTableOfLengths() {
  vector<int> ranks = {8, 12, 16, 20};
  vector<int> lengths = {20, 50, 100, 200, 500, 1000};
  for (const auto& l : lengths) {
    for (const auto& N : ranks) {
      cout << "(" << N << "," << l << ") -> ";
      TTP_Conf ttp_conf;
      ttp_conf.nBL = N / 2 - 1; // 5;     // # Generators in BL
      ttp_conf.nBR = N / 2 - 1; // 5;     // # Generators in BR
      ttp_conf.N = N; // 12;    // Group rank
      ttp_conf.nGamma = 10; // Tuple size
      ttp_conf.len_z = l; // 18;  // Conjugator's length
      ttp_conf.len_w = l;  // Word's length

      long long int total_length = 0;
      const int nExp = 100;
      for (int i = 0; i < nExp; i++) {
        BSets bs = BSets::generateEqual(ttp_conf.N);
        TTPTuple dw = AEKeyExchange::generateTuples(ttp_conf, bs);
        TTPTuple init_tuple = dw.takeModuloDeltaSQ(N);
        total_length += init_tuple.length();
      }
      cout << total_length / nExp << endl;
    }
  }
  // Output:
  //(8, 20) -> 3986
  //(12, 20) -> 6244
  //(16, 20) -> 9573
  //(20, 20) -> 12586
  //(8, 50) -> 10005
  //(12, 50) -> 15842
  //(16, 50) -> 21968
  //(20, 50) -> 28576
  //(8, 100) -> 20258
  //(12, 100) -> 31485
  //(16, 100) -> 43089
  //(20, 100) -> 56229
  //(8, 200) -> 41361
  //(12, 200) -> 64291
  //(16, 200) -> 86824
  //(20, 200) -> 110392
  //(8, 500) -> 103100
  //(12, 500) -> 159551
  //(16, 500) -> 218255
  //(20, 500) -> 275611
  //(8, 1000) -> 206835
  //(12, 1000) -> 320829
  //(16, 1000) -> 438238
  //(20, 1000) -> 556195

  // 17708
  // 31392
  // 52102
  // 72553
  // 44449
  // 79648
  // 119565
  // 164729
  // 90000
  // 158295
  // 234520
  // 324138
  // 183754
  // 323233
  // 472556
  // 636368
  // 458041
  // 802167
  // 1187895
  // 1588793
  // 918904
  // 1613018
  // 2385195
  // 3206254
}

// const vector<double> len = {
//    (log(8) + 1) / log(2) * 3986,
//    (log(12) + 1) / log(2) * 6244,
//    (log(16) + 1) / log(2) * 9573,
//    (log(20) + 1) / log(2) * 12586,
//    (log(8) + 1) / log(2) * 10005,
//    (log(12) + 1) / log(2) * 15842,
//    (log(16) + 1) / log(2) * 21968,
//    (log(20) + 1) / log(2) * 28576,
//    (log(8) + 1) / log(2) * 20258,
//    (log(12) + 1) / log(2) * 31485,
//    (log(16) + 1) / log(2) * 43089,
//    (log(20) + 1) / log(2) * 56229,
//    (log(8) + 1) / log(2) * 41361,
//    (log(12) + 1) / log(2) * 64291,
//    (log(16) + 1) / log(2) * 86824,
//    (log(20) + 1) / log(2) * 110392,
//    (log(8) + 1) / log(2) * 103100,
//    (log(12) + 1) / log(2) * 159551,
//    (log(16) + 1) / log(2) * 218255,
//    (log(20) + 1) / log(2) * 275611,
//    (log(8) + 1) / log(2) * 206835,
//    (log(12) + 1) / log(2) * 320829,
//    (log(16) + 1) / log(2) * 438238,
//    (log(20) + 1) / log(2) * 556195,
//};
// for (const auto l : len) {
//  cout << static_cast<int>(l) << endl;
//}
// cout << endl;
// return 0;
// std::ios_base::sync_with_stdio(false);

int main(int argc, char* argv[]) {
  // Declare the supported options.
  po::options_description desc("Allowed options");
  desc.add_options()
      ("help", "produces help message")
      ("n", po::value<int>()->default_value(12), "n for B_n")
      ("z", po::value<int>()->default_value(100), "length of |z|")
      ("w", po::value<int>()->default_value(100), "length of |w|");

  po::variables_map vm;

  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  po::notify(vm);

  if (vm.count("help")) {
    std::cout << desc << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Using n = " << vm["n"].as<int>()
            << ", |z| = " << vm["z"].as<int>()
            << ", |w| = " << vm["w"].as<int>()
            << std::endl;

  RandLib::ur.reset();
  long s1, s2;
  RandLib::ur.getseed(s1, s2);
  cout << "Seed : " << s1 << " " << s2 << endl;

  // test_shortenBraid2();
  // getDecompositionDetails();
  // compareGenerateWordFunctions();
  // abelinizationTest();
  // prepareTableOfLengths();
  // return 0;

  TTP_Conf ttp_conf;

  // AE suggested parameters
  // ttp_conf.nBL    =  9; // 5;     // # Generators in BL
  // ttp_conf.nBR    =  9; // 5;     // # Generators in BR
  // ttp_conf.N      = ttp_conf.nBL + ttp_conf.nBR + 2; // 12;    // Group rank
  // ttp_conf.nGamma = 10; // Tuple size
  //
  // ttp_conf.len_z  = 500; // 18;  // Conjugator's length
  // ttp_conf.len_w  = 500;  // Word's length
  auto arg_n = vm["n"].as<int>();
  ttp_conf.nBL = arg_n / 2 - 1;                 // 5;     // # Generators in BL
  ttp_conf.nBR = arg_n / 2 - 1;                 // 5;     // # Generators in BR
  ttp_conf.N = ttp_conf.nBL + ttp_conf.nBR + 2; // 12;    // Group rank
  ttp_conf.nGamma = 10;                         // Tuple size

  ttp_conf.len_z = vm["z"].as<int>();               // 18;  // Conjugator's length
  ttp_conf.len_w = vm["w"].as<int>();               // Word's length

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
    BSets bs = BSets::generateEqual(ttp_conf.N);
    // BSets bs = BSets::generateRandom( ttp_conf.N );
    cout << "==============================================" << endl;
    cout << "Experiment #" << i + 1 << endl;
    cout << "BS generated : " << bs << endl;
    TTPTuple dw = AEKeyExchange::generateTuples(ttp_conf, bs);

    of << "Experiment #" << i + 1 << endl;

    TTPAttack ttp(ttp_conf.N, bs);
    clock_t begin_t = clock();
    if (ttp.run(dw)) {
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
