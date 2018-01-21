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

void test_special_instance1() {
  // 1. Define an instance
  const int N = 8;
  const vector<Word> gens = { "x1"_w, "x1^-1"_w, "x2"_w, "x2^-1"_w, "x3"_w, "x3^-1"_w, "x4"_w, "x4^-1"_w, "x5"_w, "x5^-1"_w, "x6"_w, "x6^-1"_w, "x7"_w, "x7^-1"_w, };
  BSets BS;
  BS.BL = { "x1"_w, "x2"_w, "x3"_w, };
  BS.BR = { "x5"_w, "x6"_w, "x7"_w, };

  TTPTuple T;
  T.WL = { "x5 x6 x5^-1 x6^2 x7 x5^-1 x6^-1 x7 x5 x6^-1 x7 x5^-1 x6^-3 x5 x6^-1 x7 x5 x6^4 x7 x5^-1 x6 x5^-2 x6^2 x7 x6^-2 x7 x6^6 x7 x5^-1 x6^2 x7 x5^-2 x6^5 x7 x5 x6^-1 x7 x5^-1 x6^-1 x7 x6^2 x7 x5^-1 x6^2 x7 x5^-2 x6^2 x5^-1 x6 x7^2 x6^-3 x7 x5^-1 x6^-2 x7 x5^2 x6^-5 x7^2 x5^-3 x6^-2 x5 x6^-1 x7 x5^-1 x6 x7^2 x6^2 x7 x6^-1 x7^2 x5^-1 x6^-1 x7^3 x5^3 x6^-1 x7^2 x5^-1 x6^-1 x5 x6^-1 x7 x5^-1 x6^-1 x7 x5 x6 x5^2 x6^4 x5^-2 x6 x5^-1 x6^4 x5^-1 x6 x5 x6 x7^3 x6^-1 x5 x6^-5 x7 x5^-1 x6^-1 x5 x6^-1 x7^3 x6^-1 x7 x6 x7 x5^-1 x6^-1 x7 x5 x6^-1 x7 x5^-1 x6^-1 x5 x6^-1 x7^2 x5^-1 x6^-2 x5 x6^-1 x7 x6^-3 x7 x6^-1 x7 x5^-1 x6^-1 x7^2 x5^2 x6^-1 x7^2 x6^-5 x5^4 x6^-1 x7^2 x6 x5^-1 x6 x7^2 x5^-2 x6^-2 x7^2 x5^2 x6^-1 x7 x6 x5^-1 x6 x5^-1 x6 x5^-1 x6 x5^-2 x6 x7^2 x6^2 x7 x5^-1 x6^2 x5^-1 x6 x7 x6^-2 x5^-3 x6^-1 x5 x6^-1 x5 x6^-1 x7 x6 x5^-1 x6 x7 x5^-1 x6 x5^-1 x6 x7 x6^-1 x7 x6 x5^-1 x6^2 x7 x6^-1 x7 x6^-1 x7 x5^-1 x6 x5^-2 x6^2 x7 x6^-2 x7^2 x6 x5"_w, "x5 x6^-1 x7 x5^-1 x6^-1 x7 x5^-1 x6^-1 x7^2 x5 x6^-2 x5^2 x6^-1 x5 x6^-1 x7 x5^-3 x6^2 x7 x6^-1 x7 x6^3 x7^2 x5^-2 x6 x7 x6 x7 x6^-2 x7 x5^-1 x6^-1 x5^3 x6^-2 x5 x6^-1 x7 x5^-4 x6^-1 x5^-1 x6^-1 x7 x5^-1 x6^-1 x5 x6^-3 x5^-2 x6^-1 x7 x5^-1 x6^-1 x7 x6^-3 x7 x5 x6^-4 x7 x5^-1 x6^-1 x5 x6^-3 x7 x5^-1 x6^2 x7^2 x6^-1 x7^2 x5 x6^-1 x7 x5^-1 x6 x7 x6^-2 x7 x5 x6^-2 x7 x5^-1 x6 x7 x6^-2 x5^-1 x6^-2 x5^-2 x6^-5 x7 x5^-1 x6 x5^-1 x6 x7 x6^-1 x5 x6^-1 x7^3 x5^-1 x6 x5^-1 x6 x7 x6 x7 x5^2 x6^-2 x5^4 x6^-1 x7^2 x6^-1 x7 x6^-2 x7 x5^-1 x6^3 x7 x5^-1 x6^-2 x7^3 x6 x5 x6 x5 x6 x5 x6^6 x5 x6 x7^2 x6^-2 x7^2 x6 x7 x5^-1 x6^-3 x5 x6^-1 x7^3 x6 x7 x6 x7^2 x5^-3 x6^-2 x7 x5 x6^-1 x7^3 x6 x7 x6 x5^-1 x6 x7 x5^-1 x6^-2 x7 x5 x6^-1 x7 x6^-1 x7 x6 x7 x5^-1 x6^-1 x7 x5 x6^3 x7^3 x5^-1 x6^-1 x5 x6^-1 x7^2 x5^-1 x6^-1 x7 x5^-1 x6^-1 x5 x6^-1 x7^2 x6^-2 x7 x5^-1 x6^-6 x5^3 x6^-3 x7 x5 x6 x7 x5^-1 x6^-1 x7 x6^-1 x5 x6^-1 x7^2 x5^-1 x6^-1 x5^-1 x6^-1 x7^2 x5^2 x6^-1 x5 x6^-1 x5 x6^-1 x5 x6^-1 x7 x5 x6^-1 x7^2 x5^-1 x6^2 x7 x5^-1 x6^-1 x5 x6^-1 x7^2 x5^-2 x6 x7 x5^-1 x6 x7 x6 x7 x5^-1 x6^-2 x5 x6^-2 x7^2 x5^-1 x6^-1 x7 x6^-2 x5^2 x6^-2"_w, "x6^-1 x7^-1 x6 x5 x6^2 x7^-1 x5 x6^-1 x5 x6^-1 x7^-1 x5 x6^-2 x5^2 x6^-1 x7^-1 x6 x7^-1 x5 x6 x7^-2 x5 x6 x7^-1 x6^2 x7^-1 x5 x6^-7 x5 x6^-1 x5 x6^-1 x7^-1 x6 x7^-1 x5^3 x6^-3 x5^2 x6^-1 x7^-1 x6^5 x5^-1 x6 x7^-1 x5 x6^-2 x5 x6^-1 x7^-1 x5 x6^-1 x7^-1 x6^-1 x7^-1 x6^2 x5^-1 x6 x7^-1 x5^3 x6^-1 x7^-1 x6^-3 x7^-1 x6 x5^-1 x6 x7^-1 x5 x6^-1 x7^-1 x6^-1 x7^-1 x6^-1 x7^-1 x6^-1 x5^2 x6^-1 x7^-1 x6^-1 x7^-1 x5^-1 x6 x7^-1 x5 x6^-5 x7^-2 x5^-1 x6 x7^-1 x5 x6^-1 x7^-1 x6 x7^-1 x5^2 x6^2 x5^-1 x6 x5^-3 x6 x7^-1 x6 x7^-1 x5 x6^-1 x7^-2 x6^-2 x7^-1 x6 x5^-1 x6 x7^-1 x6 x7^-1 x6 x5^-1 x6^2 x5^-3 x6 x7^-1 x6 x5^2 x6^4 x7^-2 x6 x7^-1 x5 x6 x5^-1 x6^3 x5^2 x6 x7^-1 x6^3 x5 x6^2 x7^-1 x5 x6 x5^-1 x6 x7^-1 x5^4 x6^-1 x7^-2 x5^-1 x6 x7^-2 x5^2 x6^-1 x7^-1 x6^-1 x7^-1 x6^2 x5^-1 x6 x7^-1 x5 x6^2 x5^-4 x6^3 x5^-2 x6 x7^-4 x5^-1 x6 x7^-1 x6^2 x7^-1 x5 x6^2 x5^-1 x6 x7^-2 x5^-2 x6^2 x7^-1 x5^-2 x6 x5^-1 x6^2 x7^-1 x5 x6 x5^-1 x6 x7^-1 x5^-1 x6 x7^-1 x5 x6 x7^-2 x5^-1 x6 x5^-3 x6^3 x5^-2 x6^2 x7^-1 x6^-1 x5 x6^-1 x7^-1 x6^-1 x7^-1 x5^2 x6^-1 x7^-1 x5 x6 x7^-2 x6^2 x7^-1 x6 x7^-2 x6 x5^-1 x6^2 x7^-1 x6^-3 x7^-1 x6^-1 x7^-1 x5 x6 x7^-1 x6 x7^-3 x5^-1 x6 x7^-1 x5 x6 x5^-1 x6^2 x7^-1 x5 x6^5 x5^-1 x6 x7^-3 x5^-2 x6 x5^-2 x6 x5^-1 x6 x7^-1 x5^3 x6"_w, "x7^2 x6^-1 x7 x5 x6^-4 x7^-1 x5^2 x6^-1 x5 x6^-2 x7^2 x6^-1 x5 x6^-1 x7^-1 x5 x6^-1 x5 x6^-1 x7 x5 x6^-3 x7^-2 x5 x6 x5 x6 x7 x5 x6^-2 x7^-1 x5 x6^3 x7^-1 x5^2 x6^-1 x7 x6^-1 x5 x6^-1 x7^-2 x6^-2 x5 x6^-3 x7^-1 x5 x6 x7 x5^3 x6 x7 x6^4 x7^-1 x6 x7^-2 x5 x6^-1 x7^4 x6^-1 x7 x6^-1 x5^2 x6^-2 x7^-1 x5^4 x6^-1 x7 x6^-2 x7 x6^-2 x7^2 x6^-1 x5 x6^2 x7^-1 x5 x6^2 x7^-2 x6^3 x7^-1 x5 x6^-1 x7^-1 x5 x6^-1 x7 x6^-1 x7^-2 x5 x6 x5 x6^-1 x5 x6^-1 x7 x5 x6^-1 x7^-1 x5 x6 x5 x6^-1 x7^-1 x6^-2 x7^-1 x5 x6^3 x5 x6^-2 x7 x5 x6^-1 x7^-2 x6^-1 x7^-1 x5^2 x6^-2 x7 x5^2 x6^-1 x5 x6^-1 x7^-1 x5 x6 x5^2 x6^-1 x7 x6^-2 x5 x6^-1 x5^2 x6^-1 x7^-4 x5 x6 x5^2 x6^2 x5 x6^-5 x7 x6^-2 x5 x6^-2 x7^-1 x5 x6 x5 x6^-1 x7^-2 x6^-1 x5 x6^-2 x7^-1 x5 x6^-5 x7 x6^-1 x7^-3 x5 x6^-1 x7 x5 x6^-1 x7^-1 x5 x6 x5 x6^2 x7^-1 x5 x6 x7^-2 x5 x6 x5 x6^-3 x7^-1 x5^2 x6^-2 x7^-1 x5 x6 x7^-1 x6 x7^-1 x5 x6 x5 x6^-1 x5 x6^-1 x7^-2 x6^-1 x7^-1 x5^2 x6^-2 x5 x6^-1 x5 x6^-1 x7^-2 x5 x6^-1 x7 x5 x6^-1 x7 x6^-3 x7^-1 x5 x6 x5^2 x6^-1 x7 x6^-1 x5 x6 x7 x6 x7^-3 x6^4 x7^3 x5 x6^-2 x7^-2 x5 x6^2 x7^2 x5 x6^-1 x7^-1 x5 x6^-1 x7 x5 x6^3 x5 x6^-2 x7 x5 x6 x7 x6 x5^3 x6 x7 x6^2 x5^2 x6 x7 x6^2 x5^2 x6 x7 x6^7 x7 x5^3 x6 x7^2 x6^3 x7^2 x5 x6 x5^2 x6 x7 x5 x6"_w, "x6^-1 x7^-1 x5^-1 x6^-1 x5 x6^-1 x7^-3 x6^-1 x7^-1 x5^3 x6 x7^-2 x5^-1 x6 x7^-1 x5^3 x6 x7^-1 x5^-1 x6 x7^-3 x5 x6 x7^-1 x5 x6^-1 x7^-1 x6^-1 x7^-2 x5^3 x6 x5^-1 x6 x7^-1 x6^-1 x7^-1 x6^-1 x5^2 x6^-1 x7^-1 x6 x7^-1 x6^-2 x7^-1 x5 x6^6 x7^-1 x5^3 x6 x7^-2 x5^-1 x6 x7^-1 x5 x6^-1 x7^-2 x6^2 x7^-1 x5 x6 x7^-1 x5^-1 x6^2 x7^-1 x6^2 x7^-1 x6^-1 x7^-1 x5 x6^-1 x7^-1 x5 x6^-1 x7^-1 x6^-3 x7^-1 x6 x7^-4 x5 x6 x7^-1 x6^3 x7^-3 x5^-1 x6^2 x7^-2 x6^-1 x5^3 x6^-2 x7^-1 x6 x7^-1 x6^2 x7^-1 x5 x6^2 x7^-1 x6^-4 x7^-1 x6 x7^-2 x6 x7^-1 x6 x7^-2 x5 x6^-1 x5 x6^-1 x5 x6^-3 x5 x6^-1 x7^-3 x5^-1 x6^2 x7^-1 x5 x6 x5^-1 x6 x5^-3 x6 x7^-2 x5 x6^2 x7^-1 x5^-1 x6 x7^-1 x5 x6^-1 x7^-1 x6^2 x7^-2 x5^2 x6^-1 x5 x6^-1 x7^-1 x5^3 x6 x7^-1 x6^-3 x5 x6^-1 x7^-1 x5^-1 x6^2 x7^-1 x5 x6 x7^-1 x6^2 x7^-1 x6 x7^-1 x6^2 x7^-2 x6 x7^-1 x5 x6 x7^-1 x6^3 x5^-1 x6^2 x7^-1 x6 x7^-1 x5 x6^-1 x5 x6^-1 x7^-1 x6^-1 x7^-6 x6 x7^-2 x6 x7^-2 x6 x5^-1 x6 x7^-1 x5^3 x6^-1 x7^-1 x6 x7^-3 x6^2 x7^-1 x5 x6^-1 x7^-1 x6^2 x7^-1 x6 x7^-1 x5 x6^-1 x7^-2 x6^4 x7^-1 x5^-1 x6 x5^-1"_w, "x6^-1 x7^-1 x5^-1 x6^-1 x7^-1 x5^-2 x6^-2 x7^-1 x6^-2 x7^-1 x5^-3 x6^-1 x5^-2 x6^-1 x7^-4 x5^-2 x6^-1 x5 x6^-1 x5 x6^-1 x5^2 x6^-1 x7^-1 x5 x6^-1 x7^-1 x5^-1 x6^-1 x5^-3 x6^-1 x5^-1 x6^-1 x5^-1 x6^-1 x5^-1 x6^-3 x5^4 x6^-1 x5 x6^-2 x5^-1 x6^-1 x7^-1 x6 x7^-1 x5^-1 x6^-1 x7^-1 x5 x6 x7^-1 x6^3 x7^-1 x5 x6^4 x5^5 x6 x7^-1 x5 x6^2 x7^-1 x6^4 x7^-1 x5 x6 x7^-4 x5^-1 x6^-1 x7^-1 x5 x6 x7^-1 x5^-1 x6 x7^-1 x5 x6 x7^-1 x5^-1 x6^-1 x7^-1 x5^2 x6 x5^-1 x6 x7^-1 x6^-1 x7^-1 x5 x6^4 x7^-4 x5^-1 x6 x7^-1 x6^2 x7^-1 x5 x6 x5^2 x6 x7^-1 x6^2 x7^-1 x6 x5^-1 x6^2 x7^-1 x6^2 x5^-3 x6 x7^-2 x6 x7^-2 x5^-2 x6^-1 x5 x6^-2 x5^2 x6^-1 x7^-2 x5^-1 x6 x7^-1 x5 x6^2 x5^-1 x6^2 x7^-2 x5 x6 x5 x6^2 x5^-1 x6 x7^-1 x6 x7^-2 x5 x6 x5^-1 x6^2 x7^-1 x5 x6 x7^-1 x6 x7^-1 x6^2 x7^-1 x5^-1 x6^-1 x7^-1 x5 x6 x5^-1 x6 x7^-1 x5^3 x6 x7^-1 x6^2 x7^-2 x5^-1 x6 x7^-1 x5 x6^2 x7^-1 x5^2 x6^2 x7^-1 x5^3 x6 x5^-1 x6 x7^-1 x5^-1 x6^2 x7^-1 x6^-1 x5^2 x6^-1 x7^-1 x5^3 x6^4 x7^-1 x6^-1 x5 x6^-1 x7^-1 x5 x6^4 x5 x6^2 x5 x6 x7^-1 x6^3 x7^-1 x6^-1 x7^-1 x5^2 x6^3 x5^-1 x6 x7^-1 x6 x7^-2 x5 x6 x5 x6 x5^-1 x6 x7^-4 x6^2 x7^-1 x5 x6 x7^-1 x6^2 x7^-1 x5^-3 x6^-2 x7^-1 x5 x6^-3 x7^-1 x5^2 x6^-1 x7^-1 x6^-1 x7^-1 x6^-1 x7^-1 x6^-1 x7^-1 x6^-5 x7^-1 x6 x7^-4 x5^-6 x6 x7^-2 x5^-1"_w, "x7^4 x5 x6^-3 x7^2 x5^-1 x6^-2 x5^4 x6^-1 x7 x5^2 x6^-1 x7 x5^3 x6^4 x7 x5^-2 x6^3 x7^4 x5^-1 x6 x7 x5^-1 x6^-1 x7 x6 x5^-1 x6 x5^-1 x6^2 x7 x5^-1 x6^-1 x7^2 x5^-1 x6^-1 x5^3 x6^-3 x7^7 x5^-1 x6 x5^-2 x6 x5^-1 x6 x5 x6 x7 x5 x6^-2 x7 x6^-1 x7 x6^-5 x5^2 x6^-1 x7 x5^-2 x6^-1 x5^3 x6^-1 x5 x6^-4 x7^2 x6^-2 x7 x5^-2 x6^-1 x7 x5^4 x6^-1 x7 x6^-3 x5^2 x6^-1 x7 x5^-2 x6^-1 x7 x6^-1 x7 x5 x6^-1 x7 x5^-1 x6^-4 x5 x6^-2 x7 x6^-3 x7^3 x5 x6^-1 x7 x5^-1 x6^-1 x5 x6^-1 x7 x5^-1 x6^-1 x7^3 x6^-1 x7^3 x6^-1 x5^-1 x6^-3 x5 x6^-1 x7^3 x6^-1 x7 x5^-2 x6^-1 x7 x5 x6^-1 x5 x6^-2 x7 x5^-1 x6^-1 x7 x5^2 x6^-5 x7 x5^-1 x6^-1 x7 x5^-1 x6^-2 x7^2 x6^-3 x7 x6 x5^-2 x6 x7 x6 x7^2 x6^-1 x7 x5^-3 x6^-1 x7 x5^-1 x6^-1 x7 x6^-1 x7 x5^-1 x6^-1 x5^-1 x6^-1 x7^4 x5^-1 x6^-1 x7^2 x6^-1 x5 x6^-1 x7 x6^-1 x7 x5^-1 x6 x7^4 x5^-3 x6 x7 x5 x6^-2 x7 x6^-1 x7 x5^-1 x6 x7^2 x6^-1 x5 x6^-1 x7^2 x6^-1 x7 x5^-1 x6^-1 x7 x6^-6 x7 x5^-1 x6^-5 x7 x6^-2 x7 x5 x6^2 x5^-1 x6 x5^-1 x6 x5^-2 x6^3 x7 x6^-1 x7^2 x6^-2 x7^2 x6 x5^-1 x6 x5^-2"_w, "x5^2 x6^-2 x7 x6 x5^-1 x6 x7 x5^-1 x6^-2 x7 x5^-1 x6^-1 x7^2 x5^-1 x6^-1 x7 x6^-1 x7 x5 x6^-5 x7^2 x5^-1 x6^-1 x7 x6^-1 x5^-2 x6^-2 x7^2 x5^-1 x6^-1 x7 x5 x6 x5^-4 x6^2 x7^3 x5^-3 x6^-3 x7 x5 x6^-1 x7 x6 x5^-1 x6^3 x7 x5^-1 x6^-1 x7 x6^-1 x5^-5 x6^-2 x7 x5^-1 x6^-2 x5^-1 x6^-2 x7^2 x5^-1 x6^-2 x7 x5 x6^-2 x5^2 x6^-1 x7^2 x5 x6^-1 x7 x6^-2 x7 x5^-1 x6^-1 x7 x6^-1 x7 x5 x6^2 x5^-2 x6^3 x5^-1 x6 x7 x6^-1 x7 x6^-2 x7 x5^-1 x6 x7 x5^-1 x6^-1 x5 x6^-2 x7 x5^3 x6^-1 x7^2 x5^-1 x6^-2 x7^2 x5^2 x6^-3 x5^2 x6^-2 x7 x5^-1 x6^-1 x7 x5^-1 x6^-1 x5 x6^-3 x5^2 x6^-2 x5 x6^-1 x7^2 x5^-2 x6^-1 x7 x5 x6^-2 x7 x6^-1 x7 x5^-1 x6 x7^2 x5^6 x6^-1 x7 x5^-1 x6^-1 x5 x6^-1 x7 x5^-1 x6 x7^3 x5 x6^-2 x7^3 x5^-2 x6 x7 x5^-2 x6 x7 x5 x6^-1 x7^2 x5^-2 x6 x7 x6 x7 x5^-1 x6 x7 x6 x7^2 x5 x6^-1 x7^2 x6^-1 x7 x6^5 x5 x6 x7^2 x6^-1 x7^2 x5 x6 x5^2 x6 x5 x6^5 x7 x5 x6^3 x5^4 x6 x7 x6 x5^2 x6^3 x7^2 x5^2 x6 x5^2 x6 x7 x5 x6^2 x5^3 x6 x7 x6 x5^3 x6 x7 x6^5 x7 x6^2 x7 x5 x6 x5 x6^2 x5^-1 x6 x7 x5 x6"_w, "x7^2 x5^3 x6^-1 x7^2 x5^-1 x6 x5^-2 x6 x7 x5^-1 x6^-1 x7 x6^-1 x7 x5^3 x6^-1 x5 x6^-1 x5^-3 x6^-1 x7 x6^3 x5^-1 x6 x7 x5^-1 x6^-1 x5 x6^-2 x7 x6 x5^-3 x6^3 x7^3 x6 x7 x6 x7 x5^-1 x6^3 x5^-1 x6^4 x7 x5^-1 x6 x7 x5^-1 x6^3 x5^-1 x6 x7 x6^-2 x7^3 x6 x7 x5^-1 x6^-3 x7 x6^-1 x5 x6^-1 x7^2 x5^-1 x6^-1 x7 x6^-1 x7 x5 x6^-1 x7 x5^-1 x6^-2 x7 x5 x6^-1 x7^2 x6 x5 x6^3 x5^-3 x6^6 x7 x5^-1 x6^3 x5^-1 x6 x7^2 x6^-2 x7^2 x5^2 x6^-1 x5 x6^-1 x7^2 x5^-3 x6^-1 x5 x6^-2 x5 x6^-2 x7 x5^-1 x6^-1 x7 x5^-1 x6^-1 x7 x5 x6^-2 x5^2 x6^-1 x7 x6^2 x5^-1 x6^2 x7^2 x6^-2 x7^2 x5^3 x6^-2 x7 x5^-2 x6^-2 x7 x6^-1 x7 x5 x6^-1 x7 x5^-3 x6 x7 x6^2 x5^-1 x6 x5^-1 x6 x7 x5^-1 x6 x7 x5^-1 x6^-2 x7 x6^-1 x5 x6^-1 x5 x6^-1 x7 x5^-2 x6^-2 x7 x5^-2 x6^-1 x7 x5 x6^-1 x5 x6^-1 x5^3 x6^-3 x7^3 x5^2 x6^-1 x7 x6^-2 x7 x5^-1 x6^-1 x7^2 x6^-1 x7 x6^-1 x5 x6^-1 x7 x5^-1 x6 x7^2 x6^-2 x7^2 x5 x6^-1 x7 x5^-1 x6^-1 x7 x6^-1 x5^-2 x6^-1 x7 x5^-1 x6^-1 x5 x6^-4 x5 x6^-2 x7 x5 x6 x5^-2 x6 x7^2 x5^-2 x6^3 x7 x6^-1 x7 x5^-1 x6^-1 x5^-2 x6^-1 x7 x5^-1 x6^2 x7 x5^-2 x6^-1 x5^2 x6^-1 x7^2 x5^-1 x6^-1 x7 x6^-1 x7^2 x6^-2 x5 x6^-1 x7 x5^-1 x6 x7 x5 x6 x5^-1"_w, "x7 x6^-1 x5^-1 x6 x5^-1 x6 x7^2 x5^-1 x6^2 x7 x5^-1 x6 x7^-1 x6^2 x7^-2 x6 x7 x5^-1 x6^2 x5^-1 x6^3 x7 x6 x7 x6 x7 x6 x7 x5^-1 x6 x5^-1 x6 x7^-1 x6 x7 x5^-1 x6^4 x7^-1 x5^-3 x6^2 x5^-1 x6 x7^-3 x5^-1 x6 x7 x5^-1 x6 x7^-1 x5^-2 x6 x5^-1 x6^2 x7 x5^-1 x6^2 x7^-1 x6 x7^-1 x5^-1 x6 x7^-1 x6 x7^-2 x6^3 x7 x5^-2 x6^4 x7^-1 x6 x7^-1 x5^-1 x6^-2 x5^-2 x6^-2 x7^2 x5^-1 x6^2 x7^-1 x5^-1 x6^2 x7^-3 x6^2 x5^-1 x6^-3 x7 x6^-1 x7 x5^-2 x6 x7^-1 x6 x7^-1 x5^-2 x6^3 x7^-1 x5^-1 x6 x5^-1 x6 x7 x5^-1 x6^-2 x5^-1 x6 x7^-1 x6 x7 x5^-1 x6^2 x5^-1 x6 x7 x5^-3 x6 x5^-2 x6^-1 x7^-1 x6^-4 x7 x6^-1 x7^4 x5^-1 x6 x5^-1 x6 x5^-2 x6^-2 x7 x6^-1 x7^6 x5^-1 x6 x5^-1 x6^3 x7 x5^-1 x6 x7^-1 x6 x7 x6 x7 x5^-1 x6^2 x7^-1 x6^3 x7 x6 x7 x6 x7 x5^-1 x6^-2 x7 x6^-3 x7 x5^-1 x6 x7 x6 x7 x5^-1 x6^2 x7^-1 x6^3 x7 x5^-1 x6^3 x5^-3 x6 x7 x5^-1 x6 x7^-1 x6 x7^-1 x6 x7 x5^-1 x6^-1 x5^-1 x6 x5^-1 x6 x5^-1 x6 x5^-1 x6 x7 x6 x7 x5^-1 x6 x5^-2 x6^-1 x7^-1 x6^-3 x7 x6^-2 x7 x6^-1 x7 x6^-1 x7^4 x6^-1 x5^-1 x6 x7^-3 x6^3 x7 x5^-1 x6^2 x7^-1 x6 x7 x5^-1 x6^-1 x5^-1 x6^-1 x7 x5^-1 x6^-1 x5^-1 x6^-4 x5^-1 x6 x5^-1 x6^-1 x7 x6^-2"_w, };
  T.WR = { "x3 x1^4 x2^-1 x3^-1 x1^2 x2^-1 x3^-1 x1 x2^-1 x3^-2 x1^2 x2^-2 x3^-1 x1 x2^-1 x3^2 x2^-1 x3 x1 x2^-1 x1 x2^3 x1 x2^-1 x3^-1 x1 x2^-4 x1^2 x2^-2 x3^-1 x1 x2^-1 x3 x1 x2^-3 x3 x1^3 x2^-1 x3^-1 x1 x2^-1 x3^-1 x2^-2 x3^-3 x2^-1 x3 x2^-1 x3^-1 x1^3 x2^4 x3^-2 x1 x2^-1 x3^9 x2^-1 x3 x1 x2^2 x3^-2 x1 x2^-1 x1 x2 x3^-2 x2^2 x3^-1 x1 x2^-1 x3 x2^-1 x3^3 x2^-1 x3^-3 x1 x2^-1 x3 x2^-1 x1 x2^2 x3^-1 x2 x1 x2^-1 x3^-2 x2^-1 x3^2 x2^-1 x3^-1 x1 x2^3 x3^-1 x2^2 x1 x2 x3^-2 x1 x2^-1 x3^2 x1 x2 x3^-1 x2^2 x3^-1 x2 x1 x2^-1 x1 x2^-1 x1 x2 x3^-5 x2 x3^-1 x1 x2^2 x3^-2 x1 x2^-1 x3^-1 x1^3 x2 x3^-1 x2^3 x1 x2^-3 x3^-2 x1^3 x2^-1 x3 x1^2 x2^-1 x3^-1 x1 x2^-1 x1 x2^-3 x3^2 x2^-2 x1 x2 x3^-1 x2^2 x3^-1 x1^2 x2^-1 x1 x2^-1 x3^-1 x1 x2^-2 x3^3 x1^4 x2^-2 x3^-2 x1 x2^-1 x1 x2^2 x1 x2^-4 x3^-1 x1^4 x2^-1 x1 x2^-1 x1^5 x2^-1 x3^-1 x2^-3 x1 x2^3 x1 x2^-1 x3^-1 x1 x2^-2 x3^-1 x1 x2^-3 x3^-2 x1 x2^-1 x3 x1^-1"_w, "x2^-1 x3^-1 x1^-1 x2^-1 x3^-2 x1^-2 x2^-1 x3^-1 x1^-1 x2^-1 x3^-2 x2^-1 x1^-2 x2^-2 x3^-2 x2^-1 x3^-1 x1^-1 x2^-1 x3 x1^-1 x2^-5 x1^-1 x2 x3^-1 x1^-3 x2^-1 x3^-1 x2^-3 x3^3 x2^-1 x3^-2 x1^-1 x2 x3^-2 x1^-1 x2^-4 x3^2 x1^-3 x2 x3^-1 x1^-1 x2^6 x3 x1^-2 x2^2 x3^-1 x1^-3 x2^-1 x3^-1 x2^-6 x3^2 x2^-1 x3 x2^-1 x3 x2^-1 x3^-1 x1^-2 x2 x1^-1 x2^-1 x1^-1 x2^-1 x3 x2^-4 x1^-1 x2^6 x3 x1^-1 x2^2 x3 x2 x1^-1 x2 x1^-2 x2^-1 x3 x1^-1 x2 x3^-1 x1^-4 x2 x3 x1^-1 x2 x3^-1 x2 x1^-1 x2 x3^-1 x2 x1^-1 x2^-3 x3^2 x2^-3 x3 x2^-1 x3^2 x2^-1 x1^-1 x2 x3 x1^-1 x2 x3 x2 x3 x1^-2 x2 x3^-4 x2^2 x3^-1 x2 x3^3 x1^-1 x2 x3 x2 x3 x1^-1 x2 x3 x2 x3 x2^2 x3^2 x2^9 x3 x1^-2 x2 x3 x2^2 x3^2 x2^3 x3 x2 x1^-1 x2 x3 x1^-1 x2 x1^-4 x2^2 x3^-1 x1^-1 x2 x1^-2 x2 x1^-2 x2 x1^-1 x2^-1 x3 x1^-1 x2^-1 x3 x1^-1 x2^-1 x3 x2^-2 x1^-1 x2^5 x1^-1 x2^2 x1^-4 x2^-1 x1^-1 x2^-1 x3^2 x1^-1 x2^-1 x3^2 x1^-1 x2^-1 x3 x2^-1 x1^-3 x2^-1 x3 x1^-1 x2 x3^-1 x2 x1^-1 x2^-1 x1^-1 x2^-1 x3 x1^-1 x2 x3 x1^-1 x2 x3^-1 x2 x3^-1 x2^2 x1^-1 x2^-1 x3 x2^-1 x3^5 x1^-2 x2 x3^-2 x2 x3 x1^-1"_w, "x2^-1 x3^-1 x2^-1 x1^-1 x2^-1 x3^-2 x1^-1 x2^-1 x3^-1 x1^-1 x2^-1 x3^-1 x2^-1 x3^-1 x1^-1 x2^-1 x3^-2 x2^-1 x3^-2 x1^-2 x2^-2 x3^-2 x2^-1 x3^-1 x1^-2 x2^-2 x1^-4 x2^-1 x1^-2 x2^-1 x3^-2 x1^-2 x2^-2 x1 x2^-1 x1^3 x2^-1 x3^-1 x1^-1 x2 x1^-1 x2 x3^-1 x1 x2 x1 x2^2 x3^-1 x1 x2^-2 x3^-2 x1^-1 x2^-1 x1^4 x2^-1 x1 x2^-2 x1 x2^-1 x3^-1 x2 x3^-2 x2 x3^-2 x1 x2 x1^-1 x2 x3^-2 x1 x2 x3^-1 x1^-2 x2^3 x3^-1 x1^2 x2 x1^-1 x2 x3^-1 x2^-1 x3^-1 x1^4 x2 x3^-2 x1^-1 x2^-1 x1 x2^-1 x3^-3 x1 x2 x3^-1 x1^-2 x2^3 x3^-1 x1 x2^4 x3^-1 x2 x3^-1 x2 x3^-2 x1^-1 x2 x3^-1 x1 x2^-1 x3^-1 x1^3 x2^2 x1^2 x2 x1^-1 x2 x3^-1 x2 x3^-1 x1 x2^2 x3^-1 x1 x2 x1^-1 x2^2 x3^-2 x2^-3 x3^-1 x1 x2 x3^-1 x1^-1 x2 x3^-1 x2 x3^-1 x2^-1 x3^-1 x2^-1 x1 x2^-1 x3^-1 x1 x2 x1^-1 x2 x3^-2 x1 x2 x3^-1 x2 x3^-3 x1^-1 x2 x3^-1 x1 x2 x3^-2 x2 x3^-1 x2 x1^-1 x2 x1^-1 x2^2 x3^-2 x1 x2^-1 x3^-1 x1 x2 x3^-1 x1^-1 x2 x3^-2 x1 x2^2 x1^4 x2 x1 x2 x3^-1 x1 x2 x3^-2 x1^-1 x2^2 x3^-1 x1 x2 x3^-1 x1^-1 x2^-1 x3^-1 x1 x2 x1^-1 x2 x3^-1 x2 x3^-2 x2 x3^-1 x1^2 x2^-2 x1^2 x2^-1 x1 x2^-1 x3^-1 x2^3 x3^-1 x1 x2 x3^-2 x2 x1^2 x2^2 x1^-1 x2^2 x3^-1 x2 x1^-2 x2 x1^-1 x2^2 x3^-1 x2^-1 x3^-1 x1 x2^-1 x1 x2^-2 x3^-1 x1 x2 x1^-1 x2 x3^-1 x1 x2^2 x3^-1 x2^-2 x3^-2 x1 x2^-2 x3^-1 x2^-1 x1^-1"_w, "x3^2 x2^-1 x3 x2^-1 x3 x1^-1 x2^-1 x3^2 x1 x2^-1 x3 x1^-2 x2^-2 x3 x1^-1 x2^-1 x3 x1 x2^-1 x3 x1^-1 x2 x3 x2^-2 x1 x2^-1 x3 x1^-2 x2 x3 x1 x2^-1 x3 x1^-2 x2 x3 x2 x3 x1^-1 x2^3 x1^-2 x2^4 x1^-2 x2 x1^-1 x2 x3 x2^2 x3 x2^-3 x3 x1^-2 x2^-1 x1 x2^-1 x3^2 x1^-2 x2^-4 x3 x1^-1 x2^2 x3 x2^-2 x1 x2^-3 x3 x1^-1 x2 x3 x2^-1 x3 x1^-2 x2^2 x1^-1 x2 x1^-1 x2 x3 x2^-2 x1 x2^-1 x1 x2^-2 x3 x1^-2 x2^2 x1^-1 x2 x3 x1^2 x2^-1 x3 x2^-1 x1 x2^-1 x1 x2^-1 x1 x2^-2 x1^3 x2^-2 x1^2 x2^-1 x3^2 x1^-2 x2 x3 x2 x3 x1^-3 x2 x1^-1 x2 x3 x1^-1 x2^2 x3^5 x1 x2 x1 x2^4 x1^-1 x2 x1^-3 x2^3 x1^2 x2 x3 x2^-1 x3 x1^-1 x2^-1 x3 x2^-1 x1 x2^-2 x3 x2^-2 x3 x2 x3 x2 x3 x1^-2 x2^-1 x3 x2^-1 x3 x2^-2 x3 x1 x2^-3 x3 x2^-1 x3 x2^-1 x3^2 x1^-1 x2^-1 x3 x1^-1 x2^-2 x3 x1 x2^-2 x3 x2^4 x3 x1^-2 x2^-1 x1^2 x2^-2 x1^4 x2^-1 x3 x1 x2 x3 x2 x3 x1^-1 x2 x1^-1 x2 x3^2 x2^2 x3 x2 x1^-1 x2 x3 x1^-1 x2^-1 x3 x1 x2^-1 x3^4 x2^-1 x3 x1^-2"_w, "x7 x5 x4 x5 x6 x7 x3 x4 x5 x6 x2 x3 x4 x5 x1 x2 x3 x4 x3 x2 x1 x2 x3 x4 x5 x6 x7 x2 x3 x4 x5 x6 x2 x3 x4 x5 x6 x7 x1 x2 x3 x4 x5 x6 x3^-1 x2^-1 x3^-1 x1^-1 x2^-1 x3^-2 x1^-1 x2^-1 x3^-2 x2^-1 x1^-1 x2^-1 x3^-2 x1^-2 x2^-2 x3^-1 x2^-2 x3^-1 x2^-2 x3^-1 x2^-1 x3^-1 x1^-2 x2^-2 x3^-5 x1^-1 x2^5 x1^2 x2^3 x3^-5 x1 x2 x3^-1 x1 x2^3 x3^-2 x1^-1 x2 x3^-1 x2 x3^-4 x1^5 x2 x1^-4 x2 x1^-1 x2 x1^-2 x2^2 x3^-1 x2^-2 x1^3 x2^-1 x3^-1 x2^2 x3^-1 x1 x2 x3^-1 x1^4 x2^-1 x1 x2^-3 x3^-1 x2^3 x3^-1 x2 x3^-1 x2 x1 x2 x3^-1 x1 x2 x1^-1 x2 x3^-1 x1 x2 x3^-1 x1^-4 x2^-3 x3^-1 x1^2 x2^2 x3^-1 x2^-2 x1 x2^-1 x3^-1 x1 x2 x3^-1 x1 x2 x3^-1 x1 x2^2 x3^-4 x2 x3^-2 x1^4 x2^2 x1^-1 x2 x3^-1 x1 x2 x1^-1 x2^3 x3^-2 x1^-1 x2 x3^-1 x1 x2 x3^-1 x1^-2 x2 x3^-1 x2^2 x1^3 x2 x3^-1 x1^-1 x2 x1^-2 x2 x3^-1 x1^3 x2^4 x1^-1 x2^2 x3^-2 x2 x1^-1 x2^3 x3^-1 x2^3 x1^3 x2 x3^-1 x2^3 x3^-2 x1^2 x2^2 x1^-1 x2 x3^-1 x1 x2^2 x3^-2 x1^-4 x2 x1^-1 x2 x3^-1 x1 x2 x3^-1 x2^-1 x1 x2^-2 x3^-3 x2^-2 x1 x2^-2 x3^-3 x1^2 x2 x3^-1 x1^-2 x2^-1 x1 x2^-2 x3^-1 x1 x2^-3 x1 x2^-3 x3^-1 x1 x2^-1 x3^-1"_w, "x3^-1 x1 x2^2 x3 x1^2 x2^4 x3^-1 x1^3 x2^-1 x3 x2^-1 x3^-1 x1 x2^-1 x3 x1^2 x2^-1 x3^-1 x1 x2^-1 x3 x1 x2^-2 x3^-1 x1 x2^-1 x3 x2^-1 x3^-1 x1 x2^-2 x3 x1 x2 x3^-1 x1 x2^-2 x3 x2^-1 x1 x2 x3^-1 x1 x2^-1 x1 x2^-1 x3 x1 x2^-1 x3^-1 x1^2 x2^-5 x3 x1^4 x2^-1 x3 x2^-1 x3^-1 x1^2 x2^-2 x1 x2^-2 x3^-3 x2^-1 x1 x2^-2 x3 x1^2 x2^-1 x3 x2^-2 x3^3 x2^-1 x3 x2^-1 x3^-1 x1 x2^-3 x3 x2^-2 x3^-1 x1 x2 x3^-1 x1 x2^-1 x3 x2^-1 x3^-1 x1 x2^-2 x3^-1 x2^-3 x3^-1 x1 x2^-1 x3^-1 x2^-1 x3 x2^-1 x3 x1^2 x2^-2 x3 x2^-1 x3 x2^-1 x3^-1 x1 x2^-1 x3^-2 x2^-1 x3 x2^-2 x1 x2^-2 x3^-1 x1 x2^-1 x3 x1 x2^-5 x1^3 x2 x3^-2 x1 x2^-1 x1^3 x2^-1 x3 x2^-3 x3^-1 x1^2 x2^-1 x3 x2^-2 x3^3 x1^2 x2^-2 x1 x2^-2 x3 x1 x2^-1 x3^-1 x1 x2^-2 x1 x2^-2 x3 x2^-1 x3^4 x1^3 x2^-1 x3^-1 x1 x2^-1 x3 x2^-1 x1 x2^-1 x3^-2 x2^-2 x3^-2 x1 x2^-1 x3^-1 x1 x2^-2 x3^-1 x1^2 x2^-2 x3^-2 x1 x2^-1 x3^3 x2^-2 x1 x2^-3 x3^-1 x1 x2^-1 x3^-1 x2^-1 x1^2 x2^-3 x3^-1 x1 x2^3 x3 x1^2 x2^-1 x1 x2^-1 x3^-1 x1 x2 x3 x1 x2 x3^-3 x2^5 x3^-2 x2 x1 x2^-1 x3^5 x1 x2^4 x1 x2^-1 x1 x2^-1 x3^-1 x1 x2^2 x3^-2 x2^2 x3^-1 x1 x2^3 x3 x1^3 x2 x1 x2 x3 x2"_w, "x3 x1 x2^-1 x1^2 x2^-2 x3^-1 x2 x3^-1 x1^-1 x2 x3^-2 x1^-1 x2^2 x3^-1 x1 x2^2 x1^-2 x2 x1^-1 x2^2 x3^-2 x2 x3^-1 x1^3 x2 x1 x2 x3^-1 x2^-3 x3^-2 x1 x2^3 x3^-1 x1^-1 x2 x3^-1 x1 x2 x1 x2 x1 x2 x3^-1 x2^-2 x3^-1 x1 x2^-1 x3^-1 x1 x2^2 x3^-1 x2 x3^-2 x1^-2 x2^-1 x1 x2^-2 x3^-1 x2 x3^-1 x2 x1 x2 x3^-1 x1 x2 x1^-1 x2^2 x3^-1 x1 x2 x3^-1 x1^-2 x2 x1^-2 x2^3 x3^-1 x1 x2^5 x3^-2 x1^2 x2^3 x3^-2 x1^2 x2^2 x3^-1 x2^-4 x3^-1 x1^2 x2^-1 x3^-3 x2^-2 x3^-1 x2^-2 x3^-1 x1 x2 x3^-2 x1^-2 x2 x1^-1 x2^5 x3^-1 x1 x2^-5 x1^2 x2^-1 x3^-2 x1 x2^-2 x3^-1 x2 x3^-1 x1 x2 x3^-1 x1^-1 x2 x1^2 x2^2 x3^-1 x2^-1 x1^2 x2^-1 x1 x2^-1 x3^-1 x2^2 x1^-1 x2 x1^5 x2^2 x3^-1 x2 x3^-1 x2 x3^-2 x1 x2^2 x1^3 x2 x3^-1 x1 x2 x3^-3 x1^-1 x2^2 x3^-1 x2 x3^-1 x1^4 x2 x3^-1 x2 x3^-1 x1 x2^4 x3^-1 x1 x2^-1 x1^4 x2^-1 x3^-2 x1^2 x2^3 x3^-1 x1^-1 x2 x1^-2 x2^2 x3^-1 x2^3 x3^-1 x2^-1 x1 x2^-2 x1 x2^-1 x3^-1 x1 x2 x1^-1 x2 x3^-2"_w, "x2^-1 x1^4 x2^-4 x3 x1^-1 x2^-1 x1 x2^-1 x3 x1^3 x2 x3 x1^-1 x2^-2 x1 x2^-1 x3 x1 x2^-1 x1 x2^-3 x3^3 x1^-2 x2^-3 x3 x2 x3 x1^-1 x2 x3 x1^-2 x2^-1 x1 x2^-1 x3 x2 x3 x1^-1 x2^-1 x1 x2^-1 x3 x1^-1 x2^-2 x1^-1 x2^-3 x3 x1^-1 x2^-1 x3^2 x2 x3 x1^-1 x2^-2 x1 x2^-1 x3^2 x1 x2^-1 x3 x1^-1 x2^-1 x3 x2^4 x1^-2 x2^3 x3 x2^-1 x3^2 x2^4 x3 x1^-2 x2 x1^-2 x2 x3 x1^-1 x2^-1 x1^2 x2^-1 x1 x2^-1 x3^2 x1^-1 x2^-2 x1^-2 x2^-3 x3 x1^-2 x2 x3 x2 x3^2 x2^-1 x3^3 x1^-1 x2 x3 x2 x3 x2^-1 x3 x1 x2^-1 x3 x1^-1 x2 x3 x2^-1 x3 x2^-3 x3 x1^-1 x2^-1 x3 x2^-2 x3 x1^-1 x2 x1^-1 x2 x3 x1 x2^-1 x3^3 x1^-2 x2^-1 x3 x2^-1 x3 x1^-1 x2^2 x3^3 x1^-1 x2^-1 x3 x2^-2 x1^-2 x2^-1 x3 x1^-1 x2^-4 x3 x1^-1 x2^-4 x1^3 x2^-2 x3 x1^-3 x2^-2 x3 x1^-1 x2 x1^-2 x2 x3 x1^-1 x2 x3 x2 x3 x2^-2 x3 x1^-1 x2 x3^5 x2^-2 x3 x2^-1 x3^2 x1^-1 x2^-2 x3^2 x2^-1 x1 x2^-1 x3^3 x1^-1 x2^-2 x3 x1^-1 x2 x1^-1 x2 x3 x2^-2 x3 x2^-3 x1 x2^-5 x3 x1^-1 x2^-4 x3 x1 x2 x1^-1"_w, "x1 x2^-1 x3^-1 x2^-1 x1^-2 x2 x3^2 x1^-2 x2 x3^-1 x2 x3^-1 x2 x3 x1^-2 x2^-2 x3^-1 x1^-1 x2 x3 x1^-1 x2^-3 x1^-1 x2 x3^-1 x2 x1^-2 x2^2 x3 x1^-1 x2^-1 x3^-1 x1^-1 x2 x3^2 x1^-1 x2^-2 x3 x2^-1 x3^-1 x1^-1 x2 x3 x2 x3 x1^-1 x2 x1^-2 x2 x3 x1^-1 x2 x3^-1 x1^-2 x2^2 x3^2 x1^-1 x2^4 x3^-1 x1^-2 x2 x3^3 x1^-1 x2^-1 x1^-1 x2 x3^-1 x2^2 x3^2 x2^2 x3^3 x1^-1 x2^-1 x3^-1 x1^-1 x2^2 x3 x1^-1 x2^5 x1^-2 x2 x1^-1 x2 x3^-2 x2^3 x3 x1^-1 x2^3 x1^-2 x2 x3 x1^-1 x2^-2 x1^-1 x2 x1^-1 x2 x3 x1^-1 x2 x3^-1 x1^-2 x2^3 x3 x1^-1 x2 x3 x2 x1^-1 x2 x3 x1^-2 x2 x3^-1 x2 x3 x1^-1 x2^-1 x3 x2^-1 x3^-1 x1^-2 x2^2 x3^2 x2^6 x3 x1^-1 x2^6 x3^-2 x2 x3 x1^-1 x2^5 x3^-1 x1^-1 x2 x3 x1^-1 x2 x3^-1 x1^-2 x2^2 x3 x1^-1 x2^-1 x1^-1 x2^-1 x3 x2^-1 x3^-1 x1^-1 x2 x1^-1 x2 x1^-4 x2^2 x3 x1^-1 x2^2 x1^-1 x2^6 x3 x1^-1 x2^-2 x3^-1 x1^-1 x2 x3^3 x1^-5 x2^-3 x3^-3 x1^-1 x2 x3 x1^-2 x2^-1 x3 x2^-1 x1^-2 x2 x3^-1 x2 x3 x1^-1 x2 x1^-1 x2^-2 x1^-1 x2 x3 x1^-1 x2^2 x1^-1 x2^-3 x3 x2^-1 x3^2 x2^-1 x1^-1 x2 x3^-1 x1^-1 x2 x3^-5 x1^-1 x2^3 x3 x1^-1 x2 x3^-1 x2 x3^-1 x2 x3 x1^-1 x2 x1^-2 x2^-3 x3^-1 x1^-1 x2 x1^-1 x2 x1^-1 x2^-1 x3^2 x1^-1"_w, "x2 x3^-1 x2 x3^-1 x1^-2 x2^-1 x3^4 x2^-2 x3 x2^-1 x1^-1 x2^2 x3 x1^-1 x2 x3^3 x1^-4 x2^2 x3^-1 x1^-1 x2^-2 x3 x1^-1 x2 x3^-1 x2 x3 x1^-2 x2^2 x3^2 x1^-1 x2 x1^-1 x2^-1 x3 x2^-1 x3 x2^-1 x3 x2^-1 x3^4 x1^-1 x2 x3^-2 x1^-1 x2 x3^-1 x2^3 x3^-1 x1^-2 x2 x3 x1^-1 x2^4 x3^-3 x2^2 x3 x1^-1 x2^-1 x3 x2^-1 x1^-1 x2^2 x3^-1 x1^-1 x2 x3 x1^-1 x2^2 x3^-1 x1^-1 x2^2 x1^-1 x2 x3 x1^-1 x2^-1 x3 x1^-1 x2^3 x3^-1 x2^2 x3 x1^-1 x2^-2 x1^-1 x2 x3^-1 x1^-1 x2 x3 x1^-1 x2^-1 x1^-2 x2 x1^-1 x2 x3^-1 x2 x3 x1^-1 x2^2 x3^-1 x1^-2 x2 x1^-2 x2^2 x3 x1^-1 x2^-1 x1^-4 x2 x1^-1 x2^3 x3^-1 x2 x3 x1^-1 x2^-1 x3 x1^-4 x2 x3^-1 x2^2 x1^-1 x2^3 x1^-1 x2 x3 x1^-1 x2^-1 x3^2 x1^-1 x2 x3^-1 x2 x3^3 x1^-1 x2^-1 x1^-1 x2^-1 x3 x2^-1 x1^-1 x2 x3^-1 x2 x3 x1^-1 x2^2 x1^-1 x2 x3 x1^-1 x2^-1 x3^-1 x1^-1 x2 x1^-1 x2^-1 x3^4 x2^-4 x3 x2^-3 x1^-3 x2 x3^-1 x2 x3^2 x1^-1 x2^2 x3 x1^-1 x2 x3^-1 x1^-1 x2^-1 x3^3 x2^-3 x3^3 x1^-1 x2^-1 x1^-3 x2^-1 x3 x1^-1 x2^-4 x3 x1^-2 x2 x3 x2 x3^2 x1^-1 x2^-3 x3 x1^-1 x2^-1 x3 x2^-1 x3^3 x1^-1 x2^2 x1^-1 x2^-1 x3 x1^-1 x2 x3^-1 x2 x1^-1 x2 x3 x1^-2 x2 x3^-1 x2 x3^-1 x2 x1^-1 x2^-1"_w, };
  T.deltaSQL = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, };
  T.deltaSQR = { 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, };
  T.z = "x2^-2 x1 x5^-1 x4^-1 x3^-3 x7^2 x6^-1 x5 x4^-1 x3^-1 x6^-1 x5^-1 x6 x5^-1 x4 x3^-1 x2^-1 x1 x7 x6 x5 x4^-1 x3^-1 x2 x1 x5^-1 x4^-1 x3^-1 x2^-1 x7^-1 x6^-1 x5^-1 x4^-1 x3 x2^-1 x1^3 x3^-1 x2^-1 x7^-1 x6^-1 x5^-2 x4^-1 x3^-1 x2^-1 x1 x3^-1 x2^-1 x1 x2^-3 x1 x6 x5 x6 x5^2 x7^-1 x6^-1 x5 x7 x6^-1 x5 x4 x3 x2 x1 x3^-1 x2 x1^2 x3^-1 x2^-1 x1 x5^-2 x4 x3 x2 x1 x4 x3 x2 x1^4 x5^-1 x4^-1 x3^-1 x2 x1 x7^-1 x6^-1 x5 x6^-3 x5 x4 x3 x2 x1 x4^-2 x3^-3 x2^-1 x1^2 x7^-1 x6^-2 x5^-3 x4 x7 x6 x7^3 x6^3 x5^2 x6^-1 x5 x4 x3 x2 x1 x3^-1 x2^-2 x1 x4^-1 x3^-1 x7^-2 x6^-1 x7^-1 x6^-1 x5 x4 x3^-1 x5 x4 x3^-1 x2^-1 x7^-1 x6 x5 x4 x3 x2^-1 x1 x2^-1 x1 x3^-2 x2^-1 x1 x7 x6 x5 x4 x3 x7^2 x6^-3 x7 x6^-1 x5^-1 x4 x3 x4 x3 x7^-1 x6^-1 x5^-1 x6^-1 x5^-1 x4 x3 x2^-2 x4 x5 x4 x7^2 x6^-3 x5 x4^-1 x3^-1 x2 x3^-1 x2^-2 x7 x6^3 x5^-1 x4^-1 x3 x2 x1 x4^-1 x3^-1 x2 x4^-1 x3 x2^-1 x1 x2^-1 x1^2 x7 x6^-2 x5 x4^-1 x5^2 x4^-1 x5^-1 x4^-1 x3^-3 x4^3 x3^-1 x6 x7 x6^3 x5^-1 x4 x6^-1 x5^2 x4 x3^-1 x2^-1 x1 x5 x4^-1 x3 x2^4 x1 x7^-2 x6 x5^-2 x4 x3^-1 x2^-1 x6^-1 x5^-1 x4^-1 x3^-1 x5^-1 x4^3 x3^-3 x2^-1 x7^2 x6^-1 x5 x4 x3^-2 x2^-3 x1 x5^2 x4^-3 x7^-1 x6^2 x5^-2 x4^-3 x6^-1 x5^-1 x4^-1 x5^-1 x4^-1 x3^-2 x6^-1 x5^-1 x4^-1 x3^-1 x2^-1 x7^2 x6^2 x5 x4^-1 x3 x2^-1 x1^-1 x4 x3 x7^-2 x6^-1 x5^2 x4^-1 x7^-1 x6^-1 x5^2 x4^-2 x3^-1 x2 x5^-1 x4^-1 x3^-1 x7^2 x1^-1 x2^-1 x3^-1 x2^-1 x4^-2 x5^-2 x4^-1 x5^-1 x2 x1 x3 x2 x4 x2^-1 x1^-1 x4^-1 x3^-1 x2^-1 x7 x6^-1 x5^-1 x4^-1 x3^-1 x2 x1^2 x7^-1 x6^-1 x5^-1 x4^-1 x5 x6 x7^-1 x6 x2 x1^2 x2^-1 x1 x3 x2 x1 x3^-1 x2^-1 x1^-1 x6^-1 x1 x2 x3 x2^-1 x5 x2^-1 x3 x5^-1 x3^-1 x1^-1 x3 x5 x3^-1 x5^-1"_w;

  // 2. Run LBA minimization
  TTPLBA ttpLBA;
  TTPTuple red_T;
  bool red_res = ttpLBA.reduce(N, BS, T, gens, 3600 * 2, cout, red_T);
  cout << "!!! Huge success !!!" << endl;
}

int main() {
  // test_special_instance1();
  // return 0;

  RandLib::ur.reset();
  long s1, s2;
  RandLib::ur.getseed( s1, s2 );
  cout << "Seed : " << s1 << " " << s2 << endl;
  
  // getDecompositionDetails();
  // compareGenerateWordFunctions();
  // abelinizationTest();
  // return 0;

  TTP_Conf ttp_conf;
  

  // AE suggested parameters
  ttp_conf.nBL    =  3; // 5;     // # Generators in BL
  ttp_conf.nBR    =  3; // 5;     // # Generators in BR
  ttp_conf.N      = ttp_conf.nBL + ttp_conf.nBR + 2; // 12;    // Group rank
  ttp_conf.nGamma = 10; // Tuple size
  
  ttp_conf.len_z  = 500; // 18;  // Conjugator's length
  ttp_conf.len_w  = 500;  // Word's length
  
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
