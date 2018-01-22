// Copyright (C) 2007 Alex Myasnikov
// Contents: Implementation of TTP  Attack
//
// Principal Authors: Alex Myasnikov
//
// Revision History:
//


#include <set>
#include <list>
#include <time.h>

#include "TTPAttack.h"
#include "BraidGroup.h"
#include "ShortBraidForm.h"
#include "Permutation.h"
#include "errormsgs.h"
#include "ThLeftNormalForm.h"
#include <thread>
#include <mutex>
#include <fstream>

#define TEST_EQUIVALENCE

//
//
//    LBA
//
//


void TTPLBA::addNewElt(const TTPTuple& T, const set<NODE>& checkedElements, set<NODE>& uncheckedElements) {
  NODE new_node(T.length(), T);

  if (checkedElements.count(new_node) != 0) {
    return;
  }

  if (uncheckedElements.count(new_node) != 0) {
    return;
  }

  uncheckedElements.insert(new_node);
}


bool TTPLBA::process_conjugates(int N, const NODE& cur, const vector<Word>& gens,
                        const set<NODE>& checkedElements,
                        set<NODE>& uncheckedElements) {
  vector<TTPTuple> new_tuples(gens.size());

  vector<std::thread> threads;
  threads.reserve(gens.size());

  std::mutex mtx;

  for (int i = 0; i < gens.size(); ++i) {
    threads.emplace_back([&mtx, &cur, &gens, &new_tuples, N, i] () {
      mtx.lock();
      auto new_tuple = cur.second.conjugate(N, gens[i]);
      mtx.unlock();

      new_tuple.shorten(N);
      new_tuples[i] = new_tuple;
    });
  }

  for (auto& th : threads) {
    th.join();
  }

  bool progress = false;

  for (const auto &new_tuple : new_tuples) {
    addNewElt(new_tuple, checkedElements, uncheckedElements);
    progress |= (new_tuple.length() < cur.first);
  }

  return progress;
}

void TTPLBA::tryNode(int N, bool use_special_gens, const NODE& cur, const vector<Word>& gens,
                     const set<NODE>& checkedElements,
                     set<NODE>& uncheckedElements) {
  // 1. Conjugate by a long terminal segments of WL[0] and WR[0].
  // This dramatically reduces weight on the first iterations of the process
  cout << "a0" << endl;
  if (use_special_gens) {
    // @todo Apply special_gens only at the first 4-5 iterations. Then they are useless and can seriously slow down the program.
    vector<Word> special_gens;

    if (cur.second.WL[0].length() >= 200) {
      special_gens.push_back(-(cur.second.WL[0].terminalSegment(19 * cur.second.WL[0].length() / 20)));
    }
    if (cur.second.WR[0].length() >= 200) {
      special_gens.push_back(-(cur.second.WR[0].terminalSegment(19 * cur.second.WR[0].length() / 20)));
    }
    //for (const auto &w : cur.second.WL) {
    //  if (w.length() >= 20) {
    //    special_gens.push_back(-(w.terminalSegment(w.length() - 5)));
    //  }
    //}
    //for (const auto &w : cur.second.WR) {
    //  if (w.length() >= 20) {
    //    special_gens.push_back(-(w.terminalSegment(w.length() - 5)));
    //  }
    //}
    if (cur.second.WL[1].length() > 20) {
      special_gens.push_back(-(cur.second.WL[1].terminalSegment(cur.second.WL[1].length() - 5)));
    }
    if (cur.second.WR[1].length() > 20) {
      special_gens.push_back(-(cur.second.WR[1].terminalSegment(cur.second.WR[1].length() - 5)));
    }
    // special_gens.push_back("x1 x2 x3 x4 x5 x6 x7"_w);

    if (!special_gens.empty()) {
      if (process_conjugates(N, cur, special_gens, checkedElements,
                             uncheckedElements)) {
        cout << "a1" << endl;
        return;
      }
    }
  }

  // 2. Process all conjugates
  cout << "a2" << endl;
  if (process_conjugates(N, cur, gens, checkedElements, uncheckedElements)) {
    // return;
  }

  // 3. Try to fix Delta^2 power in WL
  cout << "a3" << endl;
  const auto new_tuple = cur.second.multiplyElementsByDeltaSQtoReduceLength(N, 3, false);
  addNewElt(new_tuple, checkedElements, uncheckedElements);
}

static void gen_distribution(int N, const Word &w) {
  cout << "[";
  vector<int> dist(N - 1, 0);
  for (const auto i : w) {
    dist[abs(i) - 1]++;
  }
  for (int i = 0; i < N - 1; ++i) {
    cout.width(3);
    cout << dist[i] << ".";
  }
  cout << "] -> " << w.length();
}

static ostream &printVectorOfWords(const vector<Word> &vec, ostream &os) {
  os << "{";
  for (const auto &w : vec) {
    os << "\"" << w << "\"_w, ";
  }
  os << "};";
  return os;
}

static ostream &printVectorOfInts(const vector<int> &vec, ostream &os) {
  os << "{";
  for (const auto &w : vec) {
    os << w << ", ";
  }
  os << "};";
  return os;
}

static void saveDifficultInstance(const int N, const vector<Word> &gens, const BSets &BS, const TTPTuple &t) {
  ofstream of("bad_example.txt", ios::app);
  of << "const int N = " << N << ";" << endl;
  printVectorOfWords(gens, of << "const vector<Word> gens = ") << endl;
  of << "BSets BS;" << endl;
  printVectorOfWords(BS.BL, of << "BS.BL = ") << endl;
  printVectorOfWords(BS.BR, of << "BS.BR = ") << endl << endl;
  of << "TTPTuple T;" << endl;
  printVectorOfWords(t.WL, of << "T.WL = ") << endl;
  printVectorOfWords(t.WR, of << "T.WR = ") << endl;
  printVectorOfInts(t.deltaSQL, of << "T.deltaSQL = ") << endl;
  printVectorOfInts(t.deltaSQR, of << "T.deltaSQR = ") << endl;
  of << "T.z = \"" << t.z << "\"_w;" << endl << endl;
}

bool TTPLBA::reduce(int N, const BSets &bs, const TTPTuple &theTuple,
                    const vector<Word> &gens, int sec, ostream &out,
                    TTPTuple &red_T) {
  int init_time = time(0);
  int maxIterations = 100000;

  set<NODE> checkedElements;
  set<NODE> uncheckedElements;

  // TTPTuple initTuple(theTuple.WL, theTuple.WR, Word());
  TTPTuple initTuple = theTuple;
  int best_result = initTuple.length();
  NODE init(best_result, initTuple);
  size_t stuck_check = 0;

  uncheckedElements.insert(init);
  out << "Initial length: " << best_result << endl;

  for (int c = 0; !uncheckedElements.empty() && c < maxIterations; ++c) {
    // Pick the best unprocessed node
    NODE cur = *uncheckedElements.begin();
    uncheckedElements.erase(uncheckedElements.begin());
    checkedElements.insert(cur);

    // Output some data
    for (const auto&w : cur.second.WL) {
      gen_distribution(N, w);
      cout << endl;
    }
    cout << endl;
    for (const auto&w : cur.second.WR) {
      gen_distribution(N, w);
      cout << endl;
    }

    int cur_time = time(0);
    if (best_result > cur.first) {
      best_result = cur.first;
      stuck_check = 0;
    } else {
      if (++stuck_check > 50) {
        // We are officially stuck. Save the instance to process later
        // I think we need 2 saves: (a) the original instance as it was originally generated and (b) the reduced one to start LBA from that point
        cout << " >>> STUCK <<< " << endl;

        if (checkedElements.size() % 50 == 0) {
          const auto &best = *checkedElements.begin();
          for (const auto &w : best.second.WL) {
            cout << "> " << endl;
            cout << w << endl;
            gen_distribution(N, w);
            cout << endl;
          }
          cout << endl;
          for (const auto &w : best.second.WR) {
            cout << "> " << endl;
            cout << w << endl;
            gen_distribution(N, w);
            cout << endl;
          }
        }
      }
    }

    out << "Current (best) length: " << cur.first << " (" << best_result << ")   ";
    cur.second.printPowers();
    cout << "   tm = " << cur_time << endl;

    if (cur_time - init_time > sec) {
      cout << "Failed example!" << endl;
      saveDifficultInstance(N, gens, bs, checkedElements.begin()->second);
      exit(1);
      return false; // TIME_EXPIRED;
    }

    // Termination condition: check that cur.second.WL and cur.second.WR are separated
    // if (cur.second.testTuples(N, false)) {
    if (cur.second.shortAndTestTuples(N)) {
      // (debug)
#ifdef TEST_EQUIVALENCE
      if (!theTuple.equivalent(N, cur.second)) {
        cout << "Internal check failure in TTPLBA::reduce" << endl;
        exit(1);
      }
#endif
      red_T = cur.second;
      return true;
    }

    // tryNode(N, c < 5, cur, gens, checkedElements, uncheckedElements);
    tryNode(N, true, cur, gens, checkedElements, uncheckedElements);
  }

  return false; // FAILED;
}

  

///////////////////////////////////////////////////////////////////////////////////////
//
//  TTP ATTACK
//
///////////////////////////////////////////////////////////////////////////////////////

bool TTPAttack::run(const TTPTuple &original_tuple) {
  BraidGroup B(N);
  typedef ThLeftNormalForm NF;
  const Word& original_z = original_tuple.origZ;

  // 1. (part of generation) Take all braids modulo Delta^2 and construct a single tuple of elements
  TTPTuple init_tuple = original_tuple.takeModuloDeltaSQ(N);

  // (debug)
#ifdef TEST_EQUIVALENCE
  if (!original_tuple.equivalent(N, init_tuple)) {
    cout << "Internal failure in takeModuloDeltaSQ" << endl;
    exit(1);
  }
#endif

  // 2. (attack) Attempt to restore the original delta powers
  cout << "Attempt to restore the original delta powers" << endl;
  const auto tuple1 = init_tuple.multiplyElementsByDeltaSQtoReduceLength(N, 0, true);
  // (debug)
#ifdef TEST_EQUIVALENCE
  if (!original_tuple.equivalent(N, tuple1)) {
    cout << "Internal failure in multiplyElementsByDeltaSQtoReduceLength" << endl;
    exit(1);
  }
#endif

  // (debug) Here we test if we found correct powers of Delta^2
  tuple1.printPowers();
  cout << endl;

  // 3. Apply length-based conjugacy minimization
  return LBA(original_tuple.WL.size(), original_tuple.WR.size(), tuple1, original_z);
}

bool TTPAttack::LBA(int NWL, int NWR, const TTPTuple &t, const Word &z) {
  // 1. Create the initial tuple
  TTPTuple T = t;
  T.shorten_parallel(N);
  T.printPowers();
  cout << endl;

  // 2. Prepare the generators
  vector<Word> gens;
  gens.reserve(2 * N);

  for (int i = 1; i < N; ++i) {
    gens.push_back(Word(i));
    gens.push_back(Word(-i));
  }

  // 3. Run LBA minimization
  TTPLBA ttpLBA;
  TTPTuple red_T;
  // bool red_res = ttpLBA.reduce(N, BS, T, gens, 3600 * 2, cout, red_T);
  bool red_res = ttpLBA.reduce(N, BS, T, gens, 3600, cout, red_T);

  // (debug) If LBA minimization is successful, then check correctness of computations and check if we got the original z
  if (red_res) {
    cout << "Same z: |z*z'^-1| = " << shortenBraid(N, z * red_T.z).length() << endl;
    // Checking correctness of computations (check if red_T.z is the conjugator)
    if (!t.equivalent(N, red_T)) {
      cout << "Internal LBA check failed" << endl;
      exit(1);
    }
  } else {
  }
  return red_res;
}

bool TTPAttack::oneOfSSSReps(int NWL, int NWR, const vector<ThLeftNormalForm> &theTuple) {

  for (int i = 0; i<theTuple.size(); i++) {
    //	  printStats(theTuple[i],cout); cout << " --> ";
    pair<ThLeftNormalForm, ThLeftNormalForm> sssR = theTuple[i].findSSSRepresentative();
    //		printStats(sssR.first,cout); cout << endl;



    // do test

    // convert and separate tuples 
    TTPTuple check_tuples;
    check_tuples.WL = vector<Word>(NWL);
    check_tuples.WR = vector<Word>(NWR);

    for (int j = 0; j<theTuple.size(); j++)
      if (j < NWL)
        check_tuples.WL[j] = (-sssR.second*theTuple[j] * sssR.second).getWord();
      else
        check_tuples.WR[j - NWL] = (-sssR.second*theTuple[j] * sssR.second).getWord();

    if (check_tuples.shortAndTestTuples(N)) {
      //      cout << "SUCC: " << endl; // << sssR.second << endl;
      return true;
    }

  }

  //  cout << "FAIL" << endl;  
  return false;
}