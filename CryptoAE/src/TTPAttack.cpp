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

static bool _shorterWordFound;
static vector<Word> _word;
static vector<int> _power;

static void run_thread_multiplyElementByPMDeltaSQtoReduceLength(const Word &w, int N, int i) {
  BraidGroup B(N);

  ThLeftNormalForm nf(B, w);
  const auto p = nf.getPower();

  nf.setPower(p - 2);
  const auto w1 = shortenBraid(N, nf.getReducedWord2());

  nf.setPower(p + 2);
  const auto w2 = shortenBraid(N, nf.getReducedWord2());

  if (w.length() <= w1.length() && w.length() <= w2.length()) {
    _word[i] = w;
    _power[i] = 0;
  } else if (w1.length() <= w2.length()) {
    _word[i] = w1;
    _power[i] = -1;
    _shorterWordFound = true;
  } else {
    _word[i] = w2;
    _power[i] = 1;
    _shorterWordFound = true;
  }
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

    if (cur.second.WL[0].length() >= 60) {
      special_gens.push_back(-(cur.second.WL[0].terminalSegment(9 * cur.second.WL[0].length() / 10)));
    }

    if (cur.second.WR[0].length() >= 60) {
      special_gens.push_back(-(cur.second.WR[0].terminalSegment(9 * cur.second.WR[0].length() / 10)));
    }
     if (cur.second.WL[0].length() > 20) {
       special_gens.push_back(-(cur.second.WL[0].terminalSegment(cur.second.WL[0].length() - 5)));
     }
     if (cur.second.WR[0].length() > 20) {
       special_gens.push_back(-(cur.second.WR[0].terminalSegment(cur.second.WR[0].length() - 5)));
     }

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
    return;
  }

  // 3. Try to fix Delta^2 power in WL
  {
    cout << "a3" << endl;
    const auto nWL = cur.second.WL.size();
    const auto nWR = cur.second.WR.size();
    _word = vector<Word>(nWL + nWR);
    _power = vector<int>(nWL + nWR, 0);
    _shorterWordFound = false;
    vector<std::thread> threads;
    cout << "a4" << endl;

    for (int i = 0; i < nWL; ++i) {
      threads.push_back(std::thread(run_thread_multiplyElementByPMDeltaSQtoReduceLength, cur.second.WL[i], N, i));
    }

    for (int i = 0; i < nWR; ++i) {
      threads.push_back(std::thread(run_thread_multiplyElementByPMDeltaSQtoReduceLength, cur.second.WR[i], N, nWL + i));
    }

    for (auto& th : threads) {
      th.join();
    }

    cout << "a5" << endl;
    if (_shorterWordFound) {
      auto new_tuple = cur.second;
      new_tuple.WL = vector<Word>(_word.begin(), _word.begin() + nWL);
      new_tuple.WR = vector<Word>(_word.begin() + nWL, _word.end());
      for (auto i = 0; i < nWL; ++i)
        new_tuple.deltaSQL[i] += _power[i];
      for (auto i = 0; i < nWR; ++i)
        new_tuple.deltaSQR[i] += _power[nWL + i];
      addNewElt(new_tuple, checkedElements, uncheckedElements);
    }
  }
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

        const auto &best = *checkedElements.begin();
        for (const auto&w : best.second.WL) {
          cout << "> " << endl;
          cout << w << endl;
          gen_distribution(N, w);
          cout << endl;
        }
        cout << endl;
        for (const auto&w : best.second.WR) {
          cout << "> " << endl;
          cout << w << endl;
          gen_distribution(N, w);
          cout << endl;
        }
      }
    }

    out << "Current (best) length: " << cur.first << " (" << best_result << ")   ";
    cur.second.printPowers();
    cout << "   tm = " << cur_time << endl;

    if (cur_time - init_time > sec)
      return false; // TIME_EXPIRED;

    // Termination condition: check that cur.second.WL and cur.second.WR are separated
    // if (cur.second.testTuples(N, false)) {
    if (cur.second.shortAndTestTuples(N)) {
      // (debug)
      if (!theTuple.equivalent(N, cur.second)) {
        cout << "Internal check failure in TTPLBA::reduce" << endl;
        exit(1);
      }
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
  // if (!original_tuple.equivalent(N, init_tuple)) {
  //  cout << "Internal failure in takeModuloDeltaSQ" << endl;
  //  exit(1);
  //}

  // 2. (attack) Attempt to restore the original delta powers
  cout << "Attempt to restore the original delta powers" << endl;
  const auto tuple1 = multiplyElementsByDeltaSQtoReduceLength(init_tuple);
  // (debug)
  // if (!original_tuple.equivalent(N, tuple1)) {
  //  cout << "Internal failure in multiplyElementsByDeltaSQtoReduceLength" <<
  //  endl; exit(1);
  //}

  // (debug) Here we test if we found correct powers of Delta^2
  tuple1.printPowers();
  cout << endl;

  // 3. Apply length-based conjugacy minimization
  return LBA(original_tuple.WL.size(), original_tuple.WR.size(), tuple1, original_z);
}

bool TTPAttack::oneOfSSSReps( int NWL, int NWR, const vector<ThLeftNormalForm>& theTuple )
{ 
  
  for ( int i=0;i<theTuple.size(); i++) {
    //	  printStats(theTuple[i],cout); cout << " --> ";
    pair<ThLeftNormalForm,ThLeftNormalForm> sssR = theTuple[i].findSSSRepresentative();
    //		printStats(sssR.first,cout); cout << endl;
    
    
    
    // do test
    
    // convert and separate tuples 
    TTPTuple check_tuples;
    check_tuples.WL = vector<Word>(NWL);
    check_tuples.WR = vector<Word>(NWR);
    
    for ( int j=0;j<theTuple.size();j++)
      if ( j < NWL )
	check_tuples.WL[j] = (-sssR.second*theTuple[j]*sssR.second).getWord();
      else
	check_tuples.WR[j-NWL] = (-sssR.second*theTuple[j]*sssR.second).getWord();
    
    if ( check_tuples.shortAndTestTuples( N ) ) {
//      cout << "SUCC: " << endl; // << sssR.second << endl;
      return true;
    }
    
  }

//  cout << "FAIL" << endl;  
  return false;	
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
  bool red_res = ttpLBA.reduce(N, BS, T, gens, 3600 * 2, cout, red_T);

  // (debug) If LBA minimization is successful, then check correctness of computations and check if we got the original z
  if (red_res) {
    cout << "Same z: |z*z'^-1| = " << shortenBraid(N, z * red_T.z).length() << endl;
    // Checking correctness of computations (check if red_T.z is the conjugator)
    if (!t.equivalent(N, red_T)) {
      cout << "Internal LBA check failed" << endl;
      exit(1);
    }
  }
  return red_res;
}

/*
//! Find a power Delta^2p s.t. |Delta^2p*nf| is minimal in the <Delta^2>-coset
static int multiplyByDeltaSQtoReduceLength(int N, Word &w) {
  typedef ThLeftNormalForm NF;
  BraidGroup B(N);

  NF nf(B, w);
  // -nf.getDecomposition().size()/4 is the anticipated value
  int best_nf = -nf.getDecomposition().size()/4;
  map<int, int> weights;
  static const int delta = 4;

  bool progress = true;
  while (progress) {
    progress = false;
    for (int i = best_nf - delta; i <= best_nf + delta; ++i) {
      if (weights.find(i) == weights.end()) {
        ThLeftNormalForm nf2 = nf;
        nf2.setPower(nf.getPower() + 2 * i);
        const auto cur_w = shortenBraid(N, nf2.getReducedWord2());
        weights[i] = cur_w.length();
        if (weights.find(best_nf) == weights.end() || weights[i] < weights[best_nf]) {
          progress = true;
          best_nf = i;
          w = cur_w;
        }
        cout << weights[i] << ", ";
      }
    }
  }
  nf.setPower(nf.getPower() + 2 * best_nf);
  for (const auto&p : weights) {
    // cout << p.second << ", ";
    // cout << "(" << p.first << "," << p.second << "), ";
  }
  cout << endl;
  return best_nf;
}
*/

//! Find a power Delta^2p s.t. |Delta^2p*nf| is minimal in the <Delta^2>-coset
static int multiplyByDeltaSQtoReduceLength(int N, Word &w) {
  typedef ThLeftNormalForm NF;
  BraidGroup B(N);

  NF nf(B, w);
  // -nf.getDecomposition().size()/4 is the anticipated value
  int best_nf = -nf.getDecomposition().size()/4;
  map<int, Word> weights;
  static const int delta = 2;
  const Permutation omega = Permutation::getHalfTwistPermutation(N);
  Word omegaWord = omega.geodesicWord();

  // initial data
  {
    ThLeftNormalForm nf2 = nf;
    nf2.setPower(nf.getPower() + 2 * best_nf);
    weights[best_nf] = shortenBraid(N, nf2.getReducedWord2());
    w = weights[best_nf];
    cout << weights[best_nf].length() << ", ";
  }

  bool progress = true;
  while (progress) {
    progress = false;

    for (int i = best_nf + 1; i <= best_nf + delta; ++i) {
      if (weights.find(i) != weights.end())
        continue;
      const auto &cur_w = weights[i - 1];
      const auto new_w = shortenBraid(N, omegaWord.power(2) * cur_w);
      weights[i] = new_w;
      if (weights[i].length() < weights[best_nf].length()) {
        progress = true;
        best_nf = i;
        w = new_w;
        cout << "!";
      }
      cout << weights[i].length() << ", ";
    }

    for (int i = best_nf - 1; i >= best_nf - delta; --i) {
      if (weights.find(i) != weights.end())
        continue;
      const auto &cur_w = weights[i + 1];
      const auto new_w = shortenBraid(N, omegaWord.power(-2) * cur_w);
      weights[i] = new_w;
      if (weights[i].length() < weights[best_nf].length()) {
        progress = true;
        best_nf = i;
        w = new_w;
        cout << "!";
      }
      cout << weights[i].length() << ", ";
    }
  }
  for (const auto&p : weights) {
    // cout << p.second << ", ";
    // cout << "(" << p.first << "," << p.second << "), ";
  }
  cout << endl;
  return best_nf;
}

TTPTuple TTPAttack::multiplyElementsByDeltaSQtoReduceLength(const TTPTuple &t) {
  vector<std::thread> threads;
  threads.reserve(t.WL.size() + t.WR.size());

  auto result = t;

  for (int i = 0; i < t.WL.size(); ++i) {
    threads.emplace_back([&result, i, this]() {
      result.deltaSQL[i] += multiplyByDeltaSQtoReduceLength(N, result.WL[i]);
    });
  }

  for (int i = 0; i < t.WR.size(); ++i) {
    threads.emplace_back([&result, i, this]() {
      result.deltaSQR[i] += multiplyByDeltaSQtoReduceLength(N, result.WR[i]);
    });
  }

  for (auto& th : threads) {
    th.join();
  }

  return result;
}
