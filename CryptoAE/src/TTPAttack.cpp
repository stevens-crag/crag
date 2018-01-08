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

//
//
//    LBA
//
//


void TTPLBA::addNewElt(const TTPTuple &T, const set<NODE> &checkedElements, set<NODE> &uncheckedElements) {
  NODE new_node(T.length(), T);
  if (checkedElements.find(new_node) != checkedElements.end())
    return;
  if (uncheckedElements.find(new_node) != uncheckedElements.end())
    return;
  uncheckedElements.insert(new_node);
}

static vector<TTPTuple> _tuple;

static void run_thread_conjugate(const TTPTuple& t, const Word& g, int N, int i) {
  auto t1 = t.conjugate(N, g);
  t1.shorten(N);
  _tuple[i] = t1;
}

static bool _shorterWordFound;
static vector<Word> _word;
static vector<int> _power;

static void run_thread_multiplyElementByPMDeltaSQtoReduceLength(const Word &w, int N, int i) {
  typedef ThLeftNormalForm NF;
  BraidGroup B(N);

  NF nf(B, w);
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

bool TTPLBA::process_conjugates(int N, NODE cur, const vector<Word> &gens,
                        const set<NODE> &checkedElements,
                        set<NODE> &uncheckedElements) {
  bool progress = false;
  _tuple = vector<TTPTuple>(gens.size());
  vector<std::thread> threads;
  for (int i = 0; i < gens.size(); ++i) {
    threads.push_back(std::thread(run_thread_conjugate, cur.second, gens[i], N, i));
  }
  for (auto &th : threads) {
    th.join();
  }
  for (const auto &new_tuple : _tuple) {
    addNewElt(new_tuple, checkedElements, uncheckedElements);
    progress |= (new_tuple.length() < cur.first);
  }
  return progress;
}

void TTPLBA::tryNode(int N, NODE cur, const vector<Word>& gens,
                     const set<NODE> &checkedElements,
                     set<NODE> &uncheckedElements) {
  typedef ThLeftNormalForm NF;
  BraidGroup B(N);

  // 1. Conjugate by a long terminal segments of WL[0] and WR[0].
  // This dramatically reduces weight on the first iterations of the process
  cout << "a0" << endl;
  vector<Word> special_gens;
  if (cur.second.WL[0].length() >= 60)
    special_gens.push_back(-(cur.second.WL[0].terminalSegment(9 * cur.second.WL[0].length() / 10)));
  if (cur.second.WR[0].length() >= 60)
    special_gens.push_back(-(cur.second.WR[0].terminalSegment(9 * cur.second.WR[0].length() / 10)));
  //if (cur.second.WL[0].length() > 6)
  //  gens.push_back(-(cur.second.WL[0].terminalSegment(cur.second.WL[0].length() - 3)));
  //if (cur.second.WR[0].length() > 6)
  //  gens.push_back(-(cur.second.WR[0].terminalSegment(cur.second.WR[0].length() - 3)));
  if (special_gens.size() > 0) {
    if (process_conjugates(N, cur, special_gens, checkedElements, uncheckedElements)) {
      cout << "a1" << endl;
      return;
    }
  }

  // 2. Process all conjugates
  cout << "a2" << endl;
  if (process_conjugates(N, cur, gens, checkedElements, uncheckedElements))
    return;

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
    for (auto& th : threads) th.join();
    cout << "a5" << endl;
    if (_shorterWordFound) {
      auto new_tuple = cur.second;
      cout << "f" << endl;
      new_tuple.WL = vector<Word>(_word.begin(), _word.begin() + nWL);
      new_tuple.WR = vector<Word>(_word.begin() + nWL, _word.end());
      for (auto i = 0; i < nWL; ++i)
        new_tuple.deltaSQL[i] += _power[i];
      for (auto i = 0; i < nWR; ++i)
        new_tuple.deltaSQR[i] += _power[nWL + i];
      addNewElt(new_tuple, checkedElements, uncheckedElements);
    }
    cout << "a6" << endl;
  }

  /*
  for (int i = 0; i < cur.second.WL.size(); ++i) {
    const auto w = cur.second.WL[i];
    NF nf(B, w);
    const auto p = nf.getPower();
    for (int d = -1; d <= 1; d += 2) {
      nf.setPower(p + 2 * d);
      const auto w1 = shortenBraid(N, nf.getReducedWord2());
      if (w1.length() < w.length()) {
        TTPTuple new_tuple = cur.second;
        new_tuple.WL[i] = w1;
        new_tuple.deltaSQL[i] += d;
        // (debug)
        //if (!cur.second.equivalent(N, new_tuple)) {
        //  cout << "Internal check failure in tryNode" << endl;
        //  exit(1);
        //} else {
        //  cout << "Internal check in tryNode works WL #" << i << ": " << cur.second.deltaSQL[i] << " -> " << new_tuple.deltaSQL[i] << endl;
        //}
        addNewElt(new_tuple, checkedElements, uncheckedElements);
      }
    }
  }
  // Try to fix Delta^2 power in WR
  for (int i = 0; i < cur.second.WR.size(); ++i) {
    const auto w = cur.second.WR[i];
    NF nf(B, w);
    const auto p = nf.getPower();
    for (int d = -1; d <= 1; d += 2) {
      nf.setPower(p + 2 * d);
      const auto w1 = shortenBraid(N, nf.getReducedWord2());
      if (w1.length() < w.length()) {
        TTPTuple new_tuple = cur.second;
        new_tuple.WR[i] = w1;
        new_tuple.deltaSQR[i] += d;
        // (debug)
        //if (!cur.second.equivalent(N, new_tuple)) {
        //  cout << "Internal check failure in tryNode" << endl;
        //  exit(1);
        //} else {
        //  cout << "Internal check in tryNode works WR #" << i << ": " << cur.second.deltaSQR[i] << " -> " << new_tuple.deltaSQR[i] << endl;
        //}
        addNewElt(new_tuple, checkedElements, uncheckedElements);
      }
    }
  }
  */
}


bool TTPLBA::reduce(int N, const BSets &bs, const TTPTuple &theTuple,
                    const vector<Word> &gens, int sec, ostream &out,
                    TTPTuple &red_T, const Word &z) {
  typedef ThLeftNormalForm NF;
  BraidGroup B(N);

  int init_time = time(0);
  int maxIterations = 100000;

  set<NODE> checkedElements;
  set<NODE> uncheckedElements;

  // TTPTuple initTuple(theTuple.WL, theTuple.WR, Word());
  TTPTuple initTuple = theTuple;
  int best_result = initTuple.length();
  NODE init(best_result, initTuple);

  uncheckedElements.insert(init);
  out << "Initial length: " << best_result << endl;

  for (int c = 0; uncheckedElements.size() && c < maxIterations; ++c) {

    cout << "q0" << endl;
    NODE cur = *uncheckedElements.begin();
    uncheckedElements.erase(uncheckedElements.begin());
    checkedElements.insert(cur);

    cout << "qa1" << endl;
    int cur_time = time(0);
    if (best_result > cur.first)
      best_result = cur.first;
    out << "Current (best) length: " << cur.first << " (" << best_result << ")   ";
    cur.second.printPowers();
    cout << "   tm = " << cur_time << endl;

    if (cur_time - init_time > sec)
      return false; // TIME_EXPIRED;

    // Termination condition: check that cur.second.WL and cur.second.WR are separated
    // if (cur.second.testTuples(N, false)) {
    cout << "q1" << endl;
    if (cur.second.shortAndTestTuples(N)) {
      // (debug)
      cout << "y" << endl;
      if (!theTuple.equivalent(N, cur.second)) {
        cout << "Internal check failure in TTPLBA::reduce" << endl;
        exit(1);
      }
      red_T = cur.second;
      return true;
    }
    cout << "q2" << endl;
    tryNode(N, cur, gens, checkedElements, uncheckedElements);
    cout << "q3" << endl;
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
  typedef ThLeftNormalForm NF;
  BraidGroup B(N);

  // 1. Create the initial tuple
  TTPTuple T = t;
  T.shorten_parallel(N);
  T.printPowers(); cout << endl;

  // 2. Prepare the generators
  vector<Word> gens;
  for (int i = 1; i < N; ++i) {
    gens.push_back(Word(i));
    gens.push_back(Word(-i));
  }

  // 3. Run LBA minimization
  TTPLBA ttpLBA;
  TTPTuple red_T;
  bool red_res = ttpLBA.reduce(N, BS, T, gens, 3600 * 12, cout, red_T, z);

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


/*
TTPTuple TTPAttack::multiplyElementsByDeltaSQtoReduceLength(const TTPTuple &t) {
TTPTuple result = t;
for (int i = 0; i < t.WL.size(); ++i) {
result.deltaSQL[i] += multiplyByDeltaSQtoReduceLength(N, result.WL[i]);
}
for (int i = 0; i < t.WR.size(); ++i) {
result.deltaSQR[i] += multiplyByDeltaSQtoReduceLength(N, result.WR[i]);
}
return result;
}
*/



static TTPTuple tmp_tuple;

static void run_thread1(int N, int i) {
  tmp_tuple.deltaSQL[i] += multiplyByDeltaSQtoReduceLength(N, tmp_tuple.WL[i]);
}
static void run_thread2(int N, int i) {
  tmp_tuple.deltaSQR[i] += multiplyByDeltaSQtoReduceLength(N, tmp_tuple.WR[i]);
}

TTPTuple TTPAttack::multiplyElementsByDeltaSQtoReduceLength(const TTPTuple &t) {
  vector<std::thread> threads;
  tmp_tuple = t;
  for (int i = 0; i < t.WL.size(); ++i) {
    threads.push_back(std::thread(run_thread1, N, i));
  }
  for (int i = 0; i < t.WR.size(); ++i) {
    threads.push_back(std::thread(run_thread2, N, i));
  }
  for (auto& th : threads) {
    th.join();
  }
  return tmp_tuple;
}
