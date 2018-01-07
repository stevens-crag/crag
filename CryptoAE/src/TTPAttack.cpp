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

//typedef pair< vector< Word > , vector< Word > > TTPTuple;





//
//
//  SOME CONDITIONS
//
//

bool equalUpToCommut( int N, const BSets& bs, const TTPTuple& ttp, const Word& z)
{
  
  Word z_diff = shortenBraid(N,z*ttp.z);


  if (z_diff.length() == 0)
    return true;

  
  // It seems that  z can be hard to recover 
  // because the difference z*z' may commute with both BL and BR
  // Then |BL| + |BR| = |(z*z')^-1 BL (z*z')| + |(z*z')^-1 BR (z*z')|
  // Check it here 
  
  bool commute = true; 
  for ( int i=0;i<bs.BL.size();i++)
    if (shortenBraid(N,-z_diff*bs.BL[i]*z_diff) != bs.BL[i])
      commute = false;
  
  for ( int i=0;i<bs.BR.size();i++)
    if (shortenBraid(N,-z_diff*bs.BR[i]*z_diff) != bs.BR[i])
      commute = false;


  if (commute)
    cout << "COMMUTE DIFF : ";
  else {
    cout << "TROUBLE DIFF : ";
    
    // OUTPUT TROUBLED THING
    
    cout << "TROUBLED  : " << endl
	 << "Z  : " << z  << endl
         << "Z' : " << ttp.z << endl
         << "DIFF : " << z_diff << endl;
    
    cout << " WL : " << endl;
    for ( int i=0;i<ttp.WL.size();i++){
      cout << ttp.WL[i] << endl;
    }
    cout << " WR : " << endl;
    for ( int i=0;i<ttp.WR.size();i++){
      cout << ttp.WR[i] << endl;
    }     
  }

  cout << z_diff.length() << endl;

  return commute;

  /*
    // CHECKS IF THE CONJ DOES NOT CHANGE THE LENGTH
    //
  int len = 0;
  int len_conj = 0;

  //  cout << "TROUBLED  : " 
  //       << "Z  : " << z  << endl
  //       << "Z' : " << ttp.z << endl
  //       << "DIFF : " << z_diff << endl;

  //  cout << " WL : " << endl;
  for ( int i=0;i<ttp.WL.size();i++){
    len += ttp.WL[i].length();
    len_conj += shortenBraid(N,-z_diff*ttp.WL[i]*z_diff).length();
    //    cout << ttp.WL[i] << endl;
  }
  //  cout << " WR : " << endl;
  for ( int i=0;i<ttp.WR.size();i++){
    len += ttp.WR[i].length();
    len_conj += shortenBraid(N,-z_diff*ttp.WR[i]*z_diff).length();
    //    cout << ttp.WR[i] << endl;
  }     
  
  

  if ( len_conj <= len)
    return true;
  else
    return false;
  */
}

//
//
//    LBA
//
//


void TTPLBA::addNewElt(const TTPTuple &T, const set<NODE> &checkedElements, set<NODE> &uncheckedElements) {
  int weight = T.length();
  NODE new_node(weight, T);
  if (checkedElements.find(new_node) != checkedElements.end())
    return;
  if (uncheckedElements.find(new_node) != uncheckedElements.end())
    return;
  uncheckedElements.insert(new_node);
}


void TTPLBA::tryNode(int N, NODE cur, const vector<Word>& gens,
                     const set<NODE> &checkedElements,
                     set<NODE> &uncheckedElements) {
  typedef ThLeftNormalForm NF;
  BraidGroup B(N);

  // Conjugate by a long terminal segments of WL[0] and WR[0].
  // This dramatically reduces weight on the first iterations of the process
  vector<Word> special_gens = {
      cur.second.WL[0].terminalSegment(9 * cur.second.WL[0].length() / 10),
      cur.second.WR[0].terminalSegment(9 * cur.second.WR[0].length() / 10),
      cur.second.WL[0].terminalSegment(cur.second.WL[0].length() - 3),
      cur.second.WR[0].terminalSegment(cur.second.WR[0].length() - 3),
  };
  for (const auto& b : special_gens) {
    auto new_tuple = cur.second.conjugate(N, -b);
    new_tuple.shorten(N);
    cout << cur.first << " -> " << new_tuple.length() << endl;
    if (new_tuple.length() < cur.first) {
      addNewElt(new_tuple, checkedElements, uncheckedElements);
      return;
    }
  }

  // Process all conjugates
  bool progress = false;
  for (const auto& b : gens) {
    auto new_tuple = cur.second.conjugate(N, b);
    new_tuple.shorten(N);
    addNewElt(new_tuple, checkedElements, uncheckedElements);
    progress |= (new_tuple.length() < cur.first);
  }
  if (progress)
    return;

  // Try to fix Delta^2 power in WL
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

  savTuple = initTuple;

  uncheckedElements.insert(init);
  out << "Initial length: " << best_result << endl;

  for (int c = 0; uncheckedElements.size() && c < maxIterations; ++c) {

    NODE cur = *uncheckedElements.begin();
    uncheckedElements.erase(uncheckedElements.begin());
    checkedElements.insert(cur);

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
    if (cur.second.shortAndTestTuples(N)) {
      // (debug)
      if (!theTuple.equivalent(N, cur.second)) {
        cout << "Internal check failure in TTPLBA::reduce" << endl;
        exit(1);
      }
      red_T = cur.second;
      return true;
    }
    tryNode(N, cur, gens, checkedElements, uncheckedElements);
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
  T.shorten(N);
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

