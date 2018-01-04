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

void TTPLBA::tryNode(int N, NODE cur, const vector<Word> &gens,
                     const set<NODE> &checkedElements,
                     set<NODE> &uncheckedElements, int &min_weight) {
  // This was used to exclude nodes that are too far from the current optimal
  // we better vary this value, depending on parameters of A and B
  int MAX_DELTA = 80;

  for (int i = 0; i < gens.size(); ++i) {

    int new_weight = 0;
    Word b = gens[i];
    for (int d = 0; d < 2; ++d) {

      // cout << "   #" << i+1 << "," << d+1;
      if (d == 1)
        b = -b;
      bool candidate = true;
      int delta = 0;
      vector<Word> WL = cur.second.WL;
      for (int t = 0; t < WL.size() && candidate; ++t) {
        delta -= WL[t].length();
        WL[t] = shortenBraid(N, -b * WL[t] * b);
        new_weight += WL[t].length();
        delta += WL[t].length();
      }

      vector<Word> WR = cur.second.WR;
      for (int t = 0; t < WR.size() && candidate; ++t) {
        delta -= WR[t].length();
        WR[t] = shortenBraid(N, -b * WR[t] * b);
        new_weight += WR[t].length();
        delta += WR[t].length();
      }

      //      if( new_weight > min_weight + MAX_DELTA ) {
      //      if( delta>= 0 ) { // MAX_DELTA ) {
      //	candidate = false;
      // cout << " stopped @ " << t << endl;
      //      }

      if (new_weight < min_weight) {
        min_weight = new_weight;
        // cout << "New min w : " << min_weight << endl;
      }

      if (candidate) {
        // cout << "Adding candidate " << b*cur.second.z << " : " <<  endl;
        Word conj = cur.second.z * b;
        for (int i = 0; i < WL.size(); i++) {
          //  			cout << ThLeftNormalForm( BraidGroup( N )  ,
          //  -savTuple.WL[i]*conj*WL[i]*-conj );
          // 			cout << WL[i] << " -----> " << endl;
          //       		cout <<
          //       shortenBraid(N,-savTuple.WL[i]*conj*WL[i]*-conj) << endl;
          // 			cout << ThLeftNormalForm( BraidGroup( N )  ,
          // -savTuple.WR[i]*conj*WR[i]*-conj ); 			cout <<
          // WR[i] << " -----> " << endl;
          //  			cout << shortenBraid(N,-savTuple.WR[i]*conj*WR[i]*-conj) <<
          //  endl;
        }
        addNewElt(TTPTuple(WL, WR, conj), checkedElements, uncheckedElements);
        // cout << " accepted with " << delta << endl;
      }
    }
  }
}

bool TTPLBA::reduce(int N, const BSets &bs, const TTPTuple &theTuple,
                    const vector<Word> &gens, int sec, ostream &out,
                    TTPTuple &red_T, const Word &z) {

  int maxIterations = 100000;

  set<NODE> checkedElements;
  set<NODE> uncheckedElements;

  int init_time = time(0);

  TTPTuple initTuple(theTuple.WL, theTuple.WR, Word());
  //  cout << "Start LBR. Orig tuple:" << endl;
  for (int i = 0; i < initTuple.WL.size(); i++) {
    //  		cout << ThLeftNormalForm( BraidGroup( N )  ,
    //  initTuple.WL[i] );
    //    cout << shortenBraid(N,initTuple.WL[i])<< endl;
    // 		cout << ThLeftNormalForm( BraidGroup( N )  , initTuple.WR[i] );
    //    cout << shortenBraid(N,initTuple.WR[i]) << endl;
  }

  //  int init_weight1 = sbgpGeneratorsWeight( A1 );
  int init_weight = initTuple.length();
  NODE init(init_weight, initTuple);

  int min_weight = init_weight;

  savTuple = initTuple;

  uncheckedElements.insert(init);
  out << "Initial length: " << init_weight << endl;
  int best_result = 999999;

  for (int c = 0; uncheckedElements.size() && c < maxIterations; ++c) {

    NODE cur = *uncheckedElements.begin();
    uncheckedElements.erase(uncheckedElements.begin());
    checkedElements.insert(cur);

    int cur_time = time(0);
    if (best_result > cur.first)
      best_result = cur.first;
    out << "Current (best) length: " << cur.first << " (" << best_result << ")"
        << ", tm = " << cur_time << endl;

    if (cur_time - init_time > sec)
      return false; // TIME_EXPIRED;

    //    if( equalUpToCommut(N,cur.second,z) ) {
    if (cur.second.shortAndTestTuples(N)) {
      red_T = cur.second;

      // CHECK THAT THE Z IS ACTUALL CONJUGATOR
      bool same_z = true;
      for (int i = 0; i < red_T.WL.size(); i++) {
        if (shortenBraid(N, cur.second.z * cur.second.WL[i] * -cur.second.z *
                                -theTuple.WL[i])
                .length() > 0) {
          same_z = false;
          break;
        }
      }
      for (int i = 0; i < red_T.WR.size(); i++) {
        if (shortenBraid(N, cur.second.z * cur.second.WR[i] * -cur.second.z *
                                -theTuple.WR[i])
                .length() > 0) {
          same_z = false;
          break;
        }
      }

      if (!same_z)
        cout << "LBA CHECK FOR Z FAILED" << endl;

      // END CHECK

      return simpleLBA(N, bs, cur.second, shortenBraid(N, z * cur.second.z));
      //      return true; // SUCCESSFULL;
    }

    tryNode(N, cur, gens, checkedElements, uncheckedElements, min_weight);
  }

  return false; // FAILED;
}

bool TTPLBA::simpleLBA( int N , const BSets& bs, const TTPTuple& theTuple, const Word& z, TTPTuple* ret_T )
{
  
  TTPTuple T = theTuple;
  
  T.shorten( N );
  
  vector< Word > gens(N-1);
  for (int i=0;i<N-1;i++)
    gens[i] = Word(i+1);
  

  // DO THE REDUCTION HERE
  TTPTuple red_T( T.WL,T.WR );
  Word     z_conj;

  int minLength = 999999999;
  while ( 1) { // !red_T.shortAndTestTuples( N ) ) {

    Word succGen;
    bool lenReduced = false;
    // int minLength = 999999999;
    for ( int i=0;i<gens.size();i++) {
      
      // switch for positive negative
      Word g = gens[i];
      for (int s=0;s<=1;s++){
	if (s) g = -g;
	
	int length=0;
	// DO LEFT TUPLE
	for ( int iL=0;iL<T.WL.size();iL++)
	  length += shortenBraid(N,-g*red_T.WL[iL]*g).length();

	// DO RIGTH TUPLE
	for ( int iR=0;iR<T.WR.size();iR++)
	  length += shortenBraid(N,-g*red_T.WR[iR]*g).length();	
	

	//	cout << "Try " << g << " Length : " << length << endl;

	if ( length < minLength ){
	  minLength = length;
	  succGen = g;
	  //	  cout << "New Length : " << minLength << " by " << g << endl;
	  lenReduced = true;
	}
      }

    }

    if (!lenReduced) {
      //      if ( !red_T.shortAndTestTuples( N ) )
      //      else 
      //	break;
      //      cout << "Z : " << shortenBraid(N,z) << endl;
      //      cout << "Z': " << shortenBraid(N,z_conj) << endl;
     

      Word z_diff = shortenBraid(N,z*z_conj);
      cout << "Z DIST : " << z_diff.length() << endl;
      //      if (z_diff.length() > 2)
      //	return false;
      //      else
      //	return true;
      
      red_T.z = z_conj;
      if ( equalUpToCommut( N, bs, red_T, z) )
	cout << "Z COMMUTE WITH WS" << endl;

//       It seems that  z can be hard to recover 
//       because the difference z*z' may commute with both BL and BR
//       Then |BL| + |BR| = |(z*z')^-1 BL (z*z')| + |(z*z')^-1 BR (z*z')|
//       Check it here 
//       int len = 0;
//       int len_conj = 0;
//       for ( int i=0;i<red_T.WL.size();i++){
// 	len += red_T.WL[i].length();
// 	len_conj += shortenBraid(N,-z_diff*red_T.WL[i]*z_diff).length();
//       }
//       for ( int i=0;i<red_T.WR.size();i++){
// 	len += red_T.WR[i].length();
// 	len_conj += shortenBraid(N,-z_diff*red_T.WR[i]*z_diff).length();
//       }     
      
//       if ( len_conj <= len)
// 	cout << "Z COMMUTE WITH WS" << endl;
      
      // copy reduced tuple 
      if ( ret_T )
	*ret_T = TTPTuple( red_T.WL, red_T.WR, z_conj );

      if ( red_T.shortAndTestTuples( N ) )
	return true;
      else
	return false;
      
    }
      
    // UPDATE THE TUPLE
  
    
    // UPDATE LEFT TUPLE
    for ( int iL=0;iL<T.WL.size();iL++){
      Word redWord = shortenBraid(N,-succGen*red_T.WL[iL]*succGen);
      red_T.WL[iL] = redWord;
    }
    // UPDATE RIGTH TUPLE
    for ( int iR=0;iR<T.WR.size();iR++){
      Word redWord = shortenBraid(N,-succGen*red_T.WR[iR]*succGen);
      red_T.WR[iR] = redWord;
    }

    z_conj = z_conj*succGen;   // z is such that w_red = z^-1 w z.
    
    int conj_dist = shortenBraid(N,z*z_conj).length();
    int stop_cond = red_T.shortAndTestTuples( N );
  
    cout << "DIST to Z: " << conj_dist  << " " << shortenBraid(N,z*z_conj)<< " COND : " << stop_cond << endl;
    if ( conj_dist == 0 && stop_cond == 0 ) { // show details since it should not happen
      red_T.shortAndTestTuples( N,true );
    }

    if ( conj_dist == 0 && stop_cond == 1 ) {

      
      //   // If success check if we have the original z by computing 
      bool same_z = true;
      //   if (red_res){
      cout << "Z DIST : " << shortenBraid(N,z*z_conj).length() << endl;
      
      for (int i=0;i<red_T.WL.size();i++){
	if (shortenBraid(N,z_conj*red_T.WL[i]*-z_conj*-T.WL[i]).length() > 0) {
	  same_z = false;
	  break;
	}
      }
      for (int i=0;i<red_T.WR.size();i++){
	if (shortenBraid(N,z_conj*red_T.WR[i]*-z_conj*-T.WR[i]).length() > 0) {
	  same_z = false;
	  break;
	}      
      }
      
      
      if (same_z)
	cout << "SOLUTION CORRECT" << endl;
      else
	cout << "SOLUTION IS NOT CORRECT" << endl;
      
      // copy reduced tuple 
      if ( ret_T )
	*ret_T = TTPTuple( red_T.WL, red_T.WR, z_conj );

      return true;
    }
  }
}

  

///////////////////////////////////////////////////////////////////////////////////////
//
//  TTP ATTACK
//
///////////////////////////////////////////////////////////////////////////////////////

bool TTPAttack::run(const TTPTuple &d) {
  Word original_z = d.z;

  // Convert the tuples
  vector<ThLeftNormalForm> theTuple(d.WL.size() + d.WR.size());
  for (int i = 0; i < d.WL.size(); i++) {
    theTuple[i] = ThLeftNormalForm(BraidGroup(N), d.WL[i]);
    theTuple[i].setPower(theTuple[i].getPower() % 2);
  }
  for (int i = 0; i < d.WR.size(); i++) {
    theTuple[i + d.WL.size()] = ThLeftNormalForm(BraidGroup(N), d.WR[i]);
    theTuple[i + d.WL.size()].setPower(theTuple[i + d.WL.size()].getPower() %
                                       2);
  }

  // execute attack

  // run LBA to restore delta values
  reduceDeltaLBA(theTuple);

  // Test LBA
  bool DelatLBA_succ = true;
  for (int i = 0; i < d.WL.size(); i++) {
    if (shortenBraid(N, theTuple[i].getWord() * -d.WL[i]).length() > 0) {
      DelatLBA_succ = false;
      cout << "F";
    } else {
      cout << "S";
    }
  }

  for (int i = 0; i < d.WR.size(); i++) {
    if (shortenBraid(N, theTuple[i + d.WL.size()].getWord() * -d.WR[i])
            .length() > 0) {
      DelatLBA_succ = false;
      cout << "F";
    } else {
      cout << "S";
    }
  }
  cout << endl;
  if (!DelatLBA_succ)
    cout << "WARNING!!! Delta LBA FAILED!" << endl;

  //  return oneOfSSSReps( d.WL.size(), d.WR.size(), theTuple );
  //  TTPTuple red_T;
  // return simpleLBA(d.WL.size(), d.WR.size(), theTuple, original_z );

  return LBA(d.WL.size(), d.WR.size(), theTuple, original_z);
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

bool TTPAttack::LBA( int NWL, int NWR, const vector<ThLeftNormalForm>& theTuple, const Word& z )
{ 
  
  TTPTuple T;
  T.WL = vector<Word>(NWL);
  T.WR = vector<Word>(NWR);
  
  for ( int j=0;j<theTuple.size();j++)
    if ( j < NWL )
      T.WL[j] = theTuple[j].getWord();
    else
      T.WR[j-NWL] = theTuple[j].getWord();
  
  T.shorten( N );
  
  vector< Word > gens(N-1);
  for (int i=0;i<N-1;i++)
    gens[i] = Word(i+1);

  TTPLBA ttpLBA;
  TTPTuple red_T;
  bool red_res = ttpLBA.reduce(N, BS, T, gens, 600, cout, red_T, z);

  // If success check if we have the original z by computing
  bool same_z = true;
  if (red_res) {
    cout << "SAME Z: " << shortenBraid(N, z * red_T.z).length() << endl;

    // Checking correctness of computations?
    for (int i = 0; i < red_T.WL.size(); i++) {
      if (shortenBraid(N, red_T.z * red_T.WL[i] * -red_T.z * -T.WL[i])
              .length() > 0) {
        same_z = false;
        break;
      }
    }
    for (int i = 0; i < red_T.WR.size(); i++) {
      if (shortenBraid(N, red_T.z * red_T.WR[i] * -red_T.z * -T.WR[i])
              .length() > 0) {
        same_z = false;
        break;
      }
    }
  }

  if (same_z)
    cout << "SOLUTION CORRECT" << endl;
  else
    cout << "SOLUTION IS NOT CORRECT" << endl;
  return red_res;
  
}

bool TTPAttack::simpleLBA( int NWL, int NWR, const vector<ThLeftNormalForm>& theTuple, const Word& z )
{ 
  
  TTPTuple T;
  T.WL = vector<Word>(NWL);
  T.WR = vector<Word>(NWR);
  
  for ( int j=0;j<theTuple.size();j++)
    if ( j < NWL )
      T.WL[j] = theTuple[j].getWord();
    else
      T.WR[j-NWL] = theTuple[j].getWord();
  
  T.shorten( N );
  
  vector< Word > gens(N-1);
  for (int i=0;i<N-1;i++)
    gens[i] = Word(i+1);
  
  TTPLBA ttpLBA;
  
  return ttpLBA.simpleLBA( N,BS,T,z );
}

void TTPAttack::reduceDeltaLBA(vector<ThLeftNormalForm> &theTuple) {
  for (int i = 0; i < theTuple.size(); i++) {
    ThLeftNormalForm nf = theTuple[i];

    // HERE IS LBA
    int len_sav = shortenBraid(N, nf.getWord()).length();
    ThLeftNormalForm nf_sav = nf;
    //	cout << len_sav;

    nf.setPower(nf.getPower() - 2);
    int new_len = shortenBraid(N, nf.getWord()).length();
    while (new_len < len_sav) {

      //		cout << " -> " << new_len;
      len_sav = new_len;
      nf_sav = nf;

      nf.setPower(nf.getPower() - 2);
      new_len = shortenBraid(N, nf.getWord()).length();
    }

    //	cout << "  ?= " << or_len << " : " <<
    //shortenBraid(ttp_conf.N,nf_sav.getWord() * -dw.first[i]) << endl;
    theTuple[i] = nf_sav;
  }
}

// ThLeftNormalForm cycleDecycle(const ThLeftNormalForm& nf ){
    
//   pair< ThLeftNormalForm , ThLeftNormalForm > cycled_nf = nf.cycle();
//   int inf = nf.getPower();
  
//   // @am How many deciclings are neccessary. 
//   // Checck german's paper

//   for (int l=0;l<????;l++)
//     ( cycled_nf.first.getPower() < inf){
//     inf = cycled_nf.first.getPower();
//     cycled_nf = cycled_nf.first.cycle();
//   }		
  
//   int sup = cycled_nf.first.getDecomposition().size() + cycled_nf.first.getPower();
//   cycled_nf = cycled_nf.first.decycle();
//   while ( cycled_nf.first.getDecomposition().size() + cycled_nf.first.getPower() > sup ){
//     sup  = cycled_nf.first.getDecomposition().size() + cycled_nf.first.getPower();
//     cycled_nf = cycled_nf.first.decycle();
//   }		
  
//   return cycled_nf.first;

// }
