// Copyright (C) 2007 Alexey Myasnikov
// Contents: Definition of classes for an attack on TTP algorithm 
//
//  This is an implementation of an attack on
//  TTP algorithm for generating public sets of generators of the Algebraic Eraser protocol described in 
//  I. Anshel, M. Anshel,  D. Goldfeld, S. Lemieux, "Key Agreement, the Algebraic Eraser, and Lightweight Cryptography", 
//  Algebraic Methods in Cryptography, CONM Vol 41 (2006), AMS, pp. 1-38

//
// Principal Authors: Alexey Myasnikov
//
// Revision History:
//

#ifndef _TTPAttack_H_
#define _TTPAttack_H_

#include "Word.h"
#include "AEProtocol.h"
#include "ThLeftNormalForm.h"
#include <vector>

using namespace std;

class TTPTuple;
typedef pair< int , TTPTuple > NODE;


//
//
//  THE LENGTH_BASED ATTACK TO REDUCE ELEMENTS AND RECOVER THE COMMON CONJUGATOR 
//
//

class TTPLBA
{
public:
	TTPLBA(){}	

	bool  reduce( int N , const BSets& bs, const TTPTuple& theTuple , const vector< Word >& gens, int sec, ostream& out, TTPTuple& red_T, const Word& z );
	bool simpleLBA( int N , const BSets& bs, const TTPTuple& theTuple, const Word& z, TTPTuple* ret_T = NULL );
	
private:												
 
	void addNewElt( const TTPTuple& T , const set< NODE >& checkedElements , set< NODE >& uncheckedElements );
	void tryNode( int N , NODE cur , const vector<Word>& gens , const set<NODE>& checkedElements ,  set<NODE>& uncheckedElements, int& min_weight );
 	
 	TTPTuple savTuple;
};


//
//
//  THE ATTACK ON TTP ALGORITHM 
//
//

//! This is an implementation of an attack on TTP algorithm for generating public sets of generators of the Algebraic Eraser protocol
class TTPAttack {
public:

//! Constructor 
/*! \param n  - group rank (number of braids
\param bs -the initial commuting sets of subgroup generators 
*/
  TTPAttack( int n, BSets bs) : N( n ), BS( bs )  {}

enum TTP_Result { 
  TTP_FAILED, 
//  TIME_EXPIRED,
  TTP_NOTSURE,
  TTP_SUCCESSFULL
};
 

//! Excutes te attack
/*!
 \param d - the output of TTP algorithm
 */
 bool run( const TTPTuple& d );
 
 
 
 
 private:

 inline void printStats(const ThLeftNormalForm& nf, ostream& out ){
   out << "inf : " << nf.getPower() << " sup : " << nf.getDecomposition().size() << flush;
 }
 
 
 ThLeftNormalForm cycleDecycle(const ThLeftNormalForm& nf );

 bool simpleLBA( int NWL, int NWR, const vector<ThLeftNormalForm>& theTuple, const Word& z );
 bool LBA( int NWL, int NWR, const vector<ThLeftNormalForm>& theTuple, const Word& z );
 bool oneOfSSSReps( int N1, int N2, const vector<ThLeftNormalForm>& theTuple );
 void reduceDeltaLBA( vector<ThLeftNormalForm>& theTuple );


 int N;
 BSets BS;
};

#endif
