// Copyright (C) 2007 Alexey Myasnikov
// Contents: Definition of classes for an attack on TTP algorithm 
//
//  This is an implementation of some element of AAGL2 (Algebraic Eraser) protocol described in 
//  I. Anshel, M. Anshel,  D. Goldfeld, S. Lemieux, "Key Agreement, the Algebraic Eraser, and Lightweight Cryptography", 
//  Algebraic Methods in Cryptography, CONM Vol 41 (2006), AMS, pp. 1-38

//
// Principal Authors: Alexey Myasnikov
//
// Revision History:
//

#ifndef _AEPROTOCOL_H_
#define _AEPROTOCOL_H_

#include "Word.h"
#include "ThLeftNormalForm.h"
#include <vector>
#include <utility>

using namespace std;

class TTPTuple;

///////////////////////////////////////////////////////////////////////////////
//
//    TTP Algorithm 
//
//////////////////////////////////////////////////////////////////////////////



//
//
//  CONFIGURATION
//
//

//! Set of parameters required to construct  the protocol instance
struct TTP_Conf{

 	int nBL;
	int nBR;
	int N;
	int nGamma;
	
	int len_z;
	int len_w;
		
  friend ostream& operator << (ostream& out, const TTP_Conf& ttp_conf){
	
	int len_z = 20;
	int len_w = 20;
    out << "        #BL :" <<  ttp_conf.nBL << endl;
    out << "        #BR :" <<  ttp_conf.nBR << endl;
    out << "          N :" <<  ttp_conf.N << endl;
    out << "      Gamma :" <<  ttp_conf.nGamma << endl;
    out << "      len z :" <<  ttp_conf.len_z << endl;
    out << "      len w :" <<  ttp_conf.len_w << endl;
		
    return out;
  }
  
};


//
//
// Generator sets
//
//

//! Implements construction of the initial commuting sets of subgroup generators
class BSets
{
 public:
  vector<Word> BL;
  vector<Word> BR;
  
  //! Generates two sets of commuting generators by splitting the original N generators in 2, generators are taken randomly
  static BSets generateRandom(int N);
  //! Generates two sets of commuting generators by separating first (N-2)/2 generators  from (N-2)/2 last ones
  static BSets generateEqual(int N);

  //! Prints B sets
  friend ostream& operator << ( ostream& out, const BSets& bs ) {
    out << "BL{ ";
    for (int i=0;i<bs.BL.size();i++)
      out << bs.BL[i] << " , ";
    out << "} BR{ ";
    for (int i=0;i<bs.BR.size();i++)
      out << bs.BR[i] << " , ";
    out << "}";    

    return out;
  }
  
};


//
//
//  TTP TUPLE 
//
//
//! Implements tuples corresponding to the putput of TTP algorithm
class TTPTuple {
public:
  TTPTuple() {}
  TTPTuple(const vector<Word> &L, const vector<Word> &R)
      : WL(L), WR(R), deltaSQL(L.size(), 0), deltaSQR(R.size(), 0) {}

  TTPTuple(const vector<Word> &L, const vector<Word> &R, const Word &conj)
      : WL(L), WR(R), deltaSQL(L.size(), 0), deltaSQR(R.size(), 0), z(conj) {}

  vector<Word> WL;
  vector<Word> WR;

  //! Check if auxiliary data is correct
  bool equivalent(int N, const TTPTuple &t) const;

  //! Trying to reconstruct the original Delta powers. Use power_reset = true if you want to reset the power of Delta to an "anticipated value".
  TTPTuple multiplyElementsByDeltaSQtoReduceLength(int N, const int delta = 3, bool power_reset = true) const;

  //! Take modulo DeltaSQ. Normal forms (obfuscation) applies.
  TTPTuple takeModuloDeltaSQ(int N) const;
  //! The sum of lengths of words in the tuples
  int length() const;
  //! Apply shortenBraid to each word in the tuples
  void shorten(int N);
  void shorten_parallel(int N);
  //! Test tuples for being "seprated", i.e. two nonintersecting sets of
  //! commuting generators
  /*!
   \param N - braid group rank
   \param deails - if true will print verbose information
   */
  TTPTuple conjugate(int N, const Word& b) const;
  
  //! Test if elements in tuples are separated
  bool testTuples(int N, bool details) const;
  bool testTuples2(int N, bool details) const;


  //! Performs shorten and then test on being "separated"
  bool shortAndTestTuples(int N, bool details = false) {
    shorten_parallel(N);
    // shorten(N);
    return testTuples(N, details);
  }

  //! Tuple ordering operator
  friend bool operator<(const TTPTuple &t1, const TTPTuple &t2) {
    return t1.WL < t2.WL && t1.WR < t2.WR;
  }

  void printPowers() const;

  //! Auxiliary member used in TTP-attack: conjugator used to get this tuple
  Word z;
  //! Auxiliary member used in TTP-attack: powers of Delta^2 used to get this tuple
  vector<int> deltaSQL;
  vector<int> deltaSQR;

  Word origZ;

private:
  vector<Word> origWL;
  vector<Word> origWR;
};

///////////////////////////////////////////////////////////////////////////////
//
//    AE PROTOCOL
//
//////////////////////////////////////////////////////////////////////////////


//
//
//  MatrixFp: matrix over a finite field
//
//

class MatrixFp
{
	public:
		MatrixFp( int n, int p ): the_n( n ), the_p( p ) { init(); }
		~MatrixFp() { clean(); }
		MatrixFp( const MatrixFp& m );
		MatrixFp& operator = (const MatrixFp& m);
		
	  	inline MatrixFp operator + ( const MatrixFp& w ) const;
	  	inline MatrixFp operator * ( const MatrixFp& w ) const;
		inline MatrixFp scalar_mult( int l ) const;  
		MatrixFp getPower( int e )const;
		
		static MatrixFp random( int n, int p );
		static MatrixFp ID( int n, int p );		
		
		void set(int i,int j, int v) { 
			theMatrix[i][j] = v;
		} 
	private:
	
	// METHODS
		void init();
		void clean();
	// DATA
		int the_n;
		int the_p;
//		int** theMatrix;
		vector< vector<int> > theMatrix;		
};

typedef pair<MatrixFp,Permutation> ProdElement;
typedef pair<int,Permutation> BurauGenerator;


class AEKeyExchange 
{
	public:
		AEKeyExchange( int n,int p, const TTPTuple& ttpt ):
			the_n( n ), 
			the_p( p ), 
			theTTPTuple( ttpt ),
			M0( MatrixFp::random(n,p) ) 
			{
			} 
			
		ProdElement alicePublicKey() {
			int r = 10;
			int m = 10;
			int wl = 10;
			generatePublicKey( theTTPTuple.WL , r , m , wl );
		}
		ProdElement bobPublicKey() {
			int r = 10;
			int m = 10;
			int wl = 10;
			generatePublicKey( theTTPTuple.WR , r , m , wl );
		}


		static TTPTuple generateTuples( const TTP_Conf& ttp_conf , const BSets& bs);

	private:
	
		ProdElement starMult( const ProdElement& pe, const BurauGenerator& bg );
		ProdElement generatePublicKey( const vector<Word>& v,  int r, int m, int wl );
		
		int the_n;
		int the_p;
		TTPTuple theTTPTuple;
		MatrixFp M0;
};
#endif
