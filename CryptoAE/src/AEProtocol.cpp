// Copyright (C) 2007 Alex Myasnikov
// Contents: Implementation of TTP  Attack
//
// Principal Authors: Alex Myasnikov
//
// Revision History:
//


#include <set>
#include <list>

#include "AEProtocol.h"
#include "braid_group.h"
#include "ShortBraidForm.h"
#include "Permutation.h"
#include "errormsgs.h"
#include "ThLeftNormalForm.h"
#include "RanlibCPP.h"
#include <thread>
#include "LinkedBraidStructure.h"

//typedef pair< vector< Word > , vector< Word > > TTPTuple;


static int abelinization(const Word& w) {
  int result = 0;
  for (const auto &g : w) {
    result += g > 0 ? 1 : -1;
  }
  return result;
}

//
//
//  BSets
//
//

BSets BSets::generateEqual(int N)
{
  // Naive implementation just takes first nBL  generators
  // to be BL and last nBR generators to be BR
  if ( N % 2 != 0 )
    msgs::error("TTPAttack::generateTuples(): N should be even.");
  
  BSets retB;
  int n = (N-2)/2;

  retB.BL = vector<Word>(n);
  retB.BR = vector<Word>(n);
  
  for ( int i=0;i<n;i++)
    retB.BL[i] = Word(i+1);
  
  for ( int i=0;i<n;i++)
    retB.BR[i] = Word(n+1+i+1);

  return retB;
}

BSets BSets::generateRandom(int N) {

  BSets retB;

  while (retB.BL.size() <= 2 || retB.BR.size() <= 2) {
    vector<int> partitionB(N - 1, 0);
    for (int i = 0; i < N - 1; i += 2) {
      if (RandLib::ur.rand() >= 0.5)
        partitionB[i] = -1;
      else
        partitionB[i] = 1;
    }

    // fill in the gaps
    for (int i = 0; i < N - 1; i += 1) {
      if (partitionB[i] == 0) {

        if (i == 0) // the leftmost (should never happen)
          partitionB[i] = partitionB[i + 1];

        else if (i == N - 2) // the rightmost
          partitionB[i] = partitionB[i - 1];

        else if (partitionB[i - 1] ==
                 partitionB[i + 1]) // somewhere in the middle
          partitionB[i] = partitionB[i - 1];
      }
    }

    //    copy(partitionB.begin(),partitionB.end(), ostream_iterator<int>(cout,"
    //    ") );

    // Construct BL and BR
    retB.BL.clear();
    retB.BR.clear();

    for (int i = 0; i < N - 1; i += 1) {
      if (partitionB[i] > 0)
        retB.BL.push_back(Word(i + 1));
      if (partitionB[i] < 0)
        retB.BR.push_back(Word(i + 1));
    }
  }

  return retB;

  // cout << endl << "BL : " << BL.size() << " " << "BR : " << BR.size() <<
  // endl;
}

//
//
//   TTP TUPLES
//
//
int TTPTuple::length() const {
  int result = 0;
  for (const auto &w : WL) {
    result += w.length();
  }
  for (const auto &w : WR) {
    result += w.length();
  }
  return result;
}

void TTPTuple::shorten_parallel(int N) {
  vector<std::thread> threads;
  threads.reserve(WL.size() + WR.size());

  for (int i = 0; i < WL.size(); ++i) {
    threads.emplace_back([this, N, i]() {
      WL[i] = shortenBraid(N, WL[i]);
    });
  }

  for (int i = 0; i < WR.size(); ++i) {
    threads.emplace_back([this, N, i]() {
      WR[i] = shortenBraid(N, WR[i]);
    });
  }

  for (auto& th : threads) {
    th.join();
  }
}

void TTPTuple::shorten(int N) {
  for (auto &w : WL) {
    w = shortenBraid(N, w);
  }
  for (auto &w : WR) {
    w = shortenBraid(N, w);
  }
}

bool TTPTuple::equivalent(int N, const TTPTuple &t) const {
  typedef ThLeftNormalForm NF;
  crag::braidgroup::BraidGroup B(N);

  vector<NF> L1, R1, L2, R2;
  for (int i = 0; i < WL.size(); ++i) {
    NF nf(B, z * WL[i] * -z);
    nf.setPower(nf.getPower() - 2 * deltaSQL[i]);
    L1.push_back(nf);
  }
  for (int i = 0; i < WR.size(); ++i) {
    NF nf(B, z * WR[i] * -z);
    nf.setPower(nf.getPower() - 2 * deltaSQR[i]);
    R1.push_back(nf);
  }
  for (int i = 0; i < t.WL.size(); ++i) {
    NF nf(B, t.z * t.WL[i] * -t.z);
    nf.setPower(nf.getPower() - 2 * t.deltaSQL[i]);
    L2.push_back(nf);
  }
  for (int i = 0; i < t.WR.size(); ++i) {
    NF nf(B, t.z * t.WR[i] * -t.z);
    nf.setPower(nf.getPower() - 2 * t.deltaSQR[i]);
    R2.push_back(nf);
  }

  for (int i = 0; i < WL.size(); ++i) {
    if (L1[i] != L2[i])
      return false;
  }
  for (int i = 0; i < WR.size(); ++i) {
    if (R1[i] != R2[i])
      return false;
  }

  return true; 
}

TTPTuple TTPTuple::takeModuloDeltaSQ(int N) const {
  typedef ThLeftNormalForm NF;
  crag::braidgroup::BraidGroup B(N);

  TTPTuple result = *this;
  for (int i = 0; i < result.WL.size(); ++i) {
    NF nf(B, result.WL[i]);
    const auto p = nf.getPower();
    result.deltaSQL[i] += ((nf.getPower() % 2) - p) / 2;
    // cout << nf.getPower() << ", " << nf.getPower() % 2 << ", " << result.deltaSQL[i] << endl;
    nf.setPower(nf.getPower() % 2);
    result.WL[i] = nf.getReducedWord2();
  }

  for (int i = 0; i < result.WR.size(); ++i) {
    NF nf(B, result.WR[i]);
    const auto p = nf.getPower();
    result.deltaSQR[i] += ((nf.getPower() % 2) - p) / 2;
    nf.setPower(nf.getPower() % 2);
    result.WR[i] = nf.getReducedWord2();
  }

  return result;
}

static int anticipated_delta(int N, Word &w) {
  const auto l = abelinization(w);
  const auto step = N * (N - 1);
  const auto to_achieve1 = (l % step + step) % step;
  const auto to_achieve2 = to_achieve1 - step;
  return abs(to_achieve1) < abs(to_achieve2) ? -(l - to_achieve1) / step : -(l - to_achieve2) / step;
}

//! Find a power Delta^2p s.t. |Delta^2p*nf| is minimal in the <Delta^2>-coset
static int multiplyByDeltaSQtoReduceLength(int N, Word &w, const int delta, bool power_reset) {
  typedef ThLeftNormalForm NF;
  crag::braidgroup::BraidGroup B(N);

  int best_nf;
  map<int, Word> weights;
  // initial data
  if (power_reset) {
    NF nf(B, w);
    // -nf.getDecomposition().size()/4 is the anticipated value
    // best_nf = -nf.getDecomposition().size() / 4;
    best_nf = anticipated_delta(N, w);
    ThLeftNormalForm nf2 = nf;
    nf2.setPower(nf.getPower() + 2 * best_nf);
    weights[best_nf] = shortenBraid(N, nf2.getReducedWord2());
    w = weights[best_nf];
    cout << weights[best_nf].length() << ", ";
  }
  else {
    best_nf = 0;
    weights[best_nf] = w;
  }

  const Permutation omega = Permutation::getHalfTwistPermutation(N);
  Word omegaWord = Word(omega.geodesicWord());
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

TTPTuple TTPTuple::multiplyElementsByDeltaSQtoReduceLength(int N, const int delta, bool power_reset) const {
  vector<std::thread> threads;
  threads.reserve(WL.size() + WR.size());

  auto result = *this;

  for (int i = 0; i < WL.size(); ++i) {
    threads.emplace_back([&result, N, i, delta, power_reset, this]() {
      result.deltaSQL[i] += multiplyByDeltaSQtoReduceLength(N, result.WL[i], delta, power_reset);
    });
  }

  for (int i = 0; i < WR.size(); ++i) {
    threads.emplace_back([&result, N, i, delta, power_reset, this]() {
      result.deltaSQR[i] += multiplyByDeltaSQtoReduceLength(N, result.WR[i], delta, power_reset);
    });
  }

  for (auto& th : threads) {
    th.join();
  }

  return result;
}

TTPTuple TTPTuple::conjugate(int N, const Word &b) const {
  TTPTuple result = *this;
  for (auto &w : result.WL)
    w = -b * w * b;
  for (auto &w : result.WR)
    w = -b * w * b;
  result.z *= b;
  return result;
}

void TTPTuple::printPowers() const {
  cout << "|";
  for (int i = 0; i < deltaSQL.size(); i++) {
    if (deltaSQL[i] != 0) {
      cout << "F" << "(" << deltaSQL[i] << ")";
    } else {
      cout << "S";
    }
  }
  cout << "|";
  for (int i = 0; i < deltaSQR.size(); i++) {
    if (deltaSQR[i] != 0) {
      cout << "F" << "(" << deltaSQR[i] << ")";
    } else {
      cout << "S";
    }
  }
  cout << "|";
}


bool TTPTuple::testTuples2(int N, bool details) const {
  auto t1 = *this;
  // 1. Dehornloy form configuration #1
  for (auto &w : t1.WL) {
    LinkedBraidStructure df(N - 1, w);
    df.removeRightHandles();
    w = df.translateIntoWord();
  }
  for (auto &w : t1.WR) {
    LinkedBraidStructure df(N - 1, w);
    df.removeLeftHandles();
    w = df.translateIntoWord();
  }
  if (t1.testTuples(N, details)) {
    return true;
  }

  // 2. Dehornloy form configuration #2
  auto t2 = *this;
  for (auto &w : t2.WL) {
    LinkedBraidStructure df(N - 1, w);
    df.removeLeftHandles();
    w = df.translateIntoWord();
  }
  for (auto &w : t2.WR) {
    LinkedBraidStructure df(N - 1, w);
    df.removeRightHandles();
    w = df.translateIntoWord();
  }
  if (t2.testTuples(N, details)) {
    return true;
  }

  return testTuples(N, details);
}


bool TTPTuple::testTuples(int N, bool details) const {
  if (details) {
    cout << "WL:" << endl;
    for (int i = 0; i < WL.size(); i++) {
      cout << WL[i] << endl;
    }
    cout << "WR:" << endl;
    for (int i = 0; i < WR.size(); i++) {
      cout << WR[i] << endl;
    }
  }

  vector<char> gensFlags(N - 1, 0);
  for (int i = 0; i < WL.size(); i++) {
    for (auto I = WL[i].begin(); I != WL[i].end(); I++) {
      gensFlags[abs(*I) - 1] = 1;
    }
  }

  for (int i = 0; i < WR.size(); i++) {
    for (auto I = WR[i].begin(); I != WR[i].end(); I++) {
      int gI = abs(*I) - 1;
      int glI = gI > 0 ? gI - 1 : gI;
      int grI = gI < N - 2 ? gI + 1 : gI;
      if (gensFlags[gI] || gensFlags[glI] || gensFlags[grI]) {
        if (details)
          cout << gI + 1 << " FAILS" << endl;
        return false;
      }
    }
  }
  return true;
}




//
//
//   MatrixFp
//
//

MatrixFp::MatrixFp( const MatrixFp& m  )
{
	if ( the_n != m.the_n )
		msgs::error("MatrixFp( MatrixFp) : Matrices are of different size.");
		
	for ( int i=0;i<the_n;i++ )
		for(int j=0;j<the_n;j++)
			theMatrix[i][j] = m.theMatrix[i][j];
			
}

MatrixFp& MatrixFp::operator = ( const MatrixFp& m )
{
	
	if ( the_n != m.the_n )
		msgs::error("MatrixFp( MatrixFp) : Matrices are of different size.");
		
	for ( int i=0;i<the_n;i++ )
		for(int j=0;j<the_n;j++)
			theMatrix[i][j] = m.theMatrix[i][j];
  return *this;
}

void MatrixFp::init()
{
	theMatrix = vector< vector<int> >(the_n,vector<int>(the_n,0));
//	new int*[the_n];
//	for (int i=0;i<the_n;i++){
//		theMatrix[i] = new int[the_n];	
}

void MatrixFp::clean()
{
//	for (int i=0;i<the_n;i++)
//		delete [] theMatrix[i];
//	delete [] theMatrix;
}

MatrixFp MatrixFp::random( int n, int p )
{
	MatrixFp m(n,p);
	for (int i=0;i<n;i++)
		for (int j=0;j<n;j++)
			m.theMatrix[i][j] = RandLib::ur.irand(0,p-1); // @am Change to satsify all conditions !!
	return m;
}

MatrixFp MatrixFp::ID( int n, int p )
{
	MatrixFp m(n,p);
	for (int i=0;i<n;i++)
			m.theMatrix[i][i] = 1;
	return m;
}


MatrixFp MatrixFp::operator * ( const MatrixFp& m ) const
{
	MatrixFp res(the_n,the_p);
	
	if ( the_n != m.the_n )
		msgs::error("MatrixFp * MatrixFp : Matrices are of different size.");

	for (int i=0;i<the_n;i++)
		for (int j=0;j<the_n;j++){
			int  theSum = 0;
			for (int k=0;k<the_n;k++)
				theSum = (theSum + theMatrix[i][k]*m.theMatrix[k][j]) % the_p;
			res.theMatrix[i][j] = theSum;
		}
	return res;
}


MatrixFp MatrixFp::operator + ( const MatrixFp& m ) const
{
	if ( the_n != m.the_n )
		msgs::error("MatrixFp * MatrixFp : Matrices are of different size.");

	MatrixFp res(the_n,the_p);

	for (int i=0;i<the_n;i++)
		for (int j=0;j<the_n;j++){
			res.theMatrix[i][j] = (theMatrix[i][j] + m.theMatrix[i][j]) % the_p;
		}
	return res;
}

MatrixFp MatrixFp::scalar_mult( int l ) const
{

	MatrixFp res(the_n,the_p);

	for (int i=0;i<the_n;i++)
		for (int j=0;j<the_n;j++){
			res.theMatrix[i][j] = (theMatrix[i][j]*l) % the_p;
		}
	return res;
}

MatrixFp MatrixFp::getPower( int e )const
{
	if (e<1)
		msgs::error("MatrixFp::gegPower() : Power has to be positive");
	MatrixFp res( *this );
	for (int i=0;i<e-1;i++)
		res = res*(*this);
		
	return res;
}

//
//
//  AEKeyExchange
//
//

ProdElement AEKeyExchange::starMult( const ProdElement& pe, const BurauGenerator& bg )
{
	// fix elements t_i in Fp to define homomorphism PI
	vector<int>  T(the_n);
	for (int i=0;i<the_n;i++)
		T[i] = RandLib::ur.irand(0,the_p-1); //@am check if tis is correct
		
	// Compute the star operation
	MatrixFp resM(the_n,the_p);
	Permutation resP = pe.second*bg.second;
	
	MatrixFp evalMatrix(MatrixFp::ID(the_n,the_p));
	// NOT SURE ABOUT INDICES
	int t_val = T[pe.second[bg.first]]; // this together with next line is equal to \pi(s^x_i(t))
	evalMatrix.set(bg.first,bg.first,t_val);
	if( bg.first > 0)
		evalMatrix.set(bg.first,bg.first-1,t_val);
	
	resM = pe.first*evalMatrix;
	
	return ProdElement(resM,resP);
}

ProdElement AEKeyExchange::generatePublicKey( const vector<Word>& v, int r, int m, int wl )
{
	// generate matrix n
	MatrixFp N(the_n,the_p);
	for (int i=0;i<r;i++){
		int e = 1; //@am how we choose the power?
		MatrixFp pM0 = M0.getPower( e );
		int l = RandLib::ur.irand(0,the_p-1); //@am can we allow 0???
		N = N + pM0.scalar_mult(l);
	}
	
	// generate random words  in the subgroup
	vector<Word> ws(m);
	Word wsAsOneWord;
	for (int i=0;i<m;i++){
		Word w = Word::randomWord(v.size(),wl);
		Word image_w;
		for(auto I=w.begin();I!=w.end();I++){
			int w_i = abs(*I)-1;
			if (*I > 0)
				image_w *= v[w_i];
			else
				image_w *= -v[w_i];
		}
		ws[i] = image_w;
		wsAsOneWord *= image_w;
	}
	
		
	// generate the public key
	
	ProdElement pubKey(N,Permutation());
	for (auto I=wsAsOneWord.begin();I!=wsAsOneWord.end();I++){
		int wi = abs(*I)-1;
		// @am not sure what to do with inverse. Pretend that all positive for now
		Permutation p(the_n-1);  //@am check if n correct!!!
		p.change(wi,wi+1);
		BurauGenerator w_bg(wi,p);
		pubKey = starMult(pubKey,w_bg);
	}	
	
	return pubKey;
}

TTPTuple AEKeyExchange::generateTuples(const TTP_Conf &ttp_conf, const BSets &bs) {
  int nBL = ttp_conf.nBL;
  int nBR = ttp_conf.nBR;
  int N = ttp_conf.N;
  int n = N - 1;
  int nGamma = ttp_conf.nGamma;

  int len_z = ttp_conf.len_z;
  int len_w = ttp_conf.len_w;

  // Choose the secret conjugator
  Word z = Word::randomWord(n, len_z);

  // choose the tuples

  vector<Word> wL(nGamma);
  vector<Word> wR(nGamma);
  vector<Word> origWL(nGamma);
  vector<Word> origWR(nGamma);

  for (int i = 0; i < nGamma; i++) {
    origWL[i] =
        (Word::randomWord(bs.BL.size(), len_w)).replaceGenerators(bs.BL);
    origWR[i] =
        (Word::randomWord(bs.BR.size(), len_w)).replaceGenerators(bs.BR);
    // cout << "BL : " << origWL[i] << endl;
    // cout << "BR : " << origWR[i] << endl;

    wL[i] = -z * origWL[i] * z;
    wR[i] = -z * origWR[i] * z;
  }

  // 	cout << "z : " << z << endl << endl;

  // 	cout << "BL: " << endl;
  // 	for (int i=0;i<wL.size();i++)
  // 		cout << wL[i] << " ";
  // 	cout << endl << endl;

  // 	cout << "BR: " << endl;
  // 	for (int i=0;i<wR.size();i++)
  // 		cout << wR[i] << " ";
  // 	cout << endl;

  TTPTuple ret(wL, wR);
  ret.origZ = z;
  //	ret.origWL = origWL;
  //	ret.origWR = origWR;
  return ret;
}
