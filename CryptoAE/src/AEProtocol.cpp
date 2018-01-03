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
#include "BraidGroup.h"
#include "ShortBraidForm.h"
#include "Permutation.h"
#include "errormsgs.h"
#include "ThLeftNormalForm.h"
#include "RanlibCPP.h"

//typedef pair< vector< Word > , vector< Word > > TTPTuple;


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

BSets BSets::generateRandom(int N)
{

 BSets retB;
  
  while ( retB.BL.size() <=2 || retB.BR.size() <= 2 ) {
    vector<int> partitionB(N-1,0);
    for (int i=0;i<N-1;i+=2){
      if ( RandLib::ur.rand() >= 0.5 )
	partitionB[i] = -1;
      else
	partitionB[i] =  1;
    }
    
    // fill in the gaps
    for (int i=0;i<N-1;i+=1){
      if ( partitionB[i] == 0 ){
	
	if ( i == 0 ) // the leftmost (should never happen)
	  partitionB[i] =  partitionB[i+1];
	
	else if ( i == N-2) // the rightmost
	  partitionB[i] =  partitionB[i-1];
	
	else if ( partitionB[i-1] ==  partitionB[i+1] )       // somewhere in the middle
	  partitionB[i] =  partitionB[i-1];
	
      }
      
    }
    
    //    copy(partitionB.begin(),partitionB.end(), ostream_iterator<int>(cout," ") );
    
    // Construct BL and BR
    retB.BL.clear();
    retB.BR.clear();
    
    for (int i=0;i<N-1;i+=1){
      if (partitionB[i] > 0 )
	retB.BL.push_back(Word(i+1));
      if (partitionB[i] < 0 )
	retB.BR.push_back(Word(i+1));    
    }
  }


  return retB;

  //cout << endl << "BL : " << BL.size() << " " << "BR : " << BR.size() << endl;
}


//
//
//   TTP TUPLES
//
//
int TTPTuple::length() const 
{
  		int result = 0;

  		for( int i=0 ; i<WL.size() ; ++i )
    		result += WL[i].length( );

  		for( int i=0 ; i<WR.size() ; ++i )
    		result += WR[i].length( );

  		return result;
}
		
		
void TTPTuple::shorten( int N ) 
{
  for( int i=0 ; i<WL.size() ; ++i )
    WL[i] = shortenBraid( N , WL[i] );
  
  for( int i=0 ; i<WR.size() ; ++i )
    WR[i] = shortenBraid( N , WR[i] );
}

bool TTPTuple::testTuples( int N, bool details ) const
{
  if (details){

    cout << "WL:" << endl;
    for (int i=0;i<WL.size();i++)
      cout << WL[i] << endl;

    cout << "WR:" << endl;
    for (int i=0;i<WR.size();i++)
      cout << WR[i] << endl;

  }
  
  vector<char> gensFlags(N-1,0);
  for ( int i=0;i<WL.size();i++){
    for (ConstWordIterator I=WL[i].begin();I!=WL[i].end();I++){
      gensFlags[abs(*I)-1] = 1;
    }
  }
  
  for ( int i=0;i<WR.size();i++){
    for (ConstWordIterator I=WR[i].begin();I!=WR[i].end();I++){
      int gI = abs(*I)-1;
      int glI = gI > 0 ? gI-1 : gI;
      int grI = gI < N-2 ? gI+1 : gI;
      if ( gensFlags[gI] || gensFlags[glI] || gensFlags[grI] ){
	if (details)
	  cout << gI+1 << " FAILS" << endl;
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
		for(ConstWordIterator I=w.begin();I!=w.end();I++){
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
	for (ConstWordIterator I=wsAsOneWord.begin();I!=wsAsOneWord;I++){
		int wi = abs(*I)-1;
		// @am not sure what to do with inverse. Pretend that all positive for now
		Permutation p(the_n-1);  //@am check if n correct!!!
		p.change(wi,wi+1);
		BurauGenerator w_bg(wi,p);
		pubKey = starMult(pubKey,w_bg);
	}	
	
	return pubKey;
}
		

TTPTuple AEKeyExchange::generateTuples( const TTP_Conf& ttp_conf, const BSets& bs )
{
	int nBL = ttp_conf.nBL;
	int nBR = ttp_conf.nBR;
	int N = ttp_conf.N;
	int n = N-1;
	int nGamma = ttp_conf.nGamma;
	
	int len_z = ttp_conf.len_z;
	int len_w = ttp_conf.len_w;


	// Choose the secret conjugator
	Word z = Word::randomWord(n,len_z);
	
	//choose the tuples

	vector<Word> wL(nGamma);
	vector<Word> wR(nGamma);
	vector<Word> origWL(nGamma);
	vector<Word> origWR(nGamma);

	for (int i=0;i<nGamma;i++){
	  origWL[i] = (Word::randomWord(bs.BL.size(),len_w)).replaceGenerators(bs.BL);
	  origWR[i] = (Word::randomWord(bs.BR.size(),len_w)).replaceGenerators(bs.BR);
	  //cout << "BL : " << origWL[i] << endl;
	  //cout << "BR : " << origWR[i] << endl;
 
	  wL[i] = -z*origWL[i]*z;
	  wR[i] = -z*origWR[i]*z;
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
	
	TTPTuple ret(wL,wR);
	ret.z = z;
	//	ret.origWL = origWL;
	//	ret.origWR = origWR;
	return ret;
}
