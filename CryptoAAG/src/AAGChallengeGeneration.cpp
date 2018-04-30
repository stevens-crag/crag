
#include "braid_group.h"

#include <iterator>
#include <iostream>
#include <fstream>
#include <algorithm>
#include "ShortBraidForm.h"
#include "RanlibCPP.h"
#include "MajorDump.h"

using namespace std;

#include "AAGChallengeGeneration.h"


bool doCommute(int N, const Word& w1, const Word& w2 )
{
  Word w = shortenBraid( N, w1*w2*-w1*-w2 );
  return ( w == Word() );
}
  
//
// RETURNS GENERATORS of Dp obtained from generators of two free factors in genComp
//                    and relators in Prelts. 
//
vector<Word> AAGChallenge::getDpSubgroup( int N, int Pgens, 
					  const vector<Word>& Prelts, 
					  const pair< vector<Word>, vector<Word> >& genComps,
					  int conjLen )
{
  
  vector<Word> ret( Pgens + Prelts.size() );
  vector<Word> gens1 = genComps.first;
  vector<Word> gens2 = genComps.second;
  
  Word rConj   = Word::randomWord( N-1,conjLen );
  
  for ( int i=0;i<Pgens;i++)
    ret[i]   = shortenBraid( N,-rConj*(gens1[i]*gens2[i])*rConj ); 
  
  
  // replace generators of the given relation by 
  // generators of the second factor group
  for ( int i=0;i<Prelts.size();i++){
    Word oldRel = Prelts[i];
    Word newRel ;
    for ( auto I = oldRel.begin(); I!=oldRel.end();I++){
      Word newGen = (*I > 0 ?  gens2[abs(*I)-1] : -gens2[abs(*I)-1]);
	  newRel *= newGen;
    }
    ret[i+Pgens] = shortenBraid(  N, -rConj*newRel*rConj );
  }
  
  return ret;
}

//
//
// RETURNS i'th delta obtained from tubes of size c
//
//
Word  AAGChallenge::getDelta( int i,int c )
{
  int ai = abs(i);
  Word ret;
  for ( int l=0;l<c;l++)
    for( int m=0;m<c;m++)
      ret *= Word( ai*c - l + m );
  
  if (i<0)
    return -ret;
  else
    return ret;
}


//
//
//  CREATES GENERATORS of the free factors
//
//
pair< vector<Word>,vector<Word> >  AAGChallenge::getSgGenComponentsRandom( int Pgens, int c, int k, int len )
{

  *(Dump::dump_out) << "Generate <w>,<u> randomly with length " << len <<  endl;


  vector<Word> ret1( Pgens );
  vector<Word> ret2( Pgens );

  for (int d=0;d<2;d++) { // d = 0 - we generate w's, otherwise u's
    
    for ( int i=0;i<Pgens;i++ ){
      Word w   = Word::randomWord( k/2 - 1,len );
      Word dw;
      
      for (auto I=w.begin();I!=w.end();I++){
	int delta_i = *I;
	if ( d ) delta_i += (delta_i > 0 ? k/2 : -k/2 ); // @am What if not even?
	Word w_tmp = getDelta( delta_i,c ); 
	dw *= w_tmp;
      }
      
      //Word sdw = shortenBraid(  c*k,dw );
      
      //      *(Dump::dump_out) << w << " : " << w.length() << endl 
      //	   << dw << " : " << dw.length() << endl << endl;
      //	   << sdw << " : " << sdw.length() << endl;
      
      // add to the list of generators
      if ( d )
	ret2[i] = dw;
      else
	ret1[i] = dw;
    }
  }
  
  return pair< vector<Word>,vector<Word> >(ret1,ret2);
}

pair< vector<Word>,vector<Word> >  AAGChallenge::getSgGenComponentsSquares( int Pgens, int c, int k )
{

  *(Dump::dump_out) << "Generate sqare Dp generators" << endl;

  vector<Word> ret1( Pgens );
  vector<Word> ret2( Pgens );
  
  for (int d=0;d<2;d++) { // d = 0 - we generate w's, otherwise u's
    
    for ( int i=0;i<Pgens;i++ ){
      
      int delta_i = i+1;
      if ( d ) delta_i += (delta_i > 0 ? k/2 : -k/2 ); // @am What if not even?
      Word w_tmp = getDelta( delta_i,c ); 
      
      // add to the list of generators
      if ( d )
	ret2[i] =  w_tmp*w_tmp;
      else
	ret1[i] =  w_tmp*w_tmp;
    }
  }
  
  return pair< vector<Word>,vector<Word> >(ret1,ret2);
}

//
//
//  GENERATES a random word from the subgroup sg in 
//            in generators of Bn
//
//
Word AAGChallenge::randomSubgroupWord( int N,const vector<Word>& sg )
{						
  Word decomp = Word::randomWord( sg.size(),50 );
  Word rWord;
  
  for ( auto I=decomp.begin();I!=decomp.end();I++){
    Word gen = (*I > 0 ?  sg[abs(*I)-1] : -sg[abs(*I)-1]);
    rWord *= gen; 
  }
  
  Word sWord = shortenBraid( N,rWord );
  return sWord;
}

//
//
//  GENERATES word Wn
//
//
Word AAGChallenge::specialSubgroupWord( const Word& a, const Word& t, int n )
{						
  Word tn;
  for (int i=0;i<n;i++)
    tn *= t;
  
  Word ret = -a*-tn*a*tn*a*-tn*-a*tn;
  
  Word sWord = ret; //shortenBraid( N,ret );
  return sWord;
}

Word AAGChallenge::specialSubgroupWordRandom( const Word& a, const Word& t, int len )
{	
  Word ret;					

  for (int i=0;i<len;i++){
    int n = RandLib::ur.irand(2,5);
    Word tmp_w = specialSubgroupWord( a,t,n );
    ret *= tmp_w;
  }
  
  return ret;
}

//---------------------------------------------------------------------------//
//---------------------------------- generateSubgroup -----------------------//
//---------------------------------------------------------------------------//

vector<Word> AAGChallenge::generateSubgroup( int c, int k, int conj_len  )
{
  // parameters
  int  Pgens  = 2;
  Word a(1);
  Word t(2);
  
  vector<Word> Prelts(1);
  Prelts[0] = -t*a*t*-a*-a;
  
  int D_gen_comp        = 2*Pgens;  // Pgens  of  w's and the same for u's
  
  

  // CONSTRUCT THE SUBGROUP
  bool first_try = true;
  pair< vector<Word>,vector<Word> >  vp;
  do {
  
    if (!first_try)
      *(Dump::dump_out) << "Warn! Gens comute. Try again ..." << endl;
    
    vp = getSgGenComponentsRandom(  Pgens, c, k, 10 );
   
    //pair< vector<Word>,vector<Word> >    vp = getSgGenComponentsSquares( Pgens, c, k );
    first_try = false;
  
  } while ( doCommute(c*k, vp.first[0], vp.first[1] ) || doCommute( c*k, vp.second[0],vp.second[1] ) );


  *(Dump::dump_out) << "F2xF2 generators:" << endl << "Firs : " << endl ;
  copy(vp.first.begin(),vp.first.end(),ostream_iterator<Word>(*(Dump::dump_out)," \n"));
  *(Dump::dump_out) << "Second" << endl;
  copy(vp.second.begin(),vp.second.end(),ostream_iterator<Word>(*(Dump::dump_out)," \n"));
  
    
  vector<Word> D = getDpSubgroup( c*k,Pgens,Prelts,vp,conj_len );
  
  return D;;
}

Word AAGChallenge::generateKeyDecomp(  int len )
{



  // parameters
  int  Pgens  = 2;
  Word a(1);
  Word t(2);
  
  
  // GENERATE A WORD
  Word w =  specialSubgroupWord(  a, t, 10 );
  //  Word w = specialSubgroupWordRandom( a,t, len );
  

  
  return w;
}
