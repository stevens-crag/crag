
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include "RanlibCPP.h"

#include "braid_group.h"
#include "ThRightNormalForm.h"

#include <sstream>
#include <iterator>
#include <iostream>
#include <fstream>
#include "ShortBraidForm.h"

//#define USE_MPI


#ifdef USE_MPI
#include <mpi.h>
#endif

using namespace std;

#include "LengthAttack.h"
#include "AAGKeyGeneration.h"


//
// RETURNS GENERATORS of Dp obtained from generators of two free factors in genComp
//                    and relators in Prelts. 
//
vector<Word> getSubgroup( int N, int Pgens, 
			  const vector<Word>& Prelts, 
			  const pair< vector<Word>, 
			  vector<Word> >& genComps )
{

  vector<Word> ret( Pgens + Prelts.size() );
  vector<Word> gens1 = genComps.first;
  vector<Word> gens2 = genComps.second;
  for ( int i=0;i<Pgens;i++)
    ret[i]   = shortenBraid( N,gens1[i]*gens2[i] ); 
  
 
  // replace generators of the given relation by 
  // generators of the second factor group
  for ( int i=0;i<Prelts.size();i++){
    Word oldRel = Prelts[i];
    Word newRel ;
    for ( ConstWordIterator I = oldRel.begin(); I!=oldRel.end();I++){
      Word newGen = (*I > 0 ?  gens2[abs(*I)-1] : -gens2[abs(*I)-1]);
      newRel *= newGen;
    }
    ret[i+Pgens] = shortenBraid(  N,newRel );
  }
  
  return ret;
}

//
//
// RETURNS i'th delta obtained from tubes of size c
//
//
Word getDelta( int i,int c )
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
pair< vector<Word>,vector<Word> > getSgGenComponents( int Pgens, int c, int k, int len )
{
  vector<Word> ret1( Pgens );
  vector<Word> ret2( Pgens );
  


  for (int d=0;d<2;d++) { // d = 0 - we generate w's, otherwise u's

    for ( int i=0;i<Pgens;i++ ){
      Word w   = Word::randomWord( k/2 - 1,len );
      Word dw;
      
      for (ConstWordIterator I=w.begin();I!=w.end();I++){
	int delta_i = *I;
	if ( d ) delta_i += (delta_i > 0 ? k/2 : -k/2 ); // @am What if not even?
	Word w_tmp = getDelta( delta_i,c ); 
	dw *= w_tmp;
      }
      
      //Word sdw = shortenBraid(  c*k,dw );
      
      //      cout << w << " : " << w.length() << endl 
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

//
//
//  GENERATES a random word from the subgroup sg in 
//            in generators of Bn
//
//
Word randomSubgroupWord( int N,const vector<Word>& sg )
{						
  Word decomp = Word::randomWord( sg.size(),50 );
  Word rWord;

  for ( ConstWordIterator I=decomp.begin();I!=decomp.end();I++){
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
Word specialSubgroupWord( int N,const Word& a, const Word& t, int n )
{						
  Word tn;
  for (int i=0;i<n;i++)
    tn *= t;

  Word ret = -a*-tn*a*tn*a*-tn*-a*tn;
 
  Word sWord = shortenBraid( N,ret );
  return sWord;
}

//---------------------------------------------------------------------------//
//---------------------------------- main -----------------------------------//
//---------------------------------------------------------------------------//

int main( int argc, char** argv )
{
  // parameters
  int  Pgens  = 2;
  Word a(1);
  Word t(2);

  vector<Word> Prelts(1);
  Prelts[0] = -t*a*t*-a*-a;

  int c = 3;
  int k = 6;

  int D_gen_comp_length = 10;
  int D_gen_comp        = 2*Pgens;  // Pgens  of  w's and the same for u's
  
  //
  //
  //
  //
  //

  for (int i=0;i<100;i++){

    
    
    // CONSTRUCT THE SUBGROUP
    pair< vector<Word>,vector<Word> >  vp = getSgGenComponents( Pgens, c, k, D_gen_comp_length );
    vector<Word> D = getSubgroup( c*k,Pgens,Prelts,vp );
    

    // GENERATE A WORD
    //Word w = randomSubgroupWord( c*k,D );
    Word w =  specialSubgroupWord( c*k, vp.second[0], vp.second[1], 50 );
    cout << w.length() << endl;

  }


  return 0;
}
