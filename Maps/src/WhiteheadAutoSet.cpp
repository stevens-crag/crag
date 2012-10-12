
// Contents: Implementation of the WhiteheadAutoSet class
//
// Principal Authors: Alexei Miasnikov (2003)
//
// Status: 
//
// Revision History:
//


#include "WhiteheadAutoSet.h"
#include "RanlibCPP.h"
#include "errormsgs.h"
//#include "dump.h"

Word WhiteheadAutoSet::reduceBy( const Word& w, const SetOfMaps& theSet )
{
  for (SetOfMaps::const_iterator I=theSet.begin();I!=theSet.end();I++){
    Word new_w = (I->imageOf( w )).cyclicallyReduce();
    if (new_w.length() < w.length())
      return new_w;
  }

  return w;
}


// -------------------------  NielsenAutoSet ---------------------------- //

NielsenAutoSet::NielsenAutoSet( int n ): nGens( n )
{

  
  int ngens = nGens;
  vector<Word> image( nGens );
  for (int i=0;i<ngens;i++)
    image[i] = Word(i+1);

  // for every i sent its generator to one of the mappings for j
  for (int i=0;i<ngens;i++){
    Word gi(i+1);
    for (int j=0;j<ngens;j++){
      if (i!=j){
	Word gj(j+1);

	// add map i->i j
	image[i] = gi * gj;
	theSet.insert( Map( nGens,nGens, image ) );
	image[i] = gi;
	
	// add map i->i j^-1
	image[i] = gi * (gj.inverse());
	theSet.insert( Map( nGens,nGens, image ) );
	image[i] = gi;

	// add map i->j i
	image[i] = gj*gi;
	theSet.insert( Map( nGens,nGens, image ) );
	image[i] = gi;

	// add map i->j^-1 i
	image[i] = (gj.inverse())*gi;
	theSet.insert( Map( nGens,nGens, image ) );
	image[i] = gi;
      }
    }
  }
}

const Map& NielsenAutoSet::getRandomAuto()const
{
  msgs::error("Inefficient implementation");
  int mapi = RandLib::ur.irand(0,theSet.size()-1);
  int i=0;
  for (SetOfMaps::const_iterator I=theSet.begin();I!=theSet.end();I++,i++)
    if ( i==mapi )
      return *I;
}


// -------------------------  RestrictedWhiteheadAutoSet ---------------------------- //

RestrictedWhiteheadAutoSet::RestrictedWhiteheadAutoSet( int n, bool use_conj ): nGens( n )
{
  int ngens = n;

  vector<Word> image( ngens );
  for (int i=0;i<ngens;i++)
    image[i] = Word(i+1);

  // for every i sent its generator to one of the mappings for j
  for (int i=0;i<ngens;i++){
    Word gi(i+1);
    for (int j=0;j<ngens;j++){
      if (i!=j){
	Word gj(j+1);

	// add map i->i j
	image[i] = gi * gj;
	theSet.insert( Map( ngens,ngens, image ) );
	image[i] = gi;
	
	// add map i->j^-1 i
	image[i] = (gj.inverse())*gi;
	theSet.insert( Map( ngens,ngens, image ) );
	image[i] = gi;

	if (use_conj){
	  // add map i->j^-1 i j
	  image[i] = (gj.inverse())*gi*gj;
	  theSet.insert( Map( ngens,ngens, image ) );
	  image[i] = gi;
	}

      }
    }
  }
}

const Map& RestrictedWhiteheadAutoSet::getRandomAuto()const
{
   msgs::error("Inefficient implementation");
 int mapi = RandLib::ur.irand(0,theSet.size()-1);
  int i=0;
  for (SetOfMaps::const_iterator I=theSet.begin();I!=theSet.end();I++,i++)
    if ( i==mapi)
      return *I;
}





// -------------------------  WhiteheadAutoSetType2 ---------------------------- //


WhiteheadAutoSetType2::WhiteheadAutoSetType2(  int n ): nGens( n )
{
  computeSet( n );
}


WhiteheadAutoSetType2::~WhiteheadAutoSetType2( ) 
{
}

void WhiteheadAutoSetType2::computeSet(  int nGens )
{
  int ngens = nGens;
  vector<int> tCounts(ngens,0);

  // for each fixed generator
  for (int g=0;g<ngens;g++){
    
    bool done = false;
    // compute all autos with g - fixed
    while ( !done ){
      
      //      cout << g << endl;
      //      cout << dump<vector<int> >(tCounts) << endl;
      //      cout << getMap( ngens, tCounts, Word( g+1 )) << endl << endl;
      
//       Map theMap1 = getMap( tCounts, Generator( g+1 ) );
//       theSet.insert( theMap1 );
//       Map theMap2 = getMap( tCounts, inv(Generator( g+1 )) );
//       theSet.insert( theMap2 );

      Map m1 = getMap( ngens, tCounts, Word( g+1 ) );
      theSet.insert( m1 );
      Map m2 = getMap( ngens, tCounts, (Word( g+1 )).inverse() );
      theSet.insert( m2 );
   
      int i=ngens-1;
      while ( tCounts[i]+1 == nElemAutos  ){
	if (i==0) done  = true;
	tCounts[i] = 0;
	i--;
      }
      if (!done)
	tCounts[i]++;
      
    }
  }
}

Map WhiteheadAutoSetType2::getMap(   int nGens, const vector<int>& tCounts, Word a )
{
  
  int ngens = nGens;
  vector<Word> image( nGens );

  for (int i=0;i<ngens;i++){
    Word g(i+1);
    if (g != a && g != (a.inverse()) ){
      switch (tCounts[i]){
      case 0:  // x -> x^-1
	image[i] = g;
	break;
      case 1:  // x -> x a
	image[i] = g * a;
	break;
      case 2:  // x -> a^-1 x 
	image[i] = (a.inverse()) * g;
	break;
      case 3:  // x -> a^-1 x a
	image[i] = (a.inverse()) * g * a;
	break;      
      default:
	msgs::error(" WhiteheadAutoSet::getMap(...):Cannot recognize index of an elementary transformation.");
      }
    } else
      image[i] = g;
  }
  
  Map retM(ngens,ngens, image);
  return retM;
}


const Map& WhiteheadAutoSetType2::getRandomAuto()const
{
   msgs::error("Inefficient implementation");
 
  int mapi = RandLib::ur.irand(0,theSet.size()-1);
  int i=0;
  for (SetOfMaps::const_iterator I=theSet.begin();I!=theSet.end();I++,i++)
    if ( i==mapi)
      return *I;
}


bool  WhiteheadMinimization::isMinimal( const Word& w) const
{
  const SetOfMaps& theSet = wSet.getSet();
  for( SetOfMaps::const_iterator I=theSet.begin(); I!=theSet.end();I++ ){
    Word w_trans = I->imageOf( w );
    if (w.length() > (w_trans.freelyReduce()).length())
      return false;
  }
  return true;
}


Word WhiteheadMinimization::findMinimal( const Word& w, ostream* out )const
{
  Word minWord(w);
  const SetOfMaps& theSet = wSet.getSet();
  bool haveSorterWord = false;

  do {
    haveSorterWord = false;
    for( SetOfMaps::const_iterator I=theSet.begin(); I!=theSet.end();I++ ){
      Word w_trans = I->imageOf( minWord ); 
      //   Word w_trans_reduced = w_trans.freelyReduce();
      if (w_trans.length() < minWord.length() ){
	if (out){
	  //const Map& m = *I;
	  *out << "{" << (*I) << "} " << flush; 
	}
	minWord = w_trans;
	haveSorterWord = true;
	break;
      }
    }
  } while ( haveSorterWord );
    
  return minWord;
}
