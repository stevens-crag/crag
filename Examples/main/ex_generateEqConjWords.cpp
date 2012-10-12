
#include "RanlibCPP.h"
#include "FPGroup.h"
#include "BalancedTree.h"
#include "iostream"
using namespace std;


int main( )
{
  // Create a random reduced word r
  Word r = Word::randomWord( 2 , 5 );
  // Create a vector of words containing 1 word r
  vector< Word > R(1,r);
  // Create a finitely presented group on two generators with relators R
  FPGroup G(2,R);

  // Output a group presentation
  cout << "G = " << G << endl;
  
  // Create a random identity
  Word w = G.randomIdentity_Baltimore( 100 , .1 );
  // Output w
  cout << "w = " << w << endl;
  
  return 0;
}
