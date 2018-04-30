
#include "time.h"
#include "stdlib.h"
#include "braid_group.h"
#include "ThRightNormalForm.h"


void test_centr( )
{
  int rank = 60;
  crag::braidgroup::BraidGroup B( rank );
  typedef ThRightNormalForm NF;
  
  int exp = 1;
  for( int e=0 ; e<exp ; ++e ) {
    Word w = Word::randomWord( rank-1 , 2*rank );
    cout << "w = " << w << endl;
    NF nf = NF( B , w );
    // set<ThRightNormalForm> centr = nf.computeNormalizer( );
    set<ThRightNormalForm> centr = nf.computeCentralizer( );
    
    
  }
}

void test_NF_computation() {
  int rank = 4;
  crag::braidgroup::BraidGroup B(rank);
  typedef ThRightNormalForm NF;

  int exp = 10;
  for (int e = 0; e < exp; ++e) {
    Word w = Word::randomWord(rank - 1, 6);
    // vector<int> gens = { 1, 2, 1 };
    // Word w(gens);
    cout << "w = " << w << endl;

    // compute a normal form
    NF nf = NF(B, w);
    cout << "NF = " << endl << nf << endl;

    // compute a normal form of the inverse
    NF nf_inv = -nf;
    cout << "NF_inv = " << endl << nf_inv << endl;

    // multiply normal forms
    NF nf_mult = nf * nf_inv;
    if (!nf_mult.isTrivial()) {
      cout << "Unexpected result A!!!" << endl;
      exit(1);
    }

    Word w1 = nf_inv.getWord();
    NF nf2(B, w * w1);
    if (!nf2.isTrivial()) {
      cout << "Unexpected result B!!!" << endl;
      exit(1);
    }
    cout << "=============================" << endl;
  }
}

//---------------------------------------------------------------------------//
//-------------------------------------- main -------------------------------//
//---------------------------------------------------------------------------//

int main() {
  // srand(time(0));
  // test_NF_computation();
  test_centr();

  return 0;
}
