
#include "ThLeftNormalForm.h"
#include "ThRightNormalForm.h"
#include "time.h"
#include "stdlib.h"
#include "braid_group.h"

#include "ShortBraidForm.h"

#include "PermutationEnumerator.h"


//---------------------------------------------------------------------------//
//-------------------------------------- main -------------------------------//
//---------------------------------------------------------------------------//


void test_conjugacy_uss( )
{
  int rank = 30;
  crag::braidgroup::BraidGroup B( rank );
  typedef ThLeftNormalForm NF;
  // typedef ThRightNormalForm NF;
  
  int exp = 100;
  for( int e=0 ; e<exp ; ++e ) {

    cout << "exp N " << e << endl;
    Word w = Word::randomWord( rank-1 , 500 );
    Word c = Word::randomWord( rank-1 , 500 );
    
    // Word w = Word(2) * Word(-1);
    // cout << "w = " << w << endl;
    NF nf  = NF( B , w );
    NF nf2 = NF( B , -c*w*c );
    
    pair< bool , NF > res = nf.areConjugate_uss( nf2 );
    if( res.first ) {
      cout << "Conjugacy - ok" << endl;
      if( nf==-res.second*nf2*res.second ) {
	cout << "   conjugator - ok" << endl;
      } else {
	cout << "   conjugator - error" << endl;
	exit(1);
      }
    } else {
      cout << "Conjugacy - error" << endl;
      exit(1);
    }
  }
}


//---------------------------------------------------------------------------//
//-------------------------------------- main -------------------------------//
//---------------------------------------------------------------------------//


void test_ultra_summit_conjugators( )
{
  int rank = 5;
  crag::braidgroup::BraidGroup B( rank );
  typedef ThLeftNormalForm NF;
  // typedef ThRightNormalForm NF;
  
  int exp = 100;
  for( int e=0 ; e<exp ; ++e ) {

    Word w = Word::randomWord( rank-1 , 12 );
    NF nf  = NF( B , w );
    triple< NF , NF , int > tr = nf.findUSSRepresentative( );
    const NF& nf2 = tr.first;
    cout << "Period = " << tr.third << endl;
    
    if( nf2==-tr.second*nf*tr.second ) {
      cout << "   conjugator - ok" << endl;
    } else {
      cout << "   conjugator - error" << endl;
      exit(1);
    }
    
    if( nf2.getDecomposition( ).size( )==0 ) {
      cout << "   power of Delta" << endl;
      continue;
    }

    Permutation s( rank );
    s.change(0,1);
    pair< Permutation , bool > uc = tr.first.getSimpleUltraConjugator( tr.third , s );
    Permutation sc = uc.first;
    
    
    NF conj = sc;
    NF conj_nf = -conj * tr.first * conj;
    
    if( tr.first.getPower( )==conj_nf.getPower( ) && 
	tr.first.getDecomposition( ).size( )==conj_nf.getDecomposition( ).size( ) ) {
      cout << "    SSS conj - ok" << endl;
    } else {
      cout << "    SSS conj - error" << endl;
      cout << "       " << tr.first.getPower( ) << " -> " << conj_nf.getPower( ) << endl;
      cout << "       " << tr.first.getDecomposition( ).size( ) 
	   << " -> " << conj_nf.getDecomposition( ).size( ) << endl;
      exit(1);
    }
    
    if( (-s*sc).length( )+s.length( )!=sc.length( ) ) {
      cout << "    Start - error" << endl;
      exit(1);
    } else {
      cout << "    Start - ok" << endl;
    }
    
    conj_nf.computePeriod( );
    cout << "Computations successful" << endl;
  }
}


//---------------------------------------------------------------------------//
//-------------------------------------- main -------------------------------//
//---------------------------------------------------------------------------//


Permutation find_pullback( const ThLeftNormalForm& nf , const Permutation& s )
{
  int rank = nf.getRank();
  crag::braidgroup::BraidGroup B( rank );
  typedef ThLeftNormalForm NF;
  
  Permutation result( rank );
  
  for( PermutationEnumerator PE( rank ) ; !PE.end( ) ; ++PE ) {

    Permutation c = PE.getPermutation( );
    NF conj = c;
    NF conj_nf = -conj * nf * conj;
    
    if( nf.getPower( )==conj_nf.getPower( ) && 
	nf.getDecomposition( ).size( )==conj_nf.getDecomposition( ).size( ) ) {
      
      // cout << c << " is a simple conjugator" << endl;
      Permutation transport = nf.getTransport( conj_nf , c );
      if( (-s * transport).length()+s.length()==transport.length() ) {
	if( result==Permutation( rank ) ) {
	  result = c;
	} else {
	  result = result.LeftGCD( c );
	}
	cout << "       Transport starts with s for " << c << " -> " << result << endl;
      }
    }
  }
  cout << "               -> " << result << endl;
  return result;
}

//---------------------------------------------------------------------------//
//-------------------------------------- main -------------------------------//
//---------------------------------------------------------------------------//


void test_pullback( )
{
  int rank = 10;
  crag::braidgroup::BraidGroup B( rank );
  typedef ThLeftNormalForm NF;
  // typedef ThRightNormalForm NF;
  
  int exp = 100;
  for( int e=0 ; e<exp ; ++e ) {
    
    Word w = Word::randomWord( rank-1 , 100 );
    NF nf  = NF( B , w );
    triple< NF , NF , int > tr = nf.findUSSRepresentative( );
    const NF& nf2 = tr.first;
    cout << "Period = " << tr.third << endl;
    
    if( nf2==-tr.second*nf*tr.second ) {
      cout << "   conjugator - ok" << endl;
    } else {
      cout << "   conjugator - error" << endl;
      exit(1);
    }
    
    if( nf2.getDecomposition( ).size( )==0 ) {
      cout << "   power of Delta" << endl;
      continue;
    }
    
    Permutation s( rank );
    s.change( 0 , 1 );
    cout << "   s = " << s << endl;
    Permutation pullback = nf2.getPullback( s );
    // Permutation pullback2 = find_pullback( nf2 , s );
    cout << "   pullback  = " << pullback << endl;
    // cout << "   pullback2 = " << pullback2 << endl;
    NF pb = pullback;
    NF nf3 = -pb * nf2 * pb;
    Permutation transport = nf2.getTransport( nf3 , pullback );
    cout << "   transport = " << transport << endl;

    cout << "   L: " << (-s * transport).length( ) << "  ***  "
	 << s.length( ) << "  ***  " 
	 << transport.length( ) << endl;

    // check that transport starts with s
    if( (-s * transport).length()+s.length()==transport.length() ) {
      cout << "       Transport starts with s" << endl;
    } else {
      cout << "       Transport does not start with s" << endl;
      exit(1);
    }
    /*
    if( pullback!=pullback2 ) {
      cout << "Pullbacks differ" << endl;
      exit(1);
    }
    */
  }
}


//---------------------------------------------------------------------------//
//-------------------------------------- main -------------------------------//
//---------------------------------------------------------------------------//


void test_transport( )
{
  int rank = 8;
  crag::braidgroup::BraidGroup B( rank );
  typedef ThLeftNormalForm NF;
  // typedef ThRightNormalForm NF;
  
  int exp = 100;
  for( int e=0 ; e<exp ; ++e ) {
    Word w = Word::randomWord( rank-1 , 40 );
    // Word w = Word(-1)*Word(2);
    cout << "w = " << w << endl;
    NF nf  = NF( B , w );
    triple< NF , NF , int > tr = nf.findUSSRepresentative( );
    const NF& nf2 = tr.first;
    cout << "Period = " << tr.third << endl;
    
    if( nf2==-tr.second*nf*tr.second ) {
      cout << "   conjugator - ok" << endl;
    } else {
      cout << "   conjugator - error" << endl;
      exit(1);
    }
    
    if( nf2.getDecomposition().size()==0 )
      continue;

    set< Permutation > conjugators = nf2.getSimpleSummitConjugators( );
    for( set< Permutation >::const_iterator c_it=conjugators.begin( ) ; c_it!=conjugators.end( ) ; ++c_it ) {
      NF conj = *c_it;
      NF nf3 = -conj*nf2*conj;
      
      NF transport = nf2.getTransport( nf3 , *c_it );
      
      NF nf2_right = nf2.cycle( ).first;
      NF nf3_right = nf3.cycle( ).first;
      NF nf3_right_conj = -transport * nf2_right * transport;
      if( nf3_right==nf3_right_conj ) {
	cout << "      transport ok" << endl;
      } else {
	cout << "      transport error for " << * c_it << endl;
	cout << "nf2            = " << shortenBraid( rank , nf2.getWord( ) ) << endl;
	cout << "nf3            = " << shortenBraid( rank , nf3.getWord( ) ) << endl;
	cout << "nf2_right      = " << shortenBraid( rank , nf2_right.getWord( ) ) << endl;
	cout << "nf3_right      = " << shortenBraid( rank , nf3_right.getWord( ) ) << endl;
	cout << "nf3_right_conj = " << shortenBraid( rank , nf3_right_conj.getWord( ) ) << endl;
	cout << "transport      = " << shortenBraid( rank , transport.getWord( ) ) << endl;
	cout << "==============================" << endl;
	cout << nf2 << endl;
	cout << "==============================" << endl;
	cout << nf3 << endl;
	cout << "==============================" << endl;
	cout << nf3_right << endl;
	cout << "==============================" << endl;
	cout << nf2_right << endl;
	cout << "==============================" << endl;
	cout << nf3_right_conj << endl;
	cout << "==============================" << endl;
	// shortenBraid( rank , conj_nf.getWord( ) )
	exit(1);
      }
    }
    
  }
}


//---------------------------------------------------------------------------//
//-------------------------------------- main -------------------------------//
//---------------------------------------------------------------------------//


void test_conjugacy( )
{
  int rank = 6;
  crag::braidgroup::BraidGroup B( rank );
  typedef ThLeftNormalForm NF;
  // typedef ThRightNormalForm NF;
  
  int exp = 100;
  for( int e=0 ; e<exp ; ++e ) {
    Word w = Word::randomWord( rank-1 , 20 );
    Word c = Word::randomWord( rank-1 , 20 );
    
    // Word w = Word(2) * Word(-1);
    // cout << "w = " << w << endl;
    NF nf  = NF( B , w );
    NF nf2 = NF( B , -c*w*c );
    
    pair< bool , NF > res = nf.areConjugate( nf2 );
    if( res.first ) {
      cout << "Conjugacy - ok" << endl;
      if( nf==-res.second*nf2*res.second ) {
	cout << "   conjugator - ok" << endl;
      } else {
	cout << "   conjugator - error" << endl;
	exit(1);
      }
    } else {
      cout << "Conjugacy - error" << endl;
      exit(1);
    }
  }
}


//---------------------------------------------------------------------------//
//-------------------------------------- main -------------------------------//
//---------------------------------------------------------------------------//

void test_cast_right( )
{
  int rank = 5;
  crag::braidgroup::BraidGroup B( rank );
  typedef ThLeftNormalForm NF;
  typedef ThRightNormalForm RNF;
  
  int exp = 100;
  for( int e=0 ; e<exp ; ++e ) {
    Word w = Word::randomWord( rank-1 , 20 );
    // Word w = Word(1) * Word(-2);
    // cout << "w = " << w << endl;

    NF nf = NF( B , w );
    RNF rnf = nf;

    Word w1 =  nf.getWord( );
    Word w2 = rnf.getWord( );
    
    RNF nf2( B , w1*-w2 );
    if( nf2.isTrivial( ) ) {
      cout << "Cast - ok" << endl;
    } else {
      cout << w << endl;
      cout << "Cast - error" << endl;
      exit(1);
    }
  }
}


//---------------------------------------------------------------------------//
//-------------------------------------- main -------------------------------//
//---------------------------------------------------------------------------//


void test_simple_summit_conjugators( )
{
  int rank = 3;
  crag::braidgroup::BraidGroup B( rank );
  typedef ThLeftNormalForm NF;
  // typedef ThRightNormalForm NF;
  
  int exp = 100;
  for( int e=0 ; e<exp ; ++e ) {
    Word w = Word::randomWord( rank-1 , 10 );
    // Word w = Word(2) * Word(-1);
    // cout << "w = " << w << endl;
    NF nf = NF( B , w );
    pair< NF , NF > pr = nf.findSSSRepresentative( );
    
    
    if( pr.first == (-pr.second) * nf * (pr.second) ) {
      cout << "SSS_REP - ok" << endl;
    } else {
      cout << "SSS_REP - error" << endl;
      exit(1);
    }
    
    
    Permutation s( rank );
    s.change(0,1);
    Permutation sc = pr.first.getSimpleSummitConjugator( s );
    

    NF conj = sc;
    NF conj_nf = -conj * pr.first * conj;
    
    NF conj_nf_w( B , 
		  -shortenBraid( rank , conj.getWord( ) ) * 
		  shortenBraid( rank , pr.first.getWord( ) ) * 
		  shortenBraid( rank , conj.getWord( ) ) );
    
    if( conj_nf == conj_nf_w )
      cout << "Conjugation - ok" << endl;
    else
      cout << "Conjugation - error" << endl;

    /*
    cout << "s  = " << s << endl;
    cout << "sc = " << sc << endl;
    cout << "   conj.getWord( ) = " << shortenBraid( rank , conj.getWord( ) ) << endl;
    cout << "   nf.getWord( ) = " << shortenBraid( rank , nf.getWord( ) ) << endl;
    cout << "   conj_nf.getWord( ) = " << shortenBraid( rank , conj_nf.getWord( ) ) << endl;
    cout << "   conj_nf_w.getWord( ) = " << shortenBraid( rank , conj_nf_w.getWord( ) ) << endl;
    cout << "==================================" << endl;
    cout << pr.first << endl;
    cout << "==================================" << endl;
    cout << conj_nf << endl;
    cout << "==================================" << endl;
    cout << conj << endl;
    cout << "==================================" << endl;
    */
    
    if( pr.first.getPower( )==conj_nf.getPower( ) && 
	pr.first.getDecomposition( ).size( )==conj_nf.getDecomposition( ).size( ) ) {
      cout << "    conj - ok" << endl;
    } else {
      cout << "    conj - error" << endl;
      cout << "       " << pr.first.getPower( ) << " -> " << conj_nf.getPower( ) << endl;
      cout << "       " << pr.first.getDecomposition( ).size( ) 
	   << " -> " << conj_nf.getDecomposition( ).size( ) << endl;
      exit(1);
    }
    
    
    if( (-s*sc).length( )+s.length( )!=sc.length( ) ) {
      cout << "    Start - error" << endl;
      exit(1);
    } else {
      cout << "    Start - ok" << endl;
    }
    
  }
}


//---------------------------------------------------------------------------//
//-------------------------------------- main -------------------------------//
//---------------------------------------------------------------------------//


void test_sss_rep( )
{
  int rank = 20;
  crag::braidgroup::BraidGroup B( rank );
  typedef ThLeftNormalForm NF;
  // typedef ThRightNormalForm NF;
  
  int exp = 100;
  for( int e=0 ; e<exp ; ++e ) {
    Word w = Word::randomWord( rank-1 , 200 );
    // Word w = Word(2) * Word(-1);
    // cout << "w = " << w << endl;
    NF nf = NF( B , w );
    pair< NF , NF > pr = nf.findSSSRepresentative( );
    
    
    if( pr.first == (-pr.second) * nf * (pr.second) ) {
      cout << "SSS_REP - ok" << endl;
    } else {
      cout << "SSS_REP - error" << endl;
      exit(1);
    }

    
    
  }
  
}


//---------------------------------------------------------------------------//
//-------------------------------------- main -------------------------------//
//---------------------------------------------------------------------------//


void test_cycling( )
{
  int rank = 5;
  crag::braidgroup::BraidGroup B( rank );
  typedef ThLeftNormalForm NF;
  // typedef ThRightNormalForm NF;
  
  int exp = 100;
  for( int e=0 ; e<exp ; ++e ) {
    Word w = Word::randomWord( rank-1 , 20 );
    // Word w = Word(2) * Word(-1);
    // cout << "w = " << w << endl;
    NF nf = NF( B , w );
    cout << "NF = " << nf << endl;
    
    pair< NF , NF > pr = nf.cycle( );
    
    // cout << "1st = " << pr.first << endl;
    // cout << "2nd = " << pr.second << endl;
    
    if( pr.first == (-pr.second) * nf * (pr.second) ) {
      cout << "Cycling - ok" << endl;
    } else {
      cout << "Cycling - error" << endl;
      exit(1);
    }
  }
}


//---------------------------------------------------------------------------//
//-------------------------------------- main -------------------------------//
//---------------------------------------------------------------------------//


void test_decycling( )
{
  int rank = 5;
  crag::braidgroup::BraidGroup B( rank );
  typedef ThLeftNormalForm NF;
  // typedef ThRightNormalForm NF;
  
  int exp = 100;
  for( int e=0 ; e<exp ; ++e ) {
    Word w = Word::randomWord( rank-1 , 20 );
    NF nf = NF( B , w );
    pair< NF , NF > pr = nf.decycle( );
    
    if( pr.first == (-pr.second) * nf * (pr.second) ) {
      cout << "DeCycling - ok" << endl;
    } else {
      cout << "DeCycling - error" << endl;
      exit(1);
    }
  }
}


//---------------------------------------------------------------------------//
//-------------------------------------- main -------------------------------//
//---------------------------------------------------------------------------//


void test_WordProblem( )
{
  int rank = 20;
  crag::braidgroup::BraidGroup B( rank );
  typedef ThLeftNormalForm NF;
  
  int exp = 100;
  for (int e = 0; e < exp; ++e) {
    Word w = Word::randomWord( rank-1 , 400 );
    // cout << "w = " << w << endl;
    
    // compute a normal form
    NF nf = NF( B , w );
    // cout << "NF = " << endl << nf << endl;
    
    // compute a normal form of the inverse
    NF nf_inv = nf.inverse( );
    // cout << "NF_inv = " << endl << nf_inv << endl;
    
    // multiply normal forms
    NF nf_mult = nf*nf_inv;
    if( !nf_mult.isTrivial( ) ) {
      cout << "Shit A!!!" << endl;
      exit( 1 );
    } else {
      cout << "Inverse ok" << endl;
    }

    
    Word w1 = nf_inv.getWord();
    NF nf2(B, w * w1);
    if( !nf2.isTrivial( ) ) {
      cout << "Shit B!!!" << endl;
      exit( 1 );
    } else {
      cout << "getWord - ok" << endl;
    }

    cout << "=============================" << endl;
  }

}


//---------------------------------------------------------------------------//
//-------------------------------------- main -------------------------------//
//---------------------------------------------------------------------------//


void test_inverse( )
{
  int rank = 10;
  crag::braidgroup::BraidGroup B( rank );
  typedef ThLeftNormalForm NF;
  // typedef ThRightNormalForm NF;
  
  Permutation Delta = Permutation::getHalfTwistPermutation( rank );
  Permutation triv = Permutation( rank );

  int exp = 100;
  for( int e=0 ; e<exp ; ++e ) {
    cout << "e = " << e << endl;
    
    Word w = Word::randomWord( rank-1 , 100 );
    NF nf  = NF( B , w );
    
    NF nf1 = NF( B , -w );
    NF nf2 = -nf;
    if( nf1==nf2 ) {
      cout << "Inversion - ok" << endl;
    } else {
      cout << "Inversion - error" << endl;
      exit(1);
    }
  }
}

//---------------------------------------------------------------------------//
//-------------------------------------- main -------------------------------//
//---------------------------------------------------------------------------//


void test_multiplication( )
{
  int rank = 10;
  crag::braidgroup::BraidGroup B( rank );
  typedef ThLeftNormalForm NF;
  // typedef ThRightNormalForm NF;
  
  Permutation Delta = Permutation::getHalfTwistPermutation( rank );
  Permutation triv = Permutation( rank );

  int exp = 100;
  for( int e=0 ; e<exp ; ++e ) {
    cout << "e = " << e << endl;
    
    Word w1 = Word::randomWord( rank-1 , 100 );
    NF nf1  = NF( B , w1 );
    Word w2 = Word::randomWord( rank-1 , 100 );
    NF nf2  = NF( B , w2 );
    
    NF nf3 = NF( B , w1*w2 );
    NF nf4 = nf1*nf2;
    if( nf3==nf4 ) {
      cout << "Multiplication - ok" << endl;
    } else {
      cout << "Multiplication - error" << endl;
      exit(1);
    }
  }
}

//---------------------------------------------------------------------------//
//-------------------------------------- main -------------------------------//
//---------------------------------------------------------------------------//


void test_nf( )
{
  int rank = 10;
  crag::braidgroup::BraidGroup B( rank );
  typedef ThLeftNormalForm NF;
  // typedef ThRightNormalForm NF;
  
  Permutation Delta = Permutation::getHalfTwistPermutation( rank );
  Permutation triv = Permutation( rank );

  int exp = 100;
  for( int e=0 ; e<exp ; ++e ) {
    
    cout << "e = " << e << endl;
    
    Word w = Word::randomWord( rank-1 , 100 );
    NF nf  = NF( B , w );
    
    list< Permutation > D = nf.getDecomposition( );
    for( list< Permutation >::const_iterator d_it=D.begin() ; d_it!=D.end( ) ; ) {
      Permutation d1 = *(d_it++);
      if( d_it==D.end( ) )
	break;
      Permutation d2 = *d_it;
      Permutation a = d2.LeftGCD( -d1*Delta );
      // cout << a << endl;
      if( a!=triv ) {
	cout << "Error" << endl;
	exit(1);
      }
      
    }
    // cout << "=================================" << endl;
  }
}


//---------------------------------------------------------------------------//
//-------------------------------------- main -------------------------------//
//---------------------------------------------------------------------------//


void test_transport_gcd_properties( )
{
  int rank = 5;
  crag::braidgroup::BraidGroup B( rank );
  typedef ThLeftNormalForm NF;
  // typedef ThRightNormalForm NF;
  
  int exp = 1;
  for( int e=0 ; e<exp ; ++e ) {
    Word w = Word::randomWord( rank-1 , 10 );
    NF nf = NF( B , w ).findUSSRepresentative( ).first;

    map< Permutation , Permutation > conjugators;
    for( PermutationEnumerator PE( rank ) ; !PE.end( ) ; ++PE ) {
      
      Permutation c = PE.getPermutation( );
      NF conj = c;
      NF conj_nf = -conj * nf * conj;
      
      if( nf.getPower( )==conj_nf.getPower( ) && 
	  nf.getDecomposition( ).size( )==conj_nf.getDecomposition( ).size( ) ) {
	Permutation tr = nf.getTransport( conj_nf , c );
	conjugators[c] = tr;
	cout << c << endl;
      }
    }
    
    for( map< Permutation , Permutation >::const_iterator c_it=conjugators.begin( ) ; c_it!=conjugators.end( ) ; ++c_it ) {
      map< Permutation , Permutation >::const_iterator d_it=c_it;
      for( ++d_it ; d_it!=conjugators.end( ) ; ++d_it ) {

	Permutation c = (*c_it).first;
	Permutation d = (*d_it).first;
	Permutation ct = (*c_it).second;
	Permutation dt = (*d_it).second;
	
	Permutation m  = c.LeftGCD( d );
	Permutation mt = ct.LeftGCD( dt );
	if( m!=Permutation(rank) ) {
	  cout << c << "  ^  " << d << " = " << m << endl;
	  cout << ct << "  ^  " << dt << " = " << mt << endl;
	  cout << endl;

	  if( conjugators.find( m )==conjugators.end( ) ) {
	    cout << "Shit 1" << endl;
	    exit(1);
	  }
	  if( conjugators[m]!=mt ) {
	    cout << "Shit 2" << endl;
	    exit(1);
	  }
	}
      }
      
    }
  }
}


//---------------------------------------------------------------------------//
//-------------------------------------- main -------------------------------//
//---------------------------------------------------------------------------//


int main( )
{
  srand( time(0) );
  test_WordProblem( );
  // test_cast_right( );
  // test_multiplication( );
  // test_inverse( );
  // test_cycling( );
  // test_decycling( );
  // test_sss_rep( );
  // test_simple_summit_conjugators( );
  // test_conjugacy( );
  // test_transport( );
  // test_transport_gcd_properties( );

  // test_pullback( );

  // test_ultra_summit_conjugators( );
  
  // test_conjugacy_uss( );
  
  return 0;
}
