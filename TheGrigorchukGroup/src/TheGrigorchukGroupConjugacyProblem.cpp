// Copyright (C) 2006 Alexander Ushakov
// Contents: Part of implementation of class TheGrigorchukGroupAlgorithms
// Conjugacy Problem algorithm
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//


#include "Word.h"
#include "iostream"
#include "TheGrigorchukGroupAlgorithms.h"


//---------------------------------------------------------------------------//
//----------------------- TheGrigorchukGroupAlgorithms ----------------------//
//---------------------------------------------------------------------------//


triple< int , int , int > TheGrigorchukGroupAlgorithms::abelianImage( const Word& w )
{
  triple< int , int , int > result(0,0,0);
  
  for( Word::const_iterator w_it=w.begin( ) ; w_it!=w.end() ; ++w_it ) {
    int g = *w_it;
    g = g>0 ? g : -g;
    switch( g ) {
    case 1:
      result.first = 1-result.first;
      break;
    case 2:
      result.second = 1-result.second;
      break;
    case 3:
      result.third  = 1-result.third;
      break;
    case 4:
      result.second = 1-result.second;
      result.third  = 1-result.third;
      break;
    }
  }
  
  return result;
}


//---------------------------------------------------------------------------//
//----------------------- TheGrigorchukGroupAlgorithms ----------------------//
//---------------------------------------------------------------------------//


int d=0;
ostream& margin( ostream& os )
{
  for( int i=0 ; i<d ; ++i )
    os << "  ";
  return os;
}


set< int > TheGrigorchukGroupAlgorithms::conjugate_Kcosets( const Word& w1 , const Word& w2 )
{
  Word r1 = reduce(w1);
  Word r2 = reduce(w2);
  
  // A. Check if words are not equal in the abelian quotient
  typedef triple< int , int , int > TIII;
  TIII a1 = abelianImage(r1);
  TIII a2 = abelianImage(r2);
  if( a1!=a2 )
    return set< int >( );

  // B. The basis case (words have lengths at most 1)
  if( r1.length( )<=1 && r2.length( )<=1 ) {
    // since abelinization is passed r1 must be equal to r2 
    set< int > result;
    if( r1.length( )==0 ) {
      for( int i=0 ; i<16 ; ++i )
	result.insert(i);
    } else if( r1==Word(1) ) {
      result.insert( 0 );  // 1
      result.insert( 3 );  // dad
      result.insert( 4 );  // adad
      result.insert( 7 );  // a
    } else if( r1==Word(2) ) {
      result.insert( 0 );  // 1
      result.insert( 1 );  // d
      result.insert( 8 );  // b
      result.insert( 9 );  // c
    } else if( r1==Word(3) ) {
      result.insert( 0 );  // 1
      result.insert( 1 );  // d
      result.insert( 8 );  // b
      result.insert( 9 );  // c
    } else if( r1==Word(4) ) {
      result.insert( 0 );  // 1
      result.insert( 1 );  // d
      result.insert( 8 );  // b
      result.insert( 9 );  // c
      result.insert( 5 );  // ada
      result.insert( 4 );  // adad
      result.insert( 13 ); // bada
      result.insert( 12 ); // badad
    }
    return result;
  }
  
  // -----------------------------------------------
  // C. If both elements belong to St(1)
  if( a1.first==0 ) {

    // C1. Split the words
    pair< Word , Word > sp1 = split( r1 );
    pair< Word , Word > sp2 = split( r2 );

    // C2. Assume the conjugator belongs to St(1)
    set< int > K1 = conjugate_Kcosets( sp1.first , sp2.first );
    set< int > K2 = conjugate_Kcosets( sp1.second , sp2.second );
    set< int > K12 = liftPairsKcosetsUP( K1 , K2 );

    // C3. Assume the conjugator does not belong to St(1)
    set< int > K3 = conjugate_Kcosets( sp1.first , sp2.second );
    set< int > K4 = conjugate_Kcosets( sp1.second , sp2.first );
    set< int > K34 = liftPairsKcosetsUP( K3 , K4 );


    // C4. Add a at the end (since we are outside of ST(1)) and merge classes together
    for( set< int >::const_iterator K_it=K34.begin( ) ; K_it!=K34.end( ) ; ++K_it )
      K12.insert( cosetRepresentativeKSbgp( cosetRepresentativeKSbgp( *K_it ).push_back( 1 ) ) );
    return K12;
  } 
  
  // -----------------------------------------------
  // D. If both elements do not belong to St(1)
  if( a1.first!=0 ) {

    // D1. Split the words
    push_back( r1 , 1 );
    push_back( r2 , 1 );
    pair< Word , Word > sp1 = split( r1 );
    pair< Word , Word > sp2 = split( r2 );
    
    // D2. Assume the conjugator belongs to St(1)
    set< int > K1 = conjugate_Kcosets( sp1.first*sp1.second , sp2.first*sp2.second );
    set< int > K12;
    for( set< int >::const_iterator K1_it=K1.begin() ; K1_it!=K1.end( ) ; ++K1_it ) {
      int L = liftPairKcosetsUP( *K1_it , cosetRepresentativeKSbgp( sp1.second * cosetRepresentativeKSbgp( *K1_it ) * -sp2.second ) );
      if( L!=-1 )
	K12.insert( L );
    }
    
    
    // D3. Assume the conjugator does not belong to St(1)
    set< int > K3 = conjugate_Kcosets( sp1.first*sp1.second , sp2.second*sp2.first );
    set< int > K34;
    for( set< int >::const_iterator K3_it=K3.begin() ; K3_it!=K3.end( ) ; ++K3_it ) {

      int L = liftPairKcosetsUP( *K3_it , cosetRepresentativeKSbgp( sp1.second * cosetRepresentativeKSbgp( *K3_it ) * -sp2.first ) );
      if( L!=-1 )
	K34.insert( L );
    }
    
    // D4. Add a at the end (since we are outside of ST(1)) and merge classes together
    for( set< int >::const_iterator K_it=K34.begin( ) ; K_it!=K34.end( ) ; ++K_it )
      K12.insert( cosetRepresentativeKSbgp( cosetRepresentativeKSbgp( *K_it ).push_back(1) ) );
    return K12;
  }

  // cannot reach this point
  cout << "Pizdetz" << endl;
  exit(1);

}


//---------------------------------------------------------------------------//
//----------------------- TheGrigorchukGroupAlgorithms ----------------------//
//---------------------------------------------------------------------------//


set< Word > liftPairsUP( const set< Word >& K1 , const set< Word >& K2 )
{
  set< Word > result;
  
  for( set< Word >::const_iterator K1_it=K1.begin( ) ;  K1_it!=K1.end( ) ; ++K1_it ) {
    for( set< Word >::const_iterator K2_it=K2.begin( ) ;  K2_it!=K2.end( ) ; ++K2_it ) {
      pair< Word , Word > L = TheGrigorchukGroupAlgorithms::liftToSTone( *K1_it , *K2_it );
      if( L.second.length( )==0 )
	result.insert( L.first );
    }
  }
  
  return result;
}


set< Word > merge_and_reduce( const set< Word >& K1 , const set< Word >& K2 )
{
  map< int , Word > representatives;
  for( set< Word >::const_iterator K_it=K1.begin( ) ;  K_it!=K1.end( ) ; ++K_it ) {
    int n = TheGrigorchukGroupAlgorithms::cosetRepresentativeKSbgp( *K_it );
    map< int , Word >::const_iterator r_it = representatives.find(n);
    if( r_it==representatives.end( ) || (*r_it).second.length( )>(*K_it).length( ) )
      representatives[n] = *K_it;
  }

  for( set< Word >::const_iterator K_it=K2.begin( ) ;  K_it!=K2.end( ) ; ++K_it ) {
    int n = TheGrigorchukGroupAlgorithms::cosetRepresentativeKSbgp( *K_it );
    map< int , Word >::const_iterator r_it = representatives.find(n);
    if( r_it==representatives.end( ) || (*r_it).second.length( )>(*K_it).length( ) )
      representatives[n] = *K_it;
  }

  // prepare the output
  set< Word > result;
  for( map< int , Word >::const_iterator r_it=representatives.begin() ; r_it!=representatives.end() ; ++r_it ) 
    result.insert( (*r_it).second );
  
  return result;
}


set< Word > TheGrigorchukGroupAlgorithms::findConjugator_Kcosets( const Word& w1 , const Word& w2 )
{
  Word r1 = reduce(w1);
  Word r2 = reduce(w2);
  
  // A. Check if words are not equal in the abelian quotient
  typedef triple< int , int , int > TIII;
  TIII a1 = abelianImage(r1);
  TIII a2 = abelianImage(r2);
  if( a1!=a2 )
    return set< Word >( );

  // B. The basis case (words have lengths at most 1)
  if( r1.length( )<=1 && r2.length( )<=1 ) {
    // since abelinization is passed r1 must be equal to r2 
    set< Word > result;
    if( r1.length( )==0 ) {
      for( int i=0 ; i<16 ; ++i )
	result.insert( cosetRepresentativeKSbgp(i) );
    } else if( r1==Word(1) ) {
      result.insert( Word( ) );                                 // 1
      result.insert( Word(4)*Word(1)*Word(4) );                 // dad
      result.insert( Word(1)*Word(4)*Word(1)*Word(4) );         // adad
      result.insert( Word(1) );                                 // a
    } else if( r1==Word(2) ) {
      result.insert( Word( ) );                                 // 1
      result.insert( Word(4) );                                 // d
      result.insert( Word(2) );                                 // b
      result.insert( Word(3) );                                 // c
    } else if( r1==Word(3) ) {
      result.insert( Word( ) );                                 // 1
      result.insert( Word(4) );                                 // d
      result.insert( Word(2) );                                 // b
      result.insert( Word(3) );                                 // c
    } else if( r1==Word(4) ) {
      result.insert( Word( ) );                                 // 1
      result.insert( Word(4) );                                 // d
      result.insert( Word(2) );                                 // b
      result.insert( Word(3) );                                 // c
      result.insert( Word(1)*Word(4)*Word(1) );                 // ada
      result.insert( Word(1)*Word(4)*Word(1)*Word(4) );         // adad
      result.insert( Word(2)*Word(1)*Word(4)*Word(1) );         // bada
      result.insert( Word(2)*Word(1)*Word(4)*Word(1)*Word(4) ); // badad
    }
    return result;
  }
  

  // -----------------------------------------------
  // C. If both elements belong to St(1)
  if( a1.first==0 ) {

    // C1. Split the words
    pair< Word , Word > sp1 = split( r1 );
    pair< Word , Word > sp2 = split( r2 );

    // C2. Assume the conjugator belongs to St(1)
    set< Word > K12 = liftPairsUP( findConjugator_Kcosets( sp1.first , sp2.first ) , findConjugator_Kcosets( sp1.second , sp2.second ) );
    
    // C3. Assume the conjugator does not belong to St(1)
    set< Word > K34 = liftPairsUP( findConjugator_Kcosets( sp1.first , sp2.second ) , findConjugator_Kcosets( sp1.second , sp2.first ) );

    // C4. Add a at the end (since we are outside of ST(1)) and merge classes together
    while( !K34.empty( ) ) {
      K12.insert( *K34.begin( ) * Word(1) );
      K34.erase( K34.begin( ) );
    }
    return merge_and_reduce( K12 , K34 );
  } 
  
  // -----------------------------------------------
  // D. If both elements do not belong to St(1)
  if( a1.first!=0 ) {

    // D1. Split the words
    push_back( r1 , 1 );
    push_back( r2 , 1 );
    pair< Word , Word > sp1 = split( r1 );
    pair< Word , Word > sp2 = split( r2 );

    // D2. Assume the conjugator belongs to St(1)
    set< Word > K1 = findConjugator_Kcosets( sp1.first*sp1.second , sp2.first*sp2.second );
    set< Word > K12;
    for( set< Word >::const_iterator K1_it=K1.begin() ; K1_it!=K1.end( ) ; ++K1_it ) {
      pair< Word , Word > L = TheGrigorchukGroupAlgorithms::liftToSTone( *K1_it , sp1.second * *K1_it * -sp2.second );
      if( L.second.length( )==0 )
	K12.insert( L.first );
    }
    
    // D3. Assume the conjugator does not belong to St(1)
    set< Word > K3 = findConjugator_Kcosets( sp1.first*sp1.second , sp2.second*sp2.first );
    set< Word > K34;
    for( set< Word >::const_iterator K3_it=K3.begin() ; K3_it!=K3.end( ) ; ++K3_it ) {

      pair< Word , Word > L = TheGrigorchukGroupAlgorithms::liftToSTone( *K3_it , sp1.second * *K3_it * -sp2.first );
      if( L.second.length( )==0 )
	K34.insert( L.first );
    }
    
    // D4. Add a at the end (since we are outside of ST(1)) and merge classes together
    while( !K34.empty( ) ) {
      K12.insert( *K34.begin( ) * Word(1) );
      K34.erase( K34.begin( ) );
    }
    return merge_and_reduce( K12 , K34 );
  }

  cout << "Jopa !!!" << endl;
  exit(1);
}
