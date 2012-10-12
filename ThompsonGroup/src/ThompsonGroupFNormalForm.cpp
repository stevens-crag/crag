// Copyright (C) 2007 Alexander Ushakov
// Contents: Implementation of class ThompsonGroupFNormalForm
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//



#include "ThompsonGroupFNormalForm.h"


//---------------------------------------------------------------------------//
//------------------------- ThompsonGroupFNormalForm ------------------------//
//---------------------------------------------------------------------------//


ThompsonGroupFNormalForm::ThompsonGroupFNormalForm( )
{
  
}


ThompsonGroupFNormalForm::ThompsonGroupFNormalForm( const Word& w ) : 
  Word(removeBadPairs( semiNormalFormFor( w ) ) )
{
  
}


ostream& operator << ( ostream& os , const ThompsonGroupFNormalForm& nf )
{
  os << (Word)nf  << endl;
}


//---------------------------------------------------------------------------//
//---------------------------- semiNormalFormFor ----------------------------//
//---------------------------------------------------------------------------//


Word ThompsonGroupFNormalForm::semiNormalFormFor( const Word& w )
{
  // partition of the word
  list< list< int > > units;

  for( Word::const_iterator w_it=w.begin( ) ; w_it!=w.end( ) ; ) {

    // 
    units.push_back( list< int >( 1 , *w_it ) );
    ++w_it;
    while( units.size( )>1 && ( (*----units.end( )).size( )<=(*--units.end( )).size( ) || w_it==w.end( ) ) ) {

      list< int > new_unit = mergeUnits( *----units.end( ) , *--units.end( ) );
      units.pop_back( );
      units.pop_back( );
      if( new_unit.size( ) )
        units.push_back( new_unit );
    }
    
    
  }
  
  if( units.size( )==0 )
    return Word( );

  return Word( *units.begin( ) );
}


//---------------------------------------------------------------------------//
//-------------------------------- nextOperation ----------------------------//
//---------------------------------------------------------------------------//


int ThompsonGroupFNormalForm::nextOperation( int n1 , int n2 , int n3 , int n4 )
{
  int n_1 = n1;
  int n_2 = n2;
  int n_3 = n3;
  int n_4 = n4;
  if( n1<0 ) n_1 =  1000000;
  if( n2>0 ) n_2 = -1000000;
  if( n3<0 ) n_3 =  1000000;
  if( n4>0 ) n_4 = -1000000;
  n_2 = -n_2;
  n_4 = -n_4;

  if( n_1<=n_2 && n_1<=n_3 && n_1<=n_4 )
    return 1;

  if( n_4<=n_1 && n_4<=n_2 && n_4<=n_3 )
    return 4;

  if( n_2<=n_1 && n_2<=n_3 && n_2<=n_4 ) {
    if( n2==n3 )
      return 5;
    return 2;
  }

  if( n_3<=n_1 && n_3<=n_2 && n_3<=n_4 ) {
    if( n2==n3 )
      return 6;
    return 3;
  }
}


//---------------------------------------------------------------------------//
//------------------------------- mergeUnits --------------------------------//
//---------------------------------------------------------------------------//  


list< int > ThompsonGroupFNormalForm::mergeUnits( list< int >& unit1 , list< int >& unit2 )
{
  int d1 = 0;
  int d2 = 0;
  list< int > w1;
  list< int > w2;

  // a -- the first symbol of the first word
  // b -- the last symbol of the first word
  // c -- the first symbol of the second word
  // d -- the last symbol of the second word
  int a, b, c, d;

  // while at least one word is not empty
  while( !unit1.empty( ) || !unit2.empty( ) ) {

    // 1. Compute a and b
    if( !unit1.empty( ) ) {
      a = *unit1.begin( );
      b = *--unit1.end( );
    } else {
      a =  1000000;
      b = -1000000;
    }
    
    // 2. Compute c and d
    if( !unit2.empty( ) ) {
      c = *unit2.begin( );
      d = *--unit2.end( );
    } else {
      c =  1000000;
      d = -1000000;
    }

    // 3. Need to adjust a,b,c,d by d1 and d2 resp.
    int na = a<0 ? a-d1 : a+d1;
    int nb = b<0 ? b-d1 : b+d1;
    int nc = c<0 ? c-d2 : c+d2;
    int nd = d<0 ? d-d2 : d+d2;
    
    // 4. Check if we have a cancellation in between
    if( nb+nc==0 && c!=1000000 ) {
      unit1.erase( --unit1.end( ) );
      unit2.erase( unit2.begin( ) );
      continue;
    }

    // 5. Perform the next operation 
    switch( nextOperation( na , nb , nc , nd ) ) {
    case 1: // a)
      w1.push_back( na );
      unit1.erase( unit1.begin( ) );
      continue;
    case 2: // c)
      w2.push_front( nb );
      unit1.erase( --unit1.end( ) );
      ++d2;
      continue;
    case 3: // d)
      w1.push_back( nc );
      unit2.erase( unit2.begin( ) );
      ++d1;
      continue;
    case 4: // b)
      w2.push_front( nd );
      unit2.erase( --unit2.end( ) );
      break;
    case 5:
      // cout << "  Op5" << endl;
      w2.push_front( nc );
      unit2.erase( unit2.begin( ) );
      continue;
    case 6:
      // cout << "  Op6" << endl;
      w1.push_back( nb );
      unit1.erase( --unit1.end( ) );
      continue;
    }
  }
  
  for( list< int >::iterator w_it=w2.begin( ) ; w_it!=w2.end( ) ; ++w_it ) 
    w1.push_back( Generator( *w_it ) );
  return w1;
}


//---------------------------------------------------------------------------//
//------------------------------ removeBadPairs -----------------------------//
//---------------------------------------------------------------------------//
      
      
Word ThompsonGroupFNormalForm::removeBadPairs( const Word& w )
{
  list< int > lst = w.getList( );

  list< int >::iterator l_it1 = lst.begin( );
  for( ; l_it1!=lst.end( ) && *l_it1>0 ; ++l_it1 );
  list< int >::iterator l_it2 = l_it1;

  list< int > l_part;
  list< int > r_part;
  list< int > l_weight;
  list< int > r_weight;

  while( l_it1!=lst.begin( ) || l_it2!=lst.end( ) ) {

    if( l_it1==lst.begin( ) ) {
      r_part.push_back( *(l_it2++) );
      r_weight.push_back( 0 );
      continue;
    }
    if( l_it2==lst.end( ) ) {
      l_part.push_front( *(--l_it1) );
      l_weight.push_front( 0 );
      continue;
    }

    list< int >::iterator l_it = l_it1;
    int a = *--l_it;
    int b = *l_it2;

    if( a<-b ) {
      r_part.push_back( *(l_it2++) );
      r_weight.push_back( 0 );
      continue;
    }
    if( a>-b ) {
      l_part.push_front( *(--l_it1) );
      l_weight.push_front( 0 );
      continue;
    }

    int c = 1000000;
    int d = 1000000;
    int p1 = 0;
    int p2 = 0;

    bool erase = false;
    if( !l_part.empty( ) ) {
      c = *l_part.begin( );
      p1 = *l_weight.begin( );
    }
    if( !r_part.empty( ) ) {
      d = *--r_part.end( );
      p2 = *--r_weight.end( );
    }
    c = c<0 ? c+p1 : c-p1;
    d = d<0 ? d+p2 : d-p2;

/*
    cout << "  c = " << c << endl;
    cout << "  d = " << d << endl;
    cout << "  p1 = " << p1 << endl;
    cout << "  p2 = " << p2 << endl;
*/

    if( c!=a+1 && d!=b-1 && c!=a && d!=b ) {
      // cout << "  erase" << endl;
      --l_it1;
      ++l_it2;
      if( !l_weight.empty( ) )
        ++(*l_weight.begin( ));
      if( !r_weight.empty( ) )
        ++(*--r_weight.end( ));
      continue;
    }
    r_part.push_back( *(l_it2++) );
    r_weight.push_back( 0 );
    l_part.push_front( *(--l_it1) );
    l_weight.push_front( 0 );
  }

/*
  cout << "l_part = " << Word( l_part ) << endl;
  cout << "r_part = " << Word( r_part ) << endl;
  cout << "l_weight = " << Word( l_weight ) << endl;
  cout << "r_weight = " << Word( r_weight ) << endl;
*/

  int d=0;
  list< int >::iterator li_it1 = l_part.begin( );
  list< int >::iterator li_it2 = l_weight.begin( );
  for( ; li_it1!=l_part.end( ) ; ++li_it1 , ++li_it2 ) {
    d += *li_it2;
    *li_it1 -= d;
  }

  d = 0;
  li_it1 = r_part.end( );
  li_it2 = r_weight.end( );
  for( ; li_it1!=r_part.begin( ) ; ) {
    --li_it1;
    --li_it2;
    d += *li_it2;
    *li_it1 += d;
  }

  // cout << "l_part = " << Word( l_part ) << endl;
  // cout << "r_part = " << Word( r_part ) << endl;

  return Word( l_part ) * Word( r_part );
}
