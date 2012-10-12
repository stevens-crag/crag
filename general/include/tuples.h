

#ifndef _tuples_H_
#define _tuples_H_


#include "iostream"
using namespace std;


//---------------------------------------------------------------------------//
//-------------------------------- triple -----------------------------------//
//---------------------------------------------------------------------------//


template< class T1 , class T2 , class T3 > 
class triple
{

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Constructors                                       //
  //                                                     //
  /////////////////////////////////////////////////////////

public:

  triple( const T1& t1=T1() , const T2& t2=T2() , const T3& t3=T3() ) : first(t1), second(t2), third(t3) { }


  /////////////////////////////////////////////////////////
  //                                                     //
  //  Accessors:                                         //
  //                                                     //
  /////////////////////////////////////////////////////////

 public:

  bool operator == ( const triple& t ) const {
    return first==t.first && second==t.second && third==t.third;
  }

  bool operator != ( const triple& t ) const {
    return first!=t.first || second!=t.second || third!=t.third;
  }

  bool operator> ( const triple& t ) const {
    if( first>t.first )
      return true;
    if( first<t.first )
      return false;

    if( second>t.second )
      return true;
    if( second<t.second )
      return false;

    if( third>t.third )
      return true;
    return false;
  }
  
  bool operator< ( const triple& t ) const {
    if( first<t.first )
      return true;
    if( first>t.first )
      return false;

    if( second<t.second )
      return true;
    if( second>t.second )
      return false;

    if( third<t.third )
      return true;
    return false;
  }

  /////////////////////////////////////////////////////////
  //                                                     //
  //  I/O:                                               //
  //                                                     //
  /////////////////////////////////////////////////////////

public:

  friend ostream& operator << ( ostream& os , const triple& t ) {
    os << "<" << t.first << "," << t.second << "," << t.third << ">";
    return os;
  }

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Data members                                       //
  //                                                     //
  /////////////////////////////////////////////////////////

public:

  T1 first;
  T2 second;
  T3 third;
};



//---------------------------------------------------------------------------//
//------------------------------- quadruple ---------------------------------//
//---------------------------------------------------------------------------//


template< class T1 , class T2 , class T3 , class T4 >
class quadruple
{

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Constructors                                       //
  //                                                     //
  /////////////////////////////////////////////////////////

public:

  quadruple( const T1& t1=T1() , const T2& t2=T2() , const T3& t3=T3() , const T4& t4=T4() ) : first(t1), second(t2), third(t3), fourth(t4) { }


  /////////////////////////////////////////////////////////
  //                                                     //
  //  Accessors:                                         //
  //                                                     //
  /////////////////////////////////////////////////////////

 public:

  bool operator< ( const quadruple& q ) const {
    if( first<q.first )
      return true;
    if( first>q.first )
      return false;

    if( second<q.second )
      return true;
    if( second>q.second )
      return false;

    if( third<q.third )
      return true;
    if( third>q.third )
      return false;

    if( fourth<q.fourth )
      return true;
    return false;
  }

  /////////////////////////////////////////////////////////
  //                                                     //
  //  I/O:                                               //
  //                                                     //
  /////////////////////////////////////////////////////////

public:

  friend ostream& operator << ( ostream& os , const quadruple& t ) {
    os << "<" << t.first << "," << t.second << "," << t.third << "," << t.fourth << ">";
    return os;
  }

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Data members                                       //
  //                                                     //
  /////////////////////////////////////////////////////////

public:

  T1 first;
  T2 second;
  T3 third;
  T4 fourth;
};


//---------------------------------------------------------------------------//
//------------------------------- quintuple ---------------------------------//
//---------------------------------------------------------------------------//


template< class T1 , class T2 , class T3 , class T4 , class T5 >
class quintuple
{

  /////////////////////////////////////////////////////////
  //                                                     //
  //  Constructors                                       //
  //                                                     //
  /////////////////////////////////////////////////////////

public:
  
  quintuple( const T1& t1=T1() , const T2& t2=T2() , const T3& t3=T3() , const T4& t4=T4() , const T5& t5=T5() ) : first(t1), second(t2), third(t3), fourth(t4), fifth(t5) { }


  /////////////////////////////////////////////////////////
  //                                                     //
  //  Accessors:                                         //
  //                                                     //
  /////////////////////////////////////////////////////////

 public:

  bool operator< ( const quintuple& q ) const {
    if( first<q.first )
      return true;
    if( first>q.first )
      return false;

    if( second<q.second )
      return true;
    if( second>q.second )
      return false;

    if( third<q.third )
      return true;
    if( third>q.third )
      return false;

    if( fourth<q.fourth )
      return true;
    if( fourth>q.fourth )
      return false;

    if( fifth<q.fifth )
      return true;
    return false;
  }
  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  I/O:                                               //
  //                                                     //
  /////////////////////////////////////////////////////////

public:

  friend ostream& operator << ( ostream& os , const quintuple& t ) {
    os << "<" << t.first << "," << t.second << "," << t.third << "," << t.fourth << "," << t.fifth << ">";
    return os;
  }

  
  /////////////////////////////////////////////////////////
  //                                                     //
  //  Data members                                       //
  //                                                     //
  /////////////////////////////////////////////////////////

public:

  T1 first;
  T2 second;
  T3 third;
  T4 fourth;
  T4 fifth;
};

#endif

