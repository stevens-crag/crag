// Copyright (C) 2005 Alexander Ushakov, Alexey Myasnikov
// Contents: Definition of length attack
//
// Principal Authors: Alexander Ushakov, Alexey Myasnikov
//
// Revision History:
//

#ifndef _LengthAttack_H_
#define _LengthAttack_H_

#include "Word.h"
#include <vector>
using namespace std;

enum findKey_LengthBasedResult { 
  FAILED , 
  TIME_EXPIRED ,
  SUCCESSFULL
};

#define AL1 1
#define AL2 2
#define AL3 3

typedef pair< int , vector< Word > > ELT;

//
//  LENGTH-BASED ATTACK CLASS INTERFACE
//
//! Basis interface for classes implementing the Length-Based attack
class LengthAttackBase
{
public:
  //! Returns the algorithm type
  virtual int type()=0;

  //! Attack on an instances of the AAG protocol
  /*! Executes a length-based attack on an instance of the AAG protocol 
    \param N  - rank of the braid group (number of strands)
    \param A1 - Alices subgroup generators
    \param A2 - Alices subgroup generators conjugated by Bob's private key
    \param B  - Bob's subgroup (Bob's private key belongs to it)
    \param sec - amount of time (in seconds) given to procedure to finish
    \return  - Returns findKey_LengthBasedResult::SUCCESSFULL if the attack succeeds, findKey_LengthBasedResult::FAILED if the attack fails and
    findKey_LengthBasedResult::TIME_EXPIRED if the time limit is exceeded
  */
	virtual findKey_LengthBasedResult findKey_LengthBased( int N , 
							       const vector< Word >& A1 , 
							       const vector< Word >& A2 , 
							       const vector< Word >& B , 
							       int sec = 9999999, ostream& out = cout )=0;
	
};

//
//  LENGTH-BASED ATTACK ALGORITHM 1
//

//! Implements Length-based attack on AAG protocol.
/*!
  This is a very basic version of the length-based attack which implements the best descend LBA where
  on each step we choose conjugator which gives the maximal decrease among all
  currently available tuples. See A.Myasnikov. A.Ushakov, "On the length-based attack" for more details.
 */
class LengthAttack_A1 : public  LengthAttackBase
{
public:
	LengthAttack_A1(){}	
	int type() { return AL1; }
	//! Attack on an instances of the AAG protocol
	/*! Executes a length-based attack on an instance of the AAG protocol 
	  \param N  - rank of the braid group (number of strands)
	  \param A1 - Alices subgroup generators
	  \param A2 - Alices subgroup generators conjugated by Bob's private key
	  \param B  - Bob's subgroup (Bob's private key belongs to it)
	  \param sec - amount of time (in seconds) given to procedure to finish
	  \return  - Returns findKey_LengthBasedResult::SUCCESSFULL if the attack succeeds, findKey_LengthBasedResult::FAILED if the attack fails and
	  findKey_LengthBasedResult::TIME_EXPIRED if the time limit is exceeded
	*/
	findKey_LengthBasedResult findKey_LengthBased( int N , 
						       const vector< Word >& A1 , 
						       const vector< Word >& A2 , 
						       const vector< Word >& B , 
						       int sec = 9999999, ostream& out = cout );
private:												
	int  sbgpGeneratorsWeight( const vector< Word >& A );
	void addNewElt( const vector< Word >& A , set< ELT >& checkedElements , set< ELT >& uncheckedElements );
	void tryElt( int N , const ELT& cur , const vector< Word >& B , set< ELT >& checkedElements , set< ELT >& uncheckedElements );
	bool check_ifVectorsEqual( int N , const vector< Word >& A1 , const vector< Word >& A2 );	
};

//
//  LENGTH-BASED ATTACK ALGORITHM 2
//
//! Implements Length-based attack on AAG protocol.
/*!
  This is an implementation of the  generalised length-based attack.
  This is an LBA with backtracking in which   the set of elements in
  Alice~s public set is extended by all conjugators.  
  See A.Myasnikov. A.Ushakov, "On the length-based attack" for more details.
 */
class LengthAttack_A2 : public  LengthAttackBase
{
public:
	LengthAttack_A2(){}	
	int type() { return AL2; }
	//! Attack on an instances of the AAG protocol
	/*! Executes a length-based attack on an instance of the AAG protocol 
	  \param N  - rank of the braid group (number of strands)
	  \param A1 - Alices subgroup generators
	  \param A2 - Alices subgroup generators conjugated by Bob's private key
	  \param B  - Bob's subgroup (Bob's private key belongs to it)
	  \param sec - amount of time (in seconds) given to procedure to finish
	  \return  - Returns findKey_LengthBasedResult::SUCCESSFULL if the attack succeeds, findKey_LengthBasedResult::FAILED if the attack fails and
	  findKey_LengthBasedResult::TIME_EXPIRED if the time limit is exceeded
	*/
	findKey_LengthBasedResult findKey_LengthBased( int N , 
						       const vector< Word >& A1 , 
						       const vector< Word >& A2 , 
						       const vector< Word >& B , 
						       int sec = 9999999, ostream& out = cout );
private:												
	int  sbgpGeneratorsWeight( const vector< Word >& A );
	void addNewElt( const vector< Word >& A , set< ELT >& checkedElements , set< ELT >& uncheckedElements );
	void tryElt( int N , const ELT& cur , const vector< Word >& B , set< ELT >& checkedElements , set< ELT >& uncheckedElements );
	void tryElt( int N , const ELT& cur , const vector< Word >& B , set< ELT >& checkedElements , set< ELT >& uncheckedElements, ostream& out );
	bool check_ifVectorsEqual( int N , const vector< Word >& A1 , const vector< Word >& A2 );	
};

//
//  LENGTH-BASED ATTACK ALGORITHM 3
//
//! Implements Length-based attack on AAG protocol.
/*!
  This is an implementation of the  generalised length-based attack.
  This is an LBA with backtracking in which  the set of elements in
  Alice~s public set on each iteration is extended by conjugators and two-products 
  of the "best" generator.  
  See A.Myasnikov. A.Ushakov, "On the length-based attack" for more details.
 */

class LengthAttack_A3 : public  LengthAttackBase
{
public:
	LengthAttack_A3(){}	
	int type() { return AL3; }
	//! Attack on an instances of the AAG protocol
	/*! Executes a length-based attack on an instance of the AAG protocol 
	  \param N  - rank of the braid group (number of strands)
	  \param A1 - Alices subgroup generators
	  \param A2 - Alices subgroup generators conjugated by Bob's private key
	  \param B  - Bob's subgroup (Bob's private key belongs to it)
	  \param sec - amount of time (in seconds) given to procedure to finish
	  \return  - Returns findKey_LengthBasedResult::SUCCESSFULL if the attack succeeds, findKey_LengthBasedResult::FAILED if the attack fails and
	  findKey_LengthBasedResult::TIME_EXPIRED if the time limit is exceeded
	*/
	findKey_LengthBasedResult findKey_LengthBased( int N , 
						       const vector< Word >& A1 , 
						       const vector< Word >& A2 , 
						       const vector< Word >& B , 
						       int sec = 9999999, ostream& out = cout );
 private:		
	void addProducts(  const vector<Word>& elem_set, vector<Word>& ext_set, vector<Word>& ext_set_sg_gens, const Word& sel_gen, int sel_gen_sg );
	void addAllProducts(  const vector<Word>& elem_set, vector<Word>& ext_set, vector<Word>& ext_set_sg_gens );
	int  sbgpGeneratorsWeight( const vector< Word >& A );
	void addNewElt( const vector< Word >& A , set< ELT >& checkedElements , set< ELT >& uncheckedElements );
	void tryElt( int N , const ELT& cur , const vector< Word >& B , set< ELT >& checkedElements , set< ELT >& uncheckedElements );
	void tryElt( int N , const ELT& cur , const vector< Word >& B , const vector<Word>& B_sg_gens,set< ELT >& checkedElements , 
		     set< ELT >& uncheckedElements,
		     bool is_B_extended,
		     ostream& out );
	
	bool check_ifVectorsEqual( int N , const vector< Word >& A1 , const vector< Word >& A2 );	
};

#endif
