// Copyright (C) 2007 Alexander Ushakov
// Contents: Definition of "Random FSA Generation" methods.
//
// Principal Authors: Alexander Ushakov
//
// Revision History:
//

#ifndef _RandomFSA_h_
#define _RandomFSA_h_

#include "FSA.h"


//---------------------------------------------------------------------------//
//------------------------------ RandomFSA ----------------------------------//
//---------------------------------------------------------------------------//


//! Generate a random (connected admissible) Finite State Automaton.
/*
  The procedure generates a random FSA. Moreover, the corresponding distribution
  tends to uniform as the size tends to infinity. \f$N\f$ is the number of vertices 
  in the FSA and \f$L\f$ is the size of the alphabet.
  
  Implementation is based on Bassino, Nicaud, Weil, "Random generation of finitely 
  generated subgroups of a free group".
*/
FSA randomFSA( int N , int L );


#endif
