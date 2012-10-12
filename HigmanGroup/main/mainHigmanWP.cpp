/*
 * main.cpp
 *
 *  Created on: 23.03.2011
 *      Author: Juern
 */

#include <string>
#include <cstring>
#include <set>
#include <vector>
#include <iostream>
#include <iomanip>
#include <assert.h>
#include <math.h>
#include <list>
#include <stdlib.h>
#include <time.h>
#include "Sign.h"
#include "SignMatrix.h"
#include "PowerCircuit.h"
#include "PowerCircuitGraph.h"
#include "PowerCircuitCompMatrix.h"
#include "BaumslagGersten.h"


using namespace std;
using namespace PC;
using namespace BG;

//-------------------------------------------
// switch between different implementations

//#define PC PowerCircuitCompMatrix
#define PC PowerCircuitGraph

int main()
{
	try
	{
		//-------------------------------------------
		// Some other small instances of the word
		// problem in BG(1,2)
		PC	pc1;
		cout << "Is taTAA=1 in BG(1,2)? " << (solveWPinBG(&pc1, "taTAA") ? "yes" : "no") << endl;

		PC	pc2;
		cout << "Is baB a bAB AA=1 in BG(1,2)? " << (solveWPinBG(&pc2, "baB a bAB AA") ? "yes" : "no") << endl;

		PC	pc3;
		cout << "Is A baB aa bAB AAA=1 in BG(1,2)? " << (solveWPinBG(&pc3, "A baB aa bAB AAA") ? "yes" : "no") << endl;

		PC	pc4;
		cout << "Is A baB aa bAB AA=1 in BG(1,2)? " << (solveWPinBG(&pc4, "A baB aa bAB AA") ? "yes" : "no") << endl;

		//-------------------------------------------
		// Now test whether [t, T(K)]=1, where T(K)
		// is a word which equals t^tower(K) in
		// BG(1,2) and which is build up via
		// T(0)=t, T(K+1)=bT(K)aT(K)^(-1)b^(-1).
		const int K = 7;
		string TK = "t", TK_inv = "T"; // T(0)=t
		for(int i = 1; i <= K; i++)
		{
			string w = "b" + TK + "a" + TK_inv + "B";
			string w_inv = "b" + TK + "A" + TK_inv + "B";
			TK = w; TK_inv = w_inv;
		}
		string comm_t_TK = "t(" + TK + ")T(" + TK_inv + ")";
		PC pc5;
		cout << "Is [t,T(" << K << ")]=" << comm_t_TK << "=1 in BG(1,2)? " << (solveWPinBG(&pc5, comm_t_TK) ? "yes" : "no") << endl;

		//-------------------------------------------
		// Finally, check whether [T(K), T(L)]=1
		const int L = 8;
		string TL = "t", TL_inv = "T";
		for(int i = 1; i <= L; i++)
		{
			string w = "b" + TL + "a" + TL_inv + "B";
			string w_inv = "b" + TL + "A" + TL_inv + "B";
			TL = w; TL_inv = w_inv;
		}
		string comm_TK_TL = "(" + TK + ")(" + TL + ")(" + TK_inv + ")(" + TL_inv + ")";
		PC pc6;
		cout << "Is [T(" << K << "),T(" << L << ")]=" << comm_TK_TL << "=1 in BG(1,2)? " << (solveWPinBG(&pc6, comm_TK_TL) ? "yes" : "no") << endl;
	}
	catch(char const* message)
	{
		cout << "Caught Exception: " << message << endl;
	}

	cout << "Hello Circuits!" << endl; // This indicates that the program has terminated without crashing.
	return 0;
}

