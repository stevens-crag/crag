/*
 * Sign.cpp
 *
 *  Created on: 06.04.2011
 *      Author: Juern
 */

#include <string>
#include <cstring>
#include <vector>
#include <iostream>
#include <iomanip>
#include <assert.h>
#include <math.h>
#include "Sign.h"

namespace PC
{

void addSigns(const Sign& op1, const Sign& op2, Sign& result, Sign& carry)
{
	if(op1 == ZERO || op2 == ZERO)
	{
		result = op1 | op2;
		carry = ZERO;
	}
	else if(op1 != op2)
	{
		result = ZERO;
		carry = ZERO;
	}
	else
	{
		result = ZERO;
		carry = op1;
	}
}

int compareSigns(const Sign& op1, const Sign& op2)
{
	if(op1 == op2) return 0;
	else if(op1 == MINUS) return -1;
	else if(op1 == PLUS) return 1;
	else return (op2 == PLUS ? -1 : 1);
}

int signToInt(const Sign& s)
{
	if(s==ZERO)
		return 0;
	else if(s==PLUS)
		return 1;
	else
		return -1;
}

Sign negateSign(const Sign& s)
{
	return (s == ZERO ? ZERO : (s == PLUS ? MINUS : PLUS));
}

std::string	signToString(const Sign& s)
{
	return (s == ZERO ? " 0" : (s == PLUS ? "+1" : "-1"));
}

}
