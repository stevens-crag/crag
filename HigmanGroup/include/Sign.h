/*
 * Sign.h
 *
 *  Created on: 06.04.2011
 *      Author: Juern
 */

#ifndef SIGN_H_
#define SIGN_H_

namespace PC
{

/**
 * Type for storing signs, i.e., members of the set
 * {-1,0,+1}.
 */

typedef char Sign;
const Sign ZERO = 0/*0b00*/;
const Sign PLUS = 1/*0b01*/;
const Sign MINUS = 2/*0b10*/;

void addSigns(const Sign& op1, const Sign& op2, Sign& result, Sign& carry);
int compareSigns(const Sign& op1, const Sign& op2);
int signToInt(const Sign& s);
Sign negateSign(const Sign& s);
std::string	signToString(const Sign& s);

}

#endif /* SIGN_H_ */
