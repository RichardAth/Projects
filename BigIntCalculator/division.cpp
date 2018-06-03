/*
This file is part of Alpertron Calculators.
Copyright 2015 Dario Alejandro Alpern
Alpertron Calculators is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
Alpertron Calculators is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with Alpertron Calculators.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <iostream>
#include <cstdlib>
#include <string>
#include "bignbr.h"
//#include "expression.h"
#include <math.h>

static limb approxInv[MAX_LEN];
static limb adjustedArgument[MAX_LEN];
static limb arrAux[MAX_LEN];
static int bitLengthCycle[20];

// This routine uses Newton iteration: if x is an approximate inverse square root of N,
// a better approximation is: x(3-Nxx)/2. After the inverse square root is computed,
// the square root is found just by multiplying by N.
// The argument is multiplied by a power of 4 so the most significant limb is
// between LIMB_RANGE/4 and LIMB_RANGE - 1 and there is an even number of limbs.
// At the end of the calculation, the result is divided by the power of 2.

// All computations are done in little-endian notation.
// Find power of 2 that divides the number.
// output: pNbrLimbs = pointer to number of limbs
//         pPower2 = reference to power of 2.
static void MultiplyBigNbrByMinPowerOf2(int &pPower2, const limb *number, int len, limb *dest)
{
	limb mostSignficLimb, oldLimb, newLimb;
	int index2, mask, shLeft;
	limb *ptrDest;

	shLeft = 0;
	mostSignficLimb.x = (number + len - 1)->x;
	for (mask = LIMB_RANGE / 2; mask > 0; mask >>= 1)
	{
		if ((mostSignficLimb.x & mask) != 0)
		{
			break;
		}
		shLeft++;
	}
	ptrDest = dest;
	// Multiply number by this power.
	oldLimb.x = 0;
	for (index2 = len; index2 > 0; index2--)
	{
		newLimb.x = ptrDest->x;
		(ptrDest++)->x = ((newLimb.x << shLeft) |
			(oldLimb.x >> (BITS_PER_GROUP - shLeft))) & MAX_VALUE_LIMB;
		oldLimb.x = newLimb.x;
	}
	ptrDest->x = oldLimb.x >> (BITS_PER_GROUP - shLeft);
	pPower2 = shLeft;
}

/* used for operator overloading */
// After computing the number of limbs of the results, this routine finds the inverse
// of the divisor and then multiplies it by the dividend using nbrLimbs+1 limbs.
// After that, the quotient is adjusted.
BigInteger BigIntDivide(const BigInteger &Dividend, const BigInteger &Divisor) {
	BigInteger Quotient;
	double inverse;
	limb oldLimb, newLimb;
	int nbrLimbs, nbrLimbsDividend, nbrLimbsDivisor;

	// Check whether the divisor is zero.
	if (Divisor.limbs[0].x == 0 && Divisor.nbrLimbs == 1)
	{  // Indicate overflow if divisor is zero.
		std::string line = std::to_string(__LINE__);
		std::string mesg = "cannot divide by zero: ";
		mesg += __func__;
		mesg += " line ";  mesg += line;
		mesg += " in file "; mesg += __FILE__;
		throw std::range_error(mesg);
	}
	// Get number of limbs of quotient.
	nbrLimbsDividend = Dividend.nbrLimbs;
	nbrLimbsDivisor = Divisor.nbrLimbs;
	nbrLimbs = nbrLimbsDividend - nbrLimbsDivisor;
	if (nbrLimbs < 0)
	{   // Absolute value of dividend is less than absolute value of divisor.
		Quotient.limbs[0].x = 0;
		Quotient.nbrLimbs = 1;
		Quotient.sign = SIGN_POSITIVE;
		return Quotient;
	}
	if (nbrLimbs == 0)
	{   // Both divisor and dividend have the same number of limbs.
		for (nbrLimbs = nbrLimbsDividend - 1; nbrLimbs > 0; nbrLimbs--) {
			if (Dividend.limbs[nbrLimbs].x != Divisor.limbs[nbrLimbs].x) {
				break;
			}
		}
		if (Dividend.limbs[nbrLimbs].x < Divisor.limbs[nbrLimbs].x)
		{   // Dividend is less than divisor, so quotient is zero.
			Quotient.limbs[0].x = 0;
			Quotient.nbrLimbs = 1;
			Quotient.sign = SIGN_POSITIVE;
			return Quotient;
		}
	}
	if (nbrLimbsDividend == 1) {   // If dividend is small, perform the division directly.
		Quotient.limbs[0].x = Dividend.limbs[0].x / Divisor.limbs[0].x;
		Quotient.nbrLimbs = 1;
	}
	else if (nbrLimbsDivisor == 1)
	{   // Divisor is small: use divide by int.
		// Sign of quotient is determined later.
		Quotient = Dividend; //CopyBigInt(Quotient, Dividend);
		Quotient /= Divisor.limbs[0].x;   //subtractdivide(Quotient, 0, Divisor.limbs[0].x);
	}
	else {
		int index;
		int bitLength;
		int bitLengthNbrCycles;
		int idx;
		int nbrLimbsQuotient;
		int power2;
		limb *ptrDest, *ptrQuotient, *ptrQuot;
		const limb *ptrDivisor, *ptrDividend;

		nbrLimbs += 3;    // Use this number of limbs for intermediate calculations.
		if (nbrLimbs > nbrLimbsDivisor) {
			memset(&adjustedArgument[0], 0, (nbrLimbs - nbrLimbsDivisor) * sizeof(limb));
			memcpy(&adjustedArgument[nbrLimbs - nbrLimbsDivisor], &Divisor.limbs[0], nbrLimbsDivisor * sizeof(limb));
		}
		else {
			memcpy(&adjustedArgument[0], &Divisor.limbs[nbrLimbsDivisor - nbrLimbs], nbrLimbs * sizeof(limb));
		}
		MultiplyBigNbrByMinPowerOf2(power2, adjustedArgument, nbrLimbs, adjustedArgument);
		// Initialize approximate inverse.
		inverse = MAX_VALUE_LIMB / ((double)adjustedArgument[nbrLimbs - 1].x + 1);
		approxInv[nbrLimbs - 1].x = 1;
		if (inverse <= 1) {
			approxInv[nbrLimbs - 2].x = 0;
		}
		else {
			approxInv[nbrLimbs - 2].x = (int)floor((inverse - 1)*MAX_VALUE_LIMB);
		}
		// Perform Newton approximation loop.
		// Get bit length of each cycle.
		bitLengthNbrCycles = 0;
		bitLength = nbrLimbs*BITS_PER_GROUP;
		while (bitLength >= BITS_PER_GROUP) {
			bitLengthCycle[bitLengthNbrCycles++] = bitLength;
			bitLength = (bitLength + 1) >> 1;
		}
		// Each loop increments precision.
		// Use Newton iteration: x_{n+1} = x_n(2 - x_n)
		while (--bitLengthNbrCycles >= 0) {
			limb *ptrArrAux;
			int limbLength;

			bitLength = bitLengthCycle[bitLengthNbrCycles];
			limbLength = (bitLength + 3 * (BITS_PER_GROUP)-1) / BITS_PER_GROUP;
			if (limbLength > nbrLimbs) {
				limbLength = nbrLimbs;
			}
			// Compute x(2-Nx).
			// Multiply by divisor.
			multiply(&approxInv[nbrLimbs - limbLength],
				&adjustedArgument[nbrLimbs - limbLength], arrAux, limbLength, NULL);
			// Subtract arrAux from 2.
			ptrArrAux = &arrAux[limbLength];
			for (idx = limbLength - 1; idx > 0; idx--) {
				ptrArrAux->x = MAX_VALUE_LIMB - ptrArrAux->x;
				ptrArrAux++;
			}
			ptrArrAux->x = 1 - ptrArrAux->x;
			// Multiply arrAux by approxInv.
			multiply(&arrAux[limbLength], &approxInv[nbrLimbs - limbLength], approxInv, limbLength, NULL);
			memmove(&approxInv[nbrLimbs - limbLength], &approxInv[limbLength - 1], limbLength * sizeof(limb));
		}
		// Multiply approxInv by argument to obtain the quotient.
		if (nbrLimbsDividend >= nbrLimbs) {
			multiply(&Dividend.limbs[nbrLimbsDividend - nbrLimbs],
				approxInv, approxInv, nbrLimbs, NULL);
		}
		else {
			memset(arrAux, 0, (nbrLimbs - nbrLimbsDividend) * sizeof(limb));
			memcpy(&arrAux[nbrLimbs - nbrLimbsDividend], Dividend.limbs, nbrLimbsDividend * sizeof(limb));
			multiply(arrAux, approxInv, approxInv, nbrLimbs, NULL);
		}             // approxInv holds the quotient.
					  // Shift left quotient power2 bits into result.
		ptrDest = &approxInv[nbrLimbs - 1];
		oldLimb.x = 0;
		for (index = nbrLimbs; index >= 0; index--) {
			newLimb.x = ptrDest->x;
			(ptrDest++)->x = ((newLimb.x << power2) |
				(oldLimb.x >> (BITS_PER_GROUP - power2))) & MAX_VALUE_LIMB;
			oldLimb.x = newLimb.x;
		}

		// Determine number of limbs of quotient.
		nbrLimbsQuotient = nbrLimbsDividend - nbrLimbsDivisor;
		ptrDivisor = &Divisor.limbs[nbrLimbsDivisor - 1];
		ptrDividend = &Dividend.limbs[nbrLimbsDividend - 1];
		for (idx = nbrLimbsDivisor - 1; idx > 0; idx--) {
			if (ptrDividend->x != ptrDivisor->x) {
				break;
			}
			ptrDividend--;
			ptrDivisor--;
		}
		if (ptrDividend->x >= ptrDivisor->x) {
			nbrLimbsQuotient++;
		}
		ptrQuotient = &approxInv[2 * nbrLimbs - nbrLimbsQuotient];
		if (approxInv[2 * nbrLimbs - 1].x == 0)
		{  // Most significant byte is zero, so it is not part of the quotient. 
			ptrQuotient--;
		}
		ptrQuot = ptrQuotient;
		if ((ptrQuotient - 1)->x > (7 << (BITS_PER_GROUP - 3))) {    
			// Increment quotient.
			for (idx = 0; idx <= nbrLimbsQuotient; idx++) {
				if ((++((ptrQuotient + idx)->x)) & MAX_INT_NBR) {
					break;
				}
				(ptrQuotient + idx)->x = 0;
			}
			if (idx >= nbrLimbsQuotient)
			{                // Roll back on overflow.
				for (idx = 0; idx < nbrLimbsQuotient; idx++) {
					if (--((ptrQuotient + idx)->x) >= 0) {
						break;
					}
					(ptrQuotient + idx)->x = MAX_VALUE_LIMB;
				}
			}
			if (approxInv[2 * nbrLimbs - 1].x != 0)
			{    // Most significant byte is not zero, so it is part of the quotient.
				ptrQuot = &approxInv[2 * nbrLimbs - nbrLimbsQuotient];
			}
			// Test whether the quotient is correct.
			// It is correct only if multiplied by the divisor, it is <= than the dividend.
			if (nbrLimbsQuotient > nbrLimbsDivisor) {
				memcpy(&approxInv[0], Divisor.limbs, nbrLimbsDivisor * sizeof(limb));
				memset(&approxInv[nbrLimbsDivisor], 0, (nbrLimbsQuotient - nbrLimbsDivisor) * sizeof(limb));
				multiply(&approxInv[0], ptrQuot, arrAux, nbrLimbsQuotient, NULL);
			}
			else {
				memset(&approxInv[2 * nbrLimbs], 0, (nbrLimbsDivisor - nbrLimbsQuotient) * sizeof(limb));
				multiply(Divisor.limbs, ptrQuot, arrAux, nbrLimbsDivisor, NULL);
			}
			ptrDividend = (limb *)&Dividend.limbs[Dividend.nbrLimbs - 1];
			ptrDest = &arrAux[Dividend.nbrLimbs - 1];
			for (idx = Dividend.nbrLimbs - 1; idx > 0; idx--) {
				if (ptrDividend->x != ptrDest->x) {
					break;
				}
				ptrDividend--;
				ptrDest--;
			}
			if (ptrDividend->x < ptrDest->x) {  // Decrement quotient.
				ptrQuotient = ptrQuot;
				for (idx = 0; idx < nbrLimbsQuotient; idx++) {
					if (--(ptrQuotient->x) >= 0) {
						break;
					}
					(ptrQuotient++)->x = MAX_VALUE_LIMB;
				}
				if (idx == nbrLimbsQuotient) {
					nbrLimbsQuotient--;
				}
			}
		}
		memcpy(&Quotient.limbs[0], ptrQuot, nbrLimbsQuotient * sizeof(limb));
		Quotient.nbrLimbs = nbrLimbsQuotient;
	}
	if (Dividend.sign == Divisor.sign || (Quotient.limbs[0].x == 0 && Quotient.nbrLimbs == 1)) {
		Quotient.sign = SIGN_POSITIVE;
	}
	else {
		Quotient.sign = SIGN_NEGATIVE;
	}

	while (Quotient.nbrLimbs > 1) {
		if (Quotient.limbs[Quotient.nbrLimbs - 1].x == 0)
			Quotient.nbrLimbs--;
		else break;
	}
	return Quotient;
}

/* used for operator overloading */
BigInteger BigIntDivideInt(const BigInteger &Dividend, const int Divisor) {
	BigInteger Quotient;
	int len = Dividend.nbrLimbs;
	DivBigNbrByInt((int *)Dividend.limbs, abs(Divisor), (int *)Quotient.limbs, len);
	if (Divisor >= 0)
		if (Dividend.sign == SIGN_POSITIVE)
			Quotient.sign = SIGN_POSITIVE;
		else 
			Quotient.sign = SIGN_NEGATIVE;
	else   // Divisor is -ve
		if (Dividend.sign == SIGN_POSITIVE)
			Quotient.sign = SIGN_NEGATIVE;
		else
			Quotient.sign = SIGN_POSITIVE;

	while (len > 1 && Quotient.limbs[len-1].x == 0)
		len--;  // remove any leading zeros
	Quotient.nbrLimbs = len;
	return Quotient;
}