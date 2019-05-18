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
#include <math.h>

static limb approxInv[MAX_LEN];
static limb adjustedArgument[MAX_LEN];
static limb arrAux[MAX_LEN];
static int bitLengthCycle[20];

// All computations are done in little-endian notation.
// Find power of 2 that divides the number.
// output: pNbrLimbs = pointer to number of limbs
//         pPower2 = reference to power of 2.
static void MultiplyBigNbrByMinPowerOf2(int &pPower2, const limb *number, int len, limb *dest)
{
	limb mostSignficLimb, oldLimb, newLimb;
	int index2, shLeft;
	limb *ptrDest;
	unsigned long long mask;

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
BigInteger BigIntDivide(const BigInteger &Dividend, const BigInteger &Divisor,
	BigInteger &Remainder) {
	BigInteger Quotient;
	double inverse;
	limb oldLimb, newLimb;
	int nbrLimbs, nbrLimbsDividend, nbrLimbsDivisor;
	long long rem;
	int quotadjust = 0;

	// Check whether the divisor is zero.
	if (Divisor == 0) {  // Indicate overflow if divisor is zero.
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
		Quotient = 0;
		Remainder = Dividend;
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
		{   // abs(Dividend) is less than abs(divisor), so quotient is zero.
			Quotient = 0;
			Remainder = Dividend;
			return Quotient;
		}
	}

	if (nbrLimbsDividend == 1) {   /* If dividend is small (and divisor is smaller), 
								   perform the division directly. */
		Quotient.limbs[0].x = Dividend.limbs[0].x / Divisor.limbs[0].x;
		Quotient.nbrLimbs = 1;
	}

	else if (nbrLimbsDivisor == 1)
	{   // Divisor is small: use divide by int.
		// Sign of quotient is determined later.
		Quotient = BigIntDivideInt(Dividend, Divisor.limbs[0].x, rem);
		Remainder = rem;   
	}
	else {
		/* code below becomes inaccurate at about Quotient size = 63 limbs 
		(about 1200 digits) */
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
			memcpy(&adjustedArgument[nbrLimbs - nbrLimbsDivisor], &Divisor.limbs[0], 
				nbrLimbsDivisor * sizeof(limb));
		}
		else {
			memcpy(&adjustedArgument[0], &Divisor.limbs[nbrLimbsDivisor - nbrLimbs], 
				nbrLimbs * sizeof(limb));
		}

		MultiplyBigNbrByMinPowerOf2(power2, adjustedArgument, nbrLimbs, adjustedArgument);
		// Initialize approximate inverse.
		inverse = MAX_VALUE_LIMB / ((double)adjustedArgument[nbrLimbs - 1].x + 1);
		approxInv[nbrLimbs - 1].x = 1;
		if (inverse <= 1) {
			approxInv[nbrLimbs - 2].x = 0;
		}
		else {
			approxInv[nbrLimbs - 2].x = (long long)floor((inverse - 1)*MAX_VALUE_LIMB);
		}
		// Perform Newton approximation loop.
		// Get bit length of each cycle.
		bitLengthNbrCycles = 0;
		bitLength = nbrLimbs*BITS_PER_GROUP;

		/* if lower limit below for bitLength is too high, the inverse will not be 
		accurate enough, and large quotients will be incorrect */
		while (bitLength >= BITS_PER_GROUP/2-5) {
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
			multiply(&arrAux[limbLength], &approxInv[nbrLimbs - limbLength], 
				approxInv, limbLength, NULL);
			memmove(&approxInv[nbrLimbs - limbLength], &approxInv[limbLength - 1], 
				limbLength * sizeof(limb));
		}

		// Multiply approxInv by argument to obtain the quotient.
		if (nbrLimbsDividend >= nbrLimbs) {
			multiply(&Dividend.limbs[nbrLimbsDividend - nbrLimbs],
				approxInv, approxInv, nbrLimbs, NULL);
		}
		else {
			memset(arrAux, 0, (nbrLimbs - nbrLimbsDividend) * sizeof(limb));
			memcpy(&arrAux[nbrLimbs - nbrLimbsDividend], Dividend.limbs, 
				nbrLimbsDividend * sizeof(limb));
			multiply(arrAux, approxInv, approxInv, nbrLimbs, NULL);
		} 

		// approxInv holds the quotient.
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
		if ((ptrQuotient - 1)->x > (7LL << (BITS_PER_GROUP - 3))) {    
			// Increment quotient.
			for (idx = 0; idx <= nbrLimbsQuotient; idx++) {
				if ((++((ptrQuotient + idx)->x)) & MAX_INT_NBR) {
					break;
				}
				(ptrQuotient + idx)->x = 0;
			}

			if (idx >= nbrLimbsQuotient) {
				// Roll back on overflow.
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
			//if (nbrLimbsQuotient > nbrLimbsDivisor) {
			//	memcpy(&approxInv[0], Divisor.limbs, nbrLimbsDivisor * sizeof(limb));
			//	memset(&approxInv[nbrLimbsDivisor], 0, (nbrLimbsQuotient - nbrLimbsDivisor) * sizeof(limb));
			//	multiply(&approxInv[0], ptrQuot, arrAux, nbrLimbsQuotient, NULL);
			//}
			//else {
			//	memset(&approxInv[2 * nbrLimbs], 0, (nbrLimbsDivisor - nbrLimbsQuotient) * sizeof(limb));
			//	/* arrAux = Divisor * Quot */
			//	multiply(Divisor.limbs, ptrQuot, arrAux, nbrLimbsDivisor, NULL);
			//}
			//ptrDividend = (limb *)&Dividend.limbs[Dividend.nbrLimbs - 1];
			//ptrDest = &arrAux[Dividend.nbrLimbs - 1];
			//for (idx = Dividend.nbrLimbs - 1; idx > 0; idx--) {
			//	if (ptrDividend->x != ptrDest->x) {
			//		break;
			//	}
			//	ptrDividend--;
			//	ptrDest--;
			//}

			//if (ptrDividend->x < ptrDest->x) {  // Decrement quotient.
			//	quotadjust--;
			//	ptrQuotient = ptrQuot;
			//	for (idx = 0; idx < nbrLimbsQuotient; idx++) {
			//		if (--(ptrQuotient->x) >= 0) {
			//			break;
			//		}
			//		(ptrQuotient++)->x = MAX_VALUE_LIMB;
			//	}
			//	if (idx == nbrLimbsQuotient) {
			//		nbrLimbsQuotient--;
			//	}
			//}
		}

		memcpy(&Quotient.limbs[0], ptrQuot, nbrLimbsQuotient * sizeof(limb));
		Quotient.nbrLimbs = nbrLimbsQuotient;
	}

	/* now set sign of quotient */
	if (Dividend.sign == Divisor.sign || (Quotient.limbs[0].x == 0 && Quotient.nbrLimbs == 1)) {
		/* either dividend & divisor have same sign, or quotient is zero*/
		Quotient.sign = SIGN_POSITIVE;  
	}
	else {
		Quotient.sign = SIGN_NEGATIVE;
	}

	/* adjust length (remove leading zeros) */
	while (Quotient.nbrLimbs > 1) {
		if (Quotient.limbs[Quotient.nbrLimbs - 1].x == 0)
			Quotient.nbrLimbs--;
		else break;
	}

	/* check quotient by calculating remainder, then checking that the remainder
	is in the correct range. It is normal that the magnitude of the quotient is
	sometimes reduced by 1. Any other adjustment indicates a problem. */
	BigInteger d2;
	int adjustctr = 0;
	do {
		d2 = Quotient * Divisor;
		Remainder = Dividend - d2;
		if (Remainder.sign != Dividend.sign) {
			/* need to decrease magnitude of Quotient */
			if (Quotient >= 0)
				Quotient--;
			else
				Quotient++;
			adjustctr--;
		} else
			if (abs(Remainder) > abs(Divisor)) {
				/* need to increase magnitude of Quotient */
				if (Quotient >= 0)
					Quotient++;
				else
					Quotient--;
				adjustctr++;
			}
			else break;     // Remainder is in correct range; quotient is correct
	} while (true);

#ifdef _DEBUG
	if (adjustctr < -1 || adjustctr > 0)
			std::cout << "Division: quotient adjusted by " << adjustctr 
				<< "\nDividend size = " << Dividend.nbrLimbs 
				<< " limbs. Divisor size = " << Divisor.nbrLimbs 
				<< " limbs. Quotient size = " << Quotient.nbrLimbs << '\n';

	/* redo division using Znums to check result */
	//Znum Zdividend, Zdivisor, Zquotient, Zq2, Zremainder, Zr2;
	//BigtoZ(Zdividend, Dividend);
	//BigtoZ(Zdivisor, Divisor);
	//Zquotient = Zdividend / Zdivisor;
	//Zremainder = Zdividend % Zdivisor;
	//BigtoZ(Zq2, Quotient);
	//BigtoZ(Zr2, Remainder);
	//if (Zq2 != Zquotient) {
	//	std::cout << "Dividend  = " << Zdividend
	//		<< "\nDivisor   = " << Zdivisor
	//		<< "\nQuotient  = " << Zquotient
	//		//<< "\nZQ2       = " << Zq2
	//		<< "\nDiff.     = " << (Zquotient - Zq2)
	//		//<< "\nRemainder = " << Zremainder
	//		//<< "\nZr2       = " << Zr2 
	//		<< "\nDividend size = " << Dividend.nbrLimbs << " limbs. Divisor size = "
	//		<< Divisor.nbrLimbs << " limbs. Quotient size = "
	//		<< Quotient.nbrLimbs << '\n';
	//}
#endif
	return Quotient;
}

/* used for operator overloading. Value returned is the quotient. The remainder 
is also calculated. */
BigInteger BigIntDivideInt(const BigInteger &Dividend, const long long Divisor, long long &rem) {
	BigInteger Quotient;
	int len = Dividend.nbrLimbs;
	//BigInteger remainder;

	DivBigNbrByInt((long long *)Dividend.limbs, abs(Divisor), (long long *)Quotient.limbs, 
		len, rem);
	// rem is +ve, fix sign 
	if (Dividend.sign == SIGN_NEGATIVE)
		rem = -rem;  // remainder always has same sign as dividend for truncation division
	if (Divisor >= 0)
		Quotient.sign = Dividend.sign;
	else   // Divisor is -ve
		if (Dividend.sign == SIGN_POSITIVE) {
			Quotient.sign = SIGN_NEGATIVE;
		}
		else {
			Quotient.sign = SIGN_POSITIVE;
		}

	while (len > 1 && Quotient.limbs[len-1].x == 0)
		len--;             // remove any leading zeros
	Quotient.nbrLimbs = len;

	return Quotient;
}

/* calculate Dividend mod Divisor (used for operator overloading)  */
BigInteger BigIntRemainder(const BigInteger &Dividend, const BigInteger &Divisor) {
	BigInteger quot, remainder, d2;
	quot = BigIntDivide(Dividend, Divisor, remainder);
	
	if (remainder.sign != Dividend.sign) {
		Znum Zdividend, Zdivisor, Zquot, Zremainder, Zq2, Zd2;
		BigtoZ(Zdivisor, Divisor);
		BigtoZ(Zdividend, Dividend);
		BigtoZ(Zquot, quot);
		BigtoZ(Zremainder, remainder);
		BigtoZ(Zd2, d2);
		Zq2 = Zdividend / Zdivisor;
		std::cout << "Remainder error quotient = " << Zquot<< '\n';
		if (Zquot != Zq2)
			std::cout << "quot diff = " << Zquot - Zq2 << '\n';
		if (Zd2 != (Zquot * Zdivisor)) {
			Znum diff = Zd2 ^ (Zquot *Zdivisor); // compare using XOR
			/*std::cout << "Zd2 = " << Zd2
				<< "\nError = " << Zd2 - (Zquot *Zdivisor) << '\n';*/
			gmp_printf("Error = %#Zx \n", diff);
		} 
		else {
			std::cout << "got remainder " << Zremainder << '\n';
			std::cout << "should be     " << Zdividend % Zdivisor << '\n';
		}

	}

	return remainder;
}