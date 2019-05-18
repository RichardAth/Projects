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

/* Perform Karatsuba multiplication in little endian order */

#include <string>
#include "bignbr.h"
//#include <stdio.h>
#include <cmath>
#include <cstdint>
#include <vector>
#include <iostream>

#define KARATSUBA_CUTOFF 16
//#define KARATSUBA_CUTOFF 2000    // turn of Karatsuba temporarily

static limb arr[4 * MAX_LEN];      /*    3 * changed to 4 * on 16/5/2019, because
							karatsuba exceeded arrray bound for large numbers */
static limb arrayAux[3 * MAX_LEN];
static int karatLength;
static void Karatsuba(int idxFactor1, const int length, int diffIndex, int level);

bool extralog = false;


/* this the entry point for external calls
result = factor1 * factor2 */
void multiply(const limb *factor1, const limb *factor2, limb *result,
	const int len, int *pResultLen)
{
	int length = len;
	int finalLength;
	// Compute length of numbers for each recursion.
	if (length > KARATSUBA_CUTOFF)
	{
		int div = 1;
		while (length > KARATSUBA_CUTOFF)
		{
			div *= 2;
			length = (length + 1) / 2;
		}
		length *= div;
	}
	if (length > MAX_LEN) {
		std::string line = std::to_string(__LINE__);
		std::string mesg = "number too big : cannot perform multiplication: ";
		mesg += __func__;
		mesg += " line ";  mesg += line;
		mesg += " in file "; mesg += __FILE__;
		throw std::range_error(mesg);
	}

	karatLength = length;
	memset(arr, 0, 2 * length * sizeof(limb));
	memcpy(&arr[0], factor1, len * sizeof(limb));
	memcpy(&arr[length], factor2, len * sizeof(limb));
	Karatsuba(0, length, 2 * length, 0);
	memcpy(result, &arr[2 * (karatLength - length)], 2 * length * sizeof(limb));
	if (pResultLen != NULL)
	{
		//memcpy(result, &arr[2 * (karatLength - length)], 2 * length * sizeof(limb));
		/*if (karatLength > length && arr[2 * (karatLength - length) - 1].x == 0)
		{
			*pResultLen = length * 2 - 1;
		}
		else
		{
			*pResultLen = length * 2;
		}*/
		for (finalLength = length * 2; finalLength > 1; finalLength--) {
			if (result[finalLength - 1].x != 0)
				break;
		}
		*pResultLen = finalLength;
	}
}

// The return value is the sign: true: negative.
// In result the absolute value of the difference is computed.
static int absSubtract(int idxMinuend, int idxSubtrahend,
	int idxResult, int nbrLen)
{
	int sign = 0;
	limb carry;
	int i;
	for (i = nbrLen - 1; i >= 0; i--)
	{
		if (arr[idxMinuend + i].x != arr[idxSubtrahend + i].x)
		{
			break;
		}
	}
	if (i >= 0 && arr[idxMinuend + i].x < arr[idxSubtrahend + i].x)
	{
		sign = 1;
		i = idxMinuend;    // Exchange minuend and subtrahend.
		idxMinuend = idxSubtrahend;
		idxSubtrahend = i;
	}
	carry.x = 0;
	for (i = 0; i < nbrLen; i++)
	{
		carry.x += arr[idxMinuend + i].x - arr[idxSubtrahend + i].x;
		arr[idxResult + i].x = carry.x & MAX_VALUE_LIMB;
		carry.x >>= BITS_PER_GROUP;
	}
	return sign;
}



/* multiply 2 numbers both stored in arr at idxFactor1 & 2. Result is
initially in arrayAux, then copied back into arr starting from idxFactor1. */
static void ClassicalMult(const int idxFactor1, const int idxFactor2, int nbrLen) {

	MultBigNbr((const long long *)arr+idxFactor1, (const long long *)arr+idxFactor2, 
		(long long *)arrayAux, nbrLen);
	memcpy(&arr[idxFactor1], &arrayAux[0], 2 * nbrLen * sizeof(limb));
	if (extralog)
		std::cout << "Classmult ix = " << idxFactor1 << " nbrLen = " << nbrLen << '\n';
	return;
}

#ifdef _DEBUG
/* check whether multiplication is correct. If used, multiplication is many times
slower but it has been very useful for debugging */
static bool Multcheck(std::vector<limb> f1, const int nbrLen, limb result[]) {
	std::vector <limb> newprod(nbrLen * 2 + 1);
	bool rv = true;
	newprod[nbrLen * 2].x = 0xdffddffddffddffd;  // 'canary' to detect buffer overrun

	/* redo the multiplication the slow way. */
	MultBigNbr((const long long *)f1.data(), (const long long *)f1.data() + nbrLen,
		(long long *)newprod.data(), nbrLen);
	/* compare the 2 results. Should be the same */
	for (int i = 0; i < nbrLen * 2; i++) {
		if (newprod[i].x != result[i].x) {
			std::cout << "multiplication error at index " << i
				<< "\n expected " << newprod[i].x << " got "
				<< result[i].x << '\n';
			rv = false;
		}
	}

	if (newprod[nbrLen * 2].x != 0xdffddffddffddffd) {
		std::cout << "buffer overrun on multiplication\n";
		rv = false;
	}

	return rv;
}
#endif

// Recursive Karatsuba function.
static void Karatsuba(int idxFactor1, const int nbrLen, int diffIndex, int level) {
	/* check that array index is within bounds */
	assert(idxFactor1 + nbrLen * 2 <= sizeof(arr) / sizeof(arr[0]));
#ifdef _DEBUG
	/* save original value. Use a vector so that recursion does't cause
	a stack overflow */
	/*std::vector <limb> Fsave(nbrLen*2);
	memcpy(Fsave.data(), arr + idxFactor1, nbrLen*2 * sizeof(limb));*/
#endif
	int idxFactor2 = idxFactor1 + nbrLen;
	int i;
	unsigned long long carry1First, carry1Second;
	unsigned long long carry2Second;
	limb *ptrResult, *ptrHigh, tmp;
	int middle;
	int sign;
	int halfLength;

	// Check if one of the factors is equal to zero.
	ptrResult = &arr[idxFactor1];
	for (i = nbrLen; i > 0; i--) {
		if ((ptrResult++)->x != 0) {
			break;
		}
	}

	if (i > 0) {     // First factor is not zero. Check second.
		ptrResult = &arr[idxFactor2];
		for (i = nbrLen; i > 0; i--) {
			if ((ptrResult++)->x != 0) 	{
				break;
			}
		}
	}

	if (i == 0) {    // One of the factors is equal to zero.
		for (i = nbrLen - 1; i >= 0; i--) {
			arr[idxFactor1 + i].x = arr[idxFactor2 + i].x = 0;
		}
		if (extralog)
			std::cout << "karatmult ix = " << idxFactor1 << " nbrLen = " << nbrLen
			<< " level = " << level << " (return 0) \n";
		return;
	}

	if (nbrLen <= KARATSUBA_CUTOFF) {
		// Below cutoff: perform standard classical multiplcation.
		ClassicalMult(idxFactor1, idxFactor2, nbrLen);
		return;
	}

	// Length > KARATSUBA_CUTOFF: Use Karatsuba multiplication.
	// It uses three half-length multiplications instead of four.
	//  x*y = (xH*b + xL)*(yH*b + yL)
	//  x*y = (b + 1)*(xH*yH*b + xL*yL) + (xH - xL)*(yL - yH)*b
	// The length of b is stored in variable halfLength.
	// Since the absolute values of (xH - xL) and (yL - yH) fit in
	// a single limb, there will be no overflow.

	// At this moment the order is: xL, xH, yL, yH.
	// Exchange high part of first factor with low part of 2nd factor.

	halfLength = nbrLen >> 1;    // nbrLen/2, rounded down
	for (i = idxFactor1 + halfLength; i<idxFactor2; i++) {
		tmp.x = arr[i].x;
		arr[i].x = arr[i + halfLength].x;
		arr[i + halfLength].x = tmp.x;
	}

	// At this moment the order is: xL, yL, xH, yH.
	// Get absolute values of (xH-xL) and (yL-yH) and the signs.
	sign = absSubtract(idxFactor1, idxFactor2, diffIndex, halfLength);
	sign ^= absSubtract(idxFactor2 + halfLength, idxFactor1 + halfLength,
		diffIndex + halfLength, halfLength);
	middle = diffIndex;
	diffIndex += nbrLen;
	Karatsuba(idxFactor1, halfLength, diffIndex, level+1); // Multiply both low parts.
	Karatsuba(idxFactor2, halfLength, diffIndex, level+1); // Multiply both high parts.
	Karatsuba(middle, halfLength, diffIndex, level+1);     // Multiply the differences.
					// Process all carries at the end.
					// Obtain (b+1)(xH*yH*b + xL*yL) = xH*yH*b^2 + (xL*yL+xH*yH)*b + xL*yL
					// The first and last terms are already in correct locations.
	ptrResult = &arr[idxFactor1 + halfLength];
	carry1First = carry1Second = carry2Second = 0;

	for (i = halfLength; i > 0; i--) {
		// The sum of three ints overflows an unsigned int variable,
		// so two adds are required. Also carries must be separated in
		// order to avoid overflow:
		// 00000001 + 7FFFFFFF + 7FFFFFFF = FFFFFFFF
		unsigned long long accum1Lo = carry1First + ptrResult->x + (ptrResult + halfLength)->x;
		unsigned long long accum2Lo;
		carry1First = accum1Lo >> BITS_PER_GROUP;
		accum2Lo = carry2Second + (accum1Lo & MAX_VALUE_LIMB) +
			(ptrResult - halfLength)->x;
		carry2Second = accum2Lo >> BITS_PER_GROUP;
		accum1Lo = carry1Second + (accum1Lo & MAX_VALUE_LIMB) +
			(ptrResult + nbrLen)->x;
		carry1Second = accum1Lo >> BITS_PER_GROUP;
		(ptrResult + halfLength)->x = accum1Lo & MAX_VALUE_LIMB;
		ptrResult->x = accum2Lo & MAX_VALUE_LIMB;
		ptrResult++;
	}

	(ptrResult + halfLength)->x += carry1First + carry1Second;
	ptrResult->x += carry1First + carry2Second;
	// Process carries.
	ptrResult = &arr[idxFactor1];
	carry1First = 0;
	for (i = 2 * nbrLen; i > 0; i--) {
		carry1First += ptrResult->x;
		(ptrResult++)->x = carry1First & MAX_VALUE_LIMB;
		carry1First >>= BITS_PER_GROUP;
	}

	// Compute final product.
	ptrHigh = &arr[middle];
	ptrResult = &arr[idxFactor1 + halfLength];

	if (sign != 0) {            // (xH-xL) * (yL-yH) is negative.
		long long borrow = 0;

		for (i = nbrLen; i > 0; i--) {
			borrow += ptrResult->x - (ptrHigh++)->x;
			(ptrResult++)->x = borrow & MAX_VALUE_LIMB;
			borrow >>= BITS_PER_GROUP;     /* note: For microsoft C/C++ right shift of a 
	        -ve number is eqivalent to dividing it by 2 but it is rounded DOWN towards 
			-infinity. In this case if borrow is -ve the result is -1, if positive the 
			result is 0 */
		}

		for (i = halfLength; i > 0; i--) {
			borrow += ptrResult->x;
			(ptrResult++)->x = borrow & MAX_VALUE_LIMB;
			borrow >>= BITS_PER_GROUP;
		}
	}

	else {            // (xH-xL) * (yL-yH) is positive or zero.
		unsigned long long carry = 0;
		for (i = nbrLen; i > 0; i--) {
			carry += (unsigned long long)ptrResult->x + (unsigned long long)(ptrHigh++)->x;
			(ptrResult++)->x = (long long)(carry & MAX_VALUE_LIMB);
			carry >>= BITS_PER_GROUP;
		}

		for (i = halfLength; i > 0; i--) {
			carry += (unsigned long long)ptrResult->x;
			(ptrResult++)->x = (long long)(carry & MAX_VALUE_LIMB);
			carry >>= BITS_PER_GROUP;
		}
	}
	if (extralog)
		std::cout << "karatmult ix = " << idxFactor1 << " nbrLen = " << nbrLen 
			<< " level = " << level << '\n';
#ifdef _DEBUG
	/*if (!Multcheck(Fsave, nbrLen, arr + idxFactor1)) {
		std::cout << "karatsuba multiplication error. F1 index = " << idxFactor1
			<< " nbrLen = " << nbrLen << '\n';
		abort();
	}*/
#endif
}
