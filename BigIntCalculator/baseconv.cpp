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


#include "pch.h"
#pragma warning(disable : 4996)

#include "bignbr.h"
#include "bigint.h"

#define DIGITS_PER_LIMB 9
#define MAX_LIMB_CONVERSION 1000000000
#define FIRST_MULT  (1 << (BITS_PER_GROUP/2))
#define SECOND_MULT (LIMB_RANGE / FIRST_MULT)

static limb power10000[MAX_LEN];

/* equivalent to:
_ui64toa((unsigned long long)nbr, *pOutput, 10));
*pOutput += (strlen(*pOutput)); */
void int2dec(char **pOutput, long long nbr)
{
	char *ptrOutput = *pOutput;
	int significantZero = 0;
	unsigned long long div = 1000000000000000000;
	unsigned long long value = (unsigned long long)nbr;
	while (div > 0)
	{
		long long digit;

		digit = value / div;
		if (digit > 0 || significantZero != 0)
		{
			significantZero = 1;
			*ptrOutput++ = (char)(digit + (int)'0');
		}
		value %= div;
		div /= 10;
	}
	if (significantZero == 0)
	{
		*ptrOutput++ = '0';
	}
	*ptrOutput = '\0';   // null terminate O/P
	*pOutput = ptrOutput;
}

// Convert little-endian number to a string with space every groupLen digits.
// In order to perform a faster conversion, use groups of DIGITS_PER_LIMB digits.
// output to char array decimal
void Bin2Dec(const limb binary[], char* decimal, int nbrLimbs, int groupLength)
{
	int len, index, index2, count;
	const limb *ptrSrc = binary + nbrLimbs - 1;
	char *ptrDest;
	int significantZero = 0;
	int groupCtr;
	int digit[DIGITS_PER_LIMB];
	int digits = 0;
	bool showDigitsText = true;

	if (groupLength <= 0)
	{
		groupLength = -groupLength;
		showDigitsText = false;
	}
	power10000[0].x = ptrSrc->x % MAX_LIMB_CONVERSION;
	power10000[1].x = ptrSrc->x / MAX_LIMB_CONVERSION;
	len = (power10000[1].x == 0 ? 1 : 2); // Initialize array length.
	for (index = nbrLimbs - 2; index >= 0; index--)
	{
		double dCarry, dQuotient;
		limb *ptrPower;

		// Multiply by FIRST_MULT and then by SECOND_MULT, so there is never
		// more than 53 bits in the product.

		ptrPower = power10000;
		dQuotient = 0;
		for (index2 = 0; index2 < len; index2++)
		{
			dCarry = dQuotient + (double)ptrPower->x * FIRST_MULT;
			dQuotient = floor(dCarry / MAX_LIMB_CONVERSION);
			(ptrPower++)->x = (int)(dCarry - dQuotient * MAX_LIMB_CONVERSION);
		}
		if (dQuotient != 0)
		{
			(ptrPower++)->x = (int)dQuotient;
			len++;
		}
		ptrPower = power10000;
		dQuotient = (double)(--ptrSrc)->x;
		for (index2 = 0; index2 < len; index2++)
		{
			dCarry = dQuotient + (double)ptrPower->x * SECOND_MULT;
			dQuotient = floor(dCarry / MAX_LIMB_CONVERSION);
			(ptrPower++)->x = (int)(dCarry - dQuotient * MAX_LIMB_CONVERSION);
		}
		if (dQuotient != 0)
		{
			(ptrPower++)->x = (int)dQuotient;
			len++;
		}
	}
	// At this moment the array power10000 has the representation
	// of the number in base 10000 in little-endian. Convert to
	// ASCII separating every groupLength characters.
	ptrDest = decimal;
	ptrSrc = &power10000[len - 1];
	groupCtr = len * DIGITS_PER_LIMB;
	if (groupLength != 0)
	{
		groupCtr %= groupLength;
		if (groupCtr == 0)
		{
			groupCtr = groupLength;
		}
	}
	for (index = len; index > 0; index--)
	{
		int value = (int)(ptrSrc--)->x;
		for (count = 0; count < DIGITS_PER_LIMB; count++)
		{
			digit[count] = value % 10;
			value /= 10;
		}
		for (count = DIGITS_PER_LIMB - 1; count >= 0; count--)
		{
			if (digit[count] != 0 || significantZero != 0)
			{
				digits++;
				*ptrDest++ = (char)(digit[count] + (int)'0');
				if (groupCtr == 1)
				{
					*ptrDest++ = ' ';
				}
				significantZero = 1;
			}
			if (--groupCtr == 0)
			{
				groupCtr = groupLength;
			}
		}
	}
	if (significantZero == 0)
	{     // Number is zero.
		*ptrDest++ = '0';
		*ptrDest = '\0';
		return;
	}
	if (digits > 30 && showDigitsText)
	{
		*ptrDest++ = '(';
		int2dec(&ptrDest, digits);
		strcpy(ptrDest, " digits)");
		ptrDest += strlen(ptrDest);
	}
	else if (ptrDest > decimal)
	{
		*(ptrDest - 1) = '\0';       // Add terminator.
	}
}

/* convert BigInteger to a text string */
void Bin2Dec(const BigInteger &BigInt, char *decimal, int groupLength) {
	if (BigInt.sign == SIGN_NEGATIVE) {
		*decimal++ = '-';
	}
	Bin2Dec(BigInt.limbs, decimal, BigInt.nbrLimbs, groupLength);
}

/* << definitions for output streams (std::cout etc)
definitions below are crude but they work. Only decimal output
is supported for direct output */
std::ostream  &operator<<(std::ostream  &s, const BigInteger  &n) {
	static std::string  buffer;   /* using std::string takes care of memory 
								  management without risking memory leaks */
	buffer.clear();
	int length = n.nbrLimbs * BITS_PER_GROUP; // get length in bits
	double buflen = (double)length*log10(2);  // get length in decimal digits
	buflen = buflen * 7 / 6 + 5; // allow for space every 6 digits + sign + trailing null
	buffer.resize((size_t)buflen);  // N.B. buffer larger than needed is OK, but 
								   // too small could be disastrous
	Bin2Dec(n, &buffer[0], -6);
	auto sLen = strnlen_s(&buffer[0], (size_t)buflen);
	buffer.resize(sLen);
	return s << buffer;
}
