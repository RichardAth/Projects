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

#define _CRT_SECURE_NO_DEPRECATE

#include <string>

#include <stdio.h>
#include <math.h>

#include "bignbr.h"

#define DIGITS_PER_LIMB 9
#define MAX_LIMB_CONVERSION 1000000000
#define FIRST_MULT  (1 << (BITS_PER_GROUP/2))
#define SECOND_MULT (LIMB_RANGE / FIRST_MULT)

extern int lang;

static limb power10000[MAX_LEN];
static limb temp[MAX_LEN];
//static void add(limb *addend1, limb *addend2, limb *sum, int length);

/* Convert number to little-endian.
decimal: pointer  to input in ascii chars base 10,
binary:  pointer to output output
digits:  number of digits in input */
//void Dec2Bin(const char *decimal, limb *binary, int digits, int *bitGroups)
//{
//	// First step: separate in groups of DIGITS_PER_LIMB digits.
//	const char *ptrSrc;
//	limb *ptrDest, *ptrBinary;
//	int digit, groupNbr, multiplier;
//	int outerGroup, innerGroup;
//	int nbrGroups = 1;
//	for (nbrGroups = 1; nbrGroups*DIGITS_PER_LIMB < digits; nbrGroups *= 2)
//	{  // this for loop does nothing!!
//	}
//	memset(binary, 0, nbrGroups * sizeof(limb));
//	memset(power10000, 0, nbrGroups * sizeof(limb));
//	power10000[0].x = MAX_LIMB_CONVERSION;
//	ptrDest = binary;
//	for (ptrSrc = decimal + digits - 1; ptrSrc >= decimal + DIGITS_PER_LIMB - 1; ptrSrc -= DIGITS_PER_LIMB)
//	{
//		int limbContents = 0;
//		for (digit = DIGITS_PER_LIMB - 1; digit >= 0; digit--)
//		{
//			limbContents = (limbContents * 10) + *(ptrSrc - digit) - '0';
//		}
//		(ptrDest++)->x = limbContents;
//	}
//	digit = 0;
//	multiplier = 1;
//	for (; ptrSrc >= decimal; ptrSrc--)
//	{
//		digit += multiplier * (*ptrSrc - '0');
//		multiplier *= 10;
//	}
//	(ptrDest++)->x = digit;
//	for (outerGroup = 1; outerGroup < nbrGroups; outerGroup += outerGroup)
//	{
//		for (innerGroup = 0; innerGroup < nbrGroups; innerGroup += 2 * outerGroup)
//		{
//			ptrBinary = binary + innerGroup;
//			multiply(power10000, ptrBinary + outerGroup, temp, outerGroup, NULL);
//			memset(ptrBinary + outerGroup, 0, outerGroup * sizeof(limb));
//			add(temp, ptrBinary, ptrBinary, 2 * outerGroup);
//		}
//		if (outerGroup * 2 < nbrGroups)
//		{    // Square power10000.
//			multiply(power10000, power10000, temp, outerGroup, NULL);
//			memcpy(power10000, temp, (outerGroup * 2) * sizeof(limb));
//		}
//	}
//	// Determine first non-significant group.
//	*bitGroups = 1;    // Initialize number of groups in advance.
//	for (groupNbr = nbrGroups - 1; groupNbr > 0; groupNbr--)
//	{
//		if ((binary + groupNbr)->x != 0)
//		{            // First non-significant group found.
//			*bitGroups = groupNbr + 1;
//			break;
//		}
//	}
//}

/* equivalent to:
_itoa((unsigned int)nbr, *ptrOutput, 10));
*ptroutput += (strlen(*ptrOutput)); */
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
