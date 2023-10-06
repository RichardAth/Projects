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


static void Bin2DecLoop(char** ppDest, int digit[], bool* pSignificantZero,
	int valueGrp, int grpLen, int* pGroupCtr, int* pDigits)
{
	int digits = *pDigits;
	int count;
	int value = valueGrp;
	char* ptrDest = *ppDest;
	int groupCtr = *pGroupCtr;
	bool significantZero = *pSignificantZero;
	for (count = 0; count < DIGITS_PER_LIMB; count++)
	{
		digit[count] = value % 10;
		value /= 10;
	}
	for (count = DIGITS_PER_LIMB - 1; count >= 0; count--)
	{
		if ((digit[count] != 0) || significantZero)
		{
			digits++;
			*ptrDest = (char)(digit[count] + '0');
			ptrDest++;
			if (groupCtr == 1)
			{
				*ptrDest = ' ';
				ptrDest++;
			}
			significantZero = true;
		}
		groupCtr--;
		if (groupCtr == 0)
		{
			groupCtr = grpLen;
		}
	}
	*ppDest = ptrDest;
	*pGroupCtr = groupCtr;
	*pSignificantZero = significantZero;
	*pDigits = digits;
}


// Convert little-endian number to a string with space every groupLen digits.
  // In order to perform a faster conversion, use groups of DIGITS_PER_LIMB digits.
void Bin2Dec(char** ppDecimal, const limb* binary, int nbrLimbs, int groupLength)
{
	int grpLen = groupLength;
	int len;
	int index;
	int index2;
	const limb* ptrSrc = binary + nbrLimbs - 1;
	char* ptrDest;
	bool significantZero = false;
	int groupCtr;
	int digit[DIGITS_PER_LIMB];
	int digits = 0;
	bool showDigitsText = true;
	unsigned int firstMult;
	unsigned int secondMult;
	firstMult = (unsigned int)BITS_PER_GROUP / 2U;
	firstMult = 1U << firstMult;
	secondMult = LIMB_RANGE / firstMult;

	if (grpLen <= 0)
	{
		grpLen = -grpLen;
		showDigitsText = false;
	}
	power10000[0] = *ptrSrc % MAX_LIMB_CONVERSION;
	power10000[1] = *ptrSrc / MAX_LIMB_CONVERSION;
	len = ((power10000[1] == 0) ? 1 : 2); // Initialize array length.
	for (index = nbrLimbs - 2; index >= 0; index--)
	{
		double dCarry;
		double dQuotient;
		limb* ptrPower;

		// Multiply by firstMult and then by secondMult, so there is never
		// more than 53 bits in the product.

		ptrPower = power10000;
		dQuotient = 0;
		for (index2 = 0; index2 < len; index2++)
		{
			dCarry = dQuotient + ((double)*ptrPower * (double)firstMult);
			dQuotient = floor(dCarry / (double)MAX_LIMB_CONVERSION);
			*ptrPower = (int)(dCarry - (dQuotient * (double)MAX_LIMB_CONVERSION));
			ptrPower++;
		}
		if (dQuotient != 0.0)
		{
			*ptrPower = (int)dQuotient;
			len++;
		}
		ptrPower = power10000;
		ptrSrc--;
		dQuotient = *ptrSrc;
		for (index2 = 0; index2 < len; index2++)
		{
			dCarry = dQuotient + ((double)*ptrPower * (double)secondMult);
			dQuotient = floor(dCarry / (double)MAX_LIMB_CONVERSION);
			*ptrPower = (int)(dCarry - (dQuotient * (double)MAX_LIMB_CONVERSION));
			ptrPower++;
		}
		if (dQuotient != 0.0)
		{
			*ptrPower = (int)dQuotient;
			ptrPower++;
			len++;
		}
	}
	// At this moment the array power10000 has the representation
	// of the number in base 10000 in little-endian. Convert to
	// ASCII separating every groupLength characters.
	ptrDest = *ppDecimal;
	ptrSrc = &power10000[len - 1];
	groupCtr = len * DIGITS_PER_LIMB;
	if (grpLen != 0)
	{
		groupCtr %= grpLen;
		if (groupCtr == 0)
		{
			groupCtr = grpLen;
		}
	}
	for (index = len; index > 0; index--)
	{
		int value = *ptrSrc;
		ptrSrc--;
		Bin2DecLoop(&ptrDest, digit, &significantZero, value, grpLen, &groupCtr, &digits);
	}
	if (!significantZero)
	{     // Number is zero.
		*ptrDest = '0';
		ptrDest++;
	}
	else if ((digits > 30) && showDigitsText)
	{
		*ptrDest = '(';
		ptrDest++;
		int2dec(&ptrDest, digits);
		copyStr(&ptrDest, (lang ? " dígitos)" : " digits)"));
	}
	else
	{
		if (ptrDest > *ppDecimal)
		{
			ptrDest--;               // Delete trailing space.
		}
	}
	*ptrDest = '\0';             // Add terminator.
	*ppDecimal = ptrDest;
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
	power10000[0] = *ptrSrc % MAX_LIMB_CONVERSION;
	power10000[1] = *ptrSrc / MAX_LIMB_CONVERSION;
	len = (power10000[1] == 0 ? 1 : 2); // Initialize array length.
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
			dCarry = dQuotient + (double)*ptrPower * FIRST_MULT;
			dQuotient = floor(dCarry / MAX_LIMB_CONVERSION);
			*(ptrPower++) = (int)(dCarry - dQuotient * MAX_LIMB_CONVERSION);
		}
		if (dQuotient != 0)
		{
			*(ptrPower++) = (int)dQuotient;
			len++;
		}
		ptrPower = power10000;
		dQuotient = (double)*(--ptrSrc);
		for (index2 = 0; index2 < len; index2++)
		{
			dCarry = dQuotient + (double)*ptrPower * SECOND_MULT;
			dQuotient = floor(dCarry / MAX_LIMB_CONVERSION);
			*(ptrPower++) = (int)(dCarry - dQuotient * MAX_LIMB_CONVERSION);
		}
		if (dQuotient != 0)
		{
			*(ptrPower++) = (int)dQuotient;
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
		int value = (int)*(ptrSrc--);
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
		std::strcpy(ptrDest, " digits)");
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


void BigInteger2Dec(char** ppDecimal, const BigInteger* pBigInt, int groupLength)
{
	char* ptrDecimal = *ppDecimal;
	if (pBigInt->sign == SIGN_NEGATIVE)
	{
		*ptrDecimal = '-';
		ptrDecimal++;
	}
	Bin2Dec(&ptrDecimal, pBigInt->limbs, pBigInt->nbrLimbs, groupLength);
	*ppDecimal = ptrDecimal;
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

/* instant print */
void PrintBigInteger(const BigInteger* pBigInt, int groupLength) {
	char buffer[5000] = { 0 };
	char* pbuffer = buffer;
	BigInteger2Dec(&pbuffer, pBigInt, groupLength);
	printf("%s", buffer);
}

void copyStr(char** pptrString, const char* stringToCopy)
{
	char* ptrString = *pptrString;
	const char* ptrStringToCopy = stringToCopy;
	while (*ptrStringToCopy != '\0')
	{
		*ptrString = *ptrStringToCopy;
		ptrString++;
		ptrStringToCopy++;
	}
	*ptrString = '\0';
	*pptrString = ptrString;
}
