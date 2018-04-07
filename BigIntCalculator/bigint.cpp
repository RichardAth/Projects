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
#include <stdexcept>
#include <limits.h>
#include <string>
#include <assert.h>
#include <math.h>
#include "bignbr.h"
#include "expression.h"

static BigInteger Temp, Temp2, Temp3, Base, Power, expon;
static char ProcessExpon[5000];
static char primes[10003];
limb Mult1[MAX_LEN];
extern limb Mult2[MAX_LEN];
static limb Mult3[MAX_LEN];
static limb Mult4[MAX_LEN];
int q[MAX_LEN];
extern limb TestNbr[MAX_LEN];
extern limb MontgomeryMultR1[MAX_LEN];

/* copy source to dest */
BigInteger &CopyBigInt(BigInteger &pDest, const BigInteger &pSrc) {
	pDest.sign = pSrc.sign;
	pDest.nbrLimbs = pSrc.nbrLimbs;
	memcpy(pDest.limbs, pSrc.limbs, (pSrc.nbrLimbs) * sizeof(limb));
	while (pDest.nbrLimbs > 1) {
		if (pDest.limbs[pDest.nbrLimbs - 1].x == 0)
			pDest.nbrLimbs--;
		else
			break;
	}
	return pDest;
}

 /* sum = addend1 + addend2 */
void BigIntAdd(const BigInteger &pAddend1, const BigInteger &pAddend2, BigInteger &pSum) {
	int ctr, nbrLimbs;
	const limb *ptrBiggerAdd, *ptrSmallerAdd;
	limb *ptrSum;
	bool A1Smaller = false;

	if (pAddend1.nbrLimbs < pAddend2.nbrLimbs) {
		A1Smaller = true;
		/* the absolute value of addend1 is less than the absolute value of addend2.*/
	}           
	
	else if (pAddend1.nbrLimbs == pAddend2.nbrLimbs) 	{
		for (ctr = pAddend1.nbrLimbs - 1; ctr >= 0; ctr--) {
			if (pAddend1.limbs[ctr].x != pAddend2.limbs[ctr].x) {
				break;
			}
		}
		if (ctr >= 0 && pAddend1.limbs[ctr].x < pAddend2.limbs[ctr].x) {
			/* the absolute value of addend1 is less than the absolute value of addend2.*/
			A1Smaller = true;
		} 
	}
	if (A1Smaller) {
		nbrLimbs = pAddend1.nbrLimbs;
		ptrBiggerAdd = pAddend2.limbs;
		ptrSmallerAdd = pAddend1.limbs;
	}
	else {
		// the absolute value of addend1 is >= the absolute value of addend2.
		nbrLimbs = pAddend2.nbrLimbs;
		ptrBiggerAdd = pAddend1.limbs;
		ptrSmallerAdd = pAddend2.limbs;
	}
	ptrSum = pSum.limbs;

	if (pAddend1.sign == pAddend2.sign)
	{             // Both addends have the same sign. Sum their absolute values.
		unsigned int carry = 0;
		for (ctr = 0; ctr < nbrLimbs; ctr++) {
			carry = (carry >> BITS_PER_GROUP) + (unsigned int)ptrBiggerAdd[ctr].x +
				(unsigned int)ptrSmallerAdd[ctr].x;
			ptrSum[ctr].x = (int)(carry & MAX_INT_NBR);
		}
		if (A1Smaller)
			nbrLimbs = pAddend2.nbrLimbs;
		else 
			nbrLimbs = pAddend1.nbrLimbs;
		for (; ctr < nbrLimbs; ctr++) {
			carry = (carry >> BITS_PER_GROUP) + (unsigned int)ptrBiggerAdd[ctr].x;
			ptrSum[ctr].x = (int)(carry & MAX_INT_NBR);
		}
		if (carry >= LIMB_RANGE)
		{
			ptrSum[ctr].x = 1;
			nbrLimbs++;
		}
	}
	else {    // addends have different signs. Subtract their absolute values.
		int borrow = 0;
		for (ctr = 0; ctr < nbrLimbs; ctr++) {
			borrow = (borrow >> BITS_PER_INT_GROUP) + ptrBiggerAdd[ctr].x - ptrSmallerAdd[ctr].x;
			ptrSum[ctr].x = borrow & MAX_INT_NBR;
		}
		if (A1Smaller)
			nbrLimbs = pAddend2.nbrLimbs;
		else
			nbrLimbs = pAddend1.nbrLimbs;
		for (; ctr < nbrLimbs; ctr++) {
			borrow = (borrow >> BITS_PER_INT_GROUP) + ptrBiggerAdd[ctr].x;
			ptrSum[ctr].x = borrow & MAX_INT_NBR;
		}
	}

	while (nbrLimbs > 1 && pSum.limbs[nbrLimbs - 1].x == 0) {    
		nbrLimbs--;  // delete leading zeros.
	}

	pSum.nbrLimbs = nbrLimbs;
	if (A1Smaller)
		pSum.sign = pAddend2.sign;  // use sign of addend with larger absolute value
	else
		pSum.sign = pAddend1.sign;

	if (pSum.nbrLimbs == 1 && pSum.limbs[0].x == 0) {          
		pSum.sign = SIGN_POSITIVE; // Result is zero.
	}
}

/* Dest = -Dest */
void BigIntNegate (BigInteger &pDest) {
	
	if (pDest.sign == SIGN_POSITIVE && (pDest.nbrLimbs != 1 || pDest.limbs[0].x != 0))
	{
		pDest.sign = SIGN_NEGATIVE;  // pDest > 0, now < 0
	}
	else
	{
		pDest.sign = SIGN_POSITIVE; // pDest <=0, now > 0
	}
}

/* Difference = Minuend - Subtrahend */
void BigIntSubt(const BigInteger &pMinuend, const BigInteger &pSubtrahend, 
	BigInteger &pDifference) {
	BigInteger temp;
	CopyBigInt(temp, pSubtrahend);
	BigIntNegate(temp);
	BigIntAdd(pMinuend, temp, pDifference);
}

/* Factor1 will be expanded to the length of Factor2 or vice versa 
Changed design from returning error code to throw an exception, because error 
codes need to be tested for after every call, passed up to the calling routine,
tested for again and so on up to the top level. In practise this wasn't done consistently. */
void BigIntMultiply(const BigInteger &pFactor1, const BigInteger &pFactor2, BigInteger &pProduct)
{
	int nbrLimbsFactor1 = pFactor1.nbrLimbs;
	int nbrLimbsFactor2 = pFactor2.nbrLimbs;
	int nbrLimbs;
	if ((pFactor1.nbrLimbs == 1 && pFactor1.limbs[0].x == 0) ||
		(pFactor2.nbrLimbs == 1 && pFactor2.limbs[0].x == 0))
	{    // one or both factors are zero.
		pProduct.nbrLimbs = 1;      // Product is zero.
		pProduct.limbs[0].x = 0;
		pProduct.sign = SIGN_POSITIVE;
		return;
	}
	if (pFactor1.nbrLimbs + pFactor2.nbrLimbs > 66438 / BITS_PER_GROUP + 1)  // 2^66438 ~ 10^20000
	{
		//return EXPR_INTERM_TOO_HIGH;
		std::string line = std::to_string(__LINE__);
		std::string mesg = "cannot multiply: product out of range ";
		mesg += __func__;
		mesg += " line ";  mesg += line;
		mesg += " in file "; mesg += __FILE__;
		throw std::range_error(mesg);
	}

	/* Factor1 will be expanded to the length of Factor2 or vice versa. It is necessary to
	override the const-ness, but the value represented does not change */
	if (nbrLimbsFactor1<nbrLimbsFactor2)
	{
		memset(&((BigInteger&)pFactor1).limbs[nbrLimbsFactor1], 0, (nbrLimbsFactor2 - nbrLimbsFactor1) * sizeof(limb));
		nbrLimbs = nbrLimbsFactor2;
	}
	else
	{
		memset(&((BigInteger&)pFactor2).limbs[nbrLimbsFactor2], 0, (nbrLimbsFactor1 - nbrLimbsFactor2) * sizeof(limb));
		nbrLimbs = nbrLimbsFactor1;
	}
	multiply(&pFactor1.limbs[0], &pFactor2.limbs[0], 
		&pProduct.limbs[0], nbrLimbs, &nbrLimbs);
	nbrLimbs = nbrLimbsFactor1 + nbrLimbsFactor2;
	if (nbrLimbs > MAX_LEN)  // limit applied earlier is probably lower, this is just insurance
	{
		std::string line = std::to_string(__LINE__);
		std::string mesg = "cannot multiply: product out of range ";
		mesg += __func__;
		mesg += " line ";  mesg += line;
		mesg += " in file "; mesg += __FILE__;
		throw std::range_error(mesg);
	}
	if (pProduct.limbs[nbrLimbs - 1].x == 0)
	{
		nbrLimbs--;
	}
	pProduct.nbrLimbs = nbrLimbs;
	if (nbrLimbs == 1 && pProduct.limbs[0].x == 0)
	{
		pProduct.sign = SIGN_POSITIVE;  // product is zero
	}
	else
	{
		if (pFactor1.sign == pFactor2.sign)
		{
			pProduct.sign = SIGN_POSITIVE;
		}
		else
		{
			pProduct.sign = SIGN_NEGATIVE;
		}
	}
	return;
}

/* calculate Dividend mod Divisor 
uses global variable base */
void BigIntRemainder(const BigInteger &pDividend, const BigInteger &pDivisor,
	BigInteger &pRemainder) {
	if (pDivisor.limbs[0].x == 0 && pDivisor.nbrLimbs == 1)
	{   // If divisor = 0, then remainder is the dividend.
		return;
	}
	CopyBigInt(Temp2, pDividend);
	BigIntDivide(pDividend, pDivisor, Base);   // Get quotient of division.
	BigIntMultiply(Base, pDivisor, Base);
	BigIntSubt(Temp2, Base, pRemainder);
	return;
}
 
/* bigint = value */
void intToBigInteger(BigInteger &bigint, int value) {
	if (value >= 0)
	{
		bigint.limbs[0].x = value;
		bigint.sign = SIGN_POSITIVE;
	}
	else
	{
		bigint.limbs[0].x = -value;
		bigint.sign = SIGN_NEGATIVE;
	}
	bigint.nbrLimbs = 1;
}

/* bigint = value */
void longToBigInteger(BigInteger &bigint, long long value) {
	int nbrLimbs = 0;
	bigint.sign = SIGN_POSITIVE;
	if (value < 0)
	{
		bigint.sign = SIGN_NEGATIVE;
		value = -value;
	}
	do
	{
		bigint.limbs[nbrLimbs++].x = (int)value & MAX_VALUE_LIMB;
		value >>= BITS_PER_GROUP;
	} while (value != 0);
	bigint.nbrLimbs = nbrLimbs;
}

/* BigInt = e^logar*/
void expBigNbr(BigInteger &bigInt, double logar)
{
	int mostSignificantLimb;
	logar /= log(2);  // convert log to base 2
	bigInt.sign = SIGN_POSITIVE;
	bigInt.nbrLimbs = (int)floor(logar / BITS_PER_GROUP);
	mostSignificantLimb = (int)floor(exp((logar - BITS_PER_GROUP*bigInt.nbrLimbs) * log(2)) + 0.5);
	if (mostSignificantLimb == LIMB_RANGE)
	{
		mostSignificantLimb = 1;
		bigInt.nbrLimbs++;
	}
	bigInt.nbrLimbs++;
	if (bigInt.nbrLimbs > MAX_LEN) {
		std::string line = std::to_string(__LINE__);
		std::string mesg = "number too big : cannot expand BigInteger: ";
		mesg += __func__;
		mesg += " line ";  mesg += line;
		mesg += " in file "; mesg += __FILE__;
		throw std::range_error(mesg);
	}
	memset(bigInt.limbs, 0, bigInt.nbrLimbs * sizeof(limb));
	bigInt.limbs[bigInt.nbrLimbs - 1].x = mostSignificantLimb;
}

/* estimate natural log of BigInt. Only the most significant 62 bits are 
taken into account because floating point numbers have limited accuracy anyway. */
double logBigNbr(const BigInteger &pBigInt)
{
	int nbrLimbs;
	double logar;
	nbrLimbs = pBigInt.nbrLimbs;
	if (nbrLimbs == 1)
	{
		logar = log((double)(pBigInt.limbs[0].x));
	}
	else
	{
		double value = pBigInt.limbs[nbrLimbs - 2].x +
			(double)pBigInt.limbs[nbrLimbs - 1].x * LIMB_RANGE;
		if (nbrLimbs == 2)
		{
			logar = log(value);
		}
		else
		{
			logar = log(value + (double)pBigInt.limbs[nbrLimbs - 3].x / LIMB_RANGE);
		}
		logar += (double)((nbrLimbs - 2)*BITS_PER_GROUP)*log(2);
	}
	return logar;
}

/* estimate natural log of BigNbr. Only the most significant 62 bits are
taken into account because floating point numbers have limited accuracy anyway. */
double logLimbs(const limb *pBigNbr, int nbrLimbs)
{
	double logar;
	if (nbrLimbs > 1)
	{
		logar = log( (double)pBigNbr[nbrLimbs - 2].x +
				(double)pBigNbr[nbrLimbs - 1].x * LIMB_RANGE) +
			(double)(nbrLimbs - 2)*log((double)LIMB_RANGE);
	}
	else
	{
		logar = log((double)pBigNbr[nbrLimbs - 1].x) +
			(double)(nbrLimbs - 1)*log((double)LIMB_RANGE);
	}
	return logar;
}

/* calculate base^exponent. Throw exception if result is out of range */
void BigIntPowerIntExp(const BigInteger &pBase, int exponent, BigInteger &pPower)
{
	int mask;
	double base;
	if (pBase.nbrLimbs == 1 && pBase.limbs[0].x == 0)
	{     // Base = 0 -> power = 0
		pPower.limbs[0].x = 0;
		pPower.nbrLimbs = 1;
		pPower.sign = SIGN_POSITIVE;
		return; // base = 0, so result is zero
	}
	base = logBigNbr(pBase);
	if (base*(double)exponent > 46051)
	{   // More than 20000 digits. 46051 = log(10^20000)
		std::string line = std::to_string(__LINE__);
		std::string mesg = "number too big : cannot create exponent: ";
		mesg += __func__;
		mesg += " line ";  mesg += line;
		mesg += " in file "; mesg += __FILE__;
		throw std::range_error(mesg);
	}
	CopyBigInt(Base, pBase);
	pPower.sign = SIGN_POSITIVE;
	pPower.nbrLimbs = 1;
	pPower.limbs[0].x = 1;
	for (mask = 1 << 30; mask != 0; mask >>= 1)
	{
		if ((exponent & mask) != 0)
		{
			for (; mask != 0; mask >>= 1)
			{
				BigIntMultiply(pPower, pPower, pPower);
				if ((exponent & mask) != 0)
				{
					BigIntMultiply(pPower, Base, pPower);
				}
			}
			break;
		}
	}
	return;
}

/* divide by 2, use right shift for speed */
void BigIntDivide2(BigInteger &pArg)
{
	int nbrLimbs = pArg.nbrLimbs;
	int ctr = nbrLimbs - 1;
	unsigned int carry;
	//limb *ptrLimb = &pArg->limbs[ctr];
	limb *ptrLimb = pArg.limbs;
	carry = 0;
	for (; ctr >= 0; ctr--)
	{
		carry = (carry << BITS_PER_GROUP) + (unsigned int)ptrLimb[ctr].x;
		ptrLimb[ctr].x = (int)(carry >> 1);
		carry &= 1;
	}
	if (nbrLimbs > 1 && pArg.limbs[nbrLimbs - 1].x == 0)
	{     // Most significant limb is zero, so reduce size by one limb.
		pArg.nbrLimbs--;
	}
}

/* arg = arg*2^power. Throw exception if product is too large */
static void BigIntMutiplyPower2(BigInteger &pArg, int power2)
{
	int ctr;
	int nbrLimbs = pArg.nbrLimbs;
	limb *ptrLimbs = pArg.limbs;

	for (; power2 > 0; power2--) {
		/*each time round the loop multiplies arg by 2 */
		unsigned int carry = 0;
		for (ctr = 0; ctr < nbrLimbs; ctr++)
		{
			carry += (unsigned int)ptrLimbs[ctr].x << 1;
			ptrLimbs[ctr].x = (int)(carry & MAX_VALUE_LIMB);
			carry >>= BITS_PER_GROUP;
		}
		if (carry != 0)
		{
			ptrLimbs[ctr].x = (int)carry;
			nbrLimbs++;
			if (nbrLimbs > MAX_LEN) {
				std::string line = std::to_string(__LINE__);
				std::string mesg = "number too big : cannot do multiplication : ";
				mesg += __func__;
				mesg += " line ";  mesg += line;
				mesg += " in file "; mesg += __FILE__;
				throw std::range_error(mesg);
			}
		}
	}
	pArg.nbrLimbs = nbrLimbs;
}

/* return true if Nbr1 == Nbr2 */
bool TestBigNbrEqual(const BigInteger &Nbr1, const BigInteger &Nbr2) {
	int ctr;
	/*const limb *ptrLimbs1 = Nbr1.limbs;
	const limb *ptrLimbs2 = Nbr2.limbs;*/
	auto N1Limbs = Nbr1.nbrLimbs;
	auto N2Limbs = Nbr2.nbrLimbs;
	while (N1Limbs > 1)
		if (Nbr1.limbs[N1Limbs - 1].x == 0)
			N1Limbs--;
		else
			break;
	while (N2Limbs > 1)
		if (Nbr2.limbs[N2Limbs - 1].x == 0)
			N2Limbs--;
		else
			break;

	if (N1Limbs != N2Limbs) {        
		return false;  // Sizes of numbers are different.
	}
	if (Nbr1.sign != Nbr2.sign) { 
	       // Sign of numbers are different.
		if (N1Limbs == 1 && Nbr1.limbs[0].x == 0 && Nbr2.limbs[0].x == 0) {              
			return true; // Both numbers are zero.
		}
		return false; // differents signs, therefore cannot be equal
	}

	// Check whether both numbers are equal.
	for (ctr = N1Limbs - 1; ctr >= 0; ctr--) {
		if (Nbr1.limbs[ctr].x != Nbr2.limbs[ctr].x) {
			return false;  // Numbers are different.
		}
	}        
	return true;  // Numbers are equal.
}

/* return true if Nbr1 < Nbr2 */
bool TestBigNbrLess(const BigInteger &Nbr1, const BigInteger &Nbr2) {
	int ctr;
	/*const limb *ptrLimbs1 = Nbr1.limbs;
	const limb *ptrLimbs2 = Nbr2.limbs;*/
	auto N1Limbs = Nbr1.nbrLimbs;
	auto N2Limbs = Nbr2.nbrLimbs;
	while (N1Limbs > 1)
		if (Nbr1.limbs[N1Limbs - 1].x == 0)
			N1Limbs--;
		else
			break;
	while (N2Limbs > 1)
		if (Nbr2.limbs[N2Limbs - 1].x == 0)
			N2Limbs--;
		else
			break;


	if (Nbr1.sign != Nbr2.sign) {
		// Sign of numbers are different.
		if (N1Limbs == 1 && Nbr1.limbs[0].x == 0 && Nbr2.limbs[0].x == 0) {
			return false; // Both numbers are zero i.e Nbr1 not less than Nbr2
		}
		else return (Nbr1.sign == SIGN_NEGATIVE);
	}

	/* numbers have same sign */
	if (N1Limbs != N2Limbs) {
		/* length of numbers is different*/
		if (Nbr1.sign == SIGN_POSITIVE)
			return N1Limbs < N2Limbs;
		else
			return N1Limbs > N2Limbs;
	}

	// numbers have same sign and length. Check whether both numbers are equal.
	for (ctr = N1Limbs - 1; ctr >= 0; ctr--) {
		if (Nbr1.limbs[ctr].x < Nbr2.limbs[ctr].x) {
			return true;  // Nbr1 < Nbr2.
		}
		else if (Nbr1.limbs[ctr].x > Nbr2.limbs[ctr].x) {
			return false;  // Nbr1 > Nbr2.
		}
	}
	return false;  // Numbers are equal.
}

/* return true if Nbr1 > Nbr2 */
bool TestBigNbrGtr(const BigInteger &Nbr1, const BigInteger &Nbr2) {
	return TestBigNbrLess(Nbr2, Nbr1);
}

/* return true if Nbr1 >= Nbr2 */
bool TestBigNbrGe(const BigInteger &Nbr1, const BigInteger &Nbr2) {
	return !TestBigNbrLess(Nbr1, Nbr2);
}

/* return true if Nbr1 <= Nbr2 */
bool TestBigNbrLe(const BigInteger &Nbr1, const BigInteger &Nbr2) {
	return !TestBigNbrLess(Nbr2, Nbr1);
}

/* calculate GCD of arg1 & arg2*/
void BigIntGcd(const BigInteger &pArg1, const BigInteger &pArg2, BigInteger &pResult)
{
	int nbrLimbs1 = pArg1.nbrLimbs;
	int nbrLimbs2 = pArg2.nbrLimbs;
	int power2;
	if (nbrLimbs1 == 1 && pArg1.limbs[0].x == 0)
	{               // First argument is zero, so the GCD is second argument.
		CopyBigInt(pResult, pArg2);
		return;
	}
	if (nbrLimbs2 == 1 && pArg2.limbs[0].x == 0)
	{               // Second argument is zero, so the GCD is first argument.
		CopyBigInt(pResult, pArg1);
		return;
	}
	// Reuse Base and Power temporary variables.
	CopyBigInt(Base, pArg1);   // Base = Arg1
	CopyBigInt(Power, pArg2);   // Power = Arg2
	Base.sign = SIGN_POSITIVE;
	Power.sign = SIGN_POSITIVE;
	power2 = 0;
	while (((Base.limbs[0].x | Power.limbs[0].x) & 1) == 0)
	{  // Both values are even
		BigIntDivide2(Base);
		BigIntDivide2(Power);
		power2++;
	}
	//while (TestBigNbrEqual(Base, Power) == 0)
	while (Base != Power)
	{    // Main GCD loop.
		if ((Base.limbs[0].x & 1) == 0)
		{          // Number is even. Divide it by 2.
			BigIntDivide2(Base);
			continue;
		}
		if ((Power.limbs[0].x & 1) == 0)
		{          // Number is even. Divide it by 2.
			BigIntDivide2(Power);
			continue;
		}
		BigIntSubt(Base, Power, pResult);
		if (pResult.sign == SIGN_POSITIVE)
		{
			CopyBigInt(Base, pResult);
			BigIntDivide2(Base);
		}
		else
		{
			CopyBigInt(Power, pResult);
			Power.sign = SIGN_POSITIVE;
			BigIntDivide2(Power);
		}
	}
	CopyBigInt(pResult, Base);
	BigIntMutiplyPower2(pResult, power2); /* pResult *= 2^power     */
}

/* add addend to big number */
static void addToAbsValue(limb *pLimbs, int *pNbrLimbs, int addend) {
	int ctr;
	int nbrLimbs = *pNbrLimbs;
	pLimbs[0].x += addend;
	if ((unsigned int)pLimbs->x < LIMB_RANGE)
	{     // No overflow. Go out of routine.
		return;
	}
	pLimbs[0].x -= LIMB_RANGE;
	for (ctr = 1; ctr < nbrLimbs; ctr++)
	{
		//pLimbs++;        // Point to next most significant limb.
		if (pLimbs[ctr].x != MAX_INT_NBR)
		{   // No overflow. Go out of routine.
			pLimbs[ctr].x++;   // Add carry.
			return;
		}
		pLimbs[ctr].x = 0;
	}
	(*pNbrLimbs)++;        // Result has an extra limb.
	pLimbs[ctr].x = 1;   // Most significant limb must be 1.
}

/* subtract from big number */
static void subtFromAbsValue(limb *pLimbs, int *pNbrLimbs, int subt) {
	int nbrLimbs = *pNbrLimbs;
	pLimbs[0].x -= subt;
	if (pLimbs[0].x < 0)
	{
		int ctr;
		for (ctr = 1; ctr < nbrLimbs; ctr++)
		{
			pLimbs[ctr - 1].x += LIMB_RANGE;
			if (--(pLimbs[ctr].x) >= 0)
			{
				break;
			}
		}
		if (nbrLimbs > 1 && pLimbs[nbrLimbs - 1].x == 0)
		{
			nbrLimbs--;
		}
	}
	*pNbrLimbs = nbrLimbs;
}

/* i = (i-subt)/divisor   */
void subtractdivide(BigInteger &pBigInt, int subt, int divisor)
{
	int nbrLimbs = pBigInt.nbrLimbs;
	// Point to most significant limb.
	limb *pLimbs;
	int ctr;
	int remainder = 0;
	double dDivisor = (double)divisor;
	double dInvDivisor = 1 / dDivisor;
	double dLimb = (double)LIMB_RANGE;

#if 0
	char *ptrOutput = output;
	*ptrOutput++ = '2';
	*ptrOutput++ = '(';
	int2dec(&ptrOutput, pBigInt->sign);
	*ptrOutput++ = ',';
	*ptrOutput++ = ' ';
	int2dec(&ptrOutput, pBigInt->nbrLimbs);
	*ptrOutput++ = ';';
	*ptrOutput++ = ' ';
	int2dec(&ptrOutput, pBigInt->limbs[0].x);
	*ptrOutput++ = ',';
	*ptrOutput++ = ' ';
	int2dec(&ptrOutput, pBigInt->limbs[1].x);
	*ptrOutput++ = ',';
	*ptrOutput++ = ' ';
	*ptrOutput++ = ')';
	*ptrOutput++ = ',';
	*ptrOutput++ = ' ';
	int2dec(&ptrOutput, subt);
	*ptrOutput++ = ',';
	*ptrOutput++ = ' ';
	int2dec(&ptrOutput, divisor);
	//  databack(output);
	if ((unsigned int)pBigInt->limbs[0].x >= LIMB_RANGE)
	{
		remainder = 1;
	}
#endif
	if (subt >= 0)
	{
		if (pBigInt.sign == SIGN_POSITIVE)
		{               // Subtract subt to absolute value.
			subtFromAbsValue(pBigInt.limbs, &nbrLimbs, subt);
		}
		else
		{               // Add subt to absolute value.
			addToAbsValue(pBigInt.limbs, &nbrLimbs, subt);
		}
	}
	else
	{
		if (pBigInt.sign == SIGN_POSITIVE)
		{               // Subtract subt to absolute value.
			addToAbsValue(pBigInt.limbs, &nbrLimbs, -subt);
		}
		else
		{               // Add subt to absolute value.
			subtFromAbsValue(pBigInt.limbs, &nbrLimbs, -subt);
		}
	}
	//pLimbs = pBigInt->limbs + nbrLimbs - 1;
	pLimbs = pBigInt.limbs;
	// Divide number by divisor.
	for (ctr = nbrLimbs - 1; ctr >= 0; ctr--)
	{
		double dDividend, dQuotient;
		unsigned int quotient, dividend;
		dividend = (remainder << BITS_PER_INT_GROUP) + pLimbs[ctr].x;
		dDividend = (double)remainder * dLimb + pLimbs[ctr].x;
		dQuotient = dDividend * dInvDivisor + 0.5;
		quotient = (unsigned int)dQuotient;   // quotient has correct value or 1 more.
		remainder = dividend - quotient * divisor;
		if (remainder < 0)
		{     // remainder not in range 0 <= remainder < divisor. Adjust.
			quotient--;
			remainder += divisor;
		}
		//(pLimbs--)->x = (int)quotient;
		pLimbs[ctr].x = (int)quotient;
	}

	while (nbrLimbs > 1 && pBigInt.limbs[nbrLimbs - 1].x == 0)
	{   // Most significant limb is now zero, so discard it.
		nbrLimbs--;
	}
	pBigInt.nbrLimbs = nbrLimbs;
}

/* calculate BigInt modulo divisor */
int getRemainder(const BigInteger &pBigInt, int divisor) {
	int ctr;
	int remainder = 0;
	int nbrLimbs = pBigInt.nbrLimbs;
	double dDivisor = (double)divisor;
	double dLimb = 0x80000000;
	const limb *pLimb = pBigInt.limbs;    // point to first limb
	for (ctr = nbrLimbs - 1; ctr >= 0; ctr--)
	{
		int quotient, dividend;
		double dQuotient, dDividend;
		dividend = (remainder << BITS_PER_INT_GROUP) + pLimb[ctr].x;
		dDividend = (double)remainder * dLimb + pLimb[ctr].x;
		dQuotient = floor(dDividend / dDivisor + 0.5);
		quotient = (int)(unsigned int)dQuotient;   // quotient has correct value or 1 more.
		remainder = dividend - quotient * divisor;
		if ((unsigned int)remainder >= (unsigned int)divisor)
		{     // remainder not in range 0 <= remainder < divisor. Adjust.
			quotient--;
			remainder += divisor;
		}
	}
	if (pBigInt.sign == SIGN_NEGATIVE && remainder != 0) 	{
		remainder = divisor - remainder;
	}
	return remainder;
}

/* result += addend */
void addbigint(BigInteger &pResult, int addend) {
	int nbrLimbs = pResult.nbrLimbs;
	limb *pResultLimbs = pResult.limbs;
	auto sign = pResult.sign;

	if (addend < 0) {
		// reverse signs of addend and result
		addend = -addend;
		if (sign == SIGN_POSITIVE) 	{
			sign = SIGN_NEGATIVE;
		}
		else {
			sign = SIGN_POSITIVE;
		}
	}

	if (sign == SIGN_POSITIVE) {   // Add addend to absolute value of pResult.
		addToAbsValue(pResultLimbs, &nbrLimbs, addend);
	}
	else {  // Subtract addend from absolute value of pResult.
		if (nbrLimbs == 1) 	{
			pResultLimbs[0].x -= addend;
			if (pResultLimbs[0].x < 0) 	{
				pResultLimbs[0].x = -pResultLimbs[0].x;  // reverse sign of result
				BigIntNegate(pResult);
			}
		}
		else {     // More than one limb.
			subtFromAbsValue(pResultLimbs, &nbrLimbs, addend);
		}
	}
	pResult.nbrLimbs = nbrLimbs;
}


// Get number of bits of given big integer.
static int bitLength(const BigInteger &pBigNbr)
{
	unsigned int mask;
	int bitCount;
	const int lastLimb = pBigNbr.nbrLimbs - 1;
	int bitLen = lastLimb*BITS_PER_GROUP;
	unsigned int limb = (unsigned int)(pBigNbr.limbs[lastLimb].x);
	mask = 1;
	for (bitCount = 0; bitCount < BITS_PER_GROUP; bitCount++)
	{
		if (limb < mask)
		{
			break;
		}
		mask *= 2;
	}
	return bitLen + bitCount;
}

/* returns nbrMod^Expon%currentPrime*/
int intModPow(int NbrMod, int Expon, int currentPrime)
{
	unsigned int power = 1;
	unsigned int square = (unsigned int)NbrMod;
	while (Expon != 0)
	{
		if ((Expon & 1) == 1)
		{
			power = (power * square) % (unsigned int)currentPrime;
		}
		square = (square * square) % (unsigned int)currentPrime;
		Expon >>= 1;
	}
	return (int)power;
}

/* copy value to global var. temp */
static void InitTempFromInt(int value)
{
	Temp.nbrLimbs = 1;
	Temp.limbs[0].x = value;
	Temp.sign = SIGN_POSITIVE;
}

/* creates a BigInteger from a list of values.   
uses global value NumberLength for number of limbs. 1st entry in list is number of values
that follow */
void UncompressBigInteger(/*@in@*/const int *ptrValues, /*@out@*/BigInteger &bigint)
{
	if (NumberLength > MAX_LEN || NumberLength < 0 || ptrValues[0] > MAX_LEN) {
		std::string line = std::to_string(__LINE__);
		std::string mesg = "number too big : cannot convert to BigInteger: ";
		mesg += __func__;
		mesg += " line ";  mesg += line;
		mesg += " in file "; mesg += __FILE__;
		throw std::range_error(mesg);
	}
	limb *destLimb = bigint.limbs;
	if (NumberLength == 1)
	{
		destLimb[0].x = ptrValues[1];
		bigint.nbrLimbs = 1;
	}
	else
	{
		int ctr, nbrLimbs;
		nbrLimbs = ptrValues[0];
		bigint.nbrLimbs = nbrLimbs;
		for (ctr = 0; ctr < nbrLimbs; ctr++)
		{
			destLimb[ctr].x = ptrValues[ctr+1];
		}
		for (; ctr < NumberLength; ctr++)
		{
			destLimb[ctr].x = 0;  // zeroise any remaining limbs
		}
	}
}

/* uses global value NumberLength for number of limbs. */
static int getNbrLimbs(const limb *bigNbr)
{
	//const limb *ptrLimb = bigNbr + NumberLength;
	auto ix = NumberLength;
	while (ix > 0)
	{
		if (bigNbr[ix-1].x != 0)
		{
			return (int)(ix);
		}
		ix--;   // reduce length because most significant limb is zero
	}
	return 1;
}

/* creates a list of values from a BigInteger, 1st entry in list is number of 
values that follow. uses global value NumberLength for number of limbs. */
void CompressBigInteger(/*@out@*/int *ptrValues, /*@in@*/const BigInteger &bigint)
{
	const limb *destLimb = bigint.limbs;
	if (NumberLength == 1) {
		ptrValues[0] = 1;
		ptrValues[1] = destLimb[0].x;
	}
	else {
		int ctr, nbrLimbs;
		nbrLimbs = getNbrLimbs(bigint.limbs);  //nbrLimbs <= NumberLength
		//*ptrValues++ = nbrLimbs;
		ptrValues[0] = nbrLimbs;
		for (ctr = 0; ctr < nbrLimbs; ctr++) {
			ptrValues[ctr+1] = destLimb[ctr].x;
		}
	}
}

/* uses global value NumberLength for number of limbs. */
void UncompressLimbsBigInteger(/*@in@*/const limb *ptrValues, /*@out@*/BigInteger &bigint)
{
	if (NumberLength == 1)
	{
		bigint.limbs[0].x = ptrValues[0].x;
		bigint.nbrLimbs = 1;
	}
	else {
		memcpy(bigint.limbs, ptrValues, NumberLength * sizeof(limb));

		int nbrLimbs;
		for (nbrLimbs = NumberLength-1; nbrLimbs > 1; nbrLimbs--)
		{
			if (ptrValues[nbrLimbs].x != 0)
			{
				break;
			}
		}
		bigint.nbrLimbs = nbrLimbs+1;
	}
}

/* uses global value NumberLength for number of limbs. */
void CompressLimbsBigInteger(/*@out@*/limb *ptrValues, /*@in@*/const BigInteger &bigint)
{
	if (NumberLength == 1)
	{
		ptrValues[0].x = bigint.limbs[0].x;
	}
	else
	{
		int nbrLimbs = bigint.nbrLimbs;
		if (nbrLimbs >= NumberLength)
		{
			memcpy(ptrValues, bigint.limbs, NumberLength * sizeof(limb));
		}
		else
		{
			memcpy(ptrValues, bigint.limbs, nbrLimbs * sizeof(limb));
			memset(ptrValues + nbrLimbs, 0, (NumberLength - nbrLimbs) * sizeof(limb));
		}
	}
}

void UncompressIntLimbs(/*@in@*/const int *ptrValues, /*@out@*/limb *bigint, int nbrLen)
{
	int nbrLimbs = *ptrValues;
	memcpy(bigint, ptrValues + 1, nbrLimbs * sizeof(limb));
	memset(bigint + nbrLimbs, 0, (nbrLen - nbrLimbs) * sizeof(limb));
}

void CompressIntLimbs(/*@out@*/int *ptrValues, /*@in@*/const limb *bigint, int nbrLen)
{
	int nbrLimbs;
	memcpy(ptrValues + 1, bigint, (nbrLen - 1) * sizeof(limb));
	for (nbrLimbs = nbrLen - 1; nbrLimbs > 1; nbrLimbs--)
	{
		if (ptrValues[nbrLimbs] != 0)
		{
			break;
		}
	}
	*ptrValues = nbrLimbs;
}

// This routine checks whether the number pointed by pNbr is a perfect power. 
// If it is not, it returns one. If it is a perfect power, it returns the  
// exponent and  it fills the buffer pointed by pBase with the base.
int PowerCheck(const BigInteger &pBigNbr, BigInteger &pBase)
{
	limb *ptrLimb;
	double dN;
	int nbrLimbs = pBigNbr.nbrLimbs;
	const int maxExpon = bitLength(pBigNbr);  // if maxExpon > 5000 throw an exception
											// this corresponds to about 161 limbs 
	int h, j;
	int modulus;
	int intLog2root;
	int primesLength, Exponent;
	double log2N, log2root;
	int prime2310x1[] =
	{ 2311, 4621, 9241, 11551, 18481, 25411, 32341, 34651, 43891, 50821 };
	// Primes of the form 2310x+1.
	bool expon2 = true, expon3 = true, expon5 = true;
	bool expon7 = true, expon11 = true;
	for (h = 0; h < sizeof(prime2310x1) / sizeof(prime2310x1[0]); h++)
	{
		int testprime = prime2310x1[h];
		int mod = getRemainder(pBigNbr, testprime);
		if (expon2 && intModPow(mod, testprime / 2, testprime) > 1)
		{
			expon2 = false;
		}
		if (expon3 && intModPow(mod, testprime / 3, testprime) > 1)
		{
			expon3 = false;
		}
		if (expon5 && intModPow(mod, testprime / 5, testprime) > 1)
		{
			expon5 = false;
		}
		if (expon7 && intModPow(mod, testprime / 7, testprime) > 1)
		{
			expon7 = false;
		}
		if (expon11 && intModPow(mod, testprime / 11, testprime) > 1)
		{
			expon11 = false;
		}
	}

	primesLength = 2 * maxExpon + 3;
	if (primesLength >= sizeof(primes) / sizeof(primes[0]) || 
		maxExpon >= sizeof(ProcessExpon) / sizeof(ProcessExpon[0])) {	
		
		std::string line = std::to_string(__LINE__);
		std::string mesg = "number too big : cannot generate prime list. function : ";
		mesg += __func__;
		mesg += " line ";  mesg += line;
		mesg += " in file "; mesg += __FILE__;
		throw std::range_error(mesg);
	}
	
	for (h = 2; h <= maxExpon; h++)
	{
		ProcessExpon[h] = true;
	}
	for (h = 2; h < primesLength; h++)
	{
		primes[h] = true;
	}
	for (h = 2; h * h < primesLength; h++)
	{ // Generation of primes
		for (j = h * h; j < primesLength; j += h)
		{ // using Eratosthenes sieve
			primes[j] = false;
		}
	}
	for (h = 13; h < primesLength; h++)
	{
		if (primes[h])
		{
			int processed = 0;
			for (j = 2 * h + 1; j < primesLength; j += 2 * h)
			{
				if (primes[j])
				{
					modulus = getRemainder(pBigNbr, j);
					if (intModPow(modulus, j / h, j) > 1)
					{
						for (j = h; j <= maxExpon; j += h)
						{
							ProcessExpon[j] = false;
						}
						break;
					}
				}
				if (++processed > 10)
				{
					break;
				}
			}
		}
	}
	log2N = logBigNbr(pBigNbr) / log(2);
	for (Exponent = maxExpon; Exponent >= 2; Exponent--)
	{
		if (Exponent % 2 == 0 && !expon2)
		{
			continue; // Not a square
		}
		if (Exponent % 3 == 0 && !expon3)
		{
			continue; // Not a cube
		}
		if (Exponent % 5 == 0 && !expon5)
		{
			continue; // Not a fifth power
		}
		if (Exponent % 7 == 0 && !expon7)
		{
			continue; // Not a 7th power
		}
		if (Exponent % 11 == 0 && !expon11)
		{
			continue; // Not an 11th power
		}
		if (!ProcessExpon[Exponent])
		{
			continue;
		}
		// Initialize approximation to n-th root (n = Exponent).
		log2root = log2N / Exponent;
		intLog2root = (int)floor(log2root / BITS_PER_GROUP);
		nbrLimbs = intLog2root + 1;
		ptrLimb = &pBase.limbs[nbrLimbs - 1];
		dN = exp((log2root - intLog2root*BITS_PER_GROUP) * log(2));
		// All approximations must be >= than true answer.
		if (nbrLimbs == 1)
		{
			ptrLimb->x = (int)(unsigned int)ceil(dN);
			if ((unsigned int)ptrLimb->x == LIMB_RANGE)
			{
				nbrLimbs = 2;
				ptrLimb->x = 0;
				(ptrLimb + 1)->x = 1;
			}
		}
		else
		{
			dN += 1 / (double)LIMB_RANGE;
			ptrLimb->x = (int)trunc(dN);
			dN -= trunc(dN);
			(ptrLimb - 1)->x = (int)trunc(dN*LIMB_RANGE);
		}
		pBase.nbrLimbs = nbrLimbs;
		// Perform Newton iteration for n-th root.
		for (;;)
		{   // Check whether the approximate root is actually exact.
			BigIntPowerIntExp(pBase, Exponent - 1, Temp3); // Temp3 <- x^(e-1)
			BigIntMultiply(Temp3, pBase, Temp2);        // Temp2 <- x^e 
														  // necessary to override const-ness of BigNbr
			BigIntSubt(pBigNbr, Temp2, Temp2);            // Compare to radicand.
			if (Temp2.nbrLimbs == 1 && Temp2.limbs[0].x == 0)
			{                     // Perfect power, so go out.
				return Exponent;
			}
			if (Temp2.sign == SIGN_POSITIVE)
			{                     // x^e > radicand -> not perfect power, so go out.
				break;
			}
			BigIntDivide(pBigNbr, Temp3, Temp);         // Temp -> N/x^(e-1)
			BigIntSubt(Temp, pBase, Temp2);             // Temp2 -> N/x^(e-1) - x
			if (Temp2.nbrLimbs == 1 && Temp2.limbs[0].x == 0)
			{     // New approximation will be the same as previous. Go out.
				break;
			}
			InitTempFromInt(Exponent - 1);
			BigIntSubt(Temp2, Temp, Temp2);
			InitTempFromInt(Exponent);
			BigIntDivide(Temp2, Temp, Temp2);
			BigIntAdd(Temp2, pBase, pBase);
		}
	}
	pBase.nbrLimbs = pBigNbr.nbrLimbs;
	memcpy(pBase.limbs, pBigNbr.limbs, pBase.nbrLimbs * sizeof(limb));
	return 1;
}

/* return false if value mod p is not 1 */
bool checkOne(const limb *value, int nbrLimbs)
{
	int idx;
	for (idx = 0; idx < nbrLimbs; idx++)
	{
		if ((value++)->x != MontgomeryMultR1[idx].x)
		{
			return false;    // Go out if value is not 1 (mod p)
		}
	}
	return true;
}

/* return false if value mod p is not -1*/
bool checkMinusOne(const limb *value, int nbrLimbs)
{
	int idx;
	unsigned int carry;
	carry = 0;
	for (idx = 0; idx < nbrLimbs; idx++)
	{
		carry += (unsigned int)(value++)->x + (unsigned int)MontgomeryMultR1[idx].x;
		if ((carry & MAX_VALUE_LIMB) != (unsigned int)TestNbr[idx].x)
		{
			return false;    // Go out if value is not -1 (mod p)
		}
		carry >>= BITS_PER_GROUP;
	}
	return true;
}

// Find power of 2 that divides the number.
// output: pNbrLimbs = pointer to number of limbs
//         pShRight = pointer to power of 2.
void DivideBigNbrByMaxPowerOf2(int *pShRight, limb *number, int *pNbrLimbs)
{
	int power2 = 0;
	long long mask;
	int index, index2,  shRg;
	int nbrLimbs = *pNbrLimbs;
	// Start from least significant limb (number zero).
	for (index = 0; index < nbrLimbs; index++)
	{
		if (number[index].x != 0)
		{
			break;
		}
		power2 += BITS_PER_GROUP;
	}
	for (mask = 0x1; mask <= MAX_VALUE_LIMB; mask *= 2)
	{
		if ((number[index].x & mask) != 0)
		{
			break;
		}
		power2++;
	}
	// Divide number by this power.
	shRg = power2 % BITS_PER_GROUP; // Shift right bit counter
	if ((number[nbrLimbs - 1].x & (-(1 << shRg))) != 0)
	{   // Most significant bits set.
		*pNbrLimbs = nbrLimbs - index;
	}
	else
	{   // Most significant bits not set.
		*pNbrLimbs = nbrLimbs - index - 1;
	}
	// Move number shRg bits to the right.
	mask = (1 << shRg) - 1;
	for (index2 = index; index2 < nbrLimbs; index2++)
	{
		number[index2].x = ((number[index2].x >> shRg) |
			(number[index2 + 1].x << (BITS_PER_GROUP - shRg))) &
			MAX_VALUE_LIMB;
	}
	if (index > 0)
	{   // Move limbs to final position.
		memmove(number, &number[index], (nbrLimbs - index) * sizeof(limb));
	}
	*pShRight = power2;
}

// Calculate Jacobi symbol by following algorithm 2.3.5 of C&P book.
int JacobiSymbol(int upper, int lower)
{
	int tmp;
	int a = upper % lower;
	int m = lower;
	int t = 1;
	while (a != 0)
	{
		while ((a & 1) == 0)
		{     // a is even.
			a >>= 1;
			if ((m & 7) == 3 || (m & 7) == 5)
			{   // m = 3 or m = 5 (mod 8)
				t = -t;
			}
		}
		tmp = a; a = m; m = tmp;   // Exchange a and m.
		if ((a & m & 3) == 3)
		{   // a = 3 and m = 3 (mod 4)
			t = -t;
		}
		a = a % m;
	}
	if (m == 1 || m == -1)
	{
		return t;
	}
	return 0;
}

/* uses global value NumberLength for number of limbs. */
static void Halve(limb *pValue)
{
	if ((pValue[0].x & 1) == 0)
	{    // Number to halve is even. Divide by 2.
		DivBigNbrByInt((int *)pValue, 2, (int *)pValue, NumberLength);
	}
	else
	{    // Number to halve is odd. Add modulus and then divide by 2.
		AddBigNbr((int *)pValue, (int *)TestNbr, (int *)pValue, NumberLength + 1);
		DivBigNbrByInt((int *)pValue, 2, (int *)pValue, NumberLength + 1);
	}
}

// BPSW primality test:
// 1) If the input number is 2-SPRP composite, indicate composite and go out.
// 2) Find the first D in the sequence 5, -7, 9, -11, 13, -15, ...
//    for which the Jacobi symbol (D/n) is −1. Set P = 1 and Q = (1 - D) / 4.
// 3) Perform a strong Lucas probable prime test on n using parameters D, P,
//    and Q. If n is not a strong Lucas probable prime, then n is composite.
//    Otherwise, n is almost certainly prime.
// Output: 0 = probable prime.
//         1 = composite: not 2-Fermat pseudoprime.
//         2 = composite: does not pass 2-SPRP test.
//         3 = composite: does not pass strong Lucas test.
int BpswPrimalityTest(/*@in@*/const BigInteger *pValue)
{
	int i, Mult3Len, ctr, D, absQ, mult, mask, index, signPowQ;
	int insidePowering = false;
	int nbrLimbs = pValue->nbrLimbs;
	limb *limbs = (limb *)pValue->limbs;  // override const-ness of pValue
	if (nbrLimbs == 1 && limbs->x <= 2)
	{
		return 0;    // Indicate prime.
	}
	if ((limbs->x & 1) == 0)
	{
		return 1;    // Number is even and different from 2. Indicate composite.
	}
	// Perform 2-SPRP test
	(limbs + nbrLimbs)->x = 0;     // doesn't change actual value of pValue
	memcpy(q, limbs, (nbrLimbs + 1) * sizeof(limb));
	q[0]--;                     // q = p - 1 (p is odd, so there is no carry).
	memcpy(Mult3, q, (nbrLimbs + 1) * sizeof(q[0]));
	Mult3Len = nbrLimbs;
	DivideBigNbrByMaxPowerOf2(&ctr, Mult3, &Mult3Len);
	memcpy(TestNbr, limbs, (nbrLimbs + 1) * sizeof(limb));
	GetMontgomeryParms(nbrLimbs);
	modPowBaseInt(2, Mult3, Mult3Len, Mult1); // Mult1 = base^Mult3.
											  // If Mult1 != 1 and Mult1 = TestNbr-1, perform full test.
	if (!checkOne(Mult1, nbrLimbs) && !checkMinusOne(Mult1, nbrLimbs))
	{
		for (i = 0; i < ctr; i++)
		{               // Loop that squares number.
			modmult(Mult1, Mult1, Mult4);
			if (checkOne(Mult4, nbrLimbs) != 0)
			{  // Current value is 1 but previous value is not 1 or -1: composite
				return 2;       // Composite. Not 2-strong probable prime.
			}
			if (checkMinusOne(Mult4, nbrLimbs) != 0)
			{
				i = -1;         // Number is strong pseudoprime.
				break;
			}
			memcpy(Mult1, Mult4, nbrLimbs * sizeof(limb));
		}
		if (i == ctr)
		{
			return 1;         // Not 2-Fermat probable prime.
		}
		if (i != -1)
		{
			return 2;         // Composite. Not 2-strong probable prime.
		}
	}
	// At this point, the number is 2-SPRP, so find value of D.
	mult = 1;
	for (D = 5; ; D += 2)
	{
		int rem = getRemainder(*pValue, D);
		if (JacobiSymbol(rem, D*mult) == -1)
		{
			break;
		}
		mult = -mult;
	}
	absQ = (D + 1) >> 2;
	// Perform strong Lucas primality test on n with parameters D, P=1, Q just found.
	// Let d*2^s = n+1 where d is odd.
	// Then U_d = 0 or v_{d*2^r} = 0 for some r < s.
	// Use the following recurrences:
	// U_0 = 0, V_0 = 2.
	// U_{2k} = U_k * V_k
	// V_{2k} = (V_k)^2 - 2*Q^K
	// U_{2k+1} = (U_{2k} + V_{2k})/2
	// V_{2k+1} = (D*U_{2k} + V_{2k})/2
	// Use the following temporary variables:
	// Mult1 for Q^n, Mult3 for U, Mult4 for V, Mult2 for temporary.
	memcpy(Mult1, MontgomeryMultR1, (nbrLimbs + 1) * sizeof(limb)); // Q^0 <- 1.
	signPowQ = 1;
	memset(Mult3, 0, (nbrLimbs + 1) * sizeof(limb));                // U_0 <- 0.
	memcpy(Mult4, MontgomeryMultR1, (nbrLimbs + 1) * sizeof(limb));
	AddBigNbrMod(Mult4, Mult4, Mult4);                              // V_0 <- 2.
	CopyBigInt(expon, *pValue);
	addbigint(expon, 1);                            // expon <- n + 1.
	Temp.limbs[nbrLimbs].x = 0;
	Temp2.limbs[nbrLimbs].x = 0;
	expon.limbs[expon.nbrLimbs].x = 0;
	DivideBigNbrByMaxPowerOf2(&ctr, expon.limbs, &expon.nbrLimbs);
	for (index = expon.nbrLimbs - 1; index >= 0; index--)
	{
		int groupExp = (int)(expon.limbs[index].x);
		for (mask = 1 << (BITS_PER_GROUP - 1); mask > 0; mask >>= 1)
		{
			if (insidePowering)
			{
				// U_{2k} = U_k * V_k
				// V_{2k} = (V_k)^2 - 2*Q^K
				modmult(Mult3, Mult4, Mult3);          // U <- U * V
				modmult(Mult4, Mult4, Mult4);          // V <- V * V
				if (signPowQ > 0)
				{
					SubtBigNbrMod(Mult4, Mult1, Mult4);  // V <- V - Q^k
					SubtBigNbrMod(Mult4, Mult1, Mult4);  // V <- V - Q^k
				}
				else
				{
					AddBigNbrMod(Mult4, Mult1, Mult4);   // V <- V - Q^k
					AddBigNbrMod(Mult4, Mult1, Mult4);   // V <- V - Q^k
				}
				signPowQ = 1;                          // Indicate it is positive. 
				modmult(Mult1, Mult1, Mult1);          // Square power of Q.
			}
			if ((groupExp & mask) != 0)
			{        // Bit of exponent is equal to 1.
					 // U_{2k+1} = (U_{2k} + V_{2k})/2
					 // V_{2k+1} = (D*U_{2k} + V_{2k})/2
				Mult3[NumberLength].x = 0;
				Mult4[NumberLength].x = 0;
				AddBigNbrMod(Mult3, Mult4, Temp.limbs);
				Halve(Temp.limbs);                     // Temp <- (U + V)/2
				MultBigNbrByIntModN((int *)Mult3, D, (int *)Temp2.limbs, (int *)TestNbr, nbrLimbs);
				if (mult > 0)
				{      // D is positive
					AddBigNbrMod(Mult4, Temp2.limbs, Mult4);
				}
				else
				{      // D is negative.
					SubtBigNbrMod(Mult4, Temp2.limbs, Mult4);
				}
				Halve(Mult4);                       // V <- (V +/- U*D)/2
				memcpy(Mult3, Temp.limbs, NumberLength * sizeof(limb));
				modmultInt(Mult1, absQ, Mult1);     // Multiply power of Q by Q.
				signPowQ = -mult;                   // Attach correct sign to power.
				insidePowering = true;
			}
		}
	}
	// If U is zero, the number passes the BPSW primality test.
	if (BigNbrIsZero(Mult3))
	{
		return 0;      // Indicate number is probable prime.
	}
	for (index = 0; index < ctr; index++)
	{
		// If V is zero, the number passes the BPSW primality test.
		if (BigNbrIsZero(Mult4))
		{
			return 0;    // Indicate number is probable prime.
		}
		modmult(Mult4, Mult4, Mult4);          // V <- V * V
		if (signPowQ > 0)
		{
			SubtBigNbrMod(Mult4, Mult1, Mult4);  // V <- V - Q^k
			SubtBigNbrMod(Mult4, Mult1, Mult4);  // V <- V - Q^k
		}
		else
		{
			AddBigNbrMod(Mult4, Mult1, Mult4);   // V <- V - Q^k
			AddBigNbrMod(Mult4, Mult1, Mult4);   // V <- V - Q^k
		}
		modmult(Mult1, Mult1, Mult1);          // Square power of Q.
		signPowQ = 1;                          // Indicate it is positive.
	}
	return 3;        // Number does not pass strong Lucas test.
}


/* returns true iff value is zero*/
bool BigNbrIsZero(const limb *value)
{
	int ctr;
	for (ctr = 0; ctr < NumberLength; ctr++)
	{
		if (value[ctr].x != 0)
		{
			return false;  // Number is not zero.
		}
	}
	return true;      // Number is zero
}

/* note that ptrLimb points AFTER last valid value in limbs.
up to 3 most significant limbs are used. */
double getMantissa(const limb *ptrLimb, int nbrLimbs)
{
	double dN = (double)(ptrLimb - 1)->x;
	double dInvLimb = 1 / (double)LIMB_RANGE;
	if (nbrLimbs > 1)
	{
		dN += (double)(ptrLimb - 2)->x * dInvLimb;
	}
	if (nbrLimbs > 2)
	{
		dN += (double)(ptrLimb - 3)->x * dInvLimb * dInvLimb;
	}
	return dN;
}
