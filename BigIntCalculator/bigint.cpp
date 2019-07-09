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
//#include <iostream>
#include <stdexcept>
#include <climits>
#include <string>
#include <cassert>
#include <cmath>
#include "bignbr.h"
#include "bigint.h"
#include "factor.h"

static BigInteger Base;     // work area
static limb approxInv[MAX_LEN];
static limb adjustedArgument[MAX_LEN];
static limb arrAux[MAX_LEN];
static int bitLengthCycle[20];


/*return addend1 + addend2 (used to overload + operator) */
BigInteger BigIntAdd(const BigInteger &Addend1, const BigInteger &Addend2) {
	int ctr, nbrLimbs;
	const limb *ptrBiggerAdd, *ptrSmallerAdd;
	limb *ptrSum;
	bool A1Smaller = false;
	BigInteger Sum;  // temporary variable 

	if (Addend1.nbrLimbs < Addend2.nbrLimbs) {
		A1Smaller = true;
		/* the absolute value of addend1 is less than the absolute value of addend2.*/
	}

	else if (Addend1.nbrLimbs == Addend2.nbrLimbs) {
		for (ctr = Addend1.nbrLimbs - 1; ctr >= 0; ctr--) {
			if (Addend1.limbs[ctr].x != Addend2.limbs[ctr].x) {
				break;
			}
		}
		if (ctr >= 0 && Addend1.limbs[ctr].x < Addend2.limbs[ctr].x) {
			/* the absolute value of addend1 is less than the absolute value of addend2.*/
			A1Smaller = true;
		}
	}
	if (A1Smaller) {
		nbrLimbs = Addend1.nbrLimbs;
		ptrBiggerAdd = Addend2.limbs;
		ptrSmallerAdd = Addend1.limbs;
	}
	else {
		// the absolute value of addend1 is >= the absolute value of addend2.
		nbrLimbs = Addend2.nbrLimbs;
		ptrBiggerAdd = Addend1.limbs;
		ptrSmallerAdd = Addend2.limbs;
	}
	ptrSum = Sum.limbs;

	if (Addend1.sign == Addend2.sign)
	{             // Both addends have the same sign. Sum their absolute values.
		unsigned int carry = 0;
		for (ctr = 0; ctr < nbrLimbs; ctr++) {
			carry = (carry >> BITS_PER_GROUP) + (unsigned int)ptrBiggerAdd[ctr].x +
				(unsigned int)ptrSmallerAdd[ctr].x;
			ptrSum[ctr].x = (int)(carry & MAX_INT_NBR);
		}
		if (A1Smaller)
			nbrLimbs = Addend2.nbrLimbs;
		else
			nbrLimbs = Addend1.nbrLimbs;
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
			nbrLimbs = Addend2.nbrLimbs;
		else
			nbrLimbs = Addend1.nbrLimbs;
		for (; ctr < nbrLimbs; ctr++) {
			borrow = (borrow >> BITS_PER_INT_GROUP) + ptrBiggerAdd[ctr].x;
			ptrSum[ctr].x = borrow & MAX_INT_NBR;
		}
	}

	while (nbrLimbs > 1 && Sum.limbs[nbrLimbs - 1].x == 0) {
		nbrLimbs--;  // delete leading zeros.
	}

	Sum.nbrLimbs = nbrLimbs;
	if (A1Smaller)
		Sum.sign = Addend2.sign;  // use sign of addend with larger absolute value
	else
		Sum.sign = Addend1.sign;

	if (Sum.nbrLimbs == 1 && Sum.limbs[0].x == 0) {
		Sum.sign = SIGN_POSITIVE; // Result is zero, so sign is +ve
	}
	return Sum;
}

/* Dest = -Dest */
static void BigIntNegate (BigInteger &pDest) {
	
	if (pDest.sign == SIGN_POSITIVE && (pDest.nbrLimbs != 1 || pDest.limbs[0].x != 0))
	{
		pDest.sign = SIGN_NEGATIVE;  // pDest > 0, now < 0
	}
	else
	{
		pDest.sign = SIGN_POSITIVE; // pDest <=0, now >= 0
	}
}


/* return minuend - subtrahend (used to overload - operator) */
BigInteger BigIntSubt(const BigInteger &Minuend, const BigInteger &Subtrahend) {
	BigInteger Difference;
	static BigInteger temp;
	temp = Subtrahend;   // copy Subtrahend to temporary variable
	BigIntNegate(temp);
	Difference = Minuend + temp;
	return Difference;
}

/* returns Factor1 * Factor2 (used to overload * operator)
Factor1 will be expanded to the length of Factor2 or vice versa 
Changed design from returning error code to throw an exception, because error 
codes need to be tested for after every call, passed up to the calling routine,
tested for again and so on up to the top level. In practise this wasn't done consistently. */
BigInteger BigIntMultiply(const BigInteger &Factor1, const BigInteger &Factor2)
{
	int nbrLimbsFactor1 = Factor1.nbrLimbs;
	int nbrLimbsFactor2 = Factor2.nbrLimbs;
	int nbrLimbs;
	BigInteger Product;                  // temporary variable 
	limb Prodl[MAX_LEN * 2] = { 0 };    // temporary variable 

	if (Factor1 == 0 || Factor2 == 0) {    // one or both factors are zero.
		Product = 0;                       // Product is zero.
		return Product;
	}

	if (Factor1.nbrLimbs + Factor2.nbrLimbs > 66438 / BITS_PER_GROUP + 1)  // 2^66438 ~ 10^20000
	{
		// product out of range; throw exception
		std::string line = std::to_string(__LINE__);
		std::string mesg = "cannot multiply: product out of range ";
		mesg += __func__;
		mesg += " line ";  mesg += line;
		mesg += " in file "; mesg += __FILE__;
		throw std::range_error(mesg);
	}

	/* Factor1 will be expanded to the length of Factor2 or vice versa. It is necessary to
	override the const-ness, but the value represented does not change */
	if (nbrLimbsFactor1 < nbrLimbsFactor2) {
		memset(&((BigInteger&)Factor1).limbs[nbrLimbsFactor1], 0, (nbrLimbsFactor2 - nbrLimbsFactor1) * sizeof(limb));
		nbrLimbs = nbrLimbsFactor2;
	}
	else {
		memset(&((BigInteger&)Factor2).limbs[nbrLimbsFactor2], 0, (nbrLimbsFactor1 - nbrLimbsFactor2) * sizeof(limb));
		nbrLimbs = nbrLimbsFactor1;
	}

	multiply(&Factor1.limbs[0], &Factor2.limbs[0], Prodl, nbrLimbs, NULL);
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
	LimbsToBigInteger(Prodl, Product, nbrLimbs);
	if (Product.limbs[nbrLimbs - 1].x == 0) {
		nbrLimbs--;  // remove leading zeros
	}
	Product.nbrLimbs = nbrLimbs;

	/* set sign of product */
	if (nbrLimbs == 1 && Product.limbs[0].x == 0) {
		Product.sign = SIGN_POSITIVE;  // product is zero
	}
	else {
		if (Factor1.sign == Factor2.sign) {
			Product.sign = SIGN_POSITIVE;
		}
		else {
			Product.sign = SIGN_NEGATIVE;
		}
	}
	return Product;
}

// m *= n (used for operator overloading)
void MultBigNbrByInt(BigInteger &m, int n) {
	long long prod, carry = 0;
	int i;
	bool pos = true;

	if (n < 0) {
		pos = false;
		n = -n;   // get abs value of n
	}
	if (n == 0 || m == 0) {
		m = 0;
		return;
	}
	for (i = 0; i < m.nbrLimbs; i++) {
		prod = carry + (long long)m.limbs[i].x * n;
		carry = prod >> BITS_PER_GROUP;
		m.limbs[i].x = prod & MAX_VALUE_LIMB;
	}
	if (carry != 0) {
		m.limbs[i].x = (int)carry;
		m.nbrLimbs++;
	}
	if (!pos) {  // if n is -ve flip sign of product 
		if (m.sign = SIGN_POSITIVE)
			m.sign = SIGN_NEGATIVE;
		else
			m.sign = SIGN_POSITIVE;
	}
}

/* calculate Dividend mod Divisor (used for operator overloading) 
uses global variable Base */
BigInteger BigIntRemainder(const BigInteger &Dividend, const BigInteger &Divisor) {
	BigInteger Remainder;  
	if (Divisor == 0) {   // If divisor = 0, then remainder is the dividend.
		Remainder = Dividend;
		return Remainder;
	}
	Base = Dividend / Divisor;    // Get quotient of division.
	Base = Base*Divisor;
	return Dividend - Base;
}

/* calculate base^exponent. */
void BigIntPowerIntExp(const BigInteger &base, int exponent, BigInteger &Power) {
	int mask;
	if (base == 0) {     // base = 0 -> power = 0
		Power = 0;
		return;       // base = 0, so result is zero
	}
	if (exponent == 0) {
		Power = 1LL;
		return;
	}
	Power = 1;
	for (mask = 1 << 30; mask != 0; mask >>= 1) {
		/* first look for most significant bit of exponent */
		if ((exponent & mask) != 0) {
			for (; mask != 0; mask >>= 1) {
				Power *= Power;  
				if ((exponent & mask) != 0) {
					Power *= base; 
				}
			}
			break;  // we are finished
		}
	}
	return;
}
 
/* BigInt = e^logar.  This is the inverse function of LogBigNbr. As the result
is an integer it is not exact, and also for large values only at best the 1st
15 significant digits of the result are accurate. */
void expBigInt(BigInteger &bigInt, double logar) {
	double e = std::exp(logar);
	auto rv = fpclassify(e);  // check for overflow
	
	if (rv != FP_INFINITE && rv != FP_NAN) {
		bigInt = e;  /* note: the assignment statement uses DoubleToBigInt to
					convert floating point to BigInteger */
	}
	else 
		/* simple approach didn't work, but more complicated method below handles
		 much larger numbers. The strategy is to split the calculation into two parts.
		 One part uses Big Integers only and so avoids floating point overflow. 
		 The other part uses floating point, but the maximum number size is limited. 
		 Mathematically this method is sound, but in practise it is less accurate. */
	{
		BigInteger base=1;   
		base <<= 128;        //  new base for logs = 2^128
		double logb = logar / (std::log(2.0)*128);  // convert log to new base
		double intpf;   // integer part of log
		double fracpf = modf(logb, &intpf);  // split into log integer and fraction
		/* the desired result = base^(intpf+fracpf)
		         = base^intpf * base^fracpf */
		int intp = (int)round(intpf);         // exact conversion
		BigIntPowerIntExp(base, intp, bigInt);  /* bigInt = base^intpf*/

		/* by using a large base we make frac large, so that conversion to integer 
		doesn't introduce errors */
		double frac = std::pow(2, fracpf*128);  //  frac = base^fracpf
		BigInteger fracBI;
		fracBI = frac;    // convert to integer      
		bigInt = bigInt * fracBI;
	}
}


/* convert double dvalue to bigInt. Conversion is only accurate to about 15 
significant digits. Used for operator overloading. */
void DoubleToBigInt(BigInteger &bigInt, double dvalue) {

	if (dvalue - 0.5 > LLONG_MIN && dvalue + 0.5 < LLONG_MAX) {
		long long vv = (long long)round(dvalue); // convert directly to long long if possible
		bigInt = vv;
		return;
	}
	Znum temp;
	mpz_set_d(ZT(temp), dvalue);  // convert double to Znum
	// this method has been tested to be at least as accurate as direct conversion,
	// and obviously it's easier to use standard library functions.
	ZtoBig(bigInt, temp);        // convert Znum to BigInt
	return;
}

/* estimate natural log of BigInt. */
double logBigNbr (const BigInteger &pBigInt) {
	int nbrLimbs;
	double logar;
	nbrLimbs = pBigInt.nbrLimbs;
	if (nbrLimbs == 1) {
		logar = log((double)(pBigInt.limbs[0].x));
	}
	else {
		double value = pBigInt.limbs[nbrLimbs - 2].x +
			(double)pBigInt.limbs[nbrLimbs - 1].x * LIMB_RANGE;
		if (nbrLimbs == 2) {
			logar = log(value);
		}
		else {
			logar = log(value + (double)pBigInt.limbs[nbrLimbs - 3].x / LIMB_RANGE);
		}
		logar += (double)((nbrLimbs - 2)*BITS_PER_GROUP)*log(2);
	}
	return logar;
}

/* divide by 2, use right shift for speed */
//void BigIntDivide2(BigInteger &pArg) {
//	int nbrLimbs = pArg.nbrLimbs;
//	int ctr = nbrLimbs - 1;
//	unsigned int carry;
//	//limb *ptrLimb = &pArg->limbs[ctr];
//	limb *ptrLimb = pArg.limbs;
//	carry = 0;
//	for (; ctr >= 0; ctr--)
//	{
//		carry = (carry << BITS_PER_GROUP) + (unsigned int)ptrLimb[ctr].x;
//		ptrLimb[ctr].x = (int)(carry >> 1);
//		carry &= 1;
//	}
//	if (nbrLimbs > 1 && pArg.limbs[nbrLimbs - 1].x == 0)
//	{     // Most significant limb is zero, so reduce size by one limb.
//		pArg.nbrLimbs--;
//	}
//}

/* arg = arg*2^power. Throw exception if product is too large */
//static void BigIntMutiplyPower2(BigInteger &pArg, int power2)
//{
//	int ctr;
//	int nbrLimbs = pArg.nbrLimbs;
//	limb *ptrLimbs = pArg.limbs;
//
//	for (; power2 > 0; power2--) {
//		/*each time round the loop multiplies arg by 2 */
//		unsigned int carry = 0;
//		for (ctr = 0; ctr < nbrLimbs; ctr++)
//		{
//			carry += (unsigned int)ptrLimbs[ctr].x << 1;
//			ptrLimbs[ctr].x = (int)(carry & MAX_VALUE_LIMB);
//			carry >>= BITS_PER_GROUP;
//		}
//		if (carry != 0)
//		{
//			ptrLimbs[ctr].x = (int)carry;
//			nbrLimbs++;
//			if (nbrLimbs > MAX_LEN) {
//				std::string line = std::to_string(__LINE__);
//				std::string mesg = "number too big : cannot do multiplication : ";
//				mesg += __func__;
//				mesg += " line ";  mesg += line;
//				mesg += " in file "; mesg += __FILE__;
//				throw std::range_error(mesg);
//			}
//		}
//	}
//	pArg.nbrLimbs = nbrLimbs;
//}

/* return true if Nbr1 == Nbr2 (used for operator overloading)*/
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

/* return true if Nbr1 < Nbr2 (used for operator overloading) */
bool TestBigNbrLess(const BigInteger &Nbr1, const BigInteger &Nbr2) {
	int ctr;
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

/* calculate GCD of arg1 & arg2*/
//void BigIntGcd(const BigInteger &Arg1, const BigInteger &Arg2, BigInteger &Result)
//{
//	int nbrLimbs1 = Arg1.nbrLimbs;
//	int nbrLimbs2 = Arg2.nbrLimbs;
//	int power2;
//	static BigInteger Power;
//	if (Arg1 == 0)
//	{               // First argument is zero, so the GCD is second argument.
//		Result = Arg2;    //CopyBigInt(Result, pArg2);
//		return;
//	}
//	if (Arg2 == 0)
//	{               // Second argument is zero, so the GCD is first argument.
//		Result = Arg1;		//CopyBigInt(Result, pArg1);
//		return;
//	}
//	// Reuse Base and Power temporary variables.
//	Base = Arg1;     // CopyBigInt(Base, Arg1);   
//	Power = Arg2;   //  CopyBigInt(Power, Arg2); 
//	Base.sign = SIGN_POSITIVE;
//	Power.sign = SIGN_POSITIVE;
//	power2 = 0;
//	while (((Base.limbs[0].x | Power.limbs[0].x) & 1) == 0)
//	{  // Both values are even
//		BigIntDivide2(Base);
//		BigIntDivide2(Power);
//		power2++;
//	}
//
//	while (Base != Power)	//while (TestBigNbrEqual(Base, Power) == 0)
//	{    // Main GCD loop.
//		if (Base.isEven()) {     // Number is even. Divide it by 2.
//			BigIntDivide2(Base);
//			continue;
//		}
//		if (Power.isEven())  {     // Number is even. Divide it by 2.
//			BigIntDivide2(Power);
//			continue;
//		}
//		Result = Base - Power; // BigIntSubt(Base, Power, Result);
//		if (Result >= 0) {
//			Base = Result; // CopyBigInt(Base, Result);
//			BigIntDivide2(Base);
//		}
//		else {
//			Power = Result; // CopyBigInt(Power, Result);
//			Power.sign = SIGN_POSITIVE;
//			BigIntDivide2(Power);
//		}
//	}
//	Result = Base; // CopyBigInt(Result, Base);
//	BigIntMutiplyPower2(Result, power2); /* Result *= 2^power     */
//}

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

/* i = (i-subt)/divisor. Assume (without checking) that divisor > 0.
Does not appear to handle -ve divisor 
used for operator overloading */
void subtractdivide(BigInteger &i, int subt, int divisor)
{
	int nbrLimbs = i.nbrLimbs;
	int remainder = 0;
	
#if 0
	char *ptrOutput = output;
	*ptrOutput++ = '2';
	*ptrOutput++ = '(';
	int2dec(&ptrOutput, i->sign);
	*ptrOutput++ = ',';
	*ptrOutput++ = ' ';
	int2dec(&ptrOutput, i->nbrLimbs);
	*ptrOutput++ = ';';
	*ptrOutput++ = ' ';
	int2dec(&ptrOutput, i->limbs[0].x);
	*ptrOutput++ = ',';
	*ptrOutput++ = ' ';
	int2dec(&ptrOutput, i->limbs[1].x);
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
	if ((unsigned int)i->limbs[0].x >= LIMB_RANGE)
	{
		remainder = 1;
	}
#endif
	/* if subt is not zero subtract it from i*/
	if (subt > 0) {
		if (i >= 0) {       // Subtract subt from absolute value.
			subtFromAbsValue(i.limbs, &nbrLimbs, subt);
		}
		else {               // Add subt to absolute value.
			addToAbsValue(i.limbs, &nbrLimbs, subt);
		}
	}
	else if (subt < 0) {  // subt < 0
		if (i >= 0) {    // Subtract subt from absolute value.
			addToAbsValue(i.limbs, &nbrLimbs, -subt);
		}
		else {               // Add subt to absolute value.
			subtFromAbsValue(i.limbs, &nbrLimbs, -subt);
		}
	}

	// Divide number by divisor.
	i = BigIntDivideInt(i, divisor);
}

/* calculate BigInt modulo divisor (used for operator overloading) */
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

/* result += addend (used for operator overloading) */
void addbigint(BigInteger &Result, int addend) {
	int nbrLimbs = Result.nbrLimbs;
	limb *pResultLimbs = Result.limbs;
	auto sign = Result.sign;

	if (addend < 0) {
		// reverse signs of addend and result
		addend = -addend;
		if (sign == SIGN_POSITIVE) 	{
			sign = SIGN_NEGATIVE;  // addend and result have opposite signs
		}
		else {
			sign = SIGN_POSITIVE;   // addend and result have same sign
		}
	}

	if (sign == SIGN_POSITIVE) {   // Add abs(addend) to absolute value of pResult.
		addToAbsValue(pResultLimbs, &nbrLimbs, addend);
	}
	else {  // Subtract abs(addend) from absolute value of pResult.
		if (nbrLimbs == 1) 	{
			pResultLimbs[0].x -= addend;
			if (pResultLimbs[0].x < 0) 	{
				pResultLimbs[0].x = -pResultLimbs[0].x;  // reverse sign of result
				BigIntNegate(Result);
			}
		}
		else {     // More than one limb.
			subtFromAbsValue(pResultLimbs, &nbrLimbs, addend);
		}
	}
	Result.nbrLimbs = nbrLimbs;
}


/* creates a BigInteger from a list of values.   
uses global value NumberLength for number of ints. 1st entry in list is number of values
that follow */
//void IntsToBigInteger(/*@in@*/const int *ptrValues, /*@out@*/BigInteger &bigint)
//{
//	if (NumberLength > MAX_LEN || NumberLength < 0 || ptrValues[0] > MAX_LEN) {
//		std::string line = std::to_string(__LINE__);
//		std::string mesg = "number too big : cannot convert to BigInteger: ";
//		mesg += __func__;
//		mesg += " line ";  mesg += line;
//		mesg += " in file "; mesg += __FILE__;
//		throw std::range_error(mesg);
//	}
//	limb *destLimb = bigint.limbs;
//	bigint.sign = SIGN_POSITIVE;
//	if (NumberLength == 1) {
//		destLimb[0].x = ptrValues[1];
//		bigint.nbrLimbs = 1;
//	}
//	else {
//		memcpy(destLimb, ptrValues+1, ptrValues[0] * sizeof(ptrValues[0]));
//		bigint.nbrLimbs = ptrValues[0];
//		if (NumberLength > ptrValues[0])  // clear most significant limbs to zero if required
//			memset(destLimb + ptrValues[0], 0, (NumberLength - ptrValues[0]) * sizeof(ptrValues[0]));
//	}
//}

/* uses global value NumberLength for starting value for number of limbs. */
//static int getNbrLimbs(const limb *bigNbr)
//{
//	//const limb *ptrLimb = bigNbr + NumberLength;
//	auto ix = NumberLength;
//	while (ix > 0)
//	{
//		if (bigNbr[ix-1].x != 0)
//		{
//			return (int)(ix);
//		}
//		ix--;   // reduce length because most significant limb is zero
//	}
//	return 1;  // BigNbr is zero
//}

/* creates a list of values from a BigInteger, 1st entry in list is number of 
values that follow. Also uses global value NumberLength for number of ints. */
//void BigIntegerToInts(/*@out@*/int *ptrValues, /*@in@*/const BigInteger &bigint) {
//	const limb *srcLimb = bigint.limbs;
//	if (NumberLength == 1) {
//		ptrValues[0] = 1;
//		ptrValues[1] = srcLimb[0].x;
//	}
//	else {
//		int nbrLimbs = getNbrLimbs(bigint.limbs); //nbrLimbs <= NumberLength
//		assert(nbrLimbs == bigint.nbrLimbs);
//		ptrValues[0] = nbrLimbs;
//		memcpy(ptrValues + 1, srcLimb, nbrLimbs * sizeof(ptrValues[0]));
//	}
//}

/* convert limbs to BigInteger. */
void LimbsToBigInteger(/*@in@*/const limb *ptrValues, 
	/*@out@*/BigInteger &bigint, int NumLen) {
	if (NumLen > MAX_LEN || NumLen < 0 ) {
		std::string line = std::to_string(__LINE__);
		std::string mesg = "number too big : cannot convert to BigInteger: ";
		mesg += __func__;
		mesg += " line ";  mesg += line;
		mesg += " in file "; mesg += __FILE__;
		throw std::range_error(mesg);
	}
	if (NumLen == 1) {
		bigint.limbs[0].x = ptrValues[0].x;
		bigint.nbrLimbs = 1;
	}
	else { 
		memcpy(bigint.limbs, ptrValues, NumLen * sizeof(limb));

		int nbrLimbs;   // remove any leading zeros
		for (nbrLimbs = NumLen-1; nbrLimbs > 1; nbrLimbs--) {
			if (ptrValues[nbrLimbs].x != 0) {
				break;
			}
		}
		bigint.nbrLimbs = nbrLimbs+1;
	}
	bigint.sign = SIGN_POSITIVE;
}

/* Convert BigInteger to limbs. uses global value NumberLength for number of limbs. */
void BigIntegerToLimbs(/*@out@*/limb *ptrValues, 
	/*@in@*/const BigInteger &bigint, int NumLen)
{
	if (NumLen == 1) {
		ptrValues[0].x = bigint.limbs[0].x;
		ptrValues[1].x = 0;
	}
	else {
		int nbrLimbs = bigint.nbrLimbs; // use lesser of bigint.nbrLimbs & NumLen
		if (nbrLimbs >= NumLen) {
			memcpy(ptrValues, bigint.limbs, NumLen * sizeof(limb));
			ptrValues[NumLen].x = 0;
		}
		else {
			memcpy(ptrValues, bigint.limbs, nbrLimbs * sizeof(limb));
			/* set any extra limbs to zero */
			memset(ptrValues + nbrLimbs, 0, (NumLen - nbrLimbs) * sizeof(limb));
		}
	}
}



//void DivideBigNbrByMaxPowerOf2(int *pShRight, limb *number, int *pNbrLimbs)
//{
//	int power2 = 0;
//	long long mask;
//	int index, index2,  shRg;
//	int nbrLimbs = *pNbrLimbs;
//	// Start from least significant limb (number zero).
//	for (index = 0; index < nbrLimbs; index++)
//	{
//		if (number[index].x != 0)
//		{
//			break;
//		}
//		power2 += BITS_PER_GROUP;
//	}
//	for (mask = 0x1; mask <= MAX_VALUE_LIMB; mask *= 2)
//	{
//		if ((number[index].x & mask) != 0)
//		{
//			break;
//		}
//		power2++;
//	}
//	// Divide number by this power.
//	shRg = power2 % BITS_PER_GROUP; // Shift right bit counter
//	if ((number[nbrLimbs - 1].x & (-(1 << shRg))) != 0)
//	{   // Most significant bits set.
//		*pNbrLimbs = nbrLimbs - index;
//	}
//	else
//	{   // Most significant bits not set.
//		*pNbrLimbs = nbrLimbs - index - 1;
//	}
//	// Move number shRg bits to the right.
//	mask = (1 << shRg) - 1;
//	for (index2 = index; index2 < nbrLimbs; index2++)
//	{
//		number[index2].x = ((number[index2].x >> shRg) |
//			(number[index2 + 1].x << (BITS_PER_GROUP - shRg))) &
//			MAX_VALUE_LIMB;
//	}
//	if (index > 0)
//	{   // Move limbs to final position.
//		memmove(number, &number[index], (nbrLimbs - index) * sizeof(limb));
//	}
//	*pShRight = power2;
//}



/* convert Znum to BigInteger. Returns false if number is too big to convert.
this function is also used to overload the assignment operator */
bool ZtoBig(BigInteger &number, Znum numberZ) {
	number.nbrLimbs = 0;
	bool neg = false;
	Znum remainder;

	if (numberZ < 0) {
		neg = true;
		numberZ = -numberZ;  // make numberZ +ve
	}
	int i = 0;
	while (numberZ > 0) {
		//mpz_fdiv_qr_ui(ZT(quot), ZT(remainder), ZT(numberZ), LIMB_RANGE);
		/* calculating quotient and remainder separately turns
		out to be faster */
		mpz_fdiv_r_2exp(ZT(remainder), ZT(numberZ), BITS_PER_GROUP);
		number.limbs[i].x = (int)MulPrToLong(remainder);
		mpz_fdiv_q_2exp(ZT(numberZ), ZT(numberZ), BITS_PER_GROUP);
		i++;
		if (i >= MAX_LEN) {
			return false;   // number too big to convert.
		}
	}
	number.nbrLimbs = i;
	if (neg) {
		number.sign = SIGN_NEGATIVE;
		numberZ = -numberZ;  // put back original value in numberZ
	}
	else
		number.sign = SIGN_POSITIVE;

	return true;
}

/* convert BigInteger to Znum */
void BigtoZ(Znum &numberZ, const BigInteger &number) {

	numberZ = 0;
	for (int i = number.nbrLimbs - 1; i >= 0; i--) {
		//numberZ *= LIMB_RANGE;
		mpz_mul_2exp(ZT(numberZ), ZT(numberZ), BITS_PER_GROUP);  // shift numberZ left
		numberZ += number.limbs[i].x;      // add next limb
	}
	if (number.sign == SIGN_NEGATIVE)
		numberZ = -numberZ;
}

/* convert num to long long. If num > max, truncate it (right shift)
and exp is > 0 and represents the number of discarded bits. */
long long BigToLL(const BigInteger &num, int &exp) {
	BigInteger temp = num;
	if (num.nbrLimbs == 1) {
		exp = 0;
		if (num.sign == SIGN_POSITIVE)
			return num.limbs[0].x;
		else
			return -num.limbs[0].x;
	}
	exp = num.bitLength() - 63;                // number of bits to truncate
	if (exp < 0)
		exp = 0;
	temp >>= exp;
	long long result = temp.lldata();
	return result;
}


/* inverse of BigToLL. set num = LL *2^exp */
void LLToBig(BigInteger &num, long long LL, int exp) {
	num = LL;
	if (exp == 0) { 
		return;
	}
	assert(exp > 0);
	num <<= exp;
}

/* shift first left by the number of bits specified in shiftCtr. A -ve value
in shiftCtr causes a right shift.
Right Shifts simulate 2s complement arithmetic right shift.
Mathematically, the shift result is equivalent to result = first * 2^shiftCtr,
whether ShiftCtr is +ve or -ve. */
void shiftBI(const BigInteger &first, const int shiftCtr, BigInteger &result)
{
	int delta, rem, ctr;
	long long prevLimb, curLimb;
	int ptrDest, ptrSrc;
	bool shiftleft = true;
	if (shiftCtr > 0) {
		delta = shiftCtr / BITS_PER_GROUP;
		rem = shiftCtr % BITS_PER_GROUP;
	}
	else {
		delta = (-shiftCtr) / BITS_PER_GROUP;
		rem = (-shiftCtr) % BITS_PER_GROUP;
		shiftleft = false;
	}
	int nbrLimbs = first.nbrLimbs;

	if (shiftleft) {     // Perform shift left.

		if ((first.nbrLimbs + delta) >= MAX_LEN) {
			// Shift too much to the left; would cause overflow
			std::string line = std::to_string(__LINE__);
			std::string mesg = "cannot shift left: result out of range ";
			mesg += __func__;
			mesg += " line ";  mesg += line;
			mesg += " in file "; mesg += __FILE__;
			throw std::range_error(mesg);
		}


		result.nbrLimbs = first.nbrLimbs + delta;
		result.sign = first.sign;
		prevLimb = 0;
		ptrSrc = nbrLimbs - 1;
		ptrDest = nbrLimbs + delta;

		for (ctr = nbrLimbs - 1; ctr >= 0; ctr--)
		{  // Process starting from most significant limb.
			curLimb = first.limbs[ctr].x;
			result.limbs[ptrDest].x = ((curLimb >> (BITS_PER_GROUP - rem))
				| (prevLimb << rem)) & MAX_INT_NBR;
			ptrDest--;
			prevLimb = curLimb;
		}

		result.limbs[ptrDest].x = (prevLimb << rem) & MAX_INT_NBR;
		if (delta > 0) {
			memset(result.limbs, 0, delta * sizeof(limb));
		}
		//result.nbrLimbs += delta;
		if (result.limbs[result.nbrLimbs].x != 0) {
			result.nbrLimbs++;
		}
	}

	else {     // Perform shift right.
		int isNegative = 0;
		if (shiftCtr > first.nbrLimbs * BITS_PER_GROUP)
		{   // Shift too much to the right. Result is zero or -1.
			if (first.sign == SIGN_POSITIVE)
				result = 0;
			else
				result = -1;

			return;
		}

		result = first;
		if (first.sign == SIGN_NEGATIVE)
		{   // If it is negative, add 1, perform shift right, and finally subtract 1 from result.
			isNegative = 1;
			result++;     
		}
		// Shift right the absolute value.
		result.limbs[nbrLimbs].x = 0;
		curLimb = result.limbs[delta].x;
		ptrDest = 0;

		for (ctr = delta; ctr <= nbrLimbs; ctr++)
		{  // Process starting from least significant limb.
			prevLimb = curLimb;
			curLimb = result.limbs[ctr + 1].x;
			result.limbs[ptrDest++].x = ((prevLimb >> rem)
				| (curLimb << (BITS_PER_GROUP - rem))) & MAX_INT_NBR;
		}

		result.nbrLimbs -= delta;
		if (result.nbrLimbs == 0 || result.limbs[result.nbrLimbs].x != 0) {
			result.nbrLimbs++;
		}
		if (isNegative) {    // Adjust negative number.
			result--;        
		}
	}
	return;
}


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
	if (Divisor == 0) {  // Indicate error if divisor is zero.
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
		Quotient = Dividend; 
		Quotient /= Divisor.limbs[0].x;   
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
		bitLength = nbrLimbs * BITS_PER_GROUP;
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

/* Return Dividend/Divisor. Used for operator overloading */
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

	while (len > 1 && Quotient.limbs[len - 1].x == 0)
		len--;  // remove any leading zeros
	Quotient.nbrLimbs = len;
	return Quotient;
}