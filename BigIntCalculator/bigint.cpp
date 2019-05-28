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
#include <stdexcept>
#include <climits>
#include <string>
#include <cassert>
#include <cmath>
#include "bignbr.h"
#include "factor.h"

extern bool *primeFlags;
extern unsigned long long *primeList;
extern unsigned int prime_list_count;
extern unsigned long long int primeListMax;
bool isPrime2(unsigned __int64 num);
BigInteger Base;

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
	BigInteger Product;    // temporary variable 
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

	multiply(&Factor1.limbs[0], &Factor2.limbs[0], Prodl, nbrLimbs, &nbrLimbs);
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
 
/* BigInt = e^logar
throw exception if result would be to large for a BigInteger */
//void expBigInt(BigInteger &bigInt, double logar)
//{
//	int mostSignificantLimb;
//	logar /= log(2);  // convert log to base 2
//	/* calculate required value for most significnt limb, initially in floating point */
//	bigInt.nbrLimbs = (int)floor(logar / BITS_PER_GROUP);  // nbrLimbs will be increased as required
//	double value = round(exp((logar - BITS_PER_GROUP*bigInt.nbrLimbs) * log(2)));
//
//	/* convert double to BigInteger */
//	bigInt.sign = SIGN_POSITIVE;
//	
//	mostSignificantLimb = (int) value;
//	if (mostSignificantLimb == LIMB_RANGE)
//	{
//		mostSignificantLimb = 1;
//		bigInt.nbrLimbs++;
//	}
//	bigInt.nbrLimbs++;
//	if (bigInt.nbrLimbs > MAX_LEN) {
//		std::string line = std::to_string(__LINE__);
//		std::string mesg = "number too big : cannot expand BigInteger: ";
//		mesg += __func__;
//		mesg += " line ";  mesg += line;
//		mesg += " in file "; mesg += __FILE__;
//		throw std::range_error(mesg);
//	}
//	memset(bigInt.limbs, 0, bigInt.nbrLimbs * sizeof(limb));
//	bigInt.limbs[bigInt.nbrLimbs - 1].x = mostSignificantLimb;
//}


/* convert double dvalue to bigInt. Conversion is only accurate to about 15 significant digits. */
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
//double logBigNbr (const BigInteger &pBigInt) {
//	int nbrLimbs;
//	double logar;
//	nbrLimbs = pBigInt.nbrLimbs;
//	if (nbrLimbs == 1) {
//		logar = log((double)(pBigInt.limbs[0].x));
//	}
//	else {
//		double value = pBigInt.limbs[nbrLimbs - 2].x +
//			(double)pBigInt.limbs[nbrLimbs - 1].x * LIMB_RANGE;
//		if (nbrLimbs == 2) {
//			logar = log(value);
//		}
//		else {
//			logar = log(value + (double)pBigInt.limbs[nbrLimbs - 3].x / LIMB_RANGE);
//		}
//		logar += (double)((nbrLimbs - 2)*BITS_PER_GROUP)*log(2);
//	}
//	return logar;
//}
double logBigNbr(const Znum &BigInt) {
	double BigId;
#ifdef __MPIR_VERSION
	long BiExp;   // changed for MPIR version 3.0.0
#else
    long BiExp;
#endif
	BigId = mpz_get_d_2exp(&BiExp, ZT(BigInt)); // BigId * 2^BiExp = BigInt 
	double logval = log(BigId)  + BiExp *  log(2);
	return logval;
}

/* estimate natural log of BigNbr. Only the most significant 62 bits are
taken into account because floating point numbers have limited accuracy anyway. */
//double logLimbs(const limb *pBigNbr, int nbrLimbs)
//{
//	double logar;
//	if (nbrLimbs > 1)
//	{
//		logar = log( (double)pBigNbr[nbrLimbs - 2].x +
//				(double)pBigNbr[nbrLimbs - 1].x * LIMB_RANGE) +
//			(double)(nbrLimbs - 2)*log((double)LIMB_RANGE);
//	}
//	else
//	{
//		logar = log((double)pBigNbr[nbrLimbs - 1].x) +
//			(double)(nbrLimbs - 1)*log((double)LIMB_RANGE);
//	}
//	return logar;
//}

/* calculate base^exponent. */
//void BigIntPowerIntExp(const BigInteger &pBase, int exponent, BigInteger &Power) {
//	int mask;
//	if (pBase == 0) {     // Base = 0 -> power = 0
//		Power = 0;
//		return; // base = 0, so result is zero
//	}
//	Power = 1;
//	for (mask = 1 << 30; mask != 0; mask >>= 1) {
//		if ((exponent & mask) != 0) {
//			for (; mask != 0; mask >>= 1) {
//				Power *= Power;// BigIntMultiply(Power, Power, Power);
//				if ((exponent & mask) != 0) {
//					Power *= pBase; //BigIntMultiply(Power, Base, Power);
//				}
//			}
//			break;
//		}
//	}
//	return;
//}
void BigIntPowerIntExp(const Znum &Base, int exponent, Znum &Power) {
	mpz_pow_ui(ZT(Power), ZT(Base), exponent);
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
				BigIntNegate(Result);
			}
		}
		else {     // More than one limb.
			subtFromAbsValue(pResultLimbs, &nbrLimbs, addend);
		}
	}
	Result.nbrLimbs = nbrLimbs;
}

/* returns nbrMod^Expon%currentPrime. 
overflow could occur if currentPrime > 2^31 
the alternative is to use intDoubleModPow */
static long long intModPow(long long NbrMod, long long Expon, long long currentPrime)
{
	unsigned long long power = 1;
	unsigned long long square = (unsigned long long)NbrMod;
	while (Expon != 0)
	{
		if ((Expon & 1) == 1)
		{
			power = (power * square) % (unsigned long long)currentPrime;
		}
		square = (square * square) % (unsigned long long)currentPrime;
		Expon >>= 1;
	}
	return (long long)power;
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

/* convert limbs to BigInteger. uses global value NumberLength for number of limbs. */
void LimbsToBigInteger(/*@in@*/const limb *ptrValues, 
	/*@out@*/BigInteger &bigint, int NumberLength) {
	if (NumberLength > MAX_LEN || NumberLength < 0 ) {
		std::string line = std::to_string(__LINE__);
		std::string mesg = "number too big : cannot convert to BigInteger: ";
		mesg += __func__;
		mesg += " line ";  mesg += line;
		mesg += " in file "; mesg += __FILE__;
		throw std::range_error(mesg);
	}
	if (NumberLength == 1) {
		bigint.limbs[0].x = ptrValues[0].x;
		bigint.nbrLimbs = 1;
	}
	else { 
		memcpy(bigint.limbs, ptrValues, NumberLength * sizeof(limb));

		int nbrLimbs;   // remove any leading zeros
		for (nbrLimbs = NumberLength-1; nbrLimbs > 1; nbrLimbs--) {
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
	/*@in@*/const BigInteger &bigint, int NumberLength)
{
	if (NumberLength == 1) {
		ptrValues[0].x = bigint.limbs[0].x;
		ptrValues[1].x = 0;
	}
	else {
		int nbrLimbs = bigint.nbrLimbs; // use lesser of bigint.nbrLimbs & NumberLength
		if (nbrLimbs >= NumberLength) {
			memcpy(ptrValues, bigint.limbs, NumberLength * sizeof(limb));
			ptrValues[NumberLength].x = 0;
		}
		else {
			memcpy(ptrValues, bigint.limbs, nbrLimbs * sizeof(limb));
			/* set any extra limbs to zero */
			memset(ptrValues + nbrLimbs, 0, (NumberLength - nbrLimbs) * sizeof(limb));
		}
	}
}
/* number = numberZ*/
int ZtoLimbs(limb *number, Znum numberZ, int NumberLength) {
// note: numberZ is a copy of the original. Its value is changed
	bool neg = false;
	Znum quot, remainder;

	if (numberZ < 0) {
		neg = true;
		numberZ = -numberZ;  // make numberZ +ve
	}
	int i = 0;
	while (numberZ > 0) {
		if (i >= MAX_LEN || NumberLength > MAX_LEN) {
			// number too big to convert.
			std::string line = std::to_string(__LINE__);
			std::string mesg = "number too big : cannot convert to limbs: ";
			mesg += __func__;
			mesg += " line ";  mesg += line;
			mesg += " in file "; mesg += __FILE__;
			throw std::range_error(mesg);
		}
		//mpz_fdiv_qr_ui(ZT(quot), ZT(remainder), ZT(numberZ), LIMB_RANGE);
		/* calculating quotient and remainder separately turns
		out to be faster */
		mpz_fdiv_r_2exp(ZT(remainder), ZT(numberZ), BITS_PER_GROUP);
		number[i].x = (int)MulPrToLong(remainder);
		mpz_fdiv_q_2exp(ZT(numberZ), ZT(numberZ), BITS_PER_GROUP);
		//numberZ = quot;
		i++;
	}
	if (i < NumberLength) {
		/* set any extra limbs to zero */
		memset(number + i, 0, (NumberLength - i) * sizeof(limb));
	}
	if (neg) {
		ChSignBigNbr((int *)number, i + 1);
	}
	return i;
}
int ZtoBigNbr(int number[], Znum numberZ) {
	// note: numberZ is a copy of the original. Its value is changed
	bool neg = false;
	Znum quot, remainder;

	if (numberZ < 0) {
		neg = true;
		numberZ = -numberZ;  // make numberZ +ve
	}
	int i = 0;
	while (numberZ > 0) {
		mpz_fdiv_qr_ui(ZT(quot), ZT(remainder), ZT(numberZ), LIMB_RANGE);
		//number[i] = (int)MulPrToLong(remainder);
		number[i] = (int)mpz_get_si(ZT(remainder));  // faster?? - no possibility of overflow here
		numberZ = quot;
		i++;
	}

	if (neg) {
		ChSignBigNbr(number, i);
	}
	return i;
}

void LimbstoZ(const limb *number, Znum &numberZ, int NumberLength) {
	numberZ = 0;
	for (int i = NumberLength - 1; i >= 0; i--) {
		mpz_mul_2exp(ZT(numberZ), ZT(numberZ), BITS_PER_GROUP);  // shift numberZ left
		numberZ += number[i].x;      // add next limb
	}
}
//void UncompressIntLimbs(/*@in@*/const int *ptrValues, /*@out@*/limb *bigNbr, int nbrLen)
//{
//	int nbrLimbs = *ptrValues;
//	memcpy(bigNbr, ptrValues + 1, nbrLimbs * sizeof(limb));
//	/* if nbrLen > nbrLimbs set extra limbs to 0*/
//	memset(bigNbr + nbrLimbs, 0, (nbrLen - nbrLimbs) * sizeof(limb));
//}

//void CompressIntLimbs(/*@out@*/int *ptrValues, /*@in@*/const limb *bigint, int nbrLen)
//{
//	int nbrLimbs;
//	memcpy(ptrValues + 1, bigint, (nbrLen - 1) * sizeof(limb));
//	for (nbrLimbs = nbrLen - 1; nbrLimbs > 1; nbrLimbs--)
//	{
//		if (ptrValues[nbrLimbs] != 0)
//		{
//			break;
//		}
//	}
//	*ptrValues = nbrLimbs;
//}

// This routine checks whether the number factor is a perfect power. 
// If it is not, it returns one. If it is a perfect power, it returns the  
// exponent and  the base such that base^exponent = factor.
long long PowerCheck(const Znum &factor, Znum &Base, long long upperBound) {
	/* upperbound is the largest number already tested as a factor by trial division
	i.e. factor has no factors < upperBound. This can be used to put a much
	smaller limit on maxExpon (max about 2000) */

	/* upperBound^maxExpon ≈ factor */
	unsigned long long maxExpon = (unsigned long long) (ceil(logBigNbr(factor) / log(upperBound)));
										  
	int h;
	long long modulus, Exponent;
	unsigned long long primesLength, j;
	int prime2310x1[] =
	{ 2311, 4621, 9241, 11551, 18481, 25411, 32341, 34651, 43891, 50821 };
	// Primes of the form 2310x+1.
	bool expon2 = true, expon3 = true, expon5 = true;
	bool expon7 = true, expon11 = true;
	Znum Zmod, Root;
	std::vector<bool>ProcessExpon(maxExpon+1);

	for (h = 0; h < sizeof(prime2310x1) / sizeof(prime2310x1[0]); h++) {
		int testprime = prime2310x1[h];
		// Zmod = mod = Bigint%testprime
		auto mod = mpz_mod_ui(ZT(Zmod),ZT(factor), testprime); // getRemainder(factor, testprime);
		if (expon2 && intModPow(mod, testprime / 2, testprime) > 1) {
			expon2 = false;
		}
		if (expon3 && intModPow(mod, testprime / 3, testprime) > 1) {
			expon3 = false;
		}
		if (expon5 && intModPow(mod, testprime / 5, testprime) > 1) {
			expon5 = false;
		}
		if (expon7 && intModPow(mod, testprime / 7, testprime) > 1) {
			expon7 = false;
		}
		if (expon11 && intModPow(mod, testprime / 11, testprime) > 1) {
			expon11 = false;
		}
	}

	primesLength = 2 * maxExpon + 3;
	if (primesLength > primeListMax) {
		std::string line = std::to_string(__LINE__);
		std::string mesg = "number too big : cannot generate prime list. function : ";
		mesg += __func__;
		mesg += " line ";  mesg += line;
		mesg += " in file "; mesg += __FILE__;
		throw std::range_error(mesg);
	}

	for (h = 2; h <= maxExpon; h++) {
		ProcessExpon[h] = true;
	}

	for (size_t ix=5, h = primeList[ix]; h < primesLength/2; ix++, h = primeList[ix]) {
		int processed = 0;
		for (j = 2 * h + 1; j < primesLength; j += 2 * h) {
			if (isPrime2(j)) {
				modulus = mpz_mod_ui(ZT(Zmod),ZT(factor), j); // getRemainder(factor, j);
				if (intModPow(modulus, j / h, j) > 1) {
					for (j = h; j <= maxExpon; j += h) {
						ProcessExpon[j] = false;
					}
					break;
				}
			}
			if (++processed > 10) {
				break;
			}
		}
	}

	/* check possible exponent values */
	for (Exponent = maxExpon; Exponent >= 2; Exponent--) {
		if (Exponent % 2 == 0 && !expon2) {
			continue; // Not a square
		}
		if (Exponent % 3 == 0 && !expon3) {
			continue; // Not a cube
		}
		if (Exponent % 5 == 0 && !expon5) {
			continue; // Not a fifth power
		}
		if (Exponent % 7 == 0 && !expon7) {
			continue; // Not a 7th power
		}
		if (Exponent % 11 == 0 && !expon11) {
			continue; // Not an 11th power
		}
		if (!ProcessExpon[Exponent]) {
			continue;
		}
		if (mpz_root(ZT(Root), ZT(factor), Exponent) != 0) {
			Base = Root;   // factor is a perfect power
			return Exponent;
		}
	}

	Base = factor;   // not perfect power
	return 1;
}

/* return true if value is 1 (mod p)*/
//bool checkOne(const limb *value, int nbrLimbs)
//{
//	int idx;
//	for (idx = 0; idx < nbrLimbs; idx++)
//	{
//		if ((value++)->x != MontgomeryMultR1[idx].x)
//		{
//			return false;    // Go out if value is not 1 (mod p)
//		}
//	}
//	return true; // value = 1 (mod p)
//}

/* return true if value is -1 (mod p) */
//bool checkMinusOne(const limb *value, int nbrLimbs)
//{
//	int idx;
//	unsigned int carry;
//	carry = 0;
//	for (idx = 0; idx < nbrLimbs; idx++)
//	{
//		carry += (unsigned int)(value++)->x + (unsigned int)MontgomeryMultR1[idx].x;
//		if ((carry & MAX_VALUE_LIMB) != (unsigned int)TestNbr[idx].x)
//		{
//			return false;    // Go out if value is not -1 (mod p)
//		}
//		carry >>= BITS_PER_GROUP;
//	}
//	return true;    // value is -1 (mod p)
//}

// Find power of 2 that divides the number.
// output: pNbrLimbs = pointer to number of limbs
//         pShRight = pointer to power of 2.

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

// Find largest power of 2 that divides the number. Divide number
// by that power and return value of power in ShRight.
void DivideBigNbrByMaxPowerOf2(int &ShRight, Znum &number) {
	Znum two = 2;
	auto shift = mpz_remove(ZT(number), ZT(number), ZT(two));
	ShRight = (int)shift;
}

// Calculate Jacobi symbol by following algorithm 2.3.5 of C&P book.
long long JacobiSymbol(long long upper, long long lower)
{
	long long tmp;
	long long a = upper % lower;
	long long m = lower;
	long long t = 1;
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
//static void Halve(limb *pValue)
//{
//	if ((pValue[0].x & 1) == 0)
//	{    // Number to halve is even. Divide by 2.
//		DivBigNbrByInt((int *)pValue, 2, (int *)pValue, NumberLength);
//	}
//	else
//	{    // Number to halve is odd. Add modulus and then divide by 2.
//		AddBigNbr((int *)pValue, (int *)TestNbr, (int *)pValue, NumberLength + 1);
//		DivBigNbrByInt((int *)pValue, 2, (int *)pValue, NumberLength + 1);
//	}
//}
//static void Halve(Znum & Value, const Znum & mod) {
//	if (!ZisEven(Value)) {
//		Value += mod; // Number to halve is odd. Add modulus and then divide by 2.
//	}
//	Value /= 2;
//	return;
//}

// BPSW primality test:
// 1) If the input number is 2-SPRP composite, indicate composite and go out.
// 2) Perform a Miller-Rabin probable prime test on n 
//    If n is not a probable prime, then n is composite.
//    Otherwise, n is almost certainly prime.
// Output: 0 = probable prime.
//         1 = composite: not 2-Fermat pseudoprime.
//         2 = composite: does not pass 2-SPRP test.
//         3 = composite: does not pass Miller-Rabin test.

int PrimalityTest(const Znum &Value, long long upperBound) {
	int i, ctr;
	Znum Mult1, Mult3, Mult4;
	const Znum two = 2;

	if (Value <= 2) {
		return 0;    // Indicate prime.
	}
	if (ZisEven(Value)) {
		return 1;    // Number is even and different from 2. Indicate composite.
	}

	// Perform base 2 Strong Probable Prime test (2-SPRP)
	//see https://en.wikipedia.org/wiki/Probable_prime
	
	Mult3 = Value - 1;
	DivideBigNbrByMaxPowerOf2(ctr, Mult3);  // Mult3 /= 2^ctr
	mpz_powm(ZT(Mult1), ZT(two), ZT(Mult3), ZT(Value));  // Mult1 = 2^Mult3 (mod Value)
	
	if(Mult1 != 1 && Mult1 != Value-1) {
		for (i = 0; i < ctr; i++)
		{               // Loop that squares number.
			mpz_powm(ZT(Mult4), ZT(Mult1), ZT(two), ZT(Value));
			if (Mult4 ==1)
			{  // Current value is 1 but previous value is not 1 or -1: composite
				return 2;       // Composite. Not 2-strong probable prime.
			}
			if (Mult4 == Value-1) {
				i = -1;         // Number is strong pseudoprime.
				break;
			}
			Mult1 = Mult4;
		}
		if (i == ctr) {
			return 1;         // composite. Not 2-Fermat probable prime.
		}
		if (i != -1) {
			return 2;         // Composite. Not 2-strong probable prime.
		}
	}
#ifdef __MPIR_VERSION
	static bool first = true;
	static gmp_randstate_t rstate;
	if (first) {
		gmp_randinit_default(rstate);
		first = false;
	}

	auto rv = mpz_likely_prime_p(ZT(Value), rstate, upperBound);
#else
	auto rv = mpz_probab_prime_p(ZT(Value), 16);
#endif
	if (rv == 0)
		return 3;			// composite - fails Miller-Rabin test
	else
		return 0;			// probably prime;
}

/* returns true iff value is zero*/
bool BigNbrIsZero(const limb *value, int NumberLength) {
	int ctr;
	for (ctr = 0; ctr < NumberLength; ctr++) {
		if (value[ctr].x != 0) {
			return false;  // Number is not zero.
		}
	}
	return true;      // Number is zero
}

/* convert value in Limb to a double (floating point)
note that ptrLimb points AFTER last valid value in limbs.
up to 3 most significant limbs are used. */
double getMantissa(const limb *ptrLimb, int nbrLimbs) {
	double dN = (double)(ptrLimb - 1)->x;
	double dInvLimb = 1 / (double)LIMB_RANGE;
	if (nbrLimbs > 1) {
		dN += (double)(ptrLimb - 2)->x * dInvLimb;
	}
	if (nbrLimbs > 2) {
		dN += (double)(ptrLimb - 3)->x * dInvLimb * dInvLimb;
	}
	return dN;
}

/* convert Znum to BigInteger. Returns false if number is too big to convert.
this function is also used to overload the assignment operator */
bool ZtoBig(BigInteger &number, Znum numberZ) {
	number.nbrLimbs = 0;
	bool neg = false;
	Znum quot, remainder;

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
		//numberZ = quot;
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

/* convert integer list to Znum. */
void ValuestoZ(Znum &numberZ, const int number[], int NumberLength) {
	numberZ = 0;
	for (int i = NumberLength-1; i >= 0; i--) {
		//numberZ *= LIMB_RANGE;
		mpz_mul_2exp(ZT(numberZ), ZT(numberZ), BITS_PER_GROUP);  // shift numberZ left
		numberZ += number[i];
	}
}
