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
#include <intrin.h>
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
		unsigned long long carry = 0;
		for (ctr = 0; ctr < nbrLimbs; ctr++) {
			carry = (carry >> BITS_PER_GROUP) + 
				(unsigned long long)ptrBiggerAdd[ctr].x +
				(unsigned long long)ptrSmallerAdd[ctr].x;
			ptrSum[ctr].x = (long long)(carry & MAX_INT_NBR);
		}
		if (A1Smaller)
			nbrLimbs = Addend2.nbrLimbs;
		else
			nbrLimbs = Addend1.nbrLimbs;
		for (; ctr < nbrLimbs; ctr++) {
			carry = (carry >> BITS_PER_GROUP) + (unsigned long long)ptrBiggerAdd[ctr].x;
			ptrSum[ctr].x = (long long)(carry & MAX_INT_NBR);
		}
		if (carry >= LIMB_RANGE)
		{
			ptrSum[ctr].x = 1;
			nbrLimbs++;
		}
	}
	else {    // addends have different signs. Subtract their absolute values.
		long long borrow = 0;
		for (ctr = 0; ctr < nbrLimbs; ctr++) {
			borrow = (borrow >> BITS_PER_INT_GROUP) + ptrBiggerAdd[ctr].x 
				- ptrSmallerAdd[ctr].x;
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
	limb Prodl[MAX_LEN * 2] = { 0 };    // temporary variable 
	BigInteger Product;
	int nbrLimbsFactor1 = Factor1.nbrLimbs;
	int nbrLimbsFactor2 = Factor2.nbrLimbs;
	int nbrLimbs;

	if (Factor1== 0 || Factor2 == 0) {    // one or both factors are zero.
		Product = 0;     // Product is zero.
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
		memset(&((BigInteger&)Factor1).limbs[nbrLimbsFactor1], 0, 
			(nbrLimbsFactor2 - nbrLimbsFactor1) * sizeof(limb));
		nbrLimbs = nbrLimbsFactor2;
	}
	else {
		memset(&((BigInteger&)Factor2).limbs[nbrLimbsFactor2], 0, 
			(nbrLimbsFactor1 - nbrLimbsFactor2) * sizeof(limb));
		nbrLimbs = nbrLimbsFactor1;
	}

	multiply(&Factor1.limbs[0], &Factor2.limbs[0], Prodl, nbrLimbs, &nbrLimbs);

	/* MultBigNbr below is the underlying function used by multiply. */
	/*MultBigNbr((const long long *)Factor1.limbs, 
		(const long long *)Factor2.limbs, 
		(long long *)Prodl,
		nbrLimbs);*/
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
BigInteger BigIntMultiply(const BigInteger &Factor, const long long num) {
	BigInteger Product;    // temporary variable 
	unsigned long long plow, phigh = 0, carry = 0;
	size_t ix = 0;
	bool neg = false;
	long long f2 = num;

	if (num < 0) {
		neg = true;
		f2 = -f2;
	}
	for (ix = 0; ix < Factor.nbrLimbs; ix++) {
		addMult(plow, phigh, Factor.limbs[ix].x, f2, carry);
		Product.limbs[ix].x = plow;
		carry = phigh;
	}

	Product.nbrLimbs = Factor.nbrLimbs;

	if (carry != 0) {
		Product.nbrLimbs = Factor.nbrLimbs + 1;
		Product.limbs[ix].x = carry;
	}
	Product.sign = Factor.sign;
	if (neg)
		BigIntNegate(Product);
	return Product;
}

/* calculate m1*m2+c. return top 63 bits in phigh, bottom 63 bits in plow. 
Top bit of m1, m2 and c must not be set, otherwise overflow may occur. */
inline void addMult(unsigned long long &plow, unsigned long long &phigh, 
	const unsigned long long m1, const unsigned long long m2, 
	const unsigned long long c) {

	plow = _umul128(m1, m2, &phigh);  // p(128 bits) = m1*m2
	/* move top (carry) bit of plow to phigh */
	phigh <<= 64 - BITS_PER_GROUP;  // assume that long long is 64 bits wide
	phigh += (plow >> BITS_PER_GROUP);  
	plow &= MAX_INT_NBR;        // clear top bit of plow

	plow += c;
	/* if carry bit is set by adding c, add carry to phigh then reset it */
	phigh += (plow >> BITS_PER_GROUP);
	plow &= MAX_INT_NBR;
}

 
/* BigInt = e^log
throw exception if result would be to large for a BigInteger */
void expBigInt(BigInteger &bigInt, double log)
{
	long long mostSignificantLimb;
	double logar;
	logar = log/std::log(2);  // convert log to base 2
	/* calculate required value for most significant limb, initially in floating point */
	
	if (logar < 1023) {
		/* if exp will not cause FP overflow, use direct method */
		bigInt = std::exp(log);    /* actual conversion is done by DoubleToBigInt */
		return;
	}

	memset(bigInt.limbs, 0, MAX_LEN * sizeof(limb));

	/* method below is not very accurate, but does avoid overflow.  */
	bigInt.nbrLimbs = (int)floor(logar / BITS_PER_GROUP);  // nbrLimbs will be increased as required
	double value = round(exp((logar - BITS_PER_GROUP*bigInt.nbrLimbs) * std::log(2)));
	double value2 = round(exp((logar - BITS_PER_GROUP * (bigInt.nbrLimbs-1)) * std::log(2)));
	value2 -= value * (double)LIMB_RANGE;
	/* value2 is approximate. Make sure it's in legal range*/
	if (value2 < 0)
		value2 = 0;
	if (value2 > MAX_INT_NBR)
		value2 = (double)MAX_INT_NBR;

	bigInt.sign = SIGN_POSITIVE;
	bigInt.nbrLimbs++;
	/* convert double to BigInteger */
	mostSignificantLimb = (long long) value;
	if (mostSignificantLimb >= LIMB_RANGE)
	{
		mostSignificantLimb = 1;
		bigInt.nbrLimbs++;
		bigInt.limbs[bigInt.nbrLimbs - 1].x = 1;
		bigInt.limbs[bigInt.nbrLimbs - 2].x = 0;
		bigInt.limbs[bigInt.nbrLimbs - 3].x = (long long)value2;
	}
	else {
		bigInt.limbs[bigInt.nbrLimbs - 1].x = mostSignificantLimb;
		bigInt.limbs[bigInt.nbrLimbs - 2].x = (long long)value2;
	}

	if (bigInt.nbrLimbs > MAX_LEN) {
		std::string line = std::to_string(__LINE__);
		std::string mesg = "number too big : cannot expand BigInteger: ";
		mesg += __func__;
		mesg += " line ";  mesg += line;
		mesg += " in file "; mesg += __FILE__;
		throw std::range_error(mesg);
	}
}


/* convert double dvalue to bigInt. Conversion is only accurate to about 15 significant digits. */
void DoubleToBigInt(BigInteger &bigInt, double dvalue) {

	if (dvalue - 0.5 > LLONG_MIN && dvalue + 0.5 < LLONG_MAX) {
		long long vv = (long long)std::round(dvalue); // convert directly to long long if possible
		bigInt = vv;
		return;
	}
	Znum temp;
	mpz_set_d(ZT(temp), dvalue);  // convert double to Znum
	// this method has been tested to be at least as accurate as direct conversion,
	// and obviously it's easier to use standard library functions.
	bool ok = ZtoBig(bigInt, temp);        // convert Znum to BigInt
	if (ok)
		return;
	else {
		std::string line = std::to_string(__LINE__);
		std::string mesg = "number too big : cannot convert to BigInteger: ";
		mesg += __func__;
		mesg += " line ";  mesg += line;
		mesg += " in file "; mesg += __FILE__;
		throw std::range_error(mesg);
	}
}

/* estimate natural log of BigInt. Sign of BigInt is ignored */
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
double logLimbs(const limb *pBigNbr, int nbrLimbs)
{
	double logar;
	
	logar = log((double)pBigNbr[nbrLimbs - 1].x) +
		(double)(nbrLimbs - 1)*log((double)LIMB_RANGE);
	
	return logar;
}

/* calculate base^exponent. */
void BigIntPowerIntExp(const BigInteger &pBase, int exponent, BigInteger &Power) {
	int mask;
	if (pBase == 0) {     // Base = 0 -> power = 0
		Power = 0;
		return; // base = 0, so result is zero
	}
	Power = 1;
	for (mask = 1 << 30; mask != 0; mask >>= 1) {
		if ((exponent & mask) != 0) {
			for (; mask != 0; mask >>= 1) {
				Power *= Power;// BigIntMultiply(Power, Power, Power);
				if ((exponent & mask) != 0) {
					Power *= pBase; //BigIntMultiply(Power, Base, Power);
				}
			}
			break;
		}
	}
	return;
}
void BigIntPowerIntExp(const Znum &Base, int exponent, Znum &Power) {
	mpz_pow_ui(ZT(Power), ZT(Base), exponent);
}

/* divide by 2, use right shift for speed */
//void BigIntDivide2(BigInteger &pArg) {
//	int nbrLimbs = pArg.nbrLimbs;
//	int ctr = nbrLimbs - 1;
//	unsigned long long carry;
//	//limb *ptrLimb = &pArg->limbs[ctr];
//	limb *ptrLimb = pArg.limbs;
//	carry = 0;
//	for (; ctr >= 0; ctr--)
//	{
//		carry = (carry << BITS_PER_GROUP) + (unsigned long long)ptrLimb[ctr].x;
//		ptrLimb[ctr].x = (long long)(carry >> 1);
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
//		unsigned long long carry = 0;
//		for (ctr = 0; ctr < nbrLimbs; ctr++)
//		{
//			carry += (unsigned long long)ptrLimbs[ctr].x << 1;
//			ptrLimbs[ctr].x = (long long)(carry & MAX_VALUE_LIMB);
//			carry >>= BITS_PER_GROUP;
//		}
//		if (carry != 0)
//		{
//			ptrLimbs[ctr].x = (long long)carry;
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
static void addToAbsValue(limb *pLimbs, int *pNbrLimbs, long long addend) {
	int ctr;
	int nbrLimbs = *pNbrLimbs;
	pLimbs[0].x += addend;
	if ((unsigned long long)pLimbs->x < LIMB_RANGE)
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
static void subtFromAbsValue(limb *pLimbs, int *pNbrLimbs, long long subt) {
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
void subtractdivide(BigInteger &i, long long subt, long long divisor)
{
	int nbrLimbs = i.nbrLimbs;
	long long remainder = 0;
	
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
	i = BigIntDivideInt(i, divisor, remainder);
}


/* result += addend (used for operator overloading) */
BigInteger BigIntAddInt(BigInteger &Result, long long addend) {
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
	return Result;
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
void IntsToBigInteger(/*@in@*/const long long *ptrValues, /*@out@*/BigInteger &bigint)
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
	bigint.sign = SIGN_POSITIVE;
	if (NumberLength == 1) {
		destLimb[0].x = ptrValues[1];
		bigint.nbrLimbs = 1;
	}
	else {
		memcpy(destLimb, ptrValues+1, ptrValues[0] * sizeof(ptrValues[0]));
		bigint.nbrLimbs = (int)ptrValues[0];
		if (NumberLength > ptrValues[0])  // clear most significant limbs to zero if required
			memset(destLimb + ptrValues[0], 0, (NumberLength - ptrValues[0]) * sizeof(ptrValues[0]));
	}
}

/* uses global value NumberLength for starting value for number of limbs. */
static int getNbrLimbs(const limb *bigNbr) {
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
	return 1;  // BigNbr is zero
}

/* creates a list of values from a BigInteger, 1st entry in list is number of 
values that follow. Also uses global value NumberLength for number of ints. 
Assumes that bigint is positive. */
void BigIntegerToInts(/*@out@*/long long *ptrValues, /*@in@*/const BigInteger &bigint) {
	const limb *srcLimb = bigint.limbs;
	if (NumberLength == 1) {
		ptrValues[0] = 1;
		ptrValues[1] = srcLimb[0].x;
	}
	else {
		int nbrLimbs = getNbrLimbs(bigint.limbs); //nbrLimbs <= NumberLength
		assert(nbrLimbs == bigint.nbrLimbs);
		ptrValues[0] = nbrLimbs;
		memcpy(ptrValues + 1, srcLimb, nbrLimbs * sizeof(ptrValues[0]));
	}
}

/* convert limbs to BigInteger.  
Assumes that Values contains a +ve number */
void LimbsToBigInteger(/*@in@*/const limb Values[], 
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
		bigint.limbs[0].x = Values[0].x;
		bigint.nbrLimbs = 1;
	}
	else { 
		memcpy(bigint.limbs, Values, NumLen * sizeof(limb));

		int nbrLimbs;   // remove any leading zeros
		for (nbrLimbs = NumLen-1; nbrLimbs > 1; nbrLimbs--) {
			if (Values[nbrLimbs].x != 0) {
				break;
			}
		}
		bigint.nbrLimbs = nbrLimbs+1;
	}
	bigint.sign = SIGN_POSITIVE;
}

/* Convert BigInteger to limbs. */
void BigIntegerToLimbs(/*@out@*/limb Values[], 
	/*@in@*/const BigInteger &bigint, int NumLen)
{
	if (NumLen == 1) {
		Values[0].x = bigint.limbs[0].x;
		Values[1].x = 0;
	}
	else {
		int nbrLimbs = bigint.nbrLimbs; // use lesser of bigint.nbrLimbs & NumLen
		if (nbrLimbs >= NumLen) {
			memcpy(Values, bigint.limbs, NumLen * sizeof(limb));
			Values[NumLen].x = 0;
		}
		else {
			memcpy(Values, bigint.limbs, nbrLimbs * sizeof(limb));
			/* set any extra limbs to zero */
			memset(Values + nbrLimbs, 0, (NumLen - nbrLimbs) * sizeof(limb));
		}
	}
	if (bigint.sign == SIGN_NEGATIVE) {
		ChSignBigNbrB((long long *)Values, std::max(NumLen, bigint.nbrLimbs));
	}
}

/* number = numberZ. returns number of limbs. 
limb after most significant is set to zero */
int ZtoLimbs(limb *number, Znum numberZ, int NumLen) {
// note: numberZ is a copy of the original. Its value is changed
	bool neg = false;
	Znum quot, remainder;

	if (numberZ < 0) {
		neg = true;
		numberZ = -numberZ;  // make numberZ +ve
	}
	int i = 0;
	while (numberZ > 0) {
		if (i >= MAX_LEN || NumLen > MAX_LEN) {
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
		number[i].x = MulPrToLong(remainder);
		mpz_fdiv_q_2exp(ZT(numberZ), ZT(numberZ), BITS_PER_GROUP);
		//numberZ = quot;
		i++;
	}
	if (i < NumLen) {
		/* set any extra limbs to zero */
		memset(number + i, 0, (NumLen - i+1) * sizeof(limb));
	}
	else
		number[i].x = 0;    
	if (neg) {
		ChSignBigNbr((long long *)number, i + 1);
	}
	return i;
}
int ZtoBigNbr(long long number[], Znum numberZ) {
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
		//number[i] = MulPrToLong(remainder);
		number[i] = mpz_get_si(ZT(remainder));  // faster?? - no possibility of overflow here
		numberZ = quot;
		i++;
	}

	if (neg) {
		ChSignBigNbr(number, i);
	}
	return i;
}

/* convert limbs to Znum. assumes number is +ve */
void LimbstoZ(const limb *number, Znum &numberZ, int NumLen) {
	numberZ = 0;
	
	for (int i = NumLen - 1; i >= 0; i--) {
		mpz_mul_2exp(ZT(numberZ), ZT(numberZ), BITS_PER_GROUP);  // shift numberZ left
		numberZ += number[i].x;      // add next limb
	}
}
//void UncompressIntLimbs(/*@in@*/const long long *ptrValues, /*@out@*/limb *bigNbr, int nbrLen)
//{
//	long long nbrLimbs = *ptrValues;
//	memcpy(bigNbr, ptrValues + 1, nbrLimbs * sizeof(limb));
//	/* if nbrLen > nbrLimbs set extra limbs to 0*/
//	memset(bigNbr + nbrLimbs, 0, (nbrLen - nbrLimbs) * sizeof(limb));
//}

//void CompressIntLimbs(/*@out@*/long long *ptrValues, /*@in@*/const limb *bigint, int nbrLen)
//{
//	int nbrLimbs;
//	memcpy(ptrValues + 1, bigint, (nbrLen - 1) * sizeof(limb));
//	for (nbrLimbs = nbrLen - 1; nbrLimbs > 1; nbrLimbs--) {
//		if (ptrValues[nbrLimbs] != 0) {
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
	smaller limit on maxExpon (max about 3600 for a 20,000 digit number) */

	/* upperBound^maxExpon ≈ factor */
	unsigned long long maxExpon = (unsigned long long) (ceil(logBigNbr(factor) / log(upperBound)));
										  
	int h;
	long long modulus, Exponent;
	unsigned long long MaxPrime, j;
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

	MaxPrime = 2 * maxExpon + 3;
	if (MaxPrime > primeListMax) {
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

	for (size_t ix=5, h = primeList[ix]; (2*h+1) < MaxPrime; ix++, h = primeList[ix]) {
		int processed = 0;
		for (j = 2 * h + 1; j < MaxPrime; j += 2 * h) {
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

	/* check possible exponent values. Note that largest found exponent value
	is returned although, unless this value is prime, any divisor of this
	value is also a valid exponent. */
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

/* This routine checks whether the number factor is a perfect power. If it is 
not, it returns one. If it is a perfect power, it returns the exponent and  the 
base such that base^exponent = factor. The smallest exponent and largest base 
are returned e.g. for 2^6 it would return exponent 2 and base 8, not exponent 6 
and base 2 */
long long PowerCheckv2(const Znum &factor, Znum &Base, long long upperBound) {
	/* upperbound is the largest number already tested as a factor by trial division
	i.e. factor has no factors < upperBound. This can be used to put a much
	smaller limit on maxExpon (max about 3600 for a 20,000 digit number) */

	/* upperBound^maxExpon ≈ factor */
	unsigned long long maxExpon = (unsigned long long) (ceil(logBigNbr(factor) / log(upperBound)));
	Znum root;

	for (size_t ix = 0, p = primeList[ix]; p <= maxExpon; p = primeList[++ix]) {
		if (mpz_root(ZT(root), ZT(factor), p) != 0) {
			Base = root;
			return p;
		}
	}

	Base = factor;   // not perfect power
	return 1;
}

/* This routine checks whether the number factor is a perfect power. If it is
not, it returns one. If it is a perfect power, it returns the exponent and  the
base such that base^exponent = factor. The largest exponent and smallest base
are returned e.g. for 2^6 it would return exponent 6 and base 2, not exponent 2
and base 6 */
long long PowerCheckNew(const Znum &factor, Znum &Base, long long upperBound) {
	/* upperbound is the largest number already tested as a factor by trial division
	i.e. factor has no factors < upperBound.  */
	long long exponent = 1, tempexp;
	Znum tempbase, tempfactor= factor;

	/* each time round the loop a prime factor of the final exponent is found. 
	the result of each loop is combined so that the exponent is maximised
	and the base is minimised. This is optimised for the common case where a random 
	number is not a perfect power (but it's still slightly slower than DA's original version) */
	do {
		tempexp = PowerCheckv2(tempfactor, tempbase, upperBound);
		if (tempexp != 1) {
			exponent *= tempexp;
			tempfactor = tempbase;
		}
	} while (tempexp != 1);

	Base = tempfactor;
	return exponent;
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
bool BigNbrIsZero(const limb *value, int NumLen) {
	int ctr;
	for (ctr = 0; ctr < NumLen; ctr++) {
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
	double dInvLimb = 1.0 / (double)LIMB_RANGE;
	if (nbrLimbs > 1) {
		dN += (double)(ptrLimb - 2)->x * dInvLimb;
	}
	if (nbrLimbs > 2) {
		dN += (double)(ptrLimb - 3)->x * dInvLimb * dInvLimb;
	}
	return dN;
}


/* convert num to long long. If num > max, truncate it (right shift)
and exp is > 0 and represents the number of discarded bits. */
long long BigToLL(const BigInteger &num, int &exp) {
	if (num.nbrLimbs == 1) {
		exp = 0;
		if (num.sign == SIGN_POSITIVE)
			return num.limbs[0].x ;
		else 
			return -num.limbs[0].x;
	}
	exp = num.bitLength() - 63;                // number of bits to truncate
	int shift = exp % BITS_PER_INT_GROUP;      // number of bits to shift
	if (shift == 0) {
		if (num.sign == SIGN_POSITIVE)
			return num.limbs[num.nbrLimbs - 1].x;
		else
			return -num.limbs[num.nbrLimbs - 1].x;
	}
	long long result = num.limbs[num.nbrLimbs - 1].x << (BITS_PER_INT_GROUP - shift);
	result += num.limbs[num.nbrLimbs - 2].x >> (shift);
	if (num.sign == SIGN_POSITIVE)
		return result;
	else
		return -result;
}


/* inverse of BigToLL. set num = long *2^exp */
void LLToBig(BigInteger &num, long long LL, int exp) {
	if (exp == 0) {
		num = LL;
		return;
	}
	assert(exp > 0);
	num = 0;
	num.nbrLimbs = (exp + BITS_PER_GROUP - 1) / BITS_PER_GROUP;
	assert(num.nbrLimbs <= MAX_LEN);
	int shift = exp % BITS_PER_GROUP;             // get shift
	num.limbs[num.nbrLimbs - 1].x = (LL >> shift);
	if (shift > 0) {
		int s2 = BITS_PER_GROUP - shift;  // get left shift for next limb
		num.limbs[num.nbrLimbs - 2].x = (LL << s2) & MAX_INT_NBR;
 	}
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
		number.limbs[i].x = MulPrToLong(remainder);
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
		//numberZ = -numberZ;  // put back original value in numberZ
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
void ValuestoZ(Znum &numberZ, const long long number[], int NumLen) {
	numberZ = 0;
	for (int i = NumLen-1; i >= 0; i--) {
		//numberZ *= LIMB_RANGE;
		mpz_mul_2exp(ZT(numberZ), ZT(numberZ), BITS_PER_GROUP);  // shift numberZ left
		numberZ += number[i];
	}
}

BigInteger abs(const BigInteger &Num) {
	if (Num.sign == SIGN_POSITIVE)
		return Num;
	else {
		BigInteger rv = Num;
		rv.sign = SIGN_POSITIVE;
		return rv;
	}
}

/* shift first left by the number of bits specified in shiftCtr. A -ve value 
in shiftCtr causes a right shift. 
Right Shifts simulate 2s complement arithmetic right shift.
Mathematically, the shift result is equivalent to result = first * 2^shiftCtr,
whether ShiftCtr is +ve or -ve. */
void shiftBI(const BigInteger &first, const int shiftCtr,  BigInteger &result)
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
		
		if ((first.nbrLimbs  + delta) >= MAX_LEN) {   
			// Shift too much to the left; would cause overflow
			std::string line = std::to_string(__LINE__);
			std::string mesg = "cannot shift left: result out of range ";
			mesg += __func__;
			mesg += " line ";  mesg += line;
			mesg += " in file "; mesg += __FILE__;
			throw std::range_error(mesg);
		}


		result.nbrLimbs = first.nbrLimbs+delta;
		result.sign = first.sign;
		prevLimb = 0;
		ptrSrc = nbrLimbs - 1;
		ptrDest = nbrLimbs + delta;

		for (ctr = nbrLimbs-1; ctr >= 0; ctr--)
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
		if ( shiftCtr > first.nbrLimbs * BITS_PER_GROUP)
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
			result++;     //addbigint(result, 1);
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
			result--;        // addbigint(result, -1);
		}
	}
	return;
}
