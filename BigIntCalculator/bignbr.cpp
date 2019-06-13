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
#include <string>
#include <cmath>
#include "bignbr.h"
#include "int128.h"

#define MAX_LIMBS_SIQS 15

/* nbr = - nbr */
void ChSignBigNbr(long long nbr[], int length)
{
	long long carry = 0;
	int ctr;
	for (ctr = 0; ctr < length; ctr++)
	{
		carry -= nbr[ctr];
		nbr[ctr] = carry & MAX_INT_NBR;
		carry >>= BITS_PER_INT_GROUP;
	}
}

/* nbr = - nbr */
void ChSignBigNbrB(long long nbr[], int length)
{
	long long carry = 0;
	int ctr;
	for (ctr = 0; ctr < length - 1; ctr++)
	{
		carry -= nbr[ctr];
		nbr[ctr] = carry & MAX_INT_NBR;
		carry >>= BITS_PER_INT_GROUP;
	}
	nbr[ctr] = carry - nbr[ctr];  /* last word does not have most significant bit masked off */
}
int BigNbrLen(const long long Nbr[], int nbrLen) {
	int ix;
	for (ix = nbrLen; ix > 0; ix--) {
		if (Nbr[ix - 1] != 0)
			break;
	}
	return ix;
}
/* Sum = Nbr1+Nbr2  */
void AddBigNbr(const long long Nbr1[], const long long Nbr2[], long long Sum[], int nbrLen)
{
	unsigned long long carry = 0;
	int i;
	/* Sum = Nbr1 + Nbr2 */
	for (i = 0; i < nbrLen; i++)
	{
		carry = (carry >> BITS_PER_INT_GROUP) + (unsigned long long)Nbr1[i] + 
			(unsigned long long)Nbr2[i];
		Sum[i] = (long long)(carry & MAX_INT_NBR);
	}

	carry = carry >> BITS_PER_INT_GROUP;
	if (carry != 0)
		Sum[i] = carry;

}
void AddBigNbr(const Znum &Nbr1, const Znum &Nbr2, Znum &Sum) {
	//Sum = Nbr1 + Nbr2;
	mpz_add(ZT(Sum), ZT(Nbr1), ZT(Nbr2));
}

/* Diff = Nbr1-Nbr2 */
void SubtractBigNbr(const long long Nbr1[], const long long Nbr2[], long long Diff[], int nbrLen)
{
	long long borrow = 0;
	int i;
	for (i = 0; i < nbrLen; i++)
	{
		borrow = (borrow >> BITS_PER_INT_GROUP) + Nbr1[i] - Nbr2[i];
		Diff[i] = borrow & MAX_INT_NBR;
	}
}
void SubtractBigNbr(const Znum &Nbr1, const Znum &Nbr2, Znum &Diff) {
	//Diff = Nbr1 - Nbr2;
	mpz_sub(ZT(Diff), ZT(Nbr1), ZT(Nbr2));
}

/* Sum = Nbr1+Nbr2 - allows for -ve numbers  */
void AddBigNbrB(const long long Nbr1[], const long long Nbr2[], long long Sum[], int nbrLen)
{
	unsigned long long carry = 0;
	int i;
	for (i = 0; i < nbrLen - 1; i++)
	{
		carry = (carry >> BITS_PER_INT_GROUP) + (unsigned long long)Nbr1[i] + (unsigned long long)Nbr2[i];
		Sum[i] = (long long)(carry & MAX_INT_NBR);
	}
	carry = (carry >> BITS_PER_INT_GROUP) + (unsigned long long)Nbr1[i] + (unsigned long long)Nbr2[i];
	Sum[i] = (long long)carry;   // last 'digit' is not masked with MAX_INT_NBR
}
void AddBigNbrB(const Znum &Nbr1, const Znum &Nbr2, Znum &Sum) {
	//Sum = Nbr1 + Nbr2;
	mpz_add(ZT(Sum), ZT(Nbr1), ZT(Nbr2));
}

/* Diff = Nbr1-Nbr2 - allows for -ve numbers */
void SubtractBigNbrB(const long long Nbr1[], const long long Nbr2[], long long Diff[], int nbrLen)
{
	long long borrow = 0;
	int i;
	for (i = 0; i < nbrLen - 1; i++)
	{
		borrow = (borrow >> BITS_PER_INT_GROUP) + Nbr1[i] - Nbr2[i];
		Diff[i] = borrow & MAX_INT_NBR;
	}
	borrow = (borrow >> BITS_PER_INT_GROUP) + Nbr1[i] - Nbr2[i];
	Diff[i] = borrow;  /* B version differs from SubtractBigNbr only in that the most significant
					   word does not have the most significant bit masked off*/
}
void SubtractBigNbrB(const Znum &Nbr1, const Znum &Nbr2, Znum &Diff) {
	//Diff = Nbr1 - Nbr2;
	mpz_sub(ZT(Diff), ZT(Nbr1), ZT(Nbr2));
}

/* Sum = Nbr1+Nbr2 (mod Mod) */
//void AddBigNbrModN(const long long Nbr1[], const long long Nbr2[], 
//  long long Sum[], const long long Mod[], int nbrLen){
//	long long borrow = 0;
//	unsigned long long carry = 0;
//	int i;
//	for (i = 0; i < nbrLen; i++)
//	{
//		carry = (carry >> BITS_PER_INT_GROUP) + 
//			(unsigned long long)Nbr1[i] + (unsigned long long)Nbr2[i];
//		Sum[i] = (int)(carry & MAX_INT_NBR);
//	}
//	carry >>= BITS_PER_INT_GROUP;
//	//pSum -= nbrLen;
//	for (i = 0; i < nbrLen; i++)
//	{
//		borrow = (borrow >> BITS_PER_INT_GROUP) + Sum[i] - Mod[i];
//		Sum[i] = borrow & MAX_INT_NBR;
//	}
//	borrow >>= BITS_PER_INT_GROUP;
//	if (borrow + (long long)carry != 0)
//	{    // Sum is less than zero. Add Mod again.
//		//pSum -= nbrLen;
//		//pMod -= nbrLen;
//		carry = 0;
//		for (i = 0; i < nbrLen; i++)
//		{
//			carry = (carry >> BITS_PER_INT_GROUP) + 
//				(unsigned long long)Sum[i] + (unsigned long long)Mod[i];
//			Sum[i] = (long long)(carry & MAX_INT_NBR);
//		}
//	}
//}
void AddBigNbrModN(const Znum &Nbr1, const Znum &Nbr2, Znum &Diff, const Znum &Mod) {
	//Diff = (Nbr1 + Nbr2) % Mod;
	mpz_add(ZT(Diff), ZT(Nbr1), ZT(Nbr2));    // Diff = Nbr1 + Nbr2
	mpz_mod(ZT(Diff), ZT(Diff), ZT(Mod));  // Diff %= Mod
}

/* Diff = Nbr1-Nbr2 (mod Mod)*/
//void SubtractBigNbrModN(const long long Nbr1[], const long long Nbr2[], 
//	long long Diff[], const long long Mod[], int nbrLen)
//{
//	long long borrow = 0;
//	int i;
//	for (i = 0; i < nbrLen; i++)
//	{
//		borrow = (borrow >> BITS_PER_INT_GROUP) + Nbr1[i] - Nbr2[i];
//		Diff[i] = borrow & MAX_INT_NBR;
//	}
//	if (borrow != 0)
//	{
//		unsigned long long carry = 0;
//		//pDiff -= nbrLen;
//		for (i = 0; i < nbrLen; i++)
//		{
//			carry = (carry >> BITS_PER_INT_GROUP) 
//				+ (unsigned long long)Diff[i] + (unsigned long long)Mod[i];
//			Diff[i] = (long long)(carry & MAX_INT_NBR);
//		}
//	}
//}
void SubtractBigNbrModN(const Znum &Nbr1, const Znum &Nbr2, Znum &Diff, const Znum &Mod) {
	//Diff = (Nbr1 - Nbr2) % Mod;
	mpz_sub(ZT(Diff), ZT(Nbr1), ZT(Nbr2));    // Diff = Nbr1 - Nbr2
	mpz_mod(ZT(Diff), ZT(Diff), ZT(Mod));  // Diff %= Mod
}

/* Prod = Factor*f2. 
On return nbrLen is incremented by 1 if length of Prod > length of Factor */
void MultBigNbrByInt(const long long Factor[], long long f2, 
	long long bigProd[], int &nbrLen) {
	unsigned long long plow, phigh = 0, carry = 0;
	size_t ix = 0;
	bool factorPositive = true;

	if (f2 < 0) {     // If factor is negative, indicate it and compute its absolute value.
		factorPositive = false;
		f2 = -f2;
	}

	for (ix = 0; ix < nbrLen; ix++) {
		addMult(plow, phigh, Factor[ix], f2, carry);
		bigProd[ix] = plow;
		carry = phigh;
	}
	bigProd[ix] = carry;

	if (carry != 0) {
		nbrLen++;
	}

	if (!factorPositive) {   // If factor is negative, change sign of product.
		ChSignBigNbr(bigProd, nbrLen);
	}
}
void MultBigNbrByInt(const Znum &bigFactor, int factor, Znum &bigProd) {
	//bigProd = bigFactor * factor;
	mpz_mul_si(ZT(bigProd), ZT(bigFactor), factor);
}

/* bigProd = bigFactor*factor */
void MultBigNbrByIntB(const long long bigFactor[], int factor, long long bigProd[], int nbrLen)
{
	long long *bigProduct = bigProd;
	long long low;
	double dFactor;
	double dVal = 1 / (double)(1ULL << BITS_PER_INT_GROUP);
	bool factorPositive = true;
	long long ctr, carry;
	if (factor < 0)
	{     // If factor is negative, indicate it and compute its absolute value.
		factorPositive = false;
		factor = -factor;
	}
	dFactor = (double)factor;
	carry = 0;
	for (ctr = 0; ctr < nbrLen - 1; ctr++) {
		low = (bigFactor[ctr] * factor + carry) & MAX_INT_NBR;
		// Subtract or add 0x20000000 so the multiplication by dVal is not nearly an integer.
		// In that case, there would be an error of +/- 1.
		if (low < HALF_INT_RANGE) {
			carry = (long long)floor(((double)bigFactor[ctr] * dFactor + (double)carry + HALF_INT_RANGE / 2)*dVal);
		}
		else {
			carry = (long long)floor(((double)bigFactor[ctr] * dFactor + (double)carry - HALF_INT_RANGE / 2)*dVal);
		}
		bigProduct[ctr] = low;
		//bigFactor++;
	}

	low = bigFactor[ctr] * factor + carry;
	if ((low & MAX_INT_NBR) < HALF_INT_RANGE) {
		carry = (long long)floor(((double)bigFactor[ctr] * dFactor + (double)carry + HALF_INT_RANGE / 2)*dVal);
	}
	else {
		carry = (long long)floor(((double)bigFactor[ctr] * dFactor + (double)carry - HALF_INT_RANGE / 2)*dVal);
	}
	bigProduct[ctr] = low;
	if (!factorPositive) {       // If factor is negative, change sign of product.
		ChSignBigNbrB(bigProd, nbrLen);
	}
}
void MultBigNbrByIntB(const Znum &bigFactor, int factor, Znum &bigProd) {
	//bigProd = bigFactor*factor;
	mpz_mul_si(ZT(bigProd), ZT(bigFactor), factor);
}

/* Quotient = Dividend/divisor. Dividend and divisor must be +ve */
void DivBigNbrByInt(const long long Dividend[], long long divisor, 
	long long Quotient[], int nbrLen, long long &rem) {
	int ctr;
	_int128 remainder = 0;
	_int128 quotient;
	_int128 dividend;

	for (ctr = nbrLen - 1; ctr >= 0; ctr--) {
		dividend = (remainder << BITS_PER_INT_GROUP) + Dividend[ctr];
		quotient = dividend / divisor;
		remainder = dividend - quotient * divisor;
		assert(HI64(quotient) == 0);  // assert that quotient does not overflow long long
		Quotient[ctr] = quotient;
	}
	
	rem = (long long)remainder;   // NB cast is essential
}
void DivBigNbrByInt(const Znum &Dividend, int divisor, Znum &Quotient) {
	Quotient = Dividend / divisor;
}

/* calculate dividend%divisor. remainder has same type as divisor 
eror occurs if divisor is zero!! */
//long long RemDivBigNbrByInt(const long long Dividend[], long long divisor, int nbrLen)
//{
//	int ctr;
//	long long remainder = 0;
//	double dDivisor = (double)divisor;
//	double dLimb = 0x80000000;
//
//	for (ctr = nbrLen - 1; ctr >= 0; ctr--) {
//		unsigned long long quotient, dividend;
//		double dQuotient, dDividend;
//		dividend = (remainder << BITS_PER_INT_GROUP) + Dividend[ctr];
//		dDividend = (double)remainder * dLimb + Dividend[ctr];
//		dQuotient = floor(dDividend / dDivisor + 0.5);
//		quotient = (unsigned int)dQuotient;   // quotient has correct value or 1 more.
//		remainder = dividend - quotient * divisor;
//		if ((unsigned long long)remainder >= (unsigned long long)divisor)
//		{     // remainder not in range 0 <= remainder < divisor. Adjust.
//			quotient--;
//			remainder += divisor;
//		}
//	}
//	return remainder;
//}
mpir_ui RemDivBigNbrByInt(const Znum &Dividend, mpir_ui divisor) {
	return mpz_fdiv_ui(ZT(Dividend), divisor);
	/* Note that using % operator with Znums returns a Znum, even though the 
	divisor is an integer. This way is much more efficient */
}

/* Prod = Fact1*Fact2 */
void MultBigNbr(const Znum &Fact1, const Znum &Fact2, Znum &Prod) {
	//Prod = Fact1*Fact2;
	mpz_mul(ZT(Prod), ZT(Fact1), ZT(Fact2));
}
long long aux[MAX_LEN+1];
void MultBigNbr(const long long Fact1[], const long long Fact2[], 
	long long Prod[], int nbrLen) {
	int prodlen = nbrLen*2;   // length of product is twice length of Fact1 & Fact2

	/* important!! Prod must not overlap with Fact1 or Fact2 !! */
	memset(Prod, 0, prodlen * sizeof(long long));  // clear out any previous stuff 

	for (int ix = 0; ix < nbrLen; ix++) {
		auto nl2 = nbrLen;
		if (Fact2[ix] != 0) {
			memset(aux, 0, sizeof(aux));
			/* note: on return from MultBigNbrByInt, adjusted length of product is in nl2.
			nl2 = nbrLen or nbrLen+1 */
			MultBigNbrByInt(Fact1, Fact2[ix], aux, nl2);     // get partial product in aux
			AddBigNbr(aux, Prod + ix, Prod + ix, nbrLen+1);  // add aux to prod
		}
	}

#ifdef _DEBUG
	/* check that result is correct */
	//Znum F1, F2, prod, prod2, diff;
	//size_t F1size, F2size, prodsize, prod2size;
	//LimbstoZ((limb *)Fact1, F1, nbrLen);
	//LimbstoZ((limb *)Fact2, F2, nbrLen);
	//prod = F1 * F2;
	//LimbstoZ((limb *)Prod, prod2, nbrLen * 2);
	//diff = prod ^ prod2;        // compare using XOR
	//if (prod != prod2) {
	//	printf_s("Check Multiplication Failure nbrlen = %d", nbrLen);
	//	gmp_printf(" diff = %#Zx \nF1 = %#Zx \nF2 = %#Zx \n", diff, F1, F2);
	//	gmp_printf("prod2 = %#Zx \n", prod2);
	//	gmp_printf("prod  = %#Zx \n", prod);
	//	F1size = mpz_sizeinbase(ZT(F1), 2);
	//	F2size = mpz_sizeinbase(ZT(F2), 2);
	//	prodsize = mpz_sizeinbase(ZT(prod), 2);
	//	prod2size = mpz_sizeinbase(ZT(prod2), 2);
	//	printf_s("F1size = %zd F2size =%zd, prodsize =%zd, prod2size = %zd \n",
	//		F1size, F2size, prodsize, prod2size);
	//}
#endif

}

/* bigNbr = value */
void IntToBigNbr(long long value, long long bigNbr[], int nbrLength)
{
	if (value >= 0)
	{     // value is positive.
		bigNbr[0] = value;
		value = 0;
	}
	else
	{     // value is negative.
		bigNbr[0] = value & MAX_INT_NBR;
		value = MAX_INT_NBR;
	}
	for (int i = 1; i < nbrLength; i++)
	{
		bigNbr[i] = value;
	}
}

/* BigNbr = BigInt. Return value is number of limbs */
int BigIntToBigNbr(const BigInteger &BigInt, long long BigNbr[])
{
	/*int nbrLenBigNbr = BigInt.nbrLimbs;
	const limb *Limbs = BigInt.limbs;*/
	memcpy(BigNbr, BigInt.limbs, BigInt.nbrLimbs * sizeof(limb));
	return BigInt.nbrLimbs;
}

/* BigInt = BigNum */
void BigNbrToBigInt(BigInteger &pBigInt, const long long BigNum[], int nbrLenBigNum)
{
	int nbrLimbs;

	limb *Limbs = pBigInt.limbs;
	pBigInt.sign = SIGN_POSITIVE;
	memcpy(Limbs, BigNum, nbrLenBigNum * sizeof(limb));

	nbrLimbs = nbrLenBigNum;
	do 	{
		if (Limbs[nbrLimbs-1].x != 0) {
			break;
		}
	} while (--nbrLimbs > 1);
	pBigInt.nbrLimbs = nbrLimbs;
}

//void GcdBigNbr(const long long *pNbr1, const long long *pNbr2, long long *pGcd, int nbrLen)
//{
//	static BigInteger BigInt1, BigInt2, BigGcd;   // follow recommendation not to use stack
//
//	BigNbrToBigInt(BigInt1, pNbr1, nbrLen);
//	BigNbrToBigInt(BigInt2, pNbr2, nbrLen);
//	BigIntGcd(BigInt1, BigInt2, BigGcd);
//	memset(pGcd, 0, NumberLength * sizeof(limb));
//	BigIntToBigNbr(BigGcd, pGcd);
//}
//
//// Compute Nbr <- Nbr mod Mod.
////static void AdjustBigNbrModN(long long *Nbr, const long long *Mod, int nbrLen) {
////	/* note use of cast to change ints to limbs! */
////	AdjustModN((limb *)Nbr, (limb *)Mod, nbrLen);  
////}

/* Prod = Nbr1*Nbr2 (mod Mod) */
//void MultBigNbrModN(const long long Nbr1[], const long long Nbr2[], long long Prod[], const long long Mod[], int nbrLen) {
//	int i;
//	long long arr[MAX_LIMBS_SIQS];
//
//	if (nbrLen >= 2 && Mod[nbrLen - 1] == 0) {
//		nbrLen--;
//	}
//	(long long)Nbr2[nbrLen] = 0;   // // 'value' of Nbr2 is not changed
//	memset(Prod, 0, nbrLen * sizeof(Prod[0]));
//	i = nbrLen;
//	do {
//		long long Nbr = Nbr1[--i];
//		int j = nbrLen;
//		do {
//			Prod[j] = Prod[j - 1];
//		} while (--j > 0);
//		Prod[0] = 0;
//		AdjustBigNbrModN(Prod, Mod, nbrLen);
//		MultBigNbrByIntModN(Nbr2, Nbr, arr, Mod, nbrLen);
//		AddBigNbrModN(arr, Prod, Prod, Mod, nbrLen);
//	} while (i > 0);
//}
void MultBigNbrModN(const Znum &Nbr1, const Znum &Nbr2, Znum &Prod, const Znum &Mod) {
	//Prod = Nbr1*Nbr2;
	mpz_mul(ZT(Prod), ZT(Nbr1), ZT(Nbr2));
	mpz_mod(ZT(Prod), ZT(Prod), ZT(Mod));
}

/* Prod = Nbr1*Nbr2 (mod Mod) */
//void MultBigNbrByIntModN(const long long Nbr1[], long long Nbr2, 
//	long long Prod[], const long long Mod[], int nbrLen)
//{
//	if (nbrLen >= 2 && *(Mod + nbrLen - 1) == 0) {
//		nbrLen--;
//	}
//	(long long)Nbr1[nbrLen] = 0;  // ´value' of Nbr1 is not changed
//	modmultIntExtended((limb *)Nbr1, Nbr2, (limb *)Prod, (limb *)Mod, nbrLen);
//}
void MultBigNbrByIntModN(const Znum &Nbr1, int Nbr2, Znum &Prod, const Znum &Mod) {
	//Prod = Nbr1*Nbr2;
	mpz_mul_si(ZT(Prod), ZT(Nbr1), Nbr2);
	mpz_mod(ZT(Prod), ZT(Prod), ZT(Mod));
}

/* calculate NbrMod^Expon%currentPrime.
calculation uses doubles to avoid risk of overflow. */
int intDoubleModPow (int NbrMod, int Expon, int currentPrime) {
	double Power = 1;
	double Square = NbrMod;
	double Modulus = currentPrime;
	while (Expon != 0) {
		if ((Expon & 1) == 1) {
			Power *= Square;
			Power -= floor(Power / Modulus)*Modulus;
		}
		Square *= Square;
		Square -= floor(Square / Modulus)*Modulus;
		Expon >>= 1;
	}
	return (int)Power;
}

/* get modular inverse of num wrt mod*/
//void ModInvBigNbr(int *num, int *inv, int *mod, int nbrLenBigInt)
//{
//	int NumberLengthBigInt;
//	int NumberLengthBak = NumberLength;
//	static BigInteger Numerator, Denominator, Modulus, Quotient;  // follow recommendation not to use stack
//	memset(inv, 0, nbrLenBigInt * sizeof(int));
//	while (nbrLenBigInt > 1) {
//		if (*(mod + nbrLenBigInt - 1) != 0) {
//			break;
//		}
//		nbrLenBigInt--;
//	}
//	NumberLength = nbrLenBigInt;
//	memcpy(TestNbr, mod, NumberLength * sizeof(limb));
//	TestNbr[NumberLength].x = 0;
//	GetMontgomeryParms(NumberLength);
//	BigNbrToBigInt(Denominator, num, nbrLenBigInt);
//	BigNbrToBigInt(Modulus, mod, nbrLenBigInt);
//	Numerator.sign = SIGN_POSITIVE;
//	Numerator.nbrLimbs = 1;
//	Numerator.limbs[0].x = 1;    // Numerator <- 1.
//	BigIntModularDivision(Numerator, Denominator, Modulus, Quotient);
//	NumberLengthBigInt = BigIntToBigNbr(Quotient, inv);
//	NumberLength = NumberLengthBak;
//	if (NumberLengthBigInt < NumberLength)
//	{
//		memset(inv + NumberLengthBigInt, 0, (NumberLength - NumberLengthBigInt) * sizeof(int));
//	}
//	return;
//}
void ModInvBigNbr(Znum &num, Znum &inv, Znum &mod) {
	auto rv = mpz_invert(ZT(inv), ZT(num), ZT(mod));
	assert(rv != 0);
}