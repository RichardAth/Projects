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


/* nbr = - nbr */
void ChSignBigNbr(int nbr[], int length)
{
	int carry = 0;
	int ctr;
	for (ctr = 0; ctr < length; ctr++)
	{
		carry -= nbr[ctr];
		nbr[ctr] = carry & MAX_INT_NBR;
		carry >>= BITS_PER_INT_GROUP;
	}
}

/* nbr = - nbr */
void ChSignBigNbrB(int nbr[], int length)
{
	int carry = 0;
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
//void AddBigNbr(const int Nbr1[], const int Nbr2[], int Sum[], int nbrLen)
//{
//	unsigned int carry = 0;
//	int i;
//	for (i = 0; i < nbrLen; i++)
//	{
//		carry = (carry >> BITS_PER_INT_GROUP) + (unsigned int)Nbr1[i] + (unsigned int)Nbr2[i];
//		Sum[i] = (int)(carry & MAX_INT_NBR);
//	}
//}
//void AddBigNbr(const Znum &Nbr1, const Znum &Nbr2, Znum &Sum) {
//	//Sum = Nbr1 + Nbr2;
//	mpz_add(ZT(Sum), ZT(Nbr1), ZT(Nbr2));
//}

/* Diff = Nbr1-Nbr2 */
//void SubtractBigNbr(const int Nbr1[], const int Nbr2[], int Diff[], int nbrLen)
//{
//	int borrow = 0;
//	int i;
//	for (i = 0; i < nbrLen; i++)
//	{
//		borrow = (borrow >> BITS_PER_INT_GROUP) + Nbr1[i] - Nbr2[i];
//		Diff[i] = borrow & MAX_INT_NBR;
//	}
//}
//void SubtractBigNbr(const Znum &Nbr1, const Znum &Nbr2, Znum &Diff) {
//	//Diff = Nbr1 - Nbr2;
//	mpz_sub(ZT(Diff), ZT(Nbr1), ZT(Nbr2));
//}

/* Sum = Nbr1+Nbr2  */
//void AddBigNbrB(const int Nbr1[], const int Nbr2[], int Sum[], int nbrLen)
//{
//	unsigned int carry = 0;
//	int i;
//	for (i = 0; i < nbrLen - 1; i++)
//	{
//		carry = (carry >> BITS_PER_INT_GROUP) + (unsigned int)Nbr1[i] + (unsigned int)Nbr2[i];
//		Sum[i] = (int)(carry & MAX_INT_NBR);
//	}
//	carry = (carry >> BITS_PER_INT_GROUP) + (unsigned int)Nbr1[i] + (unsigned int)Nbr2[i];
//	Sum[i] = (int)carry;   // last 'digit' is not masked with MAX_INT_NBR
//}
//void AddBigNbrB(const Znum &Nbr1, const Znum &Nbr2, Znum &Sum) {
//	//Sum = Nbr1 + Nbr2;
//	mpz_add(ZT(Sum), ZT(Nbr1), ZT(Nbr2));
//}

/* Diff = Nbr1-Nbr2 */
//void SubtractBigNbrB(const int Nbr1[], const int Nbr2[], int Diff[], int nbrLen)
//{
//	int borrow = 0;
//	int i;
//	for (i = 0; i < nbrLen - 1; i++)
//	{
//		borrow = (borrow >> BITS_PER_INT_GROUP) + Nbr1[i] - Nbr2[i];
//		Diff[i] = borrow & MAX_INT_NBR;
//	}
//	borrow = (borrow >> BITS_PER_INT_GROUP) + Nbr1[i] - Nbr2[i];
//	Diff[i] = borrow;  /* B version differs from SubtractBigNbr only in that the most significant
//					   word does not have the most significant bit masked off*/
//}
//void SubtractBigNbrB(const Znum &Nbr1, const Znum &Nbr2, Znum &Diff) {
//	//Diff = Nbr1 - Nbr2;
//	mpz_sub(ZT(Diff), ZT(Nbr1), ZT(Nbr2));
//}

/* Sum = Nbr1+Nbr2 (mod Mod) */
//void AddBigNbrModN(const int Nbr1[], const int Nbr2[], int Sum[], const int Mod[], int nbrLen)
//{
//	int borrow = 0;
//	unsigned int carry = 0;
//	int i;
//	for (i = 0; i < nbrLen; i++)
//	{
//		carry = (carry >> BITS_PER_INT_GROUP) + 
//			(unsigned int)Nbr1[i] + (unsigned int)Nbr2[i];
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
//	if (borrow + (int)carry != 0)
//	{    // Sum is less than zero. Add Mod again.
//		//pSum -= nbrLen;
//		//pMod -= nbrLen;
//		carry = 0;
//		for (i = 0; i < nbrLen; i++)
//		{
//			carry = (carry >> BITS_PER_INT_GROUP) + (unsigned int)Sum[i] + (unsigned int)Mod[i];
//			Sum[i] = (int)(carry & MAX_INT_NBR);
//		}
//	}
//}
void AddBigNbrModN(const Znum &Nbr1, const Znum &Nbr2, Znum &Diff, const Znum &Mod) {
	//Diff = (Nbr1 + Nbr2) % Mod;
	mpz_add(ZT(Diff), ZT(Nbr1), ZT(Nbr2));    // Diff = Nbr1 + Nbr2
	while (Diff < 0)
		Diff += Mod;
	mpz_mod(ZT(Diff), ZT(Diff), ZT(Mod));  // Diff %= Mod
}

/* Diff = Nbr1-Nbr2 (mod Mod)*/
//void SubtractBigNbrModN(const int Nbr1[], const int Nbr2[], int Diff[], const int Mod[], int nbrLen)
//{
//	int borrow = 0;
//	int i;
//	for (i = 0; i < nbrLen; i++)
//	{
//		borrow = (borrow >> BITS_PER_INT_GROUP) + Nbr1[i] - Nbr2[i];
//		Diff[i] = borrow & MAX_INT_NBR;
//	}
//	if (borrow != 0)
//	{
//		unsigned int carry = 0;
//		//pDiff -= nbrLen;
//		for (i = 0; i < nbrLen; i++)
//		{
//			carry = (carry >> BITS_PER_INT_GROUP) + (unsigned int)Diff[i] + (unsigned int)Mod[i];
//			Diff[i] = (int)(carry & MAX_INT_NBR);
//		}
//	}
//}
void SubtractBigNbrModN(const Znum &Nbr1, const Znum &Nbr2, Znum &Diff, const Znum &Mod) {
	//Diff = (Nbr1 - Nbr2) % Mod;
	mpz_sub(ZT(Diff), ZT(Nbr1), ZT(Nbr2));    // Diff = Nbr1 - Nbr2
	while (Diff < 0)
		Diff += Mod;
	mpz_mod(ZT(Diff), ZT(Diff), ZT(Mod));  // Diff %= Mod
}

/* bigProd = bigFactor*factor */
//void MultBigNbrByInt(const int bigFactor[], int factor, int bigProd[], int nbrLen)
//{
//	int *bigProduct = bigProd;
//	double dFactor;
//	double dVal = 1 / (double)(1U << BITS_PER_INT_GROUP);
//	bool factorPositive = true;
//	int ctr, carry;
//	if (factor < 0)
//	{     // If factor is negative, indicate it and compute its absolute value.
//		factorPositive = false;
//		factor = -factor;
//	}
//	dFactor = (double)factor;
//	carry = 0;
//	for (ctr = 0; ctr < nbrLen; ctr++)
//	{
//		int low = (bigFactor[ctr] * factor + carry) & MAX_INT_NBR;
//		// Subtract or add 0x20000000 so the multiplication by dVal is not nearly an integer.
//		// In that case, there would be an error of +/- 1.
//		if (low < HALF_INT_RANGE)
//		{
//			carry = (int)floor(((double)bigFactor[ctr] * dFactor + (double)carry + HALF_INT_RANGE / 2)*dVal);
//		}
//		else
//		{
//			carry = (int)floor(((double)bigFactor[ctr] * dFactor + (double)carry - HALF_INT_RANGE / 2)*dVal);
//		}
//		bigProduct[ctr] = low;
//		//bigFactor++;
//	}
//	if (!factorPositive)
//	{         // If factor is negative, change sign of product.
//		ChSignBigNbr(bigProd, nbrLen);
//	}
//}
//void MultBigNbrByInt(const Znum &bigFactor, int factor, Znum &bigProd) {
//	//bigProd = bigFactor * factor;
//	mpz_mul_si(ZT(bigProd), ZT(bigFactor), factor);
//}

/* bigProd = bigFactor*factor */
//void MultBigNbrByIntB(const int bigFactor[], int factor, int bigProd[], int nbrLen)
//{
//	int *bigProduct = bigProd;
//	int low;
//	double dFactor;
//	double dVal = 1 / (double)(1U << BITS_PER_INT_GROUP);
//	bool factorPositive = true;
//	int ctr, carry;
//	if (factor < 0)
//	{     // If factor is negative, indicate it and compute its absolute value.
//		factorPositive = false;
//		factor = -factor;
//	}
//	dFactor = (double)factor;
//	carry = 0;
//	for (ctr = 0; ctr < nbrLen - 1; ctr++)
//	{
//		low = (bigFactor[ctr] * factor + carry) & MAX_INT_NBR;
//		// Subtract or add 0x20000000 so the multiplication by dVal is not nearly an integer.
//		// In that case, there would be an error of +/- 1.
//		if (low < HALF_INT_RANGE)
//		{
//			carry = (int)floor(((double)bigFactor[ctr] * dFactor + (double)carry + HALF_INT_RANGE / 2)*dVal);
//		}
//		else
//		{
//			carry = (int)floor(((double)bigFactor[ctr] * dFactor + (double)carry - HALF_INT_RANGE / 2)*dVal);
//		}
//		bigProduct[ctr] = low;
//		//bigFactor++;
//	}
//	low = bigFactor[ctr] * factor + carry;
//	if ((low & MAX_INT_NBR) < HALF_INT_RANGE)
//	{
//		carry = (int)floor(((double)bigFactor[ctr] * dFactor + (double)carry + HALF_INT_RANGE / 2)*dVal);
//	}
//	else
//	{
//		carry = (int)floor(((double)bigFactor[ctr] * dFactor + (double)carry - HALF_INT_RANGE / 2)*dVal);
//	}
//	bigProduct[ctr] = low;
//	if (!factorPositive)
//	{         // If factor is negative, change sign of product.
//		ChSignBigNbrB(bigProd, nbrLen);
//	}
//}
//void MultBigNbrByIntB(const Znum &bigFactor, int factor, Znum &bigProd) {
//	//bigProd = bigFactor*factor;
//	mpz_mul_si(ZT(bigProd), ZT(bigFactor), factor);
//}

/* Quotient = Dividend/divisor */
void DivBigNbrByInt(const int Dividend[], int divisor, int Quotient[], int nbrLen)
{
	int ctr;
	int remainder = 0;
	double dDivisor = (double)divisor;  // assume divisor > 0
	double dLimb = 0x80000000;
	for (ctr = nbrLen - 1; ctr >= 0; ctr--)
	{
		double dDividend, dQuotient;
		int quotient, dividend;
		dividend = (remainder << BITS_PER_INT_GROUP) + Dividend[ctr];
		dDividend = (double)remainder * dLimb + Dividend[ctr];
		dQuotient = floor(dDividend / dDivisor + 0.5);
		quotient = (unsigned int)dQuotient;   // quotient has correct value or 1 more.
		remainder = dividend - quotient * divisor;
		if ((unsigned int)remainder >= (unsigned int)divisor)
		{     // remainder not in range 0 <= remainder < divisor. Adjust.
			quotient--;
			remainder += divisor;
		}
		Quotient[ctr] = quotient;
	}
}
//void DivBigNbrByInt(const Znum &Dividend, int divisor, Znum &Quotient) {
//	Quotient = Dividend / divisor;
//}

/* calculate dividend%divisor. remainder has same type as divisor 
eror occurs if divisor is zero!! */
//int RemDivBigNbrByInt(const int Dividend[], int divisor, int nbrLen)
//{
//	int ctr;
//	int remainder = 0;
//	double dDivisor = (double)divisor;
//	double dLimb = 0x80000000;
//	//pDividend += nbrLen - 1;
//	for (ctr = nbrLen - 1; ctr >= 0; ctr--)
//	{
//		unsigned int quotient, dividend;
//		double dQuotient, dDividend;
//		dividend = (remainder << BITS_PER_INT_GROUP) + Dividend[ctr];
//		dDividend = (double)remainder * dLimb + Dividend[ctr];
//		dQuotient = floor(dDividend / dDivisor + 0.5);
//		quotient = (unsigned int)dQuotient;   // quotient has correct value or 1 more.
//		remainder = dividend - quotient * divisor;
//		if ((unsigned int)remainder >= (unsigned int)divisor)
//		{     // remainder not in range 0 <= remainder < divisor. Adjust.
//			quotient--;
//			remainder += divisor;
//		}
//		//pDividend--;
//	}
//	return remainder;
//}
mpir_ui RemDivBigNbrByInt(const Znum &Dividend, mpir_ui divisor) {
	return mpz_fdiv_ui(ZT(Dividend), divisor);
	/* Note that using % operator with Znums returns a Znum, even though the 
	divisor is an integer. This way is much more efficient */
}

/* Prod = Fact1*Fact2 */
//void MultBigNbr(const int Fact1[], const int Fact2[], int Prod[], int nbrLen)
//{
//	double dRangeLimb = (double)(1U << BITS_PER_INT_GROUP);
//	double dInvRangeLimb = 1 / dRangeLimb;
//	int low = 0;
//	int i, j,k=0 ;
//	int factor1, factor2;
//	double dAccumulator = 0;
//	for (i = 0; i < nbrLen; i++)
//	{
//		for (j = 0; j <= i; j++)
//		{
//			factor1 = Fact1[j];
//			factor2 = Fact2[i - j];
//			low += factor1*factor2;
//			dAccumulator += (double)factor1 * (double)factor2;
//		}
//		low &= MAX_INT_NBR;    // Trim extra bits.
//		Prod[i] = low;
//		// Subtract or add 0x20000000 so the multiplication by dVal is not nearly an integer.
//		// In that case, there would be an error of +/- 1.
//		if (low < HALF_INT_RANGE)
//		{
//			dAccumulator = floor((dAccumulator + HALF_INT_RANGE / 2)*dInvRangeLimb);
//		}
//		else
//		{
//			dAccumulator = floor((dAccumulator - HALF_INT_RANGE / 2)*dInvRangeLimb);
//		}
//		low = (int)(dAccumulator - floor(dAccumulator * dInvRangeLimb) * dRangeLimb);
//	}
//	Prod[i] = low;
//	Prod[i + 1] = (int)floor(dAccumulator / dRangeLimb);
//}
//void MultBigNbr(const Znum &Fact1, const Znum &Fact2, Znum &Prod) {
//	//Prod = Fact1*Fact2;
//	mpz_mul(ZT(Prod), ZT(Fact1), ZT(Fact2));
//}

/* bigNbr = value */
//void IntToBigNbr(int value, int bigNbr[], int nbrLength)
//{
//	if (value >= 0)
//	{     // value is positive.
//		bigNbr[0] = value;
//		value = 0;
//	}
//	else
//	{     // value is negative.
//		bigNbr[0] = value & MAX_INT_NBR;
//		value = MAX_INT_NBR;
//	}
//	for (int i = 1; i < nbrLength; i++)
//	{
//		bigNbr[i] = value;
//	}
//}

/* BigNbr = BigInt. Return value is number of limbs */
//int BigIntToBigNbr(const BigInteger &BigInt, int BigNbr[])
//{
//	/*int nbrLenBigNbr = BigInt.nbrLimbs;
//	const limb *Limbs = BigInt.limbs;*/
//	memcpy(BigNbr, BigInt.limbs, BigInt.nbrLimbs * sizeof(int));
//	return BigInt.nbrLimbs;
//}

/* BigInt = BigNum */
//void BigNbrToBigInt(BigInteger &pBigInt, const int BigNum[], int nbrLenBigNum)
//{
//	int nbrLimbs;
//
//	limb *Limbs = pBigInt.limbs;
//	pBigInt.sign = SIGN_POSITIVE;
//	memcpy(Limbs, BigNum, nbrLenBigNum * sizeof(int));
//
//	nbrLimbs = nbrLenBigNum;
//	do 	{
//		if (Limbs[nbrLimbs-1].x != 0) {
//			break;
//		}
//	} while (--nbrLimbs > 1);
//	pBigInt.nbrLimbs = nbrLimbs;
//}

//void GcdBigNbr(const int *pNbr1, const int *pNbr2, int *pGcd, int nbrLen)
//{
//	static BigInteger BigInt1, BigInt2, BigGcd;   // follow recommendation not to use stack
//
//	BigNbrToBigInt(BigInt1, pNbr1, nbrLen);
//	BigNbrToBigInt(BigInt2, pNbr2, nbrLen);
//	BigIntGcd(BigInt1, BigInt2, BigGcd);
//	memset(pGcd, 0, NumberLength * sizeof(int));
//	BigIntToBigNbr(BigGcd, pGcd);
//}

// Compute Nbr <- Nbr mod Mod.
//static void AdjustBigNbrModN(int *Nbr, const int *Mod, int nbrLen) {
//	/* note use of cast to change ints to limbs! */
//	AdjustModN((limb *)Nbr, (limb *)Mod, nbrLen);  
//}

/* Prod = Nbr1*Nbr2 (mod Mod) */
//void MultBigNbrModN(const int Nbr1[], const int Nbr2[], int Prod[], const int Mod[], int nbrLen) {
//	int i;
//	int arr[MAX_LIMBS_SIQS];
//
//	if (nbrLen >= 2 && Mod[nbrLen - 1] == 0) {
//		nbrLen--;
//	}
//	(int)Nbr2[nbrLen] = 0;   // // 'value' of Nbr2 is not changed
//	memset(Prod, 0, nbrLen * sizeof(Prod[0]));
//	i = nbrLen;
//	do {
//		int Nbr = Nbr1[--i];
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
//void MultBigNbrByIntModN(const int Nbr1[], int Nbr2, int Prod[], const int Mod[], int nbrLen)
//{
//	if (nbrLen >= 2 && *(Mod + nbrLen - 1) == 0) {
//		nbrLen--;
//	}
//	(int)Nbr1[nbrLen] = 0;  // ´value' of Nbr1 is not changed
//	modmultIntExtended((limb *)Nbr1, Nbr2, (limb *)Prod, (limb *)Mod, nbrLen);
//}
void MultBigNbrByIntModN(const Znum &Nbr1, int Nbr2, Znum &Prod, const Znum &Mod) {
	//Prod = Nbr1*Nbr2;
	mpz_mul_si(ZT(Prod), ZT(Nbr1), Nbr2);
	mpz_mod(ZT(Prod), ZT(Prod), ZT(Mod));
}

/* calculate NbrMod^Expon%currentPrime.
calculation uses doubles to avoid risk of overflow. */
int modPower (int NbrMod, int Expon, int currentPrime) {
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
	double dInvLimb = 1 / (double)LIMB_RANGE;
	if (nbrLimbs > 1) {
		dN += (double)(ptrLimb - 2)->x * dInvLimb;
	}
	if (nbrLimbs > 2) {
		dN += (double)(ptrLimb - 3)->x * dInvLimb * dInvLimb;
	}
	return dN;
}

void LimbstoZ(const limb *number, Znum &numberZ, int NumLen) {
	numberZ = 0;
	for (int i = NumLen - 1; i >= 0; i--) {
		mpz_mul_2exp(ZT(numberZ), ZT(numberZ), BITS_PER_GROUP);  // shift numberZ left
		numberZ += number[i].x;      // add next limb
	}
}

int ZtoLimbs(limb *number, Znum numberZ, int NumLen) {
	// note: numberZ is a copy of the original. Its value is changed
	bool neg = false;
	Znum remainder;

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
		number[i].x = (int)MulPrToLong(remainder);
		mpz_fdiv_q_2exp(ZT(numberZ), ZT(numberZ), BITS_PER_GROUP);

		i++;
	}
	if (i < NumLen) {
		/* set any extra limbs to zero */
		memset(number + i, 0, (NumLen - i) * sizeof(limb));
	}
	if (neg) {
		ChSignBigNbr((int *)number, i + 1);
	}
	return i;
}

/* number = numberZ*/
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

/* convert integer list to Znum. */
void ValuestoZ(Znum &numberZ, const int number[], int NumLen) {
	numberZ = 0;
	for (int i = NumLen - 1; i >= 0; i--) {
		//numberZ *= LIMB_RANGE;
		mpz_mul_2exp(ZT(numberZ), ZT(numberZ), BITS_PER_GROUP);  // shift numberZ left
		numberZ += number[i];
	}
}

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

double logBigNbr(const Znum &BigInt) {
	double BigId;
#ifdef __MPIR_VERSION
	long BiExp;   // changed for MPIR version 3.0.0
#else
	long BiExp;
#endif
	BigId = mpz_get_d_2exp(&BiExp, ZT(BigInt)); // BigId * 2^BiExp = BigInt 
	double logval = log(BigId) + BiExp * log(2);
	return logval;
}

/* returns nbrMod^Expon%currentPrime.
overflow could occur if currentPrime > 2^31
the alternative is to use modPower */
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

// This routine checks whether the number factor is a perfect power. 
// If it is not, it returns exponent 1. If it is a perfect power, it returns the  
// exponent and  the base such that base^exponent = factor.
long long PowerCheck(const Znum &factor, Znum &Base, long long upperBound) {
	/* upperbound is the largest number already tested as a factor by trial division
	i.e. factor has no factors < upperBound. This can be used to put a much
	smaller limit on maxExpon (max about 2000) */

	/* upperBound^maxExpon ≈ factor */
	unsigned long long maxExpon = (unsigned long long) (ceil(logBigNbr(factor) / log(upperBound)));

	int h;
	long long modulus, Exponent;
	unsigned long long maxPrime, j;
	int prime2310x1[] =
	{ 2311, 4621, 9241, 11551, 18481, 25411, 32341, 34651, 43891, 50821 };
	// Primes of the form 2310x+1.
	bool expon2 = true, expon3 = true, expon5 = true;
	bool expon7 = true, expon11 = true;
	Znum Zmod, Root;
	std::vector<bool>ProcessExpon(maxExpon + 1);

	for (h = 0; h < sizeof(prime2310x1) / sizeof(prime2310x1[0]); h++) {
		int testprime = prime2310x1[h];
		// Zmod = mod = Bigint%testprime
		auto mod = mpz_mod_ui(ZT(Zmod), ZT(factor), testprime); // getRemainder(factor, testprime);
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

	maxPrime = 2 * maxExpon + 3;
	if (maxPrime > primeListMax) {
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

	for (size_t ix = 5, h = primeList[ix]; h < maxPrime / 2; ix++, h = primeList[ix]) {
		int processed = 0;
		for (j = 2 * h + 1; j < maxPrime; j += 2 * h) {
			if (isPrime2(j)) {
				modulus = mpz_mod_ui(ZT(Zmod), ZT(factor), j); // getRemainder(factor, j);
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