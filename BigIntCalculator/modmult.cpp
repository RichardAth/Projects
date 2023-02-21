/* calculators/modmult.c
c25525c on 10 Dec 2017
@alpertron alpertron Added general division modulo any composite
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

static limb aux3[MAX_LEN];
static limb aux4[MAX_LEN];
static limb aux5[MAX_LEN];
static limb aux6[MAX_LEN];

/* the functions below are now used only by the built-in ECM function. */

BigInteger TestNbrBI;   /* value set up by ecm */

/* values below are set up by calling GetMontgomeryParms*/
limb * const TestNbr = TestNbrBI.limbs;
BigInteger MontgomeryMultNBI;
limb * const MontgomeryMultN = MontgomeryMultNBI.limbs;   // used by  modmult, ComputeInversePower2, etc
BigInteger  MontgomeryMultR1BI;
limb * const MontgomeryMultR1 = MontgomeryMultR1BI.limbs;
BigInteger  MontgomeryMultR2BI;
limb *const MontgomeryMultR2 = MontgomeryMultR2BI.limbs;  // used by ModInvBigNbr

static int powerOf2Exponent;
static limb aux[MAX_LEN], aux2[MAX_LEN];

long long lModularMult = 0;    // count of number of modular multiplications used to control status display
mmCback modmultCallback = nullptr;     // function pointer
static limb U[MAX_LEN], V[MAX_LEN], R[MAX_LEN], S[MAX_LEN];
static limb Ubak[MAX_LEN], Vbak[MAX_LEN];
static BigInteger tmpDen, tmpNum;


// Find the inverse of value m 2^(NumberLength*BITS_PER_GROUP)
void ComputeInversePower2(const limb *value, limb *result, limb *tmp)
{
	int N, x, j;
	limb Cy;
	int currLen;
	x = N = (int)*value;       // 2 least significant bits of inverse correct.
	x = x * (2 - N * x);         // 4 least significant bits of inverse correct.
	x = x * (2 - N * x);         // 8 least significant bits of inverse correct.
	x = x * (2 - N * x);         // 16 least significant bits of inverse correct.
	x = x * (2 - N * x);         // 32 least significant bits of inverse correct.
	*result = x & MAX_VALUE_LIMB;
	for (currLen = 2; currLen <= NumberLength * 2; currLen <<= 1) {
		multiply(value, result, tmp, currLen, NULL);    // tmp <- N * x
		Cy = 2 - tmp[0];
		tmp[0] = Cy & MAX_VALUE_LIMB;
		for (j = 1; j < currLen; j++) {
			Cy = (Cy >> BITS_PER_GROUP) - tmp[j];
			tmp[j] = Cy & MAX_VALUE_LIMB;
		}                                               // tmp <- 2 - N * x
		multiply(result, tmp, result, currLen, NULL);   // tmp <- x * (2 - N * x)
	}
}

// Compute Nbr <- Nbr mod Modulus.
// Modulus has NumberLength limbs.
static void AdjustModN(limb *Nbr, const limb *Modulus, int nbrLen)
{
	int i, carry;
	int TrialQuotient;
	double dNbr, dModulus, dTrialQuotient;
	double dAccumulator, dDelta;
	double dVal = 1 / (double)LIMB_RANGE;
	double dSquareLimb = (double)LIMB_RANGE * (double)LIMB_RANGE;

	dModulus = getMantissa(Modulus + nbrLen, nbrLen);
	dNbr = getMantissa(Nbr + nbrLen + 1, nbrLen + 1) * LIMB_RANGE;
	TrialQuotient = (int)(unsigned int)floor(dNbr / dModulus + 0.5);
	if ((unsigned int)TrialQuotient >= LIMB_RANGE)
	{   // Maximum value for limb.
		TrialQuotient = MAX_VALUE_LIMB;
	}
	// Compute Nbr <- Nbr - TrialQuotient * Modulus
	dTrialQuotient = (double)TrialQuotient;
	carry = 0;
	dAccumulator = 0;
	dDelta = 0;
	for (i = 0; i <= nbrLen; i++)
	{
		int low = (Nbr[i] - Modulus[i] * TrialQuotient + carry) & MAX_INT_NBR;
		// Subtract or add 0x20000000 so the multiplication by dVal is not nearly an integer.
		// In that case, there would be an error of +/- 1.
		dAccumulator = Nbr[i] - Modulus[i] * dTrialQuotient + carry + dDelta;
		dDelta = 0;
		if (dAccumulator < 0) {
			dAccumulator += dSquareLimb;
			dDelta = -(double)LIMB_RANGE;
		}
		if (low < HALF_INT_RANGE) {
			carry = (int)floor((dAccumulator + HALF_INT_RANGE / 2)*dVal);
		}
		else {
			carry = (int)floor((dAccumulator - HALF_INT_RANGE / 2)*dVal);
		}
		Nbr[i] = low;
	}
	Nbr[i] = carry & MAX_INT_NBR;
	if ((Nbr[nbrLen] & MAX_VALUE_LIMB) != 0) {
		unsigned int cy = 0;
		for (i = 0; i < nbrLen; i++) {
			cy += (unsigned int)Nbr[i] + (unsigned int)Modulus[i];
			Nbr[i] = (int)(cy & MAX_VALUE_LIMB);
			cy >>= BITS_PER_GROUP;
		}
		Nbr[nbrLen] = 0;
	}
}


void GetMontgomeryParmsPowerOf2(int powerOf2)
{
	NumberLength = (powerOf2 + BITS_PER_GROUP - 1) / BITS_PER_GROUP;
	int NumberLengthBytes = NumberLength * (int)sizeof(limb);
	powerOf2Exponent = powerOf2;
	(void)memset(MontgomeryMultR1, 0, NumberLengthBytes);
	(void)memset(MontgomeryMultR2, 0, NumberLengthBytes);
	MontgomeryMultR1[0] = 1;
	MontgomeryMultR2[0] = 1;
}

// Let R be a power of 2 of at least len limbs.
// Compute R1 = MontgomeryMultR1 and N = MontgomeryN using the formulas:
// R1 = R mod M 
// N = -M^(-1) mod R
// also compute R2 = R^2 mod  M.
//   & powerOf2Exponent
/* uses global variables TestNbr, NumberLength */
void GetMontgomeryParms(int len) {
	int j;

	TestNbr[len] = 0;    // add a leading zero

	if (len == 1) {
		/* in reality, there can't be just 1 limb because numbers that fit into
		1 limb will already be factored before we get to this point, therefore
		code below is never executed. */
		MontgomeryMultR1BI = 1;
		MontgomeryMultR2BI = 1;
		return;
	}

	/* Check whether TestNbr is a power of 2. 
	In reality, it can't be, because TestNbr is the number to be factored and any
	small factors will already have been removed. */
	powerOf2Exponent = 0;    // Indicate not power of 2 in advance.
	for (j = 0; j < NumberLength - 1; j++) {
		if (TestNbr[j] != 0) {
			break;   /* j is index of least significant non-zero limb. the 
			least significant limb can't be zero because TstNbr must be odd, 
			therfore j is zero */
		}
	}
	if (j == len - 1) {  // is there only 1 non-zero limb?
		/* in reality, there can't be just 1 because numbers that fit into 
		1 limb will already be factored before we get to this point, therefore
		code below is never executed. */
		int value = TestNbr[NumberLength - 1];
		for (j = 0; j < BITS_PER_GROUP; j++) {
			if (value == 1) {
				powerOf2Exponent = (len - 1)*BITS_PER_GROUP + j;
				MontgomeryMultR1BI = 1;
				MontgomeryMultR2BI = 1;
				return;
			}
			value >>= 1;
		}
	}

	// Compute MontgomeryMultN as -1/TestNbr (mod 2^k) using Newton method,
	// which doubles the precision for each iteration.
	// In the formula above: k = BITS_PER_GROUP * NumberLength.
	if (len >= 8) {
		limb *ptrResult;
		limb Carry;
		ComputeInversePower2(TestNbr, MontgomeryMultN, aux);
		ptrResult = &MontgomeryMultN[0];
		Carry = 0;          // Change sign.
		for (j = 0; j < NumberLength; j++) {
			Carry = (Carry >> BITS_PER_GROUP) - *ptrResult;
			*ptrResult = Carry & MAX_VALUE_LIMB;
			ptrResult++;
		}
		*ptrResult = 0;
	}
	else { // NumberLength < 8
		int x, N;
		x = N = (int)TestNbr[0];   // 2 least significant bits of inverse correct.
		x = x * (2 - N * x);         // 4 least significant bits of inverse correct.
		x = x * (2 - N * x);         // 8 least significant bits of inverse correct.
		x = x * (2 - N * x);         // 16 least significant bits of inverse correct.
		x = x * (2 - N * x);         // 32 least significant bits of inverse correct.
		MontgomeryMultNBI = (long long)(-x) & MAX_VALUE_LIMB;    // Change sign
	}
	// Compute MontgomeryMultR1 as 1 in Montgomery notation,
	// this is 2^(NumberLength*BITS_PER_GROUP) % TestNbr.
	j = len;
	MontgomeryMultR1[j] = 1;
	do {
		MontgomeryMultR1[--j] = 0;
	} while (j > 0);
	MontgomeryMultR1BI.nbrLimbs = len + 1;

	MontgomeryMultR1BI %= TestNbrBI;
	MontgomeryMultR1[MontgomeryMultR1BI.nbrLimbs] = 0;

	/* adjust limb length of R1 if necessary */
	for (int NumberLengthR1 = NumberLength; NumberLengthR1 > 0; NumberLengthR1--) {
		if (MontgomeryMultR1[NumberLengthR1 - 1] != 0) {
			MontgomeryMultR1BI.nbrLimbs = NumberLengthR1;
			break;
		}
	}

	// Compute MontgomeryMultR2 as 2^(2*NumberLength*BITS_PER_GROUP) % TestNbr.
	MontgomeryMultR2BI = MontgomeryMultR1BI * MontgomeryMultR1BI;
	MontgomeryMultR2BI %= TestNbrBI;
	MontgomeryMultR2[MontgomeryMultR2BI.nbrLimbs] = 0;
	/* adjust limb length of R2 if necessary */
	for (int NumLenR2 = MontgomeryMultR2BI.nbrLimbs; NumLenR2 > 0; NumLenR2--) {
		if (MontgomeryMultR2[NumLenR2 - 1] != 0) {
			MontgomeryMultR2BI.nbrLimbs = NumLenR2;
			break;
		}
	}
}

void AddBigNbrModN(const limb* Nbr1, const limb* Nbr2, limb* Sum,
	const limb* mod, int nbrLen)
{
	unsigned int carry;
	unsigned int borrow;
	int i;

	carry = 0U;
	for (i = 0; i < nbrLen; i++)
	{
		carry = (carry >> BITS_PER_GROUP) +
			(unsigned int)*(Nbr1 + i) + (unsigned int)*(Nbr2 + i);
		Sum[i] = (int)(carry & MAX_VALUE_LIMB);
	}
	borrow = 0U;
	for (i = 0; i < nbrLen; i++)
	{
		borrow = (unsigned int)Sum[i] - (unsigned int)*(mod + i) - (borrow >> BITS_PER_GROUP);
		Sum[i] = (int)(borrow & MAX_VALUE_LIMB);
	}

	if ((carry < LIMB_RANGE) && ((int)borrow < 0))
	{
		carry = 0U;
		for (i = 0; i < nbrLen; i++)
		{
			carry = (carry >> BITS_PER_GROUP) +
				(unsigned int)*(Sum + i) + (unsigned int)*(mod + i);
			Sum[i] = (int)(carry & MAX_VALUE_LIMB);
		}
	}
}

/* Sum = Nbr1 + Nbr2 (mod m) */
void AddBigNbrModNB(const limb Nbr1[], const limb Nbr2[], limb Sum[], const limb m[], int nbrLen)
{
	unsigned int carry = 0;
	int borrow;
	int i;
	/* sum = Nbr1 + Nbr2 */
	for (i = 0; i < nbrLen; i++) {
		carry = (carry >> BITS_PER_GROUP) +
			(unsigned int)Nbr1[i] + (unsigned int)Nbr2[i];
		Sum[i] = (int)(carry & MAX_VALUE_LIMB);
	}
	/* sum -= m */
	borrow = 0;
	for (i = 0; i < nbrLen; i++)
	{
		borrow = (borrow >> BITS_PER_GROUP) +
			(unsigned int)Sum[i] - m[i];
		Sum[i] = (int)(borrow & MAX_VALUE_LIMB);
	}
	/* if Sum is -ve, Sum += m */
	if (carry < LIMB_RANGE && borrow < 0)
	{
		carry = 0;
		for (i = 0; i < nbrLen; i++)
		{
			carry = (carry >> BITS_PER_GROUP) +
				(unsigned int)Sum[i] + (unsigned int)m[i];
			Sum[i] = (int)(carry & MAX_VALUE_LIMB);
		}
	}
}


/* Diff = Nbr1-Nbr2 (mod m)*/
void SubtBigNbrModN(const limb Nbr1[], const limb Nbr2[], limb Diff[], const limb m[], int nbrLen)
{
	int i;
	int borrow = 0;
	/* Diff = Nbr1 - Nbr2 */
	for (i = 0; i < nbrLen; i++) {
		borrow = (borrow >> BITS_PER_GROUP) + Nbr1[i] - Nbr2[i];
		Diff[i] = (int)(borrow & MAX_VALUE_LIMB);
	}

	if (borrow < 0) {  /* if Diff < 0 then Diff += m */
		unsigned int carry = 0;
		for (i = 0; i < nbrLen; i++)
		{
			carry = (carry >> BITS_PER_GROUP) +
				(unsigned int)Diff[i] + (unsigned int)m[i];
			Diff[i] = (int)(carry & MAX_VALUE_LIMB);
		}
	}
}


/* product = factor1*factor2%mod */
static void smallmodmult(const int factor1, const int factor2, limb *product, int mod)
{
	if (mod < SMALL_NUMBER_BOUND)
	{
		*product = factor1 * factor2 % mod;
	}
	else
	{   // TestNbr has one limb but it is not small.
#ifdef _USING64BITS_
		*product = (int64_t)factor1 * factor2 % mod;
#else
		// Round up quotient.
		int quotient = (int)floor((double)factor1 * (double)factor2 / (double)mod + 0.5);
		int remainder = factor1 * factor2 - quotient * mod;
		if (remainder < 0)
		{    // Quotient was 1 more than expected. Adjust remainder.
			remainder += mod;
		}
		*product = remainder;
#endif
	}
}

// Multiply two numbers in Montgomery notation.
// see https://en.wikipedia.org/wiki/Montgomery_modular_multiplication
// For large numbers the REDC algorithm is:
// m <- ((T mod R)N') mod R
// t <- (T + mN) / R
// if t >= N then
//   return t - N
// else
//   return t
// end if

#ifdef _USING64BITS_
static void MontgomeryMult2(const limb pNbr1[], const limb pNbr2[], limb pProd[])
{
	int i;
	uint64_t Pr;
	int32_t borrow;
	uint32_t Nbr, MontDig;
	uint32_t Prod0, Prod1;
	Prod0 = Prod1 = 0;
	uint32_t TestNbr0 = TestNbr[0];
	uint32_t TestNbr1 = TestNbr[1];
	uint32_t Nbr2_0 = pNbr2[0];
	uint32_t Nbr2_1 = pNbr2[1];
	for (i = 0; i<2; i++)
	{
		Pr = (Nbr = pNbr1[i]) * (uint64_t)Nbr2_0 + (uint32_t)Prod0;
		MontDig = ((uint32_t)Pr * MontgomeryMultN[0]) & MAX_INT_NBR;
		Prod0 = (Pr = (((uint64_t)MontDig * TestNbr0 + Pr) >> BITS_PER_GROUP) +
			(uint64_t)MontDig * TestNbr1 + (uint64_t)Nbr * Nbr2_1 + (uint32_t)Prod1) & MAX_INT_NBR;
		Prod1 = (uint32_t)(Pr >> BITS_PER_GROUP);
	}
	if (Pr >= ((uint64_t)(TestNbr1 + 1) << BITS_PER_GROUP) || (Prod1 == TestNbr1 && Prod0 >= TestNbr0))
	{
		Prod0 = (borrow = (int32_t)Prod0 - (int32_t)TestNbr0) & MAX_INT_NBR;
		Prod1 = ((borrow >> BITS_PER_GROUP) + (int32_t)Prod1 - (int32_t)TestNbr1) & MAX_INT_NBR;
	}
	pProd[0] = Prod0;
	pProd[1] = Prod1;
}

static void MontgomeryMult3(const limb pNbr1[], const limb pNbr2[], limb pProd[])
{
	int i;
	uint64_t Pr;
	int32_t borrow;
	uint32_t Nbr, MontDig;
	uint32_t Prod0, Prod1, Prod2;
	Prod0 = Prod1 = Prod2 = 0;
	uint32_t TestNbr0 = TestNbr[0];
	uint32_t TestNbr1 = TestNbr[1];
	uint32_t TestNbr2 = TestNbr[2];
	uint32_t Nbr2_0 = pNbr2[0];
	uint32_t Nbr2_1 = pNbr2[1];
	uint32_t Nbr2_2 = pNbr2[2];
	for (i = 0; i<3; i++)
	{
		Pr = (Nbr = pNbr1[i]) * (uint64_t)Nbr2_0 + (uint32_t)Prod0;
		MontDig = ((uint32_t)Pr * MontgomeryMultN[0]) & MAX_INT_NBR;
		Prod0 = (Pr = (((uint64_t)MontDig * TestNbr0 + Pr) >> BITS_PER_GROUP) +
			(uint64_t)MontDig * TestNbr1 + (uint64_t)Nbr * Nbr2_1 + (uint32_t)Prod1) & MAX_INT_NBR;
		Prod1 = (Pr = (Pr >> BITS_PER_GROUP) +
			(uint64_t)MontDig * TestNbr2 + (uint64_t)Nbr * Nbr2_2 + (uint32_t)Prod2) & MAX_INT_NBR;
		Prod2 = (uint32_t)(Pr >> BITS_PER_GROUP);
	}
	if (Pr >= ((uint64_t)(TestNbr2 + 1) << BITS_PER_GROUP)
		|| (Prod2 == TestNbr2
			&& (Prod1 > TestNbr1 || (Prod1 == TestNbr1 && (Prod0 >= TestNbr0)))))
	{
		Prod0 = (borrow = (int32_t)Prod0 - (int32_t)TestNbr0) & MAX_INT_NBR;
		Prod1 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod1 - (int32_t)TestNbr1) & MAX_INT_NBR;
		Prod2 = ((borrow >> BITS_PER_GROUP) + (int32_t)Prod2 - (int32_t)TestNbr2) & MAX_INT_NBR;
	}
	pProd[0] = Prod0;
	pProd[1] = Prod1;
	pProd[2] = Prod2;
}

static void MontgomeryMult4(const limb pNbr1[], const limb pNbr2[], limb pProd[])
{
	int i;
	uint64_t Pr;
	int32_t borrow;
	uint32_t Nbr, MontDig;
	uint32_t Prod0, Prod1, Prod2, Prod3;
	Prod0 = Prod1 = Prod2 = Prod3 = 0;
	uint32_t TestNbr0 = TestNbr[0];
	uint32_t TestNbr1 = TestNbr[1];
	uint32_t TestNbr2 = TestNbr[2];
	uint32_t TestNbr3 = TestNbr[3];
	uint32_t Nbr2_0 = pNbr2[0];
	uint32_t Nbr2_1 = pNbr2[1];
	uint32_t Nbr2_2 = pNbr2[2];
	uint32_t Nbr2_3 = pNbr2[3];
	for (i = 0; i<4; i++)
	{
		Pr = (Nbr = pNbr1[i]) * (uint64_t)Nbr2_0 + (uint32_t)Prod0;
		MontDig = ((uint32_t)Pr * MontgomeryMultN[0]) & MAX_INT_NBR;
		Prod0 = (Pr = (((uint64_t)MontDig * TestNbr0 + Pr) >> BITS_PER_GROUP) +
			(uint64_t)MontDig * TestNbr1 + (uint64_t)Nbr * Nbr2_1 + (uint32_t)Prod1) & MAX_INT_NBR;
		Prod1 = (Pr = (Pr >> BITS_PER_GROUP) +
			(uint64_t)MontDig * TestNbr2 + (uint64_t)Nbr * Nbr2_2 + (uint32_t)Prod2) & MAX_INT_NBR;
		Prod2 = (Pr = (Pr >> BITS_PER_GROUP) +
			(uint64_t)MontDig * TestNbr3 + (uint64_t)Nbr * Nbr2_3 + (uint32_t)Prod3) & MAX_INT_NBR;
		Prod3 = (uint32_t)(Pr >> BITS_PER_GROUP);
	}
	if (Pr >= ((uint64_t)(TestNbr3 + 1) << BITS_PER_GROUP)
		|| (Prod3 == TestNbr3
			&& (Prod2 > TestNbr2
				|| (Prod2 == TestNbr2
					&& (Prod1 > TestNbr1 || (Prod1 == TestNbr1 && (Prod0 >= TestNbr0)))))))
	{
		Prod0 = (borrow = (int32_t)Prod0 - (int32_t)TestNbr0) & MAX_INT_NBR;
		Prod1 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod1 - (int32_t)TestNbr1) & MAX_INT_NBR;
		Prod2 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod2 - (int32_t)TestNbr2) & MAX_INT_NBR;
		Prod3 = ((borrow >> BITS_PER_GROUP) + (int32_t)Prod3 - (int32_t)TestNbr3) & MAX_INT_NBR;
	}
	pProd[0] = Prod0;
	pProd[1] = Prod1;
	pProd[2] = Prod2;
	pProd[3] = Prod3;
}

static void MontgomeryMult5(const limb pNbr1[], const limb pNbr2[], limb pProd[])
{
	int i;
	uint64_t Pr;
	int32_t borrow;
	uint32_t Nbr, MontDig;
	uint32_t Prod0, Prod1, Prod2, Prod3, Prod4;
	Prod0 = Prod1 = Prod2 = Prod3 = Prod4 = 0;
	uint32_t TestNbr0 = TestNbr[0];
	uint32_t TestNbr1 = TestNbr[1];
	uint32_t TestNbr2 = TestNbr[2];
	uint32_t TestNbr3 = TestNbr[3];
	uint32_t TestNbr4 = TestNbr[4];
	uint32_t Nbr2_0 = pNbr2[0];
	uint32_t Nbr2_1 = pNbr2[1];
	uint32_t Nbr2_2 = pNbr2[2];
	uint32_t Nbr2_3 = pNbr2[3];
	uint32_t Nbr2_4 = pNbr2[4];
	for (i = 0; i<5; i++)
	{
		Pr = (Nbr = pNbr1[i]) * (uint64_t)Nbr2_0 + (uint32_t)Prod0;
		MontDig = ((uint32_t)Pr * MontgomeryMultN[0]) & MAX_INT_NBR;
		Prod0 = (Pr = (((uint64_t)MontDig * TestNbr0 + Pr) >> BITS_PER_GROUP) +
			(uint64_t)MontDig * TestNbr1 + (uint64_t)Nbr * Nbr2_1 + (uint32_t)Prod1) & MAX_INT_NBR;
		Prod1 = (Pr = (Pr >> BITS_PER_GROUP) +
			(uint64_t)MontDig * TestNbr2 + (uint64_t)Nbr * Nbr2_2 + (uint32_t)Prod2) & MAX_INT_NBR;
		Prod2 = (Pr = (Pr >> BITS_PER_GROUP) +
			(uint64_t)MontDig * TestNbr3 + (uint64_t)Nbr * Nbr2_3 + (uint32_t)Prod3) & MAX_INT_NBR;
		Prod3 = (Pr = (Pr >> BITS_PER_GROUP) +
			(uint64_t)MontDig * TestNbr4 + (uint64_t)Nbr * Nbr2_4 + (uint32_t)Prod4) & MAX_INT_NBR;
		Prod4 = (uint32_t)(Pr >> BITS_PER_GROUP);
	}
	if (Pr >= ((uint64_t)(TestNbr4 + 1) << BITS_PER_GROUP)
		|| (Prod4 == TestNbr4
			&& (Prod3 > TestNbr3
				|| (Prod3 == TestNbr3
					&& (Prod2 > TestNbr2
						|| (Prod2 == TestNbr2
							&& (Prod1 > TestNbr1 || (Prod1 == TestNbr1 && (Prod0 >= TestNbr0)))))))))
	{
		Prod0 = (borrow = (int32_t)Prod0 - (int32_t)TestNbr0) & MAX_INT_NBR;
		Prod1 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod1 - (int32_t)TestNbr1) & MAX_INT_NBR;
		Prod2 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod2 - (int32_t)TestNbr2) & MAX_INT_NBR;
		Prod3 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod3 - (int32_t)TestNbr3) & MAX_INT_NBR;
		Prod4 = ((borrow >> BITS_PER_GROUP) + (int32_t)Prod4 - (int32_t)TestNbr4) & MAX_INT_NBR;
	}
	pProd[0] = Prod0;
	pProd[1] = Prod1;
	pProd[2] = Prod2;
	pProd[3] = Prod3;
	pProd[4] = Prod4;
}

static void MontgomeryMult6(const limb pNbr1[], const limb pNbr2[], limb pProd[])
{
	int i;
	uint64_t Pr;
	int32_t borrow;
	uint32_t Nbr, MontDig;
	uint32_t Prod0, Prod1, Prod2, Prod3, Prod4, Prod5;
	Prod0 = Prod1 = Prod2 = Prod3 = Prod4 = Prod5 = 0;
	uint32_t TestNbr0 = TestNbr[0];
	uint32_t TestNbr1 = TestNbr[1];
	uint32_t TestNbr2 = TestNbr[2];
	uint32_t TestNbr3 = TestNbr[3];
	uint32_t TestNbr4 = TestNbr[4];
	uint32_t TestNbr5 = TestNbr[5];
	uint32_t Nbr2_0 = pNbr2[0];
	uint32_t Nbr2_1 = pNbr2[1];
	uint32_t Nbr2_2 = pNbr2[2];
	uint32_t Nbr2_3 = pNbr2[3];
	uint32_t Nbr2_4 = pNbr2[4];
	uint32_t Nbr2_5 = pNbr2[5];
	for (i = 0; i<6; i++)
	{
		Pr = (Nbr = pNbr1[i]) * (uint64_t)Nbr2_0 + (uint32_t)Prod0;
		MontDig = ((uint32_t)Pr * MontgomeryMultN[0]) & MAX_INT_NBR;
		Prod0 = (Pr = (((uint64_t)MontDig * TestNbr0 + Pr) >> BITS_PER_GROUP) +
			(uint64_t)MontDig * TestNbr1 + (uint64_t)Nbr * Nbr2_1 + (uint32_t)Prod1) & MAX_INT_NBR;
		Prod1 = (Pr = (Pr >> BITS_PER_GROUP) +
			(uint64_t)MontDig * TestNbr2 + (uint64_t)Nbr * Nbr2_2 + (uint32_t)Prod2) & MAX_INT_NBR;
		Prod2 = (Pr = (Pr >> BITS_PER_GROUP) +
			(uint64_t)MontDig * TestNbr3 + (uint64_t)Nbr * Nbr2_3 + (uint32_t)Prod3) & MAX_INT_NBR;
		Prod3 = (Pr = (Pr >> BITS_PER_GROUP) +
			(uint64_t)MontDig * TestNbr4 + (uint64_t)Nbr * Nbr2_4 + (uint32_t)Prod4) & MAX_INT_NBR;
		Prod4 = (Pr = (Pr >> BITS_PER_GROUP) +
			(uint64_t)MontDig * TestNbr5 + (uint64_t)Nbr * Nbr2_5 + (uint32_t)Prod5) & MAX_INT_NBR;
		Prod5 = (uint32_t)(Pr >> BITS_PER_GROUP);
	}
	if (Pr >= ((uint64_t)(TestNbr5 + 1) << BITS_PER_GROUP)
		|| (Prod5 == TestNbr5
			&& (Prod4 > TestNbr4
				|| (Prod4 == TestNbr4
					&& (Prod3 > TestNbr3
						|| (Prod3 == TestNbr3
							&& (Prod2 > TestNbr2
								|| (Prod2 == TestNbr2
									&& (Prod1 > TestNbr1
										|| (Prod1 == TestNbr1
											&& (Prod0 >= TestNbr0)))))))))))
	{
		Prod0 = (borrow = (int32_t)Prod0 - (int32_t)TestNbr0) & MAX_INT_NBR;
		Prod1 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod1 - (int32_t)TestNbr1) & MAX_INT_NBR;
		Prod2 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod2 - (int32_t)TestNbr2) & MAX_INT_NBR;
		Prod3 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod3 - (int32_t)TestNbr3) & MAX_INT_NBR;
		Prod4 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod4 - (int32_t)TestNbr4) & MAX_INT_NBR;
		Prod5 = ((borrow >> BITS_PER_GROUP) + (int32_t)Prod5 - (int32_t)TestNbr5) & MAX_INT_NBR;
	}
	pProd[0] = Prod0;
	pProd[1] = Prod1;
	pProd[2] = Prod2;
	pProd[3] = Prod3;
	pProd[4] = Prod4;
	pProd[5] = Prod5;
}

static void MontgomeryMult7(const limb pNbr1[], const limb pNbr2[], limb pProd[])
{
	int i;
	uint64_t Pr;
	int32_t borrow;
	uint32_t Nbr, MontDig;
	int Prod0, Prod1, Prod2, Prod3, Prod4, Prod5, Prod6;
	Prod0 = Prod1 = Prod2 = Prod3 = Prod4 = Prod5 = Prod6 = 0;
	int TestNbr0 = (int)TestNbr[0];
	int TestNbr1 = (int)TestNbr[1];
	int TestNbr2 = (int)TestNbr[2];
	int TestNbr3 = (int)TestNbr[3];
	int TestNbr4 = (int)TestNbr[4];
	int TestNbr5 = (int)TestNbr[5];
	int TestNbr6 = (int)TestNbr[6];
	int Nbr2_0 = (int)pNbr2[0];
	int Nbr2_1 = (int)pNbr2[1];
	int Nbr2_2 = (int)pNbr2[2];
	int Nbr2_3 = (int)pNbr2[3];
	int Nbr2_4 = (int)pNbr2[4];
	int Nbr2_5 = (int)pNbr2[5];
	int Nbr2_6 = (int)pNbr2[6];
	for (i = 0; i<7; i++)
	{
		Pr = (Nbr = pNbr1[i]) * (uint64_t)Nbr2_0 + (uint32_t)Prod0;
		MontDig = ((uint32_t)Pr * MontgomeryMultN[0]) & MAX_INT_NBR;
		Prod0 = (Pr = (((uint64_t)MontDig * TestNbr0 + Pr) >> BITS_PER_GROUP) +
			(uint64_t)MontDig * TestNbr1 + (uint64_t)Nbr * Nbr2_1 + (uint32_t)Prod1) & MAX_INT_NBR;
		Prod1 = (Pr = (Pr >> BITS_PER_GROUP) +
			(uint64_t)MontDig * TestNbr2 + (uint64_t)Nbr * Nbr2_2 + (uint32_t)Prod2) & MAX_INT_NBR;
		Prod2 = (Pr = (Pr >> BITS_PER_GROUP) +
			(uint64_t)MontDig * TestNbr3 + (uint64_t)Nbr * Nbr2_3 + (uint32_t)Prod3) & MAX_INT_NBR;
		Prod3 = (Pr = (Pr >> BITS_PER_GROUP) +
			(uint64_t)MontDig * TestNbr4 + (uint64_t)Nbr * Nbr2_4 + (uint32_t)Prod4) & MAX_INT_NBR;
		Prod4 = (Pr = (Pr >> BITS_PER_GROUP) +
			(uint64_t)MontDig * TestNbr5 + (uint64_t)Nbr * Nbr2_5 + (uint32_t)Prod5) & MAX_INT_NBR;
		Prod5 = (Pr = (Pr >> BITS_PER_GROUP) +
			(uint64_t)MontDig * TestNbr6 + (uint64_t)Nbr * Nbr2_6 + (uint32_t)Prod6) & MAX_INT_NBR;
		Prod6 = (uint32_t)(Pr >> BITS_PER_GROUP);
	}
	if (Pr >= ((uint64_t)(TestNbr6 + 1) << BITS_PER_GROUP)
		|| (Prod6 == TestNbr6
			&& (Prod5 > TestNbr5
				|| (Prod5 == TestNbr5
					&& (Prod4 > TestNbr4
						|| (Prod4 == TestNbr4
							&& (Prod3 > TestNbr3
								|| (Prod3 == TestNbr3
									&& (Prod2 > TestNbr2
										|| (Prod2 == TestNbr2
											&& (Prod1 > TestNbr1
												|| (Prod1 == TestNbr1
													&& (Prod0 >= TestNbr0)))))))))))))
	{
		Prod0 = (borrow = (int32_t)Prod0 - (int32_t)TestNbr0) & MAX_INT_NBR;
		Prod1 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod1 - (int32_t)TestNbr1) & MAX_INT_NBR;
		Prod2 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod2 - (int32_t)TestNbr2) & MAX_INT_NBR;
		Prod3 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod3 - (int32_t)TestNbr3) & MAX_INT_NBR;
		Prod4 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod4 - (int32_t)TestNbr4) & MAX_INT_NBR;
		Prod5 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod5 - (int32_t)TestNbr5) & MAX_INT_NBR;
		Prod6 = ((borrow >> BITS_PER_GROUP) + (int32_t)Prod6 - (int32_t)TestNbr6) & MAX_INT_NBR;
	}
	pProd[0] = Prod0;
	pProd[1] = Prod1;
	pProd[2] = Prod2;
	pProd[3] = Prod3;
	pProd[4] = Prod4;
	pProd[5] = Prod5;
	pProd[6] = Prod6;
}

static void MontgomeryMult8(const limb pNbr1[], const limb pNbr2[], limb pProd[])
{
	int i;
	uint64_t Pr;
	int32_t borrow;
	uint32_t Nbr, MontDig;
	uint32_t Prod0, Prod1, Prod2, Prod3, Prod4, Prod5, Prod6, Prod7;
	Prod0 = Prod1 = Prod2 = Prod3 = Prod4 = Prod5 = Prod6 = Prod7 = 0;
	uint32_t TestNbr0 = TestNbr[0];
	uint32_t TestNbr1 = TestNbr[1];
	uint32_t TestNbr2 = TestNbr[2];
	uint32_t TestNbr3 = TestNbr[3];
	uint32_t TestNbr4 = TestNbr[4];
	uint32_t TestNbr5 = TestNbr[5];
	uint32_t TestNbr6 = TestNbr[6];
	uint32_t TestNbr7 = TestNbr[7];
	uint32_t Nbr2_0 = pNbr2[0];
	uint32_t Nbr2_1 = pNbr2[1];
	uint32_t Nbr2_2 = pNbr2[2];
	uint32_t Nbr2_3 = pNbr2[3];
	uint32_t Nbr2_4 = pNbr2[4];
	uint32_t Nbr2_5 = pNbr2[5];
	uint32_t Nbr2_6 = pNbr2[6];
	uint32_t Nbr2_7 = pNbr2[7];
	for (i = 0; i<8; i++)
	{
		Pr = (Nbr = pNbr1[i]) * (uint64_t)Nbr2_0 + (uint32_t)Prod0;
		MontDig = ((uint32_t)Pr * MontgomeryMultN[0]) & MAX_INT_NBR;
		Prod0 = (Pr = (((uint64_t)MontDig * TestNbr0 + Pr) >> BITS_PER_GROUP) +
			(uint64_t)MontDig * TestNbr1 + (uint64_t)Nbr * Nbr2_1 + (uint32_t)Prod1) & MAX_INT_NBR;
		Prod1 = (Pr = (Pr >> BITS_PER_GROUP) +
			(uint64_t)MontDig * TestNbr2 + (uint64_t)Nbr * Nbr2_2 + (uint32_t)Prod2) & MAX_INT_NBR;
		Prod2 = (Pr = (Pr >> BITS_PER_GROUP) +
			(uint64_t)MontDig * TestNbr3 + (uint64_t)Nbr * Nbr2_3 + (uint32_t)Prod3) & MAX_INT_NBR;
		Prod3 = (Pr = (Pr >> BITS_PER_GROUP) +
			(uint64_t)MontDig * TestNbr4 + (uint64_t)Nbr * Nbr2_4 + (uint32_t)Prod4) & MAX_INT_NBR;
		Prod4 = (Pr = (Pr >> BITS_PER_GROUP) +
			(uint64_t)MontDig * TestNbr5 + (uint64_t)Nbr * Nbr2_5 + (uint32_t)Prod5) & MAX_INT_NBR;
		Prod5 = (Pr = (Pr >> BITS_PER_GROUP) +
			(uint64_t)MontDig * TestNbr6 + (uint64_t)Nbr * Nbr2_6 + (uint32_t)Prod6) & MAX_INT_NBR;
		Prod6 = (Pr = (Pr >> BITS_PER_GROUP) +
			(uint64_t)MontDig * TestNbr7 + (uint64_t)Nbr * Nbr2_7 + (uint32_t)Prod7) & MAX_INT_NBR;
		Prod7 = (uint32_t)(Pr >> BITS_PER_GROUP);
	}
	if (Pr >= ((uint64_t)(TestNbr7 + 1) << BITS_PER_GROUP)
		|| (Prod7 == TestNbr7
			&& (Prod6 > TestNbr6
				|| (Prod6 == TestNbr6
					&& (Prod5 > TestNbr5
						|| (Prod5 == TestNbr5
							&& (Prod4 > TestNbr4
								|| (Prod4 == TestNbr4
									&& (Prod3 > TestNbr3
										|| (Prod3 == TestNbr3
											&& (Prod2 > TestNbr2
												|| (Prod2 == TestNbr2
													&& (Prod1 > TestNbr1
														|| (Prod1 == TestNbr1
															&& (Prod0 >= TestNbr0)))))))))))))))
	{
		Prod0 = (borrow = (int32_t)Prod0 - (int32_t)TestNbr0) & MAX_INT_NBR;
		Prod1 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod1 - (int32_t)TestNbr1) & MAX_INT_NBR;
		Prod2 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod2 - (int32_t)TestNbr2) & MAX_INT_NBR;
		Prod3 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod3 - (int32_t)TestNbr3) & MAX_INT_NBR;
		Prod4 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod4 - (int32_t)TestNbr4) & MAX_INT_NBR;
		Prod5 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod5 - (int32_t)TestNbr5) & MAX_INT_NBR;
		Prod6 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod6 - (int32_t)TestNbr6) & MAX_INT_NBR;
		Prod7 = ((borrow >> BITS_PER_GROUP) + (int32_t)Prod7 - (int32_t)TestNbr7) & MAX_INT_NBR;
	}
	pProd[0] = (int)Prod0;
	pProd[1] = (int)Prod1;
	pProd[2] = (int)Prod2;
	pProd[3] = (int)Prod3;
	pProd[4] = (int)Prod4;
	pProd[5] = (int)Prod5;
	pProd[6] = (int)Prod6;
	pProd[7] = (int)Prod7;
}
static void MontgomeryMult9(const limb pNbr1[], const limb pNbr2[], limb pProd[])
{
	int i;
	uint64_t Pr;
	int32_t borrow;
	uint32_t Nbr, MontDig;
	uint32_t Prod0, Prod1, Prod2, Prod3, Prod4, Prod5, Prod6, Prod7, Prod8;
	Prod0 =
		Prod1 = Prod2 = Prod3 = Prod4 = Prod5 = Prod6 = Prod7 = Prod8 = 0;
	uint32_t TestNbr0 = TestNbr[0];
	uint32_t TestNbr1 = TestNbr[1];
	uint32_t TestNbr2 = TestNbr[2];
	uint32_t TestNbr3 = TestNbr[3];
	uint32_t TestNbr4 = TestNbr[4];
	uint32_t TestNbr5 = TestNbr[5];
	uint32_t TestNbr6 = TestNbr[6];
	uint32_t TestNbr7 = TestNbr[7];
	uint32_t TestNbr8 = TestNbr[8];
	uint32_t Nbr2_0 = pNbr2[0];
	uint32_t Nbr2_1 = pNbr2[1];
	uint32_t Nbr2_2 = pNbr2[2];
	uint32_t Nbr2_3 = pNbr2[3];
	uint32_t Nbr2_4 = pNbr2[4];
	uint32_t Nbr2_5 = pNbr2[5];
	uint32_t Nbr2_6 = pNbr2[6];
	uint32_t Nbr2_7 = pNbr2[7];
	uint32_t Nbr2_8 = pNbr2[8];
	for (i = 0; i<9; i++)
	{
		Pr = (Nbr = pNbr1[i]) * (uint64_t)Nbr2_0 + (uint32_t)Prod0;
		MontDig = ((uint32_t)Pr * MontgomeryMultN[0]) & MAX_INT_NBR;
		Prod0 = (Pr = (((uint64_t)MontDig * TestNbr0 + Pr) >> BITS_PER_GROUP) +
			(uint64_t)MontDig * TestNbr1 + (uint64_t)Nbr * Nbr2_1 + (uint32_t)Prod1) & MAX_INT_NBR;
		Prod1 = (Pr = (Pr >> BITS_PER_GROUP) +
			(uint64_t)MontDig * TestNbr2 + (uint64_t)Nbr * Nbr2_2 + (uint32_t)Prod2) & MAX_INT_NBR;
		Prod2 = (Pr = (Pr >> BITS_PER_GROUP) +
			(uint64_t)MontDig * TestNbr3 + (uint64_t)Nbr * Nbr2_3 + (uint32_t)Prod3) & MAX_INT_NBR;
		Prod3 = (Pr = (Pr >> BITS_PER_GROUP) +
			(uint64_t)MontDig * TestNbr4 + (uint64_t)Nbr * Nbr2_4 + (uint32_t)Prod4) & MAX_INT_NBR;
		Prod4 = (Pr = (Pr >> BITS_PER_GROUP) +
			(uint64_t)MontDig * TestNbr5 + (uint64_t)Nbr * Nbr2_5 + (uint32_t)Prod5) & MAX_INT_NBR;
		Prod5 = (Pr = (Pr >> BITS_PER_GROUP) +
			(uint64_t)MontDig * TestNbr6 + (uint64_t)Nbr * Nbr2_6 + (uint32_t)Prod6) & MAX_INT_NBR;
		Prod6 = (Pr = (Pr >> BITS_PER_GROUP) +
			(uint64_t)MontDig * TestNbr7 + (uint64_t)Nbr * Nbr2_7 + (uint32_t)Prod7) & MAX_INT_NBR;
		Prod7 = (Pr = (Pr >> BITS_PER_GROUP) +
			(uint64_t)MontDig * TestNbr8 + (uint64_t)Nbr * Nbr2_8 + (uint32_t)Prod8) & MAX_INT_NBR;
		Prod8 = (uint32_t)(Pr >> BITS_PER_GROUP);
	}
	if (Pr >= ((uint64_t)(TestNbr8 + 1) << BITS_PER_GROUP)
		|| (Prod8 == TestNbr8
			&& (Prod7 > TestNbr7
				|| (Prod7 == TestNbr7
					&& (Prod6 > TestNbr6
						|| (Prod6 == TestNbr6
							&& (Prod5 > TestNbr5
								|| (Prod5 == TestNbr5
									&& (Prod4 > TestNbr4
										|| (Prod4 == TestNbr4
											&& (Prod3 > TestNbr3
												|| (Prod3 == TestNbr3
													&& (Prod2 > TestNbr2
														|| (Prod2 == TestNbr2
															&& (Prod1 > TestNbr1
																|| (Prod1 == TestNbr1
																	&& (Prod0 >= TestNbr0)))))))))))))))))
	{
		Prod0 = (borrow = (int32_t)Prod0 - (int32_t)TestNbr0) & MAX_INT_NBR;
		Prod1 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod1 - (int32_t)TestNbr1) & MAX_INT_NBR;
		Prod2 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod2 - (int32_t)TestNbr2) & MAX_INT_NBR;
		Prod3 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod3 - (int32_t)TestNbr3) & MAX_INT_NBR;
		Prod4 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod4 - (int32_t)TestNbr4) & MAX_INT_NBR;
		Prod5 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod5 - (int32_t)TestNbr5) & MAX_INT_NBR;
		Prod6 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod6 - (int32_t)TestNbr6) & MAX_INT_NBR;
		Prod7 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod7 - (int32_t)TestNbr7) & MAX_INT_NBR;
		Prod8 = ((borrow >> BITS_PER_GROUP) + (int32_t)Prod8 - (int32_t)TestNbr8) & MAX_INT_NBR;
	}
	pProd[0] = Prod0;
	pProd[1] = Prod1;
	pProd[2] = Prod2;
	pProd[3] = Prod3;
	pProd[4] = Prod4;
	pProd[5] = Prod5;
	pProd[6] = Prod6;
	pProd[7] = Prod7;
	pProd[8] = Prod8;
}

static void MontgomeryMult10(const limb pNbr1[], const limb pNbr2[], limb pProd[])
{
	int i;
	uint64_t Pr;
	int32_t borrow;
	uint32_t Nbr, MontDig;
	uint32_t Prod0, Prod1, Prod2, Prod3, Prod4, Prod5, Prod6, Prod7, Prod8,
		Prod9;
	Prod0 = Prod1 =
		Prod2 = Prod3 = Prod4 = Prod5 = Prod6 = Prod7 = Prod8 = Prod9 = 0;
	uint32_t TestNbr0 = TestNbr[0];
	uint32_t TestNbr1 = TestNbr[1];
	uint32_t TestNbr2 = TestNbr[2];
	uint32_t TestNbr3 = TestNbr[3];
	uint32_t TestNbr4 = TestNbr[4];
	uint32_t TestNbr5 = TestNbr[5];
	uint32_t TestNbr6 = TestNbr[6];
	uint32_t TestNbr7 = TestNbr[7];
	uint32_t TestNbr8 = TestNbr[8];
	uint32_t TestNbr9 = TestNbr[9];
	uint32_t Nbr2_0 = pNbr2[0];
	uint32_t Nbr2_1 = pNbr2[1];
	uint32_t Nbr2_2 = pNbr2[2];
	uint32_t Nbr2_3 = pNbr2[3];
	uint32_t Nbr2_4 = pNbr2[4];
	uint32_t Nbr2_5 = pNbr2[5];
	uint32_t Nbr2_6 = pNbr2[6];
	uint32_t Nbr2_7 = pNbr2[7];
	uint32_t Nbr2_8 = pNbr2[8];
	uint32_t Nbr2_9 = pNbr2[9];
	for (i = 0; i<10; i++)
	{
		Pr = (Nbr = pNbr1[i]) * (uint64_t)Nbr2_0 + (uint32_t)Prod0;
		MontDig = ((uint32_t)Pr * MontgomeryMultN[0]) & MAX_INT_NBR;
		Prod0 = (Pr = (((uint64_t)MontDig * TestNbr0 + Pr) >> BITS_PER_GROUP) +
			(uint64_t)MontDig * TestNbr1 + (uint64_t)Nbr * Nbr2_1 + (uint32_t)Prod1) & MAX_INT_NBR;
		Prod1 = (Pr = (Pr >> BITS_PER_GROUP) +
			(uint64_t)MontDig * TestNbr2 + (uint64_t)Nbr * Nbr2_2 + (uint32_t)Prod2) & MAX_INT_NBR;
		Prod2 = (Pr = (Pr >> BITS_PER_GROUP) +
			(uint64_t)MontDig * TestNbr3 + (uint64_t)Nbr * Nbr2_3 + (uint32_t)Prod3) & MAX_INT_NBR;
		Prod3 = (Pr = (Pr >> BITS_PER_GROUP) +
			(uint64_t)MontDig * TestNbr4 + (uint64_t)Nbr * Nbr2_4 + (uint32_t)Prod4) & MAX_INT_NBR;
		Prod4 = (Pr = (Pr >> BITS_PER_GROUP) +
			(uint64_t)MontDig * TestNbr5 + (uint64_t)Nbr * Nbr2_5 + (uint32_t)Prod5) & MAX_INT_NBR;
		Prod5 = (Pr = (Pr >> BITS_PER_GROUP) +
			(uint64_t)MontDig * TestNbr6 + (uint64_t)Nbr * Nbr2_6 + (uint32_t)Prod6) & MAX_INT_NBR;
		Prod6 = (Pr = (Pr >> BITS_PER_GROUP) +
			(uint64_t)MontDig * TestNbr7 + (uint64_t)Nbr * Nbr2_7 + (uint32_t)Prod7) & MAX_INT_NBR;
		Prod7 = (Pr = (Pr >> BITS_PER_GROUP) +
			(uint64_t)MontDig * TestNbr8 + (uint64_t)Nbr * Nbr2_8 + (uint32_t)Prod8) & MAX_INT_NBR;
		Prod8 = (Pr = (Pr >> BITS_PER_GROUP) +
			(uint64_t)MontDig * TestNbr9 + (uint64_t)Nbr * Nbr2_9 + (uint32_t)Prod9) & MAX_INT_NBR;
		Prod9 = (uint32_t)(Pr >> BITS_PER_GROUP);
	}
	if (Pr >= ((uint64_t)(TestNbr9 + 1) << BITS_PER_GROUP)
		|| (Prod9 == TestNbr9
			&& (Prod8 > TestNbr8
				|| (Prod8 == TestNbr8
					&& (Prod7 > TestNbr7
						|| (Prod7 == TestNbr7
							&& (Prod6 > TestNbr6
								|| (Prod6 == TestNbr6
									&& (Prod5 > TestNbr5
										|| (Prod5 == TestNbr5
											&& (Prod4 > TestNbr4
												|| (Prod4 == TestNbr4
													&& (Prod3 > TestNbr3
														|| (Prod3 == TestNbr3
															&& (Prod2 > TestNbr2
																|| (Prod2 == TestNbr2
																	&& (Prod1 > TestNbr1
																		|| (Prod1 == TestNbr1
																			&& (Prod0 >= TestNbr0)))))))))))))))))))
	{
		Prod0 = (borrow = (int32_t)Prod0 - (int32_t)TestNbr0) & MAX_INT_NBR;
		Prod1 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod1 - (int32_t)TestNbr1) & MAX_INT_NBR;
		Prod2 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod2 - (int32_t)TestNbr2) & MAX_INT_NBR;
		Prod3 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod3 - (int32_t)TestNbr3) & MAX_INT_NBR;
		Prod4 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod4 - (int32_t)TestNbr4) & MAX_INT_NBR;
		Prod5 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod5 - (int32_t)TestNbr5) & MAX_INT_NBR;
		Prod6 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod6 - (int32_t)TestNbr6) & MAX_INT_NBR;
		Prod7 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod7 - (int32_t)TestNbr7) & MAX_INT_NBR;
		Prod8 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod8 - (int32_t)TestNbr8) & MAX_INT_NBR;
		Prod9 = ((borrow >> BITS_PER_GROUP) + (int32_t)Prod9 - (int32_t)TestNbr9) & MAX_INT_NBR;
	}
	pProd[0] = Prod0;
	pProd[1] = Prod1;
	pProd[2] = Prod2;
	pProd[3] = Prod3;
	pProd[4] = Prod4;
	pProd[5] = Prod5;
	pProd[6] = Prod6;
	pProd[7] = Prod7;
	pProd[8] = Prod8;
	pProd[9] = Prod9;
}

static void MontgomeryMult11(const limb pNbr1[], const limb pNbr2[], limb pProd[])
{
	int i;
	uint64_t Pr;
	int32_t borrow;
	uint32_t Nbr, MontDig;
	uint32_t Prod0, Prod1, Prod2, Prod3, Prod4, Prod5, Prod6, Prod7, Prod8,
		Prod9, Prod10;
	Prod0 = Prod1 = Prod2 = Prod3 =
		Prod4 = Prod5 = Prod6 = Prod7 = Prod8 = Prod9 = Prod10 = 0;
	uint32_t TestNbr0 = TestNbr[0];
	uint32_t TestNbr1 = TestNbr[1];
	uint32_t TestNbr2 = TestNbr[2];
	uint32_t TestNbr3 = TestNbr[3];
	uint32_t TestNbr4 = TestNbr[4];
	uint32_t TestNbr5 = TestNbr[5];
	uint32_t TestNbr6 = TestNbr[6];
	uint32_t TestNbr7 = TestNbr[7];
	uint32_t TestNbr8 = TestNbr[8];
	uint32_t TestNbr9 = TestNbr[9];
	uint32_t TestNbr10 = TestNbr[10];
	uint32_t Nbr2_0 = pNbr2[0];
	uint32_t Nbr2_1 = pNbr2[1];
	uint32_t Nbr2_2 = pNbr2[2];
	uint32_t Nbr2_3 = pNbr2[3];
	uint32_t Nbr2_4 = pNbr2[4];
	uint32_t Nbr2_5 = pNbr2[5];
	uint32_t Nbr2_6 = pNbr2[6];
	uint32_t Nbr2_7 = pNbr2[7];
	uint32_t Nbr2_8 = pNbr2[8];
	uint32_t Nbr2_9 = pNbr2[9];
	uint32_t Nbr2_10 = pNbr2[10];
	for (i = 0; i<11; i++)
	{
		Pr = (Nbr = pNbr1[i]) * (uint64_t)Nbr2_0 + (uint32_t)Prod0;
		MontDig = ((uint32_t)Pr * MontgomeryMultN[0]) & MAX_INT_NBR;
		Prod0 = (Pr = (((uint64_t)MontDig * TestNbr0 + Pr) >> BITS_PER_GROUP) +
			(uint64_t)MontDig * TestNbr1 + (uint64_t)Nbr * Nbr2_1 + (uint32_t)Prod1) & MAX_INT_NBR;
		Prod1 = (Pr = (Pr >> BITS_PER_GROUP) +
			(uint64_t)MontDig * TestNbr2 + (uint64_t)Nbr * Nbr2_2 + (uint32_t)Prod2) & MAX_INT_NBR;
		Prod2 = (Pr = (Pr >> BITS_PER_GROUP) +
			(uint64_t)MontDig * TestNbr3 + (uint64_t)Nbr * Nbr2_3 + (uint32_t)Prod3) & MAX_INT_NBR;
		Prod3 = (Pr = (Pr >> BITS_PER_GROUP) +
			(uint64_t)MontDig * TestNbr4 + (uint64_t)Nbr * Nbr2_4 + (uint32_t)Prod4) & MAX_INT_NBR;
		Prod4 = (Pr = (Pr >> BITS_PER_GROUP) +
			(uint64_t)MontDig * TestNbr5 + (uint64_t)Nbr * Nbr2_5 + (uint32_t)Prod5) & MAX_INT_NBR;
		Prod5 = (Pr = (Pr >> BITS_PER_GROUP) +
			(uint64_t)MontDig * TestNbr6 + (uint64_t)Nbr * Nbr2_6 + (uint32_t)Prod6) & MAX_INT_NBR;
		Prod6 = (Pr = (Pr >> BITS_PER_GROUP) +
			(uint64_t)MontDig * TestNbr7 + (uint64_t)Nbr * Nbr2_7 + (uint32_t)Prod7) & MAX_INT_NBR;
		Prod7 = (Pr = (Pr >> BITS_PER_GROUP) +
			(uint64_t)MontDig * TestNbr8 + (uint64_t)Nbr * Nbr2_8 + (uint32_t)Prod8) & MAX_INT_NBR;
		Prod8 = (Pr = (Pr >> BITS_PER_GROUP) +
			(uint64_t)MontDig * TestNbr9 + (uint64_t)Nbr * Nbr2_9 + (uint32_t)Prod9) & MAX_INT_NBR;
		Prod9 = (Pr = (Pr >> BITS_PER_GROUP) +
			(uint64_t)MontDig * TestNbr10 + (uint64_t)Nbr * Nbr2_10 + (uint32_t)Prod10) & MAX_INT_NBR;
		Prod10 = (uint32_t)(Pr >> BITS_PER_GROUP);
	}
	if (Pr >= ((uint64_t)(TestNbr10 + 1) << BITS_PER_GROUP)
		|| (Prod10 == TestNbr10
			&& (Prod9 > TestNbr9
				|| (Prod9 == TestNbr9
					&& (Prod8 > TestNbr8
						|| (Prod8 == TestNbr8
							&& (Prod7 > TestNbr7
								|| (Prod7 == TestNbr7
									&& (Prod6 > TestNbr6
										|| (Prod6 == TestNbr6
											&& (Prod5 > TestNbr5
												|| (Prod5 == TestNbr5
													&& (Prod4 > TestNbr4
														|| (Prod4 == TestNbr4
															&& (Prod3 > TestNbr3
																|| (Prod3 == TestNbr3
																	&& (Prod2 > TestNbr2
																		|| (Prod2 == TestNbr2
																			&& (Prod1 > TestNbr1
																				|| (Prod1 == TestNbr1
																					&& (Prod0 >= TestNbr0)))))))))))))))))))))
	{
		Prod0 = (borrow = (int32_t)Prod0 - (int32_t)TestNbr0) & MAX_INT_NBR;
		Prod1 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod1 - (int32_t)TestNbr1) & MAX_INT_NBR;
		Prod2 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod2 - (int32_t)TestNbr2) & MAX_INT_NBR;
		Prod3 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod3 - (int32_t)TestNbr3) & MAX_INT_NBR;
		Prod4 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod4 - (int32_t)TestNbr4) & MAX_INT_NBR;
		Prod5 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod5 - (int32_t)TestNbr5) & MAX_INT_NBR;
		Prod6 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod6 - (int32_t)TestNbr6) & MAX_INT_NBR;
		Prod7 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod7 - (int32_t)TestNbr7) & MAX_INT_NBR;
		Prod8 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod8 - (int32_t)TestNbr8) & MAX_INT_NBR;
		Prod9 = (borrow = (borrow >> BITS_PER_GROUP) + (int32_t)Prod9 - (int32_t)TestNbr9) & MAX_INT_NBR;
		Prod10 = ((borrow >> BITS_PER_GROUP) + (int32_t)Prod10 - (int32_t)TestNbr10) & MAX_INT_NBR;
	}
	pProd[0] = Prod0;
	pProd[1] = Prod1;
	pProd[2] = Prod2;
	pProd[3] = Prod3;
	pProd[4] = Prod4;
	pProd[5] = Prod5;
	pProd[6] = Prod6;
	pProd[7] = Prod7;
	pProd[8] = Prod8;
	pProd[9] = Prod9;
	pProd[10] = Prod10;
}
#endif

/* product = factor1*factor2 (mod TestNbr) - values in Montgomery format
uses global variables powerOf2Exponent, NumberLength, TestNbr */
void modmult(const limb factor1[], const limb factor2[], limb product[])
{
	limb carry;
	int count;
	limb Prod[13];
	unsigned int cy;
	int index;
//#ifdef __EMSCRIPTEN__
	if (modmultCallback != nullptr)
	{
		modmultCallback();  // display status
		lModularMult++;  // increase counter used to control status display
	}
//#endif
	if (powerOf2Exponent != 0) {    // TestNbr is a power of 2; cannot use Montgomery multiplication.
		//LimbsToBigInteger(factor1, tmpNum, NumberLength);  // tmpNum = factor1
		//LimbsToBigInteger(factor2, tmpDen, NumberLength);  // tmpDen = factor2
		//tmpNum = tmpNum*tmpDen; //BigIntMultiply(tmpNum, tmpDen, tmpNum);
		multiply(factor1, factor2, product, NumberLength, NULL);
		//BigIntegerToLimbs(product, tmpNum, NumberLength);  // product = tmpNum = factor1*factor2
		*(product + powerOf2Exponent / BITS_PER_GROUP) &= 
			(1 << (powerOf2Exponent % BITS_PER_GROUP)) - 1;
		return;
	}

	if (NumberLength == 1) {
		smallmodmult(factor1[0], factor2[0], product, TestNbr[0]);
		return;
	}
	if (NumberLength <= 12) {     // Small numbers.
		int i, j;
#ifdef _USING64BITS_
		int32_t MontDig, Nbr;
		int64_t Pr;
		switch (NumberLength)
		{
		case 2:
			MontgomeryMult2(factor1, factor2, product);
			return;
		case 3:
			MontgomeryMult3(factor1, factor2, product);
			return;
		case 4:
			MontgomeryMult4(factor1, factor2, product);
			return;
		case 5:
			MontgomeryMult5(factor1, factor2, product);
			return;
		case 6:
			MontgomeryMult6(factor1, factor2, product);
			return;
		case 7:
			MontgomeryMult7(factor1, factor2, product);
			return;
		case 8:
			MontgomeryMult8(factor1, factor2, product);
			return;
		case 9:
			MontgomeryMult9(factor1, factor2, product);
			return;
		case 10:
			MontgomeryMult10(factor1, factor2, product);
			return;
		case 11:
			MontgomeryMult11(factor1, factor2, product);
			return;
		}

		/* drop through to here only if NumberLength is 12 */
		memset(Prod, 0, NumberLength * sizeof(limb));
		for (i = 0; i < NumberLength; i++) {
			Pr = (Nbr = *(factor1 + i)) * (int64_t)*factor2 + (uint32_t)Prod[0];
			MontDig = ((int32_t)Pr * MontgomeryMultN[0]) & MAX_VALUE_LIMB;
			Prod[0] = (Pr = (((int64_t)MontDig * TestNbr[0] + Pr) >> BITS_PER_GROUP) +
				(int64_t)MontDig * TestNbr[1] + (int64_t)Nbr * factor2[1] + (uint32_t)Prod[1])& MAX_VALUE_LIMB;
			for (j = 2; j < NumberLength; j++) {
				Prod[j - 1] = ((Pr = (Pr >> BITS_PER_GROUP) +
					(int64_t)MontDig * TestNbr[j] + (int64_t)Nbr * factor2[j] + (uint32_t)Prod[j]) & MAX_VALUE_LIMB);
			}
			Prod[j - 1] = (int32_t)(Pr >> BITS_PER_GROUP);
		}
#else
		double dLimbRange = (double)LIMB_RANGE;
		double dInvLimbRange = (double)1 / dLimbRange;
		memset(Prod, 0, NumberLength * sizeof(limb));
		for (i = 0; i < NumberLength; i++)
		{
			int Nbr = *(factor1 + i);
			double dNbr = (double)Nbr;
			int low = Nbr * *factor2 + Prod[0];
			double dAccum = dNbr * (double)*factor2 + (double)Prod[0];
			int MontDig = (low * MontgomeryMultN[0]) & MAX_VALUE_LIMB;
			double dMontDig = (double)MontDig;
			dAccum += dMontDig * (double)TestNbr[0];
			// At this moment dAccum is multiple of LIMB_RANGE.
			dAccum = floor(dAccum*dInvLimbRange + 0.5);
			low = ((unsigned int)dAccum + MontDig * TestNbr[1] +
				Nbr * *(factor2 + 1) + Prod[1]) & MAX_VALUE_LIMB;
			dAccum += dMontDig * TestNbr[1] + dNbr * *(factor2 + 1) + (unsigned int)Prod[1];
			Prod[0] = low;
			for (j = 2; j < NumberLength; j++)
			{
				// Subtract or add 0x20000000 so the multiplication by dVal is not nearly an integer.
				// In that case, there would be an error of +/- 1.
				if (low < HALF_INT_RANGE)
				{
					dAccum = ((dAccum + HALF_INT_RANGE / 2)*dInvLimbRange);
				}
				else
				{
					dAccum = ((dAccum - HALF_INT_RANGE / 2)*dInvLimbRange);
				}
				low = (int)(dAccum - floor(dAccum * dInvLimbRange) * dLimbRange);
				dAccum += dMontDig * TestNbr[j] + dNbr * *(factor2 + j) + (unsigned int)Prod[j];
				low = (low + MontDig * TestNbr[j] +
					Nbr * *(factor2 + j) + Prod[j]) & MAX_VALUE_LIMB;
				Prod[j - 1] = low;
			}
			if (low < HALF_INT_RANGE)
			{
				dAccum = ((dAccum + HALF_INT_RANGE / 2)*dInvLimbRange);
			}
			else
			{
				dAccum = ((dAccum - HALF_INT_RANGE / 2)*dInvLimbRange);
			}
			Prod[j - 1] = (unsigned int)dAccum;  // Most significant limb can be greater than LIMB_RANGE
		}
#endif
		for (j = NumberLength - 1; j >= 0; j--) {
			if (Prod[j] != TestNbr[j]) {
				break;
			}
		}
		if (j<0 || (unsigned int)Prod[j] >= (unsigned int)TestNbr[j])
		{        // Prod >= TestNbr, so perform Prod <- Prod - TestNbr
			carry = 0;
			for (count = 0; count < NumberLength; count++) {
				carry += Prod[count] - TestNbr[count];
				Prod[count] = carry & MAX_VALUE_LIMB;
				carry >>= BITS_PER_GROUP;
			}
		}
		memcpy(product, Prod, NumberLength * sizeof(limb));
		return;
	}

	// NumberLength > 12; 
	// Compute T
	multiply(factor1, factor2, product, NumberLength, NULL);
	// Compute m
	multiply(product, MontgomeryMultN, aux, NumberLength, NULL);
	// Compute mN
	multiply(TestNbr, aux, aux2, NumberLength, NULL);
	// Check if lowest half of mN is not zero
	for (count = NumberLength - 1; count >= 0; count--) {
		if (aux2[count] != 0)
		{    // Lowest half of mN is not zero
			break;
		}
	}
	// If lowest half of mN is zero, compute hi(T) + hi(mN)
	// else compute hi(T) + hi(mN) + 1
	// Where hi(number) is the high half of number.
	cy = (count >= 0 ? LIMB_RANGE : 0);
	index = NumberLength;
	for (count = 0; count < NumberLength; count++) {
		cy = (cy >> BITS_PER_GROUP) +
			(unsigned int)*(product + index) + (unsigned int)aux2[index];
		*(product + count) = (int)(cy & MAX_VALUE_LIMB);
		index++;
	}
	// Check whether this number is greater than TestNbr.
	if (cy < LIMB_RANGE)
	{
		for (count = NumberLength - 1; count > 0; count--) 	{
			if (*(product + count) != TestNbr[count]) {
				break;
			}
		}
	}

	if (cy >= LIMB_RANGE || *(product + count) >= TestNbr[count])
	{  // The number is greater or equal than Testnbr. Subtract it.
		int borrow = 0;
		for (count = 0; count < NumberLength; count++) {
			borrow = (borrow >> BITS_PER_GROUP) +
				*(product + count) - TestNbr[count];
			*(product + count) = (int)(borrow & MAX_VALUE_LIMB);
		}
	}
	return;
}

/* Multiply big number by integer.
result = FactorBig* factorInt (mod TestNbr)
note: factorBig is modified by adding a leading zero */
static void modmultIntExtended(const limb factorBig[], int factorInt, limb result[], 
	const limb TestNbr[], int nbrLen) {
#ifdef _USING64BITS_
	int64_t carry;
#else
	double dTrialQuotient, dAccumulator, dFactorInt;
	double dInvLimbRange = 1 / (double)LIMB_RANGE;
	int low;
#endif
	int i;
	int TrialQuotient;      // approximate value of FactorBig*FactorInt / TestNbr
	double dTestNbr, dFactorBig;

	if (nbrLen == 1) {
		smallmodmult(*factorBig, factorInt, result, TestNbr[0]);
		return;
	}
	*((limb *)factorBig + nbrLen) = 0;   // note: factorBig is modifed, but value does not change
	dTestNbr = getMantissa(&TestNbr[nbrLen], nbrLen);
	dFactorBig = getMantissa(&factorBig[nbrLen], nbrLen);
	TrialQuotient = (int)(unsigned int)floor(dFactorBig * (double)factorInt / dTestNbr + 0.5);
	if ((unsigned int)TrialQuotient >= LIMB_RANGE)
	{   // Maximum value for limb.
		TrialQuotient = MAX_VALUE_LIMB;
	}
	// Compute result <- factorBig * factorInt - TrialQuotient * TestNbr
	// If TrialQuotient is accurate this is the required value

#ifdef _USING64BITS_
	carry = 0;
	for (i = 0; i <= nbrLen; i++) {
		carry += (int64_t)factorBig[i] * factorInt -
			(int64_t)TrialQuotient * TestNbr[i];
		result[i] = (int)carry & MAX_INT_NBR;
		carry >>= BITS_PER_GROUP;

	}
#else
	dFactorInt = (double)factorInt;
	dTrialQuotient = (double)TrialQuotient;
	low = 0;
	dAccumulator = 0;
	for (i = 0; i <= nbrLen; i++)
	{
		dAccumulator += *ptrFactorBig * dFactorInt - dTrialQuotient * *ptrTestNbr;
		low += *ptrFactorBig * factorInt - TrialQuotient * *ptrTestNbr;
		low &= MAX_VALUE_LIMB;
		// Subtract or add 0x20000000 so the multiplication by dVal is not nearly an integer.
		// In that case, there would be an error of +/- 1.
		*(result + i) = low;
		if (low < HALF_INT_RANGE)
		{
			dAccumulator = floor(dAccumulator*dInvLimbRange + 0.25);
		}
		else
		{
			dAccumulator = floor(dAccumulator*dInvLimbRange - 0.25);
		}
		low = (int)dAccumulator & MAX_VALUE_LIMB;
		ptrFactorBig++;
		ptrTestNbr++;
	}
#endif
	while ((result[nbrLen] & MAX_VALUE_LIMB) != 0) {
		unsigned int cy = 0;

		for (i = 0; i <= nbrLen; i++) { 
			/* result += TestNbr */
			cy += (unsigned int)TestNbr[i] + (unsigned int)result[i];
			result[i] = (int)(cy & MAX_VALUE_LIMB);
			cy >>= BITS_PER_GROUP;
		}
	}
}

/* result = FactorBig* factorInt (mod TestNbr) 
note: result & factorBig may be the same variable */
void modmultInt(const limb factorBig[], int factorInt, limb result[]) {
	modmultIntExtended(factorBig, factorInt, result, TestNbr, NumberLength);
}


// Input: base = base in Montgomery notation.
//        exp  = exponent.
//        nbrGroupsExp = number of limbs of exponent.
// Output: power = power in Montgomery notation (mod TestNbr).
//void modPow(const limb *base, const limb *exp, int nbrGroupsExp, limb *power) {
//	int mask, index;
//	memcpy(power, MontgomeryMultR1, (NumberLength + 1) * sizeof(*power));  // power <- 1
//	for (index = nbrGroupsExp - 1; index >= 0; index--) {
//		int groupExp = (int)*(exp + index);
//		for (mask = 1 << (BITS_PER_GROUP - 1); mask > 0; mask >>= 1) {
//			modmult(power, power, power);
//			if ((groupExp & mask) != 0) {
//				modmult(power, base, power);
//			}
//		}
//	}
//}

// Output: power = power in Montgomery notation (mod TestNbr).
//void modPowBaseInt(int base, const limb *exp, int nbrGroupsExp, limb *power) {
//	int mask, index;
//	memcpy(power, MontgomeryMultR1, (NumberLength + 1) * sizeof(limb));  // power <- 1
//	for (index = nbrGroupsExp - 1; index >= 0; index--) {
//		int groupExp = (int)*(exp + index);
//		for (mask = 1 << (BITS_PER_GROUP - 1); mask > 0; mask >>= 1) {
//			modmult(power, power, power);
//			if ((groupExp & mask) != 0) {
//				modmultInt(power, base, power);
//			}
//		}
//	}
//}

/* U' <- eU + fV, V' <- gU + hV                                        */
/* U <- U', V <- V'                                                    */
static void AddMult(limb *firstBig, int e, int f, limb *secondBig, int g, int h, int nbrLen)
{
#ifdef _USING64BITS_
	int64_t carryU = 0;
	int64_t carryV = 0;
	int ctr;
	for (ctr = 0; ctr <= nbrLen; ctr++)
	{
		int u = *firstBig;
		int v = *secondBig;
		carryU += u*(int64_t)e + v*(int64_t)f;
		carryV += u*(int64_t)g + v*(int64_t)h;
		*(firstBig++) = (int)(carryU & MAX_INT_NBR);
		*(secondBig++) = (int)(carryV & MAX_INT_NBR);
		carryU >>= BITS_PER_GROUP;
		carryV >>= BITS_PER_GROUP;
	}
#else
	double dVal = 1 / (double)LIMB_RANGE;
	int ctr, carryU, carryV;
	double dFactorE = (double)e;
	double dFactorF = (double)f;
	double dFactorG = (double)g;
	double dFactorH = (double)h;
	carryU = carryV = 0;
	for (ctr = 0; ctr <= nbrLen; ctr++)
	{
		int u = *firstBig;
		int v = *secondBig;
		int lowU = (carryU + u * e + v * f) & MAX_INT_NBR;
		int lowV = (carryV + u * g + v * h) & MAX_INT_NBR;
		// Subtract or add 0.25 so the multiplication by dVal is not nearly an integer.
		// In that case, there would be an error of +/- 1.
		double dCarry = ((double)carryU + (double)u * dFactorE +
			(double)v * dFactorF)*dVal;
		if (lowU < HALF_INT_RANGE)
		{
			carryU = (int)floor(dCarry + 0.25);
		}
		else
		{
			carryU = (int)floor(dCarry - 0.25);
		}
		dCarry = ((double)carryV + (double)u * dFactorG +
			(double)v * dFactorH)*dVal;
		if (lowV < HALF_INT_RANGE)
		{
			carryV = (int)floor(dCarry + 0.25);
		}
		else
		{
			carryV = (int)floor(dCarry - 0.25);
		}
		*(firstBig++) = lowU;
		*(secondBig++) = lowV;
	}
#endif
}

// Perform first <- (first - second) / 2
// first must be greater than second.
static int HalveDifference(limb *first, limb *second, int len)
{
	int i;
	int borrow, prevLimb;
	// Perform first <- (first - second)/2.
	borrow = *first - *second;
	prevLimb = borrow & MAX_VALUE_LIMB;
	borrow >>= BITS_PER_GROUP;
	for (i = 1; i < len; i++)
	{
		int currLimb;
		borrow += *(first + i) - *(second + i);
		currLimb = borrow & MAX_VALUE_LIMB;
		borrow >>= BITS_PER_GROUP;
		*(first + i - 1) = ((prevLimb >> 1) |
			(currLimb << (BITS_PER_GROUP - 1))) & MAX_VALUE_LIMB;
		prevLimb = currLimb;
	}
	*(first + i - 1) = prevLimb >> 1;
	// Get length of result.
	for (len--; len > 0; len--)
	{
		if (*(first + len) != 0)
		{
			break;
		}
	}
	return len + 1;
}

static int modInv(int NbrMod, int currentPrime)
{
	int QQ, T1, T3;
	int V1 = 1;
	int V3 = NbrMod;
	int U1 = 0;
	int U3 = currentPrime;
	while (V3 != 0) 	{
		if (U3 < V3 + V3) {               
			// QQ = 1
			T1 = U1 - V1;
			T3 = U3 - V3;
		}
		else 	{
			QQ = U3 / V3;
			T1 = U1 - V1 * QQ;
			T3 = U3 - V3 * QQ;
		}
		U1 = V1;
		U3 = V3;
		V1 = T1;
		V3 = T3;
	}

	return U1 + (currentPrime & (U1 >> 31));
}

/***********************************************************************/
/* NAME: ModInvBigNbr                                                  */
/*                                                                     */
/* PURPOSE: Find the inverse multiplicative modulo M.                  */
/* The algorithm terminates with inv = num^(-1) mod M.                 */
/*                                                                     */
/* This routine uses Kaliski Montgomery inverse algorithm              */
/* with changes by E. Savas and C. K. Koc.                             */
/*  1. U <- M, V <- X, R <- 0, S <- 1, k <- 0                          */
/*  2. while V > 0 do                                                  */
/*  3.   if U even then U <- U / 2, S <- 2S                            */
/*  4.   elsif V even then V <- V / 2, R <- 2R                         */
/*  5.   elsif U > V  then U <- (U - V) / 2, R <- R + S, S <- 2S       */
/*  6.   else V <- (V - U) / 2, S <- S + R, R <- 2R                    */
/*  7.   k <- k + 1                                                    */
/*  8. if R >= M then R <- R - M                                       */
/*  9. R <- M - R                                                      */
/* 10. R <- MonPro(R, R2)                                              */
/* 11. return MonPro(R, 2^(m-k))                                       */
/*                                                                     */
/*  In order to reduce the calculations, several single precision      */
/*  variables are added:                                               */
/*                                                                     */
/* R' <- aR + bS, S' <-  cR + dS                                       */
/* U' <- aU - bV, V' <- -cU + dV                                       */
/***********************************************************************/
// note: both num and mod are modified (but the value is not changed)
void ModInvBigNbr(const limb num[], limb inv[], const limb mod[], int nbrLen)
{

	int len;
	int k, steps;
	int a, b, c, d;  // Coefficients used to update variables R, S, U, V.
	int size, i;
	int bitCount;
	int lenRS;
	int lenU, lenV;
	int lowU, lowV;
	double highU, highV;
	int borrow;
	if (nbrLen == 1) {
		*inv = modInv(*num, *mod);
		return;
	}
	if (powerOf2Exponent != 0) {    // TestNbr is a power of 2.
		ComputeInversePower2(num, inv, aux);
		*(inv + powerOf2Exponent / BITS_PER_GROUP) &= (1 << (powerOf2Exponent % BITS_PER_GROUP)) - 1;
		return;
	}
	//  1. U <- M, V <- X, R <- 0, S <- 1, k <- 0
	size = (nbrLen + 1) * sizeof(limb);
	*((limb *)mod + nbrLen) = 0;   // value of mod is not changed
	*((limb *)num + nbrLen) = 0;   // value of num is not changed
	memcpy(U, mod, size);
	memcpy(V, num, size);
	// Maximum value of R and S can be up to 2*M, so one more limb is needed.
	memset(R, 0, size);   // R <- 0
	memset(S, 0, size);   // S <- 1
	S[0] = 1;
	lenRS = 1;
	k = steps = 0;
	// R' <- aR + bS, S' <- cR + dS
	a = d = 1;  // R' = R, S' = S.
	b = c = 0;
	len = nbrLen;
	// Find length of U.
	for (lenU = nbrLen - 1; lenU > 0; lenU--) {
		if (U[lenU] != 0)
		{
			break;
		}
	}
	lenU++;
	// Find length of V.
	for (lenV = nbrLen - 1; lenV > 0; lenV--) {
		if (V[lenV] != 0) {
			break;
		}
	}
	lenV++;
	lowU = U[0];
	lowV = V[0];
	// Initialize highU and highV.
	if (lenU > 1 || lenV > 1) {
		if (lenV >= lenU) {
			highV = (double)V[lenV - 1] * (double)LIMB_RANGE + (double)V[lenV - 2];
			if (lenV == lenU) {
				highU = (double)U[lenV - 1] * (double)LIMB_RANGE + (double)U[lenV - 2];
			}
			else if (lenV == lenU + 1) {
				highU = (double)U[lenV - 2];
			}
			else {
				highU = 0;
			}
		}
		else {
			highU = (double)U[lenU - 1] * (double)LIMB_RANGE + (double)U[lenU - 2];
			if (lenU == lenV + 1) {
				highV = (double)V[lenU - 2];
			}
			else {
				highV = 0;
			}
		}
		//  2. while V > 0 do
		for (;;) {
			//  3.   if U even then U <- U / 2, S <- 2S
			if ((lowU & 1) == 0) {     // U is even.
				lowU >>= 1;
				highV += highV;
				// R' <- aR + bS, S' <- cR + dS
				c *= 2; d *= 2;  // Multiply S by 2.
			}
			//  4.   elsif V even then V <- V / 2, R <- 2R
			else if ((lowV & 1) == 0) {    // V is even.
				lowV >>= 1;
				highU += highU;
				// R' <- aR + bS, S' <- cR + dS
				a *= 2; b *= 2;  // Multiply R by 2.
			}
			else {
				//  5.   elsif U >= V  then U <- (U - V) / 2, R <- R + S, S <- 2S
				if (highU > highV) {     // U > V. Perform U <- (U - V) / 2
					lowU = (lowU - lowV) >> 1;
					highU -= highV;
					highV += highV;
					// R' <- aR + bS, S' <- cR + dS
					a += c; b += d;  // R <- R + S
					c *= 2; d *= 2;  // S <- 2S
				}
				//  6.   elsif V >= U then V <- (V - U) / 2, S <- S + R, R <- 2R
				else {    // V >= U. Perform V <- (V - U) / 2
					lowV = (lowV - lowU) >> 1;
					highV -= highU;
					highU += highU;
					// R' <- aR + bS, S' <- cR + dS
					c += a; d += b;  // S <- S + R
					a *= 2; b *= 2;  // R <- 2R
				}
			}
			//  7.   k <- k + 1
			// Adjust variables.
			if (++steps == BITS_PER_GROUP - 1)
			{  // compute now U and V and reset e, f, g and h.
			   // U' <- eU + fV, V' <- gU + hV
				len = (lenU > lenV ? lenU : lenV);
				memset(&U[lenU], 0, (len - lenU + 1) * sizeof(limb));
				memset(&V[lenV], 0, (len - lenV + 1) * sizeof(limb));
				memcpy(Ubak, U, (len + 1) * sizeof(limb));
				memcpy(Vbak, V, (len + 1) * sizeof(limb));
				AddMult(U, a, -b, V, -c, d, len);
				if ((U[lenU] | V[lenV]) & (1 << (BITS_PER_GROUP - 2)))
				{    // Complete expansion of U and V required for all steps.
					 //  2. while V > 0 do
					memcpy(U, Ubak, (len + 1) * sizeof(limb));
					memcpy(V, Vbak, (len + 1) * sizeof(limb));
					b = c = 0;  // U' = U, V' = V.
					a = d = 1;
					while (lenV > 1 || V[0] > 0) {
						//  3.   if U even then U <- U / 2, S <- 2S
						if ((U[0] & 1) == 0) {  // U is even.
							for (i = 0; i < lenU; i++)
							{  // Loop that divides U by 2.
								U[i] = ((U[i] >> 1) | (U[i + 1] << (BITS_PER_GROUP - 1))) & MAX_VALUE_LIMB;
							}
							if (U[lenU - 1] == 0) {
								lenU--;
							}
							// R' <- aR + bS, S' <- cR + dS
							c *= 2; d *= 2;  // Multiply S by 2.
						}
						//  4.   elsif V even then V <- V / 2, R <- 2R
						else if ((V[0] & 1) == 0) {    // V is even.
							for (i = 0; i < lenV; i++)
							{  // Loop that divides V by 2.
								V[i] = ((V[i] >> 1) | (V[i + 1] << (BITS_PER_GROUP - 1))) & MAX_VALUE_LIMB;
							}
							if (V[lenV - 1] == 0) {
								lenV--;
							}
							// R' <- aR + bS, S' <- cR + dS
							a *= 2; b *= 2;  // Multiply R by 2.
						}
						//  5.   elsif U >= V  then U <- (U - V) / 2, R <- R + S, S <- 2S
						else {
							len = (lenU > lenV ? lenU : lenV);
							for (i = len - 1; i > 0; i--) {
								if (U[i] != V[i]) {
									break;
								}
							}
							if (U[i] > V[i]) {     // U > V
								lenU = HalveDifference(U, V, len); // U <- (U - V) / 2
																   // R' <- aR + bS, S' <- cR + dS
								a += c; b += d;  // R <- R + S
								c *= 2; d *= 2;  // S <- 2S
							}
							//  6.   elsif V >= U then V <- (V - U) / 2, S <- S + R, R <- 2R
							else {    // V >= U
								lenV = HalveDifference(V, U, len); // V <- (V - U) / 2
																   // R' <- aR + bS, S' <- cR + dS
								c += a; d += b;  // S <- S + R
								a *= 2; b *= 2;  // R <- 2R
							}
						}
						//  7.   k <- k + 1
						if (++k % (BITS_PER_GROUP - 1) == 0) {
							break;
						}
					}
					if (lenV == 1 && V[0] == 0) {
						break;
					}
				}
				else {
					k += steps;
					for (i = 0; i < lenU; i++)
					{  // Loop that divides U by 2^(BITS_PER_GROUP - 1).
						U[i] = ((U[i] >> (BITS_PER_GROUP - 1)) | (U[i + 1] << 1)) & MAX_VALUE_LIMB;
					}
					U[lenU] = 0;
					while (lenU > 0 && U[lenU - 1] == 0) {
						lenU--;
					}
					for (i = 0; i < lenV; i++)
					{  // Loop that divides V by 2^(BITS_PER_GROUP - 1).
						V[i] = ((V[i] >> (BITS_PER_GROUP - 1)) | (V[i + 1] << 1)) & MAX_VALUE_LIMB;
					}
					V[lenV] = 0;
					while (lenV > 0 && V[lenV - 1] == 0) {
						lenV--;
					}
				}
				steps = 0;
				AddMult(R, a, b, S, c, d, lenRS);
				if (R[lenRS] != 0 || S[lenRS] != 0) {
					lenRS++;
				}
				lowU = U[0];
				lowV = V[0];
				b = c = 0;  // U' = U, V' = V.
				a = d = 1;
				if (lenU == 0 || lenV == 0 || (lenV == 1 && lenU == 1)) {
					break;
				}
				if (lenV >= lenU) {
					highV = (double)V[lenV - 1] * (double)LIMB_RANGE + (double)V[lenV - 2];
					if (lenV == lenU) {
						highU = (double)U[lenV - 1] * (double)LIMB_RANGE + (double)U[lenV - 2];
					}
					else if (lenV == lenU + 1) {
						highU = (double)U[lenV - 2];
						}
						else {
							highU = 0;
						}
				}
				else {
					highU = (double)U[lenU - 1] * (double)LIMB_RANGE + (double)U[lenU - 2];
					if (lenU == lenV + 1) {
						highV = (double)V[lenU - 2];
					}
					else {
						highV = 0;
					}
				}
			}
		}
	}

	if (lenU > 0) {
		//  2. while V > 0 do
		while (lowV > 0) {
			//  3.   if U even then U <- U / 2, S <- 2S
			if ((lowU & 1) == 0) {     // U is even.
				lowU >>= 1;
				// R' <- aR + bS, S' <- cR + dS
				c *= 2; d *= 2;  // Multiply S by 2.
			}
			//  4.   elsif V even then V <- V / 2, R <- 2R
			else if ((lowV & 1) == 0) {    // V is even.
				lowV >>= 1;
				// R' <- aR + bS, S' <- cR + dS
				a *= 2; b *= 2;  // Multiply R by 2.
				}
			//  5.   elsif U >= V  then U <- (U - V) / 2, R <- R + S, S <- 2S
			else if (lowU > lowV) {     // U > V. Perform U <- (U - V) / 2
				lowU = (lowU - lowV) >> 1;
				// R' <- aR + bS, S' <- cR + dS
				a += c; b += d;  // R <- R + S
				c *= 2; d *= 2;  // S <- 2S
				}
			//  6.   elsif V >= U then V <- (V - U) / 2, S <- S + R, R <- 2R
			else {    // V >= U. Perfrom V <- (V - U) / 2
				lowV = (lowV - lowU) >> 1;
				// R' <- aR + bS, S' <- cR + dS
				c += a; d += b;  // S <- S + R
				a *= 2; b *= 2;  // R <- 2R
			}

			//  7.   k <- k + 1
			if (++steps == BITS_PER_GROUP - 1) {  
			   // compute now R and S and reset a, b, c and d.
			   // R' <- aR + bS, S' <- cR + dS
			    AddMult(R, a, b, S, c, d, nbrLen + 1);
				b = c = 0;  // R' = R, S' = S.
				a = d = 1;
				k += steps;
				steps = 0;
			}
		}
	}

	AddMult(R, a, b, S, c, d, nbrLen + 1);
	k += steps;
	//  8. if R >= M then R <- R - M
	for (i = nbrLen; i > 0; i--) {
		if (R[i] != *(mod + i)) {
			break;
		}
	}

	if ((unsigned int)R[i] >= (unsigned int)*(mod + i)) {      // R >= M.
		borrow = 0;
		for (i = 0; i <= nbrLen; i++) {
			borrow += R[i] - *(mod + i);
			R[i] = borrow & MAX_VALUE_LIMB;
			borrow >>= BITS_PER_GROUP;
		}
	}

	//  9. R <- M - R
	borrow = 0;
	for (i = 0; i <= nbrLen; i++) {
		borrow += *(mod + i) - R[i];
		R[i] = borrow & MAX_VALUE_LIMB;
		borrow >>= BITS_PER_GROUP;
	}

	R[nbrLen] = 0;
	// At this moment R = x^(-1)*2^k
	// 10. R <- MonPro(R, R2)
	modmult(R, MontgomeryMultR2, R);
	R[nbrLen] = 0;
	// At this moment R = x^(-1)*2^(k+m)
	// 11. return MonPro(R, 2^(m-k))
	memset(S, 0, size);
	bitCount = nbrLen*BITS_PER_GROUP - k;
	if (bitCount < 0) {
		bitCount += nbrLen*BITS_PER_GROUP;
		S[bitCount / BITS_PER_GROUP] = 1 << (bitCount % BITS_PER_GROUP);
		modmult(R, S, inv);
	}
	else {
		S[bitCount / BITS_PER_GROUP] = 1 << (bitCount % BITS_PER_GROUP);
		modmult(R, S, inv);
		modmult(inv, MontgomeryMultR2, inv);
	}

}

// Input: base = base in Montgomery notation.
//        exp  = exponent.
//        nbrGroupsExp = number of limbs of exponent.
// Output: power = power in Montgomery notation.
void modPow(const limb* base, const limb* exp, int nbrGroupsExp, limb* power)
{
	int lenBytes = (NumberLength + 1) * (int)sizeof(*power);
	(void)memcpy(power, MontgomeryMultR1, lenBytes);  // power <- 1
	for (int index = nbrGroupsExp - 1; index >= 0; index--)
	{
		int groupExp = *(exp + index);
		for (unsigned int mask = HALF_INT_RANGE_U; mask > 0U; mask >>= 1)
		{
			modmult(power, power, power);
			if (((unsigned int)groupExp & mask) != 0U)
			{
				modmult(power, base, power);
			}
		}
	}
}

void modPowBaseInt(int base, const limb* exp, int nbrGroupsExp, limb* power)
{
	int NumberLengthBytes = (NumberLength + 1) * (int)sizeof(limb);
	(void)memcpy(power, MontgomeryMultR1, NumberLengthBytes);  // power <- 1
	for (int index = nbrGroupsExp - 1; index >= 0; index--)
	{
		int groupExp = *(exp + index);
		for (unsigned int mask = HALF_INT_RANGE_U; mask > 0U; mask >>= 1)
		{
			modmult(power, power, power);
			if (((unsigned int)groupExp & mask) != 0U)
			{
				modmultInt(power, base, power);
			}
		}
	}
}

// Compute power = base^exponent (mod modulus)
// Assumes GetMontgomeryParms routine for modulus already called.
// This works only for odd moduli.
void BigIntModularPower(const BigInteger* base, const BigInteger* exponent, BigInteger* power)
{
	int lenBytes;
	CompressLimbsBigInteger(aux5, base);
	modmult(aux5, MontgomeryMultR2, aux6);   // Convert base to Montgomery notation.
	modPow(aux6, exponent->limbs, exponent->nbrLimbs, aux5);
	lenBytes = NumberLength * (int)sizeof(limb);
	(void)memset(aux4, 0, lenBytes); // Convert power to standard notation.
	aux4[0] = 1;
	modmult(aux4, aux5, aux6);
	UncompressLimbsBigInteger(aux6, power);
}

// Compute modular division for odd moduli.
void BigIntModularDivision(const BigInteger* Num, const BigInteger* Den,
	const BigInteger* mod, BigInteger* quotient)
{
	NumberLength = mod->nbrLimbs;
	// Reduce Num modulo mod.
	tmpNum = *Num % *mod; // (void)BigIntRemainder(Num, mod, &tmpNum);
	if (tmpNum.sign == SIGN_NEGATIVE)
	{
		tmpNum += *mod;   // BigIntAdd(&tmpNum, mod, &tmpNum);
	}
	// Reduce Den modulo mod.
	tmpDen = *Den % *mod;  // (void)BigIntRemainder(Den, mod, &tmpDen);
	if (tmpDen.sign == SIGN_NEGATIVE)
	{
		tmpDen += *mod; // BigIntAdd(&tmpDen, mod, &tmpDen);
	}
	CompressLimbsBigInteger(aux3, &tmpDen);
	modmult(aux3, MontgomeryMultR2, aux3);      // aux3 <- Den in Montgomery notation
												// tmpDen.limbs <- 1 / Den in Montg notation.
	(void)ModInvBigNbr(aux3, tmpDen.limbs, TestNbr, NumberLength);
	CompressLimbsBigInteger(aux4, &tmpNum);
	modmult(tmpDen.limbs, aux4, aux3);          // aux3 <- Num / Den in standard notation.
	UncompressLimbsBigInteger(aux3, quotient);  // Get Num/Den
}