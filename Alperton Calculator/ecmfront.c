/*
This file is part of Alpertron Calculators.
Copyright 2016 Dario Alejandro Alpern
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

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include "bignbr.h"
#include "expression.h"
#include "factor.h"
#include "showtime.h"
#include "batch.h"
#ifdef __EMSCRIPTEN__
extern long long lModularMult;
#endif
extern BigInteger tofactor;
static BigInteger Quad1, Quad2, Quad3, Quad4;
extern BigInteger factorValue;
static BigInteger result;
extern int valuesProcessed;
//extern char *ptrOutput;
extern int NextEC;
static void ComputeFourSquares(struct sFactors *pstFactors);
static void GetEulerTotient(char **pptrOutput);
static void GetMobius(char **pptrOutput);
static void GetNumberOfDivisors(char **pptrOutput);
static void GetSumOfDivisors(char **pptrOutput);
static void ShowFourSquares(char **pptrOutput);
static int doFactorization;
static char *knownFactors;

#ifdef FACTORIZATION_APP
void batchCallback(char **pptrOutput)
{
	char *ptrFactorDec = tofactorDec;
	NumberLength = tofactor.nbrLimbs;
	if (tofactor.sign == SIGN_NEGATIVE)
	{
		*ptrFactorDec++ = '-';
	}
	CompressBigInteger(nbrToFactor, &tofactor);
	if (hexadecimal)
	{
		Bin2Hex(tofactor.limbs, ptrFactorDec, tofactor.nbrLimbs, groupLen);
	}
	else
	{
		Bin2Dec(tofactor.limbs, ptrFactorDec, tofactor.nbrLimbs, groupLen);
	}
	if (doFactorization)
	{
		factor(&tofactor, nbrToFactor, factorsMod, astFactorsMod, NULL);
	}
	SendFactorizationToOutput(astFactorsMod, pptrOutput, doFactorization);
}
#endif

static void ExponentToBigInteger(int exponent, BigInteger *bigint)
{
	if (exponent > MAX_VALUE_LIMB)
	{
		bigint->limbs[0].x = exponent - MAX_VALUE_LIMB;
		bigint->limbs[1].x = 1;
		bigint->nbrLimbs = 2;
	}
	else
	{
		bigint->limbs[0].x = exponent + 1;
		bigint->nbrLimbs = 1;
	}
	bigint->sign = SIGN_POSITIVE;
}

// Find number of divisors as the product of all exponents plus 1.
static void GetNumberOfDivisors(char **pptrOutput)
{
	char *ptrOutput = *pptrOutput;
	struct sFactors *pstFactor;
	int factorNumber;
	result.limbs[0].x = 1;    // Set result to 1.
	result.nbrLimbs = 1;
	result.sign = SIGN_POSITIVE;
	pstFactor = &astFactorsMod[1];
	for (factorNumber = 1; factorNumber <= astFactorsMod[0].multiplicity; factorNumber++)
	{
		int *ptrFactor = pstFactor->ptrFactor;
		if (*ptrFactor == 1 && *(ptrFactor + 1) < 2)
		{                        // Factor is 1.
			break;
		}
		ExponentToBigInteger(pstFactor->multiplicity, &factorValue);
		BigIntMultiply(&factorValue, &result, &result);
		pstFactor++;
	}
	strcpy(ptrOutput, lang ? "\nCantidad de divisores: " : "\nNumber of divisors: ");
	ptrOutput += strlen(ptrOutput);
	if (hexadecimal)
	{
		BigInteger2Hex(&result, ptrOutput, groupLen);
	}
	else
	{
		BigInteger2Dec(&result, ptrOutput, groupLen);
	}
	ptrOutput += strlen(ptrOutput);
	strcpy(ptrOutput, "\n");
	ptrOutput += strlen(ptrOutput);
	*pptrOutput = ptrOutput;
}

static void GetSumOfDivisors(char **pptrOutput)
{
	char *ptrOutput = *pptrOutput;
	SumOfDivisors(&result);
	strcpy(ptrOutput, lang ? "\nSuma de divisores: " : "\nSum of divisors: ");
	ptrOutput += strlen(ptrOutput);
	if (hexadecimal)
	{
		BigInteger2Hex(&result, ptrOutput, groupLen);
	}
	else
	{
		BigInteger2Dec(&result, ptrOutput, groupLen);
	}
	ptrOutput += strlen(ptrOutput);
	strcpy(ptrOutput, "\n");
	ptrOutput += strlen(ptrOutput);
	*pptrOutput = ptrOutput;
}

static void GetEulerTotient(char **pptrOutput)
{
	char *ptrOutput = *pptrOutput;
	Totient(&result);
	strcpy(ptrOutput, lang ? "\nPhi de Euler: " : "\nEuler's totient: ");
	ptrOutput += strlen(ptrOutput);
	if (hexadecimal)
	{
		BigInteger2Hex(&result, ptrOutput, groupLen);
	}
	else
	{
		BigInteger2Dec(&result, ptrOutput, groupLen);
	}
	ptrOutput += strlen(ptrOutput);
	strcpy(ptrOutput, "\n");
	ptrOutput += strlen(ptrOutput);
	*pptrOutput = ptrOutput;
}

// Find Mobius as zero if some exponent is > 1, 1 if the number of factors is even, 
// -1 if it is odd.
static void GetMobius(char **pptrOutput)
{
	char *ptrOutput = *pptrOutput;
	struct sFactors *pstFactor;
	int mobius = 1;
	pstFactor = &astFactorsMod[1];
	if (astFactorsMod[0].multiplicity > 1 || *pstFactor->ptrFactor != 1 ||
		*(pstFactor->ptrFactor + 1) != 1)
	{                                // Number to factor is not 1.
		int factorNumber;
		for (factorNumber = 1; factorNumber <= astFactorsMod[0].multiplicity; factorNumber++)
		{
			if (pstFactor->multiplicity == 1)
			{
				mobius = -mobius;
			}
			else
			{
				mobius = 0;
			}
			pstFactor++;
		}
	}
	strcpy(ptrOutput, "\nM�bius: ");
	ptrOutput += strlen(ptrOutput);
	if (mobius < 0)
	{
		mobius = -mobius;
		*ptrOutput++ = '-';
	}
	int2dec(&ptrOutput, mobius);
	strcpy(ptrOutput, "\n");
	ptrOutput += strlen(ptrOutput);
	*pptrOutput = ptrOutput;
}

static void ComputeFourSquares(struct sFactors *pstFactors)
{
	int indexPrimes;
	BigInteger p, q, K, Mult1, Mult2, Mult3, Mult4;
	BigInteger Tmp, Tmp1, Tmp2, Tmp3, Tmp4, M1, M2, M3, M4;
	struct sFactors *pstFactor;
	static limb minusOneMont[MAX_LEN];

	intToBigInteger(&Quad1, 1);     // 1 = 1^2 + 0^2 + 0^2 + 0^2
	intToBigInteger(&Quad2, 0);
	intToBigInteger(&Quad3, 0);
	intToBigInteger(&Quad4, 0);
	int Factorix = 1;      // Point to first factor in array of factors.
	if (pstFactors[0].multiplicity == 1 && pstFactors[1].multiplicity == 1)
	{
		if (*((pstFactors[1].ptrFactor)+1) == 1)
		{                              // Number to factor is 1.
			return;
		}
		if (*((pstFactors[1].ptrFactor)+1) == 0)
		{                             // Number to factor is 0.
			intToBigInteger(&Quad1, 0);     // 0 = 0^2 + 0^2 + 0^2 + 0^2
			return;
		}
	}
	for (Factorix = 1;  Factorix <= pstFactors[0].multiplicity;  Factorix++)
	{
		if (pstFactors[Factorix].multiplicity % 2 == 0)
		{                              // Prime factor appears twice.
			continue;
		}
		NumberLength = *pstFactors[Factorix].ptrFactor;
		UncompressBigInteger(pstFactors[Factorix].ptrFactor, &p);
		p.sign = SIGN_POSITIVE;
		CopyBigInt(&q, &p);
		addbigint(&q, -1);             // q <- p-1
		if (p.nbrLimbs == 1 && p.limbs[0].x == 2)
		{   /* Prime factor is 2 */
			intToBigInteger(&Mult1, 1); // 2 = 1^2 + 1^2 + 0^2 + 0^2
			intToBigInteger(&Mult2, 1);
			intToBigInteger(&Mult3, 0);
			intToBigInteger(&Mult4, 0);
		}
		else
		{ /* Prime factor is not 2 */
			NumberLength = p.nbrLimbs;
			memcpy(&TestNbr, p.limbs, NumberLength * sizeof(limb));
			TestNbr[NumberLength].x = 0;
			GetMontgomeryParms(NumberLength);
			memset(minusOneMont, 0, NumberLength * sizeof(limb));
			SubtBigNbrModN(minusOneMont, MontgomeryMultR1, minusOneMont, TestNbr, NumberLength);
			memset(K.limbs, 0, NumberLength * sizeof(limb));
			if ((p.limbs[0].x & 3) == 1)
			{ /* if p = 1 (mod 4) */
				CopyBigInt(&q, &p);
				subtractdivide(&q, 1, 4);     // q = (prime-1)/4
				K.limbs[0].x = 1;
				do
				{    // Loop that finds mult1 = sqrt(-1) mod prime in Montgomery notation.
					K.limbs[0].x++;
					modPow(K.limbs, q.limbs, q.nbrLimbs, Mult1.limbs);
				} while (!memcmp(Mult1.limbs, MontgomeryMultR1, NumberLength * sizeof(limb)) ||
					!memcmp(Mult1.limbs, minusOneMont, NumberLength * sizeof(limb)));
				Mult1.sign = SIGN_POSITIVE;
				memset(Mult2.limbs, 0, p.nbrLimbs * sizeof(limb));
				Mult2.limbs[0].x = 1;
				Mult2.nbrLimbs = 1;
				Mult2.sign = SIGN_POSITIVE;
				// Convert Mult1 to standard notation by multiplying by 1 in
				// Montgomery notation.
				modmult(Mult1.limbs, Mult2.limbs, Mult3.limbs);
				memcpy(Mult1.limbs, Mult3.limbs, p.nbrLimbs * sizeof(limb));
				for (Mult1.nbrLimbs = p.nbrLimbs; Mult1.nbrLimbs > 1; Mult1.nbrLimbs--)
				{  // Adjust number of limbs so the most significant limb is not zero.
					if (Mult1.limbs[Mult1.nbrLimbs - 1].x != 0)
					{
						break;
					}
				}
				for (;;)
				{
					BigIntMultiply(&Mult1, &Mult1, &Tmp);
					BigIntMultiply(&Mult2, &Mult2, &Tmp1);
					BigIntAdd(&Tmp, &Tmp1, &Tmp);
					BigIntDivide(&Tmp, &p, &K);        // K <- (mult1^2 + mult2^2) / p
					if (K.nbrLimbs == 1 && K.limbs[0].x == 1)
					{    // If K = 1...
						intToBigInteger(&Mult3, 0);
						intToBigInteger(&Mult4, 0);
						break;
					}
					BigIntRemainder(&Mult1, &K, &M1);  // M1 <- Mult1 % K
					if (M1.sign == SIGN_NEGATIVE)
					{
						BigIntAdd(&M1, &K, &M1);
					}
					BigIntRemainder(&Mult2, &K, &M2);  // M2 <- Mult2 % K
					if (M2.sign == SIGN_NEGATIVE)
					{
						BigIntAdd(&M2, &K, &M2);
					}
					CopyBigInt(&Tmp, &K);
					subtractdivide(&Tmp, -1, 2);       // Tmp <- (K+1) / 2
					BigIntSubt(&M1, &Tmp, &Tmp1);      // Tmp1 <- M1 - Tmp
					if (Tmp1.sign == SIGN_POSITIVE)    // If M1 >= K / 2 ... 
					{
						BigIntSubt(&M1, &K, &M1);        // M1 <- M1 - K
					}
					BigIntSubt(&M2, &Tmp, &Tmp1);      // Tmp1 <- M2 - Tmp
					if (Tmp1.sign == SIGN_POSITIVE)    // If M2 >= K / 2 ... 
					{
						BigIntSubt(&M2, &K, &M2);        // M1 <- M1 - K
					}
					BigIntMultiply(&Mult1, &M1, &Tmp);
					BigIntMultiply(&Mult2, &M2, &Tmp1);
					BigIntAdd(&Tmp, &Tmp1, &Tmp);
					BigIntDivide(&Tmp, &K, &Tmp2);     // Tmp2 <- (mult1*m1 + mult2*m2) / K
					BigIntMultiply(&Mult1, &M2, &Tmp);
					BigIntMultiply(&Mult2, &M1, &Tmp1);
					BigIntSubt(&Tmp, &Tmp1, &Tmp);
					BigIntDivide(&Tmp, &K, &Mult2);    // Mult2 <- (mult1*m2 - mult2*m1) / K
					CopyBigInt(&Mult1, &Tmp2);
				} /* end while */
			} /* end p = 1 (mod 4) */
			else
			{ /* if p = 3 (mod 4) */
				int mult1 = 0;
				CopyBigInt(&q, &p);
				subtractdivide(&q, 1, 2);     // q = (prime-1)/2
				memcpy(K.limbs, q.limbs, q.nbrLimbs * sizeof(limb));
				if (p.nbrLimbs > q.nbrLimbs)
				{
					K.limbs[q.nbrLimbs].x = 0;
				}
				// Compute Mult1 and Mult2 so Mult1^2 + Mult2^2 = -1 (mod p)
				memset(Mult1.limbs, 0, p.nbrLimbs * sizeof(limb));
				do
				{
					mult1++;
					// Increment Mult1 by 1 in Montgomery notation.
					AddBigNbrModN(Mult1.limbs, MontgomeryMultR1, Mult1.limbs, p.limbs, p.nbrLimbs);
					modmult(Mult1.limbs, Mult1.limbs, Tmp.limbs);
					SubtBigNbrModN(minusOneMont, Tmp.limbs, Tmp.limbs, p.limbs, p.nbrLimbs);
					modPow(Tmp.limbs, K.limbs, p.nbrLimbs, Tmp1.limbs);
					// At this moment Tmp1 = (-1 - Mult1^2)^((p-1)/2)
					// in Montgomery notation. Continue loop if it is not 1.
				} while (memcmp(Tmp1.limbs, MontgomeryMultR1, p.nbrLimbs));
				// After the loop finishes, Tmp = (-1 - Mult1^2) is a quadratic residue mod p.
				// Convert Mult1 to standard notation by multiplying by 1 in
				// Montgomery notation.
				intToBigInteger(&Mult1, mult1);
				CopyBigInt(&q, &p);
				subtractdivide(&q, -1, 4);  // q <- (p+1)/4.
											// Find Mult2 <- square root of Tmp = Tmp^q (mod p) in Montgomery notation.
				modPow(Tmp.limbs, q.limbs, p.nbrLimbs, Mult2.limbs);
				// Convert Mult2 from Montgomery notation to standard notation.
				memset(Tmp.limbs, 0, p.nbrLimbs * sizeof(limb));
				Tmp.limbs[0].x = 1;
				intToBigInteger(&Mult3, 1);
				intToBigInteger(&Mult4, 0);
				// Convert Mult2 to standard notation by multiplying by 1 in
				// Montgomery notation.
				modmult(Mult2.limbs, Tmp.limbs, Mult2.limbs);
				for (Mult2.nbrLimbs = p.nbrLimbs; Mult2.nbrLimbs > 1; Mult2.nbrLimbs--)
				{  // Adjust number of limbs so the most significant limb is not zero.
					if (Mult2.limbs[Mult2.nbrLimbs - 1].x != 0)
					{
						break;
					}
				}
				Mult2.sign = SIGN_POSITIVE;
				for (;;)
				{
					// Compute K <- (Mult1^2 + Mult2^2 + Mult3^2 + Mult4^2) / p
					BigIntMultiply(&Mult1, &Mult1, &Tmp);
					BigIntMultiply(&Mult2, &Mult2, &Tmp1);
					BigIntAdd(&Tmp, &Tmp1, &Tmp);
					BigIntMultiply(&Mult3, &Mult3, &Tmp1);
					BigIntAdd(&Tmp, &Tmp1, &Tmp);
					BigIntMultiply(&Mult4, &Mult4, &Tmp1);
					BigIntAdd(&Tmp, &Tmp1, &Tmp);
					BigIntDivide(&Tmp, &p, &K);
					if (K.nbrLimbs == 1 && K.limbs[0].x == 1)
					{   // K equals 1
						break;
					}
					if ((K.limbs[0].x & 1) == 0)
					{ // If K is even ...
						if ((Mult1.limbs[0].x + Mult2.limbs[0].x) & 1)
						{  // If Mult1 + Mult2 is odd...
							if (((Mult1.limbs[0].x + Mult3.limbs[0].x) & 1) == 0)
							{   // If Mult1 + Mult3 is even...
								CopyBigInt(&Tmp, &Mult2);
								CopyBigInt(&Mult2, &Mult3);
								CopyBigInt(&Mult3, &Tmp);
							}
							else
							{
								CopyBigInt(&Tmp, &Mult2);
								CopyBigInt(&Mult2, &Mult4);
								CopyBigInt(&Mult4, &Tmp);
							}
						} // At this moment Mult1+Mult2 = even, Mult3+Mult4 = even
						BigIntAdd(&Mult1, &Mult2, &Tmp1);
						subtractdivide(&Tmp1, 0, 2);  // Tmp1 <- (Mult1 + Mult2) / 2
						BigIntSubt(&Mult1, &Mult2, &Tmp2);
						subtractdivide(&Tmp2, 0, 2);  // Tmp2 <- (Mult1 - Mult2) / 2
						BigIntAdd(&Mult3, &Mult4, &Tmp3);
						subtractdivide(&Tmp3, 0, 2);  // Tmp3 <- (Mult3 + Mult4) / 2
						BigIntSubt(&Mult3, &Mult4, &Mult4);
						subtractdivide(&Mult4, 0, 2);  // Mult4 <- (Mult3 + Mult4) / 2
						CopyBigInt(&Mult3, &Tmp3);
						CopyBigInt(&Mult2, &Tmp2);
						CopyBigInt(&Mult1, &Tmp1);
						continue;
					} /* end if k is even */
					BigIntRemainder(&Mult1, &K, &M1);    // M1 <- Mult1 % K.
					if (M1.sign == SIGN_NEGATIVE)
					{
						BigIntAdd(&M1, &K, &M1);
					}
					BigIntRemainder(&Mult2, &K, &M2);    // M2 <- Mult2 % K.
					if (M2.sign == SIGN_NEGATIVE)
					{
						BigIntAdd(&M2, &K, &M2);
					}
					BigIntRemainder(&Mult3, &K, &M3);    // M3 <- Mult3 % K.
					if (M3.sign == SIGN_NEGATIVE)
					{
						BigIntAdd(&M3, &K, &M3);
					}
					BigIntRemainder(&Mult4, &K, &M4);    // M4 <- Mult4 % K.
					if (M4.sign == SIGN_NEGATIVE)
					{
						BigIntAdd(&M4, &K, &M4);
					}
					CopyBigInt(&Tmp, &K);
					subtractdivide(&Tmp, -1, 2);       // Tmp <- (K+1) / 2
					BigIntSubt(&M1, &Tmp, &Tmp1);      // Tmp1 <- M1 - Tmp
					if (Tmp1.sign == SIGN_POSITIVE)    // If M1 >= K / 2 ... 
					{
						BigIntSubt(&M1, &K, &M1);        // M1 <- M1 - K
					}
					BigIntSubt(&M2, &Tmp, &Tmp1);      // Tmp1 <- M2 - Tmp
					if (Tmp1.sign == SIGN_POSITIVE)    // If M2 >= K / 2 ... 
					{
						BigIntSubt(&M2, &K, &M2);        // M2 <- M2 - K
					}
					BigIntSubt(&M3, &Tmp, &Tmp1);      // Tmp1 <- M3 - Tmp
					if (Tmp1.sign == SIGN_POSITIVE)    // If M3 >= K / 2 ... 
					{
						BigIntSubt(&M3, &K, &M3);        // M3 <- M3 - K
					}
					BigIntSubt(&M4, &Tmp, &Tmp1);      // Tmp1 <- M4 - Tmp
					if (Tmp1.sign == SIGN_POSITIVE)    // If M4 >= K / 2 ... 
					{
						BigIntSubt(&M4, &K, &M4);        // M4 <- M4 - K
					}
					// Compute Tmp1 <- (Mult1*M1 + Mult2*M2 + Mult3*M3 + Mult4*M4) / K
					BigIntMultiply(&Mult1, &M1, &Tmp);
					BigIntMultiply(&Mult2, &M2, &Tmp4);
					BigIntAdd(&Tmp, &Tmp4, &Tmp);
					BigIntMultiply(&Mult3, &M3, &Tmp4);
					BigIntAdd(&Tmp, &Tmp4, &Tmp);
					BigIntMultiply(&Mult4, &M4, &Tmp4);
					BigIntAdd(&Tmp, &Tmp4, &Tmp);
					BigIntDivide(&Tmp, &K, &Tmp1);

					// Compute Tmp2 <- (Mult1*M2 - Mult2*M1 + Mult3*M4 - Mult4*M3) / K
					BigIntMultiply(&Mult1, &M2, &Tmp);
					BigIntMultiply(&Mult2, &M1, &Tmp4);
					BigIntSubt(&Tmp, &Tmp4, &Tmp);
					BigIntMultiply(&Mult3, &M4, &Tmp4);
					BigIntAdd(&Tmp, &Tmp4, &Tmp);
					BigIntMultiply(&Mult4, &M3, &Tmp4);
					BigIntSubt(&Tmp, &Tmp4, &Tmp);
					BigIntDivide(&Tmp, &K, &Tmp2);

					// Compute Tmp3 <- (Mult1*M3 - Mult3*M1 - Mult2*M4 + Mult4*M2) / K
					BigIntMultiply(&Mult1, &M3, &Tmp);
					BigIntMultiply(&Mult3, &M1, &Tmp4);
					BigIntSubt(&Tmp, &Tmp4, &Tmp);
					BigIntMultiply(&Mult2, &M4, &Tmp4);
					BigIntSubt(&Tmp, &Tmp4, &Tmp);
					BigIntMultiply(&Mult4, &M2, &Tmp4);
					BigIntAdd(&Tmp, &Tmp4, &Tmp);
					BigIntDivide(&Tmp, &K, &Tmp3);

					// Compute Mult4 <- (Mult1*M4 - Mult4*M1 + Mult2*M3 - Mult3*M2) / K
					BigIntMultiply(&Mult1, &M4, &Tmp);
					BigIntMultiply(&Mult4, &M1, &Tmp4);
					BigIntSubt(&Tmp, &Tmp4, &Tmp);
					BigIntMultiply(&Mult2, &M3, &Tmp4);
					BigIntAdd(&Tmp, &Tmp4, &Tmp);
					BigIntMultiply(&Mult3, &M2, &Tmp4);
					BigIntSubt(&Tmp, &Tmp4, &Tmp);
					BigIntDivide(&Tmp, &K, &Mult4);

					CopyBigInt(&Mult3, &Tmp3);
					CopyBigInt(&Mult2, &Tmp2);
					CopyBigInt(&Mult1, &Tmp1);
				} /* end while */
			} /* end if p = 3 (mod 4) */
		} /* end prime not 2 */

		  // Compute Tmp1 <- Mult1*Quad1 + Mult2*Quad2 + Mult3*Quad3 + Mult4*Quad4
		BigIntMultiply(&Mult1, &Quad1, &Tmp);
		BigIntMultiply(&Mult2, &Quad2, &Tmp4);
		BigIntAdd(&Tmp, &Tmp4, &Tmp);
		BigIntMultiply(&Mult3, &Quad3, &Tmp4);
		BigIntAdd(&Tmp, &Tmp4, &Tmp);
		BigIntMultiply(&Mult4, &Quad4, &Tmp4);
		BigIntAdd(&Tmp, &Tmp4, &Tmp1);

		// Compute Tmp2 <- Mult1*Quad2 - Mult2*Quad1 + Mult3*Quad4 - Mult4*Quad3
		BigIntMultiply(&Mult1, &Quad2, &Tmp);
		BigIntMultiply(&Mult2, &Quad1, &Tmp4);
		BigIntSubt(&Tmp, &Tmp4, &Tmp);
		BigIntMultiply(&Mult3, &Quad4, &Tmp4);
		BigIntAdd(&Tmp, &Tmp4, &Tmp);
		BigIntMultiply(&Mult4, &Quad3, &Tmp4);
		BigIntSubt(&Tmp, &Tmp4, &Tmp2);

		// Compute Tmp3 <- Mult1*Quad3 - Mult3*Quad1 - Mult2*Quad4 + Mult4*Quad2
		BigIntMultiply(&Mult1, &Quad3, &Tmp);
		BigIntMultiply(&Mult3, &Quad1, &Tmp4);
		BigIntSubt(&Tmp, &Tmp4, &Tmp);
		BigIntMultiply(&Mult2, &Quad4, &Tmp4);
		BigIntSubt(&Tmp, &Tmp4, &Tmp);
		BigIntMultiply(&Mult4, &Quad2, &Tmp4);
		BigIntAdd(&Tmp, &Tmp4, &Tmp3);

		// Compute Quad4 <- Mult1*Quad4 - Mult4*Quad1 + Mult2*Quad3 - Mult3*Quad2
		BigIntMultiply(&Mult1, &Quad4, &Tmp);
		BigIntMultiply(&Mult4, &Quad1, &Tmp4);
		BigIntSubt(&Tmp, &Tmp4, &Tmp);
		BigIntMultiply(&Mult2, &Quad3, &Tmp4);
		BigIntAdd(&Tmp, &Tmp4, &Tmp);
		BigIntMultiply(&Mult3, &Quad2, &Tmp4);
		BigIntSubt(&Tmp, &Tmp4, &Quad4);

		CopyBigInt(&Quad3, &Tmp3);
		CopyBigInt(&Quad2, &Tmp2);
		CopyBigInt(&Quad1, &Tmp1);
	} /* end for indexPrimes */
	pstFactor = pstFactors + 1;      // Point to first factor in array of factors.
	for (indexPrimes = pstFactors->multiplicity - 1; indexPrimes >= 0; indexPrimes--, pstFactor++)
	{
		NumberLength = *pstFactor->ptrFactor;
		UncompressBigInteger(pstFactor->ptrFactor, &p);
		BigIntPowerIntExp(&p, pstFactor->multiplicity / 2, &K);
		BigIntMultiply(&Quad1, &K, &Quad1);
		BigIntMultiply(&Quad2, &K, &Quad2);
		BigIntMultiply(&Quad3, &K, &Quad3);
		BigIntMultiply(&Quad4, &K, &Quad4);
	}
	Quad1.sign = SIGN_POSITIVE;
	Quad2.sign = SIGN_POSITIVE;
	Quad3.sign = SIGN_POSITIVE;
	Quad4.sign = SIGN_POSITIVE;
	// Sort squares
	BigIntSubt(&Quad1, &Quad2, &Tmp);
	if (Tmp.sign == SIGN_NEGATIVE)
	{   // Quad1 < Quad2, so exchange them.
		CopyBigInt(&Tmp, &Quad1);
		CopyBigInt(&Quad1, &Quad2);
		CopyBigInt(&Quad2, &Tmp);
	}
	BigIntSubt(&Quad1, &Quad3, &Tmp);
	if (Tmp.sign == SIGN_NEGATIVE)
	{   // Quad1 < Quad3, so exchange them.
		CopyBigInt(&Tmp, &Quad1);
		CopyBigInt(&Quad1, &Quad3);
		CopyBigInt(&Quad3, &Tmp);
	}
	BigIntSubt(&Quad1, &Quad4, &Tmp);
	if (Tmp.sign == SIGN_NEGATIVE)
	{   // Quad1 < Quad4, so exchange them.
		CopyBigInt(&Tmp, &Quad1);
		CopyBigInt(&Quad1, &Quad4);
		CopyBigInt(&Quad4, &Tmp);
	}
	BigIntSubt(&Quad2, &Quad3, &Tmp);
	if (Tmp.sign == SIGN_NEGATIVE)
	{   // Quad2 < Quad3, so exchange them.
		CopyBigInt(&Tmp, &Quad2);
		CopyBigInt(&Quad2, &Quad3);
		CopyBigInt(&Quad3, &Tmp);
	}
	BigIntSubt(&Quad2, &Quad4, &Tmp);
	if (Tmp.sign == SIGN_NEGATIVE)
	{   // Quad2 < Quad4, so exchange them.
		CopyBigInt(&Tmp, &Quad2);
		CopyBigInt(&Quad2, &Quad4);
		CopyBigInt(&Quad4, &Tmp);
	}
	BigIntSubt(&Quad3, &Quad4, &Tmp);
	if (Tmp.sign == SIGN_NEGATIVE)
	{   // Quad3 < Quad4, so exchange them.
		CopyBigInt(&Tmp, &Quad3);
		CopyBigInt(&Quad3, &Quad4);
		CopyBigInt(&Quad4, &Tmp);
	}
}

static void varSquared(char **pptrOutput, char letter, char sign)
{
	char *ptrOutput = *pptrOutput;
	*ptrOutput++ = ' ';
	*ptrOutput++ = letter;
	strcpy(ptrOutput, (prettyprint ? "�" : "^2"));
	ptrOutput += strlen(ptrOutput);
	*ptrOutput++ = ' ';
	*ptrOutput++ = sign;
	*pptrOutput = ptrOutput;
}

static void valueVar(char **pptrOutput, char letter, BigInteger *value)
{
	char *ptrOutput = *pptrOutput;
	strcpy(ptrOutput, "\n");
	ptrOutput += strlen(ptrOutput);
	*ptrOutput++ = letter;
	*ptrOutput++ = ' ';
	*ptrOutput++ = '=';
	*ptrOutput++ = ' ';
	if (hexadecimal)
	{
		BigInteger2Hex(value, ptrOutput, groupLen);
	}
	else
	{
		BigInteger2Dec(value, ptrOutput, groupLen);
	}
	ptrOutput += strlen(ptrOutput);
	strcpy(ptrOutput, "\n");
	ptrOutput += strlen(ptrOutput);
	*pptrOutput = ptrOutput;
}

static void ShowFourSquares(char **pptrOutput)
{
	char *ptrOutput = *pptrOutput;
	strcpy(ptrOutput, "\nn =");
	ptrOutput += strlen(ptrOutput);
	if (Quad4.nbrLimbs == 1 && Quad4.limbs[0].x == 0)
	{          // Quad4 equals zero.
		if (Quad3.nbrLimbs == 1 && Quad3.limbs[0].x == 0)
		{        // Quad3 and Quad4 equal zero.
			if (Quad2.nbrLimbs == 1 && Quad2.limbs[0].x == 0)
			{      // Quad2, Quad3 and Quad4 equal zero.
				varSquared(&ptrOutput, 'a', ' ');
				strcpy(ptrOutput, "\n");
				ptrOutput += strlen(ptrOutput);
				valueVar(&ptrOutput, 'a', &Quad1);
				*pptrOutput = ptrOutput;
				return;
			}
			varSquared(&ptrOutput, 'a', '+');
			varSquared(&ptrOutput, 'b', ' ');
			strcpy(ptrOutput, "\n");
			ptrOutput += strlen(ptrOutput);
			valueVar(&ptrOutput, 'a', &Quad1);
			valueVar(&ptrOutput, 'b', &Quad2);
			*pptrOutput = ptrOutput;
			return;
		}
		varSquared(&ptrOutput, 'a', '+');
		varSquared(&ptrOutput, 'b', '+');
		varSquared(&ptrOutput, 'c', ' ');
		strcpy(ptrOutput, "\n");
		ptrOutput += strlen(ptrOutput);
		valueVar(&ptrOutput, 'a', &Quad1);
		valueVar(&ptrOutput, 'b', &Quad2);
		valueVar(&ptrOutput, 'c', &Quad3);
		*pptrOutput = ptrOutput;
		return;
	}
	varSquared(&ptrOutput, 'a', '+');
	varSquared(&ptrOutput, 'b', '+');
	varSquared(&ptrOutput, 'c', '+');
	varSquared(&ptrOutput, 'd', ' ');
	strcpy(ptrOutput, "\n");
	ptrOutput += strlen(ptrOutput);
	valueVar(&ptrOutput, 'a', &Quad1);
	valueVar(&ptrOutput, 'b', &Quad2);
	valueVar(&ptrOutput, 'c', &Quad3);
	valueVar(&ptrOutput, 'd', &Quad4);
	*pptrOutput = ptrOutput;
}

/* factors numbers or numeric expressions using fast algorithms ECM and SIQS
tofactortext:	Expression or value to be processed, in the form of ascii text.
performFactorization: set to 1 if facorisation is to be done
				set to 0 if not (in this case only the expression is evaluated)
factors			known factors (otherwise NULL)
ptrOutput:		output buffer 
*/
void ecmFrontText(char *tofactorText, int performFactorization, char *factors, char *ptrOutput)
{

	knownFactors = factors;
	if (valuesProcessed == 0)
	{
		doFactorization = performFactorization;
	}

	/* convert toFactorText to a bigint */
	enum eExprErr rc = BatchProcessing(tofactorText, &tofactor, &ptrOutput);
	if (rc != EXPR_OK)
		return;
	if (doFactorization)
	{
		/*NumberLength = tofactor.nbrLimbs;
		CompressBigInteger(nbrToFactor, &tofactor);
		factor(&tofactor, nbrToFactor, factorsMod, astFactorsMod, knownFactors);*/
		PerformFactorization(&tofactor);
		SendFactorizationToOutput(astFactorsMod, &ptrOutput, 1);
	}
	if (valuesProcessed == 1)
	{
		if (rc == EXPR_OK && doFactorization)
		{
			if (tofactor.sign == SIGN_POSITIVE)
			{        // Number to factor is non-negative.
				ComputeFourSquares(astFactorsMod);
				ShowFourSquares(&ptrOutput);
				if (tofactor.nbrLimbs > 1 || tofactor.limbs[0].x > 0)
				{      // Number to factor is not zero.
					GetNumberOfDivisors(&ptrOutput);
					GetSumOfDivisors(&ptrOutput);
					GetEulerTotient(&ptrOutput);
					GetMobius(&ptrOutput);
				}
			}
			showElapsedTime(&ptrOutput);
		}
	}
	strcpy(ptrOutput, lang ? "\n" COPYRIGHT_SPANISH "\n" :
		"\n" COPYRIGHT_ENGLISH "\n");
}

/* process ascii text in inputstring. 
C	= continue button
otherwise inputstring consists of several null-terminated strings
format is:
grouplen,flags,?vpch etc
*/
void doWork(void)
{
	int flags;
	char *ptrData = inputString;
	char *ptrWebStorage, *ptrKnownFactors;
#ifdef __EMSCRIPTEN__
	originalTenthSecond = tenths();
#endif
	if (*ptrData == 'C')
	{    // User pressed Continue button.
		ecmFrontText(NULL, 0, NULL, output); // The 3rd parameter includes known factors.
#ifdef __EMSCRIPTEN__
		databack(output);
#endif
		return;
	}
	valuesProcessed = 0;
	groupLen = 0;
	while (*ptrData != ',')
	{
		groupLen = groupLen * 10 + (*ptrData++ - '0');
	}
	ptrData++;             // Skip comma.
	flags = *ptrData;
	if (flags == '-')
	{
		flags = -*(++ptrData);
	}
	lang = flags & 1;
	printf("flags = %d\n", flags);
	ptrData += 2;          // Skip app number and second comma.
	verbose = (*(ptrData + 1) == '1');
	prettyprint = (*(ptrData + 2) == '1');
	cunningham = (*(ptrData + 3) == '1');
	hexadecimal = (*(ptrData + 4) == '1');
	ptrData += 6;
	ptrWebStorage = ptrData + strlen(ptrData) + 1;
	ptrKnownFactors = findChar(ptrWebStorage, '=');
	if (prettyprint == 0)
	{
		groupLen = -groupLen;  // Do not show number of digts.
	}
	if (ptrKnownFactors)
	{
		ptrKnownFactors++;  // move pointer past '=' character
	}
	if (flags & 0x80)
	{
		if (ptrKnownFactors)
		{
			char *ptrText = ptrKnownFactors + strlen(ptrKnownFactors) + 1;
			NextEC = 0;
			while (*ptrText != 0)
			{
				NextEC = NextEC * 10 + (*ptrText++ & 0x0F);
			}
			flags = 2;  // do factorization.
		}
	}
	ecmFrontText(ptrData, flags & 2, ptrKnownFactors, output); // The 3rd parameter includes known factors.
#ifdef __EMSCRIPTEN__
	databack(output);
#endif
}