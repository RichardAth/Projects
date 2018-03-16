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
#include "bignbr.h"
#include "factor.h"
#include "expression.h"
#include "showtime.h"

#define PAREN_STACK_SIZE 100

static BigInteger ReValue, ImValue;
static char *ptrOutput;
static struct sFactors astFactorsNorm[1000];
static int factorsNorm[10000];
static int NbrFactorsNorm;
static limb K[MAX_LEN];
static BigInteger mult1, mult2;
static limb minusOneMont[MAX_LEN];
static void DivideGaussian(BigInteger *real, BigInteger *imag);
static BigInteger value[2];
extern BigInteger tofactor;

static void w(char *text)
{
	strcpy(ptrOutput, text);
	ptrOutput += strlen(ptrOutput);
}

/* output value of complex number */
static void showNumber(BigInteger *real, BigInteger *imag)
{
	BigInteger Tmp;
	CopyBigInt(&Tmp, imag);
	if (real->sign == SIGN_NEGATIVE)
	{
		w("-");
	}
	//Bin2Dec(real->limbs, ptrOutput, real->nbrLimbs, groupLen);
	BigInteger2Dec(real, ptrOutput, groupLen);
	ptrOutput += strlen(ptrOutput);
	if (imag->sign == SIGN_POSITIVE)
	{
		w(" + ");
	}
	else
	{
		w(" - ");
	}
	//Bin2Dec(Tmp.limbs, ptrOutput, Tmp.nbrLimbs, groupLen);
	BigInteger2Dec(&Tmp, ptrOutput, groupLen);
	ptrOutput += strlen(ptrOutput);
	w(" i");
}

static void GaussianFactorization(void)
{
	BigInteger prime, q, r, M1, M2, Tmp;
	struct sFactors *pstFactor;

	BigIntMultiply(&ReValue, &ReValue, &tofactor);
	BigIntMultiply(&ImValue, &ImValue, &Tmp);
	BigIntAdd(&tofactor, &Tmp, &tofactor);
	NbrFactorsNorm = 0;
#ifdef __EMSCRIPTEN__
	originalTenthSecond = tenths();
#endif
	if (tofactor.nbrLimbs == 1 && tofactor.limbs[0].x == 0)
	{                // Norm is zero.
		w("<ul><li>Any gaussian prime divides this number</li></ul>");
		return;
	}
	w("<ul>");
	if (tofactor.nbrLimbs > 1 || tofactor.limbs[0].x > 1)
	{           // norm greater than 1. Factor norm.
		int index, index2;
		char *ptrFactorDec = tofactorDec;
		NumberLength = tofactor.nbrLimbs;
		CompressBigInteger(nbrToFactor, &tofactor);
		strcpy(ptrFactorDec, "Re&sup2; + Im&sup2; = ");
		ptrFactorDec += strlen(ptrFactorDec);
		//Bin2Dec(ReValue.limbs, ptrFactorDec, ReValue.nbrLimbs, groupLen);
		BigInteger2Dec(&ReValue, ptrFactorDec, groupLen);
		ptrFactorDec += strlen(ptrFactorDec);
		strcpy(ptrFactorDec, "&sup2; + ");
		ptrFactorDec += strlen(ptrFactorDec);
		//Bin2Dec(ImValue.limbs, ptrFactorDec, ImValue.nbrLimbs, groupLen);
		BigInteger2Dec(&ImValue, ptrFactorDec, groupLen);
		ptrFactorDec += strlen(ptrFactorDec);
		strcpy(ptrFactorDec, "&sup2;");
		ptrFactorDec += strlen(ptrFactorDec);
		factor(&tofactor, nbrToFactor, factorsNorm, astFactorsNorm, NULL);
		NbrFactorsNorm = astFactorsNorm[0].multiplicity;
		pstFactor = &astFactorsNorm[1];
		for (index = 0; index < NbrFactorsNorm; index++)
		{
			int *ptrPrime = pstFactor->ptrFactor;
			NumberLength = *ptrPrime;
			UncompressBigInteger(ptrPrime, &prime);
			prime.sign = SIGN_POSITIVE;
			if (prime.nbrLimbs == 1 && prime.limbs[0].x == 2)
			{             // Prime factor is 2.
				for (index2 = 0; index2 < pstFactor->multiplicity; index2++)
				{
					M1.nbrLimbs = M2.nbrLimbs = 1;
					M1.limbs[0].x = M2.limbs[0].x = 1;
					M1.sign = SIGN_POSITIVE;
					M2.sign = SIGN_NEGATIVE;
					DivideGaussian(&M1, &M1);           // Divide by 1+i
					DivideGaussian(&M1, &M2);           // Divide by 1-i
				}
			}
			if ((prime.limbs[0].x & 2) == 0)
			{                               // Prime is congruent to 1 (mod 4)
				CopyBigInt(&q, &prime);
				NumberLength = prime.nbrLimbs;
				memcpy(&TestNbr, prime.limbs, NumberLength * sizeof(limb));
				TestNbr[NumberLength].x = 0;
				GetMontgomeryParms(NumberLength);
				subtractdivide(&q, 1, 4);     // q = (prime-1)/4
				memset(&K, 0, NumberLength * sizeof(limb));
				memset(minusOneMont, 0, NumberLength * sizeof(limb));
				SubtBigNbrModN(minusOneMont, MontgomeryMultR1, minusOneMont, TestNbr, NumberLength);
				K[0].x = 1;
				do
				{    // Loop that finds mult1 = sqrt(-1) mod prime in Montgomery notation.
					K[0].x++;
					modPow(K, q.limbs, q.nbrLimbs, mult1.limbs);
				} while (!memcmp(mult1.limbs, MontgomeryMultR1, NumberLength * sizeof(limb)) ||
					!memcmp(mult1.limbs, minusOneMont, NumberLength * sizeof(limb)));
				K[0].x = 1;
				modmult(mult1.limbs, K, mult1.limbs);       // Convert mult1 to standard notation.
				UncompressLimbsBigInteger(mult1.limbs, &mult1);  // Convert to Big Integer.
				mult2.nbrLimbs = 1;                // mult2 <- 1
				mult2.limbs[0].x = 1;
				mult2.sign = SIGN_POSITIVE;
				for (;;)
				{
					// norm <- (mult1^2 + mult2^2) / prime
					BigIntMultiply(&mult1, &mult1, &tofactor);
					BigIntMultiply(&mult2, &mult2, &Tmp);
					BigIntAdd(&tofactor, &Tmp, &Tmp);
					BigIntDivide(&Tmp, &prime, &tofactor);
					if (tofactor.nbrLimbs == 1 && tofactor.limbs[0].x == 1)
					{        // norm equals 1.
						break;
					}
					BigIntRemainder(&mult1, &tofactor, &M1);
					BigIntRemainder(&mult2, &tofactor, &M2);
					BigIntAdd(&M1, &M1, &Tmp);
					BigIntSubt(&tofactor, &Tmp, &Tmp);
					if (Tmp.sign == SIGN_NEGATIVE)
					{
						BigIntSubt(&M1, &tofactor, &M1);
					}
					BigIntAdd(&M2, &M2, &Tmp);
					BigIntSubt(&tofactor, &Tmp, &Tmp);
					if (Tmp.sign == SIGN_NEGATIVE)
					{
						BigIntSubt(&M2, &tofactor, &M2);
					}
					// Compute q <- (mult1*M1 + mult2*M2) / norm
					BigIntMultiply(&mult1, &M1, &q);
					BigIntMultiply(&mult2, &M2, &Tmp);
					BigIntAdd(&q, &Tmp, &Tmp);
					BigIntDivide(&Tmp, &tofactor, &q);
					// Compute Mult2 <- (mult1*M2 - mult2*M1) / tofactor
					BigIntMultiply(&mult1, &M2, &r);
					BigIntMultiply(&mult2, &M1, &Tmp);
					BigIntSubt(&r, &Tmp, &Tmp);
					BigIntDivide(&Tmp, &tofactor, &mult2);
					CopyBigInt(&mult1, &q);
					mult1.sign = SIGN_POSITIVE;    // mult1 <- abs(mult1)
					mult2.sign = SIGN_POSITIVE;    // mult2 <- abs(mult2)
				}            /* end while */
				CopyBigInt(&M1, &mult1);
				CopyBigInt(&M2, &mult2);
				BigIntSubt(&M1, &M2, &Tmp);
				if (Tmp.sign == SIGN_NEGATIVE)
				{
					CopyBigInt(&Tmp, &mult1);
					CopyBigInt(&mult1, &mult2);
					CopyBigInt(&mult2, &Tmp);
				}
				for (index2 = 0; index2 < pstFactor->multiplicity; index2++)
				{
					DivideGaussian(&mult1, &mult2);
					BigIntNegate(&mult2, &Tmp);
					DivideGaussian(&mult1, &Tmp);
				}
			}              // end p = 1 (mod 4)
			else
			{              // if p = 3 (mod 4)
				q.nbrLimbs = 1;    // q <- 0
				q.limbs[0].x = 0;
				q.sign = SIGN_POSITIVE;
				for (index2 = 0; index2 < pstFactor->multiplicity; index2++)
				{
					DivideGaussian(&prime, &q);
				}            // end p = 3 (mod 4)
			}
			pstFactor++;
		}
	}
	// Process units: 1, -1, i, -i.
	if (ReValue.nbrLimbs == 1 && ReValue.limbs[0].x == 1)
	{
		if (ReValue.sign == SIGN_POSITIVE)
		{             // Value is 1.
			if (NbrFactorsNorm == 0)
			{
				w("No gaussian prime divides this number");
			}
		}
		else
		{            // Value is -1.
			w("<li>-1</li>");
		}
	}
	else if (ImValue.sign == SIGN_POSITIVE)
	{
		w("<li>i</li>");
	}
	else
	{
		w("<li>-i</li>");
	}
	w("</ul>");
}

static void DivideGaussian(BigInteger *real, BigInteger *imag)
{
	BigInteger Tmp, norm, realNum, imagNum;
	CopyBigInt(&Tmp, real);
	Tmp.sign = SIGN_POSITIVE;
	BigIntMultiply(real, real, &norm);
	BigIntMultiply(imag, imag, &Tmp);
	BigIntAdd(&norm, &Tmp, &norm);
	BigIntMultiply(&ReValue, real, &realNum);
	BigIntMultiply(&ImValue, imag, &Tmp);
	BigIntAdd(&realNum, &Tmp, &realNum);
	BigIntMultiply(&ImValue, real, &imagNum);
	BigIntMultiply(&ReValue, imag, &Tmp);
	BigIntSubt(&imagNum, &Tmp, &imagNum);
	BigIntRemainder(&realNum, &norm, &Tmp);
	if (Tmp.nbrLimbs == 1 && Tmp.limbs[0].x == 0)
	{
		BigIntRemainder(&imagNum, &norm, &Tmp);
		if (Tmp.nbrLimbs == 1 && Tmp.limbs[0].x == 0)
		{
			BigIntDivide(&realNum, &norm, &ReValue);
			BigIntDivide(&imagNum, &norm, &ImValue);
			w("<li>");
			showNumber(real, imag);
			w("</li>");
		}
	}
}

void gaussianText(char *valueText, int doFactorization)
{
	enum eExprErr rc;
#ifndef __EMSCRIPTEN__
	groupLen = 6;
#endif
	rc = ComputeGaussianExpression(valueText, value);
	output[0] = '2';
	ptrOutput = &output[1];
	if (rc == EXPR_OK)
	{
		CopyBigInt(&ReValue, &value[0]);
		CopyBigInt(&ImValue, &value[1]);
		if (doFactorization != '0')
		{
			w(lang ? "\nFactores de " : "\nFactors of ");
		}
		else
		{
			w(lang ? "\nEl valor es " : "\nValue is equal to ");
		}
		showNumber(&ReValue, &ImValue);
		if (doFactorization != '0')
		{
			GaussianFactorization();
		}
	}
	if (rc == EXPR_OK)
	{
		showElapsedTime(&ptrOutput);
	}
	else
	{
		textError(ptrOutput, rc);
		ptrOutput = output + strlen(output);
	}
	strcat(ptrOutput, lang ? "\n" COPYRIGHT_SPANISH "\n" :
		"\n" COPYRIGHT_ENGLISH "\n");
}

#ifdef __EMSCRIPTEN__
void doWork(void)
{
	int flags;
	char *ptrData = inputString;
	groupLen = 0;
	while (*ptrData != ',')
	{
		groupLen = groupLen * 10 + (*ptrData++ - '0');
	}
	ptrData++;             // Skip comma.
	flags = *ptrData;
	lang = flags & 1;
	ptrData += 2;          // Skip flags and comma.
	gaussianText(ptrData, (flags & 2) + '0');
	databack(output);
}
#endif