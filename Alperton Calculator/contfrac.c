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

#define _CRT_SECURE_NO_DEPRECATE

#include <string.h>
#include "bignbr.h"
#include "expression.h"
#include "highlevel.h"

extern int groupLen;

static BigInteger num, den, delta;
static BigInteger Temp;
static char *ptrOutput;
static void ShowRational(BigInteger *pNum, BigInteger *pDen);
static void showText(char *text)
{
	strcpy(ptrOutput, text);
	ptrOutput += strlen(ptrOutput);
}

// the continued fraction expansion of (num + sqrt(delta))/den is calculated.
static void ContFrac(void)
{
	BigInteger intSqrt;
	BigInteger biK, biP, biM;
	BigInteger K, L, M, P, Z;

	ptrOutput = output;
	// Show formula.
	showText("2\n<var>x</var> = <span class=\"fraction\"><span class=\"offscr\">");
	strcpy(ptrOutput, lang ? " la fracci�n cuyo numerador es </span>" : 
		" the fraction whose numerator is </span>");
	ptrOutput += strlen(ptrOutput);
	showText("<span class=\"fup\">");
	BigInteger2Dec(&num, ptrOutput, groupLen);    // Show numerator.
	ptrOutput += strlen(ptrOutput);
	//  strcpy(ptrOutput, lang ? " m�s la ra�z cuadrada de " : " plus the square root of ");
	//  ptrOutput += strlen(ptrOutput);
	showText(" + <span class=\"sqrtout\"><span class=\"sqrtin\">");
	BigInteger2Dec(&delta, ptrOutput, groupLen);  // Show radicand.
	ptrOutput += strlen(ptrOutput);
	showText("</span></span></span><span class=\"bar\"> </span><span class=\"fdn\"><span class=\"offscr\">");
	strcpy(ptrOutput, lang ? " y el denominador es </span>" : " and the denominator is </span>");
	ptrOutput += strlen(ptrOutput);
	BigInteger2Dec(&den, ptrOutput, groupLen);    // Show denominator.
	ptrOutput += strlen(ptrOutput);
	showText("</span></span></span>\n");
	// Validate input.
	if (den.nbrLimbs == 1 && den.limbs[0].x == 0)
	{
		showText(lang != 0 ? "\nError: El denominador es cero.\n" : "\nError: The denominator is zero.\n");
		return;
	}
	if (delta.sign == SIGN_NEGATIVE)
	{   /* Complex number */
		showText(lang != 0 ? "\nEl n�mero no es real, por lo que no tiene desarrollo en fracciones continuas.\n" :
			"\nThe number is not real, so it does not have continued fraction expansion.\n");
		return;
	}
	showText("\n<span class=\"offscr\">");
	showText(lang ? "El desarrollo en fracci�n continua de" :
		"The expansion in continued fraction of");
	showText(" </span><var>x</var> = ");
	if (delta.nbrLimbs == 1 && delta.limbs[0].x == 0)
	{   /* Rational number */
		ShowRational(&num, &den);
		return;
	}
	(void)BigIntMultiply(&num, &num, &Temp);
	BigIntSubt(&delta, &Temp, &Temp);  // Temp <- delta - num^2.
	(void)BigIntRemainder(&Temp, &den, &Temp);
	if (Temp.nbrLimbs != 1 || Temp.limbs[0].x != 0)
	{            // If delta - num^2 is not multiple of den...
		int sign;
		(void)BigIntMultiply(&delta, &den, &delta);
		(void)BigIntMultiply(&delta, &den, &delta);   // delta <- delta*den^2
		sign = den.sign;
		(void)BigIntMultiply(&num, &den, &num);       // num <- num*den
		(void)BigIntMultiply(&den, &den, &den);       // den <- den*den
		if (sign == SIGN_NEGATIVE)
		{
			BigIntNegate(&num, &num);             // Change sign to both
			BigIntNegate(&den, &den);             // numerator and denominator.
		}
	}
	squareRoot(delta.limbs, intSqrt.limbs, delta.nbrLimbs, &intSqrt.nbrLimbs);
	intSqrt.sign = SIGN_POSITIVE;
	(void)BigIntMultiply(&intSqrt, &intSqrt, &Temp);
	if (TestBigNbrEqual(&Temp, &delta) != 0)
	{     // delta is a perfect square, so number is rational.
		BigIntAdd(&num, &intSqrt, &Temp);
		ShowRational(&Temp, &den);
	}
	else
	{
		int cont;
		char *sep;

		CopyBigInt(&biP, &den);
		CopyBigInt(&biK, &intSqrt);
		if (biP.sign == SIGN_NEGATIVE)
		{
			addbigint(&biK, 1);
		}
		BigIntAdd(&biK, &num, &biK);
		if (biK.sign == SIGN_POSITIVE)
		{
			if (den.sign == SIGN_POSITIVE)
			{
				(void)BigIntDivide(&biK, &den, &biM);
			}
			else
			{
				CopyBigInt(&biM, &den);
				addbigint(&biM, 1);
				BigIntSubt(&biM, &biK, &biM);
				BigIntNegate(&den, &den);
				(void)BigIntDivide(&biM, &den, &biM);   // biM <- ((den+1)-biK)/(-den)
				BigIntNegate(&den, &den);
			}
		}
		else
		{
			if (den.sign == SIGN_POSITIVE)
			{
				CopyBigInt(&biM, &biK);
				addbigint(&biM, 1);
				BigIntSubt(&biM, &den, &biM);
				(void)BigIntDivide(&biM, &den, &biM);   // biM <- ((biK+1)-den)/den
			}
			else
			{
				BigIntNegate(&biK, &biM);
				BigIntNegate(&den, &den);
				(void)BigIntDivide(&biM, &den, &biM);   // biM <- (-biK)/(-den)
				BigIntNegate(&den, &den);
			}
		}

		BigInteger2Dec(&biM, ptrOutput, groupLen);  // Show integer part.
		ptrOutput += strlen(ptrOutput);
		(void)BigIntMultiply(&biM, &den, &biM);
		BigIntSubt(&biM, &num, &biM);
		sep = " + //";
		cont = -1;
		K.sign = L.sign = P.sign = M.sign = SIGN_NEGATIVE;
		K.nbrLimbs = L.nbrLimbs = P.nbrLimbs = M.nbrLimbs = 1;
		K.limbs[0].x = L.limbs[0].x = P.limbs[0].x = M.limbs[0].x = 1;   // K <- L <- P <- M <- -1. 
		while ((cont < 0 || TestBigNbrEqual(&K, &P) == 0 || TestBigNbrEqual(&L, &M) == 0) && cont < 100000)
		{
			showText(sep);
			sep = ", ";
			if (cont < 0 && biP.sign == SIGN_POSITIVE && biM.sign == SIGN_POSITIVE)
			{
				BigIntSubt(&biP, &intSqrt, &Temp);
				BigIntSubt(&Temp, &biM, &Temp);
				if (Temp.sign == SIGN_NEGATIVE || (Temp.nbrLimbs == 1 && Temp.limbs[0].x == 0))
				{
					BigIntSubt(&biM, &intSqrt, &Temp);
					if (Temp.sign == SIGN_NEGATIVE || (Temp.nbrLimbs == 1 && Temp.limbs[0].x == 0))
					{
						showText("<span class=\"offscr\">");
						showText(lang ? "inicio del per�odo" : "start periodic part");
						showText("</span><span class=\"bold\">");
						CopyBigInt(&P, &biP);
						CopyBigInt(&K, &P);
						CopyBigInt(&M, &biM);
						CopyBigInt(&L, &M);
						cont = 0;
					}
				}
			}
			if (cont >= 0)
			{
				/* both numerator and denominator are positive */
				(void)BigIntMultiply(&M, &M, &Temp);
				BigIntSubt(&delta, &Temp, &Temp);
				(void)BigIntDivide(&Temp, &P, &P);     // P <- (delta - M^2) / P
				BigIntAdd(&intSqrt, &M, &Z);
				(void)BigIntDivide(&Z, &P, &Z);        // Z <- (intSqrt + M) / P
				(void)BigIntMultiply(&Z, &P, &Temp);
				BigIntSubt(&Temp, &M, &M);       // M <- Z * P - M
				cont++;
			}
			else
			{
				(void)BigIntMultiply(&biM, &biM, &Temp);
				BigIntSubt(&delta, &Temp, &Temp);
				(void)BigIntDivide(&Temp, &biP, &biP);    // biP <- (delta - biM^2) / biP
				if (biP.sign == SIGN_POSITIVE)
				{
					BigIntAdd(&intSqrt, &biM, &Z);
					(void)BigIntDivide(&Z, &biP, &Z);      // Z <- (intSqrt + biM) / biP
				}
				else
				{
					CopyBigInt(&Z, &intSqrt);
					addbigint(&Z, 1);
					BigIntAdd(&Z, &biM, &Z);
					(void)BigIntDivide(&Z, &biP, &Z);      // Z <- (intSqrt + 1 + biM) / biP
				}
				(void)BigIntMultiply(&Z, &biP, &Temp);
				BigIntSubt(&Temp, &biM, &biM);     // biM <- Z * biP - biM
			}
			BigInteger2Dec(&Z, ptrOutput, groupLen);  // Show convergent.
			ptrOutput += strlen(ptrOutput);
		}
		if (cont >= 100000)
		{  // Too many convergents.
			showText(lang != 0 ? ", ... </span>//<br />donde la parte peri�dica (truncada a partir de los 100000 convergentes) est� se�alada en negrita.\n" :
				", ... </span>//<br />where the periodic part (truncated after 100000 convergents) is marked in bold.\n");
			return;
		}
		else
		{
			showText(lang != 0 ? "</span>//<br /><span aria-hidden=\"true\">donde la parte peri�dica est� se�alada en negrita</span>" :
				"</span>//<br /><span aria-hidden=\"true\">where the periodic part is marked in bold</span>");
			if (cont > 1)
			{
				showText(lang != 0 ? " (el per�odo tiene " : " (the period has ");
				int2dec(&ptrOutput, cont);
				showText(lang != 0 ? " coeficientes)" : " coefficients)");
			}
		}
	}
	showText(".");
	showText(lang ? "\n" COPYRIGHT_SPANISH "\n" :
		"\n" COPYRIGHT_ENGLISH "\n");
}

static void ShowRational(BigInteger *pNum, BigInteger *pDen)
{
	BigInteger Tmp;
	char *sep;

	BigIntGcd(pNum, pDen, &Tmp);     // Tmp <- GCD of numerator and denominator.
	(void)BigIntDivide(pNum, &Tmp, pNum);
	(void)BigIntDivide(pDen, &Tmp, pDen);
	if (pDen->sign == SIGN_NEGATIVE)
	{
		BigIntNegate(pNum, pNum);
		BigIntNegate(pDen, pDen);
	}
	(void)BigIntDivide(pNum, pDen, &Tmp); //  Get integer part.
	(void)BigIntRemainder(pNum, pDen, pNum);
	if (pNum->sign == SIGN_NEGATIVE)
	{
		BigIntAdd(pDen, pNum, pNum);   // Convert numerator to positive.
		addbigint(&Tmp, -1);          // Adjust integer part to floor.
	}
	BigInteger2Dec(&Tmp, ptrOutput, groupLen);  // Show convergent.
	ptrOutput += strlen(ptrOutput);
	(void)BigIntRemainder(pNum, pDen, pNum);
	sep = " + //";
	while (pNum->nbrLimbs != 1 || pNum->limbs[0].x != 0)
	{      // Numerator greater than zero.
		if (pDen->nbrLimbs != 1 || pDen->limbs[0].x != 0)
		{    // Denominator greater than zero.
			(void)BigIntDivide(pDen, pNum, &Tmp);
			showText(sep);
			BigInteger2Dec(&Tmp, ptrOutput, groupLen);  // Show convergent.
			ptrOutput += strlen(ptrOutput);
			sep = ", ";
		}
		(void)BigIntRemainder(pDen, pNum, &Tmp);
		CopyBigInt(pDen, pNum);
		CopyBigInt(pNum, &Tmp);
	}
	if (sep[0] == ',')
	{         // Inside continued fraction. Close it.
		showText("//");
	}
}

/* converts text in Input to BigInteger in Number. If no error occurs return 0.
If an error occurs return 1, after outputting the title and error text
to Output. */
static int getNumber(BigInteger *pNumber, char *title, char **pptrInput)
{
	enum eExprErr rc;
	rc = ComputeExpression(*pptrInput, 1, pNumber);
	if (rc != EXPR_OK)
	{
		showText(title);
		*ptrOutput++ = ':';
		*ptrOutput++ = ' ';
		textError(ptrOutput, rc);
		return 1;
	}
	*pptrInput += strlen(*pptrInput) + 1;  // Skip terminator.
	return 0;
}

// input contains three expressions separated by 00h (null character)
// The expressions are evaluated as num, delta and den.
// the continued fraction expansion of (num + sqrt(delta))/den is calculated.
void contfracText(char *input, int GroupLen)
{
	char *ptrInput = input;
	ptrOutput = output;
	if (getNumber(&num, lang != 0 ? "Numerador" : "Numerator", &ptrInput) != 0)
	{
		return;
	}
	if (getNumber(&delta, lang != 0 ? "Argumento de la ra�z cuadrada" : "Square root argument", &ptrInput) != 0)
	{
		return;
	}
	if (getNumber(&den, lang != 0 ? "Denominador" : "Denominator", &ptrInput) != 0)
	{
		return;
	}
	groupLen = GroupLen;  // note different capitalisation; 2 different variables.
	ContFrac();
}