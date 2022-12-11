/* cyclo.c */

#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_NONSTDC_NO_DEPRECATE
#include <stdio.h>
#include <stdlib.h>
#include "integer.h"
#include "fun.h"
#define S0 65536

void NEXTPHI(USL p, USL *gptr, MPIA *b, MPIA *coef)
{
	USL cycle, grad, t, lastgrad;
	int i, j, next;
	MPI *l, *ZERO, *TEMP, *w;

	ZERO = ZEROI();
	cycle = p;
	grad = *gptr;
	next = grad;

	for (i = grad; i >= 0; i--) {
		if (cycle == p)
		{
			cycle = 1;
			TEMP = COPYI(((*coef)->A)[next]);
			ADD_TO_MPIA(*b, TEMP, i);
			FREEMPI(TEMP);
			next = next - 1;
		}
		else {
			cycle++;
			ADD_TO_MPIA(*b, ZERO, i);
		}
	}
	lastgrad = grad;
	grad = grad*(p - 1);

	t = grad / 2;
	for (i = grad; i >= t; i--) {
		/* PHI = PHI(x^p)/PHI.  First half of the coefficients. */
		l = COPYI(((*b)->A)[lastgrad]);
		ADD_TO_MPIA(*coef, l, i);
		for (j = lastgrad - 1; j >= 0; j--)
		{
			w = MULTI(l, ((*coef)->A)[j]);
			TEMP = SUBI(((*b)->A)[j], w);
			FREEMPI(w);
			ADD_TO_MPIA(*b, TEMP, j + 1);
			FREEMPI(TEMP);
		}
		FREEMPI(l);
		if (cycle == p)
		{
			cycle = 1;
			TEMP = COPYI(((*coef)->A)[next]);
			ADD_TO_MPIA(*b, TEMP, 0);
			FREEMPI(TEMP);
			next--;/*commented out by KRM on 1h26/9/02, reinstated 10th Feb 2012 by KRM */
		}
		else
		{
			cycle++;
			ADD_TO_MPIA(*b, ZERO, 0);
		}
	}
	FREEMPI(ZERO);

	t = grad / 2 - 1;
	for (j = 0; j <= t; j++) {
		TEMP = COPYI(((*coef)->A)[grad - j]);
		ADD_TO_MPIA(*coef, TEMP, j);
		FREEMPI(TEMP);
	}
	*gptr = grad;
}

int NOTYET(USL n, USL i, USL *qptr)
{
	USL q;

	q = n / i;
	*qptr = q;
	if ((n > q * i) && (i <= q))
		return 1;
	else
		return 0;
}

int LPF(USL n)
{
	USL c, i, q;

	if (n % 2 == 0) return (2);
	if (n % 3 == 0) return (3);
	else
	{
		i = 5;
		c = 0;
		while (NOTYET(n, i, &q))
		{
			i = (c == 0) ? i + 2 : i + 4;
			c = 1 - c;
		}
		if (q < i)
			return (n);
		else
			return (i);
	}
}

POLYI CYCLOTOMIC(MPI *N)
/* This algorithm is from Heinz Luneburg's book Galoisfelder, Kreisteilungs
* korper and Schieberegisterfolgen. It was originally coded in CMAT by perhaps Peter Adams.
* I imported the code to CALC on 25-26 September 2002, tightening up the code.
*/

{
	USL n, n1, socn, p, grad;
	int i;
	POLYI P;
	MPI *H, *TEMP, *MINUSONE, *ONE, *ZERO;
	MPIA b, coef;

	if (N->D) {
		printf("N >= 65536\n");
		return (POLYI)NULL;
	}
	MINUSONE = MINUS_ONEI();
	ONE = ONEI();
	ZERO = ZEROI();

	n = CONVERTI(N);
	b = BUILDMPIA();
	coef = BUILDMPIA();
	n1 = n;
	if (n1 == 1)
	{
		grad = 1;
		ADD_TO_MPIA(coef, MINUSONE, 0);
		ADD_TO_MPIA(coef, ONE, 1);
	}
	else
	{
		p = LPF(n1);
		socn = p;
		while (n1 % p == 0)
			n1 = n1 / p;
		grad = p - 1;
		for (i = 0; i <= grad; i++) {
			ADD_TO_MPIA(coef, ONE, i);
		}
		/* PHI is the pth cyclotomic polynomial. */

		while (n1 > 1)
		{
			p = LPF(n1);
			socn = socn * p;
			while (n1 % p == 0)
				n1 = n1 / p;
			NEXTPHI(p, &grad, &b, &coef);
			/* PHI is the socnth cyclotomic polynomial. */
		}
		n1 = n / socn;
		if (n1 > 1)
		{
			/* PHI = PHI(x^n). */
			grad = n1 * grad;
			for (i = grad; i >= 0; i--)
			{
				if (i % n1 == 0) {
					TEMP = COPYI((coef->A)[i / n1]);
					ADD_TO_MPIA(coef, TEMP, i);
					FREEMPI(TEMP);
				}
				else {
					ADD_TO_MPIA(coef, ZERO, i);
				}
			}
		}
	}
	FREEMPIA(b);
	P = NULL;
	for (i = 0; i <= grad; i++)
	{
		if (((coef->A)[i])->S)
		{
			H = COPYI((coef->A)[i]);
			PINSERTPI((int)i, H, &P, 0);
			FREEMPI(H);
		}
	}
	FREEMPIA(coef);
	FREEMPI(MINUSONE);
	FREEMPI(ONE);
	FREEMPI(ZERO);
	return (P);
}
