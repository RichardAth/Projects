/* program "i1.c" */
/*
*  A translation into C of the  pascal program ARITHMETIC, addition,
* subtraction, multiplication and exponentiation.
* from "SCIENTIFIC PASCAL" by H. FLANDERS pp. 342-357.
*/

#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_NONSTDC_NO_DEPRECATE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "integer.h"
#include "fun.h"

int MAX(int i, int j)
{
	return(i >= j ? i : j);
}

int MIN(int i, int j)
{
	return(i <= j ? i : j);
}

MPI *MAXMPI(MPI *I, MPI *J)
{
	return(COMPAREI(I, J) >= 0 ? COPYI(I) : COPYI(J));
}

MPI *MINMPI(MPI *I, MPI *J)
{
	return(COMPAREI(I, J) <= 0 ? COPYI(I) : COPYI(J));
}

int RSV(MPI *Aptr, MPI *Bptr)
{
	/*
	*  relative sizes of absolute values of MPI'S *Aptr and *Bptr.
	*				    1 if |*Aptr| > |*Bptr|
	*            RSV(Aptr, Bptr) =     0 if |*Aptr| = |*Bptr|
	*       			   -1 if |*Aptr| < |*Bptr|
	*/
	unsigned int j;

	if (Aptr->D > Bptr->D)
		return (1);
	else if (Aptr->D < Bptr->D)
		return (-1);
	else
	{
		j = Aptr->D;
		while ((Aptr->V[j] == Bptr->V[j]) && j)
			j--;
		if (Aptr->V[j] > Bptr->V[j])
			return (1);
		else if (Aptr->V[j] < Bptr->V[j])
			return (-1);
		else
			return (0);
	}
}

MPI *ADD0I(MPI *Aptr, MPI *Bptr)
/*
* input: Aptr and Bptr are pointers to non-negative MPI'S.
* output: *Eptr = *Aptr + *Bptr.
*/
{
	unsigned int j, m, n;
	unsigned long c, t;
	MPI *Eptr;
	MPI *tmp;

	m = (unsigned int)MIN((int)Aptr->D, (int)Bptr->D);
	n = (unsigned int)MAX((int)Aptr->D, (int)Bptr->D);
	c = 0;

	Eptr = BUILDMPI(n + 1);
	for (j = 0; j <= m; j++)
	{
		t = Aptr->V[j] + Bptr->V[j] + c;
		Eptr->V[j] = t & (R0 - 1);
		c = t >> T0;
	}
	if (Aptr->D > Bptr->D)
	{
		for (j = 1 + m; j <= Aptr->D; j++)
		{
			t = Aptr->V[j] + c;
			Eptr->V[j] = t & (R0 - 1);
			c = t >> T0;
		}
	}
	if (Aptr->D < Bptr->D)
	{
		for (j = 1 + m; j <= Bptr->D; j++)
		{
			t = Bptr->V[j] + c;
			Eptr->V[j] = t & (R0 - 1);
			c = t >> T0;
		}
	}
	if (c)
	{
		tmp = BANK_REALLOC(Eptr, n + 1);
		FREEMPI(Eptr);
		Eptr = tmp;
		Eptr->V[n + 1] = c;
	}
	Eptr->S = MAX(Aptr->S, Bptr->S);

	return  Eptr;
}

MPI *ADD0_I(MPI *Aptr, unsigned long b)
/*
* input: Aptr is a pointer to a non-negative MPI.
* b is a non-negative integer < R0.
* output: *Eptr = *Aptr + b.
*/
{
	unsigned int j, n;
	unsigned long c, t;
	MPI *Eptr, *tmp;

	if (b >= R0)
	{
		fprintf(stderr, "in ADD0_I, %lu >= R0\n", b);
		exit(1);
	}
	if (b == 0)
		return COPYI(Aptr);
	n = Aptr->D;
	Eptr = BUILDMPI(n + 1);
	c = b;
	for (j = 0; j <= n; j++)
	{
		t = Aptr->V[j] + c;
		Eptr->V[j] = t & (R0 - 1);
		c = t >> T0;
	}
	if (c)
	{
		tmp = Eptr;
		Eptr = BANK_REALLOC(Eptr, n + 1);
		FREEMPI(tmp);
		Eptr->V[n + 1] = c;
	}
	Eptr->S = b ? 1 : Aptr->S;

	return Eptr;
}

MPI *SUB0I(MPI *Aptr, MPI *Bptr)
/*
* input:  Aptr, Bptr, pointers to non-negative MPI'S, *Aptr >= *Bptr.
* output: MPI *Eptr = *Aptr - *Bptr.
*/
{
	int k;
	unsigned int d, e, j, m, n;
	unsigned long c, t;
	MPI *Eptr, *tmptr;

	if (Bptr->S == 0)
		return COPYI(Aptr);
	n = Aptr->D;
	Eptr = BUILDMPI(n + 1);
	m = Bptr->D;
	c = 1;

	/* c = 1: no borrow, c = 0: a borrow */

	for (j = 0; j <= m; j++)
	{
		t = Aptr->V[j] - Bptr->V[j] + R0 + c - 1;
		Eptr->V[j] = t & (R0 - 1);
		c = t >> T0;
	}
	if (n > m)
	{
		for (j = m + 1; j <= n; j++)
		{
			t = Aptr->V[j] + R0 + c - 1;
			Eptr->V[j] = t & (R0 - 1);
			c = t >> T0;
		}
	}
	/* finding Eptr->S and Eptr->D */

	Eptr->S = 0;
	d = 0;
	k = n;
	while (k >= 0)
	{
		if (Eptr->V[k] != 0)
		{
			Eptr->S = 1;
			d = k;
			k = -1;
		}
		else
			k--;
	}
	e = n - d;
	if (e)
	{
		tmptr = BANK_REALLOC(Eptr, d);
		FREEMPI(Eptr);
		Eptr = tmptr;
	}

	return Eptr;
}

MPI *SUB0_I(MPI *Aptr, unsigned long b)
/*
* input: Aptr a pointer to a non-negative MPI. b a non-negative integer,
* b < R0, *Aptr >= b;  output: *Eptr = *Aptr - b.
*/
{
	unsigned int e, j, n;
	unsigned long c, t;
	MPI *Eptr, *tmp;

	if (b >= R0)
	{
		fprintf(stderr, "in SUB0_I, %lu >= R0\n", b);
		exit(1);
	}
	if (b == 0)
		return COPYI(Aptr);
	n = Aptr->D;
	Eptr = BUILDMPI(n + 1);
	c = 1;
	/* c = 1: no borrow, c = 0: a borrow. */
	t = Aptr->V[0] - b + R0;
	Eptr->V[0] = t & (R0 - 1);
	c = t >> T0;
	for (j = 1; j <= n; j++)
	{
		t = Aptr->V[j] + R0 + c - 1;
		Eptr->V[j] = t & (R0 - 1);
		c = t >> T0;
	}
	if (Aptr->D)
		Eptr->S = 1;
	else
		Eptr->S = Eptr->V[0] ? 1 : 0;
	e = Eptr->V[n];
	if (e == 0 && n)
	{
		tmp = Eptr;
		Eptr = BANK_REALLOC(Eptr, n - 1);
		FREEMPI(tmp);
	}

	return Eptr;
}

MPI *ADDI(MPI *Aptr, MPI *Bptr)
/*
* input: Aptr and Bptr, pointers to MPI'S.
* output: MPI *Eptr = *Aptr + *Bptr.
*/
{
	MPI *Eptr;

	if (Aptr->S == Bptr->S)
	{
		Eptr = ADD0I(Aptr, Bptr);
		Eptr->S = Aptr->S;
	}
	else
	{
		switch (RSV(Aptr, Bptr))
		{
		case -1:
		{
			Eptr = SUB0I(Bptr, Aptr);
			Eptr->S = Bptr->S;
			break;
		}
		case 0:
		{
			Eptr = ZEROI();
			break;
		}
		default: /* case 1: */
		{
			Eptr = SUB0I(Aptr, Bptr);
			Eptr->S = Aptr->S;
		}
		}
	}
	return Eptr;
}

MPI *SUBI(MPI *Aptr, MPI *Bptr)
/*
* input: Aptr and Bptr, pointers to MPI'S.
* output: MPI *Eptr = *Aptr - *Bptr.
*/
{
	MPI *Eptr;

	Bptr->S = -Bptr->S;
	Eptr = ADDI(Aptr, Bptr);
	Bptr->S = -Bptr->S;
	return Eptr;
}

MPI *BANK_REALLOC(MPI *Aptr, unsigned int degree)
/*
* Reallocs using the memory bank. The original MPI is Aptr, and the new one
* has the given degree.
*/
{
	int size;
	MPI *Bptr;
#ifdef _WIN32
	unsigned int k;
#endif

	if (Aptr->D <= degree)
		size = Aptr->D;
	else
		size = degree;
	Bptr = BUILDMPI(1 + degree);
#ifdef _WIN32
	for (k = 0; k <= size; k++)
		Bptr->V[k] = Aptr->V[k];
#else
	/*
	bcopy((char *)(Aptr->V), (char *)(Bptr->V), (1 + size) * sizeof(unsigned long));
	*/
	Bptr->V = (USL *)memcpy((char *)(Bptr->V), (char *)(Aptr->V), (1 + size) * sizeof(unsigned long));
#endif
	Bptr->S = Aptr->S;

	return Bptr;
}

MPI *MULTI(MPI *Aptr, MPI *Bptr)
/*
* input: Aptr and Bptr, pointers to MPI'S.
* output: MPI *Eptr = *Aptr * *Bptr.
*/
{
	register unsigned int j, bk, k, m, n;
	unsigned long c = 0, t;
	int e;
	MPI *Eptr, *tmp;

	e = Aptr->S * Bptr->S;
	if (e == 0)
		return ZEROI();
	m = Aptr->D;
	n = Bptr->D;
	/*
	if (m < 1000 && n < 1000)
	{
	*/
	Eptr = BUILDMPI(m + n + 2);
	for (k = 0; k <= m; k++)
		Eptr->V[k] = 0;
	for (k = 0; k <= n; k++)
	{
		bk = Bptr->V[k];
		if (bk != 0)
		{
			c = 0;
			for (j = 0; j <= m; j++)
			{
				t = Aptr->V[j] * bk + Eptr->V[j + k] + c;
				/* t < R0 * R0 = 2 ^ O32 */
				Eptr->V[j + k] = t & (R0 - 1);
				c = t >> T0;	/* c < R0 */
			}
			Eptr->V[m + k + 1] = c;
		}
		else
			Eptr->V[m + k + 1] = 0;
	}
	if (c == 0)
	{
		tmp = BANK_REALLOC(Eptr, m + n);
		FREEMPI(Eptr);
		Eptr = tmp;
	}
	/*
	}
	*/
	/*
	else if (m < R0 && n < R0)
	Eptr = FFM(Aptr, Bptr, 3);*/ /* Fast Fourier Transform, 3 primes*/
	/*
	else
	Eptr = FFM(Aptr, Bptr, 4);*/ /* Fast Fourier Transform, 4 primes*/

	Eptr->S = e;
	return Eptr;
}

MPI *MULT_I(MPI *Aptr, long b)
/*
* input: Aptr is a pointer to an MPI.
* b is an integer, |b| < R0. output: *Eptr = *Aptr * b.
*/
{
	unsigned int j, n;
	unsigned long a, c = 0, t;
	MPI *Eptr;
	MPI *tmp;

	if (labs(b) >= R0)
	{
		fprintf(stderr, "in MULT_I, %ld >= R0\n", b);
		exit(1);
	}
	if (b == 0 || Aptr->S == 0)
		return ZEROI();
	if (b == 1)
		return COPYI(Aptr);
	if (b == -1) {
		Eptr = COPYI(Aptr);
		Eptr->S = -(Eptr->S);
		return COPYI(Eptr);
	}
	n = Aptr->D;
	Eptr = BUILDMPI(n + 1);
	Eptr->S = b > 0 ? Aptr->S : -Aptr->S;
	a = b > 0 ? b : -b;
	for (j = 0; j <= n; j++)
	{
		t = Aptr->V[j] * a + c;
		Eptr->V[j] = t & ((USL)R0 - 1);
		c = t >> T0;
	}
	if (c)
	{
		tmp = BANK_REALLOC(Eptr, n + 1);
		FREEMPI(Eptr);
		Eptr = tmp;
		Eptr->V[n + 1] = c;
	}

	return Eptr;
}

MPI *POWERI(MPI *Aptr, unsigned n)
/*
*  *Eptr = (*Aptr) ^ n, where 0 <= n < R0 * R0.
*/
{
	MPI *Btmp, *Ctmp, *Eptr;

	Eptr = ONEI();
	if (n == 0)
		return Eptr;
	Btmp = COPYI(Aptr);
	while (1)
	{
		if (n & 1)
		{
			Ctmp = Eptr;
			Eptr = MULTI(Ctmp, Btmp);
			FREEMPI(Ctmp);
			if (n == 1)
			{
				FREEMPI(Btmp);
				return Eptr;
			}
		}
		Ctmp = Btmp;
		Btmp = MULTI(Ctmp, Ctmp);
		FREEMPI(Ctmp);
		n = n >> 1;
	}
}

MPI *POWER_I(long a, unsigned n)
/*
* n is a non-negative integer, n < R0 * R0.
* a is an integer, |a| < R0.
* *Eptr = a ^ n.
*/
{
	MPI *Eptr, *Btmp, *Ctmp;

	if (labs(a) >= R0)
	{
		fprintf(stderr, "in POWER_I, %ld >= R0\n", a);
		exit(1);
	}
	if (a == 0)
		return ZEROI();
	if (n == 0)
		return ONEI();
	Eptr = ONEI();
	Btmp = BUILDMPI(1);
	Btmp->S = a > 0 ? 1 : -1;
	Btmp->V[0] = a > 0 ? a : -a;
	while (1)
	{
		if (n & 1)
		{
			Ctmp = Eptr;
			Eptr = MULTI(Ctmp, Btmp);
			FREEMPI(Ctmp);
			if (n == 1)
			{
				FREEMPI(Btmp);
				return Eptr;
			}
		}
		Ctmp = Btmp;
		Btmp = MULTI(Ctmp, Ctmp);
		FREEMPI(Ctmp);
		n = n >> 1;
	}
}

MPI *MULTI3(MPI *A, MPI *B, MPI *C)
/*
* Returns A * B * C.
*/
{
	MPI *Tmp1, *Tmp2;

	Tmp1 = MULTI(A, B);
	Tmp2 = MULTI(Tmp1, C);
	FREEMPI(Tmp1);
	return (Tmp2);
}

unsigned long MULT32(USL x, USL y)
/*
* Returns [x*y/2^32].
*/
{
	MPI *X, *Y, *Z;
	USL t;

	X = CHANGE(x);
	Y = CHANGE(y);
	FREEMPI(X);
	FREEMPI(Y);
	Z = MULTI(X, Y);
	if (Z->D <= 1)
		t = 0;
	else if (Z->D == 2)
		t = Z->V[2];
	else
		t = ((Z->V[3]) << 16) + Z->V[2];
	FREEMPI(Z);
	return (t);
}

MPI *MULTABC(MPI *A, MPI *B, MPI *C)
/*
* Returns A + B * C.
*/
{
	MPI *Temp1, *Temp2;

	Temp1 = MULTI(B, C);
	Temp2 = ADDI(A, Temp1);
	FREEMPI(Temp1);
	return (Temp2);
}

MPI *SHIFTLB(MPI *U)
/*
* Multiplies U by R0.
*/
{
	MPI *T;
	unsigned int n, i;

	n = U->D + 1;
	T = BUILDMPI(n + 1);
	T->S = U->S;
	for (i = 1; i <= n; i++)
		T->V[i] = U->V[i - 1];
	T->V[0] = 0;
	return (T);
}

MPI *MULT_II(MPI *Aptr, USL b)
/*
* output: MPI *Eptr = *Aptr * b, where b is an USL.
*/
{
	register unsigned int j, bk, k, m, n;
	unsigned long c = 0, t;
	int e;
	MPI *Eptr, *Bptr, *tmp;

	e = Aptr->S;
	if (e == 0 || b == 0)
		return ZEROI();

	if (b == 1)
		return COPYI(Aptr);
	m = Aptr->D;
	Bptr = CHANGE(b);
	n = Bptr->D;
	Eptr = BUILDMPI(m + n + 2);
	for (k = 0; k <= m; k++)
		Eptr->V[k] = 0;
	for (k = 0; k <= n; k++)
	{
		bk = Bptr->V[k];
		if (bk != 0)
		{
			c = 0;
			for (j = 0; j <= m; j++)
			{
				t = Aptr->V[j] * bk + Eptr->V[j + k] + c;
				/* t < R0 * R0 = 2 ^ O32 */
				Eptr->V[j + k] = t & (R0 - 1);
				c = t >> T0; 		/* c < R0 */
			}
			Eptr->V[m + k + 1] = c;
		}
		else
			Eptr->V[m + k + 1] = 0;
	}
	if (c == 0)
	{
		tmp = BANK_REALLOC(Eptr, m + n);
		FREEMPI(Eptr);
		Eptr = tmp;
	}
	Eptr->S = e;
	FREEMPI(Bptr);
	return Eptr;
}

int COMPAREI(MPI *Aptr, MPI *Bptr)
/*
* returns 1 if Aptr > Bptr, 0 if Aptr = Bptr, -1 if Aptr < Bptr.
*/
{
	MPI *Temp;
	int t;

	Temp = SUBI(Aptr, Bptr);
	t = Temp->S;
	FREEMPI(Temp);
	return (t);
}

int QSORTCOMPAREI(const void *Aptr, const void *Bptr)
/*
* returns 1 if Aptr > Bptr, 0 if Aptr = Bptr, -1 if Aptr < Bptr.
*/
{
	MPI *Temp;
	int t;

	Temp = SUBI(*(MPI **)Aptr, *(MPI **)Bptr);
	t = Temp->S;
	FREEMPI(Temp);
	return (t);
}

void QSORTMPI(MPI *A[], USI m, USI n)
/*
* Sorts the array of MPI's A[m],A[m+1],...,A[n] into
* increasing order.
*/
{
	USI t;

	t = n - m + 1;
	printf("before qsort...t=%d\n", t);
	qsort((char*)(A + m), t, sizeof(MPI *), QSORTCOMPAREI);
	printf("after qsort\n");
	return;
}

void QSORTMPI0(MPI *A[], USI m, USI n)
/*
* Sorts the array of non-negative MPI's A[m],A[m+1],...,A[n] into
* increasing order.
*/
{
	USI t;

	t = n - m + 1;
	qsort((char*)(A + m), t, sizeof(MPI *), QSORTRSV);
	return;
}

int QSORTRSV(const void *Aptr, const void *Bptr)
{
	/*
	*  relative sizes of absolute values of MPI'S **Aptr and **Bptr.
	*				    1 if |**Aptr| > |**Bptr|
	*            RSV(Aptr, Bptr) =     0 if |**Aptr| = |**Bptr|
	*       			   -1 if |**Aptr| < |**Bptr|
	* for use in QSORTMPI0().
	*/
	unsigned int j;

	if ((*(MPI **)Aptr)->D > (*(MPI **)Bptr)->D)
		return (1);
	else if ((*(MPI **)Aptr)->D < (*(MPI **)Bptr)->D)
		return (-1);
	else
	{
		j = (*(MPI **)Aptr)->D;
		while (((*(MPI **)Aptr)->V[j] == (*(MPI **)Bptr)->V[j]) && j)
			j--;
		if ((*(MPI **)Aptr)->V[j] >(*(MPI **)Bptr)->V[j])
			return (1);
		else if ((*(MPI **)Aptr)->V[j] < (*(MPI **)Bptr)->V[j])
			return (-1);
		else
			return (0);
	}
}

void QSORTMATI(MPMATI *Mptr, USI r, USI s)
/*
* Sorts the rows of *Mptr in increasing order of length.
* Based on Kochan, p.142-143, Programming in ANSI C.
*/
{
	USI i, j, m;
	MPI **TEMP, *temp, **R;
	m = Mptr->R;
	R = (MPI **)mmalloc(m * sizeof(MPI *));
	for (i = r; i <= s; i++)
		R[i] = DOTRI(Mptr, i, i);
	if (s > r)
	{
		for (i = r; i < s - 1; i++)
			for (j = i + 1; j <= s; j++)
			{
				if (RSV(R[i], R[j]) == 1)
				{
					temp = R[i];
					R[i] = R[j];
					R[j] = temp;
					TEMP = Mptr->V[i];
					Mptr->V[i] = Mptr->V[j];
					Mptr->V[j] = TEMP;
				}
			}
	}
	for (i = r; i <= s; i++)
		FREEMPI(R[i]);
	ffree((char *)R, m * sizeof(MPI *));
}

/* Returns respectively positive, zero, negative if MPI is positive, zero
* or negative */
int SIGNI(MPI *I)
{
	if (I == NULL)
		return 0;
	else
		return I->S;
}

MPI *MULTAB_PLUS_CDI(MPI *A, MPI *B, MPI *C, MPI *D)
/*
* Returns A*B + C*D.
*/
{
	MPI *Temp, *Temp1, *Temp2;

	Temp1 = MULTI(A, B);
	Temp2 = MULTI(C, D);
	Temp = ADDI(Temp1, Temp2);
	FREEMPI(Temp1);
	FREEMPI(Temp2);
	return (Temp);
}

MPI *MULTAB_MINUS_CDI(MPI *A, MPI *B, MPI *C, MPI *D)
/*
* Returns A*B - C*D.
*/
{
	MPI *Temp, *Temp1, *Temp2;

	Temp1 = MULTI(A, B);
	Temp2 = MULTI(C, D);
	Temp = SUBI(Temp1, Temp2);
	FREEMPI(Temp1);
	FREEMPI(Temp2);
	return (Temp);
}

MPI *ADD_ONEI(MPI *A) {
	MPI *ONE, *TEMP;

	ONE = ONEI();
	TEMP = ADDI(A, ONE);
	FREEMPI(ONE);
	return(TEMP);
}

MPI *SUB_ONEI(MPI *A) {
	MPI *ONE, *TEMP;

	ONE = ONEI();
	TEMP = SUBI(A, ONE);
	FREEMPI(ONE);
	return(TEMP);
}

/* AB +- C */
MPI *MULTAB_PLUS_MINUS_CI(MPI *A, MPI *B, int a, MPI *C) {
	MPI *TEMP, *TEMP1;

	TEMP = MULTI(A, B);
	if (a == 1) {
		TEMP1 = ADDI(TEMP, C);
	}
	else {
		TEMP1 = SUBI(TEMP, C);
	}
	FREEMPI(TEMP);
	return(TEMP1);
}

MPI *MIN_MPI_ARRAY(MPIA AA) {
	USI k, l;
	USL n;

	l = 0;
	n = AA->size;
	for (k = 1; k < n; k++) {
		if (COMPAREI(AA->A[k], AA->A[l]) < 0) {
			l = k;
		}
	}
	return(COPYI(AA->A[l]));
}
