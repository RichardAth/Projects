/* i5R.c */
/* programs for doing arithmetic with MPR's. */

#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_NONSTDC_NO_DEPRECATE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "integer.h"
#include "fun.h"

#ifdef DEBUG
extern long int nettbytes;
#endif

MPR *BUILDMPR()
/*
* mallocs space for an MPR.
*/
{
	MPR *Mptr;

	Mptr = (MPR *)mmalloc((USL)sizeof(MPR));
	return (Mptr);
}

void FREEMPR(MPR *Mptr)
/*
* deallocates the space alloted for the MPR Mptr.
*/
{
	if (Mptr == NULL)
		return;
	FREEMPI(Mptr->N);
	FREEMPI(Mptr->D);
	ffree((char *)Mptr, sizeof(MPR));
	return;
}

MPR *COPYR(MPR *Aptr)
/*
* returns a pointer to a copy of *Aptr.
*/
{
	MPR *Bptr;
	Bptr = BUILDMPR();
	Bptr->N = COPYI(Aptr->N);
	Bptr->D = COPYI(Aptr->D);
	return (Bptr);
}

unsigned int EQUALR(MPR *Aptr, MPR *Bptr)
/*
* returns 1 if *Aptr = *Bptr, 0 otherwise.
*/
{
	if (EQUALI(Aptr->N, Bptr->N) * EQUALI(Aptr->D, Bptr->D))
		return (1);
	else
		return (0);
}

MPR *MINUSR(MPR *Aptr)
/*
* *Bptr = -(*Aptr).
*/
{
	MPR *Bptr;

	Bptr = BUILDMPR();
	Bptr->N = MINUSI(Aptr->N);
	Bptr->D = COPYI(Aptr->D);
	return (Bptr);
}

MPR *RATIOI(MPI *Aptr, MPI *Bptr)
/*
* *Cptr = (*Aptr) / (*Bptr).
*/
{
	MPI *G;
	MPR *Cptr;

	Cptr = BUILDMPR();
	G = GCD(Aptr, Bptr);
	Cptr->N = INT(Aptr, G);
	Cptr->D = INT(Bptr, G);
	if ((Cptr->D)->S == -1)
	{
		(Cptr->D)->S = 1;
		(Cptr->N)->S = -((Cptr->N)->S);
	}
	FREEMPI(G);
	return (Cptr);
}

unsigned long LENGTHR(MPR *Mptr)
/*
* returns the sum of the lengths of numerator and denominator + 1,
* if Mptr is not an integer, otherwise returns the length of *Mptr.
*/
{
	if ((Mptr->D)->D || (Mptr->D)->V[0] > 1)
		return (LENGTHI(Mptr->N) + LENGTHI(Mptr->D) + 1);
	else
		return (LENGTHI(Mptr->N));
}

MPR *FINPUTR(FILE *f, unsigned int *uptr)
/*
* Converts the ratio of two decimal inputs from stream into an MPR.
* *uptr =  0 if input fails, 1 if successful.
*/
{
	int t, c;
	MPI *G, *H;
	MPR *Aptr;

	G = FINPUTI(f, uptr);
	if (*uptr == 0)
	{
		FREEMPI(G);
		return ZEROR();
	}
	t = c = fgetc(f);
	while (c == ' ')
		c = fgetc(f);
	if (c != '\n' && c != '/' && c != 'i'  && c != '+' && c != '-' && (c < '0' || c > '9'))
	{
		printf("inside INPUTSR: illegal character %c encountered:\n", c);
		FREEMPI(G);
		*uptr = 0;
		return ZEROR();
	}
	if (c == '/')
	{
		H = FINPUTI(f, uptr);
		if (*uptr == 0)
		{
			FREEMPI(G);
			FREEMPI(H);
			return ZEROR();
		}
		else if (H->S == 0)
		{
			printf("zero denominator entered\n");
			*uptr = 0;
			FREEMPI(G);
			FREEMPI(H);
			return ZEROR();
		}
		else
		{
			Aptr = RATIOI(G, H);
			FREEMPI(G);
			FREEMPI(H);
			c = fgetc(f);
			if (c != '\n')
				ungetc(c, f);
		}
	}
	else
	{
		if (c != '\n' || t == ' ')
			ungetc(c, f);
		Aptr = BUILDMPR();
		Aptr->N = G;
		Aptr->D = ONEI();
	}
	return (Aptr);
}

MPR *INPUTR(unsigned int *uptr)
{
	return FINPUTR(stdin, uptr);
}

MPR *INPUTSR(char **ptr, unsigned int *uptr)
/*
* Converts the ratio of two decimal inputs from stream into an MPR.
* *uptr = 0 if input fails, 1 if successful. For use with complex rationals.
*/
{
	char t;
	MPI *G, *H;
	MPR *Aptr;

	G = INPUTSI(ptr, uptr);
	if (*uptr == 0)
	{
		FREEMPI(G);
		return ZEROR();
	}
	t = **ptr;
	while (**ptr == ' ')
		(*ptr)++;
	if (**ptr != '\0' && **ptr != '/' && **ptr != 'i'  && **ptr != '+' && **ptr != '-' && (**ptr < '0' || **ptr > '9'))
	{
		printf("inside INPUTSR: illegal character %c entered:\n", **ptr);
		*uptr = 0;
		FREEMPI(G);
		return ZEROR();
	}
	else if (**ptr == '/')
	{
		(*ptr)++;
		H = INPUTSI(ptr, uptr);
		if (*uptr == 0)
		{
			FREEMPI(G);
			FREEMPI(H);
			return ZEROR();
		}
		else if (H->S == 0)
		{
			printf("zero denominator entered\n");
			*uptr = 0;
			FREEMPI(G);
			FREEMPI(H);
			return ZEROR();
		}
		else
		{
			Aptr = RATIOI(G, H);
			FREEMPI(G);
			FREEMPI(H);
		}
	}
	else
	{
		if (t == ' ')
			(*ptr)--;
		Aptr = BUILDMPR();
		Aptr->N = G;
		Aptr->D = ONEI();
	}
	return (Aptr);
}

void FPRINTR(FILE *outfile, MPR *Aptr)
/*
* prints the MPR *Aptr as (Aptr->N)/(Aptr->D).
*/
{
	FPRINTI(outfile, Aptr->N);
	if ((Aptr->D)->D || (Aptr->D)->V[0] > 1)
	{
		fprintf(outfile, "/");
		FPRINTI(outfile, Aptr->D);
	}
	return;
}

void PRINTR(MPR *Aptr)
{
	FPRINTR(stdout, Aptr);
	return;
}

MPR *ZEROR()
/*
* Returns the MPR ZERO.
*/
{
	MPR *Eptr;

	Eptr = BUILDMPR();
	Eptr->N = ZEROI();
	Eptr->D = ONEI();
	return Eptr;
}

MPR *ONER()
/*
* Returns the MPR ONE.
*/
{
	MPR *Eptr;

	Eptr = BUILDMPR();
	Eptr->N = ONEI();
	Eptr->D = ONEI();
	return Eptr;
}

MPR *TWOR()
/* Returns the MPR TWO */
{
	MPR *Eptr;
	Eptr = BUILDMPR();
	Eptr->N = TWOI();
	Eptr->D = ONEI();
	return Eptr;
}
MPR *MINUS_ONER()
/*
* Returns the MPR MINUS_ONE.
*/
{
	MPR *Eptr;

	Eptr = BUILDMPR();
	Eptr->N = MINUS_ONEI();
	Eptr->D = ONEI();
	return (Eptr);
}

void SPRINTR(char *buffer, MPR *Aptr)
/*
* prints the MPR *Aptr as (Aptr->N)/(Aptr->D).
*/
{
	SPRINTI(buffer, Aptr->N);
	if ((Aptr->D)->D || (Aptr->D)->V[0] > 1)
	{
		buffer += strlen(buffer);
		sprintf(buffer, "/");
		SPRINTI(buffer + 1, Aptr->D);
	}
	return;
}

void FPRINTDR(FILE *outfile, unsigned int d, MPR *Mptr)
/*
* prints *Mptr truncated to d (>= 1) decimal places to outfile.
*/
{
	int s;
	unsigned int i;
	unsigned long l;
	MPI *F, *G, *TempI;
	MPR *X, *TempR;

	s = (Mptr->N)->S;
	X = COPYR(Mptr);
	if (s == -1)
	{
		TempR = X;
		X = MINUSR(X);
		FREEMPR(TempR);
	}
	G = INT0(X->N, X->D);
	if (s == -1)
		fprintf(outfile, "-");
	FPRINTI(outfile, G);
	fprintf(outfile, ".");
	FREEMPI(G);
	G = MOD0(X->N, X->D);
	F = POWER_I(10L, d);
	TempI = F;
	F = MULTI(F, G);
	FREEMPI(TempI);
	TempI = G;
	G = INT0(F, X->D);
	FREEMPI(F);
	FREEMPI(TempI);
	l = LENGTHI(G);
	for (i = 1; i <= d - l; i++)
		fprintf(outfile, "0");
	FPRINTI(outfile, G);
	FREEMPI(G);
	FREEMPR(X);
	return;
}

void PRINTDR(unsigned int d, MPR *Mptr)
/*
* prints *Mptr truncated to d (>= 1) decimal places to stdout.
*/
{
	FPRINTDR(stdout, d, Mptr);
	return;
}

void SPRINTDR(char *buffer, unsigned int d, MPR *Mptr)
/*
* prints *Mptr truncated to d (>= 1) decimal places to buffer.
*/
{
	int s;
	unsigned int i;
	unsigned long l;
	MPI *F, *G, *TempI;
	MPR *X, *TempR;

	s = (Mptr->N)->S;
	X = COPYR(Mptr);
	if (s == -1)
	{
		TempR = X;
		X = MINUSR(X);
		FREEMPR(TempR);
	}
	G = INT0(X->N, X->D);
	if (s == -1)
	{
		sprintf(buffer, "-");
		buffer++;
	}
	SPRINTI(buffer, G);
	buffer += strlen(buffer);
	sprintf(buffer, ".");
	buffer++;
	FREEMPI(G);
	G = MOD0(X->N, X->D);
	F = POWER_I(10L, d);
	TempI = F;
	F = MULTI(F, G);
	FREEMPI(TempI);
	TempI = G;
	G = INT0(F, X->D);
	FREEMPI(TempI);
	FREEMPI(F);
	l = LENGTHI(G);
	for (i = 1; i <= d - l; i++)
	{
		sprintf(buffer, "0");
		buffer++;
	}
	SPRINTI(buffer, G);
	FREEMPI(G);
	FREEMPR(X);
}

unsigned long LENGTHDR(unsigned int d, MPR *Mptr)
/*
* returns LENGTHI(int(*Mptr)) + d + 1, if *Mptr >= 0.
* returns LENGTHI(int(abs(*Mptr))) + d + 2, if *Mptr < 0.
*/
{
	int s;
	unsigned long l;
	MPI *G;
	MPR *X, *TempR;

	s = (Mptr->N)->S;
	X = COPYR(Mptr);
	if (s == -1)
	{
		TempR = X;
		X = MINUSR(X);
		FREEMPR(TempR);
	}
	G = INT0(X->N, X->D);
	l = LENGTHI(G);
	FREEMPR(X);
	FREEMPI(G);
	if (s == -1)
		return (l + d + (USL)2);
	else
		return (l + d + (USL)1);
}

MPR *RECIPROCAL(unsigned long n)
/*
* input: unsigned int n, 0 < n < R0, output:  MPR *Aptr = 1 / n.
*/
{
	MPR *Aptr;

	Aptr = BUILDMPR();
	Aptr->N = ONEI();
	Aptr->D = CHANGE(n);
	return (Aptr);
}

MPR *ADDR(MPR *Aptr, MPR *Bptr)
/*
* *Eptr = *Aptr + *Bptr.
* P. Henrici's method: see Computer Algebra Symbolic and Algebraic
* Computation, Springer 1982,p. 200.
*/
{
	MPI *D, *E, *R, *S, *T1, *T2, *X, *Y;
	MPR *Eptr;

	if ((Aptr->N)->S == 0)
		return COPYR(Bptr);
	else if ((Bptr->N)->S == 0)
		return COPYR(Aptr);
	else
	{
		D = GCD(Aptr->D, Bptr->D);
		if (EQONEI(D))
		{
			FREEMPI(D);
			X = MULTI(Aptr->N, Bptr->D);
			Y = MULTI(Aptr->D, Bptr->N);
			Eptr = BUILDMPR();
			Eptr->N = ADDI(X, Y);
			Eptr->D = MULTI(Aptr->D, Bptr->D);
			FREEMPI(X);
			FREEMPI(Y);
			return (Eptr);
		}
		else
		{
			R = INT0(Aptr->D, D);
			S = INT0(Bptr->D, D);
			X = MULTI(Aptr->N, S);
			Y = MULTI(Bptr->N, R);
			T1 = ADDI(X, Y);
			T2 = MULTI(Aptr->D, S);
			FREEMPI(R);
			FREEMPI(S);
			FREEMPI(X);
			FREEMPI(Y);
			if (T1->S == 0)
			{
				FREEMPI(D);
				FREEMPI(T1);
				FREEMPI(T2);
				return ZEROR();
			}
			else
			{
				Eptr = BUILDMPR();
				E = GCD(T1, D);
				FREEMPI(D);
				if (EQONEI(E))
				{
					FREEMPI(E);
					Eptr->N = T1;
					Eptr->D = T2;
					return (Eptr);
				}
				else
				{
					Eptr->N = INT(T1, E);
					Eptr->D = INT(T2, E);
					FREEMPI(E);
					FREEMPI(T1);
					FREEMPI(T2);
					return (Eptr);
				}
			}
		}
	}
}

MPR *MULTR(MPR *Aptr, MPR *Bptr)
/*
* *Eptr = *Aptr * *Bptr.
* P. Henrici's method: see Computer Algebra Symbolic and Algebraic
* Computation, Springer 1982,p. 202.
*/
{
	MPI *D1 = NULL, *D2 = NULL, *R1 = NULL, *R2 = NULL, *S1 = NULL, *S2 = NULL;
	MPR *Eptr;
	if ((Aptr->N)->S == 0 || (Bptr->N)->S == 0)
		return ZEROR();
	else
	{
		D1 = GCD(Aptr->N, Bptr->D);
		D2 = GCD(Bptr->N, Aptr->D);
		if (EQONEI(D1))
		{
			R1 = COPYI(Aptr->N);
			S2 = COPYI(Bptr->D);
		}
		else
		{
			R1 = INT(Aptr->N, D1);
			S2 = INT0(Bptr->D, D1);
		}
		if (EQONEI(D2))
		{

			S1 = COPYI(Bptr->N);
			R2 = COPYI(Aptr->D);
		}
		else
		{
			S1 = INT(Bptr->N, D2);

			R2 = INT0(Aptr->D, D2);


		}
		Eptr = BUILDMPR();
		Eptr->N = MULTI(R1, S1);
		Eptr->D = MULTI(R2, S2);
		FREEMPI(D1);
		FREEMPI(D2);
		FREEMPI(R1);
		FREEMPI(R2);
		FREEMPI(S1);
		FREEMPI(S2);


		return (Eptr);
	}
}

MPR *POWERR(MPR *Aptr, unsigned int n)
/*
* *Eptr = (*Aptr) ^ n, where 0 <= n < R0 * R0.
*/
{
	MPR *B, *Temp, *Eptr;

	if (n == 0)
		return ONER();
	B = COPYR(Aptr);
	Eptr = ONER();
	while (1)
	{
		if (n & 1)
		{

			Temp = Eptr;
			Eptr = MULTR(Eptr, B);
			FREEMPR(Temp);
			if (n == 1)
			{
				FREEMPR(B);
				return (Eptr);
			}
		}
		Temp = B;
		B = MULTR(B, B);
		FREEMPR(Temp);
		n = n >> 1;
	}
}

/* As above but can find powers of negative numbers */
MPR *POWERR_2(MPR *R, int n)
{
	MPR *RETVAL;
	if (n < 0) {
		MPR *TMPR;
		RETVAL = POWERR(R, -n);
		TMPR = RETVAL;
		RETVAL = INVERSER(RETVAL);
		FREEMPR(TMPR);
	}
	else {
		RETVAL = POWERR(R, n);
	}
	return RETVAL;
}

MPR *SUBR(MPR *Aptr, MPR *Bptr)
/*
* *Eptr = *Aptr - *Bptr.
*/
{
	MPI *X, *Y, *U, *V, *G;
	MPR *Eptr;

	X = MULTI(Aptr->N, Bptr->D);
	Y = MULTI(Aptr->D, Bptr->N);
	V = MULTI(Aptr->D, Bptr->D);
	U = SUBI(X, Y);
	G = GCD(U, V);
	Eptr = BUILDMPR();
	Eptr->N = INT(U, G);
	Eptr->D = INT(V, G);
	FREEMPI(X);
	FREEMPI(Y);
	FREEMPI(U);
	FREEMPI(V);
	FREEMPI(G);
	return (Eptr);
}

MPR *INVERSER(MPR *Aptr)
/*
* *Eptr = 1 / *Aptr.
*/
{
	MPR *Eptr;
	int s;

	Eptr = BUILDMPR();
	Eptr->D = ABSI(Aptr->N);
	Eptr->N = COPYI(Aptr->D);
	s = (Aptr->N)->S;
	(Eptr->N)->S = ((Eptr->N)->S) * s;
	return (Eptr);
}

MPR *RATIOR(MPR *Aptr, MPR *Bptr)
/*
* *Eptr = *Aptr / *Bptr.
*/
{
	MPR *X, *Eptr;

	X = INVERSER(Bptr);
	Eptr = MULTR(Aptr, X);
	FREEMPR(X);
	return (Eptr);
}

MPR *FRAC_PARTR(MPR *Aptr)
/*
* *Bptr = fractional part of *Aptr.
*/
{
	MPI *X;
	MPR *Bptr;

	if (EQONEI(Aptr->D))
		Bptr = ZEROR();
	else
	{
		X = MOD(Aptr->N, Aptr->D);
		Bptr = RATIOI(X, Aptr->D);
		FREEMPI(X);
	}
	return (Bptr);
}

MPR *FRAC_PARTI(MPI *Aptr, MPI *Bptr)
/*
* *Cptr = fractional part of *Aptr/(*Bptr).
*/
{
	MPR *C, *Cptr;

	C = RATIOI(Aptr, Bptr);
	Cptr = FRAC_PARTR(C);
	FREEMPR(C);
	return (Cptr);
}

MPI *NEAREST_INTR(MPR *Aptr)
/*
* *Bptr is the nearest integer to *Aptr.
*/
{
	MPI *Y, *Z, *Temp, *Bptr;
	MPR *X;

	Y = INT(Aptr->N, Aptr->D);
	X = FRAC_PARTR(Aptr);
	Z = MULT_I(X->N, 2L);
	if (RSV(Z, X->D) <= 0)
		Bptr = Y;
	else
	{
		Temp = ONEI();
		Bptr = ADDI(Y, Temp);
		FREEMPI(Y);
		FREEMPI(Temp);
	}
	FREEMPR(X);
	FREEMPI(Z);
	return (Bptr);
}

MPI *ABS_NEAREST_INTR(MPR *Aptr)
/*
* *Bptr is the nearest integer to *Aptr, taking the integer closer to 0
* if half an odd integer.
*/
{
	MPI *Y, *Z, *Bptr, *O;
	MPR *X;
	int t;

	Y = INT(Aptr->N, Aptr->D);
	X = FRAC_PARTR(Aptr);
	Z = MULT_I(X->N, 2L);
	t = RSV(Z, X->D);
	FREEMPR(X);
	FREEMPI(Z);
	if (t < 0)
		Bptr = Y;
	else if (t == 0)
	{
		if (Y->S >= 0)
			Bptr = Y;
		else
		{
			O = ONEI();
			Bptr = ADDI(Y, O);
			FREEMPI(O);
			FREEMPI(Y);
		}
	}
	else
	{
		O = ONEI();
		Bptr = ADDI(Y, O);
		FREEMPI(Y);
		FREEMPI(O);
	}
	return (Bptr);
}

MPI *INPUTSID(char **ptr, unsigned int *uptr)
/*
* For use in converting the decimals after the decimal point to an MPI.
* *uptr is the number of decimals met after the decimal point.
*/
{
	unsigned int long n;
	char *t;
	MPI *Mptr, *Temp;

	Mptr = ZEROI();
	t = *ptr;
	while ((**ptr >= '0' && **ptr <= '9'))
	{
		Temp = Mptr;
		Mptr = MULT_I(Temp, 10L);
		FREEMPI(Temp);
		n = (unsigned long)(**ptr - '0');
		Temp = Mptr;
		Mptr = ADD0_I(Temp, n);
		FREEMPI(Temp);
		(*ptr)++;
	}
	*uptr = *ptr - t;
	return (Mptr);
}

MPR *INPUTSRD(char **ptr)
/*
* Converts a floating point decimal to an MPI. Used in inputting a file of
* matrices.
*/
{
	unsigned int n, length;
	MPI *G, *H, *K, *L, *M, *N;
	MPR *Aptr;

	G = INPUTSI(ptr, &n);
	if (n == 0)
	{
		FREEMPI(G);/* default - junk at start of number */
		return ZEROR();
	}
	else if (**ptr == '.')
	{
		(*ptr)++;
		H = INPUTSID(ptr, &length);
		K = POWER_I(10L, length);
		M = ABSI(G);
		L = MULTI(M, K);
		N = ADD0I(L, H);
		if (G->S == -1)
			N->S = -(N->S);
		Aptr = RATIOI(N, K);
		FREEMPI(H);
		FREEMPI(K);
		FREEMPI(M);
		FREEMPI(L);
		FREEMPI(N);
		FREEMPI(G);
	}
	else
	{
		Aptr = BUILDMPR();
		Aptr->N = G;
		Aptr->D = ONEI();
	}
	return (Aptr);
}

void FPRINTMATR(FILE *outfile, USI i1, USI i2, USI j1, USI j2, MPMATR *Mptr)
/*
* Modified 8/1/93 using Peter Adams' improvement from, August 1992.
*/
{
	char **str, *buff;
	unsigned int tmp, i, j, nrow, ncol, len, nstr = 0;
	int *colwidth;
	unsigned int ct;

	nrow = i2 - i1 + 1;
	ncol = j2 - j1 + 1;
	str = (char **)mmalloc(nrow * ncol * sizeof(char *));
	colwidth = (int *)mmalloc(ncol * sizeof(int));

	for (i = 0; i<ncol; i++)
		colwidth[i] = 0;
	for (i = i1; i <= i2; i++)
	{
		for (j = j1; j <= j2; j++)
		{
			tmp = 1 + LENGTHR(elt(Mptr, i, j));
			buff = (char *)mmalloc(tmp * sizeof(char));
			SPRINTR(buff, elt(Mptr, i, j));
			if ((len = strlen(buff)) > colwidth[j - j1])
				colwidth[j - j1] = len;
			str[nstr++] = strcpy((char *)mmalloc(len + 1), buff);
			ffree(buff, tmp * sizeof(char));
		}
	}
	if (outfile != stdout)
		fprintf(outfile, "%u %u\n", nrow, ncol);
	ct = 0;
	for (i = i1; i <= i2; i++)
	{
		for (j = j1; j <= j2; j++)
		{
			fprintf(outfile, "%*s ", colwidth[j - j1], str[ct]);
			ffree(str[ct], (unsigned int)(strlen(str[ct]) + 1) * sizeof(char));
			if ((ct % ncol) == (ncol - 1))
				fprintf(outfile, "\n");
			ct++;
		}
	}
	ffree((char *)str, nrow * ncol * sizeof(char *));
	ffree((char *)colwidth, ncol * sizeof(int));
	return;
}

void PRINTMATR(USI i1, USI i2, USI j1, USI j2, MPMATR *Mptr)
/*
* prints *Mptr from rows i1 to i2 and cols j1 to j2.
*/
{
	unsigned int i, j;

	i = i2 - i1 + 1;
	j = j2 - j1 + 1;
	if (MAX((int)i, (int)j) > 20)
	{
		printf("(The number of rows or columns to be printed exceeds 20;\n");
		printf("there is no point in printing this matrix on the screen.)\n");
		return;
	}
	FPRINTMATR(stdout, i1, i2, j1, j2, Mptr);
	return;
}


MPR *INTR(MPR *Aptr)
/*
*  Returns the integer part of *Aptr.
*/
{
	MPI *X;
	MPR *Bptr, *Temp;

	if (EQONEI(Aptr->D))
		Bptr = COPYR(Aptr);
	else
	{
		X = MOD(Aptr->N, Aptr->D);
		Temp = RATIOI(X, Aptr->D);
		FREEMPI(X);
		Bptr = SUBR(Aptr, Temp);
		FREEMPR(Temp);
	}
	return (Bptr);
}
