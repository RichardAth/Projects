/* i6R.c */
/* matrix operations using MPR's */

#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_NONSTDC_NO_DEPRECATE
#include <stdio.h>
/*#include <unistd.h>*/
#include <string.h>
#include <stdlib.h>
#include "integer.h"
#include "fun.h"
MPMATR *TRANSPOSER(MPMATR *Mptr);
MPMATR *MULTMATR(MPMATR *Mptr, MPMATR *Nptr);
extern USI GCDVERBOSE;
extern USI GCD3FLAG;
extern USI GCD_MAX;

MPMATR *BUILDMATR(USI m, USI n)
/*
* Allocates space for an m x n matrix of MPR's.
*/
{
	MPMATR *N;

	N = (MPMATR *)mmalloc(sizeof(MPMATR));
	N->R = m;
	N->C = n;
	N->V = (MPR **)mmalloc(m * n * sizeof(MPR *));
	return (N);
}

MPMATR *COPYMATR(MPMATR *Mptr)
/*
* a replacement for the assignment *Nptr = *Mptr.
*/
{
	unsigned int i, t;
	MPMATR *Nptr;

	t = Mptr->R * Mptr->C;
	Nptr = BUILDMATR(Mptr->R, Mptr->C);
	for (i = 0; i < t; i++)
		Nptr->V[i] = COPYR(Mptr->V[i]);
	return (Nptr);
}

void FREEMATR(MPMATR *Mptr)
/*
* frees the storage allocated to the two-dimensional array Mptr->V.
*/
{
	unsigned int i, t;

	if (Mptr == NULL)
		return;
	t = Mptr->R * Mptr->C;
	for (i = 0; i < t; i++)
		FREEMPR(Mptr->V[i]);
	ffree((char *)(Mptr->V), t * sizeof(MPR *));
	ffree((char *)Mptr, sizeof(MPMATR));
	return;
}

MPMATR *CHOLESKYR(MPMATR *A)
/*
* Input: A positive definite matrix A.
* Output: The Cholesky decomposition of A.
* See U.Finke and M. Pohst, "Improved methods for calculating vectors of
* short length in a lattice, including a complexity analysis." Math. Comp.
* 44, 463-471, 1985
*/
{
	MPMATR *Q;
	unsigned int i, j, k, l, m;
	MPR *temp1, *temp2;

	m = A->R;
	Q = ZEROMNR(m, m);
	for (i = 0; i < m; i++)
	{
		for (j = i; j < m; j++)
		{
			FREEMPR(elt(Q, i, j));
			elt(Q, i, j) = COPYR(elt(A, i, j));
		}
	}
	for (i = 0; i < m - 1; i++)
	{
		for (j = i + 1; j < m; j++)
		{
			FREEMPR(elt(Q, j, i));
			elt(Q, j, i) = COPYR(elt(Q, i, j));
			temp1 = elt(Q, i, j);
			elt(Q, i, j) = RATIOR(elt(Q, i, j), elt(Q, i, i));
			FREEMPR(temp1);
		}
		for (k = i + 1; k < m; k++)
		{
			for (l = k; l < m; l++)
			{
				temp1 = elt(Q, k, l);
				temp2 = MULTR(elt(Q, k, i), elt(Q, i, l));
				elt(Q, k, l) = SUBR(temp1, temp2);
				FREEMPR(temp1);
				FREEMPR(temp2);
			}
		}
	}
	for (j = 0; j < m - 1; j++)
	{
		for (i = j + 1; i < m; i++)
		{
			FREEMPR(elt(Q, i, j));
			elt(Q, i, j) = ZEROR();
		}
	}
	return (Q);
}

void FINCKE_POHST(MPMATR *A, MPR *C, USI filename)
/*
* Input: A matrix of integers A with LI rows spanning a lattice L.
* Output:  The integer vectors X with ||X||^2 <= C, highest nonzero coord <=0.
* See "Improved methods for calculating vectors of short length in a lattice,
* including a complexity analysis, U. Fincke and M. Pohst,
* Mathematics of Computation, 44, 1985, 463-471.
*/
{
	unsigned int i, j, n, flag, *l, *u, *t, *x, count = 0;
	unsigned int FLAG;
	MPR **X, **L, **T, **U, *Temp1, *Temp2, *Temp3, *Z, *ONE, *SUM;
	MPR *temp;
	char buff[20];
	FILE *outfile;
	MPMATR *B, *BB, *Q, *VECTOR, *MATR;
	MPMATI *MATI;
	int e;

	B = TRANSPOSER(A);
	BB = MULTMATR(A, B);
	FREEMATR(B);
	Q = CHOLESKYR(BB);
	FREEMATR(BB);
	if (filename == 1)
		strcpy(buff, "fp.out");
	else
		strcpy(buff, "slv.out");
	outfile = fopen(buff, "w");
	n = Q->R;
	i = n - 1;
	ONE = ONER();
	X = (MPR **)mmalloc(n * sizeof(MPR *));
	L = (MPR **)mmalloc(n * sizeof(MPR *));
	T = (MPR **)mmalloc(n * sizeof(MPR *));
	U = (MPR **)mmalloc(n * sizeof(MPR *));
	l = (USI *)mmalloc(n * sizeof(USI));
	t = (USI *)mmalloc(n * sizeof(USI));
	u = (USI *)mmalloc(n * sizeof(USI));
	x = (USI *)mmalloc(n * sizeof(USI));
	for (j = 0; j < n; j++)
	{
		l[j] = 0;
		t[j] = 0;
		u[j] = 0;
		x[j] = 0;
	}
	T[i] = COPYR(C); t[i] = 1;
	U[i] = ZEROR(); u[i] = 1;
bounds:
	Z = RATIOR(T[i], elt(Q, i, i));
	Temp2 = MINUSR(U[i]);
	if (l[i])
		FREEMPR(L[i]);
	else
		l[i] = 1;
	L[i] = INTROOT(Z, Temp2);
	FREEMPR(Temp2);
	Temp2 = INTROOT(Z, U[i]);
	Temp3 = MINUSR(Temp2);
	if (x[i])
		FREEMPR(X[i]);
	else
		x[i] = 1;
	X[i] = SUBR(Temp3, ONE);
	FREEMPR(Z);
	FREEMPR(Temp2);
	FREEMPR(Temp3);
loop:
	Temp1 = X[i];
	X[i] = ADDR(X[i], ONE);
	FREEMPR(Temp1);
	if (COMPAREI(X[i]->N, L[i]->N) == 1)
	{
		i++;
		goto loop;
	}
	else
	{
		if (i)
		{
			Temp1 = ADDR(X[i], U[i]);
			Temp2 = MULTR(Temp1, Temp1);
			FREEMPR(Temp1);
			Temp1 = MULTR(elt(Q, i, i), Temp2);
			FREEMPR(Temp2);
			if (t[i - 1])
				FREEMPR(T[i - 1]);
			else
				t[i - 1] = 1;
			T[i - 1] = SUBR(T[i], Temp1);
			FREEMPR(Temp1);
			i--;
			SUM = ZEROR();
			for (j = i + 1; j < n; j++)
			{
				Temp1 = SUM;
				Temp2 = MULTR(elt(Q, i, j), X[j]);
				SUM = ADDR(SUM, Temp2);
				FREEMPR(Temp1);
				FREEMPR(Temp2);
			}
			if (u[i])
				FREEMPR(U[i]);
			else
				u[i] = 1;
			U[i] = COPYR(SUM);
			FREEMPR(SUM);
			goto bounds;
		}
		else
			goto found;
	}

found:
	flag = 0;
	for (j = 0; j < n; j++)
	{
		if (!EQZEROR(X[j]))
		{
			flag = 1;
			break;
		}
	}
	if (flag == 0)
	{
		FREEMPR(ONE);
		for (j = 0; j < n; j++)
		{
			FREEMPR(X[j]);
			FREEMPR(L[j]);
			FREEMPR(T[j]);
			FREEMPR(U[j]);
		}
		ffree((char *)X, n * sizeof(MPR *));
		ffree((char *)L, n * sizeof(MPR *));
		ffree((char *)T, n * sizeof(MPR *));
		ffree((char *)U, n * sizeof(MPR *));
		ffree((char *)l, n * sizeof(USI));
		ffree((char *)t, n * sizeof(USI));
		ffree((char *)u, n * sizeof(USI));
		ffree((char *)x, n * sizeof(USI));
		fclose(outfile);
		FREEMATR(Q);
		return;
	}
	else
	{
		VECTOR = BUILDMATR(1, A->R);
		printf("LATTICE_VECTOR[%u] = ", count + 1);
		fprintf(outfile, "LATTICE_VECTOR[%u]:", count + 1);
		FLAG = 0;
		for (j = 0; j < A->R; j++)
		{
			e = (X[j]->N)->S;
			if (e == -1)
			{
				if (FLAG)
				{
					printf("+");
					fprintf(outfile, "+");
				}
				if (FLAG == 0)
					FLAG = 1;
				temp = MINUSR(X[j]);
				if (!EQONER(temp))
				{
					PRINTR(temp);
					FPRINTR(outfile, temp);
				}
				printf("b[%u]", j + 1);
				fprintf(outfile, "b[%u]", j + 1);
				FREEMPR(temp);
			}
			if (e == 1)
			{
				if (FLAG == 0)
					FLAG = 1;
				printf("-");
				fprintf(outfile, "-");
				if (!EQONER(X[j]))
				{
					PRINTR(X[j]);
					FPRINTR(outfile, X[j]);
				}
				printf("b[%u]", j + 1);
				fprintf(outfile, "b[%u]", j + 1);
			}
			elt(VECTOR, 0, j) = MINUSR(X[j]);
		}
		MATR = MULTMATR(VECTOR, A);
		FREEMATR(VECTOR);
		printf("=");
		fprintf(outfile, "=");
		MATI = BUILDMATI(1, A->C);
		for (j = 0; j < A->C; j++)
		{
			MATI->V[0][j] = COPYI((elt(MATR, 0, j))->N);
			PRINTI(MATI->V[0][j]);
			FPRINTI(outfile, MATI->V[0][j]);
			printf(" ");
			fprintf(outfile, " ");
		}
		printf(": ");
		fprintf(outfile, ": ");
		FREEMATR(MATR);
		FREEMATI(MATI);
		Temp1 = ADDR(X[0], U[0]);
		Temp2 = MULTR(Temp1, Temp1);
		FREEMPR(Temp1);
		Temp1 = MULTR(elt(Q, 0, 0), Temp2);
		FREEMPR(Temp2);
		Temp2 = SUBR(C, T[0]);
		Temp3 = ADDR(Temp2, Temp1);
		FREEMPR(Temp1);
		FREEMPR(Temp2);
		PRINTR(Temp3);
		FPRINTR(outfile, Temp3);
		printf("\n");
		fprintf(outfile, "\n");
		FREEMPR(Temp3);
		count++;
		goto loop;
	}
}

MPR *SLVECTOR(MPMATR *A, MPR *C, MPR **VALUE)
/*
* Input: A matrix of integers A with LI LLL reduced rows spanning a lattice L.
* C is the length-squared of a small lattive vector.
* Output: if 0 is returned, this means that VALUE is the shortest length.
* Then in nfunc.c, FINCKE-POHST is applied to get all the shortest vectors in
* L with highest nonzero coord > 0 are printed.
* Otherwise a shorter length is returned and VALUE will be NULL.
* See "Improved methods for calculating vectors of short length in a lattice,
* including a complexity analysis, U. Fincke and M. Pohst,
* Mathematics of Computation, 44, 1985, 463-471.
*/
{
	unsigned int i, j, n, flag, *l, *u, *t, *x, count = 0;
	MPR **X, **L, **T, **U, *Temp1, *Temp2, *Temp3, *Z, *ONE, *SUM;
	MPMATR *B, *BB, *Q;
	int s1, s2;

	*VALUE = (MPR *)NULL;
	B = TRANSPOSER(A);
	BB = MULTMATR(A, B);
	FREEMATR(B);
	Q = CHOLESKYR(BB);
	FREEMATR(BB);
	n = Q->R;
	i = n - 1;
	ONE = ONER();
	X = (MPR **)mmalloc(n * sizeof(MPR *));
	L = (MPR **)mmalloc(n * sizeof(MPR *));
	T = (MPR **)mmalloc(n * sizeof(MPR *));
	U = (MPR **)mmalloc(n * sizeof(MPR *));
	l = (USI *)mmalloc(n * sizeof(USI));
	t = (USI *)mmalloc(n * sizeof(USI));
	u = (USI *)mmalloc(n * sizeof(USI));
	x = (USI *)mmalloc(n * sizeof(USI));
	for (j = 0; j < n; j++)
	{
		l[j] = 0;
		t[j] = 0;
		u[j] = 0;
		x[j] = 0;
	}
	T[i] = COPYR(C); t[i] = 1;
	U[i] = ZEROR(); u[i] = 1;
bounds:
	Z = RATIOR(T[i], elt(Q, i, i));
	Temp2 = MINUSR(U[i]);
	if (l[i])
		FREEMPR(L[i]);
	else
		l[i] = 1;
	L[i] = INTROOT(Z, Temp2);
	FREEMPR(Temp2);
	Temp2 = INTROOT(Z, U[i]);
	Temp3 = MINUSR(Temp2);
	if (x[i])
		FREEMPR(X[i]);
	else
		x[i] = 1;
	X[i] = SUBR(Temp3, ONE);
	FREEMPR(Z);
	FREEMPR(Temp2);
	FREEMPR(Temp3);
loop:
	Temp1 = X[i];
	X[i] = ADDR(X[i], ONE);
	FREEMPR(Temp1);
	if (COMPAREI(X[i]->N, L[i]->N) == 1)
	{
		i++;
		goto loop;
	}
	else
	{
		if (i)
		{
			Temp1 = ADDR(X[i], U[i]);
			Temp2 = MULTR(Temp1, Temp1);
			FREEMPR(Temp1);
			Temp1 = MULTR(elt(Q, i, i), Temp2);
			FREEMPR(Temp2);
			if (t[i - 1])
				FREEMPR(T[i - 1]);
			else
				t[i - 1] = 1;
			T[i - 1] = SUBR(T[i], Temp1);
			FREEMPR(Temp1);
			i--;
			SUM = ZEROR();
			for (j = i + 1; j < n; j++)
			{
				Temp1 = SUM;
				Temp2 = MULTR(elt(Q, i, j), X[j]);
				SUM = ADDR(SUM, Temp2);
				FREEMPR(Temp1);
				FREEMPR(Temp2);
			}
			if (u[i])
				FREEMPR(U[i]);
			else
				u[i] = 1;
			U[i] = COPYR(SUM);
			FREEMPR(SUM);
			goto bounds;
		}
		else
			goto found;
	}

found:
	flag = 0;
	for (j = 0; j < n; j++)
	{
		if (!EQZEROR(X[j]))
		{
			flag = 1;
			break;
		}
	}
	if (flag == 0)
	{
		FREEMPR(ONE);
		for (j = 0; j < n; j++)
		{
			FREEMPR(X[j]);
			FREEMPR(L[j]);
			FREEMPR(T[j]);
			FREEMPR(U[j]);
		}
		ffree((char *)X, n * sizeof(MPR *));
		ffree((char *)L, n * sizeof(MPR *));
		ffree((char *)T, n * sizeof(MPR *));
		ffree((char *)U, n * sizeof(MPR *));
		ffree((char *)l, n * sizeof(USI));
		ffree((char *)t, n * sizeof(USI));
		ffree((char *)u, n * sizeof(USI));
		ffree((char *)x, n * sizeof(USI));
		FREEMATR(Q);
		return (ZEROR());
	}
	else
	{
		Temp1 = ADDR(X[0], U[0]);
		Temp2 = MULTR(Temp1, Temp1);
		FREEMPR(Temp1);
		Temp1 = MULTR(elt(Q, 0, 0), Temp2);
		FREEMPR(Temp2);
		Temp2 = SUBR(C, T[0]);
		Temp3 = ADDR(Temp2, Temp1);
		FREEMPR(Temp1);
		FREEMPR(Temp2);
		s1 = (Temp3->N)->S;
		if (s1)
		{
			if (*VALUE != NULL)
				FREEMPR(*VALUE);
			*VALUE = COPYR(Temp3);
		}
		s2 = RSV(Temp3->N, C->N);
		if (s1 && s2 == -1)
		{
			FREEMPR(ONE);
			for (j = 0; j < n; j++)
			{
				FREEMPR(X[j]);
				FREEMPR(L[j]);
				FREEMPR(T[j]);
				FREEMPR(U[j]);
			}
			ffree((char *)X, n * sizeof(MPR *));
			ffree((char *)L, n * sizeof(MPR *));
			ffree((char *)T, n * sizeof(MPR *));
			ffree((char *)U, n * sizeof(MPR *));
			ffree((char *)l, n * sizeof(USI));
			ffree((char *)t, n * sizeof(USI));
			ffree((char *)u, n * sizeof(USI));
			ffree((char *)x, n * sizeof(USI));
			FREEMATR(Q);
			return(Temp3);
		}
		FREEMPR(Temp3);
		count++;
		goto loop;
	}
}

MPR *INTROOT(MPR *Z, MPR *U)
/*
* Returns [sqrt(Z)+U]. First ANSWER = [sqrt(Z)] + [U]. One then
* tests if Z < ([sqrt(Z)] + 1 -{U})^2. If this does not hold, ANSWER++.
* For use in FINCKE_POHST().
*/
{
	int t;
	MPR *X, *Y, *ANSWER, *R, *S, *T1, *T2, *ONE, *Temp;

	if (EQZEROR(Z))
		return (INTR(U));
	ONE = ONER();
	X = MTHROOTR(Z, 2, 0);
	Y = INTR(U);
	ANSWER = ADDR(X, Y);
	R = SUBR(U, Y);
	FREEMPR(Y);
	S = SUBR(ONE, R);
	FREEMPR(R);
	T1 = ADDR(X, S);
	FREEMPR(S);
	FREEMPR(X);
	T2 = MULTR(T1, T1);
	FREEMPR(T1);
	t = COMPARER(T2, Z);
	FREEMPR(T2);
	if (t <= 0)
	{
		Temp = ANSWER;
		ANSWER = ADDR(ANSWER, ONE);
		FREEMPR(Temp);
	}
	FREEMPR(ONE);
	return (ANSWER);
}

int COMPARER(MPR *Aptr, MPR *Bptr)
/*
*				    1 if *Aptr > *Bptr
*            Returns               0 if *Aptr = *Bptr
*       			   -1 if *Aptr < *Bptr
*/
{
	MPI *T1, *T2;
	int t;

	T1 = MULTI(Aptr->N, Bptr->D);
	T2 = MULTI(Aptr->D, Bptr->N);
	t = COMPAREI(T1, T2);
	FREEMPI(T1);
	FREEMPI(T2);
	return (t);
}

MPMATR *ZEROMNR(USI m, USI n)
/*
* returns the zero  m x n matrix.
*/
{
	unsigned int i, j;
	MPMATR *Mptr;

	Mptr = BUILDMATR(m, n);
	for (i = 0; i <= m - 1; i++)
	{
		for (j = 0; j <= n - 1; j++)
			elt(Mptr, i, j) = ZEROR();
	}
	return (Mptr);
}

MPMATR *TRANSPOSER(MPMATR *Mptr)
/*
* returns the transpose of *Mptr.
*/
{
	MPMATR *Nptr;
	unsigned int i, j;

	Nptr = BUILDMATR(Mptr->C, Mptr->R);
	for (j = 0; j <= Nptr->R - 1; j++)
		for (i = 0; i <= Nptr->C - 1; i++)
			elt(Nptr, j, i) = COPYR(elt(Mptr, i, j));
	return (Nptr);
}

MPMATR *MULTMATR(MPMATR *Mptr, MPMATR *Nptr)
/*
* Here Mptr->C = Nptr->R.
* returns (*Mptr) * (*Nptr).
*/
{
	MPR *X, *Y, *Temp;
	unsigned int i, j, k;
	MPMATR *Lptr;

	Lptr = BUILDMATR(Mptr->R, Nptr->C);
	for (i = 0; i <= Mptr->R - 1; i++)
	{
		for (j = 0; j <= Nptr->C - 1; j++)
		{
			X = ZEROR();
			for (k = 0; k <= Mptr->C - 1; k++)
			{
				Y = MULTR(elt(Mptr, i, k), elt(Nptr, k, j));
				Temp = X;
				X = ADDR(X, Y);
				FREEMPR(Temp);
				FREEMPR(Y);
			}
			elt(Lptr, i, j) = X;
		}
	}
	return (Lptr);
}

void SHORTEST(MPMATI *AA, MPI *I[], USI filename, USI number)
/*
* Input: An m x m matrix of integers AA = Q arising from the EXTGCD
* LLLGCD or LLLGCD0. The last row of AA is assumed to be a multiplier
* of shortest length.
* Output:  All multipliers of shortest length.
* filename:1->egcdmult.out,2->lllgcdmult.out,3->lllgcd0mult.out.
* 4->gcd3.out (used in GCD3()).
* If number = 0, no printing.
* If number > 0 prints out multipliers on screen and in .out file.
* Only used in nfunc.c.
*/
{
	unsigned int i, j, k, m, n, *l, *u, *t, *x, w, count = 0;
	MPR **X, **Y, **L, **T, **U, **V, *Temp0, *Temp1, *Temp2, *Temp3;
	MPR *C, *Z, *ONE, *SUM, **N, *temp, *tempR, **XX;
	MPI *lengthi;
	char buff[25];
	FILE *outfile;
	MPMATR *A, *B, *BB, *Q, *QQ, *MATR, *VECTOR;
	MPMATI *MATI;
	int e;

	/* First get the vector [V[0],...,V[n-2], where V[i]=<A[n-1][*],A[i][*]>.*/
	m = AA->C;
	n = (AA->R) - 1;
	A = BUILDMATR(n + 1, m);
	for (i = 0; i < n + 1; i++)
		for (j = 0; j < m; j++)
		{
			elt(A, i, j) = BUILDMPR();
			elt(A, i, j)->N = COPYI(AA->V[i][j]);
			elt(A, i, j)->D = ONEI();
		}

	V = (MPR **)mmalloc(n * sizeof(MPR *));
	for (i = 0; i < n; i++)
	{
		V[i] = ZEROR();
		for (j = 0; j < m; j++)
		{
			Temp1 = V[i];
			Temp2 = MULTR(elt(A, n, j), elt(A, i, j));
			V[i] = ADDR(V[i], Temp2);
			FREEMPR(Temp1);
			FREEMPR(Temp2);
		}
	}
	B = TRANSPOSER(A);
	BB = MULTMATR(A, B);
	FREEMATR(B);
	Q = CHOLESKYR(BB);
	FREEMATR(BB);
	QQ = COPYMATR(Q);
	for (i = 0; i < n; i++)
	{
		FREEMPR(elt(QQ, i, i));
		elt(QQ, i, i) = ONER();
	}
	/* Calculate (U^{-1})^tV = Y */
	MATR = QQ;
	QQ = TRANSPOSER(QQ);
	FREEMATR(MATR);
	Y = (MPR **)mmalloc(n * sizeof(MPR *));
	Y[0] = COPYR(V[0]);
	for (w = 1; w < n; w++)
	{
		SUM = ZEROR();
		for (k = 0; k < w; k++)
		{
			Temp1 = SUM;
			Temp2 = MULTR(elt(QQ, w, k), Y[k]);
			SUM = ADDR(SUM, Temp2);
			FREEMPR(Temp1);
			FREEMPR(Temp2);
		}
		Y[w] = SUBR(V[w], SUM);
		FREEMPR(SUM);
	}
	FREEMATR(QQ);
	N = (MPR **)mmalloc(n * sizeof(MPR *));
	for (j = 0; j < n; j++)
		N[j] = RATIOR(Y[j], elt(Q, j, j));
	if (filename == 1)
		strcpy(buff, "egcdmult.out");
	else if (filename == 2)
		strcpy(buff, "lllgcdmult.out");
	else if (filename == 3)
		strcpy(buff, "lllgcd0mult.out");
	else if (filename == 4)
		strcpy(buff, "gcd3.out");
	else if (filename == 5)
		strcpy(buff, "axb.out");
	outfile = fopen(buff, "a");
	ONE = ONER();
	X = (MPR **)mmalloc(n * sizeof(MPR *));
	L = (MPR **)mmalloc(n * sizeof(MPR *));
	T = (MPR **)mmalloc(n * sizeof(MPR *));
	U = (MPR **)mmalloc(n * sizeof(MPR *));
	l = (USI *)mmalloc(n * sizeof(USI));
	t = (USI *)mmalloc(n * sizeof(USI));
	u = (USI *)mmalloc(n * sizeof(USI));
	x = (USI *)mmalloc(n * sizeof(USI));
	for (j = 0; j < n; j++)
	{
		l[j] = 0;
		t[j] = 0;
		u[j] = 0;
		x[j] = 0;
	}
	C = ZEROR();
	for (i = 0; i < n; i++)
	{
		Temp1 = C;
		Temp2 = MULTR(Y[i], Y[i]);
		Temp3 = RATIOR(Temp2, elt(Q, i, i));
		C = ADDR(C, Temp3);
		FREEMPR(Temp1);
		FREEMPR(Temp2);
		FREEMPR(Temp3);
	}
	i = n - 1;
	T[i] = COPYR(C); t[i] = 1;
	FREEMPR(C);
	U[i] = ZEROR(); u[i] = 1;
bounds:
	Z = RATIOR(T[i], elt(Q, i, i));
	Temp2 = SUBR(N[i], U[i]);
	if (l[i])
		FREEMPR(L[i]);
	else
		l[i] = 1;
	L[i] = INTROOT(Z, Temp2);
	FREEMPR(Temp2);
	Temp0 = SUBR(U[i], N[i]);
	Temp2 = INTROOT(Z, Temp0);
	FREEMPR(Temp0);
	Temp3 = MINUSR(Temp2);
	if (x[i])
		FREEMPR(X[i]);
	else
		x[i] = 1;
	X[i] = SUBR(Temp3, ONE);
	FREEMPR(Z);
	FREEMPR(Temp2);
	FREEMPR(Temp3);
loop:
	Temp1 = X[i];
	X[i] = ADDR(X[i], ONE);
	FREEMPR(Temp1);
	if (COMPAREI(X[i]->N, L[i]->N) == 1)
	{
		i++;
		if (i == n)
			goto exiting;
		else
			goto loop;
	}
	else
	{
		if (i)
		{
			Temp1 = ADDR(X[i], U[i]);
			Temp0 = SUBR(Temp1, N[i]);
			FREEMPR(Temp1);
			Temp2 = MULTR(Temp0, Temp0);
			FREEMPR(Temp0);
			Temp1 = MULTR(elt(Q, i, i), Temp2);
			FREEMPR(Temp2);
			if (t[i - 1])
				FREEMPR(T[i - 1]);
			else
				t[i - 1] = 1;
			T[i - 1] = SUBR(T[i], Temp1);
			FREEMPR(Temp1);
			i--;
			SUM = ZEROR();
			for (j = i + 1; j < n; j++)
			{
				Temp1 = SUM;
				Temp2 = MULTR(elt(Q, i, j), X[j]);
				SUM = ADDR(SUM, Temp2);
				FREEMPR(Temp1);
				FREEMPR(Temp2);
			}
			if (u[i])
				FREEMPR(U[i]);
			else
				u[i] = 1;
			U[i] = COPYR(SUM);
			FREEMPR(SUM);
			goto bounds;
		}
		else
			goto found;
	}

found:
	VECTOR = BUILDMATR(1, AA->R);
	XX = (MPR **)mmalloc((USL)(n * sizeof(MPR *)));
	for (j = 0; j < n; j++)
	{
		tempR = BUILDMPR();
		tempR->N = COPYI(I[j]);
		tempR->D = ONEI();
		XX[j] = ADDR(X[j], tempR);
		FREEMPR(tempR);
	}
	if (number)
	{
		printf("\tMULTIPLIER[%u] = ", count + 1);
		fprintf(outfile, "\tMULTIPLIER[%u]:", count + 1);
		printf("b[%u]", n + 1);
		fprintf(outfile, "b[%u]", n + 1);
		for (j = 0; j < n; j++)
		{
			e = (XX[j]->N)->S;
			if (e == -1)
			{
				printf("+");
				fprintf(outfile, "+");
				temp = MINUSR(XX[j]);
				if (!EQONER(temp))
				{
					PRINTR(temp);
					FPRINTR(outfile, temp);
				}
				printf("b[%u]", j + 1);
				fprintf(outfile, "b[%u]", j + 1);
				FREEMPR(temp);
			}
			if (e == 1)
			{
				printf("-");
				fprintf(outfile, "-");
				if (!EQONER(XX[j]))
				{
					PRINTR(XX[j]);
					FPRINTR(outfile, XX[j]);
				}
				printf("b[%u]", j + 1);
				fprintf(outfile, "b[%u]", j + 1);
			}
		}
	}
	for (j = 0; j < n; j++)
	{
		elt(VECTOR, 0, j) = MINUSR(X[j]);
		FREEMPR(XX[j]);
	}
	ffree((char *)XX, n * sizeof(MPR *));
	elt(VECTOR, 0, n) = ONER();
	MATR = MULTMATR(VECTOR, A);
	FREEMATR(VECTOR);
	if (number)
	{
		printf("=");
		fprintf(outfile, "=");
	}
	MATI = BUILDMATI(1, m);
	for (j = 0; j < m; j++)
		MATI->V[0][j] = COPYI((elt(MATR, 0, j))->N);
	if (number)
	{
		for (j = 0; j < m; j++)
		{
			PRINTI(MATI->V[0][j]);
			FPRINTI(outfile, MATI->V[0][j]);
			printf(" ");
			fprintf(outfile, " ");
		}
		fprintf(outfile, ": ");
		printf(": ");
	}
	lengthi = LENGTHSQRI(MATI, 0);
	if (number)
	{

		FPRINTI(outfile, lengthi);
		fprintf(outfile, "\n");
		PRINTI(lengthi);
		printf("\n");
	}
	FREEMPI(lengthi);
	FREEMATR(MATR);
	FREEMATI(MATI);
	count++;
	goto loop;
exiting:
	FREEMPR(ONE);
	for (j = 0; j < n; j++)
	{
		FREEMPR(N[j]);
		FREEMPR(X[j]);
		FREEMPR(L[j]);
		FREEMPR(T[j]);
		FREEMPR(U[j]);
		FREEMPR(Y[j]);
		FREEMPR(V[j]);
	}
	ffree((char *)N, n * sizeof(MPR *));
	ffree((char *)V, n * sizeof(MPR *));
	ffree((char *)Y, n * sizeof(MPR *));
	ffree((char *)X, n * sizeof(MPR *));
	ffree((char *)L, n * sizeof(MPR *));
	ffree((char *)T, n * sizeof(MPR *));
	ffree((char *)U, n * sizeof(MPR *));
	ffree((char *)l, n * sizeof(USI));
	ffree((char *)t, n * sizeof(USI));
	ffree((char *)u, n * sizeof(USI));
	ffree((char *)x, n * sizeof(USI));
	FREEMATR(Q);
	FREEMATR(A);
	if (number)
	{
		if (count > 1)
		{
			fprintf(outfile, "\tThere are exactly %u shortest multiplier vectors\n", count);
			printf("\tThere are exactly %u shortest multiplier vectors\n", count);
		}
		if (count == 1)
		{
			printf("\tThere is exactly one shortest multiplier vector\n");
			fprintf(outfile, "\tThere is exactly one shortest multiplier vector\n");
		}
	}
	fclose(outfile);
	return;
}

void SHORTESTTT(MPMATI *AA, MPI *BOUND)
/*
* This is used in INHOMFP().
* Input: An m x M matrix of integers AA whose first m - 1 rows Q[0],... * ,Q[m-2] are LI and span a lattice L.
* Output: If P is the last row, finds all x[0],...,x[m-2]
* satisfying ||P -x[0]Q[0]-...-x[m-2]Q[m-2]||^2<=BOUND.
*/
{
	unsigned int i, j, k, m, n, *l, *u, *t, *x, w, count = 0;
	unsigned int M, FLAG;
	MPR **X, **Y, **L, **T, **U, **V, *Temp0, *Temp1, *Temp2, *Temp3;
	MPR *C, *Z, *ONE, *SUM, **N, *temp;
	MPI *lengthi, *lengthj, *TempI1, *TempI2;
	MPMATR *A, *B, *BB, *Q, *QQ, *MATR, *VECTOR;
	MPMATI *MATI;
	int e;
	char buff[25];
	FILE *outfile;

	strcpy(buff, "inhomfp.out");
	outfile = fopen(buff, "w");
	/* Compute ||AA[m - 1][*]||^2. */
	m = AA->R;
	M = AA->C;
	n = m - 1;
	lengthj = ZEROI();
	for (j = 0; j < M; j++)
	{
		TempI1 = lengthj;
		TempI2 = MULTI(AA->V[n][j], AA->V[n][j]);
		lengthj = ADDI(lengthj, TempI2);
		FREEMPI(TempI1);
		FREEMPI(TempI2);
	}
	/* First get the vector [V[0],...,V[n-2], where V[i]=<A[n-1][*],A[i][*]>.*/
	A = BUILDMATR(m, M);
	for (i = 0; i < m; i++)
		for (j = 0; j < M; j++)
		{
			elt(A, i, j) = BUILDMPR();
			elt(A, i, j)->N = COPYI(AA->V[i][j]);
			elt(A, i, j)->D = ONEI();
		}

	V = (MPR **)mmalloc(n * sizeof(MPR *));
	for (i = 0; i < n; i++)
	{
		V[i] = ZEROR();
		for (j = 0; j < M; j++)
		{
			Temp1 = V[i];
			Temp2 = MULTR(elt(A, n, j), elt(A, i, j));
			V[i] = ADDR(V[i], Temp2);
			FREEMPR(Temp1);
			FREEMPR(Temp2);
		}
	}
	B = TRANSPOSER(A);
	BB = MULTMATR(A, B);
	FREEMATR(B);
	Q = CHOLESKYR(BB);
	FREEMATR(BB);
	QQ = COPYMATR(Q);
	for (i = 0; i < n; i++)
	{
		FREEMPR(elt(QQ, i, i));
		elt(QQ, i, i) = ONER();
	}
	/* Calculate (U^{-1})^tV = Y */
	MATR = QQ;
	QQ = TRANSPOSER(QQ);
	FREEMATR(MATR);
	Y = (MPR **)mmalloc(n * sizeof(MPR *));
	Y[0] = COPYR(V[0]);
	for (w = 1; w < n; w++)
	{
		SUM = ZEROR();
		for (k = 0; k < w; k++)
		{
			Temp1 = SUM;
			Temp2 = MULTR(elt(QQ, w, k), Y[k]);
			SUM = ADDR(SUM, Temp2);
			FREEMPR(Temp1);
			FREEMPR(Temp2);
		}
		Y[w] = SUBR(V[w], SUM);
		FREEMPR(SUM);
	}
	FREEMATR(QQ);
	N = (MPR **)mmalloc(n * sizeof(MPR *));
	for (j = 0; j < n; j++)
		N[j] = RATIOR(Y[j], elt(Q, j, j));
	ONE = ONER();
	X = (MPR **)mmalloc(n * sizeof(MPR *));
	L = (MPR **)mmalloc(n * sizeof(MPR *));
	T = (MPR **)mmalloc(n * sizeof(MPR *));
	U = (MPR **)mmalloc(n * sizeof(MPR *));
	l = (USI *)mmalloc(n * sizeof(USI));
	t = (USI *)mmalloc(n * sizeof(USI));
	u = (USI *)mmalloc(n * sizeof(USI));
	x = (USI *)mmalloc(n * sizeof(USI));
	for (j = 0; j < n; j++)
	{
		l[j] = 0;
		t[j] = 0;
		u[j] = 0;
		x[j] = 0;
	}
	C = ZEROR();
	for (i = 0; i < n; i++)
	{
		Temp1 = C;
		Temp2 = MULTR(Y[i], Y[i]);
		Temp3 = RATIOR(Temp2, elt(Q, i, i));
		C = ADDR(C, Temp3);
		FREEMPR(Temp1);
		FREEMPR(Temp2);
		FREEMPR(Temp3);
	}

	Temp0 = BUILDMPR();
	Temp0->N = SUBI(BOUND, lengthj);
	Temp0->D = ONEI();
	Temp1 = C;
	C = ADDR(C, Temp0);
	FREEMPR(Temp0);
	FREEMPR(Temp1);

	i = n - 1;
	T[i] = COPYR(C); t[i] = 1;
	FREEMPR(C);
	U[i] = ZEROR(); u[i] = 1;
bounds:
	Z = RATIOR(T[i], elt(Q, i, i));
	Temp2 = SUBR(N[i], U[i]);
	if (l[i])
		FREEMPR(L[i]);
	else
		l[i] = 1;
	L[i] = INTROOT(Z, Temp2);
	/*printf("at L[%u]\n", i);*/
	FREEMPR(Temp2);
	Temp0 = SUBR(U[i], N[i]);
	Temp2 = INTROOT(Z, Temp0);
	FREEMPR(Temp0);
	Temp3 = MINUSR(Temp2);
	if (x[i])
		FREEMPR(X[i]);
	else
		x[i] = 1;
	X[i] = SUBR(Temp3, ONE);
	FREEMPR(Z);
	FREEMPR(Temp2);
	FREEMPR(Temp3);
loop:
	Temp1 = X[i];
	X[i] = ADDR(X[i], ONE);
	FREEMPR(Temp1);
	if (COMPAREI(X[i]->N, L[i]->N) == 1)
	{
		i++;
		if (i == n)
			goto exiting;
		else
			goto loop;
	}
	else
	{
		if (i)
		{
			Temp1 = ADDR(X[i], U[i]);
			Temp0 = SUBR(Temp1, N[i]);
			FREEMPR(Temp1);
			Temp2 = MULTR(Temp0, Temp0);
			FREEMPR(Temp0);
			Temp1 = MULTR(elt(Q, i, i), Temp2);
			FREEMPR(Temp2);
			if (t[i - 1])
				FREEMPR(T[i - 1]);
			else
				t[i - 1] = 1;
			T[i - 1] = SUBR(T[i], Temp1);
			FREEMPR(Temp1);
			i--;
			SUM = ZEROR();
			for (j = i + 1; j < n; j++)
			{
				Temp1 = SUM;
				Temp2 = MULTR(elt(Q, i, j), X[j]);
				SUM = ADDR(SUM, Temp2);
				FREEMPR(Temp1);
				FREEMPR(Temp2);
			}
			if (u[i])
				FREEMPR(U[i]);
			else
				u[i] = 1;
			U[i] = COPYR(SUM);
			FREEMPR(SUM);
			goto bounds;
		}
		else
			goto found;
	}

found:
	VECTOR = BUILDMATR(1, AA->R);
	for (j = 0; j < n; j++)
	{
		/*
		printf("X[%u][%u] = ", count, j);
		PRINTR(X[j]);
		printf("\n");
		fprintf(outfile, "X[%u][%u] = ", count, j);
		FPRINTR(outfile, X[j]);
		fprintf(outfile,"\n");
		*/
		elt(VECTOR, 0, j) = COPYR(X[j]);
	}
	elt(VECTOR, 0, n) = ZEROR();
	MATR = MULTMATR(VECTOR, A);
	FREEMATR(VECTOR);
	printf("LV[%u] = ", count + 1);
	fprintf(outfile, "LV[%u] = ", count + 1);
	FLAG = 0;
	for (j = 0; j < n; j++)
	{
		e = (X[j]->N)->S;
		if (e == -1)
		{
			if (FLAG == 0)
				FLAG = 1;
			printf("-");
			fprintf(outfile, "-");
			temp = MINUSR(X[j]);
			if (!EQONER(temp))
			{
				PRINTR(temp);
				FPRINTR(outfile, temp);
			}
			printf("b[%u]", j + 1);
			fprintf(outfile, "b[%u]", j + 1);
			FREEMPR(temp);
		}
		if (e == 1)
		{
			if (FLAG)
			{
				printf("+");
				fprintf(outfile, "+");
			}
			if (FLAG == 0)
				FLAG = 1;
			if (!EQONER(X[j]))
			{
				PRINTR(X[j]);
				FPRINTR(outfile, X[j]);
			}
			printf("b[%u]", j + 1);
			fprintf(outfile, "b[%u]", j + 1);
		}
	}
	printf("=");
	fprintf(outfile, "=");
	MATI = BUILDMATI(1, A->C);
	for (j = 0; j < A->C; j++)
	{
		MATI->V[0][j] = COPYI((elt(MATR, 0, j))->N);
		PRINTI(MATI->V[0][j]);
		FPRINTI(outfile, MATI->V[0][j]);
		printf(" ");
		fprintf(outfile, " ");
		FREEMPI(MATI->V[0][j]);
		MATI->V[0][j] = SUBI(AA->V[n][j], ((elt(MATR, 0, j))->N));
	}
	printf(": ||b[%u]-LV[%u]||^2=", AA->R, count + 1);
	fprintf(outfile, ": ||b[%u]-LV[%u]||^2=", AA->R, count + 1);
	lengthi = LENGTHSQRI(MATI, 0);
	PRINTI(lengthi);
	FPRINTI(outfile, lengthi);
	printf("\n");
	fprintf(outfile, "\n");
	FREEMPI(lengthi);
	FREEMATR(MATR);
	FREEMATI(MATI);
	count++;
	goto loop;
exiting:
	FREEMPR(ONE);
	for (j = 0; j < n; j++)
	{
		FREEMPR(N[j]);
		FREEMPR(X[j]);
		FREEMPR(L[j]);
		FREEMPR(T[j]);
		FREEMPR(U[j]);
		FREEMPR(Y[j]);
		FREEMPR(V[j]);
	}
	ffree((char *)N, n * sizeof(MPR *));
	ffree((char *)V, n * sizeof(MPR *));
	ffree((char *)Y, n * sizeof(MPR *));
	ffree((char *)X, n * sizeof(MPR *));
	ffree((char *)L, n * sizeof(MPR *));
	ffree((char *)T, n * sizeof(MPR *));
	ffree((char *)U, n * sizeof(MPR *));
	ffree((char *)l, n * sizeof(USI));
	ffree((char *)t, n * sizeof(USI));
	ffree((char *)u, n * sizeof(USI));
	ffree((char *)x, n * sizeof(USI));
	fclose(outfile);
	FREEMATR(Q);
	FREEMATR(A);
	FREEMPI(lengthj);
	if (count > 1)
		printf("There are exactly %u lattice vectors\n", count);
	if (count == 1)
		printf("There is exactly one lattice vector\n");
	return;
}

MPMATI *SHORTESTT0(MPMATI *AA, MPI **XX[])
/*
* Used in EXTGCD(), LLLGCD() and LLLGCD0() in nfunc.c.
* Input: An m x M matrix of integers AA whose first m - 1 rows
* Q[0],...,Q[m-2] are LI and span a lattice L.
* Problem: If P is the last row, find a shorter multiplier
* P-x[0]Q[0]-...-x[m-2]Q[m-2]
*/
{
	unsigned int i, j, k, m, n, *l, *u, *t, *x, w;
	unsigned int FOUND = 0, M;
	MPR **X, **Y, **L, **T, **U, **V, *Temp0, *Temp1, *Temp2, *Temp3;
	MPR *C, *Z, *ONE, *SUM, **N;
	MPI *lengthi, *lengthj, *TempI1, *TempI2;
	MPMATR *A, *B, *BB, *Q, *QQ, *MATR, *VECTOR;
	MPMATI *MATI;
	int tt;

	MATI = NULL;
	/* Compute ||AA[m - 1][*]||^2. */
	m = AA->R;
	M = AA->C;
	n = m - 1;
	*XX = (MPI **)mmalloc((USL)(n * sizeof(MPI *)));
	lengthj = ZEROI();
	for (j = 0; j < M; j++)
	{
		TempI1 = lengthj;
		TempI2 = MULTI(AA->V[n][j], AA->V[n][j]);
		lengthj = ADDI(lengthj, TempI2);
		FREEMPI(TempI1);
		FREEMPI(TempI2);
	}
	/* First get the vector [V[0],...,V[n-2], where V[i]=<A[n-1][*],A[i][*]>.*/
	A = BUILDMATR(m, M);
	for (i = 0; i < m; i++)
		for (j = 0; j < M; j++)
		{
			elt(A, i, j) = BUILDMPR();
			elt(A, i, j)->N = COPYI(AA->V[i][j]);
			elt(A, i, j)->D = ONEI();
		}

	V = (MPR **)mmalloc(n * sizeof(MPR *));
	for (i = 0; i < n; i++)
	{
		V[i] = ZEROR();
		for (j = 0; j < M; j++)
		{
			Temp1 = V[i];
			Temp2 = MULTR(elt(A, n, j), elt(A, i, j));
			V[i] = ADDR(V[i], Temp2);
			FREEMPR(Temp1);
			FREEMPR(Temp2);
		}
	}
	B = TRANSPOSER(A);
	BB = MULTMATR(A, B);
	FREEMATR(B);
	Q = CHOLESKYR(BB);
	FREEMATR(BB);
	QQ = COPYMATR(Q);
	for (i = 0; i < n; i++)
	{
		FREEMPR(elt(QQ, i, i));
		elt(QQ, i, i) = ONER();
	}
	/* Calculate (U^{-1})^tV = Y */
	MATR = QQ;
	QQ = TRANSPOSER(QQ);
	FREEMATR(MATR);
	Y = (MPR **)mmalloc(n * sizeof(MPR *));
	Y[0] = COPYR(V[0]);
	for (w = 1; w < n; w++)
	{
		SUM = ZEROR();
		for (k = 0; k < w; k++)
		{
			Temp1 = SUM;
			Temp2 = MULTR(elt(QQ, w, k), Y[k]);
			SUM = ADDR(SUM, Temp2);
			FREEMPR(Temp1);
			FREEMPR(Temp2);
		}
		Y[w] = SUBR(V[w], SUM);
		FREEMPR(SUM);
	}
	FREEMATR(QQ);
	N = (MPR **)mmalloc(n * sizeof(MPR *));
	for (j = 0; j < n; j++) {
		N[j] = RATIOR(Y[j], elt(Q, j, j));
	}
	ONE = ONER();
	X = (MPR **)mmalloc(n * sizeof(MPR *));
	L = (MPR **)mmalloc(n * sizeof(MPR *));
	T = (MPR **)mmalloc(n * sizeof(MPR *));
	U = (MPR **)mmalloc(n * sizeof(MPR *));
	l = (USI *)mmalloc(n * sizeof(USI));
	t = (USI *)mmalloc(n * sizeof(USI));
	u = (USI *)mmalloc(n * sizeof(USI));
	x = (USI *)mmalloc(n * sizeof(USI));
	for (j = 0; j < n; j++)
	{
		l[j] = 0;
		t[j] = 0;
		u[j] = 0;
		x[j] = 0;
	}
	C = ZEROR();
	for (i = 0; i < n; i++)
	{
		Temp1 = C;
		Temp2 = MULTR(Y[i], Y[i]);
		Temp3 = RATIOR(Temp2, elt(Q, i, i));
		C = ADDR(C, Temp3);
		FREEMPR(Temp1);
		FREEMPR(Temp2);
		FREEMPR(Temp3);
	}
	i = n - 1;
	T[i] = COPYR(C); t[i] = 1;
	FREEMPR(C);
	U[i] = ZEROR(); u[i] = 1;
bounds:
	Z = RATIOR(T[i], elt(Q, i, i));
	Temp2 = SUBR(N[i], U[i]);
	if (l[i])
		FREEMPR(L[i]);
	else
		l[i] = 1;
	L[i] = INTROOT(Z, Temp2);
	FREEMPR(Temp2);
	Temp0 = SUBR(U[i], N[i]);
	Temp2 = INTROOT(Z, Temp0);
	FREEMPR(Temp0);
	Temp3 = MINUSR(Temp2);
	if (x[i])
		FREEMPR(X[i]);
	else
		x[i] = 1;
	X[i] = SUBR(Temp3, ONE);
	FREEMPR(Z);
	FREEMPR(Temp2);
	FREEMPR(Temp3);
loop:
	Temp1 = X[i];
	X[i] = ADDR(X[i], ONE);
	FREEMPR(Temp1);
	if (COMPAREI(X[i]->N, L[i]->N) == 1)
	{
		i++;
		if (i == n) {
			goto exiting;
		}
		else {
			goto loop;
		}
	}
	else
	{
		if (i)
		{
			Temp1 = ADDR(X[i], U[i]);
			Temp0 = SUBR(Temp1, N[i]);
			FREEMPR(Temp1);
			Temp2 = MULTR(Temp0, Temp0);
			FREEMPR(Temp0);
			Temp1 = MULTR(elt(Q, i, i), Temp2);
			FREEMPR(Temp2);
			if (t[i - 1])
				FREEMPR(T[i - 1]);
			else
				t[i - 1] = 1;
			T[i - 1] = SUBR(T[i], Temp1);
			FREEMPR(Temp1);
			i--;
			SUM = ZEROR();
			for (j = i + 1; j < n; j++)
			{
				Temp1 = SUM;
				Temp2 = MULTR(elt(Q, i, j), X[j]);
				SUM = ADDR(SUM, Temp2);
				FREEMPR(Temp1);
				FREEMPR(Temp2);
			}
			if (u[i])
				FREEMPR(U[i]);
			else
				u[i] = 1;
			U[i] = COPYR(SUM);
			FREEMPR(SUM);
			goto bounds;
		}
		else
			goto found;
	}

found:
	VECTOR = BUILDMATR(1, m);
	for (j = 0; j < n; j++)
		elt(VECTOR, 0, j) = MINUSR(X[j]);

	elt(VECTOR, 0, m - 1) = ONER();
	MATR = MULTMATR(VECTOR, A);
	FREEMATR(VECTOR);
	MATI = BUILDMATI(1, M);
	for (j = 0; j < M; j++)
		MATI->V[0][j] = COPYI((elt(MATR, 0, j))->N);
	lengthi = LENGTHSQRI(MATI, 0);
	tt = RSV(lengthj, lengthi);
	FREEMPI(lengthi);
	FREEMATR(MATR);
	if (tt == 1)
	{
		FREEMPI(lengthj);
		FOUND = 1;
		for (j = 0; j < n; j++)
			(*XX)[j] = COPYI(X[j]->N);
		goto exiting;
	}
	else
		FREEMATI(MATI);
	goto loop;
exiting:
	FREEMPR(ONE);
	for (j = 0; j < n; j++)
	{
		FREEMPR(N[j]);
		FREEMPR(X[j]);
		FREEMPR(L[j]);
		FREEMPR(T[j]);
		FREEMPR(U[j]);
		FREEMPR(Y[j]);
		FREEMPR(V[j]);
	}
	ffree((char *)N, n * sizeof(MPR *));
	ffree((char *)V, n * sizeof(MPR *));
	ffree((char *)Y, n * sizeof(MPR *));
	ffree((char *)X, n * sizeof(MPR *));
	ffree((char *)L, n * sizeof(MPR *));
	ffree((char *)T, n * sizeof(MPR *));
	ffree((char *)U, n * sizeof(MPR *));
	ffree((char *)l, n * sizeof(USI));
	ffree((char *)t, n * sizeof(USI));
	ffree((char *)u, n * sizeof(USI));
	ffree((char *)x, n * sizeof(USI));
	FREEMATR(Q);
	FREEMATR(A);
	if (FOUND)
	{
		if (GCDVERBOSE)
		{
			printf("found a shorter multiplier:\n");
			PRINTMATI(0, MATI->R - 1, 0, MATI->C - 1, MATI);
			lengthi = LENGTHSQRI(MATI, 0);
			printf("LENGTHSQUARED = ");
			PRINTI(lengthi);
			FREEMPI(lengthi);
			printf("\n");
			printf("search continuing:\n\n");
		}
		return(MATI);
	}
	else
	{
		FREEMPI(lengthj);
		if (GCDVERBOSE)
			printf("search ended.\n");
		for (j = 0; j < n; j++)
			(*XX)[j] = ZEROI();
		return(NULL);
	}
}

MPI ***SHORTESTX(MPMATI *AA, MPI *I[], USI *counter)
/*
* Input: An m x m matrix of integers AA = Q arising from LLLGCD.
* The last row of AA is assumed to be a multiplier of shortest length.
* Output: the XX[] defining the multipliers of shortest length.
* Only used in GCD_CONJ in nfunc.c.
* Warning: Array XX below only caters for up to GCD_MAX=100000 shortest multipliers.
* If necessary, increase it and recompile.
*/
{
	unsigned int i, j, k, m, n, *l, *u, *t, *x, w, count = 0;
	MPR **X, **Y, **L, **T, **U, **V, *Temp0, *Temp1, *Temp2, *Temp3;
	MPR *C, *Z, *ONE, *SUM, **N;
	MPI ***XX, *tempI;
	MPMATR *A, *B, *BB, *Q, *QQ, *MATR;

	XX = (MPI ***)mmalloc((USL)(GCD_MAX * sizeof(MPI **)));
	/* First get the vector [V[0],...,V[n-2], where V[i]=<A[n-1][*],A[i][*]>.*/
	m = AA->C;
	n = (AA->R) - 1;
	A = BUILDMATR(n + 1, m);
	for (i = 0; i < n + 1; i++)
		for (j = 0; j < m; j++)
		{
			elt(A, i, j) = BUILDMPR();
			elt(A, i, j)->N = COPYI(AA->V[i][j]);
			elt(A, i, j)->D = ONEI();
		}

	V = (MPR **)mmalloc(n * sizeof(MPR *));
	for (i = 0; i < n; i++)
	{
		V[i] = ZEROR();
		for (j = 0; j < m; j++)
		{
			Temp1 = V[i];
			Temp2 = MULTR(elt(A, n, j), elt(A, i, j));
			V[i] = ADDR(V[i], Temp2);
			FREEMPR(Temp1);
			FREEMPR(Temp2);
		}
	}
	B = TRANSPOSER(A);
	BB = MULTMATR(A, B);
	FREEMATR(B);
	Q = CHOLESKYR(BB);
	FREEMATR(BB);
	QQ = COPYMATR(Q);
	for (i = 0; i < n; i++)
	{
		FREEMPR(elt(QQ, i, i));
		elt(QQ, i, i) = ONER();
	}
	/* Calculate (U^{-1})^tV = Y */
	MATR = QQ;
	QQ = TRANSPOSER(QQ);
	FREEMATR(MATR);
	Y = (MPR **)mmalloc(n * sizeof(MPR *));
	Y[0] = COPYR(V[0]);
	for (w = 1; w < n; w++)
	{
		SUM = ZEROR();
		for (k = 0; k < w; k++)
		{
			Temp1 = SUM;
			Temp2 = MULTR(elt(QQ, w, k), Y[k]);
			SUM = ADDR(SUM, Temp2);
			FREEMPR(Temp1);
			FREEMPR(Temp2);
		}
		Y[w] = SUBR(V[w], SUM);
		FREEMPR(SUM);
	}
	FREEMATR(QQ);
	N = (MPR **)mmalloc(n * sizeof(MPR *));
	for (j = 0; j < n; j++)
		N[j] = RATIOR(Y[j], elt(Q, j, j));
	ONE = ONER();
	X = (MPR **)mmalloc(n * sizeof(MPR *));
	L = (MPR **)mmalloc(n * sizeof(MPR *));
	T = (MPR **)mmalloc(n * sizeof(MPR *));
	U = (MPR **)mmalloc(n * sizeof(MPR *));
	l = (USI *)mmalloc(n * sizeof(USI));
	t = (USI *)mmalloc(n * sizeof(USI));
	u = (USI *)mmalloc(n * sizeof(USI));
	x = (USI *)mmalloc(n * sizeof(USI));
	for (j = 0; j < n; j++)
	{
		l[j] = 0;
		t[j] = 0;
		u[j] = 0;
		x[j] = 0;
	}
	C = ZEROR();
	for (i = 0; i < n; i++)
	{
		Temp1 = C;
		Temp2 = MULTR(Y[i], Y[i]);
		Temp3 = RATIOR(Temp2, elt(Q, i, i));
		C = ADDR(C, Temp3);
		FREEMPR(Temp1);
		FREEMPR(Temp2);
		FREEMPR(Temp3);
	}
	i = n - 1;
	T[i] = COPYR(C); t[i] = 1;
	FREEMPR(C);
	U[i] = ZEROR(); u[i] = 1;
bounds:
	Z = RATIOR(T[i], elt(Q, i, i));
	Temp2 = SUBR(N[i], U[i]);
	if (l[i])
		FREEMPR(L[i]);
	else
		l[i] = 1;
	L[i] = INTROOT(Z, Temp2);
	FREEMPR(Temp2);
	Temp0 = SUBR(U[i], N[i]);
	Temp2 = INTROOT(Z, Temp0);
	FREEMPR(Temp0);
	Temp3 = MINUSR(Temp2);
	if (x[i])
		FREEMPR(X[i]);
	else
		x[i] = 1;
	X[i] = SUBR(Temp3, ONE);
	FREEMPR(Z);
	FREEMPR(Temp2);
	FREEMPR(Temp3);
loop:
	Temp1 = X[i];
	X[i] = ADDR(X[i], ONE);
	FREEMPR(Temp1);
	if (COMPAREI(X[i]->N, L[i]->N) == 1)
	{
		i++;
		if (i == n)
			goto exiting;
		else
			goto loop;
	}
	else
	{
		if (i)
		{
			Temp1 = ADDR(X[i], U[i]);
			Temp0 = SUBR(Temp1, N[i]);
			FREEMPR(Temp1);
			Temp2 = MULTR(Temp0, Temp0);
			FREEMPR(Temp0);
			Temp1 = MULTR(elt(Q, i, i), Temp2);
			FREEMPR(Temp2);
			if (t[i - 1])
				FREEMPR(T[i - 1]);
			else
				t[i - 1] = 1;
			T[i - 1] = SUBR(T[i], Temp1);
			FREEMPR(Temp1);
			i--;
			SUM = ZEROR();
			for (j = i + 1; j < n; j++)
			{
				Temp1 = SUM;
				Temp2 = MULTR(elt(Q, i, j), X[j]);
				SUM = ADDR(SUM, Temp2);
				FREEMPR(Temp1);
				FREEMPR(Temp2);
			}
			if (u[i])
				FREEMPR(U[i]);
			else
				u[i] = 1;
			U[i] = COPYR(SUM);
			FREEMPR(SUM);
			goto bounds;
		}
		else
			goto found;
	}

found:
	XX[count] = (MPI **)mmalloc((USL)(n * sizeof(MPI *)));
	for (j = 0; j < n; j++)
	{
		tempI = ADDI(X[j]->N, I[j]);
		XX[count][j] = MINUSI(tempI);
		FREEMPI(tempI);
	}
	count++;
	if (count > GCD_MAX)
	{
		fprintf(stderr, "count > %u\n", GCD_MAX);
		exit(1);
	}
	goto loop;
exiting:
	FREEMPR(ONE);
	for (j = 0; j < n; j++)
	{
		FREEMPR(N[j]);
		FREEMPR(X[j]);
		FREEMPR(L[j]);
		FREEMPR(T[j]);
		FREEMPR(U[j]);
		FREEMPR(Y[j]);
		FREEMPR(V[j]);
	}
	ffree((char *)N, n * sizeof(MPR *));
	ffree((char *)V, n * sizeof(MPR *));
	ffree((char *)Y, n * sizeof(MPR *));
	ffree((char *)X, n * sizeof(MPR *));
	ffree((char *)L, n * sizeof(MPR *));
	ffree((char *)T, n * sizeof(MPR *));
	ffree((char *)U, n * sizeof(MPR *));
	ffree((char *)l, n * sizeof(USI));
	ffree((char *)t, n * sizeof(USI));
	ffree((char *)u, n * sizeof(USI));
	ffree((char *)x, n * sizeof(USI));
	FREEMATR(Q);
	FREEMATR(A);
	/*
	if (count > 1)
	printf("There are exactly %u shortest multiplier vectors\n", count);
	if (count == 1)
	printf("There is exactly one shortest multiplier vector\n");
	*/
	*counter = count;
	return (XX);
}
