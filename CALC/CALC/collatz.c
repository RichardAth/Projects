/* collatz.c */
#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_NONSTDC_NO_DEPRECATE
#ifdef _WIN32
#include "unistd_DOS.h"
#else
#include <unistd.h>
#endif
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "integer.h"
#include "fun.h"
#define UPPER_BOUND 4294967294UL /* 2^32 - 2 */
#define CYCLE_MAX 150
#define BRANCH_MAX 20
#define ALG_BRANCH_MAX 20

MPI *GEN_COLLATZ(MPI *Aptr, USL d, MPI *m[], MPI *X[])
/*
* The generalized Collatz function.
* Here d is the number of branches.
* m[0],...,m[d-1] are non-zero integers (the multipliers).
* X[0],...,X[d-1] are the shifts.
* T(Aptr)=INT(m[r]*Aptr/d)+X[r], if Aptr(MOD(d))=r.
*/

{
	USL r;
	MPI *T1, *T2, *T3;

	r = MOD_(Aptr, d);
	T1 = MULTI(Aptr, m[r]);
	T2 = INT_(T1, d);
	FREEMPI(T1);
	T3 = ADDI(T2, X[r]);
	FREEMPI(T2);
	return (T3);
}

void P_CYCLE(MPI *Aptr, USL d, MPI *m[], MPI *X[], USI c)

/* the cycle starting with *Aptr is printed */
{
	MPI *B, *Temp, *MIN_ELT, *TEMP1;
	unsigned int i = 0;
	FILE *outfile;
	char buff[20];

	strcpy(buff, "cycle.out");
	outfile = fopen(buff, "a");
	B = COPYI(Aptr);
	MIN_ELT = COPYI(Aptr);
	do {
		Temp = B;
		B = GEN_COLLATZ(B, d, m, X);
		FREEMPI(Temp);
		if (RSV(B, MIN_ELT) < 0) {
			TEMP1 = MIN_ELT;
			MIN_ELT = COPYI(B);
			FREEMPI(TEMP1);
		}
	} while (!EQUALI(Aptr, B));

	Temp = B;
	B = COPYI(MIN_ELT);
	FREEMPI(Temp);
	do {
		PRINTI(B);
		FPRINTI(outfile, B);
		printf("\n");
		fprintf(outfile, "\n");
		Temp = B;
		B = GEN_COLLATZ(B, d, m, X);
		FREEMPI(Temp);
		i++;
	} while (!EQUALI(MIN_ELT, B));
	printf("length of cycle %u is %u\n", c, i);
	fprintf(outfile, "length of cycle %u is %u\n", c, i);
	fclose(outfile);
	FREEMPI(B);
	FREEMPI(MIN_ELT);
	return;
}

unsigned int IS_CYCLE(MPI *Aptr, USL d, MPI *m[], MPI *X[], MPI *Z[], USI c)
/* the cycle starting from Aptr is distinct from those starting from
* Z[0],..., Z[c]  if and only if IS_CYCLE() = 1.
*/
{
	USL k;
	MPI *B, *Temp;
	unsigned int i, t;

	B = COPYI(Aptr);
	for (k = 0; 1; k++)
	{
		for (i = 0; i < c; i++)
		{
			if (EQUALI(B, Z[i]))
			{
				FREEMPI(B);
				return (0);
			}
		}
		Temp = B;
		B = GEN_COLLATZ(B, d, m, X);
		FREEMPI(Temp);
		t = EQUALI(Aptr, B);
		if (t) 			/* found a new cycle */
		{
			FREEMPI(B);
			break;
		}
	}
	return (1);
}

void CYCLE(USL d, MPI *m[], MPI *X[], USL infty, USL RANGE)
/*
* This function searches all trajectories of the generalized Collatz
* function which start from p, |p| <= RANGE/2 (RANGE an even integer).
*/
{
	MPI *A, *B, **Z, *Temp, *B1;
	unsigned int c = 0, n, r, end, *FLAG, t = 0, i, s = 0;
	USI diverge_flag = 0;
	long int p, q, LIMIT, k;
	FILE *outfile;
	char buff[20];

	strcpy(buff, "cycle.out");
	LIMIT = CYCLE_MAX;
	FLAG = (unsigned int *)ccalloc(RANGE + 1, sizeof(unsigned int));
	Z = (MPI **)mmalloc(LIMIT * sizeof(MPI *));
	INTSETI((int *)FLAG, RANGE + 1, (int)RANGE + 1);
	for (n = 0; n <= RANGE; n++) {
		if (FLAG[n] > n) {
			if ((n & 1) == 0)
				p = (n >> 1);
			else
				p = -((n + 1) >> 1);
			q = p;
			printf("p = %ld\n", p);
			A = CHANGEI(p);
			B = COPYI(A);
			end = 0;
			for (k = 0; k <= UPPER_BOUND; k++)
			{
				Temp = A;
				A = GEN_COLLATZ(A, d, m, X);
				FREEMPI(Temp);
				if (A->D >= infty) {
					if (diverge_flag == 0) {
						outfile = fopen(buff, "a");

						fprintf(outfile, "divergent trajectory starting from q = %ld\n", q);
						printf("divergent trajectory starting from q = %ld\n", q);
						fclose(outfile);
						diverge_flag = 0;
					}
					FREEMPI(A);
					FREEMPI(B);
					break;
				}
				Temp = B;
				B1 = GEN_COLLATZ(B, d, m, X);
				B = GEN_COLLATZ(B1, d, m, X);
				FREEMPI(B1);
				FREEMPI(Temp);
				t = EQUALI(A, B);
				if (t)
				{
					/* cycle detected starting at p */
					FREEMPI(B);
					break;
				}
				if (A->D == 0 && A->V[0] <= RANGE / 2)
				{
					p = A->V[0];
					r = (A->S >= 0 ? 2 * p : 2 * p - 1);
					if (FLAG[r] < n) {
						end = 1;
						FREEMPI(A);
						FREEMPI(B);
						break;
					}
					else
						FLAG[r] = n;
				}
			}
			if (end == 1)
				continue;
			if (t) {
				if (c == 0) {
					/* store A, the start of first cycle */
					Z[0] = A;
					printf("\ncycle %u:\n", c + 1);
					outfile = fopen(buff, "a");
					fprintf(outfile, "\ncycle %u:\n", c + 1);
					fclose(outfile);
					P_CYCLE(A, d, m, X, c + 1);
					c++;
				}
				else {
					/* store start of subsequent cycles */
					s = IS_CYCLE(A, d, m, X, Z, c);
					if (s == 0)
						FREEMPI(A);
					else
					{
						if (c >= LIMIT)
						{
							LIMIT++;
							/*		ptr = Z;*/
							Z = (MPI **)rrealloc((char *)Z, (c + 1) * sizeof(MPI *), sizeof(MPI *));
							/*if (ptr != Z)
							free((char *) ptr);*/
						}
						Z[c] = A;
						printf("\ncycle %u:\n", c + 1);
						outfile = fopen(buff, "a");
						fprintf(outfile, "\ncycle %u:\n", c + 1);
						fclose(outfile);
						P_CYCLE(A, d, m, X, c + 1);
						c++;
					}
				}
			}
		}
	}
	for (i = 0; i < c; i++)
		FREEMPI(Z[i]);
	ffree((char *)Z, LIMIT * sizeof(MPI *));
	ffree((char *)FLAG, (RANGE + 1) * sizeof(USI));
	return;
}

void CYCLEX()
{
	unsigned int i, u;
	int s;
	USL d, infty, RANGE;
	MPI **m, **X;
	FILE *outfile;
	char buff[20];

	strcpy(buff, "cycle.out");
	if (access(buff, R_OK) == 0)
		unlink(buff);
	outfile = fopen(buff, "w");
	printf("enter 'infty', (the cut-off number of allowable digits (base %u) for an iterate - say <20):", R0);
	s = scanf("%lu", &infty);
	printf("infty = %lu\n", infty);
	fprintf(outfile, "infty = %lu\n", infty);
	printf("enter an even integer 'RANGE' (<= (say)500,000): (searches all trajectories of the generalized Collatz function which start from p, |p| <= RANGE)/2\n");
	s = scanf("%lu", &RANGE);
	printf("RANGE = %lu\n", RANGE);
	fprintf(outfile, "RANGE = %lu\n", RANGE);
	printf("enter the number d of branches (<= %u)\n", BRANCH_MAX);
	s = scanf("%lu", &d);
	printf("the number of branches is d = %lu\n", d);
	fprintf(outfile, "the number of branches is d = %lu\n", d);

	m = (MPI **)mmalloc(d * sizeof(MPI *));
	X = (MPI **)mmalloc(d * sizeof(MPI *));
	for (i = 0; i < d; i++) {
		printf("enter the nonzero multiplier m[%u]\n", i);
		m[i] = INPUTI(&u);;
		printf("m[%u] = ", i);
		fprintf(outfile, "m[%u] = ", i);
		PRINTI(m[i]);
		FPRINTI(outfile, m[i]);
		printf("\n");
		fprintf(outfile, "\n");
	}

	for (i = 0; i < d; i++) {
		printf("enter the shift X[%u]\n", i);
		X[i] = INPUTI(&u);;
		printf("X[%u] = ", i);
		fprintf(outfile, "X[%u] = ", i);
		PRINTI(X[i]);
		FPRINTI(outfile, X[i]);
		printf("\n");
		fprintf(outfile, "\n");
	}
	fclose(outfile);
	Flush();
	printf("\n\n");
	CYCLE(d, m, X, infty, RANGE);
	for (i = 0; i < d; i++)
		FREEMPI(m[i]);
	ffree((char *)m, d * sizeof(MPI *));
	for (i = 0; i < d; i++)
		FREEMPI(X[i]);
	ffree((char *)X, d * sizeof(MPI *));
	return;
}

void ALG_COLLATZ(MPI *Aptr, MPI *Bptr, USL d, MPI *M[], MPI *N[], MPI *X[], MPI *Y[], MPI **AAptr, MPI **BBptr)
/*
* The generalized Collatz function mod sqrt(d).
* Here d is the number of branches.
* M[0]+N[0]sqrt(d),...,M[d-1]+N[d-1]sqrt(d) are non-zero integers of Z[sqrt(d)] (the multipliers).
* X[0]+Y[0]sqrt(d),...,X[d-1]+Y[d-1]sqrt(d) are the shifts.
* M[i]i-X[i] = 0 (mod d).
* T(Aptr+Bptr.sqrt(d))=(M[i]+N[i]sqrt(d))(Aptr+Bptr.sqrt(d))-(X[i]+Y[i]sqrt(d))/sqrt(d)if Aptr(MOD(d))=i.
*/

{
	USL i;
	MPI *T1, *T2, *T3, *T4, *T5;

	i = MOD_(Aptr, d);
	T1 = MULTI(Aptr, M[i]);
	T2 = MULT_I(Bptr, d);
	T3 = MULTI(N[i], T2);
	T4 = ADDI(T1, T3);
	T5 = SUBI(T4, X[i]);
	*BBptr = INT_(T5, d);
	FREEMPI(T1);
	FREEMPI(T2);
	FREEMPI(T3);
	FREEMPI(T4);
	FREEMPI(T5);
	T1 = MULTI(Aptr, N[i]);
	T2 = MULTI(Bptr, M[i]);
	T3 = ADDI(T1, T2);
	*AAptr = SUBI(T3, Y[i]);
	FREEMPI(T1);
	FREEMPI(T2);
	FREEMPI(T3);
	return;
}

void ALG_CYCLEX()
{
	unsigned int i, u;
	int s;
	USL d, infty, RANGE, temp;
	MPI **M, **N, **X, **Y, *TEMP1, *TEMP2;
	FILE *outfile;
	char buff[20];

	strcpy(buff, "alg_cycle.out");
	if (access(buff, R_OK) == 0)
		unlink(buff);
	outfile = fopen(buff, "w");
	printf("enter 'infty', (the cut-off number of allowable digits (base %u) for an iterate - say <20):", R0);
	s = scanf("%lu", &infty);
	printf("infty = %lu\n", infty);
	fprintf(outfile, "infty = %lu\n", infty);
	printf("enter a positive integer 'RANGE' (<= (say)500,000): (searches all trajectories of the generalized Collatz function which start from spiral(x[p],y[p]) p <= RANGE\n");
	s = scanf("%lu", &RANGE);
	printf("RANGE = %lu\n", RANGE);
	fprintf(outfile, "RANGE = %lu\n", RANGE);
	printf("enter the number d of branches (<= %u)\n", ALG_BRANCH_MAX);
	s = scanf("%lu", &d);
	printf("the number of branches is d = %lu\n", d);
	fprintf(outfile, "the number of branches is d = %lu\n", d);

	M = (MPI **)mmalloc(d * sizeof(MPI *));
	N = (MPI **)mmalloc(d * sizeof(MPI *));
	X = (MPI **)mmalloc(d * sizeof(MPI *));
	Y = (MPI **)mmalloc(d * sizeof(MPI *));
	for (i = 0; i < d; i++) {
		printf("enter the first component M[%u] of non-zero multiplier\n", i);
		M[i] = INPUTI(&u);;
		printf("M[%u] = ", i);
		fprintf(outfile, "M[%u] = ", i);
		PRINTI(M[i]);
		FPRINTI(outfile, M[i]);
		printf("\n");
		fprintf(outfile, "\n");
		printf("enter the second component N[%u] of non-zero multiplier\n", i);
		N[i] = INPUTI(&u);;
		printf("N[%u] = ", i);
		fprintf(outfile, "N[%u] = ", i);
		PRINTI(N[i]);
		FPRINTI(outfile, N[i]);
		printf("\n");
		fprintf(outfile, "\n");
	}

	for (i = 0; i < d; i++) {
		printf("enter the first component X[%u] of shift\n", i);
		X[i] = INPUTI(&u);;
		printf("X[%u] = ", i);
		fprintf(outfile, "X[%u] = ", i);
		PRINTI(X[i]);
		FPRINTI(outfile, X[i]);
		printf("\n");
		fprintf(outfile, "\n");
		TEMP1 = MULT_I(M[i], i);
		TEMP2 = SUBI(TEMP1, X[i]);
		FREEMPI(TEMP1);
		temp = MOD_(TEMP2, d);
		FREEMPI(TEMP2);
		if (temp) {
			printf("M[i]i-X[i] is not divisible by d when i=%u\n", i);
			return;
		}
		printf("enter the second component Y[%u] of shift\n", i);
		Y[i] = INPUTI(&u);;
		printf("Y[%u] = ", i);
		fprintf(outfile, "Y[%u] = ", i);
		PRINTI(Y[i]);
		FPRINTI(outfile, Y[i]);
		printf("\n");
		fprintf(outfile, "\n");
	}
	fclose(outfile);
	Flush();
	printf("\n\n");
	ALG_CYCLE(d, M, N, X, Y, infty, RANGE);
	for (i = 0; i < d; i++) {
		FREEMPI(M[i]);
		FREEMPI(N[i]);
	}
	ffree((char *)M, d * sizeof(MPI *));
	ffree((char *)N, d * sizeof(MPI *));
	for (i = 0; i < d; i++) {
		FREEMPI(X[i]);
		FREEMPI(Y[i]);
	}
	ffree((char *)X, d * sizeof(MPI *));
	ffree((char *)Y, d * sizeof(MPI *));
	return;
}

void ALG_P_CYCLE(MPI *Aptr, MPI *Bptr, USL d, MPI *M[], MPI *N[], MPI *X[], MPI *Y[], USI c)

/* the cycle starting with (*Aptr, *Bptr) is printed */
{
	MPI *A, *B, *AA, *BB, *TempA, *TempB, *MIN_ELTA, *MIN_ELTB, *TEMP1, *TEMP2;
	unsigned int i = 0;
	FILE *outfile;
	char buff[20];

	strcpy(buff, "alg_cycle.out");
	outfile = fopen(buff, "a");
	A = COPYI(Aptr);
	B = COPYI(Bptr);
	MIN_ELTA = COPYI(Aptr);
	MIN_ELTB = COPYI(Bptr);
	do {
		TempA = A;
		TempB = B;
		ALG_COLLATZ(A, B, d, M, N, X, Y, &AA, &BB);
		FREEMPI(TempA);
		FREEMPI(TempB);
		A = AA;
		B = BB;
		if (RSV(A, MIN_ELTA) < 0) {
			TEMP1 = MIN_ELTA;
			TEMP2 = MIN_ELTB;
			MIN_ELTA = COPYI(A);
			MIN_ELTB = COPYI(B);
			FREEMPI(TEMP1);
			FREEMPI(TEMP2);
		}
	} while (!EQUALI(Aptr, A) || !EQUALI(Bptr, B));

	TempA = A;
	A = COPYI(MIN_ELTA);
	FREEMPI(TempA);
	TempB = B;
	B = COPYI(MIN_ELTB);
	FREEMPI(TempB);
	do {
		printf("(");
		fprintf(outfile, "(");
		PRINTI(A);
		FPRINTI(outfile, A);
		printf(",");
		fprintf(outfile, ",");
		PRINTI(B);
		FPRINTI(outfile, B);
		printf(")\n");
		fprintf(outfile, ")\n");
		TempA = A;
		TempB = B;
		ALG_COLLATZ(A, B, d, M, N, X, Y, &AA, &BB);
		FREEMPI(TempA);
		FREEMPI(TempB);
		A = AA;
		B = BB;
		i++;
	} while (!EQUALI(MIN_ELTA, A) || !EQUALI(MIN_ELTB, B));
	printf("length of cycle %u is %u\n", c, i);
	fprintf(outfile, "length of cycle %u is %u\n", c, i);
	fclose(outfile);
	FREEMPI(A);
	FREEMPI(B);
	FREEMPI(MIN_ELTA);
	FREEMPI(MIN_ELTB);
	return;
}

unsigned int ALG_IS_CYCLE(MPI *Aptr, MPI *Bptr, USL d, MPI *M[], MPI *N[], MPI *X[], MPI *Y[], MPI *Z1[], MPI *Z2[], USI c)
/* the cycle starting from (Aptr, Bptr) is distinct from those starting from
* (Z1[0],Z2[0]),..., (Z1[c],Z2[c])  if and only if ALG_IS_CYCLE() = 1.
*/
{
	USL k;
	MPI *A, *B, *TempA, *TempB, *AA, *BB;
	unsigned int i, t;

	A = COPYI(Aptr);
	B = COPYI(Bptr);
	for (k = 0; 1; k++)
	{
		for (i = 0; i < c; i++)
		{
			if (EQUALI(A, Z1[i]) && EQUALI(B, Z2[i]))
			{
				FREEMPI(A);
				FREEMPI(B);
				return (0);
			}
		}
		TempA = A;
		TempB = B;
		ALG_COLLATZ(A, B, d, M, N, X, Y, &AA, &BB);
		FREEMPI(TempA);
		FREEMPI(TempB);
		A = AA;
		B = BB;
		t = (EQUALI(Aptr, A))*(EQUALI(Bptr, B));
		if (t) 			/* found a new cycle */
		{
			FREEMPI(A);
			FREEMPI(B);
			break;
		}
	}
	return (1);
}

void SPIRAL(MPI *N, MPI **X, MPI **Y) {
	USL r, s;
	MPI *L, *M, *R, *RR, *T, *TEMP, *ONE, *TWO;
	MPI *TEMP1, *TEMP2, *TEMP3, *TEMP4, *TEMP5, *TEMP6, *TEMP7;
	if (N->S) {
		M = BIG_MTHROOT(N, (USI)2);
		s = MOD_(M, 2);
	}
	else {
		M = ZEROI();
		s = 0;
	}
	TEMP = MULT_I(N, 4);
	if (N->S) {
		T = BIG_MTHROOT(TEMP, (USI)2);
		r = MOD_(T, 2);
		FREEMPI(T);
	}
	else {
		r = 0;
	}
	FREEMPI(TEMP);
	TWO = TWOI();
	L = CEILINGI(M, TWO);
	FREEMPI(TWO);
	ONE = ONEI();
	TEMP1 = ADD0I(M, ONE);
	FREEMPI(ONE);
	TEMP2 = MULTI(M, TEMP1);
	FREEMPI(M);
	FREEMPI(TEMP1);
	TEMP3 = SUBI(N, TEMP2);
	FREEMPI(TEMP2);
	RR = CHANGE(1 - r);
	R = CHANGE(r);
	TEMP4 = MULTI(TEMP3, RR);
	FREEMPI(RR);
	TEMP5 = MULTI(TEMP3, R);
	FREEMPI(TEMP3);
	FREEMPI(R);
	TEMP6 = ADDI(TEMP4, L);
	FREEMPI(TEMP4);
	if (s == 0) {
		*X = COPYI(TEMP6);
	}
	else {
		*X = MINUSI(TEMP6);
	}
	TEMP7 = SUBI(TEMP5, L);
	FREEMPI(L);
	FREEMPI(TEMP5);
	if (s == 0) {
		*Y = COPYI(TEMP7);
	}
	else {
		*Y = MINUSI(TEMP7);
	}
	FREEMPI(TEMP6);
	FREEMPI(TEMP7);
	return;
}

void ALG_CYCLE(USL d, MPI *M[], MPI *N[], MPI *X[], MPI *Y[], USL infty, USL RANGE)
/*
* This function searches all trajectories of the generalized Collatz
* function dividing by sqrt(d)), d>1, which start from SPIRAL(n), 0 <= n <= RANGE.
*/
{
	MPI *A, *B, *A2, *B2, *AA, *BB, *AN, *BN, **Z1, **Z2, *TempA, *TempB, *TempA2, *TempB2;
	MPI *Atemp, *Btemp, *NN;
	unsigned int c = 0, r, end, *FLAG, t = 0, i, s = 0;
	USL n;
	USI diverge_flag = 0;
	long int LIMIT, k;
	FILE *outfile;
	char buff[20];

	strcpy(buff, "alg_cycle.out");
	LIMIT = CYCLE_MAX;
	FLAG = (unsigned int *)ccalloc(RANGE + 1, sizeof(unsigned int));
	Z1 = (MPI **)mmalloc(LIMIT * sizeof(MPI *));
	Z2 = (MPI **)mmalloc(LIMIT * sizeof(MPI *));
	INTSETI((int *)FLAG, RANGE + 1, RANGE + 1);
	for (n = 0; n <= RANGE; n++) {
		if (FLAG[n] > n) {
			if (n % 10000 == 0) {
				printf("n = %ld\n", n);
			}
			NN = CHANGE(n);
			SPIRAL(NN, &AN, &BN);
			FREEMPI(NN);
			A = COPYI(AN);
			B = COPYI(BN);
			A2 = COPYI(AN);
			B2 = COPYI(BN);
			FREEMPI(AN);
			FREEMPI(BN);
			end = 0;
			for (k = 0; k <= UPPER_BOUND; k++)
			{
				TempA = A;
				TempB = B;
				ALG_COLLATZ(A, B, d, M, N, X, Y, &AA, &BB);
				FREEMPI(TempA);
				FREEMPI(TempB);
				A = AA;
				B = BB;
				if (A->D >= infty || B->D >= infty) {
					if (diverge_flag == 0) {
						outfile = fopen(buff, "a");

						fprintf(outfile, "divergent trajectory starting from n = %ld\n", n);
						printf("divergent trajectory starting from n = %ld\n", n);
						fclose(outfile);
						diverge_flag = 0;
					}
					FREEMPI(A);
					FREEMPI(B);
					FREEMPI(A2);
					FREEMPI(B2);
					break;
				}
				TempA2 = A2;
				TempB2 = B2;
				ALG_COLLATZ(A2, B2, d, M, N, X, Y, &Atemp, &Btemp);
				FREEMPI(TempA2);
				FREEMPI(TempB2);
				ALG_COLLATZ(Atemp, Btemp, d, M, N, X, Y, &AA, &BB);
				FREEMPI(Atemp);
				FREEMPI(Btemp);
				A2 = AA;
				B2 = BB;
				t = EQUALI(A, A2)*EQUALI(B, B2);
				if (t)
				{
					/* cycle detected starting at n */
					printf("cycle detected starting at n=%lu\n", n);
					FREEMPI(A2);
					FREEMPI(B2);
					break;
				}
				if (A->D == 0 && A->V[0] <= RANGE && B->D == 0 && B->V[0] <= RANGE)
				{
					r = INVERSE_SPIRAL(A, B);
					if (r <= RANGE) {
						if (FLAG[r] < n) {
							end = 1;
							FREEMPI(A);
							FREEMPI(B);
							FREEMPI(A2);
							FREEMPI(B2);
							break;
						}
						else {
							FLAG[r] = n;
						}
					}
				}
			}
			if (end == 1)
				continue;
			if (t) {
				if (c == 0) {
					/* store (A,B) the start of first cycle */
					Z1[0] = A;
					Z2[0] = B;
					printf("\ncycle %u:\n", c + 1);
					outfile = fopen(buff, "a");
					fprintf(outfile, "\ncycle %u:\n", c + 1);
					fclose(outfile);
					ALG_P_CYCLE(A, B, d, M, N, X, Y, c + 1);
					c++;
				}
				else {
					/* store start of subsequent cycles */
					s = ALG_IS_CYCLE(A, B, d, M, N, X, Y, Z1, Z2, c);
					if (s == 0) {
						FREEMPI(A);
						FREEMPI(B);

					}
					else {
						if (c >= LIMIT)
						{
							LIMIT++;
							Z1 = (MPI **)rrealloc((char *)Z1, (c + 1) * sizeof(MPI *), sizeof(MPI *));
							Z2 = (MPI **)rrealloc((char *)Z2, (c + 1) * sizeof(MPI *), sizeof(MPI *));
						}
						Z1[c] = A;
						Z2[c] = B;
						printf("\ncycle %u:\n", c + 1);
						outfile = fopen(buff, "a");
						fprintf(outfile, "\ncycle %u:\n", c + 1);
						fclose(outfile);
						ALG_P_CYCLE(A, B, d, M, N, X, Y, c + 1);
						c++;
					}
				}
			}
		}
	}
	for (i = 0; i < c; i++) {
		FREEMPI(Z1[i]);
		FREEMPI(Z2[i]);
	}
	ffree((char *)Z1, LIMIT * sizeof(MPI *));
	ffree((char *)Z2, LIMIT * sizeof(MPI *));
	ffree((char *)FLAG, (RANGE + 1) * sizeof(USI));
	return;
}

USL INVERSE_SPIRAL(MPI *X, MPI *Y) {
	USL n;
	long a, b, k, m, s, x, y;
	int t;

	x = CONVERTII(X);
	y = CONVERTII(Y);
	if (x >= 0) {
		a = x;
	}
	else {
		a = -x;
	}
	if (y >= 0) {
		b = y;
	}
	else {
		b = -y;
	}
	if (a >= b) {
		k = a;
	}
	else {
		k = b;
	}
	t = COMPAREI(X, Y);
	if (t<0) {
		s = -1;
	}
	else {
		s = 1;
	}
	m = 4 * k*k + s*(2 * k + x + y);
	n = (USL)m;
	return(n);
}

MPI *INVERSE_SPIRAL1(MPI *X, MPI *Y) {
	MPI *A, *B, *K, *M, *TEMP1, *TEMP2, *TEMP3, *TEMP4;
	int t;

	A = ABSI(X);
	B = ABSI(Y);
	if (RSV(A, B) >= 0) {
		K = COPYI(A);
	}
	else {
		K = COPYI(B);
	}
	FREEMPI(A);
	FREEMPI(B);
	t = COMPAREI(X, Y);
	TEMP1 = MULT_I(K, 4);
	TEMP2 = MULTI(TEMP1, K);
	FREEMPI(TEMP1);
	TEMP1 = ADDI(X, Y);
	TEMP3 = MULT_I(K, 2);
	FREEMPI(K);
	TEMP4 = ADDI(TEMP3, TEMP1);
	FREEMPI(TEMP1);
	FREEMPI(TEMP3);
	if (t == 1) {
		M = ADDI(TEMP2, TEMP4);
	}
	else {
		M = SUBI(TEMP2, TEMP4);
	}
	FREEMPI(TEMP2);
	FREEMPI(TEMP4);
	return(M);
}
