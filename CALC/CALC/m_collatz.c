/* m_collatz.c */
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

/*
* Benoit Cloitre's 3f(x) + 1 function.
*/
MPI *M_COLLATZ(MPI *X, MPI *M)
{
	USL r;
	MPI *T, *T1, *T2, *T3, *T4, *TEMP, *ONE;

	T = ADD0_I(M, 1);
	T1 = MULTI(T, X);
	FREEMPI(T);
	T2 = INT(T1, M);
	FREEMPI(T1);
	r = MOD_(T2, 2);
	if (r == 0) {
		T3 = INT_(T2, 2);
	}
	else {
		T4 = MULT_I(T2, 3);
		TEMP = T4;
		ONE = ONEI();
		T4 = ADDI(T4, ONE);
		FREEMPI(TEMP);
		FREEMPI(ONE);
		T3 = INT_(T4, 2);
		FREEMPI(T4);
	}
	FREEMPI(T2);
	return(T3);
}

void M_P_CYCLE(MPI *Aptr, MPI *M, USI c)
/* the cycle starting with *Aptr is printed */
{
	MPI *B, *Temp, *MIN_ELT, *TEMP1;
	unsigned int i = 0;
	FILE *outfile;
	char buff[20];

	strcpy(buff, "m_cycle.out");
	outfile = fopen(buff, "a");
	B = COPYI(Aptr);
	MIN_ELT = COPYI(Aptr);
	do {
		Temp = B;
		B = M_COLLATZ(B, M);
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
		B = M_COLLATZ(B, M);
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

unsigned int M_IS_CYCLE(MPI *Aptr, MPI *M, MPI *Z[], USI c)
/* the cycle starting from Aptr is distinct from those starting from
* Z[0],..., Z[c]  if and only if M_IS_CYCLE() = 1.
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
		B = M_COLLATZ(B, M);
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

void M_CYCLE(MPI *M, USL infty, USL RANGE)
/*
* This function searches all trajectories of the M_COLLATZ
* function which start from p, |p| <= RANGE/2 (RANGE an even integer).
*/
{
	MPI *A, *B, **Z, *Temp, *B1;
	unsigned int c = 0, n, r, end, *FLAG, t = 0, i, s = 0;
	USI diverge_flag = 0;
	long int p, q, LIMIT, k;
	FILE *outfile;
	char buff[20];

	strcpy(buff, "m_cycle.out");
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
			A = CHANGEI(p);
			B = COPYI(A);
			end = 0;
			for (k = 0; k <= UPPER_BOUND; k++)
			{
				Temp = A;
				A = M_COLLATZ(A, M);
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
				B1 = M_COLLATZ(B, M);
				B = M_COLLATZ(B1, M);
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
					M_P_CYCLE(A, M, c + 1);
					c++;
				}
				else {
					/* store start of subsequent cycles */
					s = M_IS_CYCLE(A, M, Z, c);
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
						M_P_CYCLE(A, M, c + 1);
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

void M_CYCLEX()
{
	unsigned int u;
	int s;
	USL infty, RANGE;
	MPI *M;
	FILE *outfile;
	char buff[20];

	strcpy(buff, "m_cycle.out");
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

	printf("enter the nonzero M:");
	M = INPUTI(&u);
	printf("M = ");
	fprintf(outfile, "M = ");
	PRINTI(M);
	FPRINTI(outfile, M);
	printf("\n");
	fprintf(outfile, "\n");
	fclose(outfile);
	Flush();
	printf("\n\n");
	M_CYCLE(M, infty, RANGE);
	FREEMPI(M);
	return;
}

USI MM_P_CYCLE(MPI *Aptr, MPI *M, USI c)
/* the length of the cycle c starting with *Aptr is returned */
{
	MPI *B, *Temp, *MIN_ELT, *TEMP1;
	unsigned int i = 0;
	FILE *outfile;
	char buff[20];

	strcpy(buff, "mm_cycle.out");
	outfile = fopen(buff, "a");
	B = COPYI(Aptr);
	MIN_ELT = COPYI(Aptr);
	do {
		Temp = B;
		B = M_COLLATZ(B, M);
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
	PRINTI(B);
	FPRINTI(outfile, B);
	do {
		Temp = B;
		B = M_COLLATZ(B, M);
		FREEMPI(Temp);
		i++;
	} while (!EQUALI(MIN_ELT, B));
	printf(" (%u)", i);
	fprintf(outfile, " (%u)", i);
	fclose(outfile);
	FREEMPI(B);
	FREEMPI(MIN_ELT);
	return(i);
}

void MM_CYCLE(MPI *M, USL infty, USL RANGE)
/*
* This function searches all trajectories of the M_COLLATZ
* function which start from p, |p| <= RANGE/2 (RANGE an even integer).
*/
{
	MPI *A, *B, **Z, *Temp, *B1;
	unsigned int c = 0, n, r, end, *FLAG, t = 0, i, s = 0, len;
	USI diverge_flag = 0;
	long int p, q, LIMIT, k;
	FILE *outfile;
	char buff[20];

	strcpy(buff, "mm_cycle.out");
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
			A = CHANGEI(p);
			B = COPYI(A);
			end = 0;
			for (k = 0; k <= UPPER_BOUND; k++)
			{
				Temp = A;
				A = M_COLLATZ(A, M);
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
				B1 = M_COLLATZ(B, M);
				B = M_COLLATZ(B1, M);
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
					printf(" %u:", c + 1);
					outfile = fopen(buff, "a");
					fprintf(outfile, " %u:", c + 1);
					fclose(outfile);
					len = MM_P_CYCLE(A, M, c + 1);
					c++;
				}
				else {
					/* store start of subsequent cycles */
					s = M_IS_CYCLE(A, M, Z, c);
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
						printf(" %u:", c + 1);
						outfile = fopen(buff, "a");
						fprintf(outfile, " %u:", c + 1);
						fclose(outfile);
						len = MM_P_CYCLE(A, M, c + 1);
						c++;
						if (c % 6 == 0) {
							printf("\n\t");
							outfile = fopen(buff, "a");
							fprintf(outfile, "\n\t");
							fclose(outfile);
						}
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

void MM_CYCLEX()
{
	unsigned m;
	USL infty, RANGE, u, v;
	MPI *M;
	FILE *outfile;
	char buff[20];

	infty = (USL)15;
	strcpy(buff, "mm_cycle.out");
	if (access(buff, R_OK) == 0)
		unlink(buff);
	RANGE = (USL)100000;

	printf("enter the lower range u >= 7:");
	scanf("%lu", &u);
	if (u < 7) {
		printf("lower range %lu < 7\n", u);
		return;
	}
	printf("enter the upper range v >= u:");
	scanf("%lu", &v);
	if (u > v) {
		printf("lower range %lu > upper range %lu", u, v);
		return;
	}
	for (m = u; m <= v; m++) {
		M = CHANGE(m);
		printf("\nM = ");
		outfile = fopen(buff, "a");
		fprintf(outfile, "\nM = ");
		PRINTI(M);
		FPRINTI(outfile, M);
		printf(": ");
		fprintf(outfile, ": ");
		fclose(outfile);
		MM_CYCLE(M, infty, RANGE);
		FREEMPI(M);
	}
	printf("\n");
	outfile = fopen(buff, "a");
	fprintf(outfile, "\n");
	fclose(outfile);
	Flush();
	return;
}
