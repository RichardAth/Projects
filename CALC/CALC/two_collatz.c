/* two_collatz.c */
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

MPI *TWO_COLLATZ(MPI *Aptr, MPI *m0, MPI *m1, MPI *a0, MPI *a1, MPI *x0, MPI *x1)
/*
* A generalization of Benoit Cloitre's function.
* T(Aptr)=INT(ar*Aptr/mr)+xr, if Aptr(MOD(2))=r.
*/

{
	USL r;
	MPI *T1, *T2, *T3;

	r = MOD_(Aptr, (USL)2);
	if (r == 0) {
		T1 = MULTI(Aptr, a0);
		T2 = INT(T1, m0);
		FREEMPI(T1);
		T3 = ADDI(T2, x0);
		FREEMPI(T2);
	}
	else {
		T1 = MULTI(Aptr, a1);
		T2 = INT(T1, m1);
		FREEMPI(T1);
		T3 = ADDI(T2, x1);
		FREEMPI(T2);
	}
	return (T3);
}

void TWO_P_CYCLE(MPI *Aptr, MPI *m0, MPI *m1, MPI *a0, MPI *a1, MPI *x0, MPI *x1, USI c)

/* the cycle starting with *Aptr is printed */
{
	MPI *B, *Temp, *MIN_ELT, *TEMP1;
	unsigned int i = 0;
	FILE *outfile;
	char buff[20];

	strcpy(buff, "two_cycle.out");
	outfile = fopen(buff, "a");
	B = COPYI(Aptr);
	MIN_ELT = COPYI(Aptr);
	do {
		Temp = B;
		B = TWO_COLLATZ(B, m0, m1, a0, a1, x0, x1);
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
		B = TWO_COLLATZ(B, m0, m1, a0, a1, x0, x1);
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

unsigned int TWO_IS_CYCLE(MPI *Aptr, MPI *m0, MPI *m1, MPI *a0, MPI *a1, MPI *x0, MPI *x1, MPI *Z[], USI c)
/* the cycle starting from Aptr is distinct from those starting from
* Z[0],..., Z[c]  if and only if TWO_IS_CYCLE() = 1.
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
		B = TWO_COLLATZ(B, m0, m1, a0, a1, x0, x1);
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

void TWO_CYCLE(MPI *m0, MPI *m1, MPI *a0, MPI *a1, MPI *x0, MPI *x1, USL infty, USL RANGE)
/*
* This function searches all trajectories of the TWO_COLLATZ
* function which start from p, |p| <= RANGE/2 (RANGE an even integer).
*/
{
	MPI *A, *B, **Z, *Temp, *B1;
	unsigned int c = 0, n, r, end, *FLAG, t = 0, i, s = 0;
	USI diverge_flag = 0;
	long int p, q, LIMIT, k;
	FILE *outfile;
	char buff[20];

	strcpy(buff, "two_cycle.out");
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
				A = TWO_COLLATZ(A, m0, m1, a0, a1, x0, x1);
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
				B1 = TWO_COLLATZ(B, m0, m1, a0, a1, x0, x1);
				B = TWO_COLLATZ(B1, m0, m1, a0, a1, x0, x1);
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
					TWO_P_CYCLE(A, m0, m1, a0, a1, x0, x1, c + 1);
					c++;
				}
				else {
					/* store start of subsequent cycles */
					s = TWO_IS_CYCLE(A, m0, m1, a0, a1, x0, x1, Z, c);
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
						TWO_P_CYCLE(A, m0, m1, a0, a1, x0, x1, c + 1);
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

void TWO_CYCLEX()
{
	unsigned int u;
	int s;
	USL infty, RANGE;
	MPI *a0, *a1, *m0, *m1, *x0, *x1;
	FILE *outfile;
	char buff[20];

	strcpy(buff, "two_cycle.out");
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

	printf("enter the nonzero multiplier a0:");
	a0 = INPUTI(&u);
	printf("a0 = ");
	fprintf(outfile, "a0 = ");
	PRINTI(a0);
	FPRINTI(outfile, a0);
	printf("\n");
	fprintf(outfile, "\n");

	printf("enter the nonzero multiplier a1:");
	a1 = INPUTI(&u);
	printf("a1 = ");
	fprintf(outfile, "a1 = ");
	PRINTI(a1);
	FPRINTI(outfile, a1);
	printf("\n");
	fprintf(outfile, "\n");

	printf("enter the nonzero divisor m0:");
	m0 = INPUTI(&u);
	printf("m0 = ");
	fprintf(outfile, "m0 = ");
	PRINTI(m0);
	FPRINTI(outfile, m0);
	printf("\n");
	fprintf(outfile, "\n");

	printf("enter the nonzero divisor m1:");
	m1 = INPUTI(&u);
	printf("m1 = ");
	fprintf(outfile, "m1 = ");
	PRINTI(m1);
	FPRINTI(outfile, m1);
	printf("\n");
	fprintf(outfile, "\n");

	printf("enter the shift x0:");
	x0 = INPUTI(&u);
	printf("x0 = ");
	fprintf(outfile, "x0 = ");
	PRINTI(x0);
	FPRINTI(outfile, x0);
	printf("\n");
	fprintf(outfile, "\n");

	printf("enter the shift x1:");
	x1 = INPUTI(&u);
	printf("x1 = ");
	fprintf(outfile, "x1 = ");
	PRINTI(x1);
	FPRINTI(outfile, x1);
	printf("\n");
	fprintf(outfile, "\n");
	fclose(outfile);
	Flush();
	printf("\n\n");
	TWO_CYCLE(m0, m1, a0, a1, x0, x1, infty, RANGE);
	FREEMPI(a0);
	FREEMPI(a1);
	FREEMPI(m0);
	FREEMPI(m1);
	FREEMPI(x0);
	FREEMPI(x1);
	return;
}
