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

MPI *LUPEI(MPI *N, USL d)
/*
* Lu Pei's generalized 3x+1 function. (He only did the case d = 3.)
* Here d is the number of branches.
* m[0]=1 and m[1]=...=m[d-1]=d+1 are the multipliers,
* X[0]=0,X[1]=...=X[d-1]=1 are the shifts.
* T(N)=INT(m[r]*N/d)+X[r], if N(MOD(d))=r.
*/
{
	USL dplus1;
	MPI *TEMP, *COPY, *Z, *I, *D;

	D = CHANGEI((long)d);
	I = HALFMOD(N, D);
	FREEMPI(D);
	if (EQZEROI(I)) {
		Z = INT_(N, d);
	}
	else {
		dplus1 = d + 1;
		TEMP = MULT_I(N, dplus1);
		COPY = TEMP;
		TEMP = SUBI(TEMP, I);
		FREEMPI(COPY);
		Z = INT_(TEMP, d);
		FREEMPI(TEMP);
	}
	FREEMPI(I);
	return (Z);
}

void LUPEI_P_CYCLE(MPI *Aptr, USL d, USI c, USI verbose)

/* the cycle starting with *Aptr is printed */
{
	MPI *B, *Temp, *MIN_ELT, *TEMP1;
	unsigned int i = 0;
	FILE *outfile;
	char buff[20];

	strcpy(buff, "lupei_cycle.out");
	B = COPYI(Aptr);
	MIN_ELT = COPYI(Aptr);
	do {
		Temp = B;
		B = LUPEI(B, d);
		FREEMPI(Temp);
		if (RSV(B, MIN_ELT) < 0) {
			TEMP1 = MIN_ELT;
			MIN_ELT = COPYI(B);
			FREEMPI(TEMP1);
		}
		i++;
	} while (!EQUALI(Aptr, B));

	if (i>1) {

		printf("cycle %u: ", c);
		outfile = fopen(buff, "a");
		fprintf(outfile, "cycle %u: ", c);
		if (!verbose) {
			printf("minumum element ");
			fprintf(outfile, "minumum element ");
			PRINTI(MIN_ELT);
			FPRINTI(outfile, MIN_ELT);
			printf(" cycle-length %u\n", i);
			fprintf(outfile, " cycle-length %u\n", i);
			fclose(outfile);
			FREEMPI(B);
			FREEMPI(MIN_ELT);
			return;
		}
		i = 0;
		Temp = B;
		B = COPYI(MIN_ELT);
		FREEMPI(Temp);
		do {
			PRINTI(B);
			FPRINTI(outfile, B);
			printf(",");
			fprintf(outfile, ",");
			Temp = B;
			B = LUPEI(B, d);
			FREEMPI(Temp);
			i++;
		} while (!EQUALI(MIN_ELT, B));
		printf(" cycle-length %u\n", i);
		fprintf(outfile, " cycle-length %u\n", i);
		fclose(outfile);
		FREEMPI(B);
		FREEMPI(MIN_ELT);
		return;
	}
	else {
		FREEMPI(B);
		FREEMPI(MIN_ELT);
		return;
	}
	return;
}

unsigned int LUPEI_IS_CYCLE(MPI *Aptr, USL d, MPI *Z[], USI c)
/* the cycle starting from Aptr is distinct from those starting from
* Z[0],..., Z[c]  if and only if LUPEI_IS_CYCLE() = 1.
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
		B = LUPEI(B, d);
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

USI LUPEI_CYCLE(USL d, USL infty, USL RANGE, USI verbose)
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

	strcpy(buff, "lupei_cycle.out");
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
			/*	printf("p = %ld\n", p);*/
			A = CHANGEI(p);
			B = COPYI(A);
			end = 0;
			for (k = 0; k <= UPPER_BOUND; k++)
			{
				Temp = A;
				A = LUPEI(A, d);
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
				B1 = LUPEI(B, d);
				B = LUPEI(B1, d);
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
					/*printf("\ncycle %u:\n", c + 1);
					outfile = fopen(buff, "a");
					fprintf(outfile, "\ncycle %u:\n", c + 1);
					fclose(outfile);*/
					LUPEI_P_CYCLE(A, d, c + 1, verbose);
					c++;
				}
				else {
					/* store start of subsequent cycles */
					s = LUPEI_IS_CYCLE(A, d, Z, c);
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
						/*printf("\ncycle %u:\n", c + 1);
						outfile = fopen(buff, "a");
						fprintf(outfile, "\ncycle %u:\n", c + 1);
						fclose(outfile);*/
						LUPEI_P_CYCLE(A, d, c + 1, verbose);
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
	return(c);
}

void LUPEI_CYCLEX()
{
	int s;
	USI cycle_count, verbose;
	USL d, infty, RANGE, m, n;
	FILE *outfile;
	char buff[20];

	strcpy(buff, "lupei_cycle.out");
	if (access(buff, R_OK) == 0)
		unlink(buff);
	outfile = fopen(buff, "w");
	printf("enter 'infty', (the cut-off number of allowable digits (base %u) for an iterate - say <20):", R0);
	s = scanf("%lu", &infty);
	printf("U = %lu\n", infty);
	fprintf(outfile, "U = %lu\n", infty);
	printf("enter an even integer 'RANGE' (<= (say)500,000): (searches all trajectories of the generalized Collatz function which start from p, |p| <= RANGE)/2\n");
	s = scanf("%lu", &RANGE);
	printf("R = %lu\n", RANGE);
	fprintf(outfile, "R = %lu\n", RANGE);
	printf("verbose? (1) nonverbose (0)\n");
	s = scanf("%u", &verbose);
	printf("enter the number m\n");
	s = scanf("%lu", &m);
	printf("the lower value of d is m = %lu\n", m);
	fprintf(outfile, "the lower value of d is m = %lu\n", m);
	printf("enter the number n\n");
	s = scanf("%lu", &n);
	printf("the upper value of d is n = %lu\n", n);
	fprintf(outfile, "the upper value of d is n = %lu\n", n);
	fclose(outfile);

	for (d = m; d <= n; d++) {
		printf("d=%lu:", d);
		outfile = fopen(buff, "a");
		fprintf(outfile, "d=%lu:", d);
		fclose(outfile);
		cycle_count = LUPEI_CYCLE(d, infty, RANGE, verbose);
		outfile = fopen(buff, "a");
		if (cycle_count == d) {
			printf("No extra cycles\n");
			fprintf(outfile, "No extra cycles\n");
		}
		else {
			printf("\n");
			fprintf(outfile, "\n");
		}
		fclose(outfile);
	}
	return;
}
