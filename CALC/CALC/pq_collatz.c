/* pq_collatz.c */
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
* Benoit Cloitre's 3f(x) + 1 function, where f(x) = int(p*x/q).
*/
MPI *PQ_COLLATZ(MPI *X, MPI *P, MPI *Q)
{
	USL r;
	MPI *T1, *T2, *T3, *T4, *TEMP, *ONE;

	T1 = MULTI(P, X);
	T2 = INT(T1, Q);
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

void PQ_P_CYCLE(MPI *Aptr, MPI *P, MPI *Q, USI c)
/* the cycle starting with *Aptr is printed */
{
	MPI *B, *Temp, *MIN_ELT, *TEMP1;
	unsigned int i = 0;
	FILE *outfile;
	char buff[20];

	strcpy(buff, "pq_cycle.out");
	outfile = fopen(buff, "a");
	B = COPYI(Aptr);
	MIN_ELT = COPYI(Aptr);
	do {
		Temp = B;
		B = PQ_COLLATZ(B, P, Q);
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
		B = PQ_COLLATZ(B, P, Q);
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

unsigned int PQ_IS_CYCLE(MPI *Aptr, MPI *P, MPI *Q, MPI *Z[], USI c)
/* the cycle starting from Aptr is distinct from those starting from
* Z[0],..., Z[c]  if and only if PQ_IS_CYCLE() = 1.
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
		B = PQ_COLLATZ(B, P, Q);
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

void PQ_CYCLE(MPI *P, MPI *Q, USL infty, USL RANGE)
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

	strcpy(buff, "pq_cycle.out");
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
				A = PQ_COLLATZ(A, P, Q);
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
				B1 = PQ_COLLATZ(B, P, Q);
				B = PQ_COLLATZ(B1, P, Q);
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
					PQ_P_CYCLE(A, P, Q, c + 1);
					c++;
				}
				else {
					/* store start of subsequent cycles */
					s = PQ_IS_CYCLE(A, P, Q, Z, c);
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
						PQ_P_CYCLE(A, P, Q, c + 1);
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

void PQ_CYCLEX()
{
	unsigned int u;
	int s;
	USL infty, RANGE;
	MPI *P, *Q, *G, *TEMP;
	FILE *outfile;
	char buff[20];

	strcpy(buff, "pq_cycle.out");
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

	printf("enter the nonzero positive integer P:");
	P = INPUTI(&u);
	printf("P = ");
	fprintf(outfile, "P = ");
	PRINTI(P);
	FPRINTI(outfile, P);
	printf("\n");
	fprintf(outfile, "\n");
	fclose(outfile);
	Flush();
	printf("enter the nonzero positive integer Q, (P > Q > 1):");
	Q = INPUTI(&u);
	printf("Q = ");
	fprintf(outfile, "Q = ");
	PRINTI(Q);
	FPRINTI(outfile, Q);
	printf("\n");
	fprintf(outfile, "\n");
	fclose(outfile);
	Flush();
	if (RSV(P, Q) <= 0) { /* Now we have |P| > |Q| */
		printf("P <= Q\n");
		FREEMPI(P);
		FREEMPI(Q);
		return;
	}
	if (P->S <= 0 || Q->S <= 0) {
		printf("P <= 0 or  Q <= 0\n"); /* Now we have P > Q > 0 */
		FREEMPI(P);
		FREEMPI(Q);
		return;
	}
	TEMP = MOD0(P, Q);
	if (EQZEROI(TEMP)) {
		printf("Q divides P\n"); /* Now we have P > Q > 1 and Q does not divide P */
		FREEMPI(TEMP);
		FREEMPI(P);
		FREEMPI(Q);
		return;
	}
	G = GCD(P, Q);
	TEMP = P;
	P = INT0(P, G);
	FREEMPI(TEMP);
	TEMP = Q;
	Q = INT0(Q, G); /* Now gcd(P, Q) = 1 */
	FREEMPI(TEMP);
	printf("\n\n");
	PQ_CYCLE(P, Q, infty, RANGE);
	FREEMPI(P);
	FREEMPI(Q);
	return;
}
