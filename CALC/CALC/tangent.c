/* tangent.c */
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

/* This constructs the tangent function t(n).
* See slides by Richard Brent at http://maths.anu.edu.au/~brent/talks.html.
*/
MPI *TANGENT(USL n) {
	USL k, j, t;
	MPIA b;
	MPI *ONE, *TEMP, *TEMP1, *TEMP2;

	b = BUILDMPIA();
	ONE = ONEI();
	ADD_TO_MPIA(b, ONE, 1);
	for (k = 2; k <= n; k++) {
		TEMP = MULT_I(b->A[k - 1], k - 1);
		ADD_TO_MPIA(b, TEMP, k);
		FREEMPI(TEMP);
	}
	for (k = 2; k <= n; k++) {
		for (j = k; j <= n; j++) {
			t = j - k;
			TEMP1 = MULT_I(b->A[j - 1], t);
			TEMP2 = MULT_I(b->A[j], t + 2);
			TEMP = ADD0I(TEMP1, TEMP2);
			FREEMPI(TEMP1);
			FREEMPI(TEMP2);
			ADD_TO_MPIA(b, TEMP, j);
			FREEMPI(TEMP);
		}
	}
	TEMP = COPYI(b->A[n]);
	FREEMPIA(b);
	FREEMPI(ONE);
	return(TEMP);
}

MPI *TANGENTX(MPI *N) {
	MPI *NULL_MPI = NULL, *TEMP;
	USL n;

	if (N->S <= 0) {
		printf("N <= 0\n");
		return(NULL_MPI);
	}
	if (N->D > 0) {
		printf("N > R0\n");
		return(NULL_MPI);
	}
	if (N->V[0] > 2000) {
		printf("N > 2000\n");
		return(NULL_MPI);
	}
	n = CONVERTI(N);
	TEMP = TANGENT(n);
	return(TEMP);
}

void BERNOULLI(USL n, MPI **BERNOULLI_NUMERATOR, MPI **BERNOULLI_DENOMINATOR) {
	USL t;
	int s;
	MPI *TEMP, *TEMP1, *TEMP2, *K, *G, *T;

	if (n == 0) {
		*BERNOULLI_NUMERATOR = ONEI();
		*BERNOULLI_DENOMINATOR = ONEI();
		return;
	}
	if (n == 1) {
		*BERNOULLI_NUMERATOR = MINUS_ONEI();
		*BERNOULLI_DENOMINATOR = TWOI();
		return;
	}
	if (n % 2) {
		*BERNOULLI_NUMERATOR = ZEROI();
		*BERNOULLI_DENOMINATOR = ONEI();
		return;
	}
	t = n / 2;
	if (t - 1 % 2) {
		s = -1;
	}
	else {
		s = 1;
	}
	T = TANGENT(t);
	TEMP1 = MULT_I(T, n);
	FREEMPI(T);
	if (s < 0) {
		TEMP1->S = -(TEMP1->S);
	}
	K = POWER_I((long)2, (unsigned int)n);
	TEMP = SUB0_I(K, (USL)1);
	TEMP2 = MULTI(K, TEMP);
	FREEMPI(TEMP);
	FREEMPI(K);
	G = GCD(TEMP1, TEMP2);
	*BERNOULLI_NUMERATOR = INT(TEMP1, G);
	*BERNOULLI_DENOMINATOR = INT(TEMP2, G);
	FREEMPI(TEMP1);
	FREEMPI(TEMP2);
	FREEMPI(G);
}

MPI *BERNOULLIX(MPI *N, MPI **BERNOULLI_NUMERATOR, MPI**BERNOULLI_DENOMINATOR) {
	MPI *NULL_MPI = NULL, *TEMP1, *TEMP2;
	USL n;

	if (N->S < 0) {
		printf("N < 0\n");
		return(NULL_MPI);
	}
	if (N->D > 0) {
		printf("N > R0\n");
		return(NULL_MPI);
	}
	if (N->V[0] > 4000) {
		printf("N > 4000\n");
		return(NULL_MPI);
	}
	n = CONVERTI(N);
	BERNOULLI(n, &TEMP1, &TEMP2);
	*BERNOULLI_NUMERATOR = TEMP1;
	*BERNOULLI_DENOMINATOR = TEMP2;
	return(ONEI());
}

/* This algorithm is taken from http://www.numericana.com/answer/numbers.htm#partitions
* and implements a recurrence relation of Euler - see Hardy and Wright,
* An introduction to the theory of numbers (1962 edition), p. 286 .
* We calculate p(n) for 1 <= n <= 65535 = MAX_ARRAY_SIZE - 1.
*/
MPI *PARTITION(USL m) {
	MPIA P;
	MPI *S, *ONE, *TEMP;
	int i, j, k, kk, t;

	P = BUILDMPIA();
	ONE = ONEI();
	ADD_TO_MPIA(P, ONE, 0);
	FREEMPI(ONE);
	for (i = 1; i <= m; i++) {
		j = 1;
		k = 1;
		S = ZEROI();
		while (j > 0) {
			kk = k * k;
			kk = 3 * kk;
			j = i - (kk + k) / 2;
			if (k % 2) {
				t = -1;
			}
			else {
				t = 1;
			}
			if (j >= 0) {
				TEMP = S;
				if (t == 1) {
					S = SUBI(S, P->A[j]);
				}
				else {
					S = ADDI(S, P->A[j]);
				}
				FREEMPI(TEMP);
			}
			j = i - (kk - k) / 2;
			if (j >= 0) {
				TEMP = S;
				if (t == 1) {
					S = SUBI(S, P->A[j]);
				}
				else {
					S = ADDI(S, P->A[j]);
				}
				FREEMPI(TEMP);
			}
			k++;
		}
		TEMP = S;
		ADD_TO_MPIA(P, S, (USL)i);
		FREEMPI(TEMP);
	}
	TEMP = COPYI(P->A[m]);
	FREEMPIA(P);
	return(TEMP);
}

MPI *PARTITIONX(MPI *N) {
	MPI *NULL_MPI = NULL, *TEMP;
	USL n;

	if (N->S <= 0) {
		printf("N <= 0\n");
		return(NULL_MPI);
	}
	if (N->D > 1) {
		printf("N > 2^32\n");
		return(NULL_MPI);
	}
	n = CONVERTI(N);
	if (n > 65535) {
		printf("n > 65535\n");
		return(NULL_MPI);
	}
	TEMP = PARTITION(n);
	return(TEMP);
}
