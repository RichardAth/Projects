/* utility.c */
#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_NONSTDC_NO_DEPRECATE
#include <stdio.h>
#include <stdlib.h>
#include "integer.h"
#include "fun.h"
#define LIMIT 4294967295U

extern unsigned long PRIME[];

unsigned long *ERATOSTHENES(USI n)
/*
* returns an array consisting of the first  n primes.
* From "Programming in C" by S. Kochan,p.89.
*/
{
	unsigned long *prime, p, is_prime, i, prime_index = 2;

	prime = (unsigned long *)mmalloc((USL)(n * sizeof(unsigned long)));
	prime[0] = 2;
	prime[1] = 3;

	for (p = 5; p <= LIMIT; p = p + 2)
	{
		is_prime = 1;
		for (i = 1; is_prime && p >= prime[i] * prime[i]; ++i)
		{
			if (p % prime[i] == 0)
				is_prime = 0;
		}
		if (is_prime)
		{
			prime[prime_index] = p;
			++prime_index;
			if (prime_index == n)
				break;
		}
	}
	return prime;
}

void QERATOSTHENES(MPI *N, USL base[], USL soln[], USI FBASE)
/*
* returns an array consisting of the first FBASE primes p for which (N/p) = 1.
* Also returns array soln[i] a solution of x*x=N (mod base[i], 1 <= i <= FBASE.
* Based on "Programming in C" by S. Kochan,p.89.
*/
{
	unsigned long *prime, p, r;
	unsigned int t, T, is_prime, i, prime_index = 2, found = 1;
	int u;

	t = T = 2 * FBASE;
	prime = (unsigned long *)mmalloc((USL)(t * sizeof(unsigned long)));
	prime[0] = 2;
	prime[1] = 3;
	r = MOD0_(N, (USL)3);
	if (r == 1)
	{
		base[found] = 3;
		soln[found] = SQRTm(r, (USL)3);
		found++;
	}

	for (p = 5; p < R0; p = p + 2)
	{
		is_prime = 1;
		for (i = 1; is_prime && p >= prime[i] * prime[i]; ++i)
		{
			if (p % prime[i] == 0)
				is_prime = 0;
		}
		if (is_prime)
		{
			prime[prime_index] = p;
			r = MOD0_(N, prime[prime_index]);
			u = JACOBI(r, prime[prime_index]);
			if (u == 1)
			{
				base[found] = prime[prime_index];
				soln[found] = SQRTm(r, prime[prime_index]);
				/*
				printf("found base[%u] = %lu, soln[%u] = %lu\n", found, base[found] , found, soln[found]);
				*/
				/*
				fprintf(dfile,"found base[%u] = %lu, soln[%u] = %lu\n", found, base[found] , found, soln[found]);
				fflush(dfile);
				*/
				if (found == FBASE)
				{
					ffree((char *)prime, (USL)T * sizeof(unsigned long));
					return;
				}
				found++;
			}
			++prime_index;
			if (prime_index == T)
			{
				rrealloc(prime, (T + 10) * sizeof(unsigned long), (USL)(10 * sizeof(unsigned long)));
				T += 10;
			}
		}
	}
	return;
}

MPI *PROB_PRIME(MPI *M, MPI *N, USI *hptr)
/*
* returns a probable prime Q = 4n+3, with Q = M + hptr, with (N/Q) = 1.
* Also Q->D <= 1.
*/
{
	MPI *Q, *Tmp0, *Tmp1;
	unsigned int r;

	r = MOD0_(M, (USL)12);
	if (r < 7)
	{
		*hptr = 7 - r;
		Q = ADD0_I(M, (USL)(*hptr));
	}
	else if (7 < r && r < 11)
	{
		*hptr = 11 - r;
		Q = ADD0_I(M, (USL)(*hptr));
	}
	else
	{
		*hptr = 0;
		Q = COPYI(M);
	}
	while (1)
	{
		Tmp0 = SUB0_I(Q, (USL)1);
		Tmp1 = MPOWER_(2L, Tmp0, Q);
		FREEMPI(Tmp0);
		if (!EQONEI(Tmp1))
		{
			FREEMPI(Tmp1);
			Tmp1 = Q;
			Q = ADD0_I(Q, (USL)12);
			FREEMPI(Tmp1);
			*hptr += 12;
			continue;
		}
		FREEMPI(Tmp1);
		if (Q_PRIME_TEST(Q) && JACOBIB(N, Q) == 1)
			return (Q);
		Tmp1 = Q;
		Q = ADD0_I(Q, (USL)12);
		FREEMPI(Tmp1);
		*hptr += 12;
	}
}

MPI *PINCH_PRIME(MPI *M, MPI *N, long *base, USI FBASE, USI *hptr)
/*
* returns a probable prime Q = 4n+3, with Q = M + hptr, with (N/Q) = 1.
*/
{
	MPI *Q, *Tmp0, *Tmp1, *H1, *X1, *N1;
	unsigned int r, i, flag, t;

	r = MOD0_(M, (USL)12);
	if (r < 7)
	{
		*hptr = 7 - r;
		Q = ADD0_I(M, (USL)(*hptr));
	}
	else if (7 < r && r < 11)
	{
		*hptr = 11 - r;
		Q = ADD0_I(M, (USL)(*hptr));
	}
	else
	{
		*hptr = 0;
		Q = COPYI(M);
	}
	while (1)
	{
		Tmp0 = ADD0_I(Q, (USL)1);
		X1 = INT0_(Tmp0, (USL)2);
		N1 = MOD0(N, Q);
		H1 = MPOWER(N, X1, Q);
		FREEMPI(Tmp0);
		FREEMPI(X1);
		t = EQUALI(H1, N1);
		FREEMPI(H1);
		FREEMPI(N1);
		if (!t)
		{
			Tmp1 = Q;
			Q = ADD0_I(Q, (USL)12);
			FREEMPI(Tmp1);
			*hptr += 12;
			continue;
		}
		else
		{
			flag = 0;
			for (i = 0; i <= FBASE; i++)
			{
				if (MOD0_(Q, (USL)(base[i])) == 0)
				{
					flag = 1;
					Tmp1 = Q;
					Q = ADD0_I(Q, (USL)12);
					FREEMPI(Tmp1);
					*hptr += 12;
					break;
				}
			}
			if (flag)
				continue;
		}
		return (Q);
	}
}

MPI *NEXT_PRIME(MPI *M, USI *hptr)
/*
* returns a probable prime Q with Q = M + hptr.
*/
{
	MPI *P, *Q, *Tmp0, *Tmp1;
	unsigned int i, flag;
	unsigned long r, p;
	int t;

	if (M->D == 0 && M->V[0] <= PRIME[Y0 - 1])
	{
		p = M->V[0];
		for (i = 0; i < Y0; i++)
		{
			if (p <= PRIME[i])
				return(CHANGE(PRIME[i]));
		}
	}
	/* now M > PRIME[Y0 - 1]. */
	r = M->V[0] % 2;
	*hptr = 1 - r;
	/*	printf("h = %u\n", *hptr); */
	Q = ADD0_I(M, (USL)(*hptr));
	while (1)
	{
		flag = 0;
		for (i = 1; i < Y0; i++)
		{
			if (MOD0_(Q, PRIME[i]) == 0)
			{
				flag = 1;
				break;
			}
		}
		if (flag)
		{
			Tmp1 = Q;
			Q = ADD0_I(Q, (USL)2);
			FREEMPI(Tmp1);
			*hptr += 2;
			printf("h = %u\n", *hptr);
			continue;
		}
		Tmp0 = SUB0_I(Q, (USL)1);
		Tmp1 = MPOWER_(2L, Tmp0, Q);
		FREEMPI(Tmp0);
		if (!EQONEI(Tmp1))
		{
			FREEMPI(Tmp1);
			Tmp1 = Q;
			Q = ADD0_I(Q, (USL)2);
			FREEMPI(Tmp1);
			*hptr += 2;
			printf("h = %u\n", *hptr);
			/* Q fails the base 2 pseudoprime test */
			continue;
		}
		/* Q has passed the base 2 pseudoprime test */
		FREEMPI(Tmp1);
		/* if (Q_PRIME_TEST(Q))*/
		P = LUCAS(Q);
		t = P->S;
		FREEMPI(P);
		if (t)
			/* Q passes Lucas' test */
			return (Q);
		/* Q failed test*/
		Tmp1 = Q;
		Q = ADD0_I(Q, (USL)2);
		FREEMPI(Tmp1);
		*hptr += 2;
		printf("h = %u\n", *hptr);
	}
}

MPI *NEXT_PRIMEX(MPI *M)
/*
* returns the first probable prime after M.
*/
{
	unsigned int h;
	MPI *T;

	T = NEXT_PRIME(M, &h);
	return (T);
}

MPI *NEXTPRIMEAP(MPI *A, MPI *B, MPI *M)
/*
* Finds the first Lucas probable prime P, A | P - B, M <= P.
* Here A is even, B is odd, A > 1 , A > B >= 1, gcd(A,B) = 1, M > 1.
*/
{
	MPI *P, *temp, *T, *K;
	USI flag, i;
	int t;

	/* We start testing P = AK + B, where K is the least integer not less than
	* (M-B)/A. Then P is the least number congruent to B (mod A), P >= M.
	*/
	temp = SUBI(B, M);
	K = INT(temp, A);
	K->S = -(K->S);
	P = MULTABC(B, A, K);
	FREEMPI(temp);
	FREEMPI(K);
	while (1)
	{
		flag = 0;
		for (i = 1; i < Y0; i++)
		{
			if (MOD0_(P, PRIME[i]) == 0)
			{
				if (P->D == 0 && P->V[0] == PRIME[i])
					return(P);
				flag = 1; /* P is composite */
				break;
			}
		}
		if (flag)
		{
			temp = P;
			P = ADD0I(P, A);
			FREEMPI(temp);
			continue;
		}
		printf("testing "); PRINTI(P); printf("\n");
		T = LUCAS(P);
		t = T->S;
		FREEMPI(T);
		if (t)
			break;
		temp = P;
		P = ADD0I(P, A);
		FREEMPI(temp);
	}
	return (P);
}
