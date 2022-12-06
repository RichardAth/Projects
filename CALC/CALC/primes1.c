/* primes1.c */
#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_NONSTDC_NO_DEPRECATE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "integer.h"
#include "fun.h"
#include "bigprime.h"
#ifdef _WIN32
#include "unistd_DOS.h"
#else
#include <unistd.h>
#endif


extern unsigned long PRIME[];
extern unsigned int reciprocal[][Z0];
extern unsigned int VERBOSE;
extern unsigned int FREE_PRIMES;/* FREE_PRIMES = 1 destroys the Q_[i] below */
extern unsigned int ECMAX;
#define JMAX 500
USL conway_limit;


MPI *PERFECT_POWER(MPI *N)
/*
* Here N > 1.
* Returns X if N = X^k, for some X, k > 1,  NULL otherwise.
* See E. Bach and J. Sorenson, "Sieve algorithms for perfect power testing"
* Algorithmica 9 (1993) 313-328.
*/
{
	unsigned int i, t;
	int s;
	MPI *X, *Y;

	t = BINARYB(N) - 1;
	i = 0;
	while (PRIME[i] <= t)
	{
		X = BIG_MTHROOT(N, (USI)PRIME[i]);
		Y = POWERI(X, (USI)PRIME[i]);
		s = RSV(Y, N);
		FREEMPI(Y);
		if (s == 0)
			return (X);
		else
		{
			i++;
			FREEMPI(X);
		}
	}
	return ((MPI *)NULL);
}

MPI *POLLARD(MPI *Nptr)
/*
* Pollard's p-1 method of factoring *Nptr.
*/
{
	unsigned long i, b = 1;
	MPI *T, *P, *TT, *TMP, *PP;

	PP = ONEI();
	T = ADD0I(PP, PP);
	while (b <= 2) {
		b++;
		printf("b = %lu\n", b);
		i = 2;
		while (i <= 10000)
		{
			if (i % 100 == 0)
				printf("i = %lu\n", i);
			TMP = T;
			T = MPOWER_M(T, i, Nptr);
			FREEMPI(TMP);
			TT = SUB0I(T, PP);
			P = GCD(Nptr, TT);
			FREEMPI(TT);
			if (!EQONEI(P) && RSV(P, Nptr) == -1)
			{
				PRINTI(P);
				printf(" is a proper factor of ");
				PRINTI(Nptr);
				printf("\n");
				FREEMPI(PP);
				FREEMPI(T);
				return (P);
			}
			if (EQUALI(P, Nptr))
			{
				FREEMPI(P);
				b++;
				TMP = T;
				printf("GCD = *Nptr; b increased to %lu\n", b);
				T = CHANGE(b);
				FREEMPI(TMP);
				continue;
			}
			FREEMPI(P);
			i++;
		}
	}
	printf("no factors returned\n");
	FREEMPI(PP);
	FREEMPI(T);
	return ((MPI *)NULL);
}

MPI *BRENT_POLLARD(MPI *Nptr)
/*
* the Brent-Pollard method for returning a proper factor of
* a composite MPI *Nptr. (see R. Brent, BIT 20, 176 - 184).
*/
{
	MPI  *X, *Y, *Z, *Temp, *Gptr, *PRODUCT, *Temp1;
	unsigned int g, k;
	unsigned long a, r, i, s, t;

	if (VERBOSE)
	{
		printf("BRENT-POLLARD IS SEARCHING FOR A PROPER FACTOR OF ");
		PRINTI(Nptr);
		printf("\n");
	}
	a = 1;
	g = 1;
	r = 1;
	PRODUCT = ONEI();
	Gptr = ONEI();
	Y = ZEROI();
	while (g == 1)
	{
		X = COPYI(Y);
		for (i = 1; i <= r; i++)
		{
			Temp = Y;
			Y = MPOWER_M(Y, (USL)2, Nptr);
			FREEMPI(Temp);
			Temp = Y;
			Y = ADD0_I(Y, a);
			FREEMPI(Temp);
			Temp = Y;
			Y = MOD0(Y, Nptr);
			FREEMPI(Temp);
		}
		k = 0;
		while (1)
		{
			Temp = Y;
			Y = MULTI(Y, Y);
			FREEMPI(Temp);
			Temp = Y;
			Y = ADD0_I(Y, a);
			FREEMPI(Temp);
			Temp = Y;
			Y = MOD0(Y, Nptr);
			FREEMPI(Temp);
			k++;
			Z = SUBI(X, Y);
			Temp1 = PRODUCT;
			PRODUCT = MULTI(PRODUCT, Z);
			FREEMPI(Z);
			FREEMPI(Temp1);
			if (k % 10 == 0)
			{
				FREEMPI(Gptr);
				Gptr = GCD(PRODUCT, Nptr);
				FREEMPI(PRODUCT);
				PRODUCT = ONEI();
				if (Gptr->D || Gptr->V[0] > 1)
				{
					g = 0;
					break;
				}
			}
			if (k >= r)
				break;
		}
		FREEMPI(X);
		r = 2 * r;
		if (VERBOSE)
			printf("r = %lu\n", r);
		/*s = r & 8192;*/
		/*s = r & 33554432;*/
		s = r & 65536;
		t = EQUALI(Gptr, Nptr);
		if (s || t)
		{
			g = 1;
			if (t)
			{
				a++;
				if (VERBOSE)
					printf("a increased to %lu\n", a);
			}
			FREEMPI(Gptr);
			Gptr = ONEI();
			r = 1;
			FREEMPI(Y);
			Y = ZEROI();
			if (s || (a == 10))
			{
				if (VERBOSE)
				{
					printf(" BRENT_POLLARD failed to find a factor of ");
					PRINTI(Nptr);
					printf("\n");
				}
				FREEMPI(Gptr);
				FREEMPI(PRODUCT);
				FREEMPI(Y);
				return ((MPI *)NULL);
			}
		}
		/*
		if (g)
		FREEMPI(Gptr);
		*/
	}
	FREEMPI(Y);
	if (VERBOSE)
	{
		printf("--\n");
		printf("BRENT_POLLARD FINISHED: \n");
		PRINTI(Gptr);
		printf(" IS A PROPER FACTOR OF ");
		PRINTI(Nptr);
		printf("--\n");
	}
	FREEMPI(PRODUCT);
	return (Gptr);
}

unsigned long Q[JMAX];
unsigned int j, j_, K[JMAX], K_[JMAX], ltotal;
MPI *Q_[JMAX], *Q_PRIME[JMAX];

MPI *BABY_DIVIDE(MPI *Nptr)
/*
* *Eptr = 1 if *Nptr is composed of primes < 1000, otherwise *Eptr
* is the part of *Nptr composed of primes > 1000,
* The prime factors of *Nptr < 1000 are stored in the global array Q[]
* and the corresponding exponents in the global array K[].
*/
{
	MPI *X, *Temp;
	unsigned int a = Y0, i, k;
	unsigned long y;

	X = COPYI(Nptr);
	for (i = 0; i < a; i++)
	{
		k = 0;
		y = MOD0_(X, PRIME[i]);
		if (y == 0)
		{
			while (y == 0)
			{
				k++;
				Temp = X;
				X = INT0_(X, PRIME[i]);
				FREEMPI(Temp);
				y = MOD0_(X, PRIME[i]);
			}
			Q[j] = PRIME[i];
			K[j] = k;
			j++;
			if (j == JMAX)
				execerror("j = JMAX, increase JMAX", "");
			if (VERBOSE) {
				printf("%lu is a prime factor of ", PRIME[i]);
				PRINTI(Nptr);
				printf(" exponent: %u\n", k);
				printf("--\n");
			}
		}
	}
	return (X);
}

unsigned int MILLER(MPI *Nptr, USL b)
/*
* *Nptr is odd, > 1, and does not divide b, 0 < b < R0.
*  if 1 is returned, then *Nptr passes Miller's test for base b.
* if 0 is returned, then *Nptr is composite.
*/
{
	MPI *A, *C, *D, *Temp;
	unsigned int i = 0, j = 0;

	D = SUB0_I(Nptr, (USL)1);
	A = INT0_(D, (USL)2);
	while (MOD0_(A, (USL)2) == 0)
	{
		i++;
		Temp = A;
		A = INT0_(A, (USL)2);
		FREEMPI(Temp);
	}
	C = MPOWER_((long)b, A, Nptr);
	if (C->D == 0 && C->V[0] == 1)
	{
		FREEMPI(C);
		FREEMPI(D);
		FREEMPI(A);
		return (1);
	}
	while (1)
	{
		if (EQUALI(C, D))
		{
			FREEMPI(C);
			FREEMPI(D);
			FREEMPI(A);
			return (1);
		}
		j++;
		if (i < j)
		{
			FREEMPI(C);
			FREEMPI(D);
			FREEMPI(A);
			return (0);
		}
		Temp = C;
		C = MULTI(C, C);
		FREEMPI(Temp);
		Temp = C;
		C = MOD0(C, Nptr);
		FREEMPI(Temp);
		Temp = A;
		A = MULT_I(A, 2L);
		FREEMPI(Temp);
	}
}

unsigned int Q_PRIME_TEST(MPI *Nptr)
/*
* *Nptr > 1 and not equal to PRIME[0],...,PRIME[4].
* if 1 is returned, then *Nptr passes Miller's test for bases PRIME[0],
* ...,PRIME[16] and is likely to be prime.
* if 0 is returned, then *Nptr is composite.
*/
{
	unsigned int i;

	for (i = 0; i <= 16; i++)
	{
		if (MILLER(Nptr, PRIME[i]) == 0)
		{
			if (VERBOSE)
			{
				printf("MILLER'S TEST FINISHED: \n");
				PRINTI(Nptr);
				printf(" is composite \n");
				printf("--\n");
			}
			return (0);
		}
	}
	if (VERBOSE)
	{
		PRINTI(Nptr);
		printf(" passed MILLER'S TEST\n");
		printf("--\n");
	}
	return (1);

}

MPI *BIG_FACTOR(MPI *Nptr)
/*
* *Nptr > 1 is not divisible by PRIMES[0], ..., PRIMES[167].
* BIG_FACTOR returns a factor *Eptr > 1000 of *Nptr which is either
* < 1,000,000 (and hence prime) or which passes Miller's test for bases
* PRIMES[0],...,PRIMES[4] and is therefore likely to be prime.
*/
{
	MPI *B, *X, *Y, *Temp;
	unsigned int f;

	B = BUILDMPI(2);
	B->V[0] = 16960;
	B->V[1] = 15;
	B->S = 1; /* *B = 1,000,000. */

	if (RSV(Nptr, B) == -1 || Q_PRIME_TEST(Nptr))
	{
		FREEMPI(B);
		return (COPYI(Nptr));
	}
	f = 1;
	X = COPYI(Nptr);
	while (f)
	{
		Y = BRENT_POLLARD(X);
		if (Y == NULL)
			Y = PERFECT_POWER(X);
		if (Y == NULL)
			Y = MPQS1(X);
		if (Y == NULL)
		{
			printf("Switching to ECF:%u elliptic curves\n", ECMAX);
			Y = EFACTOR(X, 1279, 1);
		}
		if (Y == NULL)
		{
			printf("Switching to POLLARD p-1\n");
			Y = POLLARD(X);
		}
		if (Y == NULL) {
			printf("Trying to find a factor of composite integer ");
			PRINTI(X);
			printf("\n");
			execerror("Unfortunately no factor was found", "");
		}
		Temp = X;
		X = INT0(X, Y);
		FREEMPI(Temp);
		if (RSV(X, Y) == 1)
		{
			Temp = X;
			X = COPYI(Y);
			FREEMPI(Temp);
		}
		if (RSV(X, B) == -1)
		{
			FREEMPI(B);
			FREEMPI(Y);
			return (X);
		}
		if (Q_PRIME_TEST(X))
			f = 0;
		FREEMPI(Y);
	}
	FREEMPI(B);
	return (X);
}

unsigned int PRIME_FACTORS(MPI *Nptr)
/*
* A quasi-prime (q-prime) factor of *Nptr is > 1,000,000,
* is not divisible by PRIMES[0],...,PRIMES[167], passes Millers' test
* for bases PRIMES[0],...,PRIMES[10] and is hence likely to be prime.
* PRIME_FACTORS returns the number of q-prime factors of *Nptr, stores
* them in the global array Q_PRIME[].
* Any prime factors < 1000 of *Nptr and corresponding exponents
* are printed and placed in the arrays Q[] and K[], while any prime factors
* > 1000 and all q-prime factors and corresponding exponents of *Nptr are
* printed and placed in the arrays Q_[] and K_[].
*/
{
	MPI *B, *P, *X, *Z, *Temp;
	unsigned int k, t;

	B = BUILDMPI(2);
	B->V[0] = 16960;
	B->V[1] = 15;
	B->S = 1;
	if (VERBOSE)
	{
		printf("FACTORIZING ");
		PRINTI(Nptr);
		printf("--\n");
	}
	X = BABY_DIVIDE(Nptr);
	t = ltotal;
	while (X->D || X->V[0] != 1)
	{
		k = 0;
		P = BIG_FACTOR(X);
		while (1)
		{
			Z = MOD0(X, P);
			if (Z->S != 0)
			{
				FREEMPI(Z);
				break;
			}
			FREEMPI(Z);
			k++;
			Temp = X;
			X = INT0(X, P);
			FREEMPI(Temp);
		}
		if (VERBOSE)
		{
			if (RSV(P, B) == -1)
			{
				PRINTI(P);
				printf(" is a prime factor of ");
				PRINTI(Nptr);
				printf("\n");
			}
		}
		if (RSV(P, B) == 1)
		{
			if (VERBOSE)
			{
				PRINTI(P);
				printf(" is a q-prime factor of ");
				PRINTI(Nptr);
				printf("\n");
			}
			Q_PRIME[ltotal] = COPYI(P);
			ltotal++;
		}
		Q_[j_] = COPYI(P);
		FREEMPI(P);
		K_[j_] = k;
		j_++;
		if (j_ == JMAX)
			execerror("j_ = JMAX, increase JMAX", "");
		if (VERBOSE)
		{
			printf(" exponent: %u\n", k);
			printf("--\n");
		}
	}
	if (VERBOSE)
	{
		printf("FACTORIZATION INTO PRIMES AND Q-PRIMES COMPLETED FOR ");
		PRINTI(Nptr);
		printf("\n------------------------------------\n");
	}
	FREEMPI(B);
	FREEMPI(X);
	return (ltotal - t);
}

unsigned int SELFRIDGE(MPI *Nptr, USI *uptr)
/*
* Selfridges's test for primality - see "Prime Numbers and Computer
* Methods for Factorization" by H. Riesel, Theorem 4.4, p.106.
* *Nptr > 1 is a q-prime.
* The prime and q-prime factors of n - 1 are first found. If no q-prime
* factor is present and 1 is returned, then n is prime. However if at
* least one q-prime factor is present and 1 is returned, then n retains its
* q-prime status. If 0 is returned, then n is either composite or likely
* to be composite.
*/
{
	MPI *S, *T, *M, *N;
	unsigned int i, i_;
	unsigned long x;

	i = j;
	i_ = j_;
	M = SUB0_I(Nptr, (USL)1);
	if (VERBOSE)
	{
		printf("SELFRIDGE'S EXPONENT TEST IN PROGRESS FOR ");
		PRINTI(Nptr);
		printf("\n");
	}
	*uptr = PRIME_FACTORS(M);
	while (i <= j - 1)
	{
		for (x = 2; x < (USL)R0; x++)
		{
			S = MPOWER_((long)x, M, Nptr);
			if (S->D || S->V[0] != 1)
			{
				if (VERBOSE)
				{
					printf("SELFRIDGE'S TEST IS FINISHED:\n");
					PRINTI(Nptr);
					printf(" is not a pseudo-prime to base %lu\n", x);
					printf(" and is hence composite\n");
				}
				FREEMPI(S);
				FREEMPI(M);
				return (0);
			}
			FREEMPI(S);
			N = INT0_(M, Q[i]);
			T = MPOWER_((long)x, N, Nptr);
			FREEMPI(N);
			if (T->D || T->V[0] != 1)
			{
				i++;
				FREEMPI(T);
				break;
			}
			FREEMPI(T);
		}
		if (x == R0)
		{
			if (VERBOSE)
			{
				printf("SELFRIDGE'S TEST IS FINISHED:\n");
				PRINTI(Nptr);
				printf(" is likely to be composite\n");
			}
			FREEMPI(M);
			return (0);
		}
	}
	while (i_ <= j_ - 1)
	{
		for (x = 2; x < (USL)R0; x++)
		{
			S = MPOWER_((long)x, M, Nptr);
			if (S->D || S->V[0] != 1)
			{
				if (VERBOSE)
				{
					printf("SELFRIDGE'S TEST IS FINISHED:\n");
					PRINTI(Nptr);
					printf(" is not a pseudo-prime to base %lu\n", x);
					printf(" and is hence composite\n");
				}
				FREEMPI(M);
				FREEMPI(S);
				return (0);
			}
			FREEMPI(S);
			N = INT0(M, Q_[i_]);
			T = MPOWER_((long)x, N, Nptr);
			FREEMPI(N);
			if (T->D || T->V[0] != 1)
			{
				i_++;
				FREEMPI(T);
				break;
			}
			FREEMPI(T);
		}
		if (x == R0)
		{
			if (VERBOSE)
				printf("SELFRIDGE'S TEST IS INCONCLUSIVE:\n");
			FREEMPI(M);
			return (0);
		}
	}
	FREEMPI(M);
	if (*uptr == 0)
	{
		if (VERBOSE)
		{
			printf("SELFRIDGE'S TEST IS FINISHED:\n");
			PRINTI(Nptr);
			printf(" is prime\n");
			printf("--\n");
		}
		return (1);
	}
	else
	{
		if (VERBOSE)
		{
			printf("SELFRIDGE'S TEST IS FINISHED:\n");
			PRINTI(Nptr);
			printf(" retains its q-prime status\n");
			printf("--\n");
		}
		return (1);
	}
}

unsigned int scount; /* the number of prime factors of N < 1000 */
unsigned int s_count; /* the number of q-prime factors of N */

unsigned int FACTOR(MPI *Nptr)
/*
* a factorization of *Nptr into prime and q-prime factors is first obtained.
* Selfridge's primality test is then applied to any q-prime factors; the test
* is applied repeatedly until either a q-prime is found to be composite or
* likely to be composite (in which case the initial factorization is doubtful
* and an extra base should be used in Miller's test) or we run out of q-primes,
* in which case every q-prime factor of *Nptr is certified as prime.
*/
{
	unsigned int c, i, t, u, r, count;

	j = j_ = ltotal = 0;
	c = PRIME_FACTORS(Nptr);
	scount = j;
	s_count = j_;
	if (c == 0)
	{
		if (VERBOSE)
		{
			printf("NO Q-PRIMES:\n\n");
			PRINTI(Nptr);
			printf(" has the following factorization:\n\n");
			for (i = 0; i < scount; i++)
				printf("prime factor: %lu exponent: %u\n", Q[i], K[i]);
			for (i = 0; i < s_count; i++)
			{
				printf("prime factor: ");
				PRINTI(Q_[i]);
				printf(" exponent: %u\n", K_[i]);
			}
		}
		count = scount + s_count;
	}
	else
	{
		if (VERBOSE)
		{
			printf("TESTING Q-PRIMES FOR PRIMALITY\n");
			printf("--\n");
		}
		i = 0;
		while (c)
		{
			t = SELFRIDGE(Q_PRIME[i], &u);
			FREEMPI(Q_PRIME[i]);
			if (t == 0)
			{
				printf("do FACTOR() again with an extra base\n");
				for (r = i + 1; r < ltotal; r++)
					FREEMPI(Q_PRIME[r]);
				for (i = 0; i < j_; i++)
					FREEMPI(Q_[i]);
				return (0);
			}
			c = c + u;
			i++;
			c--;
		}
		if (VERBOSE)
		{
			printf("ALL Q-PRIMES ARE PRIMES:\n\n");
			PRINTI(Nptr);
			printf(" has the following factorization:\n\n");
			for (i = 0; i < scount; i++)
				printf("prime factor: %lu exponent: %u\n", Q[i], K[i]);
			for (i = 0; i < s_count; i++)
			{
				printf("prime factor: ");
				PRINTI(Q_[i]);
				printf(" exponent: %u\n", K_[i]);
			}
		}
		count = scount + s_count;
	}
	/*printf("s_count = %u, j_ = %u\n", s_count, j_);*/
	for (i = s_count; i < j_; i++)
		FREEMPI(Q_[i]);
	if (FREE_PRIMES)
	{
		for (i = 0; i < s_count; i++)
			FREEMPI(Q_[i]);
	}
	return (count);
}

MPI *DIVISOR(MPI *N)
/*
*  Returns the number of divisors of N,
*  returns NULL if factorization unsuccessful .
*/
{
	MPI *U, *Tmp;
	unsigned int i, s;

	if (EQONEI(N))
		return (ONEI());
	U = ONEI();
	FREE_PRIMES = 1; /* This frees the Q_[i] in FACTOR()*/
	s = FACTOR(N);
	FREE_PRIMES = 0;
	if (s == 0)
	{
		FREEMPI(U);
		return ((MPI *)NULL);
	}
	for (i = 0; i < scount; i++)
	{
		Tmp = U;
		U = MULT_I(U, 1 + (USL)K[i]);
		FREEMPI(Tmp);
	}
	for (i = 0; i < s_count; i++)
	{
		Tmp = U;
		U = MULT_I(U, 1 + (USL)K_[i]);
		FREEMPI(Tmp);
	}
	return (U);
}

MPI *MOBIUS(MPI *N)
/*
*  Returns the Mobius function mu(N),
*  returns NULL if factorization unsuccessful .
*/
{
	MPI *S;
	unsigned int i, s;

	if (EQONEI(N))
		return (ONEI());
	FREE_PRIMES = 1; /* This frees the Q_[i] in FACTOR()*/
	s = FACTOR(N);
	FREE_PRIMES = 0;
	if (s == 0)
		return ((MPI *)NULL);
	for (i = 0; i < scount; i++)
		if (K[i] >= 2)
			return (ZEROI());
	for (i = 0; i < s_count; i++)
	{
		if (K_[i] >= 2)
			return (ZEROI());
	}
	if (s % 2)
		return (MINUS_ONEI());
	else
		return (ONEI());
	return (S);
}

MPI *EULER(MPI *N)
/*
*  Returns Euler's function  phi(N),
*  returns NULL if factorization unsuccessful .
*/
{
	MPI *U, *V, *W, *Tmp;
	unsigned int i, s;

	if (EQONEI(N))
		return (ONEI());
	U = ONEI();
	/* FREE_PRIMES = 0, so we can free the Q_[i] later */
	s = FACTOR(N);
	if (s == 0)
	{
		FREEMPI(U);
		return ((MPI *)NULL);
	}
	for (i = 0; i < scount; i++)
	{
		if (K[i] == 1)
			V = ONEI();
		else
			V = POWER_I((long)(Q[i]), K[i] - 1);
		W = CHANGE(Q[i] - 1);
		Tmp = U;
		U = MULTI3(U, V, W);
		FREEMPI(Tmp);
		FREEMPI(V);
		FREEMPI(W);
	}
	for (i = 0; i < s_count; i++)
	{
		if (K_[i] == 1)
			V = ONEI();
		else
			V = POWERI(Q_[i], K_[i] - 1);
		W = SUB0_I(Q_[i], 1);
		Tmp = U;
		U = MULTI3(U, V, W);
		FREEMPI(Tmp);
		FREEMPI(V);
		FREEMPI(W);
	}
	for (i = 0; i < s_count; i++)
		FREEMPI(Q_[i]);
	return (U);
}

MPI *SIGMA(MPI *N)
/*
*  Returns sigma(N), the sum of the divisors of N,
*  returns NULL if factorization unsuccessful .
*/
{
	MPI *U, *V, *W, *Tmp;
	unsigned int i, s;

	if (EQONEI(N))
		return (ONEI());
	U = ONEI();
	/* FREE_PRIMES = 0, so we can free the Q_[i] later */
	s = FACTOR(N);
	if (s == 0)
	{
		FREEMPI(U);
		return ((MPI *)NULL);
	}
	for (i = 0; i < scount; i++)
	{
		V = POWER_I((long)(Q[i]), K[i] + 1);
		Tmp = V;
		V = SUB0_I(V, 1);
		FREEMPI(Tmp);
		Tmp = V;
		V = INT0_(V, Q[i] - 1);
		FREEMPI(Tmp);
		Tmp = U;
		U = MULTI(U, V);
		FREEMPI(Tmp);
		FREEMPI(V);
	}
	for (i = 0; i < s_count; i++)
	{
		V = POWERI(Q_[i], K_[i] + 1);
		Tmp = V;
		V = SUB0_I(V, 1);
		FREEMPI(Tmp);
		Tmp = V;
		W = SUB0_I(Q_[i], 1);
		V = INT0(V, W);
		FREEMPI(W);
		FREEMPI(Tmp);
		Tmp = U;
		U = MULTI(U, V);
		FREEMPI(Tmp);
		FREEMPI(V);
	}
	for (i = 0; i < s_count; i++)
		FREEMPI(Q_[i]);
	return (U);
}

MPI *LPRIMROOT(MPI *P)
/*
*  Returns  the least primitive root mod P, an odd prime;
*  returns NULL if factorization of P - 1 is unsuccessful .
*/
{
	MPI *Tmp, *QQ, *R, *RR, *T;
	unsigned int s, u, success = 1;
	long i;
	int r;

	T = ZEROI();
	QQ = SUB0_I(P, 1);
	s = FACTOR(QQ);
	if (s == 0)
		return (T);
	for (i = 2; i >= 2; i++)
	{
		RR = CHANGEI(i);
		r = JACOBIB(RR, P);
		FREEMPI(RR);
		if (r == -1)
		{
			for (u = 0; u < scount; u++)
			{
				Tmp = INT0_(QQ, Q[u]);
				R = MPOWER_(i, Tmp, P);
				FREEMPI(Tmp);
				if (EQONEI(R))
				{
					success = 0;
					FREEMPI(R);
					break;
				}
				FREEMPI(R);
			}
			if (success)
			{
				for (u = 0; u < s_count; u++)
				{
					Tmp = INT0(QQ, Q_[u]);
					R = MPOWER_(i, Tmp, P);
					FREEMPI(Tmp);
					if (EQONEI(R))
					{
						success = 0;
						FREEMPI(R);
						break;
					}
					FREEMPI(R);
				}
			}
			if (success)
			{
				FREEMPI(QQ);
				for (u = 0; u < s_count; u++)
					FREEMPI(Q_[u]);
				FREEMPI(T);
				T = CHANGE(i);
				break;
			}
			else
				success = 1;
		}
	}
	return (T);
}

MPI *ORDERP(MPI *A, MPI *P)
/*
* Returns the order of A mod P, a prime.
*/
{
	MPI *AA, *Tmp, *QQ, *M;
	unsigned int s, u, i;

	AA = MOD(A, P);
	Tmp = SUB0_I(AA, 1);
	if (Tmp->S == 0)
	{
		FREEMPI(Tmp);
		FREEMPI(AA);
		return (ONEI());
	}
	FREEMPI(Tmp);
	Tmp = ADD0_I(AA, 1);
	if (Tmp->S == 0)
	{
		FREEMPI(Tmp);
		FREEMPI(AA);
		return (CHANGE(2));
	}
	FREEMPI(Tmp);
	QQ = SUB0_I(P, 1);
	/* FREE_PRIMES = 0, so we can free the Q_[i] later */
	s = FACTOR(QQ);
	for (u = 0; u < scount; u++)
	{
		for (i = 1; i <= K[u]; i++)
		{
			Tmp = INT0_(QQ, Q[u]);
			M = MPOWER(AA, Tmp, P);
			FREEMPI(Tmp);
			if (!EQONEI(M))
			{
				FREEMPI(M);
				break;
			}
			FREEMPI(M);
			Tmp = QQ;
			QQ = INT0_(QQ, Q[u]);
			FREEMPI(Tmp);
		}
	}
	for (u = 0; u < s_count; u++)
	{
		for (i = 1; i <= K_[u]; i++)
		{
			Tmp = INT0(QQ, Q_[u]);
			M = MPOWER(AA, Tmp, P);
			FREEMPI(Tmp);
			if (!EQONEI(M))
			{
				FREEMPI(M);
				break;
			}
			FREEMPI(M);
			Tmp = QQ;
			QQ = INT0_(QQ, Q[u]);
			FREEMPI(Tmp);
		}
	}
	FREEMPI(AA);
	for (i = 0; i < s_count; i++)
		FREEMPI(Q_[i]);
	return (QQ);
}

MPI *ORDERQ(MPI *A, MPI *P, unsigned int n)
/*
* Returns the order of A mod P^n, P a prime.
*/
{
	MPI *D, *E, *Tmp, *Tmp1, *Q;
	unsigned int h;

	if (EQONEI(A))
		return (ONEI());
	if (EQMINUSONEI(A))
		return (CHANGE(2));
	if (n == 1)
		return (ORDERP(A, P));
	if (P->D == 0 && P->V[0] == 2)
	{
		if ((A->V[0] - 1) % 4 == 0)
			D = ONEI();
		/*if ((A->V[0] + 1) % 4 == 0)*/
		else
			D = CHANGE(2);
	}
	else
		D = ORDERP(A, P);
	Q = COPYI(P);
	E = ZEROI();
	h = 0;
	while (E->S == 0)
	{
		Tmp = Q;
		Q = MULTI(Q, P);
		FREEMPI(Tmp);
		Tmp = E;
		E = MPOWER(A, D, Q);
		FREEMPI(Tmp);
		Tmp = E;
		E = SUB0_I(E, 1);
		FREEMPI(Tmp);
		h++;
	}
	FREEMPI(E);
	FREEMPI(Q);
	if (n <= h)
		return (D);
	else
	{
		Tmp = POWERI(P, n - h);
		Tmp1 = MULTI(Tmp, D);
		FREEMPI(Tmp);
		FREEMPI(D);
		return (Tmp1);
	}
}

MPI *ORDERM(MPI *A, MPI *M)
/*
* Returns the order of A mod M.
*/
{
	MPI *O, **Y, *Tmp, *Tmp1;
	unsigned int i, s, *x;

	if (EQONEI(A))
		return (ONEI());
	if (EQMINUSONEI(A))
	{
		if (M->D == 0 && M->V[0] == 2)
			return (ONEI());
		else
			return (CHANGE(2));
	}
	/* FREE_PRIMES = 0, so we can free the Q_[i] later */
	s = FACTOR(M);
	x = (USI *)mmalloc(s * sizeof(USI));
	Y = (MPI **)mmalloc(s * sizeof(MPI *));
	for (i = 0; i < scount; i++)
	{
		x[i] = K[i];
		Y[i] = CHANGE(Q[i]);
	}
	for (i = 0; i < s_count; i++)/* bugfix here - ANU,27/9/94. */
	{
		x[i + scount] = K_[i];
		Y[i + scount] = COPYI(Q_[i]);
		FREEMPI(Q_[i]);
	}
	O = ORDERQ(A, Y[0], x[0]);
	for (i = 1; i < s; i++)
	{
		Tmp = O;
		Tmp1 = ORDERQ(A, Y[i], x[i]);
		O = LCM(O, Tmp1);
		FREEMPI(Tmp);
		FREEMPI(Tmp1);
	}
	ffree(x, s * sizeof(USI));
	for (i = 0; i < s; i++)
		FREEMPI(Y[i]);
	ffree(Y, s * sizeof(MPI *));
	return (O);
}

void FERMAT_QUOTIENT(MPI *N)
{
	USI i, n;
	MPI *P, *Q, *O, *T1, *T2, *T3, *T4, *T33, *T44;
	MPI *PP, *A, *B, *R;

	n = CONVERTI(N);
	O = ONEI();
	for (i = 2; i <= n; i++)
	{
		printf("i = %u\n", i);
		R = CHANGE(PRIME[i]);
		P = CHANGE(PRIME[i] + 1);
		Q = CHANGE(PRIME[i] - 1);
		T1 = POWERI(P, (USI)PRIME[i]);
		FREEMPI(P);
		T2 = POWERI(Q, (USI)PRIME[i]);
		FREEMPI(Q);
		T3 = SUB0I(T1, O);
		T4 = ADD0I(T2, O);
		FREEMPI(T1);
		FREEMPI(T2);
		PP = MULTI(R, R);
		FREEMPI(R);
		T33 = INT0(T3, PP);
		T44 = INT0(T4, PP);
		FREEMPI(T3);
		FREEMPI(T4);
		FREEMPI(PP);
		A = LUCAS(T33);
		FREEMPI(T33);
		if (A->S)
			printf("p = %lu gives a - prime quotient\n", PRIME[i]);
		FREEMPI(A);
		B = LUCAS(T44);
		FREEMPI(T44);
		if (B->S)
			printf("p = %lu gives a + prime quotient\n", PRIME[i]);
		FREEMPI(B);
	}
	FREEMPI(O);
	return;
}

MPI *SQROOT1(MPI *A, MPI *P, USL n)
{
	/* Solving x^2=A (mod P^n), P an odd prime not dividing A. */
	/* returns -1 if not soluble, otherwise returns x. */

	MPI *T, *K, *U, *V, *Tmp1, *Tmp2, *N;
	USI i, t;

	Tmp1 = INT_(P, 2);
	T = MPOWER(A, Tmp1, P);
	t = EQONEI(T);
	FREEMPI(T);
	FREEMPI(Tmp1);
	if (t == 0)
		return(MINUS_ONEI());
	if (n == 1)
		return(SQRTM(A, P));
	else {
		U = SQROOT1(A, P, n - 1);
		Tmp1 = MULTI(U, U);
		V = SUBI(Tmp1, A);
		FREEMPI(Tmp1);
		for (i = 0; i < n - 1; i++) {
			Tmp1 = V;
			V = INT(V, P);
			FREEMPI(Tmp1);
		}
		Tmp1 = MULT_I(U, 2);
		Tmp2 = MINUSI(V);
		FREEMPI(V);
		K = CONGR(Tmp1, Tmp2, P, &N);
		FREEMPI(Tmp1);
		FREEMPI(Tmp2);
		FREEMPI(N);
		Tmp1 = POWERI(P, (USI)(n - 1));
		Tmp2 = MULTI(K, Tmp1);
		FREEMPI(K);
		FREEMPI(Tmp1);
		Tmp1 = ADDI(U, Tmp2);
		FREEMPI(U);
		FREEMPI(Tmp2);
		return(Tmp1);
	}
}

void SQROOT1_W(Stack S)
{
	USL n;
	MPI *X;

	MPI *A = stackPop(S);
	MPI *P = stackPop(S);
	MPI *N = stackPop(S);
	n = CONVERTI(N);
	X = SQROOT1(A, P, n);
	printf("SQROOT of ");
	PRINTI(A);
	printf(" (mod ");
	PRINTI(P);
	printf("^%lu) is ", n);
	PRINTI(X);
	FREEMPI(X);
	FREEMPI(A);
	FREEMPI(P);
	FREEMPI(N);
	printf("\n");
	return;
}

USI sqroot2_number; /* global variable, used in SQROOT */

MPI *SQROOT2(MPI *A, USL n)
{
	/* Solving x^2=A (mod 2^n), A odd.
	* returns -1 if not soluble, otherwise returns x.
	* In the case of solubility, changes global variable sqroot2_number to 1,2,4
	* according as n= 1, 2  or > 2.
	*/

	MPI *U, *V, *Tmp1, *Tmp2;
	USI i, t;
	USL s;

	if (n == 1) {
		sqroot2_number = 1;
		return(ONEI());
	}
	if (n == 2) {
		s = MOD_(A, 4);
		if (s != 1)
			return(MINUS_ONEI());
		else {
			sqroot2_number = 2;
			return(ONEI());
		}
	}
	s = MOD_(A, 8);
	if (s != 1)
		return(MINUS_ONEI());
	if (n == 3) {
		sqroot2_number = 4;
		return(ONEI());
	}
	else {
		U = SQROOT2(A, n - 1);
		Tmp1 = MULTI(U, U);
		V = SUBI(Tmp1, A);
		FREEMPI(Tmp1);
		for (i = 0; i < n - 1; i++) {
			Tmp1 = V;
			V = INT_(V, 2);
			FREEMPI(Tmp1);
		}
		Tmp2 = ONEI();
		Tmp1 = ADDI(V, Tmp2);
		FREEMPI(Tmp2);
		FREEMPI(V);
		t = (Tmp1->V[0]) % 2;
		FREEMPI(Tmp1);
		if (t == 0) {
			Tmp1 = POWER_I(2, (USI)(n - 2));
			Tmp2 = U;
			U = ADDI(U, Tmp1);
			FREEMPI(Tmp1);
			FREEMPI(Tmp2);
		}
		return(U);
	}

}

void SQROOT2_W(Stack S)
{
	USL n;
	MPI *X;

	MPI *A = stackPop(S);
	MPI *N = stackPop(S);
	n = CONVERTI(N);
	X = SQROOT2(A, n);
	printf("SQROOT of ");
	PRINTI(A);
	printf(" (mod 2");
	printf("^%lu) is ", n);
	PRINTI(X);
	FREEMPI(X);
	FREEMPI(A);
	FREEMPI(N);
	printf("\n");
	return;
}

MPI *SQROOT3(MPI *A, MPI *P, USL n, MPI**EXPONENT)
/* Solving x^2=A (mod P^n), P a prime possibly dividing A.
* Returns -1 if not soluble, otherwise returns x.
* Also returns a "reduced moduli" explained as follows:
* If P does not divide A, the story is well-known.
* If P divides A, there are two cases:
* (i) P^n divides A,
* (ii) P^n does not divide A.
* In case (i) there are two cases:
* (a) n=2m+1, (b) n=2m.
* In case (a), the solution is x=0 (mod P^(m+1)).
* In case (b), the solution is x=0 (mod P^m).
* (These are called reduced moduli.)
* In case (i) gcd(A,P^n)=P^r, r<n. If r is odd, no solution.
* If r=2m, write x=(P^m)*X, A=(p^2m)*A1, P not dividing A1.
* Then x^2=A (mod P^n) becomes X^2=A1 (mod P^(n-2m)).
* If P is odd, this has 2 solutions X mod P^(n-2m) and we get two solutions
* x (mod P^(n-m)). If P=2, we get
*                            1 solution mod (2^(n-m)) if n-2m=1,
*                            2 solutions mod (2^(n-m) if n-2m=2,
*                            4 solutions mod (2^(n-m)) if n-2m>2.
*/
{
	MPI *D, *U, *V, *AA, *Tmp1, *T;
	USI r, m, s;
	int t;

	/* The next two lines are important. eg. a call sqroot(431,10,&a[],&m,&l)
	sets global variable sqroot2_number=1, whereas a subsequent call to
	sqroot(46,210,&a[],&m,&l) didn't call SQROOT2() and before the next fix
	was added, sqroot2_number remained equal to 1, causing problems as
	action is taken in SQROOT() if S[0] = 2 only if sqroot2_number > 0. */
	if ((P->D == 0) && P->V[0] == 2)
		sqroot2_number = 0;
	*EXPONENT = NULL;
	D = MOD(A, P);
	t = D->S;
	FREEMPI(D);
	if (t) {
		*EXPONENT = POWERI(P, n);
		if ((P->D == 0) && P->V[0] == 2)
			return(SQROOT2(A, n));
		else
			return(SQROOT1(A, P, n));
	}
	else {
		U = POWERI(P, n);
		V = MOD(A, U);
		t = V->S;
		FREEMPI(U);
		FREEMPI(V);
		if (t == 0) {
			if ((n % 2) == 0)
				*EXPONENT = POWERI(P, n / 2);
			else
				*EXPONENT = POWERI(P, 1 + n / 2);
			return(ZEROI());
		}
		r = 0;
		AA = COPYI(A);
		V = MOD(AA, P);
		while ((V->S == 0) && (r <= n)) {
			Tmp1 = AA;
			AA = INT(AA, P);
			FREEMPI(Tmp1);
			Tmp1 = V;
			V = MOD(AA, P);
			FREEMPI(Tmp1);
			r++;
		}
		FREEMPI(V);
		if (r % 2) {
			FREEMPI(AA);
			return(MINUS_ONEI());
		}
		else {
			m = r / 2;
			s = n - r;
			*EXPONENT = POWERI(P, n - m);
			if ((P->D == 0) && (P->V[0] == 2)) {
				T = SQROOT2(AA, s);
				if (EQMINUSONEI(T)) {
					FREEMPI(AA);
					return(T);
				}
				else {
					U = POWERI(P, m);
					V = MULTI(U, T);
					FREEMPI(AA);
					FREEMPI(T);
					FREEMPI(U);
					return(V);
				}
			}
			else {
				T = SQROOT1(AA, P, s);
				if (EQMINUSONEI(T)) {
					FREEMPI(AA);
					return(T);
				}
				else {
					U = POWERI(P, m);
					V = MULTI(U, T);
					FREEMPI(AA);
					FREEMPI(T);
					FREEMPI(U);
					return(V);
				}

			}
		}
	}
}

void SQROOT3_W(Stack S)
{
	USL n;
	MPI *X, *EXPONENT;

	MPI *A = stackPop(S);
	MPI *P = stackPop(S);
	MPI *N = stackPop(S);
	n = CONVERTI(N);
	X = SQROOT3(A, P, n, &EXPONENT);
	printf("SQROOT of ");
	PRINTI(A);
	printf(" (mod ");
	PRINTI(P);
	printf("^%lu) is ", n);
	PRINTI(X);
	FREEMPI(X);
	printf(" (mod ");
	PRINTI(EXPONENT);
	FREEMPI(A);
	FREEMPI(P);
	FREEMPI(N);
	FREEMPI(EXPONENT);
	printf(")\n");
	return;
}

MPI *SQROOT(MPI *A, MPI *N, MPIA *SOL, MPI **MODULUS, USI *lptr)
/* Returns an array SOL consisting of the solutions of x^2=A (mod N),
* with 0<=x<=N/2 where N = prod QQ[i]^KK[i].
* Returns the number of solutions (mod N) and their MODULUS if soluble,
* otherwise returns -1.
* The algorithm uses the Chinese remainder theorem after solving mod
* QQ[i]^KK[i], i=1,...,e=omega(N). Note N is even if and only if QQ[0]=2.
* Finally we join all solutions together using the CRT.
* Some care is needed to exhibit all solutions. We operate on the D[]
* array below, as follows:
* l= number of primes QQ[i] dividing N, for which P^a=QQ[i]^KK[i] does not
* divide N, while jj is the remaining number.
* D[] is the array formed by a single solution each of x^2=A mod P^a
* SS is the product of the reduced moduli for these l primes.
* OM is the product of the reduced moduli for the jj primes.
* S[] is the array formed by the moduli P^a.
* If l>1, we produce all 2^l l-tuples (e[0]*D[0],...,e[l-1]*D[l-1]),
* where e[i]=1 or -1, if s is odd or D[0]=2 and "sqroot2_number">1, but only
* 2^(l-1) * (l-1)-tuples (e[0]*D[0],...,e[l-2]*D[l-2]), if SS is even and
* "sqroot2_number"=1.
* It is convenient to use the CRT with D[l]=om, S[l]=0, even when OM=1,
* so as to produce all solutions.
* If SS is even, we place D[0] at the end of the D[] array, making it D[l-1].
* This is done so as to conveniently operate on the initial part of the D[]
* array, where the initial S[i] are only odd prime powers.
* The program aborts if the number of prime factors of N which do not divide A
* is > 31. (Not an easy program to get right. Finished 12th April 2000.)
*/
{
	USI l, e, jj, i, r, j, k, t, *KK, rr, ii;
	MPI *C, *OM, **D, **D1, **DD, **DD1, **S, *E, *X, *AA, *B, *Y, *SS;
	MPI *Tmp1, *Tmp2, *Tmp3, *SO, *Z, **QQ, *F, *X1;
	int tt;

	if (EQONEI(N)) {
		*SOL = BUILDMPIA();
		Y = ZEROI();
		ADD_TO_MPIA(*SOL, Y, 0);
		FREEMPI(Y);
		*MODULUS = ONEI();
		*lptr = 1;
		return(ONEI());
	}
	e = FACTOR(N);

	QQ = (MPI **)mmalloc(e * sizeof(MPI *));
	KK = (USI *)mmalloc(e * sizeof(USI));
	for (i = 0; i< scount; i++) {
		QQ[i] = CHANGE(Q[i]);
		KK[i] = K[i];
	}
	for (i = 0; i< s_count; i++) {
		QQ[i + scount] = Q_[i];
		KK[i + scount] = K_[i];
	}
	l = 0;
	jj = 0;
	OM = ONEI();
	SS = ONEI();
	D = (MPI **)mmalloc(e * sizeof(MPI *));
	S = (MPI **)mmalloc((e + 1) * sizeof(MPI *));
	for (i = 0; i < e; i++) {
		C = SQROOT3(A, QQ[i], KK[i], &E);
		t = EQMINUSONEI(C);
		tt = C->S;
		if (t) {
			FREEMPI(OM);
			FREEMPI(E);
			FREEMPI(SS);
			FREEMPI(C);
			for (ii = 0; ii < e; ii++)
				FREEMPI(QQ[ii]);
			ffree((char *)D, e * sizeof(MPI *));
			ffree((char *)S, (e + 1) * sizeof(MPI *));
			ffree((char *)QQ, e * sizeof(MPI *));
			ffree((char *)KK, e * sizeof(USI));
			*SOL = BUILDMPIA();
			ADD_TO_MPIA(*SOL, NULL, 0);
			*MODULUS = (MPI *)NULL;
			*lptr = 0;
			return(MINUS_ONEI());
		}
		else if (tt == 0) {
			Tmp1 = OM;
			OM = MULTI(OM, E);
			FREEMPI(Tmp1);
			FREEMPI(E); /* added 10/12.00 */
			FREEMPI(C);
		}
		else {
			D[l] = C;
			S[l] = E;
			Tmp1 = SS;
			SS = MULTI(SS, S[l]);
			FREEMPI(Tmp1);
			l++;
		}
	}
	for (i = 0; i < e; i++)
		FREEMPI(QQ[i]);
	ffree((char *)QQ, e * sizeof(MPI *));
	ffree((char *)KK, e * sizeof(USI));
	if (l)
		SO = MULTI(S[0], OM);
	else
		SO = NULL;
	if (l == 0) {
		X = INT(N, OM);
		FREEMPI(SS);
		ffree((char *)D, e * sizeof(MPI *));
		ffree((char *)S, (e + 1) * sizeof(MPI *));
		*SOL = BUILDMPIA();
		Y = ZEROI();
		ADD_TO_MPIA(*SOL, Y, 0);
		FREEMPI(Y);
		*MODULUS = OM;
		*lptr = 1;
		return(X);
	}
	if (l == 1) {
		if (!EQONEI(OM)) {
			Tmp1 = ZEROI();
			X = CHINESE(D[0], Tmp1, S[0], OM, &Tmp2);
			FREEMPI(Tmp2);
			FREEMPI(Tmp1);
		}
		else
			X = COPYI(D[0]);
		if ((S[0]->V[0]) % 2) {
			Tmp1 = MULT_I(N, 2);
			Y = INT(Tmp1, SO);
			FREEMPI(Tmp1);
			FREEMPI(SS);
			FREEMPI(OM);
			FREEMPI(D[0]);
			ffree((char *)D, e * sizeof(MPI *));
			for (i = 0; i<l; i++) /* changed 4/1/01 */
				FREEMPI(S[i]);
			ffree((char *)S, (e + 1) * sizeof(MPI *));
			*SOL = BUILDMPIA();
			ADD_TO_MPIA(*SOL, X, 0);
			FREEMPI(X);
			*MODULUS = SO;
			*lptr = 1;
			return(Y);
		}
		else {/* S[0] is even */
			t = sqroot2_number;
			if (t == 4) {
				Tmp2 = INT_(S[0], 2);
				Tmp3 = ADDI(D[0], Tmp2);
				FREEMPI(Tmp2);
				Y = MOD(Tmp3, S[0]);
				FREEMPI(Tmp3);
				if (!EQONEI(OM)) {
					Tmp1 = ZEROI();
					X1 = CHINESE(Y, Tmp1, S[0], OM, &Tmp2);
					FREEMPI(Tmp1);
					FREEMPI(Tmp2);
				}
				else
					X1 = COPYI(Y);
				FREEMPI(Y);
				FREEMPI(D[0]);
				FREEMPI(SS);
				FREEMPI(OM);
				for (i = 0; i<l; i++)/* 20/12/00 - had e not l */
					FREEMPI(S[i]);
				ffree((char *)D, e * sizeof(MPI *));
				ffree((char *)S, (e + 1) * sizeof(MPI *));
				Y = MULT_I(N, 4);
				Z = INT0(Y, SO);
				FREEMPI(Y);
				*SOL = BUILDMPIA();
				ADD_TO_MPIA(*SOL, X, 0);
				ADD_TO_MPIA(*SOL, X1, 1);
				*MODULUS = SO;
				FREEMPI(X);
				FREEMPI(X1);
				*lptr = 2; /* bugfix: used to have =1, 4/Jul/00 */
				return(Z);
			}
			else if (t == 2) {
				Y = MULT_I(N, 2);
				Z = INT0(Y, SO);
				FREEMPI(D[0]);
				FREEMPI(SS);
				FREEMPI(OM);
				FREEMPI(Y);
				ffree((char *)D, e * sizeof(MPI *));
				for (i = 0; i<l; i++)/* 20/12/00 - had e not l */
					FREEMPI(S[i]);
				ffree((char *)S, (e + 1) * sizeof(MPI *));
				*SOL = BUILDMPIA();
				ADD_TO_MPIA(*SOL, X, 0);
				FREEMPI(X);
				*MODULUS = SO;
				*lptr = 1;
				return(Z);
			}
			else {
				Z = INT0(N, SO);
				FREEMPI(D[0]);
				FREEMPI(SS);
				FREEMPI(OM);
				ffree((char *)D, e * sizeof(MPI *));
				for (i = 0; i<e; i++)
					FREEMPI(S[i]);
				ffree((char *)S, (e + 1) * sizeof(MPI *));
				*SOL = BUILDMPIA();
				ADD_TO_MPIA(*SOL, X, 0);
				FREEMPI(X);
				*MODULUS = SO;
				*lptr = 1;
				return(Z);
			}
		}
	}
	if (l > 31)
		execerror("l in SQROOT >31", "");
	F = POWER_I(2, l);
	k = (USI)CONVERTI(F);
	D1 = (MPI **)mmalloc((l + 1) * sizeof(MPI *));
	S[l] = OM;
	D1[l] = ZEROI();
	FREEMPI(SO);
	SO = MULTI(SS, OM);
	FREEMPI(SS);
	rr = 0;
	*SOL = BUILDMPIA();
	if ((N->V[0]) % 2) {
		for (i = 0; i < k; i++) {
			j = i;
			for (r = 0; r < l; r++) {
				t = j % 2;
				if (t)
					D1[r] = MINUSI(D[r]);
				else
					D1[r] = COPYI(D[r]);
				j = j / 2;
			}
			X = CHINESE_ARRAY(D1, S, &Tmp2, l + 1);
			for (r = 0; r < l; r++)
				FREEMPI(D1[r]);
			Tmp1 = MULT_I(X, 2);
			if (RSV(Tmp1, Tmp2) <= 0) {
				ADD_TO_MPIA(*SOL, X, rr);
				rr++;
			}
			FREEMPI(X);
			FREEMPI(Tmp1);
			FREEMPI(Tmp2);
		}
		for (i = 0; i < l; i++) {
			FREEMPI(D[i]);
			FREEMPI(S[i]);
		}
		FREEMPI(D1[l]);
		FREEMPI(S[l]);
		ffree((char *)D, e * sizeof(MPI *));
		ffree((char *)S, (e + 1) * sizeof(MPI *));
		ffree((char *)D1, (l + 1) * sizeof(MPI *));
		Tmp1 = MULTI(F, N);
		FREEMPI(F);
		Z = INT0(Tmp1, SO);
		FREEMPI(Tmp1);
		*MODULUS = SO;
		*lptr = rr;
		return(Z);
	}
	else {/* N is even */
		DD = (MPI **)mmalloc((l + 1) * sizeof(MPI *));
		DD1 = (MPI **)mmalloc((l + 1) * sizeof(MPI *));
		if (sqroot2_number == 4)
			DD1[l] = ZEROI();
		AA = D[0];
		B = S[0];
		for (r = 1; r < l; r++) {
			D[r - 1] = D[r];
			S[r - 1] = S[r];
		}
		D[l - 1] = AA;
		S[l - 1] = B;
		if (sqroot2_number == 4) {
			for (r = 0; r < l - 1; r++)
				DD[r] = COPYI(D[r]);
			Tmp2 = INT_(S[l - 1], 2);
			Tmp3 = ADDI(D[l - 1], Tmp2);
			FREEMPI(Tmp2);
			DD[l - 1] = MOD(Tmp3, S[l - 1]);
			FREEMPI(Tmp3);
		}
		if (sqroot2_number == 1)
			k = k / 2;

		for (i = 0; i < k; i++) {
			j = i;
			for (r = 0; r < l; r++) {
				t = j % 2;
				if (t) {
					D1[r] = MINUSI(D[r]);
					if (sqroot2_number == 4)
						DD1[r] = MINUSI(DD[r]);
				}
				else {
					D1[r] = COPYI(D[r]);
					if (sqroot2_number == 4)
						DD1[r] = COPYI(DD[r]);
				}
				j = j / 2;
			}
			X = CHINESE_ARRAY(D1, S, &Tmp2, l + 1);
			Tmp1 = MULT_I(X, 2);
			if (RSV(Tmp1, Tmp2) <= 0) {
				ADD_TO_MPIA(*SOL, X, rr);
				rr++;
			}
			FREEMPI(Tmp1);
			FREEMPI(X);
			FREEMPI(Tmp2);
			if (sqroot2_number == 4) {
				X = CHINESE_ARRAY(DD1, S, &Tmp2, l + 1);
				Tmp1 = MULT_I(X, 2);
				if (RSV(Tmp1, Tmp2) <= 0) {
					ADD_TO_MPIA(*SOL, X, rr);
					rr++;
				}
				FREEMPI(X);
				FREEMPI(Tmp1);
				FREEMPI(Tmp2);
			}
			for (r = 0; r < l; r++) {
				FREEMPI(D1[r]);
				if (sqroot2_number == 4)
					FREEMPI(DD1[r]);
			}
		}
		for (i = 0; i < l; i++) {
			FREEMPI(D[i]);
			FREEMPI(S[i]);
		}
		FREEMPI(D1[l]);
		if (sqroot2_number == 4) {
			FREEMPI(DD1[l]);
			for (i = 0; i < l; i++)
				FREEMPI(DD[i]);
		}
		FREEMPI(S[l]);
		ffree((char *)D, e * sizeof(MPI *));
		ffree((char *)S, (e + 1) * sizeof(MPI *));
		ffree((char *)D1, (l + 1) * sizeof(MPI *));
		ffree((char *)DD1, (l + 1) * sizeof(MPI *));
		ffree((char *)DD, (l + 1) * sizeof(MPI *));
		Tmp1 = MULTI(F, N);
		if (sqroot2_number == 1) {
			Tmp2 = Tmp1;
			Tmp1 = INT0_(Tmp1, 2);
			FREEMPI(Tmp2);
		}
		if (sqroot2_number == 4) {
			Tmp2 = Tmp1;
			Tmp1 = MULT_I(Tmp1, 2);
			FREEMPI(Tmp2);
		}
		FREEMPI(F);
		Z = INT0(Tmp1, SO);
		FREEMPI(Tmp1);
		*MODULUS = SO;
		*lptr = rr;
		return(Z);
	}
}

MPI *SQROOTX(MPI *A, MPI *N, MPIA *Y, MPI **M, USI *l)
{
	MPI *tmp;

	if (N->S <= 0) {
		printf("N <= 0\n");
		return(NULL);
	}
	else
	{
		tmp = SQROOT(A, N, Y, M, l);
		return(tmp);
	}
}

void CORNACCHIA(MPI *A, MPI *B, MPI *M)
/* Cornacchia's algorithm. See L'algorithme de Cornacchia, A. Nitaj,
* Expositiones Mathematicae, 358-365.
*/
{
	MPI *A1, *A2, *Tmp1, *TT, *T1, *T2, *N;
	MPI *AA, *BB, *R, *Q1, *QQ, *MA, *MB;
	MPIA ARRAY;
	USI ll, i, jj;
	int t1, t2;

	A1 = INVERSEM(A, M);
	Tmp1 = MINUSI(B);
	A2 = MULTI(A1, Tmp1);
	FREEMPI(A1);
	FREEMPI(Tmp1);
	MA = INT0(M, A);
	MB = INT0(M, B);
	N = SQROOTX(A2, M, &ARRAY, &Tmp1, &ll);
	FREEMPI(N);
	FREEMPI(A2);
	FREEMPI(Tmp1);
	for (i = 0; i < ll; i++) {
		TT = ZEROI();
		T1 = ONEI();
		AA = COPYI(M);
		BB = COPYI(ARRAY->A[i]);
		t2 = 0;
		jj = 0;
		R = COPYI(AA);
		T2 = ZEROI();
		while (t2 <= 0) {
			Tmp1 = MULTI(R, R);
			t1 = RSV(Tmp1, MA);
			FREEMPI(Tmp1);
			if (t1 <= 0) {
				printf("X=");
				PRINTI(R);
				printf(", Y=");
				if (T2->S == -1)
					T2->S = 1;
				PRINTI(T2);
				printf("\n");
				break;
			}
			jj++;
			if (jj == 1) {
				FREEMPI(R);
				R = COPYI(BB);
				FREEMPI(T2);
				T2 = ONEI();
				t2 = RSV(T2, MB);
			}
			if (jj > 1) {
				if (jj == 2) {
					FREEMPI(R);
					FREEMPI(T2);
				}
				R = MOD0(AA, BB);
				Q1 = INT0(AA, BB);
				FREEMPI(AA);
				AA = BB;
				BB = R;
				QQ = MINUSI(Q1);
				FREEMPI(Q1);
				T2 = MULTABC(TT, QQ, T1);
				FREEMPI(TT);
				FREEMPI(QQ);
				Tmp1 = MULTI(T2, T2);
				t2 = RSV(Tmp1, MB);
				FREEMPI(Tmp1);
				TT = T1;
				T1 = T2;
			}
		}
		FREEMPI(TT);
		FREEMPI(T1);
		if (jj == 1) {
			FREEMPI(T2);
			FREEMPI(R);
		}
		FREEMPI(AA);
		FREEMPI(BB);
	}
	FREEMPIA(ARRAY);
	FREEMPI(MA);
	FREEMPI(MB);
}


MPI *CORNACCHIAX(MPI *A, MPI *B, MPI *M)
{
	MPI *Tmp;
	int t;
	USI s;

	if (A->S <= 0)
	{
		printf("A <= 0\n");
		return NULL;
	}
	if (B->S <= 0)
	{
		printf("B <= 0\n");
		return NULL;
	}
	Tmp = ADD0I(A, B);
	t = RSV(Tmp, M);
	FREEMPI(Tmp);
	if (t > 0)
	{
		printf("M < A + B\n");
		return NULL;
	}
	Tmp = GCD(A, B);
	s = EQONEI(Tmp);
	FREEMPI(Tmp);
	if (s == 0)
	{
		printf("gcd(A,B) > 1\n");
		return NULL;
	}
	Tmp = GCD(A, M);
	s = EQONEI(Tmp);
	FREEMPI(Tmp);
	if (s == 0)
	{
		printf("gcd(A,M) > 1\n");
		return NULL;
	}

	CORNACCHIA(A, B, M);
	return(ONEI());
}

void GAUSS(MPI *A, MPI *B, MPI *C, MPI *N, MPI **alpha, MPI **gamma, MPI **M)
/* input: gcd(A,B,C)=1, |N|>1.
* output: (alpha,gamma), where A*alpha^2+B*alpha*gamma+C*gamma^2=M, gcd(M,N)=1.
*/
{
	MPI *ABSN, **QQ, *TMP1, *TMP2, *TMP3, *TMP4;
	MPI **XX, **YY, *MM, *NN, *G;
	USI e, i;
	int s, t;

	ABSN = ABSI(N);
	e = FACTOR(ABSN);
	FREEMPI(ABSN);
	QQ = (MPI **)mmalloc(e * sizeof(MPI *));
	for (i = 0; i< scount; i++)
		QQ[i] = CHANGE(Q[i]);
	for (i = 0; i< s_count; i++)
		QQ[i + scount] = Q_[i];
	XX = (MPI **)mmalloc(e * sizeof(MPI *));
	YY = (MPI **)mmalloc(e * sizeof(MPI *));
	for (i = 0; i < e; i++) {
		TMP1 = MOD(A, QQ[i]);
		TMP2 = MOD(C, QQ[i]);
		s = TMP1->S;
		t = TMP2->S;
		FREEMPI(TMP1);
		FREEMPI(TMP2);
		if (s) {
			XX[i] = ONEI();
			YY[i] = ZEROI();
		}
		if (!s && t) {
			XX[i] = ZEROI();
			YY[i] = ONEI();
		}
		if (!s && !t) {
			XX[i] = ONEI();
			YY[i] = ONEI();
		}
	}
	TMP1 = CHINESE_ARRAY(XX, QQ, &MM, e);
	TMP2 = CHINESE_ARRAY(YY, QQ, &NN, e);
	FREEMPI(MM);
	FREEMPI(NN);
	G = GCD(TMP1, TMP2);
	*alpha = INT(TMP1, G);
	FREEMPI(TMP1);
	*gamma = INT(TMP2, G);
	FREEMPI(TMP2);
	FREEMPI(G);
	TMP1 = MULTI3(A, *alpha, *alpha);
	TMP2 = MULTI3(B, *alpha, *gamma);
	TMP3 = MULTI3(C, *gamma, *gamma);
	TMP4 = ADDI(TMP1, TMP2);
	*M = ADDI(TMP4, TMP3);
	FREEMPI(TMP1);
	FREEMPI(TMP2);
	FREEMPI(TMP3);
	FREEMPI(TMP4);
	for (i = 0; i < e; i++) {
		FREEMPI(QQ[i]);
		FREEMPI(XX[i]);
		FREEMPI(YY[i]);
	}
	ffree((char *)QQ, e * sizeof(MPI *));
	ffree((char *)XX, e * sizeof(MPI *));
	ffree((char *)YY, e * sizeof(MPI *));
	return;
}

/*
#include <stdio.h>
#include "unistd.h"
*/

MPI *PRIME_GENERATOR(MPI *M, MPI *N)
/* lists the primes P satisfying M <= P <= N.
* Returns the number of such primes.
*/
{
	unsigned long j, k;
	MPI *C, *P, *Temp, *Temp1, *Temp2, *MM, *MMM, *NN, *ONE;
	unsigned int c = 0;
	int t;
	char buff[20];
	FILE *outfile;

	ONE = ONEI();
	strcpy(buff, "primes.out");
	if (access(buff, R_OK) == 0)
		unlink(buff);
	outfile = fopen(buff, "w");

	if (RSV(M, N) == 1) {
		Temp = COPYI(N);
		NN = COPYI(M);
		MM = Temp;
	}
	else {
		MM = COPYI(M);
		NN = COPYI(N);
	}
	if (MM->D == 0 && MM->V[0] == 2) {
		c++;
		printf("2\n");
		fprintf(outfile, "2\n");
	}
	t = (MM->V[0]) % 2;
	MMM = COPYI(MM);
	if (t == 0) {
		Temp = MM;
		MM = ADD0I(MM, ONE);
		FREEMPI(Temp);
	}
	FREEMPI(ONE);
	P = COPYI(MM);
	while (RSV(P, NN) <= 0) {
		j = 1;
		k = 1;
		while (1) {
			Temp2 = CHANGE(PRIMES[k]);
			Temp1 = MULT_I(Temp2, (long)PRIMES[k]);
			FREEMPI(Temp2);
			t = RSV(P, Temp1);
			FREEMPI(Temp1);
			if (t < 0) {
				break;
			}
			if (MOD0_(P, PRIMES[k]) == 0) {
				j = 0;
				break;
			}
			if (k >= 9591)
				break;
			k++;
		}
		if (j == 1) {
			c++;
			PRINTI(P); printf("\n");
			FPRINTI(outfile, P); fprintf(outfile, "\n");
		}
		Temp = P;
		P = ADD0_I(P, (USL)2);
		FREEMPI(Temp);
	}
	FREEMPI(P);
	printf("the number of primes in the range ");
	PRINTI(MMM);
	printf(" to ");
	PRINTI(NN);
	fprintf(outfile, "the number of primes in the range ");
	FPRINTI(outfile, MMM);
	fprintf(outfile, " to ");
	FPRINTI(outfile, NN);
	printf(" is %u\n", c);
	fprintf(outfile, " is %u\n", c);
	fclose(outfile);
	C = CHANGE(c);
	FREEMPI(MM);
	FREEMPI(NN);
	FREEMPI(MMM);
	return (C);
}

MPI *PRIME_GENERATORX(MPI *M, MPI *N)
{
	MPI *BOUND, *ONE, *NUMBER;
	int s, t;

	BOUND = BUILDMPI(3);
	BOUND->V[0] = 58368;
	BOUND->V[1] = 21515;
	BOUND->V[2] = 2;
	BOUND->S = 1;
	if (RSV(N, BOUND) >= 0) {
		printf("n>="); PRINTI(BOUND); printf("\n");
		FREEMPI(BOUND);
		return (NULL);
	}
	FREEMPI(BOUND);
	ONE = ONEI();
	s = RSV(M, ONE);
	t = RSV(N, ONE);
	FREEMPI(ONE);
	if (s <= 0) {
		printf("m <= 1\n");
		return (NULL);
	}
	if (t <= 0) {
		printf("n <= 1\n");
		return (NULL);
	}
	NUMBER = PRIME_GENERATOR(M, N);
	return (NUMBER);
}

MPI *SIGMAK(USI k, MPI *N)
/*
*  Returns the sum of the kth powers of the divisors of N,
*  returns NULL if factorization unsuccessful .
*  Here 0 < k <= 10000 and N >= 1.
*/
{
	MPI *U, *V, *W, *Tmp;
	unsigned int i, s;

	if (EQONEI(N))
		return (ONEI());
	/* FREE_PRIMES = 0, so we can free the Q_[i] later */
	s = FACTOR(N);
	if (s == 0)
	{
		return ((MPI *)NULL);
	}
	U = ONEI();
	for (i = 0; i < scount; i++)
	{
		V = POWER_I((long)(Q[i]), k*(K[i] + 1));
		Tmp = V;
		V = SUB0_I(V, 1);
		FREEMPI(Tmp);
		W = POWER_I((long)(Q[i]), k);
		Tmp = W;
		W = SUB0_I(W, 1);
		FREEMPI(Tmp);
		Tmp = V;
		V = INT0(V, W);
		FREEMPI(Tmp);
		FREEMPI(W);
		Tmp = U;
		U = MULTI(U, V);
		FREEMPI(Tmp);
		FREEMPI(V);
	}
	for (i = 0; i < s_count; i++)
	{
		V = POWERI(Q_[i], k*(K_[i] + 1));
		Tmp = V;
		V = SUB0_I(V, 1);
		FREEMPI(Tmp);
		W = POWERI(Q_[i], k);
		Tmp = W;
		W = SUB0_I(W, 1);
		FREEMPI(Tmp);
		Tmp = V;
		V = INT0(V, W);
		FREEMPI(W);
		FREEMPI(Tmp);
		Tmp = U;
		U = MULTI(U, V);
		FREEMPI(Tmp);
		FREEMPI(V);
	}
	for (i = 0; i < s_count; i++)
		FREEMPI(Q_[i]);
	return (U);
}

/*
*  Returns the sum of the kth powers of the divisors of N,
*  returns NULL if factorization unsuccessful .
*  Here 0 < k <= 10000 and N >= 1.
*/
MPI *SIGMAKX(MPI *K, MPI *N) {
	USI k;

	if (K->S <= 0) {
		printf("k <= 0\n");
		return(NULL);
	}
	if (N->S <= 0) {
		printf("n <= 0\n");
		return(NULL);
	}
	if (K->D > 0) {
		printf("k >= 2^16\n");
		return(NULL);
	}
	k = (USI)CONVERTI(K);
	if (k > 10000) {
		printf("k > 10000\n");
		return(NULL);
	}
	return(SIGMAK(k, N));
}


/* TAU(n) is Ramanujan's tau function. See T.M. Apostol, 'Modular functions
* and Dirichlet series in number theory', 20-22, 140.
* We assume n < 2^16
*/
MPI *TAU(USI n) {
	USI i;
	MPI *S, *T, *TT, *TEMP1, *TEMP2, *TEMP3, *TEMP, *TEMPS;

	if (n == 1) {
		return(ONEI());
	}
	S = ZEROI();
	for (i = 1; i<n; i++) {
		TT = CHANGE((USL)i);
		TEMP1 = SIGMAK(5, TT);
		FREEMPI(TT);
		TT = CHANGE((USL)(n - i));
		TEMP2 = SIGMAK(5, TT);
		FREEMPI(TT);
		TEMP = MULTI(TEMP1, TEMP2);
		FREEMPI(TEMP1);
		FREEMPI(TEMP2);
		TEMPS = S;
		S = ADD0I(S, TEMP);
		FREEMPI(TEMPS);
		FREEMPI(TEMP);
	}
	TT = CHANGE(n);
	TEMP1 = SIGMAK(11, TT);
	TEMP = TEMP1;
	TEMP1 = MULT_I(TEMP1, 65);
	FREEMPI(TEMP);

	TEMP2 = SIGMAK(5, TT);
	FREEMPI(TT);
	TEMP = TEMP2;
	TEMP2 = MULT_I(TEMP2, 691);
	FREEMPI(TEMP);

	TEMP3 = MULT_I(S, 252);
	TEMP = TEMP3;
	TEMP3 = MULT_I(TEMP3, 691);
	FREEMPI(TEMP);

	T = ADDI(TEMP1, TEMP2);
	FREEMPI(TEMP1);
	FREEMPI(TEMP2);
	TEMP = T;
	T = SUBI(T, TEMP3);
	FREEMPI(TEMP);
	FREEMPI(TEMP3);

	TEMP = T;
	T = INT_(T, 756);
	FREEMPI(TEMP);

	FREEMPI(S);
	return(T);
}

/* TAUX(N) is Ramanujan's tau function. See T.M. Apostol, 'Modular functions
* and Dirichlet series in number theory', 20-22, 140.
* We assume n < 2^16
*/
MPI *TAUX(MPI *N) {
	USI n;

	if (N->S <= 0) {
		printf("n <= 0\n");
		return(NULL);
	}
	if (N->D > 0) {
		printf("n >= 2^16\n");
		return(NULL);
	}
	n = (USI)CONVERTI(N);
	return(TAU(n));
}

MPI *TAU_PRIMEPOWER(USI n, USI p) {
	USI j, t, temp1, temp2;
	MPI *S, *T, *U, *TEMP, *TEMP1, *TEMP2;

	t = n / 2;
	S = ZEROI();
	for (j = 0; j <= t; j++) {
		temp1 = n - j;
		temp2 = temp1 - j;
		TEMP = BINOMIAL(temp1, temp2);
		TEMP1 = POWER_I(p, 11 * j);
		T = TAU(p);
		TEMP2 = POWERI(T, temp2);
		FREEMPI(T);
		U = MULTI3(TEMP, TEMP1, TEMP2);
		FREEMPI(TEMP);
		FREEMPI(TEMP1);
		FREEMPI(TEMP2);
		if (j % 2) {
			U->S = -(U->S);
		}
		TEMP = S;
		S = ADDI(S, U);
		FREEMPI(TEMP);
		FREEMPI(U);
	}
	return(S);
}

MPI *TAU_PRIMEPOWERX(MPI *N, MPI *P) {
	USI n, p, t;
	MPI *T;
	if (N->S <= 0) {
		printf("n <= 0\n");
		return(NULL);
	}
	if (N->D > 0) {
		printf("n >= 2^16\n");
		return(NULL);
	}
	n = (USI)CONVERTI(N);
	T = LUCAS(P);
	t = EQONEI(T);
	FREEMPI(T);
	if (t == 0) {
		printf("p is not prime\n");
		return(NULL);
	}

	p = (USI)CONVERTI(P);
	return(TAU_PRIMEPOWER(n, p));
}

MPI *TAU_COMPOSITE(USI n) {
	USI i, d, *KGLOB, *QGLOB;
	MPI *U, *N, *TEMP, *T;

	U = ONEI();
	N = CHANGE(n);
	d = FACTOR(N);
	FREEMPI(N);
	KGLOB = (USI *)mmalloc(d * sizeof(USI));
	QGLOB = (USI *)mmalloc(d * sizeof(USI));
	for (i = 0; i < scount; i++) {
		KGLOB[i] = K[i];
		QGLOB[i] = Q[i];
	}
	/*for(i = scount; i < d; i++){*/
	for (i = 0; i < s_count; i++) {
		KGLOB[i + scount] = K_[i];
		QGLOB[i + scount] = CONVERTI(Q_[i]); /* n < 2^32 here */
	}
	for (i = 0; i < d; i++) {
		TEMP = U;
		T = TAU_PRIMEPOWER(KGLOB[i], QGLOB[i]);
		U = MULTI(U, T);
		FREEMPI(TEMP);
		FREEMPI(T);
	}
	for (i = 0; i < s_count; i++)
		FREEMPI(Q_[i]);
	ffree(KGLOB, d * sizeof(USI));
	ffree(QGLOB, d * sizeof(USI));
	return(U);
}

MPI *TAU_COMPOSITEX(MPI *N) {
	USI n;

	if (N->S <= 0) {
		printf("n <= 0\n");
		return(NULL);
	}
	if (N->D > 0) {
		printf("n >= 2^16\n");
		return(NULL);
	}
	n = (USI)CONVERTI(N);
	return(TAU_COMPOSITE(n));
}

MPI *DIVISOR_PMINUS1(MPI *N, MPIA *SORTED_PRIME_FOUND) {
	/* This takes an even positive N and finds the primes such that p-1 divides N and returns an array of them. */
	USI prime_count, t, i, *kglobal, list_length, k, e;
	MPIA qglobal, list, PRIME_FOUND;
	MPI *TEMP, *L, *d, *PRIME_COUNT, *F;

	PRIME_FOUND = BUILDMPIA();
	t = FACTOR(N);
	qglobal = BUILDMPIA();
	kglobal = (USI *)mmalloc(t * sizeof(USI));
	for (i = 0; i < scount; i++) {
		TEMP = CHANGE(Q[i]);
		ADD_TO_MPIA(qglobal, TEMP, i);
		FREEMPI(TEMP);
		kglobal[i] = K[i];
	}
	for (i = 0; i < s_count; i++) {
		ADD_TO_MPIA(qglobal, Q_[i], i + scount);
		FREEMPI(Q_[i]);
		kglobal[i + scount] = K_[i];
	}
	list_length = kglobal[0] + 1;
	prime_count = 0;
	list = BUILDMPIA();
	for (k = 0; k < list_length; k++) {
		TEMP = POWERI((qglobal->A)[0], k);
		ADD_TO_MPIA(list, TEMP, k);
		FREEMPI(TEMP);
		L = ADD_ONEI((list->A)[k]);
		F = LUCAS(L);
		if (F->S) {
			ADD_TO_MPIA(PRIME_FOUND, L, prime_count);
			prime_count++;
		}
		FREEMPI(F);
		FREEMPI(L);
	}
	for (i = 1; i < t; i++) {
		e = list_length;
		for (k = 0; k < list_length; k++) {
			for (j = 1; j <= kglobal[i]; j++) {
				TEMP = POWERI((qglobal->A)[i], j);
				d = MULTI((list->A)[k], TEMP);
				FREEMPI(TEMP);
				L = ADD_ONEI(d);
				ADD_TO_MPIA(list, d, e);
				FREEMPI(d);
				F = LUCAS(L);
				if (F->S) {
					ADD_TO_MPIA(PRIME_FOUND, L, prime_count);
					prime_count++;
				}
				FREEMPI(F);
				FREEMPI(L);
				e++;
			}
		}
		list_length = list_length * (kglobal[i] + 1);
	}

	PRIME_COUNT = CHANGE(prime_count);
	SORT_ARRAY_MPI(PRIME_FOUND, SORTED_PRIME_FOUND, PRIME_COUNT);
	FREEMPIA(PRIME_FOUND);
	FREEMPIA(list);
	FREEMPIA(qglobal);
	ffree(kglobal, t * sizeof(USI));

	return(PRIME_COUNT);
}

void SORT_ARRAY_MPI(MPIA A, MPIA *B, MPI *N) {
	/* Algorithm from Stephen G. Kochan, Programming in C, 121-122. */
	/* n = number of array elements */
	USI t, s, i, j, n;
	MPI *TEMPI, *TEMPJ, *TEMP;

	*B = BUILDMPIA();
	n = CONVERTI(N);
	for (i = 0; i < n; i++) {
		TEMP = COPYI((A->A)[i]);
		ADD_TO_MPIA(*B, TEMP, i);
		FREEMPI(TEMP);
	}
	t = n - 1;
	for (i = 0; i < t; i++) {
		s = i + 1;
		for (j = s; j < n; j++) {
			if (COMPAREI(((*B)->A)[i], ((*B)->A)[j]) > 0) {
				TEMPI = COPYI(((*B)->A)[i]);
				TEMPJ = COPYI(((*B)->A)[j]);
				ADD_TO_MPIA(*B, TEMPJ, i);
				ADD_TO_MPIA(*B, TEMPI, j);
				FREEMPI(TEMPI);
				FREEMPI(TEMPJ);
				printf("\n");
			}
		}
	}
	return;
}

void CARMICHAEL(MPI *N) {
	/* Modified version of the Anthony J.Towns C program for solving phi(x)=n  26th August 2010.
	* A.J.'s program was written in April 1997, while he was in my MP206 class.
	* Also see section 3 of the paper Complexity of inverting the Euler function by
	* Scott Contini, Ernie Croot and Igor Shparlinkski, Math. Comp. 75 (2006), 983-996
	*/
	MPIA PRIME_FOUND, PHI_N, QPRIME;
	MPI *PRIME_COUNT, *RESULT, *TEMP, *TEMP0;
	USL prime_count;
	USI number, index, *a, t, u;
	int s;

	if (EQONEI(N)) {
		printf("1 and 2 are the solutions\n");
		return;
	}
	if ((N->V[0]) % 2) {
		printf("n is odd and > 1, so the number of solutions is 0\n");
		return;
	}
	PRIME_COUNT = DIVISOR_PMINUS1(N, &PRIME_FOUND);
	prime_count = CONVERTI(PRIME_COUNT);
	if (prime_count > 100) {
		printf("There are %lu primes p with p-1 dividing n. So we stop if we find two solutions.\n", prime_count);
	}
	FREEMPI(PRIME_COUNT);
	TEMP0 = ZEROI();
	ADD_TO_MPIA(PRIME_FOUND, TEMP0, prime_count);
	FREEMPI(TEMP0);

	number = 0;
	index = 0;
	QPRIME = BUILDMPIA();
	TEMP = TWOI();
	ADD_TO_MPIA(QPRIME, TEMP, 0);
	FREEMPI(TEMP);
	PHI_N = BUILDMPIA();
	ADD_TO_MPIA(PHI_N, N, 0);
	a = (USI *)mmalloc(100 * sizeof(USI));
	a[0] = 0;
	while (1) {
		t = EQONEI((PHI_N->A)[index]);
		if (t) {
			number++;
			CONSTRUCT_N(QPRIME, a, index, &RESULT);
			PRINTI(RESULT);
			FREEMPI(RESULT);
			printf(" ");
			if (number == 2 && prime_count > 100) {
				printf("are two solutions\n");
				return;
			}
			index--;
		}
		else {/* start of first else */
			TEMP = ADD_ONEI((PHI_N->A)[index]);
			t = RSV((QPRIME->A)[index], TEMP);
			FREEMPI(TEMP);
			if (((QPRIME->A)[index])->S == 0 || (index < prime_count && t == 1)) {
				if (index == 0) {
					printf("\nThe number of solutions x is %u\n", number);
					ffree(a, 100 * sizeof(USI));
					FREEMPIA(PRIME_FOUND);
					FREEMPIA(QPRIME);
					FREEMPIA(PHI_N);
					return;
				}
				TEMP = MOD0((PHI_N->A)[index], (QPRIME->A)[index - 1]);
				s = TEMP->S;
				FREEMPI(TEMP);
				if (s == 0) {
					a[index - 1] = a[index - 1] + 1;
					TEMP = INT((PHI_N->A)[index], (QPRIME->A)[index - 1]);
					ADD_TO_MPIA(PHI_N, TEMP, index);
					FREEMPI(TEMP);
					TEMP = COPYI((QPRIME->A)[index - 1]);
					ADD_TO_MPIA(QPRIME, TEMP, index);
					FREEMPI(TEMP);
				}
				else {
					index--;
				}
			}
			else {
				TEMP0 = SUB_ONEI((QPRIME->A)[index]);
				TEMP = MOD0((PHI_N->A)[index], TEMP0);
				s = TEMP->S;
				FREEMPI(TEMP);
				if (s == 0) {
					TEMP = INT((PHI_N->A)[index], TEMP0);
					FREEMPI(TEMP0);
					ADD_TO_MPIA(PHI_N, TEMP, index + 1);
					FREEMPI(TEMP);
					TEMP = COPYI((QPRIME->A)[index]);
					ADD_TO_MPIA(QPRIME, TEMP, index + 1);
					FREEMPI(TEMP);
					a[index] = 1;
					index++;
				}
				else {
					FREEMPI(TEMP0);
				}
			}
		}/* end of first else loop */
		for (u = 0; u < prime_count; u++) {
			if (EQUALI((QPRIME->A)[index], (PRIME_FOUND->A)[u])) {
				TEMP = COPYI((PRIME_FOUND->A)[u + 1]);
				ADD_TO_MPIA(QPRIME, TEMP, index);
				FREEMPI(TEMP);
				break;
			}
		}
	}/* end of big while loop */
}

void CONSTRUCT_N(MPIA QPRIME, USI *a, USI index, MPI **RESULT) {
	USI i, j;
	MPI *TEMP;

	*RESULT = ONEI();
	for (i = 0; i < index; i++) {
		for (j = 0; j < a[i]; j++) {
			TEMP = *RESULT;
			*RESULT = MULTI(*RESULT, (QPRIME->A)[i]);
			FREEMPI(TEMP);
		}
	}
	return;
}

MPI *LPRIME_FACTOR(MPI *N)
/* Returns the least prime factor of N.*/
{
	unsigned int i, s, t;
	MPIA AA;
	MPI *MINP, *QQ;

	s = FACTOR(N);
	if (s == 0) {
		return ((MPI *)NULL);
	}

	AA = BUILDMPIA();
	t = scount + s_count;
	AA->size = t;
	for (i = 0; i < t; i++) {
		if (i < scount) {
			QQ = CHANGE(Q[i]);
			ADD_TO_MPIA(AA, QQ, i);
		}
		else {
			QQ = COPYI(Q_[i - scount]);
			ADD_TO_MPIA(AA, QQ, i);
		}
		FREEMPI(QQ);
	}
	MINP = MIN_MPI_ARRAY(AA);
	for (i = 0; i < s_count; i++) {
		FREEMPI(Q_[i]);
	}
	FREEMPIA(AA);
	return(MINP);
}

MPI *CONWAY(MPI *A, MPI *B) {
	MPI *E, *Z, *F;

	E = ADD0I(A, B);
	Z = LPRIME_FACTOR(E);
	if (EQUALI(E, Z)) {
		FREEMPI(E);
		return(Z);
	}
	else {
		F = INT0(E, Z);
		FREEMPI(E);
		FREEMPI(Z);
		return(F);
	}
}

MPI *CONWAY_CYCLES(MPI *A, MPI *B) {
	USL i, j, twoi, iplus1, twoiplus1, diff;
	MPI *TEMPB, *TEMPA, *S, *T, *TEMPS, *TEMPT, *DIFF, *AA, *BB;

	i = 0;
	S = COPYI(A);
	T = COPYI(B);
	AA = COPYI(A);
	BB = COPYI(B);
	while (1) {
		i = i + 1;
		TEMPB = BB;
		TEMPA = AA;
		BB = CONWAY(AA, BB);
		FREEMPI(TEMPA);
		AA = TEMPB;

		TEMPS = S;
		TEMPT = T;
		S = CONWAY(S, T);
		FREEMPI(TEMPS);
		T = CONWAY(T, S);
		FREEMPI(TEMPT);
		if (EQUALI(AA, S) && EQUALI(BB, T)) {
			twoi = 2 * i;
			iplus1 = i + 1;
			twoiplus1 = twoi + 1;
			printf("x[%lu] = x[%lu] = ", i, twoi);
			PRINTI(AA);
			printf(", ");
			printf("x[%lu] = x[%lu] = ", iplus1, twoiplus1);
			PRINTI(BB);
			printf("\n");
			printf("cycle: ");
			PRINTI(AA);
			printf(" ");
			break;
		}
	}
	j = i;
	while (1) {
		i = i + 1;
		TEMPB = BB;
		if (!EQUALI(AA, S) || !EQUALI(BB, T)) {
			PRINTI(AA);
			printf(" ");
		}
		TEMPA = AA;
		BB = CONWAY(AA, BB);
		FREEMPI(TEMPA);
		AA = TEMPB;
		if (EQUALI(AA, S) && EQUALI(BB, T)) {
			FREEMPI(AA);
			FREEMPI(BB);
			FREEMPI(S);
			FREEMPI(T);
			diff = i - j;
			DIFF = CHANGE(diff);
			printf("\n");
			return(DIFF);
		}
	}
	if (i == conway_limit) {
		return(NULL);
	}
}

/* CONWAY_CYCLES without printing */
MPI *CONWAY_CYCLES0(MPI *A, MPI *B) {
	USL i, j, twoi, iplus1, twoiplus1, diff;
	MPI *TEMPB, *TEMPA, *S, *T, *TEMPS, *TEMPT, *AA, *BB, *DIFF;

	i = 0;
	S = COPYI(A);
	T = COPYI(B);
	AA = COPYI(A);
	BB = COPYI(B);
	while (1) {
		i = i + 1;
		TEMPB = BB;
		TEMPA = AA;
		BB = CONWAY(AA, BB);
		FREEMPI(TEMPA);
		AA = TEMPB;

		TEMPS = S;
		TEMPT = T;
		S = CONWAY(S, T);
		FREEMPI(TEMPS);
		T = CONWAY(T, S);
		FREEMPI(TEMPT);
		if (EQUALI(AA, S) && EQUALI(BB, T)) {
			twoi = 2 * i;
			iplus1 = i + 1;
			twoiplus1 = twoi + 1;
			break;
		}
	}
	j = i;
	while (1) {
		i = i + 1;
		TEMPB = BB;
		TEMPA = AA;
		BB = CONWAY(AA, BB);
		FREEMPI(TEMPA);
		AA = TEMPB;
		if (EQUALI(AA, S) && EQUALI(BB, T)) {
			FREEMPI(AA);
			FREEMPI(BB);
			FREEMPI(S);
			FREEMPI(T);
			diff = i - j;
			DIFF = CHANGE(diff);
			return(DIFF);
		}
	}
	if (i == conway_limit) {
		return(NULL);
	}
}

/* This tests Conway's conjecture for a range of initial values (a,b) with a<=m, b,=m.*/
void CONWAY_RANGE_TEST(MPI *M, MPI *N) {
	USL m, n, i, j, e;
	MPI *I, *J, *E;

	m = CONVERTI(M);
	n = CONVERTI(N);
	for (i = 1; i <= m; i++) {
		for (j = 1; j <= n; j++) {
			if (i % 10 == 0 && j % 10 == 0) {
				printf("(i,j)=(%lu,%lu)\n", i, j);
			}
			I = CHANGE(i);
			J = CHANGE(j);
			E = CONWAY_CYCLES0(I, J);
			FREEMPI(E);
			FREEMPI(I);
			FREEMPI(J);
			e = CONVERTI(E);
			if (e != 1 && e != 10 && e != 11 && e != 18 && e != 19 && e != 56 && e != 136) {
				printf("exception: (i,j)=(%lu,%lu):e=%lu\n", i, j, e);
				return;
			}
		}
	}
	printf("finished\n");
	return;
}


/* This tests Conway's conjecture for a range of initial values (a,b) with m1<=a<=m2, n1<=b<=n2.*/
void CONWAY_RANGE_TEST1(MPI *M1, MPI *M2, MPI *N1, MPI *N2) {
	USL m1, n1, m2, n2, i, j, e, temp;
	MPI *I, *J, *E;

	m1 = CONVERTI(M1);
	n1 = CONVERTI(N1);
	m2 = CONVERTI(M2);
	n2 = CONVERTI(N2);
	if (m2<m1) {
		temp = m2;
		m2 = m1;
		m1 = temp;
	}
	if (n2<n1) {
		temp = n2;
		n2 = n1;
		n1 = temp;
	}
	for (i = m1; i <= m2; i++) {
		for (j = n1; j <= n2; j++) {
			if (i % 10 == 0 && j % 10 == 0) {
				printf("(i,j)=(%lu,%lu)\n", i, j);
			}
			I = CHANGE(i);
			J = CHANGE(j);
			E = CONWAY_CYCLES0(I, J);
			FREEMPI(E);
			FREEMPI(I);
			FREEMPI(J);
			e = CONVERTI(E);
			if (e != 1 && e != 10 && e != 11 && e != 18 && e != 19 && e != 56 && e != 136) {
				printf("exception: (i,j)=(%lu,%lu):e=%lu\n", i, j, e);
				return;
			}
		}
	}
	printf("finished\n");
	return;
}

