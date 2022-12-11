/* mpqsieve.c */
/*
* Multiple Quadratic sieve method of factoring.
* Based on "The multiple polynomial quadratic sieve ", by R.D. Silverman,
* Mathematics of Computation, 48, 1987, 329-339.
*/

#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_NONSTDC_NO_DEPRECATE
#include <stdio.h>
#include <stdlib.h>
#include "integer.h"
#include "fun.h"

MPI *MPQS1(MPI *N)
/*
* The packed array version of MPQS(), done with a lot of help from Peter Adams.
*/
{
	register unsigned int z;
	unsigned int gap, h, i, j, t;
	unsigned long a, a1, b, b1, l1, l2, r, e, k;
	unsigned int TEST_ROWS;
	unsigned int factored = 0, flag, found = 0;
	unsigned int FBASE, TOLERANCE;
	unsigned int CNUF;
	unsigned long SRANGE, SSRANGE;
	unsigned int *logn, *logQ;
	unsigned long *r1, *r2, *soln, **M1, **M2, n;
	unsigned int *exponents, *CTR;
	MPI *N0, *X1, *X2, *X3, *X4, *D, *H0, *H1, *H2, *HH, *SAVE1, *Y;
	MPI *B, *A, *Q, *Tmp, **H, **cofactor, *ROOT2N, *OVER2, *ROOTN;
	MPI *P1, *P2, *T, *G, *srange, *xx, *hh;
	int *INDEX, s;
	long x, y, L, L1, L2, *base;
	unsigned int m1cols, m2cols, m1bitcols, m2bitcols, index, ctr;

	G = NULL;
	n = LENGTHI(N);
	printf("N = "); PRINTI(N); printf(" has %lu digits\n", n);
	/* FBASE = no of odd primes in factor base.  */
	if (n <= 20)
	{
		FBASE = 150; SRANGE = 5000; TOLERANCE = 10;
	}
	else if (n > 20 && n <= 25)
	{
		FBASE = 150; SRANGE = 10000; TOLERANCE = 15;
	}
	else if (n > 25 && n <= 30)
	{
		FBASE = 250; SRANGE = 15000; TOLERANCE = 20;
	}
	else if (n > 30 && n <= 35)
	{
		FBASE = 400; SRANGE = 15000; TOLERANCE = 20;
	}
	else if (n > 35 && n <= 40)
	{
		FBASE = 500; SRANGE = 25000; TOLERANCE = 20;
	}
	else if (n > 40 && n <= 45)
	{
		FBASE = 900; SRANGE = 50000; TOLERANCE = 20;
	}
	else if (n > 45 && n <= 50)
	{
		FBASE = 1200; SRANGE = 100000; TOLERANCE = 22;
	}
	else /* n > 50 */
	{
		if (n > 61)
		{
			printf("Number of digits > 60.\n");
			return(NULL);
		}
		FBASE = 3000; SRANGE = 180000; TOLERANCE = 22;
	}
	printf("FBASE = %u:", FBASE);
	printf(" sieve_range = %lu:", SRANGE);
	printf(" tolerance  = %u:\n", TOLERANCE);
	ROOTN = BIG_MTHROOT(N, 2);
	TEST_ROWS = FBASE + 2;
	INDEX = (int *)mmalloc((USL)((TEST_ROWS) * sizeof(int)));
	INTSETI(INDEX, (USL)TEST_ROWS, -1);
	CTR = (unsigned int *)ccalloc(TEST_ROWS, sizeof(unsigned int));
	m1cols = (FBASE + 2);
	m1bitcols = m1cols / O32 + 1;
	M1 = idim2(TEST_ROWS, m1bitcols);
	H = (MPI **)mmalloc((USL)((TEST_ROWS) * sizeof(MPI *)));
	cofactor = (MPI **)mmalloc((USL)((TEST_ROWS) * sizeof(MPI *)));
	base = (long *)mmalloc((USL)((FBASE + 2) * sizeof(long)));
	base[0] = 2; base[FBASE + 1] = -1;
	soln = (unsigned long *)mmalloc((USL)((FBASE + 1) * sizeof(unsigned long)));
	/* soln[i] is a solution of x*x=n (mod base[i], 1 <= i <= FBASE. */
	QERATOSTHENES(N, (USL *)base, soln, FBASE);
	CNUF = BINARYB(N) / 2 + binary(SRANGE) - (TOLERANCE * binary((USL)(base[FBASE])) / 10);
	r1 = (unsigned long *)mmalloc((USL)((FBASE + 1) * sizeof(unsigned long)));
	r2 = (unsigned long *)mmalloc((USL)((FBASE + 1) * sizeof(unsigned long)));
	/* we will sieve over -SRANGE <= x <= SRANGE. */
	X1 = MULT_I(N, 2L);
	ROOT2N = BIG_MTHROOT(X1, 2);
	FREEMPI(X1);
	SSRANGE = (USL)2 * SRANGE;
	logQ = (unsigned int *)mmalloc((USL)(((unsigned int)SSRANGE + 1) * sizeof(unsigned int)));
	logn = (unsigned int *)mmalloc((USL)((FBASE + 1) * sizeof(unsigned int)));
	/* Table lookup: logn[i] = binary((USL)(base[i])), 1 <= i <= FBASE. */
	k = N->V[0];
	if (k % 8 == 1)
		logn[0] = 3;
	if (k % 8 == 5)
		logn[0] = 2;
	if (k % 4 == 3)
		logn[0] = 1;
	for (i = 1; i <= FBASE; i++)
		logn[i] = binary((USL)(base[i]));
	h = 0;
	srange = CHANGE(SRANGE);
	X3 = INT0(ROOT2N, srange); /* bug site */
	FREEMPI(srange);
	OVER2 = BIG_MTHROOT(X3, 2);
	FREEMPI(X3);

	/* calculating the polynomial Q(x) = A*x*x + 2*B*x + C */
	Tmp = CHANGE((USL)(base[FBASE]));
	if (RSV(OVER2, Tmp) <= 0)
	{
		T = OVER2;
		OVER2 = ADD0_I(Tmp, (USL)1);
		FREEMPI(T);
	}
	FREEMPI(Tmp);
	while (factored < TEST_ROWS)
	{
		printf("changing polynomial\n");
		if (h)
		{
			Tmp = OVER2;
			hh = CHANGE(h);
			OVER2 = ADD0I(OVER2, hh); /* bugsite */
			FREEMPI(hh);
			FREEMPI(Tmp);
		}
		D = PINCH_PRIME(OVER2, N, base, FBASE, &gap);
		/*
		D = PROB_PRIME(OVER2, N, &gap);
		*/
		h += gap + 2;
		X1 = INT0_(D, (USL)4);
		N0 = MOD0(N, D);
		H0 = MPOWER(N0, X1, D);
		FREEMPI(X1);
		H1 = MULTM(N0, H0, D);
		FREEMPI(N0);
		A = MULTI(D, D);
		X1 = ADD0_I(D, (USL)1);
		X2 = INT0_(X1, (USL)2);
		FREEMPI(X1);
		Y = MULTM(X2, H0, D);
		FREEMPI(H0);
		FREEMPI(X2);
		X1 = MULTI(H1, H1);
		X2 = SUBI(N, X1);
		FREEMPI(X1);
		X3 = INT(X2, D);
		FREEMPI(X2);
		Tmp = MOD(X3, D);
		FREEMPI(X3);
		H2 = MULTM(Y, Tmp, D);
		FREEMPI(Y);
		FREEMPI(Tmp);
		X4 = MULTI(H2, D);
		FREEMPI(H2);
		B = ADD0I(H1, X4);
		FREEMPI(H1);
		FREEMPI(X4);
		SAVE1 = INVERSEM(D, N); /* 1/(D) (mod N) */
								/*
								printf("A = ");
								PRINTI(A);
								X1 = MULTI(B, B);
								printf(", B = ");
								PRINTI(B);
								X2 = SUBI(X1, N);
								C = INT(X2, A);
								X3 = MOD(X2, A);
								FREEMPI(X1);
								FREEMPI(X2);
								FREEMPI(X3);
								printf(", C = ");
								PRINTI(C);
								printf("\n");
								*/
								/* C is freed later */
								/*
								printf("calculating the two roots r1[i], r2[i] (mod base[i]) of Q(x)= 0\n");
								*/
		r = (B->V[0]) % (USL)2 ? (USL)0 : (USL)1;
		for (i = 1; i <= FBASE; i++)
		{
			/*
			printf("mod base[%u], A, B, C = :%lu, %lu, %lu\n",i, MOD_(A,(USL)(base[i])), MOD_(B,(USL)(base[i])), MOD_(C,(USL)(base[i])));
			*/
			b = MOD0_(B, (USL)(base[i]));
			b1 = MINUSm(b, (USL)(base[i]));
			a = MOD0_(A, (USL)(base[i]));
			a1 = INVERSEm(a, (USL)(base[i]));
			r1[i] = MULTm(ADDm(b1, soln[i], (USL)(base[i])), a1, (USL)(base[i]));
			r2[i] = MULTm(ADDm(b1, MINUSm(soln[i], (USL)(base[i])), (USL)(base[i])), a1, (USL)(base[i]));
			/*
			printf("a=%lu,a1 = %lu, b1 = %lu\n",a, a1, b1);
			printf("soln[%u] = %lu\n", i, soln[i]);
			printf("r1[%u] = %lu, r2[%u] = %lu ",i, r1[i], i,r2[i]);
			printf("base[%u] = %lu\n", i, base[i]);
			*/
		}
		/*
		printf("finished calculating the two roots r1, r2 (mod base[i]) of Q(x)= 0\n");
		*/

		/* sieving over the range -SRANGE <= x <= SRANGE */
		INTSETUI(logQ, (USL)SSRANGE + 1, 0);
		L = SRANGE % (USL)2 == r ? -SRANGE : -SRANGE + 1;
		for (y = L + SRANGE; y <= SSRANGE; y += 2)
			logQ[y] += logn[0];
		for (i = 1; i <= FBASE; i++)
		{
			l1 = (SRANGE + r1[i]) / (USL)(base[i]);
			l2 = (SRANGE + r2[i]) / (USL)(base[i]);
			L1 = r1[i] - (long)(l1 * (USL)(base[i]));
			L2 = r2[i] - (long)(l2 * (USL)(base[i]));
			/*
			printf("base=%ld,r1=%lu,r2=%lu,l1=%lu,l2=%lu,L1=%ld,L2=%ld\n",base[i],r1[i],r2[i],l1,l2,L1,L2);
			*/
			for (y = L1 + SRANGE; y <= SSRANGE; y += base[i])
				logQ[y] += logn[i];
			for (y = L2 + SRANGE; y <= SSRANGE; y += base[i])
				logQ[y] += logn[i];
		}
		/* Scanning the values logQ[x] to see if any exceed CNUF */
		for (y = 0; y <= SSRANGE; y++)
		{
			if (logQ[y] >= CNUF)
			{
				x = y - SRANGE;
				xx = CHANGEL(x);/* bugsite */
				X1 = MULTI(A, xx);
				FREEMPI(xx);
				X2 = ADDI(X1, B);
				X3 = MOD(X2, N);
				H[factored] = MULTM(X3, SAVE1, N);
				HH = MULTM(H[factored], H[factored], N);
				if (RSV(X2, ROOTN) <= 0)
					Q = SUBI(HH, N); /* Q->S = -1 */
				else
					Q = COPYI(HH); /* Q->S = 1 */
				FREEMPI(HH);
				FREEMPI(X1);
				FREEMPI(X2);
				FREEMPI(X3);
				cofactor[factored] = FACTORQ1(Q, factored, (USL *)base, M1, FBASE, x, r1, r2, r, INDEX, CTR);
				/* this factors Q over the factor base,
				updates M1[factored][] and returns the cofactor. */
				t = EQONEI(cofactor[factored]);
				FREEMPI(cofactor[factored]);
				if (t)
				{
					printf("Q[%u] = ", factored); PRINTI(Q); printf("\n");
					cofactor[factored] = Q;
					factored++;
				}
				else
				{
					FREEMPI(H[factored]);
					FREEMPI(Q);
				}
				if (factored == TEST_ROWS)
					break;
			}
		}
		FREEMPI(SAVE1);
		FREEMPI(A);
		FREEMPI(B);
		FREEMPI(D);
	}
	/*
	printf("INDEX = ");
	for (i = 0; i < TEST_ROWS; i++)
	printf("%d,", INDEX[i]);
	printf("\n");
	printf("CTR = ");
	for (i = 0; i < TEST_ROWS; i++)
	printf("%u,", CTR[i]);
	printf("\n");
	printf("The packed array M1 is:\n");
	printf("array M1 is:\n");
	for (i = 0; i < TEST_ROWS; i++)
	{
	for (j = 0; j <m1bitcols; j++)
	printf("%lu,", M1[i][j]);
	printf("\n");
	}
	printf("\n\n");
	print_array( M1, (int)TEST_ROWS, (int)m1cols);
	printf("\n\n");
	*/
	/* Now find the reduced form of M1. */
	m2cols = TEST_ROWS;
	m2bitcols = (m2cols / O32) + 1;
	M2 = idim2(TEST_ROWS, m2bitcols);
	ctr = 0;
	index = 0;
	for (z = 0; z < TEST_ROWS; z++)
	{
		M2[z][index] = (USL)1 << (O31 - ctr);
		ctr++;
		if (ctr == O32)
		{
			index++;
			ctr = 0;
		}
	}
	printf("doing Gauss reduction\n");
	REDUCTION1(M1, M2, TEST_ROWS, m1cols, m1bitcols, m2bitcols, INDEX, CTR);
	printf("finished Gauss reduction\n");
	printf("\n");
	/*
	printf("The packed array M1 is:\n");
	for (i = 0; i < TEST_ROWS; i++)
	{
	for (j = 0; j <m1bitcols; j++)
	printf("%lu,", M1[i][j]);
	printf("\n");
	}
	printf("\n\n");
	printf("The packed array M2 is:\n");
	for (i = 0; i < TEST_ROWS; i++)
	{
	for (j = 0; j <m2bitcols; j++)
	printf("%lu,", M2[i][j]);
	printf("\n");
	}
	printf("\n\n");
	printf("The array M1 is:\n");
	print_array( M1, (int)TEST_ROWS, (int)m1cols);
	printf("\n\n");
	printf("The array M2 is:\n");
	print_array( M2, (int)TEST_ROWS, (int)m2cols);
	printf("\n\n");
	*/

	/* Calculate P1=product of H[i]  where M1[i][] = 0; also P2 = product
	of prime powers to half-exponent sums . */
	exponents = (unsigned int *)mmalloc((USL)((FBASE + 1) * sizeof(unsigned int)));
	for (s = TEST_ROWS - 1; s >= 0; s--)
	{
		flag = 1;
		for (j = 0; j < m1bitcols; j++)
			if (M1[s][j])
			{
				flag = 0;
				break;
			}
		if (flag == 0)
			continue;
		/* found a zero row, initialize exponent sums to zero */
		INTSETUI(exponents, (USL)FBASE + 1, 0);
		P1 = ONEI();
		for (e = 0; e < TEST_ROWS; e++)
		{
			/*
			f = (M2[s][e>>5] & ((USL)1<<(31-(e&31))))?1:0;
			*/
			if (get_element(M2, s, e))
			{
				T = MULTM(P1, H[e], N);
				FREEMPI(P1);
				P1 = T;
				EXP_UPDATE(cofactor[e], (USL *)base, FBASE, exponents);
			}
		}
		P2 = ONEI();
		for (j = 0; j <= FBASE; j++)
		{
			if (exponents[j])
			{
				Y = CHANGE((USL)(base[j]));
				T = MPOWER_M(Y, (unsigned long)(exponents[j] / 2), N);
				FREEMPI(Y);
				Tmp = P2;
				P2 = MULTM(T, P2, N);
				FREEMPI(T);
				FREEMPI(Tmp);
			}
		}
		T = SUBI(P1, P2);
		G = GCD(T, N);
		FREEMPI(T);
		printf("P1 = "); PRINTI(P1); printf("\n");
		printf("P2 = "); PRINTI(P2); printf("\n");
		printf("GCD(P2 - P1, N) = "); PRINTI(G); printf("\n");
		FREEMPI(P1);
		FREEMPI(P2);
		if (!EQONEI(G) && !EQUALI(G, N))
		{
			printf("FOUND A NON_TRIVIAL FACTOR:  "); PRINTI(G); printf("\n");
			found = 1;
			break;
		}
		FREEMPI(G);
	}
	printf("N = "); PRINTI(N); printf("\n");
	FREEMPI(ROOTN);
	FREEMPI(ROOT2N);
	FREEMPI(OVER2);
	ffree((char *)base, (FBASE + 2) * sizeof(long));
	ffree((char *)soln, (FBASE + 1) * sizeof(unsigned long));
	ffree((char *)r1, (FBASE + 1) * sizeof(unsigned long));
	ffree((char *)r2, (FBASE + 1) * sizeof(unsigned long));
	ffree((char *)logn, (FBASE + 1) * sizeof(unsigned int));
	ffree((char *)logQ, (SSRANGE + 1) * sizeof(unsigned int));
	for (i = 0; i < factored; i++)
		FREEMPI(H[i]);
	ffree((char *)H, (TEST_ROWS) * sizeof(MPI *));
	for (i = 0; i < TEST_ROWS; i++)
		FREEMPI(cofactor[i]);
	ffree((char *)cofactor, (TEST_ROWS) * sizeof(MPI *));
	ffree((char *)INDEX, (TEST_ROWS) * sizeof(int));
	ffree((char *)CTR, (TEST_ROWS) * sizeof(unsigned int));
	ffree((char *)exponents, (FBASE + 1) * sizeof(unsigned int));
	ifree2(M1, TEST_ROWS, m1bitcols);
	ifree2(M2, TEST_ROWS, m2bitcols);
	if (found)
		return (G);
	else
		return (NULL);
}
