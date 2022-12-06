/* func.c */
#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_NONSTDC_NO_DEPRECATE
#include <stdio.h>
#include "integer.h"
#include "fun.h"

unsigned int VERBOSE;
unsigned int FREE_PRIMES;
extern unsigned long PRIME[];
unsigned int ECMAX;

MPI *FACTORX(MPI *N)
{
	unsigned int z;
	int s;

	if (N->S <= 0 || (N->D == 0 && N->V[0] == 1))
	{
		printf("N <= 1\n");
		return NULL;
	}
	printf("enter the number of elliptic curves to be used:");
	s = scanf("%u", &ECMAX);
	getchar();
	VERBOSE = 1;
	FREE_PRIMES = 1;
	z = FACTOR(N);
	VERBOSE = 0;
	FREE_PRIMES = 0;
	if (z == 0)
	{
		printf("Factorization unsuccessful\n");
		return NULL;
	}
	return (CHANGE(z));
}

/* Returns NULL on failure.  Wrapper handles this. */
MPI *Nextprime(MPI *N)
{
	MPI *tmp;
	if (N->S < 0)
	{
		printf("N <= 0\n");
		return NULL;
	}
	tmp = NEXT_PRIMEX(N);
	return(tmp);
}

MPI *JACOB(MPI *M, MPI *N)
/*
* Evaluates the Jacobi symbol (M/N); N odd, N > 0.
*/
{
	int t;
	MPI *MM;

	if (MOD0_(N, (USL)2) == 0)
	{
		printf("N is even\n");
		return NULL;
	}
	if (N->S < 0)
	{
		printf("N is negative\n");
		return NULL;
	}
	MM = MOD(M, N);
	if (N->D == 0)
		t = JACOBI(CONVERTI(MM), CONVERTI(N));
	else
		t = JACOBIB(MM, N);
	FREEMPI(MM);
	if (t == 0)
		return(ZEROI());
	else if (t == 1)
		return(ONEI());
	else
		return(MINUS_ONEI());
}


MPI *GCD_ARRAYVX(MPIA M, MPIA *Y)
/*
* Returns d=gcd(M[0],...,M[N-1]) and an array Y[] of MPI's such that
* d = M[0]*Y[0]+...+M[N-1]*Y[N-1]. Here N > 1.
*/
{
	MPI *Z;
	MPIA V;

	Z = GCD_ARRAYV(M, &V);
	*Y = V;
	return (Z);
}
MPI *SQRTMX(MPI *x, MPI *p)
/*
* Calculates sqrt(x) (mod p) using "A simple and fast probabilistic algorithm
* for computing square roots modulo a prime number", I.E.E.E. Trans. Inform.
* Theory, IT-32, 1986, 846-847, R. Peralta.
* Here x is a quadratic residue mod p. x can be negative.
* Returns NULL on failure.  This is taken care of by wrapper.
*/
{
	unsigned long X, P, Z;
	MPI *M, *N = NULL, *T;
	int t;

	if (p->S <= 0 || (p->D == 0 && p->V[0] <= 2))
	{
		printf("p <= 2\n");
	}
	T = LUCAS(p);
	t = T->S;
	FREEMPI(T);
	if (!t)
	{
		printf("2nd argument is not a prime\n");
		return N;
	}
	if (JACOBIB(x, p) != 1)
	{
		printf("x is not a quadratic residue mod p\n");
		return N;
	}
	if (p->D == 0)
	{
		M = MOD(x, p);
		X = CONVERTI(M);
		P = CONVERTI(p);
		Z = SQRTm(X, P);
		FREEMPI(M);
		return(CHANGE(Z));
	}
	else
	{
		N = SQRTM(x, p);
		return (N);
	}
}

MPI *CONGRX(MPI *A, MPI *B, MPI *M, MPI **N)
/*
* Returns the least solution (mod N) of the congruence AX=B(mod M), where
* N = M / gcd(A, M).
* returns NULL if function fails. Handled by wrapper CONGRX_W.
*/
{
	MPI *Z = NULL;

	if (M->S <= 0 || (M->D == 0 && M->V[0] == 1))
	{
		printf("M <= 1\n");
		*N = ZEROI();
		return Z;
	}
	Z = CONGR(A, B, M, N);
	if (Z == NULL)
	{
		printf("the congruence has no solution\n");
		*N = ZEROI();
	}
	return (Z);
}

MPI *CHINESEX(MPI *A, MPI *B, MPI *M, MPI *N, MPI **Mptr)
/*
* Returns the solution mod *Mptr=lcm[M,N] of the congruences X = A (mod M)
* and X = B (mod N), if soluble.
* Returns NULL if it failed.  This case is handled by wrapper CHINESEX_W.
*/
{
	MPI *Z = NULL;

	if (M->S <= 0 || (M->D == 0 && M->V[0] == 1))
	{
		printf("M <= 1\n");
		*Mptr = ZEROI();
		return(Z);
	}
	if (N->S <= 0 || (N->D == 0 && N->V[0] == 1))
	{
		printf("N <= 1\n");
		*Mptr = ZEROI();
		return(Z);
	}
	Z = CHINESE(A, B, M, N, Mptr);
	if (Z == NULL)
	{
		printf("the system has no solution\n");
		*Mptr = ZEROI();
	}
	return (Z);
}

MPI *CHINESE_ARRAYX(MPIA A, MPIA M, MPI **Mptr)
/*
* Returns the solution mod *Mptr=lcm[M[0],...,M[n-1] of the congruences
* X = A[i] (mod M[i]),0 <= i < N, if soluble; otherwise returns NULL.
*/
{
	unsigned int i, n;
	MPI *Z;

	n = ((A->size >= M->size) ? M->size : A->size);
	for (i = 0; i < n; i++)
	{
		if (M->A[i]->S <= 0 || (M->A[i]->D == 0 && M->A[i]->V[0] == 1))
		{
			*Mptr = ZEROI();
			printf("M[i] <= 1\n");
			return NULL;
		}
	}
	Z = CHINESE_ARRAY(A->A, M->A, Mptr, n);
	if (Z == NULL)
	{
		*Mptr = ZEROI();
		printf("the system has no solution\n");
		return NULL;

	}
	return (Z);
}

void MTHROOTX(MPI *Aptr, MPI *Bptr, MPI *M, MPI *R)
/*
* *Aptr and *Bptr are positive MPI'S.
* The mthroot of *Aptr/(*Bptr) is computed to R decimal places, R >= 0 ;
* M, R are integers, 0 < M * R < R0 * R0.
*/
{
	unsigned int m, r;

	if (Aptr->S < 0 || Bptr->S <= 0)
	{
		printf("Numerator or denominator <= 0\n");
		return;
	}
	if (M->S <= 0 || M->D >= 1 || M->V[0] == 1)
	{
		printf("m <= 1 or m >= R0\n");
		return;
	}
	if (R->S < 0 || R->D >= 1)
	{
		printf("r < 0 or r >= R0\n");
		return;
	}
	m = CONVERTI(M);
	r = CONVERTI(R);
	MTHROOT(Aptr, Bptr, m, r);
	return;
}

MPI *BIG_MTHROOTX(MPI *Aptr, MPI *M)
/*
* *Eptr, the integer part of the Mth root of the positive MPI *Aptr,
* 1 < M < R0, is obtained by Newton's method, using the integer part function.
* Returns NULL on failure.  Handled by wrapper function.
*/
{
	unsigned int m;
	if (Aptr->S < 0)
	{
		printf("argument < 0\n");
		return NULL;
	}
	if (M->S <= 0 || M->D >= 1 || M->V[0] == 1)
	{
		printf("m <= 1 or m >= R\n");
		return NULL;
	}
	m = CONVERTI(M);
	return (BIG_MTHROOT(Aptr, m));
}

MPI  *FUND_UNITX(MPI *D, MPI **Xptr, MPI **Yptr)
/*
* This is a program for finding the fundamental unit of Q(sqrt(D)).
* The algorithm is based on K. Rosen, Elementary number theory
* and its applications, p382, B.A. Venkov, Elementary Number theory, p.62
* and D. Knuth, Art of computer programming, Vol.2, p359, with Pohst's trick
* of using half the period.
* w=(1+sqrt(D))/2 if D=1 (mod 4), w=sqrt(D) otherwise.
* The norm of the fundamental unit is returned.
*
* Returns NULL on failure.  Wrapper function handles this.
*/
{
	MPI *G, *X, *tmp;
	unsigned int t;

	if (D->S <= 0)
	{
		printf("D <= 0\n");
		*Xptr = ZEROI();
		*Yptr = ZEROI();
		return NULL;
	}
	X = BIG_MTHROOT(D, 2);
	G = MULTI(X, X);
	t = EQUALI(D, G);
	FREEMPI(G);
	FREEMPI(X);
	if (t)
	{
		printf("D is a perfect square\n");
		*Xptr = ZEROI();
		*Yptr = ZEROI();
		return NULL;
	}
	tmp = FUND_UNIT(D, Xptr, Yptr);
	return (tmp);
}

MPI  *PELX(MPI *D, MPI *E, MPI **Xptr, MPI **Yptr)
/*
* This is a program for finding the least solution of Pell's equation
* x*x - D*y*y = +-1.
* The algorithm is based on K. Rosen, Elementary number theory
* and its applications, p382, B.A. Venkov, Elementary Number theory, p.62
* and D. Knuth, Art of computer programming, Vol.2, p359, with Pohst's trick
* of using half the period.
* The norm of the least solution is returned.
* The partial quotients are printed iff E->S != 0.
*
* Returns NULL on failure.  Wrapper handles this.
*/
{
	MPI *G, *X, *tmp;
	unsigned int t;

	if (D->S <= 0)
	{
		printf("D <= 0\n");
		*Xptr = ZEROI();
		*Yptr = ZEROI();
		return NULL;
	}
	X = BIG_MTHROOT(D, 2);
	G = MULTI(X, X);
	t = EQUALI(D, G);
	FREEMPI(G);
	FREEMPI(X);
	if (t)
	{
		printf("D is a perfect square\n");
		*Xptr = ZEROI();
		*Yptr = ZEROI();
		return NULL;
	}
	tmp = PEL(D, E, Xptr, Yptr);
	return (tmp);
}

MPI *SURDX(MPI *D, MPI *T, MPI *U, MPI *V, MPIA *A_ARRAY, MPIA *U_ARRAY, MPIA *V_ARRAY, MPIA *P_ARRAY, MPIA *Q_ARRAY)
/*
* This function uses the continued fraction algorithm expansion in K. Rosen,
* Elementary Number theory and its applications,p.379-381 and Knuth's
* The art of computer programming, Vol.2, p. 359. It locates the first complete
* quotient that is reduced and then uses the function REDUCED(D,U,V,i) above
* to locate and return the period of the continued fraction expansion.
*
* Returns NULL on failure.  Wrapper handles this.
*/
{
	MPI *G, *X;
	unsigned int t;
	unsigned long f;

	if (D->S <= 0)
	{
		printf("D <= 0\n");
		return NULL;
	}
	X = BIG_MTHROOT(D, 2);
	G = MULTI(X, X);
	f = EQUALI(D, G);
	FREEMPI(G);
	FREEMPI(X);
	if (f)
	{
		printf("D is a perfect square\n");
		return NULL;
	}
	if (V->S == 0)
	{
		printf("V = 0\n");
		return NULL;
	}
	t = SURD(D, T, U, V, A_ARRAY, U_ARRAY, V_ARRAY, P_ARRAY, Q_ARRAY, 0);
	return (CHANGE((USL)t));
}

MPI *MPOWERX(MPI *Aptr, MPI *Bptr, MPI *Cptr)
/*
* *Aptr, *Bptr, *Cptr, *Eptr  are MPI'S, *Cptr > 0, *Bptr >= 0,
* (*Aptr)^ *Bptr (mod *Cptr) is returned.
*
* Returns NULL on failure.  Handled by wrapper.
*/
{
	if (Cptr->S <= 0)
	{
		printf("Modulus <= 0\n");
		return NULL;
	}
	if (Cptr->D == 0 && Cptr->V[0] == 1)
	{
		printf("Modulus = 1\n");
		return NULL;
	}
	return (MPOWER(Aptr, Bptr, Cptr));
}

MPI *INVERSEMX(MPI *A, MPI *M)
/*
* Returns the inverse of A mod M.
*/
{
	MPI *temp, *Z;
	unsigned int t;

	if (M->S <= 0 || (M->D == 0 && M->V[0] == 1))
	{
		printf("M <= 1\n");
		return NULL;
	}
	temp = GCD(A, M);
	t = EQONEI(temp);
	FREEMPI(temp);
	if (!t)
	{
		printf("A is not relatively prime to M\n");
		return NULL;
	}
	temp = MOD(A, M);
	Z = INVERSEM(temp, M);
	FREEMPI(temp);
	return (Z);
}

void MILLERX(MPI *N, MPI *B)
/*
* *N is odd, > 1, and does not divide b, 0 < b < R0.
*  If *N passes Miller's test for base *B, *N is likely to be prime.
* If *N fails Miller's test for base *B, then *N is composite.
*/
{
	unsigned long b, r;
	int t;
	MPI *M;

	if (N->S <= 0 || (N->D == 0 && N->V[0] == 1))
	{
		printf("N <= 1\n");
		return;
	}
	if (B->S <= 0 || B->D >= 1 || B->V[0] == 1)
	{
		printf("BASE <= 1 or BASE >= R0\n");
		return;
	}
	if (N->V[0] % 2 == 0)
	{
		printf("N is even\n");
		return;
	}
	M = MOD0(B, N);
	t = M->S;
	FREEMPI(M);
	if (t == 0)
	{
		printf("N divides base\n");
		return;
	}
	b = CONVERTI(B);
	r = MILLER(N, b);
	PRINTI(N);
	if (r)
		printf(" passes Miller's test to base %lu\n", b);
	else
		printf(" fails Miller's test to base %lu\n", b);
	return;
}

MPI *DIVISORX(MPI *N)
/* Returns the number of divisors of N.
* Returns NULL on failure.  Handled by wrapper */
{
	MPI *M;

	if (N->S <= 0)
	{
		printf("N <= 0\n");
		return NULL;
	}
	M = DIVISOR(N);
	if (M == NULL)
	{
		printf("Factorization unsuccessful\n");
		return NULL;
	}
	return (M);
}

MPI *MOBIUSX(MPI *N)
/* Returns  Mu(N). */
{
	MPI *M;

	if (N->S <= 0)
	{
		printf("N <= 0\n");
		return NULL;
	}
	M = MOBIUS(N);
	if (M == NULL)
	{
		printf("Factorization unsuccessful\n");
		return NULL;
	}
	return (M);
}

MPI *EULERX(MPI *N)
/* Returns Euler's function phi(N). */
{
	MPI *M;

	if (N->S <= 0)
	{
		printf("N <= 0\n");
		return NULL;
	}
	M = EULER(N);
	if (M == NULL)
	{
		printf("Factorization unsuccessful\n");
		return NULL;
	}
	return (M);
}

MPI *SIGMAX(MPI *N)
/* Returns sigma(N), the sum of the divisors of N. */
{
	MPI *M;

	if (N->S <= 0)
	{
		printf("N <= 0\n");
		return NULL;
	}
	M = SIGMA(N);
	if (M == NULL)
	{
		printf("Factorization unsuccessful\n");
		return NULL;
	}
	return (M);
}

MPI *LPRIMROOTX(MPI *P)
/* Returns the least primitive root mod P. */
{
	MPI *M;

	if (P->S <= 0 || (P->D == 0 && P->V[0] == 1))
	{
		printf("P <= 1\n");
		return NULL;
	}
	M = LPRIMROOT(P);
	if (M->S == 0)
	{
		printf("factorization of P-1 unsuccessful\n");
		return NULL;
	}
	return (M);
}

MPI *ORDERPX(MPI *A, MPI *P)
/*
* Returns the order of A mod P, a prime.
*/
{
	if (P->S <= 0 || (P->D == 0 && P->V[0] == 1))
	{
		printf("P <= 1\n");
		return NULL;
	}
	return (ORDERP(A, P));
}

MPI *ORDERQX(MPI *A, MPI *P, MPI *N)
/*
* Returns the order of A mod P^N, P a prime.
*/
{
	unsigned int n;

	if (N->S <= 0 || N->D || (N->D == 0 && (N->V[0] > 1000)))
	{
		printf("N <= 0 or N > 1000\n");
		return NULL;
	}
	if (P->S <= 0 || (P->D == 0 && P->V[0] == 1))
	{
		printf("P <= 1\n");
		return NULL;
	}
	n = (USI)(CONVERTI(N));
	return (ORDERQ(A, P, n));
}

MPI *ORDERMX(MPI *A, MPI *M)
/*
* Returns the order of A mod M.
*/
{
	MPI *G;
	unsigned int t;

	if (M->S <= 0 || (M->D == 0 && M->V[0] == 1))
	{
		printf("M <= 1\n");
		return NULL;
	}
	G = GCD(A, M);
	t = EQONEI(G);
	FREEMPI(G);
	if (!t)
	{
		printf("GCD(A, M) > 1\n");
		return NULL;
	}
	return (ORDERM(A, M));
}

MPI *LUCASUX(MPI *N, MPI *Q, MPI *M)
/*
* returns U_n (mod m), where U_{k+1}=U_k-qU_{k-1}, U_0=0,U_1=1.
* D=1-4q != 0.
* We use the recurrence relations:
* U_{k+1}=U_k-qU_{k-1},
* U_{2k}=-2qU_{k-1}U_k+U_kU_k,
* U_{2k-1}=U_kU_k-qU_{k-1}U_{k-1}
* See Niven, Zuckermann and Montgomery, p.202.
* Also D. Bressoud, Factorization and Primality Testing, p.179-191 and
* C. Pomerance, Kent State Lecture Notes, 19-22.
*/
{
	if (M->S <= 0 || (M->D == 0 && M->V[0] == 1))
	{
		FREEMPI(N); FREEMPI(Q); FREEMPI(M);
		execerror("M <= 1", "");
	}
	return (LUCASU(N, Q, M));
}

MPI *LUCASBX(MPI *N)
/*
* Let N > 1 be odd. Then LUCASB(N) determines an integer D(!=-3) of the form
* 4m+1,  such that the Jacobi symbol j(D,N)= -1. Then with Q=(1-d)/4, the Lucas
* pseudoprime test examines L=U_{(N+1)/2}mod N. If L=0, N is a Lucas probable
* prime, otherwise N is composite. See "The Pseudoprimes to 25.10^9",
* Mathematics of computation, 35 (1980) 1003-1026. At the end of this paper
* it is conjectured that if N is a strong base 2 pseudoprime and a Lucas
* probable prime, then N is in fact a prime. A $30 prize is offered for a
* counterexample.
* if LUCASB(N) = 1, then N is a Lucas probable prime, else N is composite.
* Returns NULL if N is a perfect square.
*/
{
	if (N->S <= 0 || (N->D == 0 && N->V[0] == 1))
	{
		FREEMPI(N);
		execerror("N <= 1", "");
	}
	if ((N->V[0]) % 2 == 0)
	{
		FREEMPI(N);
		execerror("N is even", "");
	}
	return (LUCASB(N));
}

MPI *LUCASX(MPI *N)
/* Here N is odd and > 1.
* If LUCAS(N) returns 1, then N is a strong base 2 pseudoprime and a Lucas
* probable prime; if LUCAS(N) returns 0, then N is composite.
* See "The Pseudoprimes to 25.10^9", Mathematics of computation, 35 (1980)
* 1003-1026. At the end of this paper it is conjectured that if N is a strong
* base 2 pseudoprime and a Lucas probable prime, then N is in fact a prime.
* A $30 prize is offered for a counterexample.
*/
{
	if (N->S <= 0 || (N->D == 0 && N->V[0] == 1))
	{
		printf("N <= 1\n");
		return NULL;
	}
	if ((N->V[0]) % 2 == 0)
	{
		printf("N is even\n");
		return NULL;
	}
	return (LUCAS(N));
}

MPI *LENGTHX(MPI *N)
/* Returns the number of decimal digits of N. */
{
	unsigned long z;

	z = LENGTHI(N);
	if (N->S < 0)
		z--;
	return (CHANGE(z));
}

MPI *MULT32X(MPI *X, MPI *Y)
{
	USL x, y, z;

	x = CONVERTI(X);
	y = CONVERTI(Y);
	z = MULT32(x, y);
	return (CHANGE(z));
}

MPI *NEXTPRIMEAPX(MPI *A, MPI *B, MPI *M)
/*
* Finds the first Lucas probable prime P, A | P - B, M <= P.
* Here A is even, B is odd, A > B >= 1, gcd(A,B) = 1, M > 1.
*/
{
	MPI *tmp1;
	unsigned int s;

	if (A->V[0] % 2)
	{
		printf("A is odd\n");
		return NULL;
	}
	if (B->V[0] % 2 == 0)
	{
		printf("B is even\n");
		return NULL;
	}
	if (A->S <= 0 || (A->D == 0 && A->V[0] == 1))
	{
		printf("A <= 1\n");
		return NULL;
	}
	if (B->S <= 0)
	{
		printf("B <= 0\n");
		return NULL;
	}
	if (M->S <= 0 || (M->D == 0 && M->V[0] == 1))
	{
		printf("M <= 1\n");
		return NULL;
	}
	if (RSV(A, B) <= 0)
	{
		printf("A <= B\n");
		return NULL;
	}
	tmp1 = GCD(A, B);
	s = EQONEI(tmp1);
	FREEMPI(tmp1);
	if (!s)
	{
		printf("GCD(A,B) > 1\n");
		return NULL;
	}
	return (NEXTPRIMEAP(A, B, M));
}
void SHALLIT()
{
	MPI *K, *N, *tmp, *L, *M;
	unsigned int i;

	K = ONEI();
	L = ONEI();
	VERBOSE = 1;
	for (i = 2; i <= 20; i++)
	{
		M = CHANGE(4 * i - 2);
		tmp = MULTI(M, L);
		FREEMPI(M);
		N = ADD0I(tmp, K);
		FREEMPI(tmp);
		FACTOR(N);
		printf("completed K[%u]\n", i);
		printf("is the factorization of K[%u] = ", i);
		PRINTI(N);
		printf("\n");
		FREEMPI(K);
		K = L;
		L = N;
	}
	VERBOSE = 0;
	FREEMPI(K);
	FREEMPI(L);
}

MPI *EFACTORX(MPI *N, MPI *M, MPI *P)
/* Elliptic curve factoring.
* Warning: 2329 is the number of array elements in PRIMEPOWERS[]
* in primepow.h
*/
{
	MPI *Z;
	unsigned long m, p;

	if (N->S <= 0 || (N->D == 0 && N->V[0] == 1))
	{
		printf("N <= 1\n");
		return NULL;
	}
	if (M->S <= 0 || (M->D == 0 && M->V[0] <= 10))
	{
		printf("M <= 0 or M <= 10\n");
		return NULL;
	}
	if ((M->D == 0 && (M->V[0] > 2328)) || M->D)
	{
		printf("M > 2328\n");
		return NULL;
	}
	if (P->S <= 0 || P->D >= 2)
	{
		printf("P <= 0 or P >= 2^32\n");
		return NULL;
	}
	m = CONVERTI(M);
	p = CONVERTI(P);
	Z = EFACTOR(N, m, p);
	if (Z == NULL)
	{
		printf("Factorization unsuccessful\n");
		return NULL;
	}
	return (Z);
}
void ABS_NEAREST_INTRX()
{
	unsigned int u;
	MPI *tempI, *G, *S;
	MPR *R;

	R = BUILDMPR();
	printf("enter the numerator:");
	R->N = INPUTI(&u);
	Flush();
	printf("enter the denominator (> 0):");
	R->D = INPUTI(&u);
	Flush();
	G = GCD(R->N, R->D);
	tempI = R->N;
	R->N = INT(R->N, G);
	FREEMPI(tempI);
	tempI = R->D;
	R->D = INT(R->D, G);
	FREEMPI(tempI);
	FREEMPI(G);
	S = ABS_NEAREST_INTR(R);
	printf("abs-nearest-int of "); PRINTR(R); printf("=");  PRINTI(S);
	printf("\n");
	FREEMPI(S);
	FREEMPR(R);
	return;
}


MPI *FIBONACCI(USI n)
/*
* Returns the nth Fibonacci number.
*/
{
	MPI *A, *B, *C;
	USI i;

	C = NULL;
	if (n == 0)
		return(ZEROI());
	if (n == 1)
		return(ONEI());
	A = ZEROI();
	B = ONEI();
	for (i = 2; i <= n; i++)
	{
		C = ADD0I(A, B);
		FREEMPI(A);
		A = B;
		B = C;
	}
	FREEMPI(A);
	return(C);
}
MPI *LUCASS(USI n)
/*
* Returns the nth Lucas number.
*/
{
	MPI *A, *B, *C, *THREE;
	USI i;

	C = NULL;
	THREE = BUILDMPI(1);
	THREE->S = 1;
	THREE->V[0] = 3;
	if (n == 0)
		return(ONEI());
	if (n == 1)
		return(THREE);
	A = ZEROI();
	B = ONEI();
	for (i = 2; i <= n; i++)
	{
		C = ADD0I(A, B);
		FREEMPI(A);
		A = B;
		B = C;
	}
	FREEMPI(A);
	return(C);
}

MPI *LPRIME_FACTORX(MPI *N)
/* Returns the least prime factor of N. */
{
	MPI *M;

	if (N->S <= 0) {
		printf("N <= 0\n");
		return NULL;
	}
	if (EQONEI(N)) {
		printf("N = 1\n");
		return NULL;
	}
	M = LPRIME_FACTOR(N);
	if (M == NULL)
	{
		printf("Factorization unsuccessful\n");
		return NULL;
	}
	return (M);
}

MPI *CONWAY_CYCLESX(MPI *A, MPI *B) {
	MPI *T;

	if (A->S <= 0) {
		printf("A <= 0\n");
		return NULL;
	}
	if (B->S <= 0) {
		printf("B <= 0\n");
		return NULL;
	}
	T = CONWAY_CYCLES(A, B);
	if (T == NULL) {
		printf("no cycle found in the range i <= 100000 when initial values are ");
		PRINTI(A);
		printf(", ");
		PRINTI(B);
		printf("\n");
	}
	return(T);
}

MPI *CONWAY_CYCLES0X(MPI *A, MPI *B) {
	MPI *T;
	if (A->S <= 0) {
		printf("A <= 0\n");
		return NULL;
	}
	if (B->S <= 0) {
		printf("B <= 0\n");
		return NULL;
	}
	T = CONWAY_CYCLES0(A, B);
	if (T == NULL) {
		printf("no cycle found in the range i <= 100000 when initial values are ");
		PRINTI(A);
		printf(", ");
		PRINTI(B);
		printf("\n");
	}
	return(T);
}

MPI *CONWAY_RANGE_TESTX(MPI *A, MPI *B) {
	if (A->S <= 0) {
		printf("A <= 0\n");
		return NULL;
	}
	if (B->S <= 0) {
		printf("B <= 0\n");
		return NULL;
	}
	CONWAY_RANGE_TEST(A, B);
	return(ONEI());
}

MPI *CONWAY_RANGE_TEST1X(MPI *M1, MPI *M2, MPI *N1, MPI *N2) {
	if (M1->S <= 0) {
		printf("M1 <= 0\n");
		return NULL;
	}
	if (N1->S <= 0) {
		printf("N1 <= 0\n");
		return NULL;
	}
	if (M2->S <= 0) {
		printf("M2 <= 0\n");
		return NULL;
	}
	if (N2->S <= 0) {
		printf("N2 <= 0\n");
		return NULL;
	}
	CONWAY_RANGE_TEST1(M1, M2, N1, N2);
	return(ONEI());
}
