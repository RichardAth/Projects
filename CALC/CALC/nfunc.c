/* nfunc.c */
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

extern unsigned long PRIME[];
unsigned int MLLLVERBOSE;
unsigned int HERMITEVERBOSE;
unsigned int HERMITE1VERBOSE;
unsigned int GCDVERBOSE;
unsigned int HAVASFLAG = 0;
unsigned int FPRINTMATIFLAG = 0;
unsigned int GCD_MAX = 100000;
extern unsigned int GCDFLAG;

MPI *EUCLIDI(MPI *Pptr, MPI *Qptr, MPI **Hptr, MPI **Kptr)
/*
* gcd(Pptr, Qptr) = Hptr * Pptr + Kptr * Qptr.
*/
{
	MPI *Q, *A, *B, *C, *H1, *H2, *K1, *K2, *L1, *L2, *Tmp1, *Tmp2;
	int s;

	if (Qptr->S == 0)
	{
		s = Pptr->S;
		if (Pptr->S != 0)
		{
			if (s == 1)
				*Hptr = ONEI();
			else
				*Hptr = MINUS_ONEI();
			*Kptr = ZEROI();
			return (ABSI(Pptr));
		}
		else
		{
			*Hptr = ZEROI();
			*Kptr = ZEROI();
			return (ZEROI());
		}
	}
	A = COPYI(Pptr);
	B = ABSI(Qptr);
	C = MOD(A, B);
	s = Qptr->S;
	if (C->S == 0)
	{
		if (s == 1)
			*Kptr = ONEI();
		else
			*Kptr = MINUS_ONEI();
		*Hptr = ZEROI();
		FREEMPI(A);
		FREEMPI(C);
		return (B);
	}
	L1 = ONEI();
	K1 = ZEROI();
	L2 = ZEROI();
	K2 = ONEI();
	while (C->S != 0)
	{
		Q = INT(A, B);
		FREEMPI(A);
		A = B;
		B = C;
		C = MOD0(A, B);
		Tmp1 = MULTI(Q, K1);
		Tmp2 = MULTI(Q, K2);
		FREEMPI(Q);
		H1 = SUBI(L1, Tmp1);
		H2 = SUBI(L2, Tmp2);
		FREEMPI(Tmp1);
		FREEMPI(Tmp2);
		FREEMPI(L1);
		L1 = K1;
		FREEMPI(L2);
		L2 = K2;
		K1 = H1;
		K2 = H2;
	}
	*Hptr = K1;
	if (s == -1)
		K2->S = -(K2->S);
	*Kptr = K2;
	FREEMPI(L1);
	FREEMPI(L2);
	FREEMPI(A);
	FREEMPI(C);
	return (B);
}

void SERRET(MPI *P, MPI **Xptr, MPI **Yptr)
/*
* This program finds positive integers X, Y such that "X*X+Y*Y=P, where P is a
* prime, p=4n+1. The algorithm goes back to Serret and is from the book
* "Computational methods in number theory, part 1," edited by H.W.Lenstra and
* R.Tijdeman.
*/
{
	int i, j;
	MPI *Q, *N, *Z, *U, *V, *W, *Tmp;

	j = 0;
	Q = SUB0_I(P, (USL)1);
	Z = BIG_MTHROOT(Q, 2);
	W = MULTI(Z, Z);
	if (EQUALI(W, Q))
	{
		PRINTI(P); printf(" = "); PRINTI(Z); printf("\n");
		*Xptr = Z;
		*Yptr = ONEI();
		FREEMPI(Q); FREEMPI(W);
		return;
	}
	N = INT_(Q, (USL)4);
	FREEMPI(W); FREEMPI(Q);
	for (i = 0; i <= Y0 - 1; i++)
	{
		printf("PRIME[%d] = %lu\n", i, PRIME[i]);
		U = MPOWER_((long)PRIME[i], N, P);
		V = MULTI(U, U);
		Tmp = V; V = MOD0(V, P); FREEMPI(Tmp);
		Tmp = V; V = ADD0_I(V, (USL)1); FREEMPI(Tmp);
		if (EQUALI(V, P))
		{
			j = 1;
			printf("U = "); PRINTI(U); printf("\n");
			/* U is a square-root of -1 mod P */
			FREEMPI(V); FREEMPI(N);
			break;
		}
		FREEMPI(U); FREEMPI(V);
	}
	if (j == 1)
	{
		CONT_FRAC(P, U, Z, Xptr, Yptr);
		printf("P = X^2 + Y^2, where\n");
		printf("P = "); PRINTI(P); printf("\n");
		printf("X = "); PRINTI(*Xptr); printf("\n");
		printf("Y = "); PRINTI(*Yptr); printf("\n");
		FREEMPI(U); FREEMPI(Z);
		return;
	}
	printf("I don't think that "); PRINTI(P);
	printf(" is a prime of the form 4n+1!\n");
}

void PELL(MPI *Dptr, MPI *Eptr)
/*
* This program finds the period of the regular continued
* fraction expansion of square-root d, as well as the least
* solution x,y of the Pellian equation x*x-d*y*y=+-1.
* The algorithm is from Sierpinski's 'Theory of Numbers', p.296.
* and Davenport's 'The Higher Arithmetic', p.109.
* Here sqrt(d)=a[0]+1/a[1]+...+1/a[n-1]+1/2*a[0]+1/... ,
* The length n of the period a[1],...,a[n-1],2*a[0] is printed.
* The partial quotients are printed iff *Eptr != 0.
*/
{
	MPI *B, *C, *Y, *Z, *A0, *A1, *L, *K, *M, *N, *H, *P, *Tmp;
	unsigned int i, l, m;
	int j;

	Y = BIG_MTHROOT(Dptr, 2);
	A0 = COPYI(Y);
	B = COPYI(Y);
	Z = MULTI(Y, Y);
	if (EQUALI(Dptr, Z))
	{
		printf("You have inputted a perfect square\n");
		printf("The program is aborted\n");
		return;
	}
	C = SUB0I(Dptr, Z);
	A1 = COPYI(C);
	if (Eptr->S)
	{
		printf("\nA[0] = "); PRINTI(Y); printf("\n");
	}
	L = ONEI(); K = COPYI(A0); M = ZEROI(); N = ONEI();
	for (i = 1; 1; i++)
	{
		FREEMPI(Z); Z = ADD0I(B, A0);
		FREEMPI(Y); Y = INT0(Z, C);
		if (Eptr->S)
		{
			printf("A[%u] = ", i); PRINTI(Y); printf("\n");
		}
		FREEMPI(Z); Z = MULTI(Y, C);
		Tmp = B; B = SUB0I(Z, B); FREEMPI(Tmp);
		FREEMPI(Z); Z = MULTI(B, B);
		Tmp = Z; Z = SUB0I(Dptr, Z); FREEMPI(Tmp);
		Tmp = C; C = INT0(Z, C); FREEMPI(Tmp);
		FREEMPI(Z); Z = MULTI(K, Y);
		H = ADDI(Z, L);
		FREEMPI(Z); Z = MULTI(N, Y);
		P = ADDI(Z, M); FREEMPI(M);
		FREEMPI(L); L = K;
		M = N; K = H; N = P;
		if (EQUALI(B, A0))
		{
			if (EQUALI(C, A1))
			{
				printf("The continued fraction for sqrt(");
				PRINTI(Dptr);
				printf(") has period length equal to %u.\n", i);
				printf("Also the least solution (x, y) of x*x - ");
				PRINTI(Dptr);
				j = i % 2 ? -1 : 1;
				printf("y*y = %d\n", j);
				l = LENGTHI(L); m = LENGTHI(M);
				printf(" x has %u digits;\n", l);
				printf(" y has  %u digits;\n", m);
				if (l>500)
					return;
				if (m>500)
					return;
				printf(" and \n");
				printf("x = "); PRINTI(L); printf("\n");
				printf("y = "); PRINTI(M); printf("\n");
				FREEMPI(Z); FREEMPI(L); FREEMPI(M); FREEMPI(K);
				FREEMPI(N); FREEMPI(A0); FREEMPI(A1);
				FREEMPI(B); FREEMPI(C); FREEMPI(Y);
				return;
			}
		}

	}
}

MPI *CONGR(MPI *A, MPI *B, MPI *M, MPI **N)
/*
* Returns the least solution (mod N) of the congruence AX=B(mod M), where
* N = M / gcd(A, M), otherwise returns the null pointer.
*/
{
	MPI *D, *E, *U, *V, *X, *tmp1, *tmp2;
	int t;

	D = EUCLIDI(A, M, &U, &V);
	FREEMPI(V);
	*N = INT(M, D);
	E = MOD(B, D);
	t = E->S;
	FREEMPI(E);
	if (t)
	{
		FREEMPI(U);
		FREEMPI(D);
		return ((MPI *)NULL);
	}
	tmp1 = MULTI(U, B);
	FREEMPI(U);
	tmp2 = INT(tmp1, D);
	X = MOD(tmp2, *N);
	FREEMPI(tmp1);
	FREEMPI(tmp2);
	FREEMPI(D);
	return (X);
}

MPI *CHINESE(MPI *A, MPI *B, MPI *M, MPI *N, MPI **Mptr)
/*
* Returns the solution mod *Mptr=lcm[M,N] of the congruences X = A (mod M)
* and X = B (mod N), if soluble; otherwise returns NULL.
*/
{
	MPI * D, *E, *F, *R, *S, *T, *U, *V;
	int t;

	D = EUCLIDI(M, N, &U, &V);
	S = SUBI(A, B);
	R = MOD(S, D);
	t = R->S;
	FREEMPI(S);
	FREEMPI(R);
	*Mptr = LCM(M, N);
	if (t)
	{
		FREEMPI(D);
		FREEMPI(U);
		FREEMPI(V);
		return ((MPI *)NULL);
	}
	S = MULTI3(B, U, M);
	R = MULTI3(A, V, N);
	FREEMPI(U);
	FREEMPI(V);
	T = ADDI(S, R);
	FREEMPI(S);
	FREEMPI(R);
	E = INT(T, D);
	FREEMPI(D);
	FREEMPI(T);
	F = MOD(E, *Mptr);
	FREEMPI(E);
	return (F);
}

MPI *CHINESE_ARRAY(MPI *A[], MPI *M[], MPI **Mptr, USI n)
/*
* Returns the solution mod *Mptr=lcm[M[0],...,M[n-1] of the congruences
* X = A[i] (mod M[i]),0<=i<n, if soluble; otherwise returns NULL.
*/
{
	MPI *D, *MM, *Z, *tmpMM, *tmpZ, *S, *T, *tmp;
	unsigned int i, j;
	int t;

	for (i = 0; i < n - 1; i++)
		for (j = i + 1; j < n; j++)
		{
			D = GCD(M[i], M[j]);
			S = SUBI(A[i], A[j]);
			T = MOD(S, D);
			t = T->S;
			FREEMPI(D);
			FREEMPI(S);
			FREEMPI(T);
			if (t)
			{
				*Mptr = ZEROI();
				return ((MPI *)NULL);
			}
		}
	MM = COPYI(M[0]);
	Z = COPYI(A[0]);
	for (i = 1; i < n; i++)
	{
		tmpMM = MM;
		tmpZ = Z;
		Z = CHINESE(A[i], Z, M[i], MM, &tmp);
		MM = tmp;
		FREEMPI(tmpMM);
		FREEMPI(tmpZ);
	}
	*Mptr = MM;
	return (Z);
}

MPI *COLLATZT(MPI *Dptr)
{
	MPI *Eptr, *one, *Tmp;

	if (Dptr->V[0] & 1)
	{
		one = ONEI();
		Eptr = MULT_I(Dptr, 3L);
		Tmp = Eptr;
		Eptr = ADDI(Tmp, one);
		FREEMPI(Tmp);
		Tmp = Eptr;
		Eptr = INT_(Tmp, (USL)2);
		FREEMPI(Tmp);
		FREEMPI(one);
	}
	else
		Eptr = INT_(Dptr, (USL)2);
	return (Eptr);
}

void COLLATZ(MPI *Dptr, MPI *Eptr)
/* The iterates are printed iff *Eptr != 0. */
{
	unsigned int i;
	MPI *X, *Y, *Z, *Temp;

	Y = BUILDMPI(1);
	Y->S = -1;
	Y->V[0] = 5;
	Z = BUILDMPI(1);
	Z->S = -1;
	Z->V[0] = 17;

	X = COPYI(Dptr);
	if (EQZEROI(X))
	{
		printf("starting number = "); PRINTI(Dptr); printf("\n");
		printf("the number of iterations taken to reach 0 is %u\n", 0);
		FREEMPI(X);
		FREEMPI(Y);
		FREEMPI(Z);
		return;
	}
	for (i = 0; 1; i++)
	{
		if (EQONEI(X))
		{
			printf("starting number = "); PRINTI(Dptr); printf("\n");
			printf("the number of iterations taken to reach 1 is %u\n", i);
			break;
		}
		if (EQMINUSONEI(X))
		{
			printf("starting number = "); PRINTI(Dptr); printf("\n");
			printf("the number of iterations taken to reach -1 is %u\n", i);
			break;
		}
		if (EQUALI(X, Y))
		{
			printf("starting number = "); PRINTI(Dptr); printf("\n");
			printf("the number of iterations taken to reach -5 is %u\n", i);
			break;
		}
		if (EQUALI(X, Z))
		{
			printf("starting number = "); PRINTI(Dptr); printf("\n");
			printf("the number of iterations taken to reach -17 is %u\n", i);
			break;
		}
		Temp = X;
		X = COLLATZT(Temp);
		FREEMPI(Temp);
		if (Eptr->S)
		{
			PRINTI(X); printf("\n");
		}
	}
	FREEMPI(X);
	FREEMPI(Y);
	FREEMPI(Z);
	return;
}

MPI  *FUND_UNIT(MPI *D, MPI **Xptr, MPI **Yptr)
/*
* This is a program for finding the fundamental unit of Q(sqrt(D)).
* The algorithm is based on K. Rosen, Elementary number theory
* and its applications, p382, B.A. Venkov, Elementary Number theory, p.62
* and D. Knuth, Art of computer programming, Vol.2, p359, with Pohst's trick
* of using half the period.
* w=(1+sqrt(D))/2 if D=1 (mod 4), w=sqrt(D) otherwise.
* The norm of the fundamental unit (*Xptr) + (*Yptr)*w is returned.
*/
{
	unsigned int i;
	MPI *X, *B, *C, *H, *T, *G, *Y, *F, *E, *L, *M, *N, *U, *V, *tmp;
	MPI *tmp1, *tmp2, *K, *R, *S;
	unsigned int s, t;

	if ((D->D == 0) && (D->V[0]) == 5)
	{
		*Xptr = ZEROI();
		*Yptr = ONEI();
		return (MINUS_ONEI());
	}
	X = BIG_MTHROOT(D, 2);
	B = COPYI(X);
	C = ONEI();
	tmp = SUB0_I(X, (USL)1); H = INT0_(tmp, (USL)2); FREEMPI(tmp);
	tmp = MULT_I(H, 2L); T = ADD0_I(tmp, (USL)1); FREEMPI(tmp);
	s = (D->V[0]) % 4;
	if (s == 1)
	{
		FREEMPI(B); B = COPYI(T);
		tmp = C; C = ADD0_I(C, (USL)1); FREEMPI(tmp); /* C = 2 */
	}
	G = MULTI(X, X); tmp = ADD0_I(G, (USL)1); FREEMPI(G);
	t = EQUALI(D, tmp); FREEMPI(tmp);
	if (t) /* period 1, exceptional case */
	{
		if (s == 1)
		{
			*Xptr = SUB0_I(X, (USL)1);
			*Yptr = TWOI();
			FREEMPI(X);
		}
		else
		{
			*Xptr = X;
			*Yptr = ONEI();
		}
		FREEMPI(B); FREEMPI(C); FREEMPI(T); FREEMPI(H);
		return (MINUS_ONEI());
	}
	tmp = MULTI(T, T); tmp1 = ADD0_I(tmp, (USL)4); FREEMPI(tmp);
	t = EQUALI(D, tmp1); FREEMPI(tmp1);
	if (t) /* period 1, exceptional case */
	{
		*Xptr = H;
		*Yptr = ONEI();
		FREEMPI(B); FREEMPI(C); FREEMPI(T); FREEMPI(X);
		return (MINUS_ONEI());
	}
	FREEMPI(T);
	L = ZEROI(); K = ONEI(); M = ONEI(); N = ZEROI();
	for (i = 0; 1; i++)
	{
		tmp = ADDI(X, B); Y = INT(tmp, C); FREEMPI(tmp);
		F = COPYI(B);
		tmp = B; tmp1 = MULTI(Y, C); B = SUBI(tmp1, B);
		FREEMPI(tmp); FREEMPI(tmp1);
		E = COPYI(C); tmp = C;
		tmp1 = MULTI(B, B); tmp2 = SUBI(D, tmp1); C = INT(tmp2, C);
		FREEMPI(tmp); FREEMPI(tmp1); FREEMPI(tmp2);
		if (i == 0)
		{
			FREEMPI(Y);
			if ((D->V[0]) % 4 == 1)
				Y = COPYI(H);
			else
				Y = COPYI(X);
			FREEMPI(H);
		}
		R = L; S = M;
		tmp = MULTI(K, Y); U = ADDI(tmp, L); FREEMPI(tmp);
		tmp = MULTI(N, Y); V = ADDI(tmp, M); FREEMPI(tmp);
		FREEMPI(Y);
		L = K; K = U;
		M = N; N = V;
		/* U/V is the ith convergent to sqrt(D) or (sqrt(D)-1)/2 */
		if (i)
		{
			if (EQUALI(B, F))
				/*\alpha_H=\alpha_{H+1}, even period 2H */
			{
				tmp = ADDI(U, R); tmp1 = MULTI(M, tmp);
				FREEMPI(tmp);
				if (i % 2 == 0)
					*Xptr = ADD0_I(tmp1, (USL)1);
				else
					*Xptr = SUB0_I(tmp1, (USL)1);
				FREEMPI(tmp1);
				tmp = ADDI(V, S); *Yptr = MULTI(M, tmp);
				FREEMPI(tmp);
				FREEMPI(X);
				FREEMPI(L); FREEMPI(M);
				FREEMPI(U); FREEMPI(V);
				FREEMPI(R); FREEMPI(S);
				FREEMPI(C); FREEMPI(B);
				FREEMPI(E); FREEMPI(F);
				return (ONEI());
			}
			if (EQUALI(C, E))
				/*\beta_H=\beta_{H-1}, odd period 2H-1 */
			{
				tmp = MULTI(U, V); tmp1 = MULTI(L, M);
				*Xptr = ADDI(tmp, tmp1);
				FREEMPI(tmp); FREEMPI(tmp1);
				tmp = MULTI(V, V); tmp1 = MULTI(M, M);
				*Yptr = ADD0I(tmp, tmp1);
				FREEMPI(tmp); FREEMPI(tmp1);
				FREEMPI(X);
				FREEMPI(L); FREEMPI(M);
				FREEMPI(U); FREEMPI(V);
				FREEMPI(R); FREEMPI(S);
				FREEMPI(C); FREEMPI(B);
				FREEMPI(E); FREEMPI(F);
				return (MINUS_ONEI());
			}
		}
		FREEMPI(R); FREEMPI(S); FREEMPI(E); FREEMPI(F);
	}
}

MPI  *PEL(MPI *D, MPI *Eptr, MPI **Xptr, MPI **Yptr)
/*
* This is a program for finding the least solution of Pell's equation
* x*x - D*y*y = +-1.
* The algorithm is based on K. Rosen, Elementary number theory
* and its applications, p382, B.A. Venkov, Elementary Number theory, p.62
* and D. Knuth, Art of computer programming, Vol.2, p359, with Pohst's trick
* of using half the period.
* The norm of the least solution is returned.
* The partial quotients are printed iff *Eptr != 0.
*/
{
	unsigned int i;
	MPI *X, *B, *C, *G, *Y, *F, *E, *L, *M, *N, *U, *V, *tmp;
	MPI *tmp1, *tmp2, *K, *R, *S;
	unsigned int t;
	FILE *outfile;
	char buff[20];

	X = BIG_MTHROOT(D, 2);
	strcpy(buff, "pell.out");
	outfile = fopen(buff, "w");
	if (Eptr->S)
	{
		printf("A[0] = "); PRINTI(X); printf("\n");
		fprintf(outfile, "A[0] = "); FPRINTI(outfile, X); fprintf(outfile, "\n");
	}
	B = COPYI(X);
	C = ONEI();
	G = MULTI(X, X); tmp = ADD0_I(G, (USL)1); FREEMPI(G);
	t = EQUALI(D, tmp); FREEMPI(tmp);
	if (t) /* period 1, exceptional case */
	{
		*Xptr = X;
		*Yptr = ONEI();
		printf("period length 1\n");
		fprintf(outfile, "period length 1\n");
		fclose(outfile);
		FREEMPI(B); FREEMPI(C);
		return (MINUS_ONEI());
	}
	L = ZEROI(); K = ONEI(); M = ONEI(); N = ZEROI();
	for (i = 0; 1; i++)
	{
		tmp = ADDI(X, B); Y = INT(tmp, C); FREEMPI(tmp);
		if (i && Eptr->S)
		{
			printf("A[%u] = ", i); PRINTI(Y); printf("\n");
			fprintf(outfile, "A[%u] = ", i); FPRINTI(outfile, Y); fprintf(outfile, "\n");
		}
		F = COPYI(B);
		tmp = B; tmp1 = MULTI(Y, C); B = SUBI(tmp1, B);
		FREEMPI(tmp); FREEMPI(tmp1);
		E = COPYI(C); tmp = C;
		tmp1 = MULTI(B, B); tmp2 = SUBI(D, tmp1); C = INT(tmp2, C);
		FREEMPI(tmp); FREEMPI(tmp1); FREEMPI(tmp2);
		if (i == 0)
		{
			FREEMPI(Y); Y = COPYI(X);
		}
		R = L; S = M;
		tmp = MULTI(K, Y); U = ADDI(tmp, L); FREEMPI(tmp);
		tmp = MULTI(N, Y); V = ADDI(tmp, M); FREEMPI(tmp);
		FREEMPI(Y);
		L = K; K = U;
		M = N; N = V;
		/* U/V is the ith convergent to sqrt(D) */
		if (i)
		{
			if (EQUALI(B, F))
				/*\alpha_H=\alpha_{H+1}, even period 2H */
			{
				tmp = ADDI(U, R); tmp1 = MULTI(M, tmp);
				FREEMPI(tmp);
				if (i % 2 == 0)
					*Xptr = ADD0_I(tmp1, (USL)1);
				else
					*Xptr = SUB0_I(tmp1, (USL)1);
				FREEMPI(tmp1);
				tmp = ADDI(V, S); *Yptr = MULTI(M, tmp);
				FREEMPI(tmp);
				FREEMPI(X);
				FREEMPI(L); FREEMPI(M);
				FREEMPI(U); FREEMPI(V);
				FREEMPI(R); FREEMPI(S);
				FREEMPI(C); FREEMPI(B);
				FREEMPI(E); FREEMPI(F);
				/*printf("even period length %u\n", 2*i); */
				fprintf(outfile, "even period length %u\n", 2 * i);
				fclose(outfile);
				return (ONEI());
			}
			if (EQUALI(C, E))
				/*\beta_H=\beta_{H-1}, odd period 2H-1 */
			{
				tmp = MULTI(U, V); tmp1 = MULTI(L, M);
				*Xptr = ADDI(tmp, tmp1);
				FREEMPI(tmp); FREEMPI(tmp1);
				tmp = MULTI(V, V); tmp1 = MULTI(M, M);
				*Yptr = ADD0I(tmp, tmp1);
				FREEMPI(tmp); FREEMPI(tmp1);
				FREEMPI(X);
				FREEMPI(L); FREEMPI(M);
				FREEMPI(U); FREEMPI(V);
				FREEMPI(R); FREEMPI(S);
				FREEMPI(C); FREEMPI(B);
				FREEMPI(E); FREEMPI(F);
				/*	printf("odd period length %u\n", 2*i+1); */
				fprintf(outfile, "odd period length %u\n", 2 * i + 1);
				fclose(outfile);
				return (MINUS_ONEI());
			}
		}
		FREEMPI(R); FREEMPI(S); FREEMPI(E); FREEMPI(F);
	}
}

MPIA A_SURD; /* used in REDUCED() */
MPIA U_SURD; /* used in REDUCED() */
MPIA V_SURD; /* used in REDUCED() */
unsigned int REDUCED(MPI *D, MPI *U, MPI *V, USI i)
/*
* This is a function for finding the period of the continued fraction
* expansion of reduced quadratic irrational a=(U+sqrt(D))/V.
* Here D is non-square, 1<(U+sqrt(D))/V, -1<(U-sqrt(D))/V<0.
* The algorithm also assumes that V divides d-U*U and is based on K. Rosen,
* Elementary Number theory and its applications, p.379-381 and Knuth's The art
* of computer programming, Vol. 2, p. 359. The period length is returned if
* a is reduced.
* variable i is created by SURD(D,T,U,V,P_ARRAY,Q_ARRAY) below and indexes
* the ith convergent of (U+T*sqrt(D))/V.
*/
{
	MPI *A, *F, *R, *S, *tmp, *tmp1, *tmp2;
	unsigned int j;
	FILE *outfile;
	char buff[20];

	F = BIG_MTHROOT(D, 2);
	R = COPYI(U); S = COPYI(V);
	strcpy(buff, "surd.out");
	outfile = fopen(buff, "a");
	for (j = i; 1; j++)
	{
		tmp = ADD0I(F, R); A = INT0(tmp, S); FREEMPI(tmp);
		ADD_TO_MPIA(A_SURD, A, j);
		fprintf(outfile, ", A[%u]=", j); FPRINTI(outfile, A);
		fprintf(outfile, ", PERIOD\n");
		tmp = MULTI(A, S); tmp1 = R; R = SUB0I(tmp, R);
		FREEMPI(tmp); FREEMPI(tmp1);
		tmp = S; tmp1 = MULTI(R, R);
		tmp2 = SUB0I(D, tmp1); S = INT0(tmp2, S);
		ADD_TO_MPIA(U_SURD, R, j + 1);
		fprintf(outfile, "P[%u]=", j + 1);
		FPRINTI(outfile, R);
		ADD_TO_MPIA(V_SURD, S, j + 1);
		fprintf(outfile, ", Q[%u]=", j + 1);
		FPRINTI(outfile, S);
		FREEMPI(tmp); FREEMPI(tmp1); FREEMPI(tmp2);
		FREEMPI(A);
		if (EQUALI(U, R) && EQUALI(V, S))
		{
			fprintf(outfile, "\n");
			FREEMPI(R); FREEMPI(S); FREEMPI(F);
			fclose(outfile);
			return (j + 1 - i);
		}
		if (j == R0) {
			FREEMPI(R); FREEMPI(S); FREEMPI(F);
			execerror("j = R0", "");
		}
	}
}

unsigned int SURD(MPI *D, MPI *T, MPI *U, MPI *V, MPIA *AA_SURD, MPIA *UU_SURD, MPIA *VV_SURD, MPIA *P_SURD, MPIA *Q_SURD, USI surd_flag)
/*
* This function uses the continued fraction algorithm expansion in K. Rosen,
* Elementary Number theory and its applications,p.379-381 and Knuth's
* The art of computer programming, Vol.2, p. 359. It locates the first complete
* quotient that is reduced and then uses the function REDUCED(D,U,V,i) above
* to locate and return the period of the continued fraction expansion
* of (U+T*sqrt(D))/V.
* AA_SURD is the sequence of partial quotients up to the end of the period;
* UU_SURD and VV_SDURD are the sequences of U[i] and V[i], where the i-th
* complete quotient is (U[i]+sqrt(D))/V[i];
* P_SURD/Q_SURD give the convergents.
* Output is sent to surd.out.
* Regarding surd_flag (added 29th May 2000). This is needed in PATZ.
* If surd_flag = 1 and the period length is odd, we do an extra period.
*/
{
	MPI *A, *DD, *F, *W, *tmp, *tmp1, *tmp2, *UU, *VV, *X;
	unsigned int i, j, l;
	int z, t;
	FILE *outfile;
	char buff[20];

	A_SURD = BUILDMPIA();
	U_SURD = BUILDMPIA();
	V_SURD = BUILDMPIA();
	F = BIG_MTHROOT(D, 2);
	UU = COPYI(U); VV = COPYI(V);
	z = T->S;
	UU->S = (U->S) * z;
	VV->S = (V->S) * z;
	DD = MULTI3(D, T, T);
	W = ABSI(VV);
	tmp1 = MULTI(UU, UU); tmp2 = SUBI(DD, tmp1); FREEMPI(tmp1);
	tmp1 = MOD(tmp2, W); FREEMPI(tmp2);
	if (tmp1->S)
	{
		tmp2 = DD; DD = MULTI3(DD, VV, VV); FREEMPI(tmp2);
		tmp2 = UU; UU = MULTI(UU, W); FREEMPI(tmp2);
		tmp2 = VV; VV = MULTI(VV, W); FREEMPI(tmp2);
	}
	FREEMPI(W);
	FREEMPI(tmp1);
	FREEMPI(F);
	F = BIG_MTHROOT(DD, 2);
	strcpy(buff, "surd.out");
	outfile = fopen(buff, "w");
	if (V->S == -1)
		fprintf(outfile, "-");
	if (!EQONEI(V))
		fprintf(outfile, "(");
	if (U->S)
		FPRINTI(outfile, U);
	if (T->S >0 && U->S)
		fprintf(outfile, " + ");
	if (!EQONEI(T) && !EQMINUSONEI(T))
	{
		FPRINTI(outfile, T);
		fprintf(outfile, "*");
	}
	if (EQMINUSONEI(T))
		fprintf(outfile, "-");
	fprintf(outfile, "sqrt(");
	FPRINTI(outfile, D); fprintf(outfile, ")");
	if (!EQONEI(V))
		fprintf(outfile, ")");
	if (!EQONEI(V) && !EQMINUSONEI(V))
	{
		fprintf(outfile, "/");
		X = ABSI(V);
		FPRINTI(outfile, X);
		FREEMPI(X);
	}
	fprintf(outfile, ":\n");
	for (j = 0; 1; j++)
	{
		if (j == R0) {
			FREEMPIA(A_SURD);
			FREEMPIA(U_SURD);
			FREEMPIA(V_SURD);
			FREEMPI(DD);
			FREEMPI(F);
			FREEMPI(UU);
			FREEMPI(VV);
			fclose(outfile);
			execerror("j = R0", "");
		}
		ADD_TO_MPIA(U_SURD, UU, j);
		fprintf(outfile, "P[%u]=", j);
		FPRINTI(outfile, UU);
		ADD_TO_MPIA(V_SURD, VV, j);
		fprintf(outfile, ", Q[%u]=", j);
		FPRINTI(outfile, VV);
		if (VV->S > 0 && UU->S > 0 && (RSV(UU, F) <= 0))
		{
			tmp1 = ADD0I(UU, VV);
			z = RSV(tmp1, F);
			FREEMPI(tmp1);
			tmp1 = SUBI(VV, UU);
			t = RSV(tmp1, F);
			FREEMPI(tmp1);
			if (z == 1 && t <= 0)/* (U+sqrt(D))/V is reduced */
			{
				fclose(outfile);
				l = REDUCED(DD, UU, VV, j);
				if (surd_flag && l % 2)
					l = REDUCED(DD, UU, VV, j + l);
				outfile = fopen(buff, "a");
				fprintf(outfile, "period length =  %u\n", l);
				CONVERGENTS(A_SURD, P_SURD, Q_SURD);
				FREEMPI(DD);
				FREEMPI(F);
				FREEMPI(UU);
				FREEMPI(VV);
				fprintf(outfile, "convergents:\n");
				for (i = 0; i < A_SURD->size; i++)
				{
					fprintf(outfile, "A[%u]/B[%u]=", i, i);
					FPRINTI(outfile, (*(P_SURD))->A[i]);
					fprintf(outfile, "/");
					FPRINTI(outfile, (*(Q_SURD))->A[i]);
					fprintf(outfile, "\n");
				}
				/* Actually U_SURD and V_SURD are one unit longer than A_SURD */

				*(AA_SURD) = BUILDMPIA();
				for (i = 0; i < A_SURD->size; i++)
					ADD_TO_MPIA(*(AA_SURD), A_SURD->A[i], i);
				*(UU_SURD) = BUILDMPIA();
				*(VV_SURD) = BUILDMPIA();
				for (i = 0; i < A_SURD->size; i++)
				{
					ADD_TO_MPIA(*(UU_SURD), U_SURD->A[i], i);
					ADD_TO_MPIA(*(VV_SURD), V_SURD->A[i], i);
				}
				FREEMPIA(A_SURD);
				FREEMPIA(U_SURD);
				FREEMPIA(V_SURD);
				fclose(outfile);
				return (l); /* l is the period length */
			}
		}
		tmp1 = ADDI(F, UU);
		if (VV->S > 0)
			A = INT(tmp1, VV);
		else /* See Knuth p. 359 */
		{
			tmp = ONEI();
			tmp2 = ADDI(tmp1, tmp); A = INTI(tmp2, VV);
			FREEMPI(tmp); FREEMPI(tmp2);
		}
		FREEMPI(tmp1);
		fprintf(outfile, ", A[%u]=", j);
		FPRINTI(outfile, A);
		fprintf(outfile, "\n");
		ADD_TO_MPIA(A_SURD, A, j);
		tmp2 = MULTI(A, VV); tmp1 = UU; UU = SUBI(tmp2, UU);
		FREEMPI(A); FREEMPI(tmp2); FREEMPI(tmp1);
		tmp1 = MULTI(UU, UU); tmp2 = SUBI(DD, tmp1); FREEMPI(tmp1);
		tmp1 = VV; VV = INTI(tmp2, VV); FREEMPI(tmp1); FREEMPI(tmp2);
	}
}

MPI * JUGGLERT(MPI *Dptr)
{
	MPI *N, *M;
	unsigned long int t;

	t = MOD0_(Dptr, 3);
	printf("res class mod %lu\n", t);
	if (t == 0)
		N = COPYI(Dptr);
	else if (t == 1)
		N = POWERI(Dptr, 2);
	else
		N = POWERI(Dptr, 13);
	M = BIG_MTHROOT(N, 3);
	FREEMPI(N);
	return (M);
}

void JUGGLER(MPI *Dptr, MPI *Iptr)
{
	unsigned long i, k;
	MPI *X, *Y, *Temp;

	Y = BUILDMPI(1);
	Y->S = 1;
	Y->V[0] = 2;
	k = CONVERTI(Iptr);
	X = COPYI(Dptr);
	for (i = 0; i <= k; i++)
	{
		printf("i = %lu: size = %u\n", i, X->D);
		if (EQUALI(X, Y))
		{
			printf("starting number = "); PRINTI(Dptr); printf("\n");
			printf("the number of iterations taken to reach 2 is %lu\n", i);
			break;
		}
		else if (EQONEI(X))
		{
			printf("starting number = "); PRINTI(Dptr); printf("\n");
			printf("the number of iterations taken to reach 1 is %lu\n", i);
			break;
		}
		Temp = X;
		/*printf("i = %lu, LENGTHI(X[%lu])=%lu\n", i, i, LENGTHI(X));*/
		X = JUGGLERT(Temp);
		FREEMPI(Temp);
		/*	PRINTI(X);printf("\n"); */
	}
	FREEMPI(X);
	FREEMPI(Y);
	return;
}

void HERMITE()
{
	USI *alpha, nz, answer;
	MPMATI *MATI1, *MATI2, *MATI3, *MATI4;
	FILE *outfile;
	char buff[20];

	printf("Do you wish to use an existing matrix from a file? (Y/N)\n");
	answer = GetYN();
	if (answer)
	{
		MATI1 = FINPUTMATFILEI_I();
		if (MATI1 == NULL)
			exit(0);
	}
	else
	{
		printf("enter the matrix of integers (first row non-zero) :\n");
		MATI1 = INPUTMATI();
	}
	printf("The matrix entered is\n\n");
	PRINTMATI(0, MATI1->R - 1, 0, MATI1->C - 1, MATI1);
	printf("Do you want the unimodular transformation matrix P? (Y/N)\n");
	answer = GetYN();
	if (!answer)
	{
		alpha = KB_ROW(MATI1, &nz);
		MATI2 = HERMITE1(MATI1, alpha, nz);
		printf("The Hermite normal form = \n");
		PRINTMATI(0, MATI2->R - 1, 0, MATI2->C - 1, MATI2);
		strcpy(buff, "hermite.out");
		outfile = fopen(buff, "w");
		FPRINTMATI(outfile, 0, MATI2->R - 1, 0, MATI2->C - 1, MATI2);
		fclose(outfile);
	}
	else
	{
		alpha = KB_ROWP(MATI1, &MATI3, &nz);
		MATI2 = HERMITE1P(MATI1, MATI3, &MATI4, alpha, nz);
		FREEMATI(MATI3);
		printf("The Hermite normal form = \n");
		PRINTMATI(0, MATI2->R - 1, 0, MATI2->C - 1, MATI2);
		strcpy(buff, "hermite.out");
		outfile = fopen(buff, "w");
		FPRINTMATI(outfile, 0, MATI2->R - 1, 0, MATI2->C - 1, MATI2);
		fclose(outfile);
		printf("The unimodular transformation matrix P = \n");
		PRINTMATI(0, MATI4->R - 1, 0, MATI4->C - 1, MATI4);
		strcpy(buff, "hermitep.out");
		outfile = fopen(buff, "w");
		FPRINTMATI(outfile, 0, MATI4->R - 1, 0, MATI4->C - 1, MATI4);
		fprintf(outfile, "%u", nz);
		fclose(outfile);
		FREEMATI(MATI4);
	}
	FREEMATI(MATI2);
	FREEMATI(MATI1);
	ffree((char *)alpha, (MATI1->C) * sizeof(USI));
	return;
}

void MLLL()
{
	USI answer, m1, n1;
	int s;
	MPMATI *MATI1, *MATI2, *MATI3;
	char buff[20];
	FILE *outfile;

	printf("Do you wish to use an existing matrix from a file? (Y/N)\n");
	answer = GetYN();
	if (answer)
	{
		MATI1 = FINPUTMATFILEI_I();
		if (MATI1 == NULL)
			exit(0);
	}
	else
	{
		printf("enter the matrix of integers (first row non-zero) :\n");
		MATI1 = INPUTMATI();
	}
	printf("The matrix entered is\n\n");
	PRINTMATI(0, MATI1->R - 1, 0, MATI1->C - 1, MATI1);
	MLLLVERBOSE = 0;
	printf("VERBOSE? (Y/N)\n");
	answer = GetYN();
	if (answer)
		MLLLVERBOSE = 1;
	MATI3 = IDENTITYI(MATI1->R);
	printf("enter the parameters m1 and n1 (normally 1 and 1) :");
	s = scanf("%u %u", &m1, &n1);
	Flush();
	MATI2 = BASIS_REDUCTION(MATI1, &MATI3, 0, m1, n1);
	printf("\n\n");
	printf("The corresponding transformation matrix is\n");
	PRINTMATI(0, MATI3->R - 1, 0, MATI3->C - 1, MATI3);
	strcpy(buff, "mllltran.out");
	outfile = fopen(buff, "w");
	FPRINTMATI(outfile, 0, MATI3->R - 1, 0, MATI3->C - 1, MATI3);
	fclose(outfile);
	printf("The corresponding reduced basis is\n");
	PRINTMATI(0, MATI2->R - 1, 0, MATI2->C - 1, MATI2);
	strcpy(buff, "mlllbas.out");
	outfile = fopen(buff, "w");
	FPRINTMATI(outfile, 0, MATI2->R - 1, 0, MATI2->C - 1, MATI2);
	fclose(outfile);
	FREEMATI(MATI1);
	FREEMATI(MATI2);
	FREEMATI(MATI3);
	return;
}

void SMITH()
{
	MPI **M, *MAX_ENTRY;
	unsigned int i, j, u, r, z, answer, m1, n1;
	int s;
	MPMATI *MATI1, *MATI2, *MATI3, *Temp, *Temp1;
	char buff[20];
	FILE *outfile;

	printf("This program takes as input a rectangular matrix A of integers and computes\n");
	printf("unimodular integer matrices P, Q such that PAQ=D, where D is a diagonal matrix\n");
	printf("with positive integer diagonal elements d[1],...,d[r].\n");
	printf("Here d[i] divides d[i+1], 1<=i<=r-1.\n");
	printf("d[1],...,d[r] are called the invariant factors of A.\n\n");
	printf("Do you wish to use an existing matrix from a file? (Y/N)\n\n");
	answer = GetYN();
	if (answer)
	{
		MATI1 = FINPUTMATFILEI_I();
		if (MATI1 == NULL)
			exit(0);
	}
	else
	{
		printf("enter the matrix of integers:\n");
		MATI1 = INPUTMATI();
	}

	printf(" The matrix entered is A = \n\n");
	PRINTMATI(0, MATI1->R - 1, 0, MATI1->C - 1, MATI1);
	z = MIN((int)MATI1->R, (int)MATI1->C);
	printf("enter MAX_ENTRY:");
	MAX_ENTRY = INPUTI(&u);
	Flush();
	printf("The test matrix is %u x %u\n\n", MATI1->R, MATI1->C);
	printf("Enter alpha=m1/n1: select m1 and n1 (1/4 < alpha <= 1) :");
	s = scanf("%u %u", &m1, &n1);
	Flush();
	printf("MAX_ENTRY CUTOFF= "); PRINTI(MAX_ENTRY); printf("\n\n");
	MLLLVERBOSE = 0;
	printf("VERBOSE? (Y/N)\n\n");
	answer = GetYN();
	if (answer)
		MLLLVERBOSE = 1;
	printf("Do you want P and Q? (Y/N)\n\n");
	answer = GetYN();
	if (answer)
	{
		M = SMITHI(MATI1, &MATI2, &MATI3, &r, MAX_ENTRY, m1, n1);
		FREEMPI(MAX_ENTRY);
		printf("The row unit matrix is P = \n");
		PRINTMATI(0, MATI2->R - 1, 0, MATI2->C - 1, MATI2);
		strcpy(buff, "smithp.out");
		outfile = fopen(buff, "w");
		FPRINTMATI(outfile, 0, MATI2->R - 1, 0, MATI2->C - 1, MATI2);
		fclose(outfile);
		printf("The column unit matrix is Q = \n");
		PRINTMATI(0, MATI3->R - 1, 0, MATI3->C - 1, MATI3);
		strcpy(buff, "smithq.out");
		outfile = fopen(buff, "w");
		FPRINTMATI(outfile, 0, MATI3->R - 1, 0, MATI3->C - 1, MATI3);
		fclose(outfile);
		printf("Do you want to check that PAQ is actually equal to the quoted SNF?(Y/N)\n\n");
		answer = GetYN();
		if (answer) {
			Temp = MULTMATI(MATI2, MATI1);
			Temp1 = MULTMATI(Temp, MATI3);
			for (i = 0; i < Temp1->R; i++)
			{
				for (j = 0; j < Temp1->C; j++)
				{
					if (i != j && (Temp1->V[i][j])->S)
					{
						fprintf(stderr, "PAQ is not diagonal!\n");
						exit(1);
					}
				}
			}
			for (i = 0; i < r; i++)
			{
				if (!EQUALI(Temp1->V[i][i], M[i]))
				{
					fprintf(stderr, "PAQ is diagonal, but the diagonals are not right!\n");
					exit(1);
				}

			}
			FREEMATI(Temp);
			FREEMATI(Temp1);
		}
		FREEMATI(MATI2);
		FREEMATI(MATI3);
		printf("validated!\n\n");
	}
	else
	{
		M = SMITHI1(MATI1, &r, MAX_ENTRY, m1, n1);
		FREEMPI(MAX_ENTRY);
	}
	MLLLVERBOSE = 0;
	FREEMATI(MATI1);
	strcpy(buff, "smith.out");
	outfile = fopen(buff, "w");
	fprintf(outfile, "%u %u\n", r, 1);
	for (i = 0; i < r; i++)
	{
		printf("invariant factor d[%u] = ", i + 1);
		PRINTI(M[i]);
		FPRINTI(outfile, M[i]);
		fprintf(outfile, "\n");
		printf("\n");
		FREEMPI(M[i]);
	}
	printf("rank A = %u\n", r);
	fclose(outfile);
	ffree((char *)M, z * sizeof(MPI *));
	return;
}

void DECODEX(MPI *Eptr, MPI *Pptr, MPI *Qptr)
/*
* *Eptr is the encrytion key, *Pptr and *Qptr are the RSA primes.
* The deciphering key *Dptr is computed and file "encoded" is decoded.
* The decoded message is sent to the screen and also to the file "decoded".
*/
{
	MPI *Rptr, *Dptr;
	MPI *Tmp1, *Tmp2, *Tmp3;

	Rptr = MULTI(Pptr, Qptr);
	Tmp1 = SUB0_I(Pptr, 1);
	Tmp2 = SUB0_I(Qptr, 1);
	Tmp3 = MULTI(Tmp1, Tmp2);
	printf("(p-1)(q-1) =  :"); PRINTI(Tmp3); printf("\n");
	Dptr = INVERSEM(Eptr, Tmp3);
	printf("d =  :"); PRINTI(Dptr); printf("\n");
	FREEMPI(Tmp1);
	FREEMPI(Tmp2);
	FREEMPI(Tmp3);
	DECODE(Dptr, Rptr);
	FREEMPI(Dptr);
	FREEMPI(Rptr);
	return;
}
/*
MPI *RSAEX()
{
unsigned int u;
MPI *Pptr, *Qptr, *Eptr;

printf("enter a prime p > 355142:" );
Pptr = INPUTI(&u);
printf("enter another prime q > 355142:" );
Qptr = INPUTI(&u);
Flush();
Eptr = RSAE(Pptr, Qptr);
FREEMPI(Pptr);
FREEMPI(Qptr);
return (Eptr);
}
*/

MPI *GCD_ARRAYV(MPIA M, MPIA *Y)
/*
* n is length of array
* Returns d=gcd(M[0],...,M[n-1]) and an array Y[] of MPI's such that
* d = M[0]*Y[0]+...+M[n-1]*Y[n-1]. Here n > 1.
*/
{
	MPI *D;
	int i;
	unsigned long n = M->size;
	MPMATI *Temp, *Temp1, *Q;

	Temp = BUILDMATI(n, 1);
	for (i = 0; i < n; i++)
		Temp->V[i][0] = COPYI(M->A[i]);
	Temp1 = EXTGCD(Temp, &D, &Q, 3, 4);
	printf("The multipliers are\n");
	PRINTMATI(0, Temp1->R - 1, 0, Temp1->C - 1, Temp1);
	*Y = BUILDMPIA();
	for (i = 0; i < n; i++)
		ADD_TO_MPIA(*Y, Temp1->V[0][i], i);
	FREEMATI(Temp);
	FREEMATI(Temp1);
	FREEMATI(Q);
	return (D);

	/*
	B = (MPI **)mmalloc((USL)(n * sizeof(MPI *)));
	for (i = 0; i < n; i++)
	B[i] = COPYI(M[i]);
	for (i = 1; i < n; i++)
	{
	tmp = B[i];
	B[i] = GCD(B[i], B[i - 1]);
	FREEMPI(tmp);
	}
	D = COPYI(B[n - 1]);
	tmp = B[n - 1];
	B[n - 1] = COPYI(M[n - 1]);
	FREEMPI(tmp);
	k = ONEI();
	for (i = n - 1; i >= 1; i--)
	{
	G = EUCLIDI(B[i], B[i - 1], &U, &V);
	FREEMPI(G);
	tmp = B[i];
	B[i] = MULTI(U, k);
	FREEMPI(U);
	FREEMPI(tmp);
	if ( i == 1)
	break;
	tmp = k;
	k = MULTI(k, V);
	FREEMPI(V);
	FREEMPI(tmp);
	tmp = B[i - 1];
	B[i - 1] = COPYI(M[i - 1]);
	FREEMPI(tmp);
	}
	tmp = B[0];
	B[0] = MULTI(k, V);
	FREEMPI(tmp);
	FREEMPI(k);
	FREEMPI(V);
	*Y = B;
	return (D);
	*/
}

void EXTGCDX()
/* This is algorthm 1 of Havas, Majewski, Matthews. */
{
	USI answer, m1, n1, m, n, i, j, p;
	MPMATI *MATI1, *MATI2, *Q, *M, *MATI3;
	char buff[20];
	FILE *outfile;
	MPI *A, *T, **XX, **X, *Temp;
	int e, s;

	HAVASFLAG = 1;
	printf("Do you wish to enter your sequence of numbers from an existing column matrix? (Y/N)\n");
	answer = GetYN();
	if (answer)
	{
		MATI1 = FINPUTMATFILEI_I();
		if (MATI1 == NULL)
			exit(0);
	}
	else
	{
		printf("WARNING: Make sure the first integer in the sequence is non-zero) :\n");
		MATI1 = INPUTMATI();
	}
	printf("The matrix entered is\n\n");
	PRINTMATI(0, MATI1->R - 1, 0, MATI1->C - 1, MATI1);
	MLLLVERBOSE = 0;
	printf("VERBOSE? (Y/N)\n");
	answer = GetYN();
	printf("answer = %d\n", answer);
	if (answer)
		MLLLVERBOSE = 1;
	printf("MLLLVERBOSE = %u\n", MLLLVERBOSE);
	printf("to enter alpha=m1/n1: select m1 and n1 (normally 1 and 1) :");
	s = scanf("%u %u", &m1, &n1);
	Flush();
	MATI3 = EXTGCD(MATI1, &A, &MATI2, m1, n1);
	FREEMATI(MATI1);
	printf("The multiplier matrix found by LLL is\n");
	PRINTMATI(0, MATI2->R - 1, 0, MATI2->C - 1, MATI2);
	printf("gcd = "); PRINTI(A); printf("\n");
	FREEMPI(A);
	printf("The multipliers found by LLL are\n");
	PRINTMATI(0, MATI3->R - 1, 0, MATI3->C - 1, MATI3);
	m = MATI2->C;
	T = LENGTHSQRI(MATI3, 0);
	printf("The Euclidean norm squared = ");
	PRINTI(T);
	FREEMPI(T);
	printf("\n");
	printf("Do you want to get a shortest multiplier using Fincke_Pohst? (Y/N)\n");
	answer = GetYN();
	p = m - 1;
	XX = (MPI **)mmalloc((USL)(p * sizeof(MPI *)));
	for (j = 0; j < p; j++)
		XX[j] = ZEROI();
	if (answer)
	{
		GCDVERBOSE = 1;
		Q = MATI2;
		n = Q->R;
		while (1)
		{
			M = SHORTESTT0(Q, &X);
			for (j = 0; j < p; j++)
			{
				Temp = XX[j];
				XX[j] = ADDI(XX[j], X[j]);
				FREEMPI(Temp);
				FREEMPI(X[j]);
			}
			ffree((char *)X, p * sizeof(MPI *));
			if (M == NULL)
				break;
			else
			{
				for (j = 0; j < Q->C; j++)
				{
					FREEMPI(Q->V[n - 1][j]);
					Q->V[n - 1][j] = COPYI(M->V[0][j]);
				}
				FREEMATI(MATI3);
				MATI3 = M;
			}
		}
		printf("found a shortest multiplier vector:\n");
	}
	else
		Q = MATI2;
	strcpy(buff, "egcdmat.out");
	outfile = fopen(buff, "w");
	FPRINTMATI(outfile, 0, Q->R - 1, 0, Q->C - 1, Q);
	fclose(outfile);
	strcpy(buff, "egcdmult.out");
	outfile = fopen(buff, "w");
	if (answer)
		fprintf(outfile, "A Shortest multiplier is ");
	else
	{
		fprintf(outfile, "A not necessarily shortest multiplier is ");
		fprintf(outfile, "b[%u]=", m);
	}
	if (answer)
	{
		printf("b[%u]", m);
		fprintf(outfile, "b[%u]", m);
		for (j = 0; j < p; j++)
		{
			e = XX[j]->S;
			if (e == -1)
			{
				printf("+");
				fprintf(outfile, "+");
				Temp = MINUSI(XX[j]);
				if (!EQONEI(Temp))
				{
					PRINTI(Temp);
					FPRINTI(outfile, Temp);
				}
				printf("b[%u]", j + 1);
				fprintf(outfile, "b[%u]", j + 1);
				FREEMPI(Temp);
			}
			if (e == 1)
			{
				printf("-");
				fprintf(outfile, "-");
				if (!EQONEI(XX[j]))
				{
					PRINTI(XX[j]);
					FPRINTI(outfile, XX[j]);
				}
				printf("b[%u]", j + 1);
				fprintf(outfile, "b[%u]", j + 1);
			}
		}
	}
	if (answer) {
		printf("=");
		fprintf(outfile, "=");
		for (i = 0; i < MATI3->C; i++)
		{
			PRINTI(MATI3->V[0][i]);
			printf(" ");
			FPRINTI(outfile, MATI3->V[0][i]);
			fprintf(outfile, " ");

		}
		T = LENGTHSQRI(MATI3, 0);
		printf(": ");
		fprintf(outfile, ": ");
		PRINTI(T);
		FPRINTI(outfile, T);
		printf("\n");
		fprintf(outfile, "\n");
		FREEMPI(T);
	}
	else
	{
		for (i = 0; i < MATI3->C; i++)
		{
			FPRINTI(outfile, MATI3->V[0][i]);
			fprintf(outfile, " ");

		}
		T = LENGTHSQRI(MATI3, 0);
		fprintf(outfile, ": ");
		FPRINTI(outfile, T);
		fprintf(outfile, "\n");
		FREEMPI(T);
	}
	fclose(outfile);
	if (answer)
	{
		printf("Do you want to get all the shortest multipliers? (Y/N)\n");
		answer = GetYN();
		if (answer)
			SHORTEST(Q, XX, 1, 1);
	}
	for (j = 0; j < p; j++)
		FREEMPI(XX[j]);
	ffree((char *)XX, p * sizeof(MPI *));
	FREEMATI(MATI3);
	FREEMATI(Q);
	HAVASFLAG = 0;
	GCDVERBOSE = 0;
	return;
}

void FINCKE_POHSTX()
{
	USI answer, i, j, m, n, u;
	MPMATI *MATI1;
	MPMATR *M;
	MPR *C;

	printf("Do you wish to use an existing matrix from a file? (Y/N)\n");
	answer = GetYN();
	if (answer)
	{
		MATI1 = FINPUTMATFILEI_I();
		if (MATI1 == NULL)
			exit(0);
	}
	else
	{
		printf("enter the matrix of integers (first row non-zero) :\n");
		MATI1 = INPUTMATI();
	}
	printf("The matrix entered is\n\n");
	PRINTMATI(0, MATI1->R - 1, 0, MATI1->C - 1, MATI1);
	m = MATI1->R;
	n = MATI1->C;
	M = BUILDMATR(m, n);
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
		{
			elt(M, i, j) = BUILDMPR();
			elt(M, i, j)->N = COPYI(MATI1->V[i][j]);
			elt(M, i, j)->D = ONEI();
		}
	}
	FREEMATI(MATI1);
	printf("Enter the upper bound for Q(X):");
	C = INPUTR(&u);
	FINCKE_POHST(M, C, 1);
	FREEMATR(M);
	FREEMPR(C);
	return;
}

void IMPROVEPX()
{
	USI answer, m1, n1, norig, m, i;
	MPMATI *MATI1, *MATI2;
	char buff[20];
	FILE *outfile, *infile;
	MPI ***temp;
	int k, s;

	strcpy(buff, "hermitep.out");
	infile = fopen(buff, "r");
	MATI1 = FINPUTMATI(infile);
	s = fscanf(infile, "%u", &m); /* m = no of rows to be improved */
	fclose(infile);
	norig = MATI1->R - m;/* norig = no of zero rows in HNF */
	if (norig == 0)
	{
		printf("HNF has no non-zero rows - P cannot be improved\n");
		FREEMATI(MATI1);
		return;
	}
	temp = (MPI ***)mmalloc(m * sizeof(MPI **));
	for (i = 0; i < m; i++)
		temp[i] = MATI1->V[i];
	for (i = 0; i < norig; i++)
		MATI1->V[i] = MATI1->V[i + m];
	for (i = norig; i < MATI1->R; i++)
		MATI1->V[i] = temp[i - norig];
	printf("The matrix entered is\n\n");
	PRINTMATI(0, MATI1->R - 1, 0, MATI1->C - 1, MATI1);
	MLLLVERBOSE = 0;
	printf("VERBOSE? (Y/N)\n");
	answer = GetYN();
	if (answer)
		MLLLVERBOSE = 1;
	printf("enter the parameters m1 and n1 (normally 1 and 1) :");
	s = scanf("%u %u", &m1, &n1);
	Flush();
	GCDFLAG = 1;
	MATI2 = BASIS_REDUCTION00(MATI1, m1, n1, norig);
	GCDFLAG = 0;
	printf("\n\n");
	strcpy(buff, "improvep.out");
	outfile = fopen(buff, "w");
	for (i = norig; i < MATI2->R; i++)
		temp[i - norig] = MATI2->V[i];
	for (k = norig - 1; k >= 0; k--)
		MATI2->V[k + m] = MATI2->V[k];
	for (i = 0; i < m; i++)
		MATI2->V[i] = temp[i];
	FPRINTMATI(outfile, 0, MATI2->R - 1, 0, MATI2->C - 1, MATI2);
	fclose(outfile);
	printf("The improved transformation matrix is\n");
	PRINTMATI(0, MATI2->R - 1, 0, MATI2->C - 1, MATI2);
	FREEMATI(MATI1);
	FREEMATI(MATI2);
	ffree((char *)temp, m * sizeof(MPI **));
	return;
}

void QSORTMPIX()
/*
* Input: an array A[0],...,A[q] of MPI's.
* sorts nonnegative MPI's  A[m],...,A[n] in order of size.
*/
{
	USI m, n, p, i, u, q;
	int s;
	MPI **A;

	printf("enter 0<=m<= n<=q (q+1= no of elements, starting from A[0])  :");
	s = scanf("%u %u %u", &m, &n, &q);
	p = q + 1;
	A = (MPI **)mmalloc(p * sizeof(MPI *));
	printf("enter the array A[0],...,A[%u] of %u MPIs:", q, p);
	for (i = 0; i < p; i++)
		A[i] = INPUTI(&u);
	for (i = 0; i < p; i++)
	{
		printf("A[%u] = ", i);
		PRINTI(A[i]);
		printf("\n");
	}
	QSORTMPI0(A, m, n);
	for (i = 0; i < p; i++)
	{
		printf("A[%u] = ", i);
		PRINTI(A[i]);
		FREEMPI(A[i]);
		printf("\n");
	}
	ffree((char *)A, p * sizeof(MPI *));
	return;
}

void QSORTMATIX()
{
	USI answer, r, s;
	int ss;
	MPMATI *MATI1;

	printf("Do you wish to use an existing matrix from a file? (Y/N)\n");
	answer = GetYN();
	if (answer)
	{
		MATI1 = FINPUTMATFILEI_I();
		if (MATI1 == NULL)
			exit(0);
	}
	else
	{
		printf("enter the matrix of integers (first row non-zero) :\n");
		MATI1 = INPUTMATI();
	}
	printf("The matrix entered is\n\n");
	PRINTMATI(0, MATI1->R - 1, 0, MATI1->C - 1, MATI1);
	printf("enter the rows to be sorted :");
	ss = scanf("%u %u", &r, &s);
	QSORTMATI(MATI1, r, s);
	printf("The sorted matrix is\n");
	PRINTMATI(0, MATI1->R - 1, 0, MATI1->C - 1, MATI1);
	FREEMATI(MATI1);
	return;
}

void SCHNORRGCD(MPI *N)
{
	USI answer, m1, n1, n, t, i, j;
	int s;
	MPMATI *MATI1, *MATI2, *MATI0;
	char buff[20];
	FILE *outfile;

	HAVASFLAG = 1;
	printf("Do you wish to enter your sequence of numbers from an existing column matrix? (Y/N)\n");
	answer = GetYN();
	if (answer)
	{
		MATI1 = FINPUTMATFILEI_I();
		if (MATI1 == NULL)
			exit(0);
	}
	else
	{
		printf("WARNING: Make sure the first integer in the sequence is non-zero) :\n");
		MATI1 = INPUTMATI();
	}
	printf("The matrix entered is\n\n");
	PRINTMATI(0, MATI1->R - 1, 0, MATI1->C - 1, MATI1);
	MLLLVERBOSE = 0;
	printf("VERBOSE? (Y/N)\n");
	answer = GetYN();
	printf("answer = %d\n", answer);
	if (answer)
		MLLLVERBOSE = 1;
	printf("MLLLVERBOSE = %u\n", MLLLVERBOSE);
	printf("to enter alpha=m1/n1: select m1 and n1 (normally 1 and 1) :");
	s = scanf("%u %u", &m1, &n1);
	Flush();
	n = MATI1->R;
	t = n + 1;
	MATI0 = BUILDMATI(n, t);
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			if (i == j)
				MATI0->V[i][j] = ONEI();
			else
				MATI0->V[i][j] = ZEROI();
		}
		MATI0->V[i][n] = MULTI(MATI1->V[i][0], N);
	}
	/*
	MATI0 = BUILDMATI(t, t + 1);
	for (i = 0; i < n; i++)
	{
	for (j = 0; j < n; j++)
	{
	if ( i == j)
	MATI0->V[i][j] = ONEI();
	else
	MATI0->V[i][j] = ZEROI();
	}
	MATI0->V[i][n] = MULTI(MATI1->V[i][0], N);
	MATI0->V[i][t] = ZEROI();
	}
	for (j = 0; j < n; j++)
	MATI0->V[n][j] = ONEI();
	*/
	FREEMATI(MATI1);
	strcpy(buff, "sgcdmat.out");
	outfile = fopen(buff, "w");
	FPRINTMATI(outfile, 0, MATI0->R - 1, 0, MATI0->C - 1, MATI0);
	fclose(outfile);

	GetReturn();
	MATI2 = BASIS_REDUCTION0(MATI0, m1, n1);
	FREEMATI(MATI0);
	printf("The corresponding reduced basis is\n");
	PRINTMATI(0, MATI2->R - 1, 0, MATI2->C - 1, MATI2);
	strcpy(buff, "sgcdbas.out");
	outfile = fopen(buff, "w");
	FPRINTMATI(outfile, 0, MATI2->R - 1, 0, MATI2->C - 1, MATI2);
	fclose(outfile);
	FREEMATI(MATI2);
	HAVASFLAG = 0;
	return;
}

void SHORTESTTTX()
{
	/* This is the inhomogenous Fincke-Pohst Algorithm .
	* See Pohst-Zassenhaus .
	*/
	USI answer, u;
	MPI *BOUND;
	MPMATI *MATI1;

	printf("Do you wish to enter your matrix from an existing matrix? (Y/N)\n");
	answer = GetYN();
	if (answer)
	{
		MATI1 = FINPUTMATFILEI_I();
		if (MATI1 == NULL)
			exit(0);
	}
	else
		MATI1 = INPUTMATI();
	printf("The matrix entered is\n\n");
	PRINTMATI(0, MATI1->R - 1, 0, MATI1->C - 1, MATI1);
	MLLLVERBOSE = 0;
	printf("VERBOSE? (Y/N)\n");
	answer = GetYN();
	printf("answer = %d\n", answer);
	if (answer)
		MLLLVERBOSE = 1;
	printf("MLLLVERBOSE = %u\n", MLLLVERBOSE);
	printf("enter a positive integer (the distance bound):");
	BOUND = INPUTI(&u);
	Flush();
	if (BOUND->S <= 0)
	{
		printf("BOUND <= 0\n");
		return;
	}
	SHORTESTTT(MATI1, BOUND);
	FREEMATI(MATI1);
	FREEMPI(BOUND);
	return;
}

void FIB_MIN()
{
	USI M1, M2, N1, N2, i, m, n;
	int s;
	MPMATI *MATI1, *MATI2, *Q;
	MPI *A;
	char buff[20];
	FILE *outfile;

	printf("enter positive integers N1 <= N2, the test range:");
	s = scanf("%u%u", &N1, &N2);
	printf("enter integers (>= 2) M1 <= M2, the sequence range:");
	s = scanf("%u%u", &M1, &M2);
	Flush();
	strcpy(buff, "fibmin.out");
	outfile = fopen(buff, "w");
	outfile = fopen(buff, "a+");
	FPRINTMATIFLAG = 1;
	MLLLVERBOSE = 0;
	for (m = M1; m <= M2; m++)
	{
		fprintf(outfile, "====================\n");
		for (n = N1; n <= N2; n++)
		{
			printf("m=%u: n=%u:\n ", m, n);
			fprintf(outfile, "m=%u: n=%u:", m, n);
			MATI1 = BUILDMATI(m, 1);
			for (i = 0; i < MATI1->R; i++)
			{
				MATI1->V[i][0] = FIBONACCI(n + i);
				/*
				FPRINTI(outfile, MATI1->V[i][0]);
				fprintf(outfile, " ");
				*/
			}
			MATI2 = EXTGCD(MATI1, &A, &Q, 3, 4);
			FREEMATI(MATI1);
			FREEMATI(MATI2);
			FREEMPI(A);
			FPRINTMATI(outfile, Q->R - 1, Q->R - 1, 0, Q->C - 1, Q);
			FREEMATI(Q);
		}
	}
	fprintf(outfile, "====================\n");
	fclose(outfile);
	FPRINTMATIFLAG = 0;
	return;
}

void PRINTW1(USI n, USI r)
{
	USI x;
	x = n - r - 2;
	if (n % 2 == 0)
	{
		if (r <= n - 4)
			printf(" L[%u]       ", x);
		if (r == n - 3)
			printf(" 2           ");
		if (r == n - 2)
			printf(" 1           ");
		if (r >= n - 1)
			printf(" 0           ");
	}
	if (n % 2)
	{
		if (r <= n - 3)
			printf(" L[%u]       ", x);
		if (r == n - 2)
			printf(" 1           ");
		if (r == n - 1)
			printf(" 1           ");
		if (r >= n)
			printf(" 0           ");
	}
}

void PRINTW2(USI n, USI r, USI m)
{
	USI x, y, z;
	x = n - r - 2; y = 2 * m - n; z = r - y - 3;
	if (n % 2 == 0)
	{
		if (r < m)
		{
			if (r <= y + 1)
				printf(" L[%u]       ", x);
			if (r == y + 2)
				printf(" L[%u]+1     ", x);
			if (r == y + 3)
				printf(" L[%u]-1     ", x);
			if (r == y + 4)
				printf(" L[%u]+2     ", x);
			if (r == y + 5)
				printf(" L[%u]-2     ", x);
			if (r < m && r >= y + 6 && r % 2 == 0)
				printf(" L[%u]+L[%u] ", x, z);
			if (r < m && r >= y + 6 && r % 2)
				printf(" L[%u]-L[%u] ", x, z);
		}
		if (r == m && n == m + 3)
			printf(" 1           ");
		if (r == m && n == m + 5)
			printf(" 2           ");
		if (r == m && n >= m + 6)
			printf(" L[%u]       ", n - m - 4);
	}
	if (n % 2)
	{
		if (r < m)
		{
			if (r <= y)
				printf(" L[%u]       ", x);
			if (r == y + 1)
				printf(" L[%u]+1     ", x);
			if (r == y + 2)
				printf(" L[%u]-1     ", x);
			if (r == y + 3)
				printf(" L[%u]+1     ", x);
			if (r == y + 4)
				printf(" L[%u]-2     ", x);
			if (r < m && r >= y + 5 && r % 2 == 0)
				printf(" L[%u]+L[%u] ", x, z);
			if (r < m && r >= y + 6 && r % 2)
				printf(" L[%u]-L[%u] ", x, z);
		}
		if (r == m && n == m + 4)
			printf(" 1           ");
		if (r == m && n >= m + 5)
			printf(" L[%u]       ", n - m - 4);
	}
}

void PRINTWW()
{
	USI m, t, w, x, y, z, n, r;
	int s;

	printf("enter m:");
	s = scanf("%u", &m);
	x = m + 3; y = 2 * m, z = x - 2; w = m - 1;
	for (n = 3; n <= z; n++)
	{
		for (r = 1; r <= w; r++)
			PRINTW1(n, r);
		printf(" 0           ");
		printf("\n");
	}
	for (r = 1; r <= m; r++)
	{
		t = m - r;
		if (r < m - 1)
			printf(" L[%u]       ", t);
		else if (r == m - 1)
			printf(" L[%u]+1     ", t);
		else
			printf(" 0           ");
	}
	printf("\n");
	for (n = x; n <= y; n++)
	{
		for (r = 1; r <= m; r++)
			PRINTW2(n, r, m);
		printf("\n");
	}
	return;
}

void PRINT_DEFECT()
{
	USI M1, M2, N1, N2, i, m, n;
	int s;
	MPMATI *MATI, *MATI1, *MATI2, *MATI3, *Q1, *Q2, *Q3, *Q;
	MPI *A;
	char buff[20];
	FILE *outfile;

	printf("enter positive integers N1 <= N2, the test range:");
	s = scanf("%u%u", &N1, &N2);
	printf("enter integers (>= 2) M1 <= M2, the sequence range:");
	s = scanf("%u%u", &M1, &M2);
	Flush();
	strcpy(buff, "defect.out");
	outfile = fopen(buff, "w");
	outfile = fopen(buff, "a+");
	FPRINTMATIFLAG = 1;
	MLLLVERBOSE = 0;
	for (m = M1; m <= M2; m++)
	{
		printf("m=%u:\n ", m);
		fprintf(outfile, "====================\n");
		for (n = N1; n <= N2; n++)
		{
			printf("n=%u:\n ", n);
			fprintf(outfile, "n=%u:", n);
			MATI1 = BUILDMATI(m, 1);
			MATI2 = BUILDMATI(m, 1);
			MATI3 = BUILDMATI(m, 1);
			for (i = 0; i < MATI1->R; i++)
			{
				MATI1->V[i][0] = FIBONACCI(n + i);
				MATI2->V[i][0] = FIBONACCI(n + 1 + i);
				MATI3->V[i][0] = FIBONACCI(n + 2 + i);
			}
			MATI = EXTGCD(MATI1, &A, &Q1, 3, 4);
			FREEMATI(MATI);
			FREEMATI(MATI1);
			FREEMPI(A);
			MATI = EXTGCD(MATI2, &A, &Q2, 3, 4);
			FREEMATI(MATI);
			FREEMATI(MATI2);
			FREEMPI(A);
			MATI = EXTGCD(MATI3, &A, &Q3, 3, 4);
			FREEMATI(MATI);
			FREEMATI(MATI3);
			FREEMPI(A);
			MATI = ADDMATI(Q3, Q2);
			Q = SUBMATI(MATI, Q1);
			FPRINTMATI(outfile, Q->R - 1, Q->R - 1, 0, Q->C - 1, Q);
			FREEMATI(Q1);
			FREEMATI(Q2);
			FREEMATI(Q3);
			FREEMATI(MATI);
		}
	}
	fprintf(outfile, "====================\n");
	fclose(outfile);
	FPRINTMATIFLAG = 0;
	return;
}

void LUCAS_MIN()
{
	USI M1, M2, N1, N2, i, m, n;
	int s;
	MPMATI *MATI1, *MATI2, *Q;
	MPI *A;
	char buff[20];
	FILE *outfile;

	printf("enter positive integers N1 <= N2, the test range:");
	s = scanf("%u%u", &N1, &N2);
	printf("enter integers (>= 2) M1 <= M2, the sequence range:");
	s = scanf("%u%u", &M1, &M2);
	Flush();
	strcpy(buff, "lucasmin.out");
	outfile = fopen(buff, "w");
	outfile = fopen(buff, "a+");
	FPRINTMATIFLAG = 1;
	MLLLVERBOSE = 0;
	for (m = M1; m <= M2; m++)
	{
		fprintf(outfile, "====================\n");
		for (n = N1; n <= N2; n++)
		{
			printf("m=%u: n=%u:\n ", m, n);
			fprintf(outfile, "m=%u: n=%u:", m, n);
			MATI1 = BUILDMATI(m, 1);
			for (i = 0; i < MATI1->R; i++)
			{
				MATI1->V[i][0] = LUCASS(n + i);
				/*
				FPRINTI(outfile, MATI1->V[i][0]);
				fprintf(outfile, " ");
				*/
			}
			MATI2 = EXTGCD(MATI1, &A, &Q, 3, 4);
			FREEMATI(MATI1);
			FREEMATI(MATI2);
			FREEMPI(A);
			FPRINTMATI(outfile, Q->R - 1, Q->R - 1, 0, Q->C - 1, Q);
			FREEMATI(Q);
		}
	}
	fprintf(outfile, "====================\n");
	fclose(outfile);
	FPRINTMATIFLAG = 0;
	return;
}

void LLLGCDX()
/*
* This executes the LLLGCD algorithm of Havas, Majewski and Matthews.
* 30/1/97.
*/
{
	USI answer, m1, n1, m, n, i, j, p;
	int e, s;
	MPMATI *MATI1, *MATI2, *Q, *M, *MATI3, *MATI0;
	char buff[20];
	FILE *outfile;
	MPI *A, *T, **XX, **X, *Temp;

	HAVASFLAG = 1;
	printf("Do you wish to enter your sequence of numbers from an existing column matrix? (Y/N)\n");
	answer = GetYN();
	if (answer)
	{
		MATI1 = FINPUTMATFILEI_I();
		if (MATI1 == NULL)
			exit(0);
	}
	else
	{
		printf("WARNING: Make sure the first integer in the sequence is non-zero) :\n");
		MATI1 = INPUTMATI();
	}
	printf("The file entered is\n\n");
	PRINTMATI(0, MATI1->R - 1, 0, MATI1->C - 1, MATI1);
	MATI0 = COPYMATI(MATI1);
	GCDVERBOSE = 0;
	printf("VERBOSE? (Y/N)\n");
	answer = GetYN();
	if (answer)
		GCDVERBOSE = 1;
	printf("to enter alpha=m1/n1: select m1 and n1 (normally 1 and 1) :");
	s = scanf("%u %u", &m1, &n1);
	Flush();
	m = MATI1->R;
	MATI2 = LLLGCD(MATI1, &A, m, m1, n1);
	FREEMATI(MATI1);
	printf("The multiplier matrix found by LLL is\n");
	PRINTMATI(0, MATI2->R - 1, 0, MATI2->C - 1, MATI2);
	printf("gcd = "); PRINTI(A); printf("\n");
	FREEMPI(A);
	printf("The multiplier vector found by LLL is b[%u]=", m);
	for (i = 0; i < MATI2->R; i++)
	{
		PRINTI(MATI2->V[MATI2->R - 1][i]);
		printf(" ");
	}
	printf("\n");
	MATI3 = BUILDMATI(1, m);
	for (i = 0; i < m; i++)
		MATI3->V[0][i] = COPYI(MATI2->V[m - 1][i]);
	T = LENGTHSQRI(MATI3, 0);
	printf("The Euclidean norm squared = ");
	PRINTI(T);
	FREEMPI(T);
	printf("\n");
	printf("Do you want to get a shortest multiplier using Fincke_Pohst? (Y/N)\n");
	answer = GetYN();
	p = m - 1;
	XX = (MPI **)mmalloc((USL)(p * sizeof(MPI *)));
	for (j = 0; j < p; j++)
		XX[j] = ZEROI();
	if (answer)
	{
		GCDVERBOSE = 1;
		Q = MATI2;
		n = Q->R;
		while (1)
		{
			M = SHORTESTT0(Q, &X);
			for (j = 0; j < p; j++)
			{
				Temp = XX[j];
				XX[j] = ADDI(XX[j], X[j]);
				FREEMPI(Temp);
				FREEMPI(X[j]);
			}
			ffree((char *)X, p * sizeof(MPI *));
			if (M == NULL)
				break;
			else
			{
				for (j = 0; j < Q->C; j++)
				{
					FREEMPI(Q->V[n - 1][j]);
					Q->V[n - 1][j] = COPYI(M->V[0][j]);
				}
				FREEMATI(MATI3);
				MATI3 = M;
			}
		}
		printf("found a shortest multiplier vector:\n");
	}
	else
		Q = MATI2;
	strcpy(buff, "lllgcdmat.out");
	outfile = fopen(buff, "w");
	FPRINTMATI(outfile, 0, Q->R - 1, 0, Q->C - 1, Q);
	printf("The numbers were  \n");
	PRINTMATI(0, MATI0->R - 1, 0, MATI0->C - 1, MATI0);
	fprintf(outfile, "The numbers were  \n");
	FPRINTMATI(outfile, 0, MATI0->R - 1, 0, MATI0->C - 1, MATI0);
	FREEMATI(MATI0);
	fclose(outfile);
	strcpy(buff, "lllgcdmult.out");
	outfile = fopen(buff, "w");
	if (answer)
		fprintf(outfile, "A Shortest multiplier is ");
	else
	{
		fprintf(outfile, "A not necessarily shortest multiplier is ");
		fprintf(outfile, "b[%u]=", m);
	}
	if (answer)
	{
		printf("b[%u]", m);
		fprintf(outfile, "b[%u]", m);
		for (j = 0; j < p; j++)
		{
			e = XX[j]->S;
			if (e == -1)
			{
				printf("+");
				fprintf(outfile, "+");
				Temp = MINUSI(XX[j]);
				if (!EQONEI(Temp))
				{
					PRINTI(Temp);
					FPRINTI(outfile, Temp);
				}
				printf("b[%u]", j + 1);
				fprintf(outfile, "b[%u]", j + 1);
				FREEMPI(Temp);
			}
			if (e == 1)
			{
				printf("-");
				fprintf(outfile, "-");
				if (!EQONEI(XX[j]))
				{
					PRINTI(XX[j]);
					FPRINTI(outfile, XX[j]);
				}
				printf("b[%u]", j + 1);
				fprintf(outfile, "b[%u]", j + 1);
			}
		}
	}
	if (answer) {
		printf("=");
		fprintf(outfile, "=");
		for (i = 0; i < MATI3->C; i++)
		{
			PRINTI(MATI3->V[0][i]);
			printf(" ");
			FPRINTI(outfile, MATI3->V[0][i]);
			fprintf(outfile, " ");

		}
		T = LENGTHSQRI(MATI3, 0);
		printf(": ");
		fprintf(outfile, ": ");
		PRINTI(T);
		FPRINTI(outfile, T);
		printf("\n");
		fprintf(outfile, "\n");
		FREEMPI(T);
	}
	else
	{
		for (i = 0; i < MATI3->C; i++)
		{
			FPRINTI(outfile, MATI3->V[0][i]);
			fprintf(outfile, " ");

		}
		T = LENGTHSQRI(MATI3, 0);
		fprintf(outfile, ": ");
		FPRINTI(outfile, T);
		fprintf(outfile, "\n");
		FREEMPI(T);
	}
	fclose(outfile);
	if (answer)
	{
		printf("Do you want to get all the shortest multipliers? (Y/N)\n");
		answer = GetYN();
		if (answer)
			SHORTEST(Q, XX, 2, 1);
	}
	for (j = 0; j < p; j++)
		FREEMPI(XX[j]);
	ffree((char *)XX, p * sizeof(MPI *));
	FREEMATI(MATI3);
	FREEMATI(Q);
	GCDVERBOSE = 0;
	HAVASFLAG = 0;
	return;
}

void JACOBIGCDX()
/*
* This executes the LLLGCD algorithm of Havas, Majewski and Matthews.
* 30/1/97.
*/
{
	USI answer, m, i;
	MPMATI *MATI1, *MATI2, *MATI3;
	MPI *A, *T;

	HAVASFLAG = 1;
	printf("Do you wish to enter your sequence of numbers from an existing column matrix? (Y/N)\n");
	answer = GetYN();
	if (answer)
	{
		MATI1 = FINPUTMATFILEI_I();
		if (MATI1 == NULL)
			exit(0);
	}
	else
	{
		printf("WARNING: Make sure the first integer in the sequence is non-zero) :\n");
		MATI1 = INPUTMATI();
	}
	printf("The matrix entered is\n\n");
	PRINTMATI(0, MATI1->R - 1, 0, MATI1->C - 1, MATI1);
	m = MATI1->R;
	MLLLVERBOSE = 0;
	printf("VERBOSE? (Y/N)\n");
	answer = GetYN();
	printf("answer = %d\n", answer);
	if (answer)
		MLLLVERBOSE = 1;
	printf("MLLLVERBOSE = %u\n", MLLLVERBOSE);
	MATI2 = JACOBIGCD(MATI1, &A, m);
	MLLLVERBOSE = 0;
	FREEMATI(MATI1);
	HAVASFLAG = 0;
	printf("The multiplier matrix found by LLL is\n");
	PRINTMATI(0, MATI2->R - 1, 0, MATI2->C - 1, MATI2);
	printf("gcd = "); PRINTI(A); printf("\n");
	FREEMPI(A);
	printf("The multipliers found by JACOBI are\n");
	for (i = 0; i < MATI2->R; i++)
	{
		PRINTI(MATI2->V[m - 1][i]);
		printf(" ");
	}
	printf("\n");
	MATI3 = BUILDMATI(1, m);
	for (i = 0; i < m; i++)
		MATI3->V[0][i] = COPYI(MATI2->V[m - 1][i]);
	T = LENGTHSQRI(MATI3, 0);
	printf("The Euclidean norm squared = ");
	PRINTI(T);
	FREEMPI(T);
	FREEMATI(MATI2);
	FREEMATI(MATI3);
	printf("\n");

}

void SCHNORRHERMITE(MPI *N)
{
	USI answer, m1, n1, m, n, t, i, j;
	int s;
	MPMATI *MATI1, *MATI2, *MATI0;
	MPI *Temp;
	char buff[20];
	FILE *outfile;

	printf("Do you wish to enter your sequence of numbers from an existing column matrix? (Y/N)\n");
	answer = GetYN();
	if (answer)
	{
		MATI1 = FINPUTMATFILEI_I();
		if (MATI1 == NULL)
			exit(0);
	}
	else
	{
		printf("WARNING: Make sure the first integer in the sequence is non-zero) :\n");
		MATI1 = INPUTMATI();
	}
	printf("The matrix entered is\n\n");
	PRINTMATI(0, MATI1->R - 1, 0, MATI1->C - 1, MATI1);
	MLLLVERBOSE = 0;
	HERMITEVERBOSE = 0;
	printf("VERBOSE? (Y/N)\n");
	answer = GetYN();
	printf("answer = %d\n", answer);
	if (answer)
	{
		MLLLVERBOSE = 1;
		HERMITEVERBOSE = 1;
	}
	printf("MLLLVERBOSE = %u\n", MLLLVERBOSE);
	printf("to enter alpha=m1/n1: select m1 and n1 (normally 1 and 1) :");
	s = scanf("%u %u", &m1, &n1);
	Flush();
	n = MATI1->R;
	m = MATI1->C;
	t = n + m;
	MATI0 = BUILDMATI(n, t);
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			if (i == j)
				MATI0->V[i][j] = ONEI();
			else
				MATI0->V[i][j] = ZEROI();
		}
		for (j = n; j < t; j++)
		{
			Temp = POWERI(N, t - j);
			MATI0->V[i][j] = MULTI(MATI1->V[i][j - n], Temp);
			FREEMPI(Temp);
		}
	}
	FREEMATI(MATI1);
	strcpy(buff, "shermitemat.out");
	outfile = fopen(buff, "w");
	FPRINTMATI(outfile, 0, MATI0->R - 1, 0, MATI0->C - 1, MATI0);
	fclose(outfile);
	MATI2 = BASIS_REDUCTION000(MATI0, m1, n1, N);
	HERMITEVERBOSE = 0;
	FREEMATI(MATI0);
	for (i = 0; i < n; i++)
	{
		for (j = n; j < t; j++)
		{
			Temp = POWERI(N, t - j);
			MATI2->V[i][j] = INTI(MATI2->V[i][j], Temp);
			FREEMPI(Temp);
		}
	}
	printf("The corresponding reduced basis is\n");
	PRINTMATI(0, MATI2->R - 1, 0, MATI2->C - 1, MATI2);
	strcpy(buff, "shermitebas.out");
	outfile = fopen(buff, "w");
	FPRINTMATI(outfile, 0, MATI2->R - 1, 0, MATI2->C - 1, MATI2);
	fclose(outfile);
	FREEMATI(MATI2);
	return;
}

void LLLHERMITE1X()
{
	USI answer, rank, m1, n1;
	int s;
	MPMATI *MATI1, *MATI2, *MATI3;
	char buff[20];
	FILE *outfile;

	printf("Do you wish to use an existing matrix from a file? (Y/N)\n");
	answer = GetYN();
	if (answer)
	{
		MATI1 = FINPUTMATFILEI_I();
		if (MATI1 == NULL)
			exit(0);
	}
	else
	{
		printf("enter the matrix of integers (first row non-zero) :\n");
		MATI1 = INPUTMATI();
	}
	printf("The matrix entered is\n\n");
	PRINTMATI(0, MATI1->R - 1, 0, MATI1->C - 1, MATI1);
	HERMITE1VERBOSE = 0;
	printf("VERBOSE? (Y/N)\n");
	answer = GetYN();
	if (answer)
		HERMITE1VERBOSE = 1;
	printf("enter the parameters m1 and n1 (normally 1 and 1) :");
	s = scanf("%u %u", &m1, &n1);
	Flush();
	MATI2 = LLLHERMITE1(MATI1, &MATI3, &rank, m1, n1);
	HERMITE1VERBOSE = 0;
	printf("rank = %u\n\n", rank);
	printf("The transformation matrix is\n");
	PRINTMATI(0, MATI2->R - 1, 0, MATI2->C - 1, MATI2);
	strcpy(buff, "lllhermitetrans.out");
	outfile = fopen(buff, "w");
	FPRINTMATI(outfile, 0, MATI2->R - 1, 0, MATI2->C - 1, MATI2);
	fclose(outfile);
	printf("The row echelon Form is\n");
	PRINTMATI(0, MATI3->R - 1, 0, MATI3->C - 1, MATI3);
	strcpy(buff, "lllhermitebas.out");
	outfile = fopen(buff, "w");
	FPRINTMATI(outfile, 0, MATI3->R - 1, 0, MATI3->C - 1, MATI3);
	fclose(outfile);
	FREEMATI(MATI1);
	FREEMATI(MATI2);
	FREEMATI(MATI3);
	return;
}

void CLEMENS(MPI *N)
/* constructs Clemens' matrix and then applies LLLHERMITE with alpha=3/4. */
{
	MPI *I;
	USI rank;
	USL i, j, n;
	MPMATI *MATI1, *MATI2, *MATI3;

	n = CONVERTI(N);
	MATI1 = BUILDMATI(n, n);
	for (j = 0; j < n; j++)
		MATI1->V[0][j] = ONEI();
	for (i = 1; i < n; i++)
		MATI1->V[i][0] = ZEROI();
	for (i = 1; i < n; i++)
	{
		I = CHANGEI(i);
		for (j = 1; j < n; j++)
			MATI1->V[i][j] = MPOWER_(j, I, N);
		FREEMPI(I);
	}
	PRINTMATI(0, MATI1->R - 1, 0, MATI1->C - 1, MATI1);
	GetReturn();
	MATI2 = LLLHERMITE1(MATI1, &MATI3, &rank, 3, 4);
	printf("The Hermite normal form = \n");
	PRINTMATI(0, MATI2->R - 1, 0, MATI2->C - 1, MATI2);
	printf("The unimodular transformation matrix = \n");
	PRINTMATI(0, MATI3->R - 1, 0, MATI3->C - 1, MATI3);
	FREEMATI(MATI1);
	FREEMATI(MATI2);
	FREEMATI(MATI3);
}

void CLEMENS1(MPI *N)
/* constructs Clemens' matrix and then applies Kannan-Bachem. */
{
	MPI *I;
	USI *alpha, nz;
	USL i, j, n;
	MPMATI *MATI1, *MATI2;

	n = CONVERTI(N);
	MATI1 = BUILDMATI(n, n);
	for (j = 0; j < n; j++)
		MATI1->V[0][j] = ONEI();
	for (i = 1; i < n; i++)
		MATI1->V[i][0] = ZEROI();
	for (i = 1; i < n; i++)
	{
		I = CHANGEI(i);
		for (j = 1; j < n; j++)
			MATI1->V[i][j] = MPOWER_(j, I, N);
		FREEMPI(I);
	}
	alpha = KB_ROW(MATI1, &nz);
	printf("rank = %u\n", nz);
	MATI2 = HERMITE1(MATI1, alpha, nz);
	printf("The Hermite normal form = \n");
	PRINTMATI(0, MATI2->R - 1, 0, MATI2->C - 1, MATI2);
	FREEMATI(MATI2);
	ffree((char *)alpha, (MATI1->C) * sizeof(USI));
	FREEMATI(MATI1);
}

void CLEMENS2(MPI *N)
/* constructs Clemens' matrix and then applies Kannan-Bachem. */
/* Also produces the transformation matrix. */
{
	MPI *I;
	USI *alpha, nz;
	USL i, j, n;
	MPMATI *MATI1, *MATI2, *MATI3, *MATI4;

	n = CONVERTI(N);
	MATI1 = BUILDMATI(n, n);
	for (j = 0; j < n; j++)
		MATI1->V[0][j] = ONEI();
	for (i = 1; i < n; i++)
		MATI1->V[i][0] = ZEROI();
	for (i = 1; i < n; i++)
	{
		I = CHANGEI(i);
		for (j = 1; j < n; j++)
			MATI1->V[i][j] = MPOWER_(j, I, N);
		FREEMPI(I);
	}
	alpha = KB_ROWP(MATI1, &MATI3, &nz);
	printf("rank = %u\n", nz);
	MATI2 = HERMITE1P(MATI1, MATI3, &MATI4, alpha, nz);
	FREEMATI(MATI3);
	printf("The Hermite normal form = \n");
	PRINTMATI(0, MATI2->R - 1, 0, MATI2->C - 1, MATI2);
	FREEMATI(MATI2);
	printf("The unimodular transformation matrix = \n");
	PRINTMATI(0, MATI4->R - 1, 0, MATI4->C - 1, MATI4);
	FREEMATI(MATI4);
	ffree((char *)alpha, (MATI1->C) * sizeof(USI));
	FREEMATI(MATI1);
}

void AXB()
{
	USI answer, m1, n1, rankA, nullA, i, j, k, r, s, flag, p, m, n;
	MPMATI *MATI1, *MATI2, *MATI3, *MATI4, *MATI5, *MATI6, *MATI7;
	MPMATI *Tmp, *Q, *M;
	char buff[20];
	FILE *outfile;
	MPI **XX, **X, *Temp;
	int e, ss;

	printf("Do you wish to use an existing augmented matrix from a file? (Y/N)\n");
	answer = GetYN();
	if (answer)
	{
		Tmp = FINPUTMATFILEI_I();
		if (Tmp == NULL)
			exit(0);
	}
	else
	{
		printf("enter the augmented matrix A|B of integers (first column non-zero) :\n");
		Tmp = INPUTMATI();
	}
	MATI1 = TRANSPOSEI(Tmp);
	FREEMATI(Tmp);
	printf("The matrix entered in transposed form is\n\n");
	PRINTMATI(0, MATI1->R - 1, 0, MATI1->C - 1, MATI1);
	MLLLVERBOSE = 0;
	printf("VERBOSE? (Y/N)\n");
	answer = GetYN();
	if (answer)
		MLLLVERBOSE = 1;
	MATI3 = IDENTITYI(MATI1->R);
	printf("enter the parameters m1 and n1 (normally 1 and 1) :");
	ss = scanf("%u %u", &m1, &n1);
	Flush();
	MATI2 = BASIS_REDUCTION(MATI1, &MATI3, 0, m1, n1);
	flag = 0;
	r = MATI3->R - 1;
	s = MATI3->C - 1;
	k = 0;
	while (k < r) {
		if (!EQZEROI(MATI3->V[k][s])) {
			flag = 1;
			break;
		}
		k++;
	}
	if (flag) {
		printf("No solution of AX=B\n");
		FREEMATI(MATI1);
		FREEMATI(MATI2);
		FREEMATI(MATI3);
		MLLLVERBOSE = 0;
		return;
	}
	printf("\n\n");
	if (MLLLVERBOSE) {
		printf("The interim transformation matrix is\n");
		PRINTMATI(0, MATI3->R - 1, 0, MATI3->C - 1, MATI3);
		GetReturn();
	}
	/*
	printf("The corresponding reduced basis for Rowspace(A^t) is\n");
	PRINTMATI(0,MATI2->R-1,0,MATI2->C-1,MATI2);
	GetReturn();
	strcpy(buff, "axbimbas.out");
	outfile = fopen(buff, "w");
	FPRINTMATI(outfile,0,MATI2->R-1,0,MATI2->C-1,MATI2);
	fclose(outfile);
	*/
	rankA = MATI2->R;
	nullA = (MATI1->R) - (rankA + 1);
	if (!nullA) {
		printf("Unique solution\n");
		MATI6 = BUILDMATI(1, MATI1->R - 1);
		for (j = 0; j < MATI1->R - 1; j++)
			MATI6->V[0][j] = MINUSI(MATI3->V[MATI3->R - 1][j]);
		FREEMATI(MATI1);
		FREEMATI(MATI2);
		FREEMATI(MATI3);
		printf("The solution X of XA^t=B^t is\n");
		PRINTMATI(0, MATI6->R - 1, 0, MATI6->C - 1, MATI6);
		strcpy(buff, "axbsol.out");
		outfile = fopen(buff, "w");
		FPRINTMATI(outfile, 0, MATI6->R - 1, 0, MATI6->C - 1, MATI6);
		fclose(outfile);
		FREEMATI(MATI6);
		MLLLVERBOSE = 0;
		return;
	}
	MATI4 = BUILDMATI(1 + nullA, MATI1->R - 1);
	for (i = 0; i <= nullA; i++) {
		for (j = 0; j < MATI1->R - 1; j++) {
			MATI4->V[i][j] = COPYI(MATI3->V[rankA + i][j]);
		}
	}
	GCDFLAG = 1;
	MATI5 = BASIS_REDUCTION0(MATI4, m1, n1);
	MLLLVERBOSE = 0;
	MATI6 = BUILDMATI(1, MATI1->R - 1);
	for (j = 0; j < MATI1->R - 1; j++)
	{
		(MATI5->V[MATI5->R - 1][j])->S = -((MATI5->V[MATI5->R - 1][j])->S);
		MATI6->V[0][j] = COPYI(MATI5->V[MATI5->R - 1][j]);
	}
	GCDFLAG = 0;
	MATI7 = DELETE_ROWI(MATI5->R, MATI5);
	printf("The short basis for nullspace(A^t) is\n");
	PRINTMATI(0, MATI7->R - 1, 0, MATI7->C - 1, MATI7);
	GetReturn();
	strcpy(buff, "axbbas.out");
	outfile = fopen(buff, "w");
	FPRINTMATI(outfile, 0, MATI7->R - 1, 0, MATI7->C - 1, MATI7);
	fclose(outfile);
	p = nullA;
	m = p + 1;

	printf("A short solution X of XA^t=B^t is \n");
	printf("b[%u] = ", m); PRINTMATI(0, MATI6->R - 1, 0, MATI6->C - 1, MATI6);
	/*
	strcpy(buff, "axbsol.out");
	outfile = fopen(buff, "w");
	FPRINTMATI(outfile,0,MATI6->R-1,0,MATI6->C-1,MATI6);
	fclose(outfile);
	*/
	printf("Do you want to get the shortest multipliers using Fincke_Pohst? (Y/N)\n");
	answer = GetYN();
	XX = (MPI **)mmalloc((USL)(p * sizeof(MPI *)));
	for (j = 0; j < p; j++)
		XX[j] = ZEROI();
	if (answer)
	{
		GCDVERBOSE = 1;
		Q = MATI5;
		n = Q->R;
		while (1)
		{
			M = SHORTESTT0(Q, &X);
			for (j = 0; j < p; j++)
			{
				Temp = XX[j];
				XX[j] = ADDI(XX[j], X[j]);
				FREEMPI(Temp);
				FREEMPI(X[j]);
			}
			ffree((char *)X, p * sizeof(MPI *));
			if (M == NULL)
				break;
			else
			{
				for (j = 0; j < Q->C; j++)
				{
					FREEMPI(Q->V[n - 1][j]);
					Q->V[n - 1][j] = COPYI(M->V[0][j]);
				}
				FREEMATI(MATI6);
				MATI6 = M;
			}
		}
		printf("found a shortest solution vector:\n");
	}
	else
		Q = MATI5;
	strcpy(buff, "axb.out");
	outfile = fopen(buff, "w");
	if (answer)
		fprintf(outfile, "A shortest solution vector is ");
	else
		fprintf(outfile, "A short multiplier is ");
	if (answer)
		printf("b[%u]", m);
	fprintf(outfile, "b[%u]", m);
	for (j = 0; j < p; j++)
	{
		e = XX[j]->S;
		if (e == -1)
		{
			printf("+");
			fprintf(outfile, "+");
			Temp = MINUSI(XX[j]);
			if (!EQONEI(Temp))
			{
				PRINTI(Temp);
				FPRINTI(outfile, Temp);
			}
			printf("b[%u]", j + 1);
			fprintf(outfile, "b[%u]", j + 1);
			FREEMPI(Temp);
		}
		if (e == 1)
		{
			printf("-");
			fprintf(outfile, "-");
			if (!EQONEI(XX[j]))
			{
				PRINTI(XX[j]);
				FPRINTI(outfile, XX[j]);
			}
			printf("b[%u]", j + 1);
			fprintf(outfile, "b[%u]", j + 1);
		}
	}
	printf("\n");
	fprintf(outfile, "\n=");
	for (j = 0; j < p; j++)
		FREEMPI(XX[j]);
	ffree((char *)XX, p * sizeof(MPI *));
	for (i = 0; i < MATI6->C; i++)
	{
		FPRINTI(outfile, MATI6->V[0][i]);
		fprintf(outfile, " ");

	}
	fprintf(outfile, "\n");
	fclose(outfile);
	GCDVERBOSE = 0;
	FREEMATI(MATI1);
	FREEMATI(MATI2);
	FREEMATI(MATI3);
	FREEMATI(MATI4);
	FREEMATI(MATI5);
	FREEMATI(MATI6);
	FREEMATI(MATI7);
	return;
}

void AXB1()
{
	USI answer, m1, n1, nullA, i, j, p, m, n, flag, flag1;
	USI rank, u, v, rankminus1;
	MPMATI *MATI1, *MATI2, *MATI3, *MATI5, *MATI6, *MATI7;
	MPMATI *Tmp, *Tmp1, *Q, *M, *COEFFICIENT_MATRIX;
	char buff[20];
	FILE *outfile;
	MPI **XX, **X, *Temp;
	int e, s;

	printf("Do you wish to use an existing augmented matrix from a file? (Y/N)\n");
	answer = GetYN();
	if (answer)
	{
		Tmp = FINPUTMATFILEI_I();
		if (Tmp == NULL)
			exit(0);
	}
	else
	{
		printf("enter the augmented matrix A|B of integers (first column non-zero) :\n");
		Tmp = INPUTMATI();
	}
	COEFFICIENT_MATRIX = BUILDMATI(Tmp->R, Tmp->C - 1);
	for (i = 0; i < Tmp->R; i++) {
		for (j = 0; j < Tmp->C - 1; j++) {
			COEFFICIENT_MATRIX->V[i][j] = COPYI(Tmp->V[i][j]);
		}
	}
	if (TEST_ZEROMATI(COEFFICIENT_MATRIX, Tmp->R, Tmp->C - 1)) {
		printf("Coefficient matrix is the zero matrix\n");
		FREEMATI(Tmp);
		FREEMATI(COEFFICIENT_MATRIX);
		return;
	}
	FREEMATI(COEFFICIENT_MATRIX);
	Tmp1 = TRANSPOSEI(Tmp);
	FREEMATI(Tmp);
	printf("The matrix entered in transposed form is\n\n");
	PRINTMATI(0, Tmp1->R - 1, 0, Tmp1->C - 1, Tmp1);
	u = Tmp1->R;
	v = Tmp1->C + 1;
	MATI1 = BUILDMATI(u, v);
	for (i = 0; i < u; i++) {
		for (j = 0; j < v - 1; j++)
			MATI1->V[i][j] = COPYI(Tmp1->V[i][j]);
		MATI1->V[i][v - 1] = ZEROI();
	}
	FREEMPI(MATI1->V[u - 1][v - 1]);
	MATI1->V[u - 1][v - 1] = ONEI();
	HERMITE1VERBOSE = 0;
	printf("VERBOSE? (Y/N)\n");
	answer = GetYN();
	if (answer)
		HERMITE1VERBOSE = 1;
	printf("enter the parameters m1 and n1 (normally 1 and 1) :");
	s = scanf("%u %u", &m1, &n1);
	Flush();
	MATI3 = LLLHERMITE1(MATI1, &MATI2, &rank, m1, n1);

	if (HERMITE1VERBOSE)
	{
		printf("MATI2=\n\n");
		PRINTMATI(0, MATI2->R - 1, 0, MATI2->C - 1, MATI2);
		GetReturn();
		printf("MATI3=\n\n");
		PRINTMATI(0, MATI3->R - 1, 0, MATI3->C - 1, MATI3);
		GetReturn();
	}
	strcpy(buff, "axbtrans.out");
	outfile = fopen(buff, "w");
	FPRINTMATI(outfile, 0, MATI3->R - 1, 0, MATI3->C - 1, MATI3);
	fclose(outfile);
	HERMITE1VERBOSE = 0;
	flag = 0;
	rankminus1 = rank - 1;
	for (j = 0; j<rankminus1; j++) {
		if ((MATI2->V[j][v - 1])->S != 0) {
			flag = 1;
			break;
		}
	}
	flag1 = 0;
	for (j = 0; j<Tmp1->C - 1; j++) {
		if ((MATI2->V[rankminus1][j])->S != 0) {
			flag1 = 1;
			break;
		}
	}

	FREEMATI(Tmp1);
	if (flag == 0 && flag1 == 0 && EQONEI(MATI2->V[rankminus1][v - 1])) {
		nullA = (MATI3->R) - rank;
		if (!nullA) {
			printf("Unique solution\n");
			MATI6 = BUILDMATI(1, MATI1->R - 1);
			for (j = 0; j < MATI1->R - 1; j++)
				MATI6->V[0][j] = MINUSI(MATI3->V[rank - 1][j]);
			FREEMATI(MATI1);
			FREEMATI(MATI2);
			FREEMATI(MATI3);
			printf("The solution X of XA^t=B^t is\n");
			PRINTMATI(0, MATI6->R - 1, 0, MATI6->C - 1, MATI6);
			strcpy(buff, "axbsol.out");
			outfile = fopen(buff, "w");
			FPRINTMATI(outfile, 0, MATI6->R - 1, 0, MATI6->C - 1, MATI6);
			fclose(outfile);
			FREEMATI(MATI6);
			return;
		}
		v = MATI3->C - 1;
		MATI5 = BUILDMATI(u + 1 - rank, v);
		for (i = 0; i < u - rank; i++) {
			for (j = 0; j < v; j++) {
				MATI5->V[i][j] = COPYI(MATI3->V[i + rank][j]);
			}
		}
		for (j = 0; j < v; j++) {
			MATI5->V[u - rank][j] = MINUSI(MATI3->V[rank - 1][j]);
		}
		printf("MATI5=B=\n"); PRINTMATI(0, MATI5->R - 1, 0, MATI5->C - 1, MATI5);
		GetReturn();

		MATI6 = BUILDMATI(1, MATI5->C);
		for (j = 0; j < MATI5->C; j++)
			MATI6->V[0][j] = COPYI(MATI5->V[MATI5->R - 1][j]);
		MATI7 = DELETE_ROWI(MATI5->R, MATI5);
		if (nullA>1) {
			printf("A short basis of %u vectors for nullspace(A^t):\n", nullA);
		}
		else {
			printf("A short basis of 1 vector for nullspace(A^t):\n");
		}
		PRINTMATI(0, MATI7->R - 1, 0, MATI7->C - 1, MATI7);
		GetReturn();
		strcpy(buff, "axbbas.out");
		outfile = fopen(buff, "w");
		FPRINTMATI(outfile, 0, MATI7->R - 1, 0, MATI7->C - 1, MATI7);
		fclose(outfile);
		p = nullA;
		m = p + 1;

		printf("\n\n");
		printf("A short solution X of XA^t=b^t: \n");
		printf("b[%u] = ", m); PRINTMATI(0, MATI6->R - 1, 0, MATI6->C - 1, MATI6);

		printf("\n\n");
		printf("Do you want to get a shortest solution using Fincke_Pohst? (Y/N)\n");
		answer = GetYN();
		XX = (MPI **)mmalloc((USL)(p * sizeof(MPI *)));
		for (j = 0; j < p; j++)
			XX[j] = ZEROI();
		if (answer)
		{
			GCDVERBOSE = 1;
			Q = MATI5;
			n = Q->R;
			while (1)
			{
				M = SHORTESTT0(Q, &X);
				for (j = 0; j < p; j++)
				{
					Temp = XX[j];
					XX[j] = ADDI(XX[j], X[j]);
					FREEMPI(Temp);
					FREEMPI(X[j]);
				}
				ffree((char *)X, p * sizeof(MPI *));
				if (M == NULL)
					break;
				else
				{
					for (j = 0; j < Q->C; j++)
					{
						FREEMPI(Q->V[n - 1][j]);
						Q->V[n - 1][j] = COPYI(M->V[0][j]);
					}
					FREEMATI(MATI6);
					MATI6 = M;
				}
			}
			printf("found a shortest solution vector:\n");
		}
		else
			Q = MATI5;
		strcpy(buff, "axb.out");
		outfile = fopen(buff, "w");
		if (answer)
			fprintf(outfile, "A shortest solution is ");
		else
			fprintf(outfile, "Again A short solution is ");
		if (answer)
			printf("b[%u]", m);
		fprintf(outfile, "b[%u]", m);
		for (j = 0; j < p; j++)
		{
			e = XX[j]->S;
			if (e == -1)
			{
				printf("+");
				fprintf(outfile, "+");
				Temp = MINUSI(XX[j]);
				if (!EQONEI(Temp))
				{
					PRINTI(Temp);
					FPRINTI(outfile, Temp);
				}
				printf("b[%u]", j + 1);
				fprintf(outfile, "b[%u]", j + 1);
				FREEMPI(Temp);
			}
			if (e == 1)
			{
				printf("-");
				fprintf(outfile, "-");
				if (!EQONEI(XX[j]))
				{
					PRINTI(XX[j]);
					FPRINTI(outfile, XX[j]);
				}
				printf("b[%u]", j + 1);
				fprintf(outfile, "b[%u]", j + 1);
			}
		}
		printf("\n");
		fprintf(outfile, "\n=");
		for (i = 0; i < MATI6->C; i++)
		{
			FPRINTI(outfile, MATI6->V[0][i]);
			fprintf(outfile, " ");

		}
		fprintf(outfile, "\n");
		fclose(outfile);
		if (answer)
		{
			printf("Do you want to get all the shortest solutions? (Y/N)\n");
			answer = GetYN();
			if (answer)
				SHORTEST(Q, XX, 5, 1);
		}
		for (j = 0; j < p; j++)
			FREEMPI(XX[j]);
		ffree((char *)XX, p * sizeof(MPI *));
		GCDVERBOSE = 0;
		FREEMATI(MATI1);
		FREEMATI(MATI2);
		FREEMATI(MATI3);
		FREEMATI(MATI5);
		FREEMATI(MATI6);
		FREEMATI(MATI7);
		return;
	}
	else {
		printf("No solution of AX=B\n");
		FREEMATI(MATI1);
		FREEMATI(MATI2);
		FREEMATI(MATI3);
		return;
	}
}

void TESTAXB()
{
	USI m1, n1, nullA, g, i, j, k, p, m, n, row, col;
	USI rank, u, v, power, trial_no, answer;
	MPMATI *MATI1, *MATI2, *MATI3, *MATI5, *MATI6, *MATI7;
	MPMATI *Tmp1, *Q, *M, *TMATI1;
	char buff[20];
	FILE *outfile;
	MPI **XX, **X, *Temp, *MULTIPLIER, *SEED, *BASE, *temp, *T1, *T2;
	MPR *ratio;
	int e, equality, s;

	strcpy(buff, "testaxb.out");
	if (access(buff, R_OK) == 0)
		unlink(buff);
	outfile = fopen(buff, "a");
	printf("Enter alpha=m/n,  1/4 < m/n <= 1:  (ie enter  m  and n)");
	s = scanf("%u%u", &m1, &n1);
	printf("Enter the power of 2: ");
	s = scanf("%u", &power);
	BASE = POWER_I((long)2, power);
	printf("Enter row and column size: ");
	s = scanf("%u%u", &row, &col);
	printf("Enter number of trials: ");
	s = scanf("%u", &trial_no);
	Flush();
	printf("Enter the (odd) seed (eg 1): ");
	SEED = INPUTI(&g);
	printf("Enter the multiplier (4k+1, n odd) (eg 96069): ");
	MULTIPLIER = INPUTI(&g);
	printf("small random X? (Y/N)");
	answer = GetYN();
	for (k = 0; k < trial_no; k++)
	{
		printf("k = %u\n", k);
		if (answer)
			TMATI1 = RANDOM_MATRIXA3(row, col, SEED, MULTIPLIER, BASE);
		else
			TMATI1 = RANDOM_MATRIXA(row, col, SEED, MULTIPLIER, BASE);
		/*
		printf("TMATI1 :\n");
		PRINTMATI(0,TMATI1->R-1,0,TMATI1->C-1,TMATI1);
		GetReturn();
		*/
		temp = SEED;
		SEED = RANDOMI(SEED, MULTIPLIER, BASE);
		FREEMPI(temp);
		Tmp1 = TRANSPOSEI(TMATI1);
		FREEMATI(TMATI1);
		u = Tmp1->R;
		v = Tmp1->C + 1;
		MATI1 = BUILDMATI(u, v);
		for (i = 0; i < u; i++) {
			for (j = 0; j < v - 1; j++)
				MATI1->V[i][j] = COPYI(Tmp1->V[i][j]);
			MATI1->V[i][v - 1] = ZEROI();
		}
		FREEMPI(MATI1->V[u - 1][v - 1]);
		MATI1->V[u - 1][v - 1] = ONEI();
		FREEMATI(Tmp1);
		MATI3 = LLLHERMITE1(MATI1, &MATI2, &rank, m1, n1);
		nullA = (MATI3->R) - rank;
		v = MATI3->C - 1;
		MATI5 = BUILDMATI(u + 1 - rank, v);
		for (i = 0; i < u - rank; i++) {
			for (j = 0; j < v; j++) {
				MATI5->V[i][j] = COPYI(MATI3->V[i + rank][j]);
			}
		}
		for (j = 0; j < v; j++)
			MATI5->V[u - rank][j] = MINUSI(MATI3->V[rank - 1][j]);
		MATI6 = BUILDMATI(1, MATI5->C);
		for (j = 0; j < MATI5->C; j++)
			MATI6->V[0][j] = COPYI(MATI5->V[MATI5->R - 1][j]);
		/*
		printf("Short solution :\n");
		PRINTMATI(0,MATI6->R-1,0,MATI6->C-1,MATI6);
		GetReturn();
		*/
		T1 = LENGTHSQRI(MATI6, 0);
		MATI7 = DELETE_ROWI(MATI5->R, MATI5);
		/*
		printf("The short basis of %u vectors for nullspace(A^t):\n", nullA);
		PRINTMATI(0,MATI7->R-1,0,MATI7->C-1,MATI7);
		GetReturn();
		*/
		p = nullA;
		m = p + 1;


		XX = (MPI **)mmalloc((USL)(p * sizeof(MPI *)));
		for (j = 0; j < p; j++)
			XX[j] = ZEROI();
		Q = MATI5;
		n = Q->R;
		while (1)
		{
			M = SHORTESTT0(Q, &X);
			for (j = 0; j < p; j++)
			{
				Temp = XX[j];
				XX[j] = ADDI(XX[j], X[j]);
				FREEMPI(Temp);
				FREEMPI(X[j]);
			}
			ffree((char *)X, p * sizeof(MPI *));
			if (M == NULL)
				break;
			else
			{
				for (j = 0; j < Q->C; j++)
				{
					FREEMPI(Q->V[n - 1][j]);
					Q->V[n - 1][j] = COPYI(M->V[0][j]);
				}
				FREEMATI(MATI6);
				MATI6 = M;
			}
		}
		T2 = LENGTHSQRI(MATI6, 0);
		/*
		printf("Shortest solution :\n");
		PRINTMATI(0,MATI6->R-1,0,MATI6->C-1,MATI6);
		GetReturn();
		*/
		equality = EQUALI(T1, T2);
		if (!equality)
		{
			ratio = RATIOI(T1, T2);
			FREEMPI(T1);
			FREEMPI(T2);
			fprintf(outfile, "%u:", k);
			FPRINTDR(outfile, 3, ratio);
			FREEMPR(ratio);
			fprintf(outfile, ":b[%u]", m);
			for (j = 0; j < p; j++)
			{
				e = XX[j]->S;
				if (e == -1)
				{
					fprintf(outfile, "+");
					Temp = MINUSI(XX[j]);
					if (!EQONEI(Temp))
						FPRINTI(outfile, Temp);
					fprintf(outfile, "b[%u]", j + 1);
					FREEMPI(Temp);
				}
				if (e == 1)
				{
					fprintf(outfile, "-");
					if (!EQONEI(XX[j]))
						FPRINTI(outfile, XX[j]);
					fprintf(outfile, "b[%u]", j + 1);
				}
			}
			/*
			fprintf(outfile, "\n=");
			for (i = 0; i < MATI6->C; i++)
			{
			FPRINTI(outfile, MATI6->V[0][i]);
			fprintf(outfile," ");

			}
			*/
			fprintf(outfile, "\n");
		}
		else
		{
			FREEMPI(T1);
			FREEMPI(T2);
		}
		for (j = 0; j < p; j++)
			FREEMPI(XX[j]);
		ffree((char *)XX, p * sizeof(MPI *));
		FREEMATI(MATI1);
		FREEMATI(MATI2);
		FREEMATI(MATI3);
		FREEMATI(MATI5);
		FREEMATI(MATI6);
		FREEMATI(MATI7);
	}
	FREEMPI(BASE);
	FREEMPI(MULTIPLIER);
	FREEMPI(SEED);
	fclose(outfile);

	return;
}
void CHANGELX(MPI *M)
{
	USL n;
	MPI *N;
	int s;
	long m;

	s = M->S;
	printf("M->D = %u\n", M->D);
	n = CONVERTI(M);
	printf("n = %lu\n", n);
	m = n;
	printf("m = %ld\n", m);
	if (s == -1)
		m = -m;
	N = CHANGEL(m);
	printf("CHANGEL(n)="); PRINTI(N); printf("\n");
	FREEMPI(N);
	return;
}

void LLLGCD0X()
/*
* This executes the LLLGCD algorithm of Havas, Majewski and Matthews.
* with a search for short multipliers using the method of Matthews.
* 17/12/98.
*/
{
	USI answer, m1, n1, m, n, i, j, p;
	int e, s;
	MPMATI *MATI1, *MATI2, *Q, *M, *MATI3;
	char buff[20];
	FILE *outfile;
	MPI *A, *T, **XX, **X, *Temp;

	HAVASFLAG = 1;
	printf("Do you wish to enter your sequence of numbers from an existing column matrix? (Y/N)\n");
	answer = GetYN();
	if (answer)
	{
		MATI1 = FINPUTMATFILEI_I();
		if (MATI1 == NULL)
			exit(0);
	}
	else
	{
		printf("WARNING: Make sure the first integer in the sequence is non-zero) :\n");
		MATI1 = INPUTMATI();
	}
	printf("The file entered is\n\n");
	PRINTMATI(0, MATI1->R - 1, 0, MATI1->C - 1, MATI1);
	GCDVERBOSE = 0;
	printf("VERBOSE? (Y/N)\n");
	answer = GetYN();
	if (answer)
		GCDVERBOSE = 1;
	printf("to enter alpha=m1/n1: select m1 and n1 (normally 1 and 1) :");
	s = scanf("%u %u", &m1, &n1);
	Flush();
	m = MATI1->R;
	MATI2 = LLLGCD0(MATI1, &A, m, m1, n1);
	FREEMATI(MATI1);
	MATI3 = BUILDMATI(1, m);
	for (i = 0; i < m; i++)
		MATI3->V[0][i] = COPYI(MATI2->V[m - 1][i]);
	printf("\n");
	printf("The multiplier matrix found by LLL is\n");
	PRINTMATI(0, MATI2->R - 1, 0, MATI2->C - 1, MATI2);
	printf("gcd = "); PRINTI(A); printf("\n");
	GetReturn();
	FREEMPI(A);
	printf("Do you want to get a shortest multiplier using Fincke_Pohst? (Y/N)\n");
	answer = GetYN();
	p = m - 1;
	XX = (MPI **)mmalloc((USL)(p * sizeof(MPI *)));
	for (j = 0; j < p; j++)
		XX[j] = ZEROI();
	if (answer)
	{
		GCDVERBOSE = 1;
		Q = MATI2;
		n = Q->R;
		while (1)
		{
			M = SHORTESTT0(Q, &X);
			for (j = 0; j < p; j++)
			{
				Temp = XX[j];
				XX[j] = ADDI(XX[j], X[j]);
				FREEMPI(Temp);
				FREEMPI(X[j]);
			}
			ffree((char *)X, p * sizeof(MPI *));
			if (M == NULL)
				break;
			else
			{
				for (j = 0; j < Q->C; j++)
				{
					FREEMPI(Q->V[n - 1][j]);
					Q->V[n - 1][j] = COPYI(M->V[0][j]);
				}
				FREEMATI(MATI3);
				MATI3 = M;
			}
		}
		printf("found a shortest multiplier vector:\n");
	}
	else
		Q = MATI2;
	strcpy(buff, "lllgcd0mat.out");
	outfile = fopen(buff, "w");
	FPRINTMATI(outfile, 0, Q->R - 1, 0, Q->C - 1, Q);
	fclose(outfile);
	strcpy(buff, "lllgcd0mult.out");
	outfile = fopen(buff, "a");
	if (answer)
	{
		fprintf(outfile, "A Shortest multiplier is ");
		printf("b[%u]", m);
		fprintf(outfile, "b[%u]", m);
	}
	for (j = 0; j < p; j++)
	{
		e = XX[j]->S;
		if (e == -1)
		{
			printf("+");
			fprintf(outfile, "+");
			Temp = MINUSI(XX[j]);
			if (!EQONEI(Temp))
			{
				PRINTI(Temp);
				FPRINTI(outfile, Temp);
			}
			printf("b[%u]", j + 1);
			fprintf(outfile, "b[%u]", j + 1);
			FREEMPI(Temp);
		}
		if (e == 1)
		{
			printf("-");
			fprintf(outfile, "-");
			if (!EQONEI(XX[j]))
			{
				PRINTI(XX[j]);
				FPRINTI(outfile, XX[j]);
			}
			printf("b[%u]", j + 1);
			fprintf(outfile, "b[%u]", j + 1);
		}
	}
	if (answer)
	{
		printf("=");
		fprintf(outfile, "=");
		for (i = 0; i < MATI3->C; i++)
		{
			PRINTI(MATI3->V[0][i]);
			FPRINTI(outfile, MATI3->V[0][i]);
			fprintf(outfile, " ");
			printf(" ");
		}
		printf(": ");
		fprintf(outfile, ": ");
		T = LENGTHSQRI(MATI3, 0);
		PRINTI(T);
		printf("\n");
		FPRINTI(outfile, T);
		fprintf(outfile, "\n");
		FREEMPI(T);
	}
	fclose(outfile);
	if (answer)
	{
		printf("Do you want to get all the shortest multipliers? (Y/N)\n");
		answer = GetYN();
		if (answer)
			SHORTEST(Q, XX, 3, 1);
	}
	for (j = 0; j < p; j++)
		FREEMPI(XX[j]);
	ffree((char *)XX, p * sizeof(MPI *));
	FREEMATI(MATI3);
	FREEMATI(Q);
	GCDVERBOSE = 0;
	HAVASFLAG = 0;
	return;
}

void SLVECTORX()
{
	USI answer, i, j, m, n;
	MPMATI *MATI1;
	MPMATR *M;
	MPR *C, *D, *tempR, *VALUE;

	printf("Do you wish to use an existing matrix from a file? (Y/N)\n");
	answer = GetYN();
	if (answer)
	{
		MATI1 = FINPUTMATFILEI_I();
		if (MATI1 == NULL)
			exit(0);
	}
	else
	{
		printf("enter the matrix of integers (first row non-zero) :\n");
		MATI1 = INPUTMATI();
	}
	printf("The matrix entered is\n\n");
	PRINTMATI(0, MATI1->R - 1, 0, MATI1->C - 1, MATI1);
	C = BUILDMPR();
	C->N = LENGTHSQRI(MATI1, 0);
	C->D = ONEI();
	m = MATI1->R;
	n = MATI1->C;
	M = BUILDMATR(m, n);
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
		{
			elt(M, i, j) = BUILDMPR();
			elt(M, i, j)->N = COPYI(MATI1->V[i][j]);
			elt(M, i, j)->D = ONEI();
		}
	}
	FREEMATI(MATI1);
	D = ZEROR();
	while ((C->N)->S)
	{
		tempR = C;
		C = SLVECTOR(M, C, &VALUE);
		if (VALUE != NULL)
		{
			FREEMPR(D);
			D = COPYR(VALUE);
			FREEMPR(VALUE);
		}
		FREEMPR(tempR);
	}
	FREEMPR(C);
	FINCKE_POHST(M, D, 2);
	FREEMPR(D);
	FREEMATR(M);
	return;
}

void GCD_CONJ()
/*
* This tests a conjecture on the location of the shortest multipliers
* in the extended GCD problem. See paper off lll.html.
* Here we first get the matrices L (the mu{{ij}) and B (the LLLGCD
* trnaformation matrix. Then we test all the shortest multipliers to
* see if the associated X[0],...,X[m-2] satisfy the conjecture that
* |X[k]+SIGMA[k]| < 1 holds.
*/
{
	USI answer, m1, n1, m, n, i, j, p, count;
	int r, t, K, s;
	MPMATI *MATI1, *MATI2, *Q, *M, *MATI3;
	char buff[20];
	FILE *outfile;
	MPI *A, **XX, **X, *Temp, ***COEFF, *tempI, *T;
	MPMATR *LMATRIX;
	MPR **SIGMA, *SUMR, *tR1, *tR2, *tR3, *tempR, *tempR1;

	HAVASFLAG = 1;
	printf("Do you wish to enter your sequence of numbers from an existing column matrix? (Y/N)\n");
	answer = GetYN();
	if (answer)
	{
		MATI1 = FINPUTMATFILEI_I();
		if (MATI1 == NULL)
			exit(0);
	}
	else
	{
		printf("WARNING: Make sure the first integer in the sequence is non-zero) :\n");
		MATI1 = INPUTMATI();
	}
	printf("The matrix entered is\n\n");
	PRINTMATI(0, MATI1->R - 1, 0, MATI1->C - 1, MATI1);
	printf("to enter alpha=m1/n1: select m1 and n1 (normally 1 and 1) :");
	s = scanf("%u %u", &m1, &n1);
	Flush();
	m = MATI1->R;
	MATI2 = LLLGCDL(MATI1, &A, m, m1, n1, &LMATRIX);
	printf("The multiplier matrix found by LLL is\n");
	PRINTMATI(0, MATI2->R - 1, 0, MATI2->C - 1, MATI2);
	printf("gcd = "); PRINTI(A); printf("\n");
	FREEMPI(A);
	printf("The multiplier vector found by LLL is b[%u]=", m);
	for (i = 0; i < MATI2->R; i++)
	{
		PRINTI(MATI2->V[MATI2->R - 1][i]);
		printf(" ");
	}
	printf("\n");
	MATI3 = BUILDMATI(1, m);
	for (i = 0; i < m; i++)
		MATI3->V[0][i] = COPYI(MATI2->V[m - 1][i]);
	T = LENGTHSQRI(MATI3, 0);
	printf("The Euclidean norm squared = ");
	PRINTI(T);
	FREEMPI(T);
	printf("\n");
	p = m - 1;
	XX = (MPI **)mmalloc((USL)(p * sizeof(MPI *)));
	for (j = 0; j < p; j++)
		XX[j] = ZEROI();

	Q = MATI2;
	n = Q->R;
	while (1)
	{
		M = SHORTESTT0(Q, &X);
		for (j = 0; j < p; j++)
		{
			Temp = XX[j];
			XX[j] = ADDI(XX[j], X[j]);
			FREEMPI(Temp);
			FREEMPI(X[j]);
		}
		ffree((char *)X, p * sizeof(MPI *));
		if (M == NULL)
			break;
		else
		{
			for (j = 0; j < Q->C; j++)
			{
				FREEMPI(Q->V[n - 1][j]);
				Q->V[n - 1][j] = COPYI(M->V[0][j]);
			}
			FREEMATI(MATI3);
			MATI3 = M;
		}
	}
	printf("found a shortest multiplier lengthsquared:\n");
	T = LENGTHSQRI(MATI3, 0);
	printf(" = ");
	PRINTI(T);
	printf("\n");
	FREEMPI(T);

	strcpy(buff, "gcdconj.out");
	if (access(buff, R_OK) == 0)
		unlink(buff);
	outfile = fopen(buff, "a");


	COEFF = SHORTESTX(Q, XX, &count);

	for (j = 0; j < count; j++)
	{
		printf("COEFF[%u] = ", j + 1);
		fprintf(outfile, "COEFF[%u] = ", j + 1);
		for (i = 0; i < p; i++)
		{
			PRINTI(COEFF[j][i]);
			FPRINTI(outfile, COEFF[j][i]);
			printf(" ");
			fprintf(outfile, " ");
		}
		printf("\n");
		fprintf(outfile, "\n");
	}
	for (j = 0; j < p; j++)
		FREEMPI(XX[j]);
	ffree((char *)XX, p * sizeof(MPI *));
	SIGMA = (MPR **)mmalloc(p * sizeof(MPR *));
	for (j = 0; j < count; j++)
	{
		for (K = p - 1; K >= 0; K--)
		{
			SUMR = COPYR(elt(LMATRIX, m - 1, K));
			for (r = p - 1; r > K; r--)
			{
				tR2 = COPYR(elt(LMATRIX, r, K));
				tempR = BUILDMPR();
				tempR->N = COPYI(COEFF[j][r]);
				tempR->D = ONEI();
				tR3 = MULTR(tR2, tempR);
				FREEMPR(tempR);
				FREEMPR(tR2);
				tR1 = SUMR;
				SUMR = ADDR(SUMR, tR3);
				FREEMPR(tR1);
				FREEMPR(tR3);
			}
			SIGMA[K] = SUMR;
			tempR = BUILDMPR();
			tempR->N = COPYI(COEFF[j][K]);
			tempR->D = ONEI();
			tempR1 = ADDR(tempR, SIGMA[K]);
			tempI = ABSI(tempR1->N);
			t = RSV(tempI, tempR1->D);
			FREEMPR(tempR);
			FREEMPR(tempR1);
			FREEMPI(tempI);
			if (t >= 0)
			{
				fprintf(stderr, "conjecture false for SIGMA[%u]: K = %d\n", j + 1, K + 1);
				fprintf(outfile, "conjecture false for SIGMA[%u]: K = %d\n", j + 1, K + 1);
			}
		}
		printf("COEFF[%u], ", j + 1);
		fprintf(outfile, "COEFF[%u], ", j + 1);
		printf("SIGMA[%u]:\n", j + 1);
		fprintf(outfile, "SIGMA[%u]:\n", j + 1);
		for (K = 0; K < p; K++)
		{
			printf("COEFF[%u][%d]=", j + 1, K + 1);
			fprintf(outfile, "COEFF[%u][%d]=", j + 1, K + 1);
			PRINTI(COEFF[j][K]);
			FPRINTI(outfile, COEFF[j][K]);
			printf(", ");
			fprintf(outfile, ", ");
			printf("SIGMA[%u][%d]=", j + 1, K + 1);
			fprintf(outfile, "SIGMA[%u][%d]=", j + 1, K + 1);
			PRINTR(SIGMA[K]);
			FPRINTR(outfile, SIGMA[K]);
			printf("\n");
			fprintf(outfile, "\n");
			FREEMPR(SIGMA[K]);
		}
		printf("\n");
		fprintf(outfile, "\n");
	}
	fclose(outfile);
	ffree((char *)SIGMA, p * sizeof(MPR *));

	for (j = 0; j < count; j++)
	{
		for (i = 0; i < p; i++)
			FREEMPI(COEFF[j][i]);
		ffree((char *)COEFF[j], p * sizeof(MPI *));
	}
	ffree((char *)COEFF, GCD_MAX * sizeof(MPI **));
	FREEMATI(MATI3);
	FREEMATI(MATI1);
	FREEMATI(Q);
	FREEMATR(LMATRIX);
	HAVASFLAG = 0;
	return;
}

MPI *LEASTQNRX(MPI *P)
/*
* Returns NP, the least quadratic non-residue (mod P)
* if NP < R0, else returns 0.
*/
{
	MPI *NP, *T;
	int t;

	if (P->S <= 0 || (P->D == 0 && P->V[0] <= 2))
	{
		printf("P <= 2\n");
		return NULL;
	}
	T = LUCAS(P);
	t = T->S;
	FREEMPI(T);
	if (!t)
	{
		PRINTI(P);
		printf(" is not a prime\n");
		return NULL;
	}
	NP = LEASTQNR(P);
	if (NP->S == 0)
	{
		fprintf(stderr, "np > %lu\n", (USL)R0);
		printf("search failed\n");
		return NULL;
	}
	else
		return (NP);
}
USI TEST_ZEROMATI(MPMATI *A, USI m, USI n) {
	USI i, j;
	for (i = 0; i<m; i++) {
		for (j = 0; j<n; j++) {
			if ((A->V[i][j])->S != 0) {
				return(0);
			}
		}
	}
	return(1);
}

USI NAGELL(MPI *D, MPI *N, USI flag) {
	MPI *NN, *X1, *Y1, *ZERO, *E, *F, *T, *BOUND, *X, *T1, *TEMP, *TEMP1, *MINUSX, *S;
	MPI *TEMPX1, *TEMPY1, *U1, *U2, *V1, *V2, *U, *V, *Y;
	MPR *R, *BB;
	MPIA UU, VV;
	USL bound, y, start, flag1;
	USI t, count = 0, i;

	NN = ABSI(N);
	ZERO = ZEROI();
	S = PEL(D, ZERO, &X1, &Y1);
	if (S->S < 0) {
		TEMPX1 = X1;
		TEMPY1 = Y1;
		TEMP = MULTI(X1, X1);
		TEMP1 = MULTI3(D, Y1, Y1);
		X1 = ADD0I(TEMP, TEMP1);
		FREEMPI(TEMP);
		FREEMPI(TEMP1);
		TEMP = MULTI(TEMPX1, TEMPY1);
		Y1 = MULT_I(TEMP, 2);
		FREEMPI(TEMP);
		FREEMPI(TEMPX1);
		FREEMPI(TEMPY1);
	}
	if (flag) {
		printf("X1="); PRINTI(X1);
		printf(",");
		printf("Y1="); PRINTI(Y1);
		printf("\n");
	}
	FREEMPI(S);
	FREEMPI(ZERO);
	E = MULTI3(Y1, Y1, NN);
	if (N->S > 0) {
		TEMP = ADD0_I(X1, 1);
		F = MULT_I(TEMP, 2);
	}
	else {
		TEMP = SUB0_I(X1, 1);
		F = MULT_I(TEMP, 2);
	}
	FREEMPI(TEMP);

	T = INT0(E, F);
	BOUND = BIG_MTHROOT(T, 2);
	if (flag) {
		R = BUILDMPR();
		R->N = COPYI(E);
		R->D = COPYI(F);
		printf("Nagell upper bound:");
		BB = MTHROOTR(R, 2, 2);
		FREEMPR(R);
		PRINTDR(2, BB);
		FREEMPR(BB);
		printf("\n");
	}
	FREEMPI(E);
	FREEMPI(T);
	if (BOUND->D > 1) {
		printf("BOUND is too large\n");
		FREEMPI(BOUND);
		FREEMPI(NN);
		FREEMPI(F);
		FREEMPI(X1);
		FREEMPI(Y1);
		return(0);
	}
	bound = CONVERTI(BOUND);
	FREEMPI(BOUND);
	if (N->S > 0) {
		start = 0;
	}
	else {
		start = 1;
	}
	flag1 = 0;
	UU = BUILDMPIA();
	VV = BUILDMPIA();
	if (flag) {
		printf("Fundamental solutions:\n");
	}
	for (y = start; y <= bound; y++) {
		Y = CHANGE(y);
		TEMP = MULTI(D, Y);
		TEMP1 = MULTI(TEMP, Y);
		FREEMPI(TEMP);
		T1 = ADDI(TEMP1, N);
		FREEMPI(TEMP1);
		if (T1->S < 0) {
			FREEMPI(T1);
			FREEMPI(Y);
			continue;
		}
		else {
			X = BIG_MTHROOT(T1, 2);
			TEMP = MULTI(X, X);
			t = EQUALI(TEMP, T1);
			FREEMPI(TEMP);
			if (t) {
				flag1 = 1;
				if (flag) {
					printf("("); PRINTI(X); printf(",");
					printf("%lu)", y);
				}
				ADD_TO_MPIA(UU, X, count);
				ADD_TO_MPIA(VV, Y, count);
				count++;
				if (y == 0) {
					if (flag) {
						printf(" defines an ambiguous class");
					}
					flag1 = 0;
				}
				TEMP = MULTI(F, Y);
				TEMP1 = MULTI(TEMP, Y);
				FREEMPI(TEMP);
				TEMP = MULTI3(Y1, Y1, NN);
				if (y > 0 && (EQUALI(TEMP1, TEMP) || EQZEROI(X))) {
					if (flag) {
						printf(" defines an ambiguous class");
					}
					flag1 = 0;
				}
				if (flag) {
					printf("\n");
				}
				FREEMPI(TEMP1);
				FREEMPI(TEMP);
				if (flag1) {
					MINUSX = MINUSI(X);
					if (flag) {
						printf("("); PRINTI(MINUSX); printf(",");
						printf("%lu)\n", y);
					}
					flag1 = 0;
					FREEMPI(MINUSX);
					U1 = MULTI(X, X1);
					U2 = MULTI3(Y, Y1, D);
					U = SUBI(U1, U2);
					if (N->S < 0) {
						U->S = 1;
					}
					ADD_TO_MPIA(UU, U, count);
					FREEMPI(U1);
					FREEMPI(U2);
					FREEMPI(U);
					V1 = MULTI(X, Y1);
					V2 = MULTI(Y, X1);
					V = SUBI(V1, V2);
					if (N->S < 0) {
						V->S = 1;
					}
					ADD_TO_MPIA(VV, V, count);
					FREEMPI(V1);
					FREEMPI(V2);
					FREEMPI(V);
					count++;
				}
			}
			FREEMPI(X);
			FREEMPI(T1);
			FREEMPI(Y);
		}
	}
	if (count) {
		if (flag) {
			if (count>1) {
				printf("there are %u solution classes\n", count);
			}
			else {
				printf("there is %u solution classe\n", count);
			}
			printf("The non-negative solutions have a basis:\n");
			for (i = 0; i<count; i++) {
				printf("("); PRINTI((UU->A)[i]); printf(","); PRINTI((VV->A)[i]); printf(")\n");
			}
		}
	}
	else {
		printf("there are no solutions\n");
	}
	FREEMPIA(UU);
	FREEMPIA(VV);
	FREEMPI(NN);
	FREEMPI(F);
	FREEMPI(X1);
	FREEMPI(Y1);
	return(count);
}

MPI *NAGELLX(MPI *D, MPI *N)
{
	MPI *G, *X;
	unsigned long f;
	USI count;

	if (D->S <= 0)
	{
		printf("D <= 0\n");
		return NULL;
	}
	if (EQONEI(D))
	{
		printf("D = 1\n");
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
	if (N->S == 0)
	{
		printf("N = 0\n");
		return NULL;
	}
	if (EQONEI(N) || EQMINUSONEI(N)) {
		printf("|N| = 1\n");
		return NULL;
	}
	count = NAGELL(D, N, 1);
	return(ONEI());
}

void FRATTINI(MPI *D, MPI *N) {
	MPI *NN, *X1, *Y1, *ZERO, *E, *XX, *BOUND, *X, *TEMP, *TEMP1, *S;
	MPI *TEMPX1, *TEMPY1, *Y;
	MPR *R, *BB;
	USL bound, y;
	USI count = 0;

	NN = ABSI(N);
	ZERO = ZEROI();
	S = PEL(D, ZERO, &X1, &Y1);
	FREEMPI(ZERO);
	if (S->S < 0) {
		TEMPX1 = X1;
		TEMPY1 = Y1;
		TEMP = MULTI(X1, X1);
		TEMP1 = MULTI3(D, Y1, Y1);
		X1 = ADD0I(TEMP, TEMP1);
		FREEMPI(TEMP);
		FREEMPI(TEMP1);
		TEMP = MULTI(TEMPX1, TEMPY1);
		Y1 = MULT_I(TEMP, 2);
		FREEMPI(TEMP);
		FREEMPI(TEMPX1);
		FREEMPI(TEMPY1);
	}
	FREEMPI(S);
	printf("X1="); PRINTI(X1);
	printf(",");
	printf("Y1="); PRINTI(Y1);
	printf("\n");
	E = MULTI3(Y1, Y1, NN);
	FREEMPI(NN);
	BOUND = BIG_MTHROOT(E, 2);
	R = BUILDMPR();
	R->N = COPYI(E);
	R->D = ONEI();
	BB = MTHROOTR(R, 2, 2);
	FREEMPR(R);
	PRINTDR(2, BB);
	FREEMPR(BB);
	printf("\n");
	if (BOUND->D >= 1) {
		printf("BOUND is too large\n");
		FREEMPI(BOUND);
		FREEMPI(X1);
		FREEMPI(Y1);
		FREEMPI(E);
		return;
	}
	bound = CONVERTI(BOUND);
	TEMP1 = MULTI(BOUND, BOUND);
	FREEMPI(BOUND);
	if (EQUALI(TEMP1, E)) {
		bound--;
	}
	FREEMPI(TEMP1);
	FREEMPI(E);
	for (y = 0; y <= bound; y++) {
		Y = CHANGE(y);
		TEMP = MULTI3(D, Y, Y);
		TEMP1 = ADDI(TEMP, N);
		FREEMPI(Y);
		FREEMPI(TEMP);
		if (TEMP1->S < 0) {
			FREEMPI(TEMP1);
			continue;
		}
		else {
			X = BIG_MTHROOT(TEMP1, 2);
			XX = MULTI(X, X);
			if (EQUALI(XX, TEMP1)) {
				printf("("); PRINTI(X); printf(","); printf("%lu)", y);
				printf(" is a basic positive solution\n");
				count++;
			}
			FREEMPI(TEMP1);
			FREEMPI(X);
			FREEMPI(XX);
		}
	}
	if (count) {
		printf("there are %u solution classes\n", count);
	}
	else {
		printf("there are no solutions\n");
	}
	FREEMPI(X1);
	FREEMPI(Y1);
	return;
}

MPI *FRATTINIX(MPI *D, MPI *N)
{
	MPI *G, *X;
	unsigned long f;

	if (D->S <= 0)
	{
		printf("D <= 0\n");
		return NULL;
	}
	if (EQONEI(D))
	{
		printf("D = 1\n");
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
	if (N->S == 0)
	{
		printf("N = 0\n");
		return NULL;
	}
	if (EQONEI(N) || EQMINUSONEI(N)) {
		printf("|N| = 1\n");
		return NULL;
	}
	FRATTINI(D, N);
	return(ONEI());
}

void KASHIHARA(MPI *M, MPI *N) {
	MPI *TEMP, *Z, *TWO, *D, *NN;
	USL z, m, n;
	USI count;
	char buff[20];
	FILE *outfile;

	strcpy(buff, "kashihara.out");
	if (access(buff, R_OK) == 0)
		unlink(buff);
	outfile = fopen(buff, "a");
	TWO = TWOI();
	m = CONVERTI(M);
	n = CONVERTI(N);
	for (z = m; z <= n; z++) {
		Z = CHANGE(z);
		TEMP = MULTI(Z, Z);
		D = SUB0_I(TEMP, 1);
		NN = SUBI(TWO, TEMP);
		FREEMPI(TEMP);
		FREEMPI(Z);
		count = NAGELL(D, NN, 0);
		FREEMPI(D);
		FREEMPI(NN);
		if (z % 100000 == 0) {
			printf("z=%lu\n", z);
		}
		if (count >= 6) {
			printf("count = %u for z=%lu\n", count, z);
			fprintf(outfile, "count = %u for z=%lu\n", count, z);
		}
	}
	fclose(outfile);
	FREEMPI(TWO);
}

MPI *KASHIHARAX(MPI *M, MPI *N)
{
	MPI *TEMP;
	if (M->S <= 0 || N->S <= 0)
	{
		printf("M <= 0\n");
		return NULL;
	}
	if (M->S == 1 && M->D == 0 && M->V[0] == 1)
	{
		printf("M = 1\n");
		return NULL;
	}
	if (N->S == 1 && N->D == 0 && N->V[0] == 1)
	{
		printf("N = 1\n");
		return NULL;
	}
	if (RSV(M, N) == 1) {
		TEMP = N;
		N = M;
		M = TEMP;
	}
	KASHIHARA(M, N);
	return(ONEI());
}

void  PEL4(MPI *D, MPI *Eptr, MPI **Xptr, MPI **Yptr) {
	/* This finds the least positive solution (X1,Y1)=(*Xptr,*Yptr) of x^2-Dy^2=4.
	* where D>0 and is nonsquare.
	* We use the fact that a positive solution with gcd(x,y)=1 must be a convergent of
	* the continued fraction expansion of sqrt(d). If there is no such convergent, we use th
	* familiar midpoint method of finding the least solution of Pell's equations X^2-dY^2=1
	* and then (X1,Y1)=(2X,2Y).
	*/

	MPI *X, *G, *TMP, *TWO, *THREE, *FOUR, *SIX, *SEVEN, *EIGHT, *ELEVEN;
	MPI *TWELVE, *THIRTEEN, *FOURTEEN, *FIFTEEN, *P, *Q, *L, *M, *N, *K;
	MPI *Y, *U, *V, *OLDL, *OLDM, *MINUSP, *E, *F, *TMP1;
	MPI *TMP2, *R, *S, *FIVE, *TEN, *SEVENTEEN;
	USI t, h;
	int tt;
	long int flag = 0;

	X = BIG_MTHROOT(D, 2);
	G = MULTI(X, X);
	TMP = ADD0_I(G, (USL)1);
	t = EQUALI(D, TMP);
	FREEMPI(TMP);
	TWO = CHANGEI(2);
	FOUR = CHANGEI(4);
	if (EQUALI(D, TWO)) {
		*Xptr = CHANGEI(6);
		*Yptr = CHANGEI(4);
		flag = 1;
	}
	THREE = CHANGEI(3);
	if (EQUALI(D, THREE)) {
		*Xptr = CHANGEI(4);
		*Yptr = CHANGEI(2);
		flag = 1;
	}
	FREEMPI(THREE);
	FIVE = CHANGEI(5);
	if (EQUALI(D, FIVE)) {
		*Xptr = CHANGEI(3);
		*Yptr = CHANGEI(1);
		flag = 1;
	}
	FREEMPI(FIVE);
	SIX = CHANGEI(6);
	if (EQUALI(D, SIX)) {
		*Xptr = CHANGEI(10);
		*Yptr = CHANGEI(4);
		flag = 1;
	}
	FREEMPI(SIX);
	SEVEN = CHANGEI(7);
	if (EQUALI(D, SEVEN)) {
		*Xptr = CHANGEI(16);
		*Yptr = CHANGEI(6);
		flag = 1;
	}
	FREEMPI(SEVEN);
	EIGHT = CHANGEI(8);
	if (EQUALI(D, EIGHT)) {
		*Xptr = CHANGEI(6);
		*Yptr = CHANGEI(2);
		flag = 1;
	}
	FREEMPI(EIGHT);
	TEN = CHANGEI(10);
	if (EQUALI(D, TEN)) {
		*Xptr = CHANGEI(38);
		*Yptr = CHANGEI(12);
		flag = 1;
	}
	FREEMPI(TEN);
	ELEVEN = CHANGEI(11);
	if (EQUALI(D, ELEVEN)) {
		*Xptr = CHANGEI(20);
		*Yptr = CHANGEI(2);
		flag = 1;
	}
	FREEMPI(ELEVEN);
	TWELVE = CHANGEI(12);
	if (EQUALI(D, TWELVE)) {
		*Xptr = CHANGEI(4);
		*Yptr = ONEI();
		flag = 1;
	}
	FREEMPI(TWELVE);
	THIRTEEN = CHANGEI(13);
	if (EQUALI(D, THIRTEEN)) {
		*Xptr = CHANGEI(11);
		*Yptr = CHANGEI(3);
		flag = 1;
	}
	FREEMPI(THIRTEEN);
	FOURTEEN = CHANGEI(14);
	if (EQUALI(D, FOURTEEN)) {
		*Xptr = CHANGEI(30);
		*Yptr = CHANGEI(8);
		flag = 1;
	}
	FREEMPI(FOURTEEN);
	FIFTEEN = CHANGEI(15);
	if (EQUALI(D, FIFTEEN)) {
		*Xptr = CHANGEI(8);
		*Yptr = CHANGEI(2);
		flag = 1;
	}
	FREEMPI(FIFTEEN);
	SEVENTEEN = CHANGEI(17);
	if (t && RSV(D, SEVENTEEN)>0) /* period 1, exceptional case */
	{
		*Xptr = MULTABC(TWO, FOUR, G);
		*Yptr = MULTI(FOUR, X);
		flag = 1;
	}
	FREEMPI(SEVENTEEN);
	if (flag) {
		if (Eptr->S) {
			printf("(X1,Y1) = (");
			PRINTI(*Xptr);
			printf(",");
			PRINTI(*Yptr);
			printf(")\n");
		}
		FREEMPI(X);
		FREEMPI(G);
		FREEMPI(TWO);
		FREEMPI(FOUR);
		return;
	}

	P = ZEROI();
	Q = ONEI();
	L = ZEROI();
	K = ONEI();
	M = ONEI();
	N = ZEROI();
	for (h = 0; 1; h++) {
		TMP = ADD0I(X, P);
		Y = INT0(TMP, Q);
		FREEMPI(TMP);
		U = MULTABC(L, K, Y);
		V = MULTABC(M, N, Y);
		if (Eptr->S) {
			printf("A[%u]/B[%u] = ", h, h);
			PRINTI(U);
			printf("/");
			PRINTI(V);
			printf("\n");
		}
		OLDL = L;
		OLDM = M;
		L = K;
		M = N;
		K = U;
		N = V;
		F = COPYI(P);
		MINUSP = MINUSI(P);
		TMP = P;
		P = MULTABC(MINUSP, Y, Q);
		FREEMPI(Y);
		FREEMPI(TMP);
		FREEMPI(MINUSP);
		E = COPYI(Q);
		TMP = MULTI(P, P);
		TMP1 = SUBI(D, TMP);
		FREEMPI(TMP);
		TMP = Q;
		Q = INT(TMP1, Q);
		FREEMPI(TMP);
		FREEMPI(TMP1);
		if (h % 2) {
			tt = 1;
		}
		else {
			tt = -1;
		}
		if (EQUALI(Q, FOUR) && tt == 1) {
			printf("Q=4 & tt=1\n");
			*Xptr = COPYI(K);
			*Yptr = COPYI(N);
			flag = 4;
			break;
		}
		if (EQUALI(Q, FOUR) && tt == -1) {
			printf("Q=4 & tt = -1\n");
			TMP = MULTI3(D, N, N);
			TMP1 = MULTI(K, K);
			TMP2 = ADD0I(TMP1, TMP);
			*Xptr = INT0_(TMP2, 2);
			*Yptr = MULTI(K, N);;
			FREEMPI(TMP);
			FREEMPI(TMP1);
			FREEMPI(TMP2);
			flag = -4;
			break;
		}
		if (EQUALI(P, F)) { /* P[h] = P[h+1], even period 2h */
			if (Eptr->S) {
				printf("P[%u]=P[%u]\n", h, h + 1);
			}
			TMP = ADD0I(U, OLDL);
			TMP1 = MULTI(M, TMP);
			FREEMPI(TMP);
			if (tt == 1) {
				TMP = SUB0_I(TMP1, 1);
			}
			else {
				TMP = ADD0_I(TMP1, 1);
			}
			*Xptr = MULT_I(TMP, 2);
			FREEMPI(TMP);
			FREEMPI(TMP1);
			TMP = ADD0I(V, OLDM);
			TMP1 = MULTI(M, TMP);
			*Yptr = MULT_I(TMP1, 2);
			FREEMPI(TMP);
			FREEMPI(TMP1);
			flag = 1;
			break;
		}
		if (EQUALI(Q, E)) { /* Q[h] = Q[h+1], odd period 2h + 1 */
			if (Eptr->S) {
				printf("Q[%u]=Q[%u]\n", h, h + 1);
			}
			R = MULTAB_PLUS_CDI(K, N, L, M);
			S = MULTAB_PLUS_CDI(N, N, M, M);
			TMP = MULTI3(D, S, S);
			TMP1 = MULTI(R, R);
			TMP2 = ADD0I(TMP, TMP1);
			FREEMPI(TMP);
			FREEMPI(TMP1);
			*Xptr = MULT_I(TMP2, 2);;
			FREEMPI(TMP2);
			*Yptr = MULTI3(FOUR, R, S);
			FREEMPI(R);
			FREEMPI(S);
			flag = -1;
			break;
		}
		FREEMPI(OLDL);
		FREEMPI(OLDM);
		FREEMPI(E);
		FREEMPI(F);
	}
	FREEMPI(X);
	FREEMPI(G);
	FREEMPI(TWO);
	FREEMPI(FOUR);
	FREEMPI(OLDL);
	FREEMPI(OLDM);
	FREEMPI(L);
	FREEMPI(M);
	FREEMPI(P);
	FREEMPI(Q);
	FREEMPI(U);
	FREEMPI(V);
	FREEMPI(F);
	FREEMPI(E);
	if (Eptr->S) {
		printf("(X1,Y1) = (");
		PRINTI(*Xptr);
		printf(",");
		PRINTI(*Yptr);
		printf(")\n");
	}
	return;
}

MPI *PEL4X(MPI *D, MPI *E, MPI **Xptr, MPI **Yptr)
/*
* This is a program for finding the least solution of Pell's equation
* x*x - D*y*y = 4.
* The partial quotients are printed iff E->S != 0.
* Returns NULL on failure.  Wrapper handles this.
*/
{
	MPI *G, *X;
	unsigned int t;

	if (D->S <= 0)
	{
		printf("D <= 0\n");
		*Xptr = ZEROI();
		*Yptr = ZEROI();
		return(NULL);
	}
	X = BIG_MTHROOT(D, 2);
	G = MULTI(X, X);
	t = EQUALI(D, G);
	FREEMPI(X);
	FREEMPI(G);
	if (t)
	{
		printf("D is a perfect square\n");
		*Xptr = ZEROI();
		*Yptr = ZEROI();
		return(NULL);
	}
	PEL4(D, E, Xptr, Yptr);
	return(ONEI());
}

MPI *STOLT(MPI *A, MPI *B, MPI *C, MPI *N) {
	/* This program uses the bounds provided by Bengt Stolt, to find the fundamental solutions
	* of the diophantine equation Au^2+Buv+Cv^2=N, where N is non-zero, A>0 and D=B^2-4AC>0
	* is nonsquare.  We assume that gcd(A,B,C)=1.
	* See a paper by author, John Robertson and Anitha Srinivasan, where with small modification,
	* the Stolt bounds characterize the fundamental solutions.
	*/
	MPI *D, *TMP, *TMP1, *TWOA, *X1, *Y1, *AN, *Z, *ZZ, *MINUSAN, *G, *UB, *LB, *ROOTZZ;
	MPI *T, *DD, *T1, *SOLN, *TT, *TEMP, *TMP2, *X, *ZERO;
	USI t, soln, t1;
	USL v;
	long int lb, ub;
	MPR *RATIO;

	soln = 0;
	TMP = MULTI(A, C);
	TMP1 = MULT_I(TMP, 4);
	FREEMPI(TMP);
	TMP = MULTI(B, B);
	D = SUBI(TMP, TMP1);
	FREEMPI(TMP);
	FREEMPI(TMP1);

	printf("D = "); PRINTI(D); printf("\n");

	TWOA = MULT_I(A, 2);
	ZERO = ZEROI();
	PEL4(D, ZERO, &X1, &Y1);
	FREEMPI(ZERO);
	printf("(X1,Y1) = ("); PRINTI(X1); printf(", "); PRINTI(Y1); printf(")\n");

	AN = MULTI(A, N);
	MINUSAN = MINUSI(AN);
	if (N->S > 0) {
		TMP = SUB0_I(X1, 2);
		Z = MULTI(AN, TMP);
		FREEMPI(TMP);
		TMP = ADD0_I(X1, 2);
		ZZ = MULTI(AN, TMP);
		FREEMPI(TMP);
	}
	else {
		TMP = SUB0_I(X1, 2);
		ZZ = MULTI(MINUSAN, TMP);
		FREEMPI(TMP);
		TMP = ADD0_I(X1, 2);
		Z = MULTI(MINUSAN, TMP);
		FREEMPI(TMP);
	}
	FREEMPI(X1);
	FREEMPI(Y1);
	if (N->S > 0) {
		TMP = INT0(N, A);
		G = BIG_MTHROOT(TMP, 2);
		FREEMPI(TMP);
		TMP = MULTI3(A, G, G);
		FREEMPI(G);
		t = EQUALI(TMP, N);
		FREEMPI(TMP);
		SQROOTDR(N, A, 5);
		if (t == 0) {
			printf(" is not an integer\n");
		}
		else {
			printf(" is an integer\n");
			soln++;
		}
		LB = ONEI();
	}
	else {
		TT = MULT_I(MINUSAN, 4);
		TMP = INT0(TT, D);
		G = BIG_MTHROOT(TMP, 2);
		FREEMPI(TMP);
		TMP1 = MULTI3(D, G, G);
		t = EQUALI(TT, TMP1);
		FREEMPI(TMP1);
		SQROOTDR(TT, D, 5);
		FREEMPI(TT);
		if (t == 0) {
			printf(" is not an integer\n");
		}
		else {
			printf(" is an integer\n");
			soln++;
		}
		if (G->D > 0) {
			printf("Lower bound for v is too large\n");
			return(ZEROI());
		}
		LB = ADD0_I(G, 1);
		FREEMPI(G);
	}
	FREEMPI(MINUSAN);
	TMP = INT0(Z, D);
	UB = BIG_MTHROOT(TMP, 2);
	FREEMPI(TMP);
	if (UB->D > 0) {
		printf(" Potentinal upper bound for v is too large\n");
		return(ZEROI());
	}
	TMP = MULTI3(D, UB, UB);
	t = EQUALI(TMP, Z);
	FREEMPI(TMP);
	if (t) {
		ROOTZZ = BIG_MTHROOT(ZZ, 2);
		printf("V = "); PRINTI(UB); printf(" is an integer, U = ");
		PRINTI(ROOTZZ); printf(" is also an integer\n");
		TMP = MULTI(B, UB);
		TEMP = SUBI(ROOTZZ, TMP);
		FREEMPI(ROOTZZ);
		FREEMPI(TMP);
		TMP = MOD(TEMP, TWOA);
		t1 = EQZEROI(TMP);
		FREEMPI(TMP);
		if (t1) {
			TMP = TEMP;
			TEMP = INT(TEMP, TWOA);
			FREEMPI(TMP);
			printf("(U-BV)/2A = "); PRINTI(TEMP); printf(" is an integer and ");
			printf("("); PRINTI(TEMP); printf(","); PRINTI(UB); printf(") is a solution\n");
			soln++;
		}
		else {
			printf("(U-BV)/2A = ");
			RATIO = RATIOI(TEMP, TWOA);
			PRINTDR(5, RATIO);
			FREEMPR(RATIO);
			printf(" is not an integer\n");
		}
		FREEMPI(TEMP);
		TMP = UB;
		UB = SUB0_I(UB, 1);
		FREEMPI(TMP);
	}
	else {
		printf("V = "); SQROOTDR(Z, D, 5); printf(" is not an integer\n");
	}
	if (RSV(UB, LB) >= 0) {
		printf("The remaining solutions lie in the range for v: ["); PRINTI(LB);
		printf(","); PRINTI(UB); printf("]\n");
	}

	lb = (long)CONVERTI(LB);
	ub = (long)CONVERTI(UB);
	for (v = lb; v <= ub; v++) {
		T = MULT_I(B, v);
		TMP = MULT_I(AN, 4);
		TMP1 = CHANGEI(v*v);
		TMP2 = MULTI(D, TMP1);
		FREEMPI(TMP1);
		DD = ADDI(TMP2, TMP);
		FREEMPI(TMP);
		FREEMPI(TMP2);
		if (DD->S < 0) {
			continue;
		}
		G = ISSQUARE(DD);
		if (G->S > 0) {
			T1 = SUBI(G, T);
			TMP = MOD(T1, TWOA);
			t = EQZEROI(TMP);
			FREEMPI(TMP);
			if (t) {
				X = INT(T1, TWOA);
				printf("("); PRINTI(X); printf(",%lu)", v); printf(" is a solution\n");
				FREEMPI(X);
				soln++;
			}
			FREEMPI(T1);
			TMP = ADDI(T, G);
			T1 = MINUSI(TMP);
			FREEMPI(TMP);
			TMP = MOD(T1, TWOA);
			t = EQZEROI(TMP);
			FREEMPI(TMP);
			if (t) {
				X = INT(T1, TWOA);
				printf("("); PRINTI(X); printf(",%lu)", v); printf(" is a solution\n");
				FREEMPI(X);
				soln++;
			}
			FREEMPI(T1);
		}
		FREEMPI(T);
		FREEMPI(G);
		FREEMPI(DD);
	}
	printf("The number of fundamental solutions is %u\n", soln);
	SOLN = CHANGEI(soln);
	FREEMPI(TWOA);
	FREEMPI(AN);
	FREEMPI(D);
	FREEMPI(Z);
	FREEMPI(ZZ);
	FREEMPI(UB);
	FREEMPI(LB);
	return(SOLN);
}

void SQROOTDR(MPI *A, MPI *B, USI r) {
	MPI *T, *TMP, *TMP1, *TMP2, *TENTWOR, *TENTWO2R;

	TENTWOR = POWER_I(10, r);
	TENTWO2R = MULTI(TENTWOR, TENTWOR);
	T = MULTI(TENTWO2R, A);
	TMP1 = INT0(T, B);
	FREEMPI(T);
	TMP = BIG_MTHROOT(TMP1, 2);
	FREEMPI(TMP1);
	TMP1 = INT0(TMP, TENTWOR);
	TMP2 = MOD0(TMP, TENTWOR);
	FREEMPI(TMP);
	FREEMPI(TENTWOR);
	FREEMPI(TENTWO2R);
	printf("squareroot("); PRINTI(A); printf("/"); PRINTI(B); printf(") = "); PRINTI(TMP1); printf("."); PRINTI(TMP2);
	FREEMPI(TMP1);
	FREEMPI(TMP2);
	return;
}

MPI *ISSQUARE(MPI *D) {
	MPI *X, *TEMP;
	USI t;

	X = BIG_MTHROOT(D, 2);
	TEMP = MULTI(X, X);
	t = EQUALI(TEMP, D);
	FREEMPI(TEMP);
	if (t) {
		return(X);
	}
	else {
		FREEMPI(X);
		return(MINUS_ONEI());
	}
}

MPI *STOLTX(MPI *A, MPI *B, MPI *C, MPI *N) {
	MPI *TMP, *TMP1, *D, *G, *SOLN;

	TMP = MULTI(A, C);
	TMP1 = MULT_I(TMP, 4);
	FREEMPI(TMP);
	TMP = MULTI(B, B);
	D = SUBI(TMP, TMP1);
	FREEMPI(TMP);
	FREEMPI(TMP1);

	if (D->S <= 0) {
		printf("D <= 0\n");
		FREEMPI(D);
		return(NULL);
	}

	G = ISSQUARE(D);
	FREEMPI(D);
	if (G->S>0) {
		printf("D is a perfect square\n");
		FREEMPI(G);
		return(NULL);
	}
	FREEMPI(G);
	SOLN = STOLT(A, B, C, N);
	return(SOLN);
}
