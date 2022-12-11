/* elliptic.c */

#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_NONSTDC_NO_DEPRECATE
#include <stdio.h>
#include <stdlib.h>
#include "integer.h"
#include "fun.h"
#include "primepow.h"

extern unsigned int ECMAX;
/*extern unsigned long PRIMEPOWERS[];*/

void ADD_ELLIPTIC_Q(MPR *X1, MPR *Y1, MPR *X2, MPR *Y2, MPR **Xptr, MPR **Yptr, MPR *A, MPR *B)
/* (*Xptr, *Yptr) is the sum of the two points (X1,Y1) and (X2,Y2) on the
* rational elliptic curve y^2=x^3+A*X+B, where 4*A^3+27*b^2 != 0.
*/
{
	MPR *L, *Tmp1, *Tmp2, *Tmp3, *Tmp4, *Tmp5, *Tmp6, *Tmp7;

	Tmp1 = BUILDMPR();
	Tmp1->N = CHANGE(4);
	Tmp1->D = ONEI();
	Tmp2 = BUILDMPR();
	Tmp2->N = CHANGE(27);
	Tmp2->D = ONEI();
	Tmp3 = POWERR(A, 3);
	Tmp4 = POWERR(B, 2);
	Tmp5 = MULTR(Tmp1, Tmp3);
	Tmp6 = MULTR(Tmp2, Tmp4);
	Tmp7 = ADDR(Tmp5, Tmp6);
	if (EQZEROR(Tmp7))
	{
		fprintf(stderr, "discriminant is zero");
		exit(1);
	}
	FREEMPR(Tmp1);
	FREEMPR(Tmp2);
	FREEMPR(Tmp3);
	FREEMPR(Tmp4);
	FREEMPR(Tmp5);
	FREEMPR(Tmp6);
	FREEMPR(Tmp7);

	Tmp1 = SUBR(X1, X2);
	Tmp2 = ADDR(Y1, Y2);
	if (EQZEROR(Tmp1) * EQZEROR(Tmp2))
	{
		fprintf(stderr, "sum is the identity\n");
		exit(1);

	}
	FREEMPR(Tmp2);
	if (!EQZEROR(Tmp1))
	{
		Tmp2 = SUBR(Y1, Y2);
		L = RATIOR(Tmp2, Tmp1);
		FREEMPR(Tmp2);
	}
	else
	{
		Tmp2 = MULTR(X1, X1);
		Tmp3 = BUILDMPR();
		Tmp3->N = CHANGE(3);
		Tmp3->D = ONEI();
		Tmp4 = MULTR(Tmp3, Tmp2);
		Tmp5 = ADDR(Tmp4, A);
		Tmp6 = ADDR(Y1, Y1);
		L = RATIOR(Tmp5, Tmp6);
		FREEMPR(Tmp2);
		FREEMPR(Tmp3);
		FREEMPR(Tmp4);
		FREEMPR(Tmp5);
		FREEMPR(Tmp6);
	}
	FREEMPR(Tmp1);
	Tmp1 = ADDR(X1, X2);
	Tmp2 = MINUSR(Tmp1);
	Tmp3 = MULTR(L, L);
	*Xptr = ADDR(Tmp3, Tmp2);
	FREEMPR(Tmp1);
	FREEMPR(Tmp2);
	FREEMPR(Tmp3);
	Tmp1 = SUBR(*Xptr, X1);
	Tmp2 = MULTR(L, Tmp1);
	Tmp3 = ADDR(Tmp2, Y1);
	*Yptr = MINUSR(Tmp3);
	FREEMPR(L);
	FREEMPR(Tmp1);
	FREEMPR(Tmp2);
	FREEMPR(Tmp3);
	return;
}

/*MPI *XSUB2I();
MPI *ZSUB2I();
MPI *XSUB2IPLUS1();
MPI *ZSUB2IPLUS1(); */

void EXZPOWER(MPI *X, MPI *Z, USI k, MPI *P, MPI *Q, MPI **Aptr, MPI **Cptr, MPI *N)
{
	unsigned int i = 0, w, C[50], l;
	MPI *B, *D, *U, *V, *T, *Tmp;

	U = V = NULL;
	w = k;
	while (w)
	{
		i++;
		C[i] = w % 2;
		w = w / 2;
	}
	l = i;
	*Aptr = COPYI(X);
	*Cptr = COPYI(Z);
	B = XSUB2I(X, Z, P, Q, N);
	D = ZSUB2I(X, Z, P, Q, N);
	for (i = l - 1; i >= 1; i--)
	{
		U = XSUB2IPLUS1(X, Z, *Aptr, *Cptr, B, D, P, Q, N);
		V = ZSUB2IPLUS1(X, *Aptr, *Cptr, B, D, N);
		if (C[i] == 0)
		{
			T = XSUB2I(*Aptr, *Cptr, P, Q, N);
			Tmp = *Cptr;
			*Cptr = ZSUB2I(*Aptr, *Cptr, P, Q, N);
			FREEMPI(Tmp);
			Tmp = *Aptr;
			*Aptr = T;
			FREEMPI(Tmp);
			Tmp = B;
			B = U;
			FREEMPI(Tmp);
			Tmp = D;
			D = V;
			FREEMPI(Tmp);
		}
		else
		{
			T = XSUB2I(B, D, P, Q, N);
			Tmp = D;
			D = ZSUB2I(B, D, P, Q, N);
			FREEMPI(Tmp);
			Tmp = B;
			B = T;
			FREEMPI(Tmp);
			Tmp = *Aptr;
			*Aptr = U;
			FREEMPI(Tmp);
			Tmp = *Cptr;
			*Cptr = V;
			FREEMPI(Tmp);
		}
	}
	if (C[1] == 0)
	{
		FREEMPI(U);
		FREEMPI(V);
	}
	else
	{
		FREEMPI(B);
		FREEMPI(D);
	}
	return;
}

MPI *XSUB2I(MPI *R, MPI *S, MPI *P, MPI *Q, MPI *N)
{
	MPI *T, *G, *H, *I, *Tmp0, *Tmp, *Tmp1, *Tmp2, *Tmp3, *Tmp4, *EIGHT;

	G = MULTI(R, R);
	H = MULTI(S, S);
	I = MULTI(Q, S);
	Tmp = I;
	I = MULTI(I, H);/* I = Q * S * S * S  */
	FREEMPI(Tmp);
	Tmp = H;
	H = MULTI(H, P); /* H = P * S * S  */
	FREEMPI(Tmp);
	Tmp1 = ADDI(G, H); /* Tmp1 = R * R + P * S * S */
	Tmp2 = MULTI(R, Tmp1);
	FREEMPI(Tmp1);
	Tmp = ADDI(Tmp2, I);
	Tmp0 = Tmp;
	Tmp = MOD(Tmp, N);
	FREEMPI(Tmp0);
	FREEMPI(Tmp2);
	if (Tmp->S == 0)
	{
		FREEMPI(I);
		FREEMPI(G);
		FREEMPI(H);
		return (Tmp);
	}
	FREEMPI(Tmp);
	Tmp = SUBI(G, H);
	T = MOD(Tmp, N);
	FREEMPI(Tmp);
	Tmp1 = MULTI(T, T);
	EIGHT = CHANGE(8);
	Tmp2 = MULTI(R, I);
	FREEMPI(I);
	Tmp3 = MULTI(EIGHT, Tmp2);
	Tmp = SUBI(Tmp1, Tmp3);
	Tmp4 = MOD(Tmp, N);
	FREEMPI(T);
	FREEMPI(Tmp);
	FREEMPI(Tmp1);
	FREEMPI(Tmp2);
	FREEMPI(Tmp3);
	FREEMPI(EIGHT);
	FREEMPI(G);
	FREEMPI(H);
	return(Tmp4);
}

MPI *ZSUB2I(MPI *R, MPI *S, MPI *P, MPI *Q, MPI *N)
{
	MPI *Tmp, *Tmp1, *Tmp2, *Tmp3, *Tmp4, *Tmp5, *Tmp6, *Tmp7, *Tmp8, *FOUR;

	Tmp1 = POWERI(R, 3);
	Tmp2 = MULTI(S, S);
	Tmp3 = POWERI(S, 3);
	Tmp4 = MULTI(P, R);
	Tmp5 = MULTI(Tmp2, Tmp4);
	Tmp6 = ADDI(Tmp1, Tmp5);
	Tmp7 = MULTI(Q, Tmp3);
	Tmp = ADDI(Tmp6, Tmp7);
	Tmp8 = MOD(Tmp, N);
	FREEMPI(Tmp);
	FREEMPI(Tmp1);
	FREEMPI(Tmp2);
	FREEMPI(Tmp3);
	FREEMPI(Tmp4);
	FREEMPI(Tmp5);
	FREEMPI(Tmp6);
	FREEMPI(Tmp7);
	if (Tmp8->S == 0)
		return (Tmp8);
	else
	{
		FOUR = CHANGE(4);
		Tmp1 = MULTI(S, Tmp8);
		FREEMPI(Tmp8);
		Tmp = MULTI(FOUR, Tmp1);
		Tmp2 = MOD(Tmp, N);
		FREEMPI(Tmp);
		FREEMPI(Tmp1);
		FREEMPI(FOUR);
		return (Tmp2);

	}
}

MPI *XSUB2IPLUS1(MPI *X, MPI *Z, MPI *R, MPI *S, MPI *U, MPI *V, MPI *P, MPI *Q, MPI *N)
{
	MPI *F, *G, *Tmp1, *Tmp2, *Tmp3, *Tmp4, *Tmp5, *Tmp6, *Tmp7, *Tmp8;
	MPI *FOUR, *Tmp, *Tmp0;

	Tmp1 = MULTI(R, U); /* Tmp1 = R * U */
	Tmp3 = MULTI(V, S); /* Tmp3 = S * V */
	Tmp4 = MULTI(R, V);
	Tmp5 = MULTI(P, Tmp3); /* Tmp5 = P * S * V */
	Tmp = SUBI(Tmp1, Tmp5);
	F = MOD(Tmp, N);/* F = MOD(R * U - P * S * V, N) */
	FREEMPI(Tmp);
	Tmp6 = MULTI(U, S);
	Tmp7 = ADDI(Tmp4, Tmp6); /* Tmp7 = R * V + S * U */
	FREEMPI(Tmp4);
	FREEMPI(Tmp6);
	Tmp8 = MULTI(Q, Tmp3); /* Tmp8 = Q * S * V */
	Tmp = MULTI(Tmp8, Tmp7);
	G = MOD(Tmp, N);
	FREEMPI(Tmp);
	if (X->S != 0)
	{
		FREEMPI(Tmp1);
		FREEMPI(Tmp3);
		FREEMPI(Tmp5);
		FREEMPI(Tmp7);
		FREEMPI(Tmp8);
		Tmp1 = MULTI(F, F);
		FOUR = CHANGE(4);
		Tmp2 = MULTI(FOUR, G);
		Tmp3 = SUBI(Tmp1, Tmp2);
		Tmp = MULTI(Z, Tmp3);
		Tmp4 = MOD(Tmp, N);
		FREEMPI(Tmp);
		FREEMPI(Tmp1);
		FREEMPI(Tmp2);
		FREEMPI(Tmp3);
		FREEMPI(F);
		FREEMPI(G);
		FREEMPI(FOUR);
		return (Tmp4);
	}
	Tmp = F;
	FOUR = CHANGE(4);
	Tmp2 = MULTI(Tmp8, Tmp3);
	F = MULTI(FOUR, Tmp2);
	FREEMPI(Tmp);
	FREEMPI(FOUR);
	Tmp4 = ADDI(Tmp5, Tmp1);
	Tmp6 = ADDI(Tmp4, Tmp4);
	Tmp = G;
	G = MULTI(Tmp6, Tmp7);
	FREEMPI(Tmp);
	Tmp = ADDI(F, G);
	Tmp0 = MOD(Tmp, N);
	FREEMPI(Tmp);
	FREEMPI(F);
	FREEMPI(G);
	FREEMPI(Tmp1);
	FREEMPI(Tmp2);
	FREEMPI(Tmp3);
	FREEMPI(Tmp4);
	FREEMPI(Tmp5);
	FREEMPI(Tmp6);
	FREEMPI(Tmp7);
	FREEMPI(Tmp8);
	return (Tmp0);
}

MPI *ZSUB2IPLUS1(MPI *X, MPI *R, MPI *S, MPI *U, MPI *V, MPI *N)
{
	MPI *Tmp, *T, *Tmp1, *Tmp2, *Tmp3, *Tmp4;

	Tmp1 = MULTI(U, S);
	Tmp2 = MULTI(R, V);
	T = SUBI(Tmp1, Tmp2);
	Tmp3 = MULTI(T, T);
	FREEMPI(T);
	FREEMPI(Tmp1);
	FREEMPI(Tmp2);
	if (X->S != 0)
	{
		Tmp = MULTI(X, Tmp3);
		Tmp4 = MOD(Tmp, N);
		FREEMPI(Tmp);
		FREEMPI(Tmp3);
		return (Tmp4);
	}
	else
	{
		Tmp4 = MOD(Tmp3, N);
		FREEMPI(Tmp3);
		return (Tmp4);
	}
}

MPI *EFACTOR(MPI *N, USI m, USI p)
/*
* This is a modification of algorithm 14.5, page 216, Bressoud.
* This program uses the elliptic curve y^2=x^3-Ax+A, as in Niven,
* Zuckerman, Montgomery, page 285.
* It is assumed that m <= 2328 here.
* n is the composite number, m is the number of prime powers used in
* exponentiation on the elliptic curve, p is the starting seed for the
* random curve y^2=x^3-A*x+A mod n.
* J is the number of curves used.
*/
{
	/* NOTE: m must be > 10, p > 0. */
	MPI *FOUR, *P, *G, *A, *C, *Tmp1, *Tmp2, *Tmp3;
	MPI *ONE, *TWENTY_SEVEN, *Tmp, *AA, *CC;
	unsigned int i, j, k, w, ecmax;

	w = m - 10;
	P = Tmp = NULL;

	ONE = ONEI();
	FOUR = CHANGE(4);
	TWENTY_SEVEN = CHANGE(27);
	if (ECMAX == 0) {
		printf("ECMAX = 0\n");
		FREEMPI(ONE);
		FREEMPI(FOUR);
		FREEMPI(TWENTY_SEVEN);
		return NULL;
	}
	ecmax = ECMAX - 1;
	j = 0;
	while (p)
	{
		if (j > ecmax)
			break;
		if (p % 1 == 0)
		{
			fflush(stdout);
			printf("elliptic curve number %u, A =  %u\n", j + 1, p);
		}
		P = CHANGE(p);
		Tmp = MINUSI(P);
		Tmp1 = POWERI(P, 3);
		Tmp2 = MULTI(Tmp1, FOUR);
		Tmp3 = SUBI(Tmp2, TWENTY_SEVEN);
		G = GCD(Tmp3, N);
		FREEMPI(Tmp1);
		FREEMPI(Tmp2);
		FREEMPI(Tmp3);
		if (!EQUALI(G, ONE))
		{
			printf("discriminant not relatively prime to ");
			PRINTI(N);
			if (!EQUALI(G, N))
			{
				printf("\n found a factor ");
				PRINTI(G);
				printf(" of ");
				PRINTI(N);
				printf("\n");
				FREEMPI(P);
				FREEMPI(Tmp);
				FREEMPI(FOUR);
				FREEMPI(ONE);
				FREEMPI(TWENTY_SEVEN);
				return (G);
			}
			FREEMPI(G);
		}
		else
		{
			FREEMPI(G);
			A = ONEI();
			C = ONEI();
			k = 0;
			while (k <= m)
			{
				printf("k = %u\n", k);
				if (k <= w)
				{
					for (i = 1; i <= 10; i++)
					{
						EXZPOWER(A, C, PRIMEPOWERS[k], Tmp, P, &AA, &CC, N);
						FREEMPI(A);
						FREEMPI(C);
						A = AA;
						C = CC;
						k++;
					}
				}
				else
				{
					EXZPOWER(A, C, PRIMEPOWERS[k], Tmp, P, &AA, &CC, N);
					FREEMPI(A);
					FREEMPI(C);
					A = AA;
					C = CC;
					k++;
				}
				G = GCD(C, N);
				if (!EQUALI(G, ONE) && !EQUALI(G, N))
				{
					printf("Using elliptic curve number %u\n", j + 1);
					printf(" has found a factor ");
					PRINTI(G);
					printf(" of ");
					PRINTI(N);
					printf("\n");
					FREEMPI(P);
					FREEMPI(Tmp);
					FREEMPI(FOUR);
					FREEMPI(ONE);
					FREEMPI(TWENTY_SEVEN);
					FREEMPI(AA);
					FREEMPI(CC);
					return (G);
				}
				FREEMPI(G);
			}
			FREEMPI(A);
			FREEMPI(C);
		}
		FREEMPI(P);
		FREEMPI(Tmp);
		p = RANDOMm(p);
		j++;
	}
	FREEMPI(FOUR);
	FREEMPI(ONE);
	FREEMPI(TWENTY_SEVEN);
	return(NULL);
}

unsigned int ORDERECP(MPI *X, MPI *Z, MPI *P, MPI *Q, MPI *N)
/*
* Calculates the order of the point (X,Y,Z) on the elliptic curve
* Y^2*Z=X^3+P*X*Z^2+Q*Z^3 (mod N), N a prime.
*/
{
	unsigned int k = 1;
	MPI *A, *C;

	while (1)
	{
		EXZPOWER(X, Z, k, P, Q, &A, &C, N);
		FREEMPI(A);
		if (EQZEROI(C))
			break;
		if (k % 100 == 0)
			printf("k = %u\n", k);
		k++;
		FREEMPI(C);
	}
	FREEMPI(C);
	return (k--);
}


