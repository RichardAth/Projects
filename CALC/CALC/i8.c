/* i8.c */
#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_NONSTDC_NO_DEPRECATE
#include <stdio.h>
#include "integer.h"
#include "fun.h"

MPI *APPROXN(MPI *Aptr, unsigned int m)
/*
*  *Eptr = 2 ^ (1 + [i / m]), (where i + 1 is the number of binary digits in
the positive  MPI *Aptr) is the initial approximation to the integer part of
the m-th root of *Aptr.  (*Aptr)^(1/m) < *Eptr <= 2*(*Aptr)^(1/m).
*/
{
	MPI *Eptr;
	unsigned int i = 0;

	i = BINARYB(Aptr) - 1;
	Eptr = POWER_I(2L, 1 + (i / m));
	return (Eptr);
}

MPI *NEWTON(MPI *Aptr, unsigned int m, MPI *Fptr)
/*
* Returns the integer part of the mth root of *Aptr, where m is
* an unsigned integer, 1 < m < R0.
* The method is Newton's, using the integer part function, with initial
* approximation *Fptr > (*Aptr)^(1/m).
*/
{
	MPI *X, *Y, *Z, *Temp, *Temp1;

	X = COPYI(Fptr);
	while (1)
	{
		Y = COPYI(X);
		Z = POWERI(X, m - 1);
		Temp = INT0(Aptr, Z);
		FREEMPI(Z);
		Z = MULT_I(X, m - 1);
		Temp1 = ADD0I(Temp, Z);
		FREEMPI(X);
		FREEMPI(Temp);
		X = INT0_(Temp1, m);
		FREEMPI(Temp1);
		FREEMPI(Z);
		if (RSV(Y, X) != 1)
		{
			FREEMPI(X);
			return (Y);
		}
		FREEMPI(Y);
	}
}

MPI *BABY_MTHROOT(MPI *Aptr, unsigned int m)
/*
* *Eptr is the integer part of the mth root of *Aptr, where m is
* an unsigned integer, 1 < m < R0.
* The method is Newton's, using the integer part function, with initial
* approximation provided by APPROXN(), for "small" *Aptr.
*/
{
	MPI *X, *Eptr;

	X = APPROXN(Aptr, m);
	Eptr = NEWTON(Aptr, m, X);
	FREEMPI(X);
	return (Eptr);
}

MPI *BIG_MTHROOT(MPI *Aptr, unsigned int m)
/*
* *Eptr, the integer part of the mth root of the positive MPI *Aptr,
* 1 < m < R0, is obtained by Newton's method, using the integer part function.
*/
{
	MPI *a, *X, *Eptr, *Temp;
	unsigned int i, n, r, t;

	if (EQZEROI(Aptr)) {
		Temp = ZEROI();
		return(Temp);
	}
	n = Aptr->D;
	if (n < m)
	{
		Eptr = BABY_MTHROOT(Aptr, m);
		return (Eptr);
	}
	else
	{
		r = n / m;
		a = BUILDMPI((n % m) + 1); /* a = [*Aptr/R0^(m*r)] */
		a->S = 1;
		t = r * m;
		for (i = 0; i <= a->D; i++)
			a->V[i] = Aptr->V[i + t];
		X = BABY_MTHROOT(a, m);
		FREEMPI(a);
		Temp = X;
		X = ADD0_I(X, (USL)1);
		FREEMPI(Temp);
		Temp = X;
		X = BUILDMPI(Temp->D + r + 1); /* X = X* R0^r */
		X->S = 1;
		for (i = X->D; i >= r; i--)
			X->V[i] = Temp->V[i - r];
		FREEMPI(Temp);
		INTSETUL(X->V, (USL)r, (USL)0);
		/* sets the first r array elements to 0 */

		Eptr = NEWTON(Aptr, m, X); /* X is the initial approximation */
		FREEMPI(X);
		return (Eptr);
	}
}

void MTHROOT(MPI *Aptr, MPI *Bptr, unsigned int m, unsigned int r)
/*
* *Aptr and *Bptr are positive MPI'S.
* the mthroot of *Aptr/(*Bptr) is computed to r decimal places, r >= 0;
* m, r are integers, 0 < m * r < R0 * R0.
*/
{
	MPI *D, *E, *F, *G, *Y, *Temp;
	unsigned int i, l;

	E = POWER_I(10L, r * m);
	Temp = E;
	E = MULTI(E, Aptr);
	FREEMPI(Temp);
	Temp = E;
	E = INT0(E, Bptr);
	FREEMPI(Temp);
	if (E->S == 0)
	{
		FREEMPI(E);
		if (r == 0)
		{
			printf("%uthroot of ", m);
			PRINTI(Aptr);
			printf("/");
			PRINTI(Bptr);
			printf(" = 0\n");
			return;
		}
		else
		{
			printf("%uthroot of ", m);
			PRINTI(Aptr);
			printf("/");
			PRINTI(Bptr);
			printf(" = 0.");
			for (i = 1; i <= r; i++)
				printf("0");
			printf("\n");
			return;
		}
	}
	Y = BIG_MTHROOT(E, m);
	FREEMPI(E);
	F = POWER_I(10L, r);
	G = MOD0(Y, F);
	D = INT0(Y, F);
	FREEMPI(F);
	FREEMPI(Y);
	printf("%uthroot of ", m);
	PRINTI(Aptr);
	printf("/");
	PRINTI(Bptr);
	printf(" = ");
	PRINTI(D);
	FREEMPI(D);
	if (r == 0)
	{
		FREEMPI(G);
		printf("\n");
		return;
	}
	printf(".");
	l = LENGTHI(G); /* the number of decimal digits in G */
	for (i = 1; i <= r - l; i++)
		printf("0");
	PRINTI(G);
	printf("\n");
	FREEMPI(G);
	return;
}

MPR *MTHROOTR(MPR *Nptr, unsigned int m, unsigned int r)
/*
* *Nptr is a positive MPR.
* the mthroot of *Nptr is computed to r decimal places, r >= 0 .
* m, r are integers, 0 < m * r < R0 * R0 .
*/
{
	MPI *Tmp0I, *Tmp1I, *Tmp2I;
	MPR *Tmp1R;

	Tmp0I = POWER_I(10L, r);
	Tmp1I = POWERI(Tmp0I, m);
	Tmp2I = MULTI(Tmp1I, Nptr->N);
	Tmp1R = RATIOI(Tmp2I, Nptr->D);
	FREEMPI(Tmp1I);
	FREEMPI(Tmp2I);
	Tmp1I = INTI(Tmp1R->N, Tmp1R->D);
	FREEMPR(Tmp1R);
	if (EQZEROI(Tmp1I))
	{
		FREEMPI(Tmp1I);
		FREEMPI(Tmp0I);
		return (ZEROR());
	}
	Tmp2I = BIG_MTHROOT(Tmp1I, m);
	FREEMPI(Tmp1I);
	Tmp1R = RATIOI(Tmp2I, Tmp0I);
	FREEMPI(Tmp0I);
	FREEMPI(Tmp2I);
	return (Tmp1R);
}

MPI *PI(unsigned int r)
{
	MPI *PI, *X, *X1, *Y, *U, *V;
	MPI *Tmp, *Tmp0I, *Tmp1I, *Tmp2I, *Tmp3I, *Tmp4I, *Tmp5I;

	Tmp0I = POWER_I(10L, r);
	Tmp1I = POWERI(Tmp0I, 2);/* 10^2r */
	Tmp2I = POWERI(Tmp1I, 2);/* 10^4r */
	Tmp3I = POWERI(Tmp2I, 2);/* 10^8r */
	Tmp = MULT_I(Tmp2I, 2L);
	X = BIG_MTHROOT(Tmp, 2);
	FREEMPI(Tmp);
	Tmp4I = MULT_I(Tmp1I, 2L);
	PI = ADD0I(Tmp4I, X);
	FREEMPI(Tmp4I);
	Tmp4I = MULT_I(Tmp3I, 2L);
	Y = BIG_MTHROOT(Tmp4I, 4);
	FREEMPI(Tmp3I);
	FREEMPI(Tmp4I);
	while (1)
	{
		Tmp4I = ADD0I(X, Tmp1I);
		Tmp = X;
		X = MULTI(Tmp0I, Tmp4I);
		FREEMPI(Tmp);
		Tmp3I = BIG_MTHROOT(X, 2);
		X1 = MULT_I(Tmp3I, 2L);
		Tmp = X;
		X = INT0(X, X1);
		FREEMPI(Tmp);
		FREEMPI(X1);

		U = MULTI(X, Y);
		Tmp = U;
		U = ADD0I(U, Tmp2I);
		printf(" U = ");
		PRINTI(U);
		printf("\n");
		FREEMPI(Tmp);
		Tmp5I = ADD0I(Y, Tmp1I);
		V = MULTI(Tmp3I, Tmp5I);
		Tmp = V;
		V = MULTI(V, Tmp0I);
		FREEMPI(Tmp);
		Tmp = Y;
		Y = INT0(U, V);
		printf(" V = ");
		PRINTI(V);
		printf("\n");
		FREEMPI(U);
		FREEMPI(V);
		printf(" Y = ");
		PRINTI(Y);
		printf("\n");

		if (EQUALI(Tmp1I, Y))
		{
			FREEMPI(Tmp0I);
			FREEMPI(Tmp1I);
			FREEMPI(Tmp2I);
			FREEMPI(Tmp3I);
			FREEMPI(Tmp4I);
			FREEMPI(Tmp5I);
			return (PI);
		}
		Tmp = PI;
		PI = MULTI(PI, Tmp4I);
		FREEMPI(Tmp);
		Tmp = PI;
		PI = INT0(PI, Tmp5I);
		FREEMPI(Tmp4I);
		FREEMPI(Tmp5I);
	}
}

MPR *PII(unsigned int r)
/*
* From " A very rapidly convergent product expansion for pi," Bit 23 (1983)
* 538-540, by J.M. Borwein and P.B. Borwein. The algorithm delivers pi
* to 2^r decimals.
*/
{
	MPR * Tmp1, *Tmp2, *Tmp3, *Tmp4, *Tmp5, *Tmp6, *Tmp7, *Tmp8;
	MPR *X, *Y, *PI, *ONE, *TWO, *Tmp, *TEST, *SS, *TEN, *Tmp9;
	MPI *S, *TmpI;
	unsigned int s;

	S = POWER_I(2L, r);
	s = CONVERTI(S);
	FREEMPI(S);
	SS = BUILDMPR();
	SS->N = POWER_I(10L, s + 1);
	SS->D = ONEI();
	ONE = ONER();
	TWO = BUILDMPR();
	TWO->N = CHANGE((USL)2);
	TWO->D = ONEI();
	TEN = BUILDMPR();
	TEN->N = CHANGE((USL)10);
	TEN->D = ONEI();
	X = MTHROOTR(TWO, 2, s + 1);
	PI = ADDR(TWO, X);
	Tmp1 = MTHROOTR(X, 2, s + 1);
	Tmp2 = INVERSER(Tmp1);
	Tmp3 = ADDR(Tmp1, Tmp2);
	Tmp = X;
	X = RATIOR(Tmp3, TWO);
	printf("X = ");
	PRINTDR(s, X);
	printf("\n");
	FREEMPR(Tmp);
	FREEMPR(Tmp1);
	FREEMPR(Tmp2);
	FREEMPR(Tmp3);
	Y = MTHROOTR(TWO, 4, s + 1);
	while (1)
	{
		Tmp1 = MTHROOTR(X, 2, s + 1);
		Tmp2 = INVERSER(Tmp1);
		Tmp3 = ADDR(Tmp1, Tmp2);
		Tmp4 = RATIOR(Tmp3, TWO);
		Tmp5 = MULTR(Y, Tmp1);
		Tmp6 = ADDR(Tmp5, Tmp2);
		Tmp7 = ADDR(Y, ONE);
		Tmp8 = ADDR(X, ONE);
		Tmp = PI;
		PI = MULTR(Tmp8, PI);
		FREEMPR(Tmp);
		Tmp = PI;
		PI = RATIOR(PI, Tmp7);
		FREEMPR(Tmp);
		Tmp = X;
		X = Tmp4;
		FREEMPR(Tmp);
		printf("X = ");
		PRINTDR(s, X);
		printf("\n");
		Tmp = Y;
		Y = RATIOR(Tmp6, Tmp7);
		FREEMPR(Tmp);
		printf("Y = ");
		PRINTDR(s, Y);
		printf("\n");
		TEST = SUBR(Y, ONE);
		Tmp = TEST;
		Tmp9 = MULTR(SS, TEN);
		TEST = MULTR(Tmp9, TEST);
		FREEMPR(Tmp);
		TmpI = INT0(TEST->N, TEST->D);
		FREEMPR(TEST);
		if (TmpI->S == 0)
		{
			printf("Y - 1 = ");
			PRINTDR(s, Y);
			printf("\n");
			FREEMPR(Tmp1);
			FREEMPR(Tmp2);
			FREEMPR(Tmp3);
			FREEMPR(Tmp5);
			FREEMPR(Tmp6);
			FREEMPR(Tmp7);
			FREEMPR(Tmp8);
			FREEMPR(Tmp9);
			FREEMPI(TmpI);
			FREEMPR(ONE);
			FREEMPR(TWO);
			FREEMPR(TEN);
			FREEMPR(SS);
			FREEMPR(X);
			FREEMPR(Y);
			return (PI);
		}
		FREEMPR(Tmp1);
		FREEMPR(Tmp2);
		FREEMPR(Tmp3);
		FREEMPR(Tmp5);
		FREEMPR(Tmp6);
		FREEMPR(Tmp7);
		FREEMPR(Tmp8);
		FREEMPR(Tmp9);
		FREEMPI(TmpI);
	}
}

MPI *SQRT4TIMING(MPI *D, MPI *INTSQRTD, MPI *Y, int i)
/*
* This returns the squareroot of perfect square DY^2-i, i = 1 or -1.
* Here X=INTSQROOTD = int(sqrt(D)).
* The method is Newton's, using the integer part function, with initial
* approximation (X+1)Y if i=1, (X+1)Y+1 if i = -1.
*/
{
	MPI *X, *Z, *TEMP, *TEMP1, *TEMP2, *A;

	TEMP = ADD0_I(INTSQRTD, 1);
	TEMP1 = MULTI(TEMP, Y);
	FREEMPI(TEMP);
	TEMP2 = MULTI3(D, Y, Y);
	if (i == 1) {
		X = TEMP1;
		A = SUB0_I(TEMP2, 1);
	}
	else {
		X = ADD0_I(TEMP1, 1);
		A = ADD0_I(TEMP2, 1);
		FREEMPI(TEMP1);
	}
	FREEMPI(TEMP2);
	Z = NEWTON(A, 2, X);
	FREEMPI(X);
	FREEMPI(A);
	return (Z);
}

