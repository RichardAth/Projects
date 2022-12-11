/* cubicr.c */

#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_NONSTDC_NO_DEPRECATE
#include <stdio.h>
#include <stdlib.h>
#include "integer.h"
#include "fun.h"


void ADD_CUBICR(MPR *X1, MPR *Y1, MPR *X2, MPR *Y2, MPR **Xptr, MPR **Yptr, MPR *A1, MPR *A2, MPR *A3, MPR *A4, MPR *A6)
/*
* (*Xptr,*Yptr) is the sum of the two points (X1,Y1) and (X2,Y2) on the
* elliptic curve y^2+A1*x*y+A3*y=X^3+A2*X^2+A4*x+A6.
* See D. Husemoller, Elliptic curves, page 25.
*/
{
	MPR *Tmp, *Tmp1, *Tmp2, *Tmp3, *Tmp4, *Tmp5, *Tmp6, *Tmp7, *M;

	if (X1 == NULL && X2 != NULL)
	{
		*Xptr = COPYR(X2);
		*Yptr = COPYR(Y2);
		return;
	}
	if (X2 == NULL && X1 != NULL)
	{
		*Xptr = COPYR(X1);
		*Yptr = COPYR(Y1);
		return;
	}
	if (X2 == NULL && X1 == NULL)
	{
		*Xptr = (MPR *)NULL;
		*Yptr = (MPR *)NULL;
		return;
	}
	if (EQUALR(X1, X2))
	{
		if (!EQUALR(Y1, Y2))
		{
			*Xptr = (MPR *)NULL;
			*Yptr = (MPR *)NULL;
			return;
		}
		Tmp2 = ADDR(Y1, Y1);
		Tmp3 = MULTABCR(A3, A1, X1);
		Tmp7 = ADDR(Tmp3, Tmp2);
		FREEMPR(Tmp2);
		FREEMPR(Tmp3);
		if (EQZEROR(Tmp7))
		{
			*Xptr = (MPR *)NULL;
			*Yptr = (MPR *)NULL;
			FREEMPR(Tmp7);
			return;
		}
		else
		{
			Tmp1 = BUILDMPR();
			Tmp1->N = CHANGE(2);
			Tmp1->D = ONEI();
			Tmp2 = ADDR(X1, A2);
			Tmp3 = MULTR(Tmp1, Tmp2);
			FREEMPR(Tmp1);
			FREEMPR(Tmp2);
			Tmp4 = ADDR(X1, Tmp3);
			FREEMPR(Tmp3);
			Tmp5 = MULTR(X1, Tmp4);
			FREEMPR(Tmp4);
			Tmp = Tmp5;
			Tmp5 = ADDR(Tmp5, A4);
			FREEMPR(Tmp);
			Tmp6 = MULTR(A1, Y1);
			Tmp = SUBR(Tmp5, Tmp6);
			FREEMPR(Tmp5);
			FREEMPR(Tmp6);
			M = RATIOR(Tmp, Tmp7);
			FREEMPR(Tmp);
			FREEMPR(Tmp7);
		}
	}
	else
	{
		Tmp1 = SUBR(X1, X2);
		Tmp2 = SUBR(Y1, Y2);
		M = RATIOR(Tmp2, Tmp1);
		FREEMPR(Tmp1);
		FREEMPR(Tmp2);
	}
	Tmp1 = ADDR(M, A1);
	Tmp = Tmp1;
	Tmp1 = MULTR(Tmp1, M);
	FREEMPR(Tmp);
	Tmp2 = ADDR(X1, X2);
	Tmp = Tmp2;
	Tmp2 = ADDR(Tmp2, A2);
	FREEMPR(Tmp);
	Tmp = Tmp2;
	Tmp2 = MINUSR(Tmp2);
	FREEMPR(Tmp);
	*Xptr = ADDR(Tmp1, Tmp2);
	FREEMPR(Tmp1);
	FREEMPR(Tmp2);

	Tmp1 = SUBR(X1, *Xptr);
	Tmp = Tmp1;
	Tmp1 = MULTR(M, Tmp1);
	FREEMPR(Tmp);
	Tmp2 = MULTABCR(A3, A1, *Xptr);
	Tmp = Tmp2;
	Tmp2 = ADDR(Y1, Tmp2);
	FREEMPR(Tmp);
	Tmp = Tmp2;
	Tmp2 = MINUSR(Tmp2);
	FREEMPR(Tmp);
	*Yptr = ADDR(Tmp1, Tmp2);
	FREEMPR(Tmp1);
	FREEMPR(Tmp2);
	FREEMPR(M);
	return;
}

void ADD_CUBICRX()
/*
* Front end for ADD_CUBICR().
*/
{
	MPR  *X1, *X2, *X3, *Y1, *Y2, *Y3, *A1, *A2, *A3, *A4, *A6;
	MPR *Tmp, *Tmp1, *Tmp2, *Tmp3, *Tmp4, *Tmp5;
	MPR *TWO, *THREE, *FOUR, *NINE, *DELTA, *B2, *B4, *B6, *B8;
	unsigned int u;

	printf("Calculating the elliptic curve sum of (X1,Y1) and (X2,Y2)\n");
	printf("for the nonsingular elliptic curve y^2+A1*x*y+A3*y=X^3+A2*X^2+A4*x+A6.\n");
	printf("Enter A1: ");
	A1 = INPUTR(&u);
	printf("Enter A2: ");
	A2 = INPUTR(&u);
	printf("Enter A3: ");
	A3 = INPUTR(&u);
	printf("Enter A4: ");
	A4 = INPUTR(&u);
	printf("Enter A6: ");
	A6 = INPUTR(&u);

	TWO = BUILDMPR();
	TWO->N = CHANGE(2);
	TWO->D = ONEI();
	THREE = BUILDMPR();
	THREE->N = CHANGE(3);
	THREE->D = ONEI();
	FOUR = MULTR(TWO, TWO);
	NINE = MULTR(THREE, THREE);

	Tmp = MULTR(FOUR, A2);
	B2 = MULTABCR(Tmp, A1, A1);
	FREEMPR(Tmp);

	Tmp = MULTR(TWO, A4);
	B4 = MULTABCR(Tmp, A1, A3);
	FREEMPR(Tmp);

	Tmp = MULTR(FOUR, A6);
	B6 = MULTABCR(Tmp, A3, A3);
	FREEMPR(Tmp);

	Tmp1 = MULTR(B2, B6);
	Tmp2 = MULTR(B4, B4);
	Tmp3 = SUBR(Tmp1, Tmp2);
	B8 = RATIOR(Tmp3, FOUR);
	FREEMPR(Tmp1);
	FREEMPR(Tmp2);
	FREEMPR(Tmp3);

	/* Finally DELTA = -B2^2*B8 -8*B4^3 - 27*B6^2 + 9*B2*B4*B6.  */

	Tmp1 = MULTR3(B2, B2, B8);
	Tmp2 = MULTR(TWO, B4);
	Tmp3 = MULTR3(Tmp2, Tmp2, Tmp2);
	Tmp = Tmp3;
	Tmp3 = ADDR(Tmp1, Tmp3);
	FREEMPR(Tmp);
	FREEMPR(Tmp1);
	FREEMPR(Tmp2);

	Tmp2 = MULTR(B2, B4);
	Tmp4 = MULTR(THREE, B6);
	Tmp5 = SUBR(Tmp2, Tmp4);
	FREEMPR(Tmp2);
	FREEMPR(Tmp4);
	Tmp2 = MULTR3(NINE, B6, Tmp5);
	DELTA = SUBR(Tmp2, Tmp3);
	printf("DELTA = "); PRINTR(DELTA); printf("\n");
	FREEMPR(Tmp2);
	FREEMPR(Tmp3);
	FREEMPR(Tmp5);
	FREEMPR(TWO);
	FREEMPR(THREE);
	FREEMPR(FOUR);
	FREEMPR(NINE);
	FREEMPR(B2);
	FREEMPR(B4);
	FREEMPR(B6);
	FREEMPR(B8);

	if (EQZEROR(DELTA))
	{
		fprintf(stderr, "Discriminant is zero.\n");
		exit(1);
	}
	FREEMPR(DELTA);

	while (GetYN())
	{
		printf("Enter X1: ");
		X1 = INPUTR(&u);
		printf("Enter Y1: ");
		Y1 = INPUTR(&u);
		printf("Enter X2: ");
		X2 = INPUTR(&u);
		printf("Enter Y2: ");
		Y2 = INPUTR(&u);
		ADD_CUBICR(X1, Y1, X2, Y2, &X3, &Y3, A1, A2, A3, A4, A6);
		if (X3)
		{
			printf("(X3, Y3) = (");
			PRINTR(X3);
			printf(", ");
			PRINTR(Y3);
			printf(")\n");
			FREEMPR(X3);
			FREEMPR(Y3);
		}
		else
			printf("(X3, Y3) = IDENTITY\n");
		FREEMPR(X1);
		FREEMPR(Y1);
		FREEMPR(X2);
		FREEMPR(Y2);
	}
	FREEMPR(A1);
	FREEMPR(A2);
	FREEMPR(A3);
	FREEMPR(A4);
	FREEMPR(A6);
	return;

}

void POWER_CUBICR(MPR *X1, MPR *Y1, MPR **Xptr, MPR **Yptr, MPR *A1, MPR *A2, MPR *A3, MPR *A4, MPR *A6, unsigned int n)
/*
* (*Xptr, *Yptr)= n(X1, Y1), where 0 <= n < R0 * R0 and (X1, Y1) is on the
* elliptic curve y^2+A1*x*y+A3*y=X^3+A2*X^2+A4*x+A6.
* See D. Husemoller, Elliptic curves, page 25.
*/
{
	MPR *BXtmp, *CXtmp, *BYtmp, *CYtmp;

	*Xptr = (MPR *)NULL;
	*Yptr = (MPR *)NULL;
	if (n == 0)
	{
		printf(" about to return\n");
		return;
	}
	BXtmp = COPYR(X1);
	BYtmp = COPYR(Y1);
	while (1)
	{
		if (n & 1)
		{
			CXtmp = *Xptr;
			CYtmp = *Yptr;
			ADD_CUBICR(BXtmp, BYtmp, CXtmp, CYtmp, Xptr, Yptr, A1, A2, A3, A4, A6);
			FREEMPR(CXtmp);
			FREEMPR(CYtmp);
			if (n == 1)
			{
				FREEMPR(BXtmp);
				FREEMPR(BYtmp);
				return;
			}
		}
		CXtmp = BXtmp;
		CYtmp = BYtmp;
		ADD_CUBICR(CXtmp, CYtmp, CXtmp, CYtmp, &BXtmp, &BYtmp, A1, A2, A3, A4, A6);
		FREEMPR(CXtmp);
		FREEMPR(CYtmp);
		n = n >> 1;
	}
}

void POWER_CUBICRX()
/*
* Front end for POWER_CUBICR().
*/
{
	MPR  *X1, *X2, *Y1, *Y2, *A1, *A2, *A3, *A4, *A6;
	MPR *Tmp, *Tmp1, *Tmp2, *Tmp3, *Tmp4, *Tmp5;
	MPR *TWO, *THREE, *FOUR, *NINE, *DELTA, *B2, *B4, *B6, *B8;
	unsigned int u, n;
	int s;

	printf("Calculating n(X1,Y1) for the nonsingular elliptic curve \n");
	printf("y^2+A1*x*y+A3*y=X^3+A2*X^2+A4*x+A6.\n");
	printf("Enter A1: ");
	A1 = INPUTR(&u);
	printf("Enter A2: ");
	A2 = INPUTR(&u);
	printf("Enter A3: ");
	A3 = INPUTR(&u);
	printf("Enter A4: ");
	A4 = INPUTR(&u);
	printf("Enter A6: ");
	A6 = INPUTR(&u);

	TWO = BUILDMPR();
	TWO->N = CHANGE(2);
	TWO->D = ONEI();
	THREE = BUILDMPR();
	THREE->N = CHANGE(3);
	THREE->D = ONEI();
	FOUR = MULTR(TWO, TWO);
	NINE = MULTR(THREE, THREE);

	Tmp = MULTR(FOUR, A2);
	B2 = MULTABCR(Tmp, A1, A1);
	FREEMPR(Tmp);

	Tmp = MULTR(TWO, A4);
	B4 = MULTABCR(Tmp, A1, A3);
	FREEMPR(Tmp);

	Tmp = MULTR(FOUR, A6);
	B6 = MULTABCR(Tmp, A3, A3);
	FREEMPR(Tmp);

	Tmp1 = MULTR(B2, B6);
	Tmp2 = MULTR(B4, B4);
	Tmp3 = SUBR(Tmp1, Tmp2);
	B8 = RATIOR(Tmp3, FOUR);
	FREEMPR(Tmp1);
	FREEMPR(Tmp2);
	FREEMPR(Tmp3);

	/* Finally DELTA = -B2^2*B8 -8*B4^3 - 27*B6^2 + 9*B2*B4*B6.  */

	Tmp1 = MULTR3(B2, B2, B8);
	Tmp2 = MULTR(TWO, B4);
	Tmp3 = MULTR3(Tmp2, Tmp2, Tmp2);
	Tmp = Tmp3;
	Tmp3 = ADDR(Tmp1, Tmp3);
	FREEMPR(Tmp);
	FREEMPR(Tmp1);
	FREEMPR(Tmp2);

	Tmp2 = MULTR(B2, B4);
	Tmp4 = MULTR(THREE, B6);
	Tmp5 = SUBR(Tmp2, Tmp4);
	FREEMPR(Tmp2);
	FREEMPR(Tmp4);
	Tmp2 = MULTR3(NINE, B6, Tmp5);
	DELTA = SUBR(Tmp2, Tmp3);
	printf("DELTA = "); PRINTR(DELTA); printf("\n");
	FREEMPR(Tmp2);
	FREEMPR(Tmp3);
	FREEMPR(Tmp5);
	FREEMPR(TWO);
	FREEMPR(THREE);
	FREEMPR(FOUR);
	FREEMPR(NINE);
	FREEMPR(B2);
	FREEMPR(B4);
	FREEMPR(B6);
	FREEMPR(B8);

	if (EQZEROR(DELTA))
	{
		fprintf(stderr, "Discriminant is zero.\n");
		exit(1);
	}
	FREEMPR(DELTA);

	while (GetYN())
	{
		printf("Enter X1: ");
		X1 = INPUTR(&u);
		printf("Enter Y1: ");
		Y1 = INPUTR(&u);
		printf("Enter n: ");
		s = scanf("%u", &n);
		printf("n = %u\n", n);
		POWER_CUBICR(X1, Y1, &X2, &Y2, A1, A2, A3, A4, A6, n);
		if (X2)
		{
			printf("(X2, Y2) = (");
			PRINTR(X2);
			printf(", ");
			PRINTR(Y2);
			printf(")\n");
			FREEMPR(X2);
			FREEMPR(Y2);
		}
		else
			printf("(X2, Y2) = IDENTITY\n");
		FREEMPR(X1);
		FREEMPR(Y1);
	}
	FREEMPR(A1);
	FREEMPR(A2);
	FREEMPR(A3);
	FREEMPR(A4);
	FREEMPR(A6);
	return;
}

MPR *MULTR3(MPR *A, MPR *B, MPR *C)
/* Returns A * B * C. */
{
	MPR *Tmp1, *Tmp2;

	Tmp1 = MULTR(A, B);
	Tmp2 = MULTR(Tmp1, C);
	FREEMPR(Tmp1);
	return (Tmp2);
}

MPR *MULTABCR(MPR *A, MPR *B, MPR *C)
/* Returns A + B * C. */
{
	MPR *Temp1, *Temp2;

	Temp1 = MULTR(B, C);
	Temp2 = ADDR(A, Temp1);
	FREEMPR(Temp1);
	return (Temp2);
}

void ADD_CUBICM(MPI *X1, MPI *Y1, MPI *X2, MPI *Y2, MPI **Xptr, MPI **Yptr, MPI *A1, MPI *A2, MPI *A3, MPI *A4, MPI *A6, MPI *MODULUS)
/*
* (*Xptr,*Yptr) is the sum of the two points (X1,Y1) and (X2,Y2) on the
* elliptic curve y^2+A1*x*y+A3*y=X^3+A2*X^2+A4*x+A6 mod p=MODULUS.
* See D. Husemoller, Elliptic curves, page 25.
*/
{
	MPI *Tmp, *Tmp1, *Tmp2, *Tmp3, *Tmp4, *Tmp5, *Tmp6, *Tmp7, *M;

	if (X1 == NULL && X2 != NULL)
	{
		*Xptr = COPYI(X2);
		*Yptr = COPYI(Y2);
		return;
	}
	if (X2 == NULL && X1 != NULL)
	{
		*Xptr = COPYI(X1);
		*Yptr = COPYI(Y1);
		return;
	}
	if (X2 == NULL && X1 == NULL)
	{
		*Xptr = (MPI *)NULL;
		*Yptr = (MPI *)NULL;
		return;
	}
	if (EQUALI(X1, X2))
	{
		if (!EQUALI(Y1, Y2))
		{
			*Xptr = (MPI *)NULL;
			*Yptr = (MPI *)NULL;
			return;
		}
		Tmp2 = ADDM(Y1, Y1, MODULUS);
		Tmp3 = MULTABCM(A3, A1, X1, MODULUS);
		Tmp7 = ADDM(Tmp3, Tmp2, MODULUS);
		FREEMPI(Tmp2);
		FREEMPI(Tmp3);
		if (EQZEROI(Tmp7))
		{
			*Xptr = (MPI *)NULL;
			*Yptr = (MPI *)NULL;
			FREEMPI(Tmp7);
			return;
		}
		else
		{
			Tmp1 = CHANGE(2);
			Tmp2 = ADDM(X1, A2, MODULUS);
			Tmp3 = MULTM(Tmp1, Tmp2, MODULUS);
			FREEMPI(Tmp1);
			FREEMPI(Tmp2);
			Tmp4 = ADDM(X1, Tmp3, MODULUS);
			FREEMPI(Tmp3);
			Tmp5 = MULTM(X1, Tmp4, MODULUS);
			FREEMPI(Tmp4);
			Tmp = Tmp5;
			Tmp5 = ADDM(Tmp5, A4, MODULUS);
			FREEMPI(Tmp);
			Tmp6 = MULTM(A1, Y1, MODULUS);
			Tmp = SUBM(Tmp5, Tmp6, MODULUS);
			FREEMPI(Tmp5);
			FREEMPI(Tmp6);
			M = DIVM(Tmp, Tmp7, MODULUS);
			FREEMPI(Tmp);
			FREEMPI(Tmp7);
		}
	}
	else
	{
		Tmp1 = SUBM(X1, X2, MODULUS);
		Tmp2 = SUBM(Y1, Y2, MODULUS);
		M = DIVM(Tmp2, Tmp1, MODULUS);
		FREEMPI(Tmp1);
		FREEMPI(Tmp2);
	}
	Tmp1 = ADDM(M, A1, MODULUS);
	Tmp = Tmp1;
	Tmp1 = MULTM(Tmp1, M, MODULUS);
	FREEMPI(Tmp);
	Tmp2 = ADDM(X1, X2, MODULUS);
	Tmp = Tmp2;
	Tmp2 = ADDM(Tmp2, A2, MODULUS);
	FREEMPI(Tmp);
	Tmp = Tmp2;
	Tmp2 = MINUSM(Tmp2, MODULUS);
	FREEMPI(Tmp);
	*Xptr = ADDM(Tmp1, Tmp2, MODULUS);
	FREEMPI(Tmp1);
	FREEMPI(Tmp2);

	Tmp1 = SUBM(X1, *Xptr, MODULUS);
	Tmp = Tmp1;
	Tmp1 = MULTM(M, Tmp1, MODULUS);
	FREEMPI(Tmp);
	Tmp2 = MULTABCM(A3, A1, *Xptr, MODULUS);
	Tmp = Tmp2;
	Tmp2 = ADDM(Y1, Tmp2, MODULUS);
	FREEMPI(Tmp);
	Tmp = Tmp2;
	Tmp2 = MINUSM(Tmp2, MODULUS);
	FREEMPI(Tmp);
	*Yptr = ADDM(Tmp1, Tmp2, MODULUS);
	FREEMPI(Tmp1);
	FREEMPI(Tmp2);
	FREEMPI(M);
	return;
}

void ADD_CUBICMX()
/*
* Front end for ADD_CUBICM().
*/
{
	MPI  *X1, *X2, *X3, *Y1, *Y2, *Y3, *A1, *A2, *A3, *A4, *A6;
	MPI *Tmp, *Tmp1, *Tmp2, *Tmp3, *Tmp4, *Tmp5, *MODULUS;
	MPI *TWO, *THREE, *FOUR, *NINE, *DELTA, *B2, *B4, *B6, *B8;
	unsigned int u;

	printf("Calculating the elliptic curve sum of (X1,Y1) and (X2,Y2)\n");
	printf("mod MODULUS\n");
	printf("for the nonsingular elliptic curve y^2+A1*x*y+A3*y=X^3+A2*X^2+A4*x+A6.\n");
	printf("Enter A1: ");
	A1 = INPUTI(&u);
	printf("Enter A2: ");
	A2 = INPUTI(&u);
	printf("Enter A3: ");
	A3 = INPUTI(&u);
	printf("Enter A4: ");
	A4 = INPUTI(&u);
	printf("Enter A6: ");
	A6 = INPUTI(&u);
	printf("Enter MODULUS: ");
	MODULUS = INPUTI(&u);

	TWO = CHANGE(2);
	THREE = CHANGE(3);
	FOUR = MULTM(TWO, TWO, MODULUS);
	NINE = MULTM(THREE, THREE, MODULUS);

	Tmp = MULTM(FOUR, A2, MODULUS);
	B2 = MULTABCM(Tmp, A1, A1, MODULUS);
	FREEMPI(Tmp);

	Tmp = MULTM(TWO, A4, MODULUS);
	B4 = MULTABCM(Tmp, A1, A3, MODULUS);
	FREEMPI(Tmp);

	Tmp = MULTM(FOUR, A6, MODULUS);
	B6 = MULTABCM(Tmp, A3, A3, MODULUS);
	FREEMPI(Tmp);

	Tmp1 = MULTM(B2, B6, MODULUS);
	Tmp2 = MULTM(B4, B4, MODULUS);
	Tmp3 = SUBM(Tmp1, Tmp2, MODULUS);
	B8 = DIVM(Tmp3, FOUR, MODULUS);
	FREEMPI(Tmp1);
	FREEMPI(Tmp2);
	FREEMPI(Tmp3);

	/* Finally DELTA = -B2^2*B8 -8*B4^3 - 27*B6^2 + 9*B2*B4*B6.  */

	Tmp1 = MULTM3(B2, B2, B8, MODULUS);
	Tmp2 = MULTM(TWO, B4, MODULUS);
	Tmp3 = MULTM3(Tmp2, Tmp2, Tmp2, MODULUS);
	Tmp = Tmp3;
	Tmp3 = ADDM(Tmp1, Tmp3, MODULUS);
	FREEMPI(Tmp);
	FREEMPI(Tmp1);
	FREEMPI(Tmp2);

	Tmp2 = MULTM(B2, B4, MODULUS);
	Tmp4 = MULTM(THREE, B6, MODULUS);
	Tmp5 = SUBM(Tmp2, Tmp4, MODULUS);
	FREEMPI(Tmp2);
	FREEMPI(Tmp4);
	Tmp2 = MULTM3(NINE, B6, Tmp5, MODULUS);
	DELTA = SUBM(Tmp2, Tmp3, MODULUS);
	printf("DELTA = "); PRINTI(DELTA); printf("\n");
	FREEMPI(Tmp2);
	FREEMPI(Tmp3);
	FREEMPI(Tmp5);
	FREEMPI(TWO);
	FREEMPI(THREE);
	FREEMPI(FOUR);
	FREEMPI(NINE);
	FREEMPI(B2);
	FREEMPI(B4);
	FREEMPI(B6);
	FREEMPI(B8);

	if (EQZEROI(DELTA))
	{
		fprintf(stderr, "Discriminant is zero.\n");
		exit(1);
	}
	FREEMPI(DELTA);

	while (GetYN())
	{
		printf("Enter X1: ");
		X1 = INPUTI(&u);
		printf("Enter Y1: ");
		Y1 = INPUTI(&u);
		printf("Enter X2: ");
		X2 = INPUTI(&u);
		printf("Enter Y2: ");
		Y2 = INPUTI(&u);
		ADD_CUBICM(X1, Y1, X2, Y2, &X3, &Y3, A1, A2, A3, A4, A6, MODULUS);
		if (X3)
		{
			printf("(X3, Y3) = (");
			PRINTI(X3);
			printf(", ");
			PRINTI(Y3);
			printf(")\n");
			FREEMPI(X3);
			FREEMPI(Y3);
		}
		else
			printf("(X3, Y3) = IDENTITY\n");
		FREEMPI(X1);
		FREEMPI(Y1);
		FREEMPI(X2);
		FREEMPI(Y2);
	}
	FREEMPI(A1);
	FREEMPI(A2);
	FREEMPI(A3);
	FREEMPI(A4);
	FREEMPI(A6);
	return;

}

void POWER_CUBICM(MPI *X1, MPI *Y1, MPI **Xptr, MPI **Yptr, MPI *A1, MPI *A2, MPI *A3, MPI *A4, MPI *A6, unsigned int n, MPI *MODULUS)
/*
* (*Xptr, *Yptr)= n(X1, Y1), where 0 <= n < R0 * R0 and (X1, Y1) is on the
* elliptic curve y^2+A1*x*y+A3*y=X^3+A2*X^2+A4*x+A6 mod MODULUS.
* See D. Husemoller, Elliptic curves, page 25.
*/
{
	MPI *BXtmp, *CXtmp, *BYtmp, *CYtmp;

	*Xptr = (MPI *)NULL;
	*Yptr = (MPI *)NULL;
	if (n == 0)
	{
		printf(" about to return\n");
		return;
	}
	BXtmp = COPYI(X1);
	BYtmp = COPYI(Y1);
	while (1)
	{
		if (n & 1)
		{
			CXtmp = *Xptr;
			CYtmp = *Yptr;
			ADD_CUBICM(BXtmp, BYtmp, CXtmp, CYtmp, Xptr, Yptr, A1, A2, A3, A4, A6, MODULUS);
			FREEMPI(CXtmp);
			FREEMPI(CYtmp);
			if (n == 1)
			{
				FREEMPI(BXtmp);
				FREEMPI(BYtmp);
				return;
			}
		}
		CXtmp = BXtmp;
		CYtmp = BYtmp;
		ADD_CUBICM(CXtmp, CYtmp, CXtmp, CYtmp, &BXtmp, &BYtmp, A1, A2, A3, A4, A6, MODULUS);
		FREEMPI(CXtmp);
		FREEMPI(CYtmp);
		n = n >> 1;
	}
}

void POWER_CUBICMX()
/*
* Front end for POWER_CUBICM().
*/
{
	MPI  *X1, *X2, *Y1, *Y2, *A1, *A2, *A3, *A4, *A6;
	MPI *Tmp, *Tmp1, *Tmp2, *Tmp3, *Tmp4, *Tmp5, *MODULUS;
	MPI *TWO, *THREE, *FOUR, *NINE, *DELTA, *B2, *B4, *B6, *B8;
	unsigned int u, n;
	int s;

	printf("Calculating n(X1,Y1) for the nonsingular elliptic curve \n");
	printf("y^2+A1*x*y+A3*y=X^3+A2*X^2+A4*x+A6 mod MODULUS.\n");
	printf("Enter A1: ");
	A1 = INPUTI(&u);
	printf("Enter A2: ");
	A2 = INPUTI(&u);
	printf("Enter A3: ");
	A3 = INPUTI(&u);
	printf("Enter A4: ");
	A4 = INPUTI(&u);
	printf("Enter A6: ");
	A6 = INPUTI(&u);
	printf("Enter MODULUS: ");
	MODULUS = INPUTI(&u);

	TWO = CHANGE(2);
	THREE = CHANGE(3);
	FOUR = MULTM(TWO, TWO, MODULUS);
	NINE = MULTM(THREE, THREE, MODULUS);

	Tmp = MULTM(FOUR, A2, MODULUS);
	B2 = MULTABCM(Tmp, A1, A1, MODULUS);
	FREEMPI(Tmp);

	Tmp = MULTM(TWO, A4, MODULUS);
	B4 = MULTABCM(Tmp, A1, A3, MODULUS);
	FREEMPI(Tmp);

	Tmp = MULTM(FOUR, A6, MODULUS);
	B6 = MULTABCM(Tmp, A3, A3, MODULUS);
	FREEMPI(Tmp);

	Tmp1 = MULTM(B2, B6, MODULUS);
	Tmp2 = MULTM(B4, B4, MODULUS);
	Tmp3 = SUBM(Tmp1, Tmp2, MODULUS);
	B8 = DIVM(Tmp3, FOUR, MODULUS);
	FREEMPI(Tmp1);
	FREEMPI(Tmp2);
	FREEMPI(Tmp3);

	/* Finally DELTA = -B2^2*B8 -8*B4^3 - 27*B6^2 + 9*B2*B4*B6.  */

	Tmp1 = MULTM3(B2, B2, B8, MODULUS);
	Tmp2 = MULTM(TWO, B4, MODULUS);
	Tmp3 = MULTM3(Tmp2, Tmp2, Tmp2, MODULUS);
	Tmp = Tmp3;
	Tmp3 = ADDM(Tmp1, Tmp3, MODULUS);
	FREEMPI(Tmp);
	FREEMPI(Tmp1);
	FREEMPI(Tmp2);

	Tmp2 = MULTM(B2, B4, MODULUS);
	Tmp4 = MULTM(THREE, B6, MODULUS);
	Tmp5 = SUBM(Tmp2, Tmp4, MODULUS);
	FREEMPI(Tmp2);
	FREEMPI(Tmp4);
	Tmp2 = MULTM3(NINE, B6, Tmp5, MODULUS);
	DELTA = SUBM(Tmp2, Tmp3, MODULUS);
	printf("DELTA = "); PRINTI(DELTA); printf("\n");
	FREEMPI(Tmp2);
	FREEMPI(Tmp3);
	FREEMPI(Tmp5);
	FREEMPI(TWO);
	FREEMPI(THREE);
	FREEMPI(FOUR);
	FREEMPI(NINE);
	FREEMPI(B2);
	FREEMPI(B4);
	FREEMPI(B6);
	FREEMPI(B8);

	if (EQZEROI(DELTA))
	{
		fprintf(stderr, "Discriminant is zero.\n");
		exit(1);
	}
	FREEMPI(DELTA);

	while (GetYN())
	{
		printf("Enter X1: ");
		X1 = INPUTI(&u);
		printf("Enter Y1: ");
		Y1 = INPUTI(&u);
		printf("Enter n: ");
		s = scanf("%u", &n);
		printf("n = %u\n", n);
		POWER_CUBICM(X1, Y1, &X2, &Y2, A1, A2, A3, A4, A6, n, MODULUS);
		if (X2)
		{
			printf("(X2, Y2) = (");
			PRINTI(X2);
			printf(", ");
			PRINTI(Y2);
			printf(")\n");
			FREEMPI(X2);
			FREEMPI(Y2);
		}
		else
			printf("(X2, Y2) = IDENTITY\n");
		FREEMPI(X1);
		FREEMPI(Y1);
	}
	FREEMPI(A1);
	FREEMPI(A2);
	FREEMPI(A3);
	FREEMPI(A4);
	FREEMPI(A6);
	FREEMPI(MODULUS);
	return;
}

MPI *MULTABCM(MPI *A, MPI *B, MPI *C, MPI *M)
/* Returns A + B * C. */
{
	MPI *Temp1, *Temp2;

	Temp1 = MULTM(B, C, M);
	Temp2 = ADDM(A, Temp1, M);
	FREEMPI(Temp1);
	return (Temp2);
}

MPI *MULTM3(MPI *A, MPI *B, MPI *C, MPI *M)
/* Returns A * B * C mod M. */
{
	MPI *Tmp1, *Tmp2;

	Tmp1 = MULTM(A, B, M);
	Tmp2 = MULTM(Tmp1, C, M);
	FREEMPI(Tmp1);
	return (Tmp2);
}


unsigned int ORDER_CUBICM(MPI *X1, MPI *Y1, MPI *A1, MPI *A2, MPI *A3, MPI *A4, MPI *A6, MPI *MODULUS)
/* Returns order of (X1,Y1) on the elliptic curve
* y^2+A1*x*y+A3*y=X^3+A2*X^2+A4*x+A6 mod MODULUS.
*/
{
	unsigned int k;
	MPI *X2, *Y2, *X3, *Y3;
	k = 1;
	if (X1 == NULL)
		return 1;
	else
		X2 = COPYI(X1);
	Y2 = COPYI(Y1);
	while (1) {
		ADD_CUBICM(X1, Y1, X2, Y2, &X3, &Y3, A1, A2, A3, A4, A6, MODULUS);
		FREEMPI(X2);
		FREEMPI(Y2);
		k++;
		if (X3 == NULL)
			break;
		else
		{
			X2 = X3;
			Y2 = Y3;
		}
	}
	return (k);
}

/*unsigned int ORDER_CUBICM(MPI *X1, MPI *Y1, MPI *A1, MPI *A2, MPI *A3, MPI *A4, MPI *A6, MPI *MODULUS)*/
/* Returns order of (X1,Y1) on the elliptic curve
* y^2+A1*x*y+A3*y=X^3+A2*X^2+A4*x+A6 mod MODULUS.
*/
/*{
unsigned int k;
MPI *X2, *Y2;
k = 1;
if (X1 == NULL)
return 1;
while (1){
POWER_CUBICM(X1, Y1, &X2, &Y2, A1, A2, A3, A4, A6, k, MODULUS);
if (X2 == NULL)
break;
else
{
FREEMPI(X2);
FREEMPI(Y2);
k++;
}
}
return (k);
}
*/
void ORDER_CUBICMX()
/*
* Front end for ORDER_CUBICM().
*/
{
	MPI  *X1, *Y1, *A1, *A2, *A3, *A4, *A6;
	MPI *Tmp, *Tmp1, *Tmp2, *Tmp3, *Tmp4, *Tmp5, *MODULUS;
	MPI *TWO, *THREE, *FOUR, *NINE, *DELTA, *B2, *B4, *B6, *B8;
	unsigned int k, u;

	printf("Calculating o(X1,Y1) for the nonsingular elliptic curve \n");
	printf("y^2+A1*xy+A3*y=X^3+A2*X^2+A4*x+A6 modulo a prime.\n");
	printf("Enter A1: ");
	A1 = INPUTI(&u);
	printf("Enter A2: ");
	A2 = INPUTI(&u);
	printf("Enter A3: ");
	A3 = INPUTI(&u);
	printf("Enter A4: ");
	A4 = INPUTI(&u);
	printf("Enter A6: ");
	A6 = INPUTI(&u);
	printf("Enter MODULUS: ");
	MODULUS = INPUTI(&u);

	TWO = CHANGE(2);
	THREE = CHANGE(3);
	FOUR = MULTM(TWO, TWO, MODULUS);
	NINE = MULTM(THREE, THREE, MODULUS);

	Tmp = MULTM(FOUR, A2, MODULUS);
	B2 = MULTABCM(Tmp, A1, A1, MODULUS);
	FREEMPI(Tmp);

	Tmp = MULTM(TWO, A4, MODULUS);
	B4 = MULTABCM(Tmp, A1, A3, MODULUS);
	FREEMPI(Tmp);

	Tmp = MULTM(FOUR, A6, MODULUS);
	B6 = MULTABCM(Tmp, A3, A3, MODULUS);
	FREEMPI(Tmp);

	Tmp1 = MULTM(B2, B6, MODULUS);
	Tmp2 = MULTM(B4, B4, MODULUS);
	Tmp3 = SUBM(Tmp1, Tmp2, MODULUS);
	B8 = DIVM(Tmp3, FOUR, MODULUS);
	FREEMPI(Tmp1);
	FREEMPI(Tmp2);
	FREEMPI(Tmp3);

	/* Finally DELTA = -B2^2*B8 -8*B4^3 - 27*B6^2 + 9*B2*B4*B6.  */

	Tmp1 = MULTM3(B2, B2, B8, MODULUS);
	Tmp2 = MULTM(TWO, B4, MODULUS);
	Tmp3 = MULTM3(Tmp2, Tmp2, Tmp2, MODULUS);
	Tmp = Tmp3;
	Tmp3 = ADDM(Tmp1, Tmp3, MODULUS);
	FREEMPI(Tmp);
	FREEMPI(Tmp1);
	FREEMPI(Tmp2);

	Tmp2 = MULTM(B2, B4, MODULUS);
	Tmp4 = MULTM(THREE, B6, MODULUS);
	Tmp5 = SUBM(Tmp2, Tmp4, MODULUS);
	FREEMPI(Tmp2);
	FREEMPI(Tmp4);
	Tmp2 = MULTM3(NINE, B6, Tmp5, MODULUS);
	DELTA = SUBM(Tmp2, Tmp3, MODULUS);
	printf("DELTA = "); PRINTI(DELTA); printf("\n");
	FREEMPI(Tmp2);
	FREEMPI(Tmp3);
	FREEMPI(Tmp5);
	FREEMPI(TWO);
	FREEMPI(THREE);
	FREEMPI(FOUR);
	FREEMPI(NINE);
	FREEMPI(B2);
	FREEMPI(B4);
	FREEMPI(B6);
	FREEMPI(B8);

	if (EQZEROI(DELTA))
	{
		fprintf(stderr, "Discriminant is zero.\n");
		exit(1);
	}
	FREEMPI(DELTA);

	while (GetYN())
	{
		printf("Enter X1: ");
		X1 = INPUTI(&u);
		printf("Enter Y1: ");
		Y1 = INPUTI(&u);
		k = ORDER_CUBICM(X1, Y1, A1, A2, A3, A4, A6, MODULUS);
		printf("Order of (");
		PRINTI(X1);
		printf(",");
		PRINTI(Y1);
		printf(") on the curve\n");
		printf("y^2");
		if (A1->S)
		{
			printf("+");
			if (!EQONEI(A1))
				PRINTI(A1);
			printf("xy");
		}
		if (A3->S)
		{
			printf("+");
			if (!EQONEI(A3))
				PRINTI(A3);
			printf("y");
		}
		printf("=x^3");
		if (A2->S > 0)
		{
			printf("+");
			if (!EQONEI(A2))
				PRINTI(A2);
			printf("x^2");
		}
		if (A4->S > 0)
		{
			printf("+");
			if (!EQONEI(A4))
				PRINTI(A4);
			printf("x");
		}
		if (A6->S > 0)
		{
			printf("+");
			PRINTI(A6);
		}
		printf(" mod ");
		PRINTI(MODULUS);
		printf("\n");
		printf(" is %u\n ", k);
		FREEMPI(X1);
		FREEMPI(Y1);
	}
	FREEMPI(A1);
	FREEMPI(A2);
	FREEMPI(A3);
	FREEMPI(A4);
	FREEMPI(A6);
	FREEMPI(MODULUS);
	return;
}

unsigned int ORDER_CUBICR(MPR *X1, MPR *Y1, MPR *A1, MPR *A2, MPR *A3, MPR *A4, MPR *A6)
/* Returns order of (X1,Y1) on the elliptic curve
* y^2+A1*x*y+A3*y=X^3+A2*X^2+A4*x+A6. Uses fact that |tor(E)|<= 12.
*/
{
	unsigned int k;
	MPR *X2, *Y2;
	k = 1;
	if (X1 == NULL)
		return 1;
	while (k <= 12) {
		POWER_CUBICR(X1, Y1, &X2, &Y2, A1, A2, A3, A4, A6, k);
		if (X2 == NULL)
			break;
		else
		{
			FREEMPR(X2);
			FREEMPR(Y2);
			k++;
		}
	}
	if (k == 13)
		return (0);
	else
		return (k);
}

void ORDER_CUBICRX()
/*
* Front end for ORDER_CUBICR().
*/
{
	MPR  *X1, *Y1, *A1, *A2, *A3, *A4, *A6;
	MPR *Tmp, *Tmp1, *Tmp2, *Tmp3, *Tmp4, *Tmp5;
	MPR *TWO, *THREE, *FOUR, *NINE, *DELTA, *B2, *B4, *B6, *B8;
	unsigned int k, u;

	printf("Calculating o(X1,Y1) for the nonsingular elliptic curve \n");
	printf("y^2+A1*xy+A3y=X^3+A2*X^2+A4*x+A6.\n");
	printf("Enter A1: ");
	A1 = INPUTR(&u);
	printf("Enter A2: ");
	A2 = INPUTR(&u);
	printf("Enter A3: ");
	A3 = INPUTR(&u);
	printf("Enter A4: ");
	A4 = INPUTR(&u);
	printf("Enter A6: ");
	A6 = INPUTR(&u);
	TWO = BUILDMPR();
	TWO->N = CHANGE(2);
	TWO->D = ONEI();
	THREE = BUILDMPR();
	THREE->N = CHANGE(3);
	THREE->D = ONEI();
	FOUR = MULTR(TWO, TWO);
	NINE = MULTR(THREE, THREE);

	Tmp = MULTR(FOUR, A2);
	B2 = MULTABCR(Tmp, A1, A1);
	FREEMPR(Tmp);

	Tmp = MULTR(TWO, A4);
	B4 = MULTABCR(Tmp, A1, A3);
	FREEMPR(Tmp);

	Tmp = MULTR(FOUR, A6);
	B6 = MULTABCR(Tmp, A3, A3);
	FREEMPR(Tmp);

	Tmp1 = MULTR(B2, B6);
	Tmp2 = MULTR(B4, B4);
	Tmp3 = SUBR(Tmp1, Tmp2);
	B8 = RATIOR(Tmp3, FOUR);
	FREEMPR(Tmp1);
	FREEMPR(Tmp2);
	FREEMPR(Tmp3);

	/* Finally DELTA = -B2^2*B8 -8*B4^3 - 27*B6^2 + 9*B2*B4*B6.  */

	Tmp1 = MULTR3(B2, B2, B8);
	Tmp2 = MULTR(TWO, B4);
	Tmp3 = MULTR3(Tmp2, Tmp2, Tmp2);
	Tmp = Tmp3;
	Tmp3 = ADDR(Tmp1, Tmp3);
	FREEMPR(Tmp);
	FREEMPR(Tmp1);
	FREEMPR(Tmp2);

	Tmp2 = MULTR(B2, B4);
	Tmp4 = MULTR(THREE, B6);
	Tmp5 = SUBR(Tmp2, Tmp4);
	FREEMPR(Tmp2);
	FREEMPR(Tmp4);
	Tmp2 = MULTR3(NINE, B6, Tmp5);
	DELTA = SUBR(Tmp2, Tmp3);
	printf("DELTA = "); PRINTR(DELTA); printf("\n");
	FREEMPR(Tmp2);
	FREEMPR(Tmp3);
	FREEMPR(Tmp5);
	FREEMPR(TWO);
	FREEMPR(THREE);
	FREEMPR(FOUR);
	FREEMPR(NINE);
	FREEMPR(B2);
	FREEMPR(B4);
	FREEMPR(B6);
	FREEMPR(B8);

	if (EQZEROR(DELTA))
	{
		fprintf(stderr, "Discriminant is zero.\n");
		exit(1);
	}
	FREEMPR(DELTA);

	while (GetYN())
	{
		printf("Enter X1: ");
		X1 = INPUTR(&u);
		printf("Enter Y1: ");
		Y1 = INPUTR(&u);
		k = ORDER_CUBICR(X1, Y1, A1, A2, A3, A4, A6);
		printf("Order of (");
		PRINTR(X1);
		printf(",");
		PRINTR(Y1);
		printf(") on the curve\n");
		printf("y^2");
		if ((A1->N)->S > 0)
		{
			printf("+");
			if (!EQONER(A1))
				PRINTR(A1);
			printf("xy");
		}
		else if ((A1->N)->S < 0)
		{
			if (!EQMINUSONER(A1))
				PRINTR(A1);
			else
				printf("-");
			printf("xy");
		}
		if ((A3->N)->S > 0)
		{
			printf("+");
			if (!EQONER(A3))
				PRINTR(A3);
			printf("y");
		}
		else if ((A3->N)->S < 0)
		{
			if (!EQMINUSONER(A3))
				PRINTR(A3);
			else
				printf("-");
			printf("y");
		}
		printf("=x^3");
		if ((A2->N)->S > 0)
		{
			printf("+");
			if (!EQONER(A2))
				PRINTR(A2);
			printf("x^2");
		}
		else if ((A2->N)->S < 0)
		{
			if (!EQMINUSONER(A2))
				PRINTR(A2);
			else
				printf("-");
			printf("x^2");
		}
		if ((A4->N)->S > 0)
		{
			printf("+");
			if (!EQONER(A4))
				PRINTR(A4);
			printf("x");
		}
		else if ((A4->N)->S < 0)
		{
			if (!EQMINUSONER(A4))
				PRINTR(A4);
			else
				printf("-");
			printf("x");
		}
		if ((A6->N)->S > 0)
		{
			printf("+");
			PRINTR(A6);
		}
		else if ((A6->N)->S < 0)
			PRINTR(A6);
		printf("\n");
		if (k)
			printf(" is %u\n ", k);
		else
			printf(" is infinite\n ");
		FREEMPR(X1);
		FREEMPR(Y1);
	}
	FREEMPR(A1);
	FREEMPR(A2);
	FREEMPR(A3);
	FREEMPR(A4);
	FREEMPR(A6);
	return;
}
