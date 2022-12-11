/* The Pohst Algorithm updating only from where necessary */
#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_NONSTDC_NO_DEPRECATE
#include <stdio.h>
#include <stdlib.h>
#include "integer.h"
#include "fun.h"
extern MPI *MAXI, *PMAXI;

extern unsigned int MLLLVERBOSE;
extern unsigned int HERMITEVERBOSE;
unsigned int GCDFLAG;

MPMATI *BASIS_REDUCTION(MPMATI *Bptr, MPMATI **Eptr, USI rowstage, USI m1, USI n1)
/*
* Input: *Bptr, a matrix of MPI's, whose first row is not zero.
* Output: a pointer to an MPMATI whose rows form a reduced basis for
* the lattice spanned by the rows of *Bptr. This basis is reduced in the
* sense of the paper "Factoring polynomials with rational coefficients" by
* A. K. Lenstra, H. W. Lenstra and L. Lovasz, Math. Ann. 261, 515-534 (1982)
* using the modified version in "Solving exponential Diophantine equations
* using lattice basis reduction algorithms" by B. M. M. De Weger, J. No. Theory
* 26, 325-367 (1987). A change of basis matrix **Eptr is also returned.
* De Weger's algorithm has been changed to cater for arbitrary matrices. The
* the rows are now in general linearly dependent.
* We use the fact that the Gram Schmidt process detects the first row
* which is a linear combination of the preceding rows. We employ a modification
* of the LLL algorithm outlined by M. Pohst in J. Symbolic Computation (1987)4,
* 123-127.  We call this the MLLL algorithm.
* The last sigma rows of the matrix **Eptr are relation vectors.
* m1 / n1 is usually taken to be 3 / 4.
*/
{
	unsigned int i, k, l, n, nn, m, t, flag = 0, Flag = 0;
	unsigned int flagg, beta, K1 = 0, tau = 2, sigma = 0, rho;
	MPI **D, *X, *Y, *Z, *H, *Tmp, *R, *M1, *N1;
	MPMATI *C, *L, *B1ptr;

	M1 = CHANGE(m1);
	N1 = CHANGE(n1);
	m = Bptr->C;
	nn = Bptr->R;
	n = nn;
	/* We initial Eptr outside the function whenever we call the function. */
	/* This is because we have to do so in SMITH(). */

	/*	MAXI = MAXELTI(Bptr);
	PMAXI = MAXELTI(*Eptr); */
	B1ptr = COPYMATI(Bptr);
	D = (MPI **)mmalloc((1 + n) * sizeof(MPI *));
	D[0] = ONEI();
	for (i = 1; i <= n; i++)
		D[i] = ZEROI();
	C = ZEROMNI(n, m);
	L = ZEROMNI(n, n);

found:
	n = B1ptr->R;
	i = (K1 == 0) ? 1 : K1;
	/* K1 = no. of consecutive rows of *B1ptr that don't need updating
	for the Gram Schmidt process. */
	while (i <= B1ptr->R)
	{
		BASIS_UPDATE(i, m, &C, &L, B1ptr, D);
		flag = 1;
		for (t = 0; t < m; t++)
		{
			if (!EQZEROI(C->V[i - 1][t]))
			{
				flag = 0;
				break;
			}
		}
		if (flag)
			break;
		X = ZEROI();
		for (t = 0; t < m; t++)
		{
			H = MULTI(C->V[i - 1][t], C->V[i - 1][t]);
			Tmp = X;
			X = ADDI(X, H);
			FREEMPI(Tmp);
			FREEMPI(H);
		}
		FREEMPI(D[i]);
		D[i] = INT(X, D[i - 1]);
		FREEMPI(X);
		i++;
		if (MLLLVERBOSE)
			printf("i = %u\n", i);
	}
	beta = (flag) ? i : i - 1;
	rho = K1 = i - 1;
	if (MLLLVERBOSE)
		printf("completed updating the basis\n");
	/* Here K1 = no. of LI rows in *B1ptr found by Gram Schmidt process.
	flag = 0 means all the rho = beta rows of *B1ptr are LI;
	flag = 1 means that the first rho = beta - 1 rows of *B1ptr are LI, but the
	beta-th row is a LC of the preceding rows. So beta = number of rows of *B1ptr
	currently being examined by the LLL algorithm. */
	k = tau;
	if (MLLLVERBOSE)
		printf("beta = %u, k = %u\n", beta, k);
	while (k <= beta)
	{
		if (MLLLVERBOSE)
			printf("beta - k = %u - %u, %u\n", beta, k, beta - k);
		l = k - 1;
		Flag = STEP4(k, l, &L, &B1ptr, Eptr, D, rowstage);
		if (Flag)/* STEP 9 of POHST. */
		{
			sigma++;
			if (MLLLVERBOSE)
				printf("relation vector number %u found\n", sigma);
			tau = k++;
			goto found;
		}
		X = MULTI(D[k - 2], D[k]);
		Y = MULTI(D[k - 1], D[k - 1]);
		Tmp = Y;
		Y = MULT_I(Y, m1);
		FREEMPI(Tmp);
		R = MULTI(L->V[k - 1][k - 2], L->V[k - 1][k - 2]);
		Z = ADD0I(X, R);
		Tmp = Z;
		Z = MULT_I(Z, n1);
		FREEMPI(Tmp);
		if (RSV(Y, Z) == 1)/*& STEP 5 of POHST. */
		{
			flagg = 0;
			if (EQZEROI(D[k]) && EQZEROI(R))
			{/* CASE B=0 of STEP 7 of POHST. */
				FREEMPI(D[k - 1]);
				D[k - 1] = ZEROI();
				STEP8(k, &B1ptr, &L, Eptr, rowstage);
				if (k - 1 < K1)
					K1 = k - 1;
				/* The swap may have changed 2nd last row */
				/* of *B1ptr. */
				for (t = 0; t < m; t++)
				{
					FREEMPI(C->V[k - 2][t]);
					C->V[k - 2][t] = ZEROI();
				}
				beta--;
				flagg = 1;
				if (k > 2)
					k--;
				FREEMPI(X);
				FREEMPI(Y);
				FREEMPI(R);
				FREEMPI(Z);
				continue;
			}
			FREEMPI(Z);
			if (flagg == 0)
			{
				for (i = k + 1; i <= beta; i++)
					STEP7(i, k, &L, D);
			}
			STEP8(k, &B1ptr, &L, Eptr, rowstage);
			if (k - 2 < K1)
				K1 = k - 2;
			/* swap will change last two rows of *B1ptr. */
			FREEMPI(R);
			FREEMPI(Y);
			Y = MULTI(L->V[k - 1][k - 2], L->V[k - 1][k - 2]);
			Tmp = Y;
			Y = ADD0I(Y, X);
			FREEMPI(X);
			FREEMPI(Tmp);
			Tmp = D[k - 1];
			D[k - 1] = INT0(Y, D[k - 1]);
			FREEMPI(Tmp);
			FREEMPI(Y);
			if (k > 2)
				k--;
		}
		else
		{ /* STEP 6 of POHST. */
			FREEMPI(R);
			FREEMPI(X);
			FREEMPI(Y);
			FREEMPI(Z);
			for (l = k - 2; l >= 1; l--)
			{
				Flag = STEP4(k, l, &L, &B1ptr, Eptr, D, rowstage);
				if (Flag)
				{
					sigma++;
					if (MLLLVERBOSE)
						printf("relation vector number %u found\n", sigma);
					tau = k++;
					goto found; /* STEP 9 of POHST. */
				}
			}
			k++;
		}
	}
	FREEMPI(M1);
	FREEMPI(N1);
	FREEMATI(C);
	FREEMATI(L);
	for (i = 0; i <= nn; i++)
		FREEMPI(D[i]);
	/*ffree((char *)D, (1 + Bptr->R) * sizeof(MPI *));*/
	ffree((char *)D, (1 + nn) * sizeof(MPI *));
	if (MLLLVERBOSE)
	{
		printf("number of basis vectors found = %u ;\n", rho);
		printf("number of relation vectors found = %u .\n", sigma);
	}
	return (B1ptr);
}

unsigned int STEP4(k, l, Lptr, Bptr, Eptr, D, i)
/*
* updates *Lptr, *Bptr and *Eptr.
* returns 1 if row k of *Bptr becomes zero, returns zero otherwise.
*/
unsigned int k, l, i;
MPI *D[];
MPMATI **Lptr, **Bptr, **Eptr;
{
	unsigned int j, flag = 1, t, m, n;
	MPI *X, *Y, *R, *Tmp;
	MPMATI *TmpMATI;

	m = (*Bptr)->C;
	n = (*Eptr)->R;
	Y = MULT_I((*Lptr)->V[k - 1][l - 1], 2);
	if (RSV(Y, D[l]) == 1)
	{
		R = NEAREST_INTI((*Lptr)->V[k - 1][l - 1], D[l]);
		X = MINUSI(R);
		TmpMATI = *Bptr;
		*Bptr = ADD_MULT_ROWI(l - 1, k - 1, X, *Bptr);
		/*
		MAXI = UPDATEMAXI(MAXI, *Bptr);
		*/
		FREEMATI(TmpMATI);
		TmpMATI = *Eptr;
		*Eptr = ADD_MULT_ROWI(l + i - 1, k + i - 1, X, *Eptr);
		/*
		PMAXI = UPDATEMAXI(PMAXI, *Eptr);
		*/
		FREEMATI(TmpMATI);
		FREEMPI(X);
		for (j = 1; j < l; j++)
		{
			X = MULTI((*Lptr)->V[l - 1][j - 1], R);
			Tmp = (*Lptr)->V[k - 1][j - 1];
			(*Lptr)->V[k - 1][j - 1] = SUBI((*Lptr)->V[k - 1][j - 1], X);
			FREEMPI(Tmp);
			FREEMPI(X);
		}
		X = MULTI(D[l], R);
		Tmp = (*Lptr)->V[k - 1][l - 1];
		(*Lptr)->V[k - 1][l - 1] = SUBI((*Lptr)->V[k - 1][l - 1], X);
		FREEMPI(Tmp);
		FREEMPI(X);
		FREEMPI(R);
	}
	for (t = 0; t < m; t++)
	{
		if (!EQZEROI((*Bptr)->V[k - 1][t]))
		{
			flag = 0;
			break;
		}
	}
	if (flag)
	{
		TmpMATI = *Bptr;
		*Bptr = DELETE_ROWI(k, *Bptr);
		FREEMATI(TmpMATI);
		for (j = k - 1; j < n - i - 1; j++) {
			*Eptr = SWAP_ROWSI1(j + i, j + i + 1, *Eptr);
			/*              printf("swapping rows %u and %u\n", j+i+1,j+i+2);*/
		}
	}
	FREEMPI(Y);
	return (flag);
}

void STEP8(USI k, MPMATI **B1ptr, MPMATI **Lptr, MPMATI **Eptr, USI i)
{
	MPI *T;

	*B1ptr = SWAP_ROWSI1(k - 2, k - 1, *B1ptr);
	*Eptr = SWAP_ROWSI1(k + i - 2, k + i - 1, *Eptr);
	T = COPYI((*Lptr)->V[k - 1][k - 2]);
	*Lptr = SWAP_ROWSI1(k - 2, k - 1, *Lptr);
	FREEMPI((*Lptr)->V[k - 1][k - 2]);
	(*Lptr)->V[k - 1][k - 2] = T;
	FREEMPI((*Lptr)->V[k - 2][k - 2]);
	(*Lptr)->V[k - 2][k - 2] = ZEROI();
	return;
}

void STEP7(USI i, USI k, MPMATI **Lptr, MPI *D[])
{
	MPI *X1, *X2, *X3, *Y1, *Y2, *Tmp;

	X1 = MULTI((*Lptr)->V[i - 1][k - 2], (*Lptr)->V[k - 1][k - 2]);
	Y1 = MULTI((*Lptr)->V[i - 1][k - 1], D[k - 2]);
	Tmp = Y1;
	Y1 = ADDI(Y1, X1);
	FREEMPI(Tmp);
	FREEMPI(X1);
	X2 = MULTI((*Lptr)->V[i - 1][k - 2], D[k]);
	X3 = MINUSI((*Lptr)->V[k - 1][k - 2]);
	Y2 = MULTI((*Lptr)->V[i - 1][k - 1], X3);
	FREEMPI(X3);
	Tmp = Y2;
	Y2 = ADDI(Y2, X2);
	FREEMPI(Tmp);
	FREEMPI(X2);
	FREEMPI((*Lptr)->V[i - 1][k - 2]);
	(*Lptr)->V[i - 1][k - 2] = INT(Y1, D[k - 1]);
	FREEMPI((*Lptr)->V[i - 1][k - 1]);
	(*Lptr)->V[i - 1][k - 1] = INT(Y2, D[k - 1]);
	FREEMPI(Y1);
	FREEMPI(Y2);
	return;
}

void BASIS_UPDATE(USI i, USI m, MPMATI **Cptr, MPMATI **Lptr, MPMATI *B1ptr, MPI *D[])
{
	unsigned int j, k, t;
	MPI *X, *Tmp, *X1, *X2, *H;

	for (k = 1; k <= m; k++)
	{
		FREEMPI((*Cptr)->V[i - 1][k - 1]);
		(*Cptr)->V[i - 1][k - 1] = COPYI(B1ptr->V[i - 1][k - 1]);
	}
	for (j = 1; j < i; j++)
	{
		X = ZEROI();
		for (t = 0; t < m; t++)
		{
			if (((*Cptr)->V[j - 1][t])->S != 0 && (B1ptr->V[i - 1][t])->S != 0)
			{
				H = MULTI((*Cptr)->V[j - 1][t], B1ptr->V[i - 1][t]);
				Tmp = X;
				X = ADDI(X, H);
				FREEMPI(H);
				FREEMPI(Tmp);
			}
		}
		FREEMPI((*Lptr)->V[i - 1][j - 1]);
		(*Lptr)->V[i - 1][j - 1] = X;
		for (t = 0; t < m; t++)
		{
			if (((*Cptr)->V[i - 1][t])->S == 0)
				X1 = ZEROI();
			else
				X1 = MULTI((*Cptr)->V[i - 1][t], D[j]);
			if (((*Cptr)->V[j - 1][t])->S != 0 && ((*Lptr)->V[i - 1][j - 1])->S != 0)
				X2 = MULTI((*Cptr)->V[j - 1][t], (*Lptr)->V[i - 1][j - 1]);
			else
				X2 = ZEROI();
			Tmp = X1;
			X1 = SUBI(X1, X2);
			FREEMPI(Tmp);
			FREEMPI(X2);
			FREEMPI((*Cptr)->V[i - 1][t]);
			(*Cptr)->V[i - 1][t] = INT(X1, D[j - 1]);
			FREEMPI(X1);
		}
	}
	return;
}

void CSWAP_UPDATE(USI k, USI m, MPI *S, MPMATI **Cptr, MPI *D[])
{
	unsigned int t;
	MPI *Tmp1, *Tmp2, *Tmp3, *Tmp4;

	for (t = 0; t < m; t++)
	{
		Tmp1 = MULTI((*Cptr)->V[k - 1][t], D[k - 2]);
		Tmp2 = MULTI((*Cptr)->V[k - 2][t], S);
		Tmp3 = ADDI(Tmp1, Tmp2);
		FREEMPI(Tmp1);
		FREEMPI(Tmp2);

		Tmp1 = MULTI((*Cptr)->V[k - 2][t], D[k]);
		Tmp2 = MULTI((*Cptr)->V[k - 1][t], S);
		Tmp4 = SUBI(Tmp1, Tmp2);
		FREEMPI(Tmp1);
		FREEMPI(Tmp2);
		FREEMPI((*Cptr)->V[k - 2][t]);
		FREEMPI((*Cptr)->V[k - 1][t]);
		(*Cptr)->V[k - 2][t] = INT(Tmp3, D[k - 1]);
		(*Cptr)->V[k - 1][t] = INT(Tmp4, D[k - 1]);
		FREEMPI(Tmp3);
		FREEMPI(Tmp4);
	}
	return;
}

MPMATI *BASIS_REDUCTION0(MPMATI *Bptr, USI m1, USI n1)
/*
* Input: *Bptr, a matrix of MPI's, whose first row is not zero.
* Output: a pointer to an MPMATI whose rows form a reduced basis for
* the lattice spanned by the rows of *Bptr. This basis is reduced in the
* sense of the paper "Factoring polynomials with rational coefficients" by
* A. K. Lenstra, H. W. Lenstra and L. Lovasz, Math. Ann. 261, 515-534 (1982)
* using the modified version in "Solving exponential Diophantine equations
* using lattice basis reduction algorithms" by B. M. M. De Weger, J. No. Theory
* 26, 325-367 (1987). No change of basis matrix is returned.
* De Weger's algorithm has been changed to cater for arbitrary matrices. The
* the rows are now in general linearly dependent.
* We use the fact that the Gram Schmidt process detects the first row
* which is a linear combination of the preceding rows. We employ a modification
* of the LLL algorithm outlined by M. Pohst in J. Symbolic Computation (1987)4,
* 123-127.  We call this the MLLL algorithm.
* If we are using this algorithm to find small multipliers for the extended
* gcd problem, GCDFLAG is set in EXTGCD() and gcdflag is set below.
* m1 / n1 is usually taken to be 3 / 4.
*/
{
	unsigned int i, k, l, n, m, t, flag = 0, Flag = 0, gcdflag = 0;
	unsigned int flagg, beta, K1 = 0, tau = 2, sigma = 0, rho;
	MPI **D, *X, *Y, *Z, *H, *Tmp, *R, *M1, *N1;
	MPMATI *C, *L, *B1ptr;
	unsigned int norig;

	Z = NULL;
	m = Bptr->C;
	n = Bptr->R;
	norig = n;
	B1ptr = COPYMATI(Bptr);
	D = (MPI **)mmalloc((1 + n) * sizeof(MPI *));
	D[0] = ONEI();
	for (i = 1; i <= n; i++)
		D[i] = ZEROI();
	C = ZEROMNI(n, m);
	L = ZEROMNI(n, n);

found:
	n = B1ptr->R;
	i = (K1 == 0) ? 1 : K1;
	/* K1 = no. of consecutive rows of *B1ptr that don't need updating
	for the Gram Schmidt process. */
	while (i <= B1ptr->R)
	{
		BASIS_UPDATE(i, m, &C, &L, B1ptr, D);
		flag = 1;
		for (t = 0; t < m; t++)
		{
			if (!EQZEROI(C->V[i - 1][t]))
			{
				flag = 0;
				break;
			}
		}
		if (flag)
			break;
		X = ZEROI();
		for (t = 0; t < m; t++)
		{
			H = MULTI(C->V[i - 1][t], C->V[i - 1][t]);
			Tmp = X;
			X = ADDI(X, H);
			FREEMPI(Tmp);
			FREEMPI(H);
		}
		FREEMPI(D[i]);
		D[i] = INT(X, D[i - 1]);
		FREEMPI(X);
		i++;
		if (MLLLVERBOSE)
			printf("i = %u\n", i);
	}
	beta = (flag) ? i : i - 1;
	rho = K1 = i - 1;
	if (MLLLVERBOSE)
		printf("BASIS0 completed updating the basis\n");
	/* Here K1 = no. of LI rows in *B1ptr found by Gram Schmidt process.
	flag = 0 means all the rho = beta rows of *B1ptr are LI;
	flag = 1 means that the first rho = beta - 1 rows of *B1ptr are LI, but the
	beta-th row is a LC of the preceding rows. So beta = number of rows of *B1ptr
	currently being examined by the LLL algorithm. */
	k = tau;
	M1 = CHANGE(m1);
	N1 = CHANGE(n1);
	while (k <= beta)
	{
		if (MLLLVERBOSE)
			printf("beta - k = %u\n", beta - k);
		l = k - 1;
		Flag = STEP40(k, l, &L, &B1ptr, D);
		if (k >= norig && GCDFLAG)
		{
			gcdflag = 1;
			goto FOUND;
		}
		if (Flag)/* STEP 9 of POHST. */
		{
			sigma++;
			if (MLLLVERBOSE)
				printf("relation vector number %u found\n", sigma);
			tau = k++;
			goto found;
		}
		X = MULTI(D[k - 2], D[k]);
		Y = MULTI(D[k - 1], D[k - 1]);
		Tmp = Y;
		Y = MULT_I(Y, m1);
		FREEMPI(Tmp);
		R = MULTI(L->V[k - 1][k - 2], L->V[k - 1][k - 2]);
		Z = ADD0I(X, R);
		Tmp = Z;
		Z = MULT_I(Z, n1);
		FREEMPI(Tmp);
		if (RSV(Y, Z) == 1)/*& STEP 5 of POHST. */
		{
			flagg = 0;
			if (EQZEROI(D[k]) && EQZEROI(R))
			{/* CASE B=0 of STEP 7 of POHST. */
				FREEMPI(D[k - 1]);
				D[k - 1] = ZEROI();
				STEP80(k, &B1ptr, &L);
				if (k - 1 < K1)
					K1 = k - 1;
				/* The swap may have changed 2nd last row */
				/* of *B1ptr. */
				for (t = 0; t < m; t++)
				{
					FREEMPI(C->V[k - 2][t]);
					C->V[k - 2][t] = ZEROI();
				}
				beta--;
				flagg = 1;
				if (k > 2)
					k--;
				FREEMPI(X);
				FREEMPI(Y);
				FREEMPI(R);
				FREEMPI(Z);
				continue;
			}
			if (flagg == 0)
			{
				for (i = k + 1; i <= beta; i++)
					STEP7(i, k, &L, D);
			}
			STEP80(k, &B1ptr, &L);
			if (k - 2 < K1)
				K1 = k - 2;
			/* swap will change last two rows of *B1ptr. */
			FREEMPI(R);
			FREEMPI(Y);
			Y = MULTI(L->V[k - 1][k - 2], L->V[k - 1][k - 2]);
			Tmp = Y;
			Y = ADD0I(Y, X);
			FREEMPI(X);
			FREEMPI(Tmp);
			Tmp = D[k - 1];
			D[k - 1] = INT0(Y, D[k - 1]);
			FREEMPI(Tmp);
			FREEMPI(Y);
			if (k > 2)
				k--;
		}
		else
		{ /* STEP 6 of POHST. */
			FREEMPI(R);
			FREEMPI(X);
			FREEMPI(Y);
		FOUND:
			for (l = k - 2; l >= 1; l--)
			{
				Flag = STEP40(k, l, &L, &B1ptr, D);
				if (Flag)
				{
					FREEMPI(Z);
					sigma++;
					if (MLLLVERBOSE)
						printf("relation vector number %u found\n", sigma);
					tau = k++;
					goto found; /* STEP 9 of POHST. */
				}
			}
			k++;
		}
		if (!gcdflag)
			FREEMPI(Z);
		else
			break;
	}
	FREEMPI(M1);
	FREEMPI(N1);
	FREEMATI(C);
	/*
	printf("L = \n");
	PRINTMATI(0,L->R-1,0,L->C-1,L);
	for (i = 0; i <= norig; i++)
	{
	printf("D[%u] = ", i);PRINTI(D[i]);printf(", ");
	}
	printf("\n");
	*/
	FREEMATI(L);
	for (i = 0; i <= norig; i++)
		FREEMPI(D[i]);
	ffree((char *)D, (1 + norig) * sizeof(MPI *));
	if (MLLLVERBOSE)
	{
		printf("number of basis vectors found = %u ;\n", rho);
		printf("number of relation vectors found = %u .\n", sigma);
	}
	return (B1ptr);
}

unsigned int STEP40(k, l, Lptr, Bptr, D)
/*
* updates *Lptr, *Bptr.
* returns 1 if row k of *Bptr becomes zero, returns zero otherwise.
*/
unsigned int k, l;
MPI *D[];
MPMATI **Lptr, **Bptr;
{
	unsigned int j, flag = 1, t, m;
	MPI *X, *Y, *R, *Tmp;
	MPMATI *TmpMATI;

	m = (*Bptr)->C;
	Y = MULT_I((*Lptr)->V[k - 1][l - 1], 2);
	if (RSV(Y, D[l]) == 1)
	{
		R = NEAREST_INTI((*Lptr)->V[k - 1][l - 1], D[l]);
		X = MINUSI(R);
		TmpMATI = *Bptr;
		*Bptr = ADD_MULT_ROWI(l - 1, k - 1, X, *Bptr);
		if (MLLLVERBOSE)
		{
			printf("Row %u -> Row %u + ", k, k); PRINTI(X); printf(" x Row %u\n", l);
			PRINTMATI(0, (*Bptr)->R - 1, 0, (*Bptr)->C - 1, *Bptr);
			GetReturn();
		}
		FREEMATI(TmpMATI);
		/*
		MAXI = UPDATEMAXI(MAXI, *Bptr);
		*/
		FREEMPI(X);
		for (j = 1; j < l; j++)
		{
			X = MULTI((*Lptr)->V[l - 1][j - 1], R);
			Tmp = (*Lptr)->V[k - 1][j - 1];
			(*Lptr)->V[k - 1][j - 1] = SUBI((*Lptr)->V[k - 1][j - 1], X);
			FREEMPI(Tmp);
			FREEMPI(X);
		}
		X = MULTI(D[l], R);
		Tmp = (*Lptr)->V[k - 1][l - 1];
		(*Lptr)->V[k - 1][l - 1] = SUBI((*Lptr)->V[k - 1][l - 1], X);
		FREEMPI(Tmp);
		FREEMPI(X);
		FREEMPI(R);
	}
	for (t = 0; t < m; t++)
	{
		if (!EQZEROI((*Bptr)->V[k - 1][t]))
		{
			flag = 0;
			break;
		}
	}
	if (flag)
	{
		TmpMATI = *Bptr;
		*Bptr = DELETE_ROWI(k, *Bptr);
		FREEMATI(TmpMATI);
	}
	FREEMPI(Y);
	return (flag);
}

unsigned int STEP400(k, l, Lptr, Bptr, D, gcdflag, norig)
/*
* updates *Lptr, *Bptr.
* returns 1 if row k of *Bptr becomes zero, returns zero otherwise.
*/
unsigned int k, l, gcdflag, norig;
MPI *D[];
MPMATI **Lptr, **Bptr;
{
	unsigned int j, flag = 1, t, m;
	MPI *X, *Y, *R, *Tmp;
	MPMATI *TmpMATI;

	m = (*Bptr)->C;
	Y = MULT_I((*Lptr)->V[k - 1][l - 1], 2);
	if (RSV(Y, D[l]) == 1)
	{
		R = NEAREST_INTI((*Lptr)->V[k - 1][l - 1], D[l]);
		X = MINUSI(R);
		TmpMATI = *Bptr;
		*Bptr = ADD_MULT_ROWI(l - 1, k - 1, X, *Bptr);
		FREEMATI(TmpMATI);
		/*
		MAXI = UPDATEMAXI(MAXI, *Bptr);
		*/
		FREEMPI(X);
		for (j = 1; j < l; j++)
		{
			X = MULTI((*Lptr)->V[l - 1][j - 1], R);
			Tmp = (*Lptr)->V[k - 1][j - 1];
			(*Lptr)->V[k - 1][j - 1] = SUBI((*Lptr)->V[k - 1][j - 1], X);
			FREEMPI(Tmp);
			FREEMPI(X);
		}
		X = MULTI(D[l], R);
		Tmp = (*Lptr)->V[k - 1][l - 1];
		(*Lptr)->V[k - 1][l - 1] = SUBI((*Lptr)->V[k - 1][l - 1], X);
		FREEMPI(Tmp);
		FREEMPI(X);
		FREEMPI(R);
	}
	for (t = 0; t < m; t++)
	{
		if (!EQZEROI((*Bptr)->V[k - 1][t]))
		{
			flag = 0;
			break;
		}
	}
	if (flag)
	{
		TmpMATI = *Bptr;
		*Bptr = DELETE_ROWI(k, *Bptr);
		FREEMATI(TmpMATI);
	}
	FREEMPI(Y);
	return (flag);
}

void STEP80(k, B1ptr, Lptr)
MPMATI **B1ptr, **Lptr;
unsigned int k;
{
	MPI *T;

	*B1ptr = SWAP_ROWSI1(k - 2, k - 1, *B1ptr);
	if (MLLLVERBOSE)
	{
		printf("Swapping Rows %u and %u\n", k - 1, k);
		PRINTMATI(0, (*B1ptr)->R - 1, 0, (*B1ptr)->C - 1, *B1ptr);
		GetReturn();
	}
	T = COPYI((*Lptr)->V[k - 1][k - 2]);
	*Lptr = SWAP_ROWSI1(k - 2, k - 1, *Lptr);
	FREEMPI((*Lptr)->V[k - 1][k - 2]);
	(*Lptr)->V[k - 1][k - 2] = T;
	FREEMPI((*Lptr)->V[k - 2][k - 2]);
	(*Lptr)->V[k - 2][k - 2] = ZEROI();
	return;
}

MPMATI *EXTGCD(MPMATI *Dptr, MPI **Aptr, MPMATI **Q, USI m1, USI n1)
/*
* Input: an n x 1 vector of MPI's.
* Output: *Aptr = gcd of the vector of MPI's. Also we return a small set of
* multipliers using the LLL method of Havas and Matthews.
* parameters m1 and n1 were put in at the request of George Havas on 11/7/94.
* Normally m1/n1 = 3/4.
* matrix Q of the LLL extended gcd paper of Havas-Matthews is returned.
*/
{
	MPMATI *Temp, *P, *SM, *R;
	USI i, m, n, *alpha, nz;
	MPI **Tptr;

	n = Dptr->R;
	m = n - 1;
	alpha = KB_ROWP(Dptr, &R, &nz);
	Temp = HERMITE1P(Dptr, R, &P, alpha, nz);
	FREEMATI(R);
	ffree((char *)alpha, (Dptr->C) * sizeof(USI));
	*Aptr = COPYI(Temp->V[0][0]);
	FREEMATI(Temp);
	Tptr = P->V[0];
	for (i = 0; i < m; i++)
		P->V[i] = P->V[i + 1];
	P->V[n - 1] = Tptr;
	GCDFLAG = 1;
	*Q = BASIS_REDUCTION0(P, m1, n1);
	SM = BUILDMATI(1, n);
	for (i = 0; i < n; i++)
		SM->V[0][i] = COPYI((*Q)->V[m][i]);
	FREEMATI(P);
	GCDFLAG = 0;
	return (SM);
}

MPI *LENGTHSQRI(MPMATI *Mptr, USI i)
/*
* Returns the square of the length of row i of matrix *Mptr.
*/
{
	MPI *SUM, *T, *T1;
	USI j;

	SUM = ZEROI();
	for (j = 0; j < Mptr->C; j++)
	{
		T = SUM;
		T1 = MULTI(Mptr->V[i][j], Mptr->V[i][j]);
		SUM = ADD0I(SUM, T1);
		FREEMPI(T);
		FREEMPI(T1);
	}
	return (SUM);
}

MPI *LENGTHSQCI(MPMATI *Mptr, USI j)
/*
* Returns the square of the length of column j of matrix *Mptr.
*/
{
	MPI *SUM, *T, *T1;
	USI i;

	SUM = ZEROI();
	for (i = 0; i < Mptr->R; i++)
	{
		T = SUM;
		T1 = MULTI(Mptr->V[i][j], Mptr->V[i][j]);
		SUM = ADD0I(SUM, T1);
		FREEMPI(T);
		FREEMPI(T1);
	}
	return (SUM);
}

MPMATI *BASIS_REDUCTION00(MPMATI *Bptr, USI m1, USI n1, USI norig)
/*
* Input: *Bptr, a matrix of MPI's, whose first row is not zero.
* Output: a pointer to an MPMATI whose rows form a reduced basis for
* the lattice spanned by the rows of *Bptr. This basis is reduced in the
* sense of the paper "Factoring polynomials with rational coefficients" by
* A. K. Lenstra, H. W. Lenstra and L. Lovasz, Math. Ann. 261, 515-534 (1982)
* using the modified version in "Solving exponential Diophantine equations
* using lattice basis reduction algorithms" by B. M. M. De Weger, J. No. Theory
* 26, 325-367 (1987). No change of basis matrix is returned.
* De Weger's algorithm has been changed to cater for arbitrary matrices. The
* the rows are now in general linearly dependent.
* We use the fact that the Gram Schmidt process detects the first row
* which is a linear combination of the preceding rows. We employ a modification
* of the LLL algorithm outlined by M. Pohst in J. Symbolic Computation (1987)4,
* 123-127.  We call this the MLLL algorithm.
* If we are using this algorithm to find small multipliers for the extended
* gcd problem, GCDFLAG is set in IMPROVEP() and gcdflag is set below.
* m1 / n1 is usually taken to be 3 / 4.
* For use in IMPROVEP().
*/
{
	unsigned int i, k, n, m, t, flag = 0, Flag = 0, gcdflag = 0;
	unsigned int flagg, beta, K1 = 0, tau = 2, sigma = 0, rho, norigg;
	unsigned int K, j, REPEAT, iterate;
	int l, tt;
	MPI **D, *X, *Y, *Z, *H, *Tmp, *R, *M1, *N1;
	MPI *t1, *t2, *t3, *t4;
	MPMATI *C, *L, *B1ptr, *temp1;

	Z = NULL;
	m = Bptr->C;
	norigg = n = Bptr->R;
	B1ptr = COPYMATI(Bptr);
	D = (MPI **)mmalloc((1 + n) * sizeof(MPI *));
	D[0] = ONEI();
	for (i = 1; i <= n; i++)
		D[i] = ZEROI();
	C = ZEROMNI(n, m);
	L = ZEROMNI(n, n);

found:
	n = B1ptr->R;
	i = (K1 == 0) ? 1 : K1;
	/* K1 = no. of consecutive rows of *B1ptr that don't need updating
	for the Gram Schmidt process. */
	while (i <= B1ptr->R)
	{
		BASIS_UPDATE(i, m, &C, &L, B1ptr, D);
		flag = 1;
		for (t = 0; t < m; t++)
		{
			if (!EQZEROI(C->V[i - 1][t]))
			{
				flag = 0;
				break;
			}
		}
		if (flag)
			break;
		X = ZEROI();
		for (t = 0; t < m; t++)
		{
			H = MULTI(C->V[i - 1][t], C->V[i - 1][t]);
			Tmp = X;
			X = ADDI(X, H);
			FREEMPI(Tmp);
			FREEMPI(H);
		}
		FREEMPI(D[i]);
		D[i] = INT(X, D[i - 1]);
		FREEMPI(X);
		i++;
		if (MLLLVERBOSE) {
			printf("i = %u\n", i);
		}
	}
	/*
	strcpy(buff, "L.out");
	outfile = fopen(buff, "w");
	FPRINTMATI(outfile,0,L->R-1,0,L->C-1,L);
	fclose(outfile);
	for (t = 0;t < B1ptr->R; t++){
	printf("D[%u] = ", t);PRINTI(D[t]);printf("\n");
	}
	GetReturn();
	*/
	beta = (flag) ? i : i - 1;
	rho = K1 = i - 1;
	if (MLLLVERBOSE)
		printf("completed updating the basis\n");
	/* Here K1 = no. of LI rows in *B1ptr found by Gram Schmidt process.
	flag = 0 means all the rho = beta rows of *B1ptr are LI;
	flag = 1 means that the first rho = beta - 1 rows of *B1ptr are LI, but the
	beta-th row is a LC of the preceding rows. So beta = number of rows of *B1ptr
	currently being examined by the LLL algorithm. */
	k = tau;
	M1 = CHANGE(m1);
	N1 = CHANGE(n1);
	while (k <= beta)
	{
		if (MLLLVERBOSE)
			printf("beta - k = %u\n", beta - k);
		if (gcdflag)
			l = norig;
		else
			l = k - 1;
		Flag = STEP40(k, l, &L, &B1ptr, D);
		if (k > norig && GCDFLAG)
		{
			gcdflag = 1;
			if (MLLLVERBOSE)
				printf("improving row %u\n", k);
			goto FOUND;
		}
		if (Flag)/* STEP 9 of POHST. */
		{
			sigma++;
			if (MLLLVERBOSE)
				printf("relation vector number %u found\n", sigma);
			tau = k++;
			goto found;
		}
		X = MULTI(D[k - 2], D[k]);
		Y = MULTI(D[k - 1], D[k - 1]);
		Tmp = Y;
		Y = MULT_I(Y, m1);
		FREEMPI(Tmp);
		R = MULTI(L->V[k - 1][k - 2], L->V[k - 1][k - 2]);
		Z = ADD0I(X, R);
		Tmp = Z;
		Z = MULT_I(Z, n1);
		FREEMPI(Tmp);
		if (RSV(Y, Z) == 1)/*& STEP 5 of POHST. */
		{
			flagg = 0;
			if (EQZEROI(D[k]) && EQZEROI(R))
			{/* CASE B=0 of STEP 7 of POHST. */
				FREEMPI(D[k - 1]);
				D[k - 1] = ZEROI();
				STEP80(k, &B1ptr, &L);
				if (k - 1 < K1)
					K1 = k - 1;
				/* The swap may have changed 2nd last row */
				/* of *B1ptr. */
				for (t = 0; t < m; t++)
				{
					FREEMPI(C->V[k - 2][t]);
					C->V[k - 2][t] = ZEROI();
				}
				beta--;
				flagg = 1;
				if (k > 2)
					k--;
				FREEMPI(X);
				FREEMPI(Y);
				FREEMPI(R);
				FREEMPI(Z);
				continue;
			}
			if (flagg == 0)
			{
				for (i = k + 1; i <= beta; i++)
					STEP7(i, k, &L, D);
			}
			STEP80(k, &B1ptr, &L);
			if (k - 2 < K1)
				K1 = k - 2;
			/* swap will change last two rows of *B1ptr. */
			FREEMPI(R);
			FREEMPI(Y);
			Y = MULTI(L->V[k - 1][k - 2], L->V[k - 1][k - 2]);
			Tmp = Y;
			Y = ADD0I(Y, X);
			FREEMPI(X);
			FREEMPI(Tmp);
			Tmp = D[k - 1];
			D[k - 1] = INT0(Y, D[k - 1]);
			FREEMPI(Tmp);
			FREEMPI(Y);
			if (k > 2)
				k--;
		}
		else
		{ /* STEP 6 of POHST. */
			FREEMPI(R);
			FREEMPI(X);
			FREEMPI(Y);
		FOUND:
			if (gcdflag)
			{
				if (MLLLVERBOSE)
					printf("gcdflag set\n");
				K = norig + 1;
			}
			else
				K = k;
			for (l = K - 2; l >= 1; l--)
			{
				Flag = STEP40(k, l, &L, &B1ptr, D);
				if (Flag)
				{
					FREEMPI(Z);
					sigma++;
					if (MLLLVERBOSE)
						printf("relation vector number %u found\n", sigma);
					tau = k++;
					goto found; /* STEP 9 of POHST. */
				}
			}
			k++;
		}
		if (!gcdflag)
			FREEMPI(Z);
	}
	FREEMPI(M1);
	FREEMPI(N1);
	FREEMATI(C);
	FREEMATI(L);
	for (i = 0; i <= norigg; i++)
		FREEMPI(D[i]);
	ffree((char *)D, (1 + norigg) * sizeof(MPI *));
	/* Now to improve the last n - norig rows of *B1ptr, using Gauss
	lattice reduction - pointed out by George Havas in Sims' book - 8/11/94: */


	/*
	strcpy(buff, "progress");
	outfile = fopen(buff, "a+");
	FPRINTMATI(outfile,0,B1ptr->R-1,0,B1ptr->C-1,B1ptr);
	*/


	REPEAT = 1;
	temp1 = BUILDMATI(1, B1ptr->C);
	for (k = 0; k < B1ptr->C; k++)
		temp1->V[0][k] = ZEROI();
	iterate = 0;
	QSORTMATI(B1ptr, 0, norig - 1);
	while (REPEAT)
	{
		REPEAT = 0;
		for (j = 1; j < norig; j++)
		{
			t4 = DOTRI(B1ptr, j, j);
			for (i = 0; i < j; i++)
			{
				t1 = DOTRI(B1ptr, j, i);
				t2 = DOTRI(B1ptr, i, i);
				t3 = NEAREST_INTI(t1, t2);
				FREEMPI(t1);
				FREEMPI(t2);
				if (t3->S == 0)
				{
					FREEMPI(t3);
					continue;
				}
				t1 = t3;
				t3 = MINUSI(t3);
				FREEMPI(t1);
				for (k = 0; k < B1ptr->C; k++)
				{
					FREEMPI(temp1->V[0][k]);
					X = MULTI(t3, B1ptr->V[i][k]);
					temp1->V[0][k] = ADDI(X, B1ptr->V[j][k]);
					FREEMPI(X);
				}
				t1 = DOTRI(temp1, 0, 0);
				tt = RSV(t1, t4);
				tt = RSV(t1, t4);
				if (tt == -1)
				{
					ADD_MULT_ROWI0(i, j, t3, B1ptr);
					if (MLLLVERBOSE)
					{
						printf("Gauss-improving nullspace row %u: length squared was = ", j + 1);
						PRINTI(t4);
						printf("\n");
						printf("lengthsquared now = ");
						PRINTI(t1);
						printf("\n");
					}
					/*
					fprintf(outfile, "Gauss-improving nullspace row %u: length squared was = ", j + 1);
					FPRINTI(outfile, t4);
					fprintf(outfile, "\n");
					fprintf(outfile, "lengthsquared now = ");
					FPRINTI(outfile, t1);
					fprintf(outfile, "\n");
					fflush(outfile);
					*/
					FREEMPI(t4);
					t4 = t1;
					REPEAT = 1;
				}
				else
					FREEMPI(t1);
				FREEMPI(t3);
			}
			FREEMPI(t4);
		}
		iterate++;
		QSORTMATI(B1ptr, 0, norig - 1);
	}
	REPEAT = 1;
	iterate = 0;
	while (REPEAT)
	{
		REPEAT = 0;
		for (j = norig; j < B1ptr->R; j++)
		{
			t4 = DOTRI(B1ptr, j, j);
			for (i = 0; i < norig; i++)
			{
				t1 = DOTRI(B1ptr, j, i);
				t2 = DOTRI(B1ptr, i, i);
				t3 = NEAREST_INTI(t1, t2);
				FREEMPI(t1);
				FREEMPI(t2);
				if (t3->S == 0)
				{
					FREEMPI(t3);
					continue;
				}
				t1 = t3;
				t3 = MINUSI(t3);
				FREEMPI(t1);
				for (k = 0; k < B1ptr->C; k++)
				{
					FREEMPI(temp1->V[0][k]);
					X = MULTI(t3, B1ptr->V[i][k]);
					temp1->V[0][k] = ADDI(X, B1ptr->V[j][k]);
					FREEMPI(X);
				}
				t1 = DOTRI(temp1, 0, 0);
				tt = RSV(t1, t4);
				tt = RSV(t1, t4);
				if (tt == -1)
				{
					ADD_MULT_ROWI0(i, j, t3, B1ptr);
					if (MLLLVERBOSE)
					{
						printf("Gauss-improving row %u: length squared was = ", j + 1);
						PRINTI(t4);
						printf("\n");
						printf("lengthsquared now = ");
						PRINTI(t1);
						printf("\n");
					}
					/*
					fprintf(outfile, "Gauss-improving row %u: length squared was = ", j + 1);
					FPRINTI(outfile, t4);
					fprintf(outfile, "\n");
					fprintf(outfile, "lengthsquared now = ");
					FPRINTI(outfile, t1);
					fprintf(outfile, "\n");
					fflush(outfile);
					*/
					FREEMPI(t4);
					t4 = t1;
					REPEAT = 1;
				}
				else
					FREEMPI(t1);
				FREEMPI(t3);
			}
			FREEMPI(t4);
		}
		iterate++;
	}
	FREEMATI(temp1);
	/*
	fclose(outfile);
	*/
	return (B1ptr);
}

MPMATI *BASIS_REDUCTION000(MPMATI *Bptr, USI m1, USI n1, MPI *N)
/*
* Input: *Bptr, a matrix of MPI's, whose first row is not zero.
* Output: a pointer to an MPMATI whose rows form a reduced basis for
* the lattice spanned by the rows of *Bptr. This basis is reduced in the
* sense of the paper "Factoring polynomials with rational coefficients" by
* A. K. Lenstra, H. W. Lenstra and L. Lovasz, Math. Ann. 261, 515-534 (1982)
* using the modified version in "Solving exponential Diophantine equations
* using lattice basis reduction algorithms" by B. M. M. De Weger, J. No. Theory
* 26, 325-367 (1987). No change of basis matrix is returned.
* De Weger's algorithm has been changed to cater for arbitrary matrices. The
* the rows are now in general linearly dependent.
* We use the fact that the Gram Schmidt process detects the first row
* which is a linear combination of the preceding rows. We employ a modification
* of the LLL algorithm outlined by M. Pohst in J. Symbolic Computation (1987)4,
* 123-127.  We call this the MLLL algorithm.
* If we are using this algorithm to find small multipliers for the extended
* gcd problem, GCDFLAG is set in EXTGCD() and gcdflag is set below.
* m1 / n1 is usually taken to be 3 / 4.
*/
{
	unsigned int i, k, l, n, m, t, flag = 0, Flag = 0, gcdflag = 0;
	unsigned int flagg, beta, K1 = 0, tau = 2, sigma = 0, rho;
	MPI **D, *X, *Y, *Z, *H, *Tmp, *R, *M1, *N1;
	MPMATI *C, *L, *B1ptr;
	unsigned int norig;

	Z = NULL;
	m = Bptr->C;
	n = Bptr->R;
	norig = n;
	B1ptr = COPYMATI(Bptr);
	D = (MPI **)mmalloc((1 + n) * sizeof(MPI *));
	D[0] = ONEI();
	for (i = 1; i <= n; i++)
		D[i] = ZEROI();
	C = ZEROMNI(n, m);
	L = ZEROMNI(n, n);

found:
	n = B1ptr->R;
	i = (K1 == 0) ? 1 : K1;
	/* K1 = no. of consecutive rows of *B1ptr that don't need updating
	for the Gram Schmidt process. */
	while (i <= B1ptr->R)
	{
		BASIS_UPDATE(i, m, &C, &L, B1ptr, D);
		flag = 1;
		for (t = 0; t < m; t++)
		{
			if (!EQZEROI(C->V[i - 1][t]))
			{
				flag = 0;
				break;
			}
		}
		if (flag)
			break;
		X = ZEROI();
		for (t = 0; t < m; t++)
		{
			H = MULTI(C->V[i - 1][t], C->V[i - 1][t]);
			Tmp = X;
			X = ADDI(X, H);
			FREEMPI(Tmp);
			FREEMPI(H);
		}
		FREEMPI(D[i]);
		D[i] = INT(X, D[i - 1]);
		FREEMPI(X);
		i++;
		if (MLLLVERBOSE)
			printf("i = %u\n", i);
	}
	beta = (flag) ? i : i - 1;
	rho = K1 = i - 1;
	if (MLLLVERBOSE)
		printf("completed updating the basis\n");
	/* Here K1 = no. of LI rows in *B1ptr found by Gram Schmidt process.
	flag = 0 means all the rho = beta rows of *B1ptr are LI;
	flag = 1 means that the first rho = beta - 1 rows of *B1ptr are LI, but the
	beta-th row is a LC of the preceding rows. So beta = number of rows of *B1ptr
	currently being examined by the LLL algorithm. */
	k = tau;
	M1 = CHANGE(m1);
	N1 = CHANGE(n1);
	while (k <= beta)
	{
		if (MLLLVERBOSE)
			printf("beta - k = %u\n", beta - k);
		l = k - 1;
		Flag = STEP4000(k, l, &L, &B1ptr, D, N);
		if (k >= norig && GCDFLAG)
		{
			gcdflag = 1;
			goto FOUND;
		}
		if (Flag)/* STEP 9 of POHST. */
		{
			sigma++;
			if (MLLLVERBOSE)
				printf("relation vector number %u found\n", sigma);
			tau = k++;
			goto found;
		}
		X = MULTI(D[k - 2], D[k]);
		Y = MULTI(D[k - 1], D[k - 1]);
		Tmp = Y;
		Y = MULT_I(Y, m1);
		FREEMPI(Tmp);
		R = MULTI(L->V[k - 1][k - 2], L->V[k - 1][k - 2]);
		Z = ADD0I(X, R);
		Tmp = Z;
		Z = MULT_I(Z, n1);
		FREEMPI(Tmp);
		if (RSV(Y, Z) == 1)/*& STEP 5 of POHST. */
		{
			flagg = 0;
			if (EQZEROI(D[k]) && EQZEROI(R))
			{/* CASE B=0 of STEP 7 of POHST. */
				FREEMPI(D[k - 1]);
				D[k - 1] = ZEROI();
				STEP8000(k, &B1ptr, &L, N);
				if (k - 1 < K1)
					K1 = k - 1;
				/* The swap may have changed 2nd last row */
				/* of *B1ptr. */
				for (t = 0; t < m; t++)
				{
					FREEMPI(C->V[k - 2][t]);
					C->V[k - 2][t] = ZEROI();
				}
				beta--;
				flagg = 1;
				if (k > 2)
					k--;
				FREEMPI(X);
				FREEMPI(Y);
				FREEMPI(R);
				FREEMPI(Z);
				continue;
			}
			if (flagg == 0)
			{
				for (i = k + 1; i <= beta; i++)
					STEP7(i, k, &L, D);
			}
			STEP8000(k, &B1ptr, &L, N);
			if (k - 2 < K1)
				K1 = k - 2;
			/* swap will change last two rows of *B1ptr. */
			FREEMPI(R);
			FREEMPI(Y);
			Y = MULTI(L->V[k - 1][k - 2], L->V[k - 1][k - 2]);
			Tmp = Y;
			Y = ADD0I(Y, X);
			FREEMPI(X);
			FREEMPI(Tmp);
			Tmp = D[k - 1];
			D[k - 1] = INT0(Y, D[k - 1]);
			FREEMPI(Tmp);
			FREEMPI(Y);
			if (k > 2)
				k--;
		}
		else
		{ /* STEP 6 of POHST. */
			FREEMPI(R);
			FREEMPI(X);
			FREEMPI(Y);
		FOUND:
			for (l = k - 2; l >= 1; l--)
			{
				Flag = STEP4000(k, l, &L, &B1ptr, D, N);
				if (Flag)
				{
					FREEMPI(Z);
					sigma++;
					if (MLLLVERBOSE)
						printf("relation vector number %u found\n", sigma);
					tau = k++;
					goto found; /* STEP 9 of POHST. */
				}
			}
			k++;
		}
		if (!gcdflag)
			FREEMPI(Z);
		else
			break;
	}
	FREEMPI(M1);
	FREEMPI(N1);
	FREEMATI(C);
	FREEMATI(L);
	for (i = 0; i <= norig; i++)
		FREEMPI(D[i]);
	ffree((char *)D, (1 + norig) * sizeof(MPI *));
	if (MLLLVERBOSE)
	{
		printf("number of basis vectors found = %u ;\n", rho);
		printf("number of relation vectors found = %u .\n", sigma);
	}
	return (B1ptr);
}

unsigned int STEP4000(USI k, USI l, MPMATI **Lptr, MPMATI **Bptr, MPI *D[], MPI *N)
/*
* updates *Lptr, *Bptr.
* returns 1 if row k of *Bptr becomes zero, returns zero otherwise.
*/
{
	unsigned int i, j, flag = 1, t, m, n;
	MPI *X, *Y, *R, *Tmp, *Temp;
	MPMATI *TmpMATI;

	m = (*Bptr)->C;
	Y = MULT_I((*Lptr)->V[k - 1][l - 1], 2);
	if (RSV(Y, D[l]) == 1)
	{
		R = NEAREST_INTI((*Lptr)->V[k - 1][l - 1], D[l]);
		X = MINUSI(R);
		TmpMATI = *Bptr;
		*Bptr = ADD_MULT_ROWI(l - 1, k - 1, X, *Bptr);
		FREEMATI(TmpMATI);
		n = (*Bptr)->R;
		for (i = 0; i < n; i++)
		{
			for (j = n; j < m; j++)
			{
				Temp = POWERI(N, m - j);
				(*Bptr)->V[i][j] = INTI((*Bptr)->V[i][j], Temp);
				FREEMPI(Temp);
			}
		}
		if (HERMITEVERBOSE)
		{
			printf("Row %u -> Row %u + ", k, k); PRINTI(X); printf(" x Row %u\n", l);
			PRINTMATI(0, (*Bptr)->R - 1, 0, (*Bptr)->C - 1, *Bptr);
			GetReturn();
		}
		for (i = 0; i < n; i++)
		{
			for (j = n; j < m; j++)
			{
				Temp = POWERI(N, m - j);
				(*Bptr)->V[i][j] = MULTI((*Bptr)->V[i][j], Temp);
				FREEMPI(Temp);
			}
		}
		/*
		MAXI = UPDATEMAXI(MAXI, *Bptr);
		*/
		FREEMPI(X);
		for (j = 1; j < l; j++)
		{
			X = MULTI((*Lptr)->V[l - 1][j - 1], R);
			Tmp = (*Lptr)->V[k - 1][j - 1];
			(*Lptr)->V[k - 1][j - 1] = SUBI((*Lptr)->V[k - 1][j - 1], X);
			FREEMPI(Tmp);
			FREEMPI(X);
		}
		X = MULTI(D[l], R);
		Tmp = (*Lptr)->V[k - 1][l - 1];
		(*Lptr)->V[k - 1][l - 1] = SUBI((*Lptr)->V[k - 1][l - 1], X);

		FREEMPI(Tmp);
		FREEMPI(X);
		FREEMPI(R);
	}
	for (t = 0; t < m; t++)
	{
		if (!EQZEROI((*Bptr)->V[k - 1][t]))
		{
			flag = 0;
			break;
		}
	}
	if (flag)
	{
		TmpMATI = *Bptr;
		*Bptr = DELETE_ROWI(k, *Bptr);
		FREEMATI(TmpMATI);
	}
	FREEMPI(Y);
	return (flag);
}
void STEP8000(USI k, MPMATI **B1ptr, MPMATI **Lptr, MPI *N)
{
	MPI *T, *Temp;
	USI i, j, n, m;

	*B1ptr = SWAP_ROWSI1(k - 2, k - 1, *B1ptr);
	T = COPYI((*Lptr)->V[k - 1][k - 2]);
	*Lptr = SWAP_ROWSI1(k - 2, k - 1, *Lptr);
	FREEMPI((*Lptr)->V[k - 1][k - 2]);
	(*Lptr)->V[k - 1][k - 2] = T;
	FREEMPI((*Lptr)->V[k - 2][k - 2]);
	(*Lptr)->V[k - 2][k - 2] = ZEROI();

	n = (*B1ptr)->R;
	m = (*B1ptr)->C;
	for (i = 0; i < n; i++)
	{
		for (j = n; j < m; j++)
		{
			Temp = POWERI(N, m - j);
			(*B1ptr)->V[i][j] = INTI((*B1ptr)->V[i][j], Temp);
			FREEMPI(Temp);
		}
	}
	if (HERMITEVERBOSE)
	{
		printf("Swapping Rows %u and %u\n", k - 1, k);
		PRINTMATI(0, (*B1ptr)->R - 1, 0, (*B1ptr)->C - 1, *B1ptr);
		GetReturn();
	}
	for (i = 0; i < n; i++)
	{
		for (j = n; j < m; j++)
		{
			Temp = POWERI(N, m - j);
			(*B1ptr)->V[i][j] = MULTI((*B1ptr)->V[i][j], Temp);
			FREEMPI(Temp);
		}
	}

	return;
}

