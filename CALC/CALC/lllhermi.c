/* lllhermi.c */
#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_NONSTDC_NO_DEPRECATE
#include <stdio.h>
#include <stdlib.h>
#include "integer.h"
#include "fun.h"
#include <time.h>
extern MPI *MAXI, *PMAXI;

extern unsigned int MLLLVERBOSE;
extern unsigned int HERMITE1VERBOSE;
unsigned int keith97;

USI FLAGCOLI(MPMATI *Aptr);


MPMATI *LLLHERMITE1(MPMATI *DD, MPMATI **Aptr, USI *rank, USI m1, USI n1)
/*
* Input: an m x n matrix DD of MPI's.
* Output: *Aptr = HNF(DD). Also we return a small transformation matrix B
* if nullity > 0.
* The limiting form of  LLL applied to [I_m|N^mD_1|N^{m-1}D_2|...] is used.
* Normally m1/n1 = 3/4.
* 23/2/97.
*/
{
	unsigned int *COL;
	unsigned int i, j, k, l, m, n, col1, col2, r;
	int t;
	MPI **D, *X, *Y, *Z, *Tmp, *R, ***temp1, ***temp2;
	MPI *M1, *N1;
	MPMATI *L, *B;
	clock_t time_1, time_2;

	time_1 = clock();
	keith97 = FLAGCOLI(DD);
	m = DD->R;
	n = DD->C;
	B = IDENTITYI(m);
	*Aptr = COPYMATI(DD);
	if (keith97) {
		(B->V[m - 1][m - 1])->S = -1;
		for (j = 0; j < n; j++)
			((*Aptr)->V[m - 1][j])->S = -(((*Aptr)->V[m - 1][j])->S);
	}
	D = (MPI **)mmalloc((1 + m) * sizeof(MPI *));
	for (i = 0; i <= m; i++)
		D[i] = ONEI();

	L = ZEROMNI(m, m);

	k = 2;
	M1 = CHANGE(m1);
	N1 = CHANGE(n1);
	while (k <= m)
	{
		REDUCE1(k, k - 1, &L, &B, D, Aptr, &col1, &col2);
		X = MULTI(D[k - 2], D[k]);
		Y = MULTI(D[k - 1], D[k - 1]);
		Tmp = Y;
		Y = MULTI(Y, M1);
		FREEMPI(Tmp);
		R = MULTI(L->V[k - 1][k - 2], L->V[k - 1][k - 2]);
		Z = ADD0I(X, R);
		FREEMPI(R);
		Tmp = Z;
		Z = MULTI(Z, N1);
		FREEMPI(Tmp);
		if ((RSV(Y, Z) == 1 && col1 == n && col2 == n) || (col1 <= col2 && col1 < n))
		{
			for (i = k + 1; i <= m; i++)
				SWAP1(i, k, &L, D);
			SWAP21(k, &B, &L, Aptr);
			FREEMPI(Y);
			Y = MULTI(L->V[k - 1][k - 2], L->V[k - 1][k - 2]);
			Tmp = Y;
			Y = ADD0I(Y, X);
			FREEMPI(Tmp);
			Tmp = D[k - 1];
			D[k - 1] = INT0(Y, D[k - 1]);
			FREEMPI(Tmp);
			if (k > 2)
				k--;
		}
		else
		{
			for (l = k - 2; l >= 1; l--)
				REDUCE1(k, l, &L, &B, D, Aptr, &col1, &col2);
			k++;
		}
		FREEMPI(X);
		FREEMPI(Y);
		FREEMPI(Z);
	}
	FREEMPI(M1);
	FREEMPI(N1);
	FREEMATI(L);
	for (i = 0; i <= m; i++)
		FREEMPI(D[i]);

	COL = (USI *)mmalloc((1 + m) * sizeof(USI));
	COL[0] = 0;
	for (t = m - 1; t >= 0; t--)
	{
		r = COLSEEKI0((*Aptr), t, COL[m - 1 - t]);
		if (r == n)
			break;
		else
			COL[m - t] = r;
	}
	if (t == -1)
		*rank = m;
	else
		*rank = m - t - 1;

	temp1 = (MPI ***)mmalloc(m * sizeof(MPI **));
	temp2 = (MPI ***)mmalloc(m * sizeof(MPI **));
	for (i = 0; i <= m - 1; i++)
	{
		temp1[i] = (*Aptr)->V[m - 1 - i];
		temp2[i] = B->V[m - 1 - i];
	}
	for (i = 0; i <= m - 1; i++)
	{
		(*Aptr)->V[i] = temp1[i];
		B->V[i] = temp2[i];
	}

	ffree((char *)COL, (1 + m) * sizeof(USI));
	ffree((char *)temp1, m * sizeof(MPI **));
	ffree((char *)temp2, m * sizeof(MPI **));
	ffree((char *)D, (1 + m) * sizeof(MPI *));
	time_2 = clock();
	printf("lllhermite:CPU time taken=%f seconds\n", (float)(time_2 - time_1) / (float)CLOCKS_PER_SEC);

	return (B);
}


void REDUCE1(USI k, USI l, MPMATI **Lptr, MPMATI **Bptr, MPI *D[], MPMATI **Aptr, USI *col1, USI *col2)
/*
* updates *Lptr, *Bptr, *Aptr.
*/
{
	unsigned int j, n, t;
	MPI *X, *Y, *Q, *Tmp;
	MPMATI *TempMATI;

	n = (*Aptr)->C;
	*col1 = n;
	for (j = 0; j < n; j++)
	{
		t = ((*Aptr)->V[l - 1][j])->S;
		if (t)
		{
			*col1 = j;
			if (t == -1)
				LLLMINUS(Aptr, Bptr, Lptr, l - 1);
			break;
		}
	}
	*col2 = n;
	for (j = 0; j < n; j++)
	{
		t = ((*Aptr)->V[k - 1][j])->S;
		if (t)
		{
			*col2 = j;
			/*
			if ((k == (*Aptr)->R) && !keith97)
			{
			if (t == -1)
			LLLMINUS(Aptr, Bptr, Lptr, k - 1);
			}
			*/
			break;
		}
	}

	if (*col1 < n)
		/*		Q = NEAREST_INTI((*Aptr)->V[k - 1][*col1], (*Aptr)->V[l - 1][*col1]);*/
		Q = INTI((*Aptr)->V[k - 1][*col1], (*Aptr)->V[l - 1][*col1]);
	else
	{
		Y = MULT_I((*Lptr)->V[k - 1][l - 1], 2);
		if (RSV(Y, D[l]) == 1)
			Q = NEAREST_INTI((*Lptr)->V[k - 1][l - 1], D[l]);
		else
			Q = ZEROI();
		FREEMPI(Y);
	}
	if (Q->S == 0)
	{
		FREEMPI(Q);
		return;
	}
	X = MINUSI(Q);
	TempMATI = *Bptr;
	*Bptr = ADD_MULT_ROWI(l - 1, k - 1, X, *Bptr);
	FREEMATI(TempMATI);
	TempMATI = *Aptr;
	*Aptr = ADD_MULT_ROWI(l - 1, k - 1, X, *Aptr);
	FREEMATI(TempMATI);
	if (HERMITE1VERBOSE)
	{
		printf("Row %u -> Row %u + ", k, k); PRINTI(X); printf(" x Row %u\n", l);
		printf("P = \n");
		PRINTMATI(0, (*Bptr)->R - 1, 0, (*Bptr)->C - 1, *Bptr);
		printf("A = \n");
		PRINTMATI(0, (*Aptr)->R - 1, 0, (*Aptr)->C - 1, *Aptr);
		GetReturn();
	}
	FREEMPI(X);
	for (j = 1; j < l; j++)
	{
		X = MULTI((*Lptr)->V[l - 1][j - 1], Q);
		Tmp = (*Lptr)->V[k - 1][j - 1];
		(*Lptr)->V[k - 1][j - 1] = SUBI((*Lptr)->V[k - 1][j - 1], X);
		FREEMPI(Tmp);
		FREEMPI(X);
	}
	X = MULTI(D[l], Q);
	Tmp = (*Lptr)->V[k - 1][l - 1];
	(*Lptr)->V[k - 1][l - 1] = SUBI((*Lptr)->V[k - 1][l - 1], X);
	FREEMPI(Tmp);
	FREEMPI(X);
	FREEMPI(Q);
	return;
}


void SWAP21(USI k, MPMATI **Bptr, MPMATI **Lptr, MPMATI **Aptr)
{
	MPI *T;

	*Aptr = SWAP_ROWSI1(k - 2, k - 1, *Aptr);
	*Bptr = SWAP_ROWSI1(k - 2, k - 1, *Bptr);
	T = COPYI((*Lptr)->V[k - 1][k - 2]);
	*Lptr = SWAP_ROWSI1(k - 2, k - 1, *Lptr);
	FREEMPI((*Lptr)->V[k - 1][k - 2]);
	(*Lptr)->V[k - 1][k - 2] = T;
	FREEMPI((*Lptr)->V[k - 2][k - 2]);
	(*Lptr)->V[k - 2][k - 2] = ZEROI();
	if (HERMITE1VERBOSE)
	{
		printf("Swapping Rows %u and %u\n", k - 1, k);
		printf("P = \n");
		PRINTMATI(0, (*Bptr)->R - 1, 0, (*Bptr)->C - 1, *Bptr);
		printf("A = \n");
		PRINTMATI(0, (*Aptr)->R - 1, 0, (*Aptr)->C - 1, *Aptr);
		GetReturn();
		return;
	}
}

void LLLMINUS(MPMATI **Aptr, MPMATI **Pptr, MPMATI **Lptr, USI j)
{
	USI r, s, m;
	MPI *Temp;
	MPMATI *TempMATI;

	m = (*Aptr)->R;
	for (r = 0; r < m; r++)
		for (s = 0; s < r; s++)
		{
			if (r == j || s == j)
			{
				Temp = (*Lptr)->V[r][s];
				(*Lptr)->V[r][s] = MINUSI(Temp);
				FREEMPI(Temp);
			}
		}

	TempMATI = *Aptr;
	Temp = MINUS_ONEI();
	*Aptr = SCALAR_ROWI(j, Temp, *Aptr);
	FREEMPI(Temp);
	FREEMATI(TempMATI);
	TempMATI = *Pptr;
	Temp = MINUS_ONEI();
	*Pptr = SCALAR_ROWI(j, Temp, *Pptr);
	FREEMPI(Temp);
	FREEMATI(TempMATI);
	if (HERMITE1VERBOSE)
	{
		printf("Row %u -> - Row %u\n", j + 1, j + 1);
		printf("P = \n");
		PRINTMATI(0, (*Pptr)->R - 1, 0, (*Pptr)->C - 1, *Pptr);
		printf("A = \n");
		PRINTMATI(0, (*Aptr)->R - 1, 0, (*Aptr)->C - 1, *Aptr);
		GetReturn();
	}
	return;
}

MPMATI *SCALAR_ROWI(USI p, MPI *Aptr, MPMATI *Mptr)
/*
* multiplying the pth row of *Mptr by *Aptr.
*/
{
	unsigned int j;
	MPI *Temp;
	MPMATI *Nptr;

	Nptr = COPYMATI(Mptr);
	for (j = 0; j <= Nptr->C - 1; j++)
	{
		Temp = Nptr->V[p][j];
		Nptr->V[p][j] = MULTI(Temp, Aptr);
		FREEMPI(Temp);
	}
	return (Nptr);
}

USI COLSEEKI0(MPMATI *A, USI i, USI j)
/* Scans cols >= i, rows >= j for first nonzero column if any. Returns
* this column number if ther is one, else returns A->C;
*/
{
	USI k, l, n;

	n = A->C;
	for (k = j; k < n; k++)
	{
		for (l = 0; l <= i; l++)
		{
			if (A->V[l][k]->S)
				return k;
		}
	}
	return k;
}

USI FLAGCOLI(MPMATI *Aptr)
/* returns 0 if the first nonzero column j of A contains more than one nonzero entry,
* or contains only one nonzero entry and which is positive.
* returns 1 if the first nonzero column j of A contains only one nonzero entry, which is negative.
* This assumes A is a nonzero matrix with at least two rows.
*/
{
	USI i, j, k, m, n;

	m = Aptr->R;
	n = Aptr->C;
	for (j = 0; j < n; j++) {
		for (i = 0; i < m; i++) {
			if ((Aptr->V[i][j])->S) {
				goto found;
			}
		}
	}
found:
	for (k = i + 1; k < m; k++) {
		if ((Aptr->V[k][j])->S)
			return (0);
	}
	if ((Aptr->V[i][j])->S == 1) {
		return (0);
	}
	else {
		return (1);
	}
}
