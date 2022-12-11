/* lllgcd.c */
#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_NONSTDC_NO_DEPRECATE
#ifdef _WIN32
#include "unistd_DOS.h"
#else
#include <unistd.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "integer.h"
#include "fun.h"
extern MPI *MAXI, *PMAXI;

extern unsigned int MLLLVERBOSE;
extern unsigned int GCDVERBOSE;
USI GCD3FLAG;
extern unsigned int GCD_MAX;

MPMATI *LLLGCD(MPMATI *DD, MPI **Aptr, USI m, USI m1, USI n1)
/*
* Input: an m x 1 vector of positive MPI's.
* Output: *Aptr = gcd of the vector of MPI's. Also we return a small set of
* multipliers using the LLL method of Havas and Matthews.
* Normally m1=1,n1=1.
* matrix B of the LLL extended gcd algorithm of
* Havas-Majewski-Matthews is returned.
* 30/1/97.
*/
{
	unsigned int i, k, l;
	MPI **A, **D, *X, *Y, *Z, *Tmp, *R, *M1, *N1;
	MPMATI *L, *B;

	B = IDENTITYI(m);
	A = (MPI **)mmalloc(m * sizeof(MPI *));
	for (i = 0; i < m; i++)
		A[i] = COPYI(DD->V[i][0]);
	D = (MPI **)mmalloc((1 + m) * sizeof(MPI *));
	for (i = 0; i <= m; i++)
		D[i] = ONEI();

	L = ZEROMNI(m, m);
	M1 = CHANGE(m1);
	N1 = CHANGE(n1);

	k = 2;
	while (k <= m)
	{
		REDUCE(k, k - 1, &L, &B, D, A);
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
		if (!EQZEROI(A[k - 2]) || (EQZEROI(A[k - 1]) && (RSV(Y, Z) == 1)))
		{
			for (i = k + 1; i <= m; i++)
				SWAP1(i, k, &L, D);
			SWAP2(k, &B, &L, A);
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
		{
			FREEMPI(X);
			FREEMPI(Y);
			for (l = k - 2; l >= 1; l--)
				REDUCE(k, l, &L, &B, D, A);
			k++;
		}
		FREEMPI(Z);
	}
	FREEMPI(M1);
	FREEMPI(N1);
	*Aptr = COPYI(A[m - 1]);
	if ((*Aptr)->S == -1)
	{
		LLLGCDMINUS(&B, &L, m - 1);
		(*Aptr)->S = -((*Aptr)->S);
	}
	/*
	if (GCDVERBOSE)
	{
	printf("L = \n");
	PRINTMATI(0,L->R-1,0,L->C-1,L);
	for (i = 0; i <= m; i++)
	{
	printf("D[%u] = ", i);PRINTI(D[i]);printf(", ");
	}
	printf("\n");
	temp1  = MULTI(D[1], D[1]);
	printf("a = "); PRINTI(temp1); printf(", ");
	FREEMPI(temp1);
	temp1 = MULTI(D[1], L->V[1][0]);
	temp2 = MULT_I(temp1, 2);
	printf("2h = "); PRINTI(temp2); printf(", ");
	FREEMPI(temp1);
	FREEMPI(temp2);
	temp1 = MULTI(L->V[1][0], L->V[1][0]);
	temp2 = ADD0I(D[2], temp1);
	printf("b = "); PRINTI(temp2); printf(", ");
	FREEMPI(temp1);
	FREEMPI(temp2);
	temp1 = MULTI(D[1], L->V[2][0]);
	temp2 = MULT_I(temp1, 2);
	printf("2g = "); PRINTI(temp2); printf(", ");
	FREEMPI(temp1);
	FREEMPI(temp2);
	temp1 = MULTI(L->V[1][0], L->V[2][0]);
	FREEMPI(temp1);
	temp2 = ADDI(temp1, L->V[2][1]);
	temp1 = MULT_I(temp2, 2);
	printf("2f = "); PRINTI(temp1); printf(", ");
	FREEMPI(temp1);
	FREEMPI(temp2);
	}
	*/
	/*
	printf("L = \n");
	PRINTMATI(0,L->R-1,0,L->C-1,L);
	for (i = 0; i <= m; i++)
	{
	printf("D[%u] = ", i);PRINTI(D[i]);printf(", ");
	}
	printf("\n");
	*/
	FREEMATI(L);
	for (i = 0; i <= m; i++)
		FREEMPI(D[i]);
	for (i = 0; i < m; i++)
		FREEMPI(A[i]);
	ffree((char *)D, (1 + m) * sizeof(MPI *));
	ffree((char *)A, m * sizeof(MPI *));
	return (B);
}

void SWAP1(USI i, USI k, MPMATI **Lptr, MPI *D[])
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


void REDUCE(USI k, USI l, MPMATI **Lptr, MPMATI **Bptr, MPI *D[], MPI *A[])
/*
* updates *Lptr, *Bptr.
*/
{
	unsigned int i, j, m, n;
	MPI *X, *Y, *Q, *Tmp, *Z;
	MPMATI *TmpMATI, *TempMAT;

	m = (*Bptr)->C;
	n = (*Bptr)->R;
	if (A[l - 1]->S)
	{
		Q = INTI(A[k - 1], A[l - 1]);
		if (Q->S)
		{
			Tmp = A[k - 1];
			Z = MULTI(A[l - 1], Q);
			A[k - 1] = SUBI(A[k - 1], Z);
			FREEMPI(Tmp);
			FREEMPI(Z);
		}
	}
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
	TmpMATI = *Bptr;
	*Bptr = ADD_MULT_ROWI(l - 1, k - 1, X, *Bptr);
	FREEMATI(TmpMATI);
	if (MLLLVERBOSE)
	{
		printf("Row %u -> Row %u + ", k, k); PRINTI(X); printf(" x Row %u\n", l);
		PRINTMATI(0, (*Bptr)->R - 1, 0, (*Bptr)->C - 1, *Bptr);
	}
	if (GCDVERBOSE)
	{
		printf("Row %u -> Row %u + ", k, k); PRINTI(X); printf(" x Row %u\n", l);
		TempMAT = BUILDMATI(n, n + 1);
		for (i = 0; i < n; i++)
		{
			for (j = 0; j <= n; j++)
			{
				if (j == n)
					TempMAT->V[i][j] = COPYI(A[i]);
				else
					TempMAT->V[i][j] = COPYI((*Bptr)->V[i][j]);
			}
		}
		PRINTMATI(0, TempMAT->R - 1, 0, TempMAT->C - 1, TempMAT);
		printf("\n");
		FREEMATI(TempMAT);
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
}


void SWAP2(USI k, MPMATI **B1ptr, MPMATI **Lptr, MPI *A[])
{
	MPI *T, *S;
	unsigned int i, j, n;
	MPMATI *TempMAT;

	n = (*B1ptr)->R;
	S = A[k - 2];
	A[k - 2] = A[k - 1];
	A[k - 1] = S;

	*B1ptr = SWAP_ROWSI1(k - 2, k - 1, *B1ptr);
	T = COPYI((*Lptr)->V[k - 1][k - 2]);
	*Lptr = SWAP_ROWSI1(k - 2, k - 1, *Lptr);
	FREEMPI((*Lptr)->V[k - 1][k - 2]);
	(*Lptr)->V[k - 1][k - 2] = T;
	FREEMPI((*Lptr)->V[k - 2][k - 2]);
	(*Lptr)->V[k - 2][k - 2] = ZEROI();
	if (GCDVERBOSE)
	{
		printf("Swapping Rows %u and %u\n", k - 1, k);
		TempMAT = BUILDMATI(n, n + 1);
		for (i = 0; i < n; i++)
		{
			for (j = 0; j <= n; j++)
			{
				if (j == n)
					TempMAT->V[i][j] = COPYI(A[i]);
				else
					TempMAT->V[i][j] = COPYI((*B1ptr)->V[i][j]);
			}
		}
		PRINTMATI(0, TempMAT->R - 1, 0, TempMAT->C - 1, TempMAT);
		printf("\n");
		GetReturn();
		FREEMATI(TempMAT);
	}
	if (MLLLVERBOSE)
	{
		printf("Swapping Rows %u and %u\n", k - 1, k);
		PRINTMATI(0, (*B1ptr)->R - 1, 0, (*B1ptr)->C - 1, *B1ptr);
	}
	return;
}

MPMATI *JACOBIGCD(MPMATI *DD, MPI **Aptr, USI m)
/*
* Input: an m x 1 vector of positive MPI's.
* Output: *Aptr = gcd of the vector of MPI's. Also we return a small set of
* multipliers using a version of a method of Jacobi
* A unimodular transforming matrix B is returned.
* 31/1/97.
*/
{
	unsigned int i, z;
	int k;
	MPI **A, *Q, *X, *Tmp, **temp;
	MPMATI *L, *B;

	B = IDENTITYI(m);
	A = (MPI **)mmalloc(m * sizeof(MPI *));
	for (i = 0; i < m; i++)
		A[i] = COPYI(DD->V[i][0]);
	z = 0;
	while (1)
	{
		/* A[i] = 0 for i = 0,...,z-1, but A[z] != 0. */
		for (i = z + 1; i < m; i++)
		{
			Q = INT0(A[i], A[z]);
			X = MINUSI(Q);
			FREEMPI(Q);
			L = B;
			B = ADD_MULT_ROWI(z, i, X, B);
			FREEMPI(X);
			FREEMATI(L);
			Tmp = A[i];
			A[i] = MOD0(A[i], A[z]);
			FREEMPI(Tmp);
		}
		if (MLLLVERBOSE)
		{
			printf("After the reduction: B = \n");
			PRINTMATI(0, B->R - 1, 0, B->C - 1, B);
			printf("\n");
			for (i = 0; i < m; i++)
			{
				printf("A[%u] = ", i); PRINTI(A[i]); printf("\n");
			}
			GetReturn();
		}
		Tmp = A[z];
		temp = B->V[z];
		for (k = z; k <= m - 2; k++)
		{
			A[k] = A[k + 1];
			B->V[k] = B->V[k + 1];
		}
		A[m - 1] = Tmp;
		B->V[m - 1] = temp;
		if (MLLLVERBOSE)
		{
			printf("After the rearrangement: B = \n");
			PRINTMATI(0, B->R - 1, 0, B->C - 1, B);
			printf("\n");
			for (i = 0; i < m; i++)
			{
				printf("A[%i] = ", i); PRINTI(A[i]); printf("\n");
			}
			GetReturn();
		}
		while (A[z]->S == 0)
			z++;
		if (z == m - 1)
			break;
	}

	*Aptr = COPYI(A[m - 1]);
	if ((*Aptr)->S == -1)
	{
		(*Aptr)->S = 1;
		for (i = 0; i < m; i++)
			(B->V[m - 1][i])->S = -((B->V[m - 1][i])->S);
	}
	for (i = 0; i < m; i++)
		FREEMPI(A[i]);
	ffree((char *)A, m * sizeof(MPI *));
	return (B);
}

void GCD4()
/* We compute an optimal multipliers for a[1],a[2],a[3],a[4]
* in the range p1<=a[1]<=q1 etc. and rel prime.
*/
{
	USI l, m, p1, q1, p2, q2, p3, q3, p4, q4, m1, n1, t;
	USI i1, i2, i3, i4, p;
	int e, s;
	USL x, y, z;
	MPI *GCD, *T1, *T2, **XX, **X, *Temp, *DIFF;
	MPMATI *MATI1, *BB, *MATI3, *Q, *M;
	char buff[20];
	FILE *outfile;


	printf("Enter alpha=m/n,  1/4 < m/n <= 1:  (ie enter  m  and n)");
	s = scanf("%u%u", &m1, &n1);
	strcpy(buff, "gcd4.out");
	if (access(buff, R_OK) == 0)
		unlink(buff);
	printf("Enter  the ranges pi,qi; 2 <= pi < qi < 2^32, i = 1,...,5: ");
	s = scanf("%u%u%u%u%u%u%u%u", &p1, &q1, &p2, &q2, &p3, &q3, &p4, &q4);
	Flush();
	m = 4;
	p = m - 1;
	for (i1 = p1; i1 <= q1; i1++)
	{
		/*for (i2 = (i1+1 > p2 ? i1+1: p2); i2 <= q2; i2++)*/
		for (i2 = p2; i2 <= q2; i2++)
		{
			/*for (i3 = (i2+1 > p3 ? i2+1: p3); i3 <= q3; i3++)*/
			for (i3 = p3; i3 <= q3; i3++)
			{
				/*for (i4 = (i3+1 > p4 ? i3+1: p4); i4 <= q4; i4++)*/
				for (i4 = p4; i4 <= q4; i4++)
				{
					x = GCDm((USL)i1, (USL)i2);
					y = GCDm(x, (USL)i3);
					z = GCDm(y, (USL)i4);
					if (z == 1)
					{
						printf("(i1,i2,i3,i4)=(%u,%u,%u,%u)\n", i1, i2, i3, i4);
						MATI1 = BUILDMATI(m, 1);
						MATI1->V[0][0] = CHANGE((USL)i1);
						MATI1->V[1][0] = CHANGE((USL)i2);
						MATI1->V[2][0] = CHANGE((USL)i3);
						MATI1->V[3][0] = CHANGE((USL)i4);
						BB = LLLGCD0M(MATI1, &GCD, m, m1, n1);
						FREEMATI(MATI1);
						FREEMPI(GCD);
						MATI3 = BUILDMATI(1, m);
						for (l = 0; l < m; l++)
							MATI3->V[0][l] = COPYI(BB->V[m - 1][l]);
						T1 = LENGTHSQRI(MATI3, 0);
						Q = BB;
						XX = (MPI **)mmalloc((USL)(p * sizeof(MPI *)));
						for (t = 0; t < p; t++)
							XX[t] = ZEROI();
						while (1)
						{
							M = SHORTESTT0(Q, &X);
							for (t = 0; t < p; t++)
							{
								Temp = XX[t];
								XX[t] = ADDI(XX[t], X[t]);
								FREEMPI(Temp);
								FREEMPI(X[t]);
							}
							ffree((char *)X, p * sizeof(MPI *));
							if (M == NULL)
								break;
							else
							{
								for (l = 0; l < Q->C; l++)
								{
									FREEMPI(Q->V[m - 1][l]);
									Q->V[m - 1][l] = COPYI(M->V[0][l]);
								}
								FREEMATI(MATI3);
								MATI3 = M;
							}
						}
						T2 = LENGTHSQRI(MATI3, 0);
						if (RSV(T2, T1) == -1)
						{
							outfile = fopen(buff, "a");
							fprintf(outfile, "(%u,%u,%u,%u): ", i1, i2, i3, i4);
							fprintf(outfile, "b[%u]", m);
							for (t = 0; t < p; t++)
							{
								e = XX[t]->S;
								if (e == -1)
								{
									fprintf(outfile, "+");
									Temp = MINUSI(XX[t]);
									if (!EQONEI(Temp))
										FPRINTI(outfile, Temp);
									fprintf(outfile, "b[%u]", t + 1);
									FREEMPI(Temp);
								}
								if (e == 1)
								{
									fprintf(outfile, "-");
									if (!EQONEI(XX[t]))
										FPRINTI(outfile, XX[t]);
									fprintf(outfile, "b[%u]", t + 1);
								}
							}
							fprintf(outfile, ": ");
							DIFF = SUB0I(T1, T2);
							FPRINTI(outfile, DIFF);
							FREEMPI(DIFF);
							fprintf(outfile, "\n");
							fclose(outfile);
						}
						for (t = 0; t < p; t++)
							FREEMPI(XX[t]);
						ffree((char *)XX, p * sizeof(MPI *));
						FREEMATI(BB);
						FREEMATI(MATI3);
						FREEMPI(T1);
						FREEMPI(T2);
					}
				}
			}
		}
	}
	return;
}

void GCD5()
/* We compute an multipliers for a[1],a[2],a[3],a[4],a[5]
* in the range p1<=a[1]<=q1 etc. and rel prime.
*/
{
	USI l, m, p1, q1, p2, q2, p3, q3, p4, q4, p5, q5, m1, n1;
	USI i1, i2, i3, i4, i5, p, t;
	int e, s;
	USL x, y, z, u;
	MPI *GCD, *T1, *T2, **XX, **X, *Temp, *DIFF;
	MPMATI *MATI1, *BB, *MATI3, *Q, *M;
	char buff[20];
	FILE *outfile;

	/*for (i2 = (i1+1 > p2 ? i1+1: p2); i2 <= q2; i2++)*/
	/*for (i3 = (i2+1 > p3 ? i2+1: p3); i3 <= q3; i3++)*/
	/*for (i4 = (i3+1 > p4 ? i3+1: p4); i4 <= q4; i4++)*/
	/*for (i5 = (i4+1 > p5 ? i4+1: p5); i5 <= q5; i5++)*/

	printf("Enter alpha=m/n,  1/4 < m/n <= 1:  (ie enter  m  and n)");
	s = scanf("%u%u", &m1, &n1);
	strcpy(buff, "gcd5.out");
	if (access(buff, R_OK) == 0)
		unlink(buff);
	printf("Enter  the ranges pi,qi; 2 <= pi < qi < 2^32, i = 1,...,5: ");
	s = scanf("%u%u%u%u%u%u%u%u%u%u", &p1, &q1, &p2, &q2, &p3, &q3, &p4, &q4, &p5, &q5);
	Flush();
	m = 5;
	p = m - 1;
	for (i1 = p1; i1 <= q1; i1++)
	{
		for (i2 = p2; i2 <= q2; i2++)
		{
			for (i3 = p3; i3 <= q3; i3++)
			{
				for (i4 = p4; i4 <= q4; i4++)
				{
					for (i5 = p5; i5 <= q5; i5++)
					{
						x = GCDm((USL)i1, (USL)i2);
						y = GCDm(x, (USL)i3);
						z = GCDm(y, (USL)i4);
						u = GCDm(z, (USL)i5);
						if (u == 1)
						{
							printf("(i1,i2,i3,i4,i5)=(%u,%u,%u,%u,%u)\n", i1, i2, i3, i4, i5);
							MATI1 = BUILDMATI(m, 1);
							MATI1->V[0][0] = CHANGE((USL)i1);
							MATI1->V[1][0] = CHANGE((USL)i2);
							MATI1->V[2][0] = CHANGE((USL)i3);
							MATI1->V[3][0] = CHANGE((USL)i4);
							MATI1->V[4][0] = CHANGE((USL)i5);
							BB = LLLGCD0M(MATI1, &GCD, m, m1, n1);
							FREEMATI(MATI1);
							FREEMPI(GCD);
							MATI3 = BUILDMATI(1, m);
							for (l = 0; l < m; l++)
								MATI3->V[0][l] = COPYI(BB->V[m - 1][l]);
							T1 = LENGTHSQRI(MATI3, 0);
							Q = BB;
							XX = (MPI **)mmalloc((USL)(p * sizeof(MPI *)));
							for (t = 0; t < p; t++)
								XX[t] = ZEROI();
							while (1)
							{
								M = SHORTESTT0(Q, &X);
								for (t = 0; t < p; t++)
								{
									Temp = XX[t];
									XX[t] = ADDI(XX[t], X[t]);
									FREEMPI(Temp);
									FREEMPI(X[t]);
								}
								ffree((char *)X, p * sizeof(MPI *));
								if (M == NULL)
									break;
								else
								{
									for (l = 0; l < Q->C; l++)
									{
										FREEMPI(Q->V[m - 1][l]);
										Q->V[m - 1][l] = COPYI(M->V[0][l]);
									}
									FREEMATI(MATI3);
									MATI3 = M;
								}
							}
							T2 = LENGTHSQRI(MATI3, 0);
							outfile = fopen(buff, "a");
							if (RSV(T2, T1) == -1)
							{
								fprintf(outfile, "(%u,%u,%u,%u,%u): ", i1, i2, i3, i4, i5);
								fprintf(outfile, "b[%u]", m);
								for (t = 0; t < p; t++)
								{
									e = XX[t]->S;
									if (e == -1)
									{
										fprintf(outfile, "+");
										Temp = MINUSI(XX[t]);
										if (!EQONEI(Temp))
											FPRINTI(outfile, Temp);
										fprintf(outfile, "b[%u]", t + 1);
										FREEMPI(Temp);
									}
									if (e == 1)
									{
										fprintf(outfile, "-");
										if (!EQONEI(XX[t]))
											FPRINTI(outfile, XX[t]);
										fprintf(outfile, "b[%u]", t + 1);
									}
								}
								fprintf(outfile, ": ");
								DIFF = SUB0I(T1, T2);
								FPRINTI(outfile, DIFF);
								FREEMPI(DIFF);
								fprintf(outfile, "\n");
							}
							fclose(outfile);
							for (t = 0; t < p; t++)
								FREEMPI(XX[t]);
							ffree((char *)XX, p * sizeof(MPI *));
							FREEMATI(BB);
							FREEMATI(MATI3);
							FREEMPI(T1);
							FREEMPI(T2);
						}
					}
				}
			}
		}
	}
	return;
}

void LLLGCDMINUS(MPMATI **Pptr, MPMATI **Lptr, USI j)
{
	USI r, s, m;
	MPI *Temp;
	MPMATI *TempMATI;

	m = (*Pptr)->R;
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

	TempMATI = *Pptr;
	Temp = MINUS_ONEI();
	*Pptr = SCALAR_ROWI(j, Temp, *Pptr);
	FREEMPI(Temp);
	FREEMATI(TempMATI);
	if (GCDVERBOSE)
	{
		printf("Row %u -> - Row %u\n", j, j);
		printf("P = \n");
		PRINTMATI(0, (*Pptr)->R - 1, 0, (*Pptr)->C - 1, *Pptr);
		GetReturn();
	}
	return;
}

void GCD6()
/* We compute an optimal multipliers for a[1],a[2],a[3],a[4],a[5],a[6]
* in the range p1<=a[1]<=q1 etc. and rel prime.
*/
{
	USI l, m, p1, q1, p2, q2, p3, q3, p4, q4, p5, q5, p6, q6, m1, n1;
	USI i1, i2, i3, i4, i5, i6, p, t;
	int e, s;
	USL x, y, z, u, v, w;
	MPI *GCD, *T1, *T2, **XX, **X, *Temp;
	MPMATI *MATI1, *BB, *MATI3, *Q, *M;
	char buff[20];
	FILE *outfile;

	/*for (i2 = (i1+1 > p2 ? i1+1: p2); i2 <= q2; i2++)*/
	/*for (i3 = (i2+1 > p3 ? i2+1: p3); i3 <= q3; i3++)*/
	/*for (i4 = (i3+1 > p4 ? i3+1: p4); i4 <= q4; i4++)*/
	/*for (i5 = (i4+1 > p5 ? i4+1: p5); i5 <= q5; i5++)*/
	/*for (i6 = (i5+1 > p6 ? i5+1: p6); i6 <= q6; i6++)*/

	printf("Enter alpha=m/n,  1/4 < m/n <= 1:  (ie enter  m  and n)");
	s = scanf("%u%u", &m1, &n1);
	strcpy(buff, "gcd6.out");
	if (access(buff, R_OK) == 0)
		unlink(buff);
	printf("Enter  the ranges pi,qi; 2 <= pi < qi < 2^32, i = 1,...,6: ");
	s = scanf("%u%u%u%u%u%u%u%u%u%u%u%u", &p1, &q1, &p2, &q2, &p3, &q3, &p4, &q4, &p5, &q5, &p6, &q6);
	Flush();
	m = 6;
	p = m - 1;
	for (i1 = p1; i1 <= q1; i1++)
	{
		for (i2 = p2; i2 <= q2; i2++)
		{
			for (i3 = p3; i3 <= q3; i3++)
			{
				for (i4 = p4; i4 <= q4; i4++)
				{
					for (i5 = p5; i5 <= q5; i5++)
					{
						for (i6 = p6; i6 <= q6; i6++)
						{
							x = GCDm((USL)i1, (USL)i2);
							y = GCDm(x, (USL)i3);
							z = GCDm(y, (USL)i4);
							u = GCDm(z, (USL)i5);
							v = GCDm(u, (USL)i5);
							w = GCDm(v, (USL)i6);
							if (v == 1)
							{
								printf("(i1,i2,i3,i4,i5,i6)=(%u,%u,%u,%u,%u,%u)\n", i1, i2, i3, i4, i5, i6);
								MATI1 = BUILDMATI(m, 1);
								MATI1->V[0][0] = CHANGE((USL)i1);
								MATI1->V[1][0] = CHANGE((USL)i2);
								MATI1->V[2][0] = CHANGE((USL)i3);
								MATI1->V[3][0] = CHANGE((USL)i4);
								MATI1->V[4][0] = CHANGE((USL)i5);
								MATI1->V[5][0] = CHANGE((USL)i6);
								BB = LLLGCD(MATI1, &GCD, m, m1, n1);
								FREEMATI(MATI1);
								FREEMPI(GCD);
								MATI3 = BUILDMATI(1, m);
								for (l = 0; l < m; l++)
									MATI3->V[0][l] = COPYI(BB->V[m - 1][l]);
								T1 = LENGTHSQRI(MATI3, 0);
								Q = BB;
								XX = (MPI **)mmalloc((USL)(p * sizeof(MPI *)));
								for (t = 0; t < p; t++)
									XX[t] = ZEROI();
								while (1)
								{
									M = SHORTESTT0(Q, &X);
									for (t = 0; t < p; t++)
									{
										Temp = XX[t];
										XX[t] = ADDI(XX[t], X[t]);
										FREEMPI(Temp);
										FREEMPI(X[t]);
									}
									ffree((char *)X, p * sizeof(MPI *));
									if (M == NULL)
										break;
									else
									{
										for (l = 0; l < Q->C; l++)
										{
											FREEMPI(Q->V[m - 1][l]);
											Q->V[m - 1][l] = COPYI(M->V[0][l]);
										}
										FREEMATI(MATI3);
										MATI3 = M;
									}
								}
								T2 = LENGTHSQRI(MATI3, 0);
								if (RSV(T2, T1) == -1)
								{
									outfile = fopen(buff, "a");
									fprintf(outfile, "(%u,%u,%u,%u,%u,%u): ", i1, i2, i3, i4, i5, i6);
									fprintf(outfile, "b[%u]", m);
									for (t = 0; t < p; t++)
									{
										e = XX[t]->S;
										if (e == -1)
										{
											fprintf(outfile, "+");
											Temp = MINUSI(XX[t]);
											if (!EQONEI(Temp))
												FPRINTI(outfile, Temp);
											fprintf(outfile, "b[%u]", t + 1);
											FREEMPI(Temp);
										}
										if (e == 1)
										{
											fprintf(outfile, "-");
											if (!EQONEI(XX[t]))
												FPRINTI(outfile, XX[t]);
											fprintf(outfile, "b[%u]", t + 1);
										}
									}
									fprintf(outfile, "\n");
									fclose(outfile);
								}
								for (t = 0; t < p; t++)
									FREEMPI(XX[t]);
								ffree((char *)XX, p * sizeof(MPI *));
								FREEMATI(BB);
								FREEMATI(MATI3);
								FREEMPI(T1);
								FREEMPI(T2);
							}
						}
					}
				}
			}
		}
	}
	return;
}

MPMATI *LLLGCD0(MPMATI *DD, MPI **Aptr, USI m, USI m1, USI n1)
/*
* Input: an m x 1 vector of positive MPI's.
* Output: *Aptr = gcd of the vector of MPI's. Also we return a small set of
* multipliers using the LLL method of Havas and Matthews.
* Normally m1=1,n1=1.
* Also returns matrix B of the LLL extended gcd algorithm of
* Havas-Majewski-Matthews is returned.
* S is the shortest length squared all the m short vectors returned.
* 30/1/97.
*/
{
	unsigned int i, j, k, l, flag, p;
	MPI **A, **D, *X, *Y, *Z, *Tmp, *R, *M1, *N1;
	MPI **XX, *temp, *mu, *SUM, *SHORTER;
	MPI *temp1, *T, *tempI2, *Temp;
	MPMATI *L, *B, *MATI1, *MATI2, *MATI;
	MPR *tempR1, *tR1, *tR2, *tR3, *SUMR, **RHO;
	int e, r, kk, K, c, COUNT, count;
	char buff[20];
	FILE *outfile;

	B = IDENTITYI(m);
	A = (MPI **)mmalloc(m * sizeof(MPI *));
	for (i = 0; i < m; i++)
		A[i] = COPYI(DD->V[i][0]);
	D = (MPI **)mmalloc((1 + m) * sizeof(MPI *));
	for (i = 0; i <= m; i++)
		D[i] = ONEI();

	L = ZEROMNI(m, m);
	M1 = CHANGE(m1);
	N1 = CHANGE(n1);

	k = 2;
	while (k <= m)
	{
		REDUCE(k, k - 1, &L, &B, D, A);
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
		if (!EQZEROI(A[k - 2]) || (EQZEROI(A[k - 1]) && (RSV(Y, Z) == 1)))
		{
			for (i = k + 1; i <= m; i++)
				SWAP1(i, k, &L, D);
			SWAP2(k, &B, &L, A);
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
		{
			FREEMPI(X);
			FREEMPI(Y);
			for (l = k - 2; l >= 1; l--)
				REDUCE(k, l, &L, &B, D, A);
			k++;
		}
		FREEMPI(Z);
	}
	FREEMPI(M1);
	FREEMPI(N1);
	*Aptr = COPYI(A[m - 1]);
	if ((*Aptr)->S == -1)
	{
		LLLGCDMINUS(&B, &L, m - 1);
		(*Aptr)->S = -((*Aptr)->S);
	}
	printf("L = \n");
	PRINTMATI(0, L->R - 1, 0, L->C - 1, L);
	for (i = 0; i <= m; i++)
	{
		printf("D[%u] = ", i); PRINTI(D[i]); printf(", ");
	}
	printf("\n");
	GetReturn();

	RHO = (MPR **)mmalloc((m - 1) * sizeof(MPR *));
	for (K = m - 2; K >= 0; K--)
	{
		SUMR = RATIOI(L->V[m - 1][K], D[K + 1]);
		for (r = m - 2; r > K; r--)
		{
			tR1 = SUMR;
			tR2 = RATIOI(L->V[r][K], D[K + 1]);
			tR3 = MULTR(tR2, RHO[r]);
			FREEMPR(tR2);
			SUMR = ADDR(SUMR, tR3);
			FREEMPR(tR1);
			FREEMPR(tR3);
		}
		RHO[K] = MINUSR(SUMR);
		FREEMPR(SUMR);
	}
	strcpy(buff, "lllgcd0mult.out");
	outfile = fopen(buff, "w");
	for (K = 0; K < m - 1; K++)
	{
		fprintf(outfile, "RHO[%d]=", K + 1);
		printf("RHO[%d]=", K + 1);
		PRINTR(RHO[K]);
		FPRINTR(outfile, RHO[K]);
		printf("\n");
		fprintf(outfile, "\n");
	}

	for (K = 0; K < m - 1; K++)
	{
		fprintf(outfile, "RHO[%d]=", K + 1);
		printf("RHO[%d]=", K + 1);
		PRINTDR(2, RHO[K]);
		FPRINTDR(outfile, 2, RHO[K]);
		FREEMPR(RHO[K]);
		fprintf(outfile, "\n");
		printf("\n");
	}
	ffree((char *)RHO, (m - 1) * sizeof(MPR *));
	GetReturn();
	XX = (MPI **)mmalloc(m * sizeof(MPI *));
	for (i = 0; i < m - 1; i++)
		XX[i] = ZEROI();
	XX[m - 1] = ONEI();
	COUNT = m - 2;
	SHORTER = ZEROI();
	MATI = BUILDMATI(1, m);
	for (j = 0; j < m; j++)
		MATI->V[0][j] = ZEROI();
	for (K = m - 2; K >= 0; K--)
	{
		fprintf(outfile, "X[%d]= ", K + 1);
		printf("X[%d]=", K + 1);
		flag = 1;
		for (k = 0; k < m - 1; k++)
		{
			FREEMPI(XX[k]);
			XX[k] = ZEROI();
		}
		for (kk = K; kk >= 0; kk--)
		{
			mu = COPYI(L->V[m - 1][kk]);
			if (flag)
			{
				flag = 0;
				if (mu->S == 0)
				{
					FREEMPI(mu);
					continue;
				}
				else
				{
					/*
					tempI1 = ADDI(L->V[m - 1][kk], L->V[m - 1][kk]);
					tempI2= ABSI(tempI1);
					t = EQUALI(tempI2, D[kk]);
					FREEMPI(tempI1);
					FREEMPI(tempI2);
					if (t)
					{
					FREEMPI(mu);
					continue;
					}
					*/
					temp = XX[kk];
					if (mu->S == 1)
						XX[kk] = MINUS_ONEI();
					else
						XX[kk] = ONEI();
					FREEMPI(temp);
				}
			}
			else
			{
				SUM = COPYI(mu);
				for (r = m - 2; r >= kk + 1; r--)
				{
					temp = SUM;
					tempI2 = MULTI(XX[r], L->V[r][kk]);
					SUM = ADDI(SUM, tempI2);
					FREEMPI(tempI2);
					FREEMPI(temp);
				}
				tempR1 = RATIOI(SUM, D[kk + 1]);
				FREEMPI(SUM);
				temp1 = ABS_NEAREST_INTR(tempR1);
				/*
				if (K == 15 && kk == 14)
				{
				printf("tempR1=");PRINTR(tempR1);printf("\n");
				printf("temp1=");PRINTI(temp1);printf("\n");
				GetReturn();
				temp = temp1;
				temp1 = ZEROI();
				FREEMPI(temp);
				}
				*/
				FREEMPR(tempR1);
				temp = XX[kk];
				XX[kk] = MINUSI(temp1);
				FREEMPI(temp1);
				FREEMPI(temp);
			}
			FREEMPI(mu);
		}
		MATI1 = BUILDMATI(1, m);
		for (j = 0; j < m; j++)
			MATI1->V[0][j] = COPYI(XX[j]);
		MATI2 = MULTMATI(MATI1, B);
		for (i = 0; i < m; i++)
		{
			PRINTI(MATI2->V[0][i]); printf(" ");
			FPRINTI(outfile, MATI2->V[0][i]); fprintf(outfile, " ");

		}
		printf("=b[%u]", m);
		fprintf(outfile, "=b[%u]", m);
		p = m - 1;
		for (j = 0; j < p; j++)
		{
			e = XX[j]->S;
			if (e == -1)
			{
				printf("-");
				fprintf(outfile, "-");
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
				printf("+");
				fprintf(outfile, "+");
				if (!EQONEI(XX[j]))
				{
					PRINTI(XX[j]);
					FPRINTI(outfile, XX[j]);
				}
				printf("b[%u]", j + 1);
				fprintf(outfile, "b[%u]", j + 1);
			}
		}
		printf(": ");
		fprintf(outfile, ": ");
		T = LENGTHSQRI(MATI2, 0);
		PRINTI(T);
		FPRINTI(outfile, T);
		printf("\n");
		fprintf(outfile, "\n");
		if (K == m - 2)
		{
			FREEMPI(SHORTER);
			SHORTER = COPYI(T);
			FREEMATI(MATI);
			MATI = COPYMATI(MATI2);
		}
		else
		{
			c = RSV(T, SHORTER);
			if (c < 0)
			{
				FREEMPI(SHORTER);
				SHORTER = COPYI(T);
				FREEMATI(MATI);
				MATI = COPYMATI(MATI2);
				COUNT = K;
			}
		}
		FREEMPI(T);
		FREEMATI(MATI1);
		FREEMATI(MATI2);
	}
	fprintf(outfile, "X[0]=");
	printf("X[0]=");
	for (i = 0; i < m; i++)
	{
		PRINTI(B->V[m - 1][i]);
		printf(" ");
		FPRINTI(outfile, B->V[m - 1][i]);
		fprintf(outfile, " ");

	}
	fprintf(outfile, "=b[%u]: ", m);
	printf("=b[%u]: ", m);
	T = LENGTHSQRI(B, m - 1);
	FPRINTI(outfile, T);
	PRINTI(T);
	fprintf(outfile, "\n");
	printf("\n");
	c = RSV(T, SHORTER);
	if (c < 0)
	{
		FREEMPI(SHORTER);
		SHORTER = COPYI(T);
		for (i = 0; i < m; i++)
		{
			FREEMPI(MATI->V[0][i]);
			MATI->V[0][i] = COPYI(B->V[m - 1][i]);
		}
		COUNT = K;
	}
	FREEMPI(T);
	count = COUNT + 1;
	fprintf(outfile, "shortest of the vectors is number %d:", count);
	printf("shortest of the vectors is number %d:", count);
	for (i = 0; i < m; i++)
	{
		PRINTI(MATI->V[0][i]); printf(" ");
		FPRINTI(outfile, MATI->V[0][i]); fprintf(outfile, " ");

	}
	FREEMATI(MATI);
	fprintf(outfile, ": ");
	printf(": ");
	FPRINTI(outfile, SHORTER);
	fprintf(outfile, "\n");
	PRINTI(SHORTER);
	printf("\n");
	FREEMPI(SHORTER);
	fclose(outfile);
	FREEMATI(L);
	for (i = 0; i <= m; i++)
		FREEMPI(D[i]);
	ffree((char *)D, (1 + m) * sizeof(MPI *));
	for (i = 0; i < m; i++)
		FREEMPI(A[i]);
	ffree((char *)A, m * sizeof(MPI *));
	for (i = 0; i < m; i++)
		FREEMPI(XX[i]);
	ffree((char *)XX, m * sizeof(MPI *));
	return (B);
}

MPR *MPI_TO_MPR(MPI *N)
/* Converts an MPI to an MPR */
{
	MPR *T;

	T = BUILDMPR();
	T->N = COPYI(N);
	T->D = ONEI();
	return (T);
}

void GCD10()
/* We compute an optimal multipliers for a[1],...,a[10]
* in the range p1<=a[1]<=q1 etc. and rel prime.
*/
{
	USI l, m, p1, q1, p2, q2, p3, q3, p4, q4, p5, q5, p6, q6, m1, n1;
	USI p7, q7, p8, q8, p9, q9, p10, q10;
	USI i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, p, t;
	int e, s;
	USL x, y, z, u, v, v1, v2, v3, v4;
	MPI *GCD, *T1, *T2, **XX, **X, *Temp, *DIFF;
	MPMATI *MATI1, *BB, *MATI3, *Q, *M;
	char buff[20];
	FILE *outfile;

	printf("Enter alpha=m/n,  1/4 < m/n <= 1:  (ie enter  m  and n)");
	s = scanf("%u%u", &m1, &n1);
	strcpy(buff, "gcd10.out");
	if (access(buff, R_OK) == 0)
		unlink(buff);
	printf("Enter  the ranges pi,qi; 2 <= pi < qi < 2^32, i = 1,...,6: ");
	s = scanf("%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u", &p1, &q1, &p2, &q2, &p3, &q3, &p4, &q4, &p5, &q5, &p6, &q6, &p7, &q7, &p8, &q8, &p9, &q9, &p10, &q10);
	Flush();
	m = 10;
	p = m - 1;
	for (i1 = p1; i1 <= q1; i1++)
	{
		for (i2 = p2; i2 <= q2; i2++)
		{
			for (i3 = p3; i3 <= q3; i3++)
			{
				for (i4 = p4; i4 <= q4; i4++)
				{
					for (i5 = p5; i5 <= q5; i5++)
					{
						for (i6 = p6; i6 <= q6; i6++)
						{
							for (i7 = p7; i7 <= q7; i7++)
							{
								for (i8 = p8; i8 <= q8; i8++)
								{
									for (i9 = p9; i9 <= q9; i9++)
									{
										for (i10 = p10; i10 <= q10; i10++)
										{
											x = GCDm((USL)i1, (USL)i2);
											y = GCDm(x, (USL)i3);
											z = GCDm(y, (USL)i4);
											u = GCDm(z, (USL)i5);
											v = GCDm(u, (USL)i6);
											v1 = GCDm(v, (USL)i7);
											v2 = GCDm(v1, (USL)i8);
											v3 = GCDm(v2, (USL)i9);
											v4 = GCDm(v3, (USL)i10);
											if (v4 == 1)
											{
												printf("(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10)=(%u,%u,%u,%u,%u,%u,%u,%u,%u,%u)\n", i1, i2, i3, i4, i5, i6, i7, i8, i9, i10);
												MATI1 = BUILDMATI(m, 1);
												MATI1->V[0][0] = CHANGE((USL)i1);
												MATI1->V[1][0] = CHANGE((USL)i2);
												MATI1->V[2][0] = CHANGE((USL)i3);
												MATI1->V[3][0] = CHANGE((USL)i4);
												MATI1->V[4][0] = CHANGE((USL)i5);
												MATI1->V[5][0] = CHANGE((USL)i6);
												MATI1->V[6][0] = CHANGE((USL)i7);
												MATI1->V[7][0] = CHANGE((USL)i8);
												MATI1->V[8][0] = CHANGE((USL)i9);
												MATI1->V[9][0] = CHANGE((USL)i10);
												BB = LLLGCD0M(MATI1, &GCD, m, m1, n1);
												FREEMATI(MATI1);
												FREEMPI(GCD);
												MATI3 = BUILDMATI(1, m);
												for (l = 0; l < m; l++)
													MATI3->V[0][l] = COPYI(BB->V[m - 1][l]);
												T1 = LENGTHSQRI(MATI3, 0);
												Q = BB;
												XX = (MPI **)mmalloc((USL)(p * sizeof(MPI *)));
												for (t = 0; t < p; t++)
													XX[t] = ZEROI();
												while (1)
												{
													M = SHORTESTT0(Q, &X);
													for (t = 0; t < p; t++)
													{
														Temp = XX[t];
														XX[t] = ADDI(XX[t], X[t]);
														FREEMPI(Temp);
														FREEMPI(X[t]);
													}
													ffree((char *)X, p * sizeof(MPI *));
													if (M == NULL)
														break;
													else
													{
														for (l = 0; l < Q->C; l++)
														{
															FREEMPI(Q->V[m - 1][l]);
															Q->V[m - 1][l] = COPYI(M->V[0][l]);
														}
														FREEMATI(MATI3);
														MATI3 = M;
													}
												}
												T2 = LENGTHSQRI(MATI3, 0);
												if (RSV(T2, T1) == -1)
												{
													outfile = fopen(buff, "a");
													fprintf(outfile, "(%u,%u,%u,%u,%u,%u,%u,%u,%u,%u): ", i1, i2, i3, i4, i5, i6, i7, i8, i9, i10);
													fprintf(outfile, "b[%u]", m);
													for (t = 0; t < p; t++)
													{
														e = XX[t]->S;
														if (e == -1)
														{
															fprintf(outfile, "+");
															Temp = MINUSI(XX[t]);
															if (!EQONEI(Temp))
																FPRINTI(outfile, Temp);
															fprintf(outfile, "b[%u]", t + 1);
															FREEMPI(Temp);
														}
														if (e == 1)
														{
															fprintf(outfile, "-");
															if (!EQONEI(XX[t]))
																FPRINTI(outfile, XX[t]);
															fprintf(outfile, "b[%u]", t + 1);
														}
													}
													fprintf(outfile, ": ");
													DIFF = SUB0I(T1, T2);
													FPRINTI(outfile, DIFF);
													FREEMPI(DIFF);
													fprintf(outfile, "\n");
													fclose(outfile);
												}
												for (t = 0; t < p; t++)
													FREEMPI(XX[t]);
												ffree((char *)XX, p * sizeof(MPI *));
												FREEMATI(BB);
												FREEMATI(MATI3);
												FREEMPI(T1);
												FREEMPI(T2);
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	return;
}

void GCD11()
/* We compute an optimal multipliers for a[1],...,a[11]
* in the range p1<=a[1]<=q1 etc. and rel prime.
*/
{
	USI l, m, p1, q1, p2, q2, p3, q3, p4, q4, p5, q5, p6, q6, m1, n1;
	USI p7, q7, p8, q8, p9, q9, p10, q10, p11, q11;
	USI i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, p, t, i11;
	int e, s;
	USL x, y, z, u, v, v1, v2, v3, v4, v5;
	MPI *GCD, *T1, *T2, **XX, **X, *Temp;
	MPMATI *MATI1, *BB, *MATI3, *Q, *M;
	char buff[20];
	FILE *outfile;

	printf("Enter alpha=m/n,  1/4 < m/n <= 1:  (ie enter  m  and n)");
	s = scanf("%u%u", &m1, &n1);
	strcpy(buff, "gcd11.out");
	if (access(buff, R_OK) == 0)
		unlink(buff);
	printf("Enter  the ranges pi,qi; 2 <= pi < qi < 2^32, i = 1,...,11: ");
	s = scanf("%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u%u", &p1, &q1, &p2, &q2, &p3, &q3, &p4, &q4, &p5, &q5, &p6, &q6, &p7, &q7, &p8, &q8, &p9, &q9, &p10, &q10, &p11, &q11);
	Flush();
	m = 11;
	p = m - 1;
	for (i1 = p1; i1 <= q1; i1++)
	{
		for (i2 = p2; i2 <= q2; i2++)
		{
			for (i3 = p3; i3 <= q3; i3++)
			{
				for (i4 = p4; i4 <= q4; i4++)
				{
					for (i5 = p5; i5 <= q5; i5++)
					{
						for (i6 = p6; i6 <= q6; i6++)
						{
							for (i7 = p7; i7 <= q7; i7++)
							{
								for (i8 = p8; i8 <= q8; i8++)
								{
									for (i9 = p9; i9 <= q9; i9++)
									{
										for (i10 = p10; i10 <= q10; i10++)
										{
											for (i11 = p11; i11 <= q11; i11++)
											{
												x = GCDm((USL)i1, (USL)i2);
												y = GCDm(x, (USL)i3);
												z = GCDm(y, (USL)i4);
												u = GCDm(z, (USL)i5);
												v = GCDm(u, (USL)i6);
												v1 = GCDm(v, (USL)i7);
												v2 = GCDm(v1, (USL)i8);
												v3 = GCDm(v2, (USL)i9);
												v4 = GCDm(v3, (USL)i10);
												v5 = GCDm(v4, (USL)i11);
												if (v5 == 1)
												{
													printf("(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11)=(%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u)\n", i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11);
													MATI1 = BUILDMATI(m, 1);
													MATI1->V[0][0] = CHANGE((USL)i1);
													MATI1->V[1][0] = CHANGE((USL)i2);
													MATI1->V[2][0] = CHANGE((USL)i3);
													MATI1->V[3][0] = CHANGE((USL)i4);
													MATI1->V[4][0] = CHANGE((USL)i5);
													MATI1->V[5][0] = CHANGE((USL)i6);
													MATI1->V[6][0] = CHANGE((USL)i7);
													MATI1->V[7][0] = CHANGE((USL)i8);
													MATI1->V[8][0] = CHANGE((USL)i9);
													MATI1->V[9][0] = CHANGE((USL)i10);
													MATI1->V[10][0] = CHANGE((USL)i11);
													BB = LLLGCD(MATI1, &GCD, m, m1, n1);
													FREEMATI(MATI1);
													FREEMPI(GCD);
													MATI3 = BUILDMATI(1, m);
													for (l = 0; l < m; l++)
														MATI3->V[0][l] = COPYI(BB->V[m - 1][l]);
													T1 = LENGTHSQRI(MATI3, 0);
													Q = BB;
													XX = (MPI **)mmalloc((USL)(p * sizeof(MPI *)));
													for (t = 0; t < p; t++)
														XX[t] = ZEROI();
													while (1)
													{
														M = SHORTESTT0(Q, &X);
														for (t = 0; t < p; t++)
														{
															Temp = XX[t];
															XX[t] = ADDI(XX[t], X[t]);
															FREEMPI(Temp);
															FREEMPI(X[t]);
														}
														ffree((char *)X, p * sizeof(MPI *));
														if (M == NULL)
															break;
														else
														{
															for (l = 0; l < Q->C; l++)
															{
																FREEMPI(Q->V[m - 1][l]);
																Q->V[m - 1][l] = COPYI(M->V[0][l]);
															}
															FREEMATI(MATI3);
															MATI3 = M;
														}
													}
													T2 = LENGTHSQRI(MATI3, 0);
													if (RSV(T2, T1) == -1)
													{
														outfile = fopen(buff, "a");
														fprintf(outfile, "(%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%u): ", i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11);
														fprintf(outfile, "b[%u]", m);
														for (t = 0; t < p; t++)
														{
															e = XX[t]->S;
															if (e == -1)
															{
																fprintf(outfile, "+");
																Temp = MINUSI(XX[t]);
																if (!EQONEI(Temp))
																	FPRINTI(outfile, Temp);
																fprintf(outfile, "b[%u]", t + 1);
																FREEMPI(Temp);
															}
															if (e == 1)
															{
																fprintf(outfile, "-");
																if (!EQONEI(XX[t]))
																	FPRINTI(outfile, XX[t]);
																fprintf(outfile, "b[%u]", t + 1);
															}
														}
														fprintf(outfile, "\n");
														fclose(outfile);
													}
													for (t = 0; t < p; t++)
														FREEMPI(XX[t]);
													ffree((char *)XX, p * sizeof(MPI *));
													FREEMATI(BB);
													FREEMATI(MATI3);
													FREEMPI(T1);
													FREEMPI(T2);
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	return;
}

void GCD33()
/* We compute an optimal multipliers for a[1],a[2],a[3]
* in the range p1<=a[1]<=q1 etc. and rel prime.
*/
{
	USI l, m, p1, q1, p2, q2, p3, q3, m1, n1, t;
	USI i1, i2, i3, p;
	int e, s;
	USL x, y;
	MPI *GCD, *T1, *T2, **XX, **X, *Temp, *DIFF;
	MPMATI *MATI1, *BB, *MATI3, *Q, *M;
	char buff[20];
	FILE *outfile;


	printf("Enter alpha=m/n,  1/4 < m/n <= 1:  (ie enter  m  and n)");
	s = scanf("%u%u", &m1, &n1);
	strcpy(buff, "gcd33.out");
	if (access(buff, R_OK) == 0)
		unlink(buff);
	printf("Enter  the ranges pi,qi; 2 <= pi < qi < 2^32, i = 1,..3: ");
	s = scanf("%u%u%u%u%u%u", &p1, &q1, &p2, &q2, &p3, &q3);
	Flush();
	m = 3;
	p = m - 1;
	for (i1 = p1; i1 <= q1; i1++)
	{
		for (i2 = p2; i2 <= q2; i2++)
		{
			printf("(i1,i2)=(%u,%u)\n", i1, i2);
			for (i3 = p3; i3 <= q3; i3++)
			{
				x = GCDm((USL)i1, (USL)i2);
				y = GCDm(x, (USL)i3);
				if (y == 1)
				{
					MATI1 = BUILDMATI(m, 1);
					MATI1->V[0][0] = CHANGE((USL)i1);
					MATI1->V[1][0] = CHANGE((USL)i2);
					MATI1->V[2][0] = CHANGE((USL)i3);
					BB = LLLGCD0M(MATI1, &GCD, m, m1, n1);
					FREEMATI(MATI1);
					FREEMPI(GCD);
					MATI3 = BUILDMATI(1, m);
					for (l = 0; l < m; l++)
						MATI3->V[0][l] = COPYI(BB->V[m - 1][l]);
					T1 = LENGTHSQRI(MATI3, 0);
					Q = BB;
					XX = (MPI **)mmalloc((USL)(p * sizeof(MPI *)));
					for (t = 0; t < p; t++)
						XX[t] = ZEROI();
					while (1)
					{
						M = SHORTESTT0(Q, &X);
						for (t = 0; t < p; t++)
						{
							Temp = XX[t];
							XX[t] = ADDI(XX[t], X[t]);
							FREEMPI(Temp);
							FREEMPI(X[t]);
						}
						ffree((char *)X, p * sizeof(MPI *));
						if (M == NULL)
							break;
						else
						{
							for (l = 0; l < Q->C; l++)
							{
								FREEMPI(Q->V[m - 1][l]);
								Q->V[m - 1][l] = COPYI(M->V[0][l]);
							}
							FREEMATI(MATI3);
							MATI3 = M;
						}
					}
					T2 = LENGTHSQRI(MATI3, 0);
					if (RSV(T2, T1) == -1)
					{
						outfile = fopen(buff, "a");
						fprintf(outfile, "(%u,%u,%u): ", i1, i2, i3);
						fprintf(outfile, "b[%u]", m);
						for (t = 0; t < p; t++)
						{
							e = XX[t]->S;
							if (e == -1)
							{
								fprintf(outfile, "+");
								Temp = MINUSI(XX[t]);
								if (!EQONEI(Temp))
									FPRINTI(outfile, Temp);
								fprintf(outfile, "b[%u]", t + 1);
								FREEMPI(Temp);
							}
							if (e == 1)
							{
								fprintf(outfile, "-");
								if (!EQONEI(XX[t]))
									FPRINTI(outfile, XX[t]);
								fprintf(outfile, "b[%u]", t + 1);
							}
						}
						fprintf(outfile, ": ");
						DIFF = SUB0I(T1, T2);
						FPRINTI(outfile, DIFF);
						FREEMPI(DIFF);
						fprintf(outfile, "\n");
						fclose(outfile);
					}
					for (t = 0; t < p; t++)
						FREEMPI(XX[t]);
					ffree((char *)XX, p * sizeof(MPI *));
					FREEMATI(BB);
					FREEMATI(MATI3);
					FREEMPI(T1);
					FREEMPI(T2);
				}
			}
		}
	}
	return;
}

MPMATI *LLLGCD0M(MPMATI *DD, MPI **Aptr, USI m, USI m1, USI n1)
/*
* Input: an m x 1 vector of positive MPI's.
* Output: *Aptr = gcd of the vector of MPI's. Also we return a small set of
* multipliers using the LLL method of Havas and Matthews.
* Normally m1=1,n1=1.
* Also returns matrix B of the LLL extended gcd algorithm of
* Havas-Majewski-Matthews is returned.
* 30/1/97.
* SHORTER  is the shortest length squared all the m short vectors
* returned.
* The matrix returned is the one LLLGCD0 returns, but with the shortest
* of the X[K] as the last row. This is then used in GCD4() etc.
*/
{
	unsigned int i, j, k, l, flag, p, t;
	MPI **A, **D, *X, *Y, *Z, *Tmp, *R, *M1, *N1;
	MPI **XX, *temp, *mu, *SUM, *SHORTER;
	MPI *temp1, *T, *tempI1, *tempI2;
	MPMATI *L, *B, *MATI1, *MATI2, *MATI;
	MPR *tempR1;
	int r, kk, K, c, COUNT, count;

	B = IDENTITYI(m);
	A = (MPI **)mmalloc(m * sizeof(MPI *));
	for (i = 0; i < m; i++)
		A[i] = COPYI(DD->V[i][0]);
	D = (MPI **)mmalloc((1 + m) * sizeof(MPI *));
	for (i = 0; i <= m; i++)
		D[i] = ONEI();

	L = ZEROMNI(m, m);
	M1 = CHANGE(m1);
	N1 = CHANGE(n1);

	k = 2;
	while (k <= m)
	{
		REDUCE(k, k - 1, &L, &B, D, A);
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
		if (!EQZEROI(A[k - 2]) || (EQZEROI(A[k - 1]) && (RSV(Y, Z) == 1)))
		{
			for (i = k + 1; i <= m; i++)
				SWAP1(i, k, &L, D);
			SWAP2(k, &B, &L, A);
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
		{
			FREEMPI(X);
			FREEMPI(Y);
			for (l = k - 2; l >= 1; l--)
				REDUCE(k, l, &L, &B, D, A);
			k++;
		}
		FREEMPI(Z);
	}
	FREEMPI(M1);
	FREEMPI(N1);
	*Aptr = COPYI(A[m - 1]);
	if ((*Aptr)->S == -1)
	{
		LLLGCDMINUS(&B, &L, m - 1);
		(*Aptr)->S = -((*Aptr)->S);
	}
	/*
	printf("L = \n");
	PRINTMATI(0,L->R-1,0,L->C-1,L);
	for (i = 0; i <= m; i++)
	{
	printf("D[%u] = ", i);PRINTI(D[i]);printf(", ");
	}
	printf("\n");
	*/

	XX = (MPI **)mmalloc(m * sizeof(MPI *));
	for (i = 0; i < m - 1; i++)
		XX[i] = ZEROI();
	XX[m - 1] = ONEI();
	COUNT = m - 2;
	SHORTER = ZEROI();
	MATI = BUILDMATI(1, m);
	for (j = 0; j < m; j++)
		MATI->V[0][j] = ZEROI();
	for (K = m - 2; K >= 0; K--)
	{
		flag = 1;
		for (k = 0; k < m - 1; k++)
		{
			FREEMPI(XX[k]);
			XX[k] = ZEROI();
		}
		for (kk = K; kk >= 0; kk--)
		{
			mu = COPYI(L->V[m - 1][kk]);
			if (flag)
			{
				flag = 0;
				if (mu->S == 0)
				{
					FREEMPI(mu);
					continue;
				}
				else
				{
					tempI1 = ADDI(L->V[m - 1][kk], L->V[m - 1][kk]);
					tempI2 = ABSI(tempI1);
					t = EQUALI(tempI2, D[kk]);
					FREEMPI(tempI1);
					FREEMPI(tempI2);
					if (t)
					{
						FREEMPI(mu);
						continue;
					}
					temp = XX[kk];
					if (mu->S == 1)
						XX[kk] = MINUS_ONEI();
					else
						XX[kk] = ONEI();
					FREEMPI(temp);
				}
			}
			else
			{
				SUM = COPYI(mu);
				for (r = m - 2; r >= kk + 1; r--)
				{
					temp = SUM;
					tempI2 = MULTI(XX[r], L->V[r][kk]);
					SUM = ADDI(SUM, tempI2);
					FREEMPI(tempI2);
					FREEMPI(temp);
				}
				tempR1 = RATIOI(SUM, D[kk + 1]);
				FREEMPI(SUM);
				temp1 = ABS_NEAREST_INTR(tempR1);
				FREEMPR(tempR1);
				temp = XX[kk];
				XX[kk] = MINUSI(temp1);
				FREEMPI(temp1);
				FREEMPI(temp);
			}
			FREEMPI(mu);
		}
		MATI1 = BUILDMATI(1, m);
		for (j = 0; j < m; j++)
			MATI1->V[0][j] = COPYI(XX[j]);
		MATI2 = MULTMATI(MATI1, B);
		p = m - 1;
		T = LENGTHSQRI(MATI2, 0);
		if (K == m - 2)
		{
			FREEMPI(SHORTER);
			SHORTER = COPYI(T);
			FREEMATI(MATI);
			MATI = COPYMATI(MATI2);
		}
		else
		{
			c = RSV(T, SHORTER);
			if (c < 0)
			{
				FREEMPI(SHORTER);
				SHORTER = COPYI(T);
				FREEMATI(MATI);
				MATI = COPYMATI(MATI2);
				COUNT = K;
			}
		}
		FREEMPI(T);
		FREEMATI(MATI1);
		FREEMATI(MATI2);
	}
	T = LENGTHSQRI(B, m - 1);
	c = RSV(T, SHORTER);
	if (c < 0)
	{
		FREEMPI(SHORTER);
		SHORTER = COPYI(T);
		for (i = 0; i < m; i++)
		{
			FREEMPI(MATI->V[0][i]);
			MATI->V[0][i] = COPYI(B->V[m - 1][i]);
		}
		COUNT = K;
	}
	FREEMPI(T);
	count = COUNT + 1;
	for (i = 0; i < m; i++)
	{
		FREEMPI(B->V[m - 1][i]);
		B->V[m - 1][i] = COPYI(MATI->V[0][i]);
	}

	FREEMATI(MATI);
	FREEMPI(SHORTER);
	FREEMATI(L);
	for (i = 0; i <= m; i++)
		FREEMPI(D[i]);
	ffree((char *)D, (1 + m) * sizeof(MPI *));
	for (i = 0; i < m; i++)
		FREEMPI(A[i]);
	ffree((char *)A, m * sizeof(MPI *));
	for (i = 0; i < m; i++)
		FREEMPI(XX[i]);
	ffree((char *)XX, m * sizeof(MPI *));
	return (B);
}

MPMATI *LLLGCDL(MPMATI *DD, MPI **Aptr, USI m, USI m1, USI n1, MPMATR **LMATRIX)
/*
* Input: an m x 1 vector of positive MPI's.
* Output: *Aptr = gcd of the vector of MPI's. Also we return a small
* set of multipliers using the LLL method of Havas and Matthews.
* Also output the L=[mu[i][j]] matrix.
* matrix B of the LLL extended gcd algorithm of
* Havas-Majewski-Matthews is returned.
*/
{
	unsigned int i, j, k, l;
	MPI **A, **D, *X, *Y, *Z, *Tmp, *R, *M1, *N1;
	MPMATI *L, *B;

	B = IDENTITYI(m);
	A = (MPI **)mmalloc(m * sizeof(MPI *));
	for (i = 0; i < m; i++)
		A[i] = COPYI(DD->V[i][0]);
	D = (MPI **)mmalloc((1 + m) * sizeof(MPI *));
	for (i = 0; i <= m; i++)
		D[i] = ONEI();

	L = ZEROMNI(m, m);
	M1 = CHANGE(m1);
	N1 = CHANGE(n1);

	k = 2;
	while (k <= m)
	{
		REDUCE(k, k - 1, &L, &B, D, A);
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
		if (!EQZEROI(A[k - 2]) || (EQZEROI(A[k - 1]) && (RSV(Y, Z) == 1)))
		{
			for (i = k + 1; i <= m; i++)
				SWAP1(i, k, &L, D);
			SWAP2(k, &B, &L, A);
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
		{
			FREEMPI(X);
			FREEMPI(Y);
			for (l = k - 2; l >= 1; l--)
				REDUCE(k, l, &L, &B, D, A);
			k++;
		}
		FREEMPI(Z);
	}
	FREEMPI(M1);
	FREEMPI(N1);
	*Aptr = COPYI(A[m - 1]);
	if ((*Aptr)->S == -1)
	{
		LLLGCDMINUS(&B, &L, m - 1);
		(*Aptr)->S = -((*Aptr)->S);
	}
	*LMATRIX = ZEROMNR(m, m);
	for (i = 0; i < m; i++)
		for (j = 0; j < m; j++)
		{
			if (i > j)
			{
				FREEMPR(elt(*LMATRIX, i, j));
				elt(*LMATRIX, i, j) = RATIOI(L->V[i][j], D[j + 1]);
			}

		}
	/*
	printf("*LMATRIX = \n");
	PRINTMATR(0,(*LMATRIX)->R-1,0,(*LMATRIX)->C-1,*LMATRIX);
	printf("L = \n");
	PRINTMATI(0,L->R-1,0,L->C-1,L);
	for (i = 0; i <= m; i++)
	{
	printf("D[%u] = ", i);PRINTI(D[i]);printf(", ");
	}
	printf("\n");
	*/
	FREEMATI(L);
	for (i = 0; i <= m; i++)
		FREEMPI(D[i]);
	for (i = 0; i < m; i++)
		FREEMPI(A[i]);
	ffree((char *)D, (1 + m) * sizeof(MPI *));
	ffree((char *)A, m * sizeof(MPI *));
	return (B);
}

void GCDCONJECTURE5()
/* We test our LLLGCD conjecture for a[1],a[2],a[3],a[4],a[5]
* in the range p1<=a[1]<=q1 etc. and rel prime.
*/
{
	USI p1, q1, p2, q2, p3, q3, p4, q4, p5, q5;
	USI i1, i2, i3, i4, i5, p, record = 1;
	USL x, y, z, u;
	USI m1, n1, m, n, i, j, count;
	int r, e, K, s;
	MPMATI *MATI1, *MATI2, *Q, *M, *MATI3;
	MPI *A, **XX, **X, *Temp, ***COEFF, *tempI;
	MPMATR *LMATRIX;
	MPR **SIGMA, *SUMR, *tR1, *tR2, *tR3, *tempR, *tempR1;
	char buff[20];
	FILE *outfile;

	printf("Enter alpha=m/n,  1/4 < m/n <= 1:  (ie enter  m  and n)");
	s = scanf("%u%u", &m1, &n1);
	strcpy(buff, "conjecture5.out");
	if (access(buff, R_OK) == 0)
		unlink(buff);
	printf("Enter  the ranges pi,qi; 2 <= pi < qi < 2^32, i = 1,...,5: ");
	s = scanf("%u%u%u%u%u%u%u%u%u%u", &p1, &q1, &p2, &q2, &p3, &q3, &p4, &q4, &p5, &q5);
	Flush();
	m = 5;
	p = m - 1;
	for (i1 = p1; i1 <= q1; i1++)
	{
		for (i2 = p2; i2 <= q2; i2++)
		{
			for (i3 = p3; i3 <= q3; i3++)
			{
				printf("(i1,i2,i3)=(%u,%u,%u)\n", i1, i2, i3);
				for (i4 = p4; i4 <= q4; i4++)
				{
					for (i5 = p5; i5 <= q5; i5++)
					{
						/*printf("(i1,i2,i3,i4,i5)=(%u,%u,%u,%u,%u)\n", i1,i2,i3,i4,i5);*/
						x = GCDm((USL)i1, (USL)i2);
						y = GCDm(x, (USL)i3);
						z = GCDm(y, (USL)i4);
						u = GCDm(z, (USL)i5);
						if (u == 1)
						{
							MATI1 = BUILDMATI(m, 1);
							MATI1->V[0][0] = CHANGE((USL)i1);
							MATI1->V[1][0] = CHANGE((USL)i2);
							MATI1->V[2][0] = CHANGE((USL)i3);
							MATI1->V[3][0] = CHANGE((USL)i4);
							MATI1->V[4][0] = CHANGE((USL)i5);

							m = MATI1->R;
							MATI2 = LLLGCDL(MATI1, &A, m, m1, n1, &LMATRIX);
							FREEMATI(MATI1);
							FREEMPI(A);
							MATI3 = BUILDMATI(1, m);
							for (i = 0; i < m; i++)
								MATI3->V[0][i] = COPYI(MATI2->V[m - 1][i]);
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

							/*SHORTEST(Q, XX, 2, 0);*/
							COEFF = SHORTESTX(Q, XX, &count);
							if (count>record)
							{
								record = count;
								outfile = fopen(buff, "a");
								printf("(i1,i2,i3,i4,i5)=(%u,%u,%u,%u,%u): record=%u\n", i1, i2, i3, i4, i5, record);
								fprintf(outfile, "(i1,i2,i3,i4,i5)=(%u,%u,%u,%u,%u): record=%u\n", i1, i2, i3, i4, i5, record);
								fclose(outfile);
							}
							for (j = 0; j < p; j++)
								FREEMPI(XX[j]);
							ffree((char *)XX, p * sizeof(MPI *));
							SIGMA = (MPR **)mmalloc(p * sizeof(MPR *));
							for (j = 0; j < count; j++)
							{
								/*printf("SIGMA[%u]=", j);*/
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
									/*PRINTR(SIGMA[K]);
									printf(" ");*/
									tempR = BUILDMPR();
									tempR->N = COPYI(COEFF[j][K]);
									tempR->D = ONEI();
									tempR1 = ADDR(tempR, SIGMA[K]);
									tempI = ABSI(tempR1->N);
									e = RSV(tempI, tempR1->D);
									FREEMPR(tempR);
									FREEMPR(tempR1);
									FREEMPI(tempI);
									if (e >= 0)
									{
										fprintf(stderr, "conjecture false:j = %u, K = %d\n", j, K + 1);
										outfile = fopen(buff, "a"), count;
										fprintf(outfile, "(i1,i2,i3,i4,i5)=(%u,%u,%u,%u,%u): ", i1, i2, i3, i4, i5);
										fprintf(outfile, "conjecture false: j = %u, K = %d\n", j, K + 1);
										fclose(outfile);
										exit(1);
									}
								}
								/*printf("\n");*/
								for (K = 0; K < p; K++)
									FREEMPR(SIGMA[K]);
							}
							ffree((char *)SIGMA, p * sizeof(MPR *));

							ffree((char *)COEFF, GCD_MAX * sizeof(MPI **));
							for (j = 0; j < count; j++)
							{
								ffree((char *)COEFF[j], p * sizeof(MPI *));
								for (i = 0; i < p; i++)
									FREEMPI(COEFF[j][i]);
							}
							FREEMATI(MATI3);
							FREEMATI(Q);
							FREEMATR(LMATRIX);
						}
					}
				}
			}
		}
	}
	return;
}

void GCDCONJECTURE4()
/* We test our LLLGCD conjecture for a[1],a[2],a[3],a[4].
* in the range p1<=a[1]<=q1 etc. and rel prime.
*/
{
	USI p1, q1, p2, q2, p3, q3, p4, q4;
	USI i1, i2, i3, i4, p, record = 1;
	USL x, y, z;
	USI m1, n1, m, n, i, j, count;
	int r, e, K, s;
	MPMATI *MATI1, *MATI2, *Q, *M, *MATI3;
	MPI *A, **XX, **X, *Temp, ***COEFF, *tempI;
	MPMATR *LMATRIX;
	MPR **SIGMA, *SUMR, *tR1, *tR2, *tR3, *tempR, *tempR1;
	char buff[20];
	FILE *outfile;

	printf("Enter alpha=m/n,  1/4 < m/n <= 1:  (ie enter  m  and n)");
	s = scanf("%u%u", &m1, &n1);
	strcpy(buff, "conjecture4.out");
	if (access(buff, R_OK) == 0)
		unlink(buff);
	printf("Enter  the ranges pi,qi; 2 <= pi < qi < 2^32, i = 1,..4: ");
	s = scanf("%u%u%u%u%u%u%u%u", &p1, &q1, &p2, &q2, &p3, &q3, &p4, &q4);
	Flush();
	m = 4;
	p = m - 1;
	for (i1 = p1; i1 <= q1; i1++)
	{
		for (i2 = p2; i2 <= q2; i2++)
		{
			for (i3 = p3; i3 <= q3; i3++)
			{
				for (i4 = p4; i4 <= q4; i4++)
				{
					printf("(i1,i2,i3,i4)=(%u,%u,%u,%u)\n", i1, i2, i3, i4);
					x = GCDm((USL)i1, (USL)i2);
					y = GCDm(x, (USL)i3);
					z = GCDm(y, (USL)i4);
					if (z == 1)
					{
						printf("(i1,i2,i3,i4)=(%u,%u,%u,%u)\n", i1, i2, i3, i4);
						MATI1 = BUILDMATI(m, 1);
						MATI1->V[0][0] = CHANGE((USL)i1);
						MATI1->V[1][0] = CHANGE((USL)i2);
						MATI1->V[2][0] = CHANGE((USL)i3);
						MATI1->V[3][0] = CHANGE((USL)i4);

						m = MATI1->R;
						MATI2 = LLLGCDL(MATI1, &A, m, m1, n1, &LMATRIX);
						FREEMATI(MATI1);
						FREEMPI(A);
						MATI3 = BUILDMATI(1, m);
						for (i = 0; i < m; i++)
							MATI3->V[0][i] = COPYI(MATI2->V[m - 1][i]);
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

						COEFF = SHORTESTX(Q, XX, &count);
						if (count>record)
						{
							record = count;
							outfile = fopen(buff, "a");
							printf("(i1,i2,i3,i4)=(%u,%u,%u,%u): record=%u\n", i1, i2, i3, i4, record);
							fprintf(outfile, "(i1,i2,i3,i4)=(%u,%u,%u,%u): record=%u\n", i1, i2, i3, i4, record);
							fclose(outfile);
						}
						for (j = 0; j < p; j++)
							FREEMPI(XX[j]);
						ffree((char *)XX, p * sizeof(MPI *));
						SIGMA = (MPR **)mmalloc(p * sizeof(MPR *));
						for (j = 0; j < count; j++)
						{
							/*printf("SIGMA[%u]=", j);*/
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
								/*PRINTR(SIGMA[K]);
								printf(" ");*/
								tempR = BUILDMPR();
								tempR->N = COPYI(COEFF[j][K]);
								tempR->D = ONEI();
								tempR1 = ADDR(tempR, SIGMA[K]);
								tempI = ABSI(tempR1->N);
								e = RSV(tempI, tempR1->D);
								FREEMPR(tempR);
								FREEMPR(tempR1);
								FREEMPI(tempI);
								if (e >= 0)
								{
									fprintf(stderr, "conjecture false:j = %u, K = %d\n", j, K + 1);
									outfile = fopen(buff, "a"), count;
									fprintf(outfile, "(i1,i2,i3,i4)=(%u,%u,%u,%u): ", i1, i2, i3, i4);
									fprintf(outfile, "conjecture false: j = %u, K = %d\n", j, K + 1);
									fclose(outfile);
									exit(1);
								}
							}
							/*printf("\n");*/
							for (K = 0; K < p; K++)
								FREEMPR(SIGMA[K]);
						}
						ffree((char *)SIGMA, p * sizeof(MPR *));

						ffree((char *)COEFF, GCD_MAX * sizeof(MPI **));
						for (j = 0; j < count; j++)
						{
							ffree((char *)COEFF[j], p * sizeof(MPI *));
							for (i = 0; i < p; i++)
								FREEMPI(COEFF[j][i]);
						}
						FREEMATI(MATI3);
						FREEMATI(Q);
						FREEMATR(LMATRIX);
					}
				}
			}
		}
	}
	return;
}

void GCD3()
/* We print out all shortest multipliers for a[1],a[2],a[3].
* in the range p1<=a[1]<=q1 etc. and rel prime.
* The output is sent to gcd3.out.
*/
{
	USI p1, q1, p2, q2, p3, q3;
	USI i1, i2, i3, p;
	USL x, y;
	int s;
	USI m1, n1, m, n, i, j;
	MPMATI *MATI1, *MATI2, *Q, *M, *MATI3;
	MPI *A, **XX, **X, *Temp;
	MPMATR *LMATRIX;
	char buff[20];
	FILE *outfile;

	printf("Enter alpha=m/n,  1/4 < m/n <= 1:  (ie enter  m  and n)");
	s = scanf("%u%u", &m1, &n1);
	strcpy(buff, "gcd3.out");
	if (access(buff, R_OK) == 0)
		unlink(buff);
	printf("Enter  the ranges pi,qi; 2 <= pi < qi < 2^32, i = 1,..3: ");
	s = scanf("%u%u%u%u%u%u", &p1, &q1, &p2, &q2, &p3, &q3);
	Flush();
	m = 3;
	p = m - 1;
	for (i1 = p1; i1 <= q1; i1++)
	{
		for (i2 = p2; i2 <= q2; i2++)
		{
			for (i3 = p3; i3 <= q3; i3++)
			{
				x = GCDm((USL)i1, (USL)i2);
				y = GCDm(x, (USL)i3);
				if (y == 1)
				{
					MATI1 = BUILDMATI(m, 1);
					MATI1->V[0][0] = CHANGE((USL)i1);
					MATI1->V[1][0] = CHANGE((USL)i2);
					MATI1->V[2][0] = CHANGE((USL)i3);

					m = MATI1->R;
					MATI2 = LLLGCDL(MATI1, &A, m, m1, n1, &LMATRIX);
					FREEMATI(MATI1);
					FREEMPI(A);
					MATI3 = BUILDMATI(1, m);
					for (i = 0; i < m; i++)
						MATI3->V[0][i] = COPYI(MATI2->V[m - 1][i]);
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

					outfile = fopen(buff, "a");
					printf("(i1,i2,i3)=(%u,%u,%u):", i1, i2, i3);
					fprintf(outfile, "(i1,i2,i3)=(%u,%u,%u):", i1, i2, i3);
					printf("(MU[21],MU[31],MU[32]=");
					fprintf(outfile, "(MU[21],MU[31],MU[32]=");
					PRINTR(elt(LMATRIX, 1, 0));
					FPRINTR(outfile, elt(LMATRIX, 1, 0));
					printf(",");
					fprintf(outfile, ",");
					PRINTR(elt(LMATRIX, 2, 0));
					FPRINTR(outfile, elt(LMATRIX, 2, 0));
					printf(",");
					fprintf(outfile, ",");
					PRINTR(elt(LMATRIX, 2, 1));
					FPRINTR(outfile, elt(LMATRIX, 2, 1));
					printf(":\n");
					fprintf(outfile, ":\n");
					fclose(outfile);
					SHORTEST(Q, XX, 4, 1);
					for (j = 0; j < p; j++)
						FREEMPI(XX[j]);
					ffree((char *)XX, p * sizeof(MPI *));
					FREEMATI(MATI3);
					FREEMATI(Q);
					FREEMATR(LMATRIX);
				}
			}
		}
	}
	return;
}

void GCDCONJECTURE6()
/* We test our LLLGCD conjecture for a[1],a[2],a[3],a[4],a[5],a[6]
* in the range p1<=a[1]<=q1 etc. and rel prime.
*/
{
	USI p1, q1, p2, q2, p3, q3, p4, q4, p5, q5, p6, q6;
	USI i1, i2, i3, i4, i5, i6, p, record = 1;
	USL x, y, z, u, v;
	USI m1, n1, m, n, i, j, count;
	int r, e, K, s;
	MPMATI *MATI1, *MATI2, *Q, *M, *MATI3;
	MPI *A, **XX, **X, *Temp, ***COEFF, *tempI;
	MPMATR *LMATRIX;
	MPR **SIGMA, *SUMR, *tR1, *tR2, *tR3, *tempR, *tempR1;
	char buff[20];
	FILE *outfile;

	printf("Enter alpha=m/n,  1/4 < m/n <= 1:  (ie enter  m  and n)");
	s = scanf("%u%u", &m1, &n1);
	strcpy(buff, "conjecture6.out");
	if (access(buff, R_OK) == 0)
		unlink(buff);
	printf("Enter  the ranges pi,qi; 2 <= pi < qi < 2^32, i = 1,...,6: ");
	s = scanf("%u%u%u%u%u%u%u%u%u%u%u%u", &p1, &q1, &p2, &q2, &p3, &q3, &p4, &q4, &p5, &q5, &p6, &q6);
	Flush();
	m = 6;
	p = m - 1;
	for (i1 = p1; i1 <= q1; i1++)
	{
		for (i2 = p2; i2 <= q2; i2++)
		{
			for (i3 = p3; i3 <= q3; i3++)
			{
				printf("(i1,i2,i3)=(%u,%u,%u)\n", i1, i2, i3);
				for (i4 = p4; i4 <= q4; i4++)
				{
					for (i5 = p5; i5 <= q5; i5++)
					{
						for (i6 = p6; i6 <= q6; i6++)
						{
							/*printf("(i1,i2,i3,i4,i5,i6)=(%u,%u,%u,%u,%u,%u)\n", i1,i2,i3,i4,i5,i6);*/
							x = GCDm((USL)i1, (USL)i2);
							y = GCDm(x, (USL)i3);
							z = GCDm(y, (USL)i4);
							u = GCDm(z, (USL)i5);
							v = GCDm(u, (USL)i6);
							if (v == 1)
							{
								MATI1 = BUILDMATI(m, 1);
								MATI1->V[0][0] = CHANGE((USL)i1);
								MATI1->V[1][0] = CHANGE((USL)i2);
								MATI1->V[2][0] = CHANGE((USL)i3);
								MATI1->V[3][0] = CHANGE((USL)i4);
								MATI1->V[4][0] = CHANGE((USL)i5);
								MATI1->V[5][0] = CHANGE((USL)i6);

								m = MATI1->R;
								MATI2 = LLLGCDL(MATI1, &A, m, m1, n1, &LMATRIX);
								FREEMATI(MATI1);
								FREEMPI(A);
								MATI3 = BUILDMATI(1, m);
								for (i = 0; i < m; i++)
									MATI3->V[0][i] = COPYI(MATI2->V[m - 1][i]);
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

								COEFF = SHORTESTX(Q, XX, &count);
								if (count>record)
								{
									record = count;
									outfile = fopen(buff, "a");
									printf("(i1,i2,i3,i4,i5,i6)=(%u,%u,%u,%u,%u,%u): record=%u\n", i1, i2, i3, i4, i5, i6, record);
									fprintf(outfile, "(i1,i2,i3,i4,i5,i6)=(%u,%u,%u,%u,%u,%u): record=%u\n", i1, i2, i3, i4, i5, i6, record);
									fclose(outfile);
								}
								for (j = 0; j < p; j++)
									FREEMPI(XX[j]);
								ffree((char *)XX, p * sizeof(MPI *));
								SIGMA = (MPR **)mmalloc(p * sizeof(MPR *));
								for (j = 0; j < count; j++)
								{
									/*printf("SIGMA[%u]=", j);*/
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
										/*PRINTR(SIGMA[K]);
										printf(" ");*/
										tempR = BUILDMPR();
										tempR->N = COPYI(COEFF[j][K]);
										tempR->D = ONEI();
										tempR1 = ADDR(tempR, SIGMA[K]);
										tempI = ABSI(tempR1->N);
										e = RSV(tempI, tempR1->D);
										FREEMPR(tempR);
										FREEMPR(tempR1);
										FREEMPI(tempI);
										if (e >= 0)
										{
											fprintf(stderr, "conjecture false:j = %u, K = %d\n", j, K + 1);
											outfile = fopen(buff, "a"), count;
											fprintf(outfile, "(i1,i2,i3,i4,i5,i6)=(%u,%u,%u,%u,%u,%u): ", i1, i2, i3, i4, i5, i6);
											fprintf(outfile, "conjecture false: j = %u, K = %d\n", j, K + 1);
											fclose(outfile);
											exit(1);
										}
									}
									/*printf("\n");*/
									for (K = 0; K < p; K++)
										FREEMPR(SIGMA[K]);
								}
								ffree((char *)SIGMA, p * sizeof(MPR *));

								ffree((char *)COEFF, GCD_MAX * sizeof(MPI **));
								for (j = 0; j < count; j++)
								{
									ffree((char *)COEFF[j], p * sizeof(MPI *));
									for (i = 0; i < p; i++)
										FREEMPI(COEFF[j][i]);
								}
								FREEMATI(MATI3);
								FREEMATI(Q);
								FREEMATR(LMATRIX);
							}
						}
					}
				}
			}
		}
	}
	return;
}

void GCDCONJECTURE7()
/* We test our LLLGCD conjecture for a[1],...,a[7]
* in the range p1<=a[1]<=q1 etc. and rel prime.
*/
{
	USI p1, q1, p2, q2, p3, q3, p4, q4, p5, q5, p6, q6, p7, q7;
	USI i1, i2, i3, i4, i5, i6, i7, p, record = 1;
	USL x, y, z, u, v, w;
	USI m1, n1, m, n, i, j, count;
	int r, e, K, s;
	MPMATI *MATI1, *MATI2, *Q, *M, *MATI3;
	MPI *A, **XX, **X, *Temp, ***COEFF, *tempI;
	MPMATR *LMATRIX;
	MPR **SIGMA, *SUMR, *tR1, *tR2, *tR3, *tempR, *tempR1;
	char buff[20];
	FILE *outfile;

	printf("Enter alpha=m/n,  1/4 < m/n <= 1:  (ie enter  m  and n)");
	s = scanf("%u%u", &m1, &n1);
	strcpy(buff, "conjecture7.out");
	if (access(buff, R_OK) == 0)
		unlink(buff);
	printf("Enter  the ranges pi,qi; 2 <= pi < qi < 2^32, i = 1,...,7: ");
	s = scanf("%u%u%u%u%u%u%u%u%u%u%u%u%u%u", &p1, &q1, &p2, &q2, &p3, &q3, &p4, &q4, &p5, &q5, &p6, &q6, &p7, &q7);
	Flush();
	m = 7;
	p = m - 1;
	for (i1 = p1; i1 <= q1; i1++)
	{
		for (i2 = p2; i2 <= q2; i2++)
		{
			for (i3 = p3; i3 <= q3; i3++)
			{
				/*	printf("(i1,i2,i3)=(%u,%u,%u)\n", i1,i2,i3);*/
				for (i4 = p4; i4 <= q4; i4++)
				{
					for (i5 = p5; i5 <= q5; i5++)
					{
						for (i6 = p6; i6 <= q6; i6++)
						{
							for (i7 = p7; i7 <= q7; i7++)
							{
								printf("(i1,i2,i3,i4,i5,i6,i7)=(%u,%u,%u,%u,%u,%u,%u)\n", i1, i2, i3, i4, i5, i6, i7);
								x = GCDm((USL)i1, (USL)i2);
								y = GCDm(x, (USL)i3);
								z = GCDm(y, (USL)i4);
								u = GCDm(z, (USL)i5);
								v = GCDm(u, (USL)i6);
								w = GCDm(v, (USL)i7);
								if (w == 1)
								{
									MATI1 = BUILDMATI(m, 1);
									MATI1->V[0][0] = CHANGE((USL)i1);
									MATI1->V[1][0] = CHANGE((USL)i2);
									MATI1->V[2][0] = CHANGE((USL)i3);
									MATI1->V[3][0] = CHANGE((USL)i4);
									MATI1->V[4][0] = CHANGE((USL)i5);
									MATI1->V[5][0] = CHANGE((USL)i6);
									MATI1->V[6][0] = CHANGE((USL)i7);

									m = MATI1->R;
									MATI2 = LLLGCDL(MATI1, &A, m, m1, n1, &LMATRIX);
									FREEMATI(MATI1);
									FREEMPI(A);
									MATI3 = BUILDMATI(1, m);
									for (i = 0; i < m; i++)
										MATI3->V[0][i] = COPYI(MATI2->V[m - 1][i]);
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

									COEFF = SHORTESTX(Q, XX, &count);
									if (count>record)
									{
										record = count;
										outfile = fopen(buff, "a");
										printf("(i1,i2,i3,i4,i5,i6,i7)=(%u,%u,%u,%u,%u,%u,%u): record=%u\n", i1, i2, i3, i4, i5, i6, i7, record);
										fprintf(outfile, "(i1,i2,i3,i4,i5,i6,i7)=(%u,%u,%u,%u,%u,%u,%u): record=%u\n", i1, i2, i3, i4, i5, i6, i7, record);
										fclose(outfile);
									}
									for (j = 0; j < p; j++)
										FREEMPI(XX[j]);
									ffree((char *)XX, p * sizeof(MPI *));
									SIGMA = (MPR **)mmalloc(p * sizeof(MPR *));
									for (j = 0; j < count; j++)
									{
										/*printf("SIGMA[%u]=", j);*/
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
											/*PRINTR(SIGMA[K]);
											printf(" ");*/
											tempR = BUILDMPR();
											tempR->N = COPYI(COEFF[j][K]);
											tempR->D = ONEI();
											tempR1 = ADDR(tempR, SIGMA[K]);
											tempI = ABSI(tempR1->N);
											e = RSV(tempI, tempR1->D);
											FREEMPR(tempR);
											FREEMPR(tempR1);
											FREEMPI(tempI);
											if (e >= 0)
											{
												fprintf(stderr, "conjecture false:j = %u, K = %d\n", j, K + 1);
												outfile = fopen(buff, "a"), count;
												fprintf(outfile, "(i1,i2,i3,i4,i5,i6,i7)=(%u,%u,%u,%u,%u,%u,%u): ", i1, i2, i3, i4, i5, i6, i7);
												fprintf(outfile, "conjecture false: j = %u, K = %d\n", j, K + 1);
												fclose(outfile);
												exit(1);
											}
										}
										/*printf("\n");*/
										for (K = 0; K < p; K++)
											FREEMPR(SIGMA[K]);
									}
									ffree((char *)SIGMA, p * sizeof(MPR *));

									ffree((char *)COEFF, GCD_MAX * sizeof(MPI **));
									for (j = 0; j < count; j++)
									{
										ffree((char *)COEFF[j], p * sizeof(MPI *));
										for (i = 0; i < p; i++)
											FREEMPI(COEFF[j][i]);
									}
									FREEMATI(MATI3);
									FREEMATI(Q);
									FREEMATR(LMATRIX);
								}
							}
						}
					}
				}
			}
		}
	}
	return;
}

void GCDCONJECTUREM()
/* We test our LLLGCD conjecture for n random mtuples.
*/
{
	USI p, record = 1, power, g, trial_number;
	USI m1, n1, m, n, i, j, k, count;
	int r, e, K, s;
	MPMATI *MATI1, *MATI2, *Q, *M, *MATI3;
	MPI *A, **XX, **X, *Temp, ***COEFF, *tempI;
	MPI *BASE, *MULTIPLIER, *SEED;
	MPMATR *LMATRIX;
	MPR **SIGMA, *SUMR, *tR1, *tR2, *tR3, *tempR, *tempR1;
	char buff[20];
	FILE *outfile;

	printf("Enter alpha=m/n,  1/4 < m/n <= 1:  (ie enter  m  and n)");
	s = scanf("%u%u", &m1, &n1);
	strcpy(buff, "conjecturem.out");
	if (access(buff, R_OK) == 0)
		unlink(buff);
	printf("Enter the number m of integers: ");
	s = scanf("%u", &m);
	p = m - 1;
	printf("enter the number of trials: ");
	s = scanf("%u", &trial_number);
	printf("Enter the power of 2: ");
	s = scanf("%u", &power);
	Flush();
	BASE = POWER_I((long)2, power);
	printf("Enter the (odd) seed > 1): ");
	SEED = INPUTI(&g);
	printf("Enter the multiplier (4k+1, k odd) (eg 69069): ");
	MULTIPLIER = INPUTI(&g);
	for (i = 0; i < trial_number; i++)
	{
		printf("doing i = %u\n", i);
		MATI1 = BUILDMATI(m, 1);
		for (j = 0; j < m; j++)
		{
			MATI1->V[j][0] = SEED;
			SEED = RANDOMI(SEED, MULTIPLIER, BASE);
		}
		MATI2 = LLLGCDL(MATI1, &A, m, m1, n1, &LMATRIX);
		FREEMPI(A);
		MATI3 = BUILDMATI(1, m);
		for (j = 0; j < m; j++)
			MATI3->V[0][j] = COPYI(MATI2->V[m - 1][j]);
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

		COEFF = SHORTESTX(Q, XX, &count);
		if (count>record)
		{
			record = count;
			outfile = fopen(buff, "a");
			printf("record=%u:\n ", record);
			fprintf(outfile, "record=%u\n: ", record);
			PRINTMATI(0, MATI1->R - 1, 0, MATI1->C - 1, MATI1);
			FPRINTMATI(outfile, 0, MATI1->R - 1, 0, MATI1->C - 1, MATI1);
			fclose(outfile);
		}
		for (j = 0; j < p; j++)
			FREEMPI(XX[j]);
		ffree((char *)XX, p * sizeof(MPI *));
		SIGMA = (MPR **)mmalloc(p * sizeof(MPR *));
		for (j = 0; j < count; j++)
		{
			/*printf("SIGMA[%u]=", j);*/
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
				/*PRINTR(SIGMA[K]);
				printf(" ");*/
				tempR = BUILDMPR();
				tempR->N = COPYI(COEFF[j][K]);
				tempR->D = ONEI();
				tempR1 = ADDR(tempR, SIGMA[K]);
				tempI = ABSI(tempR1->N);
				e = RSV(tempI, tempR1->D);
				FREEMPR(tempR);
				FREEMPR(tempR1);
				FREEMPI(tempI);
				if (e >= 0)
				{
					outfile = fopen(buff, "a"), count;
					fprintf(stderr, "conjecture false:j = %u, K = %d\n", j, K + 1);
					printf("conjecture false: j = %u, K = %d\n", j, K + 1);
					fprintf(outfile, "conjecture false:j = %u, K = %d\n", j, K + 1);
					PRINTMATI(0, MATI1->R - 1, 0, MATI1->C - 1, MATI1);
					FPRINTMATI(outfile, 0, MATI1->R - 1, 0, MATI1->C - 1, MATI1);
					fclose(outfile);
					/*				exit (1);*/
				}
			}
			/*printf("\n");*/
			for (K = 0; K < p; K++)
				FREEMPR(SIGMA[K]);
		}
		FREEMATI(MATI1);
		ffree((char *)SIGMA, p * sizeof(MPR *));

		for (j = 0; j < count; j++)
		{
			for (k = 0; k < p; k++)
				FREEMPI(COEFF[j][k]);
			ffree((char *)COEFF[j], p * sizeof(MPI *));
		}
		ffree((char *)COEFF, GCD_MAX * sizeof(MPI **));
		FREEMATI(MATI3);
		FREEMATI(Q);
		FREEMATR(LMATRIX);
	}
	FREEMPI(SEED);
	FREEMPI(BASE);
	FREEMPI(MULTIPLIER);
	return;
}
