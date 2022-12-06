/* reduce.c */
#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_NONSTDC_NO_DEPRECATE
#include <stdio.h>
#include <stdlib.h>
#include "integer.h"
#include "fun.h"

/*
MPI *FACTORQ(Q, factored, base, M1, FBASE, x, r1, r2, r, SCAN)
*/
/*
* produces an exponent vector (mod 2) which updates M2[factored][*]
* and returns a cofactor.
*/
/*
unsigned int factored, base[], **M1, FBASE, r1[], r2[], r;
MPI *Q;
int x, SCAN[];
{
unsigned int i, s, R;
int t, j;
MPI *T, *Tmp;


INTSETUI(M1[factored], 2 + (USL)FBASE, 0);
Tmp = ABSI(Q);
t = x % 2;
R = (x < 0 && t) ? (unsigned int)t + 2 : (unsigned int)t;
if (R == r)
{
s = 0;
while (((Tmp->V[0]) & 1) == 0)
{
s++;
T = INT0_(Tmp, UL2);
FREEMPI(Tmp);
Tmp = T;
}
M1[factored][0] = s & 1;
if (M1[factored][0])
SCAN[factored] = 0;
}
for (i = 1; i <= FBASE; i++)
{
t = x % (int)(base[i]);
R = (x < 0 && t) ? (unsigned int)t + base[i] : (unsigned int)t;
if (R == r1[i] || R == r2[i])
{
s = 0;
while (MODINT0_(Tmp, base[i], &T) == 0)
{
s++;
FREEMPI(Tmp);
Tmp = T;
}
FREEMPI(T);
M1[factored][i] = s & 1;
if (M1[factored][i])
SCAN[factored] = i;
}
}
if (Q->S == -1)
{
M1[factored][1 + FBASE] = 1;
SCAN[factored] = 1 + FBASE;
}
return (Tmp);
}
*/

MPI *FACTORQ1(MPI *Q, USI factored, USL base[], USL **M1, USI FBASE, long x, USL r1[], USL r2[], USL r, int INDEX[], USI CTR[])
/*
* produces an exponent vector (mod 2) which updates M2[factored][*]
* and returns a cofactor. Used in Peter Adams' version MPQS1().
*/
{
	unsigned int i;
	long t;
	MPI *T, *Tmp;
	unsigned long *m1ptr, word, R, sp2, s;
	unsigned int ctr, index;

	m1ptr = M1[factored];
	ctr = 0;
	word = 0;
	index = 0;
	Tmp = ABSI(Q);
	t = x % 2L;
	R = (x < 0 && (USL)t) ? (USL)(t + 2) : (USL)t;
	if (R == r)
	{
		s = 0;
		while (((Tmp->V[0]) & (USL)1) == 0)
		{
			s++;
			T = INT0_(Tmp, (USL)2);
			FREEMPI(Tmp);
			Tmp = T;
		}
		word = (s & (USL)1);
		if (word)
		{
			INDEX[factored] = (int)index;
			CTR[factored] = ctr;
		}
		ctr++;
	}
	else
		ctr++;
	for (i = 1; i <= FBASE; i++)
	{
		t = x % (long)(base[i]);
		R = (x < 0 && (USL)t) ? (USL)(t + base[i]) : (USL)t;
		s = 0;
		if (R == r1[i] || R == r2[i])
		{
			while (MODINT0_(Tmp, base[i], &T) == 0)
			{
				s++;
				FREEMPI(Tmp);
				Tmp = T;
			}
			FREEMPI(T);
		}
		sp2 = s & 1;
		word = (word << 1) + sp2;
		if (sp2)
		{
			INDEX[factored] = (int)index;
			CTR[factored] = ctr;
		}
		ctr++;
		if (ctr == O32)
		{
			m1ptr[index] = word;
			index++;
			ctr = 0;
			word = 0;
		}
	}
	word = (word << 1) + ((Q->S == 1) ? (USL)0 : (USL)1);
	if (Q->S == -1)
	{
		INDEX[factored] = (FBASE + 1) >> O5;
		CTR[factored] = (FBASE + 1) & O31;
	}
	ctr++;
	word = word << (O32 - ctr);
	m1ptr[index] = word;
	return (Tmp);
}

unsigned long **idim2(USI row, USI col)
/*
* Builds an array of row * col unsigned longs.
* From page 182, Advanced C - tips and techniques, by P Anderson and
* G. Anderson.
*/
{
	register unsigned int i;
	unsigned long **prow, *pdata;

	pdata = (unsigned long *)ccalloc(row * col, sizeof(unsigned long));
	prow = (unsigned long **)mmalloc((USL)(row * sizeof(unsigned long *)));
	for (i = 0; i < row; i++)
	{
		prow[i] = pdata;
		pdata += col;
	}
	return (prow);
}

void ifree2(USL **pa, USI row, USI col)
/*
* Frees the 2-dimensional array of unsigned ints pa.
* From page 183, Advanced C - tips and techniques, by P Anderson and
* G. Anderson.
*/
{
	ffree((char *)*pa, row * col * sizeof(unsigned long)); /* free the data */
	ffree((char *)pa, row * sizeof(unsigned long *)); /* free the row pointers */
}
/*void REDUCTION();*/
/* this test the function REDUCTION() */
/*
main()
{
unsigned int i, j, m, n;
int s;
register int l;
register int k;
unsigned int M1[rows][cols], M2[rows][cols];

printf(" enter the number of rows and columns (<= 20) :");
s=scanf("%u%u", &m, &n);
for (i = 0; i < m; i++)
{
printf("enter row %u: ", i);
for (j = 0; j < n; j++)
s=scanf("%u", &M1[i][j]);
}
printf("matrix M1 entered = \n");
for (i = 0; i < m; i++)
{
for (j = 0; j < n; j++)
printf("%u", M1[i][j]);
printf("\n");
}
for (k = 0; k <= m; k++)
for (l = 0; l <= m; l++)
M2[k][l] = 0;
for (k = 0; k <= m; k++)
M2[k][k] = 1;
printf("identity matrix = \n");
for (i= 0; i < m; i++)
{
for (j= 0; j < n; j++)
printf("%u", M2[i][j]);
printf("\n");
}
REDUCTION(M1, M2, m, n);
printf("reduced M1  = \n");
for (i= 0; i < m; i++)
{
for (j= 0; j < n; j++)
printf("%u", M1[i][j]);
printf("\n");
}
printf("identity matrix changed to = \n");
for (i= 0; i < m; i++)
{
for (j= 0; j < n; j++)
printf("%u", M2[i][j]);
printf("\n");
}
exit (0);
}

*/
void REDUCTION(USI **M1, USI **M2, USI m, USI n, int SCAN[])
/*
* A form of right-to-left gaussian reduction on M1[m,n] over Z_2, doing the
* same operations on M2[m,m], where M2 is the unit matrix of same size
* as M1.
* From page 188, "A method of factoring and the factorization of F_7",
* M.A. Morrison and J. Brillhart, Math. Comp. 1975, 183-205.
*/
{
	register int l;
	register int k;
	register int i;
	register int j;

	for (j = n - 1; j >= 0; j--)
	{
		for (i = 0; i < m; i++)
		{
			if (SCAN[i] == j)
			{
				for (k = i + 1; k < m; k++)
				{
					if (SCAN[k] == j)
					{
						for (l = 0; l <= j; l++)
							M1[k][l] ^= M1[i][l];
						for (l = 0; l < m; l++)
							M2[k][l] ^= M2[i][l];
						for (l = n - 1; l >= 0 && M1[k][l] == 0; l--)
							;
						SCAN[k] = l;
					}
				}
			}
		}
	}
}

void REDUCTION1(USL **M1, USL **M2, USI m, USI n, USI m1bitcols, USI m2bitcols, int INDEX[], USI CTR[])
/*
* A form of right-to-left gaussian reduction on M1[m,n] over Z_2, doing the
* same operations on M2[m,m], where M2 is the unit matrix of same size
* as M1.
* From page 188, "A method of factoring and the factorization of F_7",
* Math. Comp. 1975, 183-205. Used by Peter Adams in MPQS1().
*/
{
	unsigned long *p3, *p4;
	unsigned long *p1, *p2;
	register int j, l;
	register int i, k;
	int jword, jbit;

	for (j = n - 1; j >= 0; j--)
	{
		jword = j >> O5;
		jbit = j & O31;
		for (i = 0; i < m; i++)
		{
			p3 = M1[i];
			p4 = M2[i];
			if (INDEX[i] == jword && CTR[i] == jbit)
			{
				for (k = i + 1; k < m; k++)
				{
					p1 = M1[k];
					p2 = M2[k];
					if (INDEX[k] == jword && CTR[k] == jbit)
					{
						for (l = 0; l < m1bitcols; l++)
							p1[l] ^= p3[l];
						for (l = 0; l < m2bitcols; l++)
							p2[l] ^= p4[l];
						for (l = m1bitcols - 1; l >= 0 && p1[l] == 0; l--)
							;
						if ((INDEX[k] = l) >= 0)
							CTR[k] = O31 - BINARY(p1[l]);
					}
				}
			}
		}
	}
}

void EXP_UPDATE(MPI *Q, USL base[], USI FBASE, USI exponents[])
/*
* updates the exponent vector exponent[], where Q is in S.
*/
{
	unsigned int i, s;
	MPI *Tmp, *T;

	Tmp = ABSI(Q);
	for (i = 0; i <= FBASE; i++)
	{
		s = 0;
		while (MODINT0_(Tmp, base[i], &T) == 0)
		{
			s++;
			FREEMPI(Tmp);
			Tmp = T;
		}
		FREEMPI(T);
		exponents[i] += s;
		/*
		if (s)
		printf(".%lu^%u", base[i], s);
		*/
	}
	FREEMPI(Tmp);
	return;
}

void print_array(USL **M, int rows, int cols)
/*
Prints the contents of the given array.
*/
{
	unsigned long i, j;
	int elmt;

	for (i = 0; i<rows; i++)
	{
		for (j = 0; j<cols; j++)
		{
			elmt = get_element(M, i, j);
			printf("%d", elmt);
		}
		printf("\n");
	}
	printf("\n");
}

