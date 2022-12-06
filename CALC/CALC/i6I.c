/* i6I.c */
/* matrix operations using MPI's */
#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_NONSTDC_NO_DEPRECATE
#include <stdio.h>
#ifdef _WIN32
#include "unistd_DOS.h"
#else
#include <unistd.h>
#endif
#include <string.h>
#include <stdlib.h>
#include "integer.h"
#include "fun.h"

extern unsigned int HAVASFLAG;
extern unsigned int FPRINTMATIFLAG;
extern unsigned int MLLLVERBOSE;

MPMATI *BUILDMATI(unsigned int m, unsigned int n)
/*
* Allocates space for an m x n matrix of MPI's.
*/
{
	MPMATI *N;
	unsigned int i;

	N = (MPMATI *)mmalloc(sizeof(MPMATI));
	N->R = m;
	N->C = n;
	N->V = (MPI ***)mmalloc(m * sizeof(MPI **));
	for (i = 0; i < m; i++)
		N->V[i] = (MPI **)mmalloc(n * sizeof(MPI *));
	return (N);
}

MPMATI *COPYMATI(MPMATI *Mptr)
/*
* a replacement for the assignment *Nptr = *Mptr.
*/
{
	unsigned int i, j, m, n;
	MPMATI *Nptr;

	m = Mptr->R;
	n = Mptr->C;
	Nptr = BUILDMATI(Mptr->R, Mptr->C);
	for (i = 0; i < m; i++)
		for (j = 0; j < n; j++)
			Nptr->V[i][j] = COPYI(Mptr->V[i][j]);
	return (Nptr);
}

void FREEMATI(MPMATI *Mptr)
/*
* frees the storage allocated to the two-dimensional array Mptr->V.
*/
{
	unsigned int i, j, m, n;

	if (Mptr == NULL)
		return;
	m = Mptr->R;
	n = Mptr->C;
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
			FREEMPI(Mptr->V[i][j]);
		ffree((char *)(Mptr->V[i]), n * sizeof(MPI *));
	}
	ffree((char *)(Mptr->V), m * sizeof(MPI **));
	ffree((char *)Mptr, sizeof(MPMATI));
	return;
}

MPMATI *FINPUTMATI(FILE *infile)
/*
* Inputs a matrix of MPI's from infile.
*/
{
	MPMATI *Temp;
	unsigned int u, i, j, m, n;
	char *line, *p;

	if (infile == stdin)
	{
		if (HAVASFLAG)
			printf("Enter the number of integers\n");
		else
			printf("Enter the number of rows and number of columns\n");
	}

	if (HAVASFLAG)
	{
		if (fscanf(infile, "%u", &m) != 1 || m == 0)
		{
			fprintf(stderr, "corrupt input: returning the 1 x 1 zero matrix\n");
			if (infile != stdin)
				exit(1);
			return ZEROMNI(1, 1);
		}
		/*
		if (infile != stdin)
		{
		if (fscanf(infile, "%u", &n) != 1 || n == 0)
		{
		fprintf(stderr, "corrupt input: returning the 1 x 1 zero matrix\n");
		if (infile != stdin)
		exit (1);
		return ZEROMNI(1, 1);
		}
		}
		if (infile == stdin)
		*/
		n = 1;
	}
	else
	{
		if (fscanf(infile, "%u%u", &m, &n) != 2 || m == 0 || n == 0)
		{
			fprintf(stderr, "corrupt input: returning the 1 x 1 zero matrix\n");
			if (infile != stdin)
				exit(1);
			return ZEROMNI(1, 1);
		}
	}
	fgetc(infile);
	Temp = ZEROMNI(m, n);
	for (i = 0; i < m; i++)
	{
		if (infile == stdin)
			printf("Enter row %u:\n", i + 1);
		p = line = Fgets(infile);
		for (j = 0; j < n; j++)
		{
			FREEMPI(Temp->V[i][j]);
			Temp->V[i][j] = INPUTSI(&p, &u);
			if (u == 0)
			{
				free(line);
				fprintf(stderr, "matrix input interrupted prematurely\n");
				if (infile != stdin)
					exit(1);
				GetReturn();
				return Temp;
			}
		}
		free(line);
	}
	return (Temp);
}

MPMATI *INPUTMATI()
/*
* Input a matrix of MPI's from stdin.
*/
{
	return FINPUTMATI(stdin);
}

void FPRINTMATI(FILE *outfile, USI i1, USI i2, USI j1, USI j2, MPMATI *Mptr)
/*
* Modified 8/1/93 using Peter Adams' improvement from, August 1992.
*/
{
	char **str, *buff;
	unsigned int tmp, i, j, nrow, ncol, len, nstr = 0;
	int *colwidth;
	unsigned int ct;

	nrow = i2 - i1 + 1;
	ncol = j2 - j1 + 1;
	str = (char **)mmalloc(nrow * ncol * sizeof(char *));
	colwidth = (int *)mmalloc(ncol * sizeof(int));

	for (i = 0; i<ncol; i++)
		colwidth[i] = 0;
	for (i = i1; i <= i2; i++)
	{
		for (j = j1; j <= j2; j++)
		{
			tmp = 1 + LENGTHI(Mptr->V[i][j]);
			buff = (char *)mmalloc(tmp * sizeof(char));
			SPRINTI(buff, Mptr->V[i][j]);
			if ((len = strlen(buff)) > colwidth[j - j1])
				colwidth[j - j1] = len;
			str[nstr++] = strcpy((char *)mmalloc(len + 1), buff);
			ffree(buff, tmp * sizeof(char));
		}
	}
	if (!FPRINTMATIFLAG)
	{
		if (outfile != stdout)
			fprintf(outfile, "%u %u\n", nrow, ncol);
	}
	ct = 0;
	for (i = i1; i <= i2; i++)
	{
		for (j = j1; j <= j2; j++)
		{
			fprintf(outfile, "%*s ", colwidth[j - j1], str[ct]);
			ffree(str[ct], (unsigned int)(strlen(str[ct]) + 1) * sizeof(char));
			if ((ct % ncol) == (ncol - 1))
				fprintf(outfile, "\n");
			ct++;
		}
	}
	ffree((char *)str, nrow * ncol * sizeof(char *));
	ffree((char *)colwidth, ncol * sizeof(int));
	return;
}

void PRINTMATI(USI i1, USI i2, USI j1, USI j2, MPMATI *Mptr)
/*
* prints *Mptr from rows i1 to i2 and cols j1 to j2.
*/
{
	unsigned int i, j;

	if (!Mptr) return;
	i = i2 - i1 + 1;
	j = j2 - j1 + 1;
	if (MAX((int)i, (int)j) > 50)
	{
		printf("(The number of rows or columns to be printed exceeds 50;\n");
		printf("there is no point in printing this matrix on the screen.)\n");
		return;
	}
	FPRINTMATI(stdout, i1, i2, j1, j2, Mptr);
	return;
}

MPMATI *FINPUTMATFILEI_I()
{
	char *p;
	MPMATI *Mptr;
	FILE *infile;

	printf("enter the name of your matrix file:");

	p = Gets();
	printf("p : %s\n", p);
	infile = fopen(p, "r");
	if (infile == NULL)
		execerror("does not exist or cannot be opened", "p");
	Mptr = FINPUTMATI(infile);
	fclose(infile);
	printf("matrix from file %s successfully inputted\n", p);
	GetReturn();
	return (Mptr);
}


MPMATI *ZEROMNI(USI m, USI n)
/*
* returns the zero  m x n matrix.
*/
{
	unsigned int i, j;
	MPMATI *Mptr;

	Mptr = BUILDMATI(m, n);
	for (i = 0; i <= m - 1; i++)
	{
		for (j = 0; j <= n - 1; j++)
			Mptr->V[i][j] = ZEROI();
	}
	return (Mptr);
}

MPMATI *SWAP_ROWSI(USI p, USI q, MPMATI *Mptr)
/*
* interchange rows p and q (C notation) of *Mptr.
*/
{
	MPI **temp;
	MPMATI *Nptr;

	Nptr = COPYMATI(Mptr);
	temp = Nptr->V[p];
	Nptr->V[p] = Nptr->V[q];
	Nptr->V[q] = temp;
	return (Nptr);
}

MPMATI *SWAP_COLSI(USI p, USI q, MPMATI *Mptr)
/*
* interchange cols p and q (C notation) of *Mptr.
*/
{
	unsigned int i;
	MPI *Temp;
	MPMATI *Nptr;

	Nptr = COPYMATI(Mptr);
	for (i = 0; i <= Nptr->R - 1; i++)
	{
		Temp = Nptr->V[i][p];
		Nptr->V[i][p] = Nptr->V[i][q];
		Nptr->V[i][q] = Temp;
	}
	return (Nptr);
}

MPMATI *SWAP_ROWSI1(USI p, USI q, MPMATI *Mptr)
/*
* Updates Mptr by interchange rows p and q (C notation) of *Mptr.
*/
{
	MPI **temp;

	temp = Mptr->V[p];
	Mptr->V[p] = Mptr->V[q];
	Mptr->V[q] = temp;
	return (Mptr);
}

MPMATI *SWAP_COLSI1(USI p, USI q, MPMATI *Mptr)
/*
* Updates Mptr by interchanging cols p and q (C notation) of *Mptr.
*/
{
	unsigned int i;
	MPI *Temp;

	for (i = 0; i <= Mptr->R - 1; i++)
	{
		Temp = Mptr->V[i][p];
		Mptr->V[i][p] = Mptr->V[i][q];
		Mptr->V[i][q] = Temp;
	}
	return (Mptr);
}

unsigned int *KB_ROW(MPMATI *Aptr, USI *nz)
/*
* The Kannan-Bachem reduction to column permuted Hermite normal form, for
* use in their Smith Normal form algorithm.
*/
{
	unsigned int *T, i, j, k, s, m, n, tmp, is_zero;
	MPI *D, *P, *Q, *R, *S, *Tmp0, *Tmp1, *Tmp2, *Tmp3;

	m = Aptr->R;
	n = Aptr->C;
	T = (unsigned int *)mmalloc(n * sizeof(unsigned int));
	for (i = 0; i < n; i++)
		T[i] = i;
	s = m - 1;
	/* s = index of last row of A not known to be zero. */
	i = 0;
	while (i <= s)
	{
		/*	printf("i = %u\n", i);*/
		for (k = 0; k < i; k++)
		{
			D = EUCLIDI(Aptr->V[k][T[k]], Aptr->V[i][T[k]], &P, &Q);
			Tmp0 = INT(Aptr->V[i][T[k]], D);
			R = MINUSI(Tmp0);
			FREEMPI(Tmp0);
			S = INT(Aptr->V[k][T[k]], D);
			FREEMPI(D);
			for (j = 0; j < n; j++)
			{
				Tmp0 = Aptr->V[k][j];
				Tmp1 = MULTI(P, Aptr->V[k][j]);
				Tmp2 = MULTI(Q, Aptr->V[i][j]);
				Aptr->V[k][j] = ADDI(Tmp1, Tmp2);
				FREEMPI(Tmp1);
				FREEMPI(Tmp2);
				Tmp3 = Aptr->V[i][j];
				Tmp1 = MULTI(R, Tmp0);
				FREEMPI(Tmp0);
				Tmp2 = MULTI(S, Aptr->V[i][j]);
				Aptr->V[i][j] = ADDI(Tmp1, Tmp2);
				FREEMPI(Tmp3);
				FREEMPI(Tmp1);
				FREEMPI(Tmp2);
			}
			FREEMPI(P);
			FREEMPI(Q);
			FREEMPI(R);
			FREEMPI(S);
			RODI(Aptr, T, k);
		}
		is_zero = 1;
		for (j = 0; j < n; j++)
		{
			if ((Aptr->V[i][j])->S)
			{
				is_zero = 0;
				break;
			}
		}
		if (is_zero)
		{
			if (i < s)
				SWAP_ROWSI1(i, s, Aptr);
			s--;
		}
		else
		{
			for (j = 0; j < n; j++)
			{
				if ((Aptr->V[i][T[j]])->S)
					break;
			}
			if (j != i)
			{
				tmp = T[j];
				T[j] = T[i];
				T[i] = tmp;
			}
			RODI(Aptr, T, i);
			i++;
		}
	}
	*nz = s + 1;
	return (T);
}


void RODI(MPMATI *Aptr, USI *T, USI j)
/* Kannan_Bachem subroutine: input an integer j, 0 <= j < Aptr->R such that
* elt(Aptr, j, T[j]) != 0 and elt(Aptr, j, T[k]) = 0 for 0 <= k < j.
* entries elt(Aptr, i, T[j]), 0 <= i < j, are reduced mod elt(Aptr, j, T[j]).
*/
{
	unsigned int i, k, n;
	MPI *Q;

	n = Aptr->C;
	if ((Aptr->V[j][T[j]])->S == -1)
	{
		for (k = 0; k < n; k++)
		{
			if ((Aptr->V[j][k])->S != 0)
				(Aptr->V[j][k])->S = -((Aptr->V[j][k])->S);
		}
	}
	for (i = 0; i < j; i++)
	{
		Q = INTI(Aptr->V[i][T[j]], Aptr->V[j][T[j]]);
		Aptr = ROWSUBI(j, i, Q, Aptr);
		FREEMPI(Q);
	}
}

void SORTI(USI *T, USI nz)
/*
* sorts the array T[0],...,T[nz-1] into ascending order, thereby returning
* a permutation of 0,...,nz-1 used in PERMUTE_MATI below.
*/
{
	unsigned int i, j, temp, *P;

	P = (unsigned int *)mmalloc(nz * sizeof(unsigned int));
	for (i = 0; i < nz; i++)
		P[i] = i;
	for (i = 0; i < nz - 1; i++)
	{
		for (j = i + 1; j < nz; j++)
			if (T[i] > T[j])
			{
				temp = T[i];
				T[i] = T[j];
				T[j] = temp;
				temp = P[i];
				P[i] = P[j];
				P[j] = temp;
			}
	}
	for (i = 0; i < nz; i++)
		T[i] = P[i];
	ffree((char *)P, nz * sizeof(unsigned int));
	return;
}

MPMATI *HERMITE1(MPMATI *Aptr, USI *T, USI nz)
/*
* Takes the output of KB_ROW() and permutes the rows to get the Hermite normal
* form from Aptr.
*/
{
	MPMATI *Tmp;
	unsigned int i, j;

	SORTI(T, nz);
	Tmp = COPYMATI(Aptr);
	for (i = 0; i < nz; i++)
	{
		for (j = 0; j < Aptr->C; j++)
		{
			FREEMPI(Tmp->V[i][j]);
			Tmp->V[i][j] = COPYI(Aptr->V[T[i]][j]);
		}
	}
	return (Tmp);
}

MPMATI *ROWSUBI(USI p, USI q, MPI *Aptr, MPMATI *Mptr)
/*
* Updates Mptr: subtracts *Aptr times the p-th row of *Mptr from the q-th.
* 0 <= p, q <= Mprt->R - 1.
*/
{
	MPI *Y, *Temp;
	unsigned int l;
	if (Aptr->S != 0)
	{
		for (l = 0; l <= Mptr->C - 1; l++)
		{

			if ((Mptr->V[p][l])->S != 0)
			{
				Y = MULTI(Aptr, Mptr->V[p][l]);
				Temp = Mptr->V[q][l];
				Mptr->V[q][l] = SUBI(Mptr->V[q][l], Y);
				FREEMPI(Y);
				FREEMPI(Temp);
			}
		}
	}
	return (Mptr);
}


MPMATI *COLSUBI(USI p, USI q, MPI *Aptr, MPMATI *Mptr)
/*
* subtracts *Aptr times the p-th column of *Mptr from the q-th.
* 0 <= p, q <= Mprt->C - 1.
*/
{
	MPI *Y, *Temp;
	unsigned int k;

	if (Aptr->S != 0)
	{
		for (k = 0; k <= Mptr->R - 1; k++)
		{
			if ((Mptr->V[k][p])->S != 0)
			{
				Y = MULTI(Aptr, Mptr->V[k][p]);
				Temp = Mptr->V[k][q];
				Mptr->V[k][q] = SUBI(Mptr->V[k][q], Y);
				FREEMPI(Y);
				FREEMPI(Temp);
			}
		}
	}
	return (Mptr);
}

MPMATI *IDENTITYI(USI n)
/*
* returns the identity matrix of size n.
*/
{
	unsigned int i, j;
	MPMATI *Mptr;

	Mptr = BUILDMATI(n, n);
	for (i = 0; i <= n - 1; i++)
	{
		for (j = 0; j <= n - 1; j++)
		{
			if (i == j)
				Mptr->V[i][j] = ONEI();
			else
				Mptr->V[i][j] = ZEROI();
		}
	}
	return (Mptr);
}

MPMATI *MULTMATI(MPMATI *Mptr, MPMATI *Nptr)
/*
* Here Mptr->C = Nptr->R.
* returns (*Mptr) * (*Nptr).
*/
{
	MPI *X, *Y, *Temp;
	unsigned int i, j, k;
	MPMATI *Lptr;

	Lptr = BUILDMATI(Mptr->R, Nptr->C);
	for (i = 0; i <= Mptr->R - 1; i++)
	{
		for (j = 0; j <= Nptr->C - 1; j++)
		{
			X = ZEROI();
			for (k = 0; k <= Mptr->C - 1; k++)
			{
				Y = MULTI(Mptr->V[i][k], Nptr->V[k][j]);
				Temp = X;
				X = ADDI(X, Y);
				FREEMPI(Temp);
				FREEMPI(Y);
			}
			Lptr->V[i][j] = X;
		}
	}
	return (Lptr);
}

MPMATI *ADD_MULT_ROWI(USI p, USI q, MPI *Aptr, MPMATI *Mptr)
/*
* adding *Aptr times row p to row q of *Mptr.
*/
{
	unsigned int j;
	MPI *X, *Temp;
	MPMATI *Nptr;

	Nptr = COPYMATI(Mptr);
	for (j = 0; j <= Nptr->C - 1; j++)
	{
		X = MULTI(Nptr->V[p][j], Aptr);
		Temp = Nptr->V[q][j];
		Nptr->V[q][j] = ADDI(Temp, X);
		FREEMPI(X);
		FREEMPI(Temp);
	}
	return (Nptr);
}

void *ADD_MULT_ROWI0(USI p, USI q, MPI *Aptr, MPMATI *Mptr)
/*
* adding *Aptr times row p to row q of *Mptr and overwrites *Mptr.
*/
{
	unsigned int j;
	MPI *X, *Temp;

	for (j = 0; j <= Mptr->C - 1; j++)
	{
		X = MULTI(Mptr->V[p][j], Aptr);
		Temp = Mptr->V[q][j];
		Mptr->V[q][j] = ADDI(Temp, X);
		FREEMPI(X);
		FREEMPI(Temp);
	}
	return (Mptr);
}

MPMATI *TRANSPOSEI(MPMATI *Mptr)
/*
* returns the transpose of *Mptr.
*/
{
	MPMATI *Nptr;
	unsigned int i, j;

	Nptr = BUILDMATI(Mptr->C, Mptr->R);
	for (j = 0; j < Nptr->R; j++)
		for (i = 0; i < Nptr->C; i++)
			Nptr->V[j][i] = COPYI(Mptr->V[i][j]);
	return (Nptr);
}

MPMATI *DELETE_ROWI(USI r, MPMATI *Mptr)
/*
* deletes row r of *Mptr.
*/
{
	unsigned int i, j;
	MPMATI *Nptr;

	r--;
	Nptr = BUILDMATI(Mptr->R - 1, Mptr->C);
	for (i = 0; i <= Nptr->R - 1; i++)
	{
		if (i < r)
		{
			for (j = 0; j <= Nptr->C - 1; j++)
				Nptr->V[i][j] = COPYI(Mptr->V[i][j]);
		}
		else
		{
			for (j = 0; j <= Nptr->C - 1; j++)
				Nptr->V[i][j] = COPYI(Mptr->V[i + 1][j]);
		}
	}
	return (Nptr);
}

unsigned int *KB_ROWP(MPMATI *Aptr, MPMATI **Pptr, USI *nz)
/*
* The Kannan-Bachem reduction to column permuted Hermite normal form, for
* use in their Smith Normal form algorithm. The transforming unimodular matrix
* **Pptr is returned.
*/
{
	unsigned int *T, i, j, k, s, m, n, tmp, is_zero;
	MPI *D, *P, *Q, *R, *S, *Tmp0, *Tmp1, *Tmp2, *Tmp3;
	MPI *Tmp0P, *Tmp1P, *Tmp2P, *Tmp3P;

	m = Aptr->R;
	n = Aptr->C;
	*Pptr = IDENTITYI(m);
	T = (unsigned int *)mmalloc(n * sizeof(unsigned int));
	for (i = 0; i < n; i++)
		T[i] = i;
	s = m - 1;
	/* s = index of last row of A not known to be zero. */
	i = 0;
	while (i <= s)
	{
		if (MLLLVERBOSE)
			printf("i = %u\n", i);
		for (k = 0; k < i; k++)
		{
			D = EUCLIDI(Aptr->V[k][T[k]], Aptr->V[i][T[k]], &P, &Q);
			Tmp0 = INT(Aptr->V[i][T[k]], D);
			R = MINUSI(Tmp0);
			FREEMPI(Tmp0);
			S = INT(Aptr->V[k][T[k]], D);
			FREEMPI(D);
			for (j = 0; j < n; j++)
			{
				Tmp0 = Aptr->V[k][j];
				Tmp1 = MULTI(P, Aptr->V[k][j]);
				Tmp2 = MULTI(Q, Aptr->V[i][j]);
				Aptr->V[k][j] = ADDI(Tmp1, Tmp2);
				FREEMPI(Tmp1);
				FREEMPI(Tmp2);
				Tmp3 = Aptr->V[i][j];
				Tmp1 = MULTI(R, Tmp0);
				FREEMPI(Tmp0);
				Tmp2 = MULTI(S, Aptr->V[i][j]);
				Aptr->V[i][j] = ADDI(Tmp1, Tmp2);
				FREEMPI(Tmp3);
				FREEMPI(Tmp1);
				FREEMPI(Tmp2);
			}
			for (j = 0; j < m; j++)
			{
				Tmp0P = (*Pptr)->V[k][j];
				Tmp1P = MULTI(P, (*Pptr)->V[k][j]);
				Tmp2P = MULTI(Q, (*Pptr)->V[i][j]);
				(*Pptr)->V[k][j] = ADDI(Tmp1P, Tmp2P);
				FREEMPI(Tmp1P);
				FREEMPI(Tmp2P);
				Tmp3P = (*Pptr)->V[i][j];
				Tmp1P = MULTI(R, Tmp0P);
				FREEMPI(Tmp0P);
				Tmp2P = MULTI(S, (*Pptr)->V[i][j]);
				(*Pptr)->V[i][j] = ADDI(Tmp1P, Tmp2P);
				FREEMPI(Tmp3P);
				FREEMPI(Tmp1P);
				FREEMPI(Tmp2P);
			}
			FREEMPI(P);
			FREEMPI(Q);
			FREEMPI(R);
			FREEMPI(S);
			RODIP(Aptr, *Pptr, T, k);
		}
		is_zero = 1;
		for (j = 0; j < n; j++)
		{
			if ((Aptr->V[i][j])->S)
			{
				is_zero = 0;
				break;
			}
		}
		if (is_zero)
		{
			if (i < s)
			{
				SWAP_ROWSI1(i, s, Aptr);
				SWAP_ROWSI1(i, s, *Pptr);
			}
			s--;
		}
		else
		{
			for (j = 0; j < n; j++)
			{
				if ((Aptr->V[i][T[j]])->S)
					break;
			}
			if (j != i)
			{
				tmp = T[j];
				T[j] = T[i];
				T[i] = tmp;
			}
			RODIP(Aptr, *Pptr, T, i);
			i++;
		}
	}
	*nz = s + 1;
	return (T);
}

void RODIP(MPMATI *Aptr, MPMATI *Pptr, USI *T, USI j)
/* Kannan_Bachem subroutine: input an integer j, 0 <= j < Aptr->R such that
* elt(Aptr, j, T[j]) != 0 and elt(Aptr, j, T[k]) = 0 for 0 <= k < j.
* entries elt(Aptr, i, T[j]), 0 <= i < j, are reduced mod elt(Aptr, j, T[j]).
* *Pptr is updated.
*/
{
	unsigned int i, k, m, n;
	MPI *Q;

	m = Aptr->R;
	n = Aptr->C;
	if ((Aptr->V[j][T[j]])->S == -1)
	{
		for (k = 0; k < n; k++)
		{
			if ((Aptr->V[j][k])->S != 0)
				(Aptr->V[j][k])->S = -((Aptr->V[j][k])->S);
		}
		for (k = 0; k < m; k++)
		{
			if ((Pptr->V[j][k])->S != 0)
				(Pptr->V[j][k])->S = -((Pptr->V[j][k])->S);
		}
	}
	for (i = 0; i < j; i++)
	{
		Q = INTI(Aptr->V[i][T[j]], Aptr->V[j][T[j]]);
		Aptr = ROWSUBI(j, i, Q, Aptr);
		Pptr = ROWSUBI(j, i, Q, Pptr);
		FREEMPI(Q);
	}
}

MPMATI *HERMITE1P(MPMATI *Aptr, MPMATI *Pptr, MPMATI **Qptr, USI *T, USI nz)
/*
* Takes the output of KB_ROWP() and permutes the rows to get the Hermite normal
* form from Aptr. *Pptr is updated to **Qptr.
*/
{
	MPMATI *Tmp;
	unsigned int i, j;

	SORTI(T, nz);
	Tmp = COPYMATI(Aptr);
	*Qptr = COPYMATI(Pptr);
	for (i = 0; i < nz; i++)
	{
		for (j = 0; j < Aptr->C; j++)
		{
			FREEMPI(Tmp->V[i][j]);
			Tmp->V[i][j] = COPYI(Aptr->V[T[i]][j]);
		}
		for (j = 0; j < Aptr->R; j++)
		{
			FREEMPI((*Qptr)->V[i][j]);
			(*Qptr)->V[i][j] = COPYI(Pptr->V[T[i]][j]);
		}
	}
	return (Tmp);
}

MPMATI *ADDMATI(MPMATI *Mptr, MPMATI *Nptr)
/*
* returns *Mptr + *Nptr.
*/
{
	unsigned int i, j;
	MPMATI *Lptr;

	Lptr = BUILDMATI(Mptr->R, Mptr->C);
	for (i = 0; i <= Lptr->R - 1; i++)
		for (j = 0; j <= Lptr->C - 1; j++)
			Lptr->V[i][j] = ADDI(Mptr->V[i][j], Nptr->V[i][j]);
	return (Lptr);
}

MPMATI *SUBMATI(MPMATI *Mptr, MPMATI *Nptr)
/*
* returns *Mptr - *Nptr.
*/
{
	unsigned int i, j;
	MPMATI *Lptr;

	Lptr = BUILDMATI(Mptr->R, Mptr->C);
	for (i = 0; i <= Lptr->R - 1; i++)
		for (j = 0; j <= Lptr->C - 1; j++)
			Lptr->V[i][j] = SUBI(Mptr->V[i][j], Nptr->V[i][j]);
	return (Lptr);
}

