/* i9.c */
#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_NONSTDC_NO_DEPRECATE
#include <stdio.h>
#include <stdlib.h>
#include "integer.h"
#include "fun.h"

extern long int nettbytes;
/*
MPI *MAXI, *PMAXI, *QMAXI;
*/
extern int answer;

MPI **SMITHI(MPMATI *Mptr, MPMATI **Pptr, MPMATI **Qptr, USI *ptr, MPI *Eptr, USI m1, USI n1)
/*
* returns the invariant factors of *Mptr.
* **Pptr and **Qptr are invertible matrices such that
* **Pptr x *Mptr x **Qptr is the Smith canonical form *Nptr below.
* *ptr is the number of invariant factors.
* See Rings, Modules and Linear Algebra by B. Hartley and T. O. Hawkes,
* Chapman and Hall, 1970.
* *Eptr is the cutoff above which we bring in MLLL. In any case MLLL is
* done before each block is processed, as a chance of getting small entries.
* This modification of the LLL algorithm was outlined by M. Pohst in J.
* Symbolic Computation (1987) 4, 123-127. m1/n1 is the LLL parameter alpha,
* 1/4 < alpha |= 1.
* I initially coded this I think in 1992, along the lines of suggestions by
* George Havas.
*/
{
	unsigned int i, k, l, m, n, p, q, j, t, z, size, r, s;
	unsigned int flag = 1, FLAG = 1;
	MPI *Q, **N, *Tmp;
	MPMATI *Nptr, *Temp;
	int isok;

	m = Mptr->R;
	n = Mptr->C;
	Nptr = COPYMATI(Mptr);
	*Pptr = IDENTITYI(m);
	*Qptr = IDENTITYI(n);
	/*
	MAXI = MAXELTI(Mptr);
	PMAXI = MAXELTI(*Pptr);
	QMAXI = MAXELTI(*Qptr);
	*/
	z = MIN((int)m - 1, (int)n - 1);
	N = (MPI **)mmalloc((z + 1) * sizeof(MPI *));
	for (i = 0; i <= z; i++)
	{
		ZEROTESTI(0, Nptr, &r, &s);
		if (r == Nptr->R && s == Nptr->C)
			break;
		flag = HAVAS_PIVOTI(Nptr, &p, &q, r, s);
		if (p > 0)/* must have a non-zero row 1 for MLLL */
		{
			Nptr = SWAP_ROWSI1(0, p, Nptr);
			*Pptr = SWAP_ROWSI1(i, p + i, *Pptr);
		}
		if (q > 0)
		{
			Nptr = SWAP_COLSI1(0, q, Nptr);
			*Qptr = SWAP_COLSI1(i, q + i, *Qptr);
		}
		if (flag)
			FLAG = 1;
		if (!flag && FLAG)
		{
			size = Nptr->R;
			Temp = Nptr;
			Nptr = BASIS_REDUCTION(Nptr, Pptr, i, m1, n1);
			/*Nptr = BASIS_REDUCTION0(Nptr, 3, 4);*/
			FREEMATI(Temp);
			ZEROTESTI(0, Nptr, &r, &s);
			if (r == Nptr->R && s == Nptr->C)
				break;
			flag = HAVAS_PIVOTI(Nptr, &p, &q, r, s);
			if (p > 0)
			{
				Nptr = SWAP_ROWSI1(0, p, Nptr);
				*Pptr = SWAP_ROWSI1(i, p + i, *Pptr);
			}
			if (q > 0)
			{
				Nptr = SWAP_COLSI1(0, q, Nptr);
				*Qptr = SWAP_COLSI1(i, q + i, *Qptr);
			}
			FLAG = 0;
			size = size - Nptr->R;
			m = m - size;
			z = MIN((int)m - 1, (int)n - 1);
		}
		/* Now elt(Nptr, 0, 0) != 0. */
		while (1)
		{
			while (1)
			{ /* PA bit */
				isok = 1;
				for (j = 0; j < Nptr->C; j++)
				{
					if ((Nptr->V[0][j])->S == -1)
					{
						for (t = 0; t < Nptr->R; t++)
							(Nptr->V[t][j])->S = -(Nptr->V[t][j])->S;
						for (t = 0; t < n; t++)
							((*Qptr)->V[t][j + i])->S = -((*Qptr)->V[t][j + i])->S;
					}
				}
				/* Now all non-zero elements in row 0 are positive. */
				while (1)
				{
					l = ROWSEEKI(0, Nptr);
					if (l == Nptr->C)
						break;
					else
					{
						/* Case 1, p. 112. */
						Q = INT(Nptr->V[0][l], Nptr->V[0][0]);
						Nptr = COLSUBI(0, l, Q, Nptr);
						/*
						MAXI = UPDATEMAXI(MAXI, Nptr);
						*/
						*Qptr = COLSUBI(i, l + i, Q, *Qptr);
						/*
						QMAXI = UPDATEMAXI(QMAXI, *Qptr);
						*/
						FREEMPI(Q);
						Nptr = SWAP_COLSI1(0, l, Nptr);
						*Qptr = SWAP_COLSI1(i, l + i, *Qptr);
						/* This reduces the size of elt(Nptr, 0, 0). */
					}
				}
				for (j = 1; j < Nptr->R; j++)
				{
					if ((Nptr->V[j][0])->S == -1)
					{
						for (t = 0; t < Nptr->C; t++)
							(Nptr->V[j][t])->S = -(Nptr->V[j][t])->S;
						for (t = 0; t < (*Pptr)->R; t++)
							((*Pptr)->V[j + i][t])->S = -((*Pptr)->V[j + i][t])->S;
					}
				}
				/*
				Now all non-zero elements in column 0 are positive.
				*/
				while (1)
				{
					k = COLSEEKI(0, Nptr);
					if (k == Nptr->R)
						break;
					else
					{
						isok = 0;
						/* Case 2, p. 112. */
						Q = INT0(Nptr->V[k][0], Nptr->V[0][0]);
						Nptr = ROWSUBI(0, k, Q, Nptr);
						/*
						MAXI = UPDATEMAXI(MAXI, Nptr);
						*/
						*Pptr = ROWSUBI(i, k + i, Q, *Pptr);
						/*
						PMAXI = UPDATEMAXI(PMAXI, *Pptr);
						*/
						FREEMPI(Q);
						Nptr = SWAP_ROWSI1(0, k, Nptr);
						*Pptr = SWAP_ROWSI1(i, k + i, *Pptr);
						/* This reduces the size of elt(Nptr, 0, 0). */
					}
				}
				/* now l = Nptr->C and k = Nptr->R */
				if (isok)
					break;
			} /*PA bit*/
			  /*
			  Now (Case 3, p.112) elt(Nptr, 0, 0) divides each element to its right as well as below it.
			  */
			for (l = 1; l < Nptr->C; l++)
			{
				if ((Nptr->V[0][l])->S != 0)
				{
					/* NOTE: Swapping above may have made elt(Nptr, 0, l) <0 */
					Q = INT(Nptr->V[0][l], Nptr->V[0][0]);
					Nptr = COLSUBI(0, l, Q, Nptr);
					/*
					MAXI = UPDATEMAXI(MAXI, Nptr);
					*/
					*Qptr = COLSUBI(i, l + i, Q, *Qptr);
					/*
					QMAXI = UPDATEMAXI(QMAXI, *Qptr);
					*/
					FREEMPI(Q);
				}
			}
			for (k = 1; k < Nptr->R; k++)
			{
				if ((Nptr->V[k][0])->S != 0)
				{
					Q = INT0(Nptr->V[k][0], Nptr->V[0][0]);
					FREEMPI(Nptr->V[k][0]);
					Nptr->V[k][0] = ZEROI();
					*Pptr = ROWSUBI(i, k + i, Q, *Pptr);
					/*
					PMAXI = UPDATEMAXI(PMAXI, *Pptr);
					*/
					FREEMPI(Q);
				}
			}
			/*
			Now each element to the right of and below elt(Nptr, 0, 0) is zero.
			*/
			k = MATSEEKI(0, Nptr);
			if (k != Nptr->R)
			{
				for (t = 0; t < Nptr->C; t++)
				{
					if ((Nptr->V[k][t])->S != 0)
					{
						Tmp = Nptr->V[0][t];
						Nptr->V[0][t] = ADDI(Nptr->V[0][t], Nptr->V[k][t]);
						/*
						MAXI = UPDATEMAXI(MAXI, Nptr);
						*/
						FREEMPI(Tmp);
					}
				}
				for (t = 0; t < Mptr->R; t++)
				{
					if (((*Pptr)->V[k + i][t])->S != 0)
					{
						Tmp = (*Pptr)->V[i][t];
						(*Pptr)->V[i][t] = ADDI((*Pptr)->V[i][t], (*Pptr)->V[k + i][t]);
						/*
						PMAXI = UPDATEMAXI(PMAXI, *Pptr);
						*/
						FREEMPI(Tmp);
					}
				}
			}
			else
				break;
		}
		/*
		Now elt(Nptr, 0, 0) divides all elements elt(Nptr, k, l), k,l > 0.
		*/
		N[i] = COPYI(Nptr->V[0][0]);

		if (i < z)
		{
			Temp = Nptr;
			Nptr = DELETE_ROW1I(Nptr);
			FREEMATI(Temp);
			Temp = Nptr;
			Nptr = DELETE_COL1I(Nptr);
			FREEMATI(Temp);
		}
		size = Nptr->R;
		Tmp = MAXELTI(Nptr);
		printf("found invariant factor d[%u] = ", i + 1); PRINTI(N[i]); printf("\n");
		/*
		printf("MAXI = "); PRINTI(MAXI); printf("\n");
		printf("PMAXI = "); PRINTI(PMAXI); printf("\n");
		printf("QMAXI = "); PRINTI(QMAXI); printf("\n");
		*/
		printf("MAX_ENTRY(Nptr) = "); PRINTI(Tmp); printf("\n");
		l = RSV(Tmp, Eptr);
		FREEMPI(Tmp);
		if (l == 1 && !flag)/* we don't want to use MLLL if there are 1's or -1's still present */
		{
			ZEROTESTI(0, Nptr, &r, &s);
			if (r > 0)/* must have a non-zero row 1 for MLLL */
			{
				printf("r = %u > 0\n", r);
				Nptr = SWAP_ROWSI1(0, r, Nptr);
				*Pptr = SWAP_ROWSI1(i + 1, r + i + 1, *Pptr);
			}
			Temp = Nptr;
			Nptr = BASIS_REDUCTION(Nptr, Pptr, i + 1, m1, n1);
			/*Nptr = BASIS_REDUCTION0(Nptr, 3, 4);*/
			FREEMATI(Temp);
		}
		size = size - Nptr->R;
		m = m - size;
		z = MIN((int)m - 1, (int)n - 1);
	}
	FREEMATI(Nptr);
	*ptr = i++;
	/*
	printf("MAXI = "); PRINTI(MAXI); printf("\n");
	printf("PMAXI = "); PRINTI(PMAXI); printf("\n");
	printf("QMAXI = "); PRINTI(QMAXI); printf("\n");
	FREEMPI(MAXI);
	FREEMPI(PMAXI);
	FREEMPI(QMAXI);
	*/
	return (N);
}

MPI *ROWSUMI(MPMATI *Mptr, USI i)
{
	unsigned int j;
	MPI *S, *Temp1I, *Temp2I;

	S = ZEROI();
	for (j = 0; j < Mptr->C; j++)
	{
		Temp1I = ABSI(Mptr->V[i][j]);
		Temp2I = S;
		S = ADD0I(Temp1I, S);
		FREEMPI(Temp1I);
		FREEMPI(Temp2I);
	}
	return (S);
}

MPI *COLSUMI(MPMATI *Mptr, USI j)
{
	unsigned int i;
	MPI *S, *Temp1I, *Temp2I;

	S = ZEROI();
	for (i = 0; i < Mptr->R; i++)
	{
		Temp1I = ABSI(Mptr->V[i][j]);
		Temp2I = S;
		S = ADD0I(Temp1I, S);
		FREEMPI(Temp1I);
		FREEMPI(Temp2I);
	}
	return (S);
}

MPI *MINELTI(MPMATI *Mptr)
/* returns the minimum element of *Mptr */
{
	unsigned int i, j, I, J, flag = 0;

	I = J = 0;
	for (i = 0; i < Mptr->R; i++)
	{
		for (j = 0; j < Mptr->C; j++)
		{
			if ((Mptr->V[i][j])->S)
			{
				I = i;
				J = j;
				flag = 1;
				break;
			}
		}
		if (flag)
			break;
	}
	for (i = I; i < Mptr->R; i++)
	{
		for (j = 0; j < Mptr->C; j++)
		{
			if (RSV(Mptr->V[I][J], Mptr->V[i][j]) == 1 && (Mptr->V[i][j])->S)
			{
				I = i;
				J = j;
			}
		}
	}
	return (ABSI(Mptr->V[I][J]));
}

void MINELTI1(MPMATI *Mptr, USI *iptr, USI *jptr)
/* returns the position of the minimum element of *Mptr */
{
	unsigned int i, j, I, J, flag = 0;

	I = J = 0;
	for (i = 0; i < Mptr->R; i++)
	{
		for (j = 0; j < Mptr->C; j++)
		{
			if ((Mptr->V[i][j])->S)
			{
				I = i;
				J = j;
				flag = 1;
				break;
			}
		}
		if (flag)
			break;
	}
	for (i = 0; i < Mptr->R; i++)
	{
		for (j = 0; j < Mptr->C; j++)
		{
			if (RSV(Mptr->V[I][J], Mptr->V[i][j]) == 1 && (Mptr->V[i][j])->S)
			{
				I = i;
				J = j;
			}
		}
	}
	*iptr = I;
	*jptr = J;
	return;
}

MPI *MAXELTI(MPMATI *Mptr)
/* returns the maximum element of *Mptr */
{
	unsigned int i, j, I = 0, J = 0;

	for (i = 0; i < Mptr->R; i++)
	{
		for (j = 0; j < Mptr->C; j++)
		{
			if (RSV(Mptr->V[i][j], Mptr->V[I][J]) == 1)
			{
				I = i;
				J = j;
			}
		}
	}
	return (ABSI(Mptr->V[I][J]));
}

void MAXELTI1(MPMATI *Mptr, USI *iptr, USI *jptr)
/* returns the position of the maximum element of *Mptr */
{
	unsigned int i, j, I = 0, J = 0;

	for (i = 0; i < Mptr->R; i++)
	{
		for (j = 0; j < Mptr->C; j++)
		{
			if (RSV(Mptr->V[i][j], Mptr->V[I][J]) == 1)
			{
				I = i;
				J = j;
			}
		}
	}
	*iptr = I;
	*jptr = J;
	return;
}

MPI *UPDATEMAXI(MPI *S, MPMATI *Mptr)
{
	MPI *TempI;

	TempI = MAXELTI(Mptr);
	if (RSV(TempI, S) == 1)
	{
		FREEMPI(S);
		return(TempI);
	}
	else
	{
		FREEMPI(TempI);
		return (S);
	}
}

unsigned int HAVAS_PIVOTI(MPMATI *Mptr, USI *rptr, USI *cptr, USI r, USI s)
{
	MPI **RS, **CS, *S, *Temp1I, *Temp2I, *TempI;
	unsigned int i, j, I, J, m, n, flag;

	m = Mptr->R;
	n = Mptr->C;
	RS = (MPI **)mmalloc(m * sizeof(MPI *));
	CS = (MPI **)mmalloc(n * sizeof(MPI *));
	for (i = 0; i < m; i++)
	{
		TempI = ROWSUMI(Mptr, i);
		if (TempI->S == 0)
		{
			FREEMPI(TempI);
			RS[i] = MINUS_ONEI();
		}
		else
		{
			Temp1I = ONEI();
			RS[i] = SUB0I(TempI, Temp1I);
			FREEMPI(TempI);
			FREEMPI(Temp1I);
		}
	}
	for (j = 0; j < Mptr->C; j++)
	{
		TempI = COLSUMI(Mptr, j);
		if (TempI->S == 0)
		{
			FREEMPI(TempI);
			CS[j] = MINUS_ONEI();
		}
		else
		{
			Temp1I = ONEI();
			CS[j] = SUB0I(TempI, Temp1I);
			FREEMPI(TempI);
			FREEMPI(Temp1I);
		}
	}
	S = MULTI(RS[r], CS[s]); /* elt(r,s) is the 1st non-zero elt of Nptr */
	I = r; J = s;
	flag = 0;
	/* first minimize RS[i]*CS[j] over ones or minus-ones */
	for (i = r; i < m; i++)
	{
		for (j = 0; j < n; j++)
		{
			if ((Mptr->V[i][j])->D == 0 && (Mptr->V[i][j])->V[0] == 1)
			{
				flag = 1;
				Temp1I = S;
				Temp2I = MULTI(RS[i], CS[j]);
				if (RSV(Temp1I, Temp2I) == 1)
				{
					S = Temp2I;
					FREEMPI(Temp1I);
					I = i; J = j;
				}
				else
					FREEMPI(Temp2I);
			}
		}
	}
	if (!flag) /* minimize RS[i]*CS[j] over a matrix with no ones or minus-ones */
	{
		for (i = r; i < m; i++)
		{
			for (j = 0; j < n; j++)
			{
				if ((Mptr->V[i][j])->S)
				{
					Temp1I = S;
					Temp2I = MULTI(RS[i], CS[j]);
					if (RSV(Temp1I, Temp2I) == 1)
					{
						S = Temp2I;
						FREEMPI(Temp1I);
						I = i; J = j;
					}
					else
						FREEMPI(Temp2I);

				}
			}
		}
	}
	*rptr = I;
	*cptr = J;
	for (i = 0; i < m; i++)
		FREEMPI(RS[i]);
	for (j = 0; j < n; j++)
		FREEMPI(CS[j]);
	ffree((char *)RS, m * sizeof(MPI *));
	ffree((char *)CS, n * sizeof(MPI *));
	FREEMPI(S);
	return (flag);
}

void ZEROTESTI(USI i, MPMATI *Mptr, USI *p, USI *q)
/*
* returns the least indices *p >= i, *q >= i such that elt(Mptr, *p, *q) != 0.
* otherwise *p = Mptr->R and *q = Mptr->C.
*/
{
	unsigned int k, l;

	for (k = i; k <= Mptr->R - 1; k++)
	{
		for (l = i; l <= Mptr->C - 1; l++)
		{
			if ((Mptr->V[k][l])->S != 0)
			{
				*p = k;
				*q = l;
				return;
			}
		}
	}
	*p = Mptr->R;
	*q = Mptr->C;
	return;
}

unsigned int MATSEEKI(USI i, MPMATI *Mptr)
/*
* returns the least row index k, i < k <= Mptr->R - 1, such that
* elt(Mptr, k, l) is not divisible by elt(Mptr, i, i), for some l,
* i < l <= Mptr->C - 1; otherwise returns Mptr->R.
*/
{
	int s;
	unsigned int k, l;
	MPI *P, *Q;

	Q = ABSI(Mptr->V[i][i]);
	for (k = i + 1; k <= Mptr->R - 1; k++)
	{
		for (l = i + 1; l <= Mptr->C - 1; l++)
		{
			P = MOD(Mptr->V[k][l], Q);
			s = P->S;
			FREEMPI(P);
			if (s != 0)
			{
				FREEMPI(Q);
				return (k);
			}
		}
	}
	FREEMPI(Q);
	return (Mptr->R);
}

unsigned int ROWSEEKI(USI i, MPMATI *Mptr)
/*
* returns the least integer l > i such that elt(Mptr, i, l) is not divisible
* by elt(Mptr, i, i), otherwise returns Mptr->C.
*/
{
	int s;
	unsigned int l;
	MPI *P;

	for (l = i + 1; l <= Mptr->C - 1; l++)
	{
		if ((Mptr->V[i][l])->S != 0)
		{
			P = MOD0(Mptr->V[i][l], Mptr->V[i][i]);
			s = P->S;
			FREEMPI(P);
			if (s != 0)
				return (l);
		}
	}
	return (Mptr->C);
}

unsigned int COLSEEKI(USI i, MPMATI *Mptr)
/*
* returns the least integer k > i such that Mptr->V[k][i] is not divisible
* by Mptr->V[i][i], otherwise returns Mptr->R.
*/
{
	int s;
	unsigned int k;
	MPI *P;

	for (k = i + 1; k <= Mptr->R - 1; k++)
	{
		if ((Mptr->V[k][i])->S != 0)
		{
			P = MOD0(Mptr->V[k][i], Mptr->V[i][i]);
			s = P->S;
			FREEMPI(P);
			if (s != 0)
				return (k);
		}
	}
	return (Mptr->R);
}

MPMATI *ROWSUBI0(USI p, USI q, MPI *Aptr, MPMATI *Mptr)
/*
* subtracts *Aptr times the p-th row of *Mptr from the q-th.
* 0 <= p, q <= Mprt->R - 1, where elt(Mptr, i, j) = 0 if j < p.
*/
{
	MPI *Y, *Temp;
	unsigned int l;
	if (Aptr->S != 0)
	{
		for (l = p; l <= Mptr->C - 1; l++)
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

MPMATI *COLSUBI0(USI p, USI q, MPI *Aptr, MPMATI *Mptr)
/*
* subtracts *Aptr times the p-th column of *Mptr from the q-th.
* 0 <= p, q <= Mprt->C - 1, where elt(Mptr, i, j) = 0 if i < p.
*/
{
	MPI *Y, *Temp;
	unsigned int k;

	if (Aptr->S != 0)
	{
		for (k = p; k <= Mptr->R - 1; k++)
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

MPI **SMITHI1(MPMATI *Mptr, USI *ptr, MPI *Eptr, USI m1, USI n1)
/*
* returns the invariant factors of *Mptr.
* *ptr is the number of invariant factors.
* See Rings, Modules and Linear Algebra by B. Hartley and T. O. Hawkes,
* Chapman and Hall, 1970.
* *Eptr is the cutoff above which we bring in MLLL.
*/
{
	unsigned int i, k, l, m, n, p, q, j, t, z, size, r, s;
	unsigned int flag = 1;
	MPI *Q, **N, *Tmp;
	MPMATI *Nptr, *Temp;
	int isok;

	m = Mptr->R;
	n = Mptr->C;
	Nptr = COPYMATI(Mptr);
	/*
	MAXI = MAXELTI(Mptr);
	*/
	z = MIN((int)m - 1, (int)n - 1);
	N = (MPI **)mmalloc((z + 1) * sizeof(MPI *));
	for (i = 0; i <= z; i++)
	{
		ZEROTESTI(0, Nptr, &r, &s);
		if (r == Nptr->R && s == Nptr->C)
			break;
		flag = HAVAS_PIVOTI(Nptr, &p, &q, r, s);
		if (p > 0)/* must have a non-zero row 1 for MLLL */
			Nptr = SWAP_ROWSI1(0, p, Nptr);
		if (q > 0)
			Nptr = SWAP_COLSI1(0, q, Nptr);
		/*
		if (flag)
		FLAG = 1;
		if (!flag && FLAG)
		{
		size = Nptr->R;
		Temp = Nptr;
		printf("FLAG = 1 and flag = 0: doing MLLL:\n");
		Nptr = BASIS_REDUCTION0(Nptr, m1, n1);
		FREEMATI(Temp);
		ZEROTESTI(0, Nptr, &r, &s);
		if (r == Nptr->R && s == Nptr->C)
		break;
		flag = HAVAS_PIVOTI(Nptr, &p, &q, r, s);
		if (p > 0)
		Nptr = SWAP_ROWSI1(0, p, Nptr);
		if (q > 0)
		Nptr = SWAP_COLSI1(0, q, Nptr);
		FLAG = 0;
		size = size - Nptr->R;
		m = m - size;
		z = MIN((int) m - 1, (int) n - 1);
		}
		*/
		/* Now elt(Nptr, 0, 0) != 0. */
		while (1)
		{
			while (1)
			{ /* PA bit */
				isok = 1;
				for (j = 0; j < Nptr->C; j++)
				{
					if ((Nptr->V[0][j])->S == -1)
					{
						for (t = 0; t < Nptr->R; t++)
							(Nptr->V[t][j])->S = -(Nptr->V[t][j])->S;
					}
				}
				/* Now all non-zero elements in row 0 are positive. */
				while (1)
				{
					l = ROWSEEKI(0, Nptr);
					if (l == Nptr->C)
						break;
					else
					{
						/* Case 1, p. 112. */
						Q = INT(Nptr->V[0][l], Nptr->V[0][0]);
						Nptr = COLSUBI(0, l, Q, Nptr);
						/*
						MAXI = UPDATEMAXI(MAXI, Nptr);
						*/
						FREEMPI(Q);
						Nptr = SWAP_COLSI1(0, l, Nptr);
						/* This reduces the size of elt(Nptr, 0, 0). */
					}
				}
				for (j = 1; j < Nptr->R; j++)
				{
					if ((Nptr->V[j][0])->S == -1)
					{
						for (t = 0; t < Nptr->C; t++)
							(Nptr->V[j][t])->S = -(Nptr->V[j][t])->S;
					}
				}
				/*
				Now all non-zero elements in column 0 are positive.
				*/
				while (1)
				{
					k = COLSEEKI(0, Nptr);
					if (k == Nptr->R)
						break;
					else
					{
						isok = 0;
						/* Case 2, p. 112. */
						Q = INT0(Nptr->V[k][0], Nptr->V[0][0]);
						Nptr = ROWSUBI(0, k, Q, Nptr);
						/*
						MAXI = UPDATEMAXI(MAXI, Nptr);
						*/
						FREEMPI(Q);
						Nptr = SWAP_ROWSI1(0, k, Nptr);
						/* This reduces the size of elt(Nptr, 0, 0). */
					}
				}
				/* now l = Nptr->C and k = Nptr->R */
				if (isok)
					break;
			} /*PA bit*/
			  /*
			  Now (Case 3, p.112) elt(Nptr, 0, 0) divides each element to its right as well as below it.
			  */
			for (l = 1; l < Nptr->C; l++)
			{
				if ((Nptr->V[0][l])->S != 0)
				{
					/* NOTE: Swapping above may have made elt(Nptr, 0, l) <0 */
					Q = INT(Nptr->V[0][l], Nptr->V[0][0]);
					Nptr = COLSUBI(0, l, Q, Nptr);
					/*
					MAXI = UPDATEMAXI(MAXI, Nptr);
					*/
					FREEMPI(Q);
				}
			}
			for (k = 1; k < Nptr->R; k++)
			{
				if ((Nptr->V[k][0])->S != 0)
				{
					Q = INT0(Nptr->V[k][0], Nptr->V[0][0]);
					FREEMPI(Nptr->V[k][0]);
					Nptr->V[k][0] = ZEROI();
					FREEMPI(Q);
				}
			}
			/*
			Now each element to the right of and below elt(Nptr, 0, 0) is zero.
			*/
			k = MATSEEKI(0, Nptr);
			if (k != Nptr->R)
			{
				for (t = 0; t < Nptr->C; t++)
				{
					if ((Nptr->V[k][t])->S != 0)
					{
						Tmp = Nptr->V[0][t];
						Nptr->V[0][t] = ADDI(Nptr->V[0][t], Nptr->V[k][t]);
						/*
						MAXI = UPDATEMAXI(MAXI, Nptr);
						*/
						FREEMPI(Tmp);
					}
				}
			}
			else
				break;
		}
		/*
		Now elt(Nptr, 0, 0) divides all elements elt(Nptr, k, l), k,l > 0.
		*/
		N[i] = COPYI(Nptr->V[0][0]);

		if (i < z)
		{
			Temp = Nptr;
			Nptr = DELETE_ROW1I(Nptr);
			FREEMATI(Temp);
			Temp = Nptr;
			Nptr = DELETE_COL1I(Nptr);
			FREEMATI(Temp);
		}
		size = Nptr->R;
		Tmp = MAXELTI(Nptr);
		printf("found invariant factor d[%u] = ", i + 1); PRINTI(N[i]); printf("\n");
		/*
		printf("MAXI = "); PRINTI(MAXI); printf("\n");
		printf("MAX_ENTRY(Nptr) = "); PRINTI(Tmp); printf("\n");
		printf("MAXI has %u digits\n", LENGTHI(MAXI));
		*/
		l = RSV(Tmp, Eptr);
		FREEMPI(Tmp);
		if (l == 1 && !flag)/* we don't want to use MLLL if there are
							1's or -1's still present */
		{
			ZEROTESTI(0, Nptr, &r, &s);
			if (r > 0)/* must have a non-zero row 1 for MLLL */
			{
				printf("r = %u > 0\n", r);
				Nptr = SWAP_ROWSI1(0, r, Nptr);
			}
			Temp = Nptr;
			Nptr = BASIS_REDUCTION0(Nptr, m1, n1);
			FREEMATI(Temp);
		}
		size = size - Nptr->R;
		m = m - size;
		z = MIN((int)m - 1, (int)n - 1);
	}
	FREEMATI(Nptr);
	*ptr = i++;
	/*
	printf("MAXI = "); PRINTI(MAXI); printf("\n");
	FREEMPI(MAXI);
	*/
	return (N);
}

MPMATI *DELETE_ROW1I(MPMATI *Mptr)
/*
* deletes the first row of *Mptr.
*/
{
	unsigned int i, j, m, n;
	MPMATI *Nptr;

	m = Mptr->R;
	n = Mptr->C;
	Nptr = BUILDMATI(m - 1, n);
	for (i = 0; i <= Nptr->R - 1; i++)
	{
		for (j = 0; j <= n - 1; j++)
			Nptr->V[i][j] = COPYI(Mptr->V[i + 1][j]);
	}
	return (Nptr);
}

MPMATI *DELETE_COL1I(MPMATI *Mptr)
/*
* deletes  the first column of *Mptr.
*/
{
	unsigned int i, j, m, n;
	MPMATI *Nptr;

	m = Mptr->R;
	n = Mptr->C;
	Nptr = BUILDMATI(m, n - 1);
	for (j = 0; j <= Nptr->C - 1; j++)
	{
		for (i = 0; i <= m - 1; i++)
			Nptr->V[i][j] = COPYI(Mptr->V[i][j + 1]);
	}
	return (Nptr);
}

unsigned int PIVOTI(MPMATI *Mptr, USI *rptr, USI *cptr, USI r, USI s)
{
	unsigned int i, j, m, n;

	m = Mptr->R;
	n = Mptr->C;
	*rptr = r; *cptr = s;
	for (i = r; i < m; i++)
	{
		for (j = 0; j < n; j++)
		{
			if ((Mptr->V[i][j])->D == 0 && (Mptr->V[i][j])->V[0] == 1)
			{
				*rptr = i; *cptr = j;
				return (1);
			}
		}
	}
	return (0);
}
