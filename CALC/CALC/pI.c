/* pI.c */

#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_NONSTDC_NO_DEPRECATE
#include <stdio.h>
#include <stdlib.h>
#include "integer.h"
#include "pI.h"
#include "fun.h"


#ifndef DEBUG

POLYI BANKI = NULL;

#endif

#ifndef DEBUG
void SAVEI(Ptr)
POLYI Ptr;
/*
This routine takes a POLYI (which is a linked list of TERMI's) that is
no longer needed, and links it onto the start of the free list BANKI,
which itself is just a linked list opf TERMI's.
*/
{
	POLYI tmp;

	tmp = Ptr;
	while (tmp->NEXT)
		tmp = tmp->NEXT;
	tmp->NEXT = BANKI;
	BANKI = Ptr;
	Ptr = NULL;
}
#endif

POLYI CREATEPI()
/*
* If there is anything in BANKI (the memory bank of TERMI's) then use that,
* otherwise allocate a new one.
* No initialisation is performed
*/
{
	POLYI Ptr;
#ifndef DEBUG
	if (BANKI)
	{
		Ptr = BANKI;
		BANKI = BANKI->NEXT;
	}
	else
#endif
		Ptr = (POLYI)mmalloc(sizeof(struct TERMI));
	return (Ptr);
}

int DEGREEPI(POLYI Pptr)
/*
* returns the degree of the polynomial Pptr.
*/
{
	if (Pptr == NULL)
		return (-1);
	else
		return (Pptr->DEG);
}

MPI *LEADPI(POLYI Pptr)
/*
* *Aptr is the  leading coefficient of the polynomial Pptr.
*/
{
	if (Pptr == NULL)
		return ZEROI();
	else
		return COPYI(Pptr->COEF);
}

void DELETEPI(POLYI Pptr)
/*
* returns back to the system storage dynamically allocated to the linked list
* which Pptr points to.
*/
{
	POLYI Qtmp, orig;

	if (!Pptr) return;
	orig = Pptr;
	while (Pptr != NULL)
	{
		Qtmp = Pptr->NEXT;
		FREEMPI(Pptr->COEF);
#ifdef DEBUG /* This fix was due to Peter Adams on 30/6/94 */
		SAVEPI(Pptr);
#endif
		Pptr = Qtmp;
	}
#ifndef DEBUG
	SAVEPI(orig);
#endif
	return;
}

void PINSERTPI(int DEG, MPI *COEF, POLYI *Pptr, USI OPT)
/*
* inserts a copy of a single term *Qptr into polynomial *Pptr.
* It has two options that control the insertion in case a term of *Pptr already
* has the degree of *Qptr, If OPT = 0, then the term is replaced; if OPT = 1,
* then the coefficient of *Qptr is added to that of the corresponding term of
*  *Pptr.
*/
{
	POLYI CURSOR, Ttmp;
	MPI *TempI;

	if (*Pptr == NULL)
	{
		*Pptr = CREATEPI();
		(*Pptr)->DEG = DEG;
		(*Pptr)->COEF = COPYI(COEF);
		(*Pptr)->NEXT = NULL;
		return;
	}
	else
	{
		CURSOR = *Pptr;
		while ((CURSOR->DEG > DEG) && (CURSOR->NEXT != NULL))
			CURSOR = CURSOR->NEXT;
		if (CURSOR->DEG < DEG)
		{
			/* Qptr inserted before CURSOR */
			Ttmp = CREATEPI();
			Ttmp->DEG = CURSOR->DEG;
			Ttmp->COEF = CURSOR->COEF;
			Ttmp->NEXT = CURSOR->NEXT;
			CURSOR->DEG = DEG;
			CURSOR->COEF = COPYI(COEF);
			CURSOR->NEXT = Ttmp;
			return;
		}
		else if (CURSOR->DEG == DEG)
		{
			if (OPT == 0)
				CURSOR->COEF = COPYI(COEF);
			else
			{
				TempI = CURSOR->COEF;
				CURSOR->COEF = ADDI(CURSOR->COEF, COEF);
				FREEMPI(TempI);
			}
			return;
		}
		else
		{
			/* CURSOR->DEG > Qptr->DEG and CURSOR->NEXT === NULL */
			/* Qptr inserted at end */
			Ttmp = CREATEPI();
			Ttmp->DEG = DEG;
			Ttmp->COEF = COPYI(COEF);
			Ttmp->NEXT = NULL;
			CURSOR->NEXT = Ttmp;
			return;
		}
	}
}

void PURGEPI(POLYI *Pptr)
/*
* gets rid of any zero coefficients in the polynomial Pptr.
*/
{
	POLYI R1tmp, R2tmp;

	if (*Pptr != NULL)
	{
		R1tmp = *Pptr;
		R2tmp = (*Pptr)->NEXT;
		while (R2tmp != NULL)
		{
			/* purge non-leading terms */
			if ((R2tmp->COEF)->S == 0)
			{
				R1tmp->NEXT = R2tmp->NEXT;
				R2tmp->NEXT = NULL;
				DELETEPI(R2tmp);
				R2tmp = R1tmp->NEXT;
			}
			else
			{
				R1tmp = R2tmp;
				R2tmp = R2tmp->NEXT;
			}
		}
		if (((*Pptr)->COEF)->S == 0)
		{
			/* purge leading term */
			R1tmp = *Pptr;
			*Pptr = (*Pptr)->NEXT;
			R1tmp->NEXT = NULL;
			DELETEPI(R1tmp);
		}
	}
	return;
}

POLYI COPYPI(POLYI Pptr)
/*
* makes a copy of the polynomial Pptr.
*/
{
	POLYI Rtmp, Sptr;

	Sptr = NULL;
	Rtmp = Pptr;
	while (Rtmp != NULL)
	{
		PINSERTPI(Rtmp->DEG, Rtmp->COEF, &Sptr, 0);
		Rtmp = Rtmp->NEXT;
	}
	return (Sptr);
}

void PRINTPI(POLYI Pptr)
{
	FPRINTPI(stdout, Pptr);
}

void FPRINTPI(FILE *outfile, POLYI Pptr)
/*
* outputs the terms of the polynomial Pptr.
*/
{
	POLYI CURSOR;
	int k, l;

	if (Pptr == NULL)
	{
		fprintf(outfile, "0");
		return;
	}
	else
	{
		CURSOR = Pptr;
		k = DEGREEPI(Pptr);
		while (CURSOR != NULL)
		{
			l = (CURSOR->COEF)->S;
			if (l == -1)
			{
				(CURSOR->COEF)->S = 1;
				fprintf(outfile, "- ");
			}
			else
			{
				if (CURSOR->DEG != k)
					fprintf(outfile, "+ ");
			}
			if (CURSOR->DEG == 0 || !EQONEI(CURSOR->COEF))
				FPRINTI(outfile, CURSOR->COEF);
			if (l == -1)
				(CURSOR->COEF)->S = -1;
			if (CURSOR->DEG == 1)
				fprintf(outfile, "X");
			if (CURSOR->DEG > 1)
				fprintf(outfile, "X^%d", CURSOR->DEG);
			if (CURSOR->NEXT != NULL)
				fprintf(outfile, " ");
			CURSOR = CURSOR->NEXT;
		}
		return;
	}
}



POLYI SCALARPI(MPI *Cptr, POLYI Pptr)
/*
* multiplying the polynomial Pptr by the MPI *Cptr.
*/
{
	POLYI Sptr, Rtmp;
	MPI *TempI;

	if (Pptr == NULL || Cptr->S == 0)
		Sptr = NULL;
	else
	{
		Sptr = COPYPI(Pptr);
		Rtmp = Sptr;
		while (Rtmp != NULL)
		{
			TempI = Rtmp->COEF;
			Rtmp->COEF = MULTI(Rtmp->COEF, Cptr);
			FREEMPI(TempI);
			Rtmp = Rtmp->NEXT;
		}
	}
	return (Sptr);
}

POLYI ADDPI(POLYI Pptr, POLYI Qptr)
/*
* returns the sum of two polynomials.
*/
{
	POLYI CURSOR, Sptr;

	Sptr = COPYPI(Pptr);
	CURSOR = Qptr;
	while (CURSOR != NULL)
	{
		PINSERTPI(CURSOR->DEG, CURSOR->COEF, &Sptr, 1);
		CURSOR = CURSOR->NEXT;
	}
	PURGEPI(&Sptr);
	return (Sptr);
}

POLYI SUBPI(POLYI Pptr, POLYI Qptr)
/*
* returns the difference Pptr - Qptr of two polynomials.
*/
{
	POLYI CURSOR, Sptr;
	MPI *TempI;

	Sptr = COPYPI(Pptr);
	if (Qptr == NULL)
		return (Sptr);
	CURSOR = Qptr;
	while (CURSOR != NULL)
	{
		TempI = MINUSI(CURSOR->COEF);
		PINSERTPI(CURSOR->DEG, TempI, &Sptr, 1);
		FREEMPI(TempI);
		CURSOR = CURSOR->NEXT;
	}
	PURGEPI(&Sptr);
	return (Sptr);
}

POLYI MULTPI(POLYI Pptr, POLYI Qptr)
/*
* returns the product of two polynomials.
*/
{
	POLYI Sptr, CURSOR_P, CURSOR_Q;
	MPI *TempI;
	int k;

	if (Pptr == NULL || Qptr == NULL)
		return (NULL);
	Sptr = NULL;
	CURSOR_P = Pptr;
	while (CURSOR_P != NULL)
	{
		CURSOR_Q = Qptr;
		while (CURSOR_Q != NULL)
		{
			k = CURSOR_P->DEG + CURSOR_Q->DEG;
			TempI = MULTI(CURSOR_P->COEF, CURSOR_Q->COEF);
			PINSERTPI(k, TempI, &Sptr, 1);
			FREEMPI(TempI);
			CURSOR_Q = CURSOR_Q->NEXT;
		}
		CURSOR_P = CURSOR_P->NEXT;
	}
	PURGEPI(&Sptr);
	return (Sptr);
}

MPI *VALPI(POLYI Pptr, MPI *Cptr)
/*
* Evaluating the polynomial Pptr at the MPI *Cptr.
*/
{
	int k;
	MPI *B, *Aptr, *TempI;
	POLYI CURSOR;

	if (Pptr == NULL)
		return ZEROI();
	else
		Aptr = COPYI(Pptr->COEF);
	CURSOR = Pptr;
	while (CURSOR != NULL)
	{
		if (CURSOR->NEXT == NULL)
			k = DEGREEPI(CURSOR);
		else
			k = DEGREEPI(CURSOR) - DEGREEPI(CURSOR->NEXT);
		B = POWERI(Cptr, (unsigned int)k);
		TempI = Aptr;
		Aptr = MULTI(Aptr, B);
		FREEMPI(B);
		FREEMPI(TempI);
		CURSOR = CURSOR->NEXT;
		if (CURSOR != NULL)
		{
			TempI = Aptr;
			Aptr = ADDI(Aptr, CURSOR->COEF);
			FREEMPI(TempI);
		}
	}
	return (Aptr);
}

POLYI FGPI(POLYI Pptr, POLYI Qptr)
/*
* Calculating Pptr(Qptr).
*/
{
	int k;
	MPI *Aptr;
	POLYI CURSOR, B, C, OUTPUT, TEMP;

	if (Pptr == NULL)
		return ZEROPI();
	else {
		OUTPUT = ZEROPI();
		Aptr = Pptr->COEF;
	}
	CURSOR = Pptr;
	while (CURSOR != NULL)
	{
		k = DEGREEPI(CURSOR);
		B = POWERPI(Qptr, (unsigned int)k);
		C = SCALARPI(Aptr, B);
		DELETEPI(B);
		TEMP = OUTPUT;
		OUTPUT = ADDPI(C, OUTPUT);
		DELETEPI(C);
		DELETEPI(TEMP);
		CURSOR = CURSOR->NEXT;
		if (CURSOR != NULL)
		{
			Aptr = CURSOR->COEF;
		}
	}
	return (OUTPUT);
}

POLYI DIVPI(POLYI Fptr, POLYI Gptr)
/*
* Pseudo-division of polynomials: see Knuth, Vol 2, p.407
* and H. Flanders, Scientific Pascal, p.510.
*/
{
	POLYI F, Tmp, Htmp, CURSOR, Qptr;
	int j, d;
	MPI *A, *G, *TempI;

	Qptr = NULL;
	d = DEGREEPI(Gptr);
	j = DEGREEPI(Fptr) - d;
	if (j < 0)
		return (Qptr);
	else
	{
		F = COPYPI(Fptr);
		G = LEADPI(Gptr);
		A = POWERI(G, (unsigned int)(j + 1));
		Tmp = F;
		F = SCALARPI(A, Tmp);
		FREEMPI(A);
		DELETEPI(Tmp);
		while ((j = DEGREEPI(F)) >= d)
		{
			j = j - d;
			A = LEADPI(F);
			TempI = A;
			A = INTI(A, G);
			FREEMPI(TempI);
			Htmp = COPYPI(Gptr);
			CURSOR = Htmp;
			do
			{
				CURSOR->DEG = j + CURSOR->DEG;
				TempI = CURSOR->COEF;
				CURSOR->COEF = MINUSI(CURSOR->COEF);
				FREEMPI(TempI);
				TempI = CURSOR->COEF;
				CURSOR->COEF = MULTI(CURSOR->COEF, A);
				FREEMPI(TempI);
				CURSOR = CURSOR->NEXT;
			} while (CURSOR != NULL);
			Tmp = F;
			F = ADDPI(Tmp, Htmp);
			DELETEPI(Tmp);
			DELETEPI(Htmp);
			PINSERTPI(j, A, &Qptr, 0);
			FREEMPI(A);
		}
		DELETEPI(F);
		FREEMPI(G);
		return (Qptr);
	}
}

POLYI MODPI(POLYI Fptr, POLYI Gptr)
/*
* Pseudo-division of polynomials: see Knuth, Vol 2, p.407
* and H. Flanders, Scientific Pascal, p.510.
* NOTE: This returns a polynomial that is a POSITIVE multiple of the true remainder.
* i.e. The sign of the remainder is correct.
*/
{
	POLYI F, Tmp, Htmp, CURSOR;
	int j, d, degDiff;
	MPI *A, *G, *Temp1I, *Temp2I;

	d = DEGREEPI(Gptr); /* n */
	j = DEGREEPI(Fptr) - d; /* m - n */
	degDiff = j;
	F = COPYPI(Fptr);
	if (j < 0)
		return (F);
	else
	{
		G = LEADPI(Gptr);
		A = POWERI(G, (unsigned int)(j + 1));
		Tmp = F;
		F = SCALARPI(A, Tmp);
		FREEMPI(A);
		DELETEPI(Tmp);
		while ((j = DEGREEPI(F)) >= d)
		{
			j = j - d;
			A = LEADPI(F);
			if (!EQONEI(G))
			{
				Temp2I = A;
				A = INTI(A, G);
				FREEMPI(Temp2I);
			}
			Htmp = COPYPI(Gptr);
			CURSOR = Htmp;
			do
			{
				CURSOR->DEG = j + CURSOR->DEG;
				Temp1I = CURSOR->COEF;
				CURSOR->COEF = MINUSI(CURSOR->COEF);
				FREEMPI(Temp1I);
				Temp1I = CURSOR->COEF;
				CURSOR->COEF = MULTI(CURSOR->COEF, A);
				FREEMPI(Temp1I);
				CURSOR = CURSOR->NEXT;
			} while (CURSOR != NULL);
			Tmp = F;
			F = ADDPI(Tmp, Htmp);
			DELETEPI(Tmp);
			DELETEPI(Htmp);
			FREEMPI(A);
		}

		if ((degDiff % 2 == 0) && G->S == -1) {
			POLYI TMPI;
			TMPI = F;
			F = MINUSPI(F);
			DELETEPI(TMPI);
		}
		FREEMPI(G);
		return (F);
	}
}

MPI *CONTENTPI(POLYI Pptr)
/*
* *Cptr is the (positive) content of the polynomial Pptr.
*/
{
	POLYI CURSOR;
	MPI *TempI, *Cptr;

	if (Pptr == NULL)
		return (ZEROI());
	else if (Pptr->DEG == 0)
	{
		Cptr = COPYI(Pptr->COEF);
		if (Cptr->S == -1)
			Cptr->S = 1;
		return (Cptr);
	}
	else
	{
		CURSOR = Pptr;
		Cptr = COPYI(CURSOR->COEF);
		if (CURSOR->NEXT != NULL)
		{
			do
			{
				TempI = Cptr;
				Cptr = GCD(Cptr, (CURSOR->NEXT)->COEF);
				FREEMPI(TempI);
				CURSOR = CURSOR->NEXT;
			} while (CURSOR->NEXT != NULL);
		}
		return (Cptr);
	}
}

MPI *CONTENTPI2(POLYI Pptr)
/* Returns content of Pptr such that content(Pptr)*primitive(Pptr) = Pptr */
{
	MPI *C, *C2, *L, *Z = ZEROI();
	C = CONTENTPI(Pptr);
	L = LEADPI(Pptr);
	if (COMPAREI(L, Z) < 0 && COMPAREI(C, Z) > 0) {
		C2 = MINUSI(C);
		FREEMPI(C);
		FREEMPI(L);
		FREEMPI(Z);
		return C2;
	}
	else {
		FREEMPI(L);
		FREEMPI(Z);
		return C;
	}
}


POLYI PRIMITIVEPI(POLYI Pptr)
/*
* returns the primitive part of the polynomial Pptr.
* leading coefficient is positive.
*/
{
	POLYI Tmp, Sptr, CURSOR;
	MPI *C, *TempI;

	if (Pptr == NULL)
		return (NULL);
	else
	{
		Sptr = COPYPI(Pptr);
		C = CONTENTPI(Sptr);
		if ((Sptr->COEF)->S == -1)
		{
			Tmp = Sptr;
			TempI = MINUS_ONEI();
			Sptr = SCALARPI(TempI, Tmp);
			FREEMPI(TempI);
			DELETEPI(Tmp);
		}
		CURSOR = Sptr;
		while (CURSOR != NULL)
		{
			TempI = CURSOR->COEF;
			CURSOR->COEF = INT(CURSOR->COEF, C);
			FREEMPI(TempI);
			CURSOR = CURSOR->NEXT;
		}
		FREEMPI(C);
		return (Sptr);
	}
}

POLYI GCDPI(POLYI Pptr, POLYI Qptr)
/*
* returns the gcd of polynomials Pptr and Qptr.
* see Knuth, Vol 2, p.403-408 and H. Flanders, Scientific Pascal, p.510.
*/
{
	POLYI Ptmp, Qtmp, Rtmp, Sptr;
	MPI *D, *P, *Q, *TempI;
	int j;

	if (Qptr == NULL)
	{
		if (Pptr == NULL)
			return (NULL);
		else
			return (COPYPI(Pptr));
	}
	if (Pptr == NULL)
	{
		if (Qptr == NULL)
			return (NULL);
		else
			return (COPYPI(Qptr));
	}
	P = CONTENTPI(Pptr);
	Q = CONTENTPI(Qptr);
	D = GCD(P, Q);
	FREEMPI(P);
	FREEMPI(Q);
	Ptmp = PRIMITIVEPI(Pptr);
	Qtmp = PRIMITIVEPI(Qptr);
	j = DEGREEPI(Qtmp);
	while ((Qtmp != NULL) && ((j = DEGREEPI(Qtmp)) != 0))
	{
		Rtmp = MODPI(Ptmp, Qtmp);
		DELETEPI(Ptmp);
		Ptmp = Qtmp;
		Qtmp = PRIMITIVEPI(Rtmp);
		DELETEPI(Rtmp);
	}
	if (j == 0)
	{
		TempI = Qtmp->COEF;
		Qtmp->COEF = ONEI();
		FREEMPI(TempI);
	}
	else/* Qtmp = NULL*/
		Qtmp = COPYPI(Ptmp);
	Sptr = SCALARPI(D, Qtmp);
	FREEMPI(D);
	DELETEPI(Ptmp);
	DELETEPI(Qtmp);
	return (Sptr);
}

POLYI DERIVPI(POLYI Pptr)
/*
* returns the derivative of a polynomial.
*/
{
	POLYI CURSOR, Sptr;
	MPI *TempI;

	if (DEGREEPI(Pptr) == 0)
		return (NULL);
	else
	{
		Sptr = COPYPI(Pptr);
		CURSOR = Sptr;
		while (CURSOR != NULL)
		{
			TempI = CURSOR->COEF;
			CURSOR->COEF = MULT_I(CURSOR->COEF, CURSOR->DEG);
			FREEMPI(TempI);
			CURSOR->DEG = CURSOR->DEG - 1;
			if (CURSOR->NEXT == NULL)
				break;
			else if ((CURSOR->NEXT)->DEG == 0)
			{
				DELETEPI(CURSOR->NEXT);
				CURSOR->NEXT = NULL;
				break;
			}
			else
				CURSOR = CURSOR->NEXT;
		}
		return (Sptr);
	}
}

POLYI *SQUAREFREEPI(POLYI Pptr, USI **D, USI *tptr)
/*
* returns the "squarefree decomposition" of Pptr, a non-constant polynomial:
* Pptr = G[0]^D[0]...G[*tpr - 1]^D[*tptr - 1],
* where G[i] is squarefree and is the product of the irreducible factors of
* multiplicity D[i]. See L. Childs, Introduction to Higher Algebra, p.159.
*/
{
	POLYI Ttmp, Ttmp1, Dtmp, Ftmp, Gtmp, Htmp, Ptmp, *G;
	unsigned int i = 1, j = 0, n;

	n = (unsigned int)(Pptr->DEG);
	Ptmp = PRIMITIVEPI(Pptr);
	Ftmp = DERIVPI(Ptmp);
	Dtmp = GCDPI(Ptmp, Ftmp);
	G = (POLYI *)mmalloc(n * sizeof(POLYI));
	*D = (unsigned int *)mmalloc(n * sizeof(unsigned int));
	if (DEGREEPI(Dtmp) == 0)
	{
		G[j] = Ptmp;
		(*D)[j] = i;
		j++;
		DELETEPI(Ftmp);
		DELETEPI(Dtmp);
		*tptr = j;
		return (G);
	}
	else
	{
		Gtmp = DIVPI(Ptmp, Dtmp);
		while (DEGREEPI(Dtmp))
		{
			DELETEPI(Ptmp);
			DELETEPI(Ftmp);
			Htmp = Gtmp;
			Ptmp = Dtmp;
			Ftmp = DERIVPI(Ptmp);
			Dtmp = GCDPI(Ptmp, Ftmp);
			Gtmp = DIVPI(Ptmp, Dtmp);
			Ttmp = DIVPI(Htmp, Gtmp);
			DELETEPI(Htmp);
			Ttmp1 = PRIMITIVEPI(Ttmp);
			if (Ttmp1->DEG > 0)
			{
				G[j] = COPYPI(Ttmp1);
				(*D)[j] = i;
				j++;
			}
			DELETEPI(Ttmp);
			DELETEPI(Ttmp1);
			i++;
		}
		DELETEPI(Ptmp);
		DELETEPI(Ftmp);
		G[j] = PRIMITIVEPI(Gtmp);
		(*D)[j] = i;
		j++;
		DELETEPI(Gtmp);
		DELETEPI(Dtmp);
		*tptr = j;
		return (G);
	}
}

POLYI POWERPI(POLYI Pptr, USI n)
/*
* returns Pptr ^ n, where 0 <= n < R0 * R0.
*/
{
	POLYI B, Eptr, Tmp;

	Eptr = ONEPI();
	if (n == 0)
		return (Eptr);
	B = COPYPI(Pptr);
	while (1)
	{
		if (n & 1)
		{
			Tmp = Eptr;
			Eptr = MULTPI(Tmp, B);
			DELETEPI(Tmp);
			if (n == 1)
			{
				DELETEPI(B);
				return (Eptr);
			}
		}
		Tmp = B;
		B = MULTPI(Tmp, Tmp);
		DELETEPI(Tmp);
		n = n >> 1;
	}
}

POLYI ONEPI()
/*
* returns the constant polynomial 1.
*/
{
	POLYI Pptr;

	Pptr = CREATEPI();
	Pptr->DEG = 0;
	Pptr->COEF = ONEI();
	Pptr->NEXT = NULL;
	return (Pptr);
}
POLYI ZEROPI()
/* returns the zero polynomial */
{
	POLYI Pptr;
	Pptr = CREATEPI();
	Pptr->DEG = 0;
	Pptr->COEF = ZEROI();
	Pptr->NEXT = NULL;
	return (Pptr);
}
POLYI CONSTANTPI(MPI *Aptr)
/*
* returns the constant polynomial whose constant term is *Aptr.
*/
{
	POLYI Pptr;

	if (Aptr->S == 0)
		return (ZEROPI());
	else
	{
		Pptr = CREATEPI();
		Pptr->DEG = 0;
		Pptr->COEF = COPYI(Aptr);
		Pptr->NEXT = NULL;
		return (Pptr);
	}
}

MPI *SUBRESULTANT(POLYI Pptr, POLYI Qptr)
/*
* *Aptr is the resultant of Pptr and Qptr, found using the subresultant
* algorithm, p. 130. E. Kaltofen, G. E. Collins etc, Algorithm 4.5.
* similar to Knuth, Algorithm C, p. 410.
*/
{
	unsigned delta;
	MPI *G, *G1, *H1, *H2, *Aptr, *TempI;
	POLYI Utmp, Vtmp, Rtmp, CURSOR;

	G = ONEI();
	Aptr = ONEI();
	Utmp = COPYPI(Pptr);
	Vtmp = COPYPI(Qptr);
	while (Vtmp != NULL)
	{
		Rtmp = MODPI(Utmp, Vtmp);
		if (Rtmp == NULL && DEGREEPI(Vtmp) > 0)
		{
			FREEMPI(G);
			FREEMPI(Aptr);
			DELETEPI(Utmp);
			DELETEPI(Vtmp);
			return (ZEROI());
		}
		delta = DEGREEPI(Utmp) - DEGREEPI(Vtmp);
		DELETEPI(Utmp);
		Utmp = Vtmp;
		G1 = MINUSI(G);
		H1 = MINUSI(Aptr);
		TempI = H1;
		H1 = POWERI(H1, delta);
		FREEMPI(TempI);
		TempI = H1;
		H1 = MULTI(H1, G1);
		FREEMPI(TempI);
		CURSOR = Rtmp;
		while (CURSOR != NULL)
		{
			TempI = CURSOR->COEF;
			CURSOR->COEF = INTI(CURSOR->COEF, H1);
			FREEMPI(TempI);
			CURSOR = CURSOR->NEXT;
		}
		Vtmp = Rtmp;
		TempI = G;
		G = LEADPI(Utmp);
		FREEMPI(TempI);
		TempI = G1;
		G1 = POWERI(G, delta);
		FREEMPI(TempI);
		TempI = H1;
		H1 = MULTI(Aptr, G1);
		FREEMPI(TempI);
		H2 = POWERI(Aptr, delta);
		TempI = Aptr;
		Aptr = INTI(H1, H2);
		FREEMPI(TempI);
		FREEMPI(H1);
		FREEMPI(H2);
		FREEMPI(G1);
	}
	DELETEPI(Utmp);
	FREEMPI(G);
	return (Aptr);
}

MPI *DISCRIMINANTPI(POLYI Pptr)
/*
* *Dptr = Discriminant of Pptr = (1/a[n])*(-1)^{n*(n-1)/2 * Res(Pptr, Pptr').
* See O. Perron, Algebra, Vol 1, p.212.
*/
{
	POLYI P;
	MPI *A, *Dptr, *TempI;
	int n;

	P = DERIVPI(Pptr);
	Dptr = SUBRESULTANT(Pptr, P);
	A = LEADPI(Pptr);
	TempI = Dptr;
	Dptr = INTI(Dptr, A);
	FREEMPI(TempI);
	FREEMPI(A);
	DELETEPI(P);
	n = (DEGREEPI(Pptr)) % 4;
	if (n == 2 || n == 3)
		Dptr->S = -(Dptr->S);
	return (Dptr);
}

MPI *LENGTHSQPI(POLYI Pptr)
/*
* *Cptr is the sum of the squares of the coefficients of Pptr.
*/
{
	POLYI CURSOR;
	MPI *C, *Cptr, *TempI;

	Cptr = ZEROI();
	if (Pptr == NULL)
		return (Cptr);
	CURSOR = Pptr;
	while (CURSOR != NULL)
	{
		C = MULTI(CURSOR->COEF, CURSOR->COEF);
		TempI = Cptr;
		Cptr = ADDI(C, Cptr);
		FREEMPI(C);
		FREEMPI(TempI);
		CURSOR = CURSOR->NEXT;
	}
	return (Cptr);
}

MPI *VAL0PI(POLYI Pptr)
/*
* *Mptr = Pptr(0)
*/
{
	POLYI CURSOR;

	if (Pptr == NULL)
		return ZEROI();
	CURSOR = Pptr;
	while (CURSOR->NEXT != NULL)
		CURSOR = CURSOR->NEXT;
	if (DEGREEPI(CURSOR) != 0)
		return ZEROI();
	else
		return COPYI(CURSOR->COEF);
}

/* Returns respectively positive, or negative if the sign of the polynomial
* evaluated at the supplied rational number is positive or negative.
* Returns zero if polynomial has a root at the supplied rational */
int SIGN_AT_R_PI(POLYI P, MPR *R)
{
	MPI *TMPI, *RESULT = ZEROI(), *X_TO_I, *TERM;
	MPI *A = COPYI(R->N), *D = COPYI(R->D);
	POLYI CURSOR = P;
	unsigned int n = DEGREEPI(P), i;

	if (P == NULL)
		i = 0;
	else while (CURSOR) {
		i = DEGREEPI(CURSOR);
		TERM = POWERI(A, i);

		TMPI = TERM;
		TERM = MULTI(CURSOR->COEF, TERM);
		FREEMPI(TMPI);

		X_TO_I = POWERI(D, n - i);

		TMPI = TERM;
		TERM = MULTI(TERM, X_TO_I);
		FREEMPI(TMPI);
		FREEMPI(X_TO_I);

		TMPI = RESULT;
		RESULT = ADDI(RESULT, TERM);
		FREEMPI(TMPI);
		FREEMPI(TERM);

		CURSOR = CURSOR->NEXT;
	}
	i = SIGNI(RESULT); /* shouldn't really use this as return val but it works */
	FREEMPI(RESULT);
	FREEMPI(A);
	FREEMPI(D);
	return i;
}

/* If P = p(x) then this function returns p(-x) */
POLYI P_OF_NEG_X(POLYI P)
{
	POLYI Q = COPYPI(P);
	POLYI CURSOR = Q;
	while (CURSOR) {
		if (DEGREEPI(CURSOR) % 2 == 1) /* if odd degree */
			CURSOR->COEF->S = -CURSOR->COEF->S;
		CURSOR = CURSOR->NEXT;
	}
	return Q;
}

/* Returns -1*P(x) */
POLYI MINUSPI(POLYI P)
{
	POLYI Q = COPYPI(P);
	POLYI CURSOR = Q;
	while (CURSOR) {
		CURSOR->COEF->S = -CURSOR->COEF->S;
		CURSOR = CURSOR->NEXT;
	}
	return Q;
}

/* Modifies the supplied polynomial so that p(x) = p(x+A) */

void P_OF_X_PLUS(POLYI *P, MPI *A)
{
	MPI *ZERO = ZEROI();
	MPIA COEF = BUILDMPIA();
	POLYI CURSOR = *P;
	int i, j, n = DEGREEPI(*P);

	for (i = n; i >= 0; i--) {
		ADD_TO_MPIA(COEF, ZERO, i);
	}
	while (CURSOR) {
		ADD_TO_MPIA(COEF, CURSOR->COEF, DEGREEPI(CURSOR));
		CURSOR = CURSOR->NEXT;
	}

	for (i = 0; i <= n - 1; i++)
		for (j = n - 1; j >= i; j--) {
			MPI *TMP, *MUL;
			MUL = MULTI(A, COEF->A[j + 1]);
			TMP = COEF->A[j];
			COEF->A[j] = ADDI(MUL, COEF->A[j]);
			FREEMPI(TMP);
			FREEMPI(MUL);
		}

	DELETEPI(*P);
	*P = NULL;
	for (i = 0; i <= n; i++)
		PINSERTPI(i, COEF->A[i], P, 0);
	PURGEPI(P);
	FREEMPIA(COEF);
	FREEMPI(ZERO);
}

POLYI PRIMITIVEPI_(POLYI Pptr)
/*
* returns the primitive part of the polynomial Pptr.
* retaining the sign of the leading coefficient.
*/
{
	POLYI Sptr, CURSOR;
	MPI *C, *TempI;

	if (Pptr == NULL)
		return (NULL);
	else
	{
		Sptr = COPYPI(Pptr);
		C = CONTENTPI(Sptr);
		CURSOR = Sptr;
		while (CURSOR != NULL)
		{
			TempI = CURSOR->COEF;
			CURSOR->COEF = INT(CURSOR->COEF, C);
			FREEMPI(TempI);
			CURSOR = CURSOR->NEXT;
		}
		FREEMPI(C);
		return (Sptr);
	}
}
