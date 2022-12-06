/* program:  "i2.c" */
/*
* A translation into C of the long division section of pascal program
* ARITHMETIC,  from "SCIENTIFIC PASCAL" by H. FLANDERS pp. 342-357.
*/

#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_NONSTDC_NO_DEPRECATE
#include <stdio.h>
#include <stdlib.h>
#include "integer.h"
#include "fun.h"

unsigned long l0, n0, bn, bn1;

void INIT(MPI *Aptr, MPI *Bptr, MPI **A1ptr, MPI **B1ptr)
/*
* see lemmas 2 and 3, pp. 349-350.
* input: Aptr, Bptr are pointers to MPI'S, Aptr->D >=  Bptr->D = n0 >= 1.
*  output: l0 = int[R0 / (Bptr->V[n0] + 1), *A1ptr = l0 * *Aptr,
* *B1ptr = l0 * *Bptr, bn = B1ptr->V[n0], bn1 = B1ptr->V[n0 - 1].
*/
{
	n0 = Bptr->D;
	l0 = R0 / (Bptr->V[n0] + 1);
	if (l0 > 1)
	{
		*A1ptr = MULT_I(Aptr, (long)(l0));
		*B1ptr = MULT_I(Bptr, (long)(l0));
	}
	else
	{
		*A1ptr = COPYI(Aptr);
		*B1ptr = COPYI(Bptr);
	}
	bn = (*B1ptr)->V[n0];
	bn1 = (*B1ptr)->V[n0 - 1];
	return;
}

unsigned long qp;
void FIND(unsigned int k, MPI *A1ptr, MPI *B1ptr, unsigned long R[])
/*
* finds the k-th trial quotient digit qp = Q' using lemmas 1-4.
*/
{
	unsigned long t, u, rp;

	if (k == A1ptr->D - B1ptr->D)
		R[A1ptr->D + 1] = 0; /* an extra 0 is added to the front of
							 the initial dividend to ensure qp < R0.(see p348) */

	u = (R[k + B1ptr->D + 1] << T0) + R[k + B1ptr->D];
	t = u / bn;
	if (t < R0)
		qp = t;
	else
		qp = R0 - 1;

	rp = u - qp * bn; /* Flanders erroneously has t instead of qp */

	if (rp >= R0) /* rp * R0 is then >= 2 ^ O32 > bn1 * qp. */
		return;
	while (bn1 * qp >(rp << T0) + R[k + n0 - 1])
	{	/* see lemma 4(1) */
		qp--;
		rp = rp + bn;
		if (rp >= R0) /* rp * R0 is then >= 2 ^ O32 > bn1 * qp. */
			return;
	}
	return;
}

void NEXT(int k, MPI *A1ptr, MPI *B1ptr, unsigned long R[])
/*
* the run of n0 + 1 "digits" in the current remainder, from
* the (k + n0 + 1)-th down to the k-th, is divided by *B1ptr. The quotient
* becomes the k-th "digit" of *Qptr in INT0() and the remainder R is updated.
*/
{
	long l;
	unsigned int j;
	unsigned long t, ct = 0, cr = 1;	/* initialize "carries */

	for (j = 0; j <= n0; j++)
	{	/* multiply with carry */
		t = (B1ptr->V[j]) * qp + ct;
		ct = t >> T0;
		t = t & (R0 - 1); 	/* digit to subract with carry */
		t = R[k + j] - t + R0 + cr - 1;
		R[k + j] = t & (R0 - 1);
		cr = t >> T0;
	}
	if (k == A1ptr->D - B1ptr->D)
		R[A1ptr->D + 1] = 0; /* avoids garbage value for R[k + n0 + 1] */
	l = R[k + n0 + 1] - ct + cr - 1;
	if (l < 0) { /* qp is one too large */
		qp--;
		cr = 0;
		for (j = 0; j <= n0; j++)
		{ /* add *B1ptr back */
			t = R[k + j] + B1ptr->V[j] + cr;
			R[k + j] = t & (R0 - 1);
			cr = t >> T0;
		}
		R[k + n0 + 1] = l + cr;
	}
	return;

}


MPI *INT0_(MPI *Aptr, unsigned long b)
/*
* input: Aptr is a pointer to an non-negative MPI, b is a positive integer,
* b < R0. output:  *Qptr = int(*Aptr / b).
*/
{
	unsigned int k, n;
	unsigned long t;
	MPI *Qptr, *tmp;

	if (b >= R0)
	{
		fprintf(stderr, "in INT0_I, %lu >= R0\n", b);
		exit(1);
	}
	Qptr = COPYI(Aptr);
	n = Qptr->D;
	t = Qptr->V[Qptr->D];
	if (Qptr->D != 0)
	{
		for (k = Qptr->D; k >= 1; k--)
		{
			Qptr->V[k] = t / b;
			t = Qptr->V[k - 1] + ((t % b) << T0);
		}
		if (Qptr->V[Qptr->D] == 0)
		{

			tmp = BANK_REALLOC(Qptr, n - 1);
			FREEMPI(Qptr);
			Qptr = tmp;
			/* we remark that Qptr->S = 1 here. */
		}
	}
	Qptr->V[0] = t / b;
	if (Qptr->D == 0)
		Qptr->S = Qptr->V[0] ? 1 : 0;
	return Qptr;
}

unsigned long MOD0_(MPI *Aptr, unsigned long b)
/*
* input: Aptr is a pointer to a non-negative MPI, b is a positive integer,
* b < R0.  output: *Aptr mod b.
*/
{
	MPI *R;
	unsigned int k;
	unsigned long t;
	if (b >= R0)
	{
		fprintf(stderr, "in MOD0_, %lu >= R0\n", b);
		exit(1);
	}
	R = COPYI(Aptr);
	t = R->V[R->D];
	if (R->D != 0)
	{
		for (k = R->D; k >= 1; k--)
		{
			R->V[k] = t / b;
			t = R->V[k - 1] + ((t % b) << T0);
		}
	}

	FREEMPI(R);
	return (t % b);
}

MPI *INT_(MPI *Aptr, unsigned long b)
/*
* input: Aptr is a pointer to an MPI, b is a positive integer, b < R0.
* output:  *Qptr = int(*Aptr, b).
*/
{
	MPI *Q, *Qptr;
	unsigned long r;

	if (b >= R0)
	{
		fprintf(stderr, "in INT_, %lu >= R0\n", b);
		exit(1);
	}
	if (Aptr->S >= 0)
		return (INT0_(Aptr, b));
	else
	{
		Aptr->S = 1;
		Q = INT0_(Aptr, b);
		r = MOD0_(Aptr, b);
		Aptr->S = -1;
		/* -*Aptr = *Q * b + r, 0 <= r < b,
		*Aptr = (-*Q - 1) * b + b - r */
		if (r)
		{/* add 1 to |*Q|, but not changing the sign of *Q */
			Qptr = ADD0_I(Q, (USL)1);
			Qptr->S = -1;
			FREEMPI(Q);
			return Qptr;
		}
		else
		{					/* r = 0 here */
			Q->S = -1;
			Qptr = COPYI(Q);
			FREEMPI(Q);
			return Qptr;
		}
	}
}

unsigned long MOD_(MPI *Aptr, unsigned long b)
/*
* input: Aptr is a pointer to an MPI, b is a positive integer, b < R0.
* output: *Aptr mod b.
*/
{
	unsigned long r;

	if (b >= R0)
	{
		fprintf(stderr, "in MOD_, %lu >= R0\n", b);
		exit(1);
	}
	if (Aptr->S >= 0)
		return (MOD0_(Aptr, b));
	else
	{
		Aptr->S = 1;
		r = MOD0_(Aptr, b);
		Aptr->S = -1;
		/* -*Aptr = *Q * b + r, 0 <= r < b,
		*Aptr = (-*Q - 1) * b + b - r */
		if (r)
			return (b - r);
		else
			return (0);
	}
}

MPI *INT0(MPI *Aptr, MPI *Bptr)
/*
* input: Aptr, Bptr, pointers to non-negative MPI'S, *Bptr > 0.
* output: *Qptr = int(*Aptr / *Bptr), *Rptr = *Aptr mod *Bptr.
*/
{
	MPI *A1, *B1, *Qptr, *tmp;
	int k;
	unsigned int e, f, degree;
	unsigned long *R;

	if (Bptr->D == 0)
		return (INT0_(Aptr, Bptr->V[0]));
	if (Aptr->D < Bptr->D)
		return ZEROI();
	INIT(Aptr, Bptr, &A1, &B1);
	f = 2 + A1->D;
	R = (unsigned long *)mmalloc((USL)(f * sizeof(unsigned long)));
	/* we malloc an extra place needed in FIND(). */
	for (k = 0; k <= A1->D; k++)
		R[k] = A1->V[k];
	e = A1->D - B1->D;
	Qptr = BUILDMPI(e + 1);
	for (k = e; k >= 0; k--)
	{
		FIND((unsigned int)k, A1, B1, R);
		NEXT(k, A1, B1, R);
		Qptr->V[k] = qp;
	}
	ffree((char *)R, f * sizeof(unsigned long));
	Qptr->S = 0;			/* finding Qptr->S and Qptr->D */
	degree = 0;
	k = e;
	FREEMPI(A1);
	FREEMPI(B1);
	while (k >= 0)
	{
		if (Qptr->V[k] != 0)
		{
			Qptr->S = 1;
			degree = k;
			k = -1;
		}
		else
			k--;
	}
	tmp = Qptr;
	Qptr = BANK_REALLOC(Qptr, degree);
	FREEMPI(tmp);
	return Qptr;
}

MPI *MOD0(MPI *Aptr, MPI *Bptr)
/*
* input: Aptr, Bptr, pointers to MPI'S, *Aptr >= 0, *Bptr > 0.
* output: *Rptr = *Aptr mod *Bptr.
*/
{
	MPI *A1, *B1, *Rptr, *Temp;
	int k, s;
	unsigned int d, e, f;
	unsigned long *R;

	if (Bptr->D == 0)
	{
		Rptr = BUILDMPI(1);
		Rptr->V[0] = MOD0_(Aptr, Bptr->V[0]);
		Rptr->S = Rptr->V[0] ? 1 : 0;

		return Rptr;
	}
	if (Aptr->D < Bptr->D)
		return COPYI(Aptr);
	INIT(Aptr, Bptr, &A1, &B1);
	f = 2 + A1->D;
	R = (unsigned long *)mmalloc((USL)(f * sizeof(unsigned long)));
	/* we malloc an extra place needed in FIND(). */
	for (k = 0; k <= A1->D; k++)
		R[k] = A1->V[k];
	e = A1->D - B1->D;
	for (k = e; k >= 0; k--)
	{
		FIND((unsigned int)k, A1, B1, R);
		NEXT(k, A1, B1, R);
	}
	FREEMPI(A1);
	FREEMPI(B1);

	/* finding Rptr->S and Rptr->D */

	s = 0;
	d = 0;
	k = n0;
	while (k >= 0)
	{
		if (R[k] != 0)
		{
			s = 1;
			d = k;
			k = -1;
		}
		else
			k--;
	}
	Rptr = BUILDMPI(1 + d);
	Rptr->S = s;
	for (k = 0; k <= Rptr->D; k++)
		Rptr->V[k] = R[k];
	ffree((char *)R, f * sizeof(unsigned long));
	Temp = Rptr;
	Rptr = INT0_(Temp, l0);
	FREEMPI(Temp);

	return Rptr;
}

MPI *INT(MPI *Aptr, MPI *Bptr)
/*
* input: Aptr, Bptr, pointers to MPI'S, *Bptr > 0,
* output: MPI *Qptr = int(*Aptr, *Bptr).
*/
{
	MPI *Q, *R, *Qptr;

	if (Bptr->D == 0)
		return	INT_(Aptr, Bptr->V[0]);
	if (Aptr->S >= 0)
		return INT0(Aptr, Bptr);
	else
	{
		Aptr->S = 1;
		Q = INT0(Aptr, Bptr);
		R = MOD0(Aptr, Bptr);
		Aptr->S = -1;
		/* -*Aptr = *Q * *Bptr + *R, 0 <= *R < *Bptr,
		*Aptr = (-*Q - 1) * *Bptr + *Bptr - R */
		if (R->S == 1)
		{/* add 1 to |*Q|, but not changing the sign of *Q */
			Qptr = ADD0_I(Q, (USL)1);
			Qptr->S = -1;
		}
		else
		{					/* *R = 0 here */
			Q->S = -1;
			Qptr = COPYI(Q);
		}
		FREEMPI(Q);
		FREEMPI(R);
		return Qptr;
	}
}

MPI *INTI(MPI *Aptr, MPI *Bptr)
/*
* input: Aptr, Bptr, pointers to MPI'S, *Bptr != 0,
* output: MPI *Qptr = int(*Aptr, *Bptr).
*/
{
	MPI *A, *B, *Qptr;

	if (Bptr->S > 0)
		return INT(Aptr, Bptr);
	else
	{
		A = COPYI(Aptr);
		B = COPYI(Bptr);
		A->S = -(A->S);
		B->S = -(B->S);
		Qptr = INT(A, B);
		FREEMPI(A);
		FREEMPI(B);
		return Qptr;
	}
}

MPI *MOD(MPI *Aptr, MPI *Bptr)
/*
* input: Aptr, Bptr, pointers to MPI'S, *Bptr > 0.
* output: MPI *Rptr = *Aptr mod *Bptr.
*/
{
	MPI *R, *Rptr;

	if (Bptr->D == 0)
		return CHANGE(MOD_(Aptr, Bptr->V[0]));
	if (Aptr->S >= 0)
		return 	MOD0(Aptr, Bptr);
	else
	{
		Aptr->S = 1;
		R = MOD0(Aptr, Bptr);
		Aptr->S = -1;
		/* -*Aptr = *Q * *Bptr + *R, 0 <= *R < *Bptr,
		*Aptr = (-*Q - 1) * *Bptr + *Bptr - *R */
		if (R->S == 1)
			Rptr = SUB0I(Bptr, R);
		else
			Rptr = ZEROI();
		FREEMPI(R);
		return Rptr;
	}
}

unsigned long MODINT0_(MPI *Aptr, unsigned long b, MPI **Qptr)
/*
* input: Aptr is a pointer to an non-negative MPI, b is a positive integer,
* b < R0. output: *Aptr (mod b) and *Qptr = int(*Aptr / b).
*/
{
	unsigned int k, n;
	unsigned long t;
	MPI *tmp;

	if (b >= R0)
	{
		fprintf(stderr, "in MODINT0_, %lu >= R0\n", b);
		exit(1);
	}
	*Qptr = COPYI(Aptr);
	n = (*Qptr)->D;
	t = (*Qptr)->V[(*Qptr)->D];
	if ((*Qptr)->D != 0)
	{
		for (k = (*Qptr)->D; k >= 1; k--)
		{
			(*Qptr)->V[k] = t / b;
			t = (*Qptr)->V[k - 1] + ((t % b) << T0);
		}
		if ((*Qptr)->V[(*Qptr)->D] == 0)
		{

			tmp = BANK_REALLOC(*Qptr, n - 1);
			FREEMPI(*Qptr);
			*Qptr = tmp;
			/* we remark that (*Qptr)->S = 1 here. */
		}
	}
	(*Qptr)->V[0] = t / b;
	if ((*Qptr)->D == 0)
		(*Qptr)->S = (*Qptr)->V[0] ? 1 : 0;
	return (t % b);
}

MPI *HALFMOD(MPI *Aptr, MPI *Pptr)
/* Returns R=Aptr(mod Pptr) if R <= Pptr/2, otherwise R-Pptr. */
{
	MPI *T, *U, *W;
	int t;

	T = MOD(Aptr, Pptr);
	W = MULT_II(T, 2);
	t = RSV(W, Pptr);
	FREEMPI(W);
	if (t == 1)
		U = SUBI(T, Pptr);
	else
		U = COPYI(T);
	FREEMPI(T);
	return (U);
}

MPI *HALFMODX(MPI *Aptr, MPI *Pptr)
/* Returns R=Aptr(mod Pptr) if R <= Pptr/2, otherwise R-Pptr. */
{
	if (Pptr->S <= 0) {
		printf("N <= 0\n");
		return(NULL);
	}
	else
		return(HALFMOD(Aptr, Pptr));
}

MPI *CEILINGI(MPI *A, MPI *B)
/* Returns the least integer not less than A/B. */
{
	MPI *X, *TMP, *TMP1;
	USI t;

	X = INTI(A, B);
	TMP = MULTI(X, B);
	t = EQUALI(TMP, A);
	FREEMPI(TMP);
	if (t)
		return (X);
	else {
		TMP = ONEI();
		TMP1 = ADDI(X, TMP);
		FREEMPI(TMP);
		FREEMPI(X);
		return(TMP1);
	}
}

MPI *CEILINGIX(MPI *A, MPI *B)
{
	if (B->S == 0) {
		printf("B = 0\n");
		return(NULL);
	}
	else
		return(CEILINGI(A, B));
}

MPI *NEARINT(MPI *M, MPI *N) {
	/* Returns Z, the nearest integer to M/N,
	* with Z=T if M/N=1/2+T.
	*/
	MPI *TEMP1, *TEMP2, *R;

	R = HALFMOD(M, N);
	TEMP1 = SUBI(M, R);
	TEMP2 = INT(TEMP1, N);
	FREEMPI(R);
	FREEMPI(TEMP1);
	return(TEMP2);
}

MPI *NEARINTX(MPI *M, MPI *N) {
	/* Returns Z, the nearest integer to M/N,
	* with Z=T if M/N=1/2+T.
	*/
	if (N->S <= 0) {
		printf("N <= 0\n");
		return(NULL);
	}
	else
		return(NEARINT(M, N));
}

void POWERD(MPI *A, MPI *B, MPI *D, MPI *N, MPI **AA, MPI **BB) {
	/* *A + *B\sqrt(D)=(A + \sqrt(D))^N  */
	MPI *Y, *TEMP, *TEMP1, *TEMP2, *TEMP3, *TEMP4, *TEMP5;
	MPI *X1, *X2, *T2;

	X1 = COPYI(A);
	X2 = COPYI(B);
	Y = COPYI(N);
	*AA = ONEI();
	*BB = ZEROI();
	while (Y->S) {
		while ((Y->V[0]) % 2 == 0) {
			TEMP = Y;
			Y = INT_(Y, 2);
			FREEMPI(TEMP);
			TEMP = X1;
			TEMP1 = MULTI(X2, X2);
			TEMP2 = MULTI(X1, X1);
			TEMP3 = MULTI(D, TEMP1);
			FREEMPI(TEMP1);
			X1 = ADDI(TEMP2, TEMP3);
			FREEMPI(TEMP2);
			FREEMPI(TEMP3);
			TEMP4 = MULTI(TEMP, X2);
			FREEMPI(TEMP);
			T2 = X2;
			X2 = MULT_I(TEMP4, 2);
			FREEMPI(TEMP4);
			FREEMPI(T2);
		}
		TEMP = Y;
		Y = SUB0_I(Y, 1);
		FREEMPI(TEMP);
		TEMP = *AA;
		TEMP1 = MULTI(*BB, X2);
		TEMP2 = MULTI(*AA, X1);
		TEMP3 = MULTI(D, TEMP1);
		FREEMPI(TEMP1);
		*AA = ADDI(TEMP2, TEMP3);
		FREEMPI(TEMP2);
		FREEMPI(TEMP3);
		TEMP4 = MULTI(TEMP, X2);
		FREEMPI(TEMP);
		TEMP5 = MULTI(*BB, X1);
		T2 = *BB;
		*BB = ADDI(TEMP4, TEMP5);
		FREEMPI(T2);
		FREEMPI(TEMP4);
		FREEMPI(TEMP5);
	}
	FREEMPI(Y);
	FREEMPI(X1);
	FREEMPI(X2);
	return;
}
