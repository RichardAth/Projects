/* i5I.c */
/* programs for inputting and printing MPI's. */

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


/*
This code implements a memory bank for allocated (but finished with) MPI's.
Please see the code below for more comments on what is happening.
Peter Adams, 1/4/1993.

To turn the MPIBANK off, include the following line INSIDE these comments.
To turn the MPIBANK on,  move the line outside the comments.
#define WANTMPIBANK
*/

#define WANTMPIBANK 
/*
Only maintain 4 MPIBANKS with at most 500 MPI's in each bank on a PC, due to
memory restrictions.
Maintain 20 MPIBANKS with no size limit on other machines
*/

#ifdef WANTMPIBANK
#ifdef _WIN32
#define NMPIBANKS 4
#define MAXMPIBSIZE 500
int MPI_BANKSIZE[NMPIBANKS];
#else
#define NMPIBANKS 20
#endif
#else    /* do not want the MPI bank at all */
#define NMPIBANKS 0
#endif

/* the +1 in the array dimension following is just to allow for the case when
NMPIBANKS is 0. We do not use the array in that case, but to simplify location
of the array definition inside #define statements, we just always define of
array of size at least 1 */

MPI *MPI_BANK[NMPIBANKS + 1];

/*
The MPI Bank works like this:
MPI's are of different sizes. Hence it does not work to place all free'd MPI's
into the one bank. Now (in theory) MPI's can be of ANY size. However, in most
cases they are small. Hence much malloc'ing and free'ing was being done on
MPI's of the same small sizes.

We maintain n=NMPIBANKS of free'd MPI's. Each of these banks is just a linked
list of MPI's of that size.
*/

void INITMPIBANK() {
	/*
	Initially, set the banks to being empty (a null pointer). On a PC, because we
	restrict the length allowed for each list, set a counter to 0.
	*/
	int i;

	for (i = 0; i < NMPIBANKS; i++)
	{
		MPI_BANK[i] = NULL;
#ifdef _WIN32
		MPI_BANKSIZE[i] = 0;
#endif
	}
}

/*
* mallocs space for an MPI of length n.
* If there is an MPI in the bank of this size, then use it rather than malloc.
*/
MPI *BUILDMPI(unsigned int n) {
	MPI *Mptr;
	unsigned int m;

	m = n - 1;
	if (m < NMPIBANKS && MPI_BANK[m])
	{
		Mptr = MPI_BANK[m];
		MPI_BANK[m] = MPI_BANK[m]->NEXT;
#ifdef _WIN32
		MPI_BANKSIZE[m]--;
#endif
	}
	else
	{
		Mptr = (MPI *)mmalloc(sizeof(MPI));
		Mptr->V = (unsigned long *)mmalloc(n * sizeof(unsigned long));
		Mptr->D = n - 1;
	}
	return (Mptr);
}

void FREEMPIBANK() {
	/*
	This routine should be called to clear out the MPI bank, usually at the end of
	execution. If you like, you can call it during execution, and all of the MPI's
	in the list will be free'd up.
	The only time you will probably want to call this during execution is if the
	algorithm just performed a lot of very memory intensive operations, and so there
	will be many many MPI's on the list. Freeing them up will free up some memory,
	but it will mean that subsequent calls to create MPI's will need to perform the
	(slow) malloc command until there are some more free ones.
	*/
	int i;
	MPI *tmp, *this;
	unsigned int m;

	for (i = 0; i<NMPIBANKS; i++)
	{
		this = MPI_BANK[i];
		while (this)
		{
			tmp = this->NEXT;
			m = this->D;
			ffree((char *)(this->V), (1 + m) * sizeof(unsigned long));
			ffree((char *)this, sizeof(MPI));
			this = tmp;
		}
		MPI_BANK[i] = NULL;
	}
	return;
}

void FREEMPI(MPI *Mptr) {
	/*
	* deallocates the space alloted for the MPI Mptr.
	* This involves placing it onto the appropriate MPI Bank, if there is such a
	* bank, and (on a PC) if there is room in the bank.
	* If it doesn't go into a bank, then it is free'd.
	*/
	unsigned int m;

	if (Mptr == NULL)
		return;
	m = Mptr->D;

	if (m < NMPIBANKS)
	{
#ifdef _WIN32
		if (MPI_BANKSIZE[m] == MAXMPIBSIZE)
		{
			ffree(Mptr->V, (1 + m) * sizeof(unsigned long));
			ffree(Mptr, sizeof(MPI));
		}
		else
		{
			Mptr->NEXT = MPI_BANK[m];
			MPI_BANK[m] = Mptr;
			MPI_BANKSIZE[m]++;
		}
#else
		Mptr->NEXT = MPI_BANK[m];
		MPI_BANK[m] = Mptr;
#endif
	}
	else
	{
		ffree(Mptr->V, (1 + m) * sizeof(unsigned long));
		ffree(Mptr, sizeof(MPI));
	}
}

MPI *COPYI(MPI *Aptr)
/*
* returns a pointer to a copy of *Aptr.
*/
{
	MPI *Bptr;
	if (Aptr == NULL)
		return NULL;
	Bptr = BUILDMPI(1 + Aptr->D);
	Bptr->S = Aptr->S;
	memcpy(Bptr->V, Aptr->V, (1 + Aptr->D) * sizeof(unsigned long));
	return Bptr;
}

void FPRINTI(FILE *outfile, MPI *Mptr) {
	/*
	* the MPI *Mptr is printed in decimal notation to *outfile.
	* no new-line is incorporated.
	*/
	MPI *F, *Temp;
	unsigned long *A, l;
	int i = 0;

	if (Mptr->S == 0)
	{
		fprintf(outfile, "%lu", (USL)0);
		return;
	}
	if (Mptr->S == -1)
		fprintf(outfile, "-");
	F = ABSI(Mptr);
	l = LENGTHI(F);
	A = (unsigned long *)mmalloc((USL)(l * sizeof(unsigned long)));
	do
	{
		A[i] = MOD0_(F, (USL)10);
		Temp = F;
		F = INT0_(F, (USL)10);
		FREEMPI(Temp);
		i++;
	} while (F->S);
	i--;
	FREEMPI(F);
	while (i >= 0)
	{
		fprintf(outfile, "%lu", A[i]);
		i--;
	}
	ffree((char *)A, l * sizeof(unsigned long));
	return;
}

void PRINTI(MPI *Mptr) {
	if (Mptr == NULL)
		printf("NULL");
	else
		FPRINTI(stdout, Mptr);
	return;

}

MPI *ZEROI() {
	/*
	* Returns the MPI ZERO.
	*/
	MPI *Eptr;

	Eptr = BUILDMPI(1);
	Eptr->S = 0;
	Eptr->V[0] = 0;
	return Eptr;
}

MPI *ONEI() {
	/*
	* Returns the MPI ONE.
	*/
	MPI *Eptr;

	Eptr = BUILDMPI(1);
	Eptr->S = 1;
	Eptr->V[0] = 1;
	return Eptr;
}

MPI *TWOI() {
	/*
	* Returns the MPI TWO.
	*/
	MPI *Eptr;

	Eptr = BUILDMPI(1);
	Eptr->S = 1;
	Eptr->V[0] = 2;
	return Eptr;
}

MPI *THREEI() {
	/*
	* Returns the MPI THREE.
	*/
	MPI *Eptr;

	Eptr = BUILDMPI(1);
	Eptr->S = 1;
	Eptr->V[0] = 3;
	return Eptr;
}

MPI *MINUS_ONEI() {
	/*
	* Returns the MPI MINUS_ONE.
	*/
	MPI *Eptr;

	Eptr = BUILDMPI(1);
	Eptr->S = -1;
	Eptr->V[0] = 1;
	return Eptr;
}

MPI *CHANGE(unsigned long n) {
	/*
	* converts n, 0 <= n < R0*R0 to an MPI *Mptr.
	*/
	MPI *Mptr;

	if (n == 0)
		Mptr = ZEROI();
	else if (n == 1)
		Mptr = ONEI();
	else if (n < R0)
	{
		Mptr = BUILDMPI(1);
		Mptr->S = 1;
		Mptr->V[0] = n;
	}
	else
	{
		Mptr = BUILDMPI(2);
		Mptr->S = 1;
		Mptr->V[0] = n % R0;
		Mptr->V[1] = n / R0;
	}
	return Mptr;
}

MPI *CHANGEL(long n) {
	/*
	* converts n, -R0*R0/2 <= n < R0*R0/2 to an MPI *Mptr.
	*/
	MPI *Mptr;
	unsigned long m;

	if (n >= 0)
		return (CHANGE((USL)n));
	else
	{
		m = (USL)(-n);
		Mptr = CHANGE(m);
		Mptr->S = -1;
		return Mptr;
	}
}

MPI *CHANGEI(long n) {
	/*
	* converts n, |n| < R0 to an MPI *Mptr.
	*/
	MPI *Mptr;
	long m;

	Mptr = BUILDMPI(1);
	m = n;
	if (m<0) {
		m = -n;
	}
	Mptr->V[0] = (unsigned long)m;
	if (n == 0) {
		Mptr->S = 0;
	}
	else if (n < 0) {
		Mptr->S = -1;
	}
	else {
		Mptr->S = 1;
	}
	return Mptr;
}

MPI *FINPUTI(FILE *f, unsigned int *uptr) {
	/*
	* Converts decimal input from stream into an MPI.
	* Ignores the combination '\' followed by '\n'.
	* If a rubbish character is met before decimal input, Mptr is set to 0
	* and 0 is returned. All characters up to and including the first newline
	* met are wiped.
	* If a rubbish character is met immediately after decimal input, 0 is
	* returned and all characters up to and including the first newline met
	* are wiped. Otherwise 1 is returned.
	* In any case *Mptr is set equal to any inputted decimal.
	*/
	unsigned long n;
	int s, t;
	MPI *Mptr, *Temp;

	Mptr = ZEROI();
	s = whitespace(f);
	if ((s != '-') && (s < '0' || s > '9'))
	{
		printf("illegal character1 %c entered:\n", s);
		FFlush(f);
		*uptr = 0;
		return Mptr;
	}
	else if (s == '0')
	{
		s = fgetc(f);
		if (s != ' ' && s != '\n' && s != 'i' && s != '/' && s != '+' && s != '-')
		{
			printf("illegal character2 %c entered:\n", s);
			FFlush(f);
			*uptr = 0;
			return Mptr;
		}
		else
		{
			ungetc(s, f);
			*uptr = 1;
			return Mptr;
		}
	}
	else if (s == '-')
	{
		do
			s = fgetc(f);
		while (s == ' ');
		if (s <= '0' || s > '9')
		{
			printf("illegal character3 %c entered:\n", s);
			FFlush(f);
			*uptr = 0;
			return Mptr;
		}
		t = -1;
	}
	else
		t = 1;
	while ((s == '\\') || (s >= '0' && s <= '9'))
	{
		if (s == '\\')
			/* for use with bc output, where '\' is followed by '\n' */
			s = fgetc(f);
		else
		{
			Temp = Mptr;
			Mptr = MULT_I(Temp, 10L);
			FREEMPI(Temp);
			n = (unsigned long)(s - '0');
			Temp = Mptr;
			Mptr = ADD0_I(Temp, n);
			FREEMPI(Temp);
		}
		s = fgetc(f);
	}
	if (t == -1)
		Mptr->S = -1;
	if (s != ' ' && s != '\n' && s != 'i' && s != '/' && s != '+' && s != '-')
	{
		printf("illegal character4 %c entered:\n", s);
		FFlush(f);
		*uptr = 0;
		return Mptr;
	}
	ungetc(s, f);
	*uptr = 1;
	return Mptr;
}

MPI *INPUTI(unsigned int *uptr)
{
	return  FINPUTI(stdin, uptr);
}

unsigned int EQUALI(MPI *Aptr, MPI *Bptr) {
	/*
	* EQUALI(Aptr, Bptr) = 1 if *Aptr = *Bptr, otherwise  = 0.
	*/
	int j;

	if (Aptr->S != Bptr->S)
		return (0);
	else if (Aptr->D != Bptr->D)
		return (0);
	else
	{
		j = Aptr->D;
		while ((Aptr->V[j] == Bptr->V[j]) && (j > 0))
			j--;
		if (Aptr->V[j] == Bptr->V[j])
			return (1);
		else
			return (0);
	}
}

MPI *ABSI(MPI *Aptr) {
	/*
	* *Bptr = |*Aptr|.
	*/
	MPI *Bptr;

	Bptr = COPYI(Aptr);
	Bptr->S = (Aptr->S >= 0 ? Aptr->S : -Aptr->S);
	return Bptr;
}

MPI *MINUSI(MPI *Aptr) {
	/*
	* *Bptr = -(*Aptr).
	*/
	MPI *Bptr;

	Bptr = COPYI(Aptr);
	Bptr->S = -(Aptr->S);
	return Bptr;
}

MPI *GCD(MPI *Aptr, MPI *Bptr) {
	/*
	* GCD(|*Aptr|, |*Bptr|).
	*/
	MPI *A, *B, *Cptr, *Temp;

	A = ABSI(Aptr);
	if (Bptr->S == 0)
		return A;
	B = ABSI(Bptr);
	Cptr = MOD0(A, B);
	while (Cptr->S == 1)
	{
		Temp = A;
		A = B;
		FREEMPI(Temp);
		B = Cptr;
		Cptr = MOD0(A, B);
	}
	Temp = Cptr;
	Cptr = B;
	FREEMPI(Temp);
	FREEMPI(A);
	return Cptr;
}

MPI *LCM(MPI *Aptr, MPI *Bptr) {
	/*
	* *Cptr = LCM(*Aptr, *Bptr).
	*/
	MPI *C, *N, *Cptr, *Temp;

	if (Aptr->S == 0 || Bptr->S == 0)
		return ZEROI();
	else
	{
		N = GCD(Aptr, Bptr);
		C = MULTI(Aptr, Bptr);
		Temp = C;
		C = ABSI(Temp);
		FREEMPI(Temp);
		Cptr = INT0(C, N);
		FREEMPI(C);
		FREEMPI(N);
		return Cptr;
	}
}

MPI *LCM_ARRAY(MPIA M) {
	/*
	* *Aptr = LCM(M[0], ..., M[n - 1]).
	*/
	MPI *L, *Aptr;
	unsigned int i;
	USL n = M->size;
	if (n == 0)
		return ZEROI();

	Aptr = COPYI(M->A[0]);
	if (n == 1)
		return Aptr;
	else
	{
		for (i = 1; i < n; i++)
		{
			L = Aptr;
			Aptr = LCM(L, M->A[i]);
			FREEMPI(L);
		}
		return Aptr;
	}
}

void CONT_FRAC(MPI *Aptr, MPI *Bptr, MPI *Zptr, MPI **Xptr, MPI **Yptr) {
	/*
	* Euclid's algorithm for use in 2SQUARES.
	*/
	MPI *A, *B, *C, *Temp;

	A = COPYI(Aptr);
	B = COPYI(Bptr);
	C = MOD0(A, B);
	while (C->S)
	{
		if (RSV(C, Zptr) <= 0)
		{
			*Xptr = COPYI(C);
			Temp = A;
			A = B;
			FREEMPI(Temp);
			B = C;
			*Yptr = MOD0(A, B);
			FREEMPI(A);
			FREEMPI(B);
			return;
		}
		Temp = A;
		A = B;
		FREEMPI(Temp);
		B = C;
		C = MOD0(A, B);
	}
}

void INTSETUI(unsigned int *Aptr, unsigned long n, unsigned int m) {
	/*
	* this converts the n array elements *Aptr,....*(Aptr + n - 1) to m.
	*/
	for (; n; n--)
		*(Aptr++) = m;
	return;
}

void INTSETUL(unsigned long *Aptr, unsigned long n, unsigned long m) {
	/*
	* this converts the n array elements *Aptr,....*(Aptr + n - 1) to m.
	*/
	for (; n; n--)
		*(Aptr++) = m;
	return;
}

void INTSETI(int *Aptr, unsigned long n, int m) {
	/*
	* this converts the n array elements *Aptr,....*(Aptr + n - 1) to m.
	*/
	for (; n; n--)
		*(Aptr++) = m;
	return;
}

void INTSETL(long *Aptr, unsigned long n, long m)
/*
* this converts the n array elements *Aptr,....*(Aptr + n - 1) to m.
*/
{
	for (; n; n--)
		*(Aptr++) = m;
	return;
}

unsigned long LENGTHI(MPI *Mptr)
/*
* returns the number of decimal digits in the MPI *Mptr,
* increased by 1 if *Mptr is negative.
*/
{
	MPI *F, *Temp;
	unsigned long i = 0, s;

	if (Mptr->S >= 0)
		s = 0;
	else
		s = 1;
	F = ABSI(Mptr);
	do
	{
		Temp = F;
		F = INT0_(Temp, (USL)10);
		FREEMPI(Temp);
		i++;
	} while (F->S == 1);
	FREEMPI(F);
	return (i + s);
}

unsigned long CONVERTI(MPI *N)
/*
* Returns |N| as an unsigned long, providing |N| < 2^32.
* If N's too big, you get an overflow error.
*/
{
	if (N->D > 1)
	{
		fprintf(stderr, "overflow in CONVERTI\n");
		exit(1);
	}
	if (N->S == 0)
		return (USL)0;
	if (N->D == 0)
		return N->V[0];
	else
		return ((N->V[1]) * R0 + (N->V[0]));
}

long CONVERTII(MPI *N)
/*
* Returns N as a long, providing |N| < 2^32.
* If N's too big, you get an overflow error.
*/
{
	int t;
	long int m;
	if (N->D > 1)
	{
		fprintf(stderr, "overflow in CONVERTI\n");
		exit(1);
	}
	t = N->S;
	m = (long)CONVERTI(N);
	if (t >= 0) {
		return (m);
	}
	else {
		return(-m);
	}
}

MPI *INPUTSI(char **ptr, unsigned int *uptr)
/*
* Converts decimal input from a stream into an MPI, leaving *ptr to point
* to the character following the last decimal digit.
* Ignores the combination '\' followed by '\n'.
* If a rubbish character is met before decimal input, Mptr is set to 0
* and *uptr = 0 is returned.
* If a rubbish character is met immediately after decimal input, 0 is
* returned,  Otherwise *uptr = 1 is returned.
* In any case *Mptr is set equal to any inputted decimal.
*/
{
	unsigned long n;
	int t;
	MPI *Mptr, *Temp;

	Mptr = ZEROI();
	while (**ptr == ' ' || **ptr == '\n')
		(*ptr)++;
	if (**ptr != '-' && (**ptr < '0' || **ptr > '9'))
	{
		printf("in INPUTSI1, illegal character %c entered:\n", **ptr);
		*uptr = 0;
		return (Mptr);
	}
	else if (**ptr == '0')
	{
		(*ptr)++;
		if (**ptr >= 0 || **ptr <= 9)
			t = 1;
		else if (**ptr != ' ' && **ptr != '\n' && **ptr != 'i' && **ptr != '/' && **ptr != '+' && **ptr != '-' && **ptr != '\0')
		{
			printf("in INPUTSI2, illegal character %c entered:\n", **ptr);
			*uptr = 0;
			return (Mptr);
		}
		else
		{
			*uptr = 1;
			return (Mptr);
		}
	}
	else if (**ptr == '-')
	{
		do
			(*ptr)++;
		while (**ptr == ' ');
		if (**ptr <= '0' || **ptr > '9')
		{
			printf("in INPUTSI3 illegal character %c entered:\n", **ptr);
			*uptr = 0;
			return (Mptr);
		}
		t = -1;
		while (**ptr == ' ')
			(*ptr)++;
	}
	else
		t = 1;
	while ((**ptr == '\\') || (**ptr >= '0' && **ptr <= '9'))
	{
		if (**ptr == '\\')
			/* for use with bc output, where '\' is followed by '\n'*/
			(*ptr)++;
		else
		{
			Temp = Mptr;
			Mptr = MULT_I(Temp, 10L);
			FREEMPI(Temp);
			n = (unsigned long)(**ptr - '0');
			Temp = Mptr;
			Mptr = ADD0_I(Temp, n);
			FREEMPI(Temp);
		}
		(*ptr)++;
	}
	if (t == -1)
		Mptr->S = -1;
	if (**ptr != '.' && **ptr != ' ' && **ptr != '\n' && **ptr != 'i' && **ptr != '/' && **ptr != '+' && **ptr != '-' && **ptr != '\0')
	{
		printf("in INPUTSI4 illegal character %c entered:\n", **ptr);
		*uptr = 0;
		return (Mptr);
	}
	*uptr = 1;
	return (Mptr);
}

void SPRINTI(char *buffer, MPI *Mptr)
/*
* The decimal digits of the MPI *Mptr are placed in the string buffer.
* no new-line is incorporated.
*/
{
	MPI *F, *Temp;
	unsigned long *A, l;
	int i = 0;

	if (Mptr->S == 0)
	{
		sprintf(buffer, "%lu", (USL)0);
		return;
	}
	if (Mptr->S == -1)
	{
		sprintf(buffer, "-");
		buffer++;
	}
	F = ABSI(Mptr);
	l = LENGTHI(F);
	A = (unsigned long *)mmalloc((USL)(l * sizeof(unsigned long)));
	do
	{
		A[i] = MOD0_(F, (USL)10);
		Temp = F;
		F = INT0_(Temp, (USL)10);
		FREEMPI(Temp);
		i++;
	} while (F->S);
	i--;
	while (i >= 0)
	{
		sprintf(buffer, "%lu", A[i]);
		buffer++;
		i--;
	}
	ffree((char *)A, l * sizeof(unsigned long));
	FREEMPI(F);
	return;
}

MPI *NEAREST_INTI(MPI *Aptr, MPI *Bptr)
/*
* Y is the nearest integer to *Aptr/(*Bptr), *Bptr > 0.
*/
{
	MPI *X, *Y, *Z, *Tmp1, *Tmp2;

	Y = INTI(Aptr, Bptr);
	X = MULTI(Bptr, Y);
	Z = SUBI(Aptr, X);
	FREEMPI(X);
	Tmp1 = Z;
	Z = MULT_I(Z, 2L);
	FREEMPI(Tmp1);
	if (RSV(Z, Bptr) == 1)
	{
		Tmp2 = ONEI();
		Tmp1 = Y;
		Y = ADDI(Y, Tmp2);
		FREEMPI(Tmp2);
		FREEMPI(Tmp1);
	}
	FREEMPI(Z);


	/*
	printf("Aptr = ");PRINTI(Aptr);printf("\n");
	printf("Bptr = ");PRINTI(Bptr);printf("\n");
	printf("Y = ");PRINTI(Y);printf("\n");

	*/

	return (Y);
}

MPI *GCD_ARRAY(MPIA M)
/*
* *Aptr = GCD(M[0], ..., M[n - 1]).
*/
{
	MPI *L, *Aptr;
	unsigned int i;
	unsigned n = M->size;
	Aptr = COPYI(M->A[0]);
	if (n == 1)
		return Aptr;
	else
	{
		for (i = 1; i < n; i++)
		{
			L = Aptr;
			Aptr = GCD(L, M->A[i]);
			FREEMPI(L);
		}
		return Aptr;
	}
}

void PRINTIA(MPIA Mptr)
/*
* Prints an array Mptr[0],...,Mptr[n], n=Aptr, consisting of *Aptr MPI's.
*/
{

	USI i;
	for (i = 0; i < Mptr->size; i++)
	{
		printf("[%u]:", i);
		PRINTI(Mptr->A[i]);
		printf("\n");
	}
	return;
}

/*void EUCLID(MPI *Aptr, MPI *Bptr, MPI **C[], MPI **Dptr)*/
/* Returns the partial quotients C[0],...C[n] of the continued fraction
* expansion of Aptr/Bptr. Here *Dptr = n.
*/
/*{
MPI *A, *B, *Cptr, *Temp;
USL n, i;

n = 2 * BINARYB(Bptr) + 1;
*C = (MPI **)mmalloc((USL)(n * sizeof(MPI *)));
A = COPYI(Aptr);
B = COPYI(Bptr);
(*C)[0] = INT0(A, B);
Cptr = MOD0(A, B);
i = 1;
while (Cptr->S == 1)
{
Temp = A;
A = B;
FREEMPI(Temp);
B = Cptr;
(*C)[i] = INT0(A, B);
Cptr = MOD0(A, B);
i++;
}
*C = (MPI **)rrealloc(*C, (USL)(i * sizeof(MPI *)), -((long)(n - i) * sizeof(MPI *)));
*Dptr = CHANGE(i - 1);
FREEMPI(A);
FREEMPI(B);
FREEMPI(Cptr);
return;
}
*/

void CONVERGENTS(MPIA A, MPIA *P, MPIA *Q)
/*
* P(0)=A[0],P[1)=A[0]*A[1]+1,P[i+1]=A[i+1]*P[i]+P[i-1] if i >= 1.
* Q(0)=1,Q[1)=A[1],Q[i+1]=A[i+1]*Q[i]+Q[i-1] if i >= 1.
* The arrays P[0],...,P[N] and Q[0],..., Q[N] are returned.
*/
{
	USL i, n;
	MPI *X1, *X2, *Y1, *Y2, *Z1, *Z2, *temp = NULL;
	MPI *ONE;
	char buff[20];
	FILE *outfile;
	strcpy(buff, "convergents.out");
	if (access(buff, R_OK) == 0)
		unlink(buff);
	outfile = fopen(buff, "w");

	ONE = ONEI();
	n = A->size - 1;
	Z1 = (MPI *)NULL;
	Z2 = (MPI *)NULL;

	*P = BUILDMPIA();
	*Q = BUILDMPIA();
	ADD_TO_MPIA(*P, A->A[0], 0);
	ADD_TO_MPIA(*Q, ONE, 0);

	if (n == 0)
	{
		fprintf(outfile, "P[0]=");
		FPRINTI(outfile, (*P)->A[0]);
		fprintf(outfile, ";Q[0]=");
		FPRINTI(outfile, (*Q)->A[0]);
		fprintf(outfile, "\n");
		fclose(outfile);
		FREEMPI(ONE);
		return;
	}
	temp = MULTABC(ONE, A->A[0], A->A[1]);
	ADD_TO_MPIA(*P, temp, 1);
	FREEMPI(temp);
	ADD_TO_MPIA(*Q, A->A[1], 1);
	if (n == 1)
	{
		for (i = 0; i <= n; i++)
		{
			fprintf(outfile, "P[%lu]=", i);
			FPRINTI(outfile, (*P)->A[i]);
			fprintf(outfile, ";Q[%lu]=", i);
			FPRINTI(outfile, (*Q)->A[i]);
			fprintf(outfile, "\n");
		}
		fclose(outfile);
		FREEMPI(ONE);
		return;
	}
	X1 = COPYI(A->A[0]);
	Y1 = MULTABC(ONE, A->A[0], A->A[1]);
	X2 = ONEI();
	Y2 = COPYI(A->A[1]);
	for (i = 2; i <= n; i++)
	{
		Z1 = MULTABC(X1, A->A[i], Y1);
		Z2 = MULTABC(X2, A->A[i], Y2);
		ADD_TO_MPIA(*P, Z1, i);
		ADD_TO_MPIA(*Q, Z2, i);
		FREEMPI(X1);
		FREEMPI(X2);
		X1 = Y1;
		Y1 = Z1;
		X2 = Y2;
		Y2 = Z2;
	}
	FREEMPI(X1);
	FREEMPI(X2);
	FREEMPI(Z1);
	FREEMPI(Z2);
	for (i = 0; i <= n; i++)
	{
		fprintf(outfile, "P[%lu]=", i);
		FPRINTI(outfile, (*P)->A[i]);
		fprintf(outfile, ";Q[%lu]=", i);
		FPRINTI(outfile, (*Q)->A[i]);
		fprintf(outfile, "\n");
	}
	fclose(outfile);
	FREEMPI(ONE);
	return;
}

MPI *SUMI(MPI *A[], USL n)
/* Returns A[0]+...+A[n] */
{
	MPI *S, *temp;
	USL i;

	S = ZEROI();
	for (i = 0; i <= n; i++)
	{
		temp = S;
		S = ADDI(S, A[i]);
		FREEMPI(temp);
	}
	return (S);
}

void LAGRANGE(POLYI P, MPIA *AA, MPI *M)
/*
* P(x)=A[n]x^n+...+A[0], A[n]>0, is a polynomial
* with integer coefficients, having no rational roots
* and having exactly one real positive root x, this being > 1.
* The method of Lagrange (1797) is used to find the
* the first m+1 partial quotients AA[0],...,AA[m] of x.
* (See Knuth, Art of computer programming, volume2,
* problem 13, 4.5.3. Also S. Lang and H. Trotter,
* 'Continued fractions for some algebraic numbers'
* J. fuer Math. 255 (1972) 112-134; Addendum 267 (1974) ibid. 219-220.)
* Also page 261, Number Theory with Applications, by R. Kumanduri and
* C. Romero.

*/
{
	int i, j, m, d;
	MPR *HIGH;
	MPR *ONE = ONER();
	MPI *ZERO = ZEROI();
	POLYI Q = COPYPI(P), CURSOR;
	MPI *TMP, *HI, *LO, *ON, *TEMP, *TEMP1, *ROOT, *MID;
	MPIA COEF;

	ON = ONEI();

	m = CONVERTI(M);
	d = DEGREEPI(P);
	/* I should a check in here for a single positive root.  I will get around
	* to it */

	*AA = BUILDMPIA();


	for (i = 0; i <= m; i++) {
		LO = ONEI();
		HIGH = CAUCHY(Q); /* this is a power of 2 */
		HI = COPYI(HIGH->N);
		FREEMPR(HIGH);
		TEMP = ADD0I(LO, ON);
		while (RSV(TEMP, HI)) {
			FREEMPI(TEMP);
			TEMP1 = ADD0I(LO, HI);
			MID = INT0_(TEMP1, 2);
			FREEMPI(TEMP1);
			TEMP1 = VALPI(Q, MID);
			if (TEMP1->S>0) {
				TEMP = HI;
				HI = MID;
				FREEMPI(TEMP);
			}
			else {
				TEMP = LO;
				LO = MID;
				FREEMPI(TEMP);
			}
			FREEMPI(TEMP1);
			TEMP = ADD0I(LO, ON);
		}
		FREEMPI(TEMP);
		ROOT = LO;
		FREEMPI(HI); /* HI=LOW+1 here */
		ADD_TO_MPIA(*AA, ROOT, i);
		P_OF_X_PLUS(&Q, ROOT);
		FREEMPI(ROOT);

		/* This performs transformation p(x) = -x^d*p(1/x +a) */
		COEF = BUILDMPIA();
		for (j = 0; j <= d; j++)
			ADD_TO_MPIA(COEF, ZERO, j);

		CURSOR = Q;
		while (CURSOR) {
			TMP = MINUSI(CURSOR->COEF);
			ADD_TO_MPIA(COEF, TMP, DEGREEPI(CURSOR));
			FREEMPI(TMP);
			CURSOR = CURSOR->NEXT;
		}

		DELETEPI(Q);
		Q = NULL;
		for (j = 0; j <= d; j++) {
			PINSERTPI(d - j, COEF->A[j], &Q, 0);
		}
		PURGEPI(&Q);
		FREEMPIA(COEF);
		/* end of transformation */

		if (DEGREEPI(Q)<d) {/* noticed by KRM mid Sept 2002 */
			printf("rational root\n");
			break;
		}
	}
	DELETEPI(Q);
	FREEMPI(ON);
	FREEMPR(ONE);
	FREEMPI(ZERO);
}

void EUCLID(MPI *Aptr, MPI *Bptr, MPIA *Q, MPIA *R, MPIA *S, MPIA *T, MPI **Dptr)
/*
* Returns Q[0]=NULL,Q[1],...Q[n],Q[n+1]=NULL,
* R[0],...R[n + 1],
* S[0],...S[n + 1], T[0],...T[n + 1]
* for Euclid's algorithm on R[0]=Aptr, R[1]=Bptr.
* R[0]=R[1]*Q[1]+R[2]
* R[1]=R[1]*Q[2]+R[3]
* .....
* R[n-2]=R[n-1]*Q[n-1]+R[n]
* R[n-1]=R[n]*Q[n], R[n+1]=0.
* S[0]=1,S[1]=0, S[j]=s[j-2]-Q[j-1]*S[j-1],
* T[0]=0,T[1]=1, T[j]=T[j-2]-Q[j-1]*T[j-1], j=2,...,n+1
* Here *Dptr = n.
*/
{
	MPI *QQ;
	USL n, i, j;
	char buff[20];
	FILE *outfile;
	MPI *ONE = ONEI();
	MPI *ZERO = ZEROI();
	MPI *TMPI;

	strcpy(buff, "euclid.out");
	if (access(buff, R_OK) == 0)
		unlink(buff);
	outfile = fopen(buff, "a");
	n = 2 * BINARYB(Bptr) + 1;

	(*Q) = BUILDMPIA();
	(*R) = BUILDMPIA();
	(*S) = BUILDMPIA();
	(*T) = BUILDMPIA();

	ADD_TO_MPIA((*R), Aptr, 0);
	ADD_TO_MPIA((*R), Bptr, 1);
	ADD_TO_MPIA((*S), ONE, 0);
	ADD_TO_MPIA((*S), ZERO, 1);

	ADD_TO_MPIA((*T), ZERO, 0);
	ADD_TO_MPIA((*T), ONE, 1);
	ADD_TO_MPIA(*Q, NULL, 0);
	TMPI = INT0((*R)->A[0], (*R)->A[1]);
	ADD_TO_MPIA((*Q), TMPI, 1);
	FREEMPI(TMPI);
	TMPI = MOD0((*R)->A[0], (*R)->A[1]);
	ADD_TO_MPIA((*R), TMPI, 2);
	FREEMPI(TMPI);
	i = 2;
	while ((*R)->A[i]->S == 1)
	{

		QQ = MINUSI((*Q)->A[i - 1]);

		TMPI = MULTABC((*S)->A[i - 2], QQ, (*S)->A[i - 1]);
		ADD_TO_MPIA(*S, TMPI, i);
		FREEMPI(TMPI);

		TMPI = MULTABC((*T)->A[i - 2], QQ, (*T)->A[i - 1]);
		ADD_TO_MPIA(*T, TMPI, i);
		FREEMPI(TMPI);
		FREEMPI(QQ);

		TMPI = INT0((*R)->A[i - 1], (*R)->A[i]);
		ADD_TO_MPIA(*Q, TMPI, i);
		FREEMPI(TMPI);

		TMPI = MOD0((*R)->A[i - 1], (*R)->A[i]);
		ADD_TO_MPIA(*R, TMPI, i + 1);
		FREEMPI(TMPI);
		i++;
	}

	QQ = MINUSI((*Q)->A[i - 1]);

	TMPI = MULTABC((*S)->A[i - 2], QQ, (*S)->A[i - 1]);
	ADD_TO_MPIA(*S, TMPI, i);
	FREEMPI(TMPI);
	TMPI = MULTABC((*T)->A[i - 2], QQ, (*T)->A[i - 1]);
	ADD_TO_MPIA(*T, TMPI, i);
	FREEMPI(TMPI);
	FREEMPI(QQ);

	ADD_TO_MPIA(*Q, NULL, i);
	i++;
	*Dptr = CHANGE(i - 1);
	for (j = 0; j < i; j++) {
		if (j == 0 || j == i - 1)
			fprintf(outfile, "Q[%lu]=-", j);
		else
		{
			fprintf(outfile, "Q[%lu]=", j);
			FPRINTI(outfile, (*Q)->A[j]);
		}
		fprintf(outfile, ", R[%lu]=", j); FPRINTI(outfile, (*R)->A[j]);
		fprintf(outfile, ", S[%lu]=", j); FPRINTI(outfile, (*S)->A[j]);
		fprintf(outfile, ", T[%lu]=", j); FPRINTI(outfile, (*T)->A[j]);
		fprintf(outfile, " \n");
	}
	fclose(outfile);
	FREEMPI(ZERO);
	FREEMPI(ONE);
	return;
}


/* Allocates space for an array initially of size 11 (enough to hold a[0] to
* a[10] and sets these slots to contain the zero MPI*/
MPIA BUILDMPIA()
{
	MPIA MA;
	int i;
	if (!(MA = mmalloc(sizeof(struct _MPIA)))) {
		printf("Not enough memory to allocate MPIA\n");
		exit(1);
	}
	MA->A = (MPI **)mmalloc((USL)((MIN_ARRAY_SIZE) * sizeof(MPI *)));
	MA->size = 0;
	MA->slots = MIN_ARRAY_SIZE;
	for (i = 0; i < MIN_ARRAY_SIZE; i++)
		MA->A[i] = ZEROI();
	return MA;
}


/* ADD_TO_MPIA
* Adds the supplied MPI at the subscript n
* If slot already exists MPI at that slot is freed and the new one is added.
* If the n is greater than the number of slots then the array is correctly
* resized and the new MPI is added.  Slots between the previous last slots
* and the new subscript n are initialised to zero.
*/
void ADD_TO_MPIA(MPIA MA, MPI *V, USL n)
{
	int i;

	if (n > MAX_ARRAY_SIZE - 1) {
		printf("Error: Array index must be in range 0 < n < %d\n",
			MAX_ARRAY_SIZE - 1);
		return;
	}

	if (n > MA->slots - 1) {
		MA->A = (MPI **)rrealloc(MA->A, (n + 1) * sizeof(MPI *),
			(n + 1 - (MA->slots)) * sizeof(MPI*));
		for (i = MA->slots; i < n + 1; i++)
			(MA->A)[i] = ZEROI();
		MA->slots = n + 1;
	}
	FREEMPI(MA->A[n]);
	if (V != NULL)
		MA->A[n] = COPYI(V);
	else MA->A[n] = NULL;

	if (n + 1 > MA->size)
		MA->size = n + 1;
}

/* Frees an MPIA previously returned by BUILDMPIA.  Undefined behaviour occurs
* if this is not the case.
* It will free all the slots of the MPI's they contain.
*/
void FREEMPIA(MPIA MA)
{
	int i;
	if (!MA)
		return;
	for (i = 0; i < MA->slots; i++)
		if (MA->A != NULL)
			FREEMPI(MA->A[i]);
	ffree(MA->A, (MA->slots) * sizeof(MPI *));
	ffree(MA, sizeof(struct _MPIA));

}

/* Inserts the term V at position n.  If n > size of array - 1 then this
* function behaves the same as ADD_TO_MPIA
* By insert I mean that the object is placed in at that subscript and the
* rest are shuffled up.
*/
void MPIA_INSERT(MPIA MA, MPI *V, USL n)
{
	/* add an extra slot */
	if (n > MA->size - 1)
		ADD_TO_MPIA(MA, V, n);
	else { /* n < MA->size */
		if (MA->slots == MA->size) { /* then realloc space */
			if (!(MA->A = rrealloc(MA->A, (++MA->slots) * sizeof(MPI *), \
				sizeof(MPI *)))) {
				printf("Not enough memory to increase size of MPIA\n");
				exit(1);
			}
			memmove(MA->A + n + 1, MA->A + n, sizeof(MPI *) * (MA->slots - n - 1));
		}
		else { /* no need to realloc */
			FREEMPI(MA->A[MA->size]);
			memmove(MA->A + n + 1, MA->A + n, sizeof(MPI *) * (MA->size - n));
		}
		MA->size++;
		MA->A[n] = COPYI(V);
	}
}

MPI *EUCLIDI1(MPI *M, MPI *N)
/*
Returns the length of Euclid's algorithm
*/
{
	MPI *I, *A, *B, *R, *TEMP, *TEMP1, *T1, *T2, *T3;
	USL i, k;

	k = 0;
	i = 1;
	A = COPYI(M);
	B = COPYI(N);
	R = MOD0(A, B);
	while (R->S > 0) {
		TEMP = A;
		A = B;
		FREEMPI(TEMP);
		B = R;
		R = MOD0(A, B);
		i++;
		if (i % 1000 == 0) {
			printf("i=%lu\n", i);
		}
		if (i == 2147483647) {
			i = 0;
			k++;
			if (k == 2147483647) {
				printf("exiting prematurely, k = 2147483647\n");
				return(NULL);
			}
		}
	}
	FREEMPI(A);
	FREEMPI(B);
	FREEMPI(R);
	T1 = CHANGE(k);
	T2 = CHANGE(2147483648UL);
	T3 = CHANGE(i);
	TEMP1 = MULTI(T1, T2);
	I = ADD0I(TEMP1, T3); /*I = k*2147483648 + i */
	FREEMPI(T1);
	FREEMPI(T2);
	FREEMPI(T3);
	FREEMPI(TEMP1);
	return (I);
}

MPI *CFRACN(MPI *N) {
	USL n;
	MPI *L, *Z, *ONE, *TWO, *AA, *BB, *TEMP;

	system("date");
	n = CONVERTI(N);
	Z = POWER_I(2, n);
	ONE = ONEI();
	TWO = TWOI();

	POWERD(ONE, ONE, TWO, Z, &AA, &BB);
	TEMP = BB;
	BB = INT0(BB, Z);
	FREEMPI(TEMP);
	L = EUCLIDI1(AA, BB);
	FREEMPI(Z);
	FREEMPI(ONE);
	FREEMPI(TWO);
	FREEMPI(AA);
	FREEMPI(BB);
	if ((L->V[0]) % 2) {
		TEMP = L;
		L = ADD0_I(L, 1);
		FREEMPI(TEMP);
	}
	system("date");
	return(L);
}

MPI *CFRAC_PERIOD(MPI *D)
/*
* This is a program for finding the period of the cfrac of sqrt(D).
* The algorithm is based on K. Rosen, Elementary number theory
* and its applications, p382, B.A. Venkov, Elementary Number theory, p.62
* and D. Knuth, Art of computer programming, Vol.2, p359, with Pohst's trick
* of using half the period.
*/
{
	MPI *P, *Q, *E, *F, *G, *X, *Y, *TEMP, *TEMP1, *H, *PERIOD;
	MPI *T2, *T3;
	USL h, k;
	USI t;

	X = BIG_MTHROOT(D, 2);
	G = MULTI(X, X);

	TEMP = ADD0_I(G, (USL)1);
	FREEMPI(G);
	t = EQUALI(D, TEMP);
	FREEMPI(TEMP);
	if (t) /* period 1, exceptional case */
	{
		FREEMPI(X);
		return(ONEI());
	}
	P = ZEROI();
	Q = ONEI();
	h = 0;
	k = 0;
	while (1)
	{
		TEMP = ADD0I(X, P);
		Y = INT0(TEMP, Q);
		FREEMPI(TEMP);
		F = P;
		TEMP = MULTI(Y, Q);
		FREEMPI(Y);
		P = SUBI(TEMP, P);
		FREEMPI(TEMP);
		t = EQUALI(P, F);
		FREEMPI(F);
		if (t) { /*P_h=P_{h+1}, even period 2h */
				 /*                  printf("P[h]=P[h+1]\n");*/
			T2 = CHANGE(2147483648UL);
			/*T2 = CHANGE(32768);*/
			T3 = CHANGE(h);
			TEMP1 = MULT_I(T2, k);
			H = ADD0I(TEMP1, T3); /*H = k*2147483648 + h */
			FREEMPI(T2);
			FREEMPI(T3);
			FREEMPI(TEMP1);
			PERIOD = MULT_I(H, 2);
			/*printf("H=");
			PRINTI(H);
			printf("\n");*/
			FREEMPI(H);
			FREEMPI(X);
			FREEMPI(P);
			FREEMPI(Q);
			return(PERIOD);
		}
		E = Q;
		TEMP = MULTI(P, P);
		TEMP1 = SUBI(D, TEMP);
		FREEMPI(TEMP);
		Q = INT0(TEMP1, Q);
		t = EQUALI(Q, E);
		FREEMPI(TEMP1);
		FREEMPI(E);
		if (t) { /*Q_h=Q_{h-1}, odd period 2h+1 */
				 /*printf("Q[h]=Q[h-1]\n");*/
			T2 = CHANGE(2147483648UL);
			/*T2 = CHANGE(32768);*/
			T3 = CHANGE(h);
			TEMP1 = MULT_I(T2, k);
			H = ADD0I(TEMP1, T3); /*H = k*2147483648 + h */
			FREEMPI(T2);
			FREEMPI(T3);
			FREEMPI(TEMP1);
			TEMP = MULT_I(H, 2);
			/*printf("H=");
			PRINTI(H);
			printf("\n");*/
			FREEMPI(H);
			PERIOD = ADD0_I(TEMP, 1);
			FREEMPI(TEMP);
			FREEMPI(X);
			FREEMPI(P);
			FREEMPI(Q);
			return(PERIOD);
		}
		h++;
		if (h == 2147483648UL) {
			h = 0;
			k++;
			if (k == 2147483648UL) {
				printf("exiting prematurely, k = 2147483648\n");
				return(NULL);
			}
		}
	}
}


MPI *NSCF_PERIOD(MPI *D, MPI *FLAG)
/*
* This is a program for finding the period of the NSCF of sqrt(D).
*  using the half period method of nscf_midpoint.pdf.
* FLAG=1 prints the type of midpoint - P, Q or PQ.
* FLAG=NULL suppresses printing and is used in CFRAC_COUNT below.
*/
{
	MPI *C, *P, *Q, *G, *X, *TEMP, *TEMP1, *TEMP2, *H;
	MPI *T2, *T3, *OLDP, *OLDQ, *P1, *P2, *Q1, *Q2, *PERIOD;
	USL h, k, flag;
	USI t;
	int a;

	X = BIG_MTHROOT(D, 2);
	G = MULTI(X, X);

	TEMP = ADD0_I(G, (USL)1);
	FREEMPI(G);
	t = EQUALI(D, TEMP);
	FREEMPI(TEMP);
	if (t) /* period 1, exceptional case */
	{
		FREEMPI(X);
		return(ONEI());
	}
	P = ZEROI();
	Q = ONEI();
	h = 0;
	k = 0;
	flag = 0;
	while (1)
	{
		TEMP = ADD0I(X, P);
		C = INT0(TEMP, Q);
		FREEMPI(TEMP);
		TEMP = MULTI(C, Q);
		P1 = SUBI(TEMP, P);
		FREEMPI(TEMP);
		TEMP = MULTI(P1, P1);
		Q1 = SUBI(D, TEMP);
		FREEMPI(TEMP);
		TEMP = Q1;
		Q1 = INT0(Q1, Q);
		FREEMPI(TEMP);
		P2 = ADD0I(P1, Q);
		TEMP = ADD0I(P2, P1);
		Q2 = SUB0I(TEMP, Q1);
		FREEMPI(TEMP);
		OLDP = P;
		OLDQ = Q;
		if (RSV(Q1, Q2) <0) {
			a = 1;
			P = P1;
			Q = Q1;
			FREEMPI(P2);
			FREEMPI(Q2);
		}
		else {
			a = -1;
			P = P2;
			Q = Q2;
			FREEMPI(P1);
			FREEMPI(Q1);
		}
		FREEMPI(C);
		t = EQUALI(P, OLDP);
		FREEMPI(OLDP);
		if (t) { /*P_h=P_{h+1}, even period 2h */
			if (FLAG) {
				printf("P-test, P_h=P_{h+1}\n");
			}
			T2 = CHANGE(2147483648UL);
			/*T2 = CHANGE(32768);*/
			T3 = CHANGE(h);
			TEMP1 = MULT_I(T2, k);
			H = ADD0I(TEMP1, T3); /*H = k*2147483648 + h */
			FREEMPI(T2);
			FREEMPI(T3);
			FREEMPI(TEMP1);
			PERIOD = MULT_I(H, 2);
			/*printf("H=");
			PRINTI(H);
			printf("\n");*/
			FREEMPI(H);
			FREEMPI(OLDQ);
			break;
		}
		t = EQUALI(Q, OLDQ);
		FREEMPI(OLDQ);
		if (t) { /*Q_h=Q_{h+1}, odd period 2h+1 */
			if (FLAG) {
				printf("Q-test:Q_h=Q_{h+1}\n");
			}
			T2 = CHANGE(2147483648UL);
			/*T2 = CHANGE(32768);*/
			T3 = CHANGE(h);
			TEMP1 = MULT_I(T2, k);
			H = ADD0I(TEMP1, T3); /*H = k*2147483648 + h */
			FREEMPI(T2);
			FREEMPI(T3);
			FREEMPI(TEMP1);
			TEMP = MULT_I(H, 2);
			/*printf("H=");
			PRINTI(H);
			printf("\n");*/
			FREEMPI(H);
			PERIOD = ADD0_I(TEMP, 1);
			FREEMPI(TEMP);
			break;
		}
		if (flag) {
			flag = 0;
			T2 = CHANGE(2147483648UL);
			T3 = CHANGE(h);
			TEMP1 = MULT_I(T2, k);
			H = ADD0I(TEMP1, T3); /*H = k*2147483648 + h */
			FREEMPI(T2);
			FREEMPI(T3);
			FREEMPI(TEMP1);
			PERIOD = MULT_I(H, 2);
			/*printf("H=");
			PRINTI(H);
			printf("\n");*/
			FREEMPI(H);
			break;
		}
		TEMP1 = INT0_(OLDQ, (USL)2);
		TEMP2 = ADD0I(Q, TEMP1);
		FREEMPI(TEMP1);
		t = EQUALI(P, TEMP2);
		FREEMPI(TEMP2);
		if (t && a == -1 && flag == 0) {
			if (FLAG) {
				printf("PQ-test\n");
			}
			flag = 1;
		}
		h++;
		if (h == 2147483648UL) {
			h = 0;
			k++;
			if (k == 2147483648UL) {
				printf("exiting prematurely, k = 2147483648\n");
				return(ZEROI());
			}
		}
	}
	FREEMPI(X);
	FREEMPI(P);
	FREEMPI(Q);
	return(PERIOD);
}

void CFRAC_COUNT(MPI *M, MPI *N) {
	MPI *NSCF_SUM, *RCF_SUM, *D, *TEMP, *TEMP1, *TEMP2, *ONE;
	MPR *RATIO;
	USL d, m, n, r;
	int t;

	if (N->D < 2) {
		m = CONVERTI(M);
		n = CONVERTI(N);
		NSCF_SUM = ZEROI();
		RCF_SUM = ZEROI();

		for (d = m; d <= n; d++) {
			D = CHANGE(d);
			t = SQUARETEST(D);
			if (t) {
				TEMP1 = NSCF_PERIOD(D, (MPI *)NULL);
				TEMP2 = CFRAC_PERIOD(D);
				FREEMPI(D);
				TEMP = NSCF_SUM;
				NSCF_SUM = ADD0I(NSCF_SUM, TEMP1);
				FREEMPI(TEMP);
				FREEMPI(TEMP1);
				TEMP = RCF_SUM;
				RCF_SUM = ADD0I(RCF_SUM, TEMP2);
				FREEMPI(TEMP);
				FREEMPI(TEMP2);
			}
			else {
				FREEMPI(D);
			}
			if (d % 10 == 0) {
				printf("%lu:", d);
				PRINTI(NSCF_SUM);
				printf(":");
				PRINTI(RCF_SUM);
				printf(":");
				RATIO = RATIOI(NSCF_SUM, RCF_SUM);
				PRINTDR((USI)7, RATIO);
				printf("\n");
				FREEMPR(RATIO);
			}
		}
		FREEMPI(NSCF_SUM);
		FREEMPI(RCF_SUM);
		return;
	}
	else {
		ONE = ONEI();
		NSCF_SUM = ZEROI();
		RCF_SUM = ZEROI();
		D = COPYI(M);
		while (RSV(D, N) <= 0) {
			t = SQUARETEST(D);
			if (t) {
				TEMP1 = NSCF_PERIOD(D, (MPI *)NULL);
				TEMP2 = CFRAC_PERIOD(D);
				TEMP = NSCF_SUM;
				NSCF_SUM = ADD0I(NSCF_SUM, TEMP1);
				FREEMPI(TEMP);
				FREEMPI(TEMP1);
				TEMP = RCF_SUM;
				RCF_SUM = ADD0I(RCF_SUM, TEMP2);
				FREEMPI(TEMP);
				FREEMPI(TEMP2);
			}
			r = MOD0_(D, (USL)10);
			if (r == 0) {
				PRINTI(D);
				printf(": ");
				PRINTI(NSCF_SUM);
				printf(":");
				PRINTI(RCF_SUM);
				printf(":");
				RATIO = RATIOI(NSCF_SUM, RCF_SUM);
				PRINTDR((USI)7, RATIO);
				printf("\n");
				FREEMPR(RATIO);
			}
			TEMP = D;
			D = ADD0I(D, ONE);
			FREEMPI(TEMP);
		}
		FREEMPI(NSCF_SUM);
		FREEMPI(RCF_SUM);
		FREEMPI(ONE);
		FREEMPI(D);
		return;
	}
}

USI SQUARETEST(MPI *D) {
	MPI *X, *TEMP;
	USI t;

	X = BIG_MTHROOT(D, 2);
	TEMP = MULTI(X, X);
	FREEMPI(X);
	t = EQUALI(TEMP, D);
	FREEMPI(TEMP);
	if (t) {
		return(0);
	}
	else {
		return(1);
	}
}

MPI *NSCF_PERIOD0(MPI *D, MPI *Eptr, MPI **Aptr, MPI **Bptr)
/*
* This is a program for finding the period of the NSCF of sqrt(D).
* using the half period method of nscf_midpoint.pdf.
* It also returns the least solution of Pell's \pm 1 equation.
* If Eptr = 1, the midpoint and its type are printed.
*/
{
	MPI *C, *P, *Q, *G, *X, *TEMP, *TEMP1, *TEMP2, *H;
	MPI *T2, *T3, *OLDP, *OLDQ, *P1, *P2, *Q1, *Q2, *PERIOD;
	MPI *B1, *C1, *B2, *C2, *OLDB1, *OLDC1, *OLDB2, *OLDC2, *B;
	USL h, k, flag;
	USI t;
	int a, olda;

	X = BIG_MTHROOT(D, 2);
	G = MULTI(X, X);

	TEMP = ADD0_I(G, (USL)1);
	FREEMPI(G);
	t = EQUALI(D, TEMP);
	FREEMPI(TEMP);
	if (t) /* period 1, exceptional case */
	{
		if (Eptr->S) {
			printf("d=x^2+1\n");
			PRINT_MPI(X, "x");
			printf("x^2-dy^2=1\n");
			printf("y=1\n");
		}
		*Aptr = X;
		*Bptr = ONEI();
		return(ONEI());
	}
	P = ZEROI();
	Q = ONEI();
	h = 0;
	k = 0;
	flag = 0;
	B1 = ZEROI();
	C1 = ONEI();
	B2 = ONEI();
	C2 = ZEROI();
	a = 1;
	while (1)
	{
		OLDC1 = C1;
		OLDC2 = C2;
		OLDB1 = B1;
		OLDB2 = B2;
		olda = a;
		TEMP = ADD0I(X, P);
		C = INT0(TEMP, Q);
		FREEMPI(TEMP);
		TEMP = MULTI(C, Q);
		P1 = SUBI(TEMP, P);
		FREEMPI(TEMP);
		TEMP = MULTI(P1, P1);
		Q1 = SUBI(D, TEMP);
		FREEMPI(TEMP);
		TEMP = Q1;
		Q1 = INT0(Q1, Q);
		FREEMPI(TEMP);
		P2 = ADD0I(P1, Q);
		TEMP = ADD0I(P2, P1);
		Q2 = SUB0I(TEMP, Q1);
		FREEMPI(TEMP);
		OLDP = P;
		OLDQ = Q;
		if (RSV(Q1, Q2) <0) {
			a = 1;
			P = P1;
			Q = Q1;
			FREEMPI(P2);
			FREEMPI(Q2);
			B = C;
		}
		else {
			a = -1;
			P = P2;
			Q = Q2;
			FREEMPI(P1);
			FREEMPI(Q1);
			B = ADD0_I(C, 1);
			FREEMPI(C);
		}

		if (h == 0 && k == 0) {
			C1 = B;
			C2 = ONEI();
		}
		else {
			C1 = MULTAB_PLUS_MINUS_CI(B, C1, olda, B1);
			C2 = MULTAB_PLUS_MINUS_CI(B, C2, olda, B2);
			FREEMPI(B);
		}
		B1 = OLDC1;
		B2 = OLDC2;

		t = EQUALI(P, OLDP);
		FREEMPI(OLDP);
		if (t) { /*P_h=P_{h+1}, even period 2h */
			T2 = CHANGE(2147483648UL);
			T3 = CHANGE(h);
			TEMP1 = MULT_I(T2, k);
			H = ADD0I(TEMP1, T3); /*H = k*2147483648 + h */
			FREEMPI(T2);
			FREEMPI(T3);
			FREEMPI(TEMP1);
			PERIOD = MULT_I(H, 2);
			if (Eptr->S) {
				printf("P-test, P_h=P_{h+1}\n");
				PRINT_MPI(H, "h");
				printf("x^2-dy^2=1\n");
			}
			FREEMPI(H);
			if (olda == 1) {
				*Aptr = MULTAB_PLUS_CDI(C1, OLDC2, OLDC1, OLDB2);
				TEMP = ADDI(C2, OLDB2);
				*Bptr = MULTI(OLDC2, TEMP);
				FREEMPI(TEMP);
			}
			else {
				*Aptr = MULTAB_MINUS_CDI(C1, OLDC2, OLDC1, OLDB2);
				TEMP = SUBI(C2, OLDB2);
				*Bptr = MULTI(OLDC2, TEMP);
				FREEMPI(TEMP);
			}
			FREEMPI(OLDQ);
			break;
		}
		t = EQUALI(Q, OLDQ);
		FREEMPI(OLDQ);
		if (t) { /*Q_h=Q_{h+1}, odd period 2h+1 */
			T2 = CHANGE(2147483648UL);
			T3 = CHANGE(h);
			TEMP1 = MULT_I(T2, k);
			H = ADD0I(TEMP1, T3); /*H = k*2147483648 + h */
			FREEMPI(T2);
			FREEMPI(T3);
			FREEMPI(TEMP1);
			TEMP = MULT_I(H, 2);
			if (Eptr->S) {
				printf("Q-test:Q_h=Q_{h+1}\n");
				PRINT_MPI(H, "h");
				printf("x^2-dy^2=%d\n", -a);
			}
			FREEMPI(H);
			PERIOD = ADD0_I(TEMP, 1);
			FREEMPI(TEMP);
			if (a == 1) {
				*Aptr = MULTAB_PLUS_CDI(C1, C2, OLDC1, OLDC2);
				*Bptr = MULTAB_PLUS_CDI(C2, C2, OLDC2, OLDC2);
			}
			else {
				*Aptr = MULTAB_MINUS_CDI(C1, C2, OLDC1, OLDC2);
				*Bptr = MULTAB_MINUS_CDI(C2, C2, OLDC2, OLDC2);
			}
			break;
		}
		if (flag) {
			T2 = CHANGE(2147483648UL);
			T3 = CHANGE(h);
			TEMP1 = MULT_I(T2, k);
			H = ADD0I(TEMP1, T3); /*H = k*2147483648 + h */
			FREEMPI(T2);
			FREEMPI(T3);
			FREEMPI(TEMP1);
			PERIOD = MULT_I(H, 2);
			if (Eptr->S) {
				PRINT_MPI(H, "h");
				printf("x^2-dy^2= -1\n");
			}
			FREEMPI(H);
			TEMP = SUBI(OLDC1, OLDB1);
			*Aptr = MULTAB_MINUS_CDI(C1, OLDC2, OLDB2, TEMP);
			FREEMPI(TEMP);
			TEMP = MULT_I(OLDC2, 2);
			*Bptr = MULTAB_MINUS_CDI(TEMP, OLDC2, C2, OLDB2);
			FREEMPI(TEMP);
			break;
		}
		TEMP1 = INT0_(OLDQ, (USL)2);
		TEMP2 = ADD0I(Q, TEMP1);
		FREEMPI(TEMP1);
		t = EQUALI(P, TEMP2);
		FREEMPI(TEMP2);
		if (t && a == -1 && flag == 0) {
			if (Eptr->S) {
				printf("PQ-test, P[h]=Q[h]+Q[h-1]/2\n");
			}
			flag = 1;
		}
		h++;
		if (h == 2147483648UL) {
			h = 0;
			k++;
			if (k == 2147483648UL) {
				printf("exiting prematurely, k = 2147483648\n");
				return(NULL);
			}
		}
		FREEMPI(OLDB1);
		FREEMPI(OLDB2);
	}
	FREEMPI(X);
	FREEMPI(P);
	FREEMPI(Q);
	FREEMPI(B1);
	FREEMPI(B2);
	FREEMPI(C1);
	FREEMPI(C2);
	FREEMPI(OLDB1);
	FREEMPI(OLDB2);
	return(PERIOD);
}

/* This returns the RCF period-length of sqrt(D) and the least solution of Pell's \pm 1 equation.
* If Eptr = 1, the midpoint and its type are printed. */
MPI  *RCF_PERIOD0(MPI *D, MPI *Eptr, MPI **Xptr, MPI **Yptr)
{
	MPI *P, *Q, *E, *F, *G, *X, *Y, *TEMP, *TEMP1, *H, *PERIOD;
	MPI *T2, *T3, *B1, *B2, *C1, *C2, *OLDB1, *OLDB2, *OLDC1, *OLDC2;
	USL h, k;
	USI t;

	X = BIG_MTHROOT(D, 2);
	G = MULTI(X, X);

	TEMP = ADD0_I(G, (USL)1);
	FREEMPI(G);
	t = EQUALI(D, TEMP);
	FREEMPI(TEMP);
	if (t) /* period 1, exceptional case */
	{
		if (Eptr->S) {
			printf("d=x^2+1\n");
			PRINT_MPI(X, "x");
			printf("x^2-dy^2=1\n");
			printf("y=1\n");
		}
		*Xptr = X;
		*Yptr = ONEI();
		return(ONEI());
	}
	P = ZEROI();
	Q = ONEI();
	h = 0;
	k = 0;
	B1 = ZEROI();
	C1 = ONEI();
	B2 = ONEI();
	C2 = ZEROI();
	while (1)
	{
		OLDC1 = C1;
		OLDC2 = C2;
		OLDB1 = B1;
		OLDB2 = B2;
		TEMP = ADD0I(X, P);
		Y = INT0(TEMP, Q);
		FREEMPI(TEMP);
		F = P;
		TEMP = MULTI(Y, Q);
		P = SUBI(TEMP, P);
		FREEMPI(TEMP);
		t = EQUALI(P, F);
		FREEMPI(F);
		if (h == 0 && k == 0) {
			C1 = Y;
			C2 = ONEI();
		}
		else {
			C1 = MULTABC(B1, Y, C1);
			C2 = MULTABC(B2, Y, C2);
			FREEMPI(Y);
		}
		B1 = OLDC1;
		B2 = OLDC2;

		if (t) { /*P_h=P_{h+1}, even period 2h */
			T2 = CHANGE(2147483648UL);
			T3 = CHANGE(h);
			TEMP1 = MULT_I(T2, k);
			H = ADD0I(TEMP1, T3); /*H = k*2147483648 + h */
			FREEMPI(T2);
			FREEMPI(T3);
			FREEMPI(TEMP1);
			PERIOD = MULT_I(H, 2);
			if (Eptr->S) {
				printf("P[h]=P[h+1]\n");
				PRINT_MPI(H, "h");
				printf("x^2-dy^2=1\n");
			}
			FREEMPI(H);
			*Xptr = MULTAB_PLUS_CDI(OLDC1, C2, OLDB1, OLDC2);
			*Yptr = MULTAB_PLUS_CDI(OLDC2, C2, OLDC2, OLDB2);
			break;
		}
		E = Q;
		TEMP = MULTI(P, P);
		TEMP1 = SUBI(D, TEMP);
		FREEMPI(TEMP);
		Q = INT0(TEMP1, Q);
		t = EQUALI(Q, E);
		FREEMPI(TEMP1);
		FREEMPI(E);
		if (t) { /*Q_h=Q_{h-1}, odd period 2h+1 */
			T2 = CHANGE(2147483648UL);
			T3 = CHANGE(h);
			TEMP1 = MULT_I(T2, k);
			H = ADD0I(TEMP1, T3); /*H = k*2147483648 + h */
			FREEMPI(T2);
			FREEMPI(T3);
			FREEMPI(TEMP1);
			TEMP = MULT_I(H, 2);
			if (Eptr->S) {
				printf("Q[h]=Q[h+1]\n");
				PRINT_MPI(H, "h");
				printf("x^2-dy^2= -1\n");
			}
			FREEMPI(H);
			PERIOD = ADD0_I(TEMP, 1);
			FREEMPI(TEMP);
			*Xptr = MULTAB_PLUS_CDI(C1, C2, OLDC1, OLDC2);
			*Yptr = MULTAB_PLUS_CDI(C2, C2, OLDC2, OLDC2);
			break;
		}
		h++;
		if (h == 2147483648UL) {
			h = 0;
			k++;
			if (k == 2147483648UL) {
				printf("exiting prematurely, k = 2147483648\n");
				return(NULL);
			}
		}
		FREEMPI(OLDB1);
		FREEMPI(OLDB2);
	}
	FREEMPI(X);
	FREEMPI(P);
	FREEMPI(Q);
	FREEMPI(B1);
	FREEMPI(B2);
	FREEMPI(C1);
	FREEMPI(C2);
	FREEMPI(OLDB1);
	FREEMPI(OLDB2);
	return(PERIOD);
}

/* This returns the least solution of Pell's \pm 1 equation. */
void  RCF_PERIOD00(MPI *D, MPI **Xptr, MPI **Yptr)
{
	MPI *P, *Q, *E, *F, *G, *X, *Y, *TEMP, *TEMP1;
	MPI *B1, *B2, *C1, *C2, *OLDB1, *OLDB2, *OLDC1, *OLDC2;
	USI t, flag;

	X = BIG_MTHROOT(D, 2);
	G = MULTI(X, X);
	t = EQUALI(D, G);
	if (t) /* D is a perfect square */
	{
		*Xptr = NULL;
		*Yptr = NULL;
		FREEMPI(G);
		FREEMPI(X);
		return;
	}
	TEMP = ADD0_I(G, (USL)1);
	FREEMPI(G);
	t = EQUALI(D, TEMP);
	FREEMPI(TEMP);
	if (t) /* period 1, exceptional case */
	{
		*Xptr = X;
		*Yptr = ONEI();
		return;
	}
	P = ZEROI();
	Q = ONEI();
	B1 = ZEROI();
	C1 = ONEI();
	B2 = ONEI();
	C2 = ZEROI();
	flag = 0;
	while (1)
	{
		OLDC1 = C1;
		OLDC2 = C2;
		OLDB1 = B1;
		OLDB2 = B2;
		TEMP = ADD0I(X, P);
		Y = INT0(TEMP, Q);
		FREEMPI(TEMP);
		F = P;
		TEMP = MULTI(Y, Q);
		P = SUBI(TEMP, P);
		FREEMPI(TEMP);
		t = EQUALI(P, F);
		FREEMPI(F);
		if (flag == 0) {
			C1 = Y;
			C2 = ONEI();
			flag = 1;
		}
		else {
			C1 = MULTABC(B1, Y, C1);
			C2 = MULTABC(B2, Y, C2);
			FREEMPI(Y);
		}
		B1 = OLDC1;
		B2 = OLDC2;

		if (t) { /*P_h=P_{h+1}, even period 2h */
			*Xptr = MULTAB_PLUS_CDI(OLDC1, C2, OLDB1, OLDC2);
			*Yptr = MULTAB_PLUS_CDI(OLDC2, C2, OLDC2, OLDB2);
			break;
		}
		E = Q;
		TEMP = MULTI(P, P);
		TEMP1 = SUBI(D, TEMP);
		FREEMPI(TEMP);
		Q = INT0(TEMP1, Q);
		t = EQUALI(Q, E);
		FREEMPI(TEMP1);
		FREEMPI(E);
		if (t) { /*Q_h=Q_{h-1}, odd period 2h+1 */
			*Xptr = MULTAB_PLUS_CDI(C1, C2, OLDC1, OLDC2);
			*Yptr = MULTAB_PLUS_CDI(C2, C2, OLDC2, OLDC2);
			break;
		}
		FREEMPI(OLDB1);
		FREEMPI(OLDB2);
	}
	FREEMPI(X);
	FREEMPI(P);
	FREEMPI(Q);
	FREEMPI(B1);
	FREEMPI(B2);
	FREEMPI(C1);
	FREEMPI(C2);
	FREEMPI(OLDB1);
	FREEMPI(OLDB2);
	return;
}

void NSCF_PERIOD00(MPI *D, MPI **Aptr, MPI **Bptr)
/*
* This is a program for finding the fundamenatl solution of Pell's equation using NSCF cfrac of sqrt(D).
*/
{
	MPI *C, *P, *Q, *G, *X, *TEMP, *TEMP1, *TEMP2;
	MPI *OLDP, *OLDQ, *P1, *P2, *Q1, *Q2;
	MPI *B1, *C1, *B2, *C2, *OLDB1, *OLDC1, *OLDB2, *OLDC2, *B;
	USL flag, flag1;
	USI t;
	int a, olda;

	X = BIG_MTHROOT(D, 2);
	G = MULTI(X, X);
	t = EQUALI(D, G);
	if (t) /* D is a perfect square */
	{
		*Aptr = NULL;
		*Bptr = NULL;
		FREEMPI(G);
		FREEMPI(X);
		return;
	}
	TEMP = ADD0_I(G, (USL)1);
	FREEMPI(G);
	t = EQUALI(D, TEMP);
	FREEMPI(TEMP);
	if (t) /* period 1, exceptional case */
	{
		*Aptr = X;
		*Bptr = ONEI();
		return;
	}
	P = ZEROI();
	Q = ONEI();
	flag = 0;
	flag1 = 0;
	B1 = ZEROI();
	C1 = ONEI();
	B2 = ONEI();
	C2 = ZEROI();
	a = 1;
	while (1)
	{
		OLDC1 = C1;
		OLDC2 = C2;
		OLDB1 = B1;
		OLDB2 = B2;
		olda = a;
		TEMP = ADD0I(X, P);
		C = INT0(TEMP, Q);
		FREEMPI(TEMP);
		TEMP = MULTI(C, Q);
		P1 = SUBI(TEMP, P);
		FREEMPI(TEMP);
		TEMP = MULTI(P1, P1);
		Q1 = SUBI(D, TEMP);
		FREEMPI(TEMP);
		TEMP = Q1;
		Q1 = INT0(Q1, Q);
		FREEMPI(TEMP);
		P2 = ADD0I(P1, Q);
		TEMP = ADD0I(P2, P1);
		Q2 = SUB0I(TEMP, Q1);
		FREEMPI(TEMP);
		OLDP = P;
		OLDQ = Q;
		if (RSV(Q1, Q2) <0) {
			a = 1;
			P = P1;
			Q = Q1;
			FREEMPI(P2);
			FREEMPI(Q2);
			B = C;
		}
		else {
			a = -1;
			P = P2;
			Q = Q2;
			FREEMPI(P1);
			FREEMPI(Q1);
			B = ADD0_I(C, 1);
			FREEMPI(C);
		}

		if (flag == 0) {
			C1 = B;
			C2 = ONEI();
			flag = 1;
		}
		else {
			C1 = MULTAB_PLUS_MINUS_CI(B, C1, olda, B1);
			C2 = MULTAB_PLUS_MINUS_CI(B, C2, olda, B2);
			FREEMPI(B);
		}
		B1 = OLDC1;
		B2 = OLDC2;

		t = EQUALI(P, OLDP);
		FREEMPI(OLDP);
		if (t) { /*P_h=P_{h+1}, even period 2h */
			if (olda == 1) {
				*Aptr = MULTAB_PLUS_CDI(C1, OLDC2, OLDC1, OLDB2);
				TEMP = ADDI(C2, OLDB2);
				*Bptr = MULTI(OLDC2, TEMP);
				FREEMPI(TEMP);
			}
			else {
				*Aptr = MULTAB_MINUS_CDI(C1, OLDC2, OLDC1, OLDB2);
				TEMP = SUBI(C2, OLDB2);
				*Bptr = MULTI(OLDC2, TEMP);
				FREEMPI(TEMP);
			}
			FREEMPI(OLDQ);
			break;
		}
		t = EQUALI(Q, OLDQ);
		FREEMPI(OLDQ);
		if (t) { /*Q_h=Q_{h+1}, odd period 2h+1 */
			if (a == 1) {
				*Aptr = MULTAB_PLUS_CDI(C1, C2, OLDC1, OLDC2);
				*Bptr = MULTAB_PLUS_CDI(C2, C2, OLDC2, OLDC2);
			}
			else {
				*Aptr = MULTAB_MINUS_CDI(C1, C2, OLDC1, OLDC2);
				*Bptr = MULTAB_MINUS_CDI(C2, C2, OLDC2, OLDC2);
			}
			break;
		}
		if (flag1) {
			TEMP = SUBI(OLDC1, OLDB1);
			*Aptr = MULTAB_MINUS_CDI(C1, OLDC2, OLDB2, TEMP);
			FREEMPI(TEMP);
			TEMP = MULT_I(OLDC2, 2);
			*Bptr = MULTAB_MINUS_CDI(TEMP, OLDC2, C2, OLDB2);
			FREEMPI(TEMP);
			break;
		}
		TEMP1 = INT0_(OLDQ, (USL)2);
		TEMP2 = ADD0I(Q, TEMP1);
		FREEMPI(TEMP1);
		t = EQUALI(P, TEMP2);
		FREEMPI(TEMP2);
		if (t && a == -1 && flag1 == 0) {
			flag1 = 1;
		}
		FREEMPI(OLDB1);
		FREEMPI(OLDB2);
	}
	FREEMPI(X);
	FREEMPI(P);
	FREEMPI(Q);
	FREEMPI(B1);
	FREEMPI(B2);
	FREEMPI(C1);
	FREEMPI(C2);
	FREEMPI(OLDB1);
	FREEMPI(OLDB2);
	return;
}

/* This returns the least solution of Pell's \pm 1 equation. */
void  RCF_PERIOD000(MPI *D, MPI **Xptr, MPI **Yptr)
{
	MPI *P, *Q, *E, *F, *G, *X, *Y, *TEMP, *TEMP1;
	MPI *B2, *C2, *OLDB2, *OLDC2;
	USI t, flag;

	X = BABY_MTHROOT(D, 2);
	G = MULTI(X, X);
	t = EQUALI(D, G);
	if (t) /* D is a perfect square */
	{
		*Xptr = NULL;
		*Yptr = NULL;
		FREEMPI(G);
		FREEMPI(X);
		return;
	}
	TEMP = ADD0_I(G, (USL)1);
	FREEMPI(G);
	t = EQUALI(D, TEMP);
	FREEMPI(TEMP);
	if (t) /* period 1, exceptional case */
	{
		*Xptr = X;
		*Yptr = ONEI();
		return;
	}
	P = ZEROI();
	Q = ONEI();
	B2 = ONEI();
	C2 = ZEROI();
	flag = 0;
	while (1)
	{
		OLDC2 = C2;
		OLDB2 = B2;
		TEMP = ADD0I(X, P);
		Y = INT0(TEMP, Q);
		FREEMPI(TEMP);
		F = P;
		TEMP = MULTI(Y, Q);
		P = SUBI(TEMP, P);
		FREEMPI(TEMP);
		t = EQUALI(P, F);
		FREEMPI(F);
		if (flag == 0) {
			C2 = ONEI();
			flag = 1;
		}
		else {
			C2 = MULTABC(B2, Y, C2);
		}
		FREEMPI(Y);
		B2 = OLDC2;

		if (t) { /*P_h=P_{h+1}, even period 2h */
			*Yptr = MULTAB_PLUS_CDI(OLDC2, C2, OLDC2, OLDB2);
			*Xptr = SQRT4TIMING(D, X, *Yptr, -1);
			break;
		}
		E = Q;
		TEMP = MULTI(P, P);
		TEMP1 = SUBI(D, TEMP);
		FREEMPI(TEMP);
		Q = INT0(TEMP1, Q);
		t = EQUALI(Q, E);
		FREEMPI(TEMP1);
		FREEMPI(E);
		if (t) { /*Q_h=Q_{h-1}, odd period 2h+1 */
			*Yptr = MULTAB_PLUS_CDI(C2, C2, OLDC2, OLDC2);
			*Xptr = SQRT4TIMING(D, X, *Yptr, 1);
			break;
		}
		FREEMPI(OLDB2);
	}
	FREEMPI(X);
	FREEMPI(P);
	FREEMPI(Q);
	FREEMPI(B2);
	FREEMPI(C2);
	FREEMPI(OLDB2);
	return;
}

void NSCF_PERIOD000(MPI *D, MPI **Aptr, MPI **Bptr)
/*
* This is a program for finding the fundamenatl solution of Pell's equation using NSCF cfrac of sqrt(D).
*/
{
	MPI *C, *P, *Q, *G, *X, *TEMP, *TEMP1, *TEMP2;
	MPI *OLDP, *OLDQ, *P1, *P2, *Q1, *Q2;
	MPI *B2, *C2, *OLDB2, *OLDC2, *B;
	USL flag, flag1;
	USI t;
	int a, olda;

	X = BABY_MTHROOT(D, 2);
	G = MULTI(X, X);
	t = EQUALI(D, G);
	if (t) /* D is a perfect square */
	{
		*Aptr = NULL;
		*Bptr = NULL;
		FREEMPI(G);
		FREEMPI(X);
		return;
	}
	TEMP = ADD0_I(G, (USL)1);
	FREEMPI(G);
	t = EQUALI(D, TEMP);
	FREEMPI(TEMP);
	if (t) /* period 1, exceptional case */
	{
		*Aptr = X;
		*Bptr = ONEI();
		return;
	}
	P = ZEROI();
	Q = ONEI();
	flag = 0;
	flag1 = 0;
	B2 = ONEI();
	C2 = ZEROI();
	a = 1;
	while (1)
	{
		OLDC2 = C2;
		OLDB2 = B2;
		olda = a;
		TEMP = ADD0I(X, P);
		C = INT0(TEMP, Q);
		FREEMPI(TEMP);
		TEMP = MULTI(C, Q);
		P1 = SUBI(TEMP, P);
		FREEMPI(TEMP);
		TEMP = MULTI(P1, P1);
		Q1 = SUBI(D, TEMP);
		FREEMPI(TEMP);
		TEMP = Q1;
		Q1 = INT0(Q1, Q);
		FREEMPI(TEMP);
		P2 = ADD0I(P1, Q);
		TEMP = ADD0I(P2, P1);
		Q2 = SUB0I(TEMP, Q1);
		FREEMPI(TEMP);
		OLDP = P;
		OLDQ = Q;
		if (RSV(Q1, Q2) <0) {
			a = 1;
			P = P1;
			Q = Q1;
			FREEMPI(P2);
			FREEMPI(Q2);
			B = C;
		}
		else {
			a = -1;
			P = P2;
			Q = Q2;
			FREEMPI(P1);
			FREEMPI(Q1);
			B = ADD0_I(C, 1);
			FREEMPI(C);
		}

		if (flag == 0) {
			C2 = ONEI();
			flag = 1;
		}
		else {
			C2 = MULTAB_PLUS_MINUS_CI(B, C2, olda, B2);
		}
		FREEMPI(B);
		B2 = OLDC2;

		t = EQUALI(P, OLDP);
		FREEMPI(OLDP);
		if (t) { /*P_h=P_{h+1}, even period 2h */
			if (olda == 1) {
				TEMP = ADDI(C2, OLDB2);
				*Bptr = MULTI(OLDC2, TEMP);
				FREEMPI(TEMP);
			}
			else {
				TEMP = SUBI(C2, OLDB2);
				*Bptr = MULTI(OLDC2, TEMP);
				FREEMPI(TEMP);
			}
			FREEMPI(OLDQ);
			*Aptr = SQRT4TIMING(D, X, *Bptr, -1);
			break;
		}
		t = EQUALI(Q, OLDQ);
		FREEMPI(OLDQ);
		if (t) { /*Q_h=Q_{h+1}, odd period 2h+1 */
			if (a == 1) {
				*Bptr = MULTAB_PLUS_CDI(C2, C2, OLDC2, OLDC2);
			}
			else {
				*Bptr = MULTAB_MINUS_CDI(C2, C2, OLDC2, OLDC2);
			}
			*Aptr = SQRT4TIMING(D, X, *Bptr, a);
			break;
		}
		if (flag1) {
			TEMP = MULT_I(OLDC2, 2);
			*Bptr = MULTAB_MINUS_CDI(TEMP, OLDC2, C2, OLDB2);
			FREEMPI(TEMP);
			*Aptr = SQRT4TIMING(D, X, *Bptr, 1);
			break;
		}
		TEMP1 = INT0_(OLDQ, (USL)2);
		TEMP2 = ADD0I(Q, TEMP1);
		FREEMPI(TEMP1);
		t = EQUALI(P, TEMP2);
		FREEMPI(TEMP2);
		if (t && a == -1 && flag1 == 0) {
			flag1 = 1;
		}
		FREEMPI(OLDB2);
	}
	FREEMPI(X);
	FREEMPI(P);
	FREEMPI(Q);
	FREEMPI(B2);
	FREEMPI(C2);
	FREEMPI(OLDB2);
	return;
}

MPI *CFRAC_PERIOD1(MPI *D)
/*
* This is a program for finding the period of the cfrac of sqrt(D).
* The algorithm is based on K. Rosen, Elementary number theory
* and its applications, p382, B.A. Venkov, Elementary Number theory, p.62
* and D. Knuth, Art of computer programming, Vol.2, p359, with Pohst's trick
* of using half the period.  We use single precision for the period.
*/
{
	MPI *P, *Q, *E, *F, *G, *X, *Y, *TEMP, *TEMP1, *PERIOD;
	USL h, period;
	USI t;

	X = BIG_MTHROOT(D, 2);
	G = MULTI(X, X);

	TEMP = ADD0_I(G, (USL)1);
	FREEMPI(G);
	t = EQUALI(D, TEMP);
	FREEMPI(TEMP);
	if (t) /* period 1, exceptional case */
	{
		FREEMPI(X);
		return(ONEI());
	}
	P = ZEROI();
	Q = ONEI();
	h = 0;
	while (1)
	{
		TEMP = ADD0I(X, P);
		Y = INT0(TEMP, Q);
		FREEMPI(TEMP);
		F = P;
		TEMP = MULTI(Y, Q);
		FREEMPI(Y);
		P = SUBI(TEMP, P);
		FREEMPI(TEMP);
		t = EQUALI(P, F);
		FREEMPI(F);
		if (t) { /*P_h=P_{h+1}, even period 2h */
				 /*                  printf("P[h]=P[h+1]\n");*/
			if (h<2147483648UL) {
				period = 2 * h;
				PERIOD = CHANGE(period);
			}
			else {
				printf("Exiting prematurely - period >= 2^32\n");
				PERIOD = ZEROI();
			}
			FREEMPI(X);
			FREEMPI(P);
			FREEMPI(Q);
			return(PERIOD);
		}
		E = Q;
		TEMP = MULTI(P, P);
		TEMP1 = SUBI(D, TEMP);
		FREEMPI(TEMP);
		Q = INT0(TEMP1, Q);
		t = EQUALI(Q, E);
		FREEMPI(TEMP1);
		FREEMPI(E);
		if (t) { /*Q_h=Q_{h-1}, odd period 2h+1 */
				 /*printf("Q[h]=Q[h-1]\n");*/
			if (h<2147483647UL) {
				period = 2 * h + 1;
				PERIOD = CHANGE(period);
			}
			else {
				printf("Exiting prematurely - period >= 2^32\n");
				PERIOD = ZEROI();
			}
			FREEMPI(X);
			FREEMPI(P);
			FREEMPI(Q);
			return(PERIOD);
		}
		h++;
	}
}

MPI *NSCF_PERIOD1(MPI *D)
/*
* This is a program for finding the period of the NSCF of sqrt(D).
* using the half period method of nscf_midpoint.pdf.
* Single precision is used for the period calculation.
*/
{
	MPI *C, *P, *Q, *G, *X, *TEMP, *TEMP1, *TEMP2;
	MPI *OLDP, *OLDQ, *P1, *P2, *Q1, *Q2, *PERIOD;
	USL h, flag, period;
	USI t;
	int a;

	X = BIG_MTHROOT(D, 2);
	G = MULTI(X, X);

	TEMP = ADD0_I(G, (USL)1);
	FREEMPI(G);
	t = EQUALI(D, TEMP);
	FREEMPI(TEMP);
	if (t) /* period 1, exceptional case */
	{
		FREEMPI(X);
		return(ONEI());
	}
	P = ZEROI();
	Q = ONEI();
	h = 0;
	flag = 0;
	while (1)
	{
		TEMP = ADD0I(X, P);
		C = INT0(TEMP, Q);
		FREEMPI(TEMP);
		TEMP = MULTI(C, Q);
		P1 = SUBI(TEMP, P);
		FREEMPI(TEMP);
		TEMP = MULTI(P1, P1);
		Q1 = SUBI(D, TEMP);
		FREEMPI(TEMP);
		TEMP = Q1;
		Q1 = INT0(Q1, Q);
		FREEMPI(TEMP);
		P2 = ADD0I(P1, Q);
		TEMP = ADD0I(P2, P1);
		Q2 = SUB0I(TEMP, Q1);
		FREEMPI(TEMP);
		OLDP = P;
		OLDQ = Q;
		if (RSV(Q1, Q2) <0) {
			a = 1;
			P = P1;
			Q = Q1;
			FREEMPI(P2);
			FREEMPI(Q2);
		}
		else {
			a = -1;
			P = P2;
			Q = Q2;
			FREEMPI(P1);
			FREEMPI(Q1);
		}
		FREEMPI(C);
		t = EQUALI(P, OLDP);
		FREEMPI(OLDP);
		if (t) { /*P_h=P_{h+1}, even period 2h */
			printf("P-test, P_h=P_{h+1}\n");
			if (h<2147483648UL) {
				period = 2 * h;
				PERIOD = CHANGE(period);
			}
			else {
				printf("Exiting prematurely - period >= 2^32\n");
				PERIOD = ZEROI();
			}
			FREEMPI(OLDQ);
			break;
		}
		t = EQUALI(Q, OLDQ);
		FREEMPI(OLDQ);
		if (t) { /*Q_h=Q_{h+1}, odd period 2h+1 */
			printf("Q-test:Q_h=Q_{h+1}\n");
			if (h<2147483647UL) {
				period = 2 * h + 1;
				PERIOD = CHANGE(period);
			}
			else {
				printf("Exiting prematurely - period >= 2^32\n");
				PERIOD = ZEROI();
			}
			break;
		}
		if (flag) {
			flag = 0;
			if (h<2147483647UL) {
				period = 2 * h;
				PERIOD = CHANGE(period);
			}
			else {
				printf("Exiting prematurely - period >= 2^32\n");
				PERIOD = ZEROI();
			}
			break;
		}
		TEMP1 = INT0_(OLDQ, (USL)2);
		TEMP2 = ADD0I(Q, TEMP1);
		FREEMPI(TEMP1);
		t = EQUALI(P, TEMP2);
		FREEMPI(TEMP2);
		if (t && a == -1 && flag == 0) {
			printf("PQ-test\n");
			flag = 1;
		}
		h++;
	}
	FREEMPI(X);
	FREEMPI(P);
	FREEMPI(Q);
	return(PERIOD);
}

/* This returns the least solution of Pell's \pm 1 equation.
* We use Jim's identity Q1 = Ax(P-P1)-QQ. */
void  RCF_PERIOD_QIMPROVED(MPI *D, MPI **Xptr, MPI **Yptr)
{
	MPI *P, *Q, *E, *F, *G, *X, *Y, *TEMP, *TEMP1, *QQQ;
	MPI *B1, *B2, *C1, *C2, *OLDB1, *OLDB2, *OLDC1, *OLDC2;
	USI t, flag;

	X = BIG_MTHROOT(D, 2);
	G = MULTI(X, X);
	t = EQUALI(D, G);
	if (t) /* D is a perfect square */
	{
		*Xptr = NULL;
		*Yptr = NULL;
		FREEMPI(G);
		FREEMPI(X);
		return;
	}
	TEMP = ADD0_I(G, (USL)1);
	FREEMPI(G);
	t = EQUALI(D, TEMP);
	FREEMPI(TEMP);
	if (t) /* period 1, exceptional case */
	{
		*Xptr = X;
		*Yptr = ONEI();
		return;
	}
	P = ZEROI();
	Q = ONEI();
	B1 = ZEROI();
	C1 = ONEI();
	B2 = ONEI();
	C2 = ZEROI();
	E = COPYI(D);
	flag = 0;
	while (1)
	{
		OLDC1 = C1;
		OLDC2 = C2;
		OLDB1 = B1;
		OLDB2 = B2;
		TEMP = ADD0I(X, P);
		Y = INT0(TEMP, Q);
		FREEMPI(TEMP);
		F = P;
		TEMP = MULTI(Y, Q);
		P = SUBI(TEMP, P);
		FREEMPI(TEMP);
		QQQ = Q;
		TEMP = SUBI(P, F);
		if (EQONEI(Y)) {
			Q = SUBI(E, TEMP);
		}
		else {
			TEMP1 = MULTI(Y, TEMP);
			Q = SUBI(E, TEMP1);
			FREEMPI(TEMP1);
		}
		FREEMPI(E);
		E = QQQ;
		FREEMPI(TEMP);
		t = EQUALI(P, F);
		FREEMPI(F);
		if (flag == 0) {
			C1 = Y;
			C2 = ONEI();
			flag = 1;
		}
		else {
			C1 = MULTABC(B1, Y, C1);
			C2 = MULTABC(B2, Y, C2);
			FREEMPI(Y);
		}
		B1 = OLDC1;
		B2 = OLDC2;

		if (t) { /*P_h=P_{h+1}, even period 2h */
			*Xptr = MULTAB_PLUS_CDI(OLDC1, C2, OLDB1, OLDC2);
			*Yptr = MULTAB_PLUS_CDI(OLDC2, C2, OLDC2, OLDB2);
			break;
		}
		t = EQUALI(Q, E);
		if (t) { /*Q_h=Q_{h-1}, odd period 2h+1 */
			*Xptr = MULTAB_PLUS_CDI(C1, C2, OLDC1, OLDC2);
			*Yptr = MULTAB_PLUS_CDI(C2, C2, OLDC2, OLDC2);
			break;
		}
		FREEMPI(OLDB1);
		FREEMPI(OLDB2);
	}
	FREEMPI(X);
	FREEMPI(E);
	FREEMPI(P);
	FREEMPI(Q);
	FREEMPI(B1);
	FREEMPI(B2);
	FREEMPI(C1);
	FREEMPI(C2);
	FREEMPI(OLDB1);
	FREEMPI(OLDB2);
	return;
}
