#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_NONSTDC_NO_DEPRECATE
#include <stdio.h>
#include <stdlib.h>
#include "integer.h"
#include "pI.h"
#include "fun.h"
#include "stack.h"
#include "wrappers.h"
#include "calc.h"

/* Free MPI's. Free POLYI's. DO NOT FREE ARRAYS! or VARIABLES!*/

MPI *GCD_W(Stack s)
{
	/* to free or not to free.  That is the question.  I say not! */
	MPI *Aptr = stackPop(s);
	MPI *Bptr = stackPop(s);
	MPI *Result = GCD(Aptr, Bptr);
	FREEMPI(Aptr);
	FREEMPI(Bptr);
	return Result;
}

MPI *LCM_W(Stack s)
{
	/* to free or not to free.  That is the question.  I say not! */
	MPI *Aptr = stackPop(s);
	MPI *Bptr = stackPop(s);
	MPI *Result = LCM(Aptr, Bptr);
	FREEMPI(Aptr);
	FREEMPI(Bptr);
	return Result;
}

/* MPI *LCM_ARRAY(MPI *M[]) */
MPI *LCM_ARRAY_W(Stack s)
{
	MPIA M = stackPop(s);
	MPI *Result;
	Result = LCM_ARRAY(M);
	return Result;
}


/*void SERRET(MPI *P, MPI **Xptr, MPI **Yptr)*/
void SERRET_W(Stack s)
{
	MPI *P = stackPop(s);
	MPI **Xptr = stackPop(s);
	MPI **Yptr = stackPop(s);
	SERRET(P, Xptr, Yptr);
	FREEMPI(P);
	/* Notice I did not free Xptr or Ypr. They are variables */
}

void LAGRANGE_W(Stack s)
{
	POLYI P = stackPop(s);
	MPIA *AA = stackPop(s);
	MPI *M = stackPop(s);
	MPI *ZERO, *TEMP, *ONE, *G0, *G1, *GC, *H, *K, *DIFF;
	MPR *HH;
	POLYI Q, R;
	int d;

	if (DEGREEPI(P)<2) {
		printf("P is linear\n");
		DELETEPI(P);
		FREEMPI(M);
		return;
	}
	TEMP = LEADPI(P);
	d = TEMP->S;
	FREEMPI(TEMP);
	if (d <0) {
		printf("leading coefficient is negative\n");
		DELETEPI(P);
		FREEMPI(M);
		return;
	}
	ZERO = ZEROI();
	TEMP = VALPI(P, ZERO);
	FREEMPI(ZERO);
	d = TEMP->S;
	FREEMPI(TEMP);
	if (d == 0) {
		printf("P(0)=0\n");
		DELETEPI(P);
		FREEMPI(M);
		return;
	}
	Q = DERIVPI(P);
	R = GCDPI(P, Q);
	DELETEPI(Q);
	d = DEGREEPI(R);
	DELETEPI(R);
	if (d>0) {
		printf("P has a repeated root\n");
		DELETEPI(P);
		FREEMPI(M);
		return;
	}
	ONE = ONEI();
	TEMP = VALPI(P, ONE);
	d = TEMP->S;
	FREEMPI(TEMP);
	FREEMPI(ONE);
	if (d == 0) {
		printf("P(1)=0\n");
		DELETEPI(P);
		FREEMPI(M);
		return;
	}
	/* now to check that P has but one positive root, also >1 */
	ZERO = ZEROI();
	ONE = ONEI();
	G0 = STURM_SEQUENCE(P, ZERO, ZERO);
	G1 = STURM_SEQUENCE(P, ONE, ZERO);
	HH = CAUCHY(P);
	H = COPYI(HH->N);
	FREEMPR(HH);

	K = VALPI(P, H);
	d = K->S;
	FREEMPI(K);
	if (d == 0) {
		TEMP = H;
		H = ADD0I(H, ONE);
		FREEMPI(TEMP);
	}
	GC = STURM_SEQUENCE(P, H, ZERO);
	FREEMPI(ZERO);
	FREEMPI(H);
	if (EQUALI(G0, G1)) {
		DIFF = SUB0I(G1, GC);
		if (RSV(DIFF, ONE)>0) {
			printf("P has ");
			PRINTI(DIFF);
			printf(" roots > 1\n");
			FREEMPI(DIFF);
			FREEMPI(ONE);
			DELETEPI(P);
			FREEMPI(M);
			FREEMPI(G0);
			FREEMPI(G1);
			FREEMPI(GC);
			return;
		}
		if (DIFF->S == 0) {
			printf("P has no roots > 1\n");
			FREEMPI(DIFF);
			FREEMPI(ONE);
			DELETEPI(P);
			FREEMPI(M);
			FREEMPI(G0);
			FREEMPI(G1);
			FREEMPI(GC);
			return;
		}
	}
	else {
		printf("a has a root between 0 and 1\n");
		return;
	}
	FREEMPI(G0);
	FREEMPI(G1);
	FREEMPI(GC);
	FREEMPI(ONE);
	FREEMPI(DIFF);
	/* now P has exactly one positve root and it's > 1 */

	LAGRANGE(P, AA, M);
	DELETEPI(P);
	FREEMPI(M);
	return;
}

/* void EUCLID(MPI *Aptr, MPI *Bptr, MPIA *Q, MPIA *R, MPIA *S, MPIA *T, MPI **Dptr) */
void EUCLID_W(Stack s)
{
	MPI *Aptr = stackPop(s);
	MPI *Bptr = stackPop(s);
	MPIA *Q = stackPop(s);
	MPIA *R = stackPop(s);
	MPIA *S = stackPop(s);
	MPIA *T = stackPop(s);
	MPI **Dptr = stackPop(s);
	EUCLID(Aptr, Bptr, Q, R, S, T, Dptr);
	FREEMPI(Aptr);
	FREEMPI(Bptr);
}

/*void CONVERGENTS(MPIA A, MPIA *P, MPIA *Q) */
void CONVERGENTS_W(Stack s)
{
	MPIA A = stackPop(s);
	MPIA *P = stackPop(s);
	MPIA *Q = stackPop(s);
	CONVERGENTS(A, P, Q);
}

/* MPI *EUCLIDI(MPI *Pptr, MPI *Qptr, MPI **Hptr, MPI **Kptr) */
MPI *EUCLIDI_W(Stack s)
{
	MPI *P = stackPop(s);
	MPI *Q = stackPop(s);
	MPI **H = stackPop(s);
	MPI **K = stackPop(s);
	MPI *Result = EUCLIDI(P, Q, H, K);
	FREEMPI(P);
	FREEMPI(Q);
	return Result;
}

/*MPI *JACOB(MPI *M, MPI *N) */
MPI *JACOB_W(Stack s)
{
	MPI *M = stackPop(s);
	MPI *N = stackPop(s);
	MPI *Result = JACOB(M, N);
	FREEMPI(M);
	FREEMPI(N);
	return Result;
}

/* MPI *SQRTMX(MPI *x, MPI *p) */
MPI *SQRTMX_W(Stack s)
{
	MPI *x = stackPop(s);
	MPI *p = stackPop(s);
	MPI *Result = SQRTMX(x, p);
	if (Result == NULL)
		rettype = FUNC_FAIL;
	FREEMPI(x);
	FREEMPI(p);
	return Result;
}

/* MPI *GCD_ARRAY(MPIA M) */
MPI *GCD_ARRAY_W(Stack s)
{
	MPIA M = stackPop(s);
	return GCD_ARRAY(M);
}

/* MPI *GCD_ARRAYVX(MPIA M, MPIA *Y) */
MPI *GCD_ARRAYVX_W(Stack s)
{
	MPIA M = stackPop(s);
	MPIA* Y = stackPop(s);
	return GCD_ARRAYVX(M, Y);
}
/* MPI *CONGRX(MPI *A, MPI *B, MPI *M, MPI **N) */
MPI *CONGRX_W(Stack s)
{
	MPI *A = stackPop(s);
	MPI *B = stackPop(s);
	MPI *M = stackPop(s);
	MPI **N = stackPop(s);
	MPI *Result = CONGRX(A, B, M, N);
	if (Result == NULL)
		rettype = FUNC_FAIL; /* do not print value.  see parse.y */
	FREEMPI(A);
	FREEMPI(B);
	FREEMPI(M);
	return Result;
}
/* MPI *CHINESEX(MPI *A, MPI *B, MPI *M, MPI *N, MPI **Mptr) */
MPI *CHINESEX_W(Stack s)
{
	MPI *A = stackPop(s);
	MPI *B = stackPop(s);
	MPI *M = stackPop(s);
	MPI *N = stackPop(s);
	MPI **Mptr = stackPop(s);
	MPI *Result = CHINESEX(A, B, M, N, Mptr);
	if (Result == NULL)
		rettype = FUNC_FAIL;
	FREEMPI(A);
	FREEMPI(B);
	FREEMPI(M);
	FREEMPI(N);
	return Result;
}

/* MPI *CHINESE_ARRAYX(MPIA A, MPIA M, MPI **Mptr) */
MPI *CHINESE_ARRAYX_W(Stack s)
{
	MPIA A = stackPop(s);
	MPIA M = stackPop(s);
	MPI **Mptr = stackPop(s);
	MPI *Result = CHINESE_ARRAYX(A, M, Mptr);
	return Result;
}
/* MPI *BIG_MTHROOTX(MPI *Aptr, MPI *M) */
MPI *BIG_MTHROOTX_W(Stack s)
{
	MPI *Aptr = stackPop(s);
	MPI *M = stackPop(s);
	MPI *Result = BIG_MTHROOTX(Aptr, M);
	if (Result == NULL)
		rettype = FUNC_FAIL;
	FREEMPI(Aptr);
	FREEMPI(M);
	return Result;
}
/* MPI  *FUND_UNITX(MPI *D, MPI **Xptr, MPI **Yptr) */
MPI *FUND_UNITX_W(Stack s)
{
	MPI *D = stackPop(s);
	MPI **Xptr = stackPop(s);
	MPI **Yptr = stackPop(s);
	MPI *Result = FUND_UNITX(D, Xptr, Yptr);
	if (Result == NULL)
		rettype = FUNC_FAIL;
	FREEMPI(D);
	return Result;
}
/* MPI  *PELX(MPI *D, MPI *E, MPI **Xptr, MPI **Yptr) */
MPI *PELX_W(Stack s)
{
	MPI *D = stackPop(s);
	MPI *E = stackPop(s);
	MPI **Xptr = stackPop(s);
	MPI **Yptr = stackPop(s);
	MPI *Result = PELX(D, E, Xptr, Yptr);
	if (Result == NULL)
		rettype = FUNC_FAIL;
	FREEMPI(D);
	FREEMPI(E);
	return Result;
}
/* MPI *SURDX(MPI *D, MPI *T, MPI *U, MPI *V, MPIA *A_ARRAY, MPIA *U_ARRAY, MPIA *V_ARRAY, MPIA *P_ARRAY, MPIA *Q_ARRAY) */
MPI *SURDX_W(Stack s)
{
	MPI *D = stackPop(s);
	MPI *T = stackPop(s);
	MPI *U = stackPop(s);
	MPI *V = stackPop(s);
	MPIA *A_ARRAY = stackPop(s);
	MPIA *U_ARRAY = stackPop(s);
	MPIA *V_ARRAY = stackPop(s);
	MPIA *P_ARRAY = stackPop(s);
	MPIA *Q_ARRAY = stackPop(s);

	MPI *Result = SURDX(D, T, U, V, A_ARRAY, U_ARRAY, V_ARRAY, P_ARRAY, Q_ARRAY);
	if (Result == NULL)
		rettype = FUNC_FAIL;
	FREEMPI(D);
	FREEMPI(T);
	FREEMPI(U);
	FREEMPI(V);
	return Result;
}

/* MPI *PATZX(MPI *D, MPI *N)*/

void PATZ_W(Stack s)
{
	MPI *D = stackPop(s);
	MPI *N = stackPop(s);
	MPI *Result = PATZX(D, N);
	if (Result == NULL)
		rettype = FUNC_FAIL;
	FREEMPI(Result);
	FREEMPI(D);
	FREEMPI(N);
	return;
}

/* MPI *MPOWERX(MPI *Aptr, MPI *Bptr, MPI *Cptr) */
MPI *MPOWERX_W(Stack s)
{
	MPI *Aptr = stackPop(s);
	MPI *Bptr = stackPop(s);
	MPI *Cptr = stackPop(s);
	MPI *Result = MPOWERX(Aptr, Bptr, Cptr);
	if (Result == NULL)
		rettype = FUNC_FAIL;
	FREEMPI(Aptr);
	FREEMPI(Bptr);
	FREEMPI(Cptr);
	return Result;
}
/* MPI *Nextprime(MPI *N) */
MPI *Nextprime_W(Stack s)
{
	MPI *N = stackPop(s);
	MPI *Result = Nextprime(N);
	if (Result == NULL)
		rettype = FUNC_FAIL;
	FREEMPI(N);
	return Result;
}

/* MPI *INVERSEMX(MPI *A, MPI *M) */
MPI *INVERSEMX_W(Stack s)
{
	MPI *A = stackPop(s);
	MPI *M = stackPop(s);
	MPI *Result = INVERSEMX(A, M);
	if (Result == NULL)
		rettype = FUNC_FAIL;
	FREEMPI(A);
	FREEMPI(M);
	return Result;
}
/* MPI *FACTORX(MPI *N) */
MPI *FACTORX_W(Stack s)
{
	MPI *N = stackPop(s);
	MPI *Result = FACTORX(N);
	if (Result == NULL)
		rettype = FUNC_FAIL;
	FREEMPI(N);
	return Result;
}
/* MPI *DIVISORX(MPI *N) */
MPI *DIVISORX_W(Stack s)
{
	MPI *N = stackPop(s);
	MPI *Result = DIVISORX(N);
	if (Result == NULL)
		rettype = FUNC_FAIL;
	FREEMPI(N);
	return Result;
}
/* MPI *MOBIUSX(MPI *N) */
MPI *MOBIUSX_W(Stack s)
{
	MPI *N = stackPop(s);
	MPI *Result = MOBIUSX(N);
	if (Result == NULL)
		rettype = FUNC_FAIL;
	FREEMPI(N);
	return Result;
}

/* MPI *EULERX(MPI *N) */
MPI *EULERX_W(Stack s)
{
	MPI *N = stackPop(s);
	MPI *Result = EULERX(N);
	if (Result == NULL)
		rettype = FUNC_FAIL;
	FREEMPI(N);
	return Result;
}
/* MPI *SIGMAX(MPI *N) */
MPI *SIGMAX_W(Stack s)
{
	MPI *N = stackPop(s);
	MPI *Result = SIGMAX(N);
	if (Result == NULL)
		rettype = FUNC_FAIL;
	FREEMPI(N);
	return Result;
}
/* MPI *LPRIMROOTX(MPI *P) */
MPI *LPRIMROOTX_W(Stack s)
{
	MPI *N = stackPop(s);
	MPI *Result = LPRIMROOTX(N);
	if (Result == NULL)
		rettype = FUNC_FAIL;
	FREEMPI(N);
	return Result;

}
/* MPI *ORDERMX(MPI *A, MPI *M*/
MPI *ORDERMX_W(Stack s)
{
	MPI *A = stackPop(s);
	MPI *M = stackPop(s);
	MPI *Result = ORDERMX(A, M);
	if (Result == NULL)
		rettype = FUNC_FAIL;
	FREEMPI(A);
	FREEMPI(M);
	return Result;
}
/* MPI *LUCASX(MPI *N) */
MPI *LUCASX_W(Stack s)
{
	MPI *N = stackPop(s);
	MPI *Result = LUCASX(N);
	if (Result == NULL)
		rettype = FUNC_FAIL;
	FREEMPI(N);
	return Result;
}
/* MPI *LENGTHX(MPI *N) */
MPI *LENGTHX_W(Stack s)
{
	MPI *N = stackPop(s);
	MPI *Result = LENGTHX(N);
	FREEMPI(N);
	return Result;
}
/* MPI *RSAE(MPI *Pptr, MPI *Qptr) */
MPI *RSAE_W(Stack s)
{
	MPI *P = stackPop(s);
	MPI *Q = stackPop(s);
	MPI *Result = RSAEX(P, Q);
	if (Result == NULL)
		printf("either p or q is not a prime <= 355142\n");
	rettype = FUNC_FAIL;
	FREEMPI(P);
	FREEMPI(Q);
	return Result;
}
/* MPI *NEXTPRIMEAPX(MPI *A, MPI *B, MPI *M) */
MPI *NEXTPRIMEAPX_W(Stack s)
{
	MPI *A = stackPop(s);
	MPI *B = stackPop(s);
	MPI *M = stackPop(s);
	MPI *Result = NEXTPRIMEAPX(A, B, M);
	if (Result == NULL)
		rettype = FUNC_FAIL;
	FREEMPI(A);
	FREEMPI(B);
	FREEMPI(M);
	return Result;
}
/* MPI *POLLARD(MPI *Nptr) */
MPI *POLLARD_W(Stack s)
{
	MPI *N = stackPop(s);
	MPI *Result = POLLARD(N);
	FREEMPI(N);
	return Result;
}
/* MPI *EFACTORX(MPI *N, MPI *M, MPI *P) */
MPI *EFACTORX_W(Stack s)
{
	MPI *N = stackPop(s);
	MPI *M = stackPop(s);
	MPI *P = stackPop(s);
	MPI *Result = EFACTORX(N, M, P);
	if (Result == NULL)
		rettype = FUNC_FAIL;
	FREEMPI(N);
	FREEMPI(M);
	FREEMPI(P);
	return Result;
}
/* MPI *PERFECT_POWER(MPI *N) */
MPI *PERFECT_POWER_W(Stack s)
{
	MPI *N = stackPop(s);
	MPI *Result = PERFECT_POWER(N);
	FREEMPI(N);
	return Result;
}

/*MPI *LEASTQNRX(MPI *P) */
MPI *LEASTQNRX_W(Stack s)
{
	MPI *P = stackPop(s);
	MPI *Result = LEASTQNRX(P);
	FREEMPI(P);
	if (Result == NULL)
		rettype = FUNC_FAIL;
	return Result;
}
MPI* CONTENTPI2_W(Stack s)
{
	POLYI P = stackPop(s);
	MPI *Result = CONTENTPI2(P);
	DELETEPI(P);
	return Result;
}

/* void COLLATZ(MPI *Dptr, MPI *Eptr) */
void COLLATZ_W(Stack s)
{
	MPI *D = stackPop(s);
	MPI *E = stackPop(s);
	COLLATZ(D, E);
	FREEMPI(D);
	FREEMPI(E);

}
/* void MTHROOTX(MPI *Aptr, MPI *Bptr, MPI *M, MPI *R) */
void MTHROOTX_W(Stack s)
{
	MPI *A = stackPop(s);
	MPI *B = stackPop(s);
	MPI *M = stackPop(s);
	MPI *R = stackPop(s);
	MTHROOTX(A, B, M, R);
	FREEMPI(A);
	FREEMPI(B);
	FREEMPI(M);
	FREEMPI(R);
}
/* void MILLERX(MPI *N, MPI *B) */
void MILLERX_W(Stack s)
{
	MPI *N = stackPop(s);
	MPI *B = stackPop(s);
	MILLERX(N, B);
	FREEMPI(N);
	FREEMPI(B);
}
/* void JUGGLER(MPI *Dptr, MPI *Iptr) */
void JUGGLER_W(Stack s)
{
	MPI *D = stackPop(s);
	MPI *I = stackPop(s);
	JUGGLER(D, I);
	FREEMPI(D);
	FREEMPI(I);
}

/* void ENCODE(MPI *Eptr, MPI *Rptr) */
void ENCODE_W(Stack s)
{
	MPI *E = stackPop(s);
	MPI *R = stackPop(s);
	ENCODE(E, R);
	FREEMPI(E);
	FREEMPI(R);
}
/* void DECODEX(MPI *Eptr, MPI *Pptr, MPI *Qptr) */
void DECODEX_W(Stack s)
{
	MPI *E = stackPop(s);
	MPI *P = stackPop(s);
	MPI *Q = stackPop(s);
	DECODEX(E, P, Q);
	FREEMPI(E);
	FREEMPI(P);
	FREEMPI(Q);

}
/* void SCHNORRGCD(MPI *N) */
void SCHNORRGCD_W(Stack s)
{
	MPI *N = stackPop(s);
	SCHNORRGCD(N);
	FREEMPI(N);
}

/* void SCHNORRHERMITE(MPI *N) */
void SCHNORRHERMITE_W(Stack s)
{
	MPI *N = stackPop(s);
	SCHNORRHERMITE(N);
	FREEMPI(N);
}

POLYI PRIMITIVEPI_W(Stack s)
{
	POLYI P = stackPop(s);
	POLYI Result = PRIMITIVEPI(P);
	DELETEPI(P);
	return Result;
}

/* MPI *CORNACCHIAX(MPI *A, MPI *B, MPI *M) */
void CORNACCHIA_W(Stack S)
{
	MPI *A = stackPop(S);
	MPI *B = stackPop(S);
	MPI *M = stackPop(S);

	MPI *Result = CORNACCHIAX(A, B, M);
	if (Result == NULL)
		rettype = FUNC_FAIL;
	FREEMPI(Result);
	FREEMPI(A);
	FREEMPI(B);
	FREEMPI(M);
	return;
}

MPI *SQROOT_W(Stack S)
{
	MPI *Tmp;
	USI ll;

	MPI *A = stackPop(S);
	MPI *N = stackPop(S);
	MPIA *Y = stackPop(S);
	MPI **M = stackPop(S);
	MPI **L = stackPop(S);
	MPI *R = SQROOTX(A, N, Y, M, &ll);
	FREEMPI(A);
	FREEMPI(N);
	if (R == NULL)
		rettype = FUNC_FAIL;
	else {
		*L = CHANGE((USL)ll);
		if (EQMINUSONEI(R)) {
			Tmp = R;
			R = ZEROI();
			FREEMPI(Tmp);
		}
	}
	return (R);
}

MPI *QUADRATIC_W(Stack S)
{
	MPI *A, *B, *C, *N, *R;
	MPIA *Y;

	A = stackPop(S);
	B = stackPop(S);
	C = stackPop(S);
	N = stackPop(S);
	Y = stackPop(S);
	R = QUADRATICX(A, B, C, N, Y);
	FREEMPI(A);
	FREEMPI(B);
	FREEMPI(C);
	FREEMPI(N);
	if (R == NULL) {
		R = ZEROI();
		rettype = FUNC_FAIL;
	}
	return (R);
}
/*
void BINARYFORM_W(Stack S)
{
MPI *A, *B, *C, *N;

A = stackPop(S);
B = stackPop(S);
C = stackPop(S);
N = stackPop(S);
BINARYFORM(A, B, C, N);
FREEMPI(A);
FREEMPI(B);
FREEMPI(C);
FREEMPI(N);
return;
}
*/
void GAUSS_W(Stack S)
{
	MPI *A, *B, *C, *N, **alpha, **gamma, **M;
	A = stackPop(S);
	B = stackPop(S);
	C = stackPop(S);
	N = stackPop(S);
	alpha = stackPop(S);
	gamma = stackPop(S);

	M = stackPop(S);
	GAUSS(A, B, C, N, alpha, gamma, M);
	FREEMPI(A);
	FREEMPI(B);
	FREEMPI(C);
	FREEMPI(N);
	return;
}

void BINFORM_W(Stack S)
{
	MPI *A, *B, *C, *N, *V;

	A = stackPop(S);
	B = stackPop(S);
	C = stackPop(S);
	N = stackPop(S);
	V = stackPop(S);
	BINARYFORM1(A, B, C, N, V);
	FREEMPI(A);
	FREEMPI(B);
	FREEMPI(C);
	FREEMPI(N);
	FREEMPI(V);
	return;
}

MPI *HALFMOD_W(Stack S)
{
	MPI *A, *B, *R;

	A = stackPop(S);
	B = stackPop(S);
	R = HALFMODX(A, B);
	FREEMPI(A);
	FREEMPI(B);
	if (R == NULL) {
		rettype = FUNC_FAIL;
	}
	return (R);
}

/* MPI *LOGX(MPI *A, MPI *B, MPI *D, MPI *R, MPIA *M, MPI **L) */
void LOG_W(Stack s)
{
	MPI *A = stackPop(s);
	MPI *B = stackPop(s);
	MPI *D = stackPop(s);
	MPI *R = stackPop(s);
	MPIA *M = stackPop(s);
	MPI **L = stackPop(s);
	MPI *Z;

	Z = LOGX(A, B, D, R, M, L);
	if (Z == NULL)
		rettype = FUNC_FAIL;
	else
		FREEMPI(Z);
	FREEMPI(A);
	FREEMPI(B);
	FREEMPI(D);
	FREEMPI(R);
	return;
}

/* MPI *CEILINGIX(MPI *A, MPI *B) */
MPI *CEILINGI_W(Stack s)
{
	MPI *A = stackPop(s);
	MPI *B = stackPop(s);

	MPI *Z = CEILINGIX(A, B);
	if (Z == NULL)
		rettype = FUNC_FAIL;
	FREEMPI(A);
	FREEMPI(B);
	return(Z);
}

/*void TESTLOG1(MPI *A, MPI *B, MPI *D, MPI *R)*/
void TESTLOG1_W(Stack s)
{
	MPI *A = stackPop(s);
	MPI *B = stackPop(s);
	MPI *D = stackPop(s);
	MPI *R = stackPop(s);

	TESTLOG1(A, B, D, R);
	FREEMPI(A);
	FREEMPI(B);
	FREEMPI(D);
	FREEMPI(R);
}

/*void TESTLOGX(MPI *A, MPI *B, MPI *D, MPI *M, MPI *N)*/
void TESTLOG_W(Stack s)
{
	MPI *A = stackPop(s);
	MPI *B = stackPop(s);
	MPI *D = stackPop(s);
	MPI *M = stackPop(s);
	MPI *N = stackPop(s);
	MPI *Z;

	Z = TESTLOGX(A, B, D, M, N);
	if (Z == NULL)
		rettype = FUNC_FAIL;
	else
		FREEMPI(Z);

	FREEMPI(A);
	FREEMPI(B);
	FREEMPI(D);
	FREEMPI(M);
	FREEMPI(N);
}

/* MPI *SUBRESULTANT(POLYI P, POLY Q) */
MPI *SUBRESULTANT_W(Stack s)
{
	MPI *Result;
	POLYI P = stackPop(s);
	POLYI Q = stackPop(s);
	if (DEGREEPI(P) <= 0) {
		printf("First argument is a constant\n");
		rettype = FUNC_FAIL;
	}
	if (DEGREEPI(Q) <= 0) {
		printf("Second argument is a constant\n");
		rettype = FUNC_FAIL;
	}
	Result = SUBRESULTANT(P, Q);
	DELETEPI(P);
	DELETEPI(Q);
	return Result;
}


/* MPI *DISCRIMINANTPI(POLYI Pptr)*/
MPI *DISCRIMINANTPI_W(Stack s) {
	MPI *Result;
	POLYI P = stackPop(s);
	if (DEGREEPI(P) <= 1) {
		printf("argument has form aX+b\n");
		rettype = FUNC_FAIL;
	}
	Result = DISCRIMINANTPI(P);
	DELETEPI(P);
	return Result;
}

/* POLYI DERIVPI(POLYI Pptr) */
POLYI DERIVPI_W(Stack s) {
	POLYI Result;
	POLYI P = stackPop(s);
	Result = DERIVPI(P);
	DELETEPI(P);
	return Result;
}

/* MPI *PRIME_GENERATORX(MPI *M, MPI *N) */
MPI *PRIME_GENERATOR_W(Stack s) {
	MPI *Result;

	MPI *M = stackPop(s);
	MPI *N = stackPop(s);
	Result = PRIME_GENERATORX(M, N);
	if (Result == NULL) {
		FREEMPI(M);
		FREEMPI(N);
		rettype = FUNC_FAIL;
	}
	FREEMPI(M);
	FREEMPI(N);
	return(Result);
}

/* POLYI GCDPI(POLYI Pptr, POLYI Qptr) */
POLYI GCDPI_W(Stack s) {
	POLYI Result;
	POLYI P = stackPop(s);
	POLYI Q = stackPop(s);
	Result = GCDPI(P, Q);
	DELETEPI(P);
	DELETEPI(Q);
	return Result;
}


/* MPI *STURM_SEQUENCE(POLYI P, MPI *B, MPI *E)*/
MPI *STURM_SEQUENCE_W(Stack s) {
	MPI *Result, *TEMP;
	POLYI R, Q;
	int d;
	POLYI P = stackPop(s);
	MPI *B = stackPop(s);
	MPI *E = stackPop(s);

	Q = DERIVPI(P);
	R = GCDPI(P, Q);
	d = DEGREEPI(R);
	DELETEPI(Q);
	DELETEPI(R);
	if (d) {
		printf("P has a multiple root\n");
		rettype = FUNC_FAIL;
		DELETEPI(P);
		FREEMPI(B);
		FREEMPI(E);
		return NULL;
	}
	else if (DEGREEPI(P) <2) {
		printf("P is linear\n");
		rettype = FUNC_FAIL;
		DELETEPI(P);
		FREEMPI(B);
		FREEMPI(E);
		return NULL;
	}
	else {
		TEMP = VALPI(P, B);
		d = TEMP->S;
		FREEMPI(TEMP);
		if (d == 0) {
			printf("P(b)=0\n");
			rettype = FUNC_FAIL;
			DELETEPI(P);
			FREEMPI(B);
			FREEMPI(E);
			return NULL;
		}
		else {
			Result = STURM_SEQUENCE(P, B, E);
			DELETEPI(P);
			FREEMPI(B);
			FREEMPI(E);
			return Result;
		}
	}
}

POLYI CYCLOTOMIC_W(Stack s) {
	POLYI P;
	MPI *N = stackPop(s);
	P = CYCLOTOMIC(N);
	if (P == NULL) {
		rettype = FUNC_FAIL;
	}
	FREEMPI(N);
	return(P);
}

MPI *POS_W(Stack s) {
	MPI *P;
	MPI *D = stackPop(s);
	P = POSX(D);
	if (P == NULL)
		rettype = FUNC_FAIL;
	FREEMPI(D);
	return(P);
}

MPI *NEG_W(Stack s) {
	MPI *H;
	MPI *D = stackPop(s);
	MPI *FLAG = stackPop(s);
	H = NEGX(D, FLAG);
	if (H == NULL)
		rettype = FUNC_FAIL;
	FREEMPI(D);
	FREEMPI(FLAG);
	return(H);
}

MPI *NEARINT_W(Stack s) {
	MPI *A, *B, *R;

	A = stackPop(s);
	B = stackPop(s);
	R = NEARINTX(A, B);
	FREEMPI(A);
	FREEMPI(B);
	if (R == NULL) {
		rettype = FUNC_FAIL;
	}
	return (R);
}

MPI *REDUCE_NEG_W(Stack s) {
	MPI *A, *B, *C, *I;
	A = stackPop(s);
	B = stackPop(s);
	C = stackPop(s);
	I = REDUCE_NEGX(A, B, C);
	FREEMPI(A);
	FREEMPI(B);
	FREEMPI(C);
	if (I == NULL) {
		rettype = FUNC_FAIL;
	}
	return(I);
}

MPI *REDUCE_POS_W(Stack s) {
	MPI *I;
	MPI *A = stackPop(s);
	MPI *B = stackPop(s);
	MPI *C = stackPop(s);
	I = REDUCE_POSX(A, B, C);
	FREEMPI(A);
	FREEMPI(B);
	FREEMPI(C);
	if (I == NULL) {
		rettype = FUNC_FAIL;
	}
	return(I);
}

MPI *POS0_W(Stack s) {
	MPI *P;
	MPI *D = stackPop(s);
	P = POS0X(D);
	if (P == NULL)
		rettype = FUNC_FAIL;
	FREEMPI(D);
	return(P);
}

MPI *TABLENEG_W(Stack s) {
	MPI *P;
	MPI *M = stackPop(s);
	MPI *N = stackPop(s);
	P = TABLENEGX(M, N);
	if (P == NULL)
		rettype = FUNC_FAIL;
	FREEMPI(M);
	FREEMPI(N);
	return(P);
}

MPI *TABLEPOS_W(Stack s) {
	MPI *P;
	MPI *M = stackPop(s);
	MPI *N = stackPop(s);
	P = TABLEPOSX(M, N);
	if (P == NULL)
		rettype = FUNC_FAIL;
	FREEMPI(M);
	FREEMPI(N);
	return(P);
}

MPI *DAVISON_W(Stack s) {
	MPI *T;
	MPI *L = stackPop(s);
	MPI *M = stackPop(s);
	MPI *N = stackPop(s);
	T = DAVISONX(L, M, N);
	if (T == NULL)
		rettype = FUNC_FAIL;
	FREEMPI(L);
	FREEMPI(M);
	FREEMPI(N);
	return(T);
}

MPI *RANEY1_W(Stack s) {
	MPI *T;
	MPI *P = stackPop(s);
	MPI *Q = stackPop(s);
	MPI *R = stackPop(s);
	MPI *S = stackPop(s);
	T = RANEY1X(P, Q, R, S);
	if (T == NULL)
		rettype = FUNC_FAIL;
	FREEMPI(P);
	FREEMPI(Q);
	FREEMPI(R);
	FREEMPI(S);
	return(T);
}

MPI *UNIMODULAR_W(Stack s) {
	MPI *T;
	MPI *P = stackPop(s);
	MPI *Q = stackPop(s);
	MPI *R = stackPop(s);
	MPI *S = stackPop(s);
	T = UNIMODULARX(P, Q, R, S);
	if (T == NULL)
		rettype = FUNC_FAIL;
	FREEMPI(P);
	FREEMPI(Q);
	FREEMPI(R);
	FREEMPI(S);
	return(T);
}

void TWOADICSQRT_W(Stack s) {
	MPI *T;
	MPI *A = stackPop(s);
	MPI *N = stackPop(s);
	MPIA *DIGITS = stackPop(s);
	T = TWOADICSQRTX(A, N, DIGITS);
	if (T == NULL)
		rettype = FUNC_FAIL;
	FREEMPI(A);
	FREEMPI(N);
	FREEMPI(T);
}

void PADICSQRT_W(Stack s) {
	MPI *T;
	MPI *A = stackPop(s);
	MPI *N = stackPop(s);
	MPI *P = stackPop(s);
	MPIA *DIGITS = stackPop(s);
	T = PADICSQRTX(A, N, P, DIGITS);
	if (T == NULL)
		rettype = FUNC_FAIL;
	FREEMPI(A);
	FREEMPI(N);
	FREEMPI(P);
	FREEMPI(T);
}

MPI *SIGMAK_W(Stack s) {
	MPI *T;
	MPI *K = stackPop(s);
	MPI *N = stackPop(s);
	T = SIGMAKX(K, N);
	if (T == NULL)
		rettype = FUNC_FAIL;
	FREEMPI(K);
	FREEMPI(N);
	return(T);
}

MPI *TAU_W(Stack s) {
	MPI *T;
	MPI *N = stackPop(s);
	T = TAUX(N);
	if (T == NULL)
		rettype = FUNC_FAIL;
	FREEMPI(N);
	return(T);
}

MPI *TAU_PRIMEPOWER_W(Stack s) {
	MPI *T;
	MPI *N = stackPop(s);
	MPI *P = stackPop(s);
	T = TAU_PRIMEPOWERX(N, P);
	if (T == NULL)
		rettype = FUNC_FAIL;
	FREEMPI(N);
	FREEMPI(P);
	return(T);
}

MPI *TAU_COMPOSITE_W(Stack s) {
	MPI *T;
	MPI *N = stackPop(s);
	T = TAU_COMPOSITEX(N);
	if (T == NULL)
		rettype = FUNC_FAIL;
	FREEMPI(N);
	return(T);
}

MPI *REP_DEFINITE_W(Stack s) {
	MPI *A, *B, *C, *M, *PRINT_FLAG, *I;
	A = stackPop(s);
	B = stackPop(s);
	C = stackPop(s);
	M = stackPop(s);
	PRINT_FLAG = stackPop(s);
	I = REP_DEFINITEX(A, B, C, M, PRINT_FLAG);
	FREEMPI(A);
	FREEMPI(B);
	FREEMPI(C);
	FREEMPI(M);
	FREEMPI(PRINT_FLAG);
	if (I == NULL) {
		rettype = FUNC_FAIL;
	}
	return(I);
}

void POWERD_W(Stack s)
{
	MPI *A = stackPop(s);
	MPI *B = stackPop(s);
	MPI *D = stackPop(s);
	MPI *N = stackPop(s);
	MPI **AA = stackPop(s);
	MPI **BB = stackPop(s);
	POWERD(A, B, D, N, AA, BB);
	FREEMPI(A);
	FREEMPI(B);
	FREEMPI(D);
	FREEMPI(N);
}

MPI *EUCLIDI1_W(Stack s)
{
	MPI *Aptr = stackPop(s);
	MPI *Bptr = stackPop(s);
	MPI *Result;
	if (Aptr->S <= 0) {
		Result = ZEROI();
		printf("A <= 0\n");
	}
	else if (Bptr->S <= 0) {
		Result = ZEROI();
		printf("B <= 0\n");
	}
	else {
		Result = EUCLIDI1(Aptr, Bptr);
	}
	FREEMPI(Aptr);
	FREEMPI(Bptr);
	return Result;
}

MPI *CFRACN_W(Stack s)
{
	MPI *TEMP1, *TEMP2;
	MPI *N = stackPop(s);
	if (N->S<0) {
		printf("n < 0\n");
		FREEMPI(N);
		return(ZEROI());
	}
	else {
		MPI *M = CFRACN(N);
		printf("the  continued fraction expansion of sqrt(2^");
		TEMP1 = MULT_I(N, 2);
		TEMP2 = ADD0_I(TEMP1, 1);
		PRINTI(TEMP2);
		FREEMPI(TEMP1);
		FREEMPI(TEMP2);
		printf(") has period-length ");
		PRINTI(M);
		FREEMPI(N);
		printf("\n");
		return(M);
	}
}

MPI *CFRAC_PERIOD_W(Stack s)
{
	USI t;
	MPI *P, *X, *Y;
	MPI *D = stackPop(s);
	X = BIG_MTHROOT(D, 2);
	Y = MULTI(X, X);
	t = EQUALI(D, Y);
	FREEMPI(X);
	FREEMPI(Y);
	if (t) {
		printf("D is a square\n");
		FREEMPI(D);
		return(ZEROI());
	}
	P = CFRAC_PERIOD(D);
	FREEMPI(D);
	return(P);
}
void MIDPT_COUNT_W(Stack s)
{
	MPI *Mptr = stackPop(s);
	MPI *Nptr = stackPop(s);
	MIDPT_COUNT(Mptr, Nptr);
	FREEMPI(Mptr);
	FREEMPI(Nptr);
	return;
}

MPI *NSCF_PERIOD_W(Stack s)
{
	USI t;
	MPI *P, *X, *Y, *ONE;
	MPI *D = stackPop(s);
	X = BIG_MTHROOT(D, 2);
	Y = MULTI(X, X);
	t = EQUALI(D, Y);
	FREEMPI(X);
	FREEMPI(Y);
	if (t) {
		printf("D is a square\n");
		FREEMPI(D);
		return(ZEROI());
	}
	ONE = ONEI();
	P = NSCF_PERIOD(D, ONE);
	FREEMPI(ONE);
	FREEMPI(D);
	return(P);
}

void CFRAC_COUNT_W(Stack s)
{
	MPI *Mptr = stackPop(s);
	MPI *Nptr = stackPop(s);
	CFRAC_COUNT(Mptr, Nptr);
	FREEMPI(Mptr);
	FREEMPI(Nptr);
	return;
}

/*void SPIRAL(MPI *P, MPI **Xptr, MPI **Yptr)*/
void SPIRAL_W(Stack s)
{
	MPI *P = stackPop(s);
	MPI **Xptr = stackPop(s);
	MPI **Yptr = stackPop(s);
	if (P->S < 0) {
		printf("Argument 1 is negative\n");
		FREEMPI(P);
		return;
	}
	SPIRAL(P, Xptr, Yptr);
	FREEMPI(P);
	/* Notice I did not free Xptr or Ypr. They are variables */
}

MPI *SPIRAL_INVERSE_W(Stack s)
{
	/* to free or not to free.  That is the question.  I say not! */
	MPI *Aptr = stackPop(s);
	MPI *Bptr = stackPop(s);
	MPI *Result = INVERSE_SPIRAL1(Aptr, Bptr);
	FREEMPI(Aptr);
	FREEMPI(Bptr);
	return Result;
}

MPI *NICF_PERIOD_W(Stack s)
{
	USI t;
	MPI *P, *X, *Y;
	MPI *D = stackPop(s);
	X = BIG_MTHROOT(D, 2);
	Y = MULTI(X, X);
	t = EQUALI(D, Y);
	FREEMPI(X);
	FREEMPI(Y);
	if (t) {
		printf("D is a square\n");
		FREEMPI(D);
		return(ZEROI());
	}
	P = NICF_PERIOD(D);
	FREEMPI(D);
	return(P);
}

MPI *NSCF_PERIOD0_W(Stack s)
{
	USI t;
	MPI *P, *X, *Y;
	MPI *D = stackPop(s);
	MPI *E = stackPop(s);
	MPI **Xptr = stackPop(s);
	MPI **Yptr = stackPop(s);
	X = BIG_MTHROOT(D, 2);
	Y = MULTI(X, X);
	t = EQUALI(D, Y);
	FREEMPI(X);
	FREEMPI(Y);
	if (t) {
		printf("D is a square\n");
		FREEMPI(D);
		FREEMPI(E);
		return(ZEROI());
	}
	P = NSCF_PERIOD0(D, E, Xptr, Yptr);
	FREEMPI(D);
	FREEMPI(E);
	return(P);
}

MPI *RCF_PERIOD0_W(Stack s)
{
	USI t;
	MPI *P, *X, *Y;
	MPI *D = stackPop(s);
	MPI *E = stackPop(s);
	MPI **Xptr = stackPop(s);
	MPI **Yptr = stackPop(s);
	X = BIG_MTHROOT(D, 2);
	Y = MULTI(X, X);
	t = EQUALI(D, Y);
	FREEMPI(X);
	FREEMPI(Y);
	if (t) {
		printf("D is a square\n");
		FREEMPI(D);
		FREEMPI(E);
		return(ZEROI());
	}
	P = RCF_PERIOD0(D, E, Xptr, Yptr);
	FREEMPI(D);
	FREEMPI(E);
	return(P);
}

MPI *NICF_PERIOD0_W(Stack s)
{
	USI t;
	MPI *P, *X, *Y;
	MPI *D = stackPop(s);
	MPI *E = stackPop(s);
	MPI **Xptr = stackPop(s);
	MPI **Yptr = stackPop(s);
	X = BIG_MTHROOT(D, 2);
	Y = MULTI(X, X);
	t = EQUALI(D, Y);
	FREEMPI(X);
	FREEMPI(Y);
	if (t) {
		printf("D is a square\n");
		FREEMPI(D);
		FREEMPI(E);
		return(ZEROI());
	}
	P = NICF_PERIOD0(D, E, Xptr, Yptr);
	FREEMPI(D);
	FREEMPI(E);
	return(P);
}

void TIME_COUNT_RCF_W(Stack s)
{
	MPI *Mptr = stackPop(s);
	MPI *Nptr = stackPop(s);
	if (RSV(Mptr, Nptr)> 0) {
		printf("Mptr > Nptr\n");
		rettype = FUNC_FAIL;
	}
	TIME_COUNT_RCF(Mptr, Nptr);
	FREEMPI(Mptr);
	FREEMPI(Nptr);
	return;
}

void TIME_COUNT_NICF_W(Stack s)
{
	MPI *Mptr = stackPop(s);
	MPI *Nptr = stackPop(s);
	if (RSV(Mptr, Nptr)> 0) {
		printf("Mptr > Nptr\n");
		rettype = FUNC_FAIL;
	}
	TIME_COUNT_NICF(Mptr, Nptr);
	FREEMPI(Mptr);
	FREEMPI(Nptr);
	return;
}

void TIME_COUNT_NSCF_W(Stack s)
{
	MPI *Mptr = stackPop(s);
	MPI *Nptr = stackPop(s);
	if (RSV(Mptr, Nptr)> 0) {
		printf("Mptr > Nptr\n");
		rettype = FUNC_FAIL;
	}
	TIME_COUNT_NSCF(Mptr, Nptr);
	FREEMPI(Mptr);
	FREEMPI(Nptr);
	return;
}

void PERIOD_TIME_COUNT_NSCF_W(Stack s)
{
	MPI *Mptr = stackPop(s);
	MPI *Nptr = stackPop(s);
	if (Mptr->D >= 2 || Nptr->D >= 2) {
		printf("Mptr->D >=2 or Nptr->D >= 2\n");
		rettype = FUNC_FAIL;
	}
	if (RSV(Mptr, Nptr)> 0) {
		printf("Mptr > Nptr\n");
		rettype = FUNC_FAIL;
	}
	PERIOD_TIME_COUNT_NSCF(Mptr, Nptr);
	FREEMPI(Mptr);
	FREEMPI(Nptr);
	return;
}

void PERIOD_TIME_COUNT_NICF_W(Stack s)
{
	MPI *Mptr = stackPop(s);
	MPI *Nptr = stackPop(s);
	if (Mptr->D >= 2 || Nptr->D >= 2) {
		printf("Mptr->D >=2 or Nptr->D >= 2\n");
		rettype = FUNC_FAIL;
	}
	if (RSV(Mptr, Nptr)> 0) {
		printf("Mptr > Nptr\n");
		rettype = FUNC_FAIL;
	}
	PERIOD_TIME_COUNT_NICF(Mptr, Nptr);
	FREEMPI(Mptr);
	FREEMPI(Nptr);
	return;
}

void PERIOD_TIME_COUNT_RCF_W(Stack s)
{
	MPI *Mptr = stackPop(s);
	MPI *Nptr = stackPop(s);
	if (Mptr->D >= 2 || Nptr->D >= 2) {
		printf("Mptr->D >=2 or Nptr->D >= 2\n");
		rettype = FUNC_FAIL;
	}
	if (RSV(Mptr, Nptr)> 0) {
		printf("Mptr > Nptr\n");
		rettype = FUNC_FAIL;
	}

	PERIOD_TIME_COUNT_RCF(Mptr, Nptr);
	FREEMPI(Mptr);
	FREEMPI(Nptr);
	return;
}

void TIME_COUNT_RCF00_W(Stack s)
{
	MPI *Mptr = stackPop(s);
	MPI *Nptr = stackPop(s);
	if (RSV(Mptr, Nptr)> 0) {
		printf("Mptr > Nptr\n");
		rettype = FUNC_FAIL;
	}
	TIME_COUNT_RCF00(Mptr, Nptr);
	FREEMPI(Mptr);
	FREEMPI(Nptr);
	return;
}

void TIME_COUNT_NICF00_W(Stack s)
{
	MPI *Mptr = stackPop(s);
	MPI *Nptr = stackPop(s);
	if (RSV(Mptr, Nptr)> 0) {
		printf("Mptr > Nptr\n");
		rettype = FUNC_FAIL;
	}
	TIME_COUNT_NICF00(Mptr, Nptr);
	FREEMPI(Mptr);
	FREEMPI(Nptr);
	return;
}

void TIME_COUNT_NSCF00_W(Stack s)
{
	MPI *Mptr = stackPop(s);
	MPI *Nptr = stackPop(s);
	if (RSV(Mptr, Nptr)> 0) {
		printf("Mptr > Nptr\n");
		rettype = FUNC_FAIL;
	}
	if (Nptr->D > 1) {
		printf("Mptr > 2^32\n");
		rettype = FUNC_FAIL;
	}
	TIME_COUNT_NSCF00(Mptr, Nptr);
	FREEMPI(Mptr);
	FREEMPI(Nptr);
	return;
}

void TIME_COUNT_RCF000_W(Stack s)
{
	MPI *Mptr = stackPop(s);
	MPI *Nptr = stackPop(s);
	if (RSV(Mptr, Nptr)> 0) {
		printf("Mptr > Nptr\n");
		rettype = FUNC_FAIL;
	}
	TIME_COUNT_RCF000(Mptr, Nptr);
	FREEMPI(Mptr);
	FREEMPI(Nptr);
	return;
}

void TIME_COUNT_NSCF000_W(Stack s)
{
	MPI *Mptr = stackPop(s);
	MPI *Nptr = stackPop(s);
	if (RSV(Mptr, Nptr)> 0) {
		printf("Mptr > Nptr\n");
		rettype = FUNC_FAIL;
	}
	TIME_COUNT_NSCF000(Mptr, Nptr);
	FREEMPI(Mptr);
	FREEMPI(Nptr);
	return;
}

void TIME_COUNT_NICF000_W(Stack s)
{
	MPI *Mptr = stackPop(s);
	MPI *Nptr = stackPop(s);
	if (RSV(Mptr, Nptr)> 0) {
		printf("Mptr > Nptr\n");
		rettype = FUNC_FAIL;
	}
	TIME_COUNT_NICF000(Mptr, Nptr);
	FREEMPI(Mptr);
	FREEMPI(Nptr);
	return;
}

void TIME_COUNT_RCF00_D_W(Stack s)
{
	MPI *Mptr = stackPop(s);
	MPI *Nptr = stackPop(s);
	if (RSV(Mptr, Nptr)> 0) {
		printf("Mptr > Nptr\n");
		rettype = FUNC_FAIL;
	}
	TIME_COUNT_RCF_D(Mptr, Nptr);
	FREEMPI(Mptr);
	FREEMPI(Nptr);
	return;
}

void TIME_COUNT_NSCF00_D_W(Stack s)
{
	MPI *Mptr = stackPop(s);
	MPI *Nptr = stackPop(s);
	if (RSV(Mptr, Nptr)> 0) {
		printf("Mptr > Nptr\n");
		rettype = FUNC_FAIL;
	}
	TIME_COUNT_NSCF_D(Mptr, Nptr);
	FREEMPI(Mptr);
	FREEMPI(Nptr);
	return;
}

void TIME_COUNT_NICF00_D_W(Stack s)
{
	MPI *Mptr = stackPop(s);
	MPI *Nptr = stackPop(s);
	if (RSV(Mptr, Nptr)> 0) {
		printf("Mptr > Nptr\n");
		rettype = FUNC_FAIL;
	}
	TIME_COUNT_NICF_D(Mptr, Nptr);
	FREEMPI(Mptr);
	FREEMPI(Nptr);
	return;
}

MPI *CFRAC_PERIOD1_W(Stack s)
{
	USI t;
	MPI *P, *X, *Y;
	MPI *D = stackPop(s);
	X = BIG_MTHROOT(D, 2);
	Y = MULTI(X, X);
	t = EQUALI(D, Y);
	FREEMPI(X);
	FREEMPI(Y);
	if (t) {
		printf("D is a square\n");
		FREEMPI(D);
		return(NULL);
	}
	P = CFRAC_PERIOD1(D);
	FREEMPI(D);
	return(P);
}

void PERIOD_TIME_COUNT_RCF1_W(Stack s)
{
	MPI *Mptr = stackPop(s);
	MPI *Nptr = stackPop(s);
	if (RSV(Mptr, Nptr)> 0) {
		printf("Mptr > Nptr\n");
		rettype = FUNC_FAIL;
	}
	PERIOD_TIME_COUNT_RCF1(Mptr, Nptr);
	FREEMPI(Mptr);
	FREEMPI(Nptr);
	return;
}

void PERIOD_TIME_COUNT_NSCF1_W(Stack s)
{
	MPI *Mptr = stackPop(s);
	MPI *Nptr = stackPop(s);
	if (RSV(Mptr, Nptr)> 0) {
		printf("Mptr > Nptr\n");
		rettype = FUNC_FAIL;
	}
	PERIOD_TIME_COUNT_NSCF1(Mptr, Nptr);
	FREEMPI(Mptr);
	FREEMPI(Nptr);
	return;
}

void PERIOD_TIME_COUNT_NICF1_W(Stack s)
{
	MPI *Mptr = stackPop(s);
	MPI *Nptr = stackPop(s);
	if (RSV(Mptr, Nptr)> 0) {
		printf("Mptr > Nptr\n");
		rettype = FUNC_FAIL;
	}
	PERIOD_TIME_COUNT_NICF1(Mptr, Nptr);
	FREEMPI(Mptr);
	FREEMPI(Nptr);
	return;
}

MPI *NSCF_PERIOD1_W(Stack s)
{
	USI t;
	MPI *P, *X, *Y;
	MPI *D = stackPop(s);
	X = BIG_MTHROOT(D, 2);
	Y = MULTI(X, X);
	t = EQUALI(D, Y);
	FREEMPI(X);
	FREEMPI(Y);
	if (t) {
		printf("D is a square\n");
		FREEMPI(D);
		return(NULL);
	}
	P = NSCF_PERIOD1(D);
	FREEMPI(D);
	return(P);
}

MPI *NICF_PERIOD1_W(Stack s)
{
	USI t;
	MPI *P, *X, *Y;
	MPI *D = stackPop(s);
	X = BIG_MTHROOT(D, 2);
	Y = MULTI(X, X);
	t = EQUALI(D, Y);
	FREEMPI(X);
	FREEMPI(Y);
	if (t) {
		printf("D is a square\n");
		FREEMPI(D);
		return(NULL);
	}
	P = NICF_PERIOD1(D);
	FREEMPI(D);
	return(P);
}

void TIME_COUNT_RCF_QIMPROVED_W(Stack s)
{
	MPI *Mptr = stackPop(s);
	MPI *Nptr = stackPop(s);
	if (RSV(Mptr, Nptr)> 0) {
		printf("Mptr > Nptr\n");
		rettype = FUNC_FAIL;
	}
	TIME_COUNT_RCF_QIMPROVED(Mptr, Nptr);
	FREEMPI(Mptr);
	FREEMPI(Nptr);
	return;
}

void TIME_COUNT_NICF_QIMPROVED_W(Stack s)
{
	MPI *Mptr = stackPop(s);
	MPI *Nptr = stackPop(s);
	if (RSV(Mptr, Nptr)> 0) {
		printf("Mptr > Nptr\n");
		rettype = FUNC_FAIL;
	}
	TIME_COUNT_NICF_QIMPROVED(Mptr, Nptr);
	FREEMPI(Mptr);
	FREEMPI(Nptr);
	return;
}

MPI *DIVISOR_PMINUS1_W(Stack s)
{
	MPI *N = stackPop(s);
	MPIA *PRIME_FOUND = stackPop(s);
	MPI *M;
	if (N->S <= 0) {
		printf("N <= 0\n");
		rettype = FUNC_FAIL;
	}
	if ((N->V[0]) % 2) {
		printf("N is odd\n");
		rettype = FUNC_FAIL;
	}
	M = DIVISOR_PMINUS1(N, PRIME_FOUND);
	FREEMPI(N);
	return(M);
}

void SORT_ARRAY_MPI_W(Stack s) {
	MPIA A = stackPop(s);
	MPIA *SORTED_ARRAY = stackPop(s);
	MPI *N = stackPop(s);
	SORT_ARRAY_MPI(A, SORTED_ARRAY, N);
	FREEMPI(N);
	return;
}

void CARMICHAEL_W(Stack s) {
	MPI *N = stackPop(s);
	CARMICHAEL(N);
	FREEMPI(N);
	return;
}

MPI *TANGENT_W(Stack s) {
	MPI *N = stackPop(s);
	MPI *I;
	I = TANGENTX(N);
	if (I == NULL) {
		rettype = FUNC_FAIL;
	}
	FREEMPI(N);
	return(I);
}

void BERNOULLI_W(Stack s) {
	MPI *N = stackPop(s);
	MPI **BERNOULLI_NUMERATOR = stackPop(s);
	MPI **BERNOULLI_DENOMONINATOR = stackPop(s);
	MPI *I;
	I = BERNOULLIX(N, BERNOULLI_NUMERATOR, BERNOULLI_DENOMONINATOR);
	if (I == NULL) {
		rettype = FUNC_FAIL;
	}
	else {
		FREEMPI(I);
	}
	FREEMPI(N);
	return;
}

MPI *PARTITION_W(Stack s) {
	MPI *N = stackPop(s);
	MPI *I;
	I = PARTITIONX(N);
	if (I == NULL) {
		rettype = FUNC_FAIL;
	}
	FREEMPI(N);
	return(I);
}

void EXCEPTIONALS_W(Stack s) {
	MPI *N = stackPop(s);
	EXCEPTIONALS(N);
	FREEMPI(N);
	return;
}

void EXCEPTIONALS0_W(Stack s) {
	MPI *N = stackPop(s);
	EXCEPTIONALS0(N);
	FREEMPI(N);
	return;
}

POLYI FGPI_W(Stack s)
{
	POLYI P = stackPop(s);
	POLYI Q = stackPop(s);
	POLYI Result = FGPI(P, Q);
	DELETEPI(P);
	DELETEPI(Q);
	return Result;
}

void EXCEPTIONALS10_W(Stack s) {
	MPI *N = stackPop(s);
	EXCEPTIONALS10(N);
	FREEMPI(N);
	return;
}

void GPLUS_W(Stack s) {
	POLYI k = stackPop(s);
	POLYI x = stackPop(s);
	POLYI y = stackPop(s);
	GPLUS(k, x, y);
	DELETEPI(k);
	DELETEPI(x);
	DELETEPI(y);
}

void GZERO_W(Stack s) {
	POLYI k = stackPop(s);
	POLYI x = stackPop(s);
	POLYI y = stackPop(s);
	GZERO(k, x, y);
	DELETEPI(k);
	DELETEPI(x);
	DELETEPI(y);
}

void GMINUS_W(Stack s) {
	POLYI k = stackPop(s);
	POLYI x = stackPop(s);
	POLYI y = stackPop(s);
	GMINUS(k, x, y);
	DELETEPI(k);
	DELETEPI(x);
	DELETEPI(y);
}

void NAGELL_W(Stack s) {
	MPI *D = stackPop(s);
	MPI *N = stackPop(s);
	MPI *Result = NAGELLX(D, N);
	if (Result == NULL)
		rettype = FUNC_FAIL;
	FREEMPI(Result);
	FREEMPI(D);
	FREEMPI(N);
	return;
}

void FRATTINI_W(Stack s) {
	MPI *D = stackPop(s);
	MPI *N = stackPop(s);
	MPI *Result = FRATTINIX(D, N);
	if (Result == NULL)
		rettype = FUNC_FAIL;
	FREEMPI(Result);
	FREEMPI(D);
	FREEMPI(N);
	return;
}

void KASHIHARA_W(Stack s) {
	MPI *M = stackPop(s);
	MPI *N = stackPop(s);
	MPI *Result = KASHIHARAX(M, N);
	if (Result == NULL)
		rettype = FUNC_FAIL;
	FREEMPI(Result);
	FREEMPI(M);
	FREEMPI(N);
	return;
}

void BRANCH_W(Stack s) {
	MPI *N = stackPop(s);
	MPI *TYPE = stackPop(s);
	MPI *FLAG = stackPop(s);
	MPI *Result = BRANCHX(N, TYPE, FLAG);
	if (Result == NULL)
		rettype = FUNC_FAIL;
	FREEMPI(Result);
	FREEMPI(N);
	FREEMPI(TYPE);
	FREEMPI(FLAG);
	return;
}
MPI *TERNARY_W(Stack s)
{
	MPI *M = stackPop(s);
	MPIA *Y = stackPop(s);
	MPI *Result = TERNARYX(M, Y);
	FREEMPI(M);
	if (Result == NULL)
		rettype = FUNC_FAIL;
	return(Result);
}
/* void  PEL4X(MPI *D, MPI *E, MPI **Xptr, MPI **Yptr) */
MPI *PEL4X_W(Stack s)
{
	MPI *D = stackPop(s);
	MPI *E = stackPop(s);
	MPI **Xptr = stackPop(s);
	MPI **Yptr = stackPop(s);
	MPI *Result = PEL4X(D, E, Xptr, Yptr);
	if (Result == NULL)
		rettype = FUNC_FAIL;
	FREEMPI(D);
	FREEMPI(E);
	return(Result);
}

MPI *STOLT_W(Stack s) {
	MPI *T;
	MPI *P = stackPop(s);
	MPI *Q = stackPop(s);
	MPI *R = stackPop(s);
	MPI *S = stackPop(s);
	T = STOLTX(P, Q, R, S);
	if (T == NULL)
		rettype = FUNC_FAIL;
	FREEMPI(P);
	FREEMPI(Q);
	FREEMPI(R);
	FREEMPI(S);
	return(T);
}

/* MPI *MIN_MPI_ARRAY(MPI *M[]) */
MPI *MIN_MPI_ARRAY_W(Stack s)
{
	MPIA M = stackPop(s);
	MPI *Result;
	Result = MIN_MPI_ARRAY(M);
	return Result;
}

/* MPI *LPRIME_FACTORX(MPI *N) */
MPI *LPRIME_FACTORX_W(Stack s)
{
	MPI *N = stackPop(s);
	MPI *Result = LPRIME_FACTORX(N);
	if (Result == NULL)
		rettype = FUNC_FAIL;
	FREEMPI(N);
	return Result;
}

/* MPI *CONWAY_CYCLESX(MPI *N) */
MPI *CONWAY_CYCLESX_W(Stack s) {
	MPI *A = stackPop(s);
	MPI *B = stackPop(s);
	MPI *Result = CONWAY_CYCLESX(A, B);
	if (Result == NULL)
		rettype = FUNC_FAIL;
	FREEMPI(A);
	FREEMPI(B);
	return(Result);
}

/* MPI *CONWAY_CYCLES0X(MPI *N) */
MPI *CONWAY_CYCLES0X_W(Stack s) {
	MPI *A = stackPop(s);
	MPI *B = stackPop(s);
	MPI *Result = CONWAY_CYCLES0X(A, B);
	if (Result == NULL)
		rettype = FUNC_FAIL;
	FREEMPI(A);
	FREEMPI(B);
	return(Result);
}

void CONWAY_RANGE_TEST_W(Stack s) {
	MPI *D = stackPop(s);
	MPI *N = stackPop(s);
	MPI *Result = CONWAY_RANGE_TESTX(D, N);
	if (Result == NULL)
		rettype = FUNC_FAIL;
	FREEMPI(Result);
	FREEMPI(D);
	FREEMPI(N);
	return;
}

void CONWAY_RANGE_TEST1_W(Stack s) {
	MPI *M1 = stackPop(s);
	MPI *M2 = stackPop(s);
	MPI *N1 = stackPop(s);
	MPI *N2 = stackPop(s);
	MPI *Result = CONWAY_RANGE_TEST1X(M1, M2, N1, N2);
	if (Result == NULL)
		rettype = FUNC_FAIL;
	FREEMPI(Result);
	FREEMPI(M1);
	FREEMPI(M2);
	FREEMPI(N1);
	FREEMPI(N2);
	return;
}
