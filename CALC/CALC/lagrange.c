#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_NONSTDC_NO_DEPRECATE
#include "integer.h"
#include <stdio.h>
#include <stdlib.h>
#include "fun.h"
#include "stack.h"
/*#include "unistd.h"*/
#include <string.h>

extern USI sqroot2_number;
void PATZ(MPI *D, MPI* N)
/* This solves the diophantine equations x^2-Dy^2=N and -N.
* If patz_verbose=1, all solutions in the range to be tested are printed.
* If patz_verbose=0, only the fundamental solutions are printed.
*/
{
	MPI *Q0, *S, *M, *M0, *TMP1, *TMP2, *X, *QQ, *ONE;
	MPI *X1, *X2, *Y1, *Y2, *Y, *TMP3, *TMP4, *Q0_2;
	MPIA A, AA, U, V, P, Q, FSX1, FSX2, FSY1, FSY2;
	USI l, period_length, k, r, s, flag1, flag2, kk, tmp;
	USI FLAG1, FLAG2, patz_verbose = 1, n1 = 0, n2 = 0;
	int t0, t, tt, ttt = 0 /*(to get rid of compiler complaint) */;
	USL i, j, m0;
	FILE *outfile;
	char buff[20];


	FSX1 = BUILDMPIA();
	FSX2 = BUILDMPIA();
	FSY1 = BUILDMPIA();
	FSY2 = BUILDMPIA();
	Q0 = ABSI(N);
	Q0_2 = INT0_(Q0, 2);
	tt = N->S;
	/* First solve x^2=D (mod Q_0), where Q_0=|N|. */
	S = SQROOT(D, Q0, &A, &M, &l);
	if (EQMINUSONEI(S)) { /* cleanup of nettbytes */
		FREEMPI(S);
		FREEMPIA(A);
		FREEMPI(M);
		FREEMPI(N);
		FREEMPI(Q0_2);
		FREEMPIA(FSX1);
		FREEMPIA(FSY1);
		FREEMPIA(FSX2);
		FREEMPIA(FSY2);
		printf("x^2=");
		PRINTI(D);
		FREEMPI(D);
		printf("(mod ");
		PRINTI(Q0);
		FREEMPI(Q0);
		printf(") not soluble;\n");
		execerror("Hence no primitive solutions", "");
	}
	FREEMPI(S);
	M0 = INT0(Q0, M);
	/*	FREEMPI(M); *//* bugfix 4th July 2000 */
	m0 = CONVERTI(M0);
	FREEMPI(M0);
	if (l * m0 > 100) {
		FREEMPI(Q0);
		FREEMPIA(A);
		FREEMPI(D);
		FREEMPI(N);
		FREEMPI(Q0_2);
		FREEMPIA(FSX1);
		FREEMPIA(FSY1);
		FREEMPIA(FSX2);
		FREEMPIA(FSY2);
		execerror("number of positive solutions of X^2=D(mod Q0) > 100", "");
	}
	ONE = ONEI();
	flag1 = 0;
	flag2 = 0;
	if (EQONEI(Q0))
		kk = 1;
	else
		kk = 0;
	strcpy(buff, "patz.out");
	outfile = fopen(buff, "w");
	for (i = 0; i < l; i++) {
		for (j = 0; j < m0; j++) {
			TMP2 = MULT_I(M, j); /* bugfix: formerly had M0 here */
			QQ = ADD0I(A->A[i], TMP2);
			FREEMPI(TMP2);
			/* inserted on 27th March 2002 */
			TMP2 = MULT_I(QQ, 2);
			t = RSV(TMP2, Q0);
			FREEMPI(TMP2);
			TMP2 = MULT_I(A->A[i], 2);/* inserted on 28th March 2002 */
			t0 = RSV(TMP2, M);
			FREEMPI(TMP2);
			if (t == 1) {
				if (t0) {
					TMP2 = QQ;
					QQ = SUB0I(Q0, QQ);
					FREEMPI(TMP2);
				}
				else {
					FREEMPI(QQ);
					continue;
				}
			}
			X1 = ZEROI();
			Y1 = ZEROI();
			X2 = ZEROI();
			Y2 = ZEROI();
			FLAG1 = 0;
			FLAG2 = 0;
			for (t = 1; t >= -1; t = t - 2) {
				/* finding the rcf of (QQ+sqrt(D))/Q0 and (-QQ+sqrt(D))/Q0 */
				if (t == -1)
					QQ->S = -(QQ->S);
				if (QQ->S == 0) /* (0+sqrt(D))/Q0= (-0+sqrt(D))/Q0 */
					t = -2;
				tmp = EQUALI(Q0_2, QQ);
				if (((Q0)->V[0]) % 2 == 0 && tmp)
					t = -2;
				/* (Q0/2+sqrt(D))/Q0 and (-Q0/2+sqrt(D))/Q0  give the same fundamental
				solutions */

				printf("rcf:(");
				fprintf(outfile, "rcf:(");
				PRINTI(QQ);
				FPRINTI(outfile, QQ);
				printf("+sqrt(");
				fprintf(outfile, "+sqrt(");
				PRINTI(D);
				FPRINTI(outfile, D);
				printf("))/");
				fprintf(outfile, "))/");
				PRINTI(Q0);
				FPRINTI(outfile, Q0);
				printf("\n");
				fprintf(outfile, "\n");
				period_length = SURD(D, ONE, QQ, Q0, &AA, &U, &V, &P, &Q, 1);
				r = AA->size;
				if (period_length % 2)
					period_length = 2 * period_length;
				s = r - period_length;
				for (k = kk; k < r; k++) {
					if ((V->A[k])->D == 0 && (V->A[k])->V[0] == 1) {
						TMP1 = MULTI(Q0, P->A[k - 1]);
						TMP2 = MULTI(QQ, Q->A[k - 1]);
						X = SUBI(TMP1, TMP2);
						FREEMPI(TMP1);
						FREEMPI(TMP2);
						if (patz_verbose) {
							printf("(X,Y) = (");
							fprintf(outfile, "(X,Y) = (");
							PRINTI(X);
							FPRINTI(outfile, X);
							printf(", ");
							fprintf(outfile, ", ");
							PRINTI(Q->A[k - 1]);
							FPRINTI(outfile, Q->A[k - 1]);
						}
						Y = COPYI(Q->A[k - 1]);
						if (patz_verbose) {
							printf("); X^2-");
							fprintf(outfile, "); X^2-");
							PRINTI(D);
							FPRINTI(outfile, D);
							printf("*Y^2=");
							fprintf(outfile, "*Y^2=");
						}
						if (k % 2)
							ttt = -((V->A[k])->S);
						else
							ttt = (V->A[k])->S;
						if (ttt == -1) {
							if (patz_verbose) {
								printf("-");
								fprintf(outfile, "-");
							}
							if (FLAG1) {
								if (RSV(Y, Y1) < 0) {
									FREEMPI(X1);
									FREEMPI(Y1);
									X1 = X;
									Y1 = Y;
									ADD_TO_MPIA(FSX1, X1, n1 - 1);
									ADD_TO_MPIA(FSY1, Y1, n1 - 1);
								}
								else {
									FREEMPI(X);
									FREEMPI(Y);
								}
							}
							else {
								FLAG1 = 1;
								FREEMPI(X1);
								FREEMPI(Y1);
								X1 = X;
								Y1 = Y;
								ADD_TO_MPIA(FSX1, X1, n1);
								ADD_TO_MPIA(FSY1, Y1, n1);
								n1++;
							}
						}
						else {
							if (FLAG2) {
								if (RSV(Y, Y2) < 0 && FLAG2) {
									FREEMPI(X2);
									FREEMPI(Y2);
									X2 = X;
									Y2 = Y;
									ADD_TO_MPIA(FSX2, X2, n2 - 1);
									ADD_TO_MPIA(FSY2, Y2, n2 - 1);
								}
								else {
									FREEMPI(X);
									FREEMPI(Y);
								}
							}
							else {
								FLAG2 = 1;
								FREEMPI(X2);
								FREEMPI(Y2);
								X2 = X;
								Y2 = Y;
								ADD_TO_MPIA(FSX2, X2, n2);
								ADD_TO_MPIA(FSY2, Y2, n2);
								n2++;
							}
						}
						if (patz_verbose) {
							PRINTI(Q0);
							FPRINTI(outfile, Q0);
							printf("\n");
							fprintf(outfile, "\n");
							GetReturn();
						}
						if (k >= s) {
							if (period_length % 2)
							{
								flag1 = 1;
								flag2 = 1;
							}
							else {
								if (k % 2)
									flag1 = 1;
								if ((k % 2) == 0)
									flag2 = 1;
							}
						}
					}
				}
				FREEMPIA(AA);
				FREEMPIA(U);
				FREEMPIA(V);
				FREEMPIA(P);
				FREEMPIA(Q);
			}
			if (Y1->S) {
				printf("Fundamental solution for x^2-");
				fprintf(outfile, "Fundamental solution for x^2-");
				PRINTI(D);
				FPRINTI(outfile, D);
				printf("y^2=-");
				fprintf(outfile, "y^2=-");
				PRINTI(Q0);
				FPRINTI(outfile, Q0);
				printf(": (x,y)=(");
				fprintf(outfile, ": (x,y)=(");
				if (X1->S <0)
					X1->S = -(X1->S);
				PRINTI(X1);
				FPRINTI(outfile, X1);
				printf(",");
				fprintf(outfile, ",");
				PRINTI(Y1);
				FPRINTI(outfile, Y1);
				printf(")\n");
				fprintf(outfile, ")\n");
			}
			FREEMPI(X1);
			FREEMPI(Y1);
			if (Y2->S) {
				printf("Fundamental solution for x^2-");
				fprintf(outfile, "Fundamental solution for x^2-");
				PRINTI(D);
				FPRINTI(outfile, D);
				printf("y^2=");
				fprintf(outfile, "y^2=");
				PRINTI(Q0);
				FPRINTI(outfile, Q0);
				printf(": (x,y)=(");
				fprintf(outfile, ": (x,y)=(");
				if (X2->S <0)
					X2->S = -(X2->S);
				PRINTI(X2);
				FPRINTI(outfile, X2);
				printf(",");
				fprintf(outfile, ",");
				PRINTI(Y2);
				FPRINTI(outfile, Y2);
				printf(")\n");
				fprintf(outfile, ")\n");
			}
			FREEMPI(X2);
			FREEMPI(Y2);
			FREEMPI(QQ);
		}
	}
	FREEMPIA(A);
	FREEMPI(ONE);
	printf("--------------------------\n");
	fprintf(outfile, "--------------------------\n");
	printf("X^2-");
	fprintf(outfile, "X^2-");
	PRINTI(D);
	FPRINTI(outfile, D);
	printf("*Y^2=");
	fprintf(outfile, "*Y^2=");
	if (flag1) {
		printf("-");
		fprintf(outfile, "-");
		PRINTI(Q0);
		FPRINTI(outfile, Q0);
		printf(" is soluble\n");
		fprintf(outfile, " is soluble\n");
		for (i = 0; i < n1; i++) {
			printf("Fundamental solution (");
			fprintf(outfile, "Fundamental solution (");
			if ((FSX1->A)[i]->S <0)
				(FSX1->A)[i]->S = -((FSX1->A)[i]->S);
			PRINTI((FSX1->A)[i]);
			FPRINTI(outfile, (FSX1->A)[i]);
			printf(",");
			fprintf(outfile, ",");
			PRINTI((FSY1->A)[i]);
			FPRINTI(outfile, (FSY1->A)[i]);
			printf(")\n");
			fprintf(outfile, ")\n");
			TMP1 = GCD((FSX1->A)[i], (FSY1->A)[i]);
			if (EQONEI(TMP1) == 0) {
				printf("imprimitive solution!\n");
				fprintf(outfile, "imprimitive solution!\n");
			}
			FREEMPI(TMP1);
		}
	}
	else {
		printf("-");
		fprintf(outfile, "-");
		PRINTI(Q0);
		FPRINTI(outfile, Q0);
		printf(" is insoluble\n");
		fprintf(outfile, " is insoluble\n");
	}
	printf("--------------------------\n");
	fprintf(outfile, "--------------------------\n");
	printf("X^2-");
	fprintf(outfile, "X^2-");
	PRINTI(D);
	FPRINTI(outfile, D);
	printf("*Y^2=");
	fprintf(outfile, "*Y^2=");
	if (flag2) {
		PRINTI(Q0);
		FPRINTI(outfile, Q0);
		printf(" is soluble\n");
		fprintf(outfile, " is soluble\n");
		if (EQONEI(Q0) && ttt == -1) {
			/* calculating eta^2, when Norm(eta)= -1 */
			TMP1 = MULTI((FSX1->A)[0], (FSX1->A)[0]);
			TMP2 = MULTI((FSY1->A)[0], (FSY1->A)[0]);
			TMP3 = MULTI(TMP2, D);
			FREEMPI(TMP2);
			TMP4 = ADD0I(TMP1, TMP3);
			FREEMPI(TMP1);
			FREEMPI(TMP3);
			ADD_TO_MPIA(FSX2, TMP4, 0);
			FREEMPI(TMP4);
			TMP1 = MULTI((FSX1->A)[0], (FSY1->A)[0]);
			TMP2 = MULT_I(TMP1, 2);
			FREEMPI(TMP1);
			ADD_TO_MPIA(FSY2, TMP2, 0);
			FREEMPI(TMP2);
			n2 = 1; /* need to increment n2 from 0 */
		}
		for (i = 0; i < n2; i++) {
			printf("Fundamental solution (");
			fprintf(outfile, "Fundamental solution (");
			if ((FSX2->A)[i]->S <0)
				(FSX2->A)[i]->S = -((FSX2->A)[i]->S);
			PRINTI((FSX2->A)[i]);
			FPRINTI(outfile, (FSX2->A)[i]);
			printf(",");
			fprintf(outfile, ",");
			PRINTI((FSY2->A)[i]);
			FPRINTI(outfile, (FSY2->A)[i]);
			printf(")\n");
			fprintf(outfile, ")\n");
			TMP1 = GCD((FSX2->A)[i], (FSY2->A)[i]);
			if (EQONEI(TMP1) == 0) {
				printf("imprimitive solution!\n");
				fprintf(outfile, "imprimitive solution!\n");
			}
			FREEMPI(TMP1);
		}
	}
	else {
		PRINTI(Q0);
		FPRINTI(outfile, Q0);
		printf(" is insoluble\n");
		fprintf(outfile, " is insoluble\n");
	}
	FREEMPI(Q0);
	FREEMPI(Q0_2);
	FREEMPIA(FSX1);
	FREEMPIA(FSY1);
	FREEMPIA(FSX2);
	FREEMPIA(FSY2);
	FREEMPI(M); /* placed this here on 4th July 2000 */
	fclose(outfile);

	return;
}

MPI *PATZX(MPI *D, MPI *N)
{
	MPI *G, *X;
	unsigned long f;

	if (D->S <= 0)
	{
		printf("D <= 0\n");
		return NULL;
	}
	if (EQONEI(D))
	{
		printf("D = 1\n");
		return NULL;
	}
	X = BIG_MTHROOT(D, 2);
	G = MULTI(X, X);
	f = EQUALI(D, G);
	FREEMPI(G);
	FREEMPI(X);
	if (f)
	{
		printf("D is a perfect square\n");
		return NULL;
	}
	if (N->S == 0)
	{
		printf("N = 0\n");
		return NULL;
	}
	PATZ(D, N);
	return(ONEI());
}

MPI *QUADRATIC(MPI *A, MPI *B, MPI *C, MPI *N, MPIA *SOL)
/* Solving the congruence AX^2+BX+C=0 (mod N). */
/* returns -1 if no solution, otherwise returns the number of solutions.
* The actual solutions are returned as the array SOL.
* Finished 14th December 2000.
*/
{
	MPI *D, *T, *TMP1, *TMP2, *M, *MM, *Z, *M0, *FOUR;
	MPI *NN, *QQ, *BB, *NNN, *NUMBER, *ABSN, *HALFN;
	MPIA SOL1, SOL2;
	USL i, j, k, m0, s;
	USI l, flag = 0;
	int t, tt;

	TMP1 = MULTI(B, B);
	FOUR = CHANGEI(4);
	TMP2 = MULTI3(FOUR, A, C);
	FREEMPI(FOUR);
	ABSN = ABSI(N);
	NN = MULT_I(ABSN, 4);
	NNN = MULT_I(ABSN, 2);
	D = SUBI(TMP1, TMP2);
	FREEMPI(TMP1);
	FREEMPI(TMP2);
	BB = COPYI(B);
	if (((B->V[0]) % 2) == 0) {
		TMP1 = D;
		D = INT_(D, 4);
		FREEMPI(TMP1);
		TMP1 = BB;
		BB = INT_(B, 2);
		FREEMPI(TMP1);
		TMP1 = NN;
		NN = INT_(NN, 4);
		FREEMPI(TMP1);
	}

	T = SQROOT(D, NN, &SOL1, &M, &l);
	*SOL = BUILDMPIA();
	SOL2 = BUILDMPIA();
	k = 0;
	if (EQMINUSONEI(T)) {
		FREEMPI(D);
		FREEMPI(NN);
		FREEMPI(BB);
		FREEMPI(NNN);
		FREEMPIA(SOL1);
		FREEMPI(M);
		FREEMPI(T);
		FREEMPI(ABSN);
		ADD_TO_MPIA(*SOL, NULL, k);
		return(ZEROI());
	}
	FREEMPI(T);

	M0 = INT0(NN, M);
	m0 = CONVERTI(M0);
	FREEMPI(M0);

	s = N->V[0] % 2;
	HALFN = INT_(ABSN, 2);
	for (i = 0; i < l; i++) {
		for (j = 0; j < m0; j++) {
			TMP2 = MULT_I(M, j);
			QQ = ADD0I(SOL1->A[i], TMP2);
			FREEMPI(TMP2);
			if (((B->V[0]) % 2)) {
				t = RSV(QQ, NNN);
				if (t == 1) {
					TMP1 = QQ;
					QQ = SUBI(NN, QQ);
					FREEMPI(TMP1);
				}
				ADD_TO_MPIA(SOL2, QQ, k);
				k++;
			}
			else {
				ADD_TO_MPIA(SOL2, QQ, k);
				k++;
				if (s == 0) {
					tt = RSV(QQ, HALFN);
					if (tt == 0)
						flag = 1;
				}
				if (QQ->S && (flag == 0)) {
					TMP2 = MINUSI(QQ);
					ADD_TO_MPIA(SOL2, TMP2, k);
					k++;
					FREEMPI(TMP2);
				}
			}
			FREEMPI(QQ);
		}
	}
	for (i = 0; i < k; i++) {
		TMP1 = SUBI(SOL2->A[i], BB);
		if (((B->V[0]) % 2) == 0)
			TMP2 = COPYI(TMP1);
		else
			TMP2 = INT_(TMP1, 2);
		Z = CONGR(A, TMP2, ABSN, &MM);
		ADD_TO_MPIA(*SOL, Z, i);
		FREEMPI(Z);
		FREEMPI(TMP1);
		FREEMPI(TMP2);
		FREEMPI(MM);
	}
	FREEMPIA(SOL1);
	FREEMPIA(SOL2);
	FREEMPI(NNN);
	FREEMPI(NN);
	FREEMPI(ABSN);
	FREEMPI(HALFN);
	FREEMPI(BB);
	FREEMPI(M);
	FREEMPI(D);
	NUMBER = CHANGEI(k);
	return(NUMBER);
}

MPI *QUADRATICX(MPI *A, MPI *B, MPI *C, MPI *N, MPIA *SOL)
{
	MPI *T, *X;
	USI t;

	if (A->S == 0 || N->S <= 0) {
		*SOL = BUILDMPIA();
		ADD_TO_MPIA(*SOL, NULL, 0);
		printf("A = 0 or N <= 0\n");
		return(NULL);
	}
	X = GCD(A, N);
	t = EQONEI(X);
	FREEMPI(X);
	if (!t) {
		*SOL = BUILDMPIA();
		ADD_TO_MPIA(*SOL, NULL, 0);
		printf("gcd(a,n)>1\n");
		return(NULL);
	}
	T = QUADRATIC(A, B, C, N, SOL);
	return(T);
}

/* Below we implement the algorithm in my paper
* "The diophantine equation ax^2+bxy+cy^2=N, D=b^2-4ac >0".
* gcd(A,N)=1."
*/
void BINARYFORM(MPI *A, MPI *B, MPI *C, MPI *N, MPIA *FSX1, MPIA *FSY1, USI *N1, USI FLAG, USI verbose)
/* first solve the congruence a\theta^2+b\theta+c=0 (mod |N|)
* Let \Delta = B[1]^2-ac if B=2*B[1], else D = B^2-4AC and Q = A*|N|.
* Then let n= 2A\theta + B, P = int(n/2).
* If B is even,
*       let \omega=(-P +\sqrt{\Delta})/Q,\omega* =(-P -\sqrt{\Delta})/Q,
* If B is odd.
*       let \omega=(-n +\sqrt{D})/2Q,\omega* =(-n -\sqrt{D})/2Q,
* Let l = period length of omega and omega*.
* If B is even (or odd):
*                    test up to the first period of \omega for
* Q[k]=(-1)^{k}N/|N| (Q[k]=2(-1)^{k}N/|N|) if l is even,
* but up to the second period, if l is odd.
* Similarly for \omega* but with Q[k]=(-1)^(k+1)N/|N| (Q[k]=2(-1)^(k+1)N/|N|)
* There will be no solution if the test always fails. Otherwise
* let X/y=p[k-1]/q[k-1], be the corresponding convergent.
* Then x=y\theta+|N|X will be a solution of our diophantine equation
* for the class determined by n. Choose y to be minimal.
* If FLAG = 1, we print output, as GCD(A,N)=1 in BINARYFORM1.
* If FLAG = 0, we do not print output, as GCD(A,N)>1 in BINARYFORM1.
*/
{
	MPI *L, *TMP1, *TMP2, *FOUR, *D, *ABSN, *THETA, *TMP3, *TMP4, *MODULUS;
	MPI *Q, *n, *P, *P01, *P02, *Q01, *Q02, *MINUSQ, *MINUS2Q;
	MPI *ONE, *X, *TWOQ, *Y, *TWOABSN, *x, *y, *SUM;
	MPI *Z, *DELTA, *BETA, *TWOA, *TWOC, *TMP5, *TMP6;
	MPIA SOL, AA1, U1, V1, P1, Q1;
	MPIA  AA2, U2, V2, P2, Q2, FSX, FSY;
	USL l, v1, v2;
	USI s, t, period_length, k, r1, r2, r, d1, d2, flag = 0;
	USI i, FLAGE = 0, FLAG1 = 0, FLAG2 = 0;
	int tt, ttt, u, n1 = 0;
	FILE *outfile;
	char buff[20];

	TWOA = MULT_I(A, 2);
	TWOC = MULT_I(C, 2);

	/* first solve the congruence a\theta^2+b\theta+c=0 (mod |N|) */
	L = QUADRATIC(A, B, C, N, &SOL);
	if (L->S == 0) {
		FREEMPI(L);
		FREEMPIA(SOL);
		execerror("There are no primitive solutions", "");
	}
	FSX = BUILDMPIA();
	FSY = BUILDMPIA();

	strcpy(buff, "binform.out");
	outfile = fopen(buff, "a");
	/* Let \Delta = B[1]^2-ac if B=2*B[1], else D = B^2-4AC.     */
	TMP1 = MULTI(B, B);
	FOUR = CHANGEI(4);
	TMP2 = MULTI3(FOUR, A, C);
	FREEMPI(FOUR);
	D = SUBI(TMP1, TMP2);
	FREEMPI(TMP1);
	FREEMPI(TMP2);
	t = (B->V[0]) % 2;
	if (t) {
		printf("D=");
		fprintf(outfile, "D=");
		PRINTI(D);
		FPRINTI(outfile, D);
		printf("\n");
		fprintf(outfile, "\n");
		if (((D->D == 0) && (D->V[0] == 5)) && (((A->S) * (N->S)) == -1))
			flag = 1;
	}
	else {
		TMP1 = D;
		D = INT_(D, 4);
		FREEMPI(TMP1);
		printf("Delta=");
		fprintf(outfile, "Delta=");
		PRINTI(D);
		FPRINTI(outfile, D);
		printf("\n");
		fprintf(outfile, "\n");
	}
	ABSN = ABSI(N);
	TWOABSN = MULT_I(ABSN, 2);
	tt = N->S;
	Q = MULTI(A, ABSN);
	TWOQ = MULT_I(Q, 2);
	MINUS2Q = MINUSI(TWOQ);
	MINUSQ = MINUSI(Q);
	ONE = ONEI();
	l = CONVERTI(L);
	FREEMPI(L);

	for (r = 0; r < l; r++) {
		THETA = COPYI(SOL->A[r]);
		if (verbose) {
			printf("THETA[%u]=", r);
			fprintf(outfile, "THETA[%u]=", r);
			PRINTI(THETA);
			FPRINTI(outfile, THETA);
			printf("\n");
			fprintf(outfile, "\n");
		}
		TMP1 = MULT_I(THETA, 2);
		TMP2 = MULTI(TMP1, A);
		FREEMPI(TMP1);
		n = ADDI(TMP2, B);
		FREEMPI(TMP2);
		P = INT_(n, 2);
		if (verbose) {
			printf("n=");
			fprintf(outfile, "n=");
			PRINTI(n);
			FPRINTI(outfile, n);
			printf("\n");
			fprintf(outfile, "\n");
		}
		if (t == 0) {
			P01 = MINUSI(P);
			P02 = COPYI(P);
			Q01 = COPYI(Q);
			Q02 = COPYI(MINUSQ);
		}
		else {
			P01 = MINUSI(n);
			P02 = COPYI(n);
			Q01 = COPYI(TWOQ);
			Q02 = COPYI(MINUS2Q);
		}
		period_length = SURD(D, ONE, P01, Q01, &AA1, &U1, &V1, &P1, &Q1, 1);
		period_length = SURD(D, ONE, P02, Q02, &AA2, &U2, &V2, &P2, &Q2, 1);
		if (period_length % 2)
			period_length = 2 * period_length;
		r1 = AA1->size;
		r2 = AA2->size;
		if (flag) {
			s = r1 - period_length;
			if (s > 1) {
				TMP1 = SUBI(P1->A[s - 1], P1->A[s - 2]);
				TMP2 = SUBI(Q1->A[s - 1], Q1->A[s - 2]);
				TMP3 = MULTI(ABSN, TMP1);
				TMP4 = MULTI(THETA, TMP2);
				X = ADDI(TMP3, TMP4);
				FREEMPI(TMP1);
				FREEMPI(TMP3);
				FREEMPI(TMP4);
			}
			else {
				TMP1 = ONEI();
				X = SUBI(P1->A[0], TMP1);
				FREEMPI(TMP1);
				TMP2 = COPYI(Q1->A[0]);
			}
			/* print from here */
			if (verbose) {
				printf("Exceptional solution: s=%u, (X,Y) = (", s);
				fprintf(outfile, "Exceptional solution: s=%u, (X,Y) = (", s);
				PRINTI(X);
				FPRINTI(outfile, X);
				printf(", ");
				fprintf(outfile, ", ");
				PRINTI(TMP2);
				FPRINTI(outfile, TMP2);
				printf(")\n");
				fprintf(outfile, ")\n");
				GetReturn();
			}
			/* to from here */
			ADD_TO_MPIA(FSX, X, n1);
			ADD_TO_MPIA(FSY, TMP2, n1);
			FLAGE = 1;
			FREEMPI(X);
			FREEMPI(TMP2);
		}

		for (k = 1; k < r1; k++) {
			ttt = (k % 2) ? -1 : 1;
			d1 = (V1->A[k])->D;
			v1 = (V1->A[k])->V[0];
			/*printf("V1->A[%u]=", k); PRINTI(V1->A[k]);printf("\n"); */
			if ((!t && ((d1 == 0) && (v1 == 1))) || (t && ((d1 == 0) && (v1 == 2)))) {
				if (verbose) {
					printf("processing omega:\n");
					fprintf(outfile, "processing omega:\n");
					printf("rcf:(");
					fprintf(outfile, "rcf:(");
					PRINTI(P01);
					FPRINTI(outfile, P01);
					printf("+sqrt(");
					fprintf(outfile, "+sqrt(");
					PRINTI(D);
					FPRINTI(outfile, D);
					printf("))/");
					fprintf(outfile, "))/");
					PRINTI(Q01);
					FPRINTI(outfile, Q01);
					printf("\n");
					fprintf(outfile, "\n");
					GetReturn();
				}
				if (!flag && ((V1->A[k])->S == ttt * tt)) {/* calculate x=y\theta+|N|X */
					TMP1 = MULTI(ABSN, P1->A[k - 1]);
					TMP2 = MULTI(THETA, Q1->A[k - 1]);
					X = ADDI(TMP1, TMP2);
					FREEMPI(TMP1);
					FREEMPI(TMP2);
					Y = COPYI(Q1->A[k - 1]);
					/* print from here */
					if (verbose) {
						printf("(X1,Y1) = (");
						fprintf(outfile, "(X1,Y1) = (");
						PRINTI(X);
						FPRINTI(outfile, X);
						printf(", ");
						fprintf(outfile, ", ");
						PRINTI(Y);
						FPRINTI(outfile, Y);
						printf(")\n");
						fprintf(outfile, ")\n");
						GetReturn();
					}
					/* to from here */
					ADD_TO_MPIA(FSX, X, n1);
					ADD_TO_MPIA(FSY, Y, n1);
					FREEMPI(X);
					FREEMPI(Y);
					FLAG1 = 1;
					break;
				}
			}
		}
		for (k = 1; k < r2; k++) {
			ttt = (k % 2) ? 1 : -1;
			d2 = (V2->A[k])->D;
			v2 = (V2->A[k])->V[0];
			/*	printf("V2->A[%u]=", k); PRINTI(V2->A[k]);printf("\n");*/
			if ((!t && ((d2 == 0) && (v2 == 1))) || (t && ((d2 == 0) && (v2 == 2)))) {
				if (verbose) {
					printf("processing omega_star\n");
					fprintf(outfile, "processing omega_star\n");
					printf("rcf:(");
					fprintf(outfile, "rcf:(");
					PRINTI(P02);
					FPRINTI(outfile, P02);
					printf("+sqrt(");
					fprintf(outfile, "+sqrt(");
					PRINTI(D);
					FPRINTI(outfile, D);
					printf("))/");
					fprintf(outfile, "))/");
					PRINTI(Q02);
					FPRINTI(outfile, Q02);
					printf("\n");
					fprintf(outfile, "\n");
					GetReturn();
				}
				if ((V2->A[k])->S == ttt * tt) {/* calculate x=y\theta+|N|X */
					TMP1 = MULTI(ABSN, P2->A[k - 1]);
					TMP2 = MULTI(THETA, Q2->A[k - 1]);
					X = ADDI(TMP1, TMP2);
					FREEMPI(TMP1);
					FREEMPI(TMP2);
					Y = COPYI(Q2->A[k - 1]);
					/* print from here */
					if (verbose) {
						printf("(X2,Y2) = (");
						fprintf(outfile, "(X2,Y2) = (");
						PRINTI(X);
						FPRINTI(outfile, X);
						printf(", ");
						fprintf(outfile, ", ");
						PRINTI(Y);
						FPRINTI(outfile, Y);
						printf(")\n");
						fprintf(outfile, ")\n");
						GetReturn();
					}
					/* to from here */
					if (!FLAG1 && !FLAGE) {
						ADD_TO_MPIA(FSX, X, n1);
						ADD_TO_MPIA(FSY, Y, n1);
						FLAG2 = 1;
					}
					else {
						u = RSV(FSY->A[n1], Y);
						if (u > 0) {
							TMP1 = FSX->A[n1];
							TMP2 = FSY->A[n1];
							FSX->A[n1] = COPYI(X);
							FSY->A[n1] = COPYI(Y);
							FREEMPI(TMP1);
							FREEMPI(TMP2);
						}
					}
					FREEMPI(X);
					FREEMPI(Y);
					break;
				}
			}
		}
		if (FLAGE || FLAG1 || FLAG2) {
			n1++;
			FLAGE = 0;
			FLAG1 = 0;
			FLAG2 = 0;
		}
		FREEMPIA(AA1);
		FREEMPIA(AA2);
		FREEMPIA(U1);
		FREEMPIA(U2);
		FREEMPIA(V1);
		FREEMPIA(V2);
		FREEMPIA(P1);
		FREEMPIA(P2);
		FREEMPIA(Q1);
		FREEMPIA(Q2);
		FREEMPI(THETA);
		FREEMPI(P01);
		FREEMPI(P02);
		FREEMPI(Q01);
		FREEMPI(Q02);
		FREEMPI(P);
		FREEMPI(n);
	}
	/*TMP1 = MULTI(A, N);
	t = EQONEI(TMP1);
	FREEMPI(TMP1);
	if(t && (n1 == 0)){
	ADD_TO_MPIA(FSX, ONE, n1);
	TMP1 = ZEROI();
	ADD_TO_MPIA(FSY, TMP1, n1);
	FREEMPI(TMP1);
	n1++;
	}*/
	if (FLAG) {
		if (n1 == 0)
			printf("There are no primitive solutions (x,y)!\n");
		else {
			for (i = 0; i < n1; i++) {
				printf("solution (");
				fprintf(outfile, "solution (");
				PRINTI(FSX->A[i]);
				FPRINTI(outfile, FSX->A[i]);
				printf(",");
				fprintf(outfile, ",");
				PRINTI(FSY->A[i]);
				FPRINTI(outfile, FSY->A[i]);
				printf(")");
				fprintf(outfile, ")");
				x = FSX->A[i];
				y = FSY->A[i];
				Z = EUCLIDI(x, y, &DELTA, &BETA);
				FREEMPI(Z);
				BETA->S = -(BETA->S);
				/* x(DELTA)-y(BETA)=1 */
				TMP1 = MULTI3(TWOA, x, BETA);
				TMP2 = MULTI3(TWOC, y, DELTA);
				TMP5 = MULTI3(B, x, DELTA);
				TMP6 = MULTI3(B, BETA, y);
				/* now add TMP1, TMP2, TPM5, TMP6 */
				SUM = ADDI(TMP1, TMP2);
				FREEMPI(TMP1);
				FREEMPI(TMP2);
				FREEMPI(BETA);
				FREEMPI(DELTA);
				TMP1 = SUM;
				SUM = ADDI(SUM, TMP5);
				FREEMPI(TMP1);
				FREEMPI(TMP5);
				TMP1 = SUM;
				SUM = ADDI(SUM, TMP6);
				FREEMPI(TMP1);
				FREEMPI(TMP6);
				TMP1 = SUM;
				SUM = MOD(SUM, TWOABSN);
				FREEMPI(TMP1);
				if (RSV(SUM, ABSN) > 0) {
					TMP1 = SUM;
					SUM = SUBI(SUM, TWOABSN);
					FREEMPI(TMP1);
				}
				MODULUS = MULT_I(ABSN, 2);
				printf(": n = ");
				fprintf(outfile, ": n = ");
				PRINTI(SUM);
				FPRINTI(outfile, SUM);
				printf("(mod ");
				PRINTI(MODULUS);
				printf(")\n");
				fprintf(outfile, "(mod ");
				FPRINTI(outfile, MODULUS);
				fprintf(outfile, ")\n");
				FREEMPI(SUM);
				FREEMPI(MODULUS);
				GetReturn();
			}
		}
	}
	FREEMPIA(SOL);
	FREEMPI(D);
	FREEMPI(Q);
	FREEMPI(ABSN);
	FREEMPI(TWOABSN);
	FREEMPI(TWOQ);
	FREEMPI(TWOA);
	FREEMPI(TWOC);
	FREEMPI(MINUSQ);
	FREEMPI(MINUS2Q);
	FREEMPI(ONE);
	*N1 = n1;
	*FSX1 = FSX;
	*FSY1 = FSY;
	fclose(outfile);
	return;
}

void BINARYFORM1(MPI *AA, MPI *BB, MPI *CC, MPI *NN, MPI *VERB)
/* We reduce to the case GCD(AA,NN)=1 and then use BINARYFORM.
* The n of p.3 of paper http://www.maths.uq.edu.au/~krm/binary.pdf
* is also printed.
* Verbose output if VERB = 1; otherwise 0;
* Completed 21/12/00.
*/
{
	MPI *a, *b, *c, *x, *y, *TWOA, *TWOC, *MODULUS;
	MPI *G, *alpha, *beta, *M, *gamma, *delta, *Z;
	MPI *ABSN, *TWOABSN, *BETA, *DELTA, *X, *D, *FOUR;
	MPI *TMP1, *TMP2, *TMP3, *TMP4, *TMP5, *TMP6, *SUM;
	MPI *A, *B, *C, *N;
	USI t, length, i, FLAG = 0, verb;
	MPIA FSX1, FSY1;
	int s;
	FILE *outfile;
	char buff[20];
	strcpy(buff, "binform.out");
	outfile = fopen(buff, "w");

	/* Check that NN is non-zero */
	s = NN->S;
	if (s == 0)
	{
		printf("NN = 0\n");
		return;
	}
	/* Check to see d=gcd(AA,BB,CC) divides NN and if so
	* to replace AA,BB,CC,NN by AA/d,BB/d,CC/d,NN/d.
	*/
	TMP1 = GCD(AA, BB);
	TMP2 = GCD(TMP1, CC);
	TMP3 = MOD(NN, TMP2);
	s = TMP3->S;
	FREEMPI(TMP1);
	FREEMPI(TMP3);
	if (s) {
		printf("gcd(A,B,C) does not divide N\n");
		FREEMPI(TMP2);
		return;
	}
	else {
		A = INT(AA, TMP2);
		B = INT(BB, TMP2);
		C = INT(CC, TMP2);
		N = INT(NN, TMP2);
	}
	FREEMPI(TMP2);

	/* Check that B^2-4AC > 0 */
	TMP1 = MULTI(B, B);
	FOUR = CHANGEI(4);
	TMP2 = MULTI3(FOUR, A, C);
	FREEMPI(FOUR);
	D = SUBI(TMP1, TMP2);
	s = D->S;
	FREEMPI(TMP1);
	FREEMPI(TMP2);
	if (s <= 0)
	{
		printf("D <= 0\n");
		FREEMPI(D);
		FREEMPI(A);
		FREEMPI(B);
		FREEMPI(C);
		FREEMPI(N);
		return;
	}
	/* Check that B^2-4AC is not a perfect square */
	X = BIG_MTHROOT(D, 2);
	G = MULTI(X, X);
	t = EQUALI(D, G);
	FREEMPI(G);
	FREEMPI(X);
	FREEMPI(D);
	if (t)
	{
		printf("D is a perfect square\n");
		FREEMPI(A);
		FREEMPI(B);
		FREEMPI(C);
		FREEMPI(N);
		return;
	}

	TWOA = MULT_I(A, 2);
	TWOC = MULT_I(C, 2);
	ABSN = ABSI(N);
	TWOABSN = MULT_I(ABSN, 2);
	G = GCD(A, N);
	t = EQONEI(G);
	FREEMPI(G);
	if (t) {
		a = COPYI(A);
		b = COPYI(B);
		c = COPYI(C);
		FLAG = 1;
	}
	else {
		GAUSS(A, B, C, N, &alpha, &gamma, &M);
		Z = EUCLIDI(alpha, gamma, &delta, &beta);
		FREEMPI(Z);
		beta->S = -(beta->S);
		a = M;
		strcpy(buff, "binform.out");
		outfile = fopen(buff, "w");

		/* alpha * delta - beta * gamma = 1 */
		if (VERB->S) {
			printf("alpha="); PRINTI(alpha); printf("\n");
			fprintf(outfile, "alpha="); FPRINTI(outfile, alpha); fprintf(outfile, "\n");
			printf("beta="); PRINTI(beta); printf("\n");
			fprintf(outfile, "beta="); FPRINTI(outfile, beta); fprintf(outfile, "\n");
			printf("gamma="); PRINTI(gamma); printf("\n");
			fprintf(outfile, "gamma="); FPRINTI(outfile, gamma); fprintf(outfile, "\n");
			printf("delta="); PRINTI(delta); printf("\n");
			fprintf(outfile, "delta="); FPRINTI(outfile, delta); fprintf(outfile, "\n");
			GetReturn();
		}

		/* now calculate a, b, c, when A*x^2+B*x*y+C*y^2 is
		transformed into a*X^2+b*X*Y+c*Y^2 under the
		transformation x=alpha*X+beta*Y, y=gamma*X+delta*Y */

		TMP1 = MULTI3(A, beta, beta);
		TMP2 = MULTI3(B, beta, delta);
		TMP3 = MULTI3(C, delta, delta);
		TMP4 = ADDI(TMP1, TMP2);
		c = ADDI(TMP4, TMP3);
		FREEMPI(TMP1);
		FREEMPI(TMP2);
		FREEMPI(TMP3);
		FREEMPI(TMP4);
		TMP1 = MULTI3(TWOA, alpha, beta);
		TMP2 = MULTI3(TWOC, gamma, delta);
		TMP5 = MULTI3(B, alpha, delta);
		TMP6 = MULTI3(B, beta, gamma);
		/* now add TMP1, TMP2, TPM5, TMP6 */
		SUM = ADDI(TMP1, TMP2);
		FREEMPI(TMP1);
		FREEMPI(TMP2);
		TMP1 = SUM;
		SUM = ADDI(SUM, TMP5);
		FREEMPI(TMP1);
		FREEMPI(TMP5);
		TMP1 = SUM;
		SUM = ADDI(SUM, TMP6);
		FREEMPI(TMP1);
		FREEMPI(TMP6);
		b = COPYI(SUM);
		FREEMPI(SUM);
		/* Now we have gcd(a,N)=1 */
	}
	verb = CONVERTI(VERB);
	fclose(outfile);
	BINARYFORM(a, b, c, N, &FSX1, &FSY1, &length, FLAG, verb);

	outfile = fopen(buff, "a");
	if (!t) {
		for (i = 0; i < length; i++) {
			printf("i=%u: ", i);
			fprintf(outfile, "i=%u: ", i);
			TMP1 = MULTI(alpha, FSX1->A[i]);
			TMP2 = MULTI(beta, FSY1->A[i]);
			TMP3 = MULTI(gamma, FSX1->A[i]);
			TMP4 = MULTI(delta, FSY1->A[i]);
			x = ADDI(TMP1, TMP2);
			FREEMPI(TMP1);
			FREEMPI(TMP2);
			y = ADDI(TMP3, TMP4);
			FREEMPI(TMP3);
			FREEMPI(TMP4);
			printf("Solution (");
			fprintf(outfile, "Solution (");
			PRINTI(x);
			FPRINTI(outfile, x);
			printf(",");
			fprintf(outfile, ",");
			PRINTI(y);
			FPRINTI(outfile, y);
			printf("): ");
			fprintf(outfile, "): ");
			Z = EUCLIDI(x, y, &DELTA, &BETA);
			FREEMPI(Z);
			BETA->S = -(BETA->S);
			TMP1 = MULTI3(TWOA, x, BETA);
			TMP2 = MULTI3(TWOC, y, DELTA);
			TMP5 = MULTI3(B, x, DELTA);
			TMP6 = MULTI3(B, BETA, y);
			/* now add TMP1, TMP2, TPM5, TMP6 */
			SUM = ADDI(TMP1, TMP2);
			FREEMPI(TMP1);
			FREEMPI(TMP2);
			FREEMPI(BETA);
			FREEMPI(DELTA);
			TMP1 = SUM;
			SUM = ADDI(SUM, TMP5);
			FREEMPI(TMP1);
			FREEMPI(TMP5);
			TMP1 = SUM;
			SUM = ADDI(SUM, TMP6);
			FREEMPI(TMP1);
			FREEMPI(TMP6);
			TMP1 = SUM;
			SUM = MOD(SUM, TWOABSN);
			FREEMPI(TMP1);
			if (RSV(SUM, ABSN) > 0) {
				TMP1 = SUM;
				SUM = SUBI(SUM, TWOABSN);
				FREEMPI(TMP1);
			}
			MODULUS = MULT_I(ABSN, 2);
			printf(": n = ");
			fprintf(outfile, ": n = ");
			PRINTI(SUM);
			FPRINTI(outfile, SUM);
			printf("(mod ");
			PRINTI(MODULUS);
			printf(")\n");
			fprintf(outfile, "(mod ");
			FPRINTI(outfile, MODULUS);
			fprintf(outfile, ")\n");
			FREEMPI(SUM);
			FREEMPI(MODULUS);
			FREEMPI(x);
			FREEMPI(y);
		}
		FREEMPI(alpha);
		FREEMPI(beta);
		FREEMPI(gamma);
		FREEMPI(delta);
	}
	FREEMPIA(FSX1);
	FREEMPIA(FSY1);
	FREEMPI(a);
	FREEMPI(b);
	FREEMPI(c);
	FREEMPI(TWOA);
	FREEMPI(TWOC);
	FREEMPI(ABSN);
	FREEMPI(TWOABSN);
	FREEMPI(A);
	FREEMPI(B);
	FREEMPI(C);
	FREEMPI(N);
	fclose(outfile);

	return;
}
