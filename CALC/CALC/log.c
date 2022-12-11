#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_NONSTDC_NO_DEPRECATE
#include "integer.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fun.h"
#include "stack.h"
#ifdef _WIN32
#include "unistd_DOS.h"
#else
#include <unistd.h>
#endif

void INTLOG(Stack s)
/* This performs variant 0 of Shank's log_b(a),
* working as in bc scale 2r, but doing it in integers.
* See http://www.maths.uq.edu.au/~krm/log.pdf
* We print out the partial quotients n[0],...,n[r-1] if a>b,
* n[0],...,n[r] if a<b,
* and, with c=10^{2r},  int[a[i]*c], where the a[i] are Shank's:
* b[0]=[a[0]*c],,...,b[r+1]=[a[r+1]*c] if a[0]=a < a[1]=b,
* b[0]=[a[1]*c],,...,b[r]=[a[r+1]*c] if a[0]=a > a[1]=b.
* The output is not guaranteed to be correct, but seems certain to be so.
* See D. Shanks, "A logarithm algorithm", Math. Tables and Other Aids to
* Computation 8, (1954). 60--64.
*/
{
	USI j;
	USL i, r, rr, n;
	MPI *AA, *BB, *C, *TMP1, *TMP2;

	MPI *A = stackPop(s);
	MPI *B = stackPop(s);
	MPI *R = stackPop(s);
	MPI *N = stackPop(s);
	i = 0;
	j = 2;
	rr = CONVERTI(R);
	FREEMPI(R);
	n = CONVERTI(N);
	FREEMPI(N);
	if (RSV(A, B) == -1)
		rr++;
	r = rr + rr;
	C = POWER_I(10, r);
	AA = MULTI(A, C);
	if (n) {
		printf("A[0] = ");
		PRINTI(AA);
		printf("\n");
	}
	BB = MULTI(B, C);
	if (n) {
		printf("A[1] = ");
		PRINTI(BB);
		printf("\n");
	}
	while (rr) {
		if (RSV(AA, BB) >= 0) {
			TMP1 = MULTI(AA, C);
			TMP2 = AA;
			AA = INT(TMP1, BB);
			FREEMPI(TMP1);
			FREEMPI(TMP2);
			i++;
		}
		else {
			if (n) {
				printf("n[%u]=%lu, A[%u]=", j - 2, i, j);
				PRINTI(AA);
			}
			else
				printf("n[%u]=%lu", j - 2, i);
			printf("\n");
			j++;
			TMP1 = BB;
			BB = AA;
			AA = TMP1;
			rr--;
			i = 0; /* reset the partial quotient count */

		}
	}
	FREEMPI(AA);
	FREEMPI(BB);
	FREEMPI(A);
	FREEMPI(B);
	FREEMPI(C);
	return;
}

void LOG_(Stack s)
/* This performs Shank's log a to base b algorithm, working in effect with
* scale 2r, but doing it in integers.
* We work with scale = 2r and believe that this is sufficient to correctly
* output partial quotients a[0],...,a[s], where s=r-1 if n>b, but r if n<b.
* the decimal expansion of log_b{n} is printed, truncated correct to
* as many decimal places as possible: the r-1th convergent is compared with
* rth convergent and decimal expansion is truncated where they differ.
* e=0 prints only the partial quotients, convergents and decimal expansion,
* while e !=0 is verbose.
* Initially written in by Alan Offer and improved by Sean Seefried, Dec 1999.
* Modified by Sean Seefried on 17th January 2000.
* Calc version by Keith Matthews, 13th March 2000.
* See D. Shanks, "A logarithm algorithm", Math. Tables and Other Aids to
* Computation 8, (1954). 60--64.
*/
{
	USI j;
	int tt;
	USL i, r, rr, n;
	MPI *AA, *BB, *C, *TMP1, *TMP2;
	MPI *PN, *QP, *QN, *PP, *X, *Y, *Z;
	char buff[20];
	FILE *outfile;

	MPI *A = stackPop(s);
	MPI *B = stackPop(s);
	MPI *R = stackPop(s);
	MPI *N = stackPop(s);
	i = 0;
	j = 2;
	PN = ONEI();
	QP = ONEI();
	QN = ZEROI();
	PP = ZEROI();
	rr = CONVERTI(R);
	FREEMPI(R);
	n = CONVERTI(N);
	FREEMPI(N);
	if (RSV(A, B) == -1)
		rr++;
	r = rr + rr;
	C = POWER_I(10, r);
	AA = MULTI(A, C);
	printf("Partial Quotient, Convergent\n");
	printf("----------------------------\n");
	if (n) {
		printf("A[0] = ");
		PRINTI(AA);
		printf("\n");
	}
	strcpy(buff, "log.out");
	if (access(buff, R_OK) == 0)
		unlink(buff);
	outfile = fopen(buff, "w");
	fprintf(outfile, "A[0] = ");
	FPRINTI(outfile, AA);
	fprintf(outfile, "\n");
	BB = MULTI(B, C);
	if (n) {
		printf("A[1] = ");
		PRINTI(BB);
		printf("\n");
	}
	fprintf(outfile, "A[1] = ");
	FPRINTI(outfile, BB);
	fprintf(outfile, "\n");
	while (rr) {
		if (RSV(AA, BB) >= 0) {
			TMP1 = MULTI(AA, C);
			TMP2 = AA;
			AA = INT(TMP1, BB);
			FREEMPI(TMP1);
			FREEMPI(TMP2);
			i++;
			TMP1 = PN;
			PN = ADDI(PN, PP);
			FREEMPI(TMP1);
			TMP1 = QN;
			QN = ADDI(QN, QP);
			FREEMPI(TMP1);
			if (EQUALI(BB, C)) {
				PRINTI(A);
				FPRINTI(outfile, A);
				printf("^");
				fprintf(outfile, "^");
				PRINTI(PP);
				FPRINTI(outfile, PP);
				printf("=");
				fprintf(outfile, "=");
				PRINTI(B);
				FPRINTI(outfile, B);
				printf("^");
				fprintf(outfile, "^");
				PRINTI(QP);
				FPRINTI(outfile, QP);
				printf("\n");
				fprintf(outfile, "\n");
				FREEMPI(AA);
				FREEMPI(BB);
				FREEMPI(A);
				FREEMPI(B);
				FREEMPI(C);
				FREEMPI(PN);
				FREEMPI(PP);
				FREEMPI(QP);
				FREEMPI(QN);
				return;
			}
		}
		else {
			if (n) {
				printf("n[%u]=%lu, A[%u]=", j - 2, i, j);
				PRINTI(AA);
			}
			else
				printf("n[%u]=%lu", j - 2, i);
			fprintf(outfile, "n[%u]=%lu, A[%u]=", j - 2, i, j);
			FPRINTI(outfile, AA);
			j++;
			TMP1 = BB;
			BB = AA;
			AA = TMP1;
			TMP1 = PN;
			PN = PP;
			PP = TMP1;
			TMP1 = QN;
			QN = QP;
			QP = TMP1;
			printf(", p[%u]/q[%u]=", j - 3, j - 3);
			PRINTI(QP);
			printf("/");
			PRINTI(PP);
			fprintf(outfile, ", p[%u]/q[%u]=", j - 3, j - 3);
			FPRINTI(outfile, QP);
			fprintf(outfile, "/");
			FPRINTI(outfile, PP);
			printf("\n");
			fprintf(outfile, "\n");
			rr--;
			i = 0; /* reset the partial quotient count */

		}
	}
	printf("The log of ");
	fprintf(outfile, "The log of ");
	PRINTI(A);
	FPRINTI(outfile, A);
	printf(" to base ");
	fprintf(outfile, " to base ");
	PRINTI(B);
	FPRINTI(outfile, B);
	if (j == 3) {
		printf(" has integer part ");
		fprintf(outfile, " has integer part ");
		PRINTI(QP);
		FPRINTI(outfile, QP);
		printf("\n");
		fprintf(outfile, "\n");
		FREEMPI(AA);
		FREEMPI(BB);
		FREEMPI(A);
		FREEMPI(B);
		FREEMPI(C);
		FREEMPI(PN);
		FREEMPI(PP);
		FREEMPI(QP);
		FREEMPI(QN);
		fclose(outfile);
		return;
	}
	printf(" equals ");
	fprintf(outfile, " equals ");
	if (PN->S) {
		fclose(outfile);
		tt = COMPARE_DIGITS(QP, PP, QN, PN, 10);
		outfile = fopen(buff, "a");
		if (tt == -1) {
			X = MULTI(QP, PN);
			Y = MULTI(PP, QN);
			printf(" has integer part ");
			fprintf(outfile, " has integer part ");
			if (RSV(X, Y)< 0)
				Z = INTI(QP, PP);
			else
				Z = INTI(QN, PN);
			PRINTI(Z);
			FPRINTI(outfile, Z);
			FREEMPI(Z);
			FREEMPI(X);
			FREEMPI(Y);
			printf("\n");
			fprintf(outfile, "\n");
			FREEMPI(AA);
			FREEMPI(BB);
			FREEMPI(A);
			FREEMPI(B);
			FREEMPI(C);
			FREEMPI(PN);
			FREEMPI(PP);
			FREEMPI(QP);
			FREEMPI(QN);
			fclose(outfile);
			return;
		}
		else {
			printf(" truncated to %u decimal places\n", tt);
			fprintf(outfile, " truncated to %u decimal places\n", tt);
		}
	}
	FREEMPI(AA);
	FREEMPI(BB);
	FREEMPI(A);
	FREEMPI(B);
	FREEMPI(C);
	FREEMPI(PN);
	FREEMPI(PP);
	FREEMPI(QP);
	FREEMPI(QN);
	fclose(outfile);
	return;
}

int COMPARE_DIGITS(MPI *M, MPI *N, MPI *U, MPI *V, USL b)
/* m,n,u,v are positive integers, b is the base.
* We compare the base b expansions. If the integer parts
* of m/n and u/v are equal, we go on to compare the expansions
* after the decimal point. We print out the t common decimals and
* the number t.
*/
{
	USI i, o;
	MPI *R, *S, *T, *P, *RR, *PP;
	MPI *MM, *UU, *TMP1;
	char buff[20];
	FILE *outfile;

	R = INTI(M, N);
	RR = INTI(U, V);
	strcpy(buff, "log.out");
	outfile = fopen(buff, "a");
	if (EQUALI(R, RR)) {
		/* integer parts are equal */
		PRINTI(R);
		FPRINTI(outfile, R);
	}
	else {
		FREEMPI(R);
		FREEMPI(RR);
		fclose(outfile);
		return(-1);
	}
	TMP1 = MULTI(R, N);
	MM = SUBI(M, TMP1);
	FREEMPI(TMP1);
	TMP1 = MULTI(RR, V);
	UU = SUBI(U, TMP1);
	FREEMPI(TMP1);
	i = 0;
	P = MULT_I(MM, b);
	FREEMPI(MM);
	PP = MULT_I(UU, b);
	FREEMPI(UU);
	S = INT(P, N); /* the first decimal digit of {m/n} */
	T = INT(PP, V);/* the first decimal digit of {u/v} */

	FREEMPI(R);
	FREEMPI(RR);
	R = MOD(P, N);
	RR = MOD(PP, V);
	FREEMPI(P);
	FREEMPI(PP);
	printf(".");
	fprintf(outfile, ".");
	o = EQUALI(S, T);
	while (o) {
		PRINTI(S);
		FPRINTI(outfile, S);
		FREEMPI(S);
		FREEMPI(T);
		P = MULT_I(R, b);
		PP = MULT_I(RR, b);
		FREEMPI(R);
		FREEMPI(RR);
		R = MOD(P, N);
		RR = MOD(PP, V);
		S = INT(P, N); /* the first decimal digit of {m/n} */
		T = INT(PP, V);/* the first decimal digit of {u/v} */
		o = EQUALI(S, T);
		FREEMPI(P);
		FREEMPI(PP);
		i++;
	}
	FREEMPI(S);
	FREEMPI(T);
	FREEMPI(R);
	FREEMPI(RR);
	printf("\n");
	fprintf(outfile, "\n");
	fclose(outfile);
	return (i);
	/* the no. of common decimal digits after the decimal point */
}

void INTLOG1(Stack s)
/* This performs variant 1 of Shank's log_b(a),
* See http://www.maths.uq.edu.au/~krm/log.pdf
*/
{
	USI j;
	USL i, r, rr, n;
	MPI *AA, *BB, *C, *TMP1, *TMP2, *E, *F;

	MPI *A = stackPop(s);
	MPI *B = stackPop(s);
	MPI *R = stackPop(s);
	MPI *N = stackPop(s);
	i = 0;
	j = 2;
	rr = CONVERTI(R);
	FREEMPI(R);
	n = CONVERTI(N);
	FREEMPI(N);
	if (RSV(A, B) == -1)
		rr++;
	r = rr + rr;
	C = POWER_I(10, r);
	F = POWER_I(10, rr + 2);
	E = ADD0I(C, F);
	FREEMPI(F);
	AA = MULTI(A, C);
	if (n) {
		printf("A[0] = "); PRINTI(AA); printf("\n");
	}
	BB = MULTI(B, C);
	if (n) {
		printf("A[1] = "); PRINTI(BB); printf("\n");
	}
	while (RSV(BB, E) >= 0) {
		if (RSV(AA, BB) >= 0) {
			TMP1 = MULTI(AA, C);
			TMP2 = AA;
			AA = INT(TMP1, BB);
			FREEMPI(TMP1);
			FREEMPI(TMP2);
			i++;
		}
		else {
			if (n) {
				printf("n[%u]=%lu, A[%u]=", j - 2, i, j);
				PRINTI(AA);
			}
			else
				printf("n[%u]=%lu", j - 2, i);
			printf("\n");
			j++;
			TMP1 = BB;
			BB = AA;
			AA = TMP1;
			/*rr--;*/
			i = 0; /* reset the partial quotient count */

		}
	}
	FREEMPI(AA);
	FREEMPI(BB);
	FREEMPI(A);
	FREEMPI(B);
	FREEMPI(C);
	FREEMPI(E);
	return;
}

void INTLOG2(Stack s)
/* This performs variant 2 of Shank's log_b(a),
* See http://www.maths.uq.edu.au/~krm/log.pdf
*/
{
	USI j;
	USL d, i, r, rr, n;
	MPI *AA, *BB, *C, *TMP1, *TMP2, *E, *F;
	MPI *U, *X, *Y, *Z;

	MPI *A = stackPop(s);
	MPI *B = stackPop(s);
	MPI *R = stackPop(s);
	MPI *N = stackPop(s);
	i = 0;
	j = 2;
	rr = CONVERTI(R);
	FREEMPI(R);
	n = CONVERTI(N);
	FREEMPI(N);
	if (RSV(A, B) == -1)
		rr++;
	r = rr + rr;
	C = POWER_I(10, r);
	F = POWER_I(10, rr + 2);
	E = ADD0I(C, F);
	FREEMPI(F);
	AA = MULTI(A, C);
	if (n) {
		printf("A[0] = "); PRINTI(AA); printf("\n");
	}
	BB = MULTI(B, C);
	if (n) {
		printf("A[1] = "); PRINTI(BB); printf("\n");
	}
	while (RSV(BB, E) >= 0) {
		Y = ONEI();
		X = COPYI(AA);
		while (1) {
			U = INT0(X, Y);
			TMP1 = X;
			TMP2 = Y;
			X = MULTI(C, X);
			Y = MULTI(BB, Y);
			FREEMPI(TMP1);
			FREEMPI(TMP2);
			Z = MULTI(C, Y);
			d = RSV(Z, X);
			FREEMPI(Z);
			if (d == 1)
				break;
			i++;
			FREEMPI(U);
		}
		FREEMPI(X);
		FREEMPI(Y);
		FREEMPI(AA);
		AA = COPYI(U);
		FREEMPI(U);
		if (n) {
			printf("n[%u]=%lu, A[%u]=", j - 2, i, j);
			PRINTI(AA);
		}
		else
			printf("n[%u]=%lu", j - 2, i);
		printf("\n");
		j++;
		TMP1 = BB;
		BB = AA;
		AA = TMP1;
		i = 0; /* reset the partial quotient count */
	}
	FREEMPI(AA);
	FREEMPI(BB);
	FREEMPI(A);
	FREEMPI(B);
	FREEMPI(C);
	FREEMPI(E);
	return;
}

void LOG1(Stack s)
/* This performs variant 1 of Shank's log a to base b algorithm, working in
* effect with scale 2r, but doing it in integers.
*/
{
	USI j;
	int tt;
	USL i, r, rr, n;
	MPI *AA, *BB, *C, *TMP1, *TMP2;
	MPI *PN, *QP, *QN, *PP, *X, *Y, *Z;
	MPI *E, *F;
	char buff[20];
	FILE *outfile;

	MPI *A = stackPop(s);
	MPI *B = stackPop(s);
	MPI *D = stackPop(s);
	MPI *R = stackPop(s);
	MPI *N = stackPop(s);
	i = 0;
	j = 2;
	PN = ONEI();
	QP = ONEI();
	QN = ZEROI();
	PP = ZEROI();
	rr = CONVERTI(R);
	FREEMPI(R);
	n = CONVERTI(N);
	FREEMPI(N);
	if (RSV(A, B) == -1)
		rr++;
	r = rr + rr;
	C = POWERI(D, r);
	F = POWERI(D, rr + 2);
	FREEMPI(D);
	E = ADD0I(C, F);
	FREEMPI(F);

	AA = MULTI(A, C);
	printf("Partial Quotient, Convergent\n");
	printf("----------------------------\n");
	if (n) {
		printf("A[0] = ");
		PRINTI(AA);
		printf("\n");
	}
	strcpy(buff, "log1.out");
	if (access(buff, R_OK) == 0)
		unlink(buff);

	outfile = fopen(buff, "w");
	fprintf(outfile, "A[0] = ");
	FPRINTI(outfile, AA);
	fprintf(outfile, "\n");
	BB = MULTI(B, C);
	if (n) {
		printf("A[1] = ");
		PRINTI(BB);
		printf("\n");
	}
	fprintf(outfile, "A[1] = ");
	FPRINTI(outfile, BB);
	fprintf(outfile, "\n");
	while (RSV(BB, E) >= 0) {
		if (RSV(AA, BB) >= 0) {
			TMP1 = MULTI(AA, C);
			TMP2 = AA;
			AA = INT(TMP1, BB);
			FREEMPI(TMP1);
			FREEMPI(TMP2);
			i++;
			TMP1 = PN;
			PN = ADDI(PN, PP);
			FREEMPI(TMP1);
			TMP1 = QN;
			QN = ADDI(QN, QP);
			FREEMPI(TMP1);
			if (EQUALI(BB, C)) {
				PRINTI(A);
				FPRINTI(outfile, A);
				printf("^");
				fprintf(outfile, "^");
				PRINTI(PP);
				FPRINTI(outfile, PP);
				printf("=");
				fprintf(outfile, "=");
				PRINTI(B);
				FPRINTI(outfile, B);
				printf("^");
				fprintf(outfile, "^");
				PRINTI(QP);
				FPRINTI(outfile, QP);
				printf("\n");
				fprintf(outfile, "\n");
				FREEMPI(AA);
				FREEMPI(BB);
				FREEMPI(A);
				FREEMPI(B);
				FREEMPI(C);
				FREEMPI(PN);
				FREEMPI(PP);
				FREEMPI(QP);
				FREEMPI(QN);
				return;
			}
		}
		else {
			if (n) {
				printf("n[%u]=%lu, A[%u]=", j - 2, i, j);
				PRINTI(AA);
			}
			else
				printf("n[%u]=%lu", j - 2, i);
			fprintf(outfile, "n[%u]=%lu, A[%u]=", j - 2, i, j);
			FPRINTI(outfile, AA);
			j++;
			TMP1 = BB;
			BB = AA;
			AA = TMP1;
			TMP1 = PN;
			PN = PP;
			PP = TMP1;
			TMP1 = QN;
			QN = QP;
			QP = TMP1;
			printf(", p[%u]/q[%u]=", j - 3, j - 3);
			PRINTI(QP);
			printf("/");
			PRINTI(PP);
			fprintf(outfile, ", p[%u]/q[%u]=", j - 3, j - 3);
			FPRINTI(outfile, QP);
			fprintf(outfile, "/");
			FPRINTI(outfile, PP);
			printf("\n");
			fprintf(outfile, "\n");
			/*	rr--; */
			i = 0; /* reset the partial quotient count */

		}
	}
	printf("The log of ");
	fprintf(outfile, "The log of ");
	PRINTI(A);
	FPRINTI(outfile, A);
	printf(" to base ");
	fprintf(outfile, " to base ");
	PRINTI(B);
	FPRINTI(outfile, B);
	if (j == 3) {
		printf(" has integer part ");
		fprintf(outfile, " has integer part ");
		PRINTI(QP);
		FPRINTI(outfile, QP);
		printf("\n");
		fprintf(outfile, "\n");
		FREEMPI(AA);
		FREEMPI(BB);
		FREEMPI(A);
		FREEMPI(B);
		FREEMPI(C);
		FREEMPI(PN);
		FREEMPI(PP);
		FREEMPI(QP);
		FREEMPI(QN);
		fclose(outfile);
		return;
	}
	printf(" equals ");
	fprintf(outfile, " equals ");
	if (PN->S) {
		fclose(outfile);
		tt = COMPARE_DIGITS(QP, PP, QN, PN, 10);
		outfile = fopen(buff, "a");
		if (tt == -1) {
			X = MULTI(QP, PN);
			Y = MULTI(PP, QN);
			printf(" has integer part ");
			fprintf(outfile, " has integer part ");
			if (RSV(X, Y)< 0)
				Z = INTI(QP, PP);
			else
				Z = INTI(QN, PN);
			PRINTI(Z);
			FPRINTI(outfile, Z);
			FREEMPI(Z);
			FREEMPI(X);
			FREEMPI(Y);
			printf("\n");
			fprintf(outfile, "\n");
			FREEMPI(AA);
			FREEMPI(BB);
			FREEMPI(A);
			FREEMPI(B);
			FREEMPI(C);
			FREEMPI(PN);
			FREEMPI(PP);
			FREEMPI(QP);
			FREEMPI(QN);
			fclose(outfile);
			return;
		}
		else {
			printf(" truncated to %u decimal places\n", tt);
			fprintf(outfile, " truncated to %u decimal places\n", tt);
		}
	}
	FREEMPI(AA);
	FREEMPI(BB);
	FREEMPI(A);
	FREEMPI(B);
	FREEMPI(C);
	FREEMPI(PN);
	FREEMPI(PP);
	FREEMPI(QP);
	FREEMPI(QN);
	FREEMPI(E);
	fclose(outfile);
	return;
}

void LOGG(Stack S)
/* This performs variant 1 of Shank's log_b(a),
* See http://www.maths.uq.edu.au/~krm/log.pdf
*/
{
	USL i, j, n, nn, r, rr, s, ss;
	MPI *A, *B, *C, *AA, *BB, *CC, *TMP1, *TMP2;

	MPI *X = stackPop(S);
	MPI *Y = stackPop(S);
	MPI *D = stackPop(S);
	MPI *R = stackPop(S);
	MPI *RR = stackPop(S);
	i = 0;
	j = 0;

	s = 2;
	ss = 2;
	r = CONVERTI(R);
	rr = CONVERTI(RR);
	FREEMPI(R);
	FREEMPI(RR);
	if (r >= rr) {
		FREEMPI(X);
		FREEMPI(Y);
		FREEMPI(D);
		execerror("arg 4 >= arg5", "");
	}
	C = POWERI(D, (USI)r);
	A = MULTI(X, C);
	B = MULTI(Y, C);
	CC = POWERI(D, (USI)rr);
	AA = MULTI(X, CC);
	BB = MULTI(Y, CC);
	FREEMPI(D);
	FREEMPI(X);
	FREEMPI(Y);
	while (RSV(B, C) >= 0 && RSV(BB, CC) >= 0) {
		while (RSV(A, B) >= 0) {
			TMP1 = MULTI(A, C);
			TMP2 = A;
			A = INT(TMP1, B);
			FREEMPI(TMP1);
			FREEMPI(TMP2);
			i++;
		}
		n = i;
		printf("n=%lu\n", n);
		s++;
		TMP1 = B;
		B = A;
		A = TMP1;
		i = 0; /* reset the partial quotient count */

		printf("AA=");
		PRINTI(AA);
		printf("\n");
		printf("BB=");
		PRINTI(BB);
		printf("\n");
		while (RSV(AA, BB) >= 0) {
			TMP1 = MULTI(AA, CC);
			TMP2 = AA;
			AA = INT(TMP1, BB);
			FREEMPI(TMP1);
			FREEMPI(TMP2);
			j++;
			printf("AA=");
			PRINTI(AA);
			printf("\n");
			printf("BB=");
			PRINTI(BB);
			printf("\n");
		}
		nn = j;
		printf("nn=%lu\n", nn);
		ss++;
		TMP1 = BB;
		BB = AA;
		AA = TMP1;
		j = 0; /* reset the partial quotient count */
		if (n != nn) {
			printf("%lu != %lu\n", n, nn);
			printf("breaking\n");
			break;
		}
		printf("n[%lu]=%lu\n", s - 3, n);
	}
	FREEMPI(A);
	FREEMPI(B);
	FREEMPI(AA);
	FREEMPI(BB);
	FREEMPI(C);
	FREEMPI(CC);
	return;
}

void LOGGG(Stack S)
/* This performs variant 1 of Shank's log_b(a),
* See http://www.maths.uq.edu.au/~krm/log.pdf
*/
{
	int tt;
	USL e, i, j, n, nn, r, rr, s, ss;
	MPI *A, *B, *C, *AA, *BB, *CC, *TMP1, *TMP2;
	MPI *PN, *QP, *QN, *PP, *X, *Y, *Z, *QNN, *PNN;
	char buff[20];
	FILE *outfile;

	MPI *XX = stackPop(S);
	MPI *YY = stackPop(S);
	MPI *D = stackPop(S);
	MPI *R = stackPop(S);
	MPI *RR = stackPop(S);
	MPI *E = stackPop(S);
	i = 0;
	j = 0;
	e = CONVERTI(E);
	FREEMPI(E);
	PN = ONEI();
	QP = ONEI();
	QN = ZEROI();
	PP = ZEROI();
	PNN = ONEI();
	QNN = ZEROI();

	s = 2;
	ss = 2;
	r = CONVERTI(R);
	rr = CONVERTI(RR);
	FREEMPI(R);
	FREEMPI(RR);
	if (r >= rr) {
		FREEMPI(XX);
		FREEMPI(YY);
		FREEMPI(D);
		FREEMPI(E);
		execerror("arg 4 >= arg5", "");
	}
	C = POWERI(D, (USI)r);
	A = MULTI(XX, C);
	B = MULTI(YY, C);
	CC = POWERI(D, (USI)rr);
	AA = MULTI(XX, CC);
	BB = MULTI(YY, CC);
	FREEMPI(D);
	strcpy(buff, "log.out");
	outfile = fopen(buff, "w");

	while (RSV(B, C) > 0 && RSV(BB, CC) > 0) {
		if (e) {
			FREEMPI(PNN);
			FREEMPI(QNN);
			PNN = COPYI(PN);
			QNN = COPYI(QN);
		}
		while (RSV(A, B) >= 0) {
			TMP1 = MULTI(A, C);
			TMP2 = A;
			A = INT(TMP1, B);
			FREEMPI(TMP1);
			FREEMPI(TMP2);
			i++;
			if (e) {
				TMP1 = PN;
				PN = ADDI(PN, PP);
				FREEMPI(TMP1);
				TMP1 = QN;
				QN = ADDI(QN, QP);
				FREEMPI(TMP1);
			}
		}
		n = i;
		s++;
		TMP1 = B;
		B = A;
		A = TMP1;
		i = 0; /* reset the partial quotient count */

		while (RSV(AA, BB) >= 0) {
			TMP1 = MULTI(AA, CC);
			TMP2 = AA;
			AA = INT(TMP1, BB);
			FREEMPI(TMP1);
			FREEMPI(TMP2);
			j++;
		}
		nn = j;
		ss++;
		TMP1 = BB;
		BB = AA;
		AA = TMP1;
		j = 0; /* reset the partial quotient count */
		if (n != nn) {
			/*	printf("%lu != %lu\n", n, nn);
			printf("breaking\n");*/
			break;
		}
		if (e) {
			printf("n[%lu]=%lu:", s - 3, n);
			fprintf(outfile, "n[%lu]=%lu:", s - 3, n);
		}
		else {
			printf("n[%lu]=%lu\n", s - 3, n);
			fprintf(outfile, "n[%lu]=%lu\n", s - 3, n);
		}
		if (e) {
			TMP1 = PN;
			PN = PP;
			PP = TMP1;
			TMP1 = QN;
			QN = QP;
			QP = TMP1;
			if (EQUALI(B, C)) {
				/* tricky-needed if B=C is reached and n !=nn is not encountered */
				FREEMPI(PNN);
				FREEMPI(QNN);
				PNN = COPYI(PN);
				QNN = COPYI(QN);
			}
			printf("p[%lu]/q[%lu]=", s - 3, s - 3);
			PRINTI(QP);
			printf("/");
			PRINTI(PP);
			fprintf(outfile, "p[%lu]/q[%lu]=", s - 3, s - 3);
			FPRINTI(outfile, QP);
			fprintf(outfile, "/");
			FPRINTI(outfile, PP);
			printf("\n");
			fprintf(outfile, "\n");
		}
	}
	if (e) {
		printf("The log of ");
		fprintf(outfile, "The log of ");
		PRINTI(XX);
		FPRINTI(outfile, XX);
		printf(" to base ");
		fprintf(outfile, " to base ");
		PRINTI(YY);
		FPRINTI(outfile, YY);
		FREEMPI(XX);
		FREEMPI(YY);
		if (s == 3) {
			printf(" has integer part ");
			fprintf(outfile, " has integer part ");
			PRINTI(QP);
			FPRINTI(outfile, QP);
			printf("\n");
			fprintf(outfile, "\n");
			FREEMPI(AA);
			FREEMPI(BB);
			FREEMPI(A);
			FREEMPI(B);
			FREEMPI(C);
			FREEMPI(PN);
			FREEMPI(PP);
			FREEMPI(QP);
			FREEMPI(QN);
			FREEMPI(PNN);
			FREEMPI(QNN);
			fclose(outfile);
			return;
		}
		printf(" equals ");
		fprintf(outfile, " equals ");
		if (PN->S) {
			fclose(outfile);
			tt = COMPARE_DIGITS(QP, PP, QNN, PNN, 10);
			outfile = fopen(buff, "a");
			if (tt == -1) {
				X = MULTI(QP, PN);
				Y = MULTI(PP, QN);
				printf(" has integer part ");
				fprintf(outfile, " has integer part ");
				if (RSV(X, Y)< 0)
					Z = INTI(QP, PP);
				else
					Z = INTI(QN, PN);
				PRINTI(Z);
				FPRINTI(outfile, Z);
				FREEMPI(Z);
				FREEMPI(X);
				FREEMPI(Y);
				printf("\n");
				fprintf(outfile, "\n");
				FREEMPI(AA);
				FREEMPI(BB);
				FREEMPI(A);
				FREEMPI(B);
				FREEMPI(C);
				FREEMPI(CC);
				FREEMPI(PN);
				FREEMPI(PP);
				FREEMPI(QP);
				FREEMPI(QN);
				FREEMPI(QNN);
				FREEMPI(PNN);
				fclose(outfile);
				return;
			}
			else {
				printf(" truncated to %u decimal places\n", tt);
				fprintf(outfile, " truncated to %u decimal places\n", tt);
			}
		}
	}

	FREEMPI(A);
	FREEMPI(B);
	FREEMPI(AA);
	FREEMPI(BB);
	FREEMPI(C);
	FREEMPI(CC);
	FREEMPI(PN);
	FREEMPI(PP);
	FREEMPI(QP);
	FREEMPI(QN);
	FREEMPI(XX);
	FREEMPI(YY);
	FREEMPI(QNN);
	FREEMPI(PNN);
	fclose(outfile);

	return;
}

void SHANKSLOG(Stack s)
/* This performs Shank's log_b(a),
* See http://www.maths.uq.edu.au/~krm/log.pdf
* See D. Shanks, "A logarithm algorithm", Math. Tables and Other Aids to
* Computation 8, (1954). 60--64.
* We print the a[i] for i <= n, to d decimal places.
*/
{
	USL i, n;
	USI j, d;
	MPR *AA, *BB, *TMP1;

	MPI *A = stackPop(s);
	MPI *B = stackPop(s);
	MPI *D = stackPop(s);
	MPI *N = stackPop(s);
	i = 0;
	j = 2;
	d = (USI)CONVERTI(D);
	FREEMPI(D);
	n = CONVERTI(N);
	FREEMPI(N);
	printf("a[0] = ");
	PRINTI(A);
	printf("\n");
	printf("a[1] = ");
	PRINTI(B);
	printf("\n");
	AA = BUILDMPR();
	AA->N = COPYI(A);
	AA->D = ONEI();
	FREEMPI(A);
	BB = BUILDMPR();
	BB->N = COPYI(B);
	BB->D = ONEI();
	FREEMPI(B);
	while (j <= n) {
		printf("Starting AA ");
		PRINTDR(d, AA);
		printf("; BB to ");
		PRINTDR(d, BB);
		printf("\n");
		while (COMPARER(AA, BB) >= 0) {
			TMP1 = AA;
			AA = RATIOR(AA, BB);
			printf("decreasing AA to ");
			PRINTDR(d, AA);
			printf("; BB=");
			PRINTDR(d, BB);
			printf("\n");
			FREEMPR(TMP1);
			i++;
		}
		printf("n[%u]=%lu, A[%u]=", j - 2, i, j);
		PRINTDR(d, AA);
		printf("\n");
		j++;
		TMP1 = BB;
		BB = AA;
		AA = TMP1;
		i = 0; /* reset the partial quotient count */
	}
	FREEMPR(AA);
	FREEMPR(BB);
	return;
}

void LOG(MPI *A, MPI *B, MPI *D, MPI *R, MPIA *M, MPI **L)
/* Returns an array M[] of L positive integers that are hopefully
* partial quotients of log(A)/log(B), using C=D^R.
* Here A > B > 1, D > 1, R >= 1
* Uses an algorithm in http://www.numbertheory.org/pdfs/log.pdf
*/
{
	USL i, ii, r, s;
	USI l, flag = 1;
	MPI *AA, *BB, *AAA, *BBB, *C, *TMP1, *TMP2, *I;
	MPIA P, Q;
	char buff[20];
	FILE *outfile;
	P = NULL;
	Q = NULL;

	s = 0;
	r = CONVERTI(R);
	C = POWERI(D, r);
	AA = MULTI(A, C);
	BB = MULTI(B, C);
	AAA = COPYI(AA);
	BBB = COPYI(BB);
	*M = BUILDMPIA();
	i = ii = 0;
	strcpy(buff, "log.out");
	outfile = fopen(buff, "w");
	fprintf(outfile, "Output from log(");
	FPRINTI(outfile, A);
	fprintf(outfile, ",");
	FPRINTI(outfile, B);
	fprintf(outfile, ",");
	FPRINTI(outfile, D);
	fprintf(outfile, ",");
	FPRINTI(outfile, R);
	fprintf(outfile, ")\n");

	while (RSV(BB, C) > 0 && RSV(BBB, C) > 0) {
		i = 0;
		while (RSV(AA, BB) >= 0) {
			TMP1 = MULTI(AA, C);
			TMP2 = AA;
			AA = INT(TMP1, BB);
			FREEMPI(TMP1);
			FREEMPI(TMP2);
			i++;
		}
		TMP1 = BB;
		BB = AA;
		AA = TMP1;
		ii = 0;
		while (RSV(AAA, BBB) >= 0) {
			TMP1 = MULTI(AAA, C);
			TMP2 = AAA;
			AAA = CEILINGI(TMP1, BBB);
			FREEMPI(TMP1);
			FREEMPI(TMP2);
			ii++;
		}
		TMP1 = BBB;
		BBB = AAA;
		AAA = TMP1;
		if (!EQUALI(BB, BBB))
			flag = 0;
		if (i == ii) {
			I = CHANGEI(i);
			ADD_TO_MPIA(*M, I, s);
			/*	fprintf(outfile, "m[%lu]=%lu\n", s, i);*/
			fprintf(outfile, "%lu", i); /*for use in Latexing*/
			if ((s + 1) % 20 != 0)
				fprintf(outfile, "&");
			else
				fprintf(outfile, "\\\\\n");

			FREEMPI(I);
		}
		else {
			printf("M[%lu]=%lu != %lu=M'[%lu]\n", s, i, ii, s);
			fprintf(outfile, "M[%lu]=%lu != %lu=M'[%lu]\n", s, i, ii, s);
			flag = 0;
			break;
		}
		s++;
	}
	*L = CHANGE(s);
	if (flag && EQUALI(BB, C) && EQUALI(BBB, C)) {
		printf("log("); PRINTI(A); printf(")");
		printf("/log("); PRINTI(B); printf(")");
		printf(" is likely to be the rational number ");
		l = (*M)->size;
		CONVERGENTS(*M, &P, &Q);
		PRINTI(P->A[l - 1]); printf("/"); PRINTI(Q->A[l - 1]);
		printf("\n");
		FREEMPIA(P);
		FREEMPIA(Q);
	}
	FREEMPI(C);
	FREEMPI(AA);
	FREEMPI(AAA);
	FREEMPI(BB);
	FREEMPI(BBB);
	fclose(outfile);
}

MPI *LOGX(MPI *A, MPI *B, MPI *D, MPI *R, MPIA *M, MPI **L)
{
	MPI *TMP = ONEI();
	int s;

	s = COMPAREI(B, TMP);
	if (s <= 0) {
		printf("B <= 1\n");
		FREEMPI(TMP);
		M = NULL;
		return NULL;
	}
	s = COMPAREI(D, TMP);
	if (s <= 0) {
		printf("D <= 1\n");
		FREEMPI(TMP);
		M = NULL;
		return NULL;
	}
	s = COMPAREI(A, B);
	if (s <= 0) {
		printf("A <= B\n");
		FREEMPI(TMP);
		M = NULL;
		return NULL;
	}
	s = COMPAREI(R, TMP);
	if (s < 0) {
		printf("R < 1\n");
		FREEMPI(TMP);
		M = NULL;
		return NULL;
	}
	FREEMPI(TMP);

	LOG(A, B, D, R, M, L);
	return ONEI();
}

void TESTLOG1(MPI *A, MPI *B, MPI *D, MPI *R)
/* This is similar to LOG, but is for use in TESTLOG2 */
{
	USL i, ii, r, s;
	MPI *AA, *BB, *AAA, *BBB, *C, *TMP1, *TMP2;
	char buff[20];
	FILE *outfile;

	s = 0;
	r = CONVERTI(R);
	C = POWERI(D, r);
	AA = MULTI(A, C);
	BB = MULTI(B, C);
	AAA = COPYI(AA);
	BBB = COPYI(BB);
	i = ii = 0;
	strcpy(buff, "testlog.out");
	outfile = fopen(buff, "a");
	while (RSV(BB, C) > 0 && RSV(BBB, C) > 0) {
		i = 0;
		while (RSV(AA, BB) >= 0) {
			TMP1 = MULTI(AA, C);
			TMP2 = AA;
			AA = INT(TMP1, BB);
			FREEMPI(TMP1);
			FREEMPI(TMP2);
			i++;
		}
		TMP1 = BB;
		BB = AA;
		AA = TMP1;
		ii = 0;
		while (RSV(AAA, BBB) >= 0) {
			TMP1 = MULTI(AAA, C);
			TMP2 = AAA;
			AAA = CEILINGI(TMP1, BBB);
			FREEMPI(TMP1);
			FREEMPI(TMP2);
			ii++;
		}
		TMP1 = BBB;
		BBB = AAA;
		AAA = TMP1;
		if (i == ii) {
			printf("%lu,", i);
			fflush(stdout);
			fprintf(outfile, "%lu,", i);
		}
		else
			break;
		s++;
	}
	printf("\n");
	fprintf(outfile, "\n");
	FREEMPI(C);
	FREEMPI(AA);
	FREEMPI(AAA);
	FREEMPI(BB);
	FREEMPI(BBB);
	fclose(outfile);
}

void TESTLOG(MPI *A, MPI *B, MPI *D, MPI *M, MPI *N)
/* runs TESTLOG1(A,B,D,R) for R=M,...,N. */
{
	USL r, m, n;
	MPI *R;
	char buff[20];
	FILE *outfile;

	strcpy(buff, "testlog.out");
	outfile = fopen(buff, "w");
	m = CONVERTI(M);
	n = CONVERTI(N);
	for (r = m; r <= n; r++) {
		R = CHANGE(r);
		TESTLOG1(A, B, D, R);
		FREEMPI(R);
	}
	fclose(outfile);
}

MPI *TESTLOGX(MPI *A, MPI *B, MPI *D, MPI *M, MPI *N)
{
	MPI *TMP = ONEI();
	int s;

	s = COMPAREI(B, TMP);
	if (s <= 0) {
		printf("B <= 1\n");
		FREEMPI(TMP);
		M = NULL;
		return NULL;
	}
	s = COMPAREI(D, TMP);
	if (s <= 0) {
		printf("D <= 1\n");
		FREEMPI(TMP);
		M = NULL;
		return NULL;
	}
	s = COMPAREI(A, B);
	if (s <= 0) {
		printf("A <= B\n");
		FREEMPI(TMP);
		M = NULL;
		return NULL;
	}
	s = COMPAREI(M, TMP);
	if (s < 0) {
		printf("M < 1\n");
		FREEMPI(TMP);
		M = NULL;
		return NULL;
	}
	s = COMPAREI(M, TMP);
	if (s < 0) {
		printf("M < 1\n");
		FREEMPI(TMP);
		M = NULL;
		return NULL;
	}
	s = COMPAREI(N, M);
	if (s < 0) {
		printf("N < M\n");
		FREEMPI(TMP);
		M = NULL;
		return NULL;
	}
	FREEMPI(TMP);

	TESTLOG(A, B, D, M, N);
	return ONEI();
}
