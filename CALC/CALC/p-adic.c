#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_NONSTDC_NO_DEPRECATE
#include "integer.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fun.h"
#include "stack.h"
/*#include "unistd.h"*/

/* MULT_PADIC forms the product of a=a[0]+a[1]p+...+a[m]p^m
* and b=b[0]+b[1]p+...+b[n]p^n, where a[m] and b[n] are nonzero.
* The output is prod[0]+prod[1]p+...+prod[l]p^l, where prod[l] is nonzero.
* If a or b is zero, we return 0 at the start, otherwise return l.
* The program is an adaption of one in i1.c
* from http://www.numbertheory.org/calc/
*/
void MULT_PADIC(MPIA A, MPIA B, MPI *P, MPIA *PROD, USI m, USI n, USI *l) {
	MPI *C, *T, *TEMP, *TEMP1, *TEMP2, *BK;
	USI j, k, temp;

	*PROD = BUILDMPIA();
	if ((m == 0 && EQZEROI(A->A[0])) || (n == 0 && EQZEROI(B->A[0]))) {
		TEMP = ZEROI();
		ADD_TO_MPIA(*PROD, TEMP, 0);
		FREEMPI(TEMP);
		*l = 0;
		return;
	}
	C = ZEROI();
	*l = m + n + 1;
	for (k = 0; k <= m; k++) {
		TEMP = ZEROI();
		ADD_TO_MPIA(*PROD, TEMP, k);
		FREEMPI(TEMP);
	}
	for (k = 0; k <= n; k++) {
		BK = B->A[k];
		if (EQZEROI(BK) == 0) {
			TEMP = C;
			C = ZEROI();
			FREEMPI(TEMP);
			for (j = 0; j <= m; j++) {
				temp = j + k;
				TEMP1 = MULTI(A->A[j], BK);
				TEMP2 = ADDI(TEMP1, (*PROD)->A[temp]);
				T = ADDI(TEMP2, C);
				TEMP = MOD(T, P);
				ADD_TO_MPIA(*PROD, TEMP, temp);
				FREEMPI(TEMP);
				FREEMPI(TEMP1);
				FREEMPI(TEMP2);
				TEMP = C;
				C = INT0(T, P);
				FREEMPI(TEMP);
				FREEMPI(T);
			}
			temp = m + k + 1;
			TEMP = COPYI(C);
			ADD_TO_MPIA(*PROD, TEMP, temp);
			FREEMPI(TEMP);
		}
		else {
			temp = m + k + 1;
			TEMP = ZEROI();
			ADD_TO_MPIA(*PROD, TEMP, temp);
			FREEMPI(TEMP);
		}
	}
	if (EQZEROI(C)) {
		*l = m + n;
	}
	FREEMPI(C);
	return;
}

/*
* RSV_PADIC() is an adaption of one in i1.c
* from http://www.numbertheory.org/calc/
* it returns 1,0,-1 according as A >,=,< B.
* It is assumed that A and B are to the same base.
*/
int RSV_PADIC(MPIA A, MPIA B, USI m, USI n) {
	USI j;
	int t;

	if (m > n)
		return(1);
	if (m < n)
		return(-1);
	j = m;
	while (EQUALI(A->A[j], B->A[j]) && j) {
		j--;
	}
	t = RSV(A->A[j], B->A[j]);
	if (t > 0)
		return(1);
	if (t < 0)
		return(-1);
	return(0);
}

/*
* SUB0_PADIC() is an adaption of one in i1.c
* from http://www.numbertheory.org/calc/
* Here A >= B and returns l and the array
* diff[]=A-B=diff[0]+diff[1]*P+...+diff[l]*P^l.
*/

void SUB_PADIC(MPIA A, MPIA B, MPI *P, MPIA *DIFF, USI m, USI n, USI *l) {
	USI d, j, temp;
	int k;
	MPI *C, *T, *TEMP, *TEMP1, *TEMP2, *ONE;

	*DIFF = BUILDMPIA();
	if (n == 0 && EQZEROI(B->A[0])) {
		for (j = 0; j <= m; j++) {
			TEMP = A->A[j];
			ADD_TO_MPIA(*DIFF, TEMP, j);
			FREEMPI(TEMP);
		}
		*l = m;
		return;
	}
	ONE = ONEI();
	C = ONEI();
	*l = m;
	for (j = 0; j <= n; j++) {
		TEMP = SUBI(A->A[j], B->A[j]);
		TEMP1 = ADD0I(P, C);
		TEMP2 = SUBI(TEMP1, ONE);
		T = ADDI(TEMP, TEMP2);
		FREEMPI(TEMP);
		FREEMPI(TEMP1);
		FREEMPI(TEMP2);
		TEMP = MOD0(T, P);
		ADD_TO_MPIA(*DIFF, TEMP, j);
		FREEMPI(TEMP);
		TEMP = C;
		C = INT0(T, P);
		FREEMPI(TEMP);
		FREEMPI(T);
	}
	if (m > n) {
		temp = n + 1;
		for (j = temp; j <= m; j++) {
			TEMP1 = ADD0I(P, C);
			TEMP2 = SUBI(TEMP1, ONE);
			T = ADD0I(A->A[j], TEMP2);
			FREEMPI(TEMP1);
			FREEMPI(TEMP2);
			TEMP = MOD0(T, P);
			ADD_TO_MPIA(*DIFF, TEMP, j);
			FREEMPI(TEMP);
			TEMP = C;
			C = INT0(T, P);
			FREEMPI(TEMP);
			FREEMPI(T);
		}
	}
	d = 0;
	k = (int)m;
	while (k >= 0) {
		if (EQZEROI((*DIFF)->A[k]) == 0) {
			d = k;
			break;
		}
		else
			k--;
	}
	if (m != d) {
		*l = d;
	}
	FREEMPI(ONE);
	FREEMPI(C);
	return;
}

/* base(b,n) gives the base b expansion of n > 0:
* base[]=baseb[0]+baseb[1]*b+ ...+baseb[i]*b^i.
* The integer i is returned, along with  baseb[].
*/
void BASE_PADIC(MPI *B, MPI *N, MPIA *BASE, USI *j) {
	MPI *TEMP, *T, *X, *Q;
	USI i;

	i = 0;
	X = COPYI(N);
	(*BASE) = BUILDMPIA();
	while (RSV(X, B) >= 0) {
		Q = INT0(X, B);
		TEMP = MULTI(Q, B);
		T = SUB0I(X, TEMP);
		FREEMPI(TEMP);
		TEMP = T;
		ADD_TO_MPIA(*BASE, TEMP, i);
		FREEMPI(TEMP);
		TEMP = X;
		X = Q;
		FREEMPI(TEMP);
		i++;
	}
	TEMP = X;
	ADD_TO_MPIA(*BASE, TEMP, i);
	FREEMPI(TEMP);
	*j = i;
	return;
}

/*
* twoadicsqrt(a,n) n>=3, returns the first n binary digits of a
* 2adic sqroot x of a positive integer a=8k+1. Here x=1 or 5 (mod 8).
* See http://www.numbertheory.org/courses/MP313/lectures/lecture23/page3.html
* and http://www.numbertheory.org/courses/MP313/solns/soln3/page1.html.
*/
void TWOADICSQRT(MPI *A, USI n, MPIA *DIGITS) {
	MPI *ONE, *TWO, *TEMP;
	USI i, j, k, l, x, y;
	MPIA BASE, PROD, DIFF;
	int z;
	char buff[20];
	FILE *outfile;

	strcpy(buff, "2-adic.out");
	outfile = fopen(buff, "w");

	ONE = ONEI();
	*DIGITS = BUILDMPIA();
	TEMP = ONEI();
	ADD_TO_MPIA(*DIGITS, TEMP, 0);
	FREEMPI(TEMP);
	if (n == 1) {
		printf("1=x[0] is the first digit of the 2-adic square root of ");
		fprintf(outfile, "1=x[0] is the first digit of the 2-adic square root of ");
		PRINTI(A);
		FPRINTI(outfile, A);
		printf("\n");
		fprintf(outfile, "\n");
		fflush(stdout);
		FREEMPI(ONE);
		printf("where x=x[0]+x[1]2+... is congruent to 1 (mod 4)\n");
		fprintf(outfile, "where x=x[0]+x[1]2+... is congruent to 1 (mod 4)\n");
		fflush(stdout);
		fclose(outfile);
		return;
	}
	TEMP = ZEROI();
	ADD_TO_MPIA(*DIGITS, TEMP, 1);
	FREEMPI(TEMP);
	if (n == 2) {
		printf("10=x[0]x[1] are the first two digits of the 2-adic square root of ");
		fprintf(outfile, "10=x[0]x[1] are the first two digits of the 2-adic square root of ");
		PRINTI(A);
		FPRINTI(outfile, A);
		printf("\n");
		fprintf(outfile, "\n");
		fflush(stdout);
		printf("where x=x[0]+x[1]2+... is congruent to 1 (mod 4)\n");
		fprintf(outfile, "where x=x[0]+x[1]2+... is congruent to 1 (mod 4)\n");
		fflush(stdout);
		FREEMPI(ONE);
		fclose(outfile);
		return;
	}
	for (k = 2; k < n; k++) {
		TEMP = ZEROI();
		ADD_TO_MPIA(*DIGITS, TEMP, k);
		FREEMPI(TEMP);
	}
	TWO = ADD0I(ONE, ONE);
	FREEMPI(ONE);
	BASE_PADIC(TWO, A, &BASE, &x);
	l = 0;
	for (k = 3; k <= n; k++) {
		MULT_PADIC(*DIGITS, *DIGITS, TWO, &PROD, l, l, &y);
		z = RSV_PADIC(PROD, BASE, y, x);
		if (z < 0) {
			SUB_PADIC(BASE, PROD, TWO, &DIFF, x, y, &j);
		}
		else {
			SUB_PADIC(PROD, BASE, TWO, &DIFF, y, x, &j);
		}
		FREEMPIA(PROD);
		if (EQZEROI(DIFF->A[k]) == 0) {
			TEMP = ONEI();
			ADD_TO_MPIA(*DIGITS, TEMP, k - 1);
			FREEMPI(TEMP);
			if (k == 3) {
				l = 2;
			}
			else {
				l = k - 1;
			}
		}
		FREEMPIA(DIFF);
	}
	FREEMPIA(BASE);
	for (i = 0; i < n; i++) {
		PRINTI((*DIGITS)->A[i]);
		FPRINTI(outfile, (*DIGITS)->A[i]);
		printf(" ");
		fprintf(outfile, " ");
		fflush(stdout);
	}
	if (n > 1) {
		printf("=x[0]...x[%u]\n", n - 1);
		fprintf(outfile, "=x[0]...x[%u]\n", n - 1);
		fflush(stdout);
	}
	else {
		printf("=x[0]\n");
		fprintf(outfile, "=x[0]\n");
		fflush(stdout);
	}
	if (n > 1) {
		printf("are the first %u digits of the 2-adic square root x of ", n);
		fprintf(outfile, "are the first %u digits of the 2-adic square root x of ", n);
	}
	else {
		printf("is the first digit of the 2-adic square root x of ");
		fprintf(outfile, "is the first digit of the 2-adic square root x of ");
	}
	PRINTI(A);
	FPRINTI(outfile, A);
	printf("\n");
	fprintf(outfile, "\n");
	printf("where x=x[0]+x[1]2+... is congruent to 1 (mod 4)\n");
	fprintf(outfile, "where x=x[0]+x[1]2+... is congruent to 1 (mod 4)\n");
	fclose(outfile);
	fflush(stdout);
	FREEMPI(TWO);
	return;
}

MPI *TWOADICSQRTX(MPI *A, MPI *N, MPIA *DIGITS) {
	USI n;
	USL t;
	if (A->S <= 0) {
		printf("A <= 0\n");
		return(NULL);
	}
	t = MOD0_(A, (USL)8);
	if (t != 1) {
		printf("n not of the from 8k+1\n");
		return(NULL);
	}
	if (N->S <= 0) {
		printf("n <= 0\n");
		return(NULL);
	}
	if (N->D > 1) {
		printf("n > 2^16.\n");
		return(NULL);
	}
	n = (USI)CONVERTI(N);
	TWOADICSQRT(A, n, DIGITS);
	return (ONEI());
}

/*
* PADICSQRT(a,n,p) returns the first n p-adic digits of a
* p-adic sqroot x of a positive integer a. Here x=b(mod p),
* where b^2=a(mod p) and 0<b<p.
* See http://www.numbertheory.org/courses/MP313/lectures/lecture23/page3.html
* and http://www.numbertheory.org/courses/MP313/solns/soln3/page3.html.
*/

void PADICSQRT(MPI *A, USI n, MPI *P, MPIA *DIGITS) {
	MPI *U, *PP, *TEMP, *TEMP1, *TEMP2;
	MPIA BASE, PROD, DIFF;
	USI i, j, temp, x, l, y, k;
	int z;
	long int sign;
	char buff[20];
	FILE *outfile;

	strcpy(buff, "p-adic.out");
	outfile = fopen(buff, "w");

	U = SQRTM(A, P);
	BASE_PADIC(P, A, &BASE, &x);
	sign = 1;
	l = 0;
	TEMP = COPYI(U);
	*DIGITS = BUILDMPIA();
	ADD_TO_MPIA(*DIGITS, TEMP, 0);
	FREEMPI(TEMP);
	for (k = 1; k <= n; k++) {
		MULT_PADIC(*DIGITS, *DIGITS, P, &PROD, l, l, &y);
		z = RSV_PADIC(PROD, BASE, y, x);
		if (z < 0) {
			SUB_PADIC(BASE, PROD, P, &DIFF, x, y, &j);
			sign = -1;
		}
		else {
			SUB_PADIC(PROD, BASE, P, &DIFF, y, x, &j);
		}
		FREEMPIA(PROD);
		if (EQZEROI(DIFF->A[k]) == 0) {
			TEMP = MULT_I(U, 2);
			TEMP1 = MULT_I(TEMP, sign);
			FREEMPI(TEMP);
			TEMP2 = MINUSI(DIFF->A[k]);
			TEMP = CONGR(TEMP1, TEMP2, P, &PP);
			FREEMPI(PP);
			FREEMPI(TEMP1);
			FREEMPI(TEMP2);
			ADD_TO_MPIA(*DIGITS, TEMP, k);
			FREEMPI(TEMP);
		}
		else {
			TEMP = ZEROI();
			ADD_TO_MPIA(*DIGITS, TEMP, k);
			FREEMPI(TEMP);
		}
		l = k;
		FREEMPIA(DIFF);
	}
	FREEMPIA(BASE);
	for (i = 0; i < n; i++) {
		if (i) {
			printf(",");
			fprintf(outfile, ",");
			PRINTI((*DIGITS)->A[i]);
			FPRINTI(outfile, (*DIGITS)->A[i]);
			fflush(stdout);
		}
		else {
			PRINTI((*DIGITS)->A[i]);
			FPRINTI(outfile, (*DIGITS)->A[i]);
			fflush(stdout);
		}
	}
	temp = n - 1;
	if (n > 1) {
		printf("=x[0]...x[%u]\n", n - 1);
		fprintf(outfile, "=x[0]...x[%u]\n", n - 1);
		fflush(stdout);
	}
	else {
		printf("=x[0]\n");
		fprintf(outfile, "=x[0]\n");
		fflush(stdout);
	}
	if (n > 1) {
		printf("are the first %u digits of the p-adic square root x of ", n);
		fprintf(outfile, "are the first %u digits of the p-adic square root x of ", n);
	}
	else {
		printf("is the first digit of the p-adic square root x of ");
		fprintf(outfile, "is the first digit of the p-adic square root x of ");
	}
	PRINTI(A);
	FPRINTI(outfile, A);
	printf("\n");
	fprintf(outfile, "\n");
	printf("where x=x[0]+x[1]P+... is congruent to ");
	fprintf(outfile, "where x=x[0]+x[1]P+... is congruent to ");
	PRINTI(U);
	FPRINTI(outfile, U);
	printf(" (mod ");
	fprintf(outfile, " (mod ");
	PRINTI(P);
	FPRINTI(outfile, P);
	printf(")\n");
	fprintf(outfile, ")\n");
	fclose(outfile);
	fflush(stdout);
	FREEMPI(U);
	return;
}

MPI *PADICSQRTX(MPI *A, MPI *N, MPI *P, MPIA *DIGITS) {
	int t;
	USI n;
	MPI *T;

	if (A->S <= 0) {
		printf("A <= 0\n");
		return(NULL);
	}
	if (N->S <= 0) {
		printf("n <= 0\n");
		return(NULL);
	}
	if (N->D > 1) {
		printf("n > 2^16.\n");
		return(NULL);
	}
	n = (USI)CONVERTI(N);

	if (P->S <= 0 || (P->D == 0 && P->V[0] <= 2)) {
		printf("P <= 2\n");
	}
	T = LUCAS(P);
	t = T->S;
	FREEMPI(T);
	if (!t)
	{
		printf("3rd argument is not a prime\n");
		return NULL;
	}
	if (JACOBIB(A, P) != 1)
	{
		printf("X is not a quadratic residue mod P\n");
		return NULL;
	}
	PADICSQRT(A, n, P, DIGITS);
	return(ONEI());
}
