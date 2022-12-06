/* reduction.c */
#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_NONSTDC_NO_DEPRECATE
#ifdef _WIN32
#include "unistd_DOS.h"
#else
#include <unistd.h>
#endif
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "integer.h"
#include "fun.h"

extern unsigned long PRIME[];
unsigned int llen;

unsigned int SQUAREFREE_TEST1000(USL n) {
	USI i;
	USL temp;
	for (i = 0; i <= 167; i++) {
		temp = PRIME[i] * PRIME[i];
		if (n % temp == 0)
			return(0);
	}
	return(1);
}

USI global_count;
MPIA GLOBALARRAY_U;
MPIA GLOBALARRAY_V;

unsigned int PERIOD(MPI *D, MPI *U, MPI *V) {
	USI j, k, cycle_length;
	MPI *F, *R, *S, *A, *tmp, *tmp1, *tmp2;

	F = BIG_MTHROOT(D, 2);
	k = global_count; /* initially set to 0 in POS() below */
	R = COPYI(U);
	S = COPYI(V);
	for (j = k; 1; j++) {
		tmp = ADD0I(F, R);
		A = INT0(tmp, S);
		FREEMPI(tmp);
		tmp = MULTI(A, S);
		FREEMPI(A);
		tmp1 = R;
		R = SUB0I(tmp, R);
		FREEMPI(tmp);
		FREEMPI(tmp1);
		tmp = S;
		tmp1 = MULTI(R, R);
		tmp2 = SUB0I(D, tmp1);
		S = INT0(tmp2, S);
		ADD_TO_MPIA(GLOBALARRAY_U, R, j);
		ADD_TO_MPIA(GLOBALARRAY_V, S, j);
		FREEMPI(tmp);
		FREEMPI(tmp1);
		FREEMPI(tmp2);

		if (EQUALI(U, R) && EQUALI(V, S)) {
			FREEMPI(R);
			FREEMPI(S);
			FREEMPI(F);
			break;
		}
		if (j == R0) {
			FREEMPI(R);
			FREEMPI(S);
			FREEMPI(F);
			execerror("j = R0", "");
		}
	}
	cycle_length = j + 1 - k;
	global_count = global_count + cycle_length;
	return(cycle_length);
}

unsigned int CHECK_UVARRAYS(MPI *U, MPI *V) {
	USI i;
	for (i = 0; i < global_count; i++) {
		if (EQUALI(U, (GLOBALARRAY_U)->A[i]) && EQUALI(V, (GLOBALARRAY_V)->A[i])) {
			return(0);
		}
	}
	return(1);
}

USI POS(MPI *D) {
	USI a, b, e, t, z, parity_flag, class_number, x;
	USL d, f; /* d will be restricted to d < 10^6. */
	MPI *A, *B, *C, *TEMP, *F, *G, *H, *I, *J, *K, *L;
	MPI *ONE, *R, *DD;
	int sign, sign1;

	GLOBALARRAY_U = BUILDMPIA();
	GLOBALARRAY_V = BUILDMPIA();
	class_number = 0;
	ONE = ONEI();
	global_count = 0;
	parity_flag = 0;
	d = CONVERTI(D);
	/* creates a fundamental discriminant */
	if ((d - 1) % 4 != 0) {
		d = 4 * d;
		e = 2;
	}
	else {
		e = 1;
	}
	DD = CHANGE(d);
	printf("Creating a complete list of reduced forms of discriminant %lu\n", d);
	F = BIG_MTHROOT(DD, 2);
	f = CONVERTI(F);
	G = INT0_(F, 2);
	for (a = 1; a <= f; a++) {
		A = CHANGE(a);
		I = MULT_I(A, 2);
		J = MULT_I(A, 4);
		for (b = e; b <= f; b = b + 2) {
			B = CHANGE(b);
			TEMP = MULTI(B, B);
			H = SUBI(TEMP, DD);
			FREEMPI(TEMP);
			TEMP = MOD(H, J);
			sign = TEMP->S;
			FREEMPI(TEMP);
			if (sign == 0) {
				TEMP = SUBI(F, I);
				sign1 = COMPAREI(TEMP, B);
				FREEMPI(TEMP);
				if ((RSV(A, G) <= 0 && sign1 < 0) || (RSV(A, G) > 0 && sign1 >= 0)) {
					R = INT(H, J);
					TEMP = ABSI(R);
					K = GCD(A, B);
					L = GCD(K, TEMP);
					FREEMPI(K);
					FREEMPI(TEMP);
					t = EQONEI(L);
					FREEMPI(L);
					if (t) {
						TEMP = MINUSI(H);
						C = INT0(TEMP, J); /* C > 0 */
						FREEMPI(TEMP);
						TEMP = MULT_I(C, 2);
						x = CHECK_UVARRAYS(B, TEMP);
						if (x) {
							class_number = class_number + 1;
							z = PERIOD(DD, B, TEMP);
							if (z % 2) {
								parity_flag = 1;
							}
							printf("[%u]: (", class_number);
							PRINTI(A);
							printf(", ");
							PRINTI(B);
							printf(", ");
							PRINTI(R);
							printf(")\n");
						}
						FREEMPI(TEMP);
						FREEMPI(C);
					}
					FREEMPI(R);
				}
			}
			FREEMPI(H);
			FREEMPI(B);
		}
		FREEMPI(A);
		FREEMPI(I);
		FREEMPI(J);
	}
	if (parity_flag) {
		printf("x^2-"); PRINTI(D); printf("*y^2=-4 has a solution\n");
	}
	else {
		printf("x^2-"); PRINTI(D); printf("*y^2=-4 has no solution\n");
	}
	printf("h(%lu)=%u\n", d, class_number);
	FREEMPI(ONE);
	FREEMPI(F);
	FREEMPI(G);
	FREEMPI(DD);
	FREEMPIA(GLOBALARRAY_U);
	FREEMPIA(GLOBALARRAY_V);
	return(class_number);
}

MPI *POSX(MPI *D) {
	USI class_number, w;
	USL d;
	int t;
	MPI *ONE, *N, *TEMP;

	d = CONVERTI(D);
	ONE = ONEI();
	t = COMPAREI(D, ONE);
	FREEMPI(ONE);
	if (t <= 0) {
		printf(" D <= 1\n");
		return(NULL);
	}
	TEMP = POWER_I(10, 6);
	t = RSV(D, TEMP);
	FREEMPI(TEMP);
	if (t >= 0) {
		printf("D >= 10^6\n");
		return(NULL);
	}
	w = SQUAREFREE_TEST1000(d);
	if (w == 0) {
		printf("D is not squarefree\n");
		return(NULL);
	}
	else {
		class_number = POS(D);
		N = CHANGE(class_number);
		return(N);
	}
}

USI NEG(MPI *D, MPI *FLAG, MPI *TABLE_FLAG) {
	/*
	* This is Henri Cohen's Algorithm 5.3.5, p. 228.
	* for finding the class number h(D) of binary quadratic forms
	* of discriminant D, when D<0.
	* Here D=0(mod 4) or 1(mod 4).
	* If flag=1, we print only the primitive forms.
	* h(D) is returned in each case.
	* Davenport's Higher Arithmetic has a table of forms, which
	* lists the imprimitive ones with an asterisk.
	* If d is the discriminant of an imaginary quadratic field K,
	* then the primitive forms class-number h(D) is also the class number of K.
	* We print forms only when TABLE_FLAG is nonzero.
	*/
	MPI *A, *B, *BB, *GG, *GGG, *TEMP, *TEMP1, *ABSD, *ONE;
	MPI *Q, *QQ, *T;
	int t;
	USI g, h;
	USL temp;

	ONE = ONEI();
	h = 0;
	g = 1;
	if (TABLE_FLAG->S) {
		if (FLAG->S == 1) {
			printf("determining primitive forms\n");
		}
		else {
			printf("determining primitive and imprimitive forms\n");
		}
	}
	temp = MOD_(D, 4);
	if (temp == 0) {
		B = ZEROI();
	}
	else {
		B = ONEI();
	}
	ABSD = MINUSI(D);
	TEMP = INT0_(ABSD, 3);
	BB = BIG_MTHROOT(TEMP, 2);
	FREEMPI(TEMP);
	Q = NULL;
	A = NULL;
	while (RSV(B, BB) <= 0) {
		TEMP = MULTI(B, B);
		TEMP1 = ADD0I(TEMP, ABSD);
		FREEMPI(TEMP);
		Q = INT0_(TEMP1, 4);
		FREEMPI(TEMP1);
		A = COPYI(B);
		if (RSV(A, ONE) <= 0) {
			TEMP = A;
			A = ONEI();
			FREEMPI(TEMP);
		}
		QQ = BIG_MTHROOT(Q, 2);
		while (RSV(A, QQ) <= 0) {
			TEMP = MOD0(Q, A);
			t = TEMP->S;
			FREEMPI(TEMP);
			if (t == 0) {
				T = INT0(Q, A);
				if (FLAG->S == 1) {
					GG = GCD(A, B);
					GGG = GCD(GG, T);
					FREEMPI(GG);
					if (RSV(GGG, ONE) > 0) {
						g = 0;
					}
					FREEMPI(GGG);
				}
				if (g == 1) {
					if (RSV(A, B) == 0 || RSV(A, T) == 0 || B->S == 0) {
						if (TABLE_FLAG->S) {
							printf("(");
							PRINTI(A);
							printf(",");
							PRINTI(B);
							printf(",");
							PRINTI(T);
							printf(")\n");
						}
						h = h + 1;
					}
					else {
						if (TABLE_FLAG->S) {
							printf("(");
							PRINTI(A);
							printf(",");
							PRINTI(B);
							printf(",");
							PRINTI(T);
							printf(")\n");

							printf("(");
							PRINTI(A);
							printf(",");
							TEMP = MINUSI(B);
							PRINTI(TEMP);
							FREEMPI(TEMP);
							printf(",");
							PRINTI(T);
							printf(")\n");
						}
						h = h + 2;
					}
				}
				else {
					g = 1;
				}
				FREEMPI(T);
			}
			TEMP = A;
			A = ADD0_I(A, 1);
			FREEMPI(TEMP);
		}
		FREEMPI(A);
		FREEMPI(QQ);
		FREEMPI(Q);
		TEMP = B;
		B = ADD0_I(B, 2);
		FREEMPI(TEMP);
	}
	FREEMPI(B);
	FREEMPI(BB);
	FREEMPI(ABSD);
	FREEMPI(ONE);
	return(h);
}

MPI *NEGX(MPI *D, MPI *FLAG) {
	MPI *ZERO, *H, *TEMP, *ONE;
	USL temp;
	USI h;
	int s, t;

	ZERO = ZEROI();
	ONE = ONEI();
	t = COMPAREI(D, ZERO);
	if (t >= 0) {
		printf("D >= 0\n");
		FREEMPI(ZERO);
		FREEMPI(ONE);
		return(NULL);
	}
	TEMP = POWER_I(10, 6);
	s = RSV(D, TEMP);
	FREEMPI(TEMP);
	if (s >= 0) {
		printf("|D| >= 10^6\n");
		FREEMPI(ZERO);
		FREEMPI(ONE);
		return(NULL);
	}
	t = COMPAREI(FLAG, ONE);
	if (t > 0) {
		printf("flag > 1\n");
		FREEMPI(ONE);
		FREEMPI(ZERO);
		return(NULL);
	}
	t = COMPAREI(FLAG, ZERO);
	FREEMPI(ZERO);
	if (t < 0) {
		printf("flag < 0\n");
		FREEMPI(ONE);
		return(NULL);
	}
	temp = MOD_(D, 4);
	if (temp == 2 || temp == 3) {
		printf("D is not congruent to 0 or 1 (mod 4)\n");
		FREEMPI(ONE);
		return(NULL);
	}
	h = NEG(D, FLAG, ONE);
	FREEMPI(ONE);
	H = CHANGE((USL)h);
	return(H);
}

USI REDUCE_NEG(MPI *A, MPI *B, MPI *C) {
	/* This is Gauss' algorithm for reducing a positive definite binary
	* quadratic form. See L.E. Dickson, Introduction to the theory of numbers,
	* page 69. Here d=b^2-4ac <0, a>0,c>0, while the reduced form (A,B,C)
	* satisfies -A<B<=A, C>=A, with B>=0 if C=A.
	* The number of steps taken in the reduction is returned.
	*/
	MPI *D, *X, *Y, *U, *V, *TEMP1, *TEMP2, *TEMP3;
	MPI *a, *b, *c, *TEMPX, *TEMPY, *TEMPU, *TEMPV, *MINUSD, *DELTA;
	MPI *TEMPB;
	USI i;
	int r, s, t;

	i = 0;
	TEMP1 = MULTI(B, B);
	TEMP2 = MULTI(A, C);
	TEMP3 = MULT_I(TEMP2, 4);
	D = SUBI(TEMP1, TEMP3);
	FREEMPI(TEMP1);
	FREEMPI(TEMP2);
	FREEMPI(TEMP3);
	TEMP1 = MINUSI(A);
	r = COMPAREI(TEMP1, B);
	FREEMPI(TEMP1);
	s = COMPAREI(B, A);
	t = COMPAREI(A, C);
	if (r < 0 && s <= 0 && t <= 0) {/* -A < B <= A <= C */
		r = COMPAREI(A, C);
		if (t < 0 || (r == 0 && B->S >= 0)) {/* A < C or A=C and B>=0 */
			printf("(");
			PRINTI(A);
			printf(",");
			PRINTI(B);
			printf(",");
			PRINTI(C);
			printf(") is reduced\n");
			printf("Transforming matrix:1,0,0,1\n");
			FREEMPI(D);
			return (0);
		}
	}
	X = ONEI();
	Y = ZEROI();
	U = ZEROI();
	V = ONEI();
	MINUSD = MINUSI(D);
	a = COPYI(A);
	b = COPYI(B);
	c = COPYI(C);
	while (1) {
		TEMPB = b;
		TEMP1 = MINUSI(b);
		TEMP2 = MULT_I(c, 2);
		b = HALFMOD(TEMP1, TEMP2);
		FREEMPI(TEMP1);
		TEMP1 = ADDI(TEMPB, b);
		FREEMPI(TEMPB);
		TEMP3 = MINUSI(TEMP1);
		FREEMPI(TEMP1);
		DELTA = INT(TEMP3, TEMP2);
		FREEMPI(TEMP2);
		FREEMPI(TEMP3);
		TEMPU = U;
		TEMPX = X;
		TEMPY = Y;
		TEMPV = V;
		X = MINUSI(Y);
		Y = MULTABC(TEMPX, Y, DELTA);
		FREEMPI(TEMPX);
		FREEMPI(TEMPY);
		U = MINUSI(V);
		V = MULTABC(TEMPU, V, DELTA);
		FREEMPI(TEMPU);
		FREEMPI(TEMPV);
		FREEMPI(DELTA);
		TEMP1 = a;
		a = COPYI(c);
		FREEMPI(TEMP1);
		TEMP1 = c;
		TEMP2 = MULTABC(MINUSD, b, b);
		TEMP3 = MULT_I(a, 4);
		c = INT(TEMP2, TEMP3);
		FREEMPI(TEMP1);
		FREEMPI(TEMP2);
		FREEMPI(TEMP3);
		i++;
		printf("(");
		PRINTI(a);
		printf(",");
		PRINTI(b);
		printf(",");
		PRINTI(c);
		printf(")\n");
		if (COMPAREI(c, a) >= 0) {
			break;
		}
	}
	if (COMPAREI(a, c) == 0 && b->S < 0) {
		TEMP1 = b;
		b = MINUSI(b);
		FREEMPI(TEMP1);
		TEMPX = X;
		TEMPY = Y;
		X = COPYI(Y);
		Y = MINUSI(TEMPX);
		FREEMPI(TEMPX);
		FREEMPI(TEMPY);
		TEMPU = U;
		TEMPV = V;
		U = COPYI(V);
		V = MINUSI(TEMPU);
		FREEMPI(TEMPU);
		FREEMPI(TEMPV);
	}
	printf("(");
	PRINTI(a);
	printf(",");
	PRINTI(b);
	printf(",");
	PRINTI(c);
	printf(") is reduced\n");
	printf("Transforming matrix: ");
	PRINTI(X);
	printf(",");
	PRINTI(Y);
	printf(",");
	PRINTI(U);
	printf(",");
	PRINTI(V);
	printf("\n");
	FREEMPI(MINUSD);
	FREEMPI(D);
	FREEMPI(X);
	FREEMPI(Y);
	FREEMPI(U);
	FREEMPI(V);
	FREEMPI(a);
	FREEMPI(b);
	FREEMPI(c);
	return(i);
}

MPI *REDUCE_NEGX(MPI *A, MPI *B, MPI *C) {
	MPI *D, *TEMP1, *TEMP2, *TEMP3;
	USI i;

	if (A->S <= 0) {
		printf("A <= 0\n");
		return(NULL);
	}
	if (C->S <= 0) {
		printf("C <= 0\n");
		return(NULL);
	}
	TEMP1 = MULTI(B, B);
	TEMP2 = MULTI(A, C);
	TEMP3 = MULT_I(TEMP2, 4);
	D = SUBI(TEMP1, TEMP3);
	FREEMPI(TEMP1);
	FREEMPI(TEMP2);
	FREEMPI(TEMP3);
	if (D->S >= 0) {
		printf("B^2-4*A*C >= 0\n");
		FREEMPI(D);
		return(NULL);
	}
	i = REDUCE_NEG(A, B, C);
	FREEMPI(D);
	return(CHANGE(i));
}

USI CYCLE_PERIOD(MPI *D, MPI *U, MPI *V, USI i, int sign) {
	MPI *A, *F, *UU, *VV, *X, *Y, *TEMP1, *TEMP2, *TEMP3;
	USI len, j;
	char buff[20];
	FILE *outfile;

	strcpy(buff, "reducepos.out");
	/*if(access(buff, R_OK) == 0)
	unlink(buff);*/
	outfile = fopen(buff, "a");
	F = BIG_MTHROOT(D, 2);
	UU = COPYI(U);
	VV = COPYI(V);
	printf("\ncycle:\n->");
	fprintf(outfile, "\ncycle:\n->");
	/* i is created by reduce(a,b,c) below and indexes the ith
	convergent of the (U+sqrt(D))/V created there) */

	for (j = i; 1; j++) {
		TEMP1 = ADDI(F, UU);
		A = INT(TEMP1, VV);
		FREEMPI(TEMP1);

		TEMP1 = UU;
		TEMP2 = MULTI(A, VV);
		FREEMPI(A);
		UU = SUBI(TEMP2, UU);
		FREEMPI(TEMP1);
		FREEMPI(TEMP2);

		TEMP1 = VV;
		TEMP2 = MULTI(UU, UU);
		TEMP3 = SUBI(D, TEMP2);
		VV = INT(TEMP3, VV);
		FREEMPI(TEMP2);
		FREEMPI(TEMP3);

		TEMP2 = INT_(TEMP1, 2);
		X = MULT_I(TEMP2, (long)sign);
		FREEMPI(TEMP1);
		FREEMPI(TEMP2);
		if (j>i) {
			printf("->");
			fprintf(outfile, "->");
		}
		printf("(");
		fprintf(outfile, "(");
		PRINTI(X);
		FPRINTI(outfile, X);
		printf(",");
		fprintf(outfile, ",");
		FREEMPI(X);
		PRINTI(UU);
		FPRINTI(outfile, UU);
		printf(",");
		fprintf(outfile, ",");
		TEMP1 = INT_(VV, 2);
		Y = MULT_I(TEMP1, (long)(-sign));
		FREEMPI(TEMP1);
		PRINTI(Y);
		FPRINTI(outfile, Y);
		FREEMPI(Y);
		printf(")");
		fprintf(outfile, ")");

		sign = -sign;
		if (COMPAREI(UU, U) == 0 && COMPAREI(VV, V) == 0) {
			break;
		}
	}
	len = j + 1 - i;
	llen = len;

	if (len % 2 == 0) {
		FREEMPI(UU);
		FREEMPI(VV);
		FREEMPI(F);
		fclose(outfile);
		return(len);
	}
	else {
		printf("->");
		fprintf(outfile, "->");
		for (j = i; 1; j++) {
			TEMP1 = ADDI(F, UU);
			A = INT(TEMP1, VV);
			FREEMPI(TEMP1);

			TEMP1 = UU;
			TEMP2 = MULTI(A, VV);
			FREEMPI(A);
			UU = SUBI(TEMP2, UU);
			FREEMPI(TEMP1);
			FREEMPI(TEMP2);

			TEMP1 = VV;
			TEMP2 = MULTI(UU, UU);
			TEMP3 = SUBI(D, TEMP2);
			FREEMPI(TEMP2);
			VV = INT(TEMP3, VV);
			FREEMPI(TEMP3);

			TEMP2 = INT_(TEMP1, 2);
			X = MULT_I(TEMP2, (long)sign);
			FREEMPI(TEMP1);
			FREEMPI(TEMP2);

			if (j>i) {
				printf("->");
				fprintf(outfile, "->");
			}
			printf("(");
			fprintf(outfile, "(");
			PRINTI(X);
			FPRINTI(outfile, X);
			printf(",");
			fprintf(outfile, ",");
			FREEMPI(X);
			PRINTI(UU);
			FPRINTI(outfile, UU);
			printf(",");
			fprintf(outfile, ",");
			TEMP1 = INT_(VV, 2);
			Y = MULT_I(TEMP1, (long)(-sign));
			FREEMPI(TEMP1);
			PRINTI(Y);
			FPRINTI(outfile, Y);
			FREEMPI(Y);
			printf(")");
			fprintf(outfile, ")");

			sign = -sign;
			if (COMPAREI(UU, U) == 0 && COMPAREI(VV, V) == 0) {
				FREEMPI(UU);
				FREEMPI(VV);
				FREEMPI(F);
				printf("\n");
				fprintf(outfile, "\n");
				fclose(outfile);
				return(2 * len);
			}
		}
	}

}

USI REDUCE_POS(MPI *A, MPI *B, MPI *C) {
	MPI *D, *F, *ONE, *U, *V, *TEMP, *TEMP1, *TEMP2, *TEMP3, *TEMP4, *TEMP5;
	MPI *X, *Y, *AA, *BB, *CC;
	MPI *A_ORIG, *B_ORIG, *C_ORIG, *PP, *QQ, *RR, *SS, *HH;
	MPI *A_TMP, *B_TMP, *C_TMP, *T1, *T2, *T3, *T4;
	int r, s, sign;
	USI j, len, y;
	char buff[20];
	FILE *outfile;

	A_ORIG = COPYI(A);
	B_ORIG = COPYI(B);
	C_ORIG = COPYI(C);
	A_TMP = COPYI(A);
	B_TMP = COPYI(B);
	C_TMP = COPYI(C);
	PP = ONEI();
	QQ = ZEROI();
	RR = ZEROI();
	SS = ONEI();

	strcpy(buff, "reducepos.out");
	if (access(buff, R_OK) == 0)
		unlink(buff);
	outfile = fopen(buff, "w");
	ONE = ONEI();
	AA = COPYI(A);
	BB = COPYI(B);
	CC = COPYI(C);
	TEMP1 = MULTI(BB, BB);
	TEMP2 = MULTI(AA, CC);
	TEMP3 = MULT_I(TEMP2, 4);
	FREEMPI(TEMP2);
	D = SUBI(TEMP1, TEMP3);
	FREEMPI(TEMP1);
	FREEMPI(TEMP3);
	F = BIG_MTHROOT(D, 2);
	printf("(");
	fprintf(outfile, "(");
	PRINTI(AA);
	FPRINTI(outfile, AA);
	printf(",");
	fprintf(outfile, ",");
	PRINTI(BB);
	FPRINTI(outfile, BB);
	printf(",");
	fprintf(outfile, ",");
	PRINTI(CC);
	FPRINTI(outfile, CC);
	printf(")->");
	fprintf(outfile, ")->");
	sign = 1;
	U = COPYI(BB);
	TEMP1 = MULT_I(CC, 2);
	V = MULT_I(TEMP1, (long)(-sign));
	FREEMPI(TEMP1);
	for (j = 0; 1; j++) {
		if (j) {
			TEMP1 = MULTI(U, U);
			TEMP2 = SUBI(D, TEMP1);
			TEMP3 = MULT_I(V, 2);
			TEMP4 = INTI(TEMP2, TEMP3);
			X = MULT_I(TEMP4, (long)sign);
			A_TMP = COPYI(TEMP4);
			B_TMP = COPYI(U);
			FREEMPI(TEMP1);
			FREEMPI(TEMP2);
			FREEMPI(TEMP3);
			FREEMPI(TEMP4);

			TEMP1 = INT_(V, 2);
			TEMP2 = MINUSI(TEMP1);
			Y = MULT_I(TEMP2, (long)sign);
			C_TMP = COPYI(Y);
			FREEMPI(TEMP1);
			FREEMPI(TEMP2);
			if (j>1) {
				printf("->");
				fprintf(outfile, "->");
			}
			printf("(");
			fprintf(outfile, "(");
			PRINTI(X);
			FPRINTI(outfile, X);
			printf(",");
			fprintf(outfile, ",");
			PRINTI(U);
			FPRINTI(outfile, U);
			printf(",");
			fprintf(outfile, ",");
			PRINTI(Y);
			FPRINTI(outfile, Y);
			printf(")");
			fprintf(outfile, ")");
			FREEMPI(X);
			FREEMPI(Y);
		}
		sign = -sign;
		if (V->S >0 && U->S > 0 && COMPAREI(U, F) <= 0) {
			TEMP1 = ADDI(U, V);
			r = COMPAREI(F, TEMP1);
			FREEMPI(TEMP1);

			TEMP2 = SUBI(V, U);
			s = COMPAREI(TEMP2, F);
			FREEMPI(TEMP2);

			if (r < 0 && s <= 0) {
				printf("\nUnimodular matrix (");
				fprintf(outfile, "\nUnimodular matrix (");
				PRINTI(PP);
				FPRINTI(outfile, PP);
				printf(",");
				fprintf(outfile, ",");
				PRINTI(QQ);
				FPRINTI(outfile, QQ);
				printf(",");
				fprintf(outfile, ",");
				PRINTI(RR);
				FPRINTI(outfile, RR);
				printf(",");
				fprintf(outfile, ",");
				PRINTI(SS);
				FPRINTI(outfile, SS);
				printf(") transforms (");
				fprintf(outfile, ") transforms (");
				PRINTI(A_ORIG);
				FPRINTI(outfile, A_ORIG);
				printf(",");
				fprintf(outfile, ",");
				PRINTI(B_ORIG);
				FPRINTI(outfile, B_ORIG);
				printf(",");
				fprintf(outfile, ",");
				PRINTI(C_ORIG);
				FPRINTI(outfile, C_ORIG);
				printf(") to the reduced form (");
				fprintf(outfile, ") to the reduced form (");
				if (j) {
					PRINTI(A_TMP);
					FPRINTI(outfile, A_TMP);
					printf(",");
					fprintf(outfile, ",");
					PRINTI(B_TMP);
					FPRINTI(outfile, B_TMP);
					printf(",");
					fprintf(outfile, ",");
					PRINTI(C_TMP);
					FPRINTI(outfile, C_TMP);
					printf(")\n");
					fprintf(outfile, ")\n");
				}
				else {
					PRINTI(A_ORIG);
					FPRINTI(outfile, A_ORIG);
					printf(",");
					fprintf(outfile, ",");
					PRINTI(B_ORIG);
					FPRINTI(outfile, B_ORIG);
					printf(",");
					fprintf(outfile, ",");
					PRINTI(C_ORIG);
					FPRINTI(outfile, C_ORIG);
					printf(")\n");
					fprintf(outfile, ")\n");
				}
				FREEMPI(A_ORIG);
				FREEMPI(B_ORIG);
				FREEMPI(C_ORIG);
				FREEMPI(A_TMP);
				FREEMPI(B_TMP);
				FREEMPI(C_TMP);
				FREEMPI(PP);
				FREEMPI(QQ);
				FREEMPI(RR);
				FREEMPI(SS);
				fclose(outfile);
				len = CYCLE_PERIOD(D, U, V, j, sign);
				outfile = fopen(buff, "a");
				printf("cfrac has period length %u\n", llen);
				fprintf(outfile, "cfrac has period length %u\n", llen);
				printf("cycle of reduced forms has period length %u\n", len);
				fprintf(outfile, "cycle of reduced forms has period length %u\n", len);
				FREEMPI(ONE);
				FREEMPI(D);
				FREEMPI(F);
				FREEMPI(U);
				FREEMPI(V);
				FREEMPI(AA);
				FREEMPI(BB);
				FREEMPI(CC);
				fclose(outfile);
				return(len);
			}
		}
		FREEMPI(A_TMP);
		FREEMPI(B_TMP);
		FREEMPI(C_TMP);
		if (V->S > 0) {
			TEMP1 = ADDI(F, U);
			X = INTI(TEMP1, V);
			FREEMPI(TEMP1);
		}
		else {
			TEMP1 = ADDI(F, U);
			TEMP2 = ADDI(TEMP1, ONE);
			FREEMPI(TEMP1);
			X = INTI(TEMP2, V);
			FREEMPI(TEMP2);
		}
		T1 = MINUSI(QQ);
		T3 = MINUSI(SS);
		y = j % 2;
		if (y == 0) {
			HH = COPYI(X);
		}
		else {
			HH = MINUSI(X);
		}
		TEMP5 = MULTI(QQ, HH);
		T2 = ADDI(PP, TEMP5);
		FREEMPI(TEMP5);
		TEMP5 = MULTI(SS, HH);
		T4 = ADDI(RR, TEMP5);
		FREEMPI(TEMP5);
		FREEMPI(HH);
		TEMP = PP;
		PP = T1;
		FREEMPI(TEMP);
		TEMP = QQ;
		QQ = T2;
		FREEMPI(TEMP);
		TEMP = RR;
		RR = T3;
		FREEMPI(TEMP);
		TEMP = SS;
		SS = T4;
		FREEMPI(TEMP);

		TEMP1 = U;
		TEMP2 = MULTI(X, V);
		U = SUBI(TEMP2, U);
		FREEMPI(TEMP1);
		FREEMPI(TEMP2);
		FREEMPI(X);

		TEMP1 = V;
		TEMP2 = MULTI(U, U);
		TEMP3 = SUBI(D, TEMP2);
		V = INTI(TEMP3, V);
		FREEMPI(TEMP1);
		FREEMPI(TEMP2);
		FREEMPI(TEMP3);
	}
}

MPI *REDUCE_POSX(MPI *A, MPI *B, MPI *C) {
	MPI *D, *F, *G, *TEMP, *TEMP1, *TEMP2, *TEMP3;
	USI len;
	int t;

	TEMP1 = MULTI(B, B);
	TEMP2 = MULTI(A, C);
	TEMP3 = MULT_I(TEMP2, 4);
	FREEMPI(TEMP2);
	D = SUBI(TEMP1, TEMP3);
	FREEMPI(TEMP1);
	FREEMPI(TEMP3);

	F = BIG_MTHROOT(D, 2);
	G = MULTI(F, F);
	t = COMPAREI(G, D);
	FREEMPI(F);
	FREEMPI(G);
	if (t == 0) {
		printf("B^2-4AC is a perfect square\n");
		FREEMPI(D);
		return(NULL);
	}
	t = D->S;
	if (t <= 0) {
		printf("B^2-4*A*C <= 0\n");
		FREEMPI(D);
		return(NULL);
	}
	TEMP = POWER_I(10, 6);
	t = RSV(D, TEMP);
	FREEMPI(TEMP);
	if (t >= 0) {
		printf("D >= 10^6\n");
		FREEMPI(D);
		return(NULL);
	}
	FREEMPI(D);
	len = REDUCE_POS(A, B, C);
	TEMP1 = CHANGE((USL)len);
	return(TEMP1);
}

USI POS0(MPI *D) {
	USI a, b, e, t, z, parity_flag, class_number, x;
	USL d, f; /* d will be restricted to d < 10^6. */
	MPI *A, *B, *C, *TEMP, *F, *G, *H, *I, *J, *K, *L;
	MPI *ONE, *R;
	int sign, sign1;

	GLOBALARRAY_U = BUILDMPIA();
	GLOBALARRAY_V = BUILDMPIA();
	class_number = 0;
	ONE = ONEI();
	global_count = 0;
	parity_flag = 0;
	d = CONVERTI(D);

	if ((d - 1) % 4 != 0) {
		e = 2;
	}
	else {
		e = 1;
	}
	printf("Creating a complete list of reduced forms of discriminant %lu\n", d);
	F = BIG_MTHROOT(D, 2);
	f = CONVERTI(F);
	G = INT0_(F, 2);
	for (a = 1; a <= f; a++) {
		A = CHANGE(a);
		I = MULT_I(A, 2);
		J = MULT_I(A, 4);
		for (b = e; b <= f; b = b + 2) {
			B = CHANGE(b);
			TEMP = MULTI(B, B);
			H = SUBI(TEMP, D);
			FREEMPI(TEMP);
			TEMP = MOD(H, J);
			sign = TEMP->S;
			FREEMPI(TEMP);
			if (sign == 0) {
				TEMP = SUBI(F, I);
				sign1 = COMPAREI(TEMP, B);
				FREEMPI(TEMP);
				if ((RSV(A, G) <= 0 && sign1 < 0) || (RSV(A, G) > 0 && sign1 >= 0)) {
					R = INT(H, J);
					TEMP = ABSI(R);
					K = GCD(A, B);
					L = GCD(K, TEMP);
					FREEMPI(K);
					FREEMPI(TEMP);
					t = EQONEI(L);
					FREEMPI(L);
					if (t) {
						TEMP = MINUSI(H);
						C = INT0(TEMP, J); /* C > 0 */
						FREEMPI(TEMP);
						TEMP = MULT_I(C, 2);
						x = CHECK_UVARRAYS(B, TEMP);
						if (x) {
							class_number = class_number + 1;
							z = PERIOD(D, B, TEMP);
							if (z % 2) {
								parity_flag = 1;
							}
							printf("[%u]: (", class_number);
							PRINTI(A);
							printf(", ");
							PRINTI(B);
							printf(", ");
							PRINTI(R);
							printf(")\n");
						}
						FREEMPI(TEMP);
						FREEMPI(C);
					}
					FREEMPI(R);
				}
			}
			FREEMPI(H);
			FREEMPI(B);
		}
		FREEMPI(A);
		FREEMPI(I);
		FREEMPI(J);
	}

	if (parity_flag) {
		printf("x^2-%lu*y^2=-4 has a solution\n", d);
	}
	else {
		printf("x^2-%lu*y^2=-4 has no solution\n", d);
		class_number = 2 * class_number;
		printf("There are reduced forms (-a,b,-c) as well\n");
	}
	printf("h(%lu)=%u\n", d, class_number);
	FREEMPI(ONE);
	FREEMPI(F);
	FREEMPI(G);
	FREEMPIA(GLOBALARRAY_U);
	FREEMPIA(GLOBALARRAY_V);
	return(class_number);
}

MPI *POS0X(MPI *D) {
	USI class_number;
	USL d;
	int t;
	MPI *ONE, *N, *TEMP, *F, *G;

	d = CONVERTI(D);
	ONE = ONEI();
	t = COMPAREI(D, ONE);
	FREEMPI(ONE);
	if (t <= 0) {
		printf(" D <= 1\n");
		return(NULL);
	}
	TEMP = POWER_I(10, 6);
	t = RSV(D, TEMP);
	FREEMPI(TEMP);
	if (t >= 0) {
		printf("D >= 10^6\n");
		return(NULL);
	}
	F = BIG_MTHROOT(D, 2);
	G = MULTI(F, F);
	t = COMPAREI(G, D);
	FREEMPI(F);
	FREEMPI(G);
	if (t == 0) {
		printf("B^2-4AC is a perfect square\n");
		return(NULL);
	}
	else {
		class_number = POS0(D);
		N = CHANGE(class_number);
		return(N);
	}
}

USL TABLENEG(MPI *M, MPI *N) {
	/* The following gives a table of h(-d) for all squarefree d
	* in the range M<=d<=N<10^6.
	*/
	MPI *DD, *ONE, *ZERO;
	USL count, m, n, d, r;
	USI h, t;
	long dd;
	char buff[20];
	FILE *outfile;

	count = 0;
	m = CONVERTI(M);
	n = CONVERTI(N);

	strcpy(buff, "tableneg.out");
	if (access(buff, R_OK) == 0)
		unlink(buff);
	outfile = fopen(buff, "a");

	for (d = m; d <= n; d++) {
		h = SQUAREFREE_TEST1000(d);
		if (h == 0) {
			continue;
		}
		else {
			r = d % 4;
			dd = (long)d;
			if (r != 3) {
				DD = CHANGEL((long)(-4) * d);
			}
			else {

				DD = CHANGEL(-dd);
			}
			ONE = ONEI();
			ZERO = ZEROI();
			t = NEG(DD, ONE, ZERO);
			FREEMPI(DD);
			FREEMPI(ONE);
			FREEMPI(ZERO);
			count++;
			printf("h(%ld) = %u\n", -d, t);
			fprintf(outfile, "h(%ld) = %u\n", -d, t);
		}
	}
	if (count == 0) {
		printf("no squarefree d in the range [%lu,%lu]\n", m, n);
	}
	fclose(outfile);
	return (count);
}

MPI *TABLENEGX(MPI *M, MPI *N) {
	MPI *LIMIT, *COUNT;
	USL count;

	LIMIT = POWER_I(10, 6);

	if (COMPAREI(M, N) > 0) {
		printf("M > N\n");
		FREEMPI(LIMIT);
		return (NULL);
	}
	if (COMPAREI(N, LIMIT) >= 0) {
		printf("N >= 10^6\n");
		FREEMPI(LIMIT);
		return (NULL);
	}
	FREEMPI(LIMIT);
	if (M->S <= 0 || N->S <= 0) {
		printf("M < 1 or N < 1\n");
		return (NULL);
	}
	count = TABLENEG(M, N);
	COUNT = CHANGE(count);
	return(COUNT);
}

USI POS1(MPI *D, int *norm) {
	/* D is squarefree. This function performs Lagrange's method on all
	* reduced quadratic irrationals (b+\sqrt(Disc))/2|c|, where 4*c divides
	* Disc-b^2, Disc being the discriminant. The class-number h(D) of Q(sqrt(D)
	* is calculated, as well as the norm of the fundamental unit.
	* For use in TABLEPOS(M,N) below.
	*/

	USI a, b, e, t, z, parity_flag, class_number, x;
	USL d, f; /* d will be restricted to d < 10^6. */
	MPI *A, *B, *C, *TEMP, *F, *G, *H, *I, *J, *K, *L;
	MPI *ONE, *R, *DD;
	int sign, sign1;

	GLOBALARRAY_U = BUILDMPIA();
	GLOBALARRAY_V = BUILDMPIA();
	class_number = 0;
	ONE = ONEI();
	global_count = 0;
	parity_flag = 0;
	d = CONVERTI(D);
	/* creates a fundamental discriminant */
	if ((d - 1) % 4 != 0) {
		d = 4 * d;
		e = 2;
	}
	else {
		e = 1;
	}
	DD = CHANGE(d);
	F = BIG_MTHROOT(DD, 2);
	f = CONVERTI(F);
	G = INT0_(F, 2);
	for (a = 1; a <= f; a++) {
		A = CHANGE(a);
		I = MULT_I(A, 2);
		J = MULT_I(A, 4);
		for (b = e; b <= f; b = b + 2) {
			B = CHANGE(b);
			TEMP = MULTI(B, B);
			H = SUBI(TEMP, DD);
			FREEMPI(TEMP);
			TEMP = MOD(H, J);
			sign = TEMP->S;
			FREEMPI(TEMP);
			if (sign == 0) {
				TEMP = SUBI(F, I);
				sign1 = COMPAREI(TEMP, B);
				FREEMPI(TEMP);
				if ((RSV(A, G) <= 0 && sign1 < 0) || (RSV(A, G) > 0 && sign1 >= 0)) {
					R = INT(H, J);
					TEMP = ABSI(R);
					K = GCD(A, B);
					L = GCD(K, TEMP);
					FREEMPI(K);
					FREEMPI(TEMP);
					t = EQONEI(L);
					FREEMPI(L);
					if (t) {
						TEMP = MINUSI(H);
						C = INT0(TEMP, J); /* C > 0 */
						FREEMPI(TEMP);
						TEMP = MULT_I(C, 2);
						x = CHECK_UVARRAYS(B, TEMP);
						if (x) {
							class_number = class_number + 1;
							z = PERIOD(DD, B, TEMP);
							if (z % 2) {
								parity_flag = 1;
							}
						}
						FREEMPI(TEMP);
						FREEMPI(C);
					}
					FREEMPI(R);
				}
			}
			FREEMPI(H);
			FREEMPI(B);
		}
		FREEMPI(A);
		FREEMPI(I);
		FREEMPI(J);
	}
	if (e == 2) {
		d = d / 4;
	}
	if (parity_flag) {
		*norm = -1;
	}
	else {
		*norm = 1;
	}
	/*printf("h(%lu)=%u, N(eta)=%d\n", d, class_number, *norm);*/
	FREEMPI(ONE);
	FREEMPI(F);
	FREEMPI(G);
	FREEMPI(DD);
	FREEMPIA(GLOBALARRAY_U);
	FREEMPIA(GLOBALARRAY_V);
	return(class_number);
}

USL TABLEPOS(MPI *M, MPI *N) {
	/* The following gives a table of h(d) for all squarefree d
	* in the range 2<=M<=d<=N<10^6.
	* The number of squarefree d encountered is returned.
	*/
	MPI *DD;
	USL count, m, n, d;
	USI h;
	char buff[20];
	FILE *outfile;
	int norm;

	count = 0;
	m = CONVERTI(M);
	n = CONVERTI(N);

	strcpy(buff, "tablepos.out");
	if (access(buff, R_OK) == 0)
		unlink(buff);
	outfile = fopen(buff, "a");

	for (d = m; d <= n; d++) {
		h = SQUAREFREE_TEST1000(d);
		if (h == 0) {
			continue;
		}
		else {
			DD = CHANGE(d);
			h = POS1(DD, &norm);
			FREEMPI(DD);
			count++;
			printf("h(%lu) = %u, N(eta)= %d\n", d, h, norm);
			fprintf(outfile, "h(%lu) = %u, N(eta)= %d\n", d, h, norm);
		}
	}
	if (count == 0) {
		printf("no squarefree d in the range [%lu,%lu]\n", m, n);
	}
	fclose(outfile);
	return (count);
}

MPI *TABLEPOSX(MPI *M, MPI *N) {
	MPI *LIMIT, *COUNT, *ONE;
	USL count;
	int t;

	LIMIT = POWER_I(10, 6);

	if (COMPAREI(M, N) > 0) {
		printf("M > N\n");
		FREEMPI(LIMIT);
		return (NULL);
	}
	if (COMPAREI(N, LIMIT) >= 0) {
		printf("N >= 10^6\n");
		FREEMPI(LIMIT);
		return (NULL);
	}
	FREEMPI(LIMIT);
	ONE = ONEI();
	t = COMPAREI(M, ONE);
	FREEMPI(ONE);
	if (t <= 0) {
		printf("M <= 1\n");
		return (NULL);
	}
	count = TABLEPOS(M, N);
	COUNT = CHANGE(count);
	return(COUNT);
}

USI REDUCE_NEG0(MPI *A, MPI *B, MPI *C, MPI **AA, MPI **BB, MPI **CC, MPI **alpha, MPI **beta, MPI **gamma, MPI **delta, USI print_flag) {
	/* This is Gauss' algorithm for reducing a positive definite binary
	* quadratic form (A,B,C).
	* See L.E. Dickson, Introduction to the theory of numbers,
	* page 69. Here d=B^2-4AC <0, A>0,C>0, while the reduced form (AA,BB,CC)
	* satisfies -AA<BB<=AA, CC>=AA, with BB>=0 if CC=AA.
	* The number of steps taken in the reduction is returned.
	*/
	MPI *D, *X, *Y, *U, *V, *TEMP1, *TEMP2, *TEMP3;
	MPI *a, *b, *c, *TEMPX, *TEMPY, *TEMPU, *TEMPV, *MINUSD, *DELTA;
	MPI *TEMPB;
	USI i;
	int r, s, t;

	i = 0;
	TEMP1 = MULTI(B, B);
	TEMP2 = MULTI(A, C);
	TEMP3 = MULT_I(TEMP2, 4);
	D = SUBI(TEMP1, TEMP3);
	FREEMPI(TEMP1);
	FREEMPI(TEMP2);
	FREEMPI(TEMP3);
	TEMP1 = MINUSI(A);
	r = COMPAREI(TEMP1, B);
	FREEMPI(TEMP1);
	s = COMPAREI(B, A);
	t = COMPAREI(A, C);
	if (r < 0 && s <= 0 && t <= 0) {
		r = COMPAREI(A, C);
		if (t < 0 || (r == 0 && B->S >= 0)) {
			if (print_flag) {
				printf("(");
				PRINTI(A);
				printf(",");
				PRINTI(B);
				printf(",");
				PRINTI(C);
				printf(") is reduced\n");
				printf("Transforming matrix:1,0,0,1\n");
			}
			FREEMPI(D);
			*alpha = ONEI();
			*beta = ZEROI();
			*gamma = ZEROI();
			*delta = ONEI();
			*AA = COPYI(A);
			*BB = COPYI(B);
			*CC = COPYI(C);
			return (0);
		}
	}
	X = ONEI();
	Y = ZEROI();
	U = ZEROI();
	V = ONEI();
	MINUSD = MINUSI(D);
	a = COPYI(A);
	b = COPYI(B);
	c = COPYI(C);
	while (1) {
		TEMPB = b;
		TEMP1 = MINUSI(b);
		TEMP2 = MULT_I(c, 2);
		b = HALFMOD(TEMP1, TEMP2);
		FREEMPI(TEMP1);
		TEMP1 = ADDI(TEMPB, b);
		FREEMPI(TEMPB);
		TEMP3 = MINUSI(TEMP1);
		FREEMPI(TEMP1);
		DELTA = INT(TEMP3, TEMP2);
		FREEMPI(TEMP2);
		FREEMPI(TEMP3);
		TEMPU = U;
		TEMPX = X;
		TEMPY = Y;
		TEMPV = V;
		X = MINUSI(Y);
		Y = MULTABC(TEMPX, Y, DELTA);
		FREEMPI(TEMPX);
		FREEMPI(TEMPY);
		U = MINUSI(V);
		V = MULTABC(TEMPU, V, DELTA);
		FREEMPI(TEMPU);
		FREEMPI(TEMPV);
		FREEMPI(DELTA);
		TEMP1 = a;
		a = COPYI(c);
		FREEMPI(TEMP1);
		TEMP1 = c;
		TEMP2 = MULTABC(MINUSD, b, b);
		TEMP3 = MULT_I(a, 4);
		c = INT(TEMP2, TEMP3);
		FREEMPI(TEMP1);
		FREEMPI(TEMP2);
		FREEMPI(TEMP3);
		i++;
		if (print_flag) {
			printf("(");
			PRINTI(a);
			printf(",");
			PRINTI(b);
			printf(",");
			PRINTI(c);
			printf(")\n");
		}
		if (COMPAREI(c, a) >= 0) {
			break;
		}
	}
	if (COMPAREI(a, c) == 0 && b->S < 0) {
		TEMP1 = b;
		b = MINUSI(b);
		FREEMPI(TEMP1);
		TEMPX = X;
		TEMPY = Y;
		X = COPYI(Y);
		Y = MINUSI(TEMPX);
		FREEMPI(TEMPX);
		FREEMPI(TEMPY);
		TEMPU = U;
		TEMPV = V;
		U = COPYI(V);
		V = MINUSI(TEMPU);
		FREEMPI(TEMPU);
		FREEMPI(TEMPV);
	}
	*alpha = COPYI(X);
	*beta = COPYI(Y);
	*gamma = COPYI(U);
	*delta = COPYI(V);
	*AA = COPYI(a);
	*BB = COPYI(b);
	*CC = COPYI(c);
	if (print_flag) {
		printf("(");
		PRINTI(a);
		printf(",");
		PRINTI(b);
		printf(",");
		PRINTI(c);
		printf(") is reduced\n");
		printf("Transforming matrix: ");
		PRINTI(X);
		printf(",");
		PRINTI(Y);
		printf(",");
		PRINTI(U);
		printf(",");
		PRINTI(V);
		printf("\n");
	}
	FREEMPI(MINUSD);
	FREEMPI(D);
	FREEMPI(X);
	FREEMPI(Y);
	FREEMPI(U);
	FREEMPI(V);
	FREEMPI(a);
	FREEMPI(b);
	FREEMPI(c);
	return(i);
}

USI REP_DEFINITE(MPI *A, MPI *B, MPI *C, MPI *M, USI print_flag) {
	/* Given a positive definite binary quadratic form Ax^2+Bxy+Cy^2,
	* we use an algorithm of Gauss to determine if a given positive integer M
	* is expressible as M = AX^2+BXY+CY^2, X and Y integers, gcd(X,Y) = 1.
	* Note: Here D = B^2 - 4AC < 0, A > 0, C > 0.
	* See L.E. Dickson, Introduction to the theory of numbers, pages 74-75.
	*/
	USI solution_number;
	MPI *AA1, *BB1, *CC1, *R1, *S1, *T1, *U1, *D;
	MPI *AA2, *BB2, *CC2, *R2, *S2, *T2, *U2;
	MPI *TEMP0, *TEMP1, *TEMP2, *FOURM, *S, *MODULUS, *TWOM;
	MPI *L, *X, *Y, *XX, *YY, *T, *U;
	MPI *ALPHA1, *BETA1, *GAMMA1, *DELTA1, *N, *M0, *TEMP3, *TEMP;
	MPIA SOLUTION, SOL;
	USI l, i, r, s, t, numbr, m0, j, k;
	FILE *outfile;
	char buff[20];

	solution_number = 0;
	REDUCE_NEG0(A, B, C, &AA1, &BB1, &CC1, &R1, &S1, &T1, &U1, 0);
	TEMP0 = MULTI(B, B);
	TEMP1 = MULTI(A, C);
	TEMP2 = MULT_I(TEMP1, 4);
	D = SUBI(TEMP0, TEMP2);
	FREEMPI(TEMP0);
	FREEMPI(TEMP1);
	FREEMPI(TEMP2);
	FOURM = MULT_I(M, 4);
	/* First solve x^2=D (mod FOURM). */
	S = SQROOT(D, FOURM, &SOLUTION, &MODULUS, &l);
	if (EQMINUSONEI(S)) { /* cleanup of nettbytes */
		FREEMPI(S);
		FREEMPIA(SOLUTION);
		FREEMPI(MODULUS);
		printf("x^2=");
		PRINTI(D);
		FREEMPI(D);
		printf("(mod ");
		PRINTI(FOURM);
		FREEMPI(FOURM);
		printf(") not soluble.\n");
		printf("Hence no primitive solutions.\n");
		return(0);
	}
	FREEMPI(S);
	M0 = INT0(FOURM, MODULUS);
	if (M0->D > 1) {
		execerror("M0 >= 2^32", "");
	}
	m0 = CONVERTI(M0); /* valid, as 0 < M0 < 2^32 */
	TEMP = MULT_I(M0, l);
	FREEMPI(M0);
	TEMP1 = INT0_(TEMP, 100);
	FREEMPI(TEMP);
	t = TEMP1->D;
	FREEMPI(TEMP1);
	if (t) {
		execerror("number of positive solutions of X^2=N(mod FOURM) > 100", "");
	}
	numbr = 0;
	SOL = BUILDMPIA();
	ADD_TO_MPIA(SOL, NULL, 0);
	for (j = 0; j<l; j++) {
		for (k = 0; k<m0; k++) {
			TEMP1 = MULT_I(MODULUS, k);
			TEMP2 = ADD0I(SOLUTION->A[j], TEMP1);
			TEMP3 = MULT_I(TEMP2, 2);
			if (COMPAREI(TEMP3, FOURM) <= 0) {
				ADD_TO_MPIA(SOL, TEMP2, numbr);
				numbr++;
			}
			else {
				TEMP0 = SUB0I(FOURM, TEMP2);
				ADD_TO_MPIA(SOL, TEMP0, numbr);
				FREEMPI(TEMP0);
				numbr++;
			}
			FREEMPI(TEMP2);
			FREEMPI(TEMP3);
			FREEMPI(TEMP1);
		}

	}
	FREEMPI(MODULUS);
	FREEMPIA(SOLUTION);

	strcpy(buff, "repdefinite.out");
	outfile = fopen(buff, "w");

	TWOM = MULT_I(M, 2);
	for (i = 0; i < numbr; i++) {
		N = COPYI(SOL->A[i]);
		/*if(RSV(N, TWOM) > 0){
		TEMP = N;
		N = SUBI(FOURM, N);
		FREEMPI(TEMP);
		}*/
		fprintf(outfile, "N = ");
		FPRINTI(outfile, N);
		fprintf(outfile, "\n");
		if (EQUALI(N, TWOM)) {
			FREEMPI(N);
			continue;
		}
		/* now calculate L, where N^2-4ML=D */
		TEMP1 = MULTI(N, N);
		TEMP2 = SUBI(TEMP1, D);
		L = INT0(TEMP2, FOURM);
		FREEMPI(TEMP1);
		FREEMPI(TEMP2);
		REDUCE_NEG0(M, N, L, &AA2, &BB2, &CC2, &R2, &S2, &T2, &U2, 0);
		r = EQUALI(AA1, AA2);
		s = EQUALI(BB1, BB2);
		t = EQUALI(CC1, CC2);
		FREEMPI(AA2);
		FREEMPI(BB2);
		FREEMPI(CC2);
		if (r && s && t) {
			TEMP1 = MULTI(R1, U2);
			TEMP2 = MULTI(S1, T2);
			ALPHA1 = SUBI(TEMP1, TEMP2);
			FREEMPI(TEMP1);
			FREEMPI(TEMP2);

			TEMP1 = MULTI(S1, R2);
			TEMP2 = MULTI(S2, R1);
			BETA1 = SUBI(TEMP1, TEMP2);
			FREEMPI(TEMP1);
			FREEMPI(TEMP2);

			TEMP1 = MULTI(U2, T1);
			TEMP2 = MULTI(U1, T2);
			GAMMA1 = SUBI(TEMP1, TEMP2);
			FREEMPI(TEMP1);
			FREEMPI(TEMP2);

			TEMP1 = MULTI(U1, R2);
			TEMP2 = MULTI(S2, T1);
			DELTA1 = SUBI(TEMP1, TEMP2);
			FREEMPI(TEMP1);
			FREEMPI(TEMP2);

			FREEMPI(R2);
			FREEMPI(S2);
			FREEMPI(T2);
			FREEMPI(U2);
			if (print_flag) {
				fprintf(outfile, "The unimodular transformation\n");
				fprintf(outfile, "x = ");
				FPRINTI(outfile, ALPHA1);
				fprintf(outfile, "X + ");
				FPRINTI(outfile, BETA1);
				fprintf(outfile, "Y\n");
				fprintf(outfile, "y = ");
				FPRINTI(outfile, GAMMA1);
				fprintf(outfile, "X + ");
				FPRINTI(outfile, DELTA1);
				fprintf(outfile, "Y\n");
				fprintf(outfile, "converts (");
				FPRINTI(outfile, A);
				fprintf(outfile, ",");
				FPRINTI(outfile, B);
				fprintf(outfile, ",");
				FPRINTI(outfile, C);
				fprintf(outfile, ") to ");
				fprintf(outfile, "(");
				FPRINTI(outfile, M);
				fprintf(outfile, ",");
				FPRINTI(outfile, N);
				fprintf(outfile, ",");
				FPRINTI(outfile, L);
				fprintf(outfile, ")\n");
			}
			FREEMPI(N);
			FREEMPI(L);
			fprintf(outfile, "solution: (");
			FPRINTI(outfile, ALPHA1);
			fprintf(outfile, ",");
			FPRINTI(outfile, GAMMA1);
			fprintf(outfile, ")\n");

			X = MINUSI(ALPHA1);
			Y = MINUSI(GAMMA1);
			fprintf(outfile, "solution: (");
			FPRINTI(outfile, X);
			fprintf(outfile, ",");
			FPRINTI(outfile, Y);
			fprintf(outfile, ")\n");
			solution_number = solution_number + 2;

			if (D->D == 0 && D->S == -1 && D->V[0] == 4) { /* D = -4 */
				T = ZEROI();
				U = ONEI();
				AUTOMORPH(A, B, D, T, U, X, Y, &XX, &YY);
				fprintf(outfile, "solution: (");
				FPRINTI(outfile, XX);
				fprintf(outfile, ",");
				FPRINTI(outfile, YY);
				fprintf(outfile, ")\n");
				FREEMPI(XX);
				FREEMPI(YY);
				AUTOMORPH(A, B, D, T, U, ALPHA1, GAMMA1, &XX, &YY);
				fprintf(outfile, "solution: (");
				FPRINTI(outfile, XX);
				fprintf(outfile, ",");
				FPRINTI(outfile, YY);
				fprintf(outfile, ")\n");
				FREEMPI(XX);
				FREEMPI(YY);
				FREEMPI(T);
				FREEMPI(U);
				solution_number = solution_number + 2;
			}
			if (D->D == 0 && D->S == -1 && D->V[0] == 3) { /* D = -3 */
				T = ONEI();
				U = ONEI();
				AUTOMORPH(A, B, D, T, U, X, Y, &XX, &YY);
				fprintf(outfile, "solution: (");
				FPRINTI(outfile, XX);
				fprintf(outfile, ",");
				FPRINTI(outfile, YY);
				fprintf(outfile, ")\n");
				FREEMPI(XX);
				FREEMPI(YY);
				AUTOMORPH(A, B, D, T, U, ALPHA1, GAMMA1, &XX, &YY);
				FREEMPI(T);
				fprintf(outfile, "solution: (");
				FPRINTI(outfile, XX);
				fprintf(outfile, ",");
				FPRINTI(outfile, YY);
				fprintf(outfile, ")\n");
				FREEMPI(XX);
				FREEMPI(YY);
				T = MINUS_ONEI();
				AUTOMORPH(A, B, D, T, U, X, Y, &XX, &YY);
				fprintf(outfile, "solution: (");
				FPRINTI(outfile, XX);
				fprintf(outfile, ",");
				FPRINTI(outfile, YY);
				fprintf(outfile, ")\n");
				FREEMPI(XX);
				FREEMPI(YY);
				AUTOMORPH(A, B, D, T, U, ALPHA1, GAMMA1, &XX, &YY);
				fprintf(outfile, "solution: (");
				FPRINTI(outfile, XX);
				fprintf(outfile, ",");
				FPRINTI(outfile, YY);
				fprintf(outfile, ")\n");
				FREEMPI(XX);
				FREEMPI(YY);
				FREEMPI(T);
				FREEMPI(U);
				solution_number = solution_number + 4;
			}
			FREEMPI(X);
			FREEMPI(Y);
			fprintf(outfile, "\n");
			FREEMPI(ALPHA1);
			FREEMPI(BETA1);
			FREEMPI(GAMMA1);
			FREEMPI(DELTA1);
		}
		else {
			continue;
		}
	}
	FREEMPI(R1);
	FREEMPI(S1);
	FREEMPI(T1);
	FREEMPI(U1);
	FREEMPI(AA1);
	FREEMPI(BB1);
	FREEMPI(CC1);
	FREEMPI(TWOM);
	FREEMPI(FOURM);
	FREEMPI(D);
	FREEMPIA(SOL);
	fclose(outfile);
	return(solution_number);
}

void AUTOMORPH(MPI *A, MPI *B, MPI *D, MPI *T, MPI *U, MPI *X, MPI *Y, MPI **XX, MPI **YY) {
	/* From Loo-Keng Hua, Introduction to Number Theory,
	* Theorem 4.2, pages 279-281 */
	MPI *TEMP0, *TEMP1, *TEMP2, *TEMP3, *TEMP4, *TEMP5, *TEMP6;
	MPI *T1, *T2, *T3, *T4, *T5, *T6, *TEMP;

	TEMP1 = MULT_I(A, 2);
	TEMP2 = MULTI(TEMP1, X);
	TEMP0 = MULTI(Y, B);
	TEMP3 = ADDI(TEMP2, TEMP0);
	TEMP4 = MULTI(TEMP3, U);
	TEMP5 = MULTI(Y, T);
	TEMP6 = ADDI(TEMP4, TEMP5);
	*YY = INT_(TEMP6, 2);

	T1 = MULTI(TEMP3, T);
	T2 = MULTI(D, Y);
	TEMP = T2;
	T2 = MULTI(T2, U);
	FREEMPI(TEMP);
	T3 = ADDI(T1, T2);
	TEMP = T3;
	T3 = INT_(T3, 2);
	FREEMPI(TEMP);
	T4 = MULTI(B, *YY);
	T5 = SUBI(T3, T4);
	T6 = MULT_I(A, 2);
	*XX = INT(T5, T6);

	FREEMPI(TEMP0);
	FREEMPI(TEMP1);
	FREEMPI(TEMP2);
	FREEMPI(TEMP3);
	FREEMPI(TEMP4);
	FREEMPI(TEMP5);
	FREEMPI(TEMP6);
	FREEMPI(T1);
	FREEMPI(T2);
	FREEMPI(T3);
	FREEMPI(T4);
	FREEMPI(T5);
	FREEMPI(T6);
}

MPI *REP_DEFINITEX(MPI *A, MPI *B, MPI *C, MPI *M, MPI *PRINT_FLAG) {
	MPI *D, *TEMP1, *TEMP2, *TEMP3;
	USI i, print_flag;

	if (M->S <= 0) {
		printf("M <= 0\n");
		return(NULL);
	}
	if (A->S <= 0) {
		printf("A <= 0\n");
		return(NULL);
	}
	if (C->S <= 0) {
		printf("C <= 0\n");
		return(NULL);
	}
	TEMP1 = MULTI(B, B);
	TEMP2 = MULTI(A, C);
	TEMP3 = MULT_I(TEMP2, 4);
	D = SUBI(TEMP1, TEMP3);
	FREEMPI(TEMP1);
	FREEMPI(TEMP2);
	FREEMPI(TEMP3);
	if (D->S >= 0) {
		printf("B^2-4*A*C >= 0\n");
		FREEMPI(D);
		return(NULL);
	}
	print_flag = (USI)CONVERTI(PRINT_FLAG);
	if (print_flag != 0 && print_flag != 1) {
		printf("print_flag != 0 or 1\n");
		FREEMPI(D);
		return(NULL);
	}
	i = REP_DEFINITE(A, B, C, M, print_flag);
	FREEMPI(D);
	return(CHANGE(i));
}
