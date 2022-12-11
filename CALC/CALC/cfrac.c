/* cfrac.c */
/* programs for continued fractions. */

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
#include <math.h>
#include <time.h>

static USL c1 = 0, c2 = 0, c3 = 0, c4 = 0, c5 = 0, c6 = 0;
/* This function finds the NICF of sqrt(D) and exits
when it reaches one of the 6 conditions of the HCW-PAB paper */
USI MIDPT(MPI *D, MPI *FLAG) {
	MPI *A, *X, *G, *TEMP1, *TEMP2, *TEMP3, *Z1, *Z2, *F, *P, *Q, *OP, *OQ;
	USI t, k;
	int s;

	k = 1;
	X = BIG_MTHROOT(D, 2);
	G = MULTI(X, X);
	t = EQUALI(G, D);
	FREEMPI(G);
	if (t) {
		FREEMPI(X);
		return(0);
	}
	TEMP1 = MULT_II(X, 2);
	TEMP2 = ADD0_I(TEMP1, 1);
	Z1 = MULTI(TEMP2, TEMP2);
	FREEMPI(TEMP1);
	FREEMPI(TEMP2);
	Z2 = MULT_II(D, 4);
	s = COMPAREI(Z2, Z1);
	if (s == 1) {
		F = ONEI();
	}
	else {
		F = ZEROI();
	}
	FREEMPI(Z1);
	FREEMPI(Z2);

	P = ZEROI();
	Q = ONEI();
	while (1) {
		A = NINT(P, Q, X, F);
		OP = COPYI(P);
		OQ = COPYI(Q);
		TEMP1 = MULTI(A, Q);
		FREEMPI(A);
		TEMP2 = P;
		P = SUBI(TEMP1, P);
		FREEMPI(TEMP1);
		FREEMPI(TEMP2);
		TEMP1 = MULTI(P, P);
		TEMP2 = SUBI(TEMP1, D);
		FREEMPI(TEMP1);
		TEMP3 = Q;
		Q = INTI(TEMP2, Q);
		FREEMPI(TEMP2);
		FREEMPI(TEMP3);

		if (EQUALI(OP, P)) {
			if (FLAG->S) {
				printf("P[%u]=", k);
				PRINTI(P);
				printf("\n");
			}
			c1++;
			FREEMPI(OP);
			FREEMPI(OQ);
			break;
		}
		TEMP1 = SUBI(OP, P);
		TEMP2 = ABSI(TEMP1);
		FREEMPI(TEMP1);
		TEMP1 = ABSI(OQ);
		t = EQUALI(TEMP2, TEMP1);
		FREEMPI(TEMP2);
		if (t) {
			if (FLAG->S) {
				printf("|Q[%u]|=", k);
				PRINTI(TEMP1);
				FREEMPI(TEMP1);
				printf("\n");
			}
			c2++;
			FREEMPI(OP);
			FREEMPI(OQ);
			break;
		}
		else {
			FREEMPI(TEMP1);
		}
		if (EQUALI(OQ, Q)) {
			if (FLAG->S) {
				printf("Q[%u]=", k - 1);
				PRINTI(OQ);
				printf("\n");
				printf("Q[%u]=", k);
				PRINTI(Q);
				printf("\n");
			}
			c3++;
			FREEMPI(OP);
			FREEMPI(OQ);
			break;
		}
		TEMP1 = MINUSI(Q);
		t = EQUALI(OQ, TEMP1);
		FREEMPI(TEMP1);
		if (t) {
			if (FLAG->S) {
				printf("Q[%u]=", k - 1);
				PRINTI(OQ);
				printf("\n");
				printf("Q[%u]=", k);
				PRINTI(Q);
				printf("\n");
			}
			c4++;
			FREEMPI(OP);
			FREEMPI(OQ);
			break;
		}
		if ((Q->V[0]) % 2 == 0) {
			TEMP1 = INT_(Q, 2);
			TEMP2 = ADDI(OQ, TEMP1);
			FREEMPI(TEMP1);
			TEMP3 = ABSI(TEMP2);
			FREEMPI(TEMP2);
			t = EQUALI(P, TEMP3);
			FREEMPI(TEMP3);
			if (t && (Q->S == OQ->S)) {
				if (FLAG->S) {
					printf("P[%u]=", k + 1);
					PRINTI(P);
					printf("\n");
				}
				c5++;
				FREEMPI(OP);
				FREEMPI(OQ);
				break;
			}
		}
		if ((OQ->V[0]) % 2 == 0) {
			TEMP1 = INT_(OQ, 2);
			TEMP2 = ADDI(Q, TEMP1);
			FREEMPI(TEMP1);
			TEMP3 = ABSI(TEMP2);
			FREEMPI(TEMP2);
			t = EQUALI(P, TEMP3);

			FREEMPI(TEMP3);
			if (t && (Q->S == OQ->S)) {
				if (FLAG->S) {
					printf("P[%u]=", k);
					PRINTI(P);
					printf("\n");
				}
				c6++;
				FREEMPI(OP);
				FREEMPI(OQ);
				break;
			}
		}
		FREEMPI(OP);
		FREEMPI(OQ);
		k = k + 1;
	}
	FREEMPI(F);
	FREEMPI(X);
	FREEMPI(P);
	FREEMPI(Q);
	return(1);
}


/* NINT(P,Q,X,F) finds the nearest integer A to (P+sqrt(D)/Q".
* Here X=int(sqrt(D)) and F = 1 if sqrt(D)-X > 1/2, 0 if sqrt(D)-X < 1/2.
* See http://www.numbertheory.org/notes.html.
* For use in midpt(D) above.
*/

MPI *NINT(MPI *P, MPI *Q, MPI *X, MPI *F) {
	MPI *S, *T, *A, *TEMP1, *TEMP2, *TEMP3, *TEMP4, *QQ;
	int r;

	r = Q->S;
	if (Q->S < 0) {
		QQ = MINUSI(Q);
	}
	else {
		QQ = COPYI(Q);
	}

	if (r>0) {
		TEMP1 = ADDI(QQ, F);
		TEMP2 = INT_(TEMP1, 2);
		FREEMPI(TEMP1);
		TEMP3 = ADDI(X, TEMP2);
		FREEMPI(TEMP2);
		T = ADDI(P, TEMP3);
		FREEMPI(TEMP3);
		A = INTI(T, Q);
		FREEMPI(T);
	}
	else {
		TEMP1 = ONEI();
		TEMP2 = ADDI(QQ, TEMP1);
		S = ADDI(F, TEMP2);
		FREEMPI(TEMP2);
		TEMP2 = INT_(S, 2);
		TEMP3 = ADDI(P, X);
		TEMP4 = ADDI(TEMP3, TEMP2);
		FREEMPI(TEMP2);
		FREEMPI(TEMP3);
		if ((S->V[0]) % 2) {
			T = ADDI(TEMP4, TEMP1);
		}
		else {
			T = COPYI(TEMP4);
		}
		FREEMPI(TEMP4);
		FREEMPI(S);
		TEMP2 = INTI(T, Q);
		A = ADDI(TEMP2, TEMP1);
		FREEMPI(T);
		FREEMPI(TEMP1);
		FREEMPI(TEMP2);
	}
	FREEMPI(QQ);
	return(A);
}

/* NINT0(P,Q,X,F) finds the nearest integer A to (P+sqrt(D)/Q", Q > 0.
* Here X=int(sqrt(D)) and F = 1 if sqrt(D)-X > 1/2, 0 if sqrt(D)-X < 1/2.
* See http://www.numbertheory.org/notes.html.
* For use in midpt(D) above.
*/

MPI *NINT0(MPI *P, MPI *Q, MPI *X, MPI *F) {
	MPI *T, *A, *TEMP1, *TEMP2, *TEMP3;

	TEMP1 = ADD0I(Q, F);
	TEMP2 = INT0_(TEMP1, 2);
	FREEMPI(TEMP1);
	TEMP3 = ADD0I(X, TEMP2);
	FREEMPI(TEMP2);
	T = ADD0I(P, TEMP3);
	FREEMPI(TEMP3);
	A = INT0(T, Q);
	FREEMPI(T);
	return(A);
}

void MIDPT_COUNT(MPI *M, MPI *N) {
	MPI *D, *TEMP;
	USL d, m, n;
	USI t;
	if (EQUALI(M, N)) {
		TEMP = ONEI();
		t = MIDPT(M, TEMP);
		FREEMPI(TEMP);
	}
	else {
		m = CONVERTI(M);
		n = CONVERTI(N);
		TEMP = ZEROI();
		for (d = m; d <= n; d++) {
			D = CHANGEL(d);
			t = MIDPT(D, TEMP);
			FREEMPI(D);
			if (d % 1000 == 0) {
				printf("d=%lu\n", d);
			}
		}
		FREEMPI(TEMP);
	}
	printf("c1=%lu\n", c1);
	printf("c2=%lu\n", c2);
	printf("c3=%lu\n", c3);
	printf("c4=%lu\n", c4);
	printf("c5=%lu\n", c5);
	printf("c6=%lu\n", c6);
	c1 = 0;
	c2 = 0;
	c3 = 0;
	c4 = 0;
	c5 = 0;
	c6 = 0;
	return;
}

/* This function finds the period of the NICF-P of sqrt(D) and exits
when it reaches one of the 5 conditions of the nicf.pdf paper */
MPI *NICF_PERIOD(MPI *D) {
	MPI *X, *G, *Z1, *Z2, *F, *TEMP1, *TEMP2, *C, *T;
	MPI *P1, *P2, *Q1, *Q2, *TEMP, *P, *Q, *OLDP, *OLDQ, *TEMP3, *TEMP4;
	MPI *ROLDQ, *PERIOD, *T2, *T3, *H;
	USL h, k, t1;
	USI t, olda, flag4 = 0, flag5 = 0;
	int a, s;

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
	TEMP1 = MULT_II(X, 2);
	TEMP2 = ADD0_I(TEMP1, 1);
	Z1 = MULTI(TEMP2, TEMP2);
	FREEMPI(TEMP1);
	FREEMPI(TEMP2);
	Z2 = MULT_II(D, 4);
	s = COMPAREI(Z2, Z1);
	if (s == 1) {
		F = ONEI();
	}
	else {
		F = ZEROI();
	}
	FREEMPI(Z1);
	FREEMPI(Z2);

	P = ZEROI();
	Q = ONEI();
	h = 0;
	k = 0;

	a = 1;
	while (1) {
		olda = a;
		TEMP = ADDI(P, X);
		C = INT0(TEMP, Q);
		FREEMPI(TEMP);
		TEMP = MULTI(C, Q);
		P1 = SUB0I(TEMP, P);
		FREEMPI(TEMP);
		TEMP = MULTI(P1, P1);
		Q1 = SUB0I(D, TEMP);
		FREEMPI(TEMP);
		P2 = ADD0I(P1, Q);
		TEMP = MULTI(P2, P2);
		Q2 = SUB0I(TEMP, D);
		FREEMPI(TEMP);
		OLDP = P;
		OLDQ = Q;
		T = NINT0(P, Q, X, F);
		if (EQUALI(T, C)) {
			a = 1;
			P = P1;
			FREEMPI(P2);
			Q = INT0(Q1, Q);
		}
		else {
			a = -1;
			P = P2;
			Q = INT0(Q2, Q);
			FREEMPI(P1);
		}
		FREEMPI(Q1);
		FREEMPI(Q2);
		FREEMPI(C);
		FREEMPI(T);

		if (EQUALI(P, OLDP)) { /* period-length = k = 2*h, RCF period p is even */
							   /*       printf("P[h] = P[h+1]\n");*/
			break;
		}

		TEMP = ADD0I(OLDP, OLDQ);
		t = EQUALI(P, TEMP);
		FREEMPI(TEMP);
		if (t) { /* period-length = k = 2*h, RCF period p is even */
				 /*     printf("P[h+1] = P[h] + Q[h]\n");*/
			break;
		}

		if (EQUALI(Q, OLDQ)) { /* period-length = 2h+1, p even if a = -1, p odd if a = 1 */
							   /*    printf("Q[h] = Q[h+1]\n");*/
			flag4 = 1;
			break;
		}
		t1 = MOD0_(Q, 2);
		TEMP3 = INT0_(Q, 2);
		TEMP4 = ADD0I(OLDQ, TEMP3);
		FREEMPI(TEMP3);
		t = EQUALI(P, TEMP4);
		FREEMPI(TEMP4);
		if (t1 == 0 && t && a == -1) { /* period-length = 2h+1, p odd */
									   /*printf("P[h+1] = Q[h] + Q[h+1]/2 and epsilon[h+1] = -1\n");*/
			flag5 = 1;
			break;
		}
		if (h != 0 || k != 0) {
			TEMP3 = MULTI(OLDP, OLDP);
			TEMP4 = SUBI(D, TEMP3);
			FREEMPI(TEMP3);
			ROLDQ = INT(TEMP4, OLDQ);
			FREEMPI(TEMP4);
			if (olda == -1) {
				ROLDQ->S = -(ROLDQ->S);
			}
			t1 = MOD0_(ROLDQ, 2);
			TEMP3 = INT0_(ROLDQ, 2);
			TEMP4 = ADD0I(OLDQ, TEMP3);
			FREEMPI(TEMP3);
			FREEMPI(ROLDQ);
			t = EQUALI(OLDP, TEMP4);
			FREEMPI(TEMP4);
			if (t1 == 0 && t && olda == -1) { /* period-length = 2h, p odd */
											  /*printf("P[h] = Q[h] + Q[h-1]/2 and epsilon[h] = -1\n");*/
				break;
			}
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
		FREEMPI(OLDP);
		FREEMPI(OLDQ);
	}
	FREEMPI(OLDP);
	FREEMPI(OLDQ);
	FREEMPI(P);
	FREEMPI(Q);
	FREEMPI(F);
	FREEMPI(X);
	T2 = CHANGE(2147483648UL);
	T3 = CHANGE(h);
	TEMP1 = MULT_I(T2, k);
	H = ADD0I(TEMP1, T3); /*H = k*2147483648 + h */
						  /*PRINT_MPI(H, "H");*/
	FREEMPI(T2);
	FREEMPI(T3);
	FREEMPI(TEMP1);
	PERIOD = MULT_II(H, 2);
	FREEMPI(H);
	if (flag4 || flag5) {
		TEMP = PERIOD;
		PERIOD = ADD0_I(PERIOD, 1);
		FREEMPI(TEMP);
	}
	return(PERIOD);
}

/* This function finds the period and smallest solution of Pell's equation using NICF-P of sqrt(D) and exits
when it reaches one of the 5 conditions of the nicf.pdf paper.
If Eptr = 1, the midpoint and its type are printed. */
MPI *NICF_PERIOD0(MPI *D, MPI *Eptr, MPI **Xptr, MPI **Yptr)
{
	MPI *X, *G, *Z1, *Z2, *F, *TEMP1, *TEMP2, *C, *T;
	MPI *P1, *P2, *Q1, *Q2, *TEMP, *P, *Q, *OLDP, *OLDQ, *TEMP3, *TEMP4;
	MPI *ROLDQ, *PERIOD, *T2, *T3, *H;
	MPI *B1, *C1, *B2, *C2, *OLDB1, *OLDC1, *OLDB2, *OLDC2, *B;
	USL h, k, t1;
	USI t, olda, flag4 = 0, flag5 = 0;
	int a, s;

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
	TEMP1 = MULT_II(X, 2);
	TEMP2 = ADD0_I(TEMP1, 1);
	Z1 = MULTI(TEMP2, TEMP2);
	FREEMPI(TEMP1);
	FREEMPI(TEMP2);
	Z2 = MULT_II(D, 4);
	s = COMPAREI(Z2, Z1);
	if (s == 1) {
		F = ONEI();
	}
	else {
		F = ZEROI();
	}
	FREEMPI(Z1);
	FREEMPI(Z2);

	P = ZEROI();
	Q = ONEI();
	h = 0;
	k = 0;

	B1 = ZEROI();
	C1 = ONEI();
	B2 = ONEI();
	C2 = ZEROI();
	a = 1;
	while (1) {
		OLDC1 = C1;
		OLDC2 = C2;
		OLDB1 = B1;
		OLDB2 = B2;
		olda = a;
		TEMP = ADD0I(P, X);
		C = INT0(TEMP, Q);
		FREEMPI(TEMP);
		TEMP = MULTI(C, Q);
		P1 = SUB0I(TEMP, P);
		FREEMPI(TEMP);
		TEMP = MULTI(P1, P1);
		Q1 = SUB0I(D, TEMP);
		FREEMPI(TEMP);
		P2 = ADD0I(P1, Q);
		TEMP = MULTI(P2, P2);
		Q2 = SUB0I(TEMP, D);
		FREEMPI(TEMP);
		OLDP = P;
		OLDQ = Q;
		T = NINT0(P, Q, X, F);
		if (EQUALI(T, C)) {
			a = 1;
			P = P1;
			FREEMPI(P2);
			Q = INT0(Q1, Q);
			B = C;
		}
		else {
			a = -1;
			P = P2;
			FREEMPI(P1);
			Q = INT0(Q2, Q);
			B = ADD0_I(C, 1);
			FREEMPI(C);
		}
		FREEMPI(Q1);
		FREEMPI(Q2);
		FREEMPI(T);

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
		if (EQUALI(P, OLDP)) { /* period-length = k = 2*h, RCF period p is even */
			if (Eptr->S) {
				printf("P[h] = P[h+1]\n");
				printf("x^2-dy^2=1\n");
			}
			if (olda == 1) {
				*Xptr = MULTAB_PLUS_CDI(OLDC2, C1, OLDC1, OLDB2);
				*Yptr = MULTAB_PLUS_CDI(OLDC2, C2, OLDC2, OLDB2);
			}
			else {
				*Xptr = MULTAB_MINUS_CDI(OLDC2, C1, OLDC1, OLDB2);
				*Yptr = MULTAB_MINUS_CDI(OLDC2, C2, OLDC2, OLDB2);
			}
			break;
		}

		TEMP = ADD0I(OLDP, OLDQ);
		t = EQUALI(P, TEMP);
		FREEMPI(TEMP);
		if (t) { /* period-length = k = 2*h, RCF period p is even */
			if (Eptr->S) {
				printf("P[h+1] = P[h] + Q[h]\n");
				printf("x^2-dy^2=1\n");
			}
			TEMP = SUBI(OLDB2, OLDC2);
			*Xptr = MULTAB_PLUS_CDI(OLDC2, C1, OLDC1, TEMP);
			FREEMPI(TEMP);
			TEMP = SUBI(OLDB2, OLDC2);
			TEMP1 = ADDI(C2, TEMP);
			FREEMPI(TEMP);
			*Yptr = MULTI(OLDC2, TEMP1);
			FREEMPI(TEMP1);
			break;
		}

		if (EQUALI(Q, OLDQ)) { /* period-length = 2h+1, p even if a = -1, p odd if a = 1 */
			if (Eptr->S) {
				printf("Q[h] = Q[h+1]\n");
				printf("x^2-dy^2=%d\n", -a);
			}
			flag4 = 1;
			if (a == 1) {
				*Xptr = MULTAB_PLUS_CDI(C1, C2, OLDC1, OLDC2);
				*Yptr = MULTAB_PLUS_CDI(C2, C2, OLDC2, OLDC2);
			}
			else {
				*Xptr = MULTAB_MINUS_CDI(C1, C2, OLDC1, OLDC2);
				*Yptr = MULTAB_MINUS_CDI(C2, C2, OLDC2, OLDC2);
			}
			break;
		}
		t1 = MOD0_(Q, 2);
		TEMP3 = INT0_(Q, 2);
		TEMP4 = ADD0I(OLDQ, TEMP3);
		FREEMPI(TEMP3);
		t = EQUALI(P, TEMP4);
		FREEMPI(TEMP4);
		if (t1 == 0 && t && a == -1) { /* period-length = 2h+1, p odd */
			if (Eptr->S) {
				printf("P[h+1] = Q[h] + Q[h+1]/2 and epsilon[h+1] = -1\n");
				printf("x^2-dy^2= -1\n");
			}
			flag5 = 1;
			TEMP = MULTI(C1, C2);
			TEMP1 = MULT_I(OLDC1, 2);
			TEMP2 = MULTAB_PLUS_CDI(C2, OLDC1, C1, OLDC2);
			TEMP3 = MULTI(TEMP1, OLDC2);
			FREEMPI(TEMP1);
			TEMP1 = ADDI(TEMP, TEMP3);
			FREEMPI(TEMP);
			*Xptr = SUBI(TEMP1, TEMP2);
			FREEMPI(TEMP1);
			FREEMPI(TEMP2);
			FREEMPI(TEMP3);

			TEMP = MULTI(C2, C2);
			TEMP1 = MULT_I(OLDC2, 2);
			TEMP2 = SUBI(OLDC2, C2);
			TEMP3 = MULTI(TEMP1, TEMP2);
			*Yptr = ADDI(TEMP, TEMP3);
			FREEMPI(TEMP);
			FREEMPI(TEMP1);
			FREEMPI(TEMP2);
			FREEMPI(TEMP3);
			break;
		}
		if (h != 0 || k != 0) {
			TEMP3 = MULTI(OLDP, OLDP);
			TEMP4 = SUBI(D, TEMP3);
			FREEMPI(TEMP3);
			ROLDQ = INT(TEMP4, OLDQ);
			FREEMPI(TEMP4);
			if (olda == -1) {
				ROLDQ->S = -(ROLDQ->S);
			}
			t1 = MOD0_(ROLDQ, 2);
			TEMP3 = INT0_(ROLDQ, 2);
			TEMP4 = ADD0I(OLDQ, TEMP3);
			FREEMPI(TEMP3);
			FREEMPI(ROLDQ);
			t = EQUALI(OLDP, TEMP4);
			FREEMPI(TEMP4);
			if (t1 == 0 && t && olda == -1) { /* period-length = 2h, p odd */
				if (Eptr->S) {
					printf("P[h] = Q[h] + Q[h-1]/2 and epsilon[h] = -1\n");
					printf("x^2-dy^2= -1\n");
				}
				TEMP = MULT_II(OLDC1, 2);
				TEMP2 = MULTAB_PLUS_CDI(TEMP, OLDC2, OLDB1, OLDB2);
				TEMP1 = MULTAB_PLUS_CDI(OLDC2, OLDB1, OLDC1, OLDB2);
				*Xptr = SUBI(TEMP2, TEMP1);
				FREEMPI(TEMP);
				FREEMPI(TEMP1);
				FREEMPI(TEMP2);

				TEMP = MULTI(OLDB2, OLDB2);
				TEMP1 = SUBI(OLDC2, OLDB2);
				TEMP2 = MULT_II(OLDC2, 2);
				TEMP3 = MULTI(TEMP2, TEMP1);
				*Yptr = ADD0I(TEMP, TEMP3);
				FREEMPI(TEMP);
				FREEMPI(TEMP1);
				FREEMPI(TEMP2);
				FREEMPI(TEMP3);
				break;
			}
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
		FREEMPI(OLDP);
		FREEMPI(OLDQ);
	}
	FREEMPI(OLDP);
	FREEMPI(OLDQ);
	FREEMPI(P);
	FREEMPI(Q);
	FREEMPI(F);
	FREEMPI(X);
	FREEMPI(B1);
	FREEMPI(B2);
	FREEMPI(C1);
	FREEMPI(C2);
	FREEMPI(OLDB1);
	FREEMPI(OLDB2);
	T2 = CHANGE(2147483648UL);
	T3 = CHANGE(h);
	TEMP1 = MULT_I(T2, k);
	H = ADD0I(TEMP1, T3); /*H = k*2147483648 + h */
	if (Eptr->S) {
		PRINT_MPI(H, "h");
	}
	FREEMPI(T2);
	FREEMPI(T3);
	FREEMPI(TEMP1);
	PERIOD = MULT_II(H, 2);
	FREEMPI(H);
	if (flag4 || flag5) {
		TEMP = PERIOD;
		PERIOD = ADD0_I(PERIOD, 1);
		FREEMPI(TEMP);
	}
	return(PERIOD);
}

void TIME_COUNT_NSCF(MPI *M, MPI *N) {
	MPI *D, *X, *Y, *PERIOD, *ZERO, *TEMP;
	int t;

	ZERO = ZEROI();
	D = COPYI(M);
	while (RSV(D, N) <= 0) {
		t = SQUARETEST(D);
		if (t) {
			PERIOD = NSCF_PERIOD0(D, ZERO, &X, &Y);
			FREEMPI(PERIOD);
			FREEMPI(X);
			FREEMPI(Y);
		}
		TEMP = D;
		D = ADD0_I(D, 1);
		FREEMPI(TEMP);
	}
	FREEMPI(ZERO);
	FREEMPI(D);
	return;
}

void TIME_COUNT_NICF(MPI *M, MPI *N) {
	MPI *D, *X, *Y, *PERIOD, *ZERO, *TEMP;
	int t;

	ZERO = ZEROI();
	D = COPYI(M);
	while (RSV(D, N) <= 0) {
		t = SQUARETEST(D);
		if (t) {
			PERIOD = NICF_PERIOD0(D, ZERO, &X, &Y);
			FREEMPI(PERIOD);
			FREEMPI(X);
			FREEMPI(Y);
		}
		TEMP = D;
		D = ADD0_I(D, 1);
		FREEMPI(TEMP);
	}
	FREEMPI(ZERO);
	FREEMPI(D);
	return;
}

void TIME_COUNT_RCF(MPI *M, MPI *N) {
	MPI *D, *X, *Y, *PERIOD, *ZERO, *TEMP;
	int t;

	ZERO = ZEROI();
	D = COPYI(M);
	while (RSV(D, N) <= 0) {
		t = SQUARETEST(D);
		if (t) {
			PERIOD = RCF_PERIOD0(D, ZERO, &X, &Y);
			FREEMPI(PERIOD);
			FREEMPI(X);
			FREEMPI(Y);
		}
		TEMP = D;
		D = ADD0_I(D, 1);
		FREEMPI(TEMP);
	}
	FREEMPI(ZERO);
	FREEMPI(D);
	return;
}

void PERIOD_TIME_COUNT_RCF(MPI *M, MPI *N) {
	MPI *D, *PERIOD, *TEMP;
	int t;
	clock_t time_1, time_2;

	D = COPYI(M);
	time_1 = clock();
	while (RSV(D, N) <= 0) {
		t = SQUARETEST(D);
		if (t) {
			PERIOD = CFRAC_PERIOD(D);
			PRINT_MPI(D, "D");
			PRINT_MPI(PERIOD, "PERIOD");
			if (EQZEROI(PERIOD)) {
				FREEMPI(D);
				return;
			}
			else {
				FREEMPI(PERIOD);
			}
		}
		TEMP = D;
		D = ADD0_I(D, 1);
		FREEMPI(TEMP);
	}
	FREEMPI(D);
	time_2 = clock();
	printf("RCF:CPU time taken=%f seconds\n", (float)(time_2 - time_1) / (float)CLOCKS_PER_SEC);
	return;
}

void PERIOD_TIME_COUNT_NSCF(MPI *M, MPI *N) {
	MPI *D, *PERIOD, *TEMP;
	int t;
	clock_t time_1, time_2;

	D = COPYI(M);
	time_1 = clock();
	while (RSV(D, N) <= 0) {
		t = SQUARETEST(D);
		if (t) {
			PERIOD = NSCF_PERIOD(D, NULL);
			PRINT_MPI(D, "D");
			PRINT_MPI(PERIOD, "PERIOD");
			if (EQZEROI(PERIOD)) {
				FREEMPI(D);
				return;
			}
			else {
				FREEMPI(PERIOD);
			}
		}
		TEMP = D;
		D = ADD0_I(D, 1);
		FREEMPI(TEMP);
	}
	FREEMPI(D);
	time_2 = clock();
	printf("NSCF:CPU time taken=%f seconds\n", (float)(time_2 - time_1) / (float)CLOCKS_PER_SEC);
	return;
}

void PERIOD_TIME_COUNT_NICF(MPI *M, MPI *N) {
	MPI *D, *PERIOD, *TEMP;
	int t;
	clock_t time_1, time_2;

	D = COPYI(M);
	time_1 = clock();
	while (RSV(D, N) <= 0) {
		t = SQUARETEST(D);
		if (t) {
			PERIOD = NICF_PERIOD(D);
			PRINT_MPI(D, "D");
			PRINT_MPI(PERIOD, "PERIOD");
			if (EQZEROI(PERIOD)) {
				FREEMPI(D);
				return;
			}
			else {
				FREEMPI(PERIOD);
			}
		}
		TEMP = D;
		D = ADD0_I(D, 1);
		FREEMPI(TEMP);
	}
	FREEMPI(D);
	time_2 = clock();
	printf("NICF:CPU time taken=%f seconds\n", (float)(time_2 - time_1) / (float)CLOCKS_PER_SEC);
	return;
}

void TIME_COUNT_RCF00(MPI *M, MPI *N) {
	MPI *D, *X, *Y, *TEMP;
	clock_t time_1, time_2;

	D = COPYI(M);
	time_1 = clock();
	while (RSV(D, N) <= 0) {
		RCF_PERIOD00(D, &X, &Y);
		FREEMPI(X);
		FREEMPI(Y);
		TEMP = D;
		D = ADD0_I(D, 1);
		FREEMPI(TEMP);
	}
	FREEMPI(D);
	time_2 = clock();
	printf("RCF:CPU time taken=%f seconds\n", (float)(time_2 - time_1) / (float)CLOCKS_PER_SEC);
	return;
}

void TIME_COUNT_NICF00(MPI *M, MPI *N) {
	MPI *D, *X, *Y, *TEMP;
	clock_t time_1, time_2;

	D = COPYI(M);
	time_1 = clock();
	while (RSV(D, N) <= 0) {
		NICF_PERIOD00(D, &X, &Y);
		FREEMPI(X);
		FREEMPI(Y);
		TEMP = D;
		D = ADD0_I(D, 1);
		FREEMPI(TEMP);
	}
	FREEMPI(D);
	time_2 = clock();
	printf("NICF:CPU time taken=%f seconds\n", (float)(time_2 - time_1) / (float)CLOCKS_PER_SEC);
	return;
}

/* This function finds smallest solution of Pell's equation using NICF-P of sqrt(D) and exits
when it reaches one of the 5 conditions of the nicf.pdf paper.*/
void NICF_PERIOD00(MPI *D, MPI **Xptr, MPI **Yptr)
{
	MPI *X, *G, *Z1, *Z2, *F, *TEMP1, *TEMP2, *C, *T;
	MPI *P1, *P2, *Q1, *Q2, *TEMP, *P, *Q, *OLDP, *OLDQ, *TEMP3, *TEMP4;
	MPI *ROLDQ;
	MPI *B1, *C1, *B2, *C2, *OLDB1, *OLDC1, *OLDB2, *OLDC2, *B;
	USL t1;
	USI flag, t, olda, flag4 = 0, flag5 = 0;
	int a, s;

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
	t = EQUALI(D, TEMP);
	FREEMPI(TEMP);
	FREEMPI(G);
	if (t) {/* period 1, exceptional case */
		*Xptr = X;
		*Yptr = ONEI();
		return;
	}
	TEMP1 = MULT_II(X, 2);
	TEMP2 = ADD0_I(TEMP1, 1);
	Z1 = MULTI(TEMP2, TEMP2);
	FREEMPI(TEMP1);
	FREEMPI(TEMP2);
	Z2 = MULT_II(D, 4);
	s = COMPAREI(Z2, Z1);
	if (s == 1) {
		F = ONEI();
	}
	else {
		F = ZEROI();
	}
	FREEMPI(Z1);
	FREEMPI(Z2);

	P = ZEROI();
	Q = ONEI();

	B1 = ZEROI();
	C1 = ONEI();
	B2 = ONEI();
	C2 = ZEROI();
	a = 1;
	flag = 0;
	while (1) {
		OLDC1 = C1;
		OLDC2 = C2;
		OLDB1 = B1;
		OLDB2 = B2;
		olda = a;
		TEMP = ADD0I(P, X);
		C = INT0(TEMP, Q);
		FREEMPI(TEMP);
		TEMP = MULTI(C, Q);
		P1 = SUB0I(TEMP, P);
		FREEMPI(TEMP);
		TEMP = MULTI(P1, P1);
		Q1 = SUB0I(D, TEMP);
		FREEMPI(TEMP);
		P2 = ADD0I(P1, Q);
		TEMP = MULTI(P2, P2);
		Q2 = SUB0I(TEMP, D);
		FREEMPI(TEMP);
		OLDP = P;
		OLDQ = Q;
		T = NINT0(P, Q, X, F);
		if (EQUALI(T, C)) {
			a = 1;
			P = P1;
			FREEMPI(P2);
			Q = INT0(Q1, Q);
			B = C;
		}
		else {
			a = -1;
			P = P2;
			FREEMPI(P1);
			Q = INT0(Q2, Q);
			B = ADD0_I(C, 1);
			FREEMPI(C);
		}
		FREEMPI(Q1);
		FREEMPI(Q2);
		FREEMPI(T);

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
		if (EQUALI(P, OLDP)) { /* period-length = k = 2*h, RCF period p is even */
			if (olda == 1) {
				*Xptr = MULTAB_PLUS_CDI(OLDC2, C1, OLDC1, OLDB2);
				*Yptr = MULTAB_PLUS_CDI(OLDC2, C2, OLDC2, OLDB2);
			}
			else {
				*Xptr = MULTAB_MINUS_CDI(OLDC2, C1, OLDC1, OLDB2);
				*Yptr = MULTAB_MINUS_CDI(OLDC2, C2, OLDC2, OLDB2);
			}
			break;
		}

		TEMP = ADD0I(OLDP, OLDQ);
		t = EQUALI(P, TEMP);
		FREEMPI(TEMP);
		if (t) { /* period-length = k = 2*h, RCF period p is even */
			TEMP = SUBI(OLDB2, OLDC2);
			*Xptr = MULTAB_PLUS_CDI(OLDC2, C1, OLDC1, TEMP);
			FREEMPI(TEMP);
			TEMP = SUBI(OLDB2, OLDC2);
			TEMP1 = ADDI(C2, TEMP);
			FREEMPI(TEMP);
			*Yptr = MULTI(OLDC2, TEMP1);
			FREEMPI(TEMP1);
			break;
		}

		if (EQUALI(Q, OLDQ)) { /* period-length = 2h+1, p even if a = -1, p odd if a = 1 */
			flag4 = 1;
			if (a == 1) {
				*Xptr = MULTAB_PLUS_CDI(C1, C2, OLDC1, OLDC2);
				*Yptr = MULTAB_PLUS_CDI(C2, C2, OLDC2, OLDC2);
			}
			else {
				*Xptr = MULTAB_MINUS_CDI(C1, C2, OLDC1, OLDC2);
				*Yptr = MULTAB_MINUS_CDI(C2, C2, OLDC2, OLDC2);
			}
			break;
		}
		t1 = MOD0_(Q, 2);
		TEMP3 = INT0_(Q, 2);
		TEMP4 = ADD0I(OLDQ, TEMP3);
		FREEMPI(TEMP3);
		t = EQUALI(P, TEMP4);
		FREEMPI(TEMP4);
		if (t1 == 0 && t && a == -1) { /* period-length = 2h+1, p odd */
			flag5 = 1;
			TEMP = MULTI(C1, C2);
			TEMP1 = MULT_I(OLDC1, 2);
			TEMP2 = MULTAB_PLUS_CDI(C2, OLDC1, C1, OLDC2);
			TEMP3 = MULTI(TEMP1, OLDC2);
			FREEMPI(TEMP1);
			TEMP1 = ADDI(TEMP, TEMP3);
			FREEMPI(TEMP);
			*Xptr = SUBI(TEMP1, TEMP2);
			FREEMPI(TEMP1);
			FREEMPI(TEMP2);
			FREEMPI(TEMP3);

			TEMP = MULTI(C2, C2);
			TEMP1 = MULT_I(OLDC2, 2);
			TEMP2 = SUBI(OLDC2, C2);
			TEMP3 = MULTI(TEMP1, TEMP2);
			*Yptr = ADDI(TEMP, TEMP3);
			FREEMPI(TEMP);
			FREEMPI(TEMP1);
			FREEMPI(TEMP2);
			FREEMPI(TEMP3);
			break;
		}
		if (flag) {
			TEMP3 = MULTI(OLDP, OLDP);
			TEMP4 = SUBI(D, TEMP3);
			FREEMPI(TEMP3);
			ROLDQ = INT(TEMP4, OLDQ);
			FREEMPI(TEMP4);
			if (olda == -1) {
				ROLDQ->S = -ROLDQ->S;
			}
			t1 = MOD0_(ROLDQ, 2);
			TEMP3 = INT0_(ROLDQ, 2);
			TEMP4 = ADD0I(OLDQ, TEMP3);
			FREEMPI(TEMP3);
			FREEMPI(ROLDQ);
			t = EQUALI(OLDP, TEMP4);
			FREEMPI(TEMP4);
			if (t1 == 0 && t && olda == -1) { /* period-length = 2h, p odd */
				TEMP = MULT_II(OLDC1, 2);
				TEMP2 = MULTAB_PLUS_CDI(TEMP, OLDC2, OLDB1, OLDB2);
				TEMP1 = MULTAB_PLUS_CDI(OLDC2, OLDB1, OLDC1, OLDB2);
				*Xptr = SUBI(TEMP2, TEMP1);
				FREEMPI(TEMP);
				FREEMPI(TEMP1);
				FREEMPI(TEMP2);

				TEMP = MULTI(OLDB2, OLDB2);
				TEMP1 = SUBI(OLDC2, OLDB2);
				TEMP2 = MULT_II(OLDC2, 2);
				TEMP3 = MULTI(TEMP2, TEMP1);
				*Yptr = ADDI(TEMP, TEMP3);
				FREEMPI(TEMP);
				FREEMPI(TEMP1);
				FREEMPI(TEMP2);
				FREEMPI(TEMP3);
				break;
			}
		}
		FREEMPI(OLDB1);
		FREEMPI(OLDB2);
		FREEMPI(OLDP);
		FREEMPI(OLDQ);
	}
	FREEMPI(OLDP);
	FREEMPI(OLDQ);
	FREEMPI(P);
	FREEMPI(Q);
	FREEMPI(F);
	FREEMPI(X);
	FREEMPI(B1);
	FREEMPI(B2);
	FREEMPI(C1);
	FREEMPI(C2);
	FREEMPI(OLDB1);
	FREEMPI(OLDB2);
	return;
}

void TIME_COUNT_NSCF00(MPI *M, MPI *N) {
	MPI *D, *X, *Y, *TEMP;
	clock_t time_1, time_2;

	D = COPYI(M);
	time_1 = clock();
	while (RSV(D, N) <= 0) {
		NSCF_PERIOD00(D, &X, &Y);
		FREEMPI(X);
		FREEMPI(Y);
		TEMP = D;
		D = ADD0_I(D, 1);
		FREEMPI(TEMP);
	}
	FREEMPI(D);
	time_2 = clock();
	printf("NSCF:CPU time taken=%f seconds\n", (float)(time_2 - time_1) / (float)CLOCKS_PER_SEC);
	return;
}

void TIME_COUNT_RCF000(MPI *M, MPI *N) {
	MPI *D, *X, *Y, *TEMP;
	clock_t time_1, time_2;

	D = COPYI(M);
	time_1 = clock();
	while (RSV(D, N) <= 0) {
		RCF_PERIOD000(D, &X, &Y);
		FREEMPI(X);
		FREEMPI(Y);
		TEMP = D;
		D = ADD0_I(D, 1);
		FREEMPI(TEMP);
	}
	FREEMPI(D);
	time_2 = clock();
	printf("RCF:CPU time taken=%f seconds\n", (float)(time_2 - time_1) / (float)CLOCKS_PER_SEC);
	return;
}

void TIME_COUNT_NICF000(MPI *M, MPI *N) {
	MPI *D, *X, *Y, *TEMP;
	clock_t time_1, time_2;

	D = COPYI(M);
	time_1 = clock();
	while (RSV(D, N) <= 0) {
		NICF_PERIOD000(D, &X, &Y);
		FREEMPI(X);
		FREEMPI(Y);
		FREEMPI(D);
		TEMP = D;
		D = ADD0_I(D, 1);
		FREEMPI(TEMP);
	}
	FREEMPI(D);
	time_2 = clock();
	printf("NICF:CPU time taken=%f seconds\n", (float)(time_2 - time_1) / (float)CLOCKS_PER_SEC);
	return;
}

void TIME_COUNT_NSCF000(MPI *M, MPI *N) {
	MPI *D, *X, *Y, *TEMP;
	clock_t time_1, time_2;

	D = COPYI(M);
	time_1 = clock();
	while (RSV(D, N) <= 0) {
		NSCF_PERIOD000(D, &X, &Y);
		FREEMPI(X);
		FREEMPI(Y);
		TEMP = D;
		D = ADD0_I(D, 1);
		FREEMPI(TEMP);
	}
	FREEMPI(D);
	time_2 = clock();
	printf("NSCF:CPU time taken=%f seconds\n", (float)(time_2 - time_1) / (float)CLOCKS_PER_SEC);
	return;
}

/* This function finds smallest solution of Pell's equation using NICF-P of sqrt(D) and exits
when it reaches one of the 5 conditions of the nicf.pdf paper.*/
void NICF_PERIOD000(MPI *D, MPI **Xptr, MPI **Yptr)
{
	MPI *X, *G, *Z1, *Z2, *F, *TEMP1, *TEMP2, *C, *T;
	MPI *P1, *P2, *Q1, *Q2, *TEMP, *P, *Q, *OLDP, *OLDQ, *TEMP3, *TEMP4;
	MPI *ROLDQ;
	MPI *B2, *C2, *OLDB2, *OLDC2, *B;
	USL t1;
	USI flag, t, olda, flag4 = 0, flag5 = 0;
	int a, s;

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
	t = EQUALI(D, TEMP);
	FREEMPI(TEMP);
	FREEMPI(G);
	if (t) {/* period 1, exceptional case */
		*Xptr = X;
		*Yptr = ONEI();
		return;
	}
	TEMP1 = MULT_II(X, 2);
	TEMP2 = ADD0_I(TEMP1, 1);
	Z1 = MULTI(TEMP2, TEMP2);
	FREEMPI(TEMP1);
	FREEMPI(TEMP2);
	Z2 = MULT_II(D, 4);
	s = COMPAREI(Z2, Z1);
	if (s == 1) {
		F = ONEI();
	}
	else {
		F = ZEROI();
	}
	FREEMPI(Z1);
	FREEMPI(Z2);

	P = ZEROI();
	Q = ONEI();

	B2 = ONEI();
	C2 = ZEROI();
	a = 1;
	flag = 0;
	while (1) {
		OLDC2 = C2;
		OLDB2 = B2;
		olda = a;
		TEMP = ADD0I(P, X);
		C = INT0(TEMP, Q);
		FREEMPI(TEMP);
		TEMP = MULTI(C, Q);
		P1 = SUB0I(TEMP, P);
		FREEMPI(TEMP);
		TEMP = MULTI(P1, P1);
		Q1 = SUB0I(D, TEMP);
		FREEMPI(TEMP);
		P2 = ADD0I(P1, Q);
		TEMP = MULTI(P2, P2);
		Q2 = SUB0I(TEMP, D);
		FREEMPI(TEMP);
		OLDP = P;
		OLDQ = Q;
		T = NINT0(P, Q, X, F);
		if (EQUALI(T, C)) {
			a = 1;
			P = P1;
			FREEMPI(P2);
			Q = INT0(Q1, Q);
			B = C;
		}
		else {
			a = -1;
			P = P2;
			FREEMPI(P1);
			Q = INT0(Q2, Q);
			B = ADD0_I(C, 1);
			FREEMPI(C);
		}
		FREEMPI(Q1);
		FREEMPI(Q2);
		FREEMPI(T);

		if (flag == 0) {
			C2 = ONEI();
			flag = 1;
		}
		else {
			C2 = MULTAB_PLUS_MINUS_CI(B, C2, olda, B2);
		}
		FREEMPI(B);
		B2 = OLDC2;
		if (EQUALI(P, OLDP)) { /* period-length = k = 2*h, RCF period p is even */
			if (olda == 1) {
				*Yptr = MULTAB_PLUS_CDI(OLDC2, C2, OLDC2, OLDB2);
			}
			else {
				*Yptr = MULTAB_MINUS_CDI(OLDC2, C2, OLDC2, OLDB2);
			}
			*Xptr = SQRT4TIMING(D, X, *Yptr, -1);
			break;
		}

		TEMP = ADD0I(OLDP, OLDQ);
		t = EQUALI(P, TEMP);
		FREEMPI(TEMP);
		if (t) { /* period-length = k = 2*h, RCF period p is even */
			TEMP = SUBI(OLDB2, OLDC2);
			TEMP1 = ADDI(C2, TEMP);
			FREEMPI(TEMP);
			*Yptr = MULTI(OLDC2, TEMP1);
			FREEMPI(TEMP1);
			*Xptr = SQRT4TIMING(D, X, *Yptr, -1);
			break;
		}

		if (EQUALI(Q, OLDQ)) { /* period-length = 2h+1, p even if a = -1, p odd if a = 1 */
			flag4 = 1;
			if (a == 1) {
				*Yptr = MULTAB_PLUS_CDI(C2, C2, OLDC2, OLDC2);
			}
			else {
				*Yptr = MULTAB_MINUS_CDI(C2, C2, OLDC2, OLDC2);
			}
			*Xptr = SQRT4TIMING(D, X, *Yptr, -a);
			break;
		}
		t1 = MOD0_(Q, 2);
		TEMP3 = INT0_(Q, 2);
		TEMP4 = ADD0I(OLDQ, TEMP3);
		FREEMPI(TEMP3);
		t = EQUALI(P, TEMP4);
		FREEMPI(TEMP4);
		if (t1 == 0 && t && a == -1) { /* period-length = 2h+1, p odd */
			flag5 = 1;
			TEMP = MULTI(C2, C2);
			TEMP1 = MULT_I(OLDC2, 2);
			TEMP2 = SUBI(OLDC2, C2);
			TEMP3 = MULTI(TEMP1, TEMP2);
			*Yptr = ADDI(TEMP, TEMP3);
			FREEMPI(TEMP);
			FREEMPI(TEMP1);
			FREEMPI(TEMP2);
			FREEMPI(TEMP3);
			*Xptr = SQRT4TIMING(D, X, *Yptr, 1);
			break;
		}
		if (flag) {
			TEMP3 = MULTI(OLDP, OLDP);
			TEMP4 = SUBI(D, TEMP3);
			FREEMPI(TEMP3);
			ROLDQ = INT0(TEMP4, OLDQ);
			FREEMPI(TEMP4);
			if (olda == -1) {
				ROLDQ->S = -(ROLDQ->S);
			}
			t1 = MOD0_(ROLDQ, 2);
			TEMP3 = INT0_(ROLDQ, 2);
			TEMP4 = ADD0I(OLDQ, TEMP3);
			FREEMPI(TEMP3);
			FREEMPI(ROLDQ);
			t = EQUALI(OLDP, TEMP4);
			FREEMPI(TEMP4);
			if (t1 == 0 && t && olda == -1) { /* period-length = 2h, p odd */
				TEMP = MULTI(OLDB2, OLDB2);
				TEMP1 = SUBI(OLDC2, OLDB2);
				TEMP2 = MULT_II(OLDC2, 2);
				TEMP3 = MULTI(TEMP2, TEMP1);
				*Yptr = ADD0I(TEMP, TEMP3);
				FREEMPI(TEMP);
				FREEMPI(TEMP1);
				FREEMPI(TEMP2);
				FREEMPI(TEMP3);
				*Xptr = SQRT4TIMING(D, X, *Yptr, 1);
				break;
			}
		}
		FREEMPI(OLDB2);
		FREEMPI(OLDP);
		FREEMPI(OLDQ);
	}
	FREEMPI(OLDP);
	FREEMPI(OLDQ);
	FREEMPI(P);
	FREEMPI(Q);
	FREEMPI(F);
	FREEMPI(X);
	FREEMPI(B2);
	FREEMPI(C2);
	FREEMPI(OLDB2);
	return;
}

void TIME_COUNT_RCF_D(MPI *M, MPI *N) {
	MPI *D, *X, *Y, *TEMP;
	clock_t time_1, time_2;

	time_1 = clock();
	D = COPYI(M);
	while (RSV(D, N) <= 0) {
		RCF_PERIOD00(D, &X, &Y);
		FREEMPI(X);
		FREEMPI(Y);
		TEMP = D;
		D = ADD0_I(D, 1);
		FREEMPI(TEMP);
	}
	FREEMPI(D);
	time_2 = clock();
	printf("RCF:CPU time taken=%f seconds\n", (float)(time_2 - time_1) / (float)CLOCKS_PER_SEC);
	return;
}

void TIME_COUNT_NICF_D(MPI *M, MPI *N) {
	MPI *D, *X, *Y, *TEMP;
	clock_t time_1, time_2;

	printf("doing NICF\n");
	PRINT_MPI(M, "M");
	PRINT_MPI(N, "N");
	time_1 = clock();
	D = COPYI(M);
	while (RSV(D, N) <= 0) {
		NICF_PERIOD00(D, &X, &Y);
		FREEMPI(X);
		FREEMPI(Y);
		TEMP = D;
		D = ADD0_I(D, 1);
		FREEMPI(TEMP);
	}
	FREEMPI(D);
	time_2 = clock();
	printf("NICF:CPU time taken=%f seconds\n", (float)(time_2 - time_1) / (float)CLOCKS_PER_SEC);
	return;
}

void TIME_COUNT_NSCF_D(MPI *M, MPI *N) {
	MPI *D, *X, *Y, *TEMP;
	clock_t time_1, time_2;

	printf("doing NSCF\n");
	PRINT_MPI(M, "M");
	PRINT_MPI(N, "N");
	time_1 = clock();
	D = COPYI(M);
	while (RSV(D, N) <= 0) {
		PRINT_MPI(D, "D");
		NSCF_PERIOD00(D, &X, &Y);
		FREEMPI(X);
		FREEMPI(Y);
		TEMP = D;
		D = ADD0_I(D, 1);
		FREEMPI(TEMP);
	}
	FREEMPI(D);
	time_2 = clock();
	printf("NSCF:CPU time taken=%f seconds\n", (float)(time_2 - time_1) / (float)CLOCKS_PER_SEC);
	return;
}

void PERIOD_TIME_COUNT_RCF1(MPI *M, MPI *N) {
	MPI *D, *PERIOD, *TEMP;
	int t;
	clock_t time_1, time_2;

	time_1 = clock();
	D = COPYI(M);
	while (RSV(D, N) <= 0) {
		t = SQUARETEST(D);
		if (t) {
			PERIOD = CFRAC_PERIOD1(D);
			PRINT_MPI(PERIOD, "PERIOD");
			if (EQZEROI(PERIOD)) {
				FREEMPI(D);
				return;
			}
			else {
				FREEMPI(PERIOD);
			}
		}
		TEMP = D;
		D = ADD0_I(D, 1);
		FREEMPI(TEMP);
	}
	FREEMPI(D);
	time_2 = clock();
	printf("RCF:CPU time taken=%f seconds\n", (float)(time_2 - time_1) / (float)CLOCKS_PER_SEC);
	return;
}

void PERIOD_TIME_COUNT_NSCF1(MPI *M, MPI *N) {
	MPI *D, *PERIOD, *TEMP;
	int t;
	clock_t time_1, time_2;

	time_1 = clock();
	D = COPYI(M);
	while (RSV(D, N) <= 0) {
		t = SQUARETEST(D);
		if (t) {
			PERIOD = NSCF_PERIOD1(D);
			PRINT_MPI(D, "D");
			PRINT_MPI(PERIOD, "PERIOD");
			if (EQZEROI(PERIOD)) {
				FREEMPI(D);
				return;
			}
			else {
				FREEMPI(PERIOD);
			}
		}
		TEMP = D;
		D = ADD0_I(D, 1);
		FREEMPI(TEMP);
	}
	FREEMPI(D);
	time_2 = clock();
	printf("NSCF:CPU time taken=%f seconds\n", (float)(time_2 - time_1) / (float)CLOCKS_PER_SEC);
	return;
}

void PERIOD_TIME_COUNT_NICF1(MPI *M, MPI *N) {
	MPI *D, *PERIOD, *TEMP;
	int t;
	clock_t time_1, time_2;

	time_1 = clock();
	D = COPYI(M);
	while (RSV(D, N) <= 0) {
		t = SQUARETEST(D);
		if (t) {
			PERIOD = NICF_PERIOD1(D);
			PRINT_MPI(D, "D");
			PRINT_MPI(PERIOD, "PERIOD");
			if (EQZEROI(PERIOD)) {
				FREEMPI(D);
				return;
			}
			else {
				FREEMPI(PERIOD);
			}
		}
		TEMP = D;
		D = ADD0_I(D, 1);
		FREEMPI(TEMP);
	}
	FREEMPI(D);
	time_2 = clock();
	printf("NICF:CPU time taken=%f seconds\n", (float)(time_2 - time_1) / (float)CLOCKS_PER_SEC);
	return;
}

/* This function finds the period of the NICF-P of sqrt(D) and exits
* when it reaches one of the 5 conditions of the nicf.pdf paper.
* Single precision is used in the calculation of the period.
*/
MPI *NICF_PERIOD1(MPI *D) {
	MPI *X, *G, *Z1, *Z2, *F, *TEMP1, *TEMP2, *C, *T;
	MPI *P1, *P2, *Q1, *Q2, *TEMP, *P, *Q, *OLDP, *OLDQ, *TEMP3, *TEMP4;
	MPI *ROLDQ, *PERIOD;
	USL h, t1, period;
	USI t, olda, flag4 = 0, flag5 = 0;
	int a, s;

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
	TEMP1 = MULT_II(X, 2);
	TEMP2 = ADD0_I(TEMP1, 1);
	Z1 = MULTI(TEMP2, TEMP2);
	FREEMPI(TEMP1);
	FREEMPI(TEMP2);
	Z2 = MULT_II(D, 4);
	s = COMPAREI(Z2, Z1);
	if (s == 1) {
		F = ONEI();
	}
	else {
		F = ZEROI();
	}
	FREEMPI(Z1);
	FREEMPI(Z2);

	P = ZEROI();
	Q = ONEI();
	h = 0;

	a = 1;
	while (1) {
		olda = a;
		TEMP = ADDI(P, X);
		C = INT0(TEMP, Q);
		FREEMPI(TEMP);
		TEMP = MULTI(C, Q);
		P1 = SUB0I(TEMP, P);
		FREEMPI(TEMP);
		TEMP = MULTI(P1, P1);
		Q1 = SUB0I(D, TEMP);
		FREEMPI(TEMP);
		P2 = ADD0I(P1, Q);
		TEMP = MULTI(P2, P2);
		Q2 = SUB0I(TEMP, D);
		FREEMPI(TEMP);
		OLDP = P;
		OLDQ = Q;
		T = NINT0(P, Q, X, F);
		if (EQUALI(T, C)) {
			a = 1;
			P = P1;
			FREEMPI(P2);
			Q = INT0(Q1, Q);
		}
		else {
			a = -1;
			P = P2;
			Q = INT0(Q2, Q);
			FREEMPI(P1);
		}
		FREEMPI(Q1);
		FREEMPI(Q2);
		FREEMPI(C);
		FREEMPI(T);

		if (EQUALI(P, OLDP)) { /* period-length = k = 2*h, RCF period p is even */
			printf("P[h] = P[h+1]\n");
			break;
		}

		TEMP = ADD0I(OLDP, OLDQ);
		t = EQUALI(P, TEMP);
		FREEMPI(TEMP);
		if (t) { /* period-length = k = 2*h, RCF period p is even */
			printf("P[h+1] = P[h] + Q[h]\n");
			break;
		}

		if (EQUALI(Q, OLDQ)) { /* period-length = 2h+1, p even if a = -1, p odd if a = 1 */
			printf("Q[h] = Q[h+1]\n");
			flag4 = 1;
			break;
		}
		t1 = MOD0_(Q, 2);
		TEMP3 = INT0_(Q, 2);
		TEMP4 = ADD0I(OLDQ, TEMP3);
		FREEMPI(TEMP3);
		t = EQUALI(P, TEMP4);
		FREEMPI(TEMP4);
		if (t1 == 0 && t && a == -1) { /* period-length = 2h+1, p odd */
			printf("P[h+1] = Q[h] + Q[h+1]/2 and epsilon[h+1] = -1\n");
			flag5 = 1;
			break;
		}
		if (h != 0) {
			TEMP3 = MULTI(OLDP, OLDP);
			TEMP4 = SUBI(D, TEMP3);
			FREEMPI(TEMP3);
			ROLDQ = INT(TEMP4, OLDQ);
			FREEMPI(TEMP4);
			if (olda == -1) {
				ROLDQ->S = -ROLDQ->S;
			}
			t1 = MOD0_(ROLDQ, 2);
			TEMP3 = INT0_(ROLDQ, 2);
			TEMP4 = ADD0I(OLDQ, TEMP3);
			FREEMPI(TEMP3);
			FREEMPI(ROLDQ);
			t = EQUALI(OLDP, TEMP4);
			FREEMPI(TEMP4);
			if (t1 == 0 && t && olda == -1) { /* period-length = 2h, p odd */
				printf("P[h] = Q[h] + Q[h-1]/2 and epsilon[h] = -1\n");
				break;
			}
		}
		h++;
		FREEMPI(OLDP);
		FREEMPI(OLDQ);
	}
	FREEMPI(OLDP);
	FREEMPI(OLDQ);
	FREEMPI(P);
	FREEMPI(Q);
	FREEMPI(F);
	FREEMPI(X);
	if (h<2147483648UL) {
		period = 2 * h;
	}
	else {
		printf("Exiting prematurely - period >= 2^32\n");
		period = 0;
		return(ZEROI());
	}
	if (flag4 || flag5) {
		period = period + 1;
	}
	PERIOD = CHANGE(period);
	return(PERIOD);
}

void TIME_COUNT_RCF_QIMPROVED(MPI *M, MPI *N) {
	MPI *D, *X, *Y, *TEMP;
	clock_t time_1, time_2;

	printf("doing RCF\n");
	PRINT_MPI(M, "M");
	PRINT_MPI(N, "N");
	time_1 = clock();
	D = COPYI(M);
	while (RSV(D, N) <= 0) {
		RCF_PERIOD_QIMPROVED(D, &X, &Y);
		FREEMPI(X);
		FREEMPI(Y);
		TEMP = D;
		D = ADD0_I(D, 1);
		FREEMPI(TEMP);
	}
	FREEMPI(D);
	time_2 = clock();
	printf("RCF QIMPROVED:CPU time taken=%f seconds\n", (float)(time_2 - time_1) / (float)CLOCKS_PER_SEC);
	return;
}

/* This function finds smallest solution of Pell's equation using NICF-P of sqrt(D) and exits
when it reaches one of the 5 conditions of the nicf.pdf paper.
We use Lemma 2 of l_nicf=l_nscf.tex. */
void NICF_PERIOD_QIMPROVED(MPI *D, MPI **Xptr, MPI **Yptr)
{
	MPI *X, *G, *TEMP1, *TEMP2, *C, *Z, *T, *TT, *TTT;
	MPI *P1, *P2, *Q1, *Q2, *TEMP, *P, *Q, *OLDP, *OLDQ, *TEMP3, *TEMP4;
	MPI *ROLDQ;
	MPI *B1, *C1, *B2, *C2, *OLDB1, *OLDC1, *OLDB2, *OLDC2, *B;
	USL t1;
	USI flag, t, olda, flag4 = 0, flag5 = 0;
	int a;

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
	t = EQUALI(D, TEMP);
	FREEMPI(TEMP);
	if (t) {/* period 1, exceptional case */
		*Xptr = X;
		*Yptr = ONEI();
		return;
	}
	T = MULT_II(X, 2);
	TT = ADD0_I(T, 1);
	TTT = MULT_II(D, 4);
	TEMP1 = ADD0I(G, X);  /* TEMP1 = X^2 + X */
	TEMP2 = MULT_II(TEMP1, 4); /* TEMP2 = 4(X^2 + X) */
	TEMP = ADD0_I(TEMP2, 1); /* TEMP = 4(X^2 + X) +1 = (2X+1)^2 = TT*2*/
	t = RSV(TEMP, TTT);
	if (t < 0) { /* z=int(sqrt(4d)) = 2int(sqrt(d)) if z is even, but = 2int(sqrt(d))+1 if z is odd.*/
		Z = T;
		FREEMPI(TT);
	}
	else {
		Z = TT;
		FREEMPI(T);
	}
	FREEMPI(G);
	FREEMPI(TTT);
	FREEMPI(TEMP);
	FREEMPI(TEMP1);
	FREEMPI(TEMP2);

	P = ZEROI();
	Q = ONEI();

	B1 = ZEROI();
	C1 = ONEI();
	B2 = ONEI();
	C2 = ZEROI();
	a = 1;
	flag = 0;
	while (1) {
		OLDC1 = C1;
		OLDC2 = C2;
		OLDB1 = B1;
		OLDB2 = B2;
		olda = a;
		TEMP = ADD0I(P, X);
		C = INT0(TEMP, Q);
		FREEMPI(TEMP);
		TEMP = MULTI(C, Q);
		P1 = SUB0I(TEMP, P);
		FREEMPI(TEMP);
		TEMP = MULTI(P1, P1);
		Q1 = SUB0I(D, TEMP);
		FREEMPI(TEMP);
		P2 = ADD0I(P1, Q);
		TEMP = MULTI(P2, P2);
		Q2 = SUB0I(TEMP, D);
		FREEMPI(TEMP);
		OLDP = P;
		OLDQ = Q;
		TEMP = ADD0I(Q1, Q2);
		if (RSV(TEMP, Z) >= 0) {
			a = 1;
			P = P1;
			FREEMPI(P2);
			Q = INT0(Q1, Q);
			B = C;
		}
		else {
			a = -1;
			P = P2;
			FREEMPI(P1);
			Q = INT0(Q2, Q);
			B = ADD0_I(C, 1);
			FREEMPI(C);
		}
		FREEMPI(TEMP);
		FREEMPI(Q1);
		FREEMPI(Q2);

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
		if (EQUALI(P, OLDP)) { /* period-length = k = 2*h, RCF period p is even */
			if (olda == 1) {
				*Xptr = MULTAB_PLUS_CDI(OLDC2, C1, OLDC1, OLDB2);
				*Yptr = MULTAB_PLUS_CDI(OLDC2, C2, OLDC2, OLDB2);
			}
			else {
				*Xptr = MULTAB_MINUS_CDI(OLDC2, C1, OLDC1, OLDB2);
				*Yptr = MULTAB_MINUS_CDI(OLDC2, C2, OLDC2, OLDB2);
			}
			break;
		}

		TEMP = ADD0I(OLDP, OLDQ);
		t = EQUALI(P, TEMP);
		FREEMPI(TEMP);
		if (t) { /* period-length = k = 2*h, RCF period p is even */
			TEMP = SUBI(OLDB2, OLDC2);
			*Xptr = MULTAB_PLUS_CDI(OLDC2, C1, OLDC1, TEMP);
			FREEMPI(TEMP);
			TEMP = SUBI(OLDB2, OLDC2);
			TEMP1 = ADDI(C2, TEMP);
			FREEMPI(TEMP);
			*Yptr = MULTI(OLDC2, TEMP1);
			FREEMPI(TEMP1);
			break;
		}

		if (EQUALI(Q, OLDQ)) { /* period-length = 2h+1, p even if a = -1, p odd if a = 1 */
			flag4 = 1;
			if (a == 1) {
				*Xptr = MULTAB_PLUS_CDI(C1, C2, OLDC1, OLDC2);
				*Yptr = MULTAB_PLUS_CDI(C2, C2, OLDC2, OLDC2);
			}
			else {
				*Xptr = MULTAB_MINUS_CDI(C1, C2, OLDC1, OLDC2);
				*Yptr = MULTAB_MINUS_CDI(C2, C2, OLDC2, OLDC2);
			}
			break;
		}
		t1 = MOD0_(Q, 2);
		TEMP3 = INT0_(Q, 2);
		TEMP4 = ADD0I(OLDQ, TEMP3);
		FREEMPI(TEMP3);
		t = EQUALI(P, TEMP4);
		FREEMPI(TEMP4);
		if (t1 == 0 && t && a == -1) { /* period-length = 2h+1, p odd */
			flag5 = 1;
			TEMP = MULTI(C1, C2);
			TEMP1 = MULT_I(OLDC1, 2);
			TEMP2 = MULTAB_PLUS_CDI(C2, OLDC1, C1, OLDC2);
			TEMP3 = MULTI(TEMP1, OLDC2);
			FREEMPI(TEMP1);
			TEMP1 = ADDI(TEMP, TEMP3);
			FREEMPI(TEMP);
			*Xptr = SUBI(TEMP1, TEMP2);
			FREEMPI(TEMP1);
			FREEMPI(TEMP2);
			FREEMPI(TEMP3);

			TEMP = MULTI(C2, C2);
			TEMP1 = MULT_I(OLDC2, 2);
			TEMP2 = SUBI(OLDC2, C2);
			TEMP3 = MULTI(TEMP1, TEMP2);
			*Yptr = ADDI(TEMP, TEMP3);
			FREEMPI(TEMP);
			FREEMPI(TEMP1);
			FREEMPI(TEMP2);
			FREEMPI(TEMP3);
			break;
		}
		if (flag) {
			TEMP3 = MULTI(OLDP, OLDP);
			TEMP4 = SUBI(D, TEMP3);
			FREEMPI(TEMP3);
			ROLDQ = INT(TEMP4, OLDQ);
			FREEMPI(TEMP4);
			if (olda == -1) {
				ROLDQ->S = -ROLDQ->S;
			}
			t1 = MOD0_(ROLDQ, 2);
			TEMP3 = INT0_(ROLDQ, 2);
			TEMP4 = ADD0I(OLDQ, TEMP3);
			FREEMPI(TEMP3);
			FREEMPI(ROLDQ);
			t = EQUALI(OLDP, TEMP4);
			FREEMPI(TEMP4);
			if (t1 == 0 && t && olda == -1) { /* period-length = 2h, p odd */
				TEMP = MULT_II(OLDC1, 2);
				TEMP2 = MULTAB_PLUS_CDI(TEMP, OLDC2, OLDB1, OLDB2);
				TEMP1 = MULTAB_PLUS_CDI(OLDC2, OLDB1, OLDC1, OLDB2);
				*Xptr = SUBI(TEMP2, TEMP1);
				FREEMPI(TEMP);
				FREEMPI(TEMP1);
				FREEMPI(TEMP2);

				TEMP = MULTI(OLDB2, OLDB2);
				TEMP1 = SUBI(OLDC2, OLDB2);
				TEMP2 = MULT_II(OLDC2, 2);
				TEMP3 = MULTI(TEMP2, TEMP1);
				*Yptr = ADDI(TEMP, TEMP3);
				FREEMPI(TEMP);
				FREEMPI(TEMP1);
				FREEMPI(TEMP2);
				FREEMPI(TEMP3);
				break;
			}
		}
		FREEMPI(OLDB1);
		FREEMPI(OLDB2);
		FREEMPI(OLDP);
		FREEMPI(OLDQ);
	}
	FREEMPI(Z);
	FREEMPI(OLDP);
	FREEMPI(OLDQ);
	FREEMPI(P);
	FREEMPI(Q);
	FREEMPI(X);
	FREEMPI(B1);
	FREEMPI(B2);
	FREEMPI(C1);
	FREEMPI(C2);
	FREEMPI(OLDB1);
	FREEMPI(OLDB2);
	return;
}

void TIME_COUNT_NICF_QIMPROVED(MPI *M, MPI *N) {
	MPI *D, *X, *Y, *TEMP;
	clock_t time_1, time_2;

	time_1 = clock();
	D = COPYI(M);
	while (RSV(D, N) <= 0) {
		NICF_PERIOD_QIMPROVED(D, &X, &Y);
		FREEMPI(X);
		FREEMPI(Y);
		TEMP = D;
		D = ADD0_I(D, 1);
		FREEMPI(TEMP);
	}
	FREEMPI(D);
	time_2 = clock();
	printf("NICFQIMPROVED:CPU time taken=%f seconds\n", (float)(time_2 - time_1) / (float)CLOCKS_PER_SEC);
	return;
}
