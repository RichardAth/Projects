#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_NONSTDC_NO_DEPRECATE
#include "integer.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fun.h"
#include "stack.h"
/*#include "unistd.h"*/


MPI *GLOBALA;
MPI *GLOBALB;
MPI *GLOBALC;
MPI *GLOBALD;
USL count, flagl, flagr, globalr, globall;
char buff[20];
FILE *outfile;

USL POSITIVITY(MPI *A, MPI *B, MPI *C, MPI *D) {

	if (A->S >= 0 && B->S >= 0 && C->S >= 0 && D->S >= 0) {
		return(1);
	}
	else {
		return(0);
	}
}


/* NPROD(l,m) finds the least n such that A_n and D=A_0...A_n are
* non-negative.
* The matrix of global variables D=[GLOBALA,GLOBALB,GLOBALC,GLOBALD] is
* returned along with n. See Proposition 4.1, J.L. Davison, 'An algorithm for
* the continued fraction of e^{l/m}', Proceedings of the Eighth Manitoba
* Conference on Numerical Mathematics and Computing (Univ. Manitoba, Winnipeg,
* 1978), 169--179, Congress. Numer., XXII, Utilitas Math.
*/
USL NPROD(USL l, USL m) {
	USL k, s, temp;
	MPI *TEMPM, *TEMPL, *T, *TEMP1, *TEMP2, *A1, *B1, *C1, *D1;
	MPI *TEMP, *T1, *T2;

	TEMPM = CHANGE(m);
	TEMPL = CHANGE(l);
	temp = m + l;
	GLOBALA = CHANGE(temp);
	GLOBALB = CHANGE(m);
	GLOBALC = CHANGE(m);
	GLOBALD = SUBI(TEMPM, TEMPL);

	for (k = 1; 1; k++) {
		s = POSITIVITY(GLOBALA, GLOBALB, GLOBALC, GLOBALD);
		if (s)
			break;
		T = CHANGE((2 * k + 1) * m);
		T1 = ADDI(GLOBALA, GLOBALB);
		T2 = ADDI(GLOBALC, GLOBALD);
		TEMP1 = MULTI(T, T1);
		TEMP2 = MULTI(T, T2);
		FREEMPI(T);
		FREEMPI(T1);
		FREEMPI(T2);
		TEMP = MULTI(GLOBALA, TEMPL);
		A1 = ADDI(TEMP1, TEMP);
		FREEMPI(TEMP);
		TEMP = MULTI(GLOBALB, TEMPL);
		B1 = SUBI(TEMP1, TEMP);
		FREEMPI(TEMP);
		TEMP = MULTI(GLOBALC, TEMPL);
		C1 = ADDI(TEMP2, TEMP);
		FREEMPI(TEMP);
		TEMP = MULTI(GLOBALD, TEMPL);
		D1 = SUBI(TEMP2, TEMP);
		FREEMPI(TEMP);
		FREEMPI(TEMP1);
		FREEMPI(TEMP2);
		TEMP = GLOBALA;
		GLOBALA = A1;
		FREEMPI(TEMP);
		TEMP = GLOBALB;
		GLOBALB = B1;
		FREEMPI(TEMP);
		TEMP = GLOBALC;
		GLOBALC = C1;
		FREEMPI(TEMP);
		TEMP = GLOBALD;
		GLOBALD = D1;
		FREEMPI(TEMP);
	}
	FREEMPI(TEMPL);
	FREEMPI(TEMPM);
	return(k - 1);
}

MPI *REFINE(MPI *A, MPI *B, MPI *C, MPI *D) {
	MPI *G, *TEMP;

	G = GCD(A, B);
	TEMP = G;
	G = GCD(G, C);
	FREEMPI(TEMP);
	TEMP = G;
	G = GCD(G, D);
	FREEMPI(TEMP);
	return(G);
}


/* Raney factorization - G.N. Raney, Math, Annalen 206 (1973) 265-283.
* Input: a non-singular matrix A=[p,q;r,s], p,q,r,s>=0, A!=I_2, A!=[0,1;1,0].
* With L=[1,0;1,1] and R=[1,1;0,1], we express A uniquely as
* a product of non-negative powers of L and R, (at least one is positive)
* followed by a row-balanced B.
* B=[a,b;c,d] is row-balanced if (a<c & b>d) or (c<a & d>b) and a,b,c>=0.
* The number k of powers of L and R is returned.
*/

USL RANEY(MPI *P, MPI *Q, MPI *R, MPI *S) {
	USL k, i, j;
	MPI *PP, *QQ, *RR, *SS, *TEMP;

	PP = COPYI(P);
	QQ = COPYI(Q);
	RR = COPYI(R);
	SS = COPYI(S);
	k = 0;
	while (1) {
		i = 0;
		while (COMPAREI(PP, RR) >= 0 && COMPAREI(QQ, SS) >= 0) {
			TEMP = PP;
			PP = SUBI(PP, RR);
			FREEMPI(TEMP);
			TEMP = QQ;
			QQ = SUBI(QQ, SS);
			FREEMPI(TEMP);
			i++;
		}
		if (i > 0) {
			globalr = globalr + i;
			flagr = 1;
			if (flagr * flagl) {
				flagl = 0;
				printf("a[%lu]:%lu\n", count, globall);
				fprintf(outfile, "a[%lu]:%lu\n", count, globall);
				globall = 0;
				count++;
			}
			k++;
		}
		j = 0;
		while (COMPAREI(RR, PP) >= 0 && COMPAREI(SS, QQ) >= 0) {
			TEMP = RR;
			RR = SUBI(RR, PP);
			FREEMPI(TEMP);
			TEMP = SS;
			SS = SUBI(SS, QQ);
			FREEMPI(TEMP);
			j++;
		}
		if (j > 0) {
			globall = globall + j;
			flagl = 1;
			if (flagr * flagl) {
				flagr = 0;
				printf("a[%lu]:%lu\n", count, globalr);
				fprintf(outfile, "a[%lu]:%lu\n", count, globalr);
				globalr = 0;
				count++;
			}
			k++;
		}
		if ((COMPAREI(PP, RR) < 0 && COMPAREI(QQ, SS) > 0) || (COMPAREI(PP, RR) > 0 && COMPAREI(QQ, SS) < 0)) {
			break;
		}
	}
	TEMP = GLOBALA;
	GLOBALA = PP;
	FREEMPI(TEMP);
	TEMP = GLOBALB;
	GLOBALB = QQ;
	FREEMPI(TEMP);
	TEMP = GLOBALC;
	GLOBALC = RR;
	FREEMPI(TEMP);
	TEMP = GLOBALD;
	GLOBALD = SS;
	FREEMPI(TEMP);
	/*fclose(outfile);*/
	return(k);
}

/* We perform the algorithm of J.L. Davison's paper.
* With n > =0, we first find the n* of Davison's Proposition 4.1
* and apply Raney's factorisation to A_0...A_k, for n*<=k<=n*+n.
* The number (count) of partial quotients of e^{l/m} found is returned.
* count becomes positive for all large n.
* We exit the program if count reaches 100000.
*/
USL DAVISON(USL l, USL m, USL n) {
	MPI *G, *ONE, *TEMP, *T, *T1, *T2, *TEMP1, *TEMP2;
	MPI *A1, *B1, *C1, *D1, *TEMPL;
	USL i, j, k, t;

	printf("l,m,n=%lu,%lu,%lu\n", l, m, n);
	strcpy(buff, "davison.out");
	outfile = fopen(buff, "w");
	count = 0;
	flagr = 0;
	flagl = 0;
	globalr = 0;
	globall = 0;
	k = NPROD(l, m);
	G = REFINE(GLOBALA, GLOBALB, GLOBALC, GLOBALD);
	ONE = ONEI();
	if (RSV(G, ONE)>0) {
		TEMP = GLOBALA;
		GLOBALA = INT(GLOBALA, G);
		FREEMPI(TEMP);
		TEMP = GLOBALB;
		GLOBALB = INT(GLOBALB, G);
		FREEMPI(TEMP);
		TEMP = GLOBALC;
		GLOBALC = INT(GLOBALC, G);
		FREEMPI(TEMP);
		TEMP = GLOBALD;
		GLOBALD = INT(GLOBALD, G);
		FREEMPI(TEMP);
	}
	FREEMPI(G);
	i = k;
	j = k + n;
	TEMPL = CHANGE(l);
	while (i <= j) {
		if (i > k) {
			T = CHANGE((2 * i + 1) * m);
			T1 = ADDI(GLOBALA, GLOBALB);
			T2 = ADDI(GLOBALC, GLOBALD);
			TEMP1 = MULTI(T, T1);
			TEMP2 = MULTI(T, T2);
			FREEMPI(T);
			FREEMPI(T1);
			FREEMPI(T2);
			TEMP = MULTI(GLOBALA, TEMPL);
			A1 = ADDI(TEMP1, TEMP);
			FREEMPI(TEMP);
			TEMP = MULTI(GLOBALB, TEMPL);
			B1 = SUBI(TEMP1, TEMP);
			FREEMPI(TEMP);
			TEMP = MULTI(GLOBALC, TEMPL);
			C1 = ADDI(TEMP2, TEMP);
			FREEMPI(TEMP);
			TEMP = MULTI(GLOBALD, TEMPL);
			D1 = SUBI(TEMP2, TEMP);
			FREEMPI(TEMP);
			FREEMPI(TEMP1);
			FREEMPI(TEMP2);
			TEMP = GLOBALA;
			GLOBALA = A1;
			FREEMPI(TEMP);
			TEMP = GLOBALB;
			GLOBALB = B1;
			FREEMPI(TEMP);
			TEMP = GLOBALC;
			GLOBALC = C1;
			FREEMPI(TEMP);
			TEMP = GLOBALD;
			GLOBALD = D1;
			FREEMPI(TEMP);
		}
		i++;
		t = RANEY(GLOBALA, GLOBALB, GLOBALC, GLOBALD);
		if (count == 1000000) {
			printf("1000000 partial quotients found\n");
			i = j + 1;
		}
		/* the input matrix will not be I_2 or [0,1;1,0] here. */
		G = REFINE(GLOBALA, GLOBALB, GLOBALC, GLOBALD);
		if (RSV(G, ONE)>0) {
			TEMP = GLOBALA;
			GLOBALA = INT(GLOBALA, G);
			FREEMPI(TEMP);
			TEMP = GLOBALB;
			GLOBALB = INT(GLOBALB, G);
			FREEMPI(TEMP);
			TEMP = GLOBALC;
			GLOBALC = INT(GLOBALC, G);
			FREEMPI(TEMP);
			TEMP = GLOBALD;
			GLOBALD = INT(GLOBALD, G);
			FREEMPI(TEMP);
		}
		FREEMPI(G);
	}
	FREEMPI(TEMPL);
	FREEMPI(ONE);
	FREEMPI(GLOBALA);
	FREEMPI(GLOBALB);
	FREEMPI(GLOBALC);
	FREEMPI(GLOBALD);
	fflush(stdout);
	printf("n* = %lu, n = %lu\n", k, n);
	fprintf(outfile, "n* = %lu, n = %lu\n", k, n);
	printf("The number of partial quotients found for e^(%lu/%lu) is %lu\n", l, m, count);
	fprintf(outfile, "The number of partial quotients found for e^(%lu/%lu) is %lu\n", l, m, count);
	fclose(outfile);
	return(count);
}

MPI *DAVISONX(MPI *L, MPI *M, MPI *N) {
	USL l, m, n, t;
	MPI *G, *TEMP, *TEMPL, *TEMPM, *ONE;

	if (L->S <= 0) {
		printf("l<=0\n");
		return(NULL);
	}
	if (M->S <= 0) {
		printf("m<=0\n");
		return(NULL);
	}
	if (N->S < 0) {
		printf("n<0\n");
		return(NULL);
	}
	TEMP = CHANGEI(1000);
	if (RSV(L, TEMP) > 0) {
		printf("l>1000\n");
		FREEMPI(TEMP);
		return(NULL);
	}
	if (RSV(M, TEMP) > 0) {
		printf("m>1000\n");
		FREEMPI(TEMP);
		return(NULL);
	}
	FREEMPI(TEMP);
	TEMP = CHANGEI(100000);
	if (RSV(N, TEMP) > 0) {
		printf("n>=10^6\n");
		FREEMPI(TEMP);
		return(NULL);
	}
	FREEMPI(TEMP);
	TEMPL = COPYI(L);
	TEMPM = COPYI(M);
	ONE = ONEI();
	G = GCD(L, M);
	if (RSV(G, ONE)>0) {
		TEMP = TEMPL;
		TEMPL = INT0(L, G);
		FREEMPI(TEMPL);
		TEMP = TEMPM;
		TEMPM = INT0(M, G);
		FREEMPI(TEMP);
	}
	FREEMPI(G);
	FREEMPI(ONE);
	l = CONVERTI(TEMPL);
	m = CONVERTI(TEMPM);
	FREEMPI(TEMPL);
	FREEMPI(TEMPM);
	n = CONVERTI(N);
	t = DAVISON(l, m, n);
	return (CHANGE(t));
}

/* Raney factorization - G.N. Raney, Math, Annalen 206 (1973) 265-283.
* Input: a non-singular matrix A=[p,q;r,s], p,q,r,s>=0, A!=I_2, A!=[0,1;1,0].
* With L=[1,0;1,1] and R=[1,1;0,1], we express A uniquely as
* a product of non-negative powers of L and R, (at least one is positive)
* followed by a row-balanced B.
* B=[a,b;c,d] is row-balanced if (a<c & b>d) or (c<a & d>b) and a,b,c>=0.
* The number k of powers of L and R is returned.
*/

USL RANEY1(MPI *P, MPI *Q, MPI *R, MPI *S) {
	USL k, i, j;
	MPI *PP, *QQ, *RR, *SS, *TEMP;
	char buff1[20];
	FILE *outfile1;

	strcpy(buff1, "raney.out");
	outfile1 = fopen(buff1, "w");
	printf("[");
	fprintf(outfile1, "[");
	PRINTI(P);
	FPRINTI(outfile1, P);
	printf(",");
	fprintf(outfile1, ",");
	PRINTI(Q);
	FPRINTI(outfile1, Q);
	printf(";");
	fprintf(outfile1, ";");
	PRINTI(R);
	FPRINTI(outfile1, R);
	printf(",");
	fprintf(outfile1, ",");
	PRINTI(S);
	FPRINTI(outfile1, S);
	printf("]=");
	fprintf(outfile1, "]=");
	PP = COPYI(P);
	QQ = COPYI(Q);
	RR = COPYI(R);
	SS = COPYI(S);
	k = 0;
	while (1) {
		i = 0;
		while (COMPAREI(PP, RR) >= 0 && COMPAREI(QQ, SS) >= 0) {
			TEMP = PP;
			PP = SUBI(PP, RR);
			FREEMPI(TEMP);
			TEMP = QQ;
			QQ = SUBI(QQ, SS);
			FREEMPI(TEMP);
			i++;
		}
		if (i > 0) {
			k++;
			if (i > 1) {
				printf("R^%lu", i);
				fprintf(outfile1, "R^%lu", i);
			}
			else {
				printf("R");
				fprintf(outfile1, "R");
			}
		}
		j = 0;
		while (COMPAREI(RR, PP) >= 0 && COMPAREI(SS, QQ) >= 0) {
			TEMP = RR;
			RR = SUBI(RR, PP);
			FREEMPI(TEMP);
			TEMP = SS;
			SS = SUBI(SS, QQ);
			FREEMPI(TEMP);
			j++;
		}
		if (j > 0) {
			if (j > 1) {
				printf("L^%lu", j);
				fprintf(outfile1, "L^%lu", j);
			}
			else {
				printf("L");
				fprintf(outfile1, "L");
			}
			k++;

		}
		if ((COMPAREI(PP, RR) < 0 && COMPAREI(QQ, SS) > 0) || (COMPAREI(PP, RR) > 0 && COMPAREI(QQ, SS) < 0)) {
			break;
		}
	}
	printf("D, where D=[");
	fprintf(outfile1, "D, where D=[");
	PRINTI(PP);
	FPRINTI(outfile1, PP);
	printf(",");
	fprintf(outfile1, ",");
	PRINTI(QQ);
	FPRINTI(outfile1, QQ);
	printf(";");
	fprintf(outfile1, ";");
	PRINTI(RR);
	FPRINTI(outfile1, RR);
	printf(",");
	fprintf(outfile1, ",");
	PRINTI(SS);
	FPRINTI(outfile1, SS);
	printf("] is row-balanced\n");
	fprintf(outfile1, "] is row-balanced\n");
	FREEMPI(PP);
	FREEMPI(QQ);
	FREEMPI(RR);
	FREEMPI(SS);
	fclose(outfile1);
	return(k);
}

MPI *RANEY1X(MPI *P, MPI *Q, MPI *R, MPI *S) {
	MPI *TEMP1, *TEMP2, *TEMP;
	USI t;
	USL k;

	if (P->S < 0) {
		printf("P < 0\n");
		return(NULL);
	}
	if (Q->S < 0) {
		printf("Q < 0\n");
		return(NULL);
	}
	if (R->S < 0) {
		printf("R < 0\n");
		return(NULL);
	}
	if (S->S < 0) {
		printf("S < 0\n");
		return(NULL);
	}
	if (EQONEI(P) && EQZEROI(Q) && EQZEROI(R) && EQONEI(S)) {
		printf("A = I_2\n");
		return(NULL);
	}
	if (EQZEROI(P) && EQONEI(Q) && EQONEI(R) && EQZEROI(S)) {
		printf("A = J = [0,1;1,0]\n");
		return(NULL);
	}
	TEMP1 = MULTI(P, S);
	TEMP = CHANGEI(1000);
	TEMP2 = MULTI(Q, R);
	t = EQUALI(TEMP1, TEMP2);
	FREEMPI(TEMP1);
	FREEMPI(TEMP2);
	if (t) {
		printf("P*S-Q*R=0\n");
		return(NULL);
	}
	k = RANEY1(P, Q, R, S);
	return(CHANGE(k));
}

/*
* This algorithm expresses a unimodular matrix A !=I_2 or U=[0 1]
*                                                           [1 0]
* with non-negative coefficients, as a product of one of the
* following forms:
* P, UP, PU, or UPU, where P is a product of matrices of the form
* [a 1], a>0.
* [1 0]
* The representation is unique. See Kjell Kolden, 'Continued fractions
* and linear substitutions', Arch. Math. Naturvid. 50 (1949), 141-196.
* The number n of matrices in the product [d[0]1]***[d[n-1]1]
* is returned.                            [  1 0]   [  1   0]
*/

USL UNIMODULAR(MPI *P, MPI *Q, MPI *R, MPI *S) {
	USL i;
	MPI *D, *D1, *D2, *TEMP, *TEMPP, *TEMPQ, *PP, *QQ, *RR, *SS;
	char buff1[20];
	FILE *outfile1;

	i = 0;
	strcpy(buff1, "unimodular.out");
	outfile1 = fopen(buff1, "w");
	printf("[");
	fprintf(outfile1, "[");
	PRINTI(P);
	FPRINTI(outfile1, P);
	printf(",");
	fprintf(outfile1, ",");
	PRINTI(Q);
	FPRINTI(outfile1, Q);
	printf(";");
	fprintf(outfile1, ";");
	PRINTI(R);
	FPRINTI(outfile1, R);
	printf(",");
	fprintf(outfile1, ",");
	PRINTI(S);
	FPRINTI(outfile1, S);
	printf("]=");
	fprintf(outfile1, "]=");
	PP = COPYI(P);
	QQ = COPYI(Q);
	RR = COPYI(R);
	SS = COPYI(S);
	while (1) {
		if (SS->S == 0) {
			printf("U[");
			fprintf(outfile1, "U[");
			PRINTI(PP);
			FPRINTI(outfile1, PP);
			printf("]");
			fprintf(outfile1, "]");
			fflush(stdout);
			i++;
			break;
		}
		if (RR->S == 0) {
			printf("U[");
			fprintf(outfile1, "U[");
			PRINTI(QQ);
			FPRINTI(outfile1, QQ);
			printf("]");
			fprintf(outfile1, "]");
			fflush(stdout);
			i++;
			printf("U[0]");
			fprintf(outfile1, "U[0]");
			fflush(stdout);
			i++;
			break;
		}
		D1 = INT0(QQ, SS);
		D2 = INT0(PP, RR);
		D = MINMPI(D1, D2);
		FREEMPI(D1);
		FREEMPI(D2);
		TEMPP = PP;
		PP = RR;
		TEMP = MULTI(D, RR);
		RR = SUBI(TEMPP, TEMP);
		FREEMPI(TEMPP);
		FREEMPI(TEMP);
		TEMPQ = QQ;
		QQ = SS;
		TEMP = MULTI(D, SS);
		SS = SUBI(TEMPQ, TEMP);
		FREEMPI(TEMPQ);
		FREEMPI(TEMP);
		printf("U[");
		fprintf(outfile1, "U[");
		PRINTI(D);
		FPRINTI(outfile1, D);
		printf("]");
		fprintf(outfile1, "]");
		FREEMPI(D);
		fflush(stdout);
		i++;
	}
	printf("\n");
	fprintf(outfile1, "\n");
	FREEMPI(PP);
	FREEMPI(QQ);
	FREEMPI(RR);
	FREEMPI(SS);
	fclose(outfile1);
	return(i);
}

MPI *UNIMODULARX(MPI *P, MPI *Q, MPI *R, MPI *S) {
	USI t1, t2;
	USL i;
	MPI *DET, *TEMP1, *TEMP2;

	if (P->S < 0) {
		printf("P < 0\n");
		return(NULL);
	}
	if (Q->S < 0) {
		printf("Q < 0\n");
		return(NULL);
	}
	if (R->S < 0) {
		printf("R < 0\n");
		return(NULL);
	}
	if (S->S < 0) {
		printf("S < 0\n");
		return(NULL);
	}
	TEMP1 = MULTI(P, S);
	TEMP2 = MULTI(Q, R);
	DET = SUBI(TEMP1, TEMP2);
	FREEMPI(TEMP1);
	FREEMPI(TEMP2);
	t1 = EQONEI(DET);
	t2 = EQMINUSONEI(DET);
	FREEMPI(DET);
	if (EQONEI(P) && EQZEROI(Q) && EQZEROI(R) && EQONEI(S)) {
		printf("A=I_2\n");
		return(NULL);
	}
	if (EQZEROI(P) && EQONEI(Q) && EQONEI(R) && EQZEROI(S)) {
		printf("A = U = [0,1;1,0]\n");
		return(NULL);
	}
	if (t1 == 0 && t2 == 0) {
		printf("matrix is not unimodular\n");
		return(NULL);
	}
	i = UNIMODULAR(P, Q, R, S);
	return(CHANGE(i));
}
