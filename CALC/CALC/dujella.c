/* dujella.c */
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
MPI *XGLOBAL;
MPI *YGLOBAL;
POLYI globalk1;
POLYI globalx1;
POLYI globaly1;
USL counter = 0;

FILE *outfile;  /* global variable */
USI G(MPI *K, MPI *X, MPI *Y, MPI *C) {
	MPI *D, *TEMP, *TEMP1, *R, *S;

	TEMP = MULTI(K, K);
	D = ADD0_I(TEMP, (USL)1);
	TEMP1 = MULT_I(TEMP, (USI)2);
	R = ADD0_I(TEMP1, (USL)1);
	FREEMPI(TEMP);
	FREEMPI(TEMP1);
	S = MULT_I(K, (USI)2);
	TEMP = MULTI(R, X);
	TEMP1 = MULTI3(D, S, Y);
	FREEMPI(D);
	XGLOBAL = ADDI(TEMP, TEMP1);
	FREEMPI(TEMP);
	FREEMPI(TEMP1);
	TEMP = MULTI(S, X);
	TEMP1 = MULTI(R, Y);
	YGLOBAL = ADDI(TEMP, TEMP1);
	FREEMPI(TEMP);
	FREEMPI(TEMP1);
	FREEMPI(R);
	FREEMPI(S);
	if (RSV(YGLOBAL, C) <= 0) {
		counter = counter + 1;
		PRINTI(YGLOBAL);
		printf(":");
		PRINTI(XGLOBAL);
		printf(":");
		PRINTI(K);
		printf("\n");
		FPRINTI(outfile, YGLOBAL);
		fprintf(outfile, ":");
		FPRINTI(outfile, XGLOBAL);
		fprintf(outfile, ":");
		FPRINTI(outfile, K);
		fprintf(outfile, "\n");
		return(1);
	}
	else {
		FREEMPI(YGLOBAL);
		FREEMPI(XGLOBAL);
		return(0);
	}
}

USI G0(MPI *K, MPI *X, MPI *Y, MPI *C) {
	MPI *D, *TEMP, *TEMP1, *R, *S;

	TEMP = MULTI(K, K);
	D = ADD0_I(TEMP, (USL)1);
	TEMP1 = MULT_I(TEMP, (USI)2);
	R = ADD0_I(TEMP1, (USL)1);
	FREEMPI(TEMP);
	FREEMPI(TEMP1);
	S = MULT_I(K, (USI)2);
	TEMP = MULTI(R, X);
	TEMP1 = MULTI3(D, S, Y);
	FREEMPI(D);
	XGLOBAL = ADDI(TEMP, TEMP1);
	FREEMPI(TEMP);
	FREEMPI(TEMP1);
	TEMP = MULTI(S, X);
	TEMP1 = MULTI(R, Y);
	YGLOBAL = ADDI(TEMP, TEMP1);
	FREEMPI(TEMP);
	FREEMPI(TEMP1);
	FREEMPI(R);
	FREEMPI(S);
	if (RSV(YGLOBAL, C) <= 0) {
		counter = counter + 1;
		FPRINTI(outfile, YGLOBAL);
		/*    fprintf(outfile, ":");
		FPRINTI(outfile, XGLOBAL);
		fprintf(outfile, ":");
		FPRINTI(outfile, K);*/
		fprintf(outfile, "\n");
		return(1);
	}
	else {
		FREEMPI(YGLOBAL);
		FREEMPI(XGLOBAL);
		return(0);
	}
}

/* C version of BC recursive program - prints to screen and to file */
void PM2(MPI *K, MPI *X, MPI *Y, MPI *C) {
	MPI *TEMP, *TEMPX, *TEMPY;
	USL t;
	t = G(K, X, Y, C);
	if (t) {
		TEMPY = YGLOBAL;
		TEMPX = XGLOBAL;
		PM2(YGLOBAL, XGLOBAL, K, C);
		FREEMPI(TEMPY);
		FREEMPI(TEMPX);
	}
	t = G(Y, X, K, C);
	if (t) {
		TEMPY = YGLOBAL;
		TEMPX = XGLOBAL;
		PM2(YGLOBAL, XGLOBAL, Y, C);
		FREEMPI(TEMPY);
		FREEMPI(TEMPX);
	}
	TEMP = MINUSI(Y);
	t = G(K, X, TEMP, C);
	FREEMPI(TEMP);
	if (t) {
		TEMPY = YGLOBAL;
		TEMPX = XGLOBAL;
		PM2(YGLOBAL, XGLOBAL, K, C);
		FREEMPI(TEMPY);
		FREEMPI(TEMPX);
	}
}

/* improved C version of BC program forest - prints to sdout and a file,
* all exceptional solutions (k,x,y) with k<= C */
void EXCEPTIONALS(MPI *C) {
	MPI *T, *TEMP, *TEMP1, *TEMP2, *TEMP3, *TEMP4, *TEMP5, *TEMP6, *TEMP7, *TEMP8;
	MPI *TEMP9, *TEMP10, *X, *K, *TMINUS1;
	USL ecount = 0;
	char buff[20];
	strcpy(buff, "exceptionals.out");
	outfile = fopen(buff, "w");
	T = TWOI();
	while (1) {
		counter = 0;
		TEMP1 = MULTI(T, T);
		TEMP = MULT_I(TEMP1, (USL)2);
		FREEMPI(TEMP1);
		TEMP1 = ADD0_I(TEMP, (USL)1);
		X = MULTI(T, TEMP1);
		FREEMPI(TEMP1);
		K = COPYI(TEMP);
		FREEMPI(TEMP);
		if (RSV(K, C) <= 0) {
			PRINTI(K);
			printf(":");
			PRINTI(X);
			printf(":");
			PRINTI(T);
			printf("\n");
			FPRINTI(outfile, K);
			fprintf(outfile, ":");
			FPRINTI(outfile, X);
			fprintf(outfile, ":");
			FPRINTI(outfile, T);
			fprintf(outfile, "\n");
			ecount = ecount + 1;
		}
		else {
			FREEMPI(X);
			FREEMPI(K);
			FREEMPI(T);
			break;
		}
		PM2(K, X, T, C);
		FREEMPI(X);
		FREEMPI(K);
		ecount = ecount + counter;
		TEMP = T;
		T = ADD0_I(T, (USL)1);
		FREEMPI(TEMP);
	}
	T = TWOI();
	while (1) {
		counter = 0;
		TEMP = MULTI(T, T);
		TEMP1 = ADD0_I(TEMP, (USL)1); /* TEMP +1  */
		TEMP2 = SUBI(TEMP1, T);       /* TEMP -T + 1 */
		TEMP3 = MULT_I(TEMP, (USL)2); /* 2TEMP     */
		TEMP4 = ADD0_I(TEMP3, (USL)1);/* 2TEMP + 1 */
		TEMP5 = MULTI(TEMP2, TEMP4);  /* (TEMP -T + 1)(2TEMP + 1) */
		TMINUS1 = SUB0_I(T, (USL)1);    /* T - 1     */
		TEMP6 = MULTI3(T, TMINUS1, TEMP1); /* T(T-1)(TEMP + 1)  */
		TEMP7 = MULT_I(TEMP6, (USL)2); /* 2T(T-1)(TEMP + 1)  */
		X = ADD0I(TEMP5, TEMP7);
		TEMP8 = MULTI(TEMP4, TMINUS1); /* (2TEMP + 1)(T-1)  */
		TEMP9 = MULTI(TEMP2, T);       /* (TEMP -T + 1)T    */
		TEMP10 = MULT_I(TEMP9, (USL)2);/* (TEMP -T + 1)2T   */
		K = ADD0I(TEMP8, TEMP10);
		FREEMPI(TEMP);
		FREEMPI(TEMP1);
		FREEMPI(TEMP2);
		FREEMPI(TEMP3);
		FREEMPI(TEMP4);
		FREEMPI(TEMP5);
		FREEMPI(TEMP6);
		FREEMPI(TEMP7);
		FREEMPI(TEMP8);
		FREEMPI(TEMP9);
		FREEMPI(TEMP10);
		FREEMPI(TMINUS1);
		if (RSV(K, C) <= 0) {
			PRINTI(K);
			printf(":");
			PRINTI(X);
			printf(":");
			PRINTI(T);
			printf("\n");
			FPRINTI(outfile, K);
			fprintf(outfile, ":");
			FPRINTI(outfile, X);
			fprintf(outfile, ":");
			FPRINTI(outfile, T);
			fprintf(outfile, "\n");
			ecount = ecount + 1;
		}
		else {
			FREEMPI(X);
			FREEMPI(K);
			FREEMPI(T);
			break;
		}
		PM2(K, X, T, C);
		FREEMPI(X);
		FREEMPI(K);
		ecount = ecount + counter;
		TEMP = T;
		T = ADD0_I(T, (USL)1);
		FREEMPI(TEMP);
	}
	T = ONEI();
	while (1) {
		counter = 0;
		TEMP = MULTI(T, T);
		TEMP1 = ADD0_I(TEMP, (USL)1); /* TEMP +1  */
		TEMP2 = ADDI(TEMP1, T);       /* TEMP +T + 1 */
		TEMP3 = MULT_I(TEMP, (USL)2); /* 2TEMP     */
		TEMP4 = ADD0_I(TEMP3, (USL)1);/* 2TEMP + 1 */
		TEMP5 = MULTI(TEMP2, TEMP4);  /* (TEMP -T + 1)(2TEMP + 1) */
		TMINUS1 = ADD0_I(T, (USL)1);    /* T + 1     */
		TEMP6 = MULTI3(T, TMINUS1, TEMP1); /* T(T+1)(TEMP + 1)  */
		TEMP7 = MULT_I(TEMP6, (USL)2); /* 2T(T+1)(TEMP + 1)  */
		X = ADD0I(TEMP5, TEMP7);
		TEMP8 = MULTI(TEMP4, TMINUS1); /* (2TEMP + 1)(T-1)  */
		TEMP9 = MULTI(TEMP2, T);       /* (TEMP+T + 1)T    */
		TEMP10 = MULT_I(TEMP9, (USL)2);/* (TEMP+T + 1)2T   */
		K = ADD0I(TEMP8, TEMP10);
		FREEMPI(TEMP);
		FREEMPI(TEMP1);
		FREEMPI(TEMP2);
		FREEMPI(TEMP3);
		FREEMPI(TEMP4);
		FREEMPI(TEMP5);
		FREEMPI(TEMP6);
		FREEMPI(TEMP7);
		FREEMPI(TEMP8);
		FREEMPI(TEMP9);
		FREEMPI(TEMP10);
		FREEMPI(TMINUS1);
		if (RSV(K, C) <= 0) {
			PRINTI(K);
			printf(":");
			PRINTI(X);
			printf(":");
			PRINTI(T);
			printf("\n");
			FPRINTI(outfile, K);
			fprintf(outfile, ":");
			FPRINTI(outfile, X);
			fprintf(outfile, ":");
			FPRINTI(outfile, T);
			fprintf(outfile, "\n");
			ecount = ecount + 1;
		}
		else {
			FREEMPI(X);
			FREEMPI(K);
			FREEMPI(T);
			break;
		}
		PM2(K, X, T, C);
		FREEMPI(X);
		FREEMPI(K);
		ecount = ecount + counter;
		TEMP = T;
		T = ADD0_I(T, (USL)1);
		FREEMPI(TEMP);
	}
	printf("The number of exceptional solutions found is %lu\n", ecount);
	fprintf(outfile, "The number of exceptional solutions found is %lu\n", ecount);
	fclose(outfile);
	return;
}

/* C version of BC recursive program - prints only to file */
void EXCEPTIONALS0(MPI *C) {
	MPI *T, *TEMP, *TEMP1, *TEMP2, *TEMP3, *TEMP4, *TEMP5, *TEMP6, *TEMP7, *TEMP8;
	MPI *TEMP9, *TEMP10, *X, *K, *TMINUS1;
	USL ecount = 0;
	char buff[20];
	strcpy(buff, "exceptionals0.out");
	outfile = fopen(buff, "w");
	T = TWOI();
	while (1) {
		counter = 0;
		TEMP1 = MULTI(T, T);
		TEMP = MULT_I(TEMP1, (USL)2);
		FREEMPI(TEMP1);
		TEMP1 = ADD0_I(TEMP, (USL)1);
		X = MULTI(T, TEMP1);
		FREEMPI(TEMP1);
		K = COPYI(TEMP);
		FREEMPI(TEMP);
		if (RSV(K, C) <= 0) {
			FPRINTI(outfile, K);
			/*  fprintf(outfile, ":");
			FPRINTI(outfile, X);
			fprintf(outfile, ":");
			FPRINTI(outfile, T);*/
			fprintf(outfile, "\n");
			ecount = ecount + 1;
		}
		else {
			FREEMPI(X);
			FREEMPI(K);
			FREEMPI(T);
			break;
		}
		PM20(K, X, T, C);
		FREEMPI(X);
		FREEMPI(K);
		ecount = ecount + counter;
		TEMP = T;
		T = ADD0_I(T, (USL)1);
		FREEMPI(TEMP);
	}
	T = TWOI();
	while (1) {
		counter = 0;
		TEMP = MULTI(T, T);
		TEMP1 = ADD0_I(TEMP, (USL)1); /* TEMP +1  */
		TEMP2 = SUBI(TEMP1, T);       /* TEMP -T + 1 */
		TEMP3 = MULT_I(TEMP, (USL)2); /* 2TEMP     */
		TEMP4 = ADD0_I(TEMP3, (USL)1);/* 2TEMP + 1 */
		TEMP5 = MULTI(TEMP2, TEMP4);  /* (TEMP -T + 1)(2TEMP + 1) */
		TMINUS1 = SUB0_I(T, (USL)1);    /* T - 1     */
		TEMP6 = MULTI3(T, TMINUS1, TEMP1); /* T(T-1)(TEMP + 1)  */
		TEMP7 = MULT_I(TEMP6, (USL)2); /* 2T(T-1)(TEMP + 1)  */
		X = ADD0I(TEMP5, TEMP7);
		TEMP8 = MULTI(TEMP4, TMINUS1); /* (2TEMP + 1)(T-1)  */
		TEMP9 = MULTI(TEMP2, T);       /* (TEMP -T + 1)T    */
		TEMP10 = MULT_I(TEMP9, (USL)2);/* (TEMP -T + 1)2T   */
		K = ADD0I(TEMP8, TEMP10);
		FREEMPI(TEMP);
		FREEMPI(TEMP1);
		FREEMPI(TEMP2);
		FREEMPI(TEMP3);
		FREEMPI(TEMP4);
		FREEMPI(TEMP5);
		FREEMPI(TEMP6);
		FREEMPI(TEMP7);
		FREEMPI(TEMP8);
		FREEMPI(TEMP9);
		FREEMPI(TEMP10);
		FREEMPI(TMINUS1);
		if (RSV(K, C) <= 0) {
			FPRINTI(outfile, K);
			/*  fprintf(outfile, ":");
			FPRINTI(outfile, X);
			fprintf(outfile, ":");
			FPRINTI(outfile, T);*/
			fprintf(outfile, "\n");
			ecount = ecount + 1;
		}
		else {
			FREEMPI(X);
			FREEMPI(K);
			FREEMPI(T);
			break;
		}
		PM20(K, X, T, C);
		FREEMPI(X);
		FREEMPI(K);
		ecount = ecount + counter;
		TEMP = T;
		T = ADD0_I(T, (USL)1);
		FREEMPI(TEMP);
	}
	T = ONEI();
	while (1) {
		counter = 0;
		TEMP = MULTI(T, T);
		TEMP1 = ADD0_I(TEMP, (USL)1); /* TEMP +1  */
		TEMP2 = ADDI(TEMP1, T);       /* TEMP +T + 1 */
		TEMP3 = MULT_I(TEMP, (USL)2); /* 2TEMP     */
		TEMP4 = ADD0_I(TEMP3, (USL)1);/* 2TEMP + 1 */
		TEMP5 = MULTI(TEMP2, TEMP4);  /* (TEMP -T + 1)(2TEMP + 1) */
		TMINUS1 = ADD0_I(T, (USL)1);    /* T + 1     */
		TEMP6 = MULTI3(T, TMINUS1, TEMP1); /* T(T+1)(TEMP + 1)  */
		TEMP7 = MULT_I(TEMP6, (USL)2); /* 2T(T+1)(TEMP + 1)  */
		X = ADD0I(TEMP5, TEMP7);
		TEMP8 = MULTI(TEMP4, TMINUS1); /* (2TEMP + 1)(T-1)  */
		TEMP9 = MULTI(TEMP2, T);       /* (TEMP+T + 1)T    */
		TEMP10 = MULT_I(TEMP9, (USL)2);/* (TEMP+T + 1)2T   */
		K = ADD0I(TEMP8, TEMP10);
		FREEMPI(TEMP);
		FREEMPI(TEMP1);
		FREEMPI(TEMP2);
		FREEMPI(TEMP3);
		FREEMPI(TEMP4);
		FREEMPI(TEMP5);
		FREEMPI(TEMP6);
		FREEMPI(TEMP7);
		FREEMPI(TEMP8);
		FREEMPI(TEMP9);
		FREEMPI(TEMP10);
		FREEMPI(TMINUS1);
		if (RSV(K, C) <= 0) {
			FPRINTI(outfile, K);
			/*     fprintf(outfile, ":");
			FPRINTI(outfile, X);
			fprintf(outfile, ":");
			FPRINTI(outfile, T);*/
			fprintf(outfile, "\n");
			ecount = ecount + 1;
		}
		else {
			FREEMPI(X);
			FREEMPI(K);
			FREEMPI(T);
			break;
		}
		PM20(K, X, T, C);
		FREEMPI(X);
		FREEMPI(K);
		ecount = ecount + counter;
		TEMP = T;
		T = ADD0_I(T, (USL)1);
		FREEMPI(TEMP);
	}
	fprintf(outfile, "The number of exceptional solutions found is %lu\n", ecount);
	fclose(outfile);
	return;
}

void PM20(MPI *K, MPI *X, MPI *Y, MPI *C) {
	MPI *TEMP, *TEMPX, *TEMPY;
	USL t;
	t = G0(K, X, Y, C);
	if (t) {
		TEMPY = YGLOBAL;
		TEMPX = XGLOBAL;
		PM20(YGLOBAL, XGLOBAL, K, C);
		FREEMPI(TEMPY);
		FREEMPI(TEMPX);
	}
	t = G0(Y, X, K, C);
	if (t) {
		TEMPY = YGLOBAL;
		TEMPX = XGLOBAL;
		PM20(YGLOBAL, XGLOBAL, Y, C);
		FREEMPI(TEMPY);
		FREEMPI(TEMPX);
	}
	TEMP = MINUSI(Y);
	t = G0(K, X, TEMP, C);
	FREEMPI(TEMP);
	if (t) {
		TEMPY = YGLOBAL;
		TEMPX = XGLOBAL;
		PM20(YGLOBAL, XGLOBAL, K, C);
		FREEMPI(TEMPY);
		FREEMPI(TEMPX);
	}
}

/* C version of BC recursive program - prints only to file */
/* finds all Type 1 exceptional solutions (k,x,y) with k <= c */
void EXCEPTIONALS10(MPI *C) {
	MPI *T, *TEMP, *TEMP1, *TEMP2, *TEMP3, *TEMP4, *TEMP5, *TEMP6, *TEMP7, *TEMP8;
	MPI *TEMP9, *TEMP10, *X, *K, *TMINUS1;
	USL ecount = 0;
	char buff[20];
	strcpy(buff, "exceptionals10.out");
	outfile = fopen(buff, "w");
	T = TWOI();
	while (1) {
		counter = 0;
		TEMP1 = MULTI(T, T);
		TEMP = MULT_I(TEMP1, (USL)2);
		FREEMPI(TEMP1);
		TEMP1 = ADD0_I(TEMP, (USL)1);
		X = MULTI(T, TEMP1);
		FREEMPI(TEMP1);
		K = COPYI(TEMP);
		FREEMPI(TEMP);
		if (RSV(K, C) <= 0) {
			FPRINTI(outfile, K);
			fprintf(outfile, ":");
			FPRINTI(outfile, X);
			fprintf(outfile, ":");
			FPRINTI(outfile, T);
			fprintf(outfile, "\n");
			ecount = ecount + 1;
		}
		else {
			FREEMPI(X);
			FREEMPI(K);
			FREEMPI(T);
			break;
		}
		PM30(K, X, T, C);
		FREEMPI(X);
		FREEMPI(K);
		ecount = ecount + counter;
		TEMP = T;
		T = ADD0_I(T, (USL)1);
		FREEMPI(TEMP);
	}
	T = TWOI();
	while (1) {
		counter = 0;
		TEMP = MULTI(T, T);
		TEMP1 = ADD0_I(TEMP, (USL)1); /* TEMP +1  */
		TEMP2 = SUBI(TEMP1, T);       /* TEMP -T + 1 */
		TEMP3 = MULT_I(TEMP, (USL)2); /* 2TEMP     */
		TEMP4 = ADD0_I(TEMP3, (USL)1);/* 2TEMP + 1 */
		TEMP5 = MULTI(TEMP2, TEMP4);  /* (TEMP -T + 1)(2TEMP + 1) */
		TMINUS1 = SUB0_I(T, (USL)1);    /* T - 1     */
		TEMP6 = MULTI3(T, TMINUS1, TEMP1); /* T(T-1)(TEMP + 1)  */
		TEMP7 = MULT_I(TEMP6, (USL)2); /* 2T(T-1)(TEMP + 1)  */
		X = ADD0I(TEMP5, TEMP7);
		TEMP8 = MULTI(TEMP4, TMINUS1); /* (2TEMP + 1)(T-1)  */
		TEMP9 = MULTI(TEMP2, T);       /* (TEMP -T + 1)T    */
		TEMP10 = MULT_I(TEMP9, (USL)2);/* (TEMP -T + 1)2T   */
		K = ADD0I(TEMP8, TEMP10);
		FREEMPI(TEMP);
		FREEMPI(TEMP1);
		FREEMPI(TEMP2);
		FREEMPI(TEMP3);
		FREEMPI(TEMP4);
		FREEMPI(TEMP5);
		FREEMPI(TEMP6);
		FREEMPI(TEMP7);
		FREEMPI(TEMP8);
		FREEMPI(TEMP9);
		FREEMPI(TEMP10);
		FREEMPI(TMINUS1);
		if (RSV(K, C) <= 0) {
			FPRINTI(outfile, K);
			fprintf(outfile, ":");
			FPRINTI(outfile, X);
			fprintf(outfile, ":");
			FPRINTI(outfile, T);
			fprintf(outfile, "\n");
			ecount = ecount + 1;
		}
		else {
			FREEMPI(X);
			FREEMPI(K);
			FREEMPI(T);
			break;
		}
		PM30(K, X, T, C);
		FREEMPI(X);
		FREEMPI(K);
		ecount = ecount + counter;
		TEMP = T;
		T = ADD0_I(T, (USL)1);
		FREEMPI(TEMP);
	}
	T = ONEI();
	while (1) {
		counter = 0;
		TEMP = MULTI(T, T);
		TEMP1 = ADD0_I(TEMP, (USL)1); /* TEMP +1  */
		TEMP2 = ADDI(TEMP1, T);       /* TEMP +T + 1 */
		TEMP3 = MULT_I(TEMP, (USL)2); /* 2TEMP     */
		TEMP4 = ADD0_I(TEMP3, (USL)1);/* 2TEMP + 1 */
		TEMP5 = MULTI(TEMP2, TEMP4);  /* (TEMP -T + 1)(2TEMP + 1) */
		TMINUS1 = ADD0_I(T, (USL)1);    /* T + 1     */
		TEMP6 = MULTI3(T, TMINUS1, TEMP1); /* T(T+1)(TEMP + 1)  */
		TEMP7 = MULT_I(TEMP6, (USL)2); /* 2T(T+1)(TEMP + 1)  */
		X = ADD0I(TEMP5, TEMP7);
		TEMP8 = MULTI(TEMP4, TMINUS1); /* (2TEMP + 1)(T-1)  */
		TEMP9 = MULTI(TEMP2, T);       /* (TEMP+T + 1)T    */
		TEMP10 = MULT_I(TEMP9, (USL)2);/* (TEMP+T + 1)2T   */
		K = ADD0I(TEMP8, TEMP10);
		FREEMPI(TEMP);
		FREEMPI(TEMP1);
		FREEMPI(TEMP2);
		FREEMPI(TEMP3);
		FREEMPI(TEMP4);
		FREEMPI(TEMP5);
		FREEMPI(TEMP6);
		FREEMPI(TEMP7);
		FREEMPI(TEMP8);
		FREEMPI(TEMP9);
		FREEMPI(TEMP10);
		FREEMPI(TMINUS1);
		if (RSV(K, C) <= 0) {
			FPRINTI(outfile, K);
			fprintf(outfile, ":");
			FPRINTI(outfile, X);
			fprintf(outfile, ":");
			FPRINTI(outfile, T);
			fprintf(outfile, "\n");
			ecount = ecount + 1;
		}
		else {
			FREEMPI(X);
			FREEMPI(K);
			FREEMPI(T);
			break;
		}
		PM30(K, X, T, C);
		FREEMPI(X);
		FREEMPI(K);
		ecount = ecount + counter;
		TEMP = T;
		T = ADD0_I(T, (USL)1);
		FREEMPI(TEMP);
	}
	fprintf(outfile, "The number of Type 1 exceptional solutions found is %lu\n", ecount);
	fclose(outfile);
	return;
}

void PM30(MPI *K, MPI *X, MPI *Y, MPI *C) {
	MPI *TEMPX, *TEMPY;
	USL t;
	t = G0(Y, X, K, C);
	if (t) {
		TEMPY = YGLOBAL;
		TEMPX = XGLOBAL;
		PM30(YGLOBAL, XGLOBAL, Y, C);
		FREEMPI(TEMPY);
		FREEMPI(TEMPX);
	}
}

void GPLUS(POLYI k, POLYI x, POLYI y) {
	/* 23rd November 2012.
	* (k1,x1,y1)=g+(k,x,y) */
	POLYI D, TEMP, TEMP1, R, S, ONE, TEMP2;

	TEMP = MULTPI(k, k);
	ONE = ONEPI();
	D = ADDPI(TEMP, ONE);
	TEMP1 = ADDPI(TEMP, TEMP);
	R = ADDPI(TEMP1, ONE);
	DELETEPI(TEMP);
	DELETEPI(TEMP1);
	S = ADDPI(k, k);
	TEMP = MULTPI(R, x);
	TEMP1 = MULTPI(D, S);
	TEMP2 = MULTPI(TEMP1, y);
	DELETEPI(TEMP1);
	DELETEPI(D);
	DELETEPI(ONE);
	globalx1 = ADDPI(TEMP, TEMP2);
	DELETEPI(TEMP);
	DELETEPI(TEMP2);
	TEMP = MULTPI(S, x);
	TEMP1 = MULTPI(R, y);
	globalk1 = ADDPI(TEMP, TEMP1);
	DELETEPI(TEMP);
	DELETEPI(TEMP1);
	DELETEPI(R);
	DELETEPI(S);
	globaly1 = COPYPI(k);
	return;
}

void GZERO(POLYI k, POLYI x, POLYI y) {
	/* 23rd November 2012.
	* (k1,x1,y1)=g0(k,x,y) */
	GPLUS(y, x, k);
}

void GMINUS(POLYI k, POLYI x, POLYI y) {
	/* 23rd November 2012.
	* (k1,x1,y1)=g-(k,x,y) */
	POLYI minusy;

	minusy = MINUSPI(y);
	GPLUS(k, x, minusy);
	DELETEPI(minusy);
	return;
}

/* This constructs the polynomial triple (k(t),x(t),y(t)) corresponding to a path
* given by the reverse digits of the negative 3-adic expansion of N,
* where the root node is (t,t,0) if TYPE=0, (t,t^2-t+1,t-1) if TYPE=-1,
* (t,t^2+t+1,t+1) if TYPE=1. If FLAG is nonzero, we print (k(t),x(t),y(t)), otherwise only k(t).
* Dujella's unicity conjecture is equivalent to the statement that the values of k(t) as t ranges
* over t >= 2 for TYPE =0, t >= 1 for TYPE 1 and t>=2 for TYPE -1, are distinct.
*/

void BRANCH(MPI *N, MPI *TYPE, MPI *FLAG) {
	int i, r, rminus1;
	MPIA R;
	MPI *ZERO;
	POLYI K1, X1, Y1, TEMP, TEMP1, TEMP2;

	r = (int)TERNARY(N, &R);
	rminus1 = r - 1;
	K1 = CREATEPI();
	K1->DEG = 1;
	K1->COEF = ONEI();
	K1->NEXT = NULL;
	if (EQZEROI(TYPE)) {
		X1 = COPYPI(K1);
		ZERO = ZEROI();
		Y1 = CONSTANTPI(ZERO);
		FREEMPI(ZERO);
	}
	if (EQONEI(TYPE) || EQMINUSONEI(TYPE)) {
		TEMP1 = CREATEPI();
		TEMP1->DEG = 2;
		TEMP1->COEF = ONEI();
		TEMP1->NEXT = NULL;
		TEMP2 = CREATEPI();
		TEMP2->DEG = 1;
		TEMP2->COEF = COPYI(TYPE);
		TEMP2->NEXT = NULL;
		TEMP = ADDPI(TEMP1, TEMP2);
		DELETEPI(TEMP1);
		DELETEPI(TEMP2);
		TEMP1 = ONEPI();
		X1 = ADDPI(TEMP, TEMP1);
		DELETEPI(TEMP);
		DELETEPI(TEMP1);

		TEMP2 = CONSTANTPI(TYPE);
		Y1 = ADDPI(K1, TEMP2);
		DELETEPI(TEMP2);
	}
	for (i = rminus1; i >= 0; i--) {
		if (EQMINUSONEI((R->A)[i])) {
			GMINUS(K1, X1, Y1);
		}
		if (EQZEROI((R->A)[i])) {
			GZERO(K1, X1, Y1);
		}
		if (EQONEI((R->A)[i])) {
			GPLUS(K1, X1, Y1);
		}
		DELETEPI(K1);
		DELETEPI(X1);
		DELETEPI(Y1);
		K1 = globalk1;
		X1 = globalx1;
		Y1 = globaly1;

		PRINTPI(globalk1);
		if (FLAG->S) {
			printf(",");
			PRINTPI(globalx1);
			printf(",");
			PRINTPI(globaly1);
		}
		printf("\n==============\n");
	}
	FREEMPIA(R);
	DELETEPI(K1);
	DELETEPI(X1);
	DELETEPI(Y1);
	return;
}

/* This contructs the ternary expansion of a positive integer, where the digits are 0, 1 or -1.*/
/* Mainly for use in BRANCH(). */
USI TERNARY(MPI *N, MPIA *R) {
	USI i, e;
	MPI *rr, *NN;
	int r;
	MPI *TEMP, *TEMP1;

	e = 1;
	*R = BUILDMPIA();
	NN = COPYI(N);
	i = 0;
	while (NN->S > 0) {
		r = MOD3(NN);
		if (e) {
			printf("%d", r);
		}
		/*if(r == 0){
		rr = ZEROI();
		}else{
		rr = CHANGEI((long)r);
		}*/
		rr = CHANGEI((long)r);
		ADD_TO_MPIA(*R, rr, (USL)i);
		TEMP = SUBI(NN, rr);
		FREEMPI(rr);
		TEMP1 = NN;
		NN = INT0_(TEMP, (USL)3);
		FREEMPI(TEMP);
		FREEMPI(TEMP1);
		i++;
	}
	printf("\n");
	FREEMPI(NN);
	return(i);
}

int MOD3(MPI *N) {
	int r;

	r = (int)MOD0_(N, 3);
	if (r == 2) {
		r = -1;
	}
	return(r);
}
MPI *TERNARYX(MPI *N, MPIA *R) {
	MPI *E;
	MPIA S;
	USI i;
	if (N->S < 0) {
		printf("N<0, so no array returned\n");

		*R = NULL;
		return(NULL);
	}
	else {
		i = TERNARY(N, &S);
		*R = S;
		E = CHANGE(i);
		return(E);
	}
}

MPI *BRANCHX(MPI *N, MPI *TYPE, MPI *FLAG) {

	if (N->S <= 0) {
		printf("n <= 0\n");
		return(NULL);
	}
	if (EQZEROI(TYPE) || EQONEI(TYPE) || EQMINUSONEI(TYPE)) {
		BRANCH(N, TYPE, FLAG);
		return(ONEI());
	}
	else {
		printf("TYPE !- 1, 0 or -1\n");
		return(NULL);
	}
}
