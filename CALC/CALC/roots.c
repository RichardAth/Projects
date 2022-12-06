#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_NONSTDC_NO_DEPRECATE
#include "roots.h"
#include "integer.h"
#include <stdio.h>
#include <stdlib.h>
#include "fun.h"
#include "stack.h"

/* This code was written by Sean Seefried when he was a vacation student
* with me at UQ Maths. Dept.
*/
/* Simple linked list implementation of polynomial sequence */
typedef struct _PolySeq {
	POLYI poly;
	struct _PolySeq *next;
} *PolySeq;

/*
* Adds a POLYI to the end of a PolySeq.  If supplied PolySeq is NULL then a
*  new PolySeq is created.
*/
PolySeq addPolySeqTerm(PolySeq ps, POLYI P)
{
	PolySeq CURSOR = ps;
	if (ps == NULL) {
		if ((ps = mmalloc(sizeof(struct _PolySeq))) == NULL)
			return NULL;
		else {
			ps->poly = COPYPI(P);
			ps->next = NULL;
		}
	}
	else {
		while (CURSOR->next)
			CURSOR = CURSOR->next;
		/* CURSOR now points to last term */
		CURSOR->next = mmalloc(sizeof(struct _PolySeq));
		CURSOR->next->poly = COPYPI(P);
		CURSOR->next->next = NULL;
	}
	return ps;
}

/*
* Frees a polynomial sequence previously created with createPolySeq
*/
void freePolySeq(PolySeq ps)
{
	PolySeq CURSOR = ps;
	PolySeq prev;
	while (CURSOR) {
		prev = CURSOR;
		DELETEPI(CURSOR->poly);
		CURSOR = CURSOR->next;
		ffree(prev, sizeof(struct _PolySeq));
	}
}
void printPolySeq(PolySeq ps)
{
	PolySeq CURSOR = ps;
	printf("{");
	if (ps) {
		PRINTPI(CURSOR->poly);
		while (CURSOR->next) {
			printf(", ");
			CURSOR = CURSOR->next;
			PRINTPI(CURSOR->poly);
		}
	}
	printf("}\n");
}

void printPolySeqEval(PolySeq ps, MPR *R)
{
	PolySeq CURSOR = ps;
	int i;
	printf("{");
	if (ps) {
		i = SIGN_AT_R_PI(CURSOR->poly, R);
		printf("%d", i);
		while (CURSOR->next) {
			printf(", ");
			CURSOR = CURSOR->next;
			i = SIGN_AT_R_PI(CURSOR->poly, R);
			printf("%d", i);
		}
	}
	printf("}\n");
}
/*
* Returns the number of sign variations in the Polynomial Sequence at the
* evualated at the supplied value.
* A sign variation exists between two non-zero numbers c[p] and c[q] (p<q)
* of a finite or infinite sequence of numbers, c[1], c[2], c[3],... if the
* following hold:
* i. for q = p + 1, c[p] and c[q] have opposite signs.
* ii. for q >= p +2, the numbers c[p+1], ..., c[q-1] are all zero and c[p]
* and c[q] have opposite signs.
*/


int signVar(PolySeq PS, MPR *R)
{
	int prev = SIGN_AT_R_PI(PS->poly, R);
	/* previous polynomial in sequence's value at I */
	int vars = 0; /* number of variations */
	MPR *Z = ZEROR();
	int sign;
	PolySeq CURSOR = PS;
	while (CURSOR) {
		/* skip polynomials that evaluate to zero */
		sign = SIGN_AT_R_PI(CURSOR->poly, R);
		while (sign == 0 && CURSOR->next) {
			CURSOR = CURSOR->next;
			sign = SIGN_AT_R_PI(CURSOR->poly, R);
		}
		if (prev != 0 && (prev*sign < 0))
			vars++;
		prev = sign;
		CURSOR = CURSOR->next;
	}
	FREEMPR(Z);
	return vars;
}

/*
* Returns a polynomial sequence that is the sturm sequence of the supplied
* polynomial.  See Akritas, Elements of Computer Algebra, p341
*/
PolySeq sturmSeq(POLYI Pptr)
{
	POLYI R1 = COPYPI(Pptr);
	POLYI R2 = DERIVPI(R1);
	POLYI TR, R; /* temp remainder, remainder, modified remainder */
	MPI *CONT;
	PolySeq ps = NULL;

	ps = addPolySeqTerm(ps, R1);
	/* It is possible that Pptr was a constant polynomial. Hence its derivative
	* is zero or NULL.  The Sturm sequence should consists of one term - Pptr.
	*/
	if (R2 != NULL)
		ps = addPolySeqTerm(ps, R2);
	else {
		DELETEPI(R1);
		return ps;
	}
	while ((DEGREEPI(TR = MODPI(R1, R2)) != 0) && TR != NULL) {
		CONT = CONTENTPI2(TR);
		if (SIGNI(CONT) > 0) {
			POLYI TMPI;
			R = PRIMITIVEPI(TR);
			DELETEPI(TR); /* free TR. not needed */

			TMPI = R;
			R = MINUSPI(R);
			DELETEPI(TMPI);
		}
		else {
			R = PRIMITIVEPI(TR);
			DELETEPI(TR);
		}
		ps = addPolySeqTerm(ps, R);
		TR = R1;
		R1 = R2;
		R2 = R;
		DELETEPI(TR); /* deletes R1 */
		FREEMPI(CONT);
	}

	if (TR != NULL) { /* can guarantee that this is a constant polynomial */
		MPI *TMP;
		POLYI ONE_P = ONEPI();
		TMP = LEADPI(TR);
		if (SIGNI(TMP) > 0) {
			POLYI NEGONE_P;
			MPI *CONST = MINUS_ONEI();
			NEGONE_P = CONSTANTPI(CONST);
			addPolySeqTerm(ps, NEGONE_P);
			DELETEPI(NEGONE_P);
			FREEMPI(CONST);
		}
		else {
			addPolySeqTerm(ps, ONE_P);
		}
		DELETEPI(ONE_P);
		FREEMPI(TMP);
	}
	DELETEPI(TR);
	DELETEPI(R1);
	DELETEPI(R2);
	return ps;
}

/*
* Calculates an upperbound on the positve roots of the supplied polynomial
* using Cauchy's Rule. See Akritas, Elements of Computer Algebra.
*/
MPR *CAUCHY(POLYI P)
{
	POLYI Q = COPYPI(P);
	POLYI CURSOR;
	MPR *RETVAL = NULL;
	MPI *LEAD;
	MPR *TWO = TWOR();
	int t, deg = DEGREEPI(P);
	/* t a flag, deg is degree of polynomial */
	int  k, k_, k__, lambda = 0, i, j, p, q, r;
	MPI *CI_ = NULL;
	MPR *CN_ = NULL;

	/* if the leading coefficient is less than zero then multiply polynomial by
	* -1 */

	LEAD = LEADPI(Q);
	if (SIGNI(LEAD) < 0) {
		MPI *NEGONE = MINUS_ONEI();
		POLYI TMPI;
		TMPI = Q;
		Q = SCALARPI(NEGONE, Q);
		DELETEPI(TMPI);
		FREEMPI(NEGONE);
	}
	CURSOR = Q;
	while (CURSOR) {
		if (SIGNI(CURSOR->COEF) < 0) /* count negative terms in poly */
			lambda++;
		CURSOR = CURSOR->NEXT;
	}

	k__ = 0;
	j = BINARYB(LEAD) - 1; /* j = log base 2 of LEAD */
	t = 0;
	if (!(lambda == 0 || deg == 0)) {
		CURSOR = Q;
		while (CURSOR) {
			if (SIGNI(CURSOR->COEF) < 0) {
				k = Q->DEG - CURSOR->DEG;
				CI_ = MULT_I(CURSOR->COEF, -lambda);
				i = BINARYB(CI_) - 1; /* i = log base 2 of CI_ */
				p = i - j - 1;
				q = p / k;
				r = p - k*q;
				if (r < 0) {
					r += k;
					q--;
				}
				k_ = q + 1;

				if (r == (k - 1)) {
					MPR *TMPR;
					MPR *POW;
					POW = POWERR_2(TWO, k*k_);
					TMPR = MPI_TO_MPR(Q->COEF);
					CN_ = MULTR(TMPR, POW);
					FREEMPR(POW);
					FREEMPR(TMPR);
					TMPR = MPI_TO_MPR(CI_);
					if (COMPARER(TMPR, CN_) > 0) {
						k_ = k_ + 1;
					}
					FREEMPR(CN_);
					FREEMPR(TMPR);
				}
				if (t == 0 || k_ > k__) {
					k__ = k_;
					t = 1;
				}
				FREEMPI(CI_);
			}
			CURSOR = CURSOR->NEXT;
		}
	}
	RETVAL = POWERR_2(TWO, k__);

	FREEMPR(TWO);
	FREEMPI(LEAD);
	DELETEPI(Q);
	return RETVAL;

}


/* Returns a newly created Interval if possible.
* Returns NULL on error.
*/
Interval createInterval(MPR* left, MPR *right)
{
	Interval intvl;
	if (!(intvl = malloc(sizeof(struct _Interval))))
		return NULL;
	intvl->left = COPYR(left);
	intvl->right = COPYR(right);
	return intvl;
}
/*-------------------------------------------------------------------*/
/* Frees an interval previously created with createInterval
* Undefined behaviour occurs if 'intvl' does not point to a structure of
* this type. */
void freeInterval(Interval intvl)
{
	FREEMPR(intvl->left);
	FREEMPR(intvl->right);
	free(intvl);
}
/*-------------------------------------------------------------------*/

/* Returns a stack of intervals that contain only one root each respectively
* of the supplied polynomial.  If there are no roots, it returns an empty
* stack */
Stack STURM(POLYI P)
{

	Stack Bounds = stackNew();
	Stack Results = stackNew();
	PolySeq sseq; /* sturm sequence */
	MPI *ZI, *V;
	MPR *BR, *ZERO = ZEROR(), *TWO = TWOR();
	unsigned debugCounter = 0;
	unsigned numRoots, pn_flag = 0;

	/* must find square free version of polynomial or else
	* we cannot use the theorm of Sturm.  This is done simply by dividing
	* the supplied polynomial by the gcd of it and it's derivative */
	POLYI D = DERIVPI(P);
	POLYI G = GCDPI(P, D);
	POLYI Pw = DIVPI(P, G);
	DELETEPI(D);
	DELETEPI(G);

	/* if there is a zero at zero then divide polynomial by X */
	ZI = ZEROI();
	V = VALPI(Pw, ZI);
	if (COMPAREI(V, ZI) == 0) {
		MPI *ONE = ONEI();
		POLYI X = NULL, TMPI;
		PINSERTPI(1, ONE, &X, 0);
		TMPI = Pw;
		Pw = DIVPI(Pw, X);
		DELETEPI(TMPI);
		DELETEPI(X);
		FREEMPI(ONE);
	}
	FREEMPI(V);
	FREEMPI(ZI);

	while (pn_flag < 2) {
		POLYI TMPPOLY;
		sseq = sturmSeq(Pw);

		BR = CAUCHY(Pw);
		stackPush(Bounds, createInterval(ZERO, BR));
		FREEMPR(BR);
		/* main loop */

		while (!stackEmpty(Bounds) && debugCounter < 20) {
			Interval intvl;
			intvl = stackPop(Bounds);
			numRoots = signVar(sseq, intvl->left) - signVar(sseq, intvl->right);
			debugCounter++;
			/* DEBUG INFO */

			/*	printf("The number of roots in the interval (");
			PRINTR(intvl->left);
			printf(", ");
			PRINTR(intvl->right);
			printf(") is %d\n", numRoots);
			printPolySeqEval(sseq, intvl->left);
			printPolySeqEval(sseq, intvl->right);         */
			/* DEBUG END */

			if (numRoots == 1) {

				MPR *L, *R;
				L = intvl->left;
				R = intvl->right;
				if (pn_flag == 1) {
					MPR *TMPR;
					TMPR = MINUSR(L); /* swap left and right of interval because
									  * they are now negative values */
					L = MINUSR(R);
					R = TMPR;
				}
				stackPush(Results, createInterval(L, R));
				if (pn_flag == 1) {
					FREEMPR(L);
					FREEMPR(R);
				}
			}
			else if (numRoots > 1) { /* numroots > 1 */
				MPR *T, *CENTRE;
				T = ADDR(intvl->left, intvl->right);
				CENTRE = RATIOR(T, TWO);
				FREEMPR(T);  /* New intervals are (L, (L+R)/2) and ((L+R)/2, R) */
				stackPush(Bounds, createInterval(intvl->left, CENTRE));
				stackPush(Bounds, createInterval(CENTRE, intvl->right));

				/* debug code */
				/*	   printf("This interval has been subdivided into: (");
				PRINTR(intvl->left);
				printf(", ");
				PRINTR(CENTRE);
				printf(") , (");
				PRINTR(CENTRE);
				printf(", ");
				PRINTR(intvl->right);
				printf(")\n");    */
				/* debug code end */
				FREEMPR(CENTRE);
			}
			freeInterval(intvl);
		}
		TMPPOLY = Pw;
		Pw = P_OF_NEG_X(Pw);
		DELETEPI(TMPPOLY);
		freePolySeq(sseq);
		pn_flag++;
	}
	stackFree(&Bounds);

	FREEMPR(ZERO);
	FREEMPR(TWO);
	DELETEPI(Pw);
	return Results;
}
/*-------------------------------------------------------------------*/

/* A wrapper function that takes a single polynomial as its argument.
* (Handled by parser).  Returns open rational intervals that are
* guaranteed to enclose the real roots of the supplied polynomial.
*/
void STURM_W(Stack s)
{
	POLYI P = stackPop(s);
	Stack Results = STURM(P);
	FILE *outfile;

	if (!(outfile = fopen("sturm.out", "w"))) {
		printf("Could not open file sturm.out for writing\n");
		exit(1);
	}

	fprintf(outfile, "The intervals containing only one root for the polynomial:\n");
	FPRINTPI(outfile, P);
	fprintf(outfile, " are:\n");

	DELETEPI(P);
	while (!stackEmpty(Results)) {
		Interval intvl;
		intvl = stackPop(Results);

		printf("(");
		PRINTR(intvl->left);
		printf(", ");
		PRINTR(intvl->right);
		printf(")\n");

		fprintf(outfile, "(");
		FPRINTR(outfile, intvl->left);
		fprintf(outfile, ", ");
		FPRINTR(outfile, intvl->right);
		fprintf(outfile, ")\n");

		freeInterval(intvl);
	}
	fclose(outfile);
	stackFree(&Results);
}


/* Using a bisection method this finds the integer part root of the
* supplied polynomial in the supplied _OPEN_ interval.  This function is
* only guaranteed to return the correct value if the interval
* supplied contains exactly one root, and this root does not fall at the
* end points. */
MPI *rootInIntervalI(POLYI P, MPI *LEFT, MPI *RIGHT)
{
	MPI *TMP, *LOW = COPYI(LEFT), *HIGH = COPYI(RIGHT), *MID;
	MPI *ONE = ONEI();
	MPI *TWO = TWOI();
	TMP = ADDI(LOW, ONE);

	while (COMPAREI(TMP, HIGH) != 0) {
		MPI *HIGHVAL, *MIDVAL;
		FREEMPI(TMP);
		MID = ADDI(LOW, HIGH);

		TMP = MID;
		MID = INTI(MID, TWO);
		FREEMPI(TMP);
		MIDVAL = VALPI(P, MID);
		HIGHVAL = VALPI(P, HIGH);

		if (SIGNI(MIDVAL) != SIGNI(HIGHVAL)) {
			TMP = LOW;
			LOW = MID;
			FREEMPI(TMP);
		}
		else {
			TMP = HIGH;
			HIGH = MID;
			FREEMPI(TMP);
		}

		FREEMPI(MIDVAL);
		FREEMPI(HIGHVAL);
		TMP = ADDI(LOW, ONE);
	}

	FREEMPI(ONE);
	FREEMPI(TWO);
	FREEMPI(TMP);
	FREEMPI(HIGH);
	return LOW;
}

/* Using a bisection method this finds the integer part root of the
* supplied polynomial in the supplied _OPEN_ interval.  This function is
* only guaranteed to return the correct value if the interval
* supplied contains exactly one root, and this root does not fall at the
* end points.
* If RIGHT == NULL this signifies that we are searching for a root in
* the open interval (LEFT, INFINITY).
*/
MPI *rootInIntervalR(POLYI P, MPR *LEFT, MPR *RIGHT)
{
	MPR *TMP, *LOW = COPYR(LEFT), *HIGH, *MID;
	MPR *TWO = TWOR();
	MPI *INT_LOW, *INT_HIGH; /* stores the integer part of LOW and HIGH */

	if (RIGHT != NULL)
		HIGH = COPYR(RIGHT);
	else
		HIGH = CAUCHY(P);

	INT_LOW = INTI(LOW->N, LOW->D);
	INT_HIGH = INTI(HIGH->N, HIGH->D);

	while (COMPAREI(INT_LOW, INT_HIGH) != 0) {
		FREEMPI(INT_LOW);
		FREEMPI(INT_HIGH);
		MID = ADDR(LOW, HIGH);
		TMP = MID;
		MID = RATIOR(MID, TWO);
		FREEMPR(TMP);
		if (SIGN_AT_R_PI(P, MID) != SIGN_AT_R_PI(P, HIGH)) {
			TMP = LOW;
			LOW = MID;
			FREEMPR(TMP);
		}
		else {
			TMP = HIGH;
			HIGH = MID;
			FREEMPR(TMP);
		}
		INT_LOW = INTI(LOW->N, LOW->D);
		INT_HIGH = INTI(HIGH->N, HIGH->D);
	}
	FREEMPR(TWO);
	FREEMPI(INT_HIGH);
	FREEMPR(HIGH);
	FREEMPR(LOW);
	return INT_LOW;
}


/* Finds the continued fraction expansion of all real roots of the supplied
* polynomial using Lagrange's method and methods first presented in a paper
* by Cantor, Galyean and Zimmer called A Continued Fraction Algorithm for
* Real Algebraic Number
*/
void ROOTEXPANSION(Stack s)
{
	POLYI P = stackPop(s);
	MPI *M = stackPop(s);
	POLYI Q, CURSOR;
	/* Making P square free */
	POLYI D = DERIVPI(P);
	POLYI G = GCDPI(P, D);
	POLYI TMPP;

	Stack intvlStack;
	Interval intvl;
	FILE *outfile;
	MPR *ONE = ONER();
	MPI *NUMQUO, *TMPI;
	MPI *ZERO = ZEROI();
	MPR *TMP;
	MPIA LANGFRAC; /* the lagrange partial quotients */
	int i = 0, j, d;
	/* main loop */

	TMPP = P;
	P = DIVPI(P, G);
	DELETEPI(TMPP);
	DELETEPI(G);
	DELETEPI(D);

	intvlStack = STURM(P);


	if (!(outfile = fopen("rootexp.out", "w"))) {
		printf("Error: Could not open file rootexp.out for writing\n");
		exit(1);
	}

	while (!stackEmpty(intvlStack)) {
		MPI *A = NULL; /* next partial quotient */
		MPR *AR = NULL;
		MPR *RN, *RP = NULL; /* next and previous left bound */
		MPR *SN, *SP = NULL; /* next and previous right bound */
		MPIA CONTFRAC = BUILDMPIA(); /* the zimmer partial quotients */

		MPIA COEF; /* coefficients for the transformation */
				   /* SN will equal NULL to represent INFINITY */
		intvl = stackPop(intvlStack);

		Q = COPYPI(P);
		RN = COPYR(intvl->left);
		SN = COPYR(intvl->right);
		i = 0;
		while (COMPARER(RN, ONE) != 0 || SN != NULL) {
			A = rootInIntervalR(Q, RN, SN);
			AR = MPI_TO_MPR(A);

			ADD_TO_MPIA(CONTFRAC, A, i++);

			FREEMPR(RP); /* ok even if null */
			FREEMPR(SP);

			RP = RN;
			SP = SN;

			TMP = ADDR(AR, ONE);
			if (SP != NULL && COMPARER(SP, TMP) < 0) {
				MPR *TMP2;
				RN = SUBR(SP, AR);
				TMP2 = RN;
				RN = INVERSER(RN);
				FREEMPR(TMP2);
			}
			else
				RN = ONER();
			FREEMPR(TMP);

			if (COMPARER(RP, AR) > 0) {
				MPR *TMP2;
				SN = SUBR(RP, AR);
				TMP2 = SN;
				SN = INVERSER(SN);
				FREEMPR(TMP2);
			}
			else
				SN = NULL;


			P_OF_X_PLUS(&Q, A);

			FREEMPI(A);
			FREEMPR(AR);

			/* This performs transformation p(x) = x^d*p(1/x +a) */
			COEF = BUILDMPIA();
			d = DEGREEPI(Q);
			for (j = 0; j <= d; j++)
				ADD_TO_MPIA(COEF, ZERO, j);

			CURSOR = Q;
			while (CURSOR) {
				ADD_TO_MPIA(COEF, CURSOR->COEF, DEGREEPI(CURSOR));
				CURSOR = CURSOR->NEXT;
			}

			DELETEPI(Q);
			Q = NULL;
			for (j = 0; j <= d; j++)
				PINSERTPI(d - j, COEF->A[j], &Q, 0);
			PURGEPI(&Q);
			FREEMPIA(COEF);

			/* end of transformation */

		} /* the zimmer method is finished. We can use Lagrange's method now */
		FREEMPR(RN);
		FREEMPR(RP);
		FREEMPR(SN);
		FREEMPR(SP);

		if (SIGNI(Q->COEF) == -1) {
			POLYI TMPPOLY;
			TMPPOLY = Q;
			Q = MINUSPI(Q);
			DELETEPI(TMPPOLY);
		}

		/* this bit runs lagrange */
		NUMQUO = CHANGEI(i);
		TMPI = NUMQUO;
		NUMQUO = SUBI(M, NUMQUO);
		FREEMPI(TMPI);
		LAGRANGE(Q, &LANGFRAC, NUMQUO);
		FREEMPI(NUMQUO);

		/* now to output the continued fractions */
		fprintf(outfile, "Continued fraction expansion of roots for the polynomial: \n");
		FPRINTPI(outfile, P);
		fprintf(outfile, "\n");
		printf("The root in the interval (");
		PRINTR(intvl->left);
		printf(", ");
		PRINTR(intvl->right);
		printf(") has the continued fraction expansion: \n");
		fprintf(outfile, "The root in the interval (");
		FPRINTR(outfile, intvl->left);
		fprintf(outfile, ", ");
		FPRINTR(outfile, intvl->right);
		fprintf(outfile, ") has the continued fraction expansion: \n");
		freeInterval(intvl);
		for (j = 0; j < CONTFRAC->size; j++) {
			printf("[%d]: ", j);
			PRINTI(CONTFRAC->A[j]);
			printf("\n");

			fprintf(outfile, "[%d]: ", j);
			FPRINTI(outfile, CONTFRAC->A[j]);
			fprintf(outfile, "\n");
		}
		for (j = 0; j< LANGFRAC->size; j++) {
			printf("[%d]: ", j + CONTFRAC->size);
			PRINTI(LANGFRAC->A[j]);
			printf("\n");

			fprintf(outfile, "[%d]: ", j);
			FPRINTI(outfile, LANGFRAC->A[j]);
			fprintf(outfile, "\n");

		}
		printf("--------------------------------------------------\n");
		fprintf(outfile, "--------------------------------------------------\n");


		FREEMPIA(CONTFRAC);
		FREEMPIA(LANGFRAC);
		DELETEPI(Q);
	}  /* end of main loop */
	fclose(outfile);
	stackFree(&intvlStack);
	FREEMPR(ONE);
	FREEMPI(ZERO);
	DELETEPI(P);
	FREEMPI(M);
}

USI SIGN_COUNT(MPIA A) {
	/* This counts the number of sign changes in the array A */
	USI i, j, n, changes;
	MPIA S;

	j = 0;
	n = A->size;
	S = BUILDMPIA();
	for (i = 0; i<n; i++) {
		if ((A->A[i]->S) == 0) {
			continue;
		}
		else {
			ADD_TO_MPIA(S, A->A[i], j);
			j = j + 1;
		}
	}
	changes = 0;
	for (i = 0; i<j; i++) {
		if ((S->A[i]->S)*(S->A[i + 1]->S)<0) {
			changes = changes + 1;
		}
	}
	FREEMPIA(S);
	return(changes);

}

MPI *STURM_SEQUENCE(POLYI P, MPI *B, MPI *E) {
	POLYI TEMPPOLYI, PA, PB;
	MPIA VALUES;
	int n;
	MPI *TEMP, *MINUSONE;
	USI count, i, changes;

	MINUSONE = MINUS_ONEI();
	VALUES = BUILDMPIA();
	PA = PRIMITIVEPI_(P);
	if (E->S) {
		printf("STURM[0]=");
		PRINTPI(PA);
		printf("\n");
	}
	TEMP = VALPI(PA, B);
	ADD_TO_MPIA(VALUES, TEMP, 0);
	FREEMPI(TEMP);

	TEMPPOLYI = DERIVPI(PA);
	PB = PRIMITIVEPI_(TEMPPOLYI);
	DELETEPI(TEMPPOLYI);

	if (E->S) {
		printf("STURM[1]=");
		PRINTPI(PB);
		printf("\n");
	}

	TEMP = VALPI(PB, B);
	ADD_TO_MPIA(VALUES, TEMP, 1);
	FREEMPI(TEMP);

	n = DEGREEPI(PB);
	count = 2;
	while (n > 0) {
		TEMPPOLYI = PB;
		PB = MODPI(PA, PB);
		n = DEGREEPI(PB);
		DELETEPI(PA);
		PA = TEMPPOLYI;

		TEMPPOLYI = PB;
		PB = PRIMITIVEPI_(PB);
		DELETEPI(TEMPPOLYI);
		TEMPPOLYI = PB;
		PB = SCALARPI(MINUSONE, PB);
		DELETEPI(TEMPPOLYI);
		if (E->S) {
			printf("STURM[%u]=", count);
			PRINTPI(PB);
			printf("\n");
		}
		TEMP = VALPI(PB, B);
		ADD_TO_MPIA(VALUES, TEMP, count);
		FREEMPI(TEMP);
		count++;
	}
	if (E->S) {
		for (i = 0; i < count; i++) {
			PRINTI(VALUES->A[i]);
			printf(", ");
		}
		printf("\n");
	}
	changes = SIGN_COUNT(VALUES);
	FREEMPIA(VALUES);
	FREEMPI(MINUSONE);
	DELETEPI(PA);
	DELETEPI(PB);
	return(CHANGE(changes));
}

