
#ifndef __ROOTS_H__
#define __ROOTS_H__
#include "integer.h"
#include "stack.h"



typedef struct _Interval {
	MPR *left;
	MPR *right;
} *Interval;

/* Returns a newly created Interval if there is enough memory.
* Returns NULL on error.
*/
Interval createInterval(MPR* left, MPR *right);

/* Frees an interval previously created with createInterval
* Undefined behaviour occurs if 'intvl' this is not true.
*/
void freeInterval(Interval intvl);

/*
* Calculates an upperbound on the positve roots of the supplied polynomial
* using Cauchy's Rule. See Akritas, Elements of Computer Algebra.
*/
MPR *CAUCHY(POLYI P);

/* Using a bisection method this finds the integer part root of the
* supplied polynomial in the supplied _OPEN_ interval.  This function is
* only guaranteed to return the correct value if the interval
* supplied contains exactly one root, and this root does not fall at the
* end points. */
MPI *rootInIntervalI(POLYI P, MPI *LEFT, MPI *RIGHT);
/* Using a bisection method this finds the integer part root of the
* supplied polynomial in the supplied _OPEN_ interval.  This function is
* only guaranteed to return the correct value if the interval
* supplied contains exactly one root, and this root does not fall at the
* end points. */
MPI *rootInIntervalR(POLYI P, MPR *LEFT, MPR *RIGHT);


/* A wrapper function that takes a single polynomial as its argument.
* (Handled by parser).  Returns open rational intervals that are
* guaranteed to enclose the real roots of the supplied polynomial.
*/
void STURM_W(Stack s);

/* Finds the continued fraction expansion of all real roots of the supplied
* polynomial using Lagrange's method and methods first presented in a paper
* by Cantor, Galyean and Zimmer called A Continued Fraction Algorithm for
* Real Algebraic Number
*/
void ROOTEXPANSION(Stack s);


#endif

