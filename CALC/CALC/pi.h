/* pI.h */
#include "integer.h"

#ifndef __PI_H__
#define __PI_H__

POLYI CREATEPI();
/*
If there is anything in BANKI (the memory bank of TERMI's) then use that,
otherwise allocate a new one.
*/

int DEGREEPI(POLYI Pptr);
/*
* returns the degree of the polynomial Pptr.
*/

MPI *LEADPI(POLYI Pptr);
/*
* *Aptr is the  leading coefficient of the polynomial Pptr.
*/

void DELETEPI(POLYI Pptr);
/*
* returns back to the system storage dynamically allocated to the linked list
* which Pptr points to.
*/

void PINSERTPI(int DEG, MPI *COEF, POLYI *Pptr, USI OPT);
/*
* inserts a copy of a single term *Qptr into polynomial *Pptr.  It
* has two options that control the insertion in case a term of *Pptr
* already has the degree of *Qptr, If OPT = 0, then the term is
* replaced; if OPT = 1, then the coefficient of *Qptr is added to
* that of the corresponding term of *Pptr.  */


void PURGEPI(POLYI *Pptr);
/*
* gets rid of any zero coefficients in the polynomial Pptr.
*/

POLYI COPYPI(POLYI Pptr);
/*
* makes a copy of the polynomial Pptr.
*/

void FPRINTPI(FILE *outfile, POLYI Pptr);
void PRINTPI(POLYI Pptr);
/*
* outputs the terms of the polynomial Pptr.
*/



POLYI SCALARPI(MPI *Cptr, POLYI Pptr);
/*
* multiplying the polynomial Pptr by the MPI *Cptr.
*/


POLYI ADDPI(POLYI Pptr, POLYI Qptr);
/*
* returns the sum of two polynomials.
*/


POLYI SUBPI(POLYI Pptr, POLYI Qptr);
/*
* returns the difference Pptr - Qptr of two polynomials.
*/


POLYI MULTPI(POLYI Pptr, POLYI Qptr);
/*
* returns the product of two polynomials.
*/


MPI *VALPI(POLYI Pptr, MPI *Cptr);
/*
* Evaluating the polynomial Pptr at the MPI *Cptr.
*/


POLYI DIVPI(POLYI Fptr, POLYI Gptr);
/*
* Pseudo-division of polynomials: see Knuth, Vol 2, p.407
* and H. Flanders, Scientific Pascal, p.510.
*/


POLYI MODPI(POLYI Fptr, POLYI Gptr);
/*
* Pseudo-division of polynomials: see Knuth, Vol 2, p.407
* and H. Flanders, Scientific Pascal, p.510.
*/


MPI *CONTENTPI(POLYI Pptr);
/*
* *Cptr is the content of the polynomial Pptr.
*/
MPI *CONTENTPI2(POLYI Pptr);


POLYI PRIMITIVEPI(POLYI Pptr);
/*
* returns the primitive part of the polynomial Pptr.
*/


POLYI GCDPI(POLYI Pptr, POLYI Qptr);
/*
* returns the gcd of polynomials Pptr and Qptr.
* see Knuth, Vol 2, p.403-408 and H. Flanders, Scientific Pascal, p.510.
*/


POLYI DERIVPI(POLYI Pptr);
/*
* returns the derivative of a polynomial.
*/


POLYI *SQUAREFREEPI(POLYI Pptr, USI **D, USI *tptr);
/*
* returns the "squarefree decomposition" of Pptr, a non-constant polynomial:
* Pptr = G[0]^D[0]...G[*tpr - 1]^D[*tptr - 1],
* where G[i] is squarefree and is the product of the irreducible factors of
* multiplicity D[i]. See L. Childs, Introduction to Higher Algebra, p.159.
*/


POLYI POWERPI(POLYI Pptr, USI n);
/*
* returns Pptr ^ n, where 0 <= n < R0 * R0.
*/


POLYI ONEPI();
/*
* returns the constant polynomial 1.
*/
POLYI ZEROPI();

POLYI CONSTANTPI(MPI *Aptr);
/*
* returns the constant polynomial whose constant term is *Aptr.
*/


MPI *SUBRESULTANT(POLYI Pptr, POLYI Qptr);
/*
* *Aptr is the resultant of Pptr and Qptr, found using the subresultant
* algorithm, p. 130. E. Kaltofen, G. E. Collins etc, Algoritm 4.5.
* similar to Knuth, Algorithm C, p. 410.
*/


MPI *DISCRIMINANTPI(POLYI Pptr);
/*
* *Dptr = Discriminant of Pptr = (1/a[n])*(-1)^{n*(n-1)/2 * Res(Pptr, Pptr').
* See O. Perron, Algebra, Vol 1, p.212.
*/


MPI *LENGTHSQPI(POLYI Pptr);
/*
* *Cptr is the sum of the squares of the coefficients of Pptr.
*/


MPI *VAL0PI(POLYI Pptr);
/*
* *Mptr = Pptr(0)
*/

/* Returns respectively positive, or negative if the sign of the polynomial
* evaluated at the supplied rational number is positive or negative.
* Returns zero if polynomial has a root at the supplied rational */
int SIGN_AT_R_PI(POLYI P, MPR *R);

/* If P = p(x) then this function returns p(-x) */
POLYI P_OF_NEG_X(POLYI P);

/* Returns -1*P(x) */
POLYI MINUSPI(POLYI P);
void P_OF_X_PLUS(POLYI *P, MPI *A);

#endif


