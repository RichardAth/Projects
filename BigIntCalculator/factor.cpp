/*
This file is part of Alpertron Calculators.
Copyright 2015 Dario Alejandro Alpern
Alpertron Calculators is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
Alpertron Calculators is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with Alpertron Calculators.  If not, see <http://www.gnu.org/licenses/>.
*/
#define  _CRT_SECURE_NO_DEPRECATE
//#define FACTORIZATION_APP
#define __EMSCRIPTEN__

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <stdint.h>
#include <math.h>
#include <assert.h>
#include "showtime.h"
#include "bignbr.h"
#include "expression.h"
#include "factor.h"

int yieldFreq;

char lowerText[30000];
char *ptrLowerText;
extern mmCback modmultCallback;
extern long long lModularMult;
extern char *ptrInputText;

#define TYP_SIQS   200000000

#define MAX_PRIME_SIEVE 7  // Only numbers 7 or 11 are accepted here.
#if MAX_PRIME_SIEVE == 11
#define SIEVE_SIZE (2*3*5*7*11)
#define GROUP_SIZE ((2-1)*(3-1)*(5-1)*(7-1)*(11-1))
#else
#define SIEVE_SIZE (2*3*5*7)
#define GROUP_SIZE ((2-1)*(3-1)*(5-1)*(7-1))
#endif
#define HALF_SIEVE_SIZE (SIEVE_SIZE/2)

static int nbrPrimes, indexPrimes, StepECM, DegreeAurif;
static int FactorIndex;
static BigInteger power, prime;

static int indexM, maxIndexM;
static bool foundByLehman;
static bool performLehman;
static int SmallPrime[670] = { 0 }; /* Primes < 5000 */
//static int NextEC;
static int EC;            // Elliptic Curve Number
static limb A0[MAX_LEN];
static limb A02[MAX_LEN];
static limb A03[MAX_LEN];
static limb AA[MAX_LEN];
static limb DX[MAX_LEN];
static limb DZ[MAX_LEN];
static limb GD[MAX_LEN];
static limb M[MAX_LEN];
static limb TX[MAX_LEN];
static limb TZ[MAX_LEN];
static limb UX[MAX_LEN];
static limb UZ[MAX_LEN];
static limb W1[MAX_LEN];
static limb W2[MAX_LEN];
static limb W3[MAX_LEN];
static limb W4[MAX_LEN];
static limb WX[MAX_LEN];
static limb WZ[MAX_LEN];
static limb X[MAX_LEN];
static limb Z[MAX_LEN];
static limb Aux1[MAX_LEN];
static limb Aux2[MAX_LEN];
static limb Aux3[MAX_LEN];
static limb Aux4[MAX_LEN];
static limb Xaux[MAX_LEN];
static limb Zaux[MAX_LEN];
static limb root[GROUP_SIZE][MAX_LEN];
static unsigned char sieve[10 * SIEVE_SIZE];
static unsigned char sieve2310[SIEVE_SIZE];
static int sieveidx[GROUP_SIZE];
static limb GcdAccumulated[MAX_LEN];

static limb *fieldAA, *fieldTX, *fieldTZ, *fieldUX, *fieldUZ;
static void add3(limb *x3, limb *z3, limb *x2, limb *z2, limb *x1, limb *z1, limb *x, limb *z);
static void duplicate(limb *x2, limb *z2, limb *x1, limb *z1);

static void insertBigFactor(struct sFactors *pstFactors, BigInteger &divisor);
void ValuestoZ(Znum &numberZ, const int number[]);

static BigInteger Temp1, Temp2, Temp3, Temp4;

long long Gamma[386];
long long Delta[386];
long long AurifQ[386];

extern int q[MAX_LEN];
static int nbrToFactor[MAX_LEN];
struct sFactors astFactorsMod[Max_Factors];   // maximum unique factors

static BigInteger Quad1, Quad2, Quad3, Quad4;


enum eEcmResult
{
	FACTOR_NOT_FOUND = 0,
	FACTOR_FOUND,
	CHANGE_TO_SIQS,
	ERROR
};

/* ECM limits for 30, 35, ..., 95 digits */
static int limits[] = { 10, 10, 10, 10, 10, 15, 22, 26, 35, 50, 100, 150, 250 };
/******************************************************/
/* Start of code adapted from Paul Zimmermann's ECM4C */
/******************************************************/
#define ADD 6  /* number of multiplications in an addition */
#define DUP 5  /* number of multiplications in a duplicate */

/* uses global var NumberLength, returns value in YieldFrequency */
static void GetYieldFrequency(void)
{
	yieldFreq = 1000000 / (NumberLength * NumberLength) + 1;
	if (yieldFreq > 100000) {
		yieldFreq = yieldFreq / 100000 * 100000;
	}
	else if (yieldFreq > 10000) {
		yieldFreq = yieldFreq / 10000 * 10000;
	}
	else if (yieldFreq > 1000) {
		yieldFreq = yieldFreq / 1000 * 1000;
	}
	else if (yieldFreq > 100) {
		yieldFreq = yieldFreq / 100 * 100;
	}
}

/* returns the number of modular multiplications */
static int lucas_cost(int n, double v)
{
	int c, d, e, r;

	d = n;
	r = (int)((double)d / v + 0.5);
	if (r >= n) {
		return (ADD * n);
	}
	d = n - r;
	e = 2 * r - n;
	c = DUP + ADD; /* initial duplicate and final addition */
	while (d != e) {
		if (d < e) {
			r = d;
			d = e;
			e = r;
		}
		if (4 * d <= 5 * e && ((d + e) % 3) == 0) { /* condition 1 */
			r = (2 * d - e) / 3;
			e = (2 * e - d) / 3;
			d = r;
			c += 3 * ADD; /* 3 additions */
		}
		else if (4 * d <= 5 * e && (d - e) % 6 == 0) { /* condition 2 */
			d = (d - e) / 2;
			c += ADD + DUP; /* one addition, one duplicate */
		} 
		else if (d <= (4 * e)) { /* condition 3 */
			d -= e;
			c += ADD; /* one addition */
		}
		else if ((d + e) % 2 == 0) { /* condition 4 */
			d = (d - e) / 2;
			c += ADD + DUP; /* one addition, one duplicate */
		}
		else if (d % 2 == 0) { /* condition 5 */
			d /= 2;
			c += ADD + DUP; /* one addition, one duplicate */
		}
		else if (d % 3 == 0) { /* condition 6 */
			d = d / 3 - e;
			c += 3 * ADD + DUP; /* three additions, one duplicate */
		}
		else if ((d + e) % 3 == 0) { /* condition 7 */
			d = (d - 2 * e) / 3;
			c += 3 * ADD + DUP; /* three additions, one duplicate */
		}
		else if ((d - e) % 3 == 0) { /* condition 8 */
			d = (d - e) / 3;
			c += 3 * ADD + DUP; /* three additions, one duplicate */
		}
		else if (e % 2 == 0) { /* condition 9 */
			e /= 2;
			c += ADD + DUP; /* one addition, one duplicate */
		}
	}
	return c;
}

/* computes nP from P=(x:z) and puts the result in (x:z). Assumes n>2. */
void prac(int n, limb *x, limb *z, limb *xT, limb *zT, limb *xT2, limb *zT2)
{
	int d, e, r, i;
	limb *t;
	limb *xA = x, *zA = z;
	limb *xB = Aux1, *zB = Aux2;
	limb *xC = Aux3, *zC = Aux4;
	double v[] =
	{
		1.61803398875,
		1.72360679775,
		1.618347119656,
		1.617914406529,
		1.612429949509,
		1.632839806089,
		1.620181980807,
		1.580178728295,
		1.617214616534,
		1.38196601125 };

	/* chooses the best value of v */
	r = lucas_cost(n, v[0]);
	i = 0;

	for (d = 1; d < 10; d++) {
		e = lucas_cost(n, v[d]);
		if (e < r) {
			r = e;
			i = d;
		}
	}
	d = n;
	r = (int)((double)d / v[i] + 0.5);
	/* first iteration always begins by Condition 3, then a swap */
	d = n - r;
	e = 2 * r - n;
	memcpy(xB, xA, NumberLength * sizeof(limb));   // B <- A
	memcpy(zB, zA, NumberLength * sizeof(limb));
	memcpy(xC, xA, NumberLength * sizeof(limb));   // C <- A
	memcpy(zC, zA, NumberLength * sizeof(limb));
	duplicate(xA, zA, xA, zA); /* A=2*A */
	while (d != e) {
		if (d < e) {
			r = d;
			d = e;
			e = r;
			t = xA;
			xA = xB;
			xB = t;
			t = zA;
			zA = zB;
			zB = t;
		}
		/* do the first line of Table 4 whose condition qualifies */
		if (4 * d <= 5 * e && ((d + e) % 3) == 0) { /* condition 1 */
			r = (2 * d - e) / 3;
			e = (2 * e - d) / 3;
			d = r;
			add3(xT, zT, xA, zA, xB, zB, xC, zC); /* T = f(A,B,C) */
			add3(xT2, zT2, xT, zT, xA, zA, xB, zB); /* T2 = f(T,A,B) */
			add3(xB, zB, xB, zB, xT, zT, xA, zA); /* B = f(B,T,A) */
			t = xA;
			xA = xT2;
			xT2 = t;
			t = zA;
			zA = zT2;
			zT2 = t; /* swap A and T2 */
		}
		else if (4 * d <= 5 * e && (d - e) % 6 == 0) { /* condition 2 */
			d = (d - e) / 2;
			add3(xB, zB, xA, zA, xB, zB, xC, zC); /* B = f(A,B,C) */
			duplicate(xA, zA, xA, zA); /* A = 2*A */
		}
		else if (d <= (4 * e)) { /* condition 3 */
			d -= e;
			add3(xT, zT, xB, zB, xA, zA, xC, zC); /* T = f(B,A,C) */
			t = xB;
			xB = xT;
			xT = xC;
			xC = t;
			t = zB;
			zB = zT;
			zT = zC;
			zC = t; /* circular permutation (B,T,C) */
		}
		else if ((d + e) % 2 == 0) { /* condition 4 */
			d = (d - e) / 2;
			add3(xB, zB, xB, zB, xA, zA, xC, zC); /* B = f(B,A,C) */
			duplicate(xA, zA, xA, zA); /* A = 2*A */
		}
		else if (d % 2 == 0) { /* condition 5 */
			d /= 2;
			add3(xC, zC, xC, zC, xA, zA, xB, zB); /* C = f(C,A,B) */
			duplicate(xA, zA, xA, zA); /* A = 2*A */
		}
		else if (d % 3 == 0) { /* condition 6 */
			d = d / 3 - e;
			duplicate(xT, zT, xA, zA); /* T1 = 2*A */
			add3(xT2, zT2, xA, zA, xB, zB, xC, zC); /* T2 = f(A,B,C) */
			add3(xA, zA, xT, zT, xA, zA, xA, zA); /* A = f(T1,A,A) */
			add3(xT, zT, xT, zT, xT2, zT2, xC, zC); /* T1 = f(T1,T2,C) */
			t = xC;
			xC = xB;
			xB = xT;
			xT = t;
			t = zC;
			zC = zB;
			zB = zT;
			zT = t; /* circular permutation (C,B,T) */
		}
		else if ((d + e) % 3 == 0) { /* condition 7 */
			d = (d - 2 * e) / 3;
			add3(xT, zT, xA, zA, xB, zB, xC, zC); /* T1 = f(A,B,C) */
			add3(xB, zB, xT, zT, xA, zA, xB, zB); /* B = f(T1,A,B) */
			duplicate(xT, zT, xA, zA);
			add3(xA, zA, xA, zA, xT, zT, xA, zA); /* A = 3*A */
		}
		else if ((d - e) % 3 == 0) { /* condition 8 */
			d = (d - e) / 3;
			add3(xT, zT, xA, zA, xB, zB, xC, zC); /* T1 = f(A,B,C) */
			add3(xC, zC, xC, zC, xA, zA, xB, zB); /* C = f(A,C,B) */
			t = xB;
			xB = xT;
			xT = t;
			t = zB;
			zB = zT;
			zT = t; /* swap B and T */
			duplicate(xT, zT, xA, zA);
			add3(xA, zA, xA, zA, xT, zT, xA, zA); /* A = 3*A */
		}
		else if (e % 2 == 0) { /* condition 9 */
			e /= 2;
			add3(xC, zC, xC, zC, xB, zB, xA, zA); /* C = f(C,B,A) */
			duplicate(xB, zB, xB, zB); /* B = 2*B */
		}
	}
	add3(x, z, xA, zA, xB, zB, xC, zC);
}

/* adds Q=(x2:z2) and R=(x1:z1) and puts the result in (x3:z3),
using 5/6 mul, 6 add/sub and 6 mod. One assumes that Q-R=P or R-Q=P where P=(x:z).
Uses the following global variables:
- n : number to factor
- x, z : coordinates of P
- u, v, w : auxiliary variables
Modifies: x3, z3, u, v, w.
(x3,z3) may be identical to (x2,z2) and to (x,z)
*/
static void add3(limb *x3, limb *z3, limb *x2, limb *z2,
	limb *x1, limb *z1, limb *x, limb *z)
{
	limb *t = fieldTX;
	limb *u = fieldTZ;
	limb *v = fieldUX;
	limb *w = fieldUZ;
	SubtBigNbrModN(x2, z2, v, TestNbr, NumberLength); // v = x2-z2
	AddBigNbrModN(x1, z1, w, TestNbr, NumberLength);      // w = x1+z1
	modmult(v, w, u);       // u = (x2-z2)*(x1+z1)
	AddBigNbrModN(x2, z2, w, TestNbr, NumberLength);      // w = x2+z2
	SubtBigNbrModN(x1, z1, t, TestNbr, NumberLength); // t = x1-z1
	modmult(t, w, v);       // v = (x2+z2)*(x1-z1)
	AddBigNbrModN(u, v, t, TestNbr, NumberLength);        // t = 2*(x1*x2-z1*z2)
	modmult(t, t, w);       // w = 4*(x1*x2-z1*z2)^2
	SubtBigNbrModN(u, v, t, TestNbr, NumberLength);   // t = 2*(x2*z1-x1*z2)
	modmult(t, t, v);       // v = 4*(x2*z1-x1*z2)^2
	if (!memcmp(x, x3, NumberLength * sizeof(limb))) {
		memcpy(u, x, NumberLength * sizeof(int));
		memcpy(t, w, NumberLength * sizeof(int));
		modmult(z, t, w);
		modmult(v, u, z3);
		memcpy(x3, w, NumberLength * sizeof(int));
	}
	else {
		modmult(w, z, x3); // x3 = 4*z*(x1*x2-z1*z2)^2
		modmult(x, v, z3); // z3 = 4*x*(x2*z1-x1*z2)^2
	}
}

/* computes 2P=(x2:z2) from P=(x1:z1), with 5 mul, 4 add/sub, 5 mod.
Uses the following global variables:
- n : number to factor
- b : (a+2)/4 mod n
- u, v, w : auxiliary variables
Modifies: x2, z2, u, v, w
*/
#if 0
void print(limb *w)
{
	char pepe[1000];
	Bin2Dec(w, pepe, NumberLength, 0);
	printf("%s\n", pepe);
	memset(Xaux, 0, NumberLength * sizeof(limb));
	Xaux[0].x = 1;
	modmult(Xaux, w, Zaux);
	Bin2Dec(Zaux, pepe, NumberLength, 0);
	printf("%s\n\n", pepe);
}
#endif
static void duplicate(limb *x2, limb *z2, limb *x1, limb *z1)
{
	limb *u = fieldUZ;
	limb *v = fieldTX;
	limb *w = fieldTZ;
	AddBigNbrModN(x1, z1, w, TestNbr, NumberLength);      // w = x1+z1
	modmult(w, w, u);       // u = (x1+z1)^2
	SubtBigNbrModN(x1, z1, w, TestNbr, NumberLength); // w = x1-z1
	modmult(w, w, v);       // v = (x1-z1)^2
	modmult(u, v, x2);      // x2 = u*v = (x1^2 - z1^2)^2
	SubtBigNbrModN(u, v, w, TestNbr, NumberLength);   // w = u-v = 4*x1*z1
	modmult(fieldAA, w, u);
	AddBigNbrModN(u, v, u, TestNbr, NumberLength);        // u = (v+b*w)
	modmult(w, u, z2);      // z2 = (w*u)
}
/* End of code adapted from Paul Zimmermann's ECM4C */

/* compute gcd of value & Global var TestNbr. return 0 if either
is zero, 1 if gcd is 1 , 2 if gcd > 1 
Uses global variables Temp1, Temp2, Temp3, GD
Uses global var NumberLength indirectly */
static int gcdIsOne(limb *value)
{
	UncompressLimbsBigInteger(value, Temp1);    // Temp1 = value
	UncompressLimbsBigInteger(TestNbr, Temp2);  // Temp2 = TestNbr;
	// Return zero if value is zero or both numbers are equal.
	if (Temp1 == 0) {
		return 0;
	}
	//Temp3 = Temp1 - Temp2; // BigIntSubt(Temp1, Temp2, Temp3);
	if (Temp1 == Temp2) {
		return 0;
	}
	BigIntGcd(Temp1, Temp2, Temp3);  
	CompressLimbsBigInteger(GD, Temp3);  // GD = gcd(value, TestNbr)
	if (Temp3 < 2) {
		return (int)Temp3.lldata();    // GCD is less than 2.
	}
	return 2;      // GCD is greater than one.
}

static void GenerateSieve(int initial)
{
	int i, j, Q, initModQ;
	for (i = 0; i < 10 * SIEVE_SIZE; i += SIEVE_SIZE)
	{
		memcpy(&sieve[i], sieve2310, SIEVE_SIZE);
	}
#if MAX_PRIME_SIEVE == 11
	j = 5;
	Q = 13; /* Point to prime 13 */
#else
	j = 4;
	Q = 11; /* Point to prime 11 */
#endif
	do {
		if (initial > Q * Q) {
			initModQ = initial % Q;
			if (initModQ & 1) {    // initModQ is odd
				i = (Q - initModQ) >> 1;
			}
			else if (initModQ == 0) {
				i = 0;
			}
			else {    // initModQ is even
				i = Q - (initModQ >> 1);
			}

			for (; i < 10 * SIEVE_SIZE; i += Q) {
				sieve[i] = 1; /* Composite */
			}
		}
		else {
			i = Q * Q - initial;
			if (i < 20 * SIEVE_SIZE) {
				for (i = i / 2; i < 10 * SIEVE_SIZE; i += Q) {
					sieve[i] = 1; /* Composite */
				}
			}
			else {
				break;
			}
		}
		Q = SmallPrime[++j];
#if MAX_PRIME_SIEVE == 11
	} while (Q < 5000);
#else
} while (Q < 10 * SIEVE_SIZE);
#endif
}

static int Cos(int N) {
	switch (N % 8) 	{
	case 0:
		return 1;
	case 4:
		return -1;
	}
	return 0;   // default value 
}

static int intTotient(int N) {
	int totient, q, k;

	totient = q = N;
	if (q % 2 == 0) {
		totient /= 2;

		do {
			q /= 2;
		} while (q % 2 == 0);
	}

	if (q % 3 == 0) {
		totient = totient * 2 / 3;
		do {
			q /= 3;
		} while (q % 3 == 0);
	}

	k = 5;
	while (k * k <= q) {
		if (k % 3 != 0 && q % k == 0) {
			totient = totient * (k - 1) / k;
			do {
				q /= k;
			} while (q % k == 0);
		}
		k += 2;
	}

	if (q > 1) {
		totient = totient * (q - 1) / q;
	}
	return totient;
}

static int Moebius(int N) {
	int moebius, q, k;

	moebius = 1;
	q = N;

	if (q % 2 == 0) {
		moebius = -moebius;
		q /= 2;
		if (q % 2 == 0) {
			return 0;
		}
	}

	if (q % 3 == 0) {
		moebius = -moebius;
		q /= 3;
		if (q % 3 == 0) {
			return 0;
		}
	}
	k = 5;
	while (k * k <= q) {
		if (k % 3 != 0) {
			while (q % k == 0) {
				moebius = -moebius;
				q /= k;
				if (q % k == 0) {
					return 0;
				}
			}
		}
		k += 2;
	}

	if (q > 1) {
		moebius = -moebius;
	}
	return moebius;
}

static void GetAurifeuilleFactor(struct sFactors *pstFactors, int L, const BigInteger &BigBase)
{
	static BigInteger x, Csal, Dsal, Nbr1;  // follow advice not to use stack for these
	int k;

	BigIntPowerIntExp(BigBase, L, x);   // x <- BigBase^L.
	Csal = 1; //intToBigInteger(Csal, 1);
	Dsal = 1; // intToBigInteger(Dsal, 1);
	for (k = 1; k < DegreeAurif; k++) {
		Nbr1 = Gamma[k]; // longToBigInteger(Nbr1, Gamma[k]);
		Csal *= x;       // BigIntMultiply(Csal, x, Csal);
		Csal += Nbr1;    // BigIntAdd(Csal, Nbr1, Csal);      // Csal <- Csal * x + Gamma[k]
		Nbr1 = Delta[k]; // longToBigInteger(Nbr1, Delta[k]);
		Dsal *= x;       //BigIntMultiply(Dsal, x, Dsal);
		Dsal += Nbr1;    // BigIntAdd(Dsal, Nbr1, Dsal);      // Dsal <- Dsal * x + Gamma[k]
	}
	Nbr1 = Gamma[k];    // longToBigInteger(Nbr1, Gamma[k]);
	Csal = Csal*x;      // BigIntMultiply(Csal, x, Csal);
	Csal = Csal + Nbr1; // BigIntAdd(Csal, Nbr1, Csal);        // Csal <- Csal * x + Gamma[k]
	BigIntPowerIntExp(BigBase, (L + 1) / 2, Nbr1);   // Nbr1 <- Dsal * base^((L+1)/2)
	Nbr1 = Nbr1 * Dsal; //BigIntMultiply(Dsal, Nbr1, Nbr1);
	Dsal = Csal + Nbr1; // BigIntAdd(Csal, Nbr1, Dsal);
	if (pstFactors->multiplicity >= Max_Factors) {
		// factor list is full!
		std::string line = std::to_string(__LINE__);
		std::string mesg = "cannot factorise: too many factors: ";
		mesg += __func__;                 // + function name
		mesg += " line ";  mesg += line;  // + line number
		mesg += " in file "; mesg += __FILE__;   // + file name
		throw std::range_error(mesg);
	}
	insertBigFactor(pstFactors, Dsal);
	Dsal = Csal - Nbr1; // BigIntSubt(Csal, Nbr1, Dsal);
	if (pstFactors->multiplicity >= Max_Factors) {
		// factor list is full!
		std::string line = std::to_string(__LINE__);
		std::string mesg = "cannot factorise: too many factors: ";
		mesg += __func__;                 // + function name
		mesg += " line ";  mesg += line;  // + line number
		mesg += " in file "; mesg += __FILE__;   // + file name
		throw std::range_error(mesg);
	}
	insertBigFactor(pstFactors, Dsal);
	return;
}

// Get Aurifeuille factors.
// see https://en.wikipedia.org/wiki/Aurifeuillean_factorization
static void InsertAurifFactors(struct sFactors *pstFactors, const BigInteger &BigBase, 
	int Expon, int Incre)
{
	if (BigBase >= 386) 	{
		return;    // Base is very big, so go out.
	}
	auto Base = BigBase.lldata();

	if (Expon % 2 == 0 && Incre == -1) {
		do {
			Expon /= 2;
		} while (Expon % 2 == 0);

		Incre = Base % 4 - 2;
	}

	if (Expon % Base == 0
		&& Expon / Base % 2 != 0
		&& ((Base % 4 != 1 && Incre == 1) || (Base % 4 == 1 && Incre == -1)))
	{
		int N1, q, L, k;
		int N = (int)Base;
		if (N % 4 == 1) {
			N1 = N;
		}
		else {
			N1 = 2 * N;
		}
		DegreeAurif = intTotient(N1) / 2;
		for (k = 1; k <= DegreeAurif; k += 2) {
			AurifQ[k] = JacobiSymbol(N, k);
		}
		for (k = 2; k <= DegreeAurif; k += 2) {
			int t1 = k; // Calculate t2 = gcd(k, N1)
			int t2 = N1;
			while (t1 != 0) {
				int t3 = t2 % t1;
				t2 = t1;
				t1 = t3;
			}
			AurifQ[k] = Moebius(N1 / t2) * intTotient(t2) * Cos((N - 1) * k);
		}
		Gamma[0] = Delta[0] = 1;
		for (k = 1; k <= DegreeAurif / 2; k++) {
			int j;
			Gamma[k] = Delta[k] = 0;
			for (j = 0; j < k; j++) {
				Gamma[k] =
					Gamma[k]
					+ N * AurifQ[2 * k
					- 2 * j
					- 1] * Delta[j]
					- AurifQ[2 * k
					- 2 * j] * Gamma[j];
				Delta[k] =
					Delta[k]
					+ AurifQ[2 * k
					+ 1
					- 2 * j] * Gamma[j]
					- AurifQ[2 * k
					- 2 * j] * Delta[j];
			}
			Gamma[k] /= 2 * k;
			Delta[k] = (Delta[k] + Gamma[k]) / (2 * k + 1);
		}

		for (k = DegreeAurif / 2 + 1; k <= DegreeAurif; k++) {
			Gamma[k] = Gamma[DegreeAurif - k];
		}

		for (k = (DegreeAurif + 1) / 2; k < DegreeAurif; k++) {
			Delta[k] = Delta[DegreeAurif - k - 1];
		}

		q = Expon / (int)Base;
		L = 1;
		while (L * L <= q) {
			if (q % L == 0) {
				GetAurifeuilleFactor(pstFactors, L, BigBase);
				if (q != L * L) {
					GetAurifeuilleFactor(pstFactors, q / L, BigBase);
				}
			}
			L += 2;
		}
	}
	return;
}


static void Cunningham(struct sFactors *pstFactors, const BigInteger &BigBase, int Expon,
	int increment, const BigInteger &BigOriginal)
{
//#ifdef __EMSCRIPTEN__
//	char url[200];
//	char *ptrUrl;
//#endif
	int Expon2, k;
	//char *ptrFactorsAscii;
	static BigInteger Nbr1, Nbr2;      // follow advice not to use stack for these

	//factorsAscii[0] = 0;    // Indicate no new factor found in advance.
	Expon2 = Expon;
	//if (cunningham && BigOriginal->nbrLimbs > 4)
	//{   // Enter here on numbers of more than 40 digits if the user selected
	//	 //get Cunningham factors from server.
	//	#ifdef __EMSCRIPTEN__
	//			printf ("%s\n", (lang ? "4\nObteniendo los factores primitivos conocidos del servidor Web.\n" :
	//				"4\nRequesting known primitive factors from Web server.\n");
	//			// Format URL.
	//			ptrUrl = url;
	//			strcpy(ptrUrl, "factors.pl?base=");
	//			ptrUrl += strlen(ptrUrl);
	//			int2dec(&ptrUrl, BigBase->limbs[0].x);
	//			strcpy(ptrUrl, "&expon=");
	//			ptrUrl += strlen(ptrUrl);
	//			int2dec(&ptrUrl, Expon);
	//			strcpy(ptrUrl, "&type=");
	//			ptrUrl += strlen(ptrUrl);
	//			*ptrUrl++ = (increment > 0 ? 'p' : 'm');
	//			*ptrUrl = 0;
	//			getCunn(url, factorsAscii);
	//	#endif
	//}
	//ptrFactorsAscii = factorsAscii;
	//while (*ptrFactorsAscii > ' ')
	//{ // Loop through factors found in server.
	//	int nbrDigits;
	//	char *ptrEndFactor = findChar(ptrFactorsAscii, '*');
	//	if (ptrEndFactor == NULL)
	//	{
	//		nbrDigits = (int)strlen(ptrFactorsAscii);
	//		ptrEndFactor = ptrFactorsAscii + nbrDigits;
	//	}
	//	else
	//	{
	//		nbrDigits = (int)(ptrEndFactor - ptrFactorsAscii);
	//		ptrEndFactor++;
	//	}
	//	Dec2Bin(ptrFactorsAscii, Nbr1.limbs, nbrDigits, &Nbr1.nbrLimbs);
	//	Nbr1.sign = SIGN_POSITIVE;
	//	ptrFactorsAscii = ptrEndFactor;
	//	insertBigFactor(pstFactors, &Nbr1);
	//}
	while (Expon2 % 2 == 0 && increment == -1) {
		Expon2 /= 2;
		BigIntPowerIntExp(BigBase, Expon2, Nbr1);
		Nbr1 += increment; //addbigint(Nbr1, increment);
		if (pstFactors->multiplicity >= Max_Factors) {
			// factor list is full!
			std::string line = std::to_string(__LINE__);
			std::string mesg = "cannot factorise: too many factors: ";
			mesg += __func__;                 // + function name
			mesg += " line ";  mesg += line;  // + line number
			mesg += " in file "; mesg += __FILE__;   // + file name
			throw std::range_error(mesg);
		}
		insertBigFactor(pstFactors, Nbr1);
		InsertAurifFactors(pstFactors, BigBase, Expon2, 1);
	}

	k = 1;
	while (k * k <= Expon) {
		if (Expon % k == 0) {
			if (k % 2 != 0) { /* Only for odd exponent */
				BigIntPowerIntExp(BigBase, Expon / k, Nbr1);
				Nbr1 += increment; //addbigint(Nbr1, increment);
				BigIntGcd(Nbr1, BigOriginal, Nbr2);   // Nbr2 <- gcd(Base^(Expon/k)+incre, original)
				if (pstFactors->multiplicity >= Max_Factors) {
					// factor list is full!
					std::string line = std::to_string(__LINE__);
					std::string mesg = "cannot factorise: too many factors: ";
					mesg += __func__;                 // + function name
					mesg += " line ";  mesg += line;  // + line number
					mesg += " in file "; mesg += __FILE__;   // + file name
					throw std::range_error(mesg);
				}
				insertBigFactor(pstFactors, Nbr2);
				Temp1 = BigOriginal; // CopyBigInt(Temp1, *BigOriginal);
				Nbr1 = Temp1 / Nbr2; // BigIntDivide(Temp1, Nbr2, Nbr1);
				if (pstFactors->multiplicity >= Max_Factors) {
					// factor list is full!
					std::string line = std::to_string(__LINE__);
					std::string mesg = "cannot factorise: too many factors: ";
					mesg += __func__;                 // + function name
					mesg += " line ";  mesg += line;  // + line number
					mesg += " in file "; mesg += __FILE__;   // + file name
					throw std::range_error(mesg);
				}
				insertBigFactor(pstFactors, Nbr1);
				InsertAurifFactors(pstFactors, BigBase, Expon / k, increment);
			}

			if ((Expon / k) % 2 != 0) { /* Only for odd exponent */
				BigIntPowerIntExp(BigBase, k, Nbr1);
				Nbr1 += increment; //addbigint(Nbr1, increment);
				BigIntGcd(Nbr1, BigOriginal, Nbr2);   // Nbr2 <- gcd(Base^k+incre, original)
				if (pstFactors->multiplicity >= Max_Factors) {
					// factor list is full!
					std::string line = std::to_string(__LINE__);
					std::string mesg = "cannot factorise: too many factors: ";
					mesg += __func__;                 // + function name
					mesg += " line ";  mesg += line;  // + line number
					mesg += " in file "; mesg += __FILE__;   // + file name
					throw std::range_error(mesg);
				}
				insertBigFactor(pstFactors, Nbr2);
				Temp1 = BigOriginal; // CopyBigInt(Temp1, *BigOriginal);
				Nbr1 = Temp1 / Nbr2; // BigIntDivide(Temp1, Nbr2, Nbr1);
				if (pstFactors->multiplicity >= Max_Factors) {
					// factor list is full!
					std::string line = std::to_string(__LINE__);
					std::string mesg = "cannot factorise: too many factors: ";
					mesg += __func__;                 // + function name
					mesg += " line ";  mesg += line;  // + line number
					mesg += " in file "; mesg += __FILE__;   // + file name
					throw std::range_error(mesg);
				}
				insertBigFactor(pstFactors, Nbr1);
				InsertAurifFactors(pstFactors, BigBase, k, increment);
			}
		}
		k++;
	}
	return;
}

static bool ProcessExponent(struct sFactors *pstFactors, 
	const BigInteger &nbrToFactor, int Exponent)
{
#ifdef __EMSCRIPTEN__
	char status[200] = { 0 };
	char *ptrStatus;
#endif
	static BigInteger NFp1, NFm1, nthRoot, rootN1, rootN, rootbak;   // follow advice not to use stack for these
	static BigInteger nextroot, dif;     // follow advice not to use stack for these
	double log2N;

#ifdef __EMSCRIPTEN__
	int elapsedTime = (int)(tenths() - originalTenthSecond);

	if (elapsedTime / 10 != oldTimeElapsed / 10) 
	{
		oldTimeElapsed = elapsedTime;
		ptrStatus = status;
		strcpy(ptrStatus, lang ? "Transcurrió " : "Time elapsed: ");
		ptrStatus += strlen(ptrStatus);
		GetDHMS(&ptrStatus, elapsedTime / 10);
		strcpy(ptrStatus, lang ? "Exponente potencia +/- 1: " :
			"Power +/- 1 exponent: ");
		ptrStatus += strlen(ptrStatus);
		int2dec(&ptrStatus, Exponent);
		printf ("%s\n", status);
	}
#endif

	NFp1 = nbrToFactor; 
	NFp1 ++;              // NFp1 <- NumberToFactor + 1
	NFm1 = nbrToFactor;   // CopyBigInt(NFm1, *nbrToFactor);
	NFm1 --;              // NFm1 <- NumberToFactor - 1
	log2N = logBigNbr(NFp1) / Exponent;    // Find nth root of number to factor.
	expBigInt(nthRoot, log2N);
	rootbak = nthRoot;

	for (;;) {
		BigIntPowerIntExp(nthRoot, Exponent - 1, rootN1); // rootN1 <- nthRoot ^ (Exponent-1)
		rootN = nthRoot*rootN1; //BigIntMultiply(nthRoot, rootN1, rootN);     // rootN <- nthRoot ^ Exponent
		dif = NFp1 - rootN; // BigIntSubt(NFp1, rootN, dif);    
		if (dif == 0) { // Perfect power-1
			Cunningham(pstFactors, nthRoot, Exponent, -1, nbrToFactor);
			return true;
		}
		dif ++;        // dif <- dif + 1
		Temp1 = dif / rootN1;        
		subtractdivide(Temp1, 0, Exponent);        // Temp1 <- Temp1 / Exponent
		nextroot = Temp1 + nthRoot;         // 
		nextroot --;  //   nextroot = nextroot - 1             
		nthRoot = nextroot - nthRoot;     
		if (nthRoot >= 0) {
			break; // Not a perfect power
		}
		nthRoot = nextroot; 
	}

	nthRoot = rootbak;
	for (;;) {
		BigIntPowerIntExp(nthRoot, Exponent - 1, rootN1); // rootN1 <- nthRoot ^ (Exponent-1)
		rootN = nthRoot*rootN1;      // rootN <- nthRoot ^ Exponent
		dif = NFm1 - rootN;             
		if (dif == 0)
		{ // Perfect power+1
			Cunningham(pstFactors, nthRoot, Exponent, 1, nbrToFactor);
			return true;
		}
		dif ++;                        // dif <- dif + 1
		Temp1 = dif / rootN1;   
		subtractdivide(Temp1, 0, Exponent);   // Temp1 <- Temp1 / Exponent
		nextroot = Temp1 + nthRoot;       
		nextroot --;             // nextroot <- nextroot - 1
		nthRoot = nextroot - nthRoot;   
		if (nthRoot >= 0) {
			break;               // Not a perfect power
		}
		nthRoot = nextroot;      // CopyBigInt(nthRoot, nextroot);
	}
	return false;
}

/* check whether the number +/- 1 is a perfect power*/
static void PowerPM1Check(struct sFactors *pstFactors, const BigInteger &nbrToFactor)
{
	bool plus1 = false;
	bool minus1 = false;
	int Exponent = 0;
	int i, j;
	int modulus;
	int mod9 = nbrToFactor % 9; // getRemainder(*nbrToFactor, 9);
	int maxExpon = nbrToFactor.bitLength(); 
	int numPrimes = 2 * maxExpon + 3;
	/* if numPrimes >= 66472 i.e. maxExpon > 33233 i.e. nbrLimbs > 1072*/
	if ((numPrimes >> 3) >= (2 * 33231 + 3 + 7) / 8) {
		std::string line = std::to_string(__LINE__);
		std::string mesg = "cannot factorise: number exceeds 10,000 digits "; 
		mesg += __func__;
		mesg += " line ";  mesg += line;
		mesg += " in file "; mesg += __FILE__;
		throw std::range_error(mesg);
	}
	double logar = logBigNbr(nbrToFactor);
	// 33219 = logarithm base 2 of max number supported = 10^10,000.
	unsigned char ProcessExpon[(33231 + 7) / 8];
	unsigned char primes[(2 * 33231 + 3 + 7) / 8];

	memset(ProcessExpon, 0xFF, sizeof(ProcessExpon));
	memset(primes, 0xFF, sizeof(primes));

	for (i = 2; i * i < numPrimes; i++) { // Generation of primes using sieve of Eratosthenes.
		if (primes[i >> 3] & (1 << (i & 7))) {     // Number i is prime.
			for (j = i * i; j < numPrimes; j += i)
			{   // Mark multiple of i as composite.
				primes[j >> 3] &= ~(1 << (j & 7));
			}
		}
	}

	// If the number +/- 1 is multiple of a prime but not a multiple
	// of its square then the number +/- 1 cannot be a perfect power.
	for (i = 2; i < numPrimes; i++) {
		if (primes[i >> 3] & (1 << (i & 7))) {      // i is prime according to sieve.
			unsigned long long remainder;
			int index;
			Temp1 = (long long)i*(long long)i; 
			Temp2 = nbrToFactor % Temp1;     // Temp2 <- nbrToFactor % (i*i)
			remainder = Temp2.lldata();
			if (remainder % i == 1 && remainder != 1) {
				plus1 = true; // NumberFactor cannot be a power + 1
			}
			if (remainder % i == (unsigned int)i - 1 &&
				remainder != (unsigned int)i*(unsigned int)i - 1) 
			{
				minus1 = true; // NumberFactor cannot be a power - 1
			}
			if (minus1 && plus1) {
				return;
			}
			index = i / 2;
			if (!(ProcessExpon[index >> 3] & (1 << (index & 7)))) {
				continue;
			}
			modulus = remainder % i;
			if (modulus > (plus1 ? 1 : 2) && modulus < (minus1 ? i - 1 : i - 2)) {
				for (j = i / 2; j <= maxExpon; j += i / 2) {
					ProcessExpon[j >> 3] &= ~(1 << (j & 7));
				}
			}
			else {
				if (modulus == i - 2) {
					for (j = i - 1; j <= maxExpon; j += i - 1) {
						ProcessExpon[j >> 3] &= ~(1 << (j & 7));
					}
				}
			}
		}
	}

	for (j = 2; j < 100; j++) {
		double u = logar / log(j) + .000005;
		Exponent = (int)floor(u);
		if (u - Exponent > .00001)
			continue;
		if (Exponent % 3 == 0 && mod9 > 2 && mod9 < 7)
			continue;
		if (!(ProcessExpon[Exponent >> 3] & (1 << (Exponent & 7))))
			continue;
		if (ProcessExponent(pstFactors, nbrToFactor, Exponent))
			return;  // number is a perfect power +/- 1
	}

	for (; Exponent >= 2; Exponent--) {
		if (Exponent % 3 == 0 && mod9 > 2 && mod9 < 7) {
			continue;
		}
		if (!(ProcessExpon[Exponent >> 3] & (1 << (Exponent & 7)))) {
			continue;
		}
		if (ProcessExponent(pstFactors, nbrToFactor, Exponent)) {
			return;   // number is  a perfect power
		}
	}
	return;
}

// Perform Lehman algorithm
static void Lehman(const BigInteger &nbr, int k, BigInteger &factor)
{
	unsigned int bitsSqrLow[] =  // could combine low and high into 64 bit integers
	{
		0x00000003, // 3
		0x00000013, // 5
		0x00000017, // 7
		0x0000023B, // 11
		0x0000161B, // 13
		0x0001A317, // 17
		0x00030AF3, // 19
		0x0005335F, // 23
		0x13D122F3, // 29
		0x121D47B7, // 31
		0x5E211E9B, // 37
		0x82B50737, // 41   N.B. exceeds max for signed int
		0x83A3EE53, // 43   N.B. exceeds max for signed int
		0x1B2753DF, // 47
		0x3303AED3, // 53
		0x3E7B92BB, // 59
		0x0A59F23B, // 61
	};
	unsigned int bitsSqrHigh[] =
	{
		0x00000000, // 3
		0x00000000, // 5
		0x00000000, // 7
		0x00000000, // 11
		0x00000000, // 13
		0x00000000, // 17
		0x00000000, // 19
		0x00000000, // 23
		0x00000000, // 29
		0x00000000, // 31
		0x00000016, // 37
		0x000001B3, // 41
		0x00000358, // 43
		0x00000435, // 47
		0x0012DD70, // 53
		0x022B6218, // 59
		0x1713E694, // 61
	};
	int primes[] = { 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61 };
	int nbrs[17];
	int diffs[17];
	int i, j, m, r;
	static BigInteger sqrRoot, nextroot;   // follow advice not to use stack for these
	static BigInteger a, c, sqr, val;      // follow advice not to use stack for these

	if (nbr.isEven()) { // nbr Even
		r = 0;
		m = 1;
	}
	else {
		if (k % 2 == 0) { // k Even
			r = 1;
			m = 2;
		}
		else { // k Odd
			r = (k + nbr.lldata()) & 3;   // was limbs[0].x
			m = 4;
		}
	}

	sqr = k << 2; 
	sqr *= nbr; 
	//squareRoot(sqr.limbs, sqrRoot.limbs, sqr.nbrLimbs, &sqrRoot.nbrLimbs);
	a = sqr.sqRoot();

	for (;;) {
		if ((a.lldata() & (m - 1)) == r) {
			nextroot = a*a - sqr; 
			if (nextroot >= 0) {
				break;
			}
		}
		a ++;              // a <- a + 1
	}
	nextroot = a*a; 
	c = nextroot - sqr; 

	for (i = 0; i < 17; i++) {
		//int pr = pr;
		nbrs[i] = c % primes[i]; // getRemainder(c, pr);    
		diffs[i] = m * ((a% primes[i]) * 2 + m) % primes[i];
	}

	for (j = 0; j < 10000; j++) {
		for (i = 0; i < 17; i++) {
			int shiftBits = nbrs[i];
			if (shiftBits < 32) {
				if ((bitsSqrLow[i] & (1 << shiftBits)) == 0) { 
					break;  // Not a perfect square
				}
			}
			else if ((bitsSqrHigh[i] & (1 << (shiftBits - 32))) == 0) { 
				break;  // Not a perfect square
			}
		}
		if (i == 17) { // Test for perfect square
			c = m * j;           
			val = a + c; 
			c = val*val - sqr;     
   			sqrRoot = c.sqRoot();   // sqrRoot <- sqrt(c)
			sqrRoot +=  val; 
			BigIntGcd(sqrRoot, nbr, c);         // Get GCD(sqrRoot + val, nbr)
			if (c >= LIMB_RANGE) {    // Non-trivial factor has been found.
				factor = c;     // CopyBigInt(*factor, c);
				return;
			}
		}

		for (i = 0; i < 17; i++) {
			nbrs[i] = (nbrs[i] + diffs[i]) % primes[i];
			diffs[i] = (diffs[i] + 2 * m * m) % primes[i];
		}
	}
	factor = 1;   // Factor not found.
	return;
}

static enum eEcmResult ecmCurve(BigInteger &N)
{
	BigInteger potentialFactor;
#ifdef __EMSCRIPTEN__
	//char text[20];
#endif
	EC %= 50000000;   // Convert to curve number.
	for (;;) {
#ifdef __EMSCRIPTEN__
		char *ptrText;
#endif
		int I, Pass;
		int i, j, u;
		long long L1, L2, LS, P, IP, Paux = 1;

		EC++;   // increment curve number

//#ifdef __EMSCRIPTEN__
//			text[0] = '7';
			//ptrText = &text[1];
//			int2dec(&ptrText, EC);
//			printf ("%s\n", text);
//#endif
		L1 = NumberLength * 9;        // Get number of digits.
		if (L1 > 30 && L1 <= 90)    // If between 30 and 90 digits...
		{                             // Switch to SIQS.
			int limit = limits[((int)L1 - 31) / 5];
			if (EC % 50000000 >= limit)
			{                           // Switch to SIQS.
				EC += TYP_SIQS;
				return CHANGE_TO_SIQS;
			}
		}

		// Try to factor BigInteger N using Lehman algorithm. Result in potentialFactor.
		Lehman(N, EC % 50000000, potentialFactor);

		if (potentialFactor.lldata() >= LIMB_RANGE) {                // large Factor found.
			memcpy(GD, potentialFactor.limbs, NumberLength * sizeof(limb));
			memset(&GD[potentialFactor.nbrLimbs], 0,
				(NumberLength - potentialFactor.nbrLimbs) * sizeof(limb));
			foundByLehman = true;
			return FACTOR_FOUND;
		}

		L1 = 2000;
		L2 = 200000;
		LS = 45;
		Paux = EC;
		nbrPrimes = 303; /* Number of primes less than 2000 */
		if (EC > 25) {
			if (EC < 326) {
				L1 = 50000;
				L2 = 5000000;
				LS = 224;
				Paux = EC - 24;
				nbrPrimes = 5133; /* Number of primes less than 50000 */
			}
			else {
				if (EC < 2000) {
					L1 = 1000000;
					L2 = 100000000;
					LS = 1001;
					Paux = EC - 299;
					nbrPrimes = 78498; /* Number of primes less than 1000000 */
				}
				else {
					L1 = 11000000;
					L2 = 1100000000;
					LS = 3316;
					Paux = EC - 1900;
					nbrPrimes = 726517; /* Number of primes less than 11000000 */
				}
			}
		}
//#ifdef __EMSCRIPTEN__
		ptrText = ptrLowerText;  // Point after number that is being factored.
		auto elapsedTime = (int)(tenths() - originalTenthSecond);
		GetDHMSt(&ptrText, elapsedTime);
		strcpy(ptrText, lang ? " ECM Curva " : " ECM Curve ");
		ptrText += strlen(ptrText);
		int2dec(&ptrText, EC);   // Show curve number.
		strcpy(ptrText, lang ? " usando límites B1=" : " using bounds B1=");
		ptrText += strlen(ptrText);
		int2dec(&ptrText, L1);   // Show first bound.
		strcpy(ptrText, lang ? " y B2=" : " and B2=");
		ptrText += strlen(ptrText);
		int2dec(&ptrText, L2);   // Show second bound.
		strcpy(ptrText, "\n");
		ptrText += strlen(ptrText);
		printf ("%s", ptrLowerText);
#if 0
		primalityString =
			textAreaContents
			+ StringToLabel
			+ "\nLimit (B1="
			+ L1
			+ "; B2="
			+ L2
			+ ")    Curve ";
		UpperLine = "Digits in factor:   ";
		LowerLine = "Probability:        ";
		for (I = 0; I < 6; I++)
		{
			UpperLine += "    >= " + (I * 5 + 15);
			Prob =
				(int)Math.round(
					100
					* (1 - exp(-((double)L1 * (double)Paux) / ProbArray[I])));
			if (Prob == 100)
			{
				LowerLine += "    100% ";
			}
			else
			{
				if (Prob >= 10)
				{
					LowerLine += "     " + Prob + "% ";
				}
				else
				{
					LowerLine += "      " + Prob + "% ";
				}
			}
		} /* end for */
		lowerTextArea.setText(
			primalityString + EC + "\n" + UpperLine + "\n" + LowerLine);
#endif
		//#endif

		//  Compute A0 <- 2 * (EC+1)*modinv(3 * (EC+1) ^ 2 - 1, N) mod N
		// Aux2 <- 1 in Montgomery notation.
		memcpy(Aux2, MontgomeryMultR1, NumberLength * sizeof(limb));
		modmultInt(Aux2, EC + 1, Aux2);            // Aux2 <- EC + 1.
		modmultInt(Aux2, 2, Aux1);                 // Aux1 <- 2*(EC+1)
		modmultInt(Aux2, EC + 1, Aux3);            // Aux3 <- (EC + 1)^2
		modmultInt(Aux3, 3, Aux3);                 // Aux3 <- 3*(EC + 1)^2
												   // Aux2 <- 3*(EC + 1)^2 - 1 
		SubtBigNbrModN(Aux3, MontgomeryMultR1, Aux2, TestNbr, NumberLength);
		ModInvBigNbr(Aux2, Aux2, TestNbr, NumberLength);
		modmult(Aux1, Aux2, A0);                   // A0 <- 2*(EC+1)/(3*(EC+1)^2 - 1)

												   //  if A0*(A0 ^ 2 - 1)*(9 * A0 ^ 2 - 1) mod N=0 then select another curve.
		modmult(A0, A0, A02);          // A02 <- A0^2
		modmult(A02, A0, A03);         // A03 <- A0^3
		SubtBigNbrModN(A03, A0, Aux1, TestNbr, NumberLength);  // Aux1 <- A0^3 - A0
		modmultInt(A02, 9, Aux2);      // Aux2 <- 9*A0^2
		SubtBigNbrModN(Aux2, MontgomeryMultR1, Aux2, TestNbr, NumberLength); // Aux2 <- 9*A0^2-1
		modmult(Aux1, Aux2, Aux3);
		if (BigNbrIsZero(Aux3, NumberLength)) {
			continue;
		}
		//   Z <- 4 * A0 mod N
		modmultInt(A0, 4, Z);
		//   A = (-3 * A0 ^ 4 - 6 * A0 ^ 2 + 1)*modinv(4 * A0 ^ 3, N) mod N
		modmultInt(A02, 6, Aux1);      // Aux1 <- 6*A0^2
		SubtBigNbrModN(MontgomeryMultR1, Aux1, Aux1, TestNbr, NumberLength);
		modmult(A02, A02, Aux2);       // Aux2 <- A0^4
		modmultInt(Aux2, 3, Aux2);     // Aux2 <- 3*A0^4
		SubtBigNbrModN(Aux1, Aux2, Aux1, TestNbr, NumberLength);
		modmultInt(A03, 4, Aux2);      // Aux2 <- 4*A0^3
		ModInvBigNbr(Aux2, Aux3, TestNbr, NumberLength);
		modmult(Aux1, Aux3, A0);
		//   AA <- (A + 2)*modinv(4, N) mod N
		modmultInt(MontgomeryMultR1, 2, Aux2);  // Aux2 <- 2
		AddBigNbrModN(A0, Aux2, Aux1, TestNbr, NumberLength); // Aux1 <- A0+2
		modmultInt(MontgomeryMultR1, 4, Aux2);  // Aux2 <- 4
		ModInvBigNbr(Aux2, Aux2, TestNbr, NumberLength);
		modmult(Aux1, Aux2, AA);
		//   X <- (3 * A0 ^ 2 + 1) mod N
		modmultInt(A02, 3, Aux1);    // Aux1 <- 3*A0^2
		AddBigNbrModN(Aux1, MontgomeryMultR1, X, TestNbr, NumberLength);
		/**************/
		/* First step */
		/**************/
		memcpy(Xaux, X, NumberLength * sizeof(limb));
		memcpy(Zaux, Z, NumberLength * sizeof(limb));
		memcpy(GcdAccumulated, MontgomeryMultR1, (NumberLength + 1) * sizeof(limb));
		for (Pass = 0; Pass < 2; Pass++) {
			/* For powers of 2 */
			indexPrimes = 0;
			StepECM = 1;
			for (I = 1; I <= L1; I <<= 1) {
				duplicate(X, Z, X, Z);
			}
			for (I = 3; I <= L1; I *= 3) {
				duplicate(W1, W2, X, Z);
				add3(X, Z, X, Z, W1, W2, X, Z);
			}

			if (Pass == 0) {
				modmult(GcdAccumulated, Z, Aux1);
				memcpy(GcdAccumulated, Aux1, NumberLength * sizeof(limb));
			}
			else {
				if (gcdIsOne(Z) > 1) {
					return FACTOR_FOUND;
				}
			}

			/* for powers of odd primes */
			indexM = 1;
			do {
				indexPrimes++;
				P = SmallPrime[indexM];
				for (IP = P; IP <= L1; IP *= P) {
					prac((int)P, X, Z, W1, W2, W3, W4);
				}
				indexM++;
				if (Pass == 0) {
					modmult(GcdAccumulated, Z, Aux1);
					memcpy(GcdAccumulated, Aux1, NumberLength * sizeof(limb));
				}
				else {
					if (gcdIsOne(Z) > 1) { 
						return FACTOR_FOUND;
					}
				}
			} while (SmallPrime[indexM - 1] <= LS);

			P += 2;

			/* Initialize sieve2310[n]: 1 if gcd(P+2n,2310) > 1, 0 otherwise */
			u = (int)P;
			for (i = 0; i < SIEVE_SIZE; i++) {
				sieve2310[i] =
					(u % 3 == 0
						|| u % 5 == 0
						|| u % 7 == 0
#if MAX_PRIME_SIEVE == 11
						|| u % 11 == 0
#endif
						? (unsigned char)1 : (unsigned char)0);
				u += 2;
			}

			do {
				/* Generate sieve */
				GenerateSieve((int)P);

				/* Walk through sieve */

				for (i = 0; i < 10 * SIEVE_SIZE; i++) {
					if (sieve[i] != 0) {
						continue; /* Do not process composites */
					}
					if (P + 2 * i > L1) {
						break;
					}
					indexPrimes++;
					prac((int)(P + 2 * i), X, Z, W1, W2, W3, W4);
					if (Pass == 0) {
						modmult(GcdAccumulated, Z, Aux1);
						memcpy(GcdAccumulated, Aux1, NumberLength * sizeof(limb));
					}
					else {
						if (gcdIsOne(Z) > 1) {
							return FACTOR_FOUND;
						}
					}
				}
				P += 20 * SIEVE_SIZE;
			} while (P < L1);

			if (Pass == 0) {
				if (BigNbrIsZero(GcdAccumulated, NumberLength))
				{ // If GcdAccumulated is
					memcpy(X, Xaux, NumberLength * sizeof(limb));
					memcpy(Z, Zaux, NumberLength * sizeof(limb));
					continue; // multiple of TestNbr, continue.
				}
				if (gcdIsOne(GcdAccumulated) > 1) {
					return FACTOR_FOUND;
				}
				break;
			}
		} /* end for Pass */

		  /******************************************************/
		  /* Second step (using improved standard continuation) */
		  /******************************************************/
		StepECM = 2;
		j = 0;
		for (u = 1; u < SIEVE_SIZE; u += 2) {
			if (u % 3 == 0 || u % 5 == 0 || u % 7 == 0
#if MAX_PRIME_SIEVE == 11
				|| u % 11 == 0
#endif
				)
			{
				sieve2310[u / 2] = (unsigned char)1;
			}
			else {
				sieve2310[(sieveidx[j++] = u / 2)] = (unsigned char)0;
			}
		}
		memcpy(&sieve2310[HALF_SIEVE_SIZE], &sieve2310[0], HALF_SIEVE_SIZE);
		memcpy(Xaux, X, NumberLength * sizeof(limb));  // (X:Z) -> Q (output
		memcpy(Zaux, Z, NumberLength * sizeof(limb));  //         from step 1)
		for (Pass = 0; Pass < 2; Pass++) {
			int Qaux, J;
			memcpy(GcdAccumulated, MontgomeryMultR1, NumberLength * sizeof(limb));
			memcpy(UX, X, NumberLength * sizeof(limb));
			memcpy(UZ, Z, NumberLength * sizeof(limb));  // (UX:UZ) -> Q 
			ModInvBigNbr(Z, Aux1, TestNbr, NumberLength);
			modmult(Aux1, X, root[0]); // root[0] <- X/Z (Q)
			J = 0;
			AddBigNbrModN(X, Z, Aux1, TestNbr, NumberLength);
			modmult(Aux1, Aux1, W1);
			SubtBigNbrModN(X, Z, Aux1, TestNbr, NumberLength);
			modmult(Aux1, Aux1, W2);
			modmult(W1, W2, TX);
			SubtBigNbrModN(W1, W2, Aux1, TestNbr, NumberLength);
			modmult(Aux1, AA, Aux2);
			AddBigNbrModN(Aux2, W2, Aux3, TestNbr, NumberLength);
			modmult(Aux1, Aux3, TZ); // (TX:TZ) -> 2Q
			SubtBigNbrModN(X, Z, Aux1, TestNbr, NumberLength);
			AddBigNbrModN(TX, TZ, Aux2, TestNbr, NumberLength);
			modmult(Aux1, Aux2, W1);
			AddBigNbrModN(X, Z, Aux1, TestNbr, NumberLength);
			SubtBigNbrModN(TX, TZ, Aux2, TestNbr, NumberLength);
			modmult(Aux1, Aux2, W2);
			AddBigNbrModN(W1, W2, Aux1, TestNbr, NumberLength);
			modmult(Aux1, Aux1, Aux2);
			modmult(Aux2, UZ, X);
			SubtBigNbrModN(W1, W2, Aux1, TestNbr, NumberLength);
			modmult(Aux1, Aux1, Aux2);
			modmult(Aux2, UX, Z); // (X:Z) -> 3Q
			for (I = 5; I < SIEVE_SIZE; I += 2) {
				memcpy(WX, X, NumberLength * sizeof(limb));
				memcpy(WZ, Z, NumberLength * sizeof(limb));
				SubtBigNbrModN(X, Z, Aux1, TestNbr, NumberLength);
				AddBigNbrModN(TX, TZ, Aux2, TestNbr, NumberLength);
				modmult(Aux1, Aux2, W1);
				AddBigNbrModN(X, Z, Aux1, TestNbr, NumberLength);
				SubtBigNbrModN(TX, TZ, Aux2, TestNbr, NumberLength);
				modmult(Aux1, Aux2, W2);
				AddBigNbrModN(W1, W2, Aux1, TestNbr, NumberLength);
				modmult(Aux1, Aux1, Aux2);
				modmult(Aux2, UZ, X);
				SubtBigNbrModN(W1, W2, Aux1, TestNbr, NumberLength);
				modmult(Aux1, Aux1, Aux2);
				modmult(Aux2, UX, Z); // (X:Z) -> 5Q, 7Q, ...
				if (Pass == 0) {
					modmult(GcdAccumulated, Aux1, Aux2);
					memcpy(GcdAccumulated, Aux2, NumberLength * sizeof(limb));
				}
				else {
					if (gcdIsOne(Aux1) > 1) {
						return FACTOR_FOUND;
					}
				}
				if (I == HALF_SIEVE_SIZE) {
					memcpy(DX, X, NumberLength * sizeof(limb));
					memcpy(DZ, Z, NumberLength * sizeof(limb));  // (DX:DZ) -> HALF_SIEVE_SIZE*Q
				}
				if (I % 3 != 0 && I % 5 != 0 && I % 7 != 0
#if MAX_PRIME_SIEVE == 11
					&& I % 11 != 0
#endif
					)
				{
					J++;
					ModInvBigNbr(Z, Aux1, TestNbr, NumberLength);
					modmult(Aux1, X, root[J]); // root[J] <- X/Z
				}
				memcpy(UX, WX, NumberLength * sizeof(limb));  // (UX:UZ) <-
				memcpy(UZ, WZ, NumberLength * sizeof(limb));  // Previous (X:Z)
			} /* end for I */

			AddBigNbrModN(DX, DZ, Aux1, TestNbr, NumberLength);
			modmult(Aux1, Aux1, W1);
			SubtBigNbrModN(DX, DZ, Aux1, TestNbr, NumberLength);
			modmult(Aux1, Aux1, W2);
			modmult(W1, W2, X);
			SubtBigNbrModN(W1, W2, Aux1, TestNbr, NumberLength);
			modmult(Aux1, AA, Aux2);
			AddBigNbrModN(Aux2, W2, Aux3, TestNbr, NumberLength);
			modmult(Aux1, Aux3, Z);
			memcpy(UX, X, NumberLength * sizeof(limb));
			memcpy(UZ, Z, NumberLength * sizeof(limb));    // (UX:UZ) -> SIEVE_SIZE*Q
			AddBigNbrModN(X, Z, Aux1, TestNbr, NumberLength);
			modmult(Aux1, Aux1, W1);
			SubtBigNbrModN(X, Z, Aux1, TestNbr, NumberLength);
			modmult(Aux1, Aux1, W2);
			modmult(W1, W2, TX);
			SubtBigNbrModN(W1, W2, Aux1, TestNbr, NumberLength);
			modmult(Aux1, AA, Aux2);
			AddBigNbrModN(Aux2, W2, Aux3, TestNbr, NumberLength);
			modmult(Aux1, Aux3, TZ); // (TX:TZ) -> 2*SIEVE_SIZE*Q
			SubtBigNbrModN(X, Z, Aux1, TestNbr, NumberLength);
			AddBigNbrModN(TX, TZ, Aux2, TestNbr, NumberLength);
			modmult(Aux1, Aux2, W1);
			AddBigNbrModN(X, Z, Aux1, TestNbr, NumberLength);
			SubtBigNbrModN(TX, TZ, Aux2, TestNbr, NumberLength);
			modmult(Aux1, Aux2, W2);
			AddBigNbrModN(W1, W2, Aux1, TestNbr, NumberLength);
			modmult(Aux1, Aux1, Aux2);
			modmult(Aux2, UZ, X);
			SubtBigNbrModN(W1, W2, Aux1, TestNbr, NumberLength);
			modmult(Aux1, Aux1, Aux2);
			modmult(Aux2, UX, Z); // (X:Z) -> 3*SIEVE_SIZE*Q
			Qaux = (int)(L1 / (2 * SIEVE_SIZE));
			maxIndexM = (int)(L2 / (2 * SIEVE_SIZE));
			for (indexM = 0; indexM <= maxIndexM; indexM++) {
				if (indexM >= Qaux) { // If inside step 2 range... 
					if (indexM == 0) {
						ModInvBigNbr(UZ, Aux3, TestNbr, NumberLength);
						modmult(UX, Aux3, Aux1); // Aux1 <- X/Z (SIEVE_SIZE*Q)
					}
					else {
						ModInvBigNbr(Z, Aux3, TestNbr, NumberLength);
						modmult(X, Aux3, Aux1); // Aux1 <- X/Z (3,5,*
					}                         //         SIEVE_SIZE*Q)

											  /* Generate sieve */
					if (indexM % 10 == 0 || indexM == Qaux) {
						GenerateSieve(indexM / 10 * (20 * SIEVE_SIZE) + 1);
					}
					/* Walk through sieve */
					J = HALF_SIEVE_SIZE + (indexM % 10) * SIEVE_SIZE;
					for (i = 0; i < GROUP_SIZE; i++) {
						j = sieveidx[i]; // 0 < J < HALF_SIEVE_SIZE
						if (sieve[J + j] != 0 && sieve[J - 1 - j] != 0) {
							continue; // Do not process if both are composite numbers.
						}
						SubtBigNbrModN(Aux1, root[i], M, TestNbr, NumberLength);
						modmult(GcdAccumulated, M, Aux2);
						memcpy(GcdAccumulated, Aux2, NumberLength * sizeof(limb));
					}
					if (Pass != 0) {
						if (BigNbrIsZero(GcdAccumulated, NumberLength)) {
							// This curve cannot factor the number.
							break;
						}
						if (gcdIsOne(GcdAccumulated) > 1) {
							return FACTOR_FOUND;
						}
					}
				}   // End for.

				if (indexM != 0) { // Update (X:Z)
					memcpy(WX, X, NumberLength * sizeof(limb));
					memcpy(WZ, Z, NumberLength * sizeof(limb));
					SubtBigNbrModN(X, Z, Aux1, TestNbr, NumberLength);
					AddBigNbrModN(TX, TZ, Aux2, TestNbr, NumberLength);
					modmult(Aux1, Aux2, W1);
					AddBigNbrModN(X, Z, Aux1, TestNbr, NumberLength);
					SubtBigNbrModN(TX, TZ, Aux2, TestNbr, NumberLength);
					modmult(Aux1, Aux2, W2);
					AddBigNbrModN(W1, W2, Aux1, TestNbr, NumberLength);
					modmult(Aux1, Aux1, Aux2);
					modmult(Aux2, UZ, X);
					SubtBigNbrModN(W1, W2, Aux1, TestNbr, NumberLength);
					modmult(Aux1, Aux1, Aux2);
					modmult(Aux2, UX, Z);
					memcpy(UX, WX, NumberLength * sizeof(limb));
					memcpy(UZ, WZ, NumberLength * sizeof(limb));
				}
			} // end for Q

			if (Pass == 0) {
				int rc;
				if (BigNbrIsZero(GcdAccumulated, NumberLength)) { // If GcdAccumulated is zero
					memcpy(X, Xaux, NumberLength * sizeof(limb));
					memcpy(Z, Zaux, NumberLength * sizeof(limb));
					continue; // multiple of TestNbr, continue.
				}
				rc = gcdIsOne(GcdAccumulated);
				if (rc == 1) {
					break;    // GCD is one, so this curve does not find a factor.
				}
				if (rc == 0) {
					continue;
				}
				// GD <- GCD(GcdAccumulated, TestNbr)
				if (memcmp(GD, TestNbr, NumberLength * sizeof(limb)))
				{           // GCD is not 1 or TestNbr
					return FACTOR_FOUND;
				}
			}
		} /* end for Pass */

		performLehman = true;
	}       /* End curve calculation */
}

static bool ecm(BigInteger &N, struct sFactors *pstFactors)
{
	int P, Q;
#ifndef __EMSCRIPTEN__
	(void)pstFactors;     // Ignore parameter.
#endif
	fieldTX = TX;
	fieldTZ = TZ;
	fieldUX = UX;
	fieldUZ = UZ;


	fieldAA = AA;
	NumberLength = N.nbrLimbs;
	memcpy(TestNbr, N.limbs, NumberLength * sizeof(limb));
	GetYieldFrequency();  //get yield frequency (used by showECMStatus)
	GetMontgomeryParms(NumberLength);
	memset(M, 0, NumberLength * sizeof(limb));
	memset(DX, 0, NumberLength * sizeof(limb));
	memset(DZ, 0, NumberLength * sizeof(limb));
	memset(W3, 0, NumberLength * sizeof(limb));
	memset(W4, 0, NumberLength * sizeof(limb));
	memset(GD, 0, NumberLength * sizeof(limb));
	ptrLowerText = lowerText;
#ifdef __EMSCRIPTEN__ 
//	*ptrLowerText++ = '3';
	//if (pstFactors->multiplicity > 1)
	//{    // Some factorization known.
	//	int NumberLengthBak = NumberLength;
	//	strcpy(ptrLowerText, "\n");
	//	ptrLowerText += strlen(ptrLowerText);
	//	SendFactorizationToOutput(pstFactors, ptrLowerText, N->sign == SIGN_POSITIVE);
	//	strcpy(ptrLowerText, "\n");
	//	ptrLowerText += strlen(ptrLowerText);
	//	NumberLength = NumberLengthBak;
	//}
	//strcpy(ptrLowerText, lang ? "\nFactorizando " : 
	//	"\nFactorizing using Lenstra elliptic-curve factorization (ECM) ");
	//ptrLowerText += strlen(ptrLowerText);
	//if (hexadecimal)
	//{
	//	Bin2Hex(N->limbs, ptrLowerText, N->nbrLimbs, groupLen);
	//}
	//else
	//{
	//	Bin2Dec(N->limbs, ptrLowerText, N->nbrLimbs, groupLen);
	//}
	//ptrLowerText += strlen(ptrLowerText);
	//strcpy(ptrLowerText, "\n");
	//ptrLowerText += strlen(ptrLowerText);
	//printf("%s", lowerText);
#endif
	EC--; 
	if (SmallPrime[0] != 2) {
		SmallPrime[0] = 2;   // put 1st prime into SmallPrime list
		P = 3;
		indexM = 1;
		for (indexM = 1; indexM < sizeof(SmallPrime) / sizeof(SmallPrime[0]); indexM++)
		{     // Loop that fills the SmallPrime array.
			SmallPrime[indexM] = (int)P; /* Store prime */
			do
			{
				P += 2;
				for (Q = 3; Q * Q <= P; Q += 2)
				{ /* Check if P is prime */
					if (P % Q == 0)
					{
						break;  /* Composite */
					}
				}
			} while (Q * Q <= P);
		}
	}
	foundByLehman = false;
	do
	{
		enum eEcmResult ecmResp = ecmCurve(N);
		if (ecmResp == CHANGE_TO_SIQS)
		{    // Perform SIQS
			FactoringSIQSx(TestNbr, GD);
			break;
		}
		else if (ecmResp == FACTOR_FOUND)
		{
			break;
		}
		else if (ecmResp == ERROR)
			return false;
	} while (!memcmp(GD, TestNbr, NumberLength * sizeof(limb)));
#if 0
	lowerTextArea.setText("");
#endif
	StepECM = 0; /* do not show pass number on screen */
	return true;
}

/* convert factors from integer lists to Znums */
static void ConvertFactors(const sFactors *pstFactors, std::vector <Znum> &factorlist,
	std::vector<int> &exponentlist) {
	int ix = 1;
	Znum primefactor;
	factorlist.clear();
	exponentlist.clear();

	/* 1st entry in pstFactors contains the total number of factors.
	step through list of prime factors. */
	for (ix = 1; ix <= pstFactors[0].multiplicity; ix++) {
		ValuestoZ(primefactor, pstFactors[ix].ptrFactor); // convert prime factor
		factorlist.push_back(primefactor);                // and copy it
		exponentlist.push_back(pstFactors[ix].multiplicity);  // copy exponent
	}
}

static void SortFactors(struct sFactors *pstFactors)
{
	int factorNumber, factorNumber2, ctr;
	struct sFactors *pstCurFactor = pstFactors + 1;
	struct sFactors stTempFactor;
	int *ptrNewFactor;
	struct sFactors *pstNewFactor;
	for (factorNumber = 1; factorNumber <= pstFactors->multiplicity; factorNumber++, pstCurFactor++)
	{
		pstNewFactor = pstCurFactor + 1;
		for (factorNumber2 = factorNumber + 1; factorNumber2 <= pstFactors->multiplicity; factorNumber2++, pstNewFactor++)
		{
			int *ptrFactor = pstCurFactor->ptrFactor;
			int *ptrFactor2 = pstNewFactor->ptrFactor;
			if (*ptrFactor < *ptrFactor2)
			{     // Factors already in correct order.
				continue;
			}
			if (*ptrFactor == *ptrFactor2)
			{
				for (ctr = *ptrFactor; ctr > 1; ctr--)
				{
					if (*(ptrFactor + ctr) != *(ptrFactor2 + ctr))
					{
						break;
					}
				}
				if (*(ptrFactor + ctr) < *(ptrFactor2 + ctr))
				{     // Factors already in correct order.
					continue;
				}
				if (*(ptrFactor + ctr) == *(ptrFactor2 + ctr))
				{     // Factors are the same.
					pstCurFactor->multiplicity += pstNewFactor->multiplicity;
					ctr = pstFactors->multiplicity - factorNumber2;
					if (ctr > 0)
					{
						memmove(pstNewFactor, pstNewFactor + 1, ctr * sizeof(struct sFactors));
					}
					pstFactors->multiplicity--;   // Indicate one less known factor.
					continue;
				}
			}
			// Exchange both factors.
			memcpy(&stTempFactor, pstCurFactor, sizeof(struct sFactors));
			memcpy(pstCurFactor, pstNewFactor, sizeof(struct sFactors));
			memcpy(pstNewFactor, &stTempFactor, sizeof(struct sFactors));
		}
	}
	// Find location for new factors.
	ptrNewFactor = 0;
	pstCurFactor = pstFactors + 1;
	for (factorNumber = 1; factorNumber <= pstFactors->multiplicity; factorNumber++, pstCurFactor++)
	{
		int *ptrPotentialNewFactor = pstCurFactor->ptrFactor + *(pstCurFactor->ptrFactor) + 1;
		if (ptrPotentialNewFactor > ptrNewFactor)
		{
			ptrNewFactor = ptrPotentialNewFactor;
		}
	}
	pstFactors->ptrFactor = ptrNewFactor;
}

// Insert new factor found into factor array. This factor array must be sorted.
static void insertIntFactor(struct sFactors *pstFactors, struct sFactors *pstFactorDividend, int divisor)
{
	struct sFactors *pstCurFactor;
	int factorNumber;
	int *ptrFactor = pstFactorDividend->ptrFactor;
	int nbrLimbs = *ptrFactor;
	int *ptrValue;
	pstFactorDividend->upperBound = divisor;
	// Divide number by factor just found.
	DivBigNbrByInt(ptrFactor + 1, divisor, ptrFactor + 1, nbrLimbs);
	if (*(ptrFactor + nbrLimbs) == 0)
	{
		(*ptrFactor)--;
	}
	// Check whether prime is already in factor list.
	pstCurFactor = pstFactors + 1;
	for (factorNumber = 1; factorNumber <= pstFactors->multiplicity; factorNumber++, pstCurFactor++)
	{
		ptrValue = pstCurFactor->ptrFactor;  // Point to factor in factor array.
		if (*ptrValue == 1 && *(ptrValue + 1) == divisor)
		{  // Prime already found: increment multiplicity and go out.
			pstCurFactor->multiplicity += pstFactorDividend->multiplicity;
			ptrValue = pstFactorDividend->ptrFactor;
			if (*ptrValue == 1 && *(ptrValue + 1) == 1)
			{    // Dividend is 1 now so discard it.
				*pstFactorDividend = *(pstFactors + pstFactors->multiplicity--);
			}
			SortFactors(pstFactors);
			return;
		}
		if (*ptrValue > 1 || *(ptrValue + 1) > divisor)
		{   // Factor in factor list is greater than factor to insert. Exit loop.
			break;
		}
	}
	pstFactors->multiplicity++; // Indicate new known factor.
								// Move all elements.
	ptrValue = pstFactorDividend->ptrFactor;
	if (pstFactors->multiplicity > factorNumber)
	{
		memmove(pstCurFactor + 1, pstCurFactor,
			(pstFactors->multiplicity - factorNumber) * sizeof(struct sFactors));
	}
	if (*ptrValue == 1 && *(ptrValue + 1) == 1)
	{
		pstCurFactor = pstFactorDividend;
	}
	else
	{
		ptrValue = pstFactors->ptrFactor;
		pstCurFactor->ptrFactor = ptrValue;
		pstCurFactor->multiplicity = pstFactorDividend->multiplicity;
		pstFactors->ptrFactor += 2;  // Next free memory.
	}
	pstCurFactor->upperBound = 0;
	*ptrValue = 1;  // Number of limbs.
	*(ptrValue + 1) = divisor;
	SortFactors(pstFactors);
}

// Insert new factor found into factor array. This factor array must be sorted.
// The divisor must be also sorted.
static void insertBigFactor(struct sFactors *pstFactors, BigInteger &divisor)
{
	struct sFactors *pstCurFactor;
	int factorNumber;
	int lastFactorNumber = pstFactors->multiplicity;
	struct sFactors *pstNewFactor = pstFactors + lastFactorNumber + 1;
	int *ptrNewFactorLimbs = pstFactors->ptrFactor;
	pstCurFactor = pstFactors + 1;
	for (factorNumber = 1; factorNumber <= lastFactorNumber; factorNumber++, pstCurFactor++)
	{     // For each known factor...
		int *ptrFactor = pstCurFactor->ptrFactor;
		NumberLength = *ptrFactor;
		UncompressBigInteger(ptrFactor, Temp2);    // Convert known factor to Big Integer.
		BigIntGcd(divisor, Temp2, Temp3);          // Temp3 is the GCD between known factor and divisor.
		if (Temp3 < 2) { 
			continue; // divisor is not a new factor (GCD = 0 or 1).
		}
		//if (TestBigNbrEqual(Temp2, Temp3)) 	{ 
		if (Temp2 == Temp3) {
			continue; // GCD is equal to known factor.
		}
		// At this moment both GCD and known factor / GCD are new known factors. Replace the known factor by
		// known factor / GCD and generate a new known factor entry.
		NumberLength = Temp3.nbrLimbs;
		CompressBigInteger(ptrNewFactorLimbs, Temp3);      // Append new known factor.
		Temp4 = Temp2 / Temp3; // BigIntDivide(Temp2, Temp3, Temp4);     // Divide by this factor.
		NumberLength = Temp4.nbrLimbs;
		CompressBigInteger(ptrFactor, Temp4);              // Overwrite old known factor.
		pstNewFactor->multiplicity = pstCurFactor->multiplicity;
		pstNewFactor->ptrFactor = ptrNewFactorLimbs;
		pstNewFactor->upperBound = pstCurFactor->upperBound;
		pstNewFactor++;
		pstFactors->multiplicity++;
		ptrNewFactorLimbs += 1 + Temp3.nbrLimbs;
	}
	// Sort factors in ascending order. If two factors are equal, coalesce them.
	// Divide number by factor just found.
	SortFactors(pstFactors);
	return;
}

#ifdef __EMSCRIPTEN__
static void showECMStatus(void)
{
	char status[200];
	int elapsedTime;
	char *ptrStatus;
	if ((lModularMult % yieldFreq) != 0)
	{
		return;
	}
	elapsedTime = (int)(tenths() - originalTenthSecond);
	if (elapsedTime / 10 == oldTimeElapsed / 10)
	{
		return;
	}
	oldTimeElapsed = elapsedTime;
	ptrStatus = status;
	strcpy(ptrStatus, lang ? "4\nTranscurrió " : "4\nTime elapsed: ");
	ptrStatus += strlen(ptrStatus);
	GetDHMS(&ptrStatus, elapsedTime / 10);
	switch (StepECM)
	{
	case 1:
		strcpy(ptrStatus, lang ? "Paso 1: " : "Step 1: ");
		ptrStatus += strlen(ptrStatus);
		int2dec(&ptrStatus, indexPrimes / (nbrPrimes / 100));
		*ptrStatus++ = '%';
		break;
	case 2:
		strcpy(ptrStatus, lang ? "Paso 2: " : "Step 2: ");
		ptrStatus += strlen(ptrStatus);
		int2dec(&ptrStatus, maxIndexM == 0 ? 0 : indexM / (maxIndexM / 100));
		*ptrStatus++ = '%';
		break;
	}
	strcpy(ptrStatus, "\n");
	printf ("%s\n", status);
}

#endif

/* called for Carmichael numbers that have no small factors. 
Return: false = No factors found, true = factors found.
 Use: Xaux for square root of -1.
      Zaux for square root of 1. */
static bool factorCarmichael(const BigInteger &pValue, struct sFactors *pstFactors)
{
	int randomBase = 0;  // pseudo-random number
	bool factorsFound = false;
	int nbrLimbsQ, countdown, ctr;
	int nbrLimbs = pValue.nbrLimbs;
	bool sqrtOneFound = false;
	int sqrtMinusOneFound = false;
	int Aux1Len;
	limb *pValueLimbs = (limb *)pValue.limbs;  // override const-ness of Pvalue
	(pValueLimbs + nbrLimbs)->x = 0;         // doesn't change value of pValue
	memcpy(q, pValueLimbs, (nbrLimbs + 1) * sizeof(limb)); // copy p to q
	nbrLimbsQ = nbrLimbs;
	q[0]--;                     // q = p - 1 (p is odd, so there is no carry).
	memcpy(Aux1, q, (nbrLimbsQ + 1) * sizeof(q[0])); // copy q to Aux1
	Aux1Len = nbrLimbs;
	DivideBigNbrByMaxPowerOf2(&ctr, Aux1, &Aux1Len);
	memcpy(TestNbr, pValueLimbs, nbrLimbs * sizeof(limb)); // copy p to TestNbr
	TestNbr[nbrLimbs].x = 0;
	GetMontgomeryParms(nbrLimbs);
	for (countdown = 20; countdown > 0; countdown--) {
		int i;
		NumberLength = nbrLimbs;
		randomBase = (int)((uint64_t)randomBase * 89547121 + 1762281733) & MAX_INT_NBR;
		modPowBaseInt(randomBase, Aux1, Aux1Len, Aux2); // Aux2 = base^Aux1.
		// If Mult1 = 1 or Mult1 = TestNbr-1, then try next base.
		if (checkOne(Aux2, nbrLimbs) != 0 || checkMinusOne(Aux2, nbrLimbs) != 0)
		{
			continue;    // This base cannot find a factor. Try another one.
		}
		for (i = 0; i < ctr; i++)
		{              // Loop that squares number.
			modmult(Aux2, Aux2, Aux3);
			if (checkOne(Aux3, nbrLimbs) != 0)
			{            // Non-trivial square root of 1 found.
				if (!sqrtOneFound)
				{          // Save it to perform GCD later.
					memcpy(Zaux, Aux2, nbrLimbs * sizeof(limb));
					sqrtOneFound = true;
				}
				else
				{          // Try to find non-trivial factor by doing GCD.
					SubtBigNbrMod(Aux2, Zaux, Aux4);
					UncompressLimbsBigInteger(Aux4, Temp2);
					BigIntGcd(pValue, Temp2, Temp4);
					if ((Temp4 > 1) && (Temp4 != pValue))
					{          // Non-trivial factor found.
						if (pstFactors->multiplicity >= Max_Factors) {
							// factor list is full!
							std::string line = std::to_string(__LINE__);
							std::string mesg = "cannot factorise: too many factors: ";
							mesg += __func__;                 // + function name
							mesg += " line ";  mesg += line;  // + line number
							mesg += " in file "; mesg += __FILE__;   // + file name
							throw std::range_error(mesg);
						}
						insertBigFactor(pstFactors, Temp4);
						factorsFound = true;
					}
				}
				// Try to find non-trivial factor by doing GCD.
				NumberLength = nbrLimbs;
				AddBigNbrMod(Aux2, MontgomeryMultR1, Aux4);
				UncompressLimbsBigInteger(Aux4, Temp2);
				BigIntGcd(pValue, Temp2, Temp4);
				if ((Temp4.nbrLimbs != 1 || Temp4.limbs[0].x > 1) &&
					(Temp4.nbrLimbs != NumberLength ||
						memcmp(pValue.limbs, Temp4.limbs, NumberLength * sizeof(limb))))
				{          // Non-trivial factor found.
					if (pstFactors->multiplicity >= Max_Factors) {
						// factor list is full!
						std::string line = std::to_string(__LINE__);
						std::string mesg = "cannot factorise: too many factors: ";
						mesg += __func__;                 // + function name
						mesg += " line ";  mesg += line;  // + line number
						mesg += " in file "; mesg += __FILE__;   // + file name
						throw std::range_error(mesg);
					}
					insertBigFactor(pstFactors, Temp4);
					factorsFound = true;
				}
				i = ctr;
				continue;  // Find more factors.
			}
			if (checkMinusOne(Aux3, nbrLimbs) != 0)
			{            // Square root of 1 found.
				if (!sqrtMinusOneFound)
				{          // Save it to perform GCD later.
					memcpy(Xaux, Aux2, nbrLimbs * sizeof(limb));
					sqrtOneFound = true;
				}
				else
				{          // Try to find non-trivial factor by doing GCD.
					SubtBigNbrMod(Aux3, Xaux, Aux4);
					UncompressLimbsBigInteger(Aux4, Temp2);
					BigIntGcd(pValue, Temp2, Temp4);
					if ((Temp4.nbrLimbs != 1 || Temp4.limbs[0].x > 1) &&
						(Temp4.nbrLimbs != NumberLength ||
							memcmp(pValue.limbs, Temp4.limbs, NumberLength * sizeof(limb))))
					{          // Non-trivial factor found.
						if (pstFactors->multiplicity >= Max_Factors) {
							// factor list is full!
							std::string line = std::to_string(__LINE__);
							std::string mesg = "cannot factorise: too many factors: ";
							mesg += __func__;                 // + function name
							mesg += " line ";  mesg += line;  // + line number
							mesg += " in file "; mesg += __FILE__;   // + file name
							throw std::range_error(mesg);
						}
						insertBigFactor(pstFactors, Temp4);
						factorsFound = true;
					}
				}
				i = ctr;
				continue;  // Find more factors.
			}
			memcpy(Aux2, Aux3, nbrLimbs * sizeof(limb));
		}
	}
	return factorsFound;
}

// pstFactors -> ptrFactor points to end of factors.
// pstFactors -> multiplicity indicates the number of different factors.
static bool factor (const BigInteger &toFactor, struct sFactors *pstFactors)
{
	struct sFactors *pstCurFactor;
	int factorNbr, expon;
	int remainder, nbrLimbs, ctr;
	int *ptrFactor;
	int dividend;
	bool restartFactoring = false;
	int result;
	
	EC = 1;  // start with 1st curve
	NumberLength = toFactor.nbrLimbs;
	CompressBigInteger(nbrToFactor, toFactor);  // convert Biginteger to integer list
	GetYieldFrequency();   //get yield frequency (used by showECMStatus)
#ifdef __EMSCRIPTEN__
	oldTimeElapsed = 0;
	originalTenthSecond = tenths();    // record start time
	modmultCallback = showECMStatus;   // Set callback function pointer
#endif
	/* initialise factor list */
	pstFactors->multiplicity = 1;
	pstFactors->ptrFactor = nbrToFactor + 1 + *nbrToFactor;
	pstFactors->upperBound = 0;

	pstCurFactor = pstFactors + 1;
	pstCurFactor->multiplicity = 1;
	pstCurFactor->ptrFactor = nbrToFactor;
	pstCurFactor->upperBound = 2;
	if (toFactor.nbrLimbs > 1) {
		PowerPM1Check(pstFactors, toFactor);  // check if toFactor is a perfect power +/- 1
	}
	for (factorNbr = 1; factorNbr <= pstFactors->multiplicity; factorNbr++, pstCurFactor++) {
		int upperBound = pstCurFactor->upperBound;
		restartFactoring = false;
		// If number is prime, do not process it.
		if (upperBound == 0) {    
			continue; // Factor is prime.
		}
		ptrFactor = pstCurFactor->ptrFactor;
		nbrLimbs = *ptrFactor;
		while (upperBound < 100000 && nbrLimbs > 1)
		{        // Number has at least 2 limbs: Trial division by small numbers.
			while (pstCurFactor->upperBound != 0) {   // Factor found.
				ptrFactor = pstCurFactor->ptrFactor;
				assert(*ptrFactor < 100000);
				remainder = RemDivBigNbrByInt(ptrFactor + 1, upperBound, nbrLimbs);
				if (remainder != 0) {    
					break;  // Factor not found. Use new divisor.
				}
				if (pstFactors->multiplicity >= Max_Factors) {  
					// factor list is full!
					std::string line = std::to_string(__LINE__);
					std::string mesg = "cannot factorise: too many factors: ";
					mesg += __func__;                 // + function name
					mesg += " line ";  mesg += line;  // + line number
					mesg += " in file "; mesg += __FILE__;   // + file name
					throw std::range_error(mesg);
				}
				insertIntFactor(pstFactors, pstCurFactor, upperBound);
				restartFactoring = true;
				break;
			}
			if (restartFactoring) {
				factorNbr = 0;
				pstCurFactor = pstFactors;
				break;
			}
			if (nbrLimbs == 1) {     // Number completely factored.
				break;
			}
			if (upperBound == 2) {
				upperBound++;
			}
			else if (upperBound % 6 == 1) {
				upperBound += 4;  // avoid setting upperBound to any multiple of 3
			}
			else {
				upperBound += 2;
			}
		}
		if (restartFactoring) {
			continue;
		}

		if (nbrLimbs == 1) {
			dividend = *(ptrFactor + 1);
			while ((unsigned int)upperBound*(unsigned int)upperBound <= (unsigned int)dividend)
			{              // Trial division by small numbers.
				if (dividend % upperBound == 0) {            // Factor found.
					if (pstFactors->multiplicity >= Max_Factors) {
						// factor list is full!
						std::string line = std::to_string(__LINE__);
						std::string mesg = "cannot factorise: too many factors: ";
						mesg += __func__;                 // + function name
						mesg += " line ";  mesg += line;  // + line number
						mesg += " in file "; mesg += __FILE__;   // + file name
						throw std::range_error(mesg);
					}
					insertIntFactor(pstFactors, pstCurFactor, upperBound);
					restartFactoring = true;
					break;
				}
				if (upperBound == 2) {
					upperBound++;
				}
				else {
					upperBound += 2;
				}
			}
			if (restartFactoring) {
				factorNbr = 0;
				pstCurFactor = pstFactors;
				continue;
			}
			pstCurFactor->upperBound = 0;   // Number is prime.
			continue;
		}

		// No small factor. Check whether the number is prime or prime power.
		NumberLength = *pstCurFactor->ptrFactor;
		UncompressBigInteger(pstCurFactor->ptrFactor, power);
		NumberLength = power.nbrLimbs;
		expon = PowerCheck(power, prime);
		if (expon > 1) {
			CompressBigInteger(pstCurFactor->ptrFactor, prime);
			pstCurFactor->multiplicity *= expon;
		}
		result = BpswPrimalityTest(prime);
		if (result == 0) {   // Number is prime power.
			pstCurFactor->upperBound = 0;   // Indicate that number is prime.
			continue;                       // Check next factor.
		}
		if (result > 1) {    // Number is 2-Fermat probable prime. Try to factor it.
			if (factorCarmichael(prime, pstFactors)) {  // Prime factors found.
				continue;
			}
		}
		auto rv = ecm(prime, pstFactors);          // Factor number.
		if (!rv) 
			return false;
			// Check whether GD is not one. In this case we found a proper factor.
		for (ctr = 1; ctr < NumberLength; ctr++) {
			if (GD[ctr].x != 0) {
				break;
			}
		}
		if (ctr != NumberLength || GD[0].x != 1) {
			int numLimbs;
			Temp1.sign = SIGN_POSITIVE;
			numLimbs = NumberLength;
			while (numLimbs > 1) {    // adjust count of number of limbs
				if (GD[numLimbs - 1].x != 0) {
					break;
				}
				numLimbs--;
			}
			/* copy number from GD and convert it to a BigInteger */
			memcpy(Temp1.limbs, GD, numLimbs * sizeof(limb));
			Temp1.nbrLimbs = numLimbs;
			if (pstFactors->multiplicity >= Max_Factors) {
				// factor list is full!
				std::string line = std::to_string(__LINE__);
				std::string mesg = "cannot factorise: too many factors: ";
				mesg += __func__;                 // + function name
				mesg += " line ";  mesg += line;  // + line number
				mesg += " in file "; mesg += __FILE__;   // + file name
				throw std::range_error(mesg);
			}
			insertBigFactor(pstFactors, Temp1);
			factorNbr = 0;
			pstCurFactor = pstFactors;
		}    // End if
	}      // End for
#ifdef __EMSCRIPTEN__
#ifdef FACTORIZATION_APP
	SaveFactors(pstFactors);
#endif
#endif
	return true;
}

/* show that the number is the sum of 4 or fewer squares. See
https://www.alpertron.com.ar/4SQUARES.HTM */
static void ComputeFourSquares(struct sFactors *pstFactors) {
	int indexPrimes;
	static BigInteger p, q, K, Mult1, Mult2, Mult3, Mult4;  // follow advice not to use stack for these
	static BigInteger Tmp, Tmp1, Tmp2, Tmp3, Tmp4, M1, M2, M3, M4;  // follow advice not to use stack for these
	struct sFactors *pstFactor;
	static limb minusOneMont[MAX_LEN];

	Quad1 = 1;      // 1 = 1^2 + 0^2 + 0^2 + 0^2
	Quad2 = 0; 
	Quad3 = 0; 
	Quad4 = 0; 
	int Factorix = 1;      // Point to first factor in array of factors.
	if (pstFactors[0].multiplicity == 1 && pstFactors[1].multiplicity == 1) {
		/* only 1 factor i.e. number is prime */
		if (*((pstFactors[1].ptrFactor) + 1) == 1) 	{   // Number to factor is 1.
			return;
		}
		if (*((pstFactors[1].ptrFactor) + 1) == 0) {    // Number to factor is 0.
			Quad1 = 0; // intToBigInteger(Quad1, 0);     // 0 = 0^2 + 0^2 + 0^2 + 0^2
			return;
		}
	}
	for (Factorix = 1; Factorix <= pstFactors[0].multiplicity; Factorix++) {
		if (pstFactors[Factorix].multiplicity % 2 == 0) {  // Prime factor appears twice.
			continue;
		}
		NumberLength = *pstFactors[Factorix].ptrFactor;
		UncompressBigInteger(pstFactors[Factorix].ptrFactor, p);
		p.sign = SIGN_POSITIVE;
		q = p; 
		q --;            // q <- p-1
		if (p == 2) {   /* Prime factor is 2 */
			Mult1 = 1;  // 2 = 1^2 + 1^2 + 0^2 + 0^2
			Mult2 = 1; 
			Mult3 = 0; 
			Mult4 = 0; 
		}
		else {       /* Prime factor is not 2 */
			NumberLength = p.nbrLimbs;
			memcpy(&TestNbr, p.limbs, NumberLength * sizeof(limb));
			TestNbr[NumberLength].x = 0;
			GetMontgomeryParms(NumberLength);
			memset(minusOneMont, 0, NumberLength * sizeof(limb));
			SubtBigNbrModN(minusOneMont, MontgomeryMultR1, minusOneMont, TestNbr, NumberLength);
			memset(K.limbs, 0, NumberLength * sizeof(limb));
			if ((p.lldata() & 3) == 1) { /* if p = 1 (mod 4) */
				q = p; 
				subtractdivide(q, 1, 4);     // q = (prime-1)/4
				K = 1;
				do {    // Loop that finds mult1 = sqrt(-1) mod prime in Montgomery notation.
					K++;
					modPow(K.limbs, q.limbs, q.nbrLimbs, Mult1.limbs);
				} while (!memcmp(Mult1.limbs, MontgomeryMultR1, NumberLength * sizeof(limb)) ||
					!memcmp(Mult1.limbs, minusOneMont, NumberLength * sizeof(limb)));

				Mult1.sign = SIGN_POSITIVE;
				memset(Mult2.limbs, 0, p.nbrLimbs * sizeof(limb));
				Mult2 = 1;
				//Mult2.nbrLimbs = 1;
				//Mult2.sign = SIGN_POSITIVE;
				// Convert Mult1 to standard notation by multiplying by 1 in
				// Montgomery notation.
				modmult(Mult1.limbs, Mult2.limbs, Mult3.limbs);
				memcpy(Mult1.limbs, Mult3.limbs, p.nbrLimbs * sizeof(limb));
				for (Mult1.nbrLimbs = p.nbrLimbs; Mult1.nbrLimbs > 1; Mult1.nbrLimbs--)
				{  // Adjust number of limbs so the most significant limb is not zero.
					if (Mult1.limbs[Mult1.nbrLimbs - 1].x != 0) {
						break;
					}
				}
				for (;;) {
					Tmp = Mult1 * Mult1 + Mult2 * Mult2; 
					K = Tmp / p;        // K <- (mult1^2 + mult2^2) / p
					if (K == 1) {  
						Mult3 = 0; 
						Mult4 = 0; 
						break;
					}
					M1 = Mult1 % K; //BigIntRemainder(Mult1, K, M1);   
					if (M1 < 0) {
						M1 += K; // M1 = M1 + K; 
					}
					M2 = Mult2 % K;  // BigIntRemainder(Mult2, K, M2)
					if (M2 < 0) {
						M2 += K; //M2 = M2 + K; 
					}
					Tmp = K; // CopyBigInt(Tmp, K);
					subtractdivide(Tmp, -1, 2);       // Tmp <- (K+1) / 2
					if (M1 >= Tmp) // If M1 >= K / 2 ... 
					{
						M1 -= K; //M1 = M1 - K; // BigIntSubt(M1, K, M1);        /
					}
					if (M2 >= Tmp) {     // If M2 >= K / 2 ...     
						M2 -= K; // M2 = M2 - K;  //BigIntSubt(M2, K, M2);         
					}
					Tmp = Mult1*M1 + Mult2*M2;
					Tmp2 = Tmp / K;  // Tmp2 <- (mult1*m1 + mult2*m2) / K
					Tmp = Mult1*M2 - Mult2*M1; 
					Mult2 = Tmp / K;    // Mult2 <- (mult1*m2 - mult2*m1) /K
					Mult1 = Tmp2; 
				} /* end while */
			} /* end p = 1 (mod 4) */

			else { /* if p = 3 (mod 4) */
				int mult1 = 0;
				q = p; //  CopyBigInt(q, p);
				subtractdivide(q, 1, 2);     // q = (prime-1)/2
				memcpy(K.limbs, q.limbs, q.nbrLimbs * sizeof(limb));
				if (p.nbrLimbs > q.nbrLimbs) {
					K.limbs[q.nbrLimbs].x = 0;
				}
				// Compute Mult1 and Mult2 so Mult1^2 + Mult2^2 = -1 (mod p)
				memset(Mult1.limbs, 0, p.nbrLimbs * sizeof(limb));
				do {
					mult1++;
					// Increment Mult1 by 1 in Montgomery notation.
					AddBigNbrModN(Mult1.limbs, MontgomeryMultR1, Mult1.limbs, p.limbs, p.nbrLimbs);
					modmult(Mult1.limbs, Mult1.limbs, Tmp.limbs);
					SubtBigNbrModN(minusOneMont, Tmp.limbs, Tmp.limbs, p.limbs, p.nbrLimbs);
					modPow(Tmp.limbs, K.limbs, p.nbrLimbs, Tmp1.limbs);
					// At this moment Tmp1 = (-1 - Mult1^2)^((p-1)/2)
					// in Montgomery notation. Continue loop if it is not 1.
				} while (memcmp(Tmp1.limbs, MontgomeryMultR1, p.nbrLimbs));

				// After the loop finishes, Tmp = (-1 - Mult1^2) is a quadratic residue mod p.
				// Convert Mult1 to standard notation by multiplying by 1 in
				// Montgomery notation.
				Mult1 = mult1; 
				q = p; 
				subtractdivide(q, -1, 4);  // q <- (p+1)/4.
				// Find Mult2 <- square root of Tmp = Tmp^q (mod p) in Montgomery notation.
				modPow(Tmp.limbs, q.limbs, p.nbrLimbs, Mult2.limbs);
				// Convert Mult2 from Montgomery notation to standard notation.
				memset(Tmp.limbs, 0, p.nbrLimbs * sizeof(limb));
				Tmp.limbs[0].x = 1;
				Mult3 = 1; 
				Mult4 = 0; 
				// Convert Mult2 to standard notation by multiplying by 1 in
				// Montgomery notation.
				modmult(Mult2.limbs, Tmp.limbs, Mult2.limbs);
				for (Mult2.nbrLimbs = p.nbrLimbs; Mult2.nbrLimbs > 1; Mult2.nbrLimbs--)
				{  // Adjust number of limbs so the most significant limb is not zero.
					if (Mult2.limbs[Mult2.nbrLimbs - 1].x != 0) {
						break;
					}
				}
				Mult2.sign = SIGN_POSITIVE;
				for (;;) {
					// Compute K <- (Mult1^2 + Mult2^2 + Mult3^2 + Mult4^2) / p
					Tmp = Mult1*Mult1 + Mult2*Mult2; 
					Tmp = Tmp + Mult3*Mult3 + Mult4*Mult4; 
					K = Tmp / p; 
					if (K == 1) {   
						break;  // K equals 1
					}
					if (K.isEven()) { // If K is even ...
						if (Mult1.isEven() != Mult2.isEven())
						{  // If Mult1 + Mult2 is odd...
							if (Mult1.isEven() == Mult3.isEven())
							{   // If Mult1 + Mult3 is even...
								Tmp = Mult2;  // swap mult2 and mult3
								Mult2 = Mult3; 
								Mult3 = Tmp;  
							}
							else {
								Tmp = Mult2; // swap mult2 and mult4
								Mult2 = Mult4; 
								Mult4 = Tmp; 
							}
						} // At this moment Mult1+Mult2 = even, Mult3+Mult4 = even
						Tmp1 = Mult1 + Mult2; 
						subtractdivide(Tmp1, 0, 2);  // Tmp1 <- (Mult1 + Mult2) / 2
						Tmp2 = Mult1 - Mult2;
						subtractdivide(Tmp2, 0, 2);  // Tmp2 <- (Mult1 - Mult2) / 2
						Tmp3 = Mult3 + Mult4; 
						subtractdivide(Tmp3, 0, 2);  // Tmp3 <- (Mult3 + Mult4) / 2
						Mult4 = Mult3 - Mult4; 
						subtractdivide(Mult4, 0, 2);  // Mult4 <- (Mult3 - Mult4) / 2
						Mult3 = Tmp3; 
						Mult2 = Tmp2; 
						Mult1 = Tmp1;  
						continue;
					} /* end if k is even */
					M1 = Mult1 % K;    //BigIntRemainder(Mult1, K, M1) .
					if (M1 < 0) {
						M1 += K; // M1 = M1 + K; // BigIntAdd(M1, K, M1);
					}
					M2 = Mult2 % K;    // BigIntRemainder(Mult2, K, M2).
					if (M2 < 0) {
						M2 = M2 + K; // BigIntAdd(M2, K, M2);
					}
					M3 = Mult3 % K;    // BigIntRemainder(Mult3, K, M3)
					if (M3 < 0) {
						M3 += K; // M3 = M3 + K; // BigIntAdd(M3, K, M3);
					}
					M4 = Mult4 % K;    //BigIntRemainder(Mult4, K, M4) .
					if (M4 < 0) {
						M4 += K; // M4 = M4 + K; 
					}
					Tmp = K; 
					subtractdivide(Tmp, -1, 2);       // Tmp <- (K+1) / 2
					if (M1 >= Tmp) { // If M1 >= K / 2 ... 
						M1 -= K; // M1 = M1 - K;    
					}
					if (M2 >= Tmp) { // If M2 >= K / 2 ... 
						M2 -= K; // M2 = M2 - K;         
					}
	
					if (M3 >= Tmp) {  // If M3 >= K / 2 ... 
						M3 -= K; // M3 = M3 - K;        
					}
					if (M4 >= Tmp) { // If M4 >= K / 2 ... 
						M4 -= K; //M4 = M4 - K;      // BigIntSubt(M4, K, M4);    
					}
					// Compute Tmp1 <- (Mult1*M1 + Mult2*M2 + Mult3*M3 + Mult4*M4) / K
					Tmp = Mult1*M1 + Mult2*M2 + Mult3*M3 + Mult4*M4; 
					Tmp1 = Tmp / K; // BigIntDivide(Tmp, K, Tmp1);

					// Compute Tmp2 <- (Mult1*M2 - Mult2*M1 + Mult3*M4 - Mult4*M3) / K
					Tmp = Mult1*M2 - Mult2*M1 + Mult3*M4 - Mult4*M3; 
					Tmp2 = Tmp / K; 

					// Compute Tmp3 <- (Mult1*M3 - Mult3*M1 - Mult2*M4 + Mult4*M2) / K
					Tmp = Mult1*M3 - Mult3*M1 - Mult2*M4 + Mult4*M2; 
					Tmp3 = Tmp / K;

					// Compute Mult4 <- (Mult1*M4 - Mult4*M1 + Mult2*M3 - Mult3*M2) / K
					Tmp = Mult1*M4 - Mult4*M1 + Mult2*M3 - Mult3*M2; // BigIntSubt(Tmp, Tmp4, Tmp);
					Mult4 = Tmp / K; 

					Mult3 = Tmp3; 
					Mult2 = Tmp2;  
					Mult1 = Tmp1; 
				} /* end while */
			} /* end if p = 3 (mod 4) */
		} /* end prime not 2 */

		  // Compute Tmp1 <- Mult1*Quad1 + Mult2*Quad2 + Mult3*Quad3 + Mult4*Quad4
		Tmp1 = Mult1*Quad1 + Mult2*Quad2 + Mult3*Quad3 + Mult4*Quad4; 

		// Compute Tmp2 <- Mult1*Quad2 - Mult2*Quad1 + Mult3*Quad4 - Mult4*Quad3
		Tmp2 = Mult1*Quad2 - Mult2*Quad1 + Mult3*Quad4 - Mult4*Quad3; 

		// Compute Tmp3 <- Mult1*Quad3 - Mult3*Quad1 - Mult2*Quad4 + Mult4*Quad2
		Tmp3 = Mult1*Quad3 - Mult3*Quad1 - Mult2*Quad4 + Mult4*Quad2; 

		// Compute Quad4 <- Mult1*Quad4 - Mult4*Quad1 + Mult2*Quad3 - Mult3*Quad2
		Quad4 = Mult1*Quad4 - Mult4*Quad1 + Mult2*Quad3 - Mult3*Quad2; 

		Quad3 = Tmp3; 
		Quad2 = Tmp2; 
		Quad1 = Tmp1; 
	} /* end for indexPrimes */

	pstFactor = pstFactors + 1;      // Point to first factor in array of factors.
	for (indexPrimes = pstFactors->multiplicity - 1; indexPrimes >= 0; indexPrimes--, pstFactor++)
	{
		NumberLength = *pstFactor->ptrFactor;
		UncompressBigInteger(pstFactor->ptrFactor, p);
		BigIntPowerIntExp(p, pstFactor->multiplicity / 2, K);
		Quad1 *= K;  // Quad1 = Quad1*K; 
		Quad2 *= K; 
		Quad3 *= K; 
		Quad4 *= K; 
	}

	Quad1.sign = SIGN_POSITIVE;
	Quad2.sign = SIGN_POSITIVE;
	Quad3.sign = SIGN_POSITIVE;
	Quad4.sign = SIGN_POSITIVE;

	// Sort squares: largest in Quad1, smallest in Quad4
	
	if(Quad1 < Quad2) {		// Quad1 < Quad2, so exchange them.
		Tmp = Quad1; 
		Quad1 = Quad2;  
		Quad2 = Tmp; 
	}

	if (Quad1 < Quad3) {	// Quad1 < Quad3, so exchange them.
		Tmp = Quad1; 
		Quad1 = Quad3; 
		Quad3 = Tmp; 
	}

	if (Quad1 < Quad4) {	// Quad1 < Quad4, so exchange them.
		Tmp = Quad1; 
		Quad1 = Quad4; 
		Quad4 = Tmp; 
	}

	if (Quad2 < Quad3) {	// Quad2 < Quad3, so exchange them.
		Tmp = Quad2; 
		Quad2 = Quad3; 
		Quad3 = Tmp; 
	}

	if (Quad2 < Quad4) {	// Quad2 < Quad4, so exchange them.
		Tmp = Quad2; 
		Quad2 = Quad4; 
		Quad4 = Tmp; 
	}

	if (Quad3 < Quad4) {	// Quad3 < Quad4, so exchange them.
		Tmp = Quad3; 
		Quad3 = Quad4; 
		Quad4 = Tmp; 
	}
	return;
}

/********************************************************************************
code below is to interface between DA's code that uses BigIntegers and limbs,
and the new code that uses Znums, which are really mpz_t integers from MPIR or GMP
multiprecision library, with a c++ class wrapped around them that allows them
to be used pretty much like normal integers. 
**********************************************************************************/

#include "mpir.h"
#include "boost/multiprecision/gmp.hpp" 

typedef boost::multiprecision::mpz_int Znum;
/* access underlying mpz_t inside an bigint */
#define ZT(a) a.backend().data()

/* convert Znum to long long, checks for overflow. If overflow were to occur,
the function would throw an exception. */
extern long long MulPrToLong(const Znum &x);


/* convert Znum to BigInteger. Returns false if number is too big to convert. 
this function is also used to overload the assignment operator */
bool ZtoBig(BigInteger &number, Znum numberZ) {
	number.nbrLimbs = 0;
	bool neg = false;
	Znum quot, remainder;

	if (numberZ < 0) {
		neg = true;
		numberZ = -numberZ;  // make numberZ +ve
	}
	int i = 0;
	while (numberZ > 0) {
		mpz_fdiv_qr_ui(ZT(quot), ZT(remainder), ZT(numberZ), LIMB_RANGE);
		number.limbs[i].x = (int)MulPrToLong(remainder);
		numberZ = quot;
		i++;
		if (i >= MAX_LEN) {
			return false;   // number too big to convert.
		}
	}
	number.nbrLimbs = i;
	if (neg) {
		number.sign = SIGN_NEGATIVE;
		numberZ = -numberZ;  // put back original value in numberZ
	}
	else
		number.sign = SIGN_POSITIVE;

	return true;
}

/* convert BigInteger to Znum */
void BigtoZ(Znum &numberZ, const BigInteger &number) {

	numberZ = 0;
	for (int i = number.nbrLimbs - 1; i >= 0; i--) {
		numberZ *= LIMB_RANGE;
		numberZ += number.limbs[i].x;
	}
	if (number.sign == SIGN_NEGATIVE)
		numberZ = -numberZ;
}

/* convert integer list to Znum. 1st number in list is the number of entries
in the list. */
void ValuestoZ(Znum &numberZ, const int number[]) {
	numberZ = 0;
	for (int i = number[0]; i >= 1; i--) {
		numberZ *= LIMB_RANGE;
		numberZ += number[i];
	}
}

/* convert number from Znum to BigInteger and factorise it*/
bool factorise(const Znum numberZ, std::vector <Znum> &factorlist,
	std::vector <int> &exponentlist, Znum quads[]) {
	static BigInteger NtoFactor;

	try {
		factorlist.clear();
		exponentlist.clear();
		if (numberZ == 0)
			return false;  // function factor can't factorize zero

		if (ZtoBig(NtoFactor, numberZ)) {
			auto rv = factor(NtoFactor, astFactorsMod);
			if (!rv) return false;
			if (quads != nullptr) {
				ComputeFourSquares(astFactorsMod);  // result is in Quad1-4
				BigtoZ(quads[0], Quad1);
				BigtoZ(quads[1], Quad2);
				BigtoZ(quads[2], Quad3);
				BigtoZ(quads[3], Quad4);
			}
			/*  now convert Factor list from integer lists to Znum list */
			ConvertFactors(astFactorsMod, factorlist, exponentlist);
			return true;
		}
		else
			return false;  // ZtoBig can't convert numberZ to BigInteger. It's too big
	}

	/* code below catches C++ 'throw' type exceptions */
	catch (const std::exception& e) {
		fprintf_s(stderr, "\n*** a standard exception was caught, with message\n '%s'\n", e.what());

		return false;
	}
}
