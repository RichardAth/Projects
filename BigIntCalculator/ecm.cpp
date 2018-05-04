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
#include <iostream>
#include "showtime.h"
#include "bignbr.h"
#include "factor.h"

int yieldFreq;
int EC;            // Elliptic Curve Number
static int limits[] = { 10, 10, 10, 10, 10, 15, 22, 26, 35, 50, 100, 150, 250 };

#define __EMSCRIPTEN__
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

enum eEcmResult
{
	FACTOR_NOT_FOUND = 0,
	FACTOR_FOUND,
	CHANGE_TO_SIQS,
	ERROR
};

static limb A0[MAX_LEN];
static limb A02[MAX_LEN];
static limb A03[MAX_LEN];
static limb AA[MAX_LEN];
static limb DX[MAX_LEN];
static limb DZ[MAX_LEN];
limb GD[MAX_LEN];   // used to pass result back to cller
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

static int indexM, maxIndexM;
static bool foundByLehman;
static bool performLehman;
static int SmallPrime[670] = { 0 }; /* Primes < 5000 */
static int nbrPrimes, indexPrimes, StepECM;
char lowerText[30000];
char *ptrLowerText;
BigInteger Temp1;
static BigInteger Temp2, Temp3, Temp4;

/* function declarations */
static void add3(limb *x3, limb *z3, limb *x2, limb *z2, limb *x1, limb *z1, limb *x, limb *z);
static void duplicate(limb *x2, limb *z2, limb *x1, limb *z1);
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
static void prac(int n, limb *x, limb *z, limb *xT, limb *zT, limb *xT2, limb *zT2)
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
	AddBigNbrModNB(x1, z1, w, TestNbr, NumberLength);      // w = x1+z1
	modmult(v, w, u);       // u = (x2-z2)*(x1+z1)
	AddBigNbrModNB(x2, z2, w, TestNbr, NumberLength);      // w = x2+z2
	SubtBigNbrModN(x1, z1, t, TestNbr, NumberLength); // t = x1-z1
	modmult(t, w, v);       // v = (x2+z2)*(x1-z1)
	AddBigNbrModNB(u, v, t, TestNbr, NumberLength);        // t = 2*(x1*x2-z1*z2)
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
	AddBigNbrModNB(x1, z1, w, TestNbr, NumberLength);      // w = x1+z1
	modmult(w, w, u);       // u = (x1+z1)^2
	SubtBigNbrModN(x1, z1, w, TestNbr, NumberLength); // w = x1-z1
	modmult(w, w, v);       // v = (x1-z1)^2
	modmult(u, v, x2);      // x2 = u*v = (x1^2 - z1^2)^2
	SubtBigNbrModN(u, v, w, TestNbr, NumberLength);   // w = u-v = 4*x1*z1
	modmult(fieldAA, w, u);
	AddBigNbrModNB(u, v, u, TestNbr, NumberLength);        // u = (v+b*w)
	modmult(w, u, z2);      // z2 = (w*u)
}
/* End of code adapted from Paul Zimmermann's ECM4C */

/* compute gcd of value & Global var TestNbr. return 0 if either
is zero, 1 if gcd is 1 , 2 if gcd > 1.
The value of the gcd is returned in GD (also in Temp3)
Uses global variables Temp1, Temp2, Temp3, GD, NumberLength */
static int gcdIsOne(limb *value)
{
	LimbsToBigInteger(value, Temp1, NumberLength);    // Temp1 = value
	LimbsToBigInteger(TestNbr, Temp2, NumberLength);  // Temp2 = TestNbr;
													  // Return zero if value is zero or both numbers are equal.
	if (Temp1 == 0) {
		return 0;
	}
	//Temp3 = Temp1 - Temp2; // BigIntSubt(Temp1, Temp2, Temp3);
	if (Temp1 == Temp2) {
		return 0;
	}
	BigIntGcd(Temp1, Temp2, Temp3);
	BigIntegerToLimbs(GD, Temp3, NumberLength);  // GD = gcd(value, TestNbr)
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
	a = sqr.sqRoot();

	for (;;) {
		if ((a.lldata() & (m - 1)) == r) {
			nextroot = a*a - sqr;
			if (nextroot >= 0) {
				break;
			}
		}
		a++;              // a <- a + 1
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
			sqrRoot += val;
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

#ifdef __EMSCRIPTEN__
void showECMStatus(void) {
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
	strcpy(ptrStatus, lang ? "4\nTranscurri� " : "4\nTime elapsed: ");
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
	printf("%s\n", status);
}

#endif

/* can return value:
FACTOR_NOT_FOUND   (not used??)
CHANGE_TO_SIQS
FACTOR_FOUND
ERROR              (not used??)  */
static enum eEcmResult ecmCurve(BigInteger &N) {
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
			/* copy value of potentialFactor to GD */
			BigIntegerToLimbs(GD, potentialFactor, NumberLength);
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
		strcpy(ptrText, lang ? " usando l�mites B1=" : " using bounds B1=");
		ptrText += strlen(ptrText);
		int2dec(&ptrText, L1);   // Show first bound.
		strcpy(ptrText, lang ? " y B2=" : " and B2=");
		ptrText += strlen(ptrText);
		int2dec(&ptrText, L2);   // Show second bound.
		strcpy(ptrText, "\n");
		ptrText += strlen(ptrText);
		printf("%s", ptrLowerText);
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
		AddBigNbrModNB(A0, Aux2, Aux1, TestNbr, NumberLength); // Aux1 <- A0+2
		modmultInt(MontgomeryMultR1, 4, Aux2);  // Aux2 <- 4
		ModInvBigNbr(Aux2, Aux2, TestNbr, NumberLength);
		modmult(Aux1, Aux2, AA);
		//   X <- (3 * A0 ^ 2 + 1) mod N
		modmultInt(A02, 3, Aux1);    // Aux1 <- 3*A0^2
		AddBigNbrModNB(Aux1, MontgomeryMultR1, X, TestNbr, NumberLength);
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
				{ // If GcdAccumulated is multiple of TestNbr, continue.
					memcpy(X, Xaux, NumberLength * sizeof(limb));
					memcpy(Z, Zaux, NumberLength * sizeof(limb));
					continue; // 
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
			AddBigNbrModNB(X, Z, Aux1, TestNbr, NumberLength);
			modmult(Aux1, Aux1, W1);
			SubtBigNbrModN(X, Z, Aux1, TestNbr, NumberLength);
			modmult(Aux1, Aux1, W2);
			modmult(W1, W2, TX);
			SubtBigNbrModN(W1, W2, Aux1, TestNbr, NumberLength);
			modmult(Aux1, AA, Aux2);
			AddBigNbrModNB(Aux2, W2, Aux3, TestNbr, NumberLength);
			modmult(Aux1, Aux3, TZ); // (TX:TZ) -> 2Q
			SubtBigNbrModN(X, Z, Aux1, TestNbr, NumberLength);
			AddBigNbrModNB(TX, TZ, Aux2, TestNbr, NumberLength);
			modmult(Aux1, Aux2, W1);
			AddBigNbrModNB(X, Z, Aux1, TestNbr, NumberLength);
			SubtBigNbrModN(TX, TZ, Aux2, TestNbr, NumberLength);
			modmult(Aux1, Aux2, W2);
			AddBigNbrModNB(W1, W2, Aux1, TestNbr, NumberLength);
			modmult(Aux1, Aux1, Aux2);
			modmult(Aux2, UZ, X);
			SubtBigNbrModN(W1, W2, Aux1, TestNbr, NumberLength);
			modmult(Aux1, Aux1, Aux2);
			modmult(Aux2, UX, Z); // (X:Z) -> 3Q
			for (I = 5; I < SIEVE_SIZE; I += 2) {
				memcpy(WX, X, NumberLength * sizeof(limb));
				memcpy(WZ, Z, NumberLength * sizeof(limb));
				SubtBigNbrModN(X, Z, Aux1, TestNbr, NumberLength);
				AddBigNbrModNB(TX, TZ, Aux2, TestNbr, NumberLength);
				modmult(Aux1, Aux2, W1);
				AddBigNbrModNB(X, Z, Aux1, TestNbr, NumberLength);
				SubtBigNbrModN(TX, TZ, Aux2, TestNbr, NumberLength);
				modmult(Aux1, Aux2, W2);
				AddBigNbrModNB(W1, W2, Aux1, TestNbr, NumberLength);
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

			AddBigNbrModNB(DX, DZ, Aux1, TestNbr, NumberLength);
			modmult(Aux1, Aux1, W1);
			SubtBigNbrModN(DX, DZ, Aux1, TestNbr, NumberLength);
			modmult(Aux1, Aux1, W2);
			modmult(W1, W2, X);
			SubtBigNbrModN(W1, W2, Aux1, TestNbr, NumberLength);
			modmult(Aux1, AA, Aux2);
			AddBigNbrModNB(Aux2, W2, Aux3, TestNbr, NumberLength);
			modmult(Aux1, Aux3, Z);
			memcpy(UX, X, NumberLength * sizeof(limb));
			memcpy(UZ, Z, NumberLength * sizeof(limb));    // (UX:UZ) -> SIEVE_SIZE*Q
			AddBigNbrModNB(X, Z, Aux1, TestNbr, NumberLength);
			modmult(Aux1, Aux1, W1);
			SubtBigNbrModN(X, Z, Aux1, TestNbr, NumberLength);
			modmult(Aux1, Aux1, W2);
			modmult(W1, W2, TX);
			SubtBigNbrModN(W1, W2, Aux1, TestNbr, NumberLength);
			modmult(Aux1, AA, Aux2);
			AddBigNbrModNB(Aux2, W2, Aux3, TestNbr, NumberLength);
			modmult(Aux1, Aux3, TZ); // (TX:TZ) -> 2*SIEVE_SIZE*Q
			SubtBigNbrModN(X, Z, Aux1, TestNbr, NumberLength);
			AddBigNbrModNB(TX, TZ, Aux2, TestNbr, NumberLength);
			modmult(Aux1, Aux2, W1);
			AddBigNbrModNB(X, Z, Aux1, TestNbr, NumberLength);
			SubtBigNbrModN(TX, TZ, Aux2, TestNbr, NumberLength);
			modmult(Aux1, Aux2, W2);
			AddBigNbrModNB(W1, W2, Aux1, TestNbr, NumberLength);
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
							break;  // This curve cannot factor the number.
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
					AddBigNbrModNB(TX, TZ, Aux2, TestNbr, NumberLength);
					modmult(Aux1, Aux2, W1);
					AddBigNbrModNB(X, Z, Aux1, TestNbr, NumberLength);
					SubtBigNbrModN(TX, TZ, Aux2, TestNbr, NumberLength);
					modmult(Aux1, Aux2, W2);
					AddBigNbrModNB(W1, W2, Aux1, TestNbr, NumberLength);
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

/* returns true if successful. The factor found is returned in global variable GD */
bool ecm(Znum &Nz) {
	static BigInteger N;
	ZtoBig(N, Nz);  // convert from Znum to BigInteger
#ifdef _DEBUG
	//std::cout << "ecm; N = " << Nz << '\n';
#endif
	int P, Q;
#ifndef __EMSCRIPTEN__
	(void)Factors;     // Ignore parameter.
#endif
	fieldTX = TX;
	fieldTZ = TZ;
	fieldUX = UX;
	fieldUZ = UZ;


	fieldAA = AA;
	NumberLength = N.nbrLimbs;
	BigIntegerToLimbs(TestNbr, N, N.nbrLimbs);  // copy N to TestNbr
	GetYieldFrequency();      //get yield frequency (used by showECMStatus)
	GetMontgomeryParms(NumberLength);
	/* set variables to zero */
	memset(M, 0, NumberLength * sizeof(limb));
	memset(DX, 0, NumberLength * sizeof(limb));
	memset(DZ, 0, NumberLength * sizeof(limb));
	memset(W3, 0, NumberLength * sizeof(limb));
	memset(W4, 0, NumberLength * sizeof(limb));
	memset(GD, 0, NumberLength * sizeof(limb));
	ptrLowerText = lowerText;
#ifdef __EMSCRIPTEN__ 
	//	*ptrLowerText++ = '3';
	//if (pstFactors[0].multiplicity > 1)
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
	/* fill SmallPrime array unless it's already set up */
	if (SmallPrime[0] != 2) {
		SmallPrime[0] = 2;   // put 1st prime into SmallPrime list
		P = 3;
		indexM = 1;
		for (indexM = 1; indexM < sizeof(SmallPrime) / sizeof(SmallPrime[0]); indexM++)
		{     // Loop that fills the SmallPrime array.
			SmallPrime[indexM] = (int)P; /* Store prime */
			do {
				P += 2;
				for (Q = 3; Q * Q <= P; Q += 2) { /* Check if P is prime */
					if (P % Q == 0) {
						break;  /* Composite */
					}
				}
			} while (Q * Q <= P);
		}
	}
	foundByLehman = false;
	do {
		enum eEcmResult ecmResp = ecmCurve(N);
		if (ecmResp == CHANGE_TO_SIQS) {    // Perform SIQS
			FactoringSIQSx(TestNbr, GD);
			break;  // value of factor found is in global variable GD
		}
		else if (ecmResp == FACTOR_FOUND) {
			break;  // value of factor found is in global variable GD
		}
		// statements below cannot be executed as ecmResp always = CHANGE_TO_SIQS or FACTOR_FOUND 
		else if (ecmResp == ERROR)
			return false;
	} while (!memcmp(GD, TestNbr, NumberLength * sizeof(limb))); // while GD != TestNbr

#if 0
	lowerTextArea.setText("");
#endif
	StepECM = 0; /* do not show pass number on screen */
	return true;
}