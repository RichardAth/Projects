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

/* factorise numbers or numeric expressions using fast algorithms ECM and SIQS.
ECM = Lenstra elliptic-curve factorization
see https://en.wikipedia.org/wiki/Lenstra_elliptic-curve_factorization 
Profiling indicates that about 2/3 of the CPU time is used during Modular Multiplication. */

#define  _CRT_SECURE_NO_DEPRECATE
#include <iostream>
#include <Windows.h>
#include "showtime.h"
#include "bignbr.h"
#include "factor.h"

extern HANDLE hConsole;
static bool first = true;
static COORD coordScreen = { 0, 0 };    // home for the cursor 
static CONSOLE_SCREEN_BUFFER_INFO csbi;

static int yieldFreq;
int ElipCurvNo;            // Elliptic Curve Number
static int limits[] = { 10, 10, 10, 10, 10, 15, 22, 26, 35, 50, 100, 150, 250 };


#define __EMSCRIPTEN__   // turn on status messages

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
	ECM_ERROR
};

static limb AA[MAX_LEN];
static limb DX[MAX_LEN];
static limb DZ[MAX_LEN];
static limb M[MAX_LEN];
static limb TX[MAX_LEN];
static limb TZ[MAX_LEN];
static limb UX[MAX_LEN];
static limb UZ[MAX_LEN];
static limb W3[MAX_LEN];
static limb W4[MAX_LEN];
static limb Aux1[MAX_LEN];
static limb Aux2[MAX_LEN];
static limb Aux3[MAX_LEN];
static limb Aux4[MAX_LEN];

static unsigned char sieve[10 * SIEVE_SIZE];
static unsigned char sieve2310[SIEVE_SIZE];
static int sieveidx[GROUP_SIZE];
static limb GcdAccumulated[MAX_LEN];

static int indexM, maxIndexM;
static bool foundByLehman;

static int SmallPrime[670] = { 0 }; /* Primes < 5000 */
static int nbrPrimes, indexPrimes, StepECM;
char lowerText[30000];
char *ptrLowerText;
static Znum BiGD;   // used by function gcdIsOne to return result
/* forward function declarations */
static void add3(limb *x3, limb *z3, const limb *x2, const limb *z2, 
	const limb *x1, const limb *z1, const limb *x, const limb *z);
static void duplicate(limb *x2, limb *z2, const limb *x1, const limb *z1);

/******************************************************/
/* Start of code adapted from Paul Zimmermann's ECM4C */
/******************************************************/
#define ADD 6  /* number of multiplications in an addition */
#define DUP 5  /* number of multiplications in a duplicate */

/* uses global var NumberLength, returns value in YieldFrequency */
static void GetYieldFrequency(void)
{  // changed constants to allow for change in NumberLength
	yieldFreq =250000 / (NumberLength * NumberLength) + 1;  
	/* round yieldFreq down to a power of 10 */
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
	int d, e, r, LucCostI;
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
	LucCostI = 0;

	for (d = 1; d < 10; d++) {
		e = lucas_cost(n, v[d]);
		if (e < r) {
			r = e;
			LucCostI = d;
		}
	}
	d = n;
	r = (int)((double)d / v[LucCostI] + 0.5);
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
			r = d; 	d = e; 	e = r;     // rotate r, d and e values
			t = xA;  xA = xB; xB = t;  // swap Xa and Xb (pointers)
			t = zA; zA = zB; zB = t;   // swap Za and Zb (pointers)
		}
		/* do the first line of Table 4 whose condition qualifies */
		if (4 * d <= 5 * e && ((d + e) % 3) == 0) { /* condition 1 */
			r = (2 * d - e) / 3;
			e = (2 * e - d) / 3;
			d = r;
			add3(xT, zT, xA, zA, xB, zB, xC, zC); /* T = f(A,B,C) */
			add3(xT2, zT2, xT, zT, xA, zA, xB, zB); /* T2 = f(T,A,B) */
			add3(xB, zB, xB, zB, xT, zT, xA, zA); /* B = f(B,T,A) */
			t = xA; xA = xT2; xT2 = t;  // swap xA and xT2 (pointers)
			t = zA; zA = zT2; 	zT2 = t; /* swap A and T2 */
		}
		else if (4 * d <= 5 * e && (d - e) % 6 == 0) { /* condition 2 */
			d = (d - e) / 2;
			add3(xB, zB, xA, zA, xB, zB, xC, zC); /* B = f(A,B,C) */
			duplicate(xA, zA, xA, zA); /* A = 2*A */
		}
		else if (d <= (4 * e)) { /* condition 3 */
			d -= e;
			add3(xT, zT, xB, zB, xA, zA, xC, zC); /* T = f(B,A,C) */
			t = xB; xB = xT;  xT = xC; xC = t;  // rotate B, T & C (pointers)
			t = zB; zB = zT; zT = zC; zC = t; /* circular permutation (B,T,C) */
		}
		else if ((d + e) % 2 == 0) { /* condition 4 */
			d = (d - e) / 2;
			add3(xB, zB, xB, zB, xA, zA, xC, zC);   /* B = f(B,A,C) */
			duplicate(xA, zA, xA, zA);              /* A = 2*A */
		}
		else if (d % 2 == 0) { /* condition 5 */
			d /= 2;
			add3(xC, zC, xC, zC, xA, zA, xB, zB);   /* C = f(C,A,B) */
			duplicate(xA, zA, xA, zA);              /* A = 2*A */
		}
		else if (d % 3 == 0) { /* condition 6 */
			d = d / 3 - e;
			duplicate(xT, zT, xA, zA); /* T1 = 2*A */
			add3(xT2, zT2, xA, zA, xB, zB, xC, zC); /* T2 = f(A,B,C) */
			add3(xA, zA, xT, zT, xA, zA, xA, zA);   /* A = f(T1,A,A) */
			add3(xT, zT, xT, zT, xT2, zT2, xC, zC); /* T1 = f(T1,T2,C) */
			t = xC; xC = xB; xB = xT; xT = t;     // rotate C, B & T
			t = zC; zC = zB; zB = zT; zT = t;   /* circular permutation (C,B,T) */
		}
		else if ((d + e) % 3 == 0) { /* condition 7 */
			d = (d - 2 * e) / 3;
			add3(xT, zT, xA, zA, xB, zB, xC, zC); /* T1 = f(A,B,C) */
			add3(xB, zB, xT, zT, xA, zA, xB, zB); /* B = f(T1,A,B) */
			duplicate(xT, zT, xA, zA);
			add3(xA, zA, xA, zA, xT, zT, xA, zA);   /* A = 3*A */
		}
		else if ((d - e) % 3 == 0) { /* condition 8 */
			d = (d - e) / 3;
			add3(xT, zT, xA, zA, xB, zB, xC, zC); /* T1 = f(A,B,C) */
			add3(xC, zC, xC, zC, xA, zA, xB, zB); /* C = f(A,C,B) */
			t = xB; xB = xT; xT = t;   
			t = zB; zB = zT; zT = t; /* swap B and T */
			duplicate(xT, zT, xA, zA);
			add3(xA, zA, xA, zA, xT, zT, xA, zA); /* A = 3*A */
		}
		else if (e % 2 == 0) { /* condition 9 */
			e /= 2;
			add3(xC, zC, xC, zC, xB, zB, xA, zA); /* C = f(C,B,A) */
			duplicate(xB, zB, xB, zB);           /* B = 2*B */
		}
	}
	add3(x, z, xA, zA, xB, zB, xC, zC);
}

/* adds Q=(x2:z2) and R=(x1:z1) and puts the result in (x3:z3),
using 5/6 mul, 6 add/sub and 6 mod. One assumes that Q-R=P or R-Q=P where P=(x:z).
Uses the following global variables:
- n : number to factor = TestNbr
- x, z : coordinates of P
- TX, TZ, UX, UZ : auxiliary (global) variables
Modifies: x3, z3, TX, TZ, UX, UZ.
(x3,z3) may be identical to (x2,z2) and to (x,z)
*/
static void add3(limb *x3, limb *z3, const limb *x2, const limb *z2,
	const limb *x1, const limb *z1, const limb *x, const limb *z)
{
	/* N.B. all arithmetic is mod TestNbr */
	SubtBigNbrModN(x2, z2, UX, TestNbr, NumberLength);   // UX = x2-z2
	AddBigNbrModNB(x1, z1,  UZ, TestNbr, NumberLength);   //  UZ = x1+z1
	modmult(UX,  UZ, TZ);                                  // TZ = (x2-z2)*(x1+z1)
	AddBigNbrModNB(x2, z2,  UZ, TestNbr, NumberLength);   //  UZ = x2+z2
	SubtBigNbrModN(x1, z1, TX, TestNbr, NumberLength);  // TX = x1-z1
	modmult(TX,  UZ, UX);                                  // UX = (x2+z2)*(x1-z1)
	AddBigNbrModNB(TZ, UX, TX, TestNbr, NumberLength);   // TX = 2*(x1*x2-z1*z2)
	modmult(TX, TX,  UZ);                                 //  UZ = 4*(x1*x2-z1*z2)^2
	SubtBigNbrModN(TZ, UX, TX, TestNbr, NumberLength);   // TX = 2*(x2*z1-x1*z2)
	modmult(TX, TX, UX);                                 // UX = 4*(x2*z1-x1*z2)^2
	if (!memcmp(x, x3, NumberLength * sizeof(limb))) {  // if x == x3
		memcpy(TZ, x, NumberLength * sizeof(limb));       // TZ = x
		memcpy(TX,  UZ, NumberLength * sizeof(limb));     // TX = UZ = 4*(x1*x2-z1*z2)^2
		modmult(z, TX,  UZ);                             // UZ = z*TX
		modmult(UX, TZ, z3);                             // z3 = UX*TZ
		memcpy(x3,  UZ, NumberLength * sizeof(limb));     // x3 = UZ
	}
	else {
		modmult( UZ, z, x3);            // x3 = 4*z*(x1*x2-z1*z2)^2
		modmult(x, UX, z3);             // z3 = 4*x*(x2*z1-x1*z2)^2
	}
}

/* computes 2P=(x2:z2) from P=(x1:z1), with 5 mul, 4 add/sub, 5 mod.
Uses the following global variables:
- n : number to factor
- b : (a+2)/4 mod n
- UZ,  TX,  TZ : auxiliary variables
Modifies: x2, z2, UZ,  TX,  TZ
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
static void duplicate(limb *x2, limb *z2, const limb *x1, const limb *z1)
{
	/* N.B. all arithmetic is mod TestNbr */
	AddBigNbrModNB(x1, z1, TZ, TestNbr, NumberLength);    //  TZ = x1+z1 (mod testNbr)
	modmult(TZ, TZ, UZ);                                  //  UZ = (x1+z1)^2 (mod testNbr)
	SubtBigNbrModN(x1, z1, TZ, TestNbr, NumberLength);    //  TZ = x1-z1 (mod testNbr)
	modmult(TZ, TZ,  TX);                                 //  TX = (x1-z1)^2 (mod testNbr)
	modmult(UZ, TX, x2);                                  // x2 = UZ* TX (mod testNbr)
	SubtBigNbrModN(UZ, TX, TZ, TestNbr, NumberLength);    //  TZ = UZ- TX = 4*x1*z1 (mod testNbr)
	modmult(AA, TZ, UZ);                                  // UZ =  TZ*AA (mod testNbr)
	AddBigNbrModNB(UZ, TX, UZ, TestNbr, NumberLength);    // UZ = ( TX+b* TZ) (mod testNbr)
	modmult(TZ, UZ, z2);                                   // z2 = ( TZ*u) (mod testNbr)
}
/* End of code adapted from Paul Zimmermann's ECM4C */

/* compute gcd of value & N. return 0 if value is zero or value = N, 
1 if gcd is 1 , 2 if gcd > 1.
The value of the gcd is returned in BiGD
Uses global variables Temp1, BiGD, NumberLength */
static int gcdIsOne(const limb *value, const Znum &zN) {
	static Znum Temp1;

	LimbstoZ(value, Temp1, NumberLength);    // Temp1 = value

    // Return zero if value is zero or both numbers are equal.
	if (Temp1 == 0) {
		return 0;
	}
	if (Temp1 == zN) {
		return 0;
	}
	BiGD = gcd(Temp1, zN);     // BiGD = gcd(value, N)
	if (BiGD < 2) {
		return (int)MulPrToLong(BiGD);    // GCD is less than 2.
	}
	return 2;      // GCD is greater than one.
}

static void GenerateSieve(int initial) {
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
	if (elapsedTime / 10 <= oldTimeElapsed/10 +5)
	{
		return;  // less than 5 seconds since last message
	}
	oldTimeElapsed = elapsedTime;
	ptrStatus = status;
	strcpy(ptrStatus, lang ? "Transcurrió " : "Time elapsed: ");
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
	*ptrStatus++ = '\0';  // add null terminator
	printf("%s\n", status);
}

#endif

/* Perform (modified) Lehman algorithm. Try to factorise nbr.
This produces a result when nbr is a product of two factors which are not
too far apart in value. e.g larger factor is less than k times the smaller.
The factors produced may not be primes. 
We need to find integers k, a and b such that 4kn = a^2 – b^2.
It can be shown that 1 <= k <= n^(1/3) and
sqrt(4kn) <= a <= sqrt(4kn) + (n^(1/6)/(4.sqrt(k)) 
However, searching the whole range of possible values would be much too slow. */
static void LehmanZ(const Znum &nbr, int k, Znum &factor) {
	const long long bitsSqr[] = {
		0x0000000000000003, // 3
		0x0000000000000013, // 5
		0x0000000000000017, // 7
		0x000000000000023B, // 11
		0x000000000000161B, // 13
		0x000000000001A317, // 17
		0x0000000000030AF3, // 19
		0x000000000005335F, // 23
		0x0000000013D122F3, // 29
		0x00000000121D47B7, // 31
		0x000000165E211E9B, // 37
		0x000001B382B50737, // 41  
		0x0000035883A3EE53, // 43  
		0x000004351B2753DF, // 47
		0x0012DD703303AED3, // 53
		0x022B62183E7B92BB, // 59
		0x1713E6940A59F23B, // 61
	};
	int primes[] = { 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61 };
	int nbrs[17];
	int diffs[17];
	int i, j, m, r;
	Znum sqrRoot, nextroot;
	Znum a, c, sqr, val;

	if (ZisEven(nbr))
	{ // nbr is Even (this would be strange, because small factors such as 2 are removed during trial division)
		r = 0;
		m = 1;
	}
	else {
		if (k % 2 == 0) { // k Even
			r = 1;
			m = 2;
		}
		else { // k Odd
			r = (int)MulPrToLong((k + nbr)& 3) ;
			m = 4;
		}
	}

	sqr = k << 2;      // sqr = 4k
	sqr *= nbr;      // sqr = 4kn
	squareRoot(sqr, sqrRoot);      // sqrRoot = sqrt(4kn)
	a = sqrRoot;                 // a = sqrt(4kn)

	for (;;) {
		if ((a & (m - 1)) == r) {
			nextroot = a*a;
			nextroot -= sqr;
			if (nextroot >= 0)
			{
				break;  // a*a >= 4kn
			}
		}
		a++;
	}

	nextroot = a*a;
	c = nextroot - sqr;  // c = a*a- 4kn    

	for (i = 0; i < 17; i++) {
		int pr = primes[i];
		nbrs[i] = (int)MulPrToLong(c % primes[i]);
		diffs[i] = m * (int)MulPrToLong(((a %pr) * 2 + m) % pr);
	}

	for (j = 0; j < 10000; j++) {
		for (i = 0; i < 17; i++) {
			int shiftBits = nbrs[i];
			if ((bitsSqr[i] & (1LL << shiftBits)) == 0)
				break;   // Not a perfect square
		}

		if (i == 17) {            // Test for perfect square
			c = m * j;
			val = a + c;     // a + m*j       
			c = val * val;   // (a + m*j)^2       
			c -= sqr;        //  (a + m*j)^2 - 4kn  
			squareRoot(c, sqrRoot);       // sqrRoot <- sqrt(c)
			sqrRoot += val;
			c = gcd(sqrRoot, nbr); // Get GCD(sqrRoot + val, nbr)
			if (c > 1) {       // factor has been found.
				factor = c;
				return;
			}
		}

		for (i = 0; i < 17; i++)
		{
			nbrs[i] = (nbrs[i] + diffs[i]) % primes[i];
			diffs[i] = (diffs[i] + 2 * m * m) % primes[i];
		}
	}

	/* failed to factorise in 1000 cycles, so stop trying */
	factor = 1;   // intToBigInteger(factor, 1);   // Factor not found.
}

static void upOneLine(void) {
	SetConsoleCursorPosition(hConsole, coordScreen);
}

/* can return value:
FACTOR_NOT_FOUND   (not used??)
CHANGE_TO_SIQS
FACTOR_FOUND   - value of factor is returned in global variable BiGD
ERROR              (not used??)  */
static enum eEcmResult ecmCurve(const Znum &zN, Znum &Zfactor) {
	static limb A0[MAX_LEN];
	static limb A02[MAX_LEN];
	static limb A03[MAX_LEN];

	static limb W1[MAX_LEN];
	static limb W2[MAX_LEN];
	static limb WX[MAX_LEN];
	static limb WZ[MAX_LEN];
	static limb X[MAX_LEN];
	static limb Z[MAX_LEN];
	static limb Xaux[MAX_LEN];
	static limb Zaux[MAX_LEN];
	static limb root[GROUP_SIZE][MAX_LEN];

#ifdef __EMSCRIPTEN__
	//char text[20];
#endif

	for (;;) {
#ifdef __EMSCRIPTEN__
		char *ptrText;
#endif
		int I, Pass;
		int i, j, u;
		long long L1, L2, LS, P, IP, Paux = 1;

		ElipCurvNo++;   // increment curve number

#ifdef __EMSCRIPTEN__
		//			text[0] = '7';
		//ptrText = &text[1];
		//			int2dec(&ptrText, ElipCurvNo);
		//			printf ("%s\n", text);
#endif
		L2 = mpz_sizeinbase(ZT(zN), 10);   // Get number of digits.
		if (L2 > 30 && L2 <= 90)          // If between 30 and 90 digits...
		{                                 // switch to SIQS when curve No reaches limit
			int limit = limits[((int)L2 - 26) / 5];  // e.g if L2<=50, limit=10
			if (ElipCurvNo  >= limit) {                          
   				return CHANGE_TO_SIQS;           // Switch to SIQS.
			}
		}

		/* Try to factor N using Lehman algorithm. Result in Zfactor. 
		This seldom achieves anything, but when it does it saves a lot of time */
		int kx = ElipCurvNo;
		const int mult = 5;  /* could change value of mult to use Lehman more, 
							 or less, relative to ECM. Benchmark testing
							 needed to estimate best value */
		for (int k = (kx - 1)*mult + 1; k <= kx*mult; k++) {
			/* for curve no 1, k = 1 to 5, curve no 2, k = 6 to 10, etc*/
			LehmanZ(zN, k, Zfactor);
			if (Zfactor > LIMB_RANGE) {
				foundByLehman = true;     // Factor found.
//#ifdef _DEBUG
				std::cout << "Lehman factor found. k = " << k << " N= " << zN 
					<< " factor = " << Zfactor << '\n';
//#endif
				return FACTOR_FOUND;
			}
		}

		/* set L1, L2, LS, Paux and nbrPrimes according to value of ElipCurvNo */
		L1 = 2000;
		L2 = 200000;
		LS = 45;
		Paux = ElipCurvNo;
		nbrPrimes = 303; /* Number of primes less than 2000 */
		if (ElipCurvNo > 25) {
			if (ElipCurvNo < 326) {   // 26 to 325
				L1 = 50000;
				L2 = 5000000;
				LS = 224;
				Paux = ElipCurvNo - 24;
				nbrPrimes = 5133; /* Number of primes less than 50000 */
			}
			else {
				if (ElipCurvNo < 2000) {  // 326 to 1999
					L1 = 1000000;
					L2 = 100000000;
					LS = 1001;
					Paux = ElipCurvNo - 299;
					nbrPrimes = 78498; /* Number of primes less than 1000000 */
				}
				else {   // >= 2000
					L1 = 11000000;
					L2 = 1100000000;
					LS = 3316;
					Paux = ElipCurvNo - 1900;
					nbrPrimes = 726517; /* Number of primes less than 11000000 */
				}
			}
		}
#ifdef __EMSCRIPTEN__
		/* print status message */
		ptrText = ptrLowerText;  // Point after number that is being factored.
		auto elapsedTime = (int)(tenths() - originalTenthSecond);
		GetDHMSt(&ptrText, elapsedTime);
		strcpy(ptrText, lang ? " ECM Curva " : " ECM Curve ");
		ptrText += strlen(ptrText);
		int2dec(&ptrText, ElipCurvNo);   // Show curve number.
		strcpy(ptrText, lang ? " usando límites B1=" : " using bounds B1=");
		ptrText += strlen(ptrText);
		int2dec(&ptrText, L1);   // Show first bound.
		strcpy(ptrText, lang ? " y B2=" : " and B2=");
		ptrText += strlen(ptrText);
		int2dec(&ptrText, L2);   // Show second bound.
		strcpy(ptrText, "\n");
		ptrText += strlen(ptrText);

		if (first) {
			if (!GetConsoleScreenBufferInfo(hConsole, &csbi))
			{
				fprintf(stderr, "** GetConsoleScreenBufferInfo failed with %d!\n", GetLastError());
				Beep(750, 1000);
			}
			coordScreen.X = csbi.dwCursorPosition.X;  // save cursor co-ordinates
			coordScreen.Y = csbi.dwCursorPosition.Y;
			if (csbi.dwCursorPosition.Y >= csbi.dwSize.Y - 1)
				coordScreen.Y--;  // if window is full, allow for text scrolling up
		}
		else
			upOneLine();

		printf("%s", ptrLowerText);
		first = false;
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
			primalityString + ElipCurvNo + "\n" + UpperLine + "\n" + LowerLine);
#endif
#endif

		//  Compute A0 <- 2 * (ElipCurvNo+1)/(3 * (ElipCurvNo+1) ^ 2 - 1, N) mod N
		// Aux2 <- 1 in Montgomery notation.
		memcpy(Aux2, MontgomeryMultR1, NumberLength * sizeof(limb));
		modmultInt(Aux2, ElipCurvNo + 1, Aux2);    // Aux2 <- ElipCurvNo + 1 (mod TestNbr)
		modmultInt(Aux2, 2, Aux1);                 // Aux1 <- 2*(ElipCurvNo+1) (mod TestNbr)
		modmultInt(Aux2, ElipCurvNo + 1, Aux3);    // Aux3 <- (ElipCurvNo + 1)^2 (mod TestNbr)
		modmultInt(Aux3, 3, Aux3);                 // Aux3 <- 3*(ElipCurvNo + 1)^2 (mod TestNbr)
		// Aux2 <- 3*(ElipCurvNo + 1)^2 - 1 (mod TestNbr)
		SubtBigNbrModN(Aux3, MontgomeryMultR1, Aux2, TestNbr, NumberLength);
		ModInvBigNbr(Aux2, Aux2, TestNbr, NumberLength);
		modmult(Aux1, Aux2, A0);       // A0 <- 2*(ElipCurvNo+1)/(3*(ElipCurvNo+1)^2 - 1) (mod TestNbr)
#ifdef _DEBUG
		{
			Znum temp;
			LimbstoZ(A0, temp, NumberLength);
			REDC(temp, temp);
			std::cout << "curve no = " << ElipCurvNo << " length = " << NumberLength << " A0 = "
				<< temp << '\n';
		}
#endif

		//  if A0*(A0 ^ 2 - 1)*(9 * A0 ^ 2 - 1) mod N=0 then select another curve.
		modmult(A0, A0, A02);          // A02 <- A0^2
		modmult(A02, A0, A03);         // A03 <- A0^3
		SubtBigNbrModN(A03, A0, Aux1, TestNbr, NumberLength);  // Aux1 <- A0^3 - A0
		modmultInt(A02, 9, Aux2);      // Aux2 <- 9*A0^2
		SubtBigNbrModN(Aux2, MontgomeryMultR1, Aux2, TestNbr, NumberLength); // Aux2 <- 9*A0^2-1
		modmult(Aux1, Aux2, Aux3);
		if (BigNbrIsZero(Aux3, NumberLength)) {
			continue;  // select another curve
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
		memcpy(Xaux, X, NumberLength * sizeof(limb));   // Xaux = X
		memcpy(Zaux, Z, NumberLength * sizeof(limb));   // Zaux = Z
		// GcdAccumulated = 1
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
				modmult(GcdAccumulated, Z, Aux1);    // GcdAccumulated *= Z
				memcpy(GcdAccumulated, Aux1, NumberLength * sizeof(limb));
			}
			else {
				if (gcdIsOne(Z, zN) > 1) {
					//BigtoZ(Zfactor, BiGD); // copy factor to global Znum
					Zfactor = BiGD;
					return FACTOR_FOUND;     // ** found a factor in global variable BiGD**
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
					modmult(GcdAccumulated, Z, Aux1);     // GcdAccumulated *= Z
					memcpy(GcdAccumulated, Aux1, NumberLength * sizeof(limb));
				}
				else {
					if (gcdIsOne(Z, zN) > 1) {
						//BigtoZ(Zfactor, BiGD); // copy factor to global Znum
						Zfactor = BiGD;
						return FACTOR_FOUND;     // ** found a factor in global variable BiGD**
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
						modmult(GcdAccumulated, Z, Aux1);    // GcdAccumulated *= Z
						memcpy(GcdAccumulated, Aux1, NumberLength * sizeof(limb));
					}
					else {
						if (gcdIsOne(Z, zN) > 1) {
							//BigtoZ(Zfactor, BiGD); // copy factor to global Znum
							Zfactor = BiGD;
							return FACTOR_FOUND;   // ** found a factor in global variable BiGD**
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
				if (gcdIsOne(GcdAccumulated, zN) > 1) {
					//BigtoZ(Zfactor, BiGD); // copy factor to global Znum
					Zfactor = BiGD;
					return FACTOR_FOUND;   // ** found a factor **  in global variable BiGD
				}
				break;  // exit for loop
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
			AddBigNbrModNB(X, Z, Aux1, TestNbr, NumberLength);  // Aux1 = X+Z (mod TestNbr)

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
			modmult(Aux2, UX, Z);       // (X:Z) -> 3Q
			for (I = 5; I < SIEVE_SIZE; I += 2) {
				memcpy(WX, X, NumberLength * sizeof(limb));            // WX = X
				memcpy(WZ, Z, NumberLength * sizeof(limb));            // WZ = Z
				SubtBigNbrModN(X, Z, Aux1, TestNbr, NumberLength);     // Aux1 = X - Z
				AddBigNbrModNB(TX, TZ, Aux2, TestNbr, NumberLength);   // Aux2 = TX - TZ
				modmult(Aux1, Aux2, W1);                               // W1   = Aux1 * Aux2
				AddBigNbrModNB(X, Z, Aux1, TestNbr, NumberLength);     // Aux1 = X + Z
				SubtBigNbrModN(TX, TZ, Aux2, TestNbr, NumberLength);   // Aux2 = TX - TZ
				modmult(Aux1, Aux2, W2);                               // W2   = Aux1 * Aux2
				AddBigNbrModNB(W1, W2, Aux1, TestNbr, NumberLength);
				modmult(Aux1, Aux1, Aux2);                             // Aux2 = Aux1 * Aux1
				modmult(Aux2, UZ, X);                                  // X    = Aux2 * UZ
				SubtBigNbrModN(W1, W2, Aux1, TestNbr, NumberLength);   // Aux1 = W1-W2
				modmult(Aux1, Aux1, Aux2);
				modmult(Aux2, UX, Z); // (X:Z) -> 5Q, 7Q, ...
				if (Pass == 0) {
					modmult(GcdAccumulated, Aux1, Aux2);
					memcpy(GcdAccumulated, Aux2, NumberLength * sizeof(limb));
				}
				else {
					if (gcdIsOne(Aux1, zN) > 1) {
						//BigtoZ(Zfactor, BiGD); // copy factor to global Znum
						Zfactor = BiGD;
						return FACTOR_FOUND;   // ** found a factor in global variable BiGD
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
						if (gcdIsOne(GcdAccumulated, zN) > 1) {
							//BigtoZ(Zfactor, BiGD); // copy factor to global Znum
							Zfactor = BiGD;
							return FACTOR_FOUND;    // ** found a factor in global variable BiGD
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
				rc = gcdIsOne(GcdAccumulated, zN);
				if (rc == 1) {
					break;    // GCD is one, so this curve does not find a factor.
				}
				if (rc == 0) {
					continue;  // GcdAccumulated = N or 0
				}

				//if (memcmp(BNgcd, TestNbr, NumberLength * sizeof(limb)))
				if (BiGD != zN)
				{           // GCD is not 1 or TestNbr
					//BigtoZ(Zfactor, BiGD); // copy factor to global Znum
					Zfactor = BiGD;
					return FACTOR_FOUND;     // ** found a factor in global variable BiGD
				}
			}
		} /* end for Pass */

	}       /* End curve calculation */
}

/* initialise variables. zN = number to be factored */
static void ecminit(Znum zN) {
	int P, Q;

	// calculate number of limbs
	//NumberLength = (int)(mpz_sizeinbase(ZT(zN), 2) + BITS_PER_GROUP - 1) / BITS_PER_GROUP;
	ZtoBig(TestNbrBI, zN);   /* copy zN to TestNbr. NB throw exception
				             if zN is too large! (more than about 23,000 digits) */
	GetYieldFrequency();      //get yield frequency (used by showECMStatus)
	GetMontgomeryParms(NumberLength);
#ifdef _DEBUG
	GetMontgomeryParms(zN);
#endif
	first = true;

	/* set variables to zero */
	memset(M, 0, NumberLength * sizeof(limb));
	memset(DX, 0, NumberLength * sizeof(limb));
	memset(DZ, 0, NumberLength * sizeof(limb));
	memset(W3, 0, NumberLength * sizeof(limb));
	memset(W4, 0, NumberLength * sizeof(limb));
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
	ElipCurvNo--;  // decrement curve number
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
}

/* returns true if successful. The factor found is returned in global Znum Zfactor */
bool ecm(Znum &zN, long long maxdivisor) {

#ifdef _DEBUG
	//std::cout << "ecm; N = " << zN << '\n';
#endif

	ecminit(zN);  // initialise values

	foundByLehman = false;
	do {
		enum eEcmResult ecmResp = ecmCurve(zN, Zfactor);
		if (ecmResp == CHANGE_TO_SIQS) {    // Perform SIQS
			FactoringSIQS(zN, Zfactor); // factor found is returned in Zfactor2
			break;
		}
		else if (ecmResp == FACTOR_FOUND) {
			break;
		}
		// statements below cannot be executed as ecmResp always = CHANGE_TO_SIQS or FACTOR_FOUND 
		else if (ecmResp == ERROR)
			return false;
	} while (zN == BiGD);

#if 0
	lowerTextArea.setText("");
#endif
	StepECM = 0; /* do not show pass number on screen */
	return true;
}