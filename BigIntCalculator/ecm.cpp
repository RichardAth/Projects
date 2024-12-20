﻿/*
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

#include "pch.h"

//#define log 1              // remove this line to generate code without logging

extern HANDLE hConsole;   /* STD_OUTPUT_HANDLE */
typedef void(*mmCback)(void);
extern mmCback modmultCallback;   // function pointer
extern void FactoringSIQS(const Znum& NbrToFactor, Znum& Factor);


static bool first = true;
static COORD coordScreen = { 0, 0 };    // home for the cursor 
static CONSOLE_SCREEN_BUFFER_INFO csbi;

static int yieldFreq;   // used to control output of status messages 
static int ElipCurvNo;            // Elliptic Curve Number


/*
Threshold for switching to SIQS 
Digits	26-50   51-55   56-60	61-65	66-70	71-75	76-80	81-85	86-90	91-95   96-100
Curve	10		15		22		26		35		50		100		150		250     300     350 */
static int limits[] = { 10, 10, 10, 10, 10, 15, 22, 26, 35, 50, 100, 150, 250, 300, 350 };


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

static limb AA[MAX_LEN];               //   AA <- (A + 2)*modinv(4, N) mod N
static limb DX[MAX_LEN];               // (DX:DZ) -> HALF_SIEVE_SIZE*Q
static limb DZ[MAX_LEN];
static limb M[MAX_LEN];               // used by ecmCurve
static limb TX[MAX_LEN];              // used by add3, ecmCurve
static limb TZ[MAX_LEN];              // used by add3, ecmCurve
static limb UX[MAX_LEN];              // used by add3, ecmCurve
static limb UZ[MAX_LEN];              // used by add3, ecmCurve
static limb W3[MAX_LEN];              // used by ecmCurve  ecminit
static limb W4[MAX_LEN];              // used by ecmCurve, ecminit
static limb Aux1[MAX_LEN], Aux2[MAX_LEN], Aux3[MAX_LEN], Aux4[MAX_LEN];      // work variables

static unsigned char sieve[10 * SIEVE_SIZE];
static unsigned char sieve2310[SIEVE_SIZE];
static int sieveidx[GROUP_SIZE];
static limb GcdAccumulated[MAX_LEN];

static int indexM, maxIndexM;
static bool foundByLehman;

static int SmallPrime[670] = { 0 }; /* Primes < 5000 */
static int nbrPrimes, indexPrimes, StepECM=0;
char lowerText[30000];
char *ptrLowerText;
static Znum BiGD;   // used by function gcdIsOne to return result
/* forward function declarations */
static void add3(limb *x3, limb *z3, const limb *x2, const limb *z2, 
	const limb *x1, const limb *z1, const limb *x, const limb *z);
static void duplicate(limb *x2, limb *z2, const limb *x1, const limb *z1);

#ifdef log
static char logbuffer[2000];     // big enough for numbers up to about 190 limbs
//static limb temp1[MAX_LEN], temp2[MAX_LEN]; 
FILE *logfile;

//void print(limb *w) {
//	Bin2Dec(w, logbuffer, NumberLength, 0);
//	fprintf_s(logfile, "%s\n", logbuffer);
//	memset(temp1, 0, NumberLength * sizeof(limb));
//	temp1[0].x = 1;
//	modmult(temp1, w, temp2);
//	Bin2Dec(temp2, logbuffer, NumberLength, 0);
//	fprintf_s(logfile, " (= %s)\n", logbuffer);
//}


#define logf(a) {														  \
	Bin2Dec(##a, logbuffer, NumberLength, -6);							  \
	fprintf_s (logfile, #a " = %s at line %d \n", logbuffer, __LINE__);  \
}
#endif 

/******************************************************/
/* Start of code adapted from Paul Zimmermann's ECM4C */
/******************************************************/
#define ADD 6  /* number of multiplications in an addition */
#define DUP 5  /* number of multiplications in a duplicate */

/* returns value in yieldFreq, based on NumLen which is the number of limbs
in TestNbr.  yieldFreq is used to control output of status messages */
static void GetYieldFrequency(int NumLen)
{
	yieldFreq = 1000000 / (NumLen * NumLen) + 1;
	if (yieldFreq > 100000) {
		yieldFreq = yieldFreq / 100000 * 100000;   // round down to multiple of 100000
	}
	else if (yieldFreq > 10000) {
		yieldFreq = yieldFreq / 10000 * 10000;     // round down to multiple of 10000
	}
	else if (yieldFreq > 1000) {
		yieldFreq = yieldFreq / 1000 * 1000;       // round down to multiple of 1000
	}
	else if (yieldFreq > 100) {
		yieldFreq = yieldFreq / 100 * 100;         // round down to multiple of 100
	}
}

/* returns the number of modular multiplications */
static int lucas_cost(int n, double v)
{
	int cost, d, e, r;

	d = n;
	r = (int)((double)d / v + 0.5);
	if (r >= n) {
		return (ADD * n);
	}
	d = n - r;
	e = 2 * r - n;
	cost = DUP + ADD; /* initial duplicate and final addition */
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
			cost += 3 * ADD; /* 3 additions */
		}
		else if (4 * d <= 5 * e && (d - e) % 6 == 0) { /* condition 2 */
			d = (d - e) / 2;
			cost += ADD + DUP; /* one addition, one duplicate */
		}
		else if (d <= (4 * e)) { /* condition 3 */
			d -= e;
			cost += ADD; /* one addition */
		}
		else if ((d + e) % 2 == 0) { /* condition 4 */
			d = (d - e) / 2;
			cost += ADD + DUP; /* one addition, one duplicate */
		}
		else if (d % 2 == 0) { /* condition 5 */
			d /= 2;
			cost += ADD + DUP; /* one addition, one duplicate */
		}
		else if (d % 3 == 0) { /* condition 6 */
			d = d / 3 - e;
			cost += 3 * ADD + DUP; /* three additions, one duplicate */
		}
		else if ((d + e) % 3 == 0) { /* condition 7 */
			d = (d - 2 * e) / 3;
			cost += 3 * ADD + DUP; /* three additions, one duplicate */
		}
		else if ((d - e) % 3 == 0) { /* condition 8 */
			d = (d - e) / 3;
			cost += 3 * ADD + DUP; /* three additions, one duplicate */
		}
		else if (e % 2 == 0) { /* condition 9 */
			e /= 2;
			cost += ADD + DUP; /* one addition, one duplicate */
		}
	}
	return cost;
}

/* computes nP from P=(x:z) and puts the result in (x:z). Assumes n>2.
uses global variables Aux1, Aux2, Aux3, Aux4 */
static void prac(int n, limb *x, limb *z, limb *xT, limb *zT, limb *xT2, limb *zT2)
{
	int d, e, r, i;
	limb *temp;
	limb *xA = x,    *zA = z;
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
	std::memcpy(xB, xA, NumberLength * sizeof(limb));   // B <- A
	std::memcpy(zB, zA, NumberLength * sizeof(limb));
	std::memcpy(xC, xA, NumberLength * sizeof(limb));   // C <- A
	std::memcpy(zC, zA, NumberLength * sizeof(limb));
	duplicate(xA, zA, xA, zA);              /* A=2*A */
	while (d != e) {
		if (d < e) {
			r = d;   d = e;  e = r;   // swap d & e
			temp = xA; 	xA = xB; xB = temp;  // swap A and B
			temp = zA;  zA = zB; zB = temp;
		}
		/* do the first line of Table 4 whose condition qualifies */
		if (4 * d <= 5 * e && ((d + e) % 3) == 0) { /* condition 1 */
			r = (2 * d - e) / 3;
			e = (2 * e - d) / 3;
			d = r;
			add3(xT, zT, xA, zA, xB, zB, xC, zC); /* T = f(A,B,C) */
			add3(xT2, zT2, xT, zT, xA, zA, xB, zB); /* T2 = f(T,A,B) */
			add3(xB, zB, xB, zB, xT, zT, xA, zA); /* B = f(B,T,A) */
			temp = xA; 	xA = xT2; xT2 = temp;
			temp = zA; 	zA = zT2; zT2 = temp; /* swap A and T2 */
		}
		else if (4 * d <= 5 * e && (d - e) % 6 == 0) { /* condition 2 */
			d = (d - e) / 2;
			add3(xB, zB, xA, zA, xB, zB, xC, zC); /* B = f(A,B,C) */
			duplicate(xA, zA, xA, zA);       /* A = 2*A */
		}
		else if (d <= (4 * e)) {      /* condition 3 */
			d -= e;
			add3(xT, zT, xB, zB, xA, zA, xC, zC); /* T = f(B,A,C) */
			temp = xB; 	xB = xT; xT = xC; xC = temp;
			temp = zB;  zB = zT; zT = zC; zC = temp; /* circular permutation (B,T,C) */
		}
		else if ((d + e) % 2 == 0) {  /* condition 4 */
			d = (d - e) / 2;
			add3(xB, zB, xB, zB, xA, zA, xC, zC); /* B = f(B,A,C) */
			duplicate(xA, zA, xA, zA);        /* A = 2*A */
		}
		else if (d % 2 == 0) {         /* condition 5 */
			d /= 2;
			add3(xC, zC, xC, zC, xA, zA, xB, zB); /* C = f(C,A,B) */
			duplicate(xA, zA, xA, zA);        /* A = 2*A */
		}
		else if (d % 3 == 0) {        /* condition 6 */
			d = d / 3 - e;
			duplicate(xT, zT, xA, zA); /* T1 = 2*A */
			add3(xT2, zT2, xA, zA, xB, zB, xC, zC); /* T2 = f(A,B,C) */
			add3(xA, zA, xT, zT, xA, zA, xA, zA); /* A = f(T1,A,A) */
			add3(xT, zT, xT, zT, xT2, zT2, xC, zC); /* T1 = f(T1,T2,C) */
			temp = xC; xC = xB; xB = xT; xT = temp;
			temp = zC; zC = zB;	zB = zT; zT = temp; /* circular permutation (C,B,T) */
		}
		else if ((d + e) % 3 == 0) {     /* condition 7 */
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
			temp = xB; xB = xT; xT = temp;
			temp = zB; zB = zT; zT = temp; /* swap B and T */
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
- n : number to factor = TestNbr
- x, z : coordinates of P

Modifies: x3, z3,
work variables: TX, TZ, UX, UZ.
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
	if (!std::memcmp(x, x3, NumberLength * sizeof(limb))) {  // if x != x3
		std::memcpy(TZ, x, NumberLength * sizeof(limb));       // TZ = x
		std::memcpy(TX,  UZ, NumberLength * sizeof(limb));     // TX = UZ = 4*(x1*x2-z1*z2)^2
		modmult(z, TX,  UZ);                              // UZ = z*TX
		modmult(UX, TZ, z3);                              // z3 = UX*TZ
		std::memcpy(x3,  UZ, NumberLength * sizeof(limb));     // x3 = UZ
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
Modifies: x2, z2, 
work variables: UZ,  TX,  TZ
Uses global variable AA
*/
static void duplicate(limb *x2, limb *z2, const limb *x1, const limb *z1)
{
	/* N.B. all arithmetic is mod TestNbr */
	AddBigNbrModNB(x1, z1, TZ, TestNbr, NumberLength);    //  TZ = x1+z1 (mod testNbr)
	modmult(TZ, TZ, UZ);                                  //  UZ = (x1+z1)^2 (mod testNbr)
	SubtBigNbrModN(x1, z1, TZ, TestNbr, NumberLength);    //  TZ = x1-z1 (mod testNbr)
	modmult(TZ, TZ,  TX);                                 //  TX = (x1-z1)^2 (mod testNbr)
	modmult(UZ, TX, x2);                                  //  x2 = UZ* TX (mod testNbr)
	SubtBigNbrModN(UZ, TX, TZ, TestNbr, NumberLength);    //  TZ = UZ- TX = 4*x1*z1 (mod testNbr)
	modmult(AA, TZ, UZ);                                  // UZ =  TZ*AA (mod testNbr)
	AddBigNbrModNB(UZ, TX, UZ, TestNbr, NumberLength);    // UZ = ( TX+b* TZ) (mod testNbr)
	modmult(TZ, UZ, z2);                                   // z2 = ( TZ*u) (mod testNbr)
}
/* End of code adapted from Paul Zimmermann's ECM4C */

/* compute gcd of value & N. return 0 if value is zero or value = N, 
1 if gcd is 1 , 2 if gcd > 1.
The value of the gcd is returned in BiGD
Uses global variables BiGD, NumberLength */
static int gcdIsOne(const limb value[], const Znum& zN, int line) {
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
#ifdef log
	gmp_fprintf(logfile, "value = %Zd N = %Zd gcd = %Zd, line %d \n", 
		Temp1, zN, BiGD, line);
#endif
	if (BiGD < 2) {
		return (int)ZnumToLong(BiGD);    // GCD is less than 2.
	}
	return 2;      // GCD is greater than one.
}

static void GenerateSieve(int initial) {
	int i, j, Q, initModQ;
	for (i = 0; i < 10 * SIEVE_SIZE; i += SIEVE_SIZE)
	{
		std::memcpy(&sieve[i], sieve2310, SIEVE_SIZE);
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



#pragma warning(disable : 4996)
#ifdef __EMSCRIPTEN__
void showECMStatus(void) {
	char status[200];
	int elapsedTime;
	char *ptrStatus;

	/* check whether or not it's time for another status message */
	if ((lModularMult % yieldFreq) != 0) {
		return;
	}
	elapsedTime = (int)(tenths() - originalTenthSecond);
	if (elapsedTime / 10 <= oldTimeElapsed/10 +5) 	{
		return;  // less than 5 seconds since last message
	}

	oldTimeElapsed = elapsedTime;
	ptrStatus = status;
	std::strcpy(ptrStatus, lang ? "Transcurrió " : "Time elapsed: ");
	ptrStatus += std::strlen(ptrStatus);
	GetDHMS(&ptrStatus, elapsedTime / 10); // days, hours, mins, secs to o/p buffer
	switch (StepECM)
	{
	case 1:
		std::strcpy(ptrStatus, lang ? "Paso 1: " : "Step 1: ");
		ptrStatus += std::strlen(ptrStatus);
		ptrStatus += std::sprintf(ptrStatus, "%2d%%", indexPrimes / (nbrPrimes / 100));
		break;

	case 2:
		std::strcpy(ptrStatus, lang ? "Paso 2: " : "Step 2: ");
		ptrStatus += std::strlen(ptrStatus);
		ptrStatus += std::sprintf(ptrStatus, "%2d%%", maxIndexM == 0 ? 0 : indexM / (maxIndexM / 100));
		break;
	}

	*ptrStatus++ = '\0';                 // add null terminator
	printf_s("%s\n", status);            // send status to screen             
#ifdef log
	fprintf_s(logfile, "%s\n", status);  // send status to log file
#endif
}

#endif

/* see https://en.wikipedia.org/wiki/Fermat%27s_factorization_method */
/* Perform (modified) Lehman algorithm. Try to factorise nbr.
This produces a result when nbr is a product of two factors which are not
too far apart in value. e.g larger factor is less than k times the smaller.
The factors produced may not be primes. 
We need to find integers k, a and b such that 4kn = a^2 – b^2.
It can be shown that 1 <= k <= n^(1/3) and
sqrt(4kn) <= a <= sqrt(4kn) + (n^(1/6)/(4.sqrt(k)) 
However, searching the whole range of possible values would be much too slow. */
void LehmanZ(const Znum &nbr, int k, Znum &factor) {
	const long long bitsSqr[] = {
		0x0000000000000003, // 3  -> 3
		0x0000000000000013, // 5  -> 19
		0x0000000000000017, // 7  -> 23
		0x000000000000023B, // 11 -> 571
		0x000000000000161B, // 13 -> 5659
		0x000000000001A317, // 17 -> 107287 = 17*6311
		0x0000000000030AF3, // 19 -> 19941
		0x000000000005335F, // 23 -> 304831 = 569*599
		0x0000000013D122F3, // 29 -> 332473075 = 5^2 11 *107*11299
		0x00000000121D47B7, // 31 -> 303908791 = 151 * 2012641
		0x000000165E211E9B, // 37 -> 96068509339 = 11 * 47 * 109 * 1704763
		0x000001B382B50737, // 41 -> 1870503675703 = 13^2 * 67 * 71 * 131 * 17761
		0x0000035883A3EE53, // 43 -> 3678700564051 = 97 * 557 * 68087519 
		0x000004351B2753DF, // 47 -> 4626135339999 = 3 * 131 * 227 * 51856109
		0x0012DD703303AED3, // 53 -> 5310023542746835 = 5 * 11 * 19 * 757 * 6712499659
		0x022B62183E7B92BB, // 59 -> 156326468341437115 = 5 * 246773 * 126696574051
		0x1713E6940A59F23B, // 61 -> 1662926210933060155  = 5 * 283 * 769 * 887 * 13183 * 130693
	};
	int primes[] = { 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61 };
	int nbrs[17];
	int diffs[17];
	int i, j, m, r;
	Znum sqrRoot, nextroot;
	Znum a, c, sqr, val;

	if (isEven(nbr))
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
			r = (int)ZnumToLong((k + nbr)& 3) ;  // r = (k+nbr)%4 i.e. 0 or 2 ?
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
		nbrs[i] = (int)ZnumToLong(c % primes[i]);
		diffs[i] = m * (int)ZnumToLong(((a %pr) * 2 + m) % pr);
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
				if (verbose > 0) {
					std::cout << "Lehman factor found. k = " << k << " N= " << nbr
						<< " factor = " << factor << '\n';
				}
				return;
			}
		}

		for (i = 0; i < 17; i++)
		{
			nbrs[i] = (nbrs[i] + diffs[i]) % primes[i];
			diffs[i] = (diffs[i] + 2 * m * m) % primes[i];
		}
	}

	/* failed to factorise in 10000 cycles, so stop trying */
	factor = 1;     // Factor not found.
}

static void upOneLine(void) {
	SetConsoleCursorPosition(hConsole, coordScreen);
}

/* can return value:
FACTOR_NOT_FOUND   (not used??)
CHANGE_TO_SIQS
FACTOR_FOUND   - value of factor is returned in global variable BiGD
                 and in Zfactor
ERROR              (not used??)  
uses work variables UX, UZ, TX, TZ */
static enum eEcmResult ecmCurve(const Znum &zN, Znum &Zfactor) {
	/* arrays below are static to avoid risk of stack overflow */
	static limb A0[MAX_LEN];     // A0 <- 2*(ElipCurvNo+1)/(3*(ElipCurvNo+1)^2 - 1) (mod TestNbr)
	static limb A02[MAX_LEN];    // A02 <- A0^2
	static limb A03[MAX_LEN];    // A03 <- A0^3

	static limb W1[MAX_LEN];
	static limb W2[MAX_LEN];
	static limb WX[MAX_LEN];
	static limb WZ[MAX_LEN];
	static limb X[MAX_LEN];
	static limb Z[MAX_LEN];
	static limb Xaux[MAX_LEN];     // copy of X
	static limb Zaux[MAX_LEN];     // copy of Z
	static limb root[GROUP_SIZE][MAX_LEN];

	for (;;) {  // only exit from this loop is by a return statement
#ifdef __EMSCRIPTEN__
		char *ptrText;
#endif
		int I, Pass;
		int i, j, u;
		long long L1, L2, LS, P, IP, Paux = 1;

		ElipCurvNo++;   // increment curve number
		oldTimeElapsed = (int)(tenths() - originalTenthSecond);

		L2 = mpz_sizeinbase(ZT(zN),10);   // Get number of digits.
		if (L2 > 30 && L2 <= 95)          // If between 31 and 95 digits...
		{                                 // switch to SIQS when curve No reaches limit
			int limit = limits[((int)L2 - 26) / 5];  // e.g if L2<=50, limit=10
			if (ElipCurvNo  >= limit) { 
#ifdef log
				fprintf_s(logfile, "Change to SIQS\n");
#endif
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
			//if (Zfactor > LIMB_RANGE) {
			if (Zfactor > 1) {
				foundByLehman = true;     // Factor found.
#ifdef log
				gmp_fprintf(logfile, "Lehman factor found. k = %d  N=  %Zd factor = %Zd \n",
					k, ZT(zN), ZT(Zfactor));
#endif
				if (verbose > 1)
					gmp_printf("Lehman factor found. k = %d  N=  %Zd factor = %Zd \n",
						k, ZT(zN), ZT(Zfactor));
				return FACTOR_FOUND;
			}
			/*else if (Zfactor > 1)
				std::cout << "Ignored small Lehman factor found. k = " << k << " N= " << zN
				<< " factor = " << Zfactor << '\n';*/
		}

		/* set L1, L2, LS, Paux and nbrPrimes according to value of ElipCurvNo */
		L1 = 2000;
		L2 = 200000;
		LS = 45;       // used as prime limit in pass 1
		Paux = ElipCurvNo;
		nbrPrimes = 303; /* Number of primes less than 2000 */
		if (ElipCurvNo > 25) {
			if (ElipCurvNo < 326) {   // curves 26 to 325
				L1 = 50000;
				L2 = 5000000;
				LS = 224;
				Paux = ElipCurvNo - 24;
				nbrPrimes = 5133; /* Number of primes less than 50000 */
			}
			else {
				if (ElipCurvNo < 2000) {  // curves 326 to 1999
					L1 = 1000000;
					L2 = 100000000;
					LS = 1001;
					Paux = ElipCurvNo - 299;
					nbrPrimes = 78498; /* Number of primes less than 1000000 */
				}
				else {                  //  curve >= 2000
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
		if (ElipCurvNo > 1 || verbose > 0) {
			ptrText = ptrLowerText;  // Point after number that is being factored.
			auto elapsedTime = (int)(tenths() - originalTenthSecond);
			GetDHMSt(&ptrText, elapsedTime);
			std::strcpy(ptrText, lang ? " ECM Curva " : " ECM Curve ");
			ptrText += std::strlen(ptrText);
			ptrText += std::sprintf(ptrText, "%4d", ElipCurvNo);   // Show curve number.

			std::strcpy(ptrText, lang ? " usando límites B1=" : " using bounds B1=");
			ptrText += std::strlen(ptrText);
 			ptrText += std::sprintf(ptrText, "%5lld", L1);       // Show first bound.

			std::strcpy(ptrText, lang ? " y B2=" : " and B2=");
			ptrText += std::strlen(ptrText);
			ptrText += std::sprintf(ptrText, "%7lld \n", L2); // Show second bound.

			if (first) {
				if (!GetConsoleScreenBufferInfo(hConsole, &csbi))
				{
					ErrorDisp(__FUNCTION__);
					Beep(750, 1000);
				}
				coordScreen.X = csbi.dwCursorPosition.X;  // save cursor co-ordinates
				coordScreen.Y = csbi.dwCursorPosition.Y;
				if (csbi.dwCursorPosition.Y >= csbi.dwSize.Y - 1)
					coordScreen.Y--;  // if window is full, allow for text scrolling up
			}
			else
				upOneLine();

			printf_s("%s", ptrLowerText);        // send status to stdout (screen)
			first = false;
#ifdef log
			fprintf_s(logfile, "%s", ptrLowerText);   // send status to l,og file
			std::fflush(logfile);
#endif
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
		}
#endif

		//  Compute A0 <- 2 * (ElipCurvNo+1)*modinv(3 * (ElipCurvNo+1) ^ 2 - 1, N) mod N
		// Aux2 <- 1 in Montgomery notation.
		std::memcpy(Aux2, MontgomeryMultR1, NumberLength * sizeof(limb));
		modmultInt(Aux2, ElipCurvNo + 1, Aux2);    // Aux2 <- ElipCurvNo + 1 (mod TestNbr)
		modmultInt(Aux2, 2, Aux1);                 // Aux1 <- 2*(ElipCurvNo+1) (mod TestNbr)
		modmultInt(Aux2, ElipCurvNo + 1, Aux3);    // Aux3 <- (ElipCurvNo + 1)^2 (mod TestNbr)
		modmultInt(Aux3, 3, Aux3);                 // Aux3 <- 3*(ElipCurvNo + 1)^2 (mod TestNbr)
		// Aux2 <- 3*(ElipCurvNo + 1)^2 - 1 (mod TestNbr)
		SubtBigNbrModN(Aux3, MontgomeryMultR1, Aux2, TestNbr, NumberLength);
		ModInvBigNbr(Aux2, Aux2, TestNbr, NumberLength);
		modmult(Aux1, Aux2, A0);       // A0 <- 2*(ElipCurvNo+1)/(3*(ElipCurvNo+1)^2 - 1) (mod TestNbr)
#ifdef log
		logf(A0);
#endif
		//  if A0*(A0 ^ 2 - 1)*(9 * A0 ^ 2 - 1) mod N=0 then select another curve.
		modmult(A0, A0, A02);          // A02 <- A0^2
		modmult(A02, A0, A03);         // A03 <- A0^3
		SubtBigNbrModN(A03, A0, Aux1, TestNbr, NumberLength);  // Aux1 <- A0^3 - A0
		modmultInt(A02, 9, Aux2);      // Aux2 <- 9*A0^2
		SubtBigNbrModN(Aux2, MontgomeryMultR1, Aux2, TestNbr, NumberLength); // Aux2 <- 9*A0^2-1
		modmult(Aux1, Aux2, Aux3);
#ifdef log
		logf(Aux3);                     // Aux3 =  (A0^3 - A0)* (9*A0^2)
#endif
		if (BigNbrIsZero(Aux3, NumberLength)) {
			continue;  // select another curve
		}

		//   Z <- 4 * A0 mod N
		modmultInt(A0, 4, Z);
		//   A = (-3 * A0 ^ 4 - 6 * A0 ^ 2 + 1)*modinv(4 * A0 ^ 3, N) mod N
		modmultInt(A02, 6, Aux1);      // Aux1 <- 6*A0^2
		SubtBigNbrModN(MontgomeryMultR1, Aux1, Aux1, TestNbr, NumberLength); // Aux1 = 1-6*A0^2
		modmult(A02, A02, Aux2);       // Aux2 <- A0^4
		modmultInt(Aux2, 3, Aux2);     // Aux2 <- 3*A0^4
		SubtBigNbrModN(Aux1, Aux2, Aux1, TestNbr, NumberLength);   // Aux1 = (-3 * A0 ^ 4 - 6 * A0 ^ 2 + 1)
		modmultInt(A03, 4, Aux2);      // Aux2 <- 4*A0^3
		ModInvBigNbr(Aux2, Aux3, TestNbr, NumberLength);
		modmult(Aux1, Aux3, A0);     // A0 = (-3 * A0 ^ 4 - 6 * A0 ^ 2 + 1)*modinv(4 * A0 ^ 3, N) mod N
#ifdef log
		logf(Z);
		logf(A0);
#endif

		//   AA <- (A + 2)*modinv(4, N) mod N
		modmultInt(MontgomeryMultR1, 2, Aux2);           // Aux2 <- 2
		AddBigNbrModNB(A0, Aux2, Aux1, TestNbr, NumberLength); // Aux1 <- A0+2
		modmultInt(MontgomeryMultR1, 4, Aux2);            // Aux2 <- 4
		ModInvBigNbr(Aux2, Aux2, TestNbr, NumberLength);   // Aux2 = 1/4
		modmult(Aux1, Aux2, AA);                // AA = Aux1 * Aux2 = Aux1/4
		//   X <- (3 * A0 ^ 2 + 1) mod N
		modmultInt(A02, 3, Aux1);    // Aux1 <- 3*A0^2
		AddBigNbrModNB(Aux1, MontgomeryMultR1, X, TestNbr, NumberLength);
		/**************/
		/* First step */
		/**************/
		std::memcpy(Xaux, X, NumberLength * sizeof(limb));   // Xaux = X
		std::memcpy(Zaux, Z, NumberLength * sizeof(limb));   // Zaux = Z
#ifdef log
		fprintf_s(logfile, " Start first step\n");
		logf(Xaux);
#endif
		// GcdAccumulated = 1
		std::memcpy(GcdAccumulated, MontgomeryMultR1, (NumberLength + 1) * sizeof(limb));
		for (Pass = 0; Pass < 2; Pass++) {
#ifdef log
			std::fprintf(logfile, "starting pass %d \n", Pass);
#endif
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
				// GcdAccumulated *= Z
				modmult(GcdAccumulated, Z, Aux1);
				std::memcpy(GcdAccumulated, Aux1, NumberLength * sizeof(limb));
#ifdef log
				logf(GcdAccumulated);
#endif
			}
			else {
				if (gcdIsOne(Z, zN, __LINE__) > 1) {
					Zfactor = BiGD;    // copy factor to global Znum
#ifdef log
					logf(Z);
					fprintf_s(logfile, "factor found Pass 1 \n");
#endif
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
					// GcdAccumulated *= Z;
					modmult(GcdAccumulated, Z, Aux1);
					std::memcpy(GcdAccumulated, Aux1, NumberLength * sizeof(limb));
#ifdef log
					fprintf_s(logfile, "indexM = %d ", indexM);
					logf(GcdAccumulated);
#endif
				}
				else {
					if (gcdIsOne(Z, zN, __LINE__) > 1) {
						Zfactor = BiGD;      // copy factor to global Znum
#ifdef log
						logf(Z);
						fprintf_s(logfile, "Pass = 1  - factor found \n");
#endif
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
						// GcdAccumulated *= Z;
						modmult(GcdAccumulated, Z, Aux1);
						std::memcpy(GcdAccumulated, Aux1, NumberLength * sizeof(limb));
#ifdef log
						fprintf_s(logfile, "%d of %d  ", i, 10 * SIEVE_SIZE);
						logf(GcdAccumulated);
#endif
					}
					else {
						if (gcdIsOne(Z, zN, __LINE__) > 1) {
							Zfactor = BiGD;    // copy factor to global Znum
#ifdef log
							logf(Z);
							fprintf_s(logfile, "Pass = 1  - factor found \n");
#endif
							return FACTOR_FOUND;   // ** found a factor in global variable BiGD**
						}
					}
				}
				P += 20 * SIEVE_SIZE;
			} while (P < L1);

			if (Pass == 0) {
				if (BigNbrIsZero(GcdAccumulated, NumberLength))
				{ // If GcdAccumulated is multiple of TestNbr, continue.
					std::memcpy(X, Xaux, NumberLength * sizeof(limb));
					std::memcpy(Z, Zaux, NumberLength * sizeof(limb));
#ifdef log
					logf(X);
					logf(Z);
#endif
					continue; 
				}
				if (gcdIsOne(GcdAccumulated, zN, __LINE__) > 1) {
					Zfactor = BiGD;      // copy factor to global Znum
#ifdef log
					fprintf_s(logfile, "Pass = 0  - factor found \n");
#endif
					return FACTOR_FOUND;   // ** found a factor **  in global variable BiGD
				}
				break;  // exit for loop
			}
		} /* end for (Pass = 0; Pass < 2; Pass++) */

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
		std::memcpy(&sieve2310[HALF_SIEVE_SIZE], &sieve2310[0], HALF_SIEVE_SIZE);
		std::memcpy(Xaux, X, NumberLength * sizeof(limb));  // (X:Z) -> Q (output
		std::memcpy(Zaux, Z, NumberLength * sizeof(limb));  //         from step 1)
#ifdef log
		fprintf_s(logfile, " Start second step\n");
		logf(Xaux);
		logf(Zaux);
#endif
		for (Pass = 0; Pass < 2; Pass++) {
			int Qaux, J;
			std::memcpy(GcdAccumulated, MontgomeryMultR1, NumberLength * sizeof(limb));
			std::memcpy(UX, X, NumberLength * sizeof(limb));
			std::memcpy(UZ, Z, NumberLength * sizeof(limb));  // (UX:UZ) -> Q 
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
			modmult(Aux2, UX, Z); // (X:Z) -> 3Q
			for (I = 5; I < SIEVE_SIZE; I += 2) {
				std::memcpy(WX, X, NumberLength * sizeof(limb));
				std::memcpy(WZ, Z, NumberLength * sizeof(limb));
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
					// GcdAccumulated *= Aux1 (Aux1 = W1 -W2)
					modmult(GcdAccumulated, Aux1, Aux2);
					std::memcpy(GcdAccumulated, Aux2, NumberLength * sizeof(limb));
#ifdef log
					fprintf_s(logfile, "%d of %d  ", I, SIEVE_SIZE);
					logf(GcdAccumulated);
#endif
				}
				else {
					if (gcdIsOne(Aux1, zN, __LINE__) > 1) {
						Zfactor = BiGD;     // copy factor to global Znum
#ifdef log
						logf(Aux1);
						fprintf_s(logfile, "Pass 1, I = %d factor found \n", I);
#endif
						return FACTOR_FOUND;   // ** found a factor in global variable BiGD
					}
				}
				if (I == HALF_SIEVE_SIZE) {
					std::memcpy(DX, X, NumberLength * sizeof(limb));
					std::memcpy(DZ, Z, NumberLength * sizeof(limb));  // (DX:DZ) -> HALF_SIEVE_SIZE*Q
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
				std::memcpy(UX, WX, NumberLength * sizeof(limb));  // (UX:UZ) <-
				std::memcpy(UZ, WZ, NumberLength * sizeof(limb));  // Previous (X:Z)
			} /* end for I */
#ifdef log
			fprintf_s(logfile, "end 'for I' I=%d Pass = %d \n", I, Pass);
#endif

			AddBigNbrModNB(DX, DZ, Aux1, TestNbr, NumberLength);
			modmult(Aux1, Aux1, W1);
			SubtBigNbrModN(DX, DZ, Aux1, TestNbr, NumberLength);
			modmult(Aux1, Aux1, W2);
			modmult(W1, W2, X);
			SubtBigNbrModN(W1, W2, Aux1, TestNbr, NumberLength);
			modmult(Aux1, AA, Aux2);
			AddBigNbrModNB(Aux2, W2, Aux3, TestNbr, NumberLength);
			modmult(Aux1, Aux3, Z);
			std::memcpy(UX, X, NumberLength * sizeof(limb));
			std::memcpy(UZ, Z, NumberLength * sizeof(limb));    // (UX:UZ) -> SIEVE_SIZE*Q
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
						j = sieveidx[i];     // 0 < J < HALF_SIEVE_SIZE
						if (sieve[J + j] != 0 && sieve[J - 1 - j] != 0) {
							continue; // Do not process if both are composite numbers.
						}

						SubtBigNbrModN(Aux1, root[i], M, TestNbr, NumberLength);
						// GcdAccumulated *= M
						modmult(GcdAccumulated, M, Aux2);
						std::memcpy(GcdAccumulated, Aux2, NumberLength * sizeof(limb));
#ifdef log
						fprintf_s(logfile, "%d of %d  (%d of %d)",
							indexM, maxIndexM, i, GROUP_SIZE);
						logf(GcdAccumulated);
#endif
					}
					if (Pass != 0) {
						if (BigNbrIsZero(GcdAccumulated, NumberLength)) {
#ifdef  log
							fprintf_s(logfile, "break. This curve cannot factor the number \n");
#endif //  log
							break;  // This curve cannot factor the number.
						}
						if (gcdIsOne(GcdAccumulated, zN, __LINE__) > 1) {
							Zfactor = BiGD;   // copy factor to global Znum
#ifdef log
							fprintf_s(logfile, "Pass = 1, indexM = %d factor found \n", indexM);
#endif
							return FACTOR_FOUND;    // ** found a factor in global variable BiGD
						}
					}  // end if (Pass != 0)
				}   // End if (indexM >= Qaux)

				if (indexM != 0) { // Update (X:Z)
					std::memcpy(WX, X, NumberLength * sizeof(limb));
					std::memcpy(WZ, Z, NumberLength * sizeof(limb));
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
					std::memcpy(UX, WX, NumberLength * sizeof(limb));
					std::memcpy(UZ, WZ, NumberLength * sizeof(limb));
				}
			} // end for (indexM = 0; indexM <= maxIndexM; indexM++)

			if (Pass == 0) {
				int rc;
				if (BigNbrIsZero(GcdAccumulated, NumberLength)) { // If GcdAccumulated is zero
					std::memcpy(X, Xaux, NumberLength * sizeof(limb));  // X = Xaux
					std::memcpy(Z, Zaux, NumberLength * sizeof(limb));  // Z = Zaux
#ifdef log
					fprintf_s(logfile, "GcdAccumulated = 0 ");
					logf(X);
					logf(Z);
#endif
					continue; // multiple of TestNbr, continue with next pass 
				}
				rc = gcdIsOne(GcdAccumulated, zN, __LINE__);
				if (rc == 1) {
#ifdef log
					logf(GcdAccumulated);
					fprintf_s(logfile, "exit from curve %d  \n", ElipCurvNo);
#endif
					break;    // GCD is one, so this curve does not find a factor.
				}
				if (rc == 0) {
					continue;  // GcdAccumulated = N or 0
				}

				/* rc = 2 */
				if (BiGD != zN)
				{           // GCD is not 1 or TestNbr
					Zfactor = BiGD;    // copy factor to global Znum
#ifdef log
					logf(GcdAccumulated);
					fprintf_s(logfile, "factor found pass %d \n", Pass);
#endif
					return FACTOR_FOUND;     // ** found a factor in global variable BiGD
				}
			}
		} /* end for Pass */
#ifdef log
		fprintf_s(logfile, "end of processing for curve %d", ElipCurvNo);
#endif
	}       /* End curve calculation. loop continues with next curve */
	/* only exit from for loop above is by return statement */
}

/* initialise variables. zN = number to be factored */
static void ecminit(const Znum &zN) {
	int P, Q;

	ZtoBig(TestNbrBI, zN);   /* copy zN to TestNbr. NB throw exception
				             if zN is too large! (more than about 23,000 digits) */
#ifdef log
	logf(TestNbr);
#endif
	GetYieldFrequency(NumberLength);      //get yield frequency (used by showECMStatus)
	GetMontgomeryParms(NumberLength);
#ifdef log
	logf(MontgomeryMultR1);
	logf(MontgomeryMultR2);
#endif
	first = true;

	/* set variables to zero */
	std::memset(M, 0, NumberLength * sizeof(limb));
	std::memset(DX, 0, NumberLength * sizeof(limb));
	std::memset(DZ, 0, NumberLength * sizeof(limb));
	std::memset(W3, 0, NumberLength * sizeof(limb));
	std::memset(W4, 0, NumberLength * sizeof(limb));
	ptrLowerText = lowerText;

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

	modmultCallback = showECMStatus;   // Set callback function pointer
	lModularMult = 0;   // reset counter used to control status display
}

/* find factor of zN. returns true if successful. The factor found is returned 
   in Zfactor */
bool ecm(const Znum &zN, fList &Factors, Znum &Zfactor) {

	ElipCurvNo = 1;  // start with 1st curve
#ifdef log
	std::string name1 = std::tmpnam(nullptr);
	auto dotpos = name1.find('.');
	if (dotpos != std::string::npos)
		name1.resize(dotpos);
	name1 += "ECMlog.txt";
	logfile = std::fopen(name1.c_str(), "a");
	std::cout << "log file name = " << name1 << '\n';
#endif

	ecminit(zN);  // initialise values

	foundByLehman = false;
	do {
		enum eEcmResult ecmResp = ecmCurve(zN, Zfactor);
		if (ecmResp == CHANGE_TO_SIQS) {    // Perform SIQS
			FactoringSIQS(zN, Zfactor); // factor found is returned in Zfactor
			Factors.ct.siqs++;
			break;
		}
		else if (ecmResp == FACTOR_FOUND) {
			if (foundByLehman)
				Factors.ct.leh++;
			else
				Factors.ct.ecm++;
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
#ifdef log
	std::fclose(logfile);
#endif
	return true;
}