/*
This file is part of Alpertron Calculators.
Copyright 2016 Dario Alejandro Alpern
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

SIQS = Self-initializing quadratic sieve
see http://www.mersennewiki.org/index.php/Self-Initializing_Quadratic_Sieve
or https://en.wikipedia.org/wiki/Quadratic_sieve
*/
 
/* ***************************************************************
Some results from profiling using Visual Studio (20/6/2018):
PerformSiqsSeiveStage	72.53%
SeiveThread				 7.63%
__gmp_qdiv_qrnnd		 3.78%    (extended precision division, from SIQS)
MontgomeryMult9			 3.52%    (modular multiplication, from ECM)
DoTrialDivision	 3.10%
(this top 5 totals 90%, the rest is scattered over many functions)
For numbers > 100 digits, SIQS may not be used at all, in which case
the results would look completely different.
**************************************************************** */

#include "pch.h"

#pragma warning(disable : 4996)

extern HANDLE hConsole;
static bool first = true;   /* set true for 1st message of set, then false for the 
							rest. used to control cursor repositioning so that 
							subsequent messages overwrite the 1st one */
static COORD coordScreen = { 0, 0 };    // home position for the cursor 
static CONSOLE_SCREEN_BUFFER_INFO csbi;

// These defines are valid for factoring up to 10^110.
#define MAX_NBR_FACTORS         13
#define MAX_PRIMES          150000
#define MAX_LIMBS_SIQS          15
#define MAX_FACTORS_RELATION    80    // increased 25/6/19 in line with DA's calculator
#define LENGTH_OFFSET            0
#define MAX_SIEVE_LIMIT     200000   // increased 25/6/19 in line with DA's calculator
									 // for >= 108 digits this limit would be reached

// change DEBUG_SIQS to 0 to turn debug messages off, 1 to turn on
#define DEBUG_SIQS               0
#define __EMSCRIPTEN__


#ifdef __EMSCRIPTEN__
extern char lowerText[], *ptrLowerText;
static char *ptrSIQSStrings;
static int startSieveTenths;
#endif

typedef struct {
	int value;
	int modsqrt;
	int Bainv2[MAX_NBR_FACTORS];
	int Bainv2_0;
	int soln1;
	int difsoln;
} PrimeSieveData;

typedef struct {
	int value;
	int exp[6];
} PrimeTrialDivisionData;

static unsigned char SIQSInfoText[300];
#if DEBUG_SIQS
static char output[1000];
static int debug_ctr = 0;
const static int debugLimit = 1000;
#endif
const static int numberThreads = 1;
static int matrixBLength;
static long trialDivisions;
static long smoothsFound;     // number of full congruences.
static long totalPartials;    // number of partial congruences
static long partialsFound;    // number of built congruences
static int nbrPartials;        // number of partials
static long ValuesSieved;
static int nbrFactorBasePrimes;
static int congruencesFound = 0;
static long polynomialsSieved;
static int multiplier = 1;
//extern int multiplier;
static int nbrFactorsA;
static int afact[MAX_NBR_FACTORS];
//static long startTime;
//static int biModulus[MAX_LIMBS_SIQS] = { 0 };
//static int biTestNbr2[MAX_LIMBS_SIQS] = { 0 };
//static int biQuadrCoeff[MAX_LIMBS_SIQS] = { 0 };
//static int biLinearDelta[MAX_LIMBS_SIQS][MAX_LIMBS_SIQS] = { 0 };
static Znum zModulus;
static Znum zTestNbr2;
static Znum biQuadrCoeff;
static Znum biLinearDelta[MAX_LIMBS_SIQS];

static long largePrimeUpperBound;
static unsigned char logar2;
static int aindex[MAX_NBR_FACTORS] = { 0 };
//  static Thread threadArray[];
//static PrimeSieveData primeSieveData[MAX_PRIMES + 3] = { 0 };
static std::vector < PrimeSieveData>  primeSieveData;
//static PrimeTrialDivisionData primeTrialDivisionData[MAX_PRIMES] = { 0 };
static int span;
static int indexMinFactorA;
const static int threadNumber = 0;
static int nbrThreadFinishedPolySet;
static unsigned int oldSeed;
static unsigned int newSeed;
static int NbrPolynomials;
static int SieveLimit;
//static int matrixPartial[MAX_PRIMES * 8][MAX_LIMBS_SIQS / 2 + 4] = { 0 };
struct MatPart {
	int nbrPartials;
	int Divid;
	int hashindex;
	Znum value;
	int seed;
};
static std::vector <MatPart> matrixPartial;
static Znum vectLeftHandSide[MAX_PRIMES + 50] = { 0 };
static int matrixPartialHashIndex[2048] = { 0 };
static std::vector <std::vector<int>> matrixB;
static int amodq[MAX_NBR_FACTORS] = { 0 };
static int tmodqq[MAX_NBR_FACTORS] = { 0 };
static char threshold;
static int smallPrimeUpperLimit;
static int firstLimit;
static int secondLimit;
static int thirdLimit;
static int nbrPrimesUsed = 0;
/* this struct is defined as an easier way to manage memory allocation and release.
A single malloc & free are sufficient. The disadvantage is that
all arrays in the structure are accessed via a pointer e.g.
ag->vectExpParity[i] rather than vectExpParity[i] */
struct AG{
	int vectExpParity[MAX_PRIMES + 50];
	int matrixAV[MAX_PRIMES];
	int matrixV[MAX_PRIMES];
	int matrixV1[MAX_PRIMES];
	int matrixV2[MAX_PRIMES];
	int matrixXmY[MAX_PRIMES];
	int newColumns[MAX_PRIMES];
	// Matrix that holds temporary data
	int matrixCalc3[MAX_PRIMES];
	int matrixTemp2[MAX_PRIMES];
	PrimeTrialDivisionData primeTrialDivisionData[MAX_PRIMES];
	char primesUsed[MAX_PRIMES];
};
static AG *ag = nullptr;      // storage allocated using malloc

static int nbrPrimes2;
//static BigInteger factorSiqs;
const static unsigned char onlyFactoring = FALSE;
static int matrixRows, matrixCols;
static PrimeSieveData *firstPrimeSieveData;

/* function forward declarations */
static bool InsertNewRelation(
	int *rowMatrixB,
	Znum &biT, Znum &biU, Znum &biR);
static void BlockLanczos(void);
static void ShowSIQSStatus(int *rowMatrixB);
static unsigned int getFactorsOfA(unsigned int seed, int *indexA);
static void sieveThread(Znum &result);

/* Sum = Nbr1+Nbr2 (mod Mod) */
static void AddZnumModN(const Znum& Nbr1, const Znum& Nbr2, Znum& Diff, const Znum& Mod) {
	//Diff = (Nbr1 + Nbr2) % Mod;
	mpz_add(ZT(Diff), ZT(Nbr1), ZT(Nbr2));    // Diff = Nbr1 + Nbr2
	while (Diff < 0)
		Diff += Mod;
	mpz_mod(ZT(Diff), ZT(Diff), ZT(Mod));  // Diff %= Mod
}

/* Diff = Nbr1-Nbr2 (mod Mod)*/
static void SubtractZnumModN(const Znum& Nbr1, const Znum& Nbr2, Znum& Diff, const Znum& Mod) {
	//Diff = (Nbr1 - Nbr2) % Mod;
	mpz_sub(ZT(Diff), ZT(Nbr1), ZT(Nbr2));    // Diff = Nbr1 - Nbr2
	while (Diff < 0)
		Diff += Mod;
	mpz_mod(ZT(Diff), ZT(Diff), ZT(Mod));  // Diff %= Mod
}

/* calculate dividend%divisor. remainder has same type as divisor
eror occurs if divisor is zero!! */
static mpir_ui RemDivZnumByInt(const Znum& Dividend, mpir_ui divisor) {
	return mpz_fdiv_ui(ZT(Dividend), divisor);
	/* Note that using % operator with Znums returns a Znum, even though the
	divisor is an integer. This way is much more efficient */
}


/* Prod = Nbr1*Nbr2 (mod Mod) */
static void MultZnumModN(const Znum& Nbr1, const Znum& Nbr2, Znum& Prod, const Znum& Mod) {
	//Prod = Nbr1*Nbr2;
	mpz_mul(ZT(Prod), ZT(Nbr1), ZT(Nbr2));
	mpz_mod(ZT(Prod), ZT(Prod), ZT(Mod));
}

/* Prod = Nbr1*Nbr2 (mod Mod) */
static void MultZnumByIntModN(const Znum& Nbr1, int Nbr2, Znum& Prod, const Znum& Mod) {
	//Prod = Nbr1*Nbr2;
	mpz_mul_si(ZT(Prod), ZT(Nbr1), Nbr2);
	mpz_mod(ZT(Prod), ZT(Prod), ZT(Mod));
}


#ifdef __EMSCRIPTEN__
static void showMatrixSize(char *SIQSInfoText, int rows, int cols)
{
	char *ptrText = ptrLowerText;  // Point after number that is being factored.
	strcpy(ptrText, lang ? "Resolviendo la matriz de congruencias de " : "Solving ");
	ptrText += strlen(ptrText);
	int2dec(&ptrText, rows);   // Show number of rows.
	strcpy(ptrText, " X ");
	ptrText += strlen(ptrText);
	int2dec(&ptrText, cols);   // Show number of columns.
	strcpy(ptrText, lang ? " usando el algoritmo de Lanczos en bloques.\n" :
		" congruence matrix using Block Lanczos algorithm.\n");

	printf_s("%s", ptrLowerText);
}

static void upOneLine(void) {
#if DEBUG_SIQS
	/* in debug, don't mess with cursor position */
#else
	SetConsoleCursorPosition(hConsole, coordScreen);
#endif
}

static void InitSIQSStrings(int SieveLimit, int nbrFactorBasePrimes)
{
	char *ptrText = ptrLowerText;  // Point after number that is being factored.
	strcpy(ptrText, lang ? "Parámetros de SIQS: " : "SIQS parameters: ");
	ptrText += strlen(ptrText);
	int2dec(&ptrText, nbrFactorBasePrimes);   // Show number of primes in factor base.
	strcpy(ptrText, lang ? " primos, límite de la criba: " : " primes, sieve limit: ");
	ptrText += strlen(ptrText);
	int2dec(&ptrText, SieveLimit);  // Show sieve limit.
	ptrSIQSStrings = ptrText;
	strcpy(ptrText, lang ? "\nBuscando el mejor multiplicador de Knuth-Schroeppel...\n" :
		"\nSearching for Knuth-Schroeppel multiplier...\n");
	printf_s("%s", ptrLowerText);
}

// Append multiplier and factor base to SIQS string.
static void getMultAndFactorBase(int multiplier, int FactorBase)
{
	char *ptrText = ptrSIQSStrings;
	strcpy(ptrText, lang ? "\nMultiplicador: " : "\nMultiplier: ");
	ptrText += strlen(ptrText);
	int2dec(&ptrText, multiplier);  // Show Knuth-Schroeppel multiplier.
	strcpy(ptrText, lang ? ", base de factores: " : ", factor base: ");
	ptrText += strlen(ptrText);
	int2dec(&ptrText, FactorBase);  // Show factor base.
	strcpy(ptrText, "\n");
	ptrText += strlen(ptrText);
	ptrSIQSStrings = ptrText;
}

/* print status: timeSieve = Sieve time so far in seconds*/
static void ShowSIQSInfo(int timeSieve, int congruencesFound, int matrixBLength, 
	int elapsedTime, int *rowMatrixB)
{
	char SIQSInfo[500] = { 0 };
	char *ptrText = SIQSInfo;
		// calculate %age of congruences found so far
	int percentage = (int)((congruencesFound * 100) / matrixBLength); 
	/* estimate time to completion (in seconds) based on time used so far and the ratio of 
	congruences still to be found vs congruences already found */
	int u = (int)((double)timeSieve * (double)(matrixBLength - congruencesFound) / (double)congruencesFound);
	u /= 2;

		GetDHMS(&ptrText, elapsedTime);
	if (lang) {  // español
		ptrText += sprintf(ptrText, "%5d congruencias halladas con %5d primos diferentes.\n"
			"Relaciones: %5d completas y %5d obtenidas de %6d parciales.\n"
			"progreso = %2d%%. ",
			congruencesFound, nbrPrimesUsed, smoothsFound, partialsFound, 
			totalPartials, percentage);
	}
	else {  // english
		ptrText += sprintf(ptrText, "%5d congruences found with %5d different primes \n"
			"Relations: %5d full and %5d found from %6d partials.\n"
			"progress = %2d%%. ",
			congruencesFound, nbrPrimesUsed, smoothsFound, partialsFound, 
			totalPartials, percentage);
	}

	if (timeSieve > 1 && congruencesFound > 10) {
		strcpy(ptrText, lang ? "Fin de la criba en " : "End sieve in ");
		ptrText += strlen(ptrText);
		GetDHMS(&ptrText, u);      // print estimated time to completion
	}

	strcpy(ptrText, "\n");

	if (first) {
		if (!GetConsoleScreenBufferInfo(hConsole, &csbi))
		{
			fprintf(stderr, "** GetConsoleScreenBufferInfo failed with %d!\n", GetLastError());
			Beep(750, 1000);
		}
		coordScreen.X = csbi.dwCursorPosition.X;  // save cursor co-ordinates
		coordScreen.Y = csbi.dwCursorPosition.Y;
		if (csbi.dwCursorPosition.Y >= csbi.dwSize.Y - 1)
			coordScreen.Y--;   // if window is full, allow for text scrolling up

	}
	else
		upOneLine();

	printf_s("%s", SIQSInfo);
	first = false;
}

#endif

/* profiling indicates that about 70% of CPU time during SIQS factoring is used within
this function, so any attempts to improve performance should probably focus on this
function. 
Uses global variable biLinearDelta*/
static void PerformSiqsSieveStage(PrimeSieveData primeSieveData[],
	short SieveArray[], int PolynomialIndex, Znum &zLinearCoeff)
{
	short logPrimeEvenPoly, logPrimeOddPoly;
	int currentPrime, F1, F2, F3, F4, X1, X2;
	int index, index2, indexFactorA;
	int mask;
	bool polyadd;
	int S1, G0, G1, G2, G3;
	int H0, H1, H2, H3, I0, I1, I2, I3;
	PrimeSieveData *rowPrimeSieveData;  // N.B. element soln1 is modified
	// above pointer starts at row 1 and is incremented through the array.
	// N.B. only element soln1 is modified

	F1 = PolynomialIndex;
	indexFactorA = 0;
	while ((F1 & 1) == 0) {
		F1 >>= 1;
		indexFactorA++;
		assert(indexFactorA < MAX_LIMBS_SIQS);
	}
	polyadd = (F1 & 2) != 0;
	if (polyadd)   // Adjust value of B as appropriate
	{              // according to the Gray code.
		zLinearCoeff += biLinearDelta[indexFactorA];
		zLinearCoeff += biLinearDelta[indexFactorA];
	}
	else {
		zLinearCoeff -= biLinearDelta[indexFactorA];
		zLinearCoeff -= biLinearDelta[indexFactorA];
	}
	indexFactorA--;
	X1 = SieveLimit << 1;
	rowPrimeSieveData = primeSieveData + 1;
	F1 = polyadd ? -rowPrimeSieveData->Bainv2[indexFactorA] :
		rowPrimeSieveData->Bainv2[indexFactorA];

#if DEBUG_SIQS
	if (debug_ctr <= debugLimit) {
		std::cout << "X1 = " << X1
			<< " F1 = " << F1
			<< " polyadd = " << polyadd
			<< " soln1 = " << rowPrimeSieveData->soln1
			<< " zLinearCoeff = " << zLinearCoeff
			<< " indexFactorA = " <<  indexFactorA << '\n';
		debug_ctr++;
	}
#endif

	if (((rowPrimeSieveData->soln1 += F1) & 1) == 0)  // soln1 is modified
	{
		*(SieveArray + 0) = (short)(logar2 - threshold);
		*(SieveArray + 1) = (short)(-threshold);
	}
	else {
		*(SieveArray + 0) = (short)(-threshold);
		*(SieveArray + 1) = (short)(logar2 - threshold);
	}

	if (((rowPrimeSieveData->soln1 + rowPrimeSieveData->Bainv2_0) & 1) == 0)
	{
		*(SieveArray + 0) += (short)((logar2 - threshold) << 8);
		*(SieveArray + 1) += (short)((-threshold) << 8);
	}
	else {
		*(SieveArray + 0) += (short)((-threshold) << 8);
		*(SieveArray + 1) += (short)((logar2 - threshold) << 8);
	}
#if DEBUG_SIQS
	if (debug_ctr < debugLimit) {
		std::cout << "sieveArray ";
		for (int i = 0; i < 10; i++) {
			std::cout << SieveArray[i] << ", ";
		}
		std::cout << '\n';
	}
#endif

	F2 = 2;
	index = 2;
	for (;;) {
		rowPrimeSieveData = primeSieveData + index;
		currentPrime = rowPrimeSieveData->value;
		F3 = F2 * currentPrime;
		if (X1 + 1 < F3) {
			F3 = X1 + 1;
		}
		F4 = F2;
		while (F4 * 2 <= F3) {
			memcpy(SieveArray + F4, SieveArray, F4 * sizeof(*SieveArray));
			F4 *= 2;
		}

		memcpy(SieveArray + F4, SieveArray, (F3 - F4) * sizeof(*SieveArray));

		if (F3 == X1 + 1) {
			break;  /* <<<<<<<<<<< only exit from for loop */
		}

		F1 = currentPrime;
		logPrimeEvenPoly = 1;
		while (F1 >= 5) {
			F1 /= 3;
			logPrimeEvenPoly++;
		}
		logPrimeOddPoly = (short)(logPrimeEvenPoly << 8);
		F1 = polyadd ? -rowPrimeSieveData->Bainv2[indexFactorA] :
			rowPrimeSieveData->Bainv2[indexFactorA];
		index2 = (rowPrimeSieveData->soln1 + F1) % currentPrime;
		rowPrimeSieveData->soln1 = index2 += currentPrime & (index2 >> BITS_PER_INT_GROUP);   // soln1 is modified

		for (; index2 < F3; index2 += currentPrime)
		{
			*(SieveArray + index2) += logPrimeEvenPoly;
		}

		for (index2 = (rowPrimeSieveData->soln1 + currentPrime -
			rowPrimeSieveData->Bainv2_0) % currentPrime;
			index2 < F3;
			index2 += currentPrime)
		{
			*(SieveArray + index2) += logPrimeOddPoly;
		}

		if (currentPrime != multiplier) {
			for (F1 = index2 = (rowPrimeSieveData->soln1 + currentPrime -
				rowPrimeSieveData->difsoln) % currentPrime;
				index2 < F3;
				index2 += currentPrime)
			{
				*(SieveArray + index2) += logPrimeEvenPoly;
			}
			for (index2 = (F1 + currentPrime -
				rowPrimeSieveData->Bainv2_0) % currentPrime;
				index2 < F3;
				index2 += currentPrime)
			{
				*(SieveArray + index2) += logPrimeOddPoly;
			}
		}
		index++;
		F2 *= currentPrime;
	}

	F1 = (primeSieveData + smallPrimeUpperLimit)->value;
	logPrimeEvenPoly = 1;
	logPrimeOddPoly = 0x100;
	mask = 5;
	while (F1 >= 5) {
		F1 /= 3;
		logPrimeEvenPoly++;
		logPrimeOddPoly += 0x100;
		mask *= 3;
	}
	if (polyadd) {
		for (; index < smallPrimeUpperLimit; index++) 	{
			rowPrimeSieveData = primeSieveData + index;
			currentPrime = rowPrimeSieveData->value;
			if ((S1 = rowPrimeSieveData->soln1 -
				rowPrimeSieveData->Bainv2[indexFactorA]) < 0)
			{
				S1 += currentPrime;
			}
			rowPrimeSieveData->soln1 = S1;   // soln1 is modified
		}

		for (index = smallPrimeUpperLimit; index < firstLimit; index++) {
			rowPrimeSieveData = primeSieveData + index;
			currentPrime = rowPrimeSieveData->value;
			if (currentPrime >= mask) {
				mask *= 3;
				logPrimeEvenPoly++;
				logPrimeOddPoly += 0x100;
			}
			F2 = currentPrime + currentPrime;
			F3 = F2 + currentPrime;
			F4 = F3 + currentPrime;
			S1 = rowPrimeSieveData->soln1 -
				rowPrimeSieveData->Bainv2[indexFactorA];
			rowPrimeSieveData->soln1 = S1 += (S1 >> BITS_PER_INT_GROUP) & currentPrime;
			index2 = X1 / F4 * F4 + S1;
			G0 = -rowPrimeSieveData->difsoln;
			if (S1 + G0 < 0) {
				G0 += currentPrime;
			}
			G1 = G0 + currentPrime;
			G2 = G1 + currentPrime;
			G3 = G2 + currentPrime;
			H0 = -rowPrimeSieveData->Bainv2_0;
			if (S1 + H0 < 0) {
				H0 += currentPrime;
			}
			H1 = H0 + currentPrime;
			H2 = H1 + currentPrime;
			H3 = H2 + currentPrime;
			I0 = H0 - rowPrimeSieveData->difsoln;
			if (S1 + I0 < 0) {
				I0 += currentPrime;
			}
			I1 = I0 + currentPrime;
			I2 = I1 + currentPrime;
			I3 = I2 + currentPrime;

			do 	{
				*(SieveArray + index2) += logPrimeEvenPoly;
				*(SieveArray + index2 + currentPrime) += logPrimeEvenPoly;
				*(SieveArray + index2 + F2) += logPrimeEvenPoly;
				*(SieveArray + index2 + F3) += logPrimeEvenPoly;
				*(SieveArray + index2 + G0) += logPrimeEvenPoly;
				*(SieveArray + index2 + G1) += logPrimeEvenPoly;
				*(SieveArray + index2 + G2) += logPrimeEvenPoly;
				*(SieveArray + index2 + G3) += logPrimeEvenPoly;
				*(SieveArray + index2 + H0) += logPrimeOddPoly;
				*(SieveArray + index2 + H1) += logPrimeOddPoly;
				*(SieveArray + index2 + H2) += logPrimeOddPoly;
				*(SieveArray + index2 + H3) += logPrimeOddPoly;
				*(SieveArray + index2 + I0) += logPrimeOddPoly;
				*(SieveArray + index2 + I1) += logPrimeOddPoly;
				*(SieveArray + index2 + I2) += logPrimeOddPoly;
				*(SieveArray + index2 + I3) += logPrimeOddPoly;
			} while ((index2 -= F4) >= 0);
		}

		for (; index < secondLimit; index++) {
			rowPrimeSieveData = primeSieveData + index;
			currentPrime = rowPrimeSieveData->value;
			F2 = currentPrime + currentPrime;
			F3 = F2 + currentPrime;
			F4 = F2 + F2;
			X2 = X1 - F4;
			if (currentPrime >= mask) {
				mask *= 3;
				logPrimeEvenPoly++;
				logPrimeOddPoly += 0x100;
			}
			if (rowPrimeSieveData->difsoln >= 0) {
				F1 = rowPrimeSieveData->soln1 -
					rowPrimeSieveData->Bainv2[indexFactorA];
				F1 += (F1 >> BITS_PER_INT_GROUP) & currentPrime;
				index2 = (rowPrimeSieveData->soln1 = F1);
				do {
					*(SieveArray + index2) += logPrimeEvenPoly;
					*(SieveArray + index2 + currentPrime) += logPrimeEvenPoly;
					*(SieveArray + index2 + F2) += logPrimeEvenPoly;
					*(SieveArray + index2 + F3) += logPrimeEvenPoly;
				} while ((index2 += F4) <= X2);

				for (; index2 <= X1; index2 += currentPrime) {
					*(SieveArray + index2) += logPrimeEvenPoly;
				}
				index2 = F1 - rowPrimeSieveData->Bainv2_0;
				index2 += (index2 >> BITS_PER_INT_GROUP) & currentPrime;

				do 	{
					*(SieveArray + index2) += logPrimeOddPoly;
					*(SieveArray + index2 + currentPrime) += logPrimeOddPoly;
					*(SieveArray + index2 + F2) += logPrimeOddPoly;
					*(SieveArray + index2 + F3) += logPrimeOddPoly;
				} while ((index2 += F4) <= X2);

				for (; index2 <= X1; index2 += currentPrime) {
					*(SieveArray + index2) += logPrimeOddPoly;
				}

				F1 -= rowPrimeSieveData->difsoln;
				F1 += (F1 >> BITS_PER_INT_GROUP) & currentPrime;
				index2 = F1;

				do {
					*(SieveArray + index2) += logPrimeEvenPoly;
					*(SieveArray + index2 + currentPrime) += logPrimeEvenPoly;
					*(SieveArray + index2 + F2) += logPrimeEvenPoly;
					*(SieveArray + index2 + F3) += logPrimeEvenPoly;
				} while ((index2 += F4) <= X2);

				for (; index2 <= X1; index2 += currentPrime) {
					*(SieveArray + index2) += logPrimeEvenPoly;
				}

				index2 = F1 - rowPrimeSieveData->Bainv2_0;
				index2 += (index2 >> BITS_PER_INT_GROUP) & currentPrime;

				do {
					*(SieveArray + index2) += logPrimeOddPoly;
					*(SieveArray + index2 + currentPrime) += logPrimeOddPoly;
					*(SieveArray + index2 + F2) += logPrimeOddPoly;
					*(SieveArray + index2 + F3) += logPrimeOddPoly;
				} while ((index2 += F4) <= X2);

				for (; index2 <= X1; index2 += currentPrime) {
					*(SieveArray + index2) += logPrimeOddPoly;
				}
			}
		}

		for (; index < thirdLimit; index++) {
			rowPrimeSieveData = primeSieveData + index;
			currentPrime = rowPrimeSieveData->value;
			if (currentPrime >= mask) {
				mask *= 3;
				logPrimeEvenPoly++;
				logPrimeOddPoly += 0x100;
			}
			F2 = rowPrimeSieveData->soln1 - rowPrimeSieveData->Bainv2[indexFactorA];
			F2 += currentPrime & (F2 >> BITS_PER_INT_GROUP);
			index2 = (rowPrimeSieveData->soln1 = F2);

			do {
				*(SieveArray + index2) += logPrimeEvenPoly;
			} while ((index2 += currentPrime) <= X1);

			F1 = F2 - (F3 = rowPrimeSieveData->Bainv2_0);
			F1 += currentPrime & (F1 >> BITS_PER_INT_GROUP);

			do {
				*(SieveArray + F1) += logPrimeOddPoly;
			} while ((F1 += currentPrime) <= X1);

			F2 -= rowPrimeSieveData->difsoln;
			index2 = F2 += currentPrime & (F2 >> BITS_PER_INT_GROUP);

			do {
				*(SieveArray + index2) += logPrimeEvenPoly;
			} while ((index2 += currentPrime) <= X1);

			F2 += (currentPrime & ((F2 - F3) >> BITS_PER_INT_GROUP)) - F3;

			do {
				*(SieveArray + F2) += logPrimeOddPoly;
			} while ((F2 += currentPrime) <= X1);
		}

		for (; index < nbrPrimes2; index++) {
			rowPrimeSieveData = primeSieveData + index;
			currentPrime = rowPrimeSieveData->value;
			if (currentPrime >= mask) {
				mask *= 3;
				logPrimeEvenPoly++;
				logPrimeOddPoly += 0x100;
			}
			F2 = rowPrimeSieveData->soln1 - rowPrimeSieveData->Bainv2[indexFactorA];
			if ((rowPrimeSieveData->soln1 = 
				(F2 += currentPrime & (F2 >> BITS_PER_INT_GROUP))) < X1) {
				*(SieveArray + F2) += logPrimeEvenPoly;
			}

			F1 = F2 - (F3 = rowPrimeSieveData->Bainv2_0);

			if ((F1 += currentPrime & (F1 >> BITS_PER_INT_GROUP)) < X1) {
				*(SieveArray + F1) += logPrimeOddPoly;
			}
			F2 -= rowPrimeSieveData->difsoln;

			if ((F2 += currentPrime & (F2 >> BITS_PER_INT_GROUP)) < X1) {
				*(SieveArray + F2) += logPrimeEvenPoly;
			}

			if ((F2 += (currentPrime & ((F2 - F3) >> BITS_PER_INT_GROUP)) - F3) < X1) {
				*(SieveArray + F2) += logPrimeOddPoly;
			}

			rowPrimeSieveData = primeSieveData + ++index;
			currentPrime = rowPrimeSieveData->value;
			F2 = rowPrimeSieveData->soln1 - rowPrimeSieveData->Bainv2[indexFactorA];

			if ((rowPrimeSieveData->soln1 = (F2 += currentPrime & (F2 >> BITS_PER_INT_GROUP))) < X1) {
				*(SieveArray + F2) += logPrimeEvenPoly;
			}

			F1 = F2 - (F3 = rowPrimeSieveData->Bainv2_0);

			if ((F1 += currentPrime & (F1 >> BITS_PER_INT_GROUP)) < X1) {
				*(SieveArray + F1) += logPrimeOddPoly;
			}

			F2 -= rowPrimeSieveData->difsoln;

			if ((F2 += currentPrime & (F2 >> BITS_PER_INT_GROUP)) < X1) {
				*(SieveArray + F2) += logPrimeEvenPoly;
			}

			if ((F2 += (currentPrime & ((F2 - F3) >> BITS_PER_INT_GROUP)) - F3) < X1) {
				*(SieveArray + F2) += logPrimeOddPoly;
			}

			rowPrimeSieveData = primeSieveData + ++index;
			currentPrime = rowPrimeSieveData->value;
			F2 = rowPrimeSieveData->soln1 - rowPrimeSieveData->Bainv2[indexFactorA];

			if ((rowPrimeSieveData->soln1 = (F2 += currentPrime & (F2 >> BITS_PER_INT_GROUP))) < X1) {
				*(SieveArray + F2) += logPrimeEvenPoly;
			}

			F1 = F2 - (F3 = rowPrimeSieveData->Bainv2_0);

			if ((F1 += currentPrime & (F1 >> BITS_PER_INT_GROUP)) < X1) {
				*(SieveArray + F1) += logPrimeOddPoly;
			}

			F2 -= rowPrimeSieveData->difsoln;

			if ((F2 += currentPrime & (F2 >> BITS_PER_INT_GROUP)) < X1)
			{
				*(SieveArray + F2) += logPrimeEvenPoly;
			}

			F2 -= F3; 

			if ((F2 += currentPrime & (F2 >> BITS_PER_INT_GROUP)) < X1) {
				*(SieveArray + F2) += logPrimeOddPoly;
			}

			rowPrimeSieveData = primeSieveData + ++index;
			currentPrime = rowPrimeSieveData->value;
			F2 = rowPrimeSieveData->soln1 - rowPrimeSieveData->Bainv2[indexFactorA];

			if ((rowPrimeSieveData->soln1 = (F2 += currentPrime & (F2 >> BITS_PER_INT_GROUP))) < X1) {
				*(SieveArray + F2) += logPrimeEvenPoly;
			}

			F1 = F2 - (F3 = rowPrimeSieveData->Bainv2_0);

			if ((F1 += currentPrime & (F1 >> BITS_PER_INT_GROUP)) < X1) {
				*(SieveArray + F1) += logPrimeOddPoly;
			}

			F2 -= rowPrimeSieveData->difsoln;

			if ((F2 += currentPrime & (F2 >> BITS_PER_INT_GROUP)) < X1) {
				*(SieveArray + F2) += logPrimeEvenPoly;
			}

			F2 -= F3;

			if ((F2 += currentPrime & (F2 >> BITS_PER_INT_GROUP)) < X1) {
				*(SieveArray + F2) += logPrimeOddPoly;
			}
		}

		for (; index < nbrFactorBasePrimes; index++) {
			rowPrimeSieveData = primeSieveData + index;
			currentPrime = rowPrimeSieveData->value;
			if (currentPrime >= mask) {
				mask *= 3;
				logPrimeEvenPoly++;
				logPrimeOddPoly += 0x100;
			}

			F2 = rowPrimeSieveData->soln1 - rowPrimeSieveData->Bainv2[indexFactorA];

			if ((rowPrimeSieveData->soln1 = (F2 += currentPrime & (F2 >> BITS_PER_INT_GROUP))) < X1) {
				*(SieveArray + F2) += logPrimeEvenPoly;
			}
			F1 = F2 - (F3 = rowPrimeSieveData->Bainv2_0);

			if ((F1 += currentPrime & (F1 >> BITS_PER_INT_GROUP)) < X1) {
				*(SieveArray + F1) += logPrimeOddPoly;
			}
			F2 -= rowPrimeSieveData->difsoln;

			if ((F2 += currentPrime & (F2 >> BITS_PER_INT_GROUP)) < X1) {
				*(SieveArray + F2) += logPrimeEvenPoly;
			}
			F2 -= F3;

			if ((F2 += currentPrime & (F2 >> BITS_PER_INT_GROUP)) < X1) {
				*(SieveArray + F2) += logPrimeOddPoly;
			}
		}
	}

	/* polyadd is false */
	else {
		for (; index < smallPrimeUpperLimit; index++) {
			rowPrimeSieveData = primeSieveData + index;
			currentPrime = rowPrimeSieveData->value;
			S1 = rowPrimeSieveData->soln1 +
				rowPrimeSieveData->Bainv2[indexFactorA] - currentPrime;
			S1 += currentPrime & (S1 >> BITS_PER_INT_GROUP);
			rowPrimeSieveData->soln1 = S1;
		}

		for (index = smallPrimeUpperLimit; index < firstLimit; index++) {
			rowPrimeSieveData = primeSieveData + index;
			currentPrime = rowPrimeSieveData->value;

			if (currentPrime >= mask) {
				mask *= 3;
				logPrimeEvenPoly++;
				logPrimeOddPoly += 0x100;
			}
			F2 = currentPrime + currentPrime;
			F3 = F2 + currentPrime;
			F4 = F3 + currentPrime;
			S1 = rowPrimeSieveData->soln1 +
				rowPrimeSieveData->Bainv2[indexFactorA] - currentPrime;
			rowPrimeSieveData->soln1 = S1 += (S1 >> BITS_PER_INT_GROUP) & currentPrime;
			index2 = X1 / F4 * F4 + S1;
			G0 = -rowPrimeSieveData->difsoln;

			if (S1 + G0 < 0) {
				G0 += currentPrime;
			}
			G1 = G0 + currentPrime;
			G2 = G1 + currentPrime;
			G3 = G2 + currentPrime;
			H0 = -rowPrimeSieveData->Bainv2_0;

			if (S1 + H0 < 0) {
				H0 += currentPrime;
			}
			H1 = H0 + currentPrime;
			H2 = H1 + currentPrime;
			H3 = H2 + currentPrime;
			I0 = H0 - rowPrimeSieveData->difsoln;

			if (S1 + I0 < 0) {
				I0 += currentPrime;
			}
			I1 = I0 + currentPrime;
			I2 = I1 + currentPrime;
			I3 = I2 + currentPrime;

			do {
				*(SieveArray + index2) += logPrimeEvenPoly;
				*(SieveArray + index2 + currentPrime) += logPrimeEvenPoly;
				*(SieveArray + index2 + F2) += logPrimeEvenPoly;
				*(SieveArray + index2 + F3) += logPrimeEvenPoly;
				*(SieveArray + index2 + G0) += logPrimeEvenPoly;
				*(SieveArray + index2 + G1) += logPrimeEvenPoly;
				*(SieveArray + index2 + G2) += logPrimeEvenPoly;
				*(SieveArray + index2 + G3) += logPrimeEvenPoly;
				*(SieveArray + index2 + H0) += logPrimeOddPoly;
				*(SieveArray + index2 + H1) += logPrimeOddPoly;
				*(SieveArray + index2 + H2) += logPrimeOddPoly;
				*(SieveArray + index2 + H3) += logPrimeOddPoly;
				*(SieveArray + index2 + I0) += logPrimeOddPoly;
				*(SieveArray + index2 + I1) += logPrimeOddPoly;
				*(SieveArray + index2 + I2) += logPrimeOddPoly;
				*(SieveArray + index2 + I3) += logPrimeOddPoly;
			} while ((index2 -= F4) >= 0);
		}

		for (; index < secondLimit; index++) {
			rowPrimeSieveData = primeSieveData + index;
			currentPrime = rowPrimeSieveData->value;
			F2 = currentPrime + currentPrime;
			F3 = F2 + currentPrime;
			F4 = F2 + F2;
			X2 = X1 - F4;
			if (currentPrime >= mask)
			{
				mask *= 3;
				logPrimeEvenPoly++;
				logPrimeOddPoly += 0x100;
			}
			if (rowPrimeSieveData->difsoln >= 0)
			{
				F1 = rowPrimeSieveData->soln1 +
					rowPrimeSieveData->Bainv2[indexFactorA] - currentPrime;
				F1 += currentPrime & (F1 >> BITS_PER_INT_GROUP);
				index2 = (rowPrimeSieveData->soln1 = F1);
				do
				{
					*(SieveArray + index2) += logPrimeEvenPoly;
					*(SieveArray + index2 + currentPrime) += logPrimeEvenPoly;
					*(SieveArray + index2 + F2) += logPrimeEvenPoly;
					*(SieveArray + index2 + F3) += logPrimeEvenPoly;
				} while ((index2 += F4) <= X2);
				for (; index2 <= X1; index2 += currentPrime)
				{
					*(SieveArray + index2) += logPrimeEvenPoly;
				}
				index2 = F1 - rowPrimeSieveData->Bainv2_0;
				index2 += (index2 >> BITS_PER_INT_GROUP) & currentPrime;
				do
				{
					*(SieveArray + index2) += logPrimeOddPoly;
					*(SieveArray + index2 + currentPrime) += logPrimeOddPoly;
					*(SieveArray + index2 + F2) += logPrimeOddPoly;
					*(SieveArray + index2 + F3) += logPrimeOddPoly;
				} while ((index2 += F4) <= X2);
				for (; index2 <= X1; index2 += currentPrime)
				{
					*(SieveArray + index2) += logPrimeOddPoly;
				}
				F1 -= rowPrimeSieveData->difsoln;
				F1 += (F1 >> BITS_PER_INT_GROUP) & currentPrime;
				index2 = F1;
				do
				{
					*(SieveArray + index2) += logPrimeEvenPoly;
					*(SieveArray + index2 + currentPrime) += logPrimeEvenPoly;
					*(SieveArray + index2 + F2) += logPrimeEvenPoly;
					*(SieveArray + index2 + F3) += logPrimeEvenPoly;
				} while ((index2 += F4) <= X2);
				for (; index2 <= X1; index2 += currentPrime)
				{
					*(SieveArray + index2) += logPrimeEvenPoly;
				}
				index2 = F1 - rowPrimeSieveData->Bainv2_0;
				index2 += (index2 >> BITS_PER_INT_GROUP) & currentPrime;
				do
				{
					*(SieveArray + index2) += logPrimeOddPoly;
					*(SieveArray + index2 + currentPrime) += logPrimeOddPoly;
					*(SieveArray + index2 + F2) += logPrimeOddPoly;
					*(SieveArray + index2 + F3) += logPrimeOddPoly;
				} while ((index2 += F4) <= X2);
				for (; index2 <= X1; index2 += currentPrime)
				{
					*(SieveArray + index2) += logPrimeOddPoly;
				}
			}
		}

		for (; index < thirdLimit; index++) {
			rowPrimeSieveData = primeSieveData + index;
			currentPrime = rowPrimeSieveData->value;
			if (currentPrime >= mask)
			{
				mask *= 3;
				logPrimeEvenPoly++;
				logPrimeOddPoly += 0x100;
			}
			F2 = rowPrimeSieveData->soln1 +
				rowPrimeSieveData->Bainv2[indexFactorA] - currentPrime;
			index2 = (rowPrimeSieveData->soln1 = (F2 += currentPrime & (F2 >> BITS_PER_INT_GROUP)));
			do
			{
				*(SieveArray + index2) += logPrimeEvenPoly;
			} while ((index2 += currentPrime) <= X1);
			F1 = F2 - (F3 = rowPrimeSieveData->Bainv2_0);
			F1 += currentPrime & (F1 >> BITS_PER_INT_GROUP);
			do
			{
				*(SieveArray + F1) += logPrimeOddPoly;
			} while ((F1 += currentPrime) <= X1);
			F2 -= rowPrimeSieveData->difsoln;
			F1 = F2 += currentPrime & (F2 >> BITS_PER_INT_GROUP);
			do
			{
				*(SieveArray + F2) += logPrimeEvenPoly;
			} while ((F2 += currentPrime) <= X1);
			F1 -= F3;
			F1 += currentPrime & (F1 >> BITS_PER_INT_GROUP);
			do
			{
				*(SieveArray + F1) += logPrimeOddPoly;
			} while ((F1 += currentPrime) <= X1);
		}

		for (; index < nbrPrimes2; index++) {
			rowPrimeSieveData = primeSieveData + index;
			currentPrime = rowPrimeSieveData->value;
			if (currentPrime >= mask)
			{
				mask *= 3;
				logPrimeEvenPoly++;
				logPrimeOddPoly += 0x100;
			}
			F2 = rowPrimeSieveData->soln1 +
				rowPrimeSieveData->Bainv2[indexFactorA] - currentPrime;
			if ((rowPrimeSieveData->soln1 = (F2 += currentPrime & (F2 >> BITS_PER_INT_GROUP))) < X1)
			{
				*(SieveArray + F2) += logPrimeEvenPoly;
			}
			F1 = F2 - (F3 = rowPrimeSieveData->Bainv2_0);
			if ((F1 += currentPrime & (F1 >> BITS_PER_INT_GROUP)) < X1)
			{
				*(SieveArray + F1) += logPrimeOddPoly;
			}
			F2 -= rowPrimeSieveData->difsoln;
			if ((F2 += currentPrime & (F2 >> BITS_PER_INT_GROUP)) < X1)
			{
				*(SieveArray + F2) += logPrimeEvenPoly;
			}
			F2 -= F3;
			if ((F2 += currentPrime & (F2 >> BITS_PER_INT_GROUP)) < X1)
			{
				*(SieveArray + F2) += logPrimeOddPoly;
			}
			rowPrimeSieveData = primeSieveData + ++index;
			currentPrime = rowPrimeSieveData->value;
			F2 = rowPrimeSieveData->soln1 +
				rowPrimeSieveData->Bainv2[indexFactorA] - currentPrime;
			if ((rowPrimeSieveData->soln1 = (F2 += currentPrime & (F2 >> BITS_PER_INT_GROUP))) < X1)
			{
				*(SieveArray + F2) += logPrimeEvenPoly;
			}
			F1 = F2 - (F3 = rowPrimeSieveData->Bainv2_0);
			if ((F1 += currentPrime & (F1 >> BITS_PER_INT_GROUP)) < X1)
			{
				*(SieveArray + F1) += logPrimeOddPoly;
			}
			F2 -= rowPrimeSieveData->difsoln;
			if ((F2 += currentPrime & (F2 >> BITS_PER_INT_GROUP)) < X1)
			{
				*(SieveArray + F2) += logPrimeEvenPoly;
			}
			F2 -= F3;
			if ((F2 += currentPrime & (F2 >> BITS_PER_INT_GROUP)) < X1)
			{
				*(SieveArray + F2) += logPrimeOddPoly;
			}
			rowPrimeSieveData = primeSieveData + ++index;
			currentPrime = rowPrimeSieveData->value;
			F2 = rowPrimeSieveData->soln1 +
				rowPrimeSieveData->Bainv2[indexFactorA] - currentPrime;
			if ((rowPrimeSieveData->soln1 = (F2 += currentPrime & (F2 >> BITS_PER_INT_GROUP))) < X1)
			{
				*(SieveArray + F2) += logPrimeEvenPoly;
			}
			F1 = F2 - (F3 = rowPrimeSieveData->Bainv2_0);
			if ((F1 += currentPrime & (F1 >> BITS_PER_INT_GROUP)) < X1)
			{
				*(SieveArray + F1) += logPrimeOddPoly;
			}
			F2 -= rowPrimeSieveData->difsoln;
			if ((F2 += currentPrime & (F2 >> BITS_PER_INT_GROUP)) < X1)
			{
				*(SieveArray + F2) += logPrimeEvenPoly;
			}
			F2 -= F3;
			if ((F2 += currentPrime & (F2 >> BITS_PER_INT_GROUP)) < X1)
			{
				*(SieveArray + F2) += logPrimeOddPoly;
			}
			rowPrimeSieveData = primeSieveData + ++index;
			currentPrime = rowPrimeSieveData->value;
			F2 = rowPrimeSieveData->soln1 +
				rowPrimeSieveData->Bainv2[indexFactorA] - currentPrime;
			if ((rowPrimeSieveData->soln1 = (F2 += currentPrime & (F2 >> BITS_PER_INT_GROUP))) < X1)
			{
				*(SieveArray + F2) += logPrimeEvenPoly;
			}
			F1 = F2 - (F3 = rowPrimeSieveData->Bainv2_0);
			if ((F1 += currentPrime & (F1 >> BITS_PER_INT_GROUP)) < X1)
			{
				*(SieveArray + F1) += logPrimeOddPoly;
			}
			F2 -= rowPrimeSieveData->difsoln;
			if ((F2 += currentPrime & (F2 >> BITS_PER_INT_GROUP)) < X1)
			{
				*(SieveArray + F2) += logPrimeEvenPoly;
			}
			F2 -= F3;
			if ((F2 += currentPrime & (F2 >> BITS_PER_INT_GROUP)) < X1)
			{
				*(SieveArray + F2) += logPrimeOddPoly;
			}
		}

		for (; index < nbrFactorBasePrimes; index++) {
			rowPrimeSieveData = primeSieveData + index;
			currentPrime = rowPrimeSieveData->value;
			if (currentPrime >= mask)
			{
				mask *= 3;
				logPrimeEvenPoly++;
				logPrimeOddPoly += 0x100;
			}
			F2 = rowPrimeSieveData->soln1 +
				rowPrimeSieveData->Bainv2[indexFactorA] - currentPrime;
			if ((rowPrimeSieveData->soln1 = (F2 += currentPrime & (F2 >> BITS_PER_INT_GROUP))) < X1)
			{
				*(SieveArray + F2) += logPrimeEvenPoly;
			}
			F1 = F2 - (F3 = rowPrimeSieveData->Bainv2_0);
			if ((F1 += currentPrime & (F1 >> BITS_PER_INT_GROUP)) < X1)
			{
				*(SieveArray + F1) += logPrimeOddPoly;
			}
			F2 -= rowPrimeSieveData->difsoln;
			if ((F2 += currentPrime & (F2 >> BITS_PER_INT_GROUP)) < X1)
			{
				*(SieveArray + F2) += logPrimeEvenPoly;
			}
			F2 -= F3;
			if ((F2 += currentPrime & (F2 >> BITS_PER_INT_GROUP)) < X1)
			{
				*(SieveArray + F2) += logPrimeOddPoly;
			}
		}
	}
#if DEBUG_SIQS
	if (debug_ctr < debugLimit) {
		std::cout << "primeSieveData[0]" << primeSieveData[0].value << ", " << primeSieveData[0].soln1 << '\n';
		std::cout << "primeSieveData[1]" << primeSieveData[1].value << ", " << primeSieveData[1].soln1 << '\n';
		std::cout << "sievearray ";
		for (int i = 0; i < 20; i++) {
			std::cout << SieveArray[i] << ", ";
		}
		std::cout << '\n';
	}
#endif
}

/* divide biR by Divisor, return remainder in biR, setup flag MostSignficantLimbZero */
static void TrialDivisionSub(int Divisor, const int NumberLengthDividend,
	int biR[], bool &mostSignificantLimbZero) {

	int Dividend;
	double dRem = 0;
	double dLimbMult = (double)(1U << BITS_PER_INT_GROUP);
	double dCurrentPrime = (double)Divisor;
	double dDivid, dQuot;


	// Perform division
	switch (NumberLengthDividend)
	{
	case 7:     // {biR6 - biR0} <- {biR6 - biR0} / divis
		Dividend = biR[6];
		dRem = (double)(Dividend - (biR[6] = Dividend / Divisor) * Divisor);
		/* no break */
	case 6:     // {biR5 - biR0} <- {biR5 - biR0} / divis
		dDivid = (double)biR[5] + dRem*dLimbMult;
		dQuot = floor(dDivid / dCurrentPrime);
		dRem = dDivid - dQuot * dCurrentPrime;
		biR[5] = (int)dQuot;
		/* no break */
	case 5:     // {biR4 - biR0} <- {biR4 - biR0} / divis
		dDivid = (double)biR[4] + dRem*dLimbMult;
		dQuot = floor(dDivid / dCurrentPrime);
		dRem = dDivid - dQuot * dCurrentPrime;
		biR[4] = (int)dQuot;
		/* no break */
	case 4:     // {biR3 - biR0} <- {biR3 - biR0} / divis
		dDivid = (double)biR[3] + dRem*dLimbMult;
		dQuot = floor(dDivid / dCurrentPrime);
		dRem = dDivid - dQuot * dCurrentPrime;
		biR[3] = (int)dQuot;
		/* no break */
	case 3:     // {biR2 - biR0} <- {biR2 - biR0} / divis
		dDivid = (double)biR[2] + dRem*dLimbMult;
		dQuot = floor(dDivid / dCurrentPrime);
		dRem = dDivid - dQuot * dCurrentPrime;
		biR[2] = (int)dQuot;
		/* no break */
	case 2:     // {biR1 - biR0} <- {biR1 - biR0} / divis
		dDivid = (double)biR[1] + dRem*dLimbMult;
		dQuot = floor(dDivid / dCurrentPrime);
		dRem = dDivid - dQuot * dCurrentPrime;
		biR[1] = (int)dQuot;
		dDivid = (double)biR[0] + dRem*dLimbMult;
		biR[0] = (int)floor(dDivid / dCurrentPrime);
	}

	/* check whether most significant limb is zero */
	if (NumberLengthDividend > 2)
		mostSignificantLimbZero = (biR[NumberLengthDividend - 1] == 0);
	else {
		// Criteria is to fit in a double (52 bits).
		mostSignificantLimbZero = (biR[1] < (1 << (52 - BITS_PER_INT_GROUP)));
	}
}

/* set value of Divisor, return true if R is not a multiple of Divisor */
static bool TrialDivisionSubA(const int index,
	const int biR[],
	bool &fullRemainder,
	const bool oddPolynomial,
	const int index2,
	const int NumberLengthDividend,
	int rowMatrixBbeforeMerge[],
	int &Divisor,
	int &expParity,
	int &nbrColumns) {

	int iRem;
	const PrimeSieveData *rowPrimeSieveData = &primeSieveData[index];
	const PrimeTrialDivisionData *rowPrimeTrialDivisionData = &ag->primeTrialDivisionData[index];

	if (fullRemainder == false) {
		Divisor = rowPrimeSieveData->value;
		if (oddPolynomial)
		{
			iRem = index2 - rowPrimeSieveData->soln1 +
				rowPrimeSieveData->Bainv2_0;
		}
		else { // not odd polynomial
			iRem = index2 - rowPrimeSieveData->soln1;
		}

		if (iRem >= Divisor) {  // get remainder into range 0 to Divisor-1
			if ((iRem -= Divisor) >= Divisor)
			{
				if ((iRem -= Divisor) >= Divisor)
				{
					iRem %= Divisor; // if 2 subtractions are not enough use modulus
				}
			}
		}
		else {
			iRem += (iRem >> BITS_PER_INT_GROUP) & Divisor;  // if iRem is -ve add value of Divisor
		}

		if (iRem != 0 && iRem != Divisor - rowPrimeSieveData->difsoln) {
			if (expParity != 0)
			{
				rowMatrixBbeforeMerge[nbrColumns++] = index;
				expParity = 0;
			}
			return true;              // Process next prime.
		}
		fullRemainder = true;
	}
	
	else {  // fullRemainder is true
		Divisor = rowPrimeTrialDivisionData->value;
		double dRem = 0;

		for (int ix = 1; ix < NumberLengthDividend; ix++) {
			dRem += (double)biR[ix] * (double)rowPrimeTrialDivisionData->exp[ix - 1];
		}
		dRem += biR[0];

		double dCurrentPrime = (double)Divisor;
		if (dRem != floor(dRem / dCurrentPrime)*dCurrentPrime)
		{                     // Number is not a multiple of prime.
			if (expParity != 0)
			{
				rowMatrixBbeforeMerge[nbrColumns++] = index;
				expParity = 0;
			}
			return true;              // Process next prime.
		}
	}
	return false;
}

/* note: DoTrialDivision still uses DA-style BigNums */
static int DoTrialDivision(PrimeSieveData *primeSieveData,
	int rowMatrixBbeforeMerge[],		// is altered within this function
	int index2,
	Znum biDividend,					// copy, which is altered, of original value
	int rowSquares[],					// is altered within this function				
	const bool oddPolynomial)
{
	int biR[7] = { 0 }; /* contains value of biDividend, converted 
						Dividend maximum about 65 digits! */
	Znum NextDiv;
	int nbrSquares = rowSquares[0];
	int Divisor;
	int NumberLengthDividend;
	int index, rem;
	int expParity;
	int nbrColumns = rowMatrixBbeforeMerge[0];
	const PrimeSieveData *rowPrimeSieveData;

	expParity = 0;

#if DEBUG_SIQS
	if (debug_ctr < debugLimit) {
		std::cout << "TrialDivision: nbrColumns = " << nbrColumns 
			<< " NumberLength = " << NumberLengthDividend ;
		Bin2Dec((limb *)biR, output, NumberLengthDividend, 0);
		std::cout << " biR =" << output << '\n';
	}
#endif

/*************** new code ****************************************/

	if (biDividend < LLONG_MAX && biDividend > LLONG_MIN) {
		long long Dividend = mpz_get_si(ZT(biDividend));
		for (index = 1; index < nbrFactorBasePrimes; index++) {
			Divisor = ag->primeTrialDivisionData[index].value;
			for (;;) {
				rem = Dividend%Divisor;
				if (rem == 0) {
					Dividend /= Divisor;
					expParity = 1 - expParity;  // flip 0 to 1, or vice versa
					if (expParity == 0) {
						rowSquares[nbrSquares++] = Divisor;  // record squared factor
					}
				}
				else {
					if (expParity != 0) {
						rowMatrixBbeforeMerge[nbrColumns++] = index;
						expParity = 0;
					}
					break;
				}
			}
		}
		rowSquares[0] = nbrSquares;  // store new total number of squares
		rowMatrixBbeforeMerge[0] = nbrColumns;

		if (Dividend > MAX_INT_NBR)
			return 0;  // result is too large to use
		else
			return (int)Dividend;
	}

/*   really simple, works perfectly, but it's much slower than the old version! */
	// same logic as above, but using Znum instead of long long because value of 
	// biDividend won't fit into a 64 bit integer
	//for (index = 1; index < nbrFactorBasePrimes; index++) {
	//	Divisor = ag->primeTrialDivisionData[index].value;
	//	for (;;) {
	//		/* profiling shows that all the CPU time is used by the function call
	//		below, so there is no easy way to speed it up */
	//		rem = (int)mpz_tdiv_q_ui(ZT(NextDiv), ZT(biDividend), Divisor);
	//		/* note: rem is less than Divisor,so overflow cannot occur */
	//		if (rem == 0) {
	//			biDividend = NextDiv;
	//			expParity = 1 - expParity;  // flip 0 to 1, or vice versa
	//			if (expParity == 0) {
	//				rowSquares[nbrSquares++] = Divisor;
	//			}
	//		}
	//		else {
	//			if (expParity != 0) {
	//				rowMatrixBbeforeMerge[nbrColumns++] = index;
	//				expParity = 0;
	//			}
	//			break;
	//		}
	//	}
	//}
	//rowSquares[0] = nbrSquares;
	//rowMatrixBbeforeMerge[0] = nbrColumns;

	//if (biDividend > MAX_INT_NBR)
	//	return 0;  // result is too large to use
	//else
	//	return (int)mpz_get_si(ZT(biDividend));

/*************** old code *****************************************/
	//if (NumberLengthDividend <= 1) { // Dividend has one limb.
	//	for (index = 1; index < nbrFactorBasePrimes; index++) {
	//		Divisor = ag->primeTrialDivisionData[index].value;
	//		while (biR[0] % Divisor == 0)
	//		{
	//			biR[0] /= Divisor;
	//			expParity = 1 - expParity;  // flip 0 to 1, or vice versa
	//			if (expParity == 0)
	//			{
	//				rowSquares[nbrSquares++] = Divisor;
	//			}
	//		}
	//		if (expParity != 0)
	//		{
	//			rowMatrixBbeforeMerge[nbrColumns++] = index;
	//			expParity = 0;
	//		}
	//	}
	//}

	 // Dividend has at least two limbs.
	else {
		bool mostSignificantLimbZero = false;
		int nbr, iRem;
		int left, right, median;
		bool fullRemainder;
		int indexFactorA = 0;
		int newFactorAIndex;
		bool testFactorA = true;
		newFactorAIndex = aindex[0];
		NumberLengthDividend = ZtoBigNbr(biR, biDividend); // convert Dividend to BigNum
		assert(NumberLengthDividend <= 7);

		/* exit for loop when testFactorA set to false*/
		for (index = 1; testFactorA; index++) {
			fullRemainder = false;
			if (index < 3) 	{
				fullRemainder = true;
			}
			else if (index == newFactorAIndex) 	{
				fullRemainder = true;
				if (++indexFactorA == nbrFactorsA) 	{
					testFactorA = false;   // All factors of A were tested.
				}
				else {
					newFactorAIndex = aindex[indexFactorA];
				}
			}

			for (;;) {
				if (TrialDivisionSubA(index, biR, fullRemainder, oddPolynomial, index2,
					NumberLengthDividend, rowMatrixBbeforeMerge, Divisor, 
					expParity, nbrColumns))
					break; // Process next prime.

				expParity = 1 - expParity;
				if (expParity == 0) {
					rowSquares[nbrSquares++] = (int)Divisor;
				}

				// Perform division
				TrialDivisionSub(Divisor, NumberLengthDividend, biR, mostSignificantLimbZero);

				if (mostSignificantLimbZero) {
					NumberLengthDividend--;
					if (NumberLengthDividend == 1) {       // Number fits in a double
						double dDivid = (double)biR[1] * (double)(1U << BITS_PER_INT_GROUP)
							+ (double)biR[0];
						int sqrtDivid = (int)(floor(sqrt(dDivid)));
						fullRemainder = true;
						for (; index < nbrFactorBasePrimes; index++) {
							rowPrimeSieveData = primeSieveData + index;
							Divisor = rowPrimeSieveData->value;

							if (testFactorA && index == newFactorAIndex) {
								fullRemainder = true;
								if (++indexFactorA == nbrFactorsA) {
									testFactorA = false;   // All factors of A were tested.
								}
								else {
									newFactorAIndex = aindex[indexFactorA];
								}
							}

							for (;;) {
								if (fullRemainder == false) {

									if (oddPolynomial) {
										iRem = index2 - rowPrimeSieveData->soln1 +
											rowPrimeSieveData->Bainv2_0;
									}
									else {
										iRem = index2 - rowPrimeSieveData->soln1;
									}

									if (iRem >= Divisor) {
										if ((iRem -= Divisor) >= Divisor) {
											if ((iRem -= Divisor) >= Divisor) {
												iRem %= Divisor;
											}
										}
									}
									else {
										iRem += (iRem >> BITS_PER_INT_GROUP) & Divisor;
									}

									if (iRem != 0 && iRem != Divisor - rowPrimeSieveData->difsoln) {
										break;
									}
									fullRemainder = true;
								}
								else {
									if (dDivid != floor(dDivid / Divisor) * Divisor) {
										break;
									}
								}

								dDivid /= Divisor;
								sqrtDivid = (int)(floor(sqrt(dDivid))) + 1;
								expParity = 1 - expParity;
								if (expParity == 0) {
									rowSquares[nbrSquares++] = (int)Divisor;
								}
							}

							if (expParity != 0) {
								rowMatrixBbeforeMerge[nbrColumns++] = index;
								expParity = 0;
							}
							if (Divisor > sqrtDivid) {             // End of trial division.?
								rowSquares[0] = nbrSquares;
								index = nbrFactorBasePrimes - 1;
								if (dDivid <= ag->primeTrialDivisionData[index].value &&
									dDivid > 1)
								{          // Perform binary search to find the index.
									left = -1;
									median = right = nbrFactorBasePrimes;

									while (left != right) {
										median = ((right - left) >> 1) + left;
										nbr = ag->primeTrialDivisionData[median].value;
										if (nbr < dDivid) {
											if (median == left &&
												congruencesFound >= matrixBLength) {
												return 0;
											}
											left = median;
										}
										else if (nbr > dDivid) {
											right = median;
										}
										else {
											break;
										}
									}

									rowMatrixBbeforeMerge[nbrColumns++] = median;
									rowMatrixBbeforeMerge[0] = nbrColumns;
									return 1;
								}
								rowMatrixBbeforeMerge[0] = nbrColumns;
								if (dDivid > 2.147e9) {
									return 0;  // Discard this congruence because of large cofactor.
								}
								return (int)dDivid;
							}
							fullRemainder = false;
						}
						break;
					}
				}
			}             /* end inner for */
		}               /* end for */

		for (; index < nbrFactorBasePrimes; index++) {
			fullRemainder = false;
			for (;;) {
				if (TrialDivisionSubA(index, biR, fullRemainder, oddPolynomial, index2,
					NumberLengthDividend, rowMatrixBbeforeMerge, Divisor, 
					expParity, nbrColumns))
					break; // Process next prime.

				expParity = 1 - expParity;
				if (expParity == 0) {
					rowSquares[nbrSquares++] = (int)Divisor;
				}

				// Perform division
				TrialDivisionSub(Divisor, NumberLengthDividend, biR, mostSignificantLimbZero);

				if (mostSignificantLimbZero) {
					NumberLengthDividend--;
					if (NumberLengthDividend == 1) {       // Number fits in a double
						double dDivid = (double)biR[1] * (double)(1U << BITS_PER_INT_GROUP)
							+ (double)biR[0];
						int sqrtDivid = (int)(floor(sqrt(dDivid)));
						fullRemainder = true;
						for (; index < nbrFactorBasePrimes; index++) {
							rowPrimeSieveData = primeSieveData + index;
							if (rowPrimeSieveData->value == 41893) {
								// divis = 5;
							}
							Divisor = rowPrimeSieveData->value;
							if (testFactorA && index == newFactorAIndex) {
								fullRemainder = true;
								if (++indexFactorA == nbrFactorsA) {
									testFactorA = false;   // All factors of A were tested.
								}
								else {
									newFactorAIndex = aindex[indexFactorA];
								}
							}

							for (;;) {
								if (fullRemainder == false) {

									if (oddPolynomial) {
										iRem = index2 - rowPrimeSieveData->soln1 +
											rowPrimeSieveData->Bainv2_0;
									}
									else {
										iRem = index2 - rowPrimeSieveData->soln1;
									}
									if (iRem >= Divisor) {
										if ((iRem -= Divisor) >= Divisor) {
											if ((iRem -= Divisor) >= Divisor) {
												iRem %= Divisor;
											}
										}
									}
									else {
										iRem += (iRem >> BITS_PER_INT_GROUP) & Divisor;
									}

									if (iRem != 0 && iRem != Divisor - rowPrimeSieveData->difsoln)
									{
										break;
									}
									fullRemainder = true;
								}
								else {
									if (dDivid != floor(dDivid / Divisor) * Divisor) {
										break;
									}
								}

								dDivid /= Divisor;
								sqrtDivid = (int)floor(sqrt(dDivid));
								expParity = 1 - expParity;
								if (expParity == 0) {
									rowSquares[nbrSquares++] = (int)Divisor;
								}
							}  /* end for */

							if (expParity != 0) {
								rowMatrixBbeforeMerge[nbrColumns++] = index;
								expParity = 0;
							}
							if (Divisor > sqrtDivid) {        // End of trial division?
								if (dDivid >= (double)(1U << BITS_PER_INT_GROUP))
								{                   // Dividend is too large.
									return 0;
								}
								Divisor = (int)dDivid;
								rowSquares[0] = nbrSquares;
								index = nbrFactorBasePrimes - 1;
								if (Divisor <= ag->primeTrialDivisionData[index].value &&
									Divisor > 1)
								{          // Perform binary search to find the index.
									left = -1;
									median = right = nbrFactorBasePrimes;
									while (left != right) {
										median = ((right - left) >> 1) + left;
										nbr = ag->primeTrialDivisionData[median].value;
										if (nbr < Divisor) {
											if (median == left &&
												congruencesFound >= matrixBLength)
											{
												return 0;
											}
											left = median;
										}
										else if (nbr > Divisor) {
											right = median;
										}
										else {
											break;
										}
									}

									rowMatrixBbeforeMerge[nbrColumns++] = median;
									rowMatrixBbeforeMerge[0] = nbrColumns;
									return 1;
								}

								rowMatrixBbeforeMerge[0] = nbrColumns;
								return Divisor;
							}

							fullRemainder = false;
						}  /* end for */

						if (dDivid >= (double)(1U << BITS_PER_INT_GROUP))
						{                   // Dividend is too large.
							return 0;
						}
						biR[0] = (int)dDivid;
						break;
					}
				}
			}             /* end for */
		}               /* end for */
	}

	rowSquares[0] = nbrSquares;	  // store new total number of square
	rowMatrixBbeforeMerge[0] = nbrColumns;
	if (NumberLengthDividend > 1) {
		return 0;           // Very large quotient.
	}
	return biR[0];
}

static void mergeArrays(int aindex[], int nbrFactorsA, int rowMatrixB[], int rowMatrixBeforeMerge[],
	int rowSquares[])
{
	int indexAindex = 0;
	int indexRMBBM = 1;
	int indexRMB = 1;
	int nbrColumns = rowMatrixBeforeMerge[0];

	while (indexAindex < nbrFactorsA && indexRMBBM < nbrColumns)
	{
		if (aindex[indexAindex] < rowMatrixBeforeMerge[indexRMBBM])
		{
			rowMatrixB[indexRMB++] = aindex[indexAindex++];
		}
		else if (aindex[indexAindex] > rowMatrixBeforeMerge[indexRMBBM])
		{
			rowMatrixB[indexRMB++] = rowMatrixBeforeMerge[indexRMBBM++];
		}
		else
		{
			rowSquares[rowSquares[0]++] =
				ag->primeTrialDivisionData[aindex[indexAindex++]].value;
			indexRMBBM++;
		}
	}
	while (indexAindex < nbrFactorsA)
	{
		rowMatrixB[indexRMB++] = aindex[indexAindex++];
	}
	while (indexRMBBM < nbrColumns)
	{
		rowMatrixB[indexRMB++] = rowMatrixBeforeMerge[indexRMBBM++];
	}
	rowMatrixB[LENGTH_OFFSET] = indexRMB;
}

static bool SmoothRelationFound(
	unsigned char positive,
	int *rowMatrixB, int *rowMatrixBbeforeMerge,
	int index2,
	int *rowSquares, Znum &biLinearCoeff,
	Znum &biT, Znum &biU,
	Znum &biR, const bool oddPolynomial)
{
	int index;
	int nbrSquares;
	if (congruencesFound >= matrixBLength)
	{
		return true;            // All congruences already found.
	}
	// Add all elements of aindex array to the rowMatrixB array discarding
	// duplicates.
	mergeArrays(aindex, nbrFactorsA, rowMatrixB, rowMatrixBbeforeMerge, rowSquares);
	nbrSquares = rowSquares[0];
	biR = 1; 
	biT = positive ? 1 : -1;
	biU = biQuadrCoeff * (index2 - SieveLimit);        // Ax
	biU += biLinearCoeff;                              // Ax+B
	if (oddPolynomial) 	{
		biU -= biLinearDelta[0];    // Ax+B (odd)
		biU -= biLinearDelta[0];    // Ax+B (odd)
	}
	
	biU = abs(biU);   // If number is negative make it positive.

	for (index = 1; index < nbrSquares; index++)
	{
		int D = rowSquares[index];
		if (D == multiplier) 	{
			biU += zModulus; 
			biU /= D; 
		}
		else {
			biR *= D; 
		}
	}
	if (InsertNewRelation(rowMatrixB, biT, biU, biR))
	{
		smoothsFound++;
		ShowSIQSStatus(rowMatrixB);
		return true;
	}
	return false;
}

/* return true if relation added to matrixB,
or if matrixB is full, otherwise false*/
static bool PartialRelationFound(
	unsigned char positive,
	int *rowMatrixB, int *rowMatrixBbeforeMerge,
	int index2,
	int Divid, int *rowPartials,
	int *rowSquares,
	Znum &biLinearCoeff, Znum &biT,
	Znum &biR, Znum &biU, Znum &biV,
	int *indexFactorsA, const bool oddPolynomial)
{
	int index;
	int expParity;
	int D, Divisor;
	int nbrFactorsPartial;
	int prev;
	unsigned int seed;
	int hashIndex;
	//int *rowPartial;
	struct MatPart *rowPartial;
	int newDivid = (int)Divid;    // This number is greater than zero.
	int indexFactorA = 0;
	int nbrSquares;
	int nbrColumns;

	if (congruencesFound >= matrixBLength) 	{
		return true;
	}
	// Partial relation found.
	totalPartials++;
	// Check if there is already another relation with the same
	// factor outside the prime base.
	// Calculate hash index
	hashIndex = matrixPartialHashIndex[(int)(Divid & 0xFFE) >> 1];
	prev = -1;
	while (hashIndex >= 0) {
		int oldDivid;

		rowPartial = &matrixPartial[hashIndex];
		oldDivid = rowPartial->Divid;
		if (newDivid == oldDivid || newDivid == -oldDivid)
		{   // Match of partials.
			long long Rem;

			biV = rowPartial->value; // biV = Old positive square root (Ax+B).
			seed = rowPartial->seed;
			getFactorsOfA(seed, indexFactorsA);
			biR = newDivid;
			nbrFactorsPartial = 0;
			
			biT = biV * biV;     // biT = old (Ax+B)^2.
			
			biT -= zModulus;     // biT = old (Ax+B)^2 - N.
			if (oldDivid < 0) {
				rowPartials[nbrFactorsPartial++] = 0; // Insert -1 as a factor.
			}
			biT = abs(biT);   // Make it positive

			// The number is multiple of the big prime, so divide by it.
			biT /= newDivid; 
			for (index = 0; index < nbrFactorsA; index++) {
				biT /= ag->primeTrialDivisionData[indexFactorsA[index]].value;
			}

			for (index = 1; index < nbrFactorBasePrimes; index++) {
				expParity = 0;
				if (index >= indexMinFactorA && indexFactorA < nbrFactorsA) {
					if (index == indexFactorsA[indexFactorA]) {
						expParity = 1;
						indexFactorA++;
					}
				}

				const auto rowPrimeTrialDivisionData = &ag->primeTrialDivisionData[index];
				Divisor = rowPrimeTrialDivisionData->value;

				for (;;) {
					Rem = mpz_tdiv_ui(ZT(biT), Divisor);  // faster??
					if (Rem != 0) {
						break;
					}
					expParity = 1 - expParity;
					biT /= Divisor; 

					if (expParity == 0) {
						rowSquares[rowSquares[0]++] = (int)Divisor;
					}
					if (biT == 1)  	{    
						break; // biT = 1, so division has ended.
					}
				}
				if (expParity != 0)
				{
					rowPartials[nbrFactorsPartial++] = index;
				}
			}
			biT = biQuadrCoeff * (index2 - SieveLimit);
			biT += biLinearCoeff;        // biT = Ax+B
			if (oddPolynomial) {         // Ax+B (odd)
				biT -= biLinearDelta[0];
				biT -= biLinearDelta[0];
			}
			
			biT = abs(biT); // If number is negative make it positive.
			// biU = Product of old Ax+B times new Ax+B
			MultZnumModN(biV, biT, biU, zModulus);
			// Add all elements of aindex array to the rowMatrixB array discarding
			// duplicates.
			mergeArrays(aindex, nbrFactorsA, rowMatrixB, rowMatrixBbeforeMerge, rowSquares);
			rowMatrixBbeforeMerge[0] = nbrColumns = rowMatrixB[LENGTH_OFFSET];
			memcpy(&rowMatrixBbeforeMerge[1], &rowMatrixB[1], nbrColumns * sizeof(int));
			mergeArrays(rowPartials, nbrFactorsPartial, rowMatrixB, rowMatrixBbeforeMerge, rowSquares);
			nbrSquares = rowSquares[0];
			for (index = 1; index < nbrSquares; index++)
			{
				D = rowSquares[index];
				if (D != multiplier)  {
					biR *= D; 
				}
				else	{
					biU += zModulus;   
					biU /= multiplier; 
				}
			}
			if (rowMatrixB[0] > 1 &&
				InsertNewRelation(rowMatrixB, biT, biU, biR))
			{
				partialsFound++;
				ShowSIQSStatus(rowMatrixB);
				return true;
			}
			return false;
		}
		else
		{
			prev = hashIndex;
			hashIndex = rowPartial->hashindex; // Get next index for same hash.
		}
	} /* end while */

	  //  synchronized(firstPrimeSieveData)
	{
		if (hashIndex == -1 && nbrPartials < MAX_PRIMES * 8)
		{ // No match and partials table is not full.
		  // Add partial to table of partials.
			if (prev >= 0)
			{
				matrixPartial[prev].hashindex = nbrPartials;
			}
			else
			{
				matrixPartialHashIndex[(newDivid & 0xFFE) >> 1] = nbrPartials;
			}
			rowPartial = &matrixPartial[nbrPartials];
			// Add all elements of aindex array to the rowMatrixB array discarding
			// duplicates.
			mergeArrays(aindex, nbrFactorsA, rowMatrixB, rowMatrixBbeforeMerge, rowSquares);
			biR = Divid; 
			nbrSquares = rowSquares[0];
			for (index = 1; index < nbrSquares; index++)
			{
				D = rowSquares[index];
				MultZnumByIntModN(biR, D, biR, zModulus);
				if (D == multiplier)	{
					biU /= D;
				}
			}
			rowPartial->Divid = (positive ? newDivid : -newDivid);
			// Indicate last index with this hash.
			rowPartial->hashindex = -1;
			biT = biQuadrCoeff * (index2 - SieveLimit);
			biT += biLinearCoeff;     // biT = Ax+B
			if (oddPolynomial)	{      // Ax+B (odd)
				biT -= biLinearDelta[0];
				biT -= biLinearDelta[0];
			}
			
			biT = abs(biT);
			
			rowPartial->value = biT;
			rowPartial->seed = (int)oldSeed;
			nbrPartials++;
		}
	}               // End synchronized block.
	return false;
}

static void SieveLocationHit(int rowMatrixB[], int rowMatrixBbeforeMerge[],
	int index2,
	PrimeSieveData *primeSieveData,
	int rowPartials[],
	int rowSquares[], Znum &biDividend,
	Znum &biT, Znum &biLinearCoeff,
	Znum &biR, Znum &biU, Znum &biV,
	int indexFactorsA[], const bool oddPolynomial)
{
	unsigned char positive;
	int index;
	int Divid;
	int nbrColumns;

	trialDivisions++;
	if (trialDivisions >= 496) 	{
		trialDivisions = 496;
	}
	// calculate Ax
	biT = biQuadrCoeff * (index2 - SieveLimit);
	biT += biLinearCoeff;          // Ax+B

	if (oddPolynomial) {                // Ax+B (odd)
		biT -= biLinearDelta[0];
		biT -= biLinearDelta[0];
	}
#if DEBUG_SIQS
	if (debug_ctr < debugLimit) {
		std::cout << "SieveLocationHit Bit = " << biT;
		debug_ctr++;
	}
#endif
	biDividend = biT * biT;       // (Ax+B)^2
			// To factor: (Ax+B)^2-N
	biDividend -= zModulus; 
	/* factor biDividend */

	positive = true;
	if (biDividend < 0) {
		positive = false;
		biDividend = - biDividend;   // change sign to +ve
	}

#if DEBUG_SIQS
	if (debug_ctr < debugLimit) {
		std::cout << " biDividend = " << biDividend << "\n";
	}
#endif

	rowSquares[0] = 1;
	for (index = 0; index < nbrFactorsA; index++) {
		biDividend /= afact[index]; 
	}
	nbrColumns = 1;
	if (!positive) 	{                      // Insert -1 as a factor.
		rowMatrixBbeforeMerge[nbrColumns++] = 0;
	}
	rowMatrixBbeforeMerge[0] = nbrColumns;
	Divid = DoTrialDivision(primeSieveData,
		rowMatrixBbeforeMerge,
		index2, biDividend,
		rowSquares, oddPolynomial);

#if DEBUG_SIQS
	if (debug_ctr < debugLimit) {
		std::cout << " Divid = " << Divid << "\n";
	}
#endif

	if (Divid == 1) { // Smooth relation found.
		SmoothRelationFound(positive, rowMatrixB,
			rowMatrixBbeforeMerge,
			index2,
			rowSquares,
			biLinearCoeff, biT, biU, biR,
			oddPolynomial);
	}

	else {
		if (Divid > 0 && Divid < largePrimeUpperBound) {
			PartialRelationFound(positive, rowMatrixB,
				rowMatrixBbeforeMerge,
				index2,
				Divid, rowPartials,
				rowSquares, biLinearCoeff,
				biT, biR, biU, biV,
				indexFactorsA, oddPolynomial);
		}
	}
	return;
}

static unsigned int getFactorsOfA(unsigned int seed, int *indexA)
{
	int index, index2, i, tmp;
	for (index = 0; index < nbrFactorsA; index++) {
		do {
			seed = 1141592621 * seed + 321435;
			i = (int)(((double)seed * (double)span) / (double)0x100000000ll + indexMinFactorA);
			for (index2 = 0; index2 < index; index2++) {
				if (indexA[index2] == i || indexA[index2] == i + 1) {
					break;
				}
			}
		} while (index2 < index);

		indexA[index] = i;
	}

	for (index = 0; index<nbrFactorsA; index++) {   // Sort factors of A.
		for (index2 = index + 1; index2<nbrFactorsA; index2++) {
			if (indexA[index] > indexA[index2]) {
				tmp = indexA[index];
				indexA[index] = indexA[index2];
				indexA[index2] = tmp;
			}
		}
	}
	return seed;
}

/************************************************************************/
/* Multithread procedure:                                               */
/*                                                                      */
/* 1) Main thread generates factor base and other parameters.           */
/* 2) Start M threads where the number M is specified by the user in a  */
/*    box beneath the applet.                                           */
/* 3) For each polynomial:                                              */
/*    3a) Main thread generates the data for the set of 2^n polynomials.*/
/*    3b) Each child thread computes a range of polynomials             */
/*        (u*2^n/M to (u+1)*2^n/M exclusive).                           */
/* Partial and full relation routines must be synchronized.             */           
/* Note this is the only function in this file that is not 'static' i.e.*/
/* it is the only entry point                                           */
/************************************************************************/
void FactoringSIQS(const Znum &zN, Znum &Factor) {
	int FactorBase;
	int currentPrime;
	int NbrMod;

	PrimeSieveData *rowPrimeSieveData;  // elements value and modsqrt are modified via this pointer
	PrimeTrialDivisionData *rowPrimeTrialDivisionData;  // elements value and exp are modified via this pointer
	int Power2, SqrRootMod, fact;
	int D, E, Q, V, W, X, Y, Z, T1, V1, W1, Y1;
	double Temp, Prod;
	double bestadjust;
	int i, j;
	const int arrmult[] = { 1, 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43,
		47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97 };  // 1st 26 primes
	double adjustment[sizeof(arrmult) / sizeof(arrmult[0])];
	double dNumberToFactor, dlogNumberToFactor;

	first = true;
	if (!GetConsoleScreenBufferInfo(hConsole, &csbi))
	{
		fprintf(stderr, "** GetConsoleScreenBufferInfo failed with %d!\n", GetLastError());
		Beep(750, 1000);
	}

	coordScreen.X = csbi.dwCursorPosition.X;  // save cursor co-ordinates
	coordScreen.Y = csbi.dwCursorPosition.Y;

	/* initialise variables */
	nbrThreadFinishedPolySet = 0;
	trialDivisions = 0;
	smoothsFound = 0;
	totalPartials = 0;
	partialsFound = 0;
	ValuesSieved = 0;
	congruencesFound = 0;
	polynomialsSieved = 0;
	nbrPartials = 0;
	newSeed = 0;
	nbrPrimesUsed = 0;

	/* allocate storage for arrays defined as vectors */
	matrixB.resize(MAX_PRIMES + 50);
	for (int i = 0; i < MAX_PRIMES + 50; i++)
		matrixB[i].resize(MAX_FACTORS_RELATION);
	matrixPartial.resize(MAX_PRIMES * 8);
	primeSieveData.resize(MAX_PRIMES + 3);
	/* allocate storage for ther arrays */
	ag = (AG*)calloc(1, sizeof(AG));
	assert(ag != nullptr);
	memset(ag->primesUsed, 0, MAX_PRIMES * sizeof(char));

	//  threadArray = new Thread[numberThreads];
	Temp = logZnum(zN);
	/* increased number of primes on 17/6/2019 in line with DA's calculator */
	nbrFactorBasePrimes = (int)exp(sqrt(Temp * log(Temp)) * 0.363 - 1);
	if (nbrFactorBasePrimes > MAX_PRIMES)
		nbrFactorBasePrimes = MAX_PRIMES;

	/* SieveLimit is rounded down to a multiple of 8.
	if zN = 10^95 SieveLimit is about 130,768. */
	SieveLimit = (int)exp(8.5 + 0.015 * Temp) & 0xFFFFFFF8;
	if (SieveLimit > MAX_SIEVE_LIMIT)
	{
		SieveLimit = MAX_SIEVE_LIMIT;
	}
	nbrFactorsA = (int)(Temp*0.051 + 1);
	NbrPolynomials = (1 << (nbrFactorsA - 1)) - 1;

	zModulus = zN;
	zTestNbr2 = zModulus; // = N
	memset(matrixPartialHashIndex, 0xFF, sizeof(matrixPartialHashIndex));

#ifdef __EMSCRIPTEN__
	InitSIQSStrings(SieveLimit, nbrFactorBasePrimes);
	startSieveTenths = (int)(tenths() - originalTenthSecond);
#endif

	/************************/
	/* Compute startup data */
	/************************/
	/* search for best Knuth-Schroeppel multiplier */
	bestadjust = -10.0e0;
	primeSieveData[0].value = 1;
	ag->primeTrialDivisionData[0].value = 1;
	rowPrimeSieveData = &primeSieveData[1];
	rowPrimeTrialDivisionData = &ag->primeTrialDivisionData[1];
	rowPrimeSieveData->value = 2;
	rowPrimeTrialDivisionData->value = 2;
	// (2^BITS_PER_INT_GROUP)^(j+1) mod 2
	rowPrimeTrialDivisionData->exp[0] = rowPrimeTrialDivisionData->exp[1] =
		rowPrimeTrialDivisionData->exp[2] = rowPrimeTrialDivisionData->exp[3] =
		rowPrimeTrialDivisionData->exp[4] = rowPrimeTrialDivisionData->exp[5] = 0;

	NbrMod = (int)mpz_fdiv_ui(ZT(zN), 8);  // get last 3 bits (fdiv_ui returns remainder of division)
	for (j = 0; j<sizeof(arrmult) / sizeof(arrmult[0]); j++) {
		int mod = (NbrMod * arrmult[j]) & 7;  // mod = N*arrmult[j] (mod 8)
		adjustment[j] = 0.34657359;           /*  (ln 2)/2  i.e. ln(sqrt(2)) */
		if (mod == 1)
			adjustment[j] *= (4.0e0);        // ln(4)
		if (mod == 5)
			adjustment[j] *= (2.0e0);         // ln(2)
		adjustment[j] -= log((double)arrmult[j]) / (2.0e0);
	}

	/* set up adjustment array*/
	currentPrime = 3;
	while (currentPrime < 10000) {
		int halfCurrentPrime;

		/* NbrMod = Modulus % currentPrime */
		NbrMod = (int)RemDivZnumByInt(zModulus, currentPrime);
		halfCurrentPrime = (currentPrime - 1) / 2;
		/* jacobi = NbrMod ^ HalfCurrentPrime % currentPrime */
		int jacobi = (int)modPower(NbrMod, halfCurrentPrime, currentPrime);
		double dp = (double)currentPrime;
		double logp = log(dp) / dp;

		for (j = 0; j<sizeof(arrmult) / sizeof(arrmult[0]); j++) {
			if (arrmult[j] == currentPrime) {
				adjustment[j] += logp;
			}
			else if (jacobi * modPower(arrmult[j], halfCurrentPrime,
				currentPrime) % currentPrime == 1) {
				adjustment[j] += 2 * logp;
			}
		}

		do {
			currentPrime += 2;
			for (Q = 3; Q * Q <= currentPrime; Q += 2) { /* Check if currentPrime is prime */
				if (currentPrime % Q == 0) {
					break;  /* Composite */
				}
			}
		} while (Q * Q <= currentPrime);

	}  /* end while */

	/* set value of multiplier */
	for (j = 0; j<sizeof(arrmult) / sizeof(arrmult[0]); j++) {
		if (adjustment[j] > bestadjust) { /* find biggest adjustment */
			bestadjust = adjustment[j];
			multiplier = arrmult[j]; /* multiplier =  value corresponding to biggest adjustment*/
		}
	} /* end for */

	  /* Modulus = TestNbr2 * multiplier (= N * Multiplier) */
	zModulus = zTestNbr2 * multiplier; 
	FactorBase = currentPrime;
	matrixBLength = nbrFactorBasePrimes + 50;
	rowPrimeSieveData->modsqrt = (isEven(zN)) ? 0 : 1;

	switch (mpz_get_si(ZT(zModulus)) & 0x07) {  // switch on last 3 bits
	case 1:
		logar2 = 3;
		break;
	case 5:
		logar2 = 1;
		break;
	default:
		logar2 = 1;
		break;
	}

	if (multiplier != 1 && multiplier != 2) {
		rowPrimeSieveData = &primeSieveData[2];
		rowPrimeTrialDivisionData = &ag->primeTrialDivisionData[2];
		rowPrimeSieveData->value = multiplier;
		rowPrimeTrialDivisionData->value = multiplier;
		rowPrimeSieveData->modsqrt = 0;
		// The following works because multiplier has less than 16 significant bits.
		E = (int)((1U << BITS_PER_INT_GROUP) % multiplier);
		rowPrimeTrialDivisionData->exp[0] = E;  // (2^BITS_PER_INT_GROUP) mod multiplier
		D = E * E % multiplier;
		rowPrimeTrialDivisionData->exp[1] = D;  // (2^BITS_PER_INT_GROUP)^2 mod multiplier
		D = D * E % multiplier;
		rowPrimeTrialDivisionData->exp[2] = D;  // (2^BITS_PER_INT_GROUP)^3 mod multiplier
		D = D * E % multiplier;
		rowPrimeTrialDivisionData->exp[3] = D;  // (2^BITS_PER_INT_GROUP)^4 mod multiplier
		D = D * E % multiplier;
		rowPrimeTrialDivisionData->exp[4] = D;  // (2^BITS_PER_INT_GROUP)^5 mod multiplier
		D = D * E % multiplier;
		rowPrimeTrialDivisionData->exp[5] = D;  // (2^BITS_PER_INT_GROUP)^6 mod multiplier
		j = 3;
	}
	else {
		j = 2;
	}

	currentPrime = 3;
	while (j < nbrFactorBasePrimes) { /* select small primes */
		NbrMod = (int)RemDivZnumByInt(zModulus, currentPrime);

		if (currentPrime != multiplier &&
			modPowerLL(NbrMod, (currentPrime - 1) / 2, currentPrime) == 1)
		{
			double dBase, dPower, dCurrentPrime, dRem;
			/* use only if Jacobi symbol = 0 or 1 */
			rowPrimeSieveData = &primeSieveData[j];
			rowPrimeTrialDivisionData = &ag->primeTrialDivisionData[j];
			rowPrimeSieveData->value = (int)currentPrime;
			rowPrimeTrialDivisionData->value = (int)currentPrime;
			// The following works because multiplier has less than 26 significant bits.
			dBase = (double)((1U << BITS_PER_INT_GROUP) % currentPrime);
			rowPrimeTrialDivisionData->exp[0] = (int)dBase;  // (2^BITS_PER_INT_GROUP) mod currentPrime
			dCurrentPrime = (double)currentPrime;
			dPower = dBase * dBase;
			dPower -= floor(dPower / dCurrentPrime)*dCurrentPrime;
			rowPrimeTrialDivisionData->exp[1] = (int)dPower; // (2^BITS_PER_INT_GROUP)^2 mod currentPrime
			dPower *= dBase;
			dPower -= floor(dPower / dCurrentPrime)*dCurrentPrime;
			rowPrimeTrialDivisionData->exp[2] = (int)dPower; // (2^BITS_PER_INT_GROUP)^3 mod currentPrime
			dPower *= dBase;
			dPower -= floor(dPower / dCurrentPrime)*dCurrentPrime;
			rowPrimeTrialDivisionData->exp[3] = (int)dPower; // (2^BITS_PER_INT_GROUP)^4 mod currentPrime
			dPower *= dBase;
			dPower -= floor(dPower / dCurrentPrime)*dCurrentPrime;
			rowPrimeTrialDivisionData->exp[4] = (int)dPower; // (2^BITS_PER_INT_GROUP)^5 mod currentPrime
			dPower *= dBase;
			dPower -= floor(dPower / dCurrentPrime)*dCurrentPrime;
			rowPrimeTrialDivisionData->exp[5] = (int)dPower; // (2^BITS_PER_INT_GROUP)^6 mod currentPrime

			if ((currentPrime & 3) == 3) {
				SqrRootMod = (int)modPowerLL(NbrMod, (currentPrime + 1) / 4, currentPrime);
			}
			else if ((currentPrime & 7) == 5) {   // currentPrime = 5 (mod 8)
				SqrRootMod =
					(int)modPowerLL(NbrMod * 2, (currentPrime - 5) / 8, currentPrime);
				dRem = (double)2 * NbrMod * (double)SqrRootMod;
				dRem -= floor(dRem / dCurrentPrime)*dCurrentPrime;
				dRem = dRem * (double)SqrRootMod - 1;
				dRem -= floor(dRem / dCurrentPrime)*dCurrentPrime;
				dRem = dRem * (double)NbrMod;
				dRem -= floor(dRem / dCurrentPrime)*dCurrentPrime;
				dRem = dRem * (double)SqrRootMod;
				dRem -= floor(dRem / dCurrentPrime)*dCurrentPrime;
				SqrRootMod = (int)dRem;
			}
			else {
				Q = currentPrime - 1;
				E = 0;
				Power2 = 1;

				do {
					E++;
					Q /= 2;
					Power2 *= 2;
				} while ((Q & 1) == 0); /* E >= 3 */

				Power2 /= 2;
				X = 1;

				do {
					X++;
					Z = (int)modPowerLL(X, Q, currentPrime);
				} while (modPowerLL(Z, Power2, currentPrime) == 1);

				Y = Z;
				X = (int)modPowerLL(NbrMod, (Q - 1) / 2, currentPrime);
				dBase = (double)NbrMod * (double)X;
				V = (int)(dBase - floor(dBase / dCurrentPrime)*dCurrentPrime);
				dBase = (double)V * (double)X;
				W = (int)(dBase - floor(dBase / dCurrentPrime)*dCurrentPrime);

				while (W != 1) {
					T1 = 0;
					D = W;
					while (D != 1) {
						dBase = (double)D * (double)D;
						D = (int)(dBase - floor(dBase / dCurrentPrime)*dCurrentPrime);
						T1++;
					}
					D = (int)modPowerLL(Y, 1LL << (E - T1 - 1), currentPrime);
					dBase = (double)D * (double)D;
					Y1 = (int)(dBase - floor(dBase / dCurrentPrime)*dCurrentPrime);
					E = T1;
					dBase = (double)V * (double)D;
					V1 = (int)(dBase - floor(dBase / dCurrentPrime)*dCurrentPrime);
					dBase = (double)W * (double)Y1;
					W1 = (int)(dBase - floor(dBase / dCurrentPrime)*dCurrentPrime);
					Y = Y1;
					V = V1;
					W = W1;
				} /* end while */
				SqrRootMod = V;
			} /* end if */

			rowPrimeSieveData->modsqrt = (int)SqrRootMod;
			j++;
		} /* end while */

		do {
			currentPrime += 2;
			for (Q = 3; Q * Q <= currentPrime; Q += 2) { /* Check if currentPrime is prime */
				if (currentPrime % Q == 0) {
					break;  /* Composite */
				}
			}
		} while (Q * Q <= currentPrime);

	} /* End while */

	FactorBase = currentPrime;
	largePrimeUpperBound = 100 * FactorBase;

	dlogNumberToFactor = logZnum(zN); 	// find logarithm of number to factor.
	dNumberToFactor = exp(dlogNumberToFactor);   // convert NbrToFactor to floating point

#ifdef __EMSCRIPTEN__
	getMultAndFactorBase(multiplier, FactorBase);  // append Mult & Factor base to SIQS string
	printf_s("%s", lowerText);  // print Mult and Factor Base
#endif

	firstLimit = 2;
	for (j = 2; j < nbrFactorBasePrimes; j++) {
		firstLimit *= (int)(primeSieveData[j].value);
		if (firstLimit > 2 * SieveLimit) {
			break;
		}
	}

	dNumberToFactor *= multiplier;
	smallPrimeUpperLimit = j + 1;
	threshold = (int) (log(sqrt(dNumberToFactor) * SieveLimit /
				(FactorBase * 64) /
					primeSieveData[j + 1].value
				) / log(3) + 0x81
			);
	firstLimit = (int)(log(dNumberToFactor) / 3);
	for (secondLimit = firstLimit; secondLimit < nbrFactorBasePrimes; secondLimit++) {
		if (primeSieveData[secondLimit].value * 2 > SieveLimit) {
			break;
		}
	}

	for (thirdLimit = secondLimit; thirdLimit < nbrFactorBasePrimes; thirdLimit++) {
		if (primeSieveData[thirdLimit].value > 2 * SieveLimit) {
			break;
		}
	}

	nbrPrimes2 = nbrFactorBasePrimes - 4;
	Prod = sqrt(2 * dNumberToFactor) / (double)SieveLimit;

	fact = (int)pow(Prod, 1 / (float)nbrFactorsA);
	for (i = 2;; i++) {
		if (primeSieveData[i].value > fact) {
			break;
		}
	}

	span = nbrFactorBasePrimes / (2 * nbrFactorsA*nbrFactorsA);
	if (nbrFactorBasePrimes < 500) {
		span *= 2;
	}
	indexMinFactorA = i - span / 2;
	/*********************************************/
	/* Generate sieve threads                    */
	/*********************************************/
	sieveThread(Factor);			   /* return factor found in Factor */

#if 0
	for (threadNumber = 0; threadNumber<numberThreads; threadNumber++)
	{
		//new Thread(this).start();                // Start new thread.
		//synchronized(amodq)
		{
			//      while (threadArray[threadNumber] == null &&
			//        getTerminateThread() == false)
			{
				//try
				//{
				//  amodq.wait();
				//}
				//catch (InterruptedException ie) {}
			}
		}
	}
	//synchronized(matrixB)
	{
		while (factorSiqs == null && getTerminateThread() == false)
		{
			try
			{
				matrixB.wait();
			}
			catch (InterruptedException ie) {}
		}
	}
#endif

	matrixB.clear();   /* release memory */
	matrixPartial.clear();
	primeSieveData.clear();
	free(ag);


	if (/*getTerminateThread() ||*/ (Factor == 0))
	{
		//throw new ArithmeticException();
	}
#if 0
	for (threadNumber = 0; threadNumber<numberThreads; threadNumber++)
	{                 // Wake up all sieve threads so they can terminate.
		if (threadArray[threadNumber].isAlive())
		{
			//try
			{
				threadArray[threadNumber].interrupt();
			}
			//catch (Exception e) {}
		}
	}
#endif
#if 0
	synchronized(this)
	{
		saveSIQSStatistics(polynomialsSieved, trialDivisions,
			smoothsFound, totalPartials,
			partialsFound, ValuesSieved);
	}
	return factorSiqs;
#endif
}

static void ShowSIQSStatus(int *rowMatrixB) {
#ifdef __EMSCRIPTEN__
	int elapsedTime = (int)(tenths() - originalTenthSecond);
	if (elapsedTime  >= oldTimeElapsed +50)
	{   // send no more than 1 message per 5 seconds
		oldTimeElapsed = elapsedTime;
		ShowSIQSInfo((elapsedTime - startSieveTenths) / 10, congruencesFound, matrixBLength,
			elapsedTime / 10, rowMatrixB);
#if DEBUG_SIQS
		if (debug_ctr < debugLimit)
			std::cout << "Partials = " << partialsFound << " Smooth = " << smoothsFound << '\n';
#endif
	}
#endif
}

static int EraseSingletons(int nbrFactorBasePrimes) {
	int row, column, delta;
	int *rowMatrixB;
	int matrixBlength = matrixBLength;

	memset(ag->newColumns, 0, matrixBlength * sizeof(int));
	// Find singletons in matrixB storing in array vectExpParity the number
	// of primes in each column.
	do {   // The singleton removal phase must run until there are no more
		   // singletons to erase.
		memset(ag->vectExpParity, 0, matrixBLength * sizeof(limb));
		for (row = matrixBlength - 1; row >= 0; row--)
		{                  // Traverse all rows of the matrix.
			rowMatrixB = &matrixB[row][0];
			for (column = rowMatrixB[LENGTH_OFFSET] - 1; column >= 1; column--)
			{                // A prime appeared in that column.
				ag->vectExpParity[rowMatrixB[column]]++;
			}
		}
		row = 0;
		for (column = 0; column<nbrFactorBasePrimes; column++) {
			if (ag->vectExpParity[column] > 1)
			{                // Useful column found with at least 2 primes.
				ag->newColumns[column] = row;
				ag->primeTrialDivisionData[row++].value =
					ag->primeTrialDivisionData[column].value;
			}
		}
		nbrFactorBasePrimes = row;
		delta = 0;
		// Erase singletons from matrixB. The rows to be erased are those where the
		// the corresponding element of the array vectExpParity equals 1.
		for (row = 0; row < matrixBlength; row++)
		{                  // Traverse all rows of the matrix.
			rowMatrixB = &matrixB[row][0];
			for (column = rowMatrixB[LENGTH_OFFSET] - 1; column >= 1; column--)
			{                // Traverse all columns.
				if (ag->vectExpParity[rowMatrixB[column]] == 1)
				{              // Singleton found: erase this row.
					delta++;
					break;
				}
			}
			if (column == 0)
			{                // Singleton not found: move row upwards.
				memcpy(&matrixB[row - delta][0], &matrixB[row][0], sizeof(int) * MAX_FACTORS_RELATION);
				vectLeftHandSide[row - delta] = vectLeftHandSide[row];
				//memcpy(vectLeftHandSide[row - delta], vectLeftHandSide[row], sizeof(vectLeftHandSide[0]));
			}
		}
		matrixBlength -= delta;      // Update number of rows of the matrix.
		for (row = 0; row < matrixBlength; row++)
		{                  // Traverse all rows of the matrix.
			rowMatrixB = &matrixB[row][0];
			for (column = rowMatrixB[LENGTH_OFFSET]; column >= 1; column--)
			{                // Change all column indexes in this row.
				rowMatrixB[column] = ag->newColumns[rowMatrixB[column]];
			}
		}
	} while (delta > 0);           // End loop if number of rows did not change.

	ag->primeTrialDivisionData[0].exp[1] = nbrFactorBasePrimes;
	return matrixBlength;
}

/****************************************************************/
/* Linear algebra phase. returns true &factor found in biT if   */
/* factor is found, otherwise returns false                     */
/* Uses global variables zModulus, nbrFactorBasePrimes,         */
/* primeTrialDivisionData, matrixRows, matrixCols, vectExpParity */
/****************************************************************/
static bool LinearAlgebraPhase(
	Znum &biT, Znum &biR, Znum &biU)
{
	int mask, row, col, j;
	int *rowMatrixB;
	int primeIndex;

#if DEBUG_SIQS
	{
		int i;
		printf("******* Linear Algebra Phase ************\n");
		for (j = 0; j < matrixBLength; j++)
		{
			char *ptrOutput = output;

			std::cout << vectLeftHandSide[j] << ", ";
			for (i = 1; i < matrixB[j][0]; i++)
			{
				int2dec(&ptrOutput, matrixB[j][i]);
				*ptrOutput++ = ',';
			}
			*ptrOutput = 0;
			printf("%s\n", output);
		}
		//    exit(0);
	}
#endif
	// Get new number of rows after erasing singletons.
	int matrixBlen = EraseSingletons(nbrFactorBasePrimes);
	matrixBLength = matrixBlen;  // note upper case L; not same variable
	matrixRows = matrixBlen;
	matrixCols = ag->primeTrialDivisionData[0].exp[1];
	ag->primeTrialDivisionData[0].exp[1] = 0;         // Restore correct value.
	BlockLanczos();
	// The rows of matrixV indicate which rows must be multiplied so no
	// primes are multiplied an odd number of times.
	mask = 1;
	for (col = BITS_PER_INT_GROUP; col >= 0; col--) {

		biT = 1; 
		biR = 1; 
		memset(ag->vectExpParity, 0, matrixBlen * sizeof(ag->vectExpParity[0]));

		for (row = matrixBlen - 1; row >= 0; row--) {
			if ((ag->matrixV[row] & mask) != 0) {
				MultZnumModN(vectLeftHandSide[row], biR, biU, zModulus);  // U = LHS * R (mod Modulus)
				biR = biU; 
				rowMatrixB = &matrixB[row][0];
				for (j = rowMatrixB[LENGTH_OFFSET] - 1; j >= 1; j--) {
					primeIndex = rowMatrixB[j];
					ag->vectExpParity[primeIndex] ^= 1;
					if (ag->vectExpParity[primeIndex] == 0) {
						if (primeIndex == 0) {
							biT = zModulus - biT; //SubtractBigNbr(zModulus, biT, biT); // Multiply biT by -1.
						}
						else {
							MultZnumByIntModN(biT, ag->primeTrialDivisionData[primeIndex].value,
								biT, zModulus);
						}
					}
				}
			}
		}

		SubtractZnumModN(biR, biT, biR, zModulus);
		biT = gcd(biR, zTestNbr2);
#if DEBUG_SIQS
		std::cout << "col = " << col
			<< " biT = " << biT
			<< " biR = " << biR
			<< " zTestNbr2 = " << zTestNbr2 << '\n';
#endif

		if (biT > 1) {   // GCD is not zero or 1.
			if (biT != zTestNbr2) {
#if DEBUG_SIQS
				std::cout << "LinearAlgebra returns true. \n";
#endif
				return true;
			}
		}
		mask *= 2;
	}
#if DEBUG_SIQS
	std::cout << "LinearAlgebra returns false. biT = " << biT
		<< " zTestNbr2 = " << zTestNbr2 
		<< " biR = " << biR<< '\n';
#endif
	return false;
}


/* return true if either relation added to matrixB (also increment global 
var congruencesFound), or matrixB is full.
return false if relation already in matrixB
uses global variables biModulus, biTestNbr2, */
static bool InsertNewRelation(
	int *rowMatrixB,
	Znum &biT, Znum &biU, Znum &biR)
{
	int i, k;
	//int lenDivisor;
	int nbrColumns = rowMatrixB[LENGTH_OFFSET];
	// Insert it only if it is different from previous relations.
	if (congruencesFound >= matrixBLength) {  // Discard excess congruences.
		return true;
	}
#if DEBUG_SIQS
	{
		char *ptrOutput = output;
		static int nn;
		std::cout << "biT = " << biT;
		std::cout << "biU = " << biU;
		std::cout << "biR = " << biR ;
		*ptrOutput++ = ',';
		for (i = 1; i < *rowMatrixB; i++)
		{
			int2dec(&ptrOutput, *(rowMatrixB + i));
			*ptrOutput++ = ',';
		}
		*ptrOutput = 0;
		printf(" rowMatrixB = %s\n", output);
	}
	if (++nn == 3018)
	{
		nn = 3018;
	}
#endif
	// Check whether this relation is already in the matrix.
	size_t row = 0;
	for (i = 0; i < congruencesFound; i++) {
		int *curRowMatrixB = &matrixB[row][0];
		if (nbrColumns == *(curRowMatrixB + LENGTH_OFFSET)) {
			for (k = 1; k < nbrColumns; k++) {
				if (*(rowMatrixB + k) != curRowMatrixB[k]) {
					break;
				}
			}
			if (k == nbrColumns) {
				return false; // Do not insert same relation.
			}
		}
		//curRowMatrixB += MAX_FACTORS_RELATION;
		row++;
	}

	/* Convert negative numbers to the range 0 <= n < biModulus */
	if (isEven(zModulus)) { // Even modulus.
		zTestNbr2 = zModulus / 2;
		//DivBigNbrByInt(zModulus, 2, zTestNbr2);  // divide modulus by 2

		// If biR >= biModulus perform biR = biR - biModulus.
		biT = 0;
		AddZnumModN(biR, biT, biR, zTestNbr2);

		ModInvZnum(biR, biT, zTestNbr2);   // biT = Mod Inv of biR
	}

	else {             // Odd modulus
		ModInvZnum(biR, biT, zModulus);  // biT = Mod Inv 
	}

	if (biU < 0) {
		biU += zModulus; //AddBigNbr(biU, zModulus, biU);  
	}

	mpz_mod(ZT(biU), ZT(biU), ZT(zModulus));   /* biU %= zModulus */
	// Compute biU / biR  (mod biModulus)
	MultZnumModN(biU, biT, biR, zModulus);  // biR = biU * biT

	// Add relation to matrix B.
	memcpy(&matrixB[congruencesFound][0], &rowMatrixB[0], nbrColumns * sizeof(int));

	vectLeftHandSide[congruencesFound] = biR;
	congruencesFound++;
	
	nbrColumns = rowMatrixB[LENGTH_OFFSET];
	for (int k = 1; k < nbrColumns; k++) {
		if (ag->primesUsed[rowMatrixB[k]] == 0) {
			ag->primesUsed[rowMatrixB[k]] = 1;   // show that prime is used
			nbrPrimesUsed++;                  // update count of primes used
		}
	}

#if DEBUG_SIQS
		std::cout << congruencesFound << " Congruences found: biR = " << biR << '\n';
#endif
	return true;
}

static int intModInv(int NbrMod, int currentPrime) {
	int QQ, T1, T3;
	int V1 = 1;
	int V3 = NbrMod;
	int U1 = 0;
	int U3 = currentPrime;

	while (V3 != 0) {
		if (U3 < V3 + V3) {           // QQ = 1
			T1 = U1 - V1;
			T3 = U3 - V3;
		}
		else {
			QQ = U3 / V3;
			T1 = U1 - V1 * QQ;
			T3 = U3 - V3 * QQ;
		}
		U1 = V1;
		U3 = V3;
		V1 = T1;
		V3 = T3;
	}
	return U1 + (currentPrime & (U1 >> BITS_PER_INT_GROUP));
}

/* Multiply binary matrices of length m x 32 by 32 x 32 */
/* The product matrix has size m x 32. Then add it to a m x 32 matrix. */
static void MatrixMultAdd(int *LeftMatr, int *RightMatr, int *ProdMatr) {
	int matrLength = matrixBLength;
	int row;

	for (row = 0; row < matrLength; row++) {
		int col;
		int prodMatr = *(ProdMatr + row);
		int leftMatr = *(LeftMatr + row);
		for (col = 0; col<32; col++) {
			if (leftMatr < 0) {
				prodMatr ^= *(RightMatr + col);
			}
			leftMatr *= 2;
		}
		*(ProdMatr + row) = prodMatr;
	}
}

/* Multiply binary matrices of length m x 32 by 32 x 32 */
/* The product matrix has size m x 32 */
static void MatrixMultiplication(int *LeftMatr, int *RightMatr, int *ProdMatr) {
	int matrLength = 32;
	int row;

	for (row = 0; row < matrLength; row++) {
		int col;
		int prodMatr = 0;
		int leftMatr = *(LeftMatr + row);
		for (col = 0; col < 32; col++) {
			if (leftMatr < 0) {
				prodMatr ^= *(RightMatr + col);
			}
			leftMatr *= 2;
		}
		*(ProdMatr + row) = prodMatr;
	}
}

/* Multiply the transpose of a binary matrix of length n x 32 by */
/* another binary matrix of length n x 32 */
/* The product matrix has size 32 x 32 */
static void MatrTranspMult(int matrLength, int *LeftMatr, int *RightMatr, int *ProdMatr) {
	int row, col;
	int iMask = 1;

	for (col = 31; col >= 0; col--) {
		int prodMatr = 0;
		for (row = 0; row < matrLength; row++) {
			if ((*(LeftMatr + row) & iMask) != 0) {
				prodMatr ^= *(RightMatr + row);
			}
		}
		*(ProdMatr + col) = prodMatr;
		iMask *= 2;
	}
}

static void MatrixAddition(int *leftMatr, int *rightMatr, int *sumMatr) {
	int row;

	for (row = 32 - 1; row >= 0; row--) {
		*(sumMatr + row) = *(leftMatr + row) ^ *(rightMatr + row);
	}
}

static void MatrMultBySSt(int length, int *Matr, int diagS, int *Prod) {
	int row;

	for (row = length - 1; row >= 0; row--) {
		*(Prod + row) = diagS & *(Matr + row);
	}
}

/* Compute Bt * B * input matrix where B is the matrix that holds the */
/* factorization relations */
static void MultiplyAByMatrix(int *Matr, int *TempMatr, int *ProdMatr) {
	int index;
	int row;
	int *rowMatrixB;

	/* Compute TempMatr = B * Matr */
	memset(TempMatr, 0, matrixBLength * sizeof(int));
	for (row = matrixBLength - 1; row >= 0; row--) {
		int rowValue;
		rowMatrixB = &matrixB[row][0];
		rowValue = *(Matr + row);
		for (index = *(rowMatrixB + LENGTH_OFFSET) - 1; index >= 1; index--) {
			*(TempMatr + *(rowMatrixB + index)) ^= rowValue;
		}
	}

	/* Compute ProdMatr = Bt * TempMatr */
	for (row = matrixBLength - 1; row >= 0; row--) {
		int prodMatr = 0;
		rowMatrixB = &matrixB[row][0];
		for (index = *(rowMatrixB + LENGTH_OFFSET) - 1; index >= 1; index--) {
			prodMatr ^= *(TempMatr + *(rowMatrixB + index));
		}
		*(ProdMatr + row) = prodMatr;
	}
}

static void colexchange(int *XmY, int *V, int *V1, int *V2, int col1, int col2) {
	int row;
	int mask1, mask2;
	int *matr1, *matr2;

	if (col1 == col2) {          // Cannot exchange the same column.
		return;
	}          // Exchange columns col1 and col2 of V1:V2
	mask1 = 0x80000000 >> (col1 & 31);
	mask2 = 0x80000000 >> (col2 & 31);
	matr1 = (col1 >= 32 ? V1 : V2);
	matr2 = (col2 >= 32 ? V1 : V2);

	for (row = matrixBLength - 1; row >= 0; row--) {
		// If both bits are different toggle them.
		if (((matr1[row] & mask1) == 0) != ((matr2[row] & mask2) == 0)) {
			// If both bits are different toggle them.
			matr1[row] ^= mask1;
			matr2[row] ^= mask2;
		}
	}
	// Exchange columns col1 and col2 of XmY:V
	matr1 = (col1 >= 32 ? XmY : V);
	matr2 = (col2 >= 32 ? XmY : V);
	for (row = matrixBLength - 1; row >= 0; row--) {
		// If both bits are different toggle them.
		if (((matr1[row] & mask1) == 0) != ((matr2[row] & mask2) == 0)) {
			matr1[row] ^= mask1;
			matr2[row] ^= mask2;
		}
	}
}

static void coladd(int *XmY, int *V, int *V1, int *V2, int col1, int col2) {
	int row;
	int mask1, mask2;
	int *matr1, *matr2;

	if (col1 == col2) {
		return;
	}
	// Add column col1 to column col2 of V1:V2
	mask1 = 0x80000000 >> (col1 & 31);
	mask2 = 0x80000000 >> (col2 & 31);
	matr1 = (col1 >= 32 ? V1 : V2);
	matr2 = (col2 >= 32 ? V1 : V2);

	for (row = matrixBLength - 1; row >= 0; row--) {
		// If bit to add is '1'...
		if ((matr1[row] & mask1) != 0) {    // Toggle bit in destination.
			matr2[row] ^= mask2;
		}
	}

	// Add column col1 to column col2 of XmY:V
	matr1 = (col1 >= 32 ? XmY : V);
	matr2 = (col2 >= 32 ? XmY : V);
	for (row = matrixBLength - 1; row >= 0; row--) {
		// If bit to add is '1'...
		if ((matr1[row] & mask1) != 0) {
			matr2[row] ^= mask2; // Toggle bit in destination.
		}
	}
}

/* uses global variables matrixV, matrixXmY
*/
static void BlockLanczos(void)
{
	int i, j, k;
	int oldDiagonalSSt, newDiagonalSSt;
	int index, mask;
	int matrixD[32] = { 0 };
	int matrixE[32] = { 0 };
	int matrixF[32] = { 0 };
	int matrixWinv[32] = { 0 };
	int matrixWinv1[32] = { 0 };
	int matrixWinv2[32] = { 0 };
	int matrixVtV0[32] = { 0 };
	int matrixVt1V0[32] = { 0 };
	int matrixVt2V0[32] = { 0 };
	int matrixVtAV[32] = { 0 };
	int matrixVt1AV1[32] = { 0 };
	int matrixCalcParenD[32] = { 0 };
	int vectorIndex[64] = { 0 };
	int matrixTemp[32] = { 0 };
	int matrixCalc1[32] = { 0 }; // Matrix that holds temporary data
	int matrixCalc2[32] = { 0 }; // Matrix that holds temporary data
	int *matr;
	double dSeed, dMult, dDivisor, dAdd;
	int Temp, Temp1;
	int stepNbr = 0;
	int currentOrder, currentMask;
	int row, col;
	int leftCol, rightCol;
	int minind, min, minanswer;
	int *rowMatrixB, *ptrMatrixV, *ptrMatrixXmY;

	newDiagonalSSt = oldDiagonalSSt = -1;
	memset(matrixWinv, 0, sizeof(matrixWinv));
	memset(matrixWinv1, 0, sizeof(matrixWinv1));
	memset(matrixWinv2, 0, sizeof(matrixWinv2));
	memset(matrixVtV0, 0, sizeof(matrixVtV0));
	memset(matrixVt1V0, 0, sizeof(matrixVt1V0));
	memset(matrixVt2V0, 0, sizeof(matrixVt2V0));
	memset(matrixVt1AV1, 0, sizeof(matrixVt1AV1));

	/* Initialize matrix X-Y and matrix V_0 with random data */
	dSeed = (double)123456789;
	dMult = (double)62089911;
	dAdd = (double)54325442;
	dDivisor = (double)0x7FFFFFFF;
	ptrMatrixXmY = &ag->matrixXmY[matrixBLength - 1];
	for (ptrMatrixV = &ag->matrixV[matrixBLength - 1]; ptrMatrixV >= ag->matrixV; ptrMatrixV--)
	{
		double dSeed2 = (dSeed * dMult + dAdd);
		dSeed2 -= floor(dSeed2 / dDivisor) * dDivisor;
		*ptrMatrixXmY-- = (int)dSeed + (int)dSeed2;
		dSeed = (dSeed2 * dMult + dAdd);
		dSeed -= floor(dSeed / dDivisor) * dDivisor;
		dSeed2 = (dSeed * dMult + dAdd);
		dSeed2 -= floor(dSeed2 / dDivisor) * dDivisor;
		*ptrMatrixV = (int)dSeed + (int)dSeed2;
		dSeed = (dSeed2 * dMult + dAdd);
		dSeed -= floor(dSeed / dDivisor) * dDivisor;
	}
	// Compute matrix Vt(0) * V(0)
	MatrTranspMult(matrixBLength, ag->matrixV, ag->matrixV, matrixVtV0);
#ifdef __EMSCRIPTEN__
	showMatrixSize((char *)SIQSInfoText, matrixRows, matrixCols);
	first = true; 
#endif
#if DEBUG_SIQS
	{
		char *ptrOutput = output;
		strcpy(ptrOutput, "MatrixBLength = ");
		ptrOutput += strlen(ptrOutput);
		int2dec(&ptrOutput, matrixBLength);
		*ptrOutput = 0;
		printf("%s\n", output);
	}
#endif
	for (;;) {
		int indexC;
#ifdef __EMSCRIPTEN__
		int elapsedTime = (int)(tenths() - originalTenthSecond);
		if (elapsedTime / 10 != oldTimeElapsed / 10) {
			/* send message once per second */
			char SIQSInfo[200];
			char *ptrText = SIQSInfo;
			oldTimeElapsed = elapsedTime;
			//strcpy(ptrText, "4\n");
			//ptrText += strlen(ptrText);
			GetDHMS(&ptrText, elapsedTime / 10);
			strcpy(ptrText, "  ");
			ptrText += strlen(ptrText);
			strcpy(ptrText, lang ? "Progreso del álgebra lineal: " : "Linear algebra progress: ");
			ptrText += strlen(ptrText);
			int2dec(&ptrText, stepNbr * 3200 / matrixRows);
			strcpy(ptrText, "%\n");

			if (first) {
				if (!GetConsoleScreenBufferInfo(hConsole, &csbi))
				{
					fprintf(stderr, "** GetConsoleScreenBufferInfo failed with %d!\n", GetLastError());
					Beep(750, 1000);
				}
				coordScreen.X = csbi.dwCursorPosition.X;  // save cursor co-ordinates
				coordScreen.Y = csbi.dwCursorPosition.Y;
				if (csbi.dwCursorPosition.Y >= csbi.dwSize.Y - 1)
					coordScreen.Y--;   // if window is full, allow for text scrolling up

			}
			else
				upOneLine();

			printf_s("%s", SIQSInfo);
			first = false;
		}
#endif
		//if (getTerminateThread())
		//{
		//  throw new ArithmeticException();
		// }
		oldDiagonalSSt = newDiagonalSSt;
		stepNbr++;
		// Compute matrix A * V(i)
		MultiplyAByMatrix(ag->matrixV, ag->matrixCalc3, ag->matrixAV);
		// Compute matrix Vt(i) * A * V(i)
		MatrTranspMult(matrixBLength, ag->matrixV, ag->matrixAV, matrixVtAV);

		/* If Vt(i) * A * V(i) = 0, end of loop */
		for (i = sizeof(matrixVtAV) / sizeof(matrixVtAV[0]) - 1; i >= 0; i--) {
			if (matrixVtAV[i] != 0) {
				break;
			}
		}
		if (i < 0) {
			break;
		} /* End X-Y calculation loop */

		  /* Selection of S(i) and W(i) */

		memcpy(matrixTemp, matrixWinv2, sizeof(matrixTemp));
		memcpy(matrixWinv2, matrixWinv1, sizeof(matrixTemp));
		memcpy(matrixWinv1, matrixWinv, sizeof(matrixTemp));
		memcpy(matrixWinv, matrixTemp, sizeof(matrixTemp));

		mask = 1;
		for (j = 31; j >= 0; j--) {
			matrixD[j] = matrixVtAV[j]; /*  D = VtAV    */
			matrixWinv[j] = mask; /*  Winv = I    */
			mask *= 2;
		}

		index = 31;
		mask = 1;

		for (indexC = 31; indexC >= 0; indexC--) {
			if ((oldDiagonalSSt & mask) != 0) {
				matrixE[index] = indexC;
				matrixF[index] = mask;
				index--;
			}
			mask *= 2;
		}

		mask = 1;
		for (indexC = 31; indexC >= 0; indexC--) {
			if ((oldDiagonalSSt & mask) == 0) {
				matrixE[index] = indexC;
				matrixF[index] = mask;
				index--;
			}
			mask *= 2;
		}

		newDiagonalSSt = 0;
		for (j = 0; j < 32; j++) {
			currentOrder = matrixE[j];
			currentMask = matrixF[j];
			for (k = j; k < 32; k++) {
				if ((matrixD[matrixE[k]] & currentMask) != 0) {
					break;
				}
			}

			if (k < 32) {
				i = matrixE[k];
				Temp = matrixWinv[i];
				matrixWinv[i] = matrixWinv[currentOrder];
				matrixWinv[currentOrder] = Temp;
				Temp1 = matrixD[i];
				matrixD[i] = matrixD[currentOrder];
				matrixD[currentOrder] = Temp1;
				newDiagonalSSt |= currentMask;
				for (k = 31; k >= 0; k--) {
					if (k != currentOrder && ((matrixD[k] & currentMask) != 0)) {
						matrixWinv[k] ^= Temp;
						matrixD[k] ^= Temp1;
					}
				} /* end for k */
			}
			else {
				for (k = j; k < 31; k++) {
					if ((matrixWinv[matrixE[k]] & currentMask) != 0) {
						break;
					}
				}
				i = matrixE[k];
				Temp = matrixWinv[i];
				matrixWinv[i] = matrixWinv[currentOrder];
				matrixWinv[currentOrder] = Temp;
				Temp1 = matrixD[i];
				matrixD[i] = matrixD[currentOrder];
				matrixD[currentOrder] = Temp1;
				for (k = 31; k >= 0; k--) {
					if ((matrixWinv[k] & currentMask) != 0) {
						matrixWinv[k] ^= Temp;
						matrixD[k] ^= Temp1;
					}
				} /* end for k */
			} /* end if */
		} /* end for j */

		  /* Compute D(i), E(i) and F(i) */
#if DEBUG_SIQS
		char *ptrOutput = output;
		if (stepNbr < 200)
		{
			strcpy(ptrOutput, "Step #");
			ptrOutput += strlen(ptrOutput);
			int2dec(&ptrOutput, stepNbr);
			strcpy(ptrOutput, ": matrixWinv1[0] = ");
			ptrOutput += strlen(ptrOutput);
			int2dec(&ptrOutput, matrixWinv1[0]);
			*ptrOutput = 0;
			printf("%s\n", output);
		}
#endif
		if (stepNbr >= 3) {
			// F = -Winv(i-2) * (I - Vt(i-1)*A*V(i-1)*Winv(i-1)) * ParenD * S*St
			MatrixMultiplication(matrixVt1AV1, matrixWinv1, matrixCalc2);
			mask = 1; /* Add identity matrix */
			for (index = 31; index >= 0; index--) {
				matrixCalc2[index] ^= mask;
				mask *= 2;
			}
			MatrixMultiplication(matrixWinv2, matrixCalc2, matrixCalc1);
			MatrixMultiplication(matrixCalc1, matrixCalcParenD, matrixF);
			MatrMultBySSt(32, matrixF, newDiagonalSSt, matrixF);
		}

		// E = -Winv(i-1) * Vt(i)*A*V(i) * S*St
		if (stepNbr >= 2) {
			MatrixMultiplication(matrixWinv1, matrixVtAV, matrixE);
			MatrMultBySSt(32, matrixE, newDiagonalSSt, matrixE);
		}

		// ParenD = Vt(i)*A*A*V(i) * S*St + Vt(i)*A*V(i)
		// D = I - Winv(i) * ParenD
		MatrTranspMult(matrixBLength, ag->matrixAV, ag->matrixAV, matrixCalc1); // Vt(i)*A*A*V(i)
		MatrMultBySSt(32, matrixCalc1, newDiagonalSSt, matrixCalc1);
		MatrixAddition(matrixCalc1, matrixVtAV, matrixCalcParenD);
		MatrixMultiplication(matrixWinv, matrixCalcParenD, matrixD);
		mask = 1; /* Add identity matrix */
		for (index = 31; index >= 0; index--) {
			matrixD[index] ^= mask;
			mask *= 2;
		}

		/* Update value of X - Y */
		MatrixMultiplication(matrixWinv, matrixVtV0, matrixCalc1);
		MatrixMultAdd(ag->matrixV, matrixCalc1, ag->matrixXmY);

		/* Compute value of new matrix V(i) */
		// V(i+1) = A * V(i) * S * St + V(i) * D + V(i-1) * E + V(i-2) * F
		MatrMultBySSt(matrixBLength, ag->matrixAV, newDiagonalSSt, ag->matrixCalc3);
		MatrixMultAdd(ag->matrixV, matrixD, ag->matrixCalc3);
		if (stepNbr >= 2) {
			MatrixMultAdd(ag->matrixV1, matrixE, ag->matrixCalc3);
			if (stepNbr >= 3) {
				MatrixMultAdd(ag->matrixV2, matrixF, ag->matrixCalc3);
			}
		}

		/* Compute value of new matrix Vt(i)V0 */
		// Vt(i+1)V(0) = Dt * Vt(i)V(0) + Et * Vt(i-1)V(0) + Ft * Vt(i-2)V(0)
		MatrTranspMult(32, matrixD, matrixVtV0, matrixCalc2);
		if (stepNbr >= 2) {
			MatrTranspMult(32, matrixE, matrixVt1V0, matrixCalc1);
			MatrixAddition(matrixCalc1, matrixCalc2, matrixCalc2);
			if (stepNbr >= 3) {
				MatrTranspMult(32, matrixF, matrixVt2V0, matrixCalc1);
				MatrixAddition(matrixCalc1, matrixCalc2, matrixCalc2);
			}
		}
		/* cycle arrays:
		  V2 <- V1, V1 <- V, V <- Calc3, Calc3 <- V2*/
		memcpy(ag->matrixTemp2, ag->matrixV2, sizeof(ag->matrixTemp2));
		memcpy(ag->matrixV2, ag->matrixV1, sizeof(ag->matrixV2));
		memcpy(ag->matrixV1, ag->matrixV, sizeof(ag->matrixV1));
		memcpy(ag->matrixV, ag->matrixCalc3, sizeof(ag->matrixV));
		memcpy(ag->matrixCalc3, ag->matrixTemp2, sizeof(ag->matrixCalc3));
		/* Vt2V0 <- Vt1V0, Vt1V0 <- VtV0, VtV0 <- Calc2, Calc2 <- Vt2V0*/
		memcpy(matrixTemp, matrixVt2V0, sizeof(matrixTemp));
		memcpy(matrixVt2V0, matrixVt1V0, sizeof(matrixVt2V0));
		memcpy(matrixVt1V0, matrixVtV0, sizeof(matrixVt1V0));
		memcpy(matrixVtV0, matrixCalc2, sizeof(matrixVtV0));
		memcpy(matrixCalc2, matrixTemp, sizeof(matrixCalc2));
		/* swap Vt1AV1 and VtAV */
		memcpy(matrixTemp, matrixVt1AV1, sizeof(matrixTemp));
		memcpy(matrixVt1AV1, matrixVtAV, sizeof(matrixVt1AV1));
		memcpy(matrixVtAV, matrixTemp, sizeof(matrixVtAV));
	} /* end while */

	  /* Find matrix V1:V2 = B * (X-Y:V) */
	for (row = matrixBLength - 1; row >= 0; row--) {
		ag->matrixV1[row] = ag->matrixV2[row] = 0;
	}
	for (row = matrixBLength - 1; row >= 0; row--) {
		int rowMatrixXmY, rowMatrixV;

		rowMatrixB = &matrixB[row][0];
		rowMatrixXmY = ag->matrixXmY[row];
		rowMatrixV = ag->matrixV[row];
		// The vector rowMatrixB includes the indexes of the columns set to '1'.
		for (index = rowMatrixB[LENGTH_OFFSET] - 1; index >= 1; index--) {
			col = rowMatrixB[index];
			ag->matrixV1[col] ^= rowMatrixXmY;
			ag->matrixV2[col] ^= rowMatrixV;
		}
	}

	rightCol = 64;
	leftCol = 0;
	while (leftCol < rightCol) {
		for (col = leftCol; col < rightCol; col++)
		{       // For each column find the first row which has a '1'.
				// Columns outside this range must have '0' in all rows.
			matr = (col >= 32 ? ag->matrixV1 : ag->matrixV2);
			mask = 0x80000000 >> (col & 31);
			vectorIndex[col] = -1;    // indicate all rows in zero in advance.
			for (row = 0; row < matrixBLength; row++) {
				if ((matr[row] & mask) != 0)
				{               // First row for this mask is found. Store it.
					vectorIndex[col] = row;
					break;
				}
			}
		}

		for (col = leftCol; col < rightCol; col++) {
			if (vectorIndex[col] < 0)
			{  // If all zeros in col 'col', exchange it with first column with
			   // data different from zero (leftCol).
				colexchange(ag->matrixXmY, ag->matrixV, ag->matrixV1, ag->matrixV2, leftCol, col);
				vectorIndex[col] = vectorIndex[leftCol];
				vectorIndex[leftCol] = -1;  // This column now has zeros.
				leftCol++;                  // Update leftCol to exclude that column.
			}
		}

		if (leftCol == rightCol) {
			break;
		}

		// At this moment all columns from leftCol to rightCol are non-zero.
		// Get the first row that includes a '1'.
		min = vectorIndex[leftCol];
		minind = leftCol;
		for (col = leftCol + 1; col < rightCol; col++) {
			if (vectorIndex[col] < min) {
				min = vectorIndex[col];
				minind = col;
			}
		}

		minanswer = 0;
		for (col = leftCol; col < rightCol; col++) {
			if (vectorIndex[col] == min) {
				minanswer++;
			}
		}

		if (minanswer > 1) {    // Two columns with the same first row to '1'.
			for (col = minind + 1; col < rightCol; col++) {
				if (vectorIndex[col] == min)
				{        // Add first column which has '1' in the same row to
						 // the other columns so they have '0' in this row after
						 // this operation.
					coladd(ag->matrixXmY, ag->matrixV, ag->matrixV1, ag->matrixV2, minind, col);
				}
			}
		}
		else {
			rightCol--;
			colexchange(ag->matrixXmY, ag->matrixV, ag->matrixV1, ag->matrixV2, minind, rightCol);
		}
	}

	leftCol = 0; /* find linear independent solutions */
	while (leftCol < rightCol) {
		for (col = leftCol; col < rightCol; col++)
		{         // For each column find the first row which has a '1'.
			matr = (col >= 32 ? ag->matrixXmY : ag->matrixV);
			mask = 0x80000000 >> (col & 31);
			vectorIndex[col] = -1;    // indicate all rows in zero in advance.
			for (row = 0; row < matrixBLength; row++) {
				if ((matr[row] & mask) != 0)
				{         // First row for this mask is found. Store it.
					vectorIndex[col] = row;
					break;
				}
			}
		}

		for (col = leftCol; col < rightCol; col++)
		{  // If all zeros in col 'col', exchange it with last column with
		   // data different from zero (rightCol).
			if (vectorIndex[col] < 0) {
				rightCol--;                 // Update rightCol to exclude that column.
				colexchange(ag->matrixXmY, ag->matrixV, ag->matrixV1, ag->matrixV2, rightCol, col);
				vectorIndex[col] = vectorIndex[rightCol];
				vectorIndex[rightCol] = -1; // This column now has zeros.
			}
		}

		if (leftCol == rightCol) {
			break;  // exit from while loop
		}

		// At this moment all columns from leftCol to rightCol are non-zero.
		// Get the first row that includes a '1'.
		min = vectorIndex[leftCol];
		minind = leftCol;
		for (col = leftCol + 1; col < rightCol; col++) {
			if (vectorIndex[col] < min) {
				min = vectorIndex[col];
				minind = col;
			}
		}

		minanswer = 0;
		for (col = leftCol; col < rightCol; col++) {
			if (vectorIndex[col] == min) {
				minanswer++;
			}
		}

		if (minanswer > 1)
		{            // At least two columns with the same first row to '1'.
			for (col = minind + 1; col < rightCol; col++) {
				if (vectorIndex[col] == min)
				{        // Add first column which has '1' in the same row to
						 // the other columns so they have '0' in this row after
						 // this operation.
					coladd(ag->matrixXmY, ag->matrixV, ag->matrixV1, ag->matrixV2, minind, col);
				}
			}
		}
		else {
			colexchange(ag->matrixXmY, ag->matrixV, ag->matrixV1, ag->matrixV2, minind, leftCol);
			leftCol++;
		}
	}
}


/* get inverse of biQuadrCoeff (mod currentPrime) wrt currentPrime*/
static int CalcInverseA(int currentPrime) {

	long long Rem = mpz_tdiv_ui(ZT(biQuadrCoeff), currentPrime);  // faster??
	int inverseA = intModInv((int)Rem, currentPrime);
	return inverseA;
}

/*********************************/
/* Sieve thread                  */
/* return factor found in result */
/*********************************/
static void sieveThread(Znum &result) {
	int polySet;
	Znum biT, biU, biV, biR;
	PrimeSieveData *rowPrimeSieveData;
	const PrimeSieveData *rowPrimeSieveData0;
	PrimeTrialDivisionData *rowPrimeTrialDivisionData;
	short SieveArray[2 * MAX_SIEVE_LIMIT];
	int rowPartials[200];
	Znum biLinearCoeff, biDividend, biAbsLinearCoeff;

	int indexFactorsA[50];
	int rowSquares[200];
	int polynomialsPerThread = ((NbrPolynomials - 1) / numberThreads) & 0xFFFFFFFE;
	int firstPolynomial = threadNumber*polynomialsPerThread;
	int lastPolynomial = firstPolynomial + polynomialsPerThread;
	firstPolynomial |= 1;
	int grayCode = firstPolynomial ^ (firstPolynomial >> 1);
	firstPolynomial++;
	int i, PolynomialIndex, index, index2;
	int currentPrime;
	int RemB, D, Q;
	Znum *PtrLinearDelta;
	int rowMatrixBbeforeMerge[200];
	int rowMatrixB[200];
	unsigned char positive;
	int inverseA, twiceInverseA;

	biLinearCoeff = 0; //memset(biLinearCoeff, 0, sizeof(biLinearCoeff));
					   //  synchronized(amodq)
	{
		if (threadNumber == 0) {
			//      primeSieveData = this.primeSieveData;
		}
		else {
			for (i = 0; i<nbrFactorBasePrimes; i++) {
				rowPrimeSieveData = &primeSieveData[i];
				rowPrimeSieveData0 = &primeSieveData[i];
				rowPrimeSieveData->value = rowPrimeSieveData0->value;
				rowPrimeSieveData->modsqrt = rowPrimeSieveData0->modsqrt;
				memcpy(rowPrimeSieveData->Bainv2, rowPrimeSieveData0->Bainv2, sizeof(rowPrimeSieveData0->Bainv2));
			}
		}
		//    threadArray[threadNumber] = Thread.currentThread();
		//    amodq.notifyAll();
	}               // End synchronized block.
	std::cout << "firstPolynomial = " << firstPolynomial;
	std::cout << "  lastPolynomial = " << lastPolynomial << '\n';
	for (polySet = 1;; polySet++) {  // For each polynomial set...
		//if (getTerminateThread())
		//{
		//  throw new ArithmeticException();
		//}
		//synchronized(amodq)
		{
			nbrThreadFinishedPolySet++;
			if (congruencesFound >= matrixBLength/* || factorSiqs != null*/)
			{
				if (nbrThreadFinishedPolySet < polySet * numberThreads) {
					return;
				}
				if (1 /* factorSiqs == null */) {
					/* LinearAlgebraPhase called repeatedly until it returns true */
					while (!LinearAlgebraPhase(biT, biR, biU));

					result = biT;    // Factor found. copy its value to result
#if 0
					synchronized(matrixB)
					{
						factorSiqs = result;
						matrixB.notify();
					}
#endif
				}
				else {
#if 0
					synchronized(matrixB)
					{
						matrixB.notify();
					}
#endif
				}

				return;
			}

			if (nbrThreadFinishedPolySet == polySet * numberThreads) {
				/*********************************************/
				/* Initialization stage for first polynomial */
				/*********************************************/
				firstPrimeSieveData = &primeSieveData[0];
				oldSeed = newSeed;
				newSeed = getFactorsOfA(oldSeed, aindex);
				for (index = 0; index<nbrFactorsA; index++)
				{                        // Get the values of the factors of A.
					afact[index] = primeSieveData[aindex[index]].value;
				}

				// Compute the leading coefficient in biQuadrCoeff.
				biQuadrCoeff = afact[0];

				for (index = 1; index < nbrFactorsA; index++) {
					biQuadrCoeff *= afact[index];
					//MultBigNbrByInt(biQuadrCoeff, afact[index], biQuadrCoeff);
				}

#if DEBUG_SIQS
				if (debug_ctr < 100){
					std::cout << "biQuadrCoeff = " << biQuadrCoeff << '\n';
					debug_ctr++;
				}
#endif
				for (index = 0; index < nbrFactorsA; index++) {
					currentPrime = afact[index];
					// D = (biQuadrCoeff%(currentPrime*currentPrime))/currentPrime
					D = (int)RemDivZnumByInt(biQuadrCoeff,
						currentPrime*currentPrime) / currentPrime;
					Q = primeSieveData[aindex[index]].modsqrt *
						intModInv(D, currentPrime) % currentPrime;
					amodq[index] = D << 1;
					tmodqq[index] = (int)RemDivZnumByInt(zModulus,
						currentPrime*currentPrime);
					if (Q + Q > currentPrime) {
						Q = currentPrime - Q;
					}
					 biDividend = biQuadrCoeff/currentPrime; // DivBigNbrByInt(biQuadrCoeff, currentPrime, biDividend); 
					biLinearDelta[index] = biDividend * Q;
					// MultBigNbrByInt(biDividend, Q, biLinearDelta[index]);  // biLinearDelta[index] = biDividend * Q
#if 0 // DEBUG_SIQS
					std::cout << "index = " << index << " Delta = " << biLinearDelta[index] << '\n';
#endif
				}
				for (index = 1; index < nbrFactorBasePrimes; index++) {
					double dRem, dCurrentPrime;
					long long Rem;
					rowPrimeTrialDivisionData = &ag->primeTrialDivisionData[index];
					rowPrimeSieveData = &primeSieveData[index];


					currentPrime = rowPrimeTrialDivisionData->value;     // Get current prime.
					dCurrentPrime = (double)currentPrime;

					// Get A mod current prime.
					inverseA = CalcInverseA(currentPrime);

					twiceInverseA = inverseA << 1;       // and twice this value.
					dRem = (double)twiceInverseA * (double)rowPrimeSieveData->modsqrt;
					dRem -= floor(dRem / dCurrentPrime)*dCurrentPrime;
					rowPrimeSieveData->difsoln = (int)dRem;
#if 0 // DEBUG_SIQS
					if (debug_ctr < 500) {
						std::cout << "index = " << index;
						std::cout << " dRem = " << dRem << '\n';
						debug_ctr++;
					}
#endif
					for (index2 = nbrFactorsA - 1; index2 > 0; index2--) {
						PtrLinearDelta = &biLinearDelta[index2];
						Rem = mpz_tdiv_ui(ZT(biLinearDelta[index2]), currentPrime); // faster??
						Rem *= twiceInverseA;
						Rem %= currentPrime;
						rowPrimeSieveData->Bainv2[index2 - 1] = (int)Rem;
#if 0 // DEBUG_SIQS
						if (debug_ctr < 100) {
							std::cout << "index2 = " << index2;
							std::cout << " Rem = " << Rem << '\n';
							debug_ctr++;
						}
#endif
					}

					Rem = mpz_tdiv_ui(ZT(biLinearDelta[0]), currentPrime); // faster?? 
					Rem *= twiceInverseA;
					Rem %= currentPrime;
					rowPrimeSieveData->Bainv2_0 = (int)Rem;
					if (rowPrimeSieveData->Bainv2_0 != 0) {
						rowPrimeSieveData->Bainv2_0 =
							currentPrime - rowPrimeSieveData->Bainv2_0;
					}
#if 0 // DEBUG_SIQS
					if (debug_ctr < 500) {
						std::cout << " Rem = " << Rem << '\n';
						debug_ctr++;
					}
#endif
				}

				for (index2 = 0; index2 < nbrFactorsA; index2++)
				{
					primeSieveData[aindex[index2]].difsoln = -1; // Do not sieve.
				}
				//synchronized(TestNbr2)
				{
					//TestNbr2.notifyAll();
				}
			}           // End initializing first polynomial
		}             // End synchronized
#if 0
		synchronized(TestNbr2)
		{
			while (nbrThreadFinishedPolySet < polySet * numberThreads)
			{
				try
				{
					TestNbr2.wait();
				}
				catch (InterruptedException ie) {}
			}
		}
#endif
		if (/*factorSiqs != null ||*/ congruencesFound >= matrixBLength) {
			if (nbrThreadFinishedPolySet > numberThreads*polySet) {
				continue;
			}
			//synchronized(amodq)
			{
				nbrThreadFinishedPolySet++;
			}
			return;
		}

		PolynomialIndex = firstPolynomial;
		// Compute first polynomial parameters.
		biLinearCoeff = biLinearDelta[0];

		for (i = 1; i<nbrFactorsA; i++) {
			if ((grayCode & (1 << i)) == 0) {
				 biLinearCoeff += biLinearDelta[i];
				//AddBigNbrB(biLinearCoeff, biLinearDelta[i], biLinearCoeff);
			}
			else {
				 biLinearCoeff -= biLinearDelta[i];
				//SubtractBigNbrB(biLinearCoeff, biLinearDelta[i], biLinearCoeff);
			}
		}

		
		if (biLinearCoeff < 0) 	{   // Number is negative.
			positive = false;
			biAbsLinearCoeff = abs(biLinearCoeff);
		}
		else {
			positive = true;                    // B is positive.
			biAbsLinearCoeff = biLinearCoeff;   // Get B mod current prime. 
		}

#if 0 // DEBUG_SIQS
		if (debug_ctr < 1000) {
			std::cout << "polynomialIndex = " << PolynomialIndex;
			std::cout << " biAbsLinearCoeff = " << biAbsLinearCoeff << '\n';
			debug_ctr++;
		}
#endif
		for (i = nbrFactorBasePrimes - 1; i>0; i--) {
			double dRem, dCurrentPrime;
			rowPrimeSieveData = &primeSieveData[i];
			rowPrimeSieveData0 = firstPrimeSieveData + i;
			rowPrimeSieveData->difsoln = rowPrimeSieveData0->difsoln;
			rowPrimeSieveData->Bainv2_0 = rowPrimeSieveData0->Bainv2_0;
			rowPrimeTrialDivisionData = &ag->primeTrialDivisionData[i];
			currentPrime = rowPrimeTrialDivisionData->value;     // Get current prime.
			dCurrentPrime = (double)currentPrime;

			inverseA = CalcInverseA(currentPrime);

			RemB = (int)mpz_tdiv_ui(ZT(biAbsLinearCoeff), currentPrime);  // faster??
			if (positive) {
				RemB = currentPrime - RemB;
			}
			dRem = (double)inverseA * (double)(rowPrimeSieveData0->modsqrt + RemB) +
				(double)SieveLimit;
			dRem -= floor(dRem / dCurrentPrime)*dCurrentPrime;
			rowPrimeSieveData->soln1 = (int)dRem;
#if 0 // DEBUG_SIQS
			if (debug_ctr < 1000) {
				std::cout << "i = " << i;
				std::cout << " dRem = " << dRem << '\n';
				debug_ctr++;
			}
#endif
		}

		do {                       // For each polynomial...
			if (congruencesFound >= matrixBLength /*|| factorSiqs != null*/)
			{
				if (nbrThreadFinishedPolySet > numberThreads*polySet) {
					break;
				}
				//synchronized(amodq)
				{
					nbrThreadFinishedPolySet++;
				}
				return;             // Another thread finished factoring.
			}
			if (onlyFactoring) {
				polynomialsSieved += 2;
			}
			/***************/
			/* Sieve stage */
			/***************/
			PerformSiqsSieveStage(&primeSieveData[0], SieveArray,
				PolynomialIndex,
				biLinearCoeff);
			ValuesSieved += 2 * SieveLimit;
			/************************/
			/* Trial division stage */
			/************************/
			index2 = 2 * SieveLimit - 1;
			do {
				index2 -= 16;
				if (((SieveArray[index2 + 1] | SieveArray[index2 + 2] |
					SieveArray[index2 + 3] | SieveArray[index2 + 4] |
					SieveArray[index2 + 5] | SieveArray[index2 + 6] |
					SieveArray[index2 + 7] | SieveArray[index2 + 8] |
					SieveArray[index2 + 9] | SieveArray[index2 + 10] |
					SieveArray[index2 + 11] | SieveArray[index2 + 12] |
					SieveArray[index2 + 13] | SieveArray[index2 + 14] |
					SieveArray[index2 + 15] | SieveArray[index2 + 16]) & 0x8080) != 0)
				{
					for (i = 16; i>0; i--) {
						if ((SieveArray[index2 + i] & 0x80) != 0) {
							if (congruencesFound >= matrixBLength) 	{       
								index2 = 0;
								break;  // All congruences were found: stop sieving.
							}
							SieveLocationHit(rowMatrixB,
								rowMatrixBbeforeMerge,
								index2 + i, &primeSieveData[0],
								rowPartials,
								rowSquares,
								biDividend, biT,
								biLinearCoeff, biR, biU, biV,
								indexFactorsA, false);
							if (congruencesFound >= matrixBLength)	{              
								index2 = 0;
								break;   // All congruences were found: stop sieving.
							}
						}
						if (SieveArray[index2 + i] < 0) {
							if (congruencesFound >= matrixBLength) 	{       
								index2 = 0;
								break;  // All congruences were found: stop sieving.
							}
							SieveLocationHit(rowMatrixB,
								rowMatrixBbeforeMerge,
								index2 + i, &primeSieveData[0],
								rowPartials,
								rowSquares,
								biDividend, biT,
								biLinearCoeff, biR, biU, biV,
								indexFactorsA, true);
							if (congruencesFound >= matrixBLength) {               
								index2 = 0;
								break;  // All congruences were found: stop sieving.
							}
						}
					}
				}
			} while (index2 > 0);
			/*******************/
			/* Next polynomial */
			/*******************/
			PolynomialIndex += 2;
		} while (PolynomialIndex <= lastPolynomial &&
			congruencesFound < matrixBLength);
	}
#if 0


catch (ArithmeticException ae)
{
	synchronized(matrixB)
	{
		factorSiqs = null;
		matrixB.notify();
	}
}
#endif
}

#if 0
/* Implementation of algorithm explained in Gower and Wagstaff paper */
/* The variables with suffix 3 correspond to multiplier = 3 */
static int SQUFOF(long N, int queue[])
{
	double sqrt;
	int Q, Q1, P, P1, L, S;
	int i, j, r, s, t, q;
	int queueHead, queueTail, queueIndex;
	long N3;
	int Q3, Q13, P3, P13, L3, S3;
	int r3, s3, t3, q3;
	int queueHead3, queueTail3, queueIndex3;
	int QRev, Q1Rev, PRev, P1Rev;
	int tRev, qRev, uRev;
	/* Step 1: Initialize */
	N3 = 3 * N;
	if ((N & 3) == 1)
	{
		N <<= 1;
	}
	if ((N3 & 3) == 1)
	{
		N3 <<= 1;
	}
	sqrt = sqrt(N);
	S = (int)sqrt;
	if ((long)(S + 1)*(long)(S + 1) <= N)
	{
		S++;
	}
	if ((long)S*(long)S > N)
	{
		S--;
	}
	if ((long)S*(long)S == N)
	{
		return S;
	}
	Q1 = 1;
	P = S;
	Q = (int)N - P*P;
	L = (int)(2 * sqrt(2 * sqrt));
	queueHead = 0;
	queueTail = 0;

	sqrt = sqrt(N3);
	S3 = (int)sqrt;
	if ((long)(S3 + 1)*(long)(S3 + 1) <= N3)
	{
		S3++;
	}
	if ((long)S3*(long)S3 > N3)
	{
		S3--;
	}
	if ((long)S3*(long)S3 == N3)
	{
		return S3;
	}
	Q13 = 1;
	P3 = S3;
	Q3 = (int)N3 - P3*P3;
	L3 = (int)(2 * sqrt(2 * sqrt));
	queueHead3 = 100;
	queueTail3 = 100;

	/* Step 2: Cycle forward to find a proper square form */
	for (i = 0; i <= L; i++)
	{
		/* Multiplier == 1 */
		q = (S + P) / Q;
		P1 = q*Q - P;
		if (Q <= L)
		{
			if ((Q & 1) == 0)
			{
				queue[queueHead++] = Q >> 1;
				queue[queueHead++] = P % (Q >> 1);
				if (queueHead == 100)
				{
					queueHead = 0;
				}
			}
			else if (Q + Q <= L)
			{
				queue[queueHead++] = Q;
				queue[queueHead++] = P % Q;
				if (queueHead == 100)
				{
					queueHead = 0;
				}
			}
		}
		t = Q1 + q*(P - P1);
		Q1 = Q;
		Q = t;
		P = P1;
		{
			r = (int)sqrt(Q);
			if (r*r == Q)
			{
				queueIndex = queueTail;
				for (;;)
				{
					if (queueIndex == queueHead)
					{
						/* Step 3: Compute inverse square root of the square form */
						PRev = P;
						Q1Rev = r;
						uRev = (S - PRev) % r;
						uRev += (uRev >> 31) & r;
						PRev = S - uRev;
						QRev = (int)((N - (long)PRev*(long)PRev) / Q1Rev);
						/* Step 4: Cycle in the reverse direction to find a factor of N */
						for (j = i; j >= 0; j--)
						{
							qRev = (S + PRev) / QRev;
							P1Rev = qRev*QRev - PRev;
							if (PRev == P1Rev)
							{
								/* Step 5: Get the factor of N */
								if ((QRev & 1) == 0)
								{
									return QRev >> 1;
								}
								return QRev;
							}
							tRev = Q1Rev + qRev*(PRev - P1Rev);
							Q1Rev = QRev;
							QRev = tRev;
							PRev = P1Rev;
							qRev = (S + PRev) / QRev;
							P1Rev = qRev*QRev - PRev;
							if (PRev == P1Rev)
							{
								/* Step 5: Get the factor of N */
								if ((QRev & 1) == 0)
								{
									return QRev >> 1;
								}
								return QRev;
							}
							tRev = Q1Rev + qRev*(PRev - P1Rev);
							Q1Rev = QRev;
							QRev = tRev;
							PRev = P1Rev;
						}
						break;
					}
					s = queue[queueIndex++];
					t = queue[queueIndex++];
					if (queueIndex == 100)
					{
						queueIndex = 0;
					}
					if ((P - t) % s == 0)
					{
						break;
					}
				}
				if (r > 1)
				{
					queueTail = queueIndex;
				}
				if (r == 1)
				{
					queueIndex = queueTail;
					for (;;)
					{
						if (queueIndex == queueHead)
						{
							break;
						}
						if (queue[queueIndex] == 1)
						{
							return 0;
						}
						queueIndex += 2;
						if (queueIndex == 100)
						{
							queueIndex = 0;
						}
					}
				}
			}
		}
		q = (S + P) / Q;
		P1 = q*Q - P;
		if (Q <= L)
		{
			if ((Q & 1) == 0)
			{
				queue[queueHead++] = Q >> 1;
				queue[queueHead++] = P % (Q >> 1);
				if (queueHead == 100)
				{
					queueHead = 0;
				}
			}
			else if (Q + Q <= L)
			{
				queue[queueHead++] = Q;
				queue[queueHead++] = P % Q;
				if (queueHead == 100)
				{
					queueHead = 0;
				}
			}
		}
		t = Q1 + q*(P - P1);
		Q1 = Q;
		Q = t;
		P = P1;

		/* Multiplier == 3 */
		q3 = (S3 + P3) / Q3;
		P13 = q3*Q3 - P3;
		if (Q3 <= L3)
		{
			if ((Q3 & 1) == 0)
			{
				queue[queueHead3++] = Q3 >> 1;
				queue[queueHead3++] = P3 % (Q3 >> 1);
				if (queueHead3 == 200)
				{
					queueHead3 = 100;
				}
			}
			else if (Q3 + Q3 <= L3)
			{
				queue[queueHead3++] = Q3;
				queue[queueHead3++] = P3 % Q3;
				if (queueHead3 == 200)
				{
					queueHead3 = 100;
				}
			}
		}
		t3 = Q13 + q3*(P3 - P13);
		Q13 = Q3;
		Q3 = t3;
		P3 = P13;
		{
			r3 = (int)sqrt(Q3);
			if (r3*r3 == Q3)
			{
				queueIndex3 = queueTail3;
				for (;;)
				{
					if (queueIndex3 == queueHead3)
					{
						/* Step 3: Compute inverse square root of the square form */
						PRev = P3;
						Q1Rev = r3;
						uRev = (S3 - PRev) % r3;
						uRev += (uRev >> 31) & r3;
						PRev = S3 - uRev;
						QRev = (int)((N3 - (long)PRev*(long)PRev) / Q1Rev);
						/* Step 4: Cycle in the reverse direction to find a factor of N */
						for (j = i; j >= 0; j--)
						{
							qRev = (S3 + PRev) / QRev;
							P1Rev = qRev*QRev - PRev;
							if (PRev == P1Rev)
							{
								/* Step 5: Get the factor of N */
								if ((QRev & 1) == 0)
								{
									return QRev >> 1;
								}
								return QRev;
							}
							tRev = Q1Rev + qRev*(PRev - P1Rev);
							Q1Rev = QRev;
							QRev = tRev;
							PRev = P1Rev;
							qRev = (S3 + PRev) / QRev;
							P1Rev = qRev*QRev - PRev;
							if (PRev == P1Rev)
							{
								/* Step 5: Get the factor of N */
								if ((QRev & 1) == 0)
								{
									return QRev >> 1;
								}
								return QRev;
							}
							tRev = Q1Rev + qRev*(PRev - P1Rev);
							Q1Rev = QRev;
							QRev = tRev;
							PRev = P1Rev;
						}
						break;
					}
					s3 = queue[queueIndex3++];
					t3 = queue[queueIndex3++];
					if (queueIndex3 == 200)
					{
						queueIndex3 = 100;
					}
					if ((P3 - t3) % s3 == 0)
					{
						break;
					}
				}
				if (r3 > 1)
				{
					queueTail3 = queueIndex3;
				}
				if (r3 == 1)
				{
					queueIndex3 = queueTail3;
					for (;;)
					{
						if (queueIndex3 == queueHead3)
						{
							break;
						}
						if (queue[queueIndex3] == 1)
						{
							return 0;
						}
						queueIndex3 += 2;
						if (queueIndex3 == 200)
						{
							queueIndex3 = 100;
						}
					}
				}
			}
		}
		q3 = (S3 + P3) / Q3;
		P13 = q3*Q3 - P3;
		if (Q3 <= L3)
		{
			if ((Q3 & 1) == 0)
			{
				queue[queueHead3++] = Q3 >> 1;
				queue[queueHead3++] = P3 % (Q3 >> 1);
				if (queueHead3 == 200)
				{
					queueHead3 = 100;
				}
			}
			else if (Q3 + Q3 <= L3)
			{
				queue[queueHead3++] = Q3;
				queue[queueHead3++] = P3 % Q3;
				if (queueHead3 == 200)
				{
					queueHead3 = 100;
				}
			}
		}
		t3 = Q13 + q3*(P3 - P13);
		Q13 = Q3;
		Q3 = t3;
		P3 = P13;
	}
	return 0;
}
/* If 2^value mod value = 2, then the value is a probable prime (value odd) */
static unsigned char isProbablePrime(long value)
{
	long mask, montgomery2, Pr, Prod0, Prod1, MontDig;
	int N, MontgomeryMultN;
	long BaseHI, BaseLO, valueHI, valueLO;
	int x = N = (int)value; // 2 least significant bits of inverse correct.
	x = x * (2 - N * x);     // 4 least significant bits of inverse correct.
	x = x * (2 - N * x);     // 8 least significant bits of inverse correct.
	x = x * (2 - N * x);     // 16 least significant bits of inverse correct.
	x = x * (2 - N * x);     // 32 least significant bits of inverse correct.
	MontgomeryMultN = (-x) & 0x7FFFFFFF;
	mask = 1L << 62;
	montgomery2 = 2 * (mask % value);
	if (montgomery2 >= value)
	{
		montgomery2 -= value;
	}
	BaseHI = (int)(montgomery2 >> 31);
	BaseLO = (int)montgomery2 & 0x7FFFFFFF;
	valueHI = (int)(value >> 31);
	valueLO = (int)value & 0x7FFFFFFF;
	while ((mask & value) == 0)
	{
		mask >>= 1;
	}
	mask >>= 1;
	while (mask > 0)
	{
		/* Square the base */
		Pr = BaseLO * BaseLO;
		MontDig = ((int)Pr * MontgomeryMultN) & 0x7FFFFFFFL;
		Prod0 = (Pr = ((MontDig * valueLO + Pr) >> 31) +
			MontDig * valueHI + BaseLO * BaseHI) & 0x7FFFFFFFL;
		Prod1 = Pr >> 31;
		Pr = BaseHI * BaseLO + Prod0;
		MontDig = ((int)Pr * MontgomeryMultN) & 0x7FFFFFFFL;
		Prod0 = (Pr = ((MontDig * valueLO + Pr) >> 31) +
			MontDig * valueHI + BaseHI * BaseHI + Prod1) & 0x7FFFFFFFL;
		Prod1 = Pr >> 31;
		if (Prod1 > valueHI || (Prod1 == valueHI && Prod0 >= valueLO))
		{
			Prod0 = (Pr = Prod0 - valueLO) & 0x7FFFFFFFL;
			Prod1 = ((Pr >> 31) + Prod1 - valueHI) & 0x7FFFFFFFL;
		}
		BaseLO = Prod0;
		BaseHI = Prod1;

		if ((mask & value) != 0)
		{
			/* Multiply by 2 */
			Pr = 2 * ((BaseHI << 31) + BaseLO);
			if (Pr >= value)
			{
				Pr -= value;
			}
			BaseHI = (int)(Pr >> 31);
			BaseLO = (int)Pr & 0x7FFFFFFF;
		}
		mask >>= 1;
	}
	Pr = (BaseHI << 31) + BaseLO;
	return Pr == montgomery2;
}

#endif

