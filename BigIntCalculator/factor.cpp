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

#define __EMSCRIPTEN__

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cstdint>
#include <cmath>
#include <cassert>
#include "showtime.h"
#include "bignbr.h"

#include "factor.h"
extern int EC;            // Elliptic Curve Number
//int yieldFreq;

//char lowerText[30000];
//char *ptrLowerText;
extern mmCback modmultCallback;


extern bool *primeFlags;
extern unsigned long long *primeList;
extern unsigned int prime_list_count;
void generatePrimes(unsigned long long int max_val);

static int FactorIndex;

static void insertBigFactor(std::vector<zFactors> &Factors, Znum &divisor);
void ValuestoZ(Znum &numberZ, const int number[]);


long long Gamma[386];
long long Delta[386];
long long AurifQ[386];


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

void printfactors(const std::vector <zFactors> Factors) {
	for (int i = 0; i < Factors.size(); i++) {
		std::cout << Factors[i].Factor << "^" << Factors[i].exponent << " ("
			<< Factors[i].upperBound << ")  * " ;
	}
	std::cout << '\n';
}

static void GetAurifeuilleFactor(std::vector<zFactors> &pstFactors, int L, const Znum &BigBase,
	const int DegreeAurif) {
	Znum x, Csal, Dsal, Nbr1;  
	int k;

	BigIntPowerIntExp(BigBase, L, x);   // x <- BigBase^L.
	Csal = 1;      //intToBigInteger(Csal, 1);
	Dsal = 1;       // intToBigInteger(Dsal, 1);
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
	Csal += Nbr1; // BigIntAdd(Csal, Nbr1, Csal);        // Csal <- Csal * x + Gamma[k]
	BigIntPowerIntExp(BigBase, (L + 1) / 2, Nbr1);   // Nbr1 <- Dsal * base^((L+1)/2)
	Nbr1 = Nbr1 * Dsal; //BigIntMultiply(Dsal, Nbr1, Nbr1);
	Dsal = Csal + Nbr1; // BigIntAdd(Csal, Nbr1, Dsal);
 
	insertBigFactor(pstFactors, Dsal);
	Dsal = Csal - Nbr1; // BigIntSubt(Csal, Nbr1, Dsal);

	insertBigFactor(pstFactors, Dsal);
	return;
}

/* Get Aurifeuille factors.
 see https://en.wikipedia.org/wiki/Aurifeuillean_factorization */
static void InsertAurifFactors(std::vector<zFactors> &pstFactors, const Znum &BigBase,
	int Expon, int Incre)
{
	int DegreeAurif;
	if (BigBase >= 386) {
		return;    // Base is very big, so go out.
	}
	auto Base = MulPrToLong(BigBase);

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
				GetAurifeuilleFactor(pstFactors, L, BigBase, DegreeAurif);
				if (q != L * L) {
					GetAurifeuilleFactor(pstFactors, q / L, BigBase, DegreeAurif);
				}
			}
			L += 2;
		}
	}
	return;
}


static void Cunningham(std::vector<zFactors> &pstFactors, const Znum &BigBase, int Expon,
	int increment, const Znum &BigOriginal) {
	int Expon2, k;
	Znum Nbr1, Nbr2, Temp1;      

	Expon2 = Expon;

	while (Expon2 % 2 == 0 && increment == -1) {
		Expon2 /= 2;
		BigIntPowerIntExp(BigBase, Expon2, Nbr1);
		Nbr1 += increment; //addbigint(Nbr1, increment);
		insertBigFactor(pstFactors, Nbr1);
		InsertAurifFactors(pstFactors, BigBase, Expon2, 1);
	}

	k = 1;
	while (k * k <= Expon) {
		if (Expon % k == 0) {
			if (k % 2 != 0) { /* Only for odd exponent */
				BigIntPowerIntExp(BigBase, Expon / k, Nbr1);
				Nbr1 += increment;      //addbigint(Nbr1, increment);
				Nbr2 = gcd(Nbr1, BigOriginal);   // Nbr2 <- gcd(Base^(Expon/k)+incre, original)
				insertBigFactor(pstFactors, Nbr2);
				Temp1 = BigOriginal; // CopyBigInt(Temp1, *BigOriginal);
				Nbr1 = Temp1 / Nbr2; // BigIntDivide(Temp1, Nbr2, Nbr1);
				insertBigFactor(pstFactors, Nbr1);
				InsertAurifFactors(pstFactors, BigBase, Expon / k, increment);
			}

			if ((Expon / k) % 2 != 0) { /* Only for odd exponent */
				BigIntPowerIntExp(BigBase, k, Nbr1);
				Nbr1 += increment; //addbigint(Nbr1, increment);
				Nbr2 = gcd(Nbr1, BigOriginal);   // Nbr2 <- gcd(Base^k+incre, original)
				insertBigFactor(pstFactors, Nbr2);
				Temp1 = BigOriginal; // CopyBigInt(Temp1, *BigOriginal);
				Nbr1 = Temp1 / Nbr2; // BigIntDivide(Temp1, Nbr2, Nbr1);
				insertBigFactor(pstFactors, Nbr1);
				InsertAurifFactors(pstFactors, BigBase, k, increment);
			}
		}
		k++;
	}
	return;
}

static bool ProcessExponent(std::vector<zFactors> &pstFactors, const Znum &nbrToFactor, int Exponent)
{
#ifdef __EMSCRIPTEN__
	char status[200] = { 0 };
	char *ptrStatus;
#endif
	Znum NFp1, NFm1, nthRoot, rootN1, rootN, rootbak;   // follow advice not to use stack for these
	Znum nextroot, dif;     // follow advice not to use stack for these
	Znum Temp1;
	//double log2N;

#ifdef __EMSCRIPTEN__
	int elapsedTime = (int)(tenths() - originalTenthSecond);
#ifndef _DEBUG
	if (elapsedTime / 10 != oldTimeElapsed / 10)
#endif
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
		printf("%s\n", status);
	}
#endif

	NFp1 = nbrToFactor+1;// NFp1 <- NumberToFactor + 1
	NFm1 = nbrToFactor-1;  // NFm1 <- NumberToFactor - 1
           
	mpz_nthroot(ZT(nthRoot), ZT(NFp1), Exponent);  // Find nth root of number to factor.
	rootbak = nthRoot;

	for (;;) {
		BigIntPowerIntExp(nthRoot, Exponent - 1, rootN1); // rootN1 <- nthRoot ^ (Exponent-1)
		rootN = nthRoot*rootN1; //BigIntMultiply(nthRoot, rootN1, rootN);     // rootN <- nthRoot ^ Exponent
		dif = NFp1 - rootN;     // BigIntSubt(NFp1, rootN, dif);    
		if (dif == 0) {         // Perfect power-1
			Cunningham(pstFactors, nthRoot, Exponent, -1, nbrToFactor);
			return true;
		}
		dif++;        // dif <- dif + 1
		Temp1 = dif / rootN1;
		Temp1 /= Exponent; //subtractdivide(Temp1, 0, Exponent);        // Temp1 <- Temp1 / Exponent
		nextroot = Temp1 + nthRoot;         // 
		nextroot--;  //   nextroot = nextroot - 1             
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
		dif++;                        // dif <- dif + 1
		Temp1 = dif / rootN1;
		Temp1 /= Exponent; //subtractdivide(Temp1, 0, Exponent);   // Temp1 <- Temp1 / Exponent
		nextroot = Temp1 + nthRoot;
		nextroot--;             // nextroot <- nextroot - 1
		nthRoot = nextroot - nthRoot;
		if (nthRoot >= 0) {
			break;               // Not a perfect power
		}
		nthRoot = nextroot;      // CopyBigInt(nthRoot, nextroot);
	}
	return false;
}


/* check whether the number +/- 1 is a perfect power*/
static void PowerPM1Check(std::vector<zFactors> &pstFactors, const Znum &nbrToFactor)
{
	bool plus1 = false;
	bool minus1 = false;
	int Exponent = 0;
	int i, j;
	int modulus;
	int mod9 = (int)MulPrToLong(nbrToFactor % 9); // getRemainder(*nbrToFactor, 9);
	int maxExpon = (int)mpz_sizeinbase(ZT(nbrToFactor), 2);
	int numPrimes = 2 * maxExpon + 3;
	/* if numPrimes >= 66472 i.e. maxExpon > 33233 */
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
	Znum Temp1, Temp2;

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
			remainder = MulPrToLong(Temp2);
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


/* sort factors into ascending order. If two factors are equal, merge them by adding
the 2nd exponent to the 1st then removing the second entry by moving any following entries
up 1 position and adjusting the count of the total mumber of entries. */
static void SortFactors(std::vector<zFactors> &Factors) {
	auto lastfactor = Factors.size();
	zFactors temp;
	for (int i = 0; i < lastfactor-1; i++)
		for (int j = i + 1; j < lastfactor; j++) {
			if (Factors[i].Factor > Factors[j].Factor) {
			/* factors out of sequence so swap them*/
				temp = Factors[i];        // shallow copy is OK for swap
				Factors[i] = Factors[j];
				Factors[j] = temp;
			}
			else if (Factors[i].Factor == Factors[j].Factor) {
				/* factors are equal so merge them */
				Factors[i].exponent += Factors[j].exponent;
				if (Factors[j].upperBound == 0)
					Factors[i].upperBound = 0;
				else if (Factors[i].upperBound != 0
					&& Factors[i].upperBound < Factors[j].upperBound)
					Factors[i].upperBound = Factors[j].upperBound;
				/* now move index entries higher than j down 1 */
				for (int k = j; k < lastfactor-1; k++) {
					Factors[k].exponent = Factors[k + 1].exponent;
					Factors[k].Factor = Factors[k + 1].Factor;
					Factors[k].upperBound = Factors[k + 1].upperBound;
				}
				/* now adjust counters */
				j--;
				lastfactor--;
				Factors.pop_back();  // remove last entry
			}
		}
#ifdef _DEBUG
	std::cout << "result after sort" << '\n';
	printfactors(Factors);
#endif
}

/* Insert new factor found into factor array. This factor array must be sorted. 
Dividend is the known factor that divisor is a factor of */
/* assume each small factor is prime.
ix is index of non-prime factor which is a multiple of divisor */
static void insertIntFactor(std::vector<zFactors> &Factors, int pi, int ix) {
	auto lastfactor = Factors.size();
	Znum quot, qnew;
	long long rem;
	long long divisor = primeList[pi];
	int exp = 0;

	quot = Factors[ix].Factor;

	/* divide quot by divisor as many times as possible */
	while (true) {
		rem = mpz_fdiv_q_ui(ZT(qnew), ZT(quot), divisor);
		if (rem != 0)
			break;
		else {
			exp++;
			quot = qnew;
		}
	}
	if (quot != 1) {
		/* add new factor */
		Factors.resize(lastfactor + 1);
		Factors[lastfactor].exponent = exp;
		Factors[lastfactor].Factor = divisor;
		Factors[lastfactor].upperBound = -1; // show that new factor is prime

		Factors[ix].Factor = quot;
		Factors[ix].upperBound = pi;
	}
	else {
		/* replace residue of 1 with new factor */
		Factors[ix].Factor = divisor;
		Factors[ix].exponent *= exp;
		Factors[ix].upperBound = -1;
	}
	SortFactors(Factors);
#ifdef _DEBUG
	/*std::cout << "result after adding factor " << divisor << '\n';
	printfactors(Factors);*/
#endif
}

/* Insert new factor found into factor array. This factor array must be sorted.
 The divisor must be also sorted. */
static void insertBigFactor(std::vector<zFactors> &Factors, Znum &divisor) {
	auto lastfactor = Factors.size();
	auto ipoint = lastfactor;
	zFactors temp;
	Znum g;
#ifdef _DEBUG
	/*std::cout << "InsertBigFactor Divisor =" << divisor << '\n';
	printfactors(Factors);*/
#endif
	for (int i = 0; i < lastfactor; i++) {
		if (Factors[i].Factor == divisor)
			break;  // factor already found
		g = gcd(Factors[i].Factor, divisor);
		if (g != 1 && g < Factors[i].Factor) {
			Factors.resize(lastfactor + 1);
			/* we can replace Factor with 2 factors, Factor/g and g 
			(if Factor is a multiple of divisor g = divisor) */
			Factors[i].Factor /= g;
			Factors[ipoint].Factor = g;
			Factors[ipoint].exponent = Factors[i].exponent;
			Factors[ipoint].upperBound = Factors[i].upperBound;
			ipoint++;
			lastfactor++;
		}
	}
#ifdef _DEBUG
	/*std::cout << "result before sort" << '\n';
	printfactors(Factors);*/
#endif
	SortFactors(Factors);
}


/* called for Carmichael numbers that have no small factors. 
Return: false = No factors found, true = factors found.
 Use: Xaux for square root of -1.
      Zaux for square root of 1. */
static bool factorCarmichael(const Znum &pValue, std::vector<zFactors> &Factors)
{
	Znum randomBase = 0;  // pseudo-random number
	bool factorsFound = false;
	int countdown, ctr;
	bool sqrtOneFound = false;
	bool sqrtMinusOneFound = false;
	Znum Aux1, Aux2, Aux3, Aux4, Xaux, Zaux, Temp4;
	const Znum two = 2;  // for squareing

	//BigIntegerToLimbs(Aux1, pValue, pValue.nbrLimbs);  // copy p to Aux1
	//Aux1[0].x--;       // Aux1 = p - 1 (p is odd, so there is no carry).
	//Aux1Len = nbrLimbsP;
	Aux1 = pValue - 1;
	DivideBigNbrByMaxPowerOf2(ctr, Aux1);  // Aux1 /= 2^ctr
	//BigIntegerToLimbs(TestNbr, pValue, pValue.nbrLimbs);  // copy p to TestNbr
	//memcpy(TestNbr, pValueLimbs, nbrLimbsP * sizeof(limb)); // copy p to TestNbr
	//TestNbr[nbrLimbsP].x = 0;
	//GetMontgomeryParms(nbrLimbsP);
	for (countdown = 20; countdown > 0; countdown--) {
		int i;
		//NumberLength = nbrLimbsP;
		randomBase = ((uint64_t)randomBase * 89547121 + 1762281733) & MAX_INT_NBR;
		//modPowBaseInt(randomBase, Aux1, Aux1Len, Aux2); // Aux2 = base^Aux1 (mod TestNbr)
		mpz_powm(ZT(Aux2), ZT(randomBase), ZT(Aux1), ZT(pValue));
		// If Mult1 = 1 or Mult1 = TestNbr-1, then try next base.
		if (Aux2 ==1 || Aux2 == pValue-1) {
			continue;    // This base cannot find a factor. Try another one.
		}
		for (i = 0; i < ctr; i++)
		{              // Loop that squares number.
			//modmult(Aux2, Aux2, Aux3);  // Aux3 = Aux2^2 (mod p)
			mpz_powm(ZT(Aux3), ZT(Aux2), ZT(two), ZT(pValue));
			if (Aux3 == 1)
			{            // Non-trivial square root of 1 found.
				if (!sqrtOneFound)
				{          // Save it to perform GCD later.
					//memcpy(Zaux, Aux2, nbrLimbsP * sizeof(limb));
					Zaux = Aux2;
					sqrtOneFound = true;
				}
				else
				{          // Try to find non-trivial factor by doing GCD.
					//SubtBigNbrMod(Aux2, Zaux, Aux4); // Aux4 = Aux2 - Zaux (mod p)
					//LimbsToBigInteger(Aux4, Temp2, NumberLength);  // Temp2 = Aux4
					//BigIntGcd(pValue, Temp2, Temp4);  // Temp4 = gcd(p, Temp2)
					Aux4 = Aux2 - Zaux;
					mpz_mod(ZT(Aux4), ZT(Aux4), ZT(pValue));
					Temp4 = gcd(pValue, Aux4);
					if ((Temp4 > 1) && (Temp4 != pValue))
					{          // Non-trivial factor found.
						insertBigFactor(Factors, Temp4);
						factorsFound = true;
					}
				}
				// Try to find non-trivial factor by doing GCD.
				//NumberLength = nbrLimbsP;
				//AddBigNbrMod(Aux2, MontgomeryMultR1, Aux4); // Aux4 = Aux2+1 (mod p)
				//LimbsToBigInteger(Aux4, Temp2, NumberLength);
				//BigIntGcd(pValue, Temp2, Temp4);  // Temp4 = gcd(p, Temp2)
				Aux4 = Aux2 + 1;
				Aux4 %= pValue;
				Temp4 = gcd(pValue, Aux4);
				if ((Temp4 > 1) && (Temp4 != pValue))
				{          // Non-trivial factor found.
					insertBigFactor(Factors, Temp4);
					factorsFound = true;
				}
				i = ctr;
				continue;  // Find more factors.
			}
			if (Aux3 == pValue-1)
			{            // Square root of 1 found.
				if (!sqrtMinusOneFound)
				{          // Save it to perform GCD later.
					//memcpy(Xaux, Aux2, nbrLimbsP * sizeof(limb)); // Xaux = Aux2
					Xaux = Aux2;
					sqrtOneFound = true;
				}
				else
				{          // Try to find non-trivial factor by doing GCD.
					//SubtBigNbrMod(Aux3, Xaux, Aux4);  // Aux4 = Aux3 - Xaux (mod p)
					//LimbsToBigInteger(Aux4, Temp2, NumberLength);  // Temp2 = Aux4
					//BigIntGcd(pValue, Temp2, Temp4);  // Temp4 = gcd(p, Temp2)
					Aux4 = Aux3 - Xaux;
					Aux4 %= pValue;
					Temp4 = gcd(pValue, Aux4);
					if ((Temp4 > 1) && (Temp4 != pValue))
					{          // Non-trivial factor found.

						insertBigFactor(Factors, Temp4);
						factorsFound = true;
					}
				}
				i = ctr;
				continue;  // Find more factors.
			}
			//memcpy(Aux2, Aux3, nbrLimbsP * sizeof(limb)); // Aux2 = Aux3
			Aux2 = Aux3;
		}
	}
	return factorsFound;
}


static bool factor(const Znum &toFactor, std::vector<zFactors> &Factors) {
	int upperBound;
	long long testP;
	bool restart = false;  // set true if trial division has to restart
	/* initialise factor list */
	Factors[0].exponent = 1;
	Factors[0].Factor = toFactor;
	Factors[0].upperBound = 0;
	if (primeFlags == NULL) {
		generatePrimes(300000);  // takes a while, but only needed on 1st call
	}
	
	oldTimeElapsed = 0;
	originalTenthSecond = tenths();    // record start time

	if (toFactor >= 300000LL* 300000LL) {
		/* may not be able to factorise entirely by trial division, so try this first */
		PowerPM1Check(Factors, toFactor);  // check if toFactor is a perfect power +/- 1
#ifdef _DEBUG
		if (Factors.size() > 1) {
			std::cout << "PowerPM1Check result: ";
			printfactors(Factors);
		}
#endif
	}

	do {  /* use trial division */
		restart = false;
		for (int i = 0; i < Factors.size(); i++) {
			upperBound = Factors[i].upperBound;  // resume from where we left off
			if (upperBound == -1)
				continue;  // factor is prime
			/* trial division */
			while (upperBound < (int)prime_list_count) {
				testP = primeList[upperBound];
				if (testP*testP > Factors[i].Factor) {
					Factors[i].upperBound = -1; // show that residue is prime
					break;
				}
				if (Factors[i].Factor%testP == 0) {
					insertIntFactor(Factors, upperBound, i);
					restart = true;
					break;
				}
				Factors[i].upperBound = upperBound;
				upperBound++;
			}
		}
		SortFactors(Factors);  // tidy up factor list
	} while (restart);
#ifdef _DEBUG
	std::cout << "result after trial division " ;
	printfactors(Factors);
#endif

	/* Any small factors (up to 300,000) have now been found by trial division */
	/* Check whether the residue is prime or prime power.
	Given that the residue is less than about 10^10,000 the maximum exponent is 
	less than 2000.  e.g. 3000007^1826 has 10,002 digits */
	Znum Zpower, Zprime;
	int expon;
	EC = 1;  // start with 1st curve
	modmultCallback = showECMStatus;   // Set callback function pointer
	//NumberLength = (int)(mpz_sizeinbase(ZT(toFactor), 2) + BITS_PER_GROUP - 1) / BITS_PER_GROUP;
	//GetYieldFrequency();   //get yield frequency based on NumberLength (used by showECMStatus)

	for (int i = 0; i < Factors.size(); i++) {
		if (Factors[i].upperBound == -1)
			continue;  // skip if factor is known to be prime
		Zprime = Factors[i].Factor;
		testP = primeList[Factors[i].upperBound]; // get largest number used in trial division
		expon = (int)PowerCheck(Zprime, Zpower,  testP - 1); // on return; Zpower^expon = Zprime
		if (expon > 1) { /* if power is a perfect power*/
			Factors[i].Factor = Zpower;
			Factors[i].exponent *= expon;
		}
		
		int result = BpswPrimalityTestNew(Zpower, testP- 1);
		if (result == 0) {   // Number is prime.
			Factors[i].upperBound = -1; // Indicate that number is prime.
			continue;
		}
		if (result > 1) {
			if (factorCarmichael(Zpower, Factors))
				continue;
		}
		auto rv = ecm(Zpower);          // Factor number.
		if (!rv)
			return false;  // failed to factorise number
		// Check whether GD is not one. In this case we found a proper factor.
		int ctr;
		for (ctr = 1; ctr < NumberLength; ctr++) {
			if (GD[ctr].x != 0) {
				break;
			}
		}
		if (ctr != 1 || GD[0].x != 1) {
			/* GD is not 1 */
			int numLimbs;
			numLimbs = NumberLength;
			while (numLimbs > 1) {    // adjust count of number of limbs
				if (GD[numLimbs - 1].x != 0) {
					break;
				}
				numLimbs--;
			}
			LimbsToBigInteger(GD, Temp1, numLimbs); /* copy number from GD and convert it to a Znum */
			Znum Zgd;
			BigtoZ(Zgd, Temp1);
			insertBigFactor(Factors, Zgd);
			i = 0;			// restart loop at beginning!!
		}
	}
	SortFactors(Factors);  // tidy up factor list 
	return true;
}

/* compute 4 or less values the squares of which add up to prime p,  
return values in Mult1, Mult2, Mult3 and Mult4 */
static void ComputeFourSquaresNew(const Znum &p, Znum &Mult1, Znum &Mult2,
	Znum &Mult3, Znum &Mult4) {
	Znum a, q, K, Tmp, Tmp1, Tmp2, Tmp3, Tmp4, M1, M2, M3, M4; 
	Znum TestNbr;

	if (p == 2) {   /* Prime factor is 2 */
		Mult1 = 1;  // 2 = 1^2 + 1^2 + 0^2 + 0^2
		Mult2 = 1;
		Mult3 = 0;
		Mult4 = 0;
	}
	else {       /* Prime factor p is not 2 */
		if ((p & 3) == 1) { /* if p = 1 (mod 4) */
			q = (p-1)/4; // q = (prime-1)/4
	  
			a = 1;
			do {    // Loop that finds mult1^2 = (-1) mod p
				a++; 
				assert(a < p);
				/* Mult1 = a^q(mod p) */
				mpz_powm(ZT(Mult1), ZT(a), ZT(q), ZT(p));
				mpz_powm_ui(ZT(TestNbr), ZT(Mult1), 2, ZT(p));
			} while (TestNbr != p-1 && TestNbr != -1);

			Mult1 = abs(Mult1);
			Mult2 = 1;

			for (;;) {  
				Tmp = Mult1 * Mult1 + Mult2 * Mult2; // K <- (mult1^2 + mult2^2) / p
				K = Tmp / p;    // in this loop K is smaller each time round   
				if (K == 1) {  // are we there yet?
					Mult3 = 0;
					Mult4 = 0;
					break;     // we are finished
				}
				M1 = Mult1 % K;
				if (M1 < 0) {
					M1 += K;
				}
				M2 = Mult2 % K;
				if (M2 < 0) {
					M2 += K;
				}
				Tmp = (K+1)/2; // Tmp <- (K+1) / 2
				//subtractdivide(Tmp, -1, 2);       
				if (M1 >= Tmp) // If M1 >= K / 2 ... 
				{
					M1 -= K;
				}
				if (M2 >= Tmp) {     // If M2 >= K / 2 ...     
					M2 -= K;
				}
				Tmp = Mult1*M1 + Mult2*M2;
				Tmp2 = Tmp / K;  // Tmp2 <- (mult1*m1 + mult2*m2) / K
				Tmp = Mult1*M2 - Mult2*M1;
				Mult2 = Tmp / K;    // Mult2 <- (mult1*m2 - mult2*m1) /K
				Mult1 = Tmp2;       // Mult1 <- (mult1*m1 + mult2*m2) / K
			} /* end while */
		} /* end p = 1 (mod 4) */

		else { /* if p = 3 (mod 4) */
			q = (p-1)/2; //  q = (prime-1)/2
			Mult1 = 0;
			do {
				Mult1++;
				Tmp = Mult1*Mult1 + 1;
				Tmp = -Tmp;
				while (Tmp < 0) Tmp += p;
				mpz_powm(ZT(Tmp1), ZT(Tmp), ZT(q), ZT(p));
				       // At this moment Tmp1 = (-1 - Mult1^2)^((p-1)/2)(Mod p)
			} while (Tmp1 != 1);  // Continue loop if it is not 1.

			// After the loop finishes, Tmp1 = (-1 - Mult1^2) is a quadratic residue mod p.
			q = (p+1)/4;

			// Find Mult2 <- square root of Tmp1 = Tmp^q (mod p) 
			mpz_powm(ZT(Mult2), ZT(Tmp), ZT(q), ZT(p));

			Mult3 = 1;
			Mult4 = 0;

			for (;;) {
				// Compute K <- (Mult1^2 + Mult2^2 + Mult3^2 + Mult4^2) / p
				Tmp = Mult1*Mult1 + Mult2*Mult2;
				Tmp = Tmp + Mult3*Mult3 + Mult4*Mult4;
				assert(Tmp%p == 0);
				K = Tmp / p;    // in this loop K is smaller each time round
				assert(K > 0);
				if (K == 1) {
					break;  // we are done when K equals 1
				}


				if (ZisEven(K)) { // If K is even ...
					if (ZisEven(Mult1) != ZisEven(Mult2))
					{  // If Mult1 + Mult2 is odd...
						if (ZisEven(Mult1) == ZisEven(Mult3))
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
					Tmp1 = (Mult1 + Mult2)/2;   // Tmp1 <- (Mult1 + Mult2) / 2
					Tmp2 = (Mult1 - Mult2)/2;   // Tmp2 <- (Mult1 - Mult2) / 2
					Tmp3 = (Mult3 + Mult4)/2;   // Tmp3 <- (Mult3 + Mult4) / 2
 					Mult4 = (Mult3 - Mult4)/2 ; // Mult4 <- (Mult3 - Mult4) / 2
					Mult3 = Tmp3;
					Mult2 = Tmp2;
					Mult1 = Tmp1;
					continue;
				} /* end if K is even */

				M1 = Mult1 % K;
				if (M1 < 0) {
					M1 += K;
				}
				M2 = Mult2 % K;
				if (M2 < 0) {
					M2 = M2 + K;
				}
				M3 = Mult3 % K;
				if (M3 < 0) {
					M3 += K;
				}
				M4 = Mult4 % K;
				if (M4 < 0) {
					M4 += K;
				}
				Tmp = (K+1)/2;  // Tmp <- (K+1) / 2
				if (M1 >= Tmp) { // If M1 >= K / 2 ... 
					M1 -= K;     // M1 = M1 - K;    
				}
				if (M2 >= Tmp) { // If M2 >= K / 2 ... 
					M2 -= K;     // M2 = M2 - K;         
				}

				if (M3 >= Tmp) {  // If M3 >= K / 2 ... 
					M3 -= K;      // M3 = M3 - K;        
				}
				if (M4 >= Tmp) { // If M4 >= K / 2 ... 
					M4 -= K;     //M4 = M4 - K;        
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
}

/* show that the number is the sum of 4 or fewer squares. See
https://www.alpertron.com.ar/4SQUARES.HTM */
/* uses the identity:
(a^2+b^2+c^2+d^2)*(A^2+B^2+C^2+D^2) = (aA+bB+cC+dD)^2 + (aB-bA+cD-dC)^2
                                    + (aC-bD-cA+dB)^2 + (aD-dA+bC-cB)^2 
This allows us to find the sum of squares for each factor separately then combine them */
static void ComputeFourSquaresNew(std::vector <zFactors> &factorlist, 
	Znum quads[4]) {
	Znum Mult1, Mult2, Mult3, Mult4, Tmp1, Tmp2, Tmp3;
	Znum p;

	quads[0] = 1;      // 1 = 1^2 + 0^2 + 0^2 + 0^2
	quads[1] = 0;
	quads[2] = 0;
	quads[3] = 0;

	if (factorlist.size() == 1) {/* only 1 factor? */
		if (factorlist[0].Factor == 1) {   // Number to factor is 1.
			return;
		}
		if (factorlist[0].Factor == 0) {    // Number to factor is 0.
			quads[0] = 0;      // 0 = 0^2 + 0^2 + 0^2 + 0^2
			return;
		}
	}

	for (int FactorIx = 0; FactorIx < factorlist.size(); FactorIx++) {
		if (factorlist[FactorIx].exponent % 2 == 0) {
			continue; /* if Prime factor appears an even number of times, no need to
					  process it in this for loop */
		}
		
		p = factorlist[FactorIx].Factor;

		/* compute 4 or less values the squares of which add up to prime p,
		return values in Mult1, Mult2, Mult3 and Mult4 */
		ComputeFourSquaresNew(p, Mult1, Mult2, Mult3, Mult4);
		assert(p == Mult1*Mult1 + Mult2*Mult2 + Mult3*Mult3 + Mult4*Mult4);

		/* use the identity:
		(a^2+b^2+c^2+d^2)*(A^2+B^2+C^2+D^2) = (aA+bB+cC+dD)^2 + (aB-bA+cD-dC)^2
		+ (aC-bD-cA+dB)^2 + (aD-dA+bC-cB)^2 */

		Tmp1 = Mult1*quads[0] + Mult2*quads[1] + Mult3*quads[2] + Mult4*quads[3];
		Tmp2 = Mult1*quads[1] - Mult2*quads[0] + Mult3*quads[3] - Mult4*quads[2];
		Tmp3 = Mult1*quads[2] - Mult3*quads[0] - Mult2*quads[3] + Mult4*quads[1];
		quads[3] = Mult1*quads[3] - Mult4*quads[0] + Mult2*quads[2] - Mult3*quads[1];

		quads[2] = Tmp3;
		quads[1] = Tmp2;
		quads[0] = Tmp1;
	} /* end for indexPrimes */

	 /* for factors that are perfect squares, multiply quads[0]-3 by sqrt(factor) */
	for (int FactorIx = 0; FactorIx < factorlist.size(); FactorIx++) {
		if (factorlist[FactorIx].Factor >= 2) {
			mpz_pow_ui(ZT(Tmp1), ZT(factorlist[FactorIx].Factor), 
				factorlist[FactorIx].exponent / 2);
			quads[0] *= Tmp1;
			quads[1] *= Tmp1;
			quads[2] *= Tmp1;
			quads[3] *= Tmp1;
		}
	}

	quads[0] = abs(quads[0]);  // ensure results are +ve
	quads[1] = abs(quads[1]);
	quads[2] = abs(quads[2]);
	quads[3] = abs(quads[3]);

	/* Sort squares: largest in quads[0], smallest in quads[3]. This
	is equivalent to a 'bubble sort' with the loops unrolled. There are
	only 6 comparisons & exchanges for 4 items. */

	if (quads[0] < quads[1]) {		// quads[0] < quads[1], so exchange them.
		Tmp1 = quads[0];
		quads[0] = quads[1];
		quads[1] = Tmp1;
	}

	if (quads[0] < quads[2]) {	// quads[0] < quads[2], so exchange them.
		Tmp1 = quads[0];
		quads[0] = quads[2];
		quads[2] = Tmp1;
	}

	if (quads[0] < quads[3]) {	// quads[0] < quads[3], so exchange them.
		Tmp1 = quads[0];
		quads[0] = quads[3];
		quads[3] = Tmp1;
	}

	if (quads[1] < quads[2]) {	// quads[1] < quads[2], so exchange them.
		Tmp1 = quads[1];
		quads[1] = quads[2];
		quads[2] = Tmp1;
	}

	if (quads[1] < quads[3]) {	// quads[1] < quads[3], so exchange them.
		Tmp1 = quads[1];
		quads[1] = quads[3];
		quads[3] = Tmp1;
	}

	if (quads[2] < quads[3]) {	// quads[2] < quads[3], so exchange them.
		Tmp1 = quads[2];
		quads[2] = quads[3];
		quads[3] = Tmp1;
	}
	return;
}

/********************************************************************************
code below is to interface between DA's code
and the new code that uses Znums, which are really mpz_t integers from MPIR or GMP
multiprecision library, with a c++ class wrapped around them that allows them
to be used pretty much like normal integers. 
**********************************************************************************/


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

/* convert number from Znum to BigInteger and factorise it. Returns false
if unable to factorise number */
bool factorise(const Znum numberZ, std::vector <zFactors> &vfactors,
	 Znum quads[]) {

	try {
		if (numberZ == 0)
			return false;  // function factor can't factorize zero

			vfactors.resize(1);
			auto rv = factor(numberZ, vfactors);
			if (quads != nullptr) {
				ComputeFourSquaresNew(vfactors, quads);
			}
			return true;
		
	}

	/* code below catches C++ 'throw' type exceptions */
	catch (const std::exception& e) {
		fprintf_s(stderr, "\n*** a standard exception was caught, with message\n '%s'\n", e.what());

		return false;
	}
}
