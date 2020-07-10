#pragma once

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

#ifdef __GNUC__
#include "gmp.h"
#else
#include "mpir.h"
#endif
#include "boost/multiprecision/gmp.hpp" 
#define ZT(a) a.backend().data()
#define ZisEven(a) (mpz_even_p(ZT(a)) != 0)  /* true iff a is even (works for -ve a as well) */
typedef boost::multiprecision::mpz_int Znum;


class zFactors {
public:
	Znum Factor;
	int exponent;
	int upperBound;
};


void showECMStatus(void);
bool ecm(Znum &Nz, long long maxdivisor);
extern int lang;
extern Znum Zfactor, Zfactor2;

/* access underlying mpz_t inside an bigint */
#define ZT(a) a.backend().data()

bool factorise(const Znum numberZ, std::vector <zFactors> &factorlist,
	Znum Quad[]);

long long MulPrToLong(const Znum &x);


unsigned long long int gcd(unsigned long long int u, unsigned long long int v);
long long int PollardRho(long long int n);
extern int ElipCurvNo;            // Elliptic Curve Number
struct ctrs{
	int pm1;           // power + / -1
	int tdiv;          // trial division
	int prho;          // Pollard-rho
	int ecm;           // elliptic curve
	int siqs;          // SIQS
	int msieve;        // Msieve
	int yafu;          // YAFU
	int carm;          // Carmichael
	int leh;           // Lehman:
	int power;		   // perfect power
};
extern struct ctrs counters;

extern unsigned long long *primeList;
extern unsigned int prime_list_count;
extern unsigned long long int primeListMax;
extern bool msieve;
extern bool yafu;
extern int verbose;
bool callMsieve(const Znum &num, std::vector<zFactors>&Factors);
bool callYafu(const Znum &num, std::vector<zFactors>&Factors);
void insertBigFactor(std::vector<zFactors> &Factors, Znum &divisor);
void generatePrimes(unsigned long long int max_val);
void LehmanZ(const Znum &nbr, int k, Znum &factor);