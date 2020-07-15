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

class fList {
private:
	std::vector <zFactors> f;
public:
	/* counters showing how the factors were found */
	int pm1 = 0;           // power + / -1
	int tdiv = 0;          // trial division
	int prho = 0;          // Pollard-rho
	int ecm = 0;           // elliptic curve
	int siqs = 0;          // SIQS
	int msieve = 0;        // Msieve
	int yafu = 0;          // YAFU
	int carm = 0;          // Carmichael
	int leh = 0;           // Lehman:
	int power = 0;		   // perfect power

	friend void insertIntFactor(fList &Factors, int pi, long long div, ptrdiff_t ix);
	friend void insertBigFactor(fList &Factors, const Znum &divisor);
	friend void SortFactors(fList &Factors);
	friend void TrialDiv(fList &Factors, long long LehmanLimit);
	friend bool factor(const Znum &toFactor, fList &Factors);
	friend bool factorise(Znum numberZ, fList &vfactors, Znum quads[]);
	friend void squareFree(Znum &num, Znum &sq, std::vector<zFactors> &sqf);
	friend void doFactors(const Znum &Result, bool test);
	friend bool factortest(const Znum &x3, const int method = 0);

	/* methods that are in the class */

// Find Euler's Totient as the product of p^(e-1)*(p-1) where p=prime and e=exponent.
	Znum totient() const {
    // this only works if factorisation is complete!
		Znum result = 1, term;
		if (this->f.empty())
			return 0;
		for (auto i : this->f) {
			mpz_pow_ui(ZT(term), ZT(i.Factor), i.exponent - 1);  // p^(e-1)
			term = term * (i.Factor - 1);	                        // (p^(e-1)-1)*(p-1)
			result *= term;
		}
		return result;
	}

/* For any positive integer n, define μ(n) as the sum of the primitive nth roots of unity.
It has values in {−1, 0, 1} depending on the factorization of n into prime factors:
μ(n) = 1 if n is a square-free positive integer with an even number of prime factors.
μ(n) = −1 if n is a square-free positive integer with an odd number of prime factors.
μ(n) = 0 if n has a squared prime factor. */
	int mob() const {
		// this only works if factorisation is complete!
		for (auto ex : this->f) {
			if (ex.exponent > 1)
				return 0;		// n is not square-free
		}
		auto result = this->f.size();
		if ((result & 1) == 1)
			return -1;  // odd nuber of  prime factors
		else
			return 1;    // even number of prime factors 
	}
	/* initialise factor list */
	void set(const Znum &toFactor) {
		this->f.resize(1);          // change size of factor list to 1
		this->f[0].exponent = 1;
		this->f[0].Factor = toFactor; 
		this->f[0].upperBound = 0;  // assume toFactor's not prime
	}
	void print() const {
		for (auto i : this->f) {
			std::cout << i.Factor << "^" << i.exponent << " ("
				<< i.upperBound << ")  * ";
		}
		std::cout << '\n';
	}

/* indicate how the number's factors were found. No detailed breakdown
for factors found by YAFU or Msieve */
	void prCounts() const {
		if (this->f.size() == 1 && this->f[0].exponent == 1)
			return;  // number is prime or has not been factored
		std::cout << "found by";
		if (this->tdiv > 0)
			std::cout << " trial division: " << this->tdiv;
		if (this->prho > 0)
			std::cout << " Pollard-rho: " << this->prho;
		if (this->pm1 > 0)
			std::cout << " power +/- 1: " << this->pm1;
		if (this->ecm > 0)
			std::cout << " elliptic curve: " << this->ecm;
		if (this->siqs > 0)
			std::cout << " SIQS: " << this->siqs;
		if (this->msieve > 0)
			std::cout << " Msieve: " << this->msieve;
		if (this->yafu > 0)
			std::cout << " YAFU:   " << this->yafu;
		if (this->carm > 0)
			std::cout << " Carmichael: " << this->carm;
		if (this->leh > 0)
			std::cout << " Lehman: " << this->leh;
		if (this->power > 0)
			std::cout << " Perfect Power: " << this->power;
		std::cout << '\n';
	}

/* calculate the number of ways an integer n can be expressed as the sum of 2
squares x^2 and y^2. The order of the squares and the sign of x and y is significant
see http://mathworld.wolfram.com/SumofSquaresFunction.html,
also http://oeis.org/A004018 */
	Znum R2() const {
		// this only works if factorisation is complete!
		Znum b = 1;
		for (auto i : this->f) {
			if (i.Factor <= 2)
				continue; // ignore factor 1 or 2
			if (i.Factor % 4 == 3) { /* p = 4k+3? */
				if (i.exponent % 2 == 1) /* exponent is odd?*/
					return 0;
			}
			else { 		/* p = 4k + 1 */
				b *= i.exponent + 1;
			}
		}
		return b * 4;
	}

/* calculate number of divisors of n, given its list of prime factors.
NB for n=1 function returns 2, but correct value would be 1.
From just the exponents we can't distinguish n=1 from n=prime */
	Znum NoOfDivs() const {
		Znum result = 1;
		if (this->f.empty())
			return 0;
		for (auto i : this->f) {
			result *= i.exponent + 1;
		}
		return result;
	}

// sum of divisors is the product of (p^(e+1)-1)/(p-1) where p=prime factor 
// and e=exponent.
	Znum DivisorSum() const {
		Znum result = 1, term;
		if (this->f.empty())
			return 0;
		if (this->f[0].Factor == 1)
			return 1;		// special case: if original number is 1
		for (auto i: this->f) {
			mpz_pow_ui(ZT(term), ZT(i.Factor), i.exponent + 1);  // p^(e+1)
			term = (term - 1) / (i.Factor - 1);	                 // (p^(e+1)-1)/(p-1)
			result *= term;
		}
		return result;
	}


/* Concatenates the prime factors (base 10) of num according to the mode
Order of factors: Ascending or 	Descending
Repeated factors: No or Yes
*/
	Znum ConcatFact(bool descending, bool repeat) const {
		std::string result;
		Znum rvalue = 0;
		char *buffer = NULL;

		if (this->f.empty())
			return 0;
		if (descending)   /* start with largest factor */
			//for (ptrdiff_t i = this->f.size()-1; i >=0; i--) {
			for (auto i = this->f.rbegin(); i != this->f.rend(); i++) {
				buffer = mpz_get_str(NULL, 10, ZT((*i).Factor));
				if (!repeat)
					result += buffer; // concatenate factor
				else
					for (int j = 1; j <= (*i).exponent; j++)
						result += buffer; // concatenate factor
				free(buffer);
				buffer = NULL;
			}
		else  /* start with smallest factor */
			for (auto i: this->f) {
				buffer = mpz_get_str(NULL, 10, ZT(i.Factor));
				if (!repeat)
					result += buffer;  // concatenate factor
				else
					for (int j = 1; j <= i.exponent; j++)
						result += buffer;  // concatenate factor
				free(buffer);
				buffer = NULL;
			}
		mpz_set_str(ZT(rvalue), result.data(), 10); /* convert back from a string to a number */
		return rvalue;
	}
};


void showECMStatus(void);
bool ecm(Znum &Nz, long long maxdivisor, fList &factors);
extern int lang;
extern Znum Zfactor, Zfactor2;

/* access underlying mpz_t inside an bigint */
#define ZT(a) a.backend().data()

long long MulPrToLong(const Znum &x);

unsigned long long int gcd(unsigned long long int u, unsigned long long int v);
long long int PollardRho(long long int n);
extern int ElipCurvNo;            // Elliptic Curve Number

extern unsigned long long *primeList;
extern unsigned int prime_list_count;
extern unsigned long long int primeListMax;
extern bool msieve;
extern bool yafu;
extern int verbose;
bool callMsieve(const Znum &num, fList&Factors);
bool callYafu(const Znum &num, fList&Factors);
void generatePrimes(unsigned long long int max_val);
void LehmanZ(const Znum &nbr, int k, Znum &factor);