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

#include "boost/multiprecision/gmp.hpp" 
#define ZT(a) a.backend().data()
#define ZisEven(a) (mpz_even_p(ZT(a)) != 0)  /* true iff a is even (works for -ve a as well) */
typedef boost::multiprecision::mpz_int Znum;

/* ComputeNumDigits(n,r): Number of digits of n in base r. */
long long ComputeNumDigits(const Znum &n, const Znum &radix);
void ShowLargeNumber(const Znum &Bi_Nbr, int digitsInGroup, bool size, bool hex);
class zFactors {
public:
	Znum Factor;
	int exponent;
	int upperBound;    /* used during trial division to show how far we've got. 
						-1 indicats that Factor is prime */
	/* define all comparison operators */
	bool operator == (const zFactors &b) const {
		return this->Factor == b.Factor;
	}
	bool operator != (const zFactors &b) const {
		return this->Factor != b.Factor;
	}
	bool operator < (const zFactors &b) const {
		return this->Factor < b.Factor;
	}
	bool operator <= (const zFactors &b) const {
		return this->Factor <= b.Factor;
	}
	bool operator > (const zFactors &b) const {
		return this->Factor > b.Factor;
	}
	bool operator >= (const zFactors &b) const {
		return this->Factor >= b.Factor;
	}
};

class fList {
private:
	std::vector <zFactors> f;

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
public:
	friend void insertIntFactor(fList &Factors, int pi, long long div, ptrdiff_t ix);
	friend void insertBigFactor(fList &Factors, const Znum &divisor);
	friend void SortFactors(fList &Factors);
	friend void TrialDiv(fList &Factors, long long LehmanLimit);
	friend bool factor(const Znum &toFactor, fList &Factors);
	friend void ComputeFourSquares(const fList &factorlist, Znum quads[4], Znum num);

	friend bool ecm(Znum &zN, fList &Factors, Znum &Zfactor);
	friend std::vector <Znum> ModSqrt(const Znum &aa, const Znum &m);

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

	/* print factors */
	void Xprint() const {
		for (auto i : this->f) {
			std::cout << i.Factor << "^" << i.exponent << " ("
				<< i.upperBound << ")  * ";
		}
		std::cout << '\n';
	}

	/* print factors */
	void print(bool neg) const {
		if (!this->isPrime()) {
			/* print factor list */
			std::cout << " = ";
			if (neg)
				std::cout << "-";
			for (size_t i = 0; i < this->fsize(); i++) {
				if (i > 0)
					std::cout << " * ";
				ShowLargeNumber(this->f[i].Factor, 6, false, false);
				if (this->f[i].exponent > 1)
					std::cout << "^" << this->f[i].exponent;
			}
		}
		else
			std::cout << " is prime";  //number has only 1 factor
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

	/* remove any square factors. Return removed factors in sqf*/
	void sqfree(std::vector<zFactors> &sqf) const {
		zFactors temp;
		sqf.clear();
		for (auto i : this->f) {
			if (i.exponent >= 2) {
				/* copy exponent value, rounded down to multiple of 2*/
				temp.exponent = i.exponent - (i.exponent & 1);
				temp.Factor = i.Factor;
				sqf.push_back(temp);
				i.exponent -= temp.exponent;
			}
		}
	}

	/* checks whether the original number was prime */
	bool isPrime() const {
		// this only works if factorisation is complete!
		return (this->f.size() == 1 && this->f[0].exponent == 1);
	}

	/* return number of unique factors */
	size_t fsize() const {
		return this->f.size();
	}

	/* get number of digits in 2nd largest factor */
	int sndFac() const {
		if (this->f.size() > 1)
			return (int)ComputeNumDigits((this->f.end() - 2)->Factor, 10);
		else return 0;  /* return 0 if only 1 factor */
	}

	/* recheck value of product of factors, and also get total number of factors */
	Znum recheck(int &totalfacs) const {
		Znum result = 1;
		totalfacs = 0;
		for (auto i : this->f) {
			totalfacs += i.exponent;
			for (int j = 1; j <= i.exponent; j++)
				result *= i.Factor;
		}
		return result;
	}

	/* return value of smallest factor */
	Znum fmin ()const{
		if (!this->f.empty())
			return (this->f.front().Factor);
		else return -1;  /* error - factor list is empty */
	}

	/* return value of largest factor */
	Znum fmax()const {
		if (!this->f.empty())
			return this->f.back().Factor ;
		else return -1;  /* error - factor list is empty */
	}

};


extern int lang;    // 0 English, 1 = Spanish
extern std::string YafuPath;
extern std::string yafuprog;
extern std::string MsievePath;
extern std::string MsieveProg;
extern bool breakSignal;
extern std::vector <Znum> roots;   /* used by functions that return multiple values */
/* access underlying mpz_t inside an bigint */
#define ZT(a) a.backend().data()

extern int ElipCurvNo;            // Elliptic Curve Number
extern unsigned long long *primeList;
extern unsigned int prime_list_count;
extern unsigned long long int primeListMax;
extern bool msieve;
extern bool yafu;
extern int verbose;

void showECMStatus(void);
constexpr unsigned long long int gcd(unsigned long long int u, unsigned long long int v);
long long int PollardRho(long long int n, int depth = 0);
long long MulPrToLong(const Znum &x);
bool factorise(Znum numberZ, fList &vfactors, Znum quads[]);
void delfile(const std::string &path, const char * FileName);
void writeIni(void);
bool callMsieve(const Znum &num, fList&Factors);
bool callYafu(const Znum &num, fList&Factors);
void generatePrimes(unsigned long long int max_val);
void LehmanZ(const Znum &nbr, int k, Znum &factor);

int mpz_bpsw_prp(const mpz_t n); /* Baillie-Pomerance-Selfridge-Wagstaff probablistic primality test*/
int mpz_aprtcle(const mpz_t N, const int verbose);  /* APR-CL prime testing */

// returns 2^exp. exp must be less than 64
constexpr unsigned __int64 pow2(unsigned int exp) {
	//assert(exp < 64);
	return 1ULL << exp;  // exp must be less than 64
}
unsigned __int64 R3(__int64 n);
std::vector <Znum> primeModSqrt(const Znum &aa, const Znum &p);
unsigned __int64 modMult(unsigned __int64 a, unsigned __int64 b, unsigned __int64 mod);
Znum             modMult(const Znum &a, const Znum &b, Znum mod);
// calculate a^n%mod using 'bigints'   
unsigned __int64 modPower(unsigned __int64 a, unsigned __int64 n,
	unsigned __int64 mod);
Znum             modPower(const Znum &a, const Znum &n, const Znum &mod);
unsigned __int64 modPowerBi(const Znum &a, const Znum &n, unsigned __int64 mod);
constexpr __int64 power(const __int64 x, unsigned int n);
Znum              power(const Znum &x, unsigned long long n);
int jacobi(__int64 k, unsigned __int64 n);
int jacobi(const Znum &k, const Znum &n);


/* error and return codes, errors are -ve, OK is 0, FAIL is +1 */
enum class retCode
{
	//EXPR_NUMBER_TOO_LOW,
	EXPR_NUMBER_TOO_HIGH = -100,
	EXPR_INTERM_TOO_HIGH,
	EXPR_DIVIDE_BY_ZERO,
	EXPR_PAREN_MISMATCH,
	EXPR_SYNTAX_ERROR,
	EXPR_TOO_MANY_PAREN,
	EXPR_INVALID_PARAM,
	EXPR_ARGUMENTS_NOT_RELATIVELY_PRIME,
	//EXPR_BREAK,
	//EXPR_OUT_OF_MEMORY,
	//EXPR_CANNOT_USE_X_IN_EXPONENT,
	//EXPR_DEGREE_TOO_HIGH,
	EXPR_EXPONENT_TOO_LARGE,
	EXPR_EXPONENT_NEGATIVE,
	//EXPR_LEADING_COFF_MULTIPLE_OF_PRIME,
	//EXPR_CANNOT_LIFT,
	//EXPR_MODULUS_MUST_BE_GREATER_THAN_ONE,
	//EXPR_MODULUS_MUST_BE_PRIME_EXP,
	//EXPR_BASE_MUST_BE_POSITIVE,
	//EXPR_POWER_MUST_BE_POSITIVE,
	//EXPR_MODULUS_MUST_BE_NONNEGATIVE,
	//EXPR_VAR_OR_COUNTER_REQUIRED,
	EXPR_OK = 0,
	EXPR_FAIL = 1
};

struct factorsS {
	int factorcount;           // number of unique prime factors
	__int64 factorlist[19][2]; // prime factor, exponent
};

unsigned int primeFactors(unsigned __int64 tnum, factorsS &f);

/* throw exception. For a list of exception classes derived from std:: exception
see https://en.cppreference.com/w/cpp/error/exception
range_error can be used to report range errors (that is, situations where a
result of a computation cannot be represented by the destination type)
a is a text string e.g. "integer overflow". The message generated includes
the text string a, function name, line number and source file name */
#define ThrowExc(a)                              \
{                                                \
	std::string line = std::to_string(__LINE__); \
	std::string mesg = a;                        \
    mesg += " " ;                                \
	mesg += __func__;                            \
	mesg += " line ";  mesg += line;             \
	mesg += " in file "; mesg += __FILE__;       \
	throw std::range_error(mesg);                \
}

unsigned long long llSqrt(const unsigned long long n);
bool isPerfectSquare(const Znum &num);