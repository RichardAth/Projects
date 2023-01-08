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

struct factorsS {
	int factorcount;           // number of unique prime factors
	__int64 factorlist[19][2]; // prime factor, exponent
};

struct sFactors
{
	int *ptrFactor;
	int multiplicity;
	int upperBound;
	int Type;
};

#define ZT(a) a.backend().data()  /* access mpz_t within a Znum (Boost mpz_int)*/

bool isEven(const Znum& a); /* true iff a is even (works for -ve a as well) */
//#define ZisEven(a) (mpz_even_p(ZT(a)) != 0)  

/*  get approximate size (1 limb = 64 bits) */
#define numLimbs(a) abs(ZT(a)->_mp_size)

/* counters showing how the factors were found */
struct counters {
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

};

class zFactors {
public:
	Znum Factor;
	int exponent;
	int upperBound=0;    /* used during trial division to show how far we've got. 
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
	int pari = 0;          // found by Parilib
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

	friend bool ecm(const Znum &zN, fList &Factors, Znum &Zfactor);
	friend std::vector <Znum> ModSqrt(const Znum &aa, const Znum &m);
	friend size_t DivisorList(const Znum &tnum, std::vector <Znum> &divlist);
	friend Znum primRoot(const Znum &num);
	friend int classify(Znum n);
	friend void factor(const BigInteger* pValN, int factorsMod[], sFactors astFactorsMod[]);

	/* methods that are in the class */

/* Find Euler's Totient as the product of p^(e-1)*(p-1) where p=prime and e=exponent. 
  Euler's totient function counts the positive integers up to a given integer n that 
  are relatively prime to n */
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

/* find  Carmichael function λ(n) See Wikpedia
https://en.wikipedia.org/wiki/Carmichael_function#Computing_%CE%BB(n)_with_Carmichael's_theorem
AKA reduced totient function. 
See also https://oeis.org/A002322 */
	Znum carmichael() const {
		// this only works if factorisation is complete!
		if (this->f.empty())
			return 0;  /* unable to calculate */
		Znum result = 1, term;
		for (auto i : this->f) {
			if (i.Factor == 2 && i.exponent > 2)
				mpz_pow_ui(ZT(term), ZT(i.Factor), i.exponent - 2);  // p^(e-2
			else {
				mpz_pow_ui(ZT(term), ZT(i.Factor), i.exponent - 1);  // p^(e-1)
				term = term * (i.Factor - 1);	                     // (p^(e-1)-1)*(p-1)
			}
			mpz_lcm(ZT(result), ZT(result), ZT(term));
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
		std::cout << (lang? "encontrado por" : "found by");
		if (this->tdiv > 0)
			std::cout << (lang? " la división de prueba: " : " trial division: ") << this->tdiv;
		if (this->prho > 0)
			std::cout << " Pollard-rho: " << this->prho;
		if (this->pm1 > 0)
			std::cout << " power +/- 1: " << this->pm1;
		if (this->ecm > 0)
			std::cout << (lang? "curvas elípticas: " : " elliptic curve: ") << this->ecm;
		if (this->siqs > 0)
			std::cout << " SIQS: " << this->siqs;
		if (this->msieve > 0)
			std::cout << " Msieve: " << this->msieve;
		if (this->yafu > 0)
			std::cout << " YAFU:   " << this->yafu;
		if (this->pari > 0)
			std::cout << " PARI:   " << this->pari;
		if (this->carm > 0)
			std::cout << " Carmichael: " << this->carm;
		if (this->leh > 0)
			std::cout << " Lehman: " << this->leh;
		if (this->power > 0)
			std::cout << " Perfect Power: " << this->power;
		std::cout << '\n';
	}


/* return true if the number can be expressed as the sum of 2 squares. 
Either square can be zero. */
	bool twosq() const {
		// this only works if factorisation is complete!
		for (auto& i : this->f) {
			if (i.Factor % 4 == 3) { /* p = 4k+3? */
				if (i.exponent % 2 == 1) /* exponent is odd?*/
					return false;
			}
			else { 		/* p = 4k + 1, or p=2 */
				continue;
			}
		}
		return true;
	}

/* return true if the number can be expressed as the sum of a square plus 
twice a square. Either square can be zero. 
This occurs when all odd prime factors of the form 8i+5 or 8i+7
have even exponents*/
	bool sqplustwosq() const {
		// this only works if factorisation is complete!
		for (auto& i : this->f) {
			if (i.Factor % 8 == 5 || i.Factor % 8 == 7) { 
				if (i.exponent % 2 == 1) /* exponent is odd?*/
					return false;
				else continue; /* exponent is even */
			}
		}

		return true;
	}

/* calculate the number of ways an integer n can be expressed as the sum of 2
squares x^2 and y^2. The order of the squares and the sign of x and y is significant
see http://mathworld.wolfram.com/SumofSquaresFunction.html,
also http://oeis.org/A004018 */
	Znum R2() const {
		// this only works if factorisation is complete!
		Znum b = 1;
		for (auto &i : this->f) {
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

/* The number of representations of n as the sum of two squares ignoring order 
and signs. see http://oeis.org/A000161 
e.g. 325 = 18²+1 = 17²+6² = 10²+15² so R2P(325) = 3
	  25 = 5²+0  = 3²+4² so R2P(25) = 2 
	  8 = 2² + 2² so R2P(8) =1 */
	Znum R2p() const {
		// this only works if factorisation is complete!
		Znum b = 1;
		int a0 = 0;  /* exponent of prime factor 2 */
		for (auto &i : this->f) {
			if (i.Factor < 2)
				continue; // ignore factor 1
			if (i.Factor == 2) {
				a0 = i.exponent;  /* exponent of prime factor 2 */
				continue;
			}
			if (i.Factor % 4 == 3) { /* p = 4k+3? */
				if (i.exponent % 2 == 1) /* exponent is odd?*/
					return 0;
			}
			else { 		/* p = 4k + 1 */
				b *= i.exponent + 1;
			}
		}
		if (isEven(b))
			return (b / 2);
		//else return (b + 1) / 2;
		else {
		/* if b is odd the the exponents of ALL prime factors of n other than 2 are even,
		therefore n is a perfect square, or 2*perfect square */
			if ((a0 & 1) == 0)
				/* mathworld.wolfram suggests using b-1 rather than b+1 here. In effect
				zero would be disallowed as 1 of the squares. */
				return (b + 1) / 2;  /* a0 is even i.e. n is a perfect square */
			else
				return (b + 1) / 2;  /* a0 is odd i.e. n is a 2*perfect square */
		}
	}

/* calculate number of divisors of n (including 1 and n itself), given its list 
   of prime factors. */
	Znum NoOfDivs() const {
		Znum result = 1;
		if (this->f.empty())
			return 0;
		/* n=1 is a special case */
		if (this->f.size() == 1 && this->f[0].Factor == 1)
			return 1;
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

	/* if n is 2, 4, an odd prime power, or twice an odd prime power it has
       primitive roots, otherwise it does not. An odd prime is any prime
       other than 2. Returns true or false, requires list of prime factors */
	bool hasPrimitiveRoot() const {
		if (this->f.size() == 1) {
			/* only 1 factor */
			if (this->f[0].Factor > 2)
				return true;  /* n is an odd prime power */
			else if (this->f[0].exponent <= 2)
					return true;  /* n = 2 or 4 */
				else return false;  /* n is higher power of 2 */
		}
		else if (this->f.size() == 2) {
			if (this->f[0].Factor == 2 && this->f[0].exponent == 1)
				return true;  /* n is twice an odd prime */
			else
				return false;  /* not twice an odd prime power */
		}
		else
			return false;  /* n is not in any category that has primitive roots */
	}

	/* copy counters */
	struct counters getCtrs() const {
		struct counters temp;
		temp.carm   = this->carm;    // Carmichael
		temp.ecm    = this->ecm;     // elliptic curve
		temp.leh    = this->leh;     // Lehman
		temp.msieve = this->msieve;  // Msieve
		temp.pm1    = this->pm1;     // power + / -1
		temp.power  = this->power;   // perfect power
		temp.prho   = this->prho;    // Pollard-rho
		temp.siqs   = this->siqs;    // SIQS
		temp.tdiv   = this->tdiv;    // trial division
		temp.yafu   = this->yafu;    // YAFU
		return temp;
	}

};

//void squareFree(Znum& num, Znum& sq, std::vector<zFactors>& sqf);

std::vector <Znum> primeModSqrt(const Znum& aa, const Znum& p);
std::vector <Znum> ModSqrtQE(const Znum& aa, const Znum& m); /* modular square root*/

extern int ElipCurvNo;            // Elliptic Curve Number

void showECMStatus(void);
constexpr unsigned long long int gcd(unsigned long long int u, unsigned long long int v);
long long int PollardRho(long long int n, int depth = 0);

bool factorise(Znum numberZ, fList &vfactors, Znum quads[]);
bool callMsieve(const Znum& num, fList& Factors);
bool callYafu(const Znum& num, fList& Factors);
void parifactor(const Znum& n, fList& factors);

void LehmanZ(const Znum &nbr, int k, Znum &factor);

/* return a factor of N, using Shanks's square forms factorization method. */
uint64_t SQUFOF(const uint64_t N);

int mpz_bpsw_prp(const mpz_t n); /* Baillie-Pomerance-Selfridge-Wagstaff probablistic primality test*/
int mpz_aprtcle(const mpz_t N, const int verbose);  /* APR-CL prime testing */

// returns 2^exp. exp must be less than 64
constexpr unsigned __int64 pow2(unsigned int exp) {
	assert(exp < 64);
	return 1ULL << exp;  // exp must be less than 64
}


/* error and return codes, errors are -ve, OK is 0, FAIL is +1 */
enum class retCode
{
	NUMBER_TOO_LOW = -100,
	NUMBER_TOO_HIGH,
	INTERM_TOO_HIGH,
	DIVIDE_BY_ZERO,
	PAREN_MISMATCH,
	SYNTAX_ERROR,
	TOO_MANY_PAREN,
	INVALID_PARAM,
	ARGUMENTS_NOT_RELATIVELY_PRIME,
	//EXPR_BREAK,
	//EXPR_OUT_OF_MEMORY,
	//EXPR_CANNOT_USE_X_IN_EXPONENT,
	//EXPR_DEGREE_TOO_HIGH,
	EXPONENT_TOO_LARGE,
	EXPONENT_NEGATIVE,
	//EXPR_LEADING_COFF_MULTIPLE_OF_PRIME,
	//EXPR_CANNOT_LIFT,
	//EXPR_MODULUS_MUST_BE_GREATER_THAN_ONE,
	//EXPR_MODULUS_MUST_BE_PRIME_EXP,
	//EXPR_BASE_MUST_BE_POSITIVE,
	//EXPR_POWER_MUST_BE_POSITIVE,
	EXPR_MODULUS_MUST_BE_NONNEGATIVE,
	//EXPR_VAR_OR_COUNTER_REQUIRED,
	EXPR_OK = 0,
	EXPR_FAIL = 1
};


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
/* similar to throwExc but is compatible with constexpr functions */
#define ThrowExs(a)                              \
{												 \
	char mesg[180] = { 0 };					     \
	sprintf_s(mesg, sizeof(mesg), "%s %s line  %d in file %s ", a,  __func__, __LINE__, __FILE__); \
	throw std::range_error(mesg);                \
}

