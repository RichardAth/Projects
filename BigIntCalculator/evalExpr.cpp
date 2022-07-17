/*
This file is part of Alpertron Calculators.
Copyright 2015 Dario Alejandro Alpern
Alpertron Calculators is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
Alpertron Calculators is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with Alpertron Calculators.If not, see < http://www.gnu.org/licenses/>.
*/

#include "pch.h"
#include "bignbr.h"   /* access PrimalityTest function*/
#define WIN32_LEAN_AND_MEAN
#include <Windows.h>
#include "factor.h"
#include <stack>


// declarations for external function

void biperm(int n, Znum &result);  
Znum llt(const Znum &p);
size_t inverseTotient(__int64 n, std::vector<unsigned __int64> **result, bool debug,
	int level, bool dump);

/* forward references */
int new_uvar(const char *name);
int set_uvar(const char *name, const Znum &data);
int get_uvar(const char *name, Znum data);
static void free_uvars();

const unsigned long long max_prime = 1000000007;  // arbitrary limit 10^9,
std::vector <Znum> roots;   /* used by functions that return multiple values */
bool multiValue = false;

typedef struct        // used for user variables
{
	char name[40];    // variable name
	Znum data;        // and value
} uvar_t;

struct
{
	std::vector <uvar_t> vars;
	int num = 0;              /* number of user variables*/
	int alloc = 0;            /* space allocated for user variables */
} uvars;

//user variables
//uvars_t uvars ;

/* list of operators, arranged in order of priority, order is not exactly the
same as C or Python. Followed by list of function codes*/
enum class opCode {
	fact        = 21,	// !   factorial 
	prim        = 23,	// #   primorial
	unary_minus = 1,    // C and Python put unary minus above multiply, divide & modulus
	not         = 16,   // C and Python put bitwise not with unary minus
	power       = 0,
	multiply    = 2,
	divide      = 3,
	remainder   = 4,   // AKA modulus
	comb        = 5,   // nCk, also known as binomial coefficient
	plus        = 6,
	minus       = 7,
	shr         = 8,   // right shift
	shl         = 9,   // left shift
	not_greater = 10,
	not_less    = 11,
	greater     = 12,
	less        = 13,
	not_equal   = 14,
	equal       = 15,
	and         = 17,      // C and Python put AND before XOR before OR
	xor         = 18,
	or          = 19,
	leftb       = 20,
	assign      = 24,       /* assignment operator */
	rightb      = 25,       // right bracket 
/* functions */
	fn_gcd      = 100,
	fn_lcm,
	fn_modpow,
	fn_modinv,
	fn_totient,
	fn_numdivs,
	fn_sumdivs,
	fn_sumdigits,
	fn_numdigits,
	fn_revdigits,
	fn_isprime,
	fn_concatfact,
	fn_fib,
	fn_luc,
	fn_primePi,
	fn_part,
	fn_np,                /* next prime*/
	fn_pp,                /* previous prime */
	fn_r2,
	fn_r2p,
	fn_r3,
	fn_r3h,
	fn_hurwitz,
	fn_legendre,
	fn_jacobi,
	fn_kronecker,
	fn_llt,
	fn_sqrt,              /* square root */
	fn_nroot,             /* nth root */
	fn_bpsw,              /* primality test */
	fn_aprcl,             /* primality test */
	fn_numfact,           /* number of factors */
	fn_minfact,           /* smallest factor */
	fn_maxfact,           /* largest factor */
	fn_ispow,             /* check if perfect power */
	fn_modsqrt,           /* modular square root */
	fn_invtot,            /* inverse totient */
	fn_divisors,          /* list of divisors */
	fn_primroot,          /* lowest primitive root */
	fn_invalid = -1,
};


struct  functions {
	char fname[11];        // maximum name length 10 chars (allow for null terminator)
	int  NoOfParams;       // number of parameters 
	opCode  fCode;        // integer code for function
};

/* list of function names. No function name can begin with C, SHL, SHR, NOT, 
 AND, OR, XOR because this would conflict with the operator. 
 Longer names must come before short ones that start with the same letters to 
 avoid mismatches */
const static struct functions functionList[]{
	"MODPOW",    3,  opCode::fn_modpow,  		// name, number of parameters, code
	"MODINV",    2,  opCode::fn_modinv,
	"TOTIENT",   1,  opCode::fn_totient,
	"SUMDIVS",   1,  opCode::fn_sumdivs,
	"SUMDIGITS", 2,  opCode::fn_sumdigits,
	"SQRT",      1,  opCode::fn_sqrt,
	"NUMDIGITS", 2,  opCode::fn_numdigits,
	"NUMDIVS",   1,  opCode::fn_numdivs,
	"NROOT",     2,  opCode::fn_nroot,
	"REVDIGITS", 2,  opCode::fn_revdigits,
	"ISPRIME",   1,	 opCode::fn_isprime,
	"NUMFACT",   1,  opCode::fn_numfact,
	"MINFACT",   1,  opCode::fn_minfact,
	"MAXFACT",   1,  opCode::fn_maxfact,
	"MODSQRT",   2,  opCode::fn_modsqrt,    // find x such that x^2 = a mod p
	"DIVISORS",  1,  opCode::fn_divisors,   // list of divisors    
	"FactConcat",2,  opCode::fn_concatfact,     // FactConcat must come before F
	"InvTot",    1,  opCode::fn_invtot,         // inverse totient
	"PrimRoot",  1,  opCode::fn_primroot,       /* smallest primitive root */
	"GCD",       SHORT_MAX,  opCode::fn_gcd,    /* gcd, variable no of parameters */
	"LCM",       SHORT_MAX,  opCode::fn_lcm,    /* lcm, variable no of parameters */
	"F",         1,  opCode::fn_fib,			// fibonacci
	"LLT",	     1,  opCode::fn_llt,           // lucas-Lehmer test
	"LE",		 2,  opCode::fn_legendre,
	"L",         1,  opCode::fn_luc,			// Lucas Number
	"PI",		 1,  opCode::fn_primePi,		// prime-counting function. PI must come before P
	"P",         1,  opCode::fn_part,			// number of partitions
	"N",         1,  opCode::fn_np,			// next prime
	"BPSW",      1,  opCode::fn_bpsw,       // Baillie-Pomerance-Selfridge-Wagstaff
	"B",         1,  opCode::fn_pp,			// previous prime
	"R2P",       1,  opCode::fn_r2p,        // number of ways n can be expressed as sum of 2 squares ignoring order and signs
	"R2",		 1,  opCode::fn_r2,			// number of ways n can be expressed as sum of 2 squares
	"R3h",       1,  opCode::fn_r3h,        // number of ways n can be expressed as sum of 3 squares
											// calculated using hurwitz class number
	"R3",        1,  opCode::fn_r3,         // number of ways n can be expressed as sum of 3 squares
	"JA",		 2,  opCode::fn_jacobi,
	"KR",		 2,  opCode::fn_kronecker,
	"APRCL",     1,  opCode::fn_aprcl,      // APR-CL prime test
	"ISPOW",     1,  opCode::fn_ispow,
	"HCLASS",    1,  opCode::fn_hurwitz,
};

/* list of operators.  */
struct oper_list{
	char oper[4];       /* operator as ascii string*/
	opCode operCode;
	int pri;     /* operator precedence, 0 is highest, 99 is lowest */
	bool left;   /* associativity ; true = left to right, false = right to left
		operators with the same precedence must have the same associativity. */
	bool pre;    /* true if unary operator e.g. - precedes expression, false if
					it follows e.g. !, otherwise not used */
	int numOps;  /* number of operands; 1 = unary, or 2 normally, or 0 for bracket)*/
};
/* list of operators with corresponding codes and priority. For the search to
work correctly
** must precede *,
<< and <= must precede <
>> and >= must precede > in this list.
unary plus is NOT in this list. 
unary minus is last entry */
const static struct oper_list operators[]{
		{"C",	 opCode::comb,	      4,  true,  false, 2},  // combination nCk, also known as binomial coefficient
		{ "^",   opCode::power,       2,  false, false, 2},
		{ "**",  opCode::power,       2,  false, false, 2},   // can use ^ or ** for exponent
		{ "*",   opCode::multiply,    3,  true,  false, 2},
		{ "/",   opCode::divide,      3,  true,  false, 2},
		{ "%",   opCode::remainder,   3,  true,  false, 2},
		{ "+",   opCode::plus,        5,  true,  false, 2},
		{ "-",   opCode::minus,       5,  true,  false, 2},
		{ "SHL", opCode::shl,         6,  true,  false, 2},
		{ "<<",  opCode::shl,         6,  true,  false, 2},     // can use << or SHL for left shift
		{ "SHR", opCode::shr,         6,  true,  false, 2},
		{ ">>",  opCode::shr,         6,  true,  false, 2},     // can use SHR or >> for right shift
		{ "<=",  opCode::not_greater, 7,  true,  false, 2},
		{ ">=",  opCode::not_less,    7,  true,  false, 2},
		{ ">",   opCode::greater,     7,  true,  false, 2},	  // to avoid mismatches > and < must come after >> and <<
		{ "<",   opCode::less,        7,  true,  false, 2},
		{ "!=",  opCode::not_equal,   8,  true,  false, 2},
		{ "==",  opCode::equal,       8,  true,  false, 2},
		{ "NOT", opCode::not,         1,  false, true,  1},      // bitwise NOT
		{ "AND", opCode::and,         9,  true,  false, 2},      // bitwise AND
		{ "OR",  opCode:: or,        11,  true,  false, 2},      // bitwise OR
		{ "XOR", opCode::xor,        10,  true,  false, 2},      // bitwise exclusive or
	  //{ "!!",  opCode::dfact,       0,  true,  false, 1},      // double factorial
		{ "!",   opCode::fact,        0,  true,  false, 1},      // multi-factorial
		{ "#",   opCode::prim,        0,  true,  false, 1},      // primorial
		{ "(",   opCode::leftb,      99, true,  false, 0},      // left bracket
		{ ")",   opCode::rightb,     -1, true,  false, 0},      // right bracket
		{"=",    opCode::assign,     12, false, false, 2},      // assignment
		/* unary - must be the last entry */
        {"U-",   opCode::unary_minus, 1, false, true,  1 },     // unary -
};

enum class types { Operator, func, number, comma, error, uservar, end };

struct token {
	types typecode;
	int function;   /* contains function code index, when typecode = func 
					   contains operator index when typecode = Operator */
	opCode oper;    /* contains operator value, only when typecode = Operator or func */
	Znum value;     /* contains numeric value,  only when typecode = number */
	size_t userIx;  /* index into user variable list (only when typecode = uservar) */
	short numops;   /* number of operands/parameters */
};

static retCode tokenise(const std::string expr, std::vector <token> &tokens, int &asgCt);

// returns true if num is a perfect square.
bool isPerfectSquare(const Znum &num) {
	return (mpz_perfect_square_p(ZT(num)) != 0); /* true if num is a perfect square */
}

/* calculate Euler's totient for n as the product of p^(e-1)*(p-1)
where p=prime factor and e=exponent.*/
static Znum ComputeTotient(const Znum &n) {
	fList factorlist;

	if (n == 1)
		return 1;
	auto rv = factorise(n, factorlist, nullptr);
	if (rv) {
		auto tot = factorlist.totient();
		return tot;
	}
	else return 0;
}

/* calculate number of divisors of n, from its list of prime factors. */
static Znum ComputeNumDivs(const Znum &n) {
	fList factorlist;

	auto rv = factorise(n, factorlist, nullptr);
	if (rv) {
		auto divisors = factorlist.NoOfDivs();
		return divisors;
	}
	else return 0;
}

/* generate sorted list of all divisors of tnum
N.B. includes non-prime factors e.g. 1, 2, 4, 8 and 16 are divisors of 16.
the value returned is the number of divisors.
*/
static size_t DivisorList(const Znum &tnum, std::vector <Znum> &divlist) {

	fList f;          /* prime factors of tnum */
	Znum divisors;    /* number of divisors */
	long long numdiv;
	size_t ctr = 0, ct2, ccpy;
	long long numfactors;  /* number of unique factors */

	divlist.clear();

	auto rv = factorise(tnum, f, nullptr);

	if (!rv) {
		return 0;
	}

	divisors = f.NoOfDivs();
	numdiv = MulPrToLong(divisors);   /* number of divisors of tnum */
	numfactors = f.fsize();           /* number of unique prime factors of tnum */
	divlist.resize(numdiv);

	divlist[ctr++] = 1;  // 1st divisor

	if (tnum <= 1)
		return 1;

	for (long long i = 0; i < numfactors; i++) {
		ct2 = 0;
		ccpy = ctr;
		for (int x = 1; x <= f.f[i].exponent; x++) {
			for (size_t j = 0; j < ctr; j++)
				divlist[j + ccpy] = divlist[j + ct2] * f.f[i].Factor;
			ct2 += ctr;
			ccpy += ctr;
		}
		ctr = ccpy;
	}

	std::sort(divlist.begin(), divlist.end());
#ifdef _DEBUG
	//if (verbose > 0) {
	//	gmp_printf("divisors of %Zd are: ", tnum);
	//	for (int i = 0; i < numdiv; i++) {
	//		gmp_printf("%Zd ", divlist[i]);
	//	}
	//	putchar('\n');
	//}
#endif
	return numdiv;
}

/* sum of divisors is the product of (p^(e+1)-1)/(p-1)
 where p=prime factor and e=exponent. */
static Znum ComputeSumDivs(const Znum &n) {
	fList factorlist;

	if (n == 1)
		return 1;   // 1 only has 1 divisor. 
	auto rv = factorise(n, factorlist, nullptr);
	if (rv) {
		auto divisors = factorlist.DivisorSum();
		return divisors;
	}
	else return 0;
}

/* number of distinct prime factors of n */
static Znum ComputeNumFact(const Znum &n) {
	if (n >= -1 && n <= 1)
		return 1;
	fList factorlist;
	auto rv = factorise(n, factorlist, nullptr);
	if (rv)
		return factorlist.fsize();
	else
		return 0;
}

/* minimum prime factor of n. */
static Znum ComputeMinFact(const Znum &n) {
	if (n >= -1 && n <= 1)
		return n;
	fList factorlist;
	auto rv = factorise(n, factorlist, nullptr);
	if (rv)
		return factorlist.fmin();
	else
		return 0;
}

/* maximum prime factor of n.*/
static Znum ComputeMaxFact(const Znum &n) {
	if (n >= -1 && n <= 1)
		return n;
	fList factorlist;
	auto rv = factorise(n, factorlist, nullptr);
	if (rv)
		return factorlist.fmax();
	else
		return 0;
}

/* SumDigits(n,r): Sum of digits of n in base r. */
static Znum ComputeSumDigits(const Znum &n, const Znum &radix) {

	Znum argum = n, Temp;
	Znum result = 0;

	while (argum > 0)
	{
		Temp = argum % radix;
		result += Temp;
		argum /= radix;
	}
	return result;
}

/* RevDigits(n,r): finds the value obtained by writing backwards the digits of n in base r */
static Znum ComputeRevDigits(const Znum &n, const Znum &radix) {

	Znum argum = n, Temp;
	Znum result = 0;

	while (argum > 0)
	{
		result *= radix;
		Temp = argum % radix;
		result += Temp;
		argum /= radix;
	}
	return result;
}

/* get 'width' of n in bits */
static long long NoOfBits(const Znum &n) {
	auto result = mpz_sizeinbase(ZT(n), 2);  // calculate number of bits
	return result;
}

/* SHL: Shift left the number of bits specified on the right operand. If the
number of bits to shift is negative, this is actually a right shift. If result would
be too large an error is reported. */
static retCode ShiftLeft(const Znum &first, const Znum &bits, Znum &result) {
	if (bits > LLONG_MAX || bits < LLONG_MIN)
		return retCode::INVALID_PARAM;

	long long shift = MulPrToLong(bits);

	if (shift == 0) {
		result = first;
		return retCode::EXPR_OK;  // shift zero bits; nothing to do
	}

	// there is no built-in shift operator for Znums, so it is simulated
	// using multiplication or division
	if (shift > 0) {
		if (NoOfBits(first) + shift > 66439) // more than 66439 bits -> more than 20,000 decimal digits
			return retCode::INTERM_TOO_HIGH;
		mpz_mul_2exp(ZT(result), ZT(first), shift);
		return retCode::EXPR_OK;
	}
	else {   // note use of floor division
		mpz_fdiv_q_2exp(ZT(result), ZT(first), -shift);
		return retCode::EXPR_OK;
	}
}


/* IsPrime(n): returns zero if n is definitely composite,  -1 if it is a probable prime, */
static int PrimalityTest(const Znum &Value) {
	int rv;
	assert(Value >= 0);

#ifdef __MPIR_VERSION
	static gmp_randstate_t state;
	static bool first = true;
	if (first) {
		gmp_randinit_default(state);
		first = false;
	}
	rv = mpz_probable_prime_p(ZT(Value), state, 16, 0);
	//rv = mpz_likely_prime_p(ZT(Value), state, 0);
#else
	rv = mpz_probab_prime_p(ZT(Value), 16);
#endif 
	/* rv is 1 if n is probably prime, or 0 if n is definitely composite.*/
	if (rv == 1)
		return -1;  // number is probably prime (less than 1 in 2^16 chance it's not prime)
	else
		return 0;  // number is definately composite
}
/* same purpose as PrimalityTest but optimised for smaller numbers. 1st call can be very slow,
but subsequent calls are very quick. -1 means definately prime*/
static int PrimalityTestSmall(const long long Value) {
	assert(Value >= 0);
	/* deal with even values first, because getBit below only handles odd values*/
	if (Value == 2) return -1;		//  2 is only even prime
	if ((Value & 1) == 0) return 0;   // even numbers > 2 are not prime

	if (Value <= max_prime) {
		if (primeFlags == NULL || (long long)primeListMax < Value) {
			generatePrimes(max_prime);  // takes a while, but only needed on 1st call
		}
		if (!getBit(Value, primeFlags))
			return -1;		// number is prime
		else
			return 0;		// number is not prime
	}
	else
		return PrimalityTest(Value);  // should work, but should never be needed
}


/*  B(n): Previous probable prime before n */
static retCode ComputeBack(const Znum &n, Znum &p) {
	if (n < 3)
		return retCode::NUMBER_TOO_LOW;  // 2 is the smallest prime
	if (n == 3) {
		p = 2;
		return retCode::EXPR_OK;
	}
	Znum i;		// largest odd number below n
	i = ((n & 1) == 1) ? n - 2 : n - 1;  // i >= 3

	for (; i > 0; i -= 2) {
		if (PrimalityTest(i) == -1) {
			p = i;
			return retCode::EXPR_OK;
		}
	}
	abort();   //should never get here!!
}

/* calculate the number of primes <= value */
static long long primePi(const Znum &n) {
	if (n < 2) return 0;
	if (n == 2) return 1;     // 2 is the smallest prime number
	long long rv = 1;		  // include 2 in the count
	long long i;			// largest odd number <= n
	long long ncopy = MulPrToLong(n);

	i = ((ncopy & 1) == 1) ? ncopy : ncopy - 1;  // i >= 3
	/* this method is crude. takes about 5 seconds for n = 10^9 */
	for (; i >= 2; i -= 2) {
		if (PrimalityTestSmall(i) == -1) {  // test all odd numbers <= Value
			rv++;
		}
	}
	return rv;
}

/* FactConcat(m,n): Concatenates the prime factors (base 10) of num according to the mode
mode	Order of factors	Repeated factors
0		Ascending			No
1		Descending			No
2		Ascending			Yes
3		Descending			Yes
*/
static Znum  FactConcat(const Znum &mode, const Znum &num) {
	fList factorlist;
	const bool descending = ((mode & 1) == 1);
	const bool repeat = ((mode & 2) == 2);

	/* get factors of num */
	auto rv = factorise(num, factorlist, nullptr);
	if (rv)
		return factorlist.ConcatFact(descending, repeat);
	else
		return 0;  // unable to factorise number
}

/* calculate the number of ways an integer n can be expressed as the sum of 2
squares x^2 and y^2. The order of the squares and the sign of x and y is significant
*/
static Znum R2(const Znum &num) {
	if (num < 0)
		return 0;
	if (num == 0)
		return 1;
	if (num % 4 == 3)
		return 0;   // at least 1 4k+3 factor has an odd exponent

	fList factorlist;
	Znum b = 1;

	/* get factors of num */
	auto rv = factorise(num, factorlist, nullptr);
	return factorlist.R2();
}

/* The number of representations of n as the sum of two squares ignoring order and signs*/
static Znum R2p(const Znum& num) {
	if (num < 0)
		return 0;
	if (num == 0)
		return 1;
	if (num % 4 == 3)
		return 0;   // at least 1 4k+3 factor has an odd exponent

	fList factorlist;
	Znum b = 1;

	/* get factors of num */
	auto rv = factorise(num, factorlist, nullptr);
	return factorlist.R2p();
}

/* return x^n */
static Znum powerBi(const __int64 x, unsigned __int64 n) {
	Znum result;
	mpz_ui_pow_ui(ZT(result), x, n);
	return result;
}
static Znum powerBi(const Znum &x, unsigned __int64 n) {
	Znum result;
	mpz_pow_ui(ZT(result), ZT(x), n);
	return result;
}

/* Calculate maximum square divisor of num and divide num by this.
return adjusted num, square divisor, and factor list of divisor.*/
void squareFree(Znum &num, Znum &sq, std::vector<zFactors> &sqf) {
	fList factorlist;
	zFactors temp;
	auto rv = factorise(num, factorlist, nullptr);

	factorlist.sqfree(sqf);  // get factors of square in sqf

	sq = 1;
	Znum x;
	for (auto i : sqf) {  /* calculate value of square */
		mpz_pow_ui(ZT(x), ZT(i.Factor), i.exponent);
		sq *= x;
	}
	num /= sq;  // adjust value of num
}


extern unsigned __int64 R3(__int64 n);
/* calculate the number of ways an integer n can be expressed as the sum of 3
squares x^2, y^2 and z^2. The order of the squares is significant. x, y and z can
be +ve, 0 or -ve See https://oeis.org/A005875 & https://oeis.org/A005875/b005875.txt*/
static Znum R3(Znum num) {

	if (num < 20'000'000'000'000'000) {
		__int64 llnum = MulPrToLong(num);
		return R3(llnum);
	}
	Znum sum = 0, sq, multiplier = 1;
	std::vector <zFactors> sqf;   // factor list for square factor of num

	if (num < 0)
		return 0;
	if (num == 0)     // test here necessary to avoid infinite loop
		return 1;
	while ((num & 3) == 0)
		num >>= 2;      // remove even factors. note that R3(4n) =R3(n)
	if (num % 8 == 7)
		return 0;     // take short cut if possible
	squareFree(num, sq, sqf);

	for (Znum k = 1; k*k <= num; k++) {
		sum += 2 * R2(num - k * k);
	}
	sum += R2(num);  // note: this time (for k=0) we DON'T multiply R2 by 2
	/* we now have sum = R3(n) */

	if (sq > 1) {
		/* add back factors that were removed to make n squarefree, one prime at a time.
		see web.maths.unsw.edu.au/~mikeh/webpapers/paper63.ps specifically equation (3) */
		for (int i = 0; i < sqf.size(); i++) {
			auto p = sqf[i].Factor;      // prime 
			auto λ = sqf[i].exponent / 2;  // exponent/2 (original exponent must be even)
			Znum mnum = -num;
			multiplier *= (powerBi(p, λ + 1) - 1) / (p - 1)
				- mpz_jacobi(ZT(mnum), ZT(p))*(powerBi(p, λ) - 1) / (p - 1);
			num *= powerBi(p, λ * 2);
		}
		sum *= multiplier;
	}

	return sum;
}

/* find smallest primitive root of num. return -1 for error 
see https://en.wikipedia.org/wiki/Primitive_root_modulo_n */
static Znum primRoot(const Znum &num) {
	fList factorlist, totF;
	std::vector <Znum> powers;  /* values of exponent that we need to use in testing */
	Znum mp;    /* mp = a^p mod num*/
	bool skip = false;

	if (num <= 1)
		return -1;    /* don't allow -ve, 0, or 1 values */
	/* get factors of num */
	auto rv = factorise(num, factorlist, nullptr);
	assert(rv);
	if (!factorlist.hasPrimitiveRoot()) {
		return -1;  /* there are no primitive roots */
	 }
	Znum tot = factorlist.totient();   /* get totient of num */
	rv = factorise(tot, totF, nullptr); /* get factors of totient */
	assert(rv);
	for (int x = 0; x < totF.fsize(); x++)
		powers.push_back(tot / totF.f[x].Factor);
	/* powers contains 1 value for each unique prime factor of the totient. */

	/* test numbers 2 to num-1 */
	for (Znum a = 2; a < num; a++) {
		if (gcd(a, num) != 1 || isPerfectSquare(a))
			continue; /* maybe save time by perfect square test? */
		skip = false;
		/* note that we only need to test certain values of p. We don't need to
		 test all values from 1 to num. */
		for (auto p : powers) {
			mp = modPower(a, p, num);   /* mp = a^p mod num*/
			if (mp == 1) {
				/* a is not a primitive root */
				skip = true;
				break;
			}
		}
		if (!skip)
			return a;  /* a passes tests; it is a primitive root */
	}

	return -1;  /* should never reach this point */
}

/* process one operator with 1 or 2 operands.
NOT, unary minus and primorial  have 1 operand.
All the others have two. Some operators can genererate an error condition
e.g. DIVIDE_BY_ZERO otherwise return EXPR_OK. 
For functions, do any further checks needed on the parameter values, then 
evaluate the function. gcd and lcm functions have a variable number of parameters. 
Some functions can generate error codes. */
static retCode ComputeSubExpr(const opCode stackOper, const std::vector <Znum> &p, Znum &result) {

	int rv;
	Znum temp;
	retCode retcode = retCode::EXPR_OK;

	switch (stackOper) {

	case opCode::comb: {  // calculate nCk AKA binomial coefficient
		if (p[1] > INT_MAX)
			return retCode::NUMBER_TOO_HIGH;
		if (p[1] < 1)
			return retCode::NUMBER_TOO_LOW;
		long long k = MulPrToLong(p[1]);
		mpz_bin_ui(ZT(result), ZT(p[0]), k);
		return retCode::EXPR_OK;
	}
	case opCode::plus: {
		result = p[0] + p[1];
		return retCode::EXPR_OK;
	}
	case opCode::minus: {
		result = p[0] - p[1];
		return retCode::EXPR_OK;
	}
	case opCode::unary_minus: {
		result = -p[0];
		return retCode::EXPR_OK;
	}
	case opCode::divide: {
		if (p[1] == 0)
			return retCode::DIVIDE_BY_ZERO;  // result would be infinity
		/* use truncation division */
		result = p[0] / p[1];
		return retCode::EXPR_OK;
	}
	case opCode::multiply: {
		auto resultsize = NoOfBits(p[0]) + NoOfBits(p[1]);
		if (resultsize > 99960)  // more than 99960 bits -> more than 30,000 decimal digits
			return retCode::INTERM_TOO_HIGH;
		result = p[0] * p[1];
		return retCode::EXPR_OK;
	}
	case opCode::remainder: {
		if (p[1] == 0)
			return retCode::DIVIDE_BY_ZERO;  // result would be infinity
		//result = p[0] % p[1];
		/* use truncation division */
		mpz_tdiv_r(ZT(result), ZT(p[0]), ZT(p[1]));
		return retCode::EXPR_OK;
	}
	case opCode::power: {
		if (p[1] > INT_MAX)
			return retCode::EXPONENT_TOO_LARGE;
		if (p[1] < 0)
			return retCode::EXPONENT_NEGATIVE;
		long long exp = MulPrToLong(p[1]);
		auto resultsize = (NoOfBits(p[0]) - 1)* exp;  // estimate number of bits for result
		if (resultsize > 99960)  // more than 99960 bits -> more than 30,000 decimal digits
			return retCode::INTERM_TOO_HIGH;
		mpz_pow_ui(ZT(result), ZT(p[0]), exp);
		return retCode::EXPR_OK;
	}
	case opCode::equal: {
		if (p[0] == p[1])
			result = -1;
		else
			result = 0;
		return retCode::EXPR_OK;
	}
	case opCode::not_equal: {
		if (p[0] != p[1])
			result = -1;
		else
			result = 0;
		return retCode::EXPR_OK;
	}
	case opCode::greater: {
		if (p[0] > p[1])
			result = -1;
		else
			result = 0;
		return retCode::EXPR_OK;
	}
	case opCode::not_greater: {
		if (p[0] <= p[1])
			result = -1;
		else
			result = 0;
		return retCode::EXPR_OK;
	}
	case opCode::less: {
		if (p[0] < p[1])
			result = -1;
		else
			result = 0;
		return retCode::EXPR_OK;
	}
	case opCode::not_less: {
		if (p[0] <= p[1])
			result = -1;
		else
			result = 0;
		return retCode::EXPR_OK;
	}
	case opCode::shl: {
		return ShiftLeft(p[0], p[1], result);
	}
	case opCode::shr: {
		// invert sign of shift
		return ShiftLeft(p[0], -p[1], result);
	}
	case opCode::not: /* Perform binary NOT */ {   
		//result = -1 - p[0];  // assumes 2s complement binary numbers
		mpz_com(ZT(result), ZT(p[0]));
		return retCode::EXPR_OK;
	}
	case opCode::and: /* Perform binary AND. */ {  
		mpz_and(ZT(result), ZT(p[0]), ZT(p[1]));
		return retCode::EXPR_OK;
	}
	case opCode::or:  /* Perform binary OR. */ {   
		mpz_ior(ZT(result), ZT(p[0]), ZT(p[1]));
		return retCode::EXPR_OK;
	}
	case opCode::xor: /*  Perform binary XOR. */ {   
		mpz_xor(ZT(result), ZT(p[0]), ZT(p[1]));
		return retCode::EXPR_OK;
	}
	case opCode::fact: {
		/* hard-coded limits allow size limit check before calculating the factorial */
		int limits[] = { 0, 5983, 11079, 15923, 20617, 25204, 29710, 34150,
			 38536, 42873, 47172 };
		if (p[0] < 0)
			return retCode::NUMBER_TOO_LOW;
		if (p[0] > LLONG_MAX)
			return retCode::NUMBER_TOO_HIGH;

		long long temp = llabs(MulPrToLong(p[0]));
		long long t2 = MulPrToLong(p[1]);
		if (t2 < sizeof(limits) / sizeof(limits[0]) && temp > limits[t2])
			/* more than 20,000 digits in base 10 */
			return retCode::INTERM_TOO_HIGH;

		mpz_mfac_uiui(ZT(result), temp, t2);  // get multi-factorial
		if (ZT(result)->_mp_size > 1039)
			/* more than 20,000 digits in base 10 */
			return retCode::INTERM_TOO_HIGH;

		return retCode::EXPR_OK;
	}
	case opCode::prim: {
		if (p[0] > 46340)
			return retCode::INTERM_TOO_HIGH;
		if (p[0] < 0)
			return retCode::NUMBER_TOO_LOW;

		long long temp = llabs(MulPrToLong(p[0]));
		mpz_primorial_ui(ZT(result), temp);  // get primorial
		return retCode::EXPR_OK;
	}

	case opCode::fn_gcd: /* GCD */ {
		//mpz_gcd(ZT(result), ZT(p[0]), ZT(p[1]));
		result = p[0];
		/* gcd has 1 or mor parameters */
		for (int ix = 1; ix < p.size(); ix++) {
			result = gcd(result, p[ix]);
		}
		break;
	}
	case opCode::fn_lcm: /* Least Common Multiple */ {
		result = p[0];
		/* lcm has 1 or more parameters */
		for (int ix = 1; ix < p.size(); ix++) {
			result = lcm(result, p[ix]);
		}
		break;
	}
	case opCode::fn_modpow: {						// MODPOW
		if (p[2] == 0)
			return retCode::DIVIDE_BY_ZERO;
		if (p[1] < 0) {
			if (gcd(p[0], p[2]) != 1)
				return retCode::ARGUMENTS_NOT_RELATIVELY_PRIME;  // p[0] and p[2] not mutually prime
		}
		/* note: negative exponent is only allowed if p[0] & p[2] are mutually prime,
		i.e modular inverse of p[0] wrt p[2] exists. */
		mpz_powm(ZT(result), ZT(p[0]), ZT(p[1]), ZT(p[2]));
		break;
	}
	case opCode::fn_modinv: /* modular inverse */ {
		/* if an inverse doesn’t exist the return value is zero and rop is undefined*/
		rv = mpz_invert(ZT(result), ZT(p[0]), ZT(p[1]));
		if (rv == 0) {
			return retCode::ARGUMENTS_NOT_RELATIVELY_PRIME;
		}
		break;
	}

	case opCode::fn_totient: {			// totient
		if (p[0] < 1) 
			return retCode::NUMBER_TOO_LOW;;
		result = ComputeTotient(p[0]);
		break;
	}
	case opCode::fn_numdivs: {		// NUMDIVS
		if (p[0] < 1) {
			return retCode::NUMBER_TOO_LOW;
		}
		result = ComputeNumDivs(p[0]);
		break;
	}
	case opCode::fn_sumdivs: {		// SUMDIVS
		result = ComputeSumDivs(p[0]);
		break;
	}

	case opCode::fn_sumdigits: /* Sum of digits of p[0] in base r.*/ {	// SumDigits(n, r) : 
		result = ComputeSumDigits(p[0], p[1]);
		break;
	}
	case opCode::fn_numdigits: /* number of digits of p[0] in base p[1] */ {
		result = ComputeNumDigits(p[0], p[1]);
		break;
	}
	case opCode::fn_revdigits: {	// revdigits
		result = ComputeRevDigits(p[0], p[1]);
		break;
	}

	case opCode::fn_isprime: {  // isprime
		/* -1 indicates probably prime, 0 = composite */
		int rv = PrimalityTest(abs(p[0]), 0);
		/* rv: 0 = probable prime.
               1 = composite: not 2-Fermat pseudoprime.
               2 = composite: does not pass 2-SPRP test.
               3 = composite: does not pass BPSW test, but passes other tests.*/
		if (verbose > 0) {
			if (rv == 2)
				std::cout << "pseudo-prime; passes 2-fermat test but not 2-SPRP test \n";
			if (rv == 3)
				std::cout << "pseudo-prime; passes 2-SPRP test but not BPSW test \n";
		}
		if (rv == 0)
			result = -1;   /* prime */
		else
			result = 0;    /* composite */
		break;
	}
	case opCode::fn_fib: /* fibonacci */ {	
		if (p[0] > 95700 || p[0] < -95700)
		{  /* result would exceed 20,000 digits */
			return retCode::INTERM_TOO_HIGH;
		}
		long long temp = MulPrToLong(ZT(p[0]));
		bool neg = false;
		if (temp < 0) {
			/* is it a "negafibonacci" number? */
			if ((temp & 1) == 0) { neg = true; } /* set for even -ve number */
			temp = -temp;
		}
		mpz_fib_ui(ZT(result), temp);  // calculate fibonacci number
		if (neg) { result = -result; } /* flip sign for even -ve number */
		break;
	}
	case opCode::fn_luc: /* lucas number */ {
		if (p[0] < 0) {
			return retCode::NUMBER_TOO_LOW;
		}
		if (p[0] > 95700)
		{
			return retCode::INTERM_TOO_HIGH;
		}
		long long temp = MulPrToLong(p[0]);
		mpz_lucnum_ui(ZT(result), temp);  // calculate lucas number
		break;
	}

	case opCode::fn_part: /* number of partitions */ {
		if (p[0] > 1000000) {
			return retCode::NUMBER_TOO_HIGH;
			// note: biperm is limited to values <= 1,000,000
		}
		if (p[0] < 0) {
			return retCode::NUMBER_TOO_LOW;
			// note: biperm is limited to values <= 1,000,000
		}
		int temp = (int)MulPrToLong(p[0]);
		biperm(temp, result);   // calculate number of partitions
		break;
	}
	case opCode::fn_np: /* next prime */ {
		mpz_nextprime(ZT(result), ZT(p[0]));  // get next prime
		break;
	}
	case opCode::fn_pp: /* previous prime */ {
		retcode = ComputeBack(p[0], result);  // get previous prime
		if (retcode != retCode::EXPR_OK) {
			return retcode;   // error: number < 3
		}
		break;
	}

	case opCode::fn_primePi: /* count primes <= n */ {
		if (p[0] > max_prime) {
			return retCode::NUMBER_TOO_HIGH;
		}
		result = primePi(p[0]);
		break;
	}
	case opCode::fn_concatfact:  /*Concatenates the prime factors of n according to
						  the mode in m */ {
		if (p[0] < 0 || p[0] > 3) {
			return retCode::INVALID_PARAM;  // mode value invalid
		}
		result = FactConcat(p[0], p[1]);
		break;
	}
	case opCode::fn_r2: {
		result = R2(p[0]);
		break;
	}
	case opCode::fn_r2p: {
		result = R2p(p[0]);
		break;
	}
	case opCode::fn_r3: {
		result = R3(p[0]);
		break;
	}
	case opCode::fn_r3h: {
		if (p[0] < 0 || p[0] > LONGLONG_MAX) {
			return retCode::INVALID_PARAM;  // parameter out of range
		}
		result = R3h(MulPrToLong(p[0]));
		break;
	}
	case opCode::fn_hurwitz: {
		result = Hclassno12(p[0]);  /* returns 12 x hurwitz class number */
		break;
	}

	/* legendre & kronecker are in fact implemented in MPIR as aliases of jacobi */
	case opCode::fn_legendre:
	case opCode::fn_jacobi:
	case opCode::fn_kronecker: {
		result = mpz_jacobi(ZT(p[0]), ZT(p[1]));
		break;
	}

	case opCode::fn_llt: /* lucas-lehmer primality test */ {
		/* see https://en.wikipedia.org/wiki/Lucas%E2%80%93Lehmer_primality_test */
		if (p[0] > INT_MAX)
			return retCode::NUMBER_TOO_HIGH;
		if (p[0] < 0)
			return retCode::NUMBER_TOO_LOW;

		result = llt(p[0]);
		if (verbose > 0 || p[0] >= 216091) {
			if (result == 1)
				std::cout << "*** 2^" << p[0] << " -1 is prime! ***\n";
			else
				std::cout << "2^" << p[0] << " -1 is not prime \n";
		}
		break;
	}
	case opCode::fn_sqrt: {
		if (p[0] < 0) return retCode::INVALID_PARAM;
		mpz_sqrt(ZT(result), ZT(p[0]));  /* result = square root of p[0]*/
		break;
	}
	case opCode::fn_nroot: {
		/* for real numbers nroot (x, n)  = x^(1/n) has a discontinuity at n=0,
		so the function is considered to be undefined for n=0 */
		if (p[1] == 0) return retCode::INVALID_PARAM;

		/* odd root of a -ve number is -ve so allow it. Even root would be a complex
		number so don't allow it */
		if (p[0] < 0 && mpz_even_p(ZT(p[1])))
			return retCode::INVALID_PARAM;

		if (p[1] < 0) {
			result = 0; /* for real numbers nroot(x, n) with n -ve = 1/(x^(-1/n))
						   If x > 1 the floor of this value is 0 */
			break;
		}
		if (p[1] > INT_MAX) {
			/* if p[1] is very large the nth root would be close to 1  for +ve p[0]
			or -1 for -ve p[0] and odd p[1] */
			return retCode::NUMBER_TOO_HIGH;
		}

		long long temp = MulPrToLong(ZT(p[1]));
		mpz_root(ZT(result), ZT(p[0]), temp);  /* result = nth root of p[0]*/
		break;
	}
	case opCode::fn_bpsw: /* Baillie-Pomerance-Selfridge-Wagstaff probabilistic
		primality test */ {
		if (p[0] <= 1) 
			return retCode::NUMBER_TOO_LOW;

		result = mpz_bpsw_prp(ZT(p[0]));
		if (result == 0)
			std::cout << "composite \n";
		else if (result == 1)
			std::cout << "probable prime \n";
		else if (result == 2)
			std::cout << "prime \n";
		break;
	}
	case opCode::fn_aprcl: /* Adleman–Pomerance–Rumely primality test  */ {

		/*  the Adleman–Pomerance–Rumely primality test is an algorithm for 
		determining whether a number is prime. Unlike other, more efficient 
		algorithms for this purpose, it avoids the use of random numbers, so it 
		is a deterministic primality test. It is named after its inventors, 
		Leonard Adleman, Carl Pomerance, and Robert Rumely. 
		It was later improved by Henri Cohen and Hendrik Willem Lenstra, commonly 
		referred to as APR-CL*/

		if (p[0] <= 1) 
			return retCode::NUMBER_TOO_LOW;
		result = mpz_aprtcle(ZT(p[0]), verbose);
		if (result == 0)
			printf_s("composite \n");
		else if (result == 1)
			printf_s("probable prime \n");
		else if (result == 2)
			printf_s("prime \n");
		break;
	}
	case opCode::fn_numfact: /* number of factors */ {
		result = ComputeNumFact(p[0]);
		break;
	}
	case opCode::fn_minfact: /* smallest factor*/ {
		result = ComputeMinFact(p[0]);
		break;
	}
	case opCode::fn_maxfact: /* largest factor */ {
		result = ComputeMaxFact(p[0]);
		break;
	}
	case opCode::fn_ispow: /* check whether or not p[0] is a perfect power */ {
		/* return -1 if p[0] is a perfect power, otherwise 0 */
		long long MaxP = 393'203;  // use 1st  33333 primes
		if ((long long)primeListMax < MaxP) {  // get primes
			generatePrimes(MaxP);  // takes a while, but only needed on 1st call
		}
		Znum base;
		long long exp = PowerCheck(p[0], base, 2);
		if (exp > 1) {
			std::cout << p[0] << " = " << base << "^" << exp << '\n';
			result = -1;
		}
		else {
			result = 0;
		}
		break;
	}
	case opCode::fn_modsqrt: /* modular square root */ {
		/* Solve the equation given p[0] and p[1].  x^2 ≡ p[0] (mod p[1]) */
		if (p[1] <= 0)
			return retCode::NUMBER_TOO_LOW;  /* modulus must be +ve */
		roots = ModSqrt(p[0], p[1]);
		if (verbose > 0) {
			if (roots.empty())
				gmp_printf("modsqrt(%Zd, %Zd) has no roots \n", p[0], p[1]);
			else {
				gmp_printf("modsqrt(%Zd, %Zd) = ", p[0], p[1]);
				for (Znum r : roots)
					gmp_printf("%Zd, ", r);
				putchar('\n');
			}
		}
		/* the result can be: no solution: roots is empty
							  or one  or more solutions */
		if (roots.empty())
			return retCode::INVALID_PARAM;  /* no solution exists */

		result = roots[0];     /* ignore 2nd solution, if any */
		multiValue = true;     /* indicate multiple return values */

		break;
	}
	case opCode::fn_invtot: /* inverse totient */ {
		std::vector<unsigned long long> *resultsP;

		/* have to limit p[0] to small values, otherwise risk running out of memory.
		   this also eliminates any possibility of integer overflow. */
		if (p[0] > 1'000'000'000'000'000)
			return retCode::NUMBER_TOO_HIGH;
		if (p[0] <= 0)
			return retCode::NUMBER_TOO_LOW;
		/* get list of numbers x1, x2, ... such that totient(x) = p[0].
		if p[0] is zero InverseTotient just clears its cache. */
		auto size = inverseTotient(MulPrToLong(p[0]), &resultsP, false, 0, false);
		if (size == 0) {
			if (verbose > 0)
				std::cout << "Inverse Totient has no solutions \n";
			return retCode::INVALID_PARAM;
		}
		roots.clear();
		for (size_t i = 0; i < size; i++) {
			roots.push_back((*resultsP)[i]); /* copy results */
		}
		result = roots[0];
		multiValue = true;     /* indicate multiple return values */
		break;
	}
	case opCode::fn_divisors: /* list of divisors */ {
		if (p[0] < 1) {
			return retCode::NUMBER_TOO_LOW;
		}
		result = DivisorList(p[0], roots);
		multiValue = true;     /* indicate multiple return values */
		break;
	}
	case opCode::fn_primroot: /* lowest primitive root */ {
		result = primRoot(p[0]);
		if (result <= 0)
			return retCode::INVALID_PARAM;
		else break;
	}

	default:
		abort();	// should never get here
	}

	return retcode;
}


/* find next , or ) but anything enclosed in nested brackets () is ignored */
static void nextsep(token expr[], int &ix) {
	int brackets = 0;
	for (ix = 0; ; ix++) {
		if (expr[ix].typecode == types::Operator
			&& expr[ix].oper == opCode::leftb)
			brackets++;
		if (expr[ix].typecode == types::Operator
			&& expr[ix].oper == opCode::rightb)
			if (brackets > 0)
				brackets--;
			else break;
		if (brackets > 0)
			continue; /* ignore anything enclosed in brackets */
		if (expr[ix].typecode == types::comma)
			break;
		if (expr[ix].typecode == types::end)
			break; /* safety check; if brackets match up this will never be reached */
	}
}

/* print tokens in readable form */
static void printTokens(const std::vector <token> expr) {
	if (expr.size() == 0)
		return;   /* if no tokens to print do nothing but exit*/
	for (int ix = 0; ix < expr.size(); ix++) {
		switch (expr[ix].typecode) {
		case types::number:
			std::cout << expr[ix].value << ' ';
			break;
		case types::func:
			std::cout << functionList[expr[ix].function].fname << ' ';
			break;
		case types::comma:
			std::cout << ',';
			break;

		case types::end:
			std::cout << " **END** ";
			break;

		case types::error:
			std::cout << " **ERROR** ";
			break;

		case types::Operator: {
			if (expr[ix].oper == opCode::fact) {
				for (int i = 1; i <= expr[ix].value; i++)
					std::cout << '!';
				std::cout << ' ';
			}
			else {
				int ixx = expr[ix].function; 
				/* search for oper code in list of operators */
				std::cout << operators[ixx].oper << ' ';
			}
			break;
		}

		case types::uservar: {
			size_t Userix = expr[ix].userIx;  /* get index value from user variable */
			std::string name = uvars.vars[Userix].name; /* get name of user variable */
			std::cout  << name << " ";
			break;
		}

		default:
			abort(); /* unrecognised token */
		}
	}
	std::cout << '\n';
}

/* evaluate expression and return value in Result
operators and functions include:
symbol		meaning					operator Priority
n!			factorial               0
n!..!	    multi-factorial         0
p#			primorial               0
			(product of all primes
			less or equal than p)

-			unary minus				1    (right to left)
NOT			bitwise operators		1    (right to left)
** or ^		Exponentiate			2    (right to left)
*			multiply				3
/			divide					3
%			modulus (remainder)		3
C			binomial coefficient	4
+			add						5
-			subtract				5
n SHL b		shift n left b bits		6
n SHR b		shift n right b bits	6

comparison operators
>			return zero for false 	7
>=			and -1 for true			7
<									7
<=									7
==									8
!=			not equal				8

AND			bitwise operators		9
OR									11
XOR									10


________________Functions___________________________________
GCD(a, b)
MODPOW(m,n,r): finds m^n modulo r
MODINV(m,n): inverse of m modulo n, only valid when gcd(m,n)=1
SUMDIGITS(n,r): Sum of digits of n in base r.
NUMDIGITS(n,r): Number of digits of n in base r
REVDIGITS(n,r): finds the value obtained by writing backwards the digits of n in base r.
ISPRIME(n)
F(n):		Fibonacci
L(n):		Lucas
P(n):		Unrestricted Partition Number (number of decompositions of n
			into sums of integers without regard to order)
PI(n):		Number of primes <= n
N(n):		Next probable prime after n
B(n):		Previous probable prime before n
Totient(n): finds the number of positive integers less than n which are relatively prime to n.
NumDivs(n): Number of positive divisors of n either prime or composite.
SumDivs(n): Sum of positive divisors of n either prime or composite.
FactConcat(m,n): Concatenates the prime factors of n according to the mode expressed in m

Hexadecimal numbers are preceded by 0X or 0x
*/

/* search for operators e.g. '+', '*'  '!', 'NOT'
return index of operator or -1 */
static void operSearch(const std::string &expr, int &opcode) {
	opcode = -1;
	for (int ix = 0; ix < sizeof(operators) / sizeof(operators[0]); ix++) {
		/* use case-insensitive compare */
		if (_strnicmp(expr.data(), operators[ix].oper, strlen(operators[ix].oper)) == 0) {
			opcode = ix;
			return;
		}
	}
	opcode = -1;  // does not match any operator in list
	return;
}


/* convert an expression which is already tokenised to reverse polish.
Syntax checking is done, but it is not guaranteed that all syntax errors will
be detected. If an error is detected return EXIT_FAIL, otherwise EXIT_SUCCESS.
The well-known shunting algorith is used, but for function parameters a recursive
call to reversePolish is made, which also takes care of nested function calls.
Also, the initial tokenisation does not distinguish unary - from normal -, so
that is deduced from the sequence of tokens and the op code is changed if necessary.
Unary + is recognised as a special case which requires no output */
static int reversePolish(token expr[], const int exprLen, std::vector <token> &rPolish) {
	int exprIndex = 0;

	std::vector <token> operStack;  /* operator stack */
	bool leftNumber = false;      /* used for syntax checking & to recognise unary + or - */

	while (exprIndex < exprLen) {
		if (expr[exprIndex].typecode == types::end)
			break;			/* reached end of tokens so stop */

		switch (expr[exprIndex].typecode) {

		case types::number: {
			if (leftNumber)
				return EXIT_FAILURE;  /* syntax error */

			rPolish.push_back(expr[exprIndex]);
			leftNumber = true;
			break;
		}

		/* a user variable is treated the same as a number */
		case types::uservar: {
			if (leftNumber)
				return EXIT_FAILURE;  /* syntax error */

			rPolish.push_back(expr[exprIndex]);
			leftNumber = true;
			break;
		}

		case types::func: {
			/* process function. 1st get number of parameters */
			if (leftNumber)
				return EXIT_FAILURE; /* syntax error */
			int numparams = expr[exprIndex].numops;

			if (expr[exprIndex + 1].typecode != types::Operator ||
				expr[exprIndex + 1].oper != opCode::leftb) {
				return EXIT_FAILURE; /* function name not followed by (*/
			}
			int paramLen = 0;
			int ix3 = 0;
			int pcount = 1;
			for (; ; pcount++) {

				/* find delimiter marking end of parameter*/
				nextsep(expr + exprIndex + 2 + ix3, paramLen); // get next , or )
				if (expr[exprIndex + 2 + ix3 + paramLen].typecode == types::end)
					return EXIT_FAILURE; /* parameter sep. not found */
				if (paramLen == 0)
					return EXIT_FAILURE;  /* syntax error */
				int rv = reversePolish(expr + exprIndex + 2 + ix3, paramLen, rPolish);
				if (rv != EXIT_SUCCESS)
					return rv;  /* syntax error? */
				ix3 += (paramLen + 1); // move ix3 past , or )

				if (expr[exprIndex + 1 + ix3].typecode == types::Operator &&
					expr[exprIndex + 1 + ix3].oper == opCode::rightb)
					break; /* ) found */
			}
			if (numparams == SHORT_MAX) {
				/* variable number of parameters; use actual number */
				expr[exprIndex].numops = pcount;  
			}
			else if (pcount != numparams)
				return EXIT_FAILURE;  /* wrong number of parameters */

			rPolish.push_back(expr[exprIndex]); /* copy function token to output */
			exprIndex += (ix3 + 1); /* move past ) after function name */
			leftNumber = true;
			break;
		}

		case types::Operator: {
			if (expr[exprIndex].oper != opCode::leftb
				&& expr[exprIndex].oper != opCode::rightb) {

				if (expr[exprIndex].oper == opCode::minus && !leftNumber) {
					expr[exprIndex].oper = opCode::unary_minus;  /* adjust op code */
					expr[exprIndex].numops = 1;            /* adjust number of operands */
					/* unary - is last operator in list, change index in token from - to unary - */
					expr[exprIndex].function = sizeof(operators) / sizeof(operators[0]) - 1; 
				}

				bool left = operators[expr[exprIndex].function].left; /* assocativity*/
				bool pre = operators[expr[exprIndex].function].pre;  /* true when unary operator precedes expr*/
				bool unary = operators[expr[exprIndex].function].numOps == 1; /* true for unary operator */

				/* get priority of current operator */
				int expOpPri = operators[expr[exprIndex].function].pri;

				/* check for unary operator - or + or NOT */
				if (!leftNumber) {  /* if operator is not preceded by an expression */
					if (expr[exprIndex].oper == opCode::plus) {
						break; /* step past unary plus; no output required*/
					}
					else if (!unary || !pre)
						return EXIT_FAILURE;  /* non-unary operator not preceded by a number or expression */
				}
				else {
					/* operator follows expression */
					if (unary && pre)
						return EXIT_FAILURE; /* unary operator preceded by a number or expression */
				}

				/* Transfer higher priority operators from stack to output.
					Stop when lower priority operator or ( is found on stack.
					For equal priority operators the left-associative flag is checked. */
				while (operStack.size() > 0
					&& operStack.back().typecode == types::Operator
					&& operStack.back().oper != opCode::leftb) {
					/* get priority of top operator on stack */
					int stkOpPri = operators[operStack.back().function].pri;
					/* N.B. lower priority value; higher priority operator */
					if ((stkOpPri < expOpPri)
						|| (stkOpPri == expOpPri && left)) {
						/* transfer high priority stack operator to output.
						exponent operator is special because it is right-associative */
						rPolish.push_back(operStack.back());
						operStack.pop_back();
					}
					else
						break; /* stop when lower priority operator found on stack */
				}

				operStack.push_back(expr[exprIndex]); /* put current operator onto stack */

				if (unary && !pre)
					leftNumber = true;  /* factorial, primorial operators follow the
										expression they apply to*/
				else
					leftNumber = false;
			}

			else if (expr[exprIndex].oper == opCode::leftb) {
				if (leftNumber)
					return EXIT_FAILURE; /* syntax error */
				operStack.push_back(expr[exprIndex]);
			}

			else if (expr[exprIndex].oper == opCode::rightb) {
				if (!leftNumber)
					return EXIT_FAILURE; /* syntax error */
				/* right bracket; remove stacked operators up to left bracket */
				while (operStack.size() > 0
					&& operStack.back().typecode == types::Operator
					&& operStack.back().oper != opCode::leftb) {
					rPolish.push_back(operStack.back());
					operStack.pop_back();
				};
				if (operStack.empty())
					return EXIT_FAILURE;  /* missing ( */
				operStack.pop_back(); /* discard left bracket */
				leftNumber = true;
			}

			break;
		}

		default:
			return EXIT_FAILURE;  /* unkown token in input */
		}

		exprIndex++;  /* move index to next token */
	}

	/* transfer any remaining operators from stack to output */
	while (operStack.size() > 0) {
		rPolish.push_back(operStack.back());
		operStack.pop_back();
	}
	if (leftNumber)
		return EXIT_SUCCESS;  /* expression appears to be syntactically valid */
	else
		return EXIT_FAILURE;  /* syntax error e.g.  expression "2+" would trigger this */
}


/* evaluate an expression in reverse polish form. Returns EXPR_OK or error code.
If there is more than one number on the stack at the end, or at any time there
are not enough numbers on the stack to perform an operation an error is reported.
(this would indicate a syntax error not detected earlier) 
If the final operation is a function call that returns multiple values,
multiValue is set to true, otherwise it is false*/
static retCode evalExpr(const std::vector<token> &rPolish, Znum & result, bool *multiV) {
	std::stack <token> nums;   /* this stack holds both numbers and user variables */
	Znum val;
	std::vector <Znum> args;
	int index = 0;
	int rPlen = (int)rPolish.size();
	retCode retcode;
	token temp;

	temp.typecode = types::number;

	for (index = 0; index < rPlen; index++) {
		multiValue = false;     /* changed to true if function returns multiple values */

		switch (rPolish[index].typecode) {
		case types::number:
		case types::uservar:
			{  	/* push number onto stack */
				nums.push(rPolish[index]);
				break;
			}
	
		/* operators and functions are processed by taking the operand values
			from the stack, executing the operation or function and putting
			the returned value onto the stack. If there are insuffficient values
			on the stack or if the operator or function returns an error code
			exit immediately. */
		case types::func: 
		case types::Operator:
			{
				opCode oper = rPolish[index].oper;

				int NoOfArgs = rPolish[index].numops;
				if (NoOfArgs > nums.size())
					return retCode::SYNTAX_ERROR;  /* not enough operands on stack*/

				if (oper == opCode::assign) {
					/* assignment operator is fully processed here */
					temp = nums.top();  /* remove top token from stack */
					nums.pop();

					if (nums.top().typecode != types::uservar)
						return retCode::SYNTAX_ERROR;
					size_t Userix = nums.top().userIx;

					/* store new value in user variable */
					if (temp.typecode == types::number)
						uvars.vars[Userix].data = temp.value;
					else if (temp.typecode == types::uservar)
						uvars.vars[Userix].data = uvars.vars[temp.userIx].data;
					else
						return retCode::SYNTAX_ERROR;  /* wrong type of token on stack */

					nums.pop();  /* remove variable from stack */
					nums.push(temp);  /* put value back on stack */
				}
				else {
					args.clear();
					for (; NoOfArgs > 0; NoOfArgs--) {
						/* copy operand(s) from number stack to args */
						if (nums.top().typecode == types::number)
							args.insert(args.begin(), nums.top().value); /* use top value from stack*/
						else {
							size_t Userix = nums.top().userIx;      /* use top variable from stack's value */
							args.insert (args.begin(), uvars.vars[Userix].data);
						}
						nums.pop();  /* remove top value from stack */
					}

					if (oper == opCode::fact) {
						/* get 2nd operand for multifactorial. In this special case
						 the 2nd operand is in the operator token, not a number token */
						args.push_back(rPolish[index].value);
					}
					retcode = ComputeSubExpr(oper, args, val);
					if (retcode != retCode::EXPR_OK)
						return retcode;
					temp.typecode = types::number;
					temp.value = val;  /* put value returned by function or operator */
					nums.push(temp);  /*  onto stack */
				}
				break;
			}
		default:
			abort(); /* unrecognised token */
		}
	}


	if (nums.size() == 1) {
		if (nums.top().typecode == types::number)
			result = nums.top().value;
		else {
			size_t Userix = nums.top().userIx;  /* get value from user variable */
			result = uvars.vars[Userix].data;
		}
		if (multiV != nullptr)
			*multiV = multiValue;
		return retCode::EXPR_OK;
	}
	else
		return retCode::SYNTAX_ERROR;  /* too many operands on stack*/
}

/*
Added 5/6/2021

   The 'engine'at the heart of the calculator was largely rewritten;
   It was divided into 3 parts:
   1.   'Tokenise' all terms in the expression i.e. each number, operator,
		Function name, bracket & comma is turned into a token. Also check that
		the opening and closing brackets pair up correctly.
	2.  Convert to Reverse Polish. This uses the well-known 'shunting' algorithm,
		but a recursive call to the reverse polish function is made for each
		function parameter (this also takes care of nested function calls).
		Also there is a tweak for the factorial, double factorial and primorial
		functions because the operator follows the number rather than precedes it.
		Some syntax checks are made but there is no guarantee that all syntax
		errors will be detected.
	3.  Calculate the value of the reverse polish sequence. If there is more than
		one number on the stack at the end, or at any time there are not enough
		numbers on the stack to perform an operation an error is reported.
		(this would indicate a syntax error not detected earlier)*/
retCode ComputeExpr(const std::string &expr, Znum &Result, int &asgCt, bool *multiV = nullptr) {
	retCode rv;
	std::vector <token> tokens;
	std::vector <token> rPolish;

	rv = tokenise(expr, tokens, asgCt); /* 'tokenise' the expression */
	if (rv == retCode::EXPR_OK) {
		rPolish.clear();
		int ircode = reversePolish(tokens.data(), (int)tokens.size(), rPolish);
		/* convert expression to reverse polish */
		if (ircode == EXIT_SUCCESS) {
			rv = evalExpr(rPolish, Result, multiV);
			if (rv != retCode::EXPR_OK && verbose > 0) {
				std::cout << "Expression could not be evaluated \n";
				printTokens(rPolish);
			}
		}
		else {
			if (verbose > 0) {
				std::cout << "** error: could not convert to reverse polish \n";
				printTokens(rPolish);
			}
			return retCode::SYNTAX_ERROR;
		}
	}
	else {
		if (verbose > 0) {
			std::cout << "** error: could not tokenise expression «"
				<< expr << "»\n";
			printTokens(tokens);
			std::cout << "expr contains \"" << expr << "\" \n";
			//printf_s("expr contains: \"%s\" \n", expr.c_str());
			//for (auto c : expr) {
			//	printf_s("%2x", c);
			//}
			//putchar('\n');
		}
	}

	return rv;
}

/* convert expression to ´tokens´. A token is basically either a number, operator,
function name or comma. Brackets are classed as operators.
Syntax is not checked properly at this stage but if there is anything that cannot
be tokenised SYNTAX_ERROR is returned. If left & right brackets don't match up
PAREN_MISMATCH is returned.
The normal return value is EXPR_OK.
*/
static retCode tokenise(const std::string expr, std::vector <token> &tokens, int &asgCt) {
	int exprIndex = 0;
	int opIndex;
	token nxtToken;

	tokens.clear();
	asgCt = 0;

	while (exprIndex < expr.length()) {
		nxtToken.typecode = types::error;   /* should be replaced later by valid code */
		int charValue = toupper(expr[exprIndex]);

		/* skip blanks */
		while (exprIndex <= (expr.length()) && isspace(charValue)) {
			exprIndex++;
			if (exprIndex >= expr.length())
				break;
			charValue = toupper(expr[exprIndex]);
		}
		if (exprIndex >= expr.length())
			break;

		/* first look for operators. */
		operSearch(expr.substr(exprIndex), opIndex);
		if (opIndex != -1) {
			/* found operator e.g. + - etc. */
			nxtToken.function = opIndex;
			nxtToken.typecode = types::Operator;
			nxtToken.oper = operators[opIndex].operCode;
			//nxtToken.numops = opr[(int)nxtToken.oper].numOps;
			nxtToken.numops = operators[opIndex].numOps;
			nxtToken.value = 0;
			exprIndex += (int)strlen(operators[opIndex].oper);  // move index to next char after operator
			if (operators[opIndex].operCode == opCode::fact) {
				nxtToken.value = 1;
				while (expr.substr(exprIndex, 1) == "!") {
					nxtToken.value++; /* get 2nd operand for multi-factorial*/
					exprIndex++;
				}
			}
			if (operators[opIndex].operCode == opCode::assign)
				asgCt++;
		}

		/* check for , */
		else if (charValue == ',') {
			nxtToken.typecode = types::comma;
			nxtToken.value = 0;                /* not used */
			nxtToken.oper = opCode::power;  /* not used */
			exprIndex++;
		}

		/* check for number */
		else if (charValue >= '0' && charValue <= '9') {
			/* convert number from ascii to Znum */
			int exprIndexAux = exprIndex;
			if (charValue == '0' && exprIndexAux < expr.length() - 2 &&
				toupper(expr[exprIndexAux + 1]) == 'X')
			{  // hexadecimal
				std::vector<char> digits;
				exprIndexAux++;
				while (exprIndexAux < expr.length() - 1) {
					charValue = expr[exprIndexAux + 1];
					if ((charValue >= '0' && charValue <= '9') ||
						(charValue >= 'A' && charValue <= 'F') ||
						(charValue >= 'a' && charValue <= 'f')) {
						exprIndexAux++;
						digits.push_back(charValue);
					}
					else {
						break;   // jump out of inner while loop
					}
				}
				digits.push_back('\0');   // null terminate string
				/* convert hex string to bigint */
				mpz_set_str(ZT(nxtToken.value), digits.data(), 16);
				nxtToken.typecode = types::number;
				nxtToken.oper = opCode::power;  /* not used */
				exprIndex = exprIndexAux + 1;
			}

			else {                   // Decimal number.
				std::vector<char> digits;
				while (exprIndexAux < expr.length()) {
					// find position of last digit
					charValue = expr[exprIndexAux];
					if (charValue >= '0' && charValue <= '9') {
						exprIndexAux++;
						digits.push_back(charValue);
					}
					else {
						break;    // jump out of inner while loop
					}
				}
				// Generate big integer from decimal number
				digits.push_back('\0');   // null terminate string
				mpz_set_str(ZT(nxtToken.value), digits.data(), 10);
				nxtToken.typecode = types::number;
				nxtToken.oper = opCode::power;  /* not used */
				exprIndex = exprIndexAux;
			}
		}

		else {
			/* try to match function name. Names are not case-sensitive */
			for (ptrdiff_t ix = 0; ix < sizeof(functionList)/sizeof(functionList[0]); ix++) {
				if (_strnicmp(&expr[exprIndex],
					functionList[ix].fname, strlen(functionList[ix].fname)) == 0) {
					/* we have a match to a function name */
					nxtToken.typecode = types::func;
					nxtToken.function = int(ix);        /* save index into functionList*/
					nxtToken.numops = functionList[ix].NoOfParams;     
					nxtToken.oper = functionList[ix].fCode;  
					exprIndex += (int)strlen(functionList[ix].fname); 
					            /* move exprIndex past function name */
					goto next;  /* go to finish processing this token, then 
								process rest of expr */
				}
			}

			/* if we drop through to here the only possibility left is a
			   user variable */
			int exprIndexAux = exprIndex;
			if (charValue == '_' || charValue == '$') {
				while (charValue == '_' || charValue == '$' || isalnum(charValue)) {
					exprIndex++;
					if (exprIndex >= expr.size())
						break;
					charValue = expr[exprIndex];
				}

				int length = exprIndex - exprIndexAux;
				char uid[40];  /* buffer contains name of user-defined identifier */
				Znum temp;

				if (length >= sizeof(uid))  /* user id cannot exceed 39 characters */
					return retCode::SYNTAX_ERROR;  /* unable to tokenise expression*/

				strcpy_s(uid, sizeof(uid), expr.substr(exprIndexAux, length).data());
				int rv = get_uvar(uid, temp); /* does variable exist already? */
				if (rv < 0) {
					rv = new_uvar(uid);   /* set up new variable */
				}
				nxtToken.typecode = types::uservar;
				nxtToken.userIx = rv;
			}
		}
	
next:
		if (nxtToken.typecode == types::error)
			return retCode::SYNTAX_ERROR;  /* unable to tokenise expression*/
		else
			tokens.push_back(nxtToken);
	}

	int brackets = 0; /* count depth of brackets*/
	for (auto t : tokens) {
		if (t.typecode == types::Operator && t.oper == opCode::leftb)
			brackets++;
		if (t.typecode == types::Operator && t.oper == opCode::rightb) {
			brackets--;
			if (brackets < 0)
				return retCode::PAREN_MISMATCH;
		}
	}
	if (brackets > 0)
		return retCode::PAREN_MISMATCH;

	nxtToken.typecode = types::end;
	tokens.push_back(nxtToken);
	return retCode::EXPR_OK;
}

/* create a new user variable with name 'name', initialise it
and return its location in the global uvars structure */
static int new_uvar(const char *name) {
	if (uvars.num == uvars.alloc) {
		//need more room for variables
		if (uvars.num == 0) 
			uvars.alloc = 4;
		else 
			uvars.alloc *= 2;

		uvars.vars.resize(uvars.alloc);
	}
	/* copy variable name, then set variable's value to 0 */
	strcpy_s(uvars.vars[uvars.num].name, sizeof(uvars.vars[uvars.num].name),
		name);
	uvars.vars[uvars.num].data = 0;
	uvars.num++;   /* increase number of user variables */
	return uvars.num - 1;
}

//look for 'name' in the global uvars structure
//if found, copy in data and return 0 else return 1
static int set_uvar(const char *name, const Znum &data) {
	int i;

	for (i = 0; i < uvars.num; i++) {
		if (strcmp(uvars.vars[i].name, name) == 0) {
			uvars.vars[i].data =  data;
			return 0;
		}
	}  

	return 1;  /* name not found */
}

/* look for 'name' in the global uvars structure
   if found, copy out data and return index else return -1 if not found */
static int get_uvar(const char *name, Znum data)
{
	int i;

	for (i = 0; i < uvars.num; i++) {
		if (strcmp(uvars.vars[i].name, name) == 0) 	{
			data=  uvars.vars[i].data;
			return i;
		}
	}

	return -1;  /* name not found */
}

static void free_uvars() {
	uvars.vars.clear();           /* reset size to 0 */
	uvars.vars.shrink_to_fit();   /* try to release memory */
	uvars.num = 0;  
	uvars.alloc = 0;
}
void printvars(std::string name) {
	while (name.size() > 0 && isblank(name[0])) {
		name.erase(0, 1);  /* remove leading spaces */
	}
	if (uvars.num == 0)
		printf_s("No user variables defined \n");
	else
		for (int i = 0; i < uvars.num; i++)
			if (name.size() == 0 || strcmp(uvars.vars[i].name, name.c_str()) == 0)
			gmp_printf("%-16s  %Zd\n", uvars.vars[i].name, uvars.vars[i].data);
 }