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

#include "pch.h"
#include "bignbr.h"
#define WIN32_LEAN_AND_MEAN
#include <Windows.h>
#include <Mmsystem.h >   // for sound effects
#include "factor.h"
#include <stack>

//#define BIGNBR       // define to include bignbr tests 
#ifdef BIGNBR
#include "bignbr.h"
#include "bigint.h"
extern Znum zR, zR2, zNI, zN;
#endif

/* external function declaration */
void msieveParam(const std::string &expupper);   /*process Msieve commands */
void yafuParam(const std::string &command);      /*process YAFU commands */
void biperm(int n, Znum &result);   // declaration for external function
void VersionInfo(const LPCSTR path, int ver[4], std::string &modified);
char * getFileName(const char *filter, HWND owner);


/* get time in format hh:mm:ss */
const char * myTime(void) {
	static char timestamp[10];   // time in format hh:mm:ss
	struct tm newtime;

	const time_t current = time(NULL);  // time as seconds elapsed since midnight, January 1, 1970
	localtime_s(&newtime, &current);    // convert time to tm structure
	/* convert time to hh:mm:ss */
	strftime(timestamp, sizeof(timestamp), "%H:%M:%S", &newtime);
	return timestamp;
}

int lang = 0;             // 0 English, 1 = Spanish

bool hex = false;		// set true if output is in hex
bool factorFlag = true;
/* verbose value is used to turn off or on optional messages; 
higher value = more messages */
#ifdef _DEBUG
int verbose = 1;
#else
int verbose = 0;
#endif

HANDLE hConsole;

/* list of operators, arranged in order of priority. Order is not exactly the
same as C or Python. */
enum class opCode {
	fact        = 21,	// !   factorial 
	dfact       = 22,	// !!  double factorial
	prim        = 23,	// #   primorial
	unary_minus =  1,   // C and Python put unary minus above multiply, divide & modulus
	not         = 16,   // C and Python put bitwise not with unary minus
	power       =  0, 
	multiply    =  2,
	divide      =  3,
	remainder   =  4,   // AKA modulus
	comb        =  5,   // nCk, also known as binomial coefficient
	plus        =  6,
	minus       =  7,
	shr         =  8,
	shl         =  9,
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
	rightb		= 24,              // right bracket (must be highest value)
};

/* list of operator priority values. lower value = higher priority. Note: this 
order is not the same as C or Python. The priority value is obtained by using 
the opCode (cast to an integer) as the index. */


struct attrs {
	int pri;
	bool left;   /* associativity ; true = left to right, false = right to left 
		    operators with the same precedence must have the same associativity. */
	bool pre;    /* true if unary operator e.g. - precedes expression, false if
				    it follows e.g. !, otherwise not used */
	int numOps;  /* number of operands; 1 = unary, or 2, or 0 for bracket)*/
};

const static attrs opr[] = {
	{2,  false, false, 2},  // 0 power (right to left)
	{1,  false, true,  1},  // 1 unary minus (right to left)
	{3,  true,  false, 2},  // 2 multiply
	{3,  true,  false, 2},  // 3 divide
	{3,  true,  false, 2},  // 4 remainder AKA modulus
	{4,  true,  false, 2},  // 5 combination nCk, also known as binomial coefficient
	{5,  true,  false, 2},  // 6 plus
	{5,  true,  false, 2},  // 7 minus
	{6,  true,  false, 2},  // 8 shift right
	{6,  true,  false, 2},  // 9 shift left
	{7,  true,  false, 2},  // 10 compare less or equal (not greater)
	{7,  true,  false, 2},  // 11 compare greater or equal (not less)
	{7,  true,  false, 2},  // 12 greater
	{7,  true,  false, 2},  // 13 less
	{8,  true,  false, 2},  // 14 not equal
	{8,  true,  false, 2},  // 15 equal
	{1,  false, true,  1},  // 16 NOT (unary operator ,right to left))
	{9,  true,  false, 2},  // 17 AND
	{10, true,  false, 2},  // 18 XOR
	{11, true,  false, 2},  // 19 OR
	{99, true,  false, 0},  // 20 left bracket
	{0,  true,  false, 1},  // 21 ! factorial (unary operator)
	{0,  true,  false, 1},  // 22 !! double factorial (unary operator)
	{0,  true,  false, 1},  // 23 # primorial (unary operator)
	{-1, true,  false, 0},  // 24 right bracket
};

/* error and return codes, errors are -ve, OK is 0, FAIL is +1 */
enum class retCode
{
	//EXPR_NUMBER_TOO_LOW,
	EXPR_NUMBER_TOO_HIGH= -100,
	EXPR_INTERM_TOO_HIGH ,
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

enum class types { Operator, func, number, comma, error, end };

struct token {
	types typecode;
	long long function;   /* contains function code index,  only when typecode = func */
	opCode oper;    /* contains operator value, only when typecode = Operator*/
	Znum value;     /* contains numeric value,  only when typecode = number*/
};


bool *primeFlags = NULL;
unsigned long long *primeList = NULL;
unsigned int prime_list_count = 0;
unsigned long long int primeListMax = 0;
const unsigned long long max_prime = 1000000007;  // arbitrary limit 10^9,


/* function declarations, only for functions that have forward references */
static retCode ComputeExpr(const std::string &expr, Znum &ExpressionResult);
long long MulPrToLong(const Znum &x);
static bool getBit(const unsigned long long int x, const bool array[]);
void generatePrimes(unsigned long long int max_val);
retCode tokenise(const std::string expr, std::vector <token> &tokens);
static void textError(retCode rc);
static retCode evalExpr(const std::vector<token> &rPolish, Znum & result);
static int reversePolish(token expr[], const int exprLen, std::vector<token> &rPolish);
static void printTokens(const token expr[], const int exprLen);

// find leftmost 1 bit of number. Bits are numbered from 63 to 0
// An intrinsic is used that gets the bit number directly.
// If n = zero the value returned is 0, although bit 0 is zero
static int leftBit(unsigned __int64 n) {
	int bit = 0, nz = 1;

	nz = _BitScanReverse64((DWORD*)&bit, n);  // this intrinsic should be expanded inline
			// to use the fastest available machine instructions
	if (nz) return (int)bit;
	else return 0;   // n actually contains no 1-bits.
	
}

/* Convert number to hexdecimal. Ensure that if number is negative the leftmost
bit of the most significant digit is set, and conversely, if the number is positive
the leftmost bit of the most significant digit is not set. This is done by 
prefixing the output with '0' or 'f' when necessary. */
static char* getHex(Znum Bi_Nbr) {
	static char *hexbuffer = NULL;  // IMPORTANT. This must be a static variable!
	std::string obuff;
	if (Bi_Nbr >= 0) {
		hexbuffer = mpz_get_str(NULL, 16, ZT(Bi_Nbr));
		if (!(hexbuffer[0] >= '0' && hexbuffer[0] <= '7')) {
			/* we need to insert a zero. move everything down 1 byte*/
			hexbuffer = (char *)realloc(hexbuffer, strlen(hexbuffer) + 2);
			assert(hexbuffer != NULL);
			memmove(&hexbuffer[1], &hexbuffer[0], strlen(hexbuffer) + 1); 
			hexbuffer[0] = '0';
		}
	}
	else {
		if (Bi_Nbr == -1)
			obuff = 'f';		// -1 has to be a special case		
		else while (Bi_Nbr != -1) {
			Znum q, r;
			long long rll;
			mpz_fdiv_q_2exp(ZT(q), ZT(Bi_Nbr), 4);	// q = Bi_Nbr/16 rounded towards -inf
			mpz_fdiv_r_2exp(ZT(r), ZT(Bi_Nbr), 4);	// r = Bi_Nbr - q*16 (0 <= r <= 15)
			rll = MulPrToLong(r);
			obuff += (rll <= 9) ? '0' + (char)rll : 'a' + (char)rll -10;  // get char '0' to 'f'
			Bi_Nbr = q;
		}

		std::reverse(std::begin(obuff), std::end(obuff));  // get digits into correct order
		obuff += '\0';				// null-terminate string
		if (obuff[0] >= '0' && obuff[0] <= '7') {
			/* we need to insert an f */
			obuff.insert(obuff.begin(), 'f');
		}
		/* copy text from local STL string to global C-style string */
		hexbuffer = (char *)malloc(obuff.size() + 1);
		assert(hexbuffer != NULL);
		strncpy_s(hexbuffer, obuff.size(), obuff.c_str(), _TRUNCATE);
	}
	return hexbuffer;  // hexbuffer must be a static variable
}

/* output value of Nbr as ascii text to stdout */
void ShowLargeNumber(const Znum &Bi_Nbr, int digitsInGroup, bool size, bool hex) {
	std::string nbrOutput = "";
	char* buffer = NULL;
	size_t msglen, index = 0;

	if (!hex) {
		// convert to null-terminated ascii string, base 10, with sign if -ve
		// corrrectly-sized buffer is allocated automatically.
		buffer = mpz_get_str(NULL, 10, ZT(Bi_Nbr));
	}
	else {
		buffer = getHex(Bi_Nbr);
	}
	msglen = strnlen(buffer, 500000);  // arbitrary limit of 500000
	if (buffer[0] == '-') {      // if number is -ve
		nbrOutput = "-";
		index = 1;
	}
	if (hex) nbrOutput += "0x";
	for (; index < msglen; index++) {
		if ((msglen - index) % digitsInGroup == 0) {
			// divide digits into groups, if there are more than 6 digits
			if ((index > 0 && buffer[0] != '-') || (index > 1)) {
				nbrOutput += " ";  // put space after group of digits 
			}
		}
		nbrOutput += buffer[index];
	}
	if (buffer[0] == '-')   // if number is -ve
		msglen--;
	free(buffer);		// avoid memory leakage
	std::cout << nbrOutput;
	if (msglen > 6 && size)
		std::cout << " (" << msglen << " digits)";
}

/* convert biginteger to normal. Checks for overflow */
long long MulPrToLong(const Znum &x) {
	long long rv;
	// note: do not use mpz_fits_slong_p because it checks whether x fits a 32 bit integer, rather than a 64 bit integer.
	if (x >= LLONG_MIN && x <= LLONG_MAX) { // is x value OK for normal integer?
		//rv = mpz_get_si(ZT(x)); // convert to normal integer
		rv = ZT(x)->_mp_d[0];     // accessing the limb directly seems to be a lot faster
		if (ZT(x)->_mp_size < 0)  // than calling mpz_get_si
			rv = -rv;
		if (ZT(x)->_mp_size == 0)
			rv = 0;
		return rv;
	}
	else
		throw std::range_error("big number cannot be converted to 64-bit integer");
	return 0;
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

	if (n == 1)
		return 1;  // 1 only has one divisor. NoOfDivs can't handle that case
	auto rv = factorise(n, factorlist, nullptr);
	if (rv) {
		auto divisors = factorlist.NoOfDivs();
		return divisors;
	}
	else return 0;
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

	Znum argum= n, Temp;
	Znum result = 0;

	while (argum > 0)
	{
		Temp = argum%radix;
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
		Temp = argum%radix;
		result += Temp;
		argum /= radix;
	}
	return result;
}

/* NumDigits(n,r): Number of digits of n in base r. Leading zeros are not counted  */
long long ComputeNumDigits(const Znum &n, const Znum &radix)
{
	Znum result = n;
	long long digits = 0;
	
	while (result > 0) {
		result /= radix;
		digits++;
	}
	return digits;
}

/* get 'width' of n in bits */
static long long NoOfBits(const Znum &n) {
	auto result = mpz_sizeinbase(ZT(n), 2);  // calculate number of bits
	return result;
}

/* IsPrime(n): returns zero if n is definately composite,  -1 if it is a probable prime, */
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
		return retCode::EXPR_INVALID_PARAM;  // 2 is the smallest prime
	if (n == 3) {
		p = 2;
		return retCode::EXPR_OK;
	}
	Znum i;		// largest odd number below n
	i = ((n & 1) == 1) ? n - 2 : n - 1;  // i >= 3

	for (; i > 0; i-=2) {
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
	for (; i >=2; i -= 2) {
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
	std::string result;
	Znum rvalue = 0;
	const bool descending = ((mode & 1) == 1);
	const bool repeat = ((mode & 2) == 2);
	char *buffer = NULL;

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
static void squareFree(Znum &num, Znum &sq, std::vector<zFactors> &sqf) {
	fList factorlist;
	zFactors temp;
	auto rv = factorise(num, factorlist, nullptr);

	factorlist.sqfree(sqf);  // get factors of square in sqf
	
	sq = 1;
	Znum x;
	for (auto i:  sqf) {  /* calculate value of square */
		mpz_pow_ui(ZT(x), ZT(i.Factor), i.exponent);
		sq *= x;
	}
	num /= sq;  // adjust value of num
}


extern unsigned __int64 R3(__int64 n);
/* calculate the number of ways an integer n can be expressed as the sum of 3
squares x^2, y^2 and z^2. The order of the squares is significant. x, y and z can
be +ve, 0 or -ve See https://oeis.org/A005875 */
static Znum R3(Znum num) {
	
	if (num < 200000000000000) {
		__int64 llnum = MulPrToLong(num);
		return R3(llnum);
	}
	Znum sum = 0, sq, multiplier = 1;
	std::vector <zFactors> sqf;   // factor list for square factor of num

	if (num < 0)
		return 0;
	if (num == 0)     // test here necessary to avoid infinite loop
		return 1;
	if (num % 8 == 7)
		return 0;     // take short cut if possible
	while ((num & 3) == 0)
		num >>= 2;      // remove even factors. note that R3(4n) =R3(n)
	squareFree(num, sq, sqf);

	for (Znum k = 1; k*k <= num; k++) {
		sum += 2* R2(num - k * k);
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

long long lltTdivCnt = 0;  /* count no of Mersenne numbers (partly) factored by
						   trial division */
long long lltCmpCnt = 0;   /* count no of Mersenne numbers determined to be composite
						   using Lucas-Lehmer test */
long long lltPriCnt = 0;   /* count no of Mersenne numbers determined to be prime
						   using Lucas-Lehmer test */

/* n = 2^p -1. Increment lltTdivCnt & return 0 if trial division finds a factor, 
   otherwise return 1 */
static int lltTrialDiv(const Znum &n, const long long p) {
	Znum tmp, ncopy, limit;
	const int maxlimit = 1000000;  /* Limits the number of iterations for trial
									 division. Experiments indicate that increasing
									 this limit slows llt down overall. */
	bool composite = false;

	/* do trial division */
	ncopy = n;

	/*  n is a mersenne number calculated as (2^p)-1 and p is a prime number.
		The factors of n must be in the form of q = 2ip+1, where q < n.
		2ip+1 < n therefore
		i < (n-1)/2p */
		//limit = (sqrt(n) - 1) / (2 * p) + 1;   /* mpz_sqrt is probably faster */
	Znum nm1 = n - 1;
	mpz_sqrt(ZT(limit), ZT(nm1));
	/* if we find all factors < sqrt(n) the residue is the one remaining factor */

	if (limit > maxlimit) {
		/* hit this limit if p >= 41 */
		if (verbose > 1)
			gmp_printf("limit reduced from %Zd to %lld\n", limit, maxlimit);
		limit = maxlimit;   // larger limit would take too long (also encounter > 64-bit integers)
	}

	long long ll_limit = MulPrToLong(limit);
	Znum rem;
	for (long long i = 1; i <= ll_limit; i++) {
		long long d = 2 * i*p + 1;
		//if (n%d == 0)
		long long r = mpz_tdiv_r_ui(ZT(rem), ZT(n), d);
		if (r == 0) {
			composite = true;
			ncopy /= d;
			if (verbose > 0)
				printf_s("2*%lld*p+1 (= %lld) is a factor\n", i, d);
			else {
				lltTdivCnt++;  /* increment counter */
				return 0;  /* not prime; return as soon as we know this */
			}
		}
	}

	if (composite) {
		if (ncopy > 1 && verbose > 0)
			gmp_printf("%Zd is a factor\n", ncopy);
		lltTdivCnt++;   /* increment counter */
		return 0; /* not prime */
	}
	else {
		if (verbose > 0)
		printf_s("%s No factors found by trial division to %d bits, start Lucas-Lehmer test \n",
			myTime(), leftBit(2 * ll_limit * p + 1) + 1);
		return 1;  /* no factor found */
	}
}

/* get remainder when num is divided by m = 2^p-1. The obscure method used 
avoids division by 2^p-1, and divides by 2^p instead which is much faster,
as it is equivalent to a right shift by p bits. 
For efficiency both p and m are given as parameter values although a cleaner
interface would just supply p and let m be recalculated. */
static void getrem(Znum &rem, const Znum &num, const long long p, const Znum &m) {
	static Znum mod2p, q;  /* static for efficiency */
	mod2p = num & m;     // mod2p = num%(2^p) = num -q*2^p
	mpz_tdiv_q_2exp(ZT(q), ZT(num), p);  // q = num/(2^p)
	rem = q + mod2p;                    //rem = num-q*(2^p-1)
	while (rem >= m)
		rem -= m;               // reduce rem if necessary until 0 <= rem < (2^p-1)
#ifdef _DEBUG
	/* demonstrate that the result is correct */
	//Znum  temp;
	//mpz_tdiv_r(ZT(temp), ZT(num), ZT(m));
	//assert(temp == rem);
#endif
}

/* perform Lucas-Lehmer test. Return 0 if 2^p-1 is composite, 1 if prime  
see https://en.wikipedia.org/wiki/Lucas%E2%80%93Lehmer_primality_test 
also https://oeis.org/A000043 
and https://www.mersenne.org/primes/ 
the largest mersenne prime found so far is p=82589933 
p must not be too large, otherwise it would take centuries to get a result.
This code is based on the llt function in YAFU, but has been modified to make 
it much faster. 
if verbose is > 0 some extra messages are output.
If verbose is > 1 trial division is used to try to find factors */
static Znum llt(const Znum &p) {
	int rv;
	Znum n, tmp; 
	long long exp = MulPrToLong(p); /* exp = p */
	clock_t t1, t2, t3;

	/* 1st check if p is prime. Maybe better to use mpz_bpsw_prp instead?*/
#ifdef __MPIR_VERSION
	static gmp_randstate_t state;
	static bool first = true;
	if (first) {
		gmp_randinit_default(state);
		first = false;
	}
	rv = mpz_probable_prime_p(ZT(p), state, 16, 0);
#else
	rv = mpz_probab_prime_p(ZT(p), 16); /* returns 0, 1 or 2*/
#endif 
	/* rv is 1 if p is probably prime, or 0 if p is definitely composite.*/
	if (rv == 0) {
		lltTdivCnt++;  // count this result as found by trial division
		return 0;      // if p is composite, 2^p -1 is composite.
	}
	if (p <= 7) {
		lltPriCnt++;  // count as found by Lucas-Lehmer test
		return 1;    /* 1st mersenne number that is not prime is p=11 */
	}

	n = 1;
	mpz_mul_2exp(ZT(n), ZT(n), exp);   // n = 2^p
	n -= 1;                            // n = 2^p-1

	/* experiment - is mpz-probab_prime_p faster overall than trial division? 
	 Results indicate yes, definitely. Maybe because it is >99% effective in screening
	 out composites with just 2 tests so that llt is only neeeded as confirmation. 
	 Update 11/7/2021; using neither trial division nor  mpz-probab_prime_p is fastest.
	 LL test is in fact faster than miller-rabin even when MR is limitied to 2 rounds. */
	if (verbose > 1) {
		t1 = clock();
		/* use trial division if verbose > 1 because this can give some factors, 
		whereas the other tests only show whether n is prime or composite */
		int rv = lltTrialDiv(n, exp);
		t2 = clock();
		printf_s("time used by trial division = %.2f sec \n", (double)(t2 - t1) / CLOCKS_PER_SEC);
		if (rv == 0) {
			std::cout << myTime() << " 2^" << exp << " -1 is not prime *** \n";
			return 0;   /* factor found -> n is composite */
		}
	}
	
	/*do the ll test */
	tmp = 4;
	int nchars = 0;
	if (verbose > 0) 
		t2 = clock();
	for (long long i = 0; i < exp - 2; i++) {
		tmp *= tmp;
		tmp -= 2;
		//tmp %= n;
		getrem(tmp, tmp, exp, n);  /* get remainder of division by n */
		if (verbose > 0 && (i > 0) && (i & 0x7ff) == 0) { /* 7ff = 2047, print msg every 2048 iterations */
			for (int j = 0; j < nchars; j++)
				printf("\b");    // erase previous output
			nchars = printf("%s llt iteration %lld %.2f%% complete", myTime(), i,
				(double)i*100.0/(exp-2) );
			fflush(stdout);
		}
	}
	if (verbose > 0) { 
		t3 = clock();
		printf_s("\ntime used by llt = %.2f sec \n", (double)(t3 - t2) / CLOCKS_PER_SEC);
	}

	if (tmp == 0) {
		lltPriCnt++;  /* increment counter */
		return 1;             // prime
	}
	else {
		lltCmpCnt++;   /* increment counter */
		return 0;            // composite
	}
}

enum class fn_Code {
	fn_gcd,
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
	fn_np,
	fn_pp,
	fn_r2,
	fn_r3,
	fn_legendre,
	fn_jacobi,
	fn_kronecker,
	fn_llt,
	fn_sqrt,
	fn_nroot,
	fn_bpsw,
	fn_aprcl,
	fn_numfact,
	fn_minfact,
	fn_maxfact,
	fn_ispow,
	fn_invalid = -1,
} ;

struct  functions {
	char fname[11];        // maximum name length 10 chars (allow for null terminator)
	int  NoOfParams;       // number of parameters 
	fn_Code  fCode;        // integer code for function
};

/* list of function names. No function name can begin with C because this would 
 conflict with the C operator. Longer names must come before short ones 
 that start with the same letters to avoid mismatches */
const static std::array <struct functions, 31> functionList{
	"GCD",       2,  fn_Code::fn_gcd,			// name, number of parameters, code
	"MODPOW",    3,  fn_Code::fn_modpow,
	"MODINV",    2,  fn_Code::fn_modinv,
	"TOTIENT",   1,  fn_Code::fn_totient,
	"SUMDIVS",   1,  fn_Code::fn_sumdivs,
	"SUMDIGITS", 2,  fn_Code::fn_sumdigits,
	"SQRT",      1,  fn_Code::fn_sqrt,
	"NUMDIGITS", 2,  fn_Code::fn_numdigits,
	"NUMDIVS",   1,  fn_Code::fn_numdivs,
	"NROOT",     2,  fn_Code::fn_nroot,
	"REVDIGITS", 2,  fn_Code::fn_revdigits,
	"ISPRIME",   1,	 fn_Code::fn_isprime,
	"NUMFACT",   1,  fn_Code::fn_numfact,
	"MINFACT",   1,  fn_Code::fn_minfact,
	"MAXFACT",   1,  fn_Code::fn_maxfact,
	"FactConcat",2,  fn_Code::fn_concatfact,     // FactConcat must come before F
	"F",         1,  fn_Code::fn_fib,			// fibonacci
	"LLT",	     1,  fn_Code::fn_llt,           // lucas-Lehmer test
	"LE",		 2,  fn_Code::fn_legendre,
	"L",         1,  fn_Code::fn_luc,			// Lucas Number
	"PI",		 1,  fn_Code::fn_primePi,		// prime-counting function. PI must come before P
	"P",         1,  fn_Code::fn_part,			// number of partitions
	"N",         1,  fn_Code::fn_np,			// next prime
	"BPSW",      1,  fn_Code::fn_bpsw,          // Baillie-Pomerance-Selfridge-Wagstaff
	"B",         1,  fn_Code::fn_pp,			// previous prime
	"R2",		 1,  fn_Code::fn_r2,			// number of ways n can be expressed as sum of 2 primes
	"R3",        1,  fn_Code::fn_r3,
	"JA",		 2,  fn_Code:: fn_jacobi,
	"KR",		 2,  fn_Code::fn_kronecker,
	"APRCL",     1,  fn_Code::fn_aprcl,          // APR-CL prime test
	"ISPOW",     1,  fn_Code::fn_ispow,
};

/* Do any further checks needed on the parameter values, then evaluate the function. 
Only ModPow uses all 3 parameters. Some functions can generate error codes. */
static retCode ComputeFunc(fn_Code fcode, const Znum &p1, const Znum &p2, 
	const Znum &p3, Znum &result) {
	int rv;
	Znum temp;
	retCode retcode = retCode::EXPR_OK;

	/* evaluate function value using parameter values  */
	switch (fcode) {
	case fn_Code::fn_gcd: {			// GCD	
		//mpz_gcd(ZT(result), ZT(p1), ZT(p2));
		result = gcd(p1, p2);
		break;
	}
	case fn_Code::fn_modpow: {						// MODPOW
		if (p3 == 0)
			return retCode::EXPR_DIVIDE_BY_ZERO;
		if (p2 < 0) {
			if (gcd(p1, p3) != 1)
				return retCode::EXPR_EXPONENT_NEGATIVE;  // p1 and p3 not mutually prime
		}
		/* note: negative exponent is only allowed if p1 & p3 are mutually prime,
		i.e modular inverse of p1 wrt p3 exists. */
		mpz_powm(ZT(result), ZT(p1), ZT(p2), ZT(p3));
		break;
	}
	case fn_Code::fn_modinv: {						// MODINV
		/* if an inverse doesn’t exist the return value is zero and rop is undefined*/
		rv = mpz_invert(ZT(result), ZT(p1), ZT(p2));
		if (rv == 0) {
			return retCode::EXPR_ARGUMENTS_NOT_RELATIVELY_PRIME;
		}
		break;
	}

	case fn_Code::fn_totient: {			// totient
		if (p1 < 1) return retCode::EXPR_INVALID_PARAM;
		result = ComputeTotient(p1);
		break;
	}
	case fn_Code::fn_numdivs: {		// NUMDIVS
		if (p1 < 1) {
			return retCode::EXPR_INVALID_PARAM;
		}
		result = ComputeNumDivs(p1);
		break;
	}
	case fn_Code::fn_sumdivs: {		// SUMDIVS
		result = ComputeSumDivs(p1);
		break;
	}

	case fn_Code::fn_sumdigits: {		// SumDigits(n, r) : Sum of digits of n in base r.
		result = ComputeSumDigits(p1, p2);
		break;
	}
	case fn_Code::fn_numdigits: {		// numdigits
		result = ComputeNumDigits(p1, p2);
		break;
	}
	case fn_Code::fn_revdigits: {	// revdigits
		result = ComputeRevDigits(p1, p2);
		break;
	}

	case fn_Code::fn_isprime: {  // isprime
						/* -1 indicates probably prime, 0 = composite */
		result = PrimalityTest(abs(p1));
		break;
	}
	case fn_Code::fn_fib: {		// fibonacci
		if (p1 < 0) {
			return retCode::EXPR_INVALID_PARAM;
		}
		if (p1 > 95700)
		{
			return retCode::EXPR_INTERM_TOO_HIGH;
		}
		long long temp = MulPrToLong(ZT(p1));
		mpz_fib_ui(ZT(result), temp);  // calculate fibonacci number
		break;
	}
	case fn_Code::fn_luc: {		// lucas number
		if (p1 < 0) {
			return retCode::EXPR_INVALID_PARAM;
		}
		if (p1 > 95700)
		{
			return retCode::EXPR_INTERM_TOO_HIGH;
		}
		long long temp = MulPrToLong(p1);
		mpz_lucnum_ui(ZT(result), temp);  // calculate lucas number
		break;
	}

	case fn_Code::fn_part: {              // number of partitions
		if (p1 < 0 || p1 > 1000000) {
			return retCode::EXPR_INVALID_PARAM;  
			// note: biperm is limited to values <= 1,000,000
		}
		int temp = (int)MulPrToLong(p1);
		biperm(temp, result);   // calculate number of partitions
		break;
	}
	case fn_Code::fn_np: {  // next prime;
		mpz_nextprime(ZT(result), ZT(p1));  // get next prime
		break;
	}
	case fn_Code::fn_pp: {			// previous prime
		retcode = ComputeBack(p1, result);  // get previous prime
		if (retcode != retCode::EXPR_OK) {
			return retcode;   // error: number < 3
		}
		break;
	}

	case fn_Code::fn_primePi: {  // count primes <= n
		if (p1 > max_prime) {
			return retCode::EXPR_INVALID_PARAM;
		}
		result = primePi(p1);
		break;
	}
	case fn_Code::fn_concatfact: { /*Concatenates the prime factors of n according to
						  the mode in m */
		if (p1 < 0 || p1 > 3) {
			return retCode::EXPR_INVALID_PARAM;  // mode value invalid
		}
		result = FactConcat(p1, p2);
		break;
	}
	case fn_Code::fn_r2: {
		result = R2(p1);
		break;
	}
	case fn_Code::fn_r3: {
		result = R3(p1);
		break;
	}
	case fn_Code::fn_legendre: {
		/* p2 must be an odd positive prime */
		result = mpz_legendre(ZT(p1), ZT(p2));
		break;
	}
	case fn_Code::fn_jacobi: {
		/*p2 must be odd */
		result = mpz_jacobi(ZT(p1), ZT(p2));
		break;
	}
	case fn_Code::fn_kronecker: {
		result = mpz_kronecker(ZT(p1), ZT(p2));
		break;
	}
	case fn_Code::fn_llt: {
		/* see https://en.wikipedia.org/wiki/Lucas%E2%80%93Lehmer_primality_test */
		if (p1 >= 0 && p1 <= INT_MAX) {
			result = llt(p1);
			if (verbose > 0 || p1 >= 216091) {
				if (result == 1)
					std::cout << "*** 2^" << p1 << " -1 is prime! ***\n";
				else
					std::cout << "2^" << p1 << " -1 is not prime \n";
			}
		}
		 else /* for large numbers LLT takes a very long time!! */
			 return retCode::EXPR_NUMBER_TOO_HIGH;
		break;
	}
	case fn_Code::fn_sqrt: {
		if (p1 < 0) return retCode::EXPR_INVALID_PARAM;
		mpz_sqrt(ZT(result), ZT(p1));  /* result = square root of p1*/
		break;
	}
	case fn_Code::fn_nroot: {
		/* for real numbers nroot (x, n)  = x^(1/n) has a discontinuity at n=0, 
		so the function is considered to be undefined for n=0 */
		if (p2 == 0) return retCode::EXPR_INVALID_PARAM;

		/* odd root of a -ve number is -ve so allow it. Even root would be a complex
		number so don't allow it */
		if (p1 < 0 && mpz_even_p(ZT(p2)))
			return retCode::EXPR_INVALID_PARAM;

		if (p2 < 0) {
			result = 0; /* for real numbers nroot(x, n) with n -ve would return
						a number between 0 and 1. The floor of this value is 0 */
			break; 
		}
		if (p2 > INT_MAX) {
			/* if p2 is very large the nth root would be close to 1  for +ve p1 
			or -1 for -ve p1 and odd p2 */
			return retCode::EXPR_NUMBER_TOO_HIGH;
		}
		
		long long temp = MulPrToLong(ZT(p2));
		mpz_root(ZT(result), ZT(p1), temp);  /* result = nth root of p1*/
		break;
	}
	case fn_Code::fn_bpsw: {
		if (p1 <= 1) return retCode::EXPR_INVALID_PARAM;
		/* Baillie-Pomerance-Selfridge-Wagstaff probablistic primality test */
		result = mpz_bpsw_prp(ZT(p1));
		if (result == 0)
			std::cout << "composite \n";
		else if (result == 1)
			std::cout << "probable prime \n";
		else if (result == 2)
			std::cout << "prime \n";
		break;
	}
	case fn_Code::fn_aprcl: {
		if (p1 <= 1) return retCode::EXPR_INVALID_PARAM;
		result = mpz_aprtcle(ZT(p1), verbose);
		if (result == 0)
			printf_s ( "composite \n");
		else if (result == 1)
			printf_s("probable prime \n");
		else if (result == 2)
			printf_s("prime \n") ;
		break;
	}
	case fn_Code::fn_numfact: {
		result = ComputeNumFact(p1);
		break;
	}
	case fn_Code::fn_minfact: {
		result = ComputeMinFact(p1);
		break;
	}
	case fn_Code::fn_maxfact: {
		result = ComputeMaxFact(p1);
		break;
	}
	case fn_Code::fn_ispow: {
		/* return -1 if p1 is a perfect power, otherwise 0 */
			long long MaxP = 393'203;  // use 1st  33333 primes
			if ((long long)primeListMax < MaxP) {  // get primes
				generatePrimes(MaxP);  // takes a while, but only needed on 1st call
			}
			Znum base;
			long long exp = PowerCheck(p1, base, 2);
			if (exp > 1) {
				std::cout << p1 << " = " << base << "^" << exp << '\n';
				result = -1;
			}
			else {
				result = 0;
			}
			break;
		}

	default:
		abort();		// if we ever get here we have a problem
	}
	return retcode;
}

/* SHL: Shift left the number of bits specified on the right operand. If the
number of bits to shift is negative, this is actually a right shift. If result would
be too large an error is reported. */
static retCode ShiftLeft(const Znum &first, const Znum &bits, Znum &result) {
	if (bits > LLONG_MAX || bits < LLONG_MIN)
		return retCode::EXPR_INVALID_PARAM;

	long long shift = MulPrToLong(bits);

	if (shift == 0) {
		result = first;
		return retCode::EXPR_OK;  // shift zero bits; nothing to do
	}

	// there is no built-in shift operator for Znums, so it is simulated
	// using multiplication or division
	if (shift > 0) {
		if (NoOfBits(first) + shift > 66439) // more than 66439 bits -> more than 20,000 decimal digits
			return retCode::EXPR_INTERM_TOO_HIGH;
		mpz_mul_2exp(ZT(result), ZT(first), shift);
		return retCode::EXPR_OK;
	}
	else {   // note use of floor division
		mpz_fdiv_q_2exp(ZT(result), ZT(first), -shift);
		return retCode::EXPR_OK;
	}
}


/* get floor(sqrt(n))*/
unsigned long long llSqrt(const unsigned long long n) {
	double root = sqrt((double)n);
	unsigned long long  iroot = llround(n);
	while (iroot*iroot > n)
		iroot--;
	return iroot;
}

// return value of bit in array corresponding to x. false (0) indicates a prime number
static bool getBit(const unsigned long long int x, const bool array[])
{
	return array[x/2];
}

/* sets bit in 'array' corresponding to 'x' */
/* assume 'num' is odd. no bits are stored for even numbers */
static void setBit(const unsigned long long int x, bool array[]) {
	array[x / 2] = true;
	return;
}

void generatePrimes(unsigned long long int max_val) {
	unsigned long long int numsave=0, count = 1, num = 3;
	size_t plist_size = 0;  // size of prime list
	unsigned long long int sqrt_max_val;

	if (primeListMax >= max_val)
		return;  // nothing to do if prime list already generated

	max_val += 63 - (max_val + 63) % 64;  // round up to next multiple of 64
	sqrt_max_val = llSqrt(max_val);

	/* initialise flags */
	
	if (primeFlags != NULL) delete []primeFlags;	// if generatePrimes has been called before
												// clear out old prime list
	//primeFlags = (uint32_t*)calloc((max_val / 16) + 1, 1);
	primeFlags = new bool[max_val / 2 + 1]{ 0 };
	assert(primeFlags != NULL);
	
	// allocate storage for primeList if required
	{
		if (verbose > 0)
			fprintf(stdout, "Expected no of primes is %.0f\n",
				(double)max_val / (log((double)max_val) - 1));
		if (primeList != NULL) free(primeList);
		plist_size = (size_t)((double)max_val / (log((double)max_val) - 1)) * 102 / 100;
		// add 2% for safety
		primeList = (unsigned long long *)malloc(plist_size * sizeof(long long));
		assert(primeList != NULL);
		prime_list_count = 1;
		primeList[0] = 2;  // put 1st prime into list
	}

	while (num <= max_val)
	{
		if (!getBit(num, primeFlags))
		{   /* we have found a prime number */
			primeList[count] = num;
			count++;
			numsave = num;

			if (num <= sqrt_max_val)  /* check whether starting value for i is already
									  more than max_val, while avoiding an integer overflow*/
				for (unsigned long long int i = num*num; i <= max_val; i += (num * 2))
				{ /* mark all odd multiples of num in range as not prime */
				  /* starting value is num*num rather than num*3 because the flag bits for
				  multiples of num less than num*num are already set and using num*num
				  makes the program run significantly faster */
					setBit(i, primeFlags);
				}
		}
		num += 2;	// advance to next possible prime
	}

	// after completing the for loop we have found all the primes < max_val
	if (verbose > 0)
		printf("  prime %9lld is %11lld\n", count, numsave);
	primeList[count] = ULLONG_MAX;		// set end marker
	prime_list_count = (unsigned int)count;
	primeListMax = primeList[count - 1];
	return;
}

/* process one operator with 1 or 2 operands. 
NOT, unary minus, factorial, double factorial and primorial  have 1 operand. 
All the others have two. Some operators can genererate an error condition 
e.g. EXPR_DIVIDE_BY_ZERO otherwise return EXPR_OK. */
static retCode ComputeSubExpr(const opCode stackOper, const Znum &firstArg,
	const Znum &secondArg, Znum &result) {
	
	switch (stackOper) 	{

	case opCode::comb: {  // calculate nCk AKA binomial coefficient
		if (secondArg > INT_MAX)
		return retCode::EXPR_NUMBER_TOO_HIGH;
		if (secondArg < 1)
			return retCode::EXPR_INVALID_PARAM;
		long long k = MulPrToLong(secondArg);
		mpz_bin_ui(ZT(result), ZT(firstArg), k);
		return retCode::EXPR_OK;
	}
	case opCode::plus: {
		result = firstArg + secondArg; 
		return retCode::EXPR_OK;
	}
	case opCode::minus: {
		result = firstArg - secondArg; 
		return retCode::EXPR_OK;
	}
	case opCode::unary_minus: {
		result = -firstArg; 
		return retCode::EXPR_OK;
	}
	case opCode::divide: {
		if (secondArg == 0)
			return retCode::EXPR_DIVIDE_BY_ZERO;  // result would be infinity
		result = firstArg / secondArg; 
		return retCode::EXPR_OK;
	}
	case opCode::multiply: {
		auto resultsize = NoOfBits(firstArg) + NoOfBits(secondArg);
		if (resultsize > 99960)  // more than 99960 bits -> more than 30,000 decimal digits
			return retCode::EXPR_INTERM_TOO_HIGH;
		result = firstArg * secondArg;
		return retCode::EXPR_OK;
	}
	case opCode::remainder: {
		if (secondArg == 0)
			return retCode::EXPR_DIVIDE_BY_ZERO;  // result would be infinity
			result = firstArg % secondArg;   
		return retCode::EXPR_OK;
	}
	case opCode::power: {
		if (secondArg > INT_MAX)
			return retCode::EXPR_EXPONENT_TOO_LARGE;
		if (secondArg < 0)
			return retCode::EXPR_EXPONENT_NEGATIVE;
		long long exp = MulPrToLong(secondArg);
		auto resultsize = (NoOfBits(firstArg)-1)* exp;  // estimate number of bits for result
		if (resultsize > 99960)  // more than 99960 bits -> more than 30,000 decimal digits
			return retCode::EXPR_INTERM_TOO_HIGH;
		mpz_pow_ui(ZT(result), ZT(firstArg), exp);
		return retCode::EXPR_OK;
	}
	case opCode::equal: {
		if (firstArg == secondArg)
			result = -1;
		else
			result = 0;
		return retCode::EXPR_OK;
	}
	case opCode::not_equal: {
		if (firstArg != secondArg)
			result = -1;
		else
			result = 0;
		return retCode::EXPR_OK;
	}
	case opCode::greater: {
		if (firstArg > secondArg)
			result = -1;
		else
			result = 0;
		return retCode::EXPR_OK;
	}
	case opCode::not_greater: {
		if (firstArg <= secondArg)
			result = -1;
		else
			result = 0;
		return retCode::EXPR_OK;
	}
	case opCode::less: {
		if (firstArg < secondArg)
			result = -1;
		else
			result = 0;
		return retCode::EXPR_OK;
	}
	case opCode::not_less: {
		if (firstArg <= secondArg)
			result = -1;
		else
			result = 0;
		return retCode::EXPR_OK;
	}
	case opCode::shl: {
		return ShiftLeft(firstArg, secondArg, result);
	}
	case opCode::shr: {
		// invert sign of shift
		return ShiftLeft(firstArg, -secondArg, result);
	}
	case opCode::not: {   // Perform binary NOT 
		//result = -1 - firstArg;  // assumes 2s complement binary numbers
		mpz_com(ZT(result), ZT(firstArg));
		return retCode::EXPR_OK;
	}
	case opCode::and: {  // Perform binary AND.
		mpz_and(ZT(result), ZT(firstArg), ZT(secondArg));
		return retCode::EXPR_OK;
	}
	case opCode::or: {   // Perform binary OR.
		mpz_ior(ZT(result), ZT(firstArg), ZT(secondArg));
		return retCode::EXPR_OK;
	}
	case opCode::xor: {   // Perform binary XOR.
		mpz_xor(ZT(result), ZT(firstArg), ZT(secondArg));
		return retCode::EXPR_OK;
	}
	case opCode::fact: {
		if(firstArg > 5984)
			return retCode::EXPR_INTERM_TOO_HIGH;
		if (firstArg < 0)
			return retCode::EXPR_INVALID_PARAM;
		long long temp = llabs(MulPrToLong(firstArg));
		mpz_fac_ui(ZT(result), temp);  // get factorial
		return retCode::EXPR_OK;
	}
	case opCode::dfact: {
		if (firstArg > 11081)
			return retCode::EXPR_INTERM_TOO_HIGH;
		if (firstArg < 0)
			return retCode::EXPR_INVALID_PARAM;
		long long temp = llabs(MulPrToLong(firstArg));
		mpz_2fac_ui(ZT(result), temp);  // get double factorial
		return retCode::EXPR_OK;
	}
	case opCode::prim: {
		if (firstArg > 46340)
			return retCode::EXPR_INTERM_TOO_HIGH;
		if (firstArg < 0)
			return retCode::EXPR_INVALID_PARAM;
		long long temp = llabs(MulPrToLong(firstArg));
		mpz_primorial_ui(ZT(result), temp);  // get primorial
		return retCode::EXPR_OK;
	}

	default:
		abort();	// should never get here
	}
}

struct oper_list{
	char oper[4];
	opCode operCode;
	int operPrio;
} ;
/* list of operators with corresponding codes and priority. For the search to 
work correctly 
!! must precede !, 
** must precede *, 
<< and <= must precede < 
>> and >= must precede > in this list.
unary minus and unary plus are NOT in this list. */
const static struct oper_list operators[]  {
		{"C",	 opCode::comb,	      3},
		{ "^",   opCode::power,       0},
		{ "**",  opCode::power,       0},     // can use ^ or ** for exponent
		{ "*",   opCode::multiply,    2},
		{ "/",   opCode::divide,      2},
		{ "%",   opCode::remainder,   2},
		{ "+",   opCode::plus,        4},
		{ "-",   opCode::minus,       4},
		{ "SHL", opCode::shl,         5},
		{ "<<",  opCode::shl,         5},     // can use << or SHL for left shift
		{ "SHR", opCode::shr,         5},
		{ ">>",  opCode::shr,         5},     // can use SHR or >> for right shift
		{ "<=",  opCode::not_greater, 6},
		{ ">=",  opCode::not_less,    6},
		{ ">",   opCode::greater,     6},	  // to avoid mismatches > and < must come after >> and <<
		{ "<",   opCode::less,        6},
		{ "!=",  opCode::not_equal,   7},
		{ "==",  opCode::equal,       7},
		{ "NOT", opCode::not,         1},      // bitwise NOT
		{ "AND", opCode::and,         9},      // bitwise AND
		{ "OR",  opCode::or,         11},      // bitwise OR
		{ "XOR", opCode::xor,        10},      // bitwise exclusive or
		{ "!!",  opCode::dfact,       1},      // double factorial
		{ "!",   opCode::fact,        1},      // factorial
		{ "#",   opCode::prim,        1},      // primorial
		{ "(",   opCode::leftb,      12},      // left bracket
		{ ")",   opCode::rightb,      0},      // left bracket
	 } ;

/* search for operators e.g. '+', '*'  that have format <expression> <operator> <expression>
or <operator> <expression>. return index of operator or -1 */
static void operSearch(const std::string &expr, int &opcode) {
	opcode = -1;
	for (int ix = 0; ix < sizeof(operators)/sizeof(operators[0]); ix++) {
		/* use case-insensitive compare */
		if (_strnicmp(expr.data(), operators[ix].oper, strlen(operators[ix].oper)) == 0) {
			opcode = ix;
			return;
		}
	}
	opcode = -1;  // does not match any operator in list
	return;
}

/* evaluate expression and return value in ExpressionResult
operators and functions include:
symbol		meaning					operator Priority
n!			factorial               0
n!!			double factorial        0
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
static retCode ComputeExpr(const std::string &expr, Znum &Result) {
	retCode rv;
	std::vector <token> tokens;
	std::vector <token> rPolish;

	rv = tokenise(expr, tokens); /* 'tokenise' the expression */
	if (rv == retCode::EXPR_OK) {
		rPolish.clear();
		int ircode = reversePolish(tokens.data(), (int)tokens.size(), rPolish); 
		/* convert expression to reverse polish */
		if (ircode == EXIT_SUCCESS) {
			rv = evalExpr(rPolish, Result);
			if (rv != retCode::EXPR_OK && verbose > 0) {
				std::cout << "Expression could not be evaluated \n";
				printTokens(rPolish.data(), (int)rPolish.size());
			}
		}
		else {
			if (verbose > 0) {
				std::cout << "** error: could not convert to reverse polish \n";
				printTokens(rPolish.data(), (int)rPolish.size());
			}
			return retCode::EXPR_SYNTAX_ERROR;
		}
	}
	else {
		if (verbose > 0) {
			std::cout << "** error: could not tokenise expression \n";
			printTokens(tokens.data(), (int)tokens.size());
		}
	}
	return rv;
}

/* convert expression to ´tokens´. A token is basically either a number, operator,
function name or comma. Brackets are classed as operators.
Syntax is not checked properly at this stage but if there is anything that cannot
be tokenised EXPR_SYNTAX_ERROR is returned. If left & right brackets don't match up
EXPR_PAREN_MISMATCH is returned.
The normal return value is EXPR_OK.
*/
static retCode tokenise(const std::string expr, std::vector <token> &tokens) {
	int exprIndex = 0;
	int opIndex;
	token nxtToken;

	tokens.clear();

	while (exprIndex < expr.length()) {
		nxtToken.typecode = types::error;   /* should be replaced later by valid code */
		int charValue = toupper(expr[exprIndex]);
		/* first look for operators. */
		operSearch(expr.substr(exprIndex), opIndex);
		if (opIndex != -1) {
			/* found operator e.g. + - etc. */
			nxtToken.typecode = types::Operator;
			nxtToken.oper = operators[opIndex].operCode;
			nxtToken.value = 0;
			exprIndex += (int)strlen(operators[opIndex].oper);  // move index to next char after operator
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
			for (ptrdiff_t ix = 0; ix < (ptrdiff_t)functionList.size(); ix++) {
				if (_strnicmp(&expr[exprIndex],
					functionList[ix].fname, strlen(functionList[ix].fname)) == 0) {
					/* we have a match */
					nxtToken.typecode = types::func;
					nxtToken.function = ix;
					nxtToken.value = 0;                  /* not used */
					nxtToken.oper = opCode::power;  /* not used */
					exprIndex += (int)strlen(functionList[ix].fname); // move exprIndex past function name
					break;
				}
			}
		}

		if (nxtToken.typecode == types::error)
			return retCode::EXPR_SYNTAX_ERROR;  /* unable to tokenise expression*/
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
				return retCode::EXPR_PAREN_MISMATCH;
		}
	}
	if (brackets > 0)
		return retCode::EXPR_PAREN_MISMATCH;

	nxtToken.typecode = types::end;
	tokens.push_back(nxtToken);
	return retCode::EXPR_OK;
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
static void printTokens(const token expr[], const int exprLen) {
	if (exprLen == 0)
		return;   /* if no tokens to print do nothing but exit*/
	for (int ix = 0; ix < exprLen; ix++) {
		switch (expr[ix].typecode) {
		case types::number:
			std::cout << expr[ix].value << ' ';
			break;
		case types::func:
			std::cout << functionList[(int)expr[ix].function].fname << ' ';
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

		case types::Operator:
			if (expr[ix].oper == opCode::leftb)
				std::cout << '(';
			else if (expr[ix].oper == opCode::fact)
				std::cout << '!';
			else if (expr[ix].oper == opCode::dfact)
				std::cout << "!!";
			else if (expr[ix].oper == opCode::prim)
				std::cout << '#';
			else if (expr[ix].oper == opCode::rightb)
				std::cout << ')';
			else if (expr[ix].oper == opCode::unary_minus)
				std::cout << " Unary - ";
			else {
				for (int ixx = 0 ; ixx < sizeof(operators)/sizeof(operators[0]); ixx++)
					if (operators[ixx].operCode == expr[ix].oper) {
						std::cout << operators[ixx].oper << ' ';
						break;
					}
			}
			break;

		default:
			abort(); /* unrecognised token */
		}
	}
	std::cout << '\n';
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

		case types::number:	{
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
			int numparams = functionList[expr[exprIndex].function].NoOfParams;

			if (expr[exprIndex + 1].typecode != types::Operator ||
				expr[exprIndex + 1].oper != opCode::leftb) {
				return EXIT_FAILURE; /* function name not followed by (*/
			}
			int paramLen = 0;
			int ix3 = 0;
			for (int pcount = 1; pcount <= numparams; pcount++) {

				/* find delimiter marking end of parameter*/
				nextsep(expr + exprIndex + 2 + ix3, paramLen); // get next , or )
				if (expr[exprIndex + 2 + ix3 + paramLen].typecode == types::end)
					return EXIT_FAILURE; /* parameter sep. not found */
				int rv = reversePolish(expr + exprIndex + 2 + ix3, paramLen, rPolish);
				if (rv != EXIT_SUCCESS)
					return rv;  /* syntax error? */
				ix3 += (paramLen + 1); // move ix3 past , or )
			}
			if (expr[exprIndex + 1 + ix3].typecode != types::Operator ||
				expr[exprIndex + 1 + ix3].oper != opCode::rightb)
				return EXIT_FAILURE; /* no ) when required */

			rPolish.push_back(expr[exprIndex]); /* copy function token to output */
			exprIndex += (ix3 + 1); /* move past ) after function name */
			leftNumber = true;
			break;
		}

		case types::Operator: {
			if (expr[exprIndex].oper != opCode::leftb
				&& expr[exprIndex].oper != opCode::rightb) {

				if (expr[exprIndex].oper == opCode::minus && !leftNumber)
					expr[exprIndex].oper = opCode::unary_minus;  /* adjust op code */

				bool left = opr[(int)expr[exprIndex].oper].left; /* assocativity*/
				bool pre = opr[(int)expr[exprIndex].oper].pre; /* unary operator precedes expr*/
				bool unary = (opr[(int)expr[exprIndex].oper].numOps == 1);
				/* get priority of current operator */
				int expOpPri = opr[(int)expr[exprIndex].oper].pri;
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
					int stkOpPri = opr[(int)operStack.back().oper].pri;
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
(this would indicate a syntax error not detected earlier) */
static retCode evalExpr(const std::vector<token> &rPolish, Znum & result) {
	std::stack <Znum> nums;
	Znum val;
	Znum args[4];
	int index = 0;
	int rPlen = (int)rPolish.size();
	retCode retcode;

	for (index = 0; index < rPlen; index++) {
		/* push numbers onto stack */
		if (rPolish[index].typecode == types::number) {
			nums.push(rPolish[index].value);
		}

	 /* operators and functions are processed by taking the operand values
		from the stack, executing the operation or function and putting
		the returned value onto the stack. If there are insuffficient values
		on the stack or if the operator or function returns an error code
		exit immediately. */
		else if (rPolish[index].typecode == types::func) {
			fn_Code fnCode = functionList[rPolish[index].function].fCode;
			int NoOfArgs = functionList[rPolish[index].function].NoOfParams;

			if (NoOfArgs > nums.size())
				return retCode::EXPR_SYNTAX_ERROR;

			for (; NoOfArgs > 0; NoOfArgs--) {
				/* copy function parameters from number stack to args */
				args[NoOfArgs - 1] = nums.top(); /* use top value from stack*/
				nums.pop();  /* remove top value from stack */
			}
			retcode = ComputeFunc(fnCode, args[0], args[1], args[2], val);
			if (retcode != retCode::EXPR_OK)
				return retcode;
			nums.push(val);  /* put value returned by function onto stack */
		}

		else if (rPolish[index].typecode == types::Operator) {
			opCode oper = rPolish[index].oper;
			int NoOfArgs = opr[(int)oper].numOps;
			if (NoOfArgs > nums.size())
				return retCode::EXPR_SYNTAX_ERROR;  /* not enough operands on stack*/
			for (; NoOfArgs > 0; NoOfArgs--) {
				/* copy operand(s) from number stack to args */
				args[NoOfArgs - 1] = nums.top(); /* use top value from stack*/
				nums.pop();  /* remove top value from stack */
			}
			retCode rc = ComputeSubExpr(oper, args[0], args[1], val);
			if (rc != retCode::EXPR_OK)
				return rc;
			nums.push(val);  /* put generated value onto stack */
		}
		else
			abort(); /* unrecognised token */
	}


	if (nums.size() == 1) {
		result = nums.top();
		return retCode::EXPR_OK;
	}
	else
		return retCode::EXPR_SYNTAX_ERROR;  /* too many operands on stack*/
}

/* translate error code to text and output it*/
static void textError(retCode rc) {
	/*
	error codes currently used:
	EXPR_NUMBER_TOO_HIGH
	EXPR_INTERM_TOO_HIGH,		
	EXPR_DIVIDE_BY_ZERO,        
	EXPR_PAREN_MISMATCH,		
	EXPR_SYNTAX_ERROR,			
	EXPR_TOO_MANY_PAREN,		
	EXPR_INVALID_PARAM,	
	EXPR_ARGUMENTS_NOT_RELATIVELY_PRIME
	EXPR_EXPONENT_TOO_LARGE,
	EXPR_EXPONENT_NEGATIVE,				
	EXPR_OK = 0
	*/
	switch (rc)
	{
	//case EXPR_NUMBER_TOO_LOW:
	//	std::cout << (lang ? "Número muy pequeño\n" : "Number too low\n");
	//	break;
	case retCode::EXPR_NUMBER_TOO_HIGH:
		std::cout << (lang ? "Número muy grande \n" :
			"Number too high \n");
		break;
	case retCode::EXPR_INTERM_TOO_HIGH:
		std::cout << (lang ? "Número intermedio muy grande (más de 20000 dígitos\n" :
			"Intermediate number too high (more than 20000 digits)\n");
		break;
	case retCode::EXPR_DIVIDE_BY_ZERO:
		std::cout << (lang ? "División por cero\n" : "Division by zero\n");
		break;
	case retCode::EXPR_PAREN_MISMATCH:
		std::cout << (lang ? "Error de paréntesis\n" : "Parenthesis mismatch\n");
		break;
	case retCode::EXPR_SYNTAX_ERROR:
		if (lang) 	{
			std::cout << ( "Error de sintaxis\n");
		}
		else {
			std::cout << ("Syntax error\n");
		}
		break;
	case retCode::EXPR_TOO_MANY_PAREN:
		std::cout << (lang ? "Demasiados paréntesis\n" : "Too many parenthesis\n");
		break;
	case retCode::EXPR_INVALID_PARAM:
		std::cout << (lang ? "Parámetro inválido\n" : "Invalid parameter\n");
		break;
	case retCode::EXPR_ARGUMENTS_NOT_RELATIVELY_PRIME:
		std::cout << (lang ? "MCD de los argumentos no es 1\n" : "GCD of arguments is not 1\n");
		break;
	/*case EXPR_BREAK:
		std::cout << (lang ? "Detenido por el usuario\n" : "Stopped by use\nr");
		break;*/
	case retCode::EXPR_EXPONENT_NEGATIVE: {
		std::cout << "Exponent is negative\n";
		break;
	}
	case retCode::EXPR_EXPONENT_TOO_LARGE: {
		std::cout << "Exponent exceeds 2^31-1\n";
		break;
	}
	/*case EXPR_VAR_OR_COUNTER_REQUIRED:
		if (lang)
		{
			std::cout <<  "La expresión \n";
		}
		else
		{
			std::cout << "Expression #\n";
		}
		break;*/
	//case EXPR_BASE_MUST_BE_POSITIVE:
	//	std::cout << (lang ? "La base debe ser mayor que cero\n" :
	//		"Base must be greater than zero\n");
	//	break;
	//case EXPR_POWER_MUST_BE_POSITIVE:
	//	std::cout << (lang ? "La potencia debe ser mayor que cero\n" :
	//		"Power must be greater than zero\n");
	//	break;
	//case EXPR_MODULUS_MUST_BE_GREATER_THAN_ONE:
	//	std::cout << (lang ? "El módulo debe ser mayor que 1\n" : "Modulus must be greater than one\n");
	//	break;
	//case EXPR_MODULUS_MUST_BE_NONNEGATIVE:
	//	std::cout << (lang ? "El módulo no debe ser negativo\n" :
	//		"Modulus must not be negative\n");
	//	break;
	default:
		printf( "unknown error code: %d\n", (int)rc);
		break;
	}
}

/* convert s to UPPER CASE in d.  d could be the same string as s 
similar to _strupr in c */
static void strToUpper(const std::string &s, std::string &d) {
	if (&d != &s)
		d.resize(s.size());  // resize d unless it's the same string as s
	for (size_t ix = 0; ix < s.size(); ix++)
		d[ix] = toupper(s[ix]);
}

/* print elapsed time. If > 60 seconds print in hour min sec format */
static void PrintTimeUsed(double elapsed, const std::string &msg = "") {

	if (msg.size() > 1)
		std::cout << myTime() << ' ' << msg;
	auto elSec = elapsed / CLOCKS_PER_SEC; // convert ticks to seconds

	if (elSec > 10.0)
		PlaySound(TEXT("c:/Windows/Media/Alarm09.wav"), NULL,
			SND_FILENAME | SND_NODEFAULT | SND_ASYNC | SND_NOSTOP);

	if (elSec <= 60.0) {
		std::cout << elSec << " seconds\n";
	}
	else {
		/* round to nearest second */
		long long sec = (long long)std::floor(elSec); // convert to an integer
		long long min = sec / 60;  // min >= 1
		elSec -= min * 60;         // get seconds
		long long hour = min / 60; // hour may be zero
		min %= 60;                 // get minutes
		if (hour > 0)
			std::cout << ' ' << hour << "h " << min << " min " << elSec << "sec \n";
		else 
			std::cout << ' ' << min << " min " << elSec << " sec \n";
	}
}

/* remove spaces, tabs, etc  from msg */
static void removeBlanks(std::string &msg) {
	for (size_t ix = 0; ix < msg.size(); ix++) {
		if ((unsigned char)msg[ix] <= 0x7f && isspace(msg[ix])) {     // look for spaces, tabs, etc
			msg.erase(ix, 1);      // remove space character
			ix--;  // adjust index to take account of removed blank
		}
	}
}

struct summary {
	int numsize;		// number of decimal digits in number
	double time;		// time used in seconds
	int NumFacs;		// number of unique factors
	int totalFacs;		// total number of factors (=1 for a prime number)
	int testNum;		// test number (if applicable)
	int sndFac;		    // number of decimal digits in 2nd largest factor
};

std::vector <summary> results;

/* print summary - 1 line per test */
static void printSummary(void) {
	long long sec, min, hour;
	double elSec;
	printf_s("Test Num Size   time      Unique Factors Total Factors     2nd Fac\n");
	for (auto res : results) {
		/* round to nearest second */
		sec = (long long)std::floor(res.time); // convert to an integer
		min = sec / 60;  // min may be 0
		elSec = res.time - min * 60;         // get seconds including fractional part
		hour = min / 60; // hour may be zero
		min %= 60;       // get minutes
		printf_s("%4d %8d %2lld:%02lld:", res.testNum, res.numsize,
			hour, min);
		if (elSec < 10.0)
			printf_s("0%.2f ", elSec);  // insert leading zero
		else
			printf_s("%.2f ", elSec);

		printf_s("     %8d   %8d     %8d\n", res.NumFacs, res.totalFacs, res.sndFac);
	}
}

/* factorise Result, calculate number of divisors etc and print results */
static void doFactors(const Znum &Result, bool test) {
	fList factorlist;
	Znum Quad[4];
	clock_t start;
	summary sum;    // save summary of factorisation

	if (test)
		start = clock();	// used to measure execution time

	/* call DA´s magic function to factorise Result */
	bool rv = factorise(Result, factorlist, Quad);
	if (rv && factorlist.fsize() > 0) {
		factorlist.print(Result < 0);   /* print factors */

		if (abs(Result) != 1) {
			auto divisors = factorlist.NoOfDivs();
			std::cout << "\nNumber of Divisors = ";
			ShowLargeNumber(divisors, 6, false, false);
		}
		else 
			std::cout << "\nNumber of Divisors = 1";   // treat n=1 as special case

		auto divisors = factorlist.DivisorSum();
		std::cout << "\nSum of Divisors    = ";
		ShowLargeNumber(divisors, 6, false, false);
		divisors = factorlist.totient();
		std::cout << "\nTotient            = ";
		ShowLargeNumber(divisors, 6, false, false);
		if (Result > 0) {
			auto mob = factorlist.mob();  // mobius only defined for +ve integers
			std::cout << "\nMöbius             = " << mob;
		}

		/* show that the number is the sum of 4 or fewer squares. See
		https://www.alpertron.com.ar/4SQUARES.HTM */
		if (Result >= 0)
			std::cout << "\nn = ";
		else
			std::cout << "\n-n = ";
		char c = 'a';
		for (auto q : Quad) {  // print "n = a² + b² + c² + d²"
			if (q == 0)
				break;
			if (c > 'a') std::cout << "+ ";  // precede number with + unless its the 1st number
			std::cout << c << "² ";
			c++;    // change a to b, b to c, etc
		}
		c = 'a';
		for (auto q : Quad) {
			if (q == 0)
				break;
			std::cout << "\n" << c << "= " << q;
			c++;  // change a to b, b to c, etc
		}
		std::cout << "\n";
		factorlist.prCounts();  // print counts

		if (test) {
			/* recalculate result & get total number of factors */
			Znum result = factorlist.recheck(sum.totalFacs);
			if (result != Result) {
				std::cout << "Factors expected value " << Result << " actual value " << result << '\n';
				Beep(750, 1000);
			}
			result = Quad[0] * Quad[0] + Quad[1] * Quad[1] + Quad[2] * Quad[2] + Quad[3] * Quad[3];
			if (result != Result) {
				std::cout << "Quad expected value " << Result << " actual value " << result << '\n';
				Beep(750, 1000);
			}
			auto end = clock(); 
			double elapsed = (double)end - start;
			PrintTimeUsed(elapsed, "time used = ");

			/* store info for summary */
			sum.time = elapsed / CLOCKS_PER_SEC;
			sum.numsize = (int)ComputeNumDigits(result, 10);
			sum.NumFacs = (int)factorlist.fsize();
			/* get number of digits in 2nd largest factor */
			sum.sndFac = factorlist.sndFac();
			results.push_back(sum);
		}
	}
	else
		std::cout << " cannot be factorised\n";
}

/* perform some simple tests. Returns true if x3 is prime 
method = 0 for standard factorisation, != 0 to use only YAFU for factorisation */
static bool factortest(const Znum &x3, const int method=0) {
	fList factorlist;
	Znum Quad[4], result;
	long long totalFactors = 0;
	summary sum;    // save summary of test

	auto start = clock();	// used to measure execution time
	double end, elapsed;

	sum.numsize = (int)ComputeNumDigits(x3, 10);

	std::cout << "\nTest: factorise ";
	ShowLargeNumber(x3, 6, true, false);
	std::cout << '\n';
	if (method == 0) {
		factorise(x3, factorlist, Quad);  // get factors
	}
	if (method != 0) {
		factorlist.set(x3);
		callYafu(x3, factorlist);
	}

	/* recalculate result & get total number of factors */
	result = factorlist.recheck(sum.totalFacs);
	if (result != x3) {
		std::cout << "Factors expected value " << x3 << " actual value " << result << '\n';
		Beep(750, 1000);
	}
	
	/* get number of digits in 2nd largest factor */
	sum.sndFac = factorlist.sndFac();
	if (method == 0) {
		/* check that sum of squares is correct */
		result = Quad[0] * Quad[0] + Quad[1] * Quad[1] + Quad[2] * Quad[2] + Quad[3] * Quad[3];
		if (result != x3) {
			std::cout << "Quad expected value " << x3 << " actual value " << result << '\n';
			Beep(750, 1000);
		}
	}

	if (!factorlist.isPrime() ) {
		/* x3 is not prime */

		std::cout << "found " << factorlist.fsize() << " unique factors, total "
			<< sum.totalFacs << " factors\n";

		/* get number of digits in 2nd largest factor */
		sum.sndFac = factorlist.sndFac();
		if (method == 0)
			factorlist.prCounts();   // print counts

		end = clock();              // measure amount of time used
		elapsed = (double)end - start;
		PrintTimeUsed(elapsed, "time used = ");
		sum.time = elapsed / CLOCKS_PER_SEC;
		sum.NumFacs = (int)factorlist.fsize();
		results.push_back(sum);
		return false;   // not prime
	}
	else {
		std::cout << "is prime \n";
		sum.NumFacs = 1;
		sum.totalFacs = 1;
		sum.sndFac = 0;
		end = clock();              // measure amount of time used
		elapsed = (double)end - start;
		PrintTimeUsed(elapsed, "time used = ");
		sum.time = elapsed / CLOCKS_PER_SEC;
		results.push_back(sum);
		return true;    // is prime
	}
}

static void doTests(void) {
	Znum x3, x4, result, result1;
	/*std::vector <token> exprTokens;
	std::vector <token> rPolish;
	int rpLen;*/
	int i;
	int testcnt = 0;
	struct test {
		std::string text;        // text of expression to be evaluated
		long long expected_result;   // can only do tests that return a value <2^63
	};

	static test testvalues [] 
	{
		"2 - 3 + 4",                        3,
		"2 - (3+4)",                       -5,  // + and - have same priority, left-to-right evaluation
		"-5/2",                            -2,  // uses truncation division
		"-5%2",                            -1,
		"2-5/2",                            0,  // division has higher priority than -
		"(2-5)/2",                         -1,
		"7 * 11 / 2",                      38,  // * and / have same priority, left-to-right evaluation
		"7 * (11/2)",                      35,
		"-2^6",                            64,   // unary - has higher priority than ^
		"0-2^6",                          -64,   // ^ has higher priority than -
		"gcd (12,30)",                      6,
		"modinv (17,21)",                   5,
		"20 + 32^2/4*7",                 1812,    
		"20 + 32^2/(4*7)",                 56,
		"20 - 32^2/4*7",                -1772,
		"(20-32)^2 / 4 * 7",              252,   
		"(20-32)^2 / (4*7)",                5,
		"19!",             121645100408832000,   // factorial
		"19!!",  	                654729075,   // double factorial
		"37#",                  7420738134810,   // primorial
		"b(123456789)",	            123456761,   // prime before n
		"f(60)",	            1548008755920,   // fibonacci
		"l(60)",                3461452808002,   // lucas number
		"n(123456789)",             123456791,   // prime after n
		"p(150)",	              40853235313,   // number of partitions
		"modinv(7, 19)",                   11,   // modular inverse
		"7* 11 % 19",                       1,   // verfy modinv in previous test
		"modinv(11,19)",                    7,
		"modpow(8, 7, 6)",                  2,   // modular power
		"modpow(8, -2, 5)",                 4,   // negative exponent OK in some cases
		"totient(201)",		              132,
		"numdivs(7!)",                     60,
		"sumdivs(7!)",                  19344,
		"numdigits(123456789, 6)",         11,
		"sumdigits(123456789, 6)",         19,
		"revdigits(1234567890, 10)",  987654321, 
		"factconcat(2, 11!)", 22222222333355711,
		"17c7",                         19448,    // binomial coefficient
		"(17!) / ((17-7)!*7!)",         19448,
		"4 ^ 3 ^ 2",                   262144,    // NB expoentiation is right to left evaluation
		"4 ^ (3^2)",                   262144,
		"(4^3) ^ 2",                     4096,
		"4 ^ (3*2)",                     4096,
		"gcd(12, gcd(30, 40))",             2,   // nested function calls
		"7 * (11%19)",                     77,   // override normal left to right evaluation
		"(0x12345678) AND 0x018",        0x18,   // brackets round 1st number are needed so that A of AND is not considered part of number
		"0x12345678 OR 0x1",       0x12345679,
		"0x12345678 XOR 0x10",     0x12345668,
		"(NOT 0x0f23 4567 89ab cde1) * -1",  0x0f23456789abcde2,
		"5 < 6 == 7 < 8",                -1,   // returns true (== has lower priority)
		"5 < 6 != 7 < 8",                 0,   // returns false (!= has lower priority)
		"5 < (6 != 7) < 8",              -1,   // returns true; expr evaluated from left to right
		"R3(49)",                        54,
		"R2(585)",                       16,
		"SQRT(1234320)",               1110,
		"NROOT(2861381721051424,5)",   1234,
		"LLT(3217)",                      1,  // 2^3217-1 is prime
		"BPSW(2^99-1)",                   0,  //not a prime number
		"BPSW(2^127-1)",                  1,  //a prime number
		"-not1",                          2,  // operators are processed from right to left
		"not-1",                          0,  // operators are processed from right to left
		"not5#",                        -31,  // # operator evaluated before not
		"4^5#",         1152921504606846976,  // # operator evluated before exponent
		"5!!#",                       30030,  // !! operator evaluated before #
		"5#!!",           42849873690624000,  // # operator evaluated before !!
	};

	results.clear();

	auto start = clock();	// used to measure execution time
	for (i = 0; i < sizeof(testvalues) / sizeof(testvalues[0]); i++) {
		removeBlanks(testvalues[i].text);  // it is necessary to remove spaces
		/* but it is not necessary to convert to upper case */
		//stackIndex = 0;         // ensure stack is empty 

		auto  rv =ComputeExpr(testvalues[i].text, result);
		if (rv != retCode::EXPR_OK || result != testvalues[i].expected_result) {
			std::cout << "test " << i + 1 << " failed\n" <<
				"expected " << testvalues[i].text << " = " << testvalues[i].expected_result << '\n';
			std::cout << "got " << result << '\n';
			Beep(750, 1000);
		}
	}
	std::cout << i << " tests completed\n";

	for (Znum i = 1000; i <= 100000000000000000; ) {
		Znum x1, x2;
		mpz_nextprime(ZT(x1), ZT(i));  // get next prime
		i *= 10;
		mpz_nextprime(ZT(x2), ZT(i));  // get next prime
		x3 = x1*x2;
		factortest(x3);
		results.back().testNum = ++testcnt;
		x3++;
		factortest(x3);
		results.back().testNum = ++testcnt;

	}

	/* tests below have shown a problem with pollard-rho for certain numbers */
	factortest(99999999973789); // = 6478429 * 15435841
	results.back().testNum = ++testcnt;
	factortest(183038861417);   // =  408229 *   448373
	results.back().testNum = ++testcnt;
	factortest(183475587821);   // =  409477 *   448073
	results.back().testNum = ++testcnt;
	factortest(181919916457);   // =  400307 *   454451
	results.back().testNum = ++testcnt;
	factortest(199996999457);   // =  441361 *   453137
	results.back().testNum = ++testcnt;
	factortest(204493418837);   // =  401209 *   509693
	results.back().testNum = ++testcnt;

	/* exercise code specifically for power +/-1 */
	mpz_ui_pow_ui(ZT(x3), 10, 20);  // x3 = 10^20
	x3 -= 1;                        // x3 = 10^20-1
	factortest(x3);
	results.back().testNum = ++testcnt;
	x3 += 2;
	factortest(x3);
	results.back().testNum = ++testcnt;
	ComputeExpr("n(10^10)^2-1", x3);  // test power of large number-1
	factortest(x3);
	results.back().testNum = ++testcnt;

	ComputeExpr("120#-1", x3);
	factortest(x3);
	results.back().testNum = ++testcnt;
	ComputeExpr("n(10^15)^2", x3);  // test power of large number
	factortest(x3);
	results.back().testNum = ++testcnt;
	ComputeExpr("n(10^6+20)^1667", x3);  // test power of large prime number
	factortest(x3);
	results.back().testNum = ++testcnt;
	ComputeExpr("n(10^7)^3*n(10^8)^2", x3);  // test powers of two large prime number
	factortest(x3);
	results.back().testNum = ++testcnt;

	ComputeExpr("n(3*10^9+50)*n(3*10^10+500)", x3);  // test Lehman factorisation
	factortest(x3);
	results.back().testNum = ++testcnt;
	x3 = 0;
	ComputeExpr("n(10^15)^3*n(10^14)", x3);  // test Lehman factorisation
	factortest(x3);
	results.back().testNum = ++testcnt;

	/* test using carmichael numbers. note that 1st example has no small factors  */
	long long int carmichael[] = { 90256390764228001, 7156857700403137441,  1436697831295441,
		60977817398996785 };
	for (int i = 0; i < sizeof(carmichael) / sizeof(carmichael[0]); i++) {
		factortest(carmichael[i]);
		results.back().testNum = ++testcnt;
	}

	/* set x3 to large prime. see https://en.wikipedia.org/wiki/Carmichael_number */
	ComputeExpr("2967449566868551055015417464290533273077199179985304335099507"
		"5531276838753171770199594238596428121188033664754218345562493168782883", x3);
	x4 = x3 * (313 * (x3 - 1) + 1) * (353 * (x3 - 1) + 1);
	/* in general numbers > about 110 digits cannot be factorised in a reasonable time 
	but this one can, because a special algorithm just for Carmichael numbers is used. */
	factortest(x4);   // 397-digit Carmichael number
	results.back().testNum = ++testcnt;
	std::cout << "factorised 397-digit Carmichael number \n";

	ComputeExpr("n(10^24)*n(10^25)*n(10^26)*n(10^27)", x3);  
	factortest(x3);
	results.back().testNum = ++testcnt;

	ComputeExpr("n(2^97)*n(2^105)", x3);
	factortest(x3);
	results.back().testNum = ++testcnt;

	auto end = clock();   // measure amount of time used
	double elapsed = (double)end - start;
	PrintTimeUsed(elapsed, "Factorisation tests completed. Time used= ");
	printSummary();
}

/* generate large random number, up to 128 bits */
//static void largeRand(Znum &a) {
//	a = ((long long)rand() << 32) + rand();
//	a <<= 64;
//	a += ((long long)rand() << 32) + rand();
//}

//find a random 'strong' prime of size 'bits'
//follows the Handbook of applied cryptography
static void gordon(Znum &p, gmp_randstate_t &state, const long long bits) {
	/*
	SUMMARY: a strong prime p is generated.
	1. Generate two large random primes s and t of roughly equal bitlength (see Note 4.54).
	2. Select an integer i0. Find the first prime in the sequence 2it + 1,
	for i = i0; i0 + 1; i0 + 2; : : : (see Note 4.54). Denote this prime by r = 2it+ 1.
	3. Compute p0 = 2(sr-2 mod r)s - 1.
	4. Select an integer j0. Find the first prime in the sequence p0 +2jrs,
	for j = j0; j0 + 1; j0 + 2; : : : (see Note 4.54). Denote this prime by
	p = p0 + 2jrs.
	5. Return(p).

  4.54 Note (implementing Gordon’s algorithm)
	(i) The primes s and t required in step 1 can be probable primes generated by
	Algorithm 4.44. The Miller-Rabin test (Algorithm 4.24) can be used to test
	each candidate for primality in steps 2 and 4, after ruling out candidates
	that are divisible by a small prime less than some bound B. See Note 4.45 for
	guidance on selecting B. Since the Miller-Rabin test is a probabilistic
	primality test, the output of this implementation of Gordon’s algorithm is
	a probable prime.
	(ii) By carefully choosing the sizes of primes s, t and parameters i0, j0,
	one can control the exact bitlength of the resulting prime p. Note that the
	bitlengths of r and s will be about half that of p, while the bitlength of
	t will be slightly less than that of r.
	*/
	Znum r, s, t, t2, p0;

	//1. s and t should be about half the bitlength of p
	mpz_urandomb(ZT(s), state, bits/2-2);
	mpz_nextprime(ZT(s), ZT(s));        // make s prime

	mpz_urandomb(ZT(t), state, bits/2 -2);
	mpz_nextprime(ZT(t), ZT(t));        // make t prime


	// 2 Find the first prime r in the sequence 2t + 1, 4t+1 6t+1 ...
	t2 = t * 2;
	r = t2 + 1;
	while (!mpz_likely_prime_p(ZT(r), state, 0))
		r += t2;


	// 3. Compute p0 = 2(sr-2 mod r)s - 1.
	p0 = ((s*r - 2) % r)*s * 2 - 1;

	// 4. Find the first prime p in the sequence p0, p0 +2rs p0+4rs ....,
	p = p0;
	while (!mpz_likely_prime_p(ZT(p), state, 0))
		p += 2 * r*s;
	return;
}

/* generate RSA-style difficult to factor number, size = bits +/- 1 */
static void get_RSA(Znum &x, gmp_randstate_t &state, const long long bits) {
	Znum p, q;
	x = 0;

	/* keep generating random strong primes till we get a pair that produces
	a product of about the right size. */
	while (abs(NoOfBits(x)-bits) >1) {
		gordon(p, state, bits / 2);
		gordon(q, state, bits / 2);
		x = p * q;
	}
	return;
}

/* test factorisation using pseudo-random numbers of a specified size.
Command format is TEST2 [num1[,num[,num3]]] where 
num1 is the number of tests, 
num2 is the size of the numbers to be factored in bits,
if num3  NE 0 the number to be factored consists of 2 approximately same-sized prime
factors, otherwise it is a random number that can contain any number of factors. */
static void doTests2(const std::string &params) {
	long long p1;  // number of tests; must be greater than 0
	long long p2;  // size of numbers to be factored in bits (>= 48)
	long long p3=0;  // optional, if non-zero generate RSA-style difficult to factor number
	gmp_randstate_t state;
	Znum x;

	auto numParams = sscanf_s(params.data(), "%lld,%lld,%lld", &p1, &p2, &p3);
	if (p1 <= 0 || numParams < 1) {
		std::cout << "Use default 2 for number of tests \n";
		p1 = 2;
	}
	if (p2 < 48 || numParams < 2) {
		std::cout << "Use default 48 for number size in bits \n";
		p2 = 48;
	}
	if (p2 >= 500)
		std::cout << "**warning: factoring such large numbers could take weeks! \n";
	gmp_randinit_mt(state);  // use Mersenne Twister to generate pseudo-random numbers
	gmp_randseed_ui(state, 756128234);
	/* fixed seed means that the exact same tests can be repeated provided
	   that the same size of number is used each time */

	auto start = clock();	// used to measure execution time

	results.clear();

	for (int i = 1; i <= p1; i++) {
		std::cout << "\nTest # " << i << " of " << p1 << '\n';
		if (p3 == 0)
			mpz_urandomb(ZT(x), state, p2);  // get random number, size=p2 bits
		else
			get_RSA(x, state, p2);  // get RSA type number size p2 bits
		ShowLargeNumber(x, 6, true, false);
		std::cout << '\n';
		doFactors(x, true); /* factorise x, calculate number of divisors etc */
		results.back().testNum = i;
	}

	auto end = clock();   // measure amount of time used
	double elapsed = (double)end - start;
	PrintTimeUsed(elapsed, "All tests completed. Time used = ");
	printSummary();
	gmp_randclear(state);  // clear state - avoid memory leakage
}

#ifdef BIGNBR
/*  1. check basic arithmetic operators for BigIntegers
	2. test BigInteger multiplication with larger numbers
	3. BigInteger division with larger numbers
	4. Modular Multiplication using Mongomery Encoding (REDC)
*/
static void doTests3(void) {
	Znum a, a1, am, b, b1, bm, mod, p, p2, pm;
	limb aL[MAX_LEN], modL[MAX_LEN], alM[MAX_LEN], al2[MAX_LEN];
	limb bL[MAX_LEN], blM[MAX_LEN], bl2[MAX_LEN], pl[MAX_LEN], plm[MAX_LEN];
	limb one[MAX_LEN];
	int numLen = MAX_LEN-2, l2;

	BigInteger aBI, a1BI, bBI, amBI, pBI;

	auto start = clock();	// used to measure execution time

	memset(one, 0, MAX_LEN * sizeof(limb));
	one[0].x = 1;                   /* set value of one to 1 */

	srand(421040034);               // seed random number generator 
	
	/* check basic arithmetic operators for BigIntegers */
	a = 1291;        // set starting values (increased each time round loop)
	b = 131;
	for (int ctr = 1; ctr <= 1700; ctr++) {
		/* with the values currently used for a and b, ctr cannot go much beyond 1700
		or overflow would occur. */
		aBI = a;                      // copy value of a to aBI (BigInteger)
		bBI = b;
		BigtoZ(a1, aBI);              // copy value to a1 (Znum)
		assert(a1 == a);              // verify conversion to & from Biginteger

		pBI = aBI + bBI;
		BigtoZ(p, pBI);
		assert(p == (a + b));         // verify addition        

		pBI = aBI;
		pBI += bBI;
		assert(p == (a + b));         // verify addition

		pBI = aBI - bBI;
		BigtoZ(p, pBI);
		assert(p == (a - b));         // verify subtraction

		pBI = aBI;
		pBI -= bBI;
		assert(p == (a - b));         // verify subtraction

		pBI = aBI * bBI;
		BigtoZ(p, pBI);
		Znum prod;
		Znum diff;
		prod = a * b;
		diff = p ^ prod;   // p XOR prod (0 if they are equal)
		if (diff != 0) {              // verify multiplication
			std::cout << "a = " << a << " b = " << b << '\n';
			gmp_printf("p     = %Zx \na * b = %Zx \n", p, prod);
			gmp_printf("diffference = %Zx \n", diff);
		}

		pBI = aBI;
		pBI *= bBI;
		BigtoZ(p, pBI);               // p = pBI = aBI * bBI
		assert(p == (a * b));         // verify multiplication

		pBI = aBI;
		pBI *= INT_MAX;
		BigtoZ(p, pBI);               // p = pBI = aBI * INT_MAX
		assert(p == a * INT_MAX);	  // verify multiplication by int

		pBI = aBI / bBI;
		BigtoZ(p, pBI);                // p = pBI = aBI / bBI
		if (p != (a / b)) {              // verify division
			std::cout << "a = " << a << "\nb = " << b << '\n';
			std::cout << "p = " << p << "\na/b = " << a / b << '\n';
		}

		pBI = aBI;
		pBI /= bBI;            
		BigtoZ(p, pBI);                  // p = pBI = aBI / bBI
		if (p != (a / b)) {              // verify division
			std::cout << "a = " << a << " b = " << b << '\n';
			std::cout << "p = " << p << "\na/b = " << a / b << '\n';
		}

		pBI = aBI % bBI;
		BigtoZ(p, pBI);               // p = pBI = aBI % bBI
		if (p != (a%b))               // verify modulus
			std::cout << "p = " << p << "\na%b = " << a % b << '\n';

		pBI = aBI;
		pBI %= bBI;
		BigtoZ(p, pBI);                 // p = pBI = aBI % bBI
		if (p != (a%b))               // verify modulus
			std::cout << "p = " << p << "\na%b = " << a % b << '\n';

		/* check BigtoLL function */
		int exp;
		a1 = BigToLL(aBI, exp);  // a1 * 2^exp = aBI
		mpz_div_2exp(ZT(b1), ZT(a), exp);   // b1 = a/(2^exp)
		if (b1 != a1)
			gmp_printf("a = %Zx a1 = %Zx, exp = %d \n", a, a1, exp);
		assert(b1 == a1);

		pBI = aBI << 4;      // check left shift
		BigtoZ(p, pBI);      // p = pBI = aBI << 4
		p2 = a << 4;
		if (p != p2) {
			gmp_printf("a = %#Zx p = %#Zx expected %#Zx \n", a, p, p2);
		}


		pBI = aBI << 56;    // check left shift
		BigtoZ(p, pBI);     // p = pBI = aBI << 56
		p2 = a << 56;
		if (p != p2) {
			gmp_printf("a = %#Zx p = %#Zx expected %#Zx \n", a, p, p2);
		}

		pBI = aBI << 567;      // check left shift
		BigtoZ(p, pBI);        // p = pBI = aBI << 567
		p2 = a << 567;
		if (p != p2) {
			gmp_printf("a = %#Zx p = %#Zx expected %#Zx \n", a, p, p2);
		}

		pBI = aBI + INT_MAX;
		BigtoZ(p, pBI);              // p = pBI = aBI + max;
		assert(p == a + INT_MAX);    // veryify addition of int

		pBI = aBI - INT_MAX;
		BigtoZ(p, pBI);               // p = pBI = aBI - max;
		assert(p == a - INT_MAX);;    // veryify subtraction of int

		a *= 1237953;                 // increase a & b, then repeat
		b *= 129218;
	}

	auto end = clock();              // measure amount of time used
	double elapsed = (double)end - start;
	std::cout << "test stage 1 completed  time used= " << elapsed / CLOCKS_PER_SEC << " seconds\n";
	
	/* test multiplication with larger numbers */
	aBI = 0x7fffffffffffffff;
	a1BI = aBI;
	a = 0x7fffffffffffffff;
	for (bBI = 1; bBI <= 256; bBI++) {
		aBI = aBI * a1BI;
		a = a * 0x7fffffffffffffffLL;
		BigtoZ(a1, aBI);

		if (a1 != a) {
			gmp_printf("aBI = %Zx \n", ZT(a1));
			gmp_printf("  a = %Zx \n", ZT(a));
		}
	}

	/* check multiplication of large numbers */
	aBI = 17341;
	a = 17341;
	for (int i = 1; i < 13; i++) {
		if (aBI.bitLength() > MAX_LEN*BITS_PER_GROUP / 2)
			break;  // exit loop if a is too large to square it
		aBI = aBI * aBI;
		a = a * a;
		BigtoZ(a1, aBI);

		/*gmp_printf("aBI = %Zx \n", ZT(a1));
		gmp_printf("  a = %Zx \n", ZT(a));*/
		assert(a1 == a);
	}

	end = clock();              // measure amount of time used
	elapsed = (double)end - start;
	std::cout << "test stage 2 completed  time used= " << elapsed / CLOCKS_PER_SEC << " seconds\n";

	
	/* check conversion to & from floating point */
	/*for (int i = 1; i <= 30; i++) {
		expBigInt(amBI, i);
		std::cout << "e^=" << amBI << " expected " << std::exp(i) << '\n';
	}*/

	p = 10;

	Znum error, relerror;
	for (int i = 1; i < 99;  i++, p *= 1000000 ) {
		if (!ZtoBig(pBI, p))
			break;            // p too large to convert so stop
		double pdb = pBI.log();  // get natural log of p
		/*if (pdb > 708)
			break;*/
		expBigInt(amBI, pdb);   // convert back from log
		BigtoZ(am, amBI);       // convert back tp Znum
		error = p - am;         // get error
		if (error != 0) {
			double e1, e2, relErrf;
			long e1l, e2l;
			e1 = mpz_get_d_2exp(&e1l, ZT(p));
			e2 = mpz_get_d_2exp(&e2l, ZT(am));
			assert(e1l == e2l);   // check power of 2 exponent is correct
			relErrf = abs(e1 - e2) / e1;   /* should be less than 10E-14 */
			relerror = (10000000000000000LL * error) / p;
			/* print log base 10 of p, then error ratio */
			std::cout << " pdb = " << pdb / std::log(10)
				<< " error = " << relErrf <<  " relerror = " << relerror << '\n';
		}
	}

	/* check division with large numbers */
	p = 12345678901;
	p2 = 11;
	for (l2 = 1; l2 < 4560; l2++) {
		if (!ZtoBig(pBI, p))
			break;                 // exit loop if p is too big to fit BigInteger (approximately 20000 digits)
		if (!ZtoBig(amBI, p2))
			break;
		bBI = pBI / amBI;         // calculçate p/p2 using Bigintegers
		BigtoZ(b1, bBI);          // convert quotient back to Znum in b1
		if (b1 != p / p2) {
			std::cout << "p      = " << p
				<< "\np2     = " << p2
				<< "\nexpected " << p / p2
				<< "\ngot      " << b1;
			break;
		}
		bBI = pBI % amBI;
		BigtoZ(b1, bBI);
		if (b1 != p % p2) {
			std::cout << "p      = " << p
				<< "\np2     = " << p2
				<< "\nexpected " << p % p2
				<< "\ngot      " << b1;
			break;
		}
		p *= (long long)rand() * 2 + 3;
		p2 *= rand();
	}
	std::cout << "division test " << l2 << " cycles\n";
	end = clock();              // measure amount of time used
	elapsed = (double)end - start;
	std::cout << "test stage 3 completed  time used= " << elapsed / CLOCKS_PER_SEC << " seconds\n";

	/* test to compare modular multiplication using DA's code against GMP/MPIR
	BigIntegers. Both use Mongomery notation for the integers to avoid slow
	division operations. The conclusion is that GMP takes about twice as long
	as DA's code. */

	/* set up modulus and Mongomery parameters */
	largeRand(mod);					// get large random number
	mod |= 1;                       // set lowest bit (make sure mod is odd)
	GetMontgomeryParms(mod);
	ZtoLimbs(modL, mod, MAX_LEN);    // copy value of mod to modL
	while (modL[numLen - 1].x == 0)
		numLen--;                    // adjust length i.e. remove leading zeros
	memcpy(TestNbr, modL, numLen * sizeof(limb));  // set up for GetMontgomeryParms
	NumberLength = numLen;
	GetMontgomeryParms(numLen);
	modmultCallback = nullptr;      // turn off status messages from modmult

	largeRand(a);				     // get large random number a
	a %= mod;						 // ensure a < mod
	modmult(a, zR2, am);             // convert a to Montgomery (Znum) in am
	ZtoLimbs(aL, a,numLen);		     // copy value of a to aL (limbs)
	modmult(aL, MontgomeryMultR2, alM);  // convert a to Mongomery (limbs)
	modmult(alM, one, al2);          // convert a from Mongomery (limbs) 
	LimbstoZ(al2, a1, numLen);       // copy value to a1 (Znum)
	assert(a == a1);                 // check that all thes conversions work properly

	largeRand(b);				     // get large random number  b
	modmult(b%mod, zR2, bm);         // convert b to Montgomery (Znum)
	ZtoLimbs(bL, b, numLen);		 // copy value of b to bL
	modmult(bL, MontgomeryMultR2, blM);  // convert b to Mongomery (limbs)

	for (int i = 1; i < 200000000; i++) {
		modmult(alM, blM, plm);          // p = a*b mod m (limbs)
		memcpy(alM, plm, numLen * sizeof(limb));  // a = p (limbs)
		modmult(am, bm, pm);             // p = a*b mod m (Znum)
		am = pm;						 //a = p (Znum)
	}

	REDC(p, pm);					 // convert p from montgomery (Znum)
	modmult(plm, one, pl);           // convert p from Mongomery (limbs)
	LimbstoZ(pl, p2, numLen);        // copy value to p2 (Znum)
	assert(p2 == p);

	end = clock();              // measure amount of time used
	elapsed = (double)end - start;
	std::cout << "tests completed  time used= " << elapsed / CLOCKS_PER_SEC << " seconds\n";
}
#endif

/* tests for r3 function */
/* see http://oeis.org/A002102 */

extern unsigned __int64 R2(const unsigned __int64 n);
//static void doTests4(void) {
//	generatePrimes(2000);
//	for (int i = 0; i <= 9992; i++) {
//		auto rv2 = R2(i);
//		auto rv = R3(i);
//		printf_s("%4d %4lld %4lld ", i, rv2, rv);
//		std::cout << '\n';
//	}
//}

/* test factorisation of mersenne numbers see https://en.wikipedia.org/wiki/Mersenne_prime
this test will take about 2.5 hours 
*/
static void doTests4(void) {
	int px = 0;
	const int pmax = 75;
	Znum m;

	results.clear();
	auto start = clock();	// used to measure execution time
	if (primeListMax < 1000)
		generatePrimes(393203);

	for (px = 0; px <= pmax; px++) {
		if (px >= 70 && !yafu) {
			std::cout << "remaining tests skipped \n";
			break;
		}
		std::cout << "\ntest " << px + 1 << " of " << pmax+1 ;
		mpz_ui_pow_ui(ZT(m), 2, primeList[px]);  // get  m= 2^p
		m--;                // get 2^p -1
		if (factortest(m)) /* factorise m, calculate number of divisors etc */
			std::cout << "2^" << primeList[px] << " - 1 is prime \n";
		results.back().testNum = px + 1;
	}
	auto end = clock();              // measure amount of time used
	auto elapsed = (double)end - start;
	PrintTimeUsed(elapsed, "tests completed time used = ");
	printSummary();           // print summary 1 line per test
}

/* tests using only YAFU for factorisation */
static void doTests5(void) {
	results.clear();
	int testcnt = 0;

	Znum num = 49728103; // 7001 * 7103
	factortest(num, 1);
	results.back().testNum = ++testcnt;
	std::cout << "test " << testcnt << " completed \n";  // test 1

	// factorise 127 digit number
	/*  =                               280673 
	*                              2756 163353 
	*                            598990 818061 
	*                          4 527716 228491 
	*                     248158 049830 971629 
	* 33637 310674 071348 724927 955857 253537  
	*117445 937227 520353 139789 517076 610399  
	(7 factors)   */
	ComputeExpr("2056802480868100646375721251575555494408897387375737955882170045672576386016591560879707933101909539325829251496440620798637813", num);
	factortest(num, 1);
	results.back().testNum = ++testcnt;
	std::cout << "test " << testcnt << " completed \n";  // test 2

	//factorise 57 digit number
	/* P6 = 280673
	  P12 = 598990818061
	  P10 = 2756163353
	  P13 = 4527716228491
	  P18 = 248158049830971629
	*/
	ComputeExpr("520634955263678254286743265052100815100679789130966891851", num);
	factortest(num, 1);
	results.back().testNum = ++testcnt;
	std::cout << "test " << testcnt << " completed \n";  // test 3

	//factorise 80 digit number (about 3 minutes)
	/* P49 = 2412329883909990626764837681505523221799246895133
       P32 = 18138544144503826077310252140817
    */
	ComputeExpr("43756152090407155008788902702412144383525640641502974083054213255054353547943661", num);
	factortest(num, 1);
	results.back().testNum = ++testcnt;
	std::cout << "test " << testcnt << " completed \n";   // test 4

	//factorise 85 digit number (about 7 mins)
	/* factors are 1485325304578290487744454354798448608807999 and 
                   1263789702211268559063981919736415575710439 */
	ComputeExpr("1877138824359859508015524119652506869600959721781289179190693027302028679377371001561", num);
	factortest(num, 1);
	results.back().testNum = ++testcnt;
	std::cout << "test " << testcnt << " completed \n";  // test 5

	// factorise 94 digit number (about 60 mins)
	/* factors are 10910042366770069935194172108679294830357131379375349 and 
                   859735020008609871428759089831769060699941 */
	ComputeExpr("9379745492489847473195765085744210645855556204246905462578925932774371960871599319713301154409", num);
	factortest(num, 1);
	results.back().testNum = ++testcnt;
	std::cout << "test " << testcnt << " completed \n";  // test 6

	//factorise 100 digit number - takes many hours
	/* factor are 618162834186865969389336374155487198277265679 and
	              4660648728983566373964395375209529291596595400646068307 */
	ComputeExpr("2881039827457895971881627053137530734638790825166127496066674320241571446494762386620442953820735453", num);
	factortest(num, 1);
	results.back().testNum = ++testcnt;
	std::cout << "test " << testcnt << " completed \n";  // test 7

	printSummary();    // print summary - 1 line per test

	return;
}

/* tests using only Msieve for factorisation. Factorise selected Mersenne numbers */
static void doTests6(void) {
	bool yafusave = yafu;
	bool msievesave = msieve;
	Znum m;
	msieve = true;
	yafu = false;
	int testcnt = 0;

	results.clear();

	mpz_ui_pow_ui(ZT(m), 2, 277);  // get  m= 2^p
	m--;                // get 2^p -1
	factortest(m);      // 84 digits
	results.back().testNum = ++testcnt;
	/* p7  factor: 1121297
	   p38 factor: 31133636305610209482201109050392404721
	   p40 factor: 6955979459776540052280934851589652278783
	*/

	mpz_ui_pow_ui(ZT(m), 2, 293);  // get  m= 2^p
	m--;                // get 2^p -1
	factortest(m);      // 89 digits
	results.back().testNum = ++testcnt;
	/* p26 factor: 40122362455616221971122353
       p63 factor: 396645227028138890415611220710757921643910743103031701971222447 */

	mpz_ui_pow_ui(ZT(m), 2, 311);  // get  m= 2^p
	m--;                // get 2^p -1
	factortest(m);      // 94 digits
	results.back().testNum = ++testcnt;
	/* p7  factor: 5344847
	   p31 factor: 2647649373910205158468946067671
	   p57 factor: 294803681348959296477194164064643062187559537539328375831
	*/

	mpz_ui_pow_ui(ZT(m), 2, 313);  // get  m= 2^p
	m--;                // get 2^p -1
	factortest(m);      // 95 digits
	results.back().testNum = ++testcnt;
	/* p8  factor: 10960009
	   p17 factor: 14787970697180273
	   p25 factor: 3857194764289141165278097
	   p47 factor: 26693012026551688286164949958620483258358551879
	*/

	mpz_ui_pow_ui(ZT(m), 2, 353);  // get  m= 2^p
	m--;                // get 2^p -1
	factortest(m);      // 107 digits
	results.back().testNum = ++testcnt;
	/* p34 factor: 2927455476800301964116805545194017
	   p67 factor: 6725414756111955781503880188940925566051960039574573675843402666863
    */
	
	yafu = yafusave;
	msieve = msievesave;
	printSummary();           // print 1 line per test
}

/* Lucas-Lehmer test*/
static void doTests7(const std::string &params) {
	std::vector <long long> mPrimes;
	int p1, i=0;
	lltPriCnt = lltCmpCnt = lltTdivCnt = 0;  /* reset counters */
#ifdef _DEBUG
	int limit = 2000;  /* find 1st 15 Mersenne  primes, test 303 primes */
#else
	int limit = 12000; /* find 1st 23 Mersenne primes, test 1438 primes */
#endif
	auto start = clock();	// used to measure execution time

	auto numParams = sscanf_s(params.data(), "%d", &p1);
	if (numParams >= 1)
		limit = p1;    /* replace default value with user-specified value */
	generatePrimes(limit);    /* make prime list if not already done */

	for (i = 0; primeList[i] < limit; i++) {
		Znum p = primeList[i];
		Znum rv = llt(p); /* Return 0 if 2^p-1 is composite, 1 if prime  */
		if (rv == 1) {
			std::cout << myTime() << " 2^" << primeList[i] << " -1 is prime *** \n";
			mPrimes.push_back(primeList[i]);
		}
		else if (verbose > 0 || (i & 0x3f) == 0)
			/* \r instead of usual \n means that each messsage overwrites the 
			previous one */
			std::cout << myTime() << " 2^" << primeList[i] << " -1 is NOT PRIME \r";
	}

	/* print the results */
	/* see https://oeis.org/A000043 */
	std::cout << "Found " << mPrimes.size() << " Mersenne primes  out of " << i 
		<< " numbers tested \n" ;
	for (auto p : mPrimes) {
		std::cout << p << ", ";
	}
	putchar('\n');
	long long other = (long long)i - (lltPriCnt+ lltTdivCnt+ lltCmpCnt);
	printf_s("%5.2f%% primes found by llt \n", 100.0 *lltPriCnt / i);
	if (lltTdivCnt != 0)
		printf_s("%5.2f%% composites found by trial division \n", 100.*lltTdivCnt / i);
	if (lltCmpCnt != 0)
		printf_s("%5.2f%% composites found by llt \n", 100.0*lltCmpCnt / i);
	if (other != 0)
		printf_s("%5.2f%% other \n", 100.0*other / i);

	auto end = clock();              // measure amount of time used
	auto elapsed = (double)end - start;
	PrintTimeUsed(elapsed, "test 7 completed time used = ");
}


std::vector<std::string> inifile;  // copy contents of .ini here
std::string iniPath;          // path to .ini file
std::string helpFilePath = "docfile.txt";  /* can be overwritten from the .ini file */

/* (re)write the BigIntCalculator.ini file
initially a .new file is created, then any .old file is deleted
then the current .ini if any is renamed as .old,
then the .new file is renamed as .ini */
void writeIni(void) {
	std::string newFname = iniPath + "BigIntCalculator.new";
	std::string iniFname = iniPath + "BigIntCalculator.ini";
	std::string oldFname = iniPath + "BigIntCalculator.old";

	std::ofstream newStr(newFname, std::ios::out);  // open .new file for output
	if (!newStr.is_open()) {
		std::cerr << "cannot open BigIntCalculator.new \n";
		return;
	}
	/* copy from current ini file to new one */
	for (auto text : inifile) {
		newStr << text << '\n';
	}
	newStr << "yafu-path=" << YafuPath << '\n';
	newStr << "yafu-prog=" << yafuprog << '\n';
	newStr << "msieve-path=" << MsievePath << '\n';
	newStr << "msieve-prog=" << MsieveProg << '\n';
	newStr << "helpfile=" << helpFilePath << '\n';
	newStr.close();

	  // delete any previous .old
	int rc = remove(oldFname.c_str());
	if (rc != 0 && errno != ENOENT) {
		perror("could not remove BigIntCalculator.old file ");
	}
    // rename .ini as .old
	int rv = rename(iniFname.c_str(), oldFname.c_str());   
	if (rv == 0 || errno == ENOENT)
		rv = rename(newFname.c_str(), iniFname.c_str());   // .new -> .ini
	else
		perror("unable to rename BigIntCalculator.ini as BigIntCalculator.old");
}

/* read the .ini file and update paths. 
path definitions begin with yafu-path=, yafu-prog=, msieve-path= or msieve-prog= 
Paths are not case-sensitive. 
Anything else is saved and is copied if the .ini file is updated 
arg is arg[0] of the call to main, which conveniently contains the full path
for the .exe file. We use the same path for the .ini file. 
If the .ini file can't be opened a new one is created */
static void processIni(const char * arg) {
	std::string iniFname;
	std::string buffer;

	iniPath = arg;  // save path in a global variable
	auto b = iniPath.find_last_of("\\/"); // find last / character
	if (b != std::string::npos) {
		iniPath.erase(b + 1);  // remove everything after /
		iniFname = iniPath + "BigIntCalculator.ini";
	}
	else 
		iniFname = "BigIntCalculator.ini";

	std::ifstream iniStr(iniFname, std::ios::in);  // open .ini file
	if (!iniStr.is_open()) {
		/* can't open .ini file - make one from scratch */
		const time_t currtime = time(NULL);  // time as seconds elapsed since midnight, January 1, 1970
		struct tm mytime;
		char timestamp[23];   // date & time in format "dd/mm/yyyy at hh:mm:ss"

		localtime_s(&mytime, &currtime);  // convert to tm format
		/* convert to dd/mm/yyyy hh:mm:ss */
		strftime(timestamp, sizeof(timestamp), "%d/%m/%C%y at %H:%M:%S", &mytime);
		buffer = "%file originally created on ";
		buffer += timestamp;     //  "dd/mm/yyyy at hh:mm:ss"
		inifile.push_back(buffer);
		writeIni();  // create a new .ini file using hard coded initial values

		std::cout << " cannot open BigIntCalculator.ini - create new file \n";
	}
	else {    /* read in .ini file */
		while (std::getline(iniStr, buffer)) {

			if (_strnicmp("yafu-path=", buffer.c_str(), 10) == 0) {
				YafuPath = buffer.substr(10); // copy path following '=' character
			}
			else if (_strnicmp("yafu-prog=", buffer.c_str(), 10) == 0) {
				yafuprog = buffer.substr(10); // copy path following '=' character
			}
			else if (_strnicmp("msieve-path=", buffer.c_str(), 12) == 0) {
				MsievePath = buffer.substr(12); // copy path following '=' character
			}
			else if (_strnicmp("msieve-prog=", buffer.c_str(), 12) == 0) {
				MsieveProg = buffer.substr(12); // copy path following '=' character
			}
			else if (_strnicmp("helpfile=", buffer.c_str(), 9) == 0) {
				helpFilePath = buffer.substr(9);
			}

			else inifile.push_back(buffer);  // save anything not recognised
		}
		iniStr.close();
	}
}


//search the docfile for the right entry specified by s
//just search for the heading, and print everything until
//the next heading is found
static void helpfunc(const std::string &s)
{
	FILE *doc;
	char str[1024];
	bool printtopic = false;
	std::string expr = "";
	char * newpathC;
	std::string newpath;

retry:
	//open the doc file and search for a matching topic
	errno_t ecode = fopen_s(&doc, helpFilePath.data(), "r");
	if (ecode != 0) {
		/* failed to open the help file*/
		char buffer[80];
		_strerror_s(buffer, sizeof(buffer), NULL); /* convert errno to a text messsage */
		fprintf_s(stderr, "fopen error: %s\n", buffer);
		fprintf_s(stderr, "help file not found\n");

		while (std::toupper(expr[0]) != 'Y') {
			std::cout << "Do you want to search for the help file? (Y/N) \n";
			getline(std::cin, expr);
			if (std::toupper(expr[0]) == 'N')
				return;
		}
		newpathC = getFileName("Text\0*.TXT\0\0", NULL);
		if (newpathC ==NULL) {
			std::cout << "command cancelled \n";
			return;
		}
		else {
			helpFilePath = newpathC; /* copy new path for doc file to permanent storage */
			writeIni();       /* update the .ini file*/
			goto retry;       /* let's try to open the file again */
		}
	}

	/* doc file has been opened successfully */
	if (verbose > 0)
		printf_s("searching for help on '%s'\n", s.data());

	/* exit this loop when reached EOF or the next topic after the one required 
	is reached */
	while (!feof(doc)) {
		
		char *rv = fgets(str, sizeof(str), doc);   //read a line
		if (rv == NULL)
			if (feof(doc)) {
				break;
			}
			else {
				char buffer[80];
				_strerror_s(buffer, sizeof(buffer), NULL); /* convert errno to a text messsage */
				fprintf_s(stderr, "fgets error: %s\n", buffer);
				break;
			}

		//is this a header?
		if ((str[0] == '[') && (str[strlen(str) - 2] == ']')) 	{

			if (printtopic)
				break;  /* we have reached the start of the next topic, so exit 
						   Only print 1 topic per help command */

			//does it match our topic?
			str[strlen(str) - 2] = '\0'; /* overwrite ']' with null */
			if (strstr(s.data(), str + 1) != NULL)
				/* we get a match if the topic between [ and ] is contained 
				anywhere in string s */
				printtopic = true;   /* we have found the required topic*/
		}
		else {  /* not a header line */
			if (printtopic)
				printf_s("%s", str);  /* print only if within the required topic */
		}
	}

	if (feof(doc))
		putchar('\n');   /* contrary to the POSIX standard, the last line of the file 
						    may not end with newline */
	fclose(doc);
	return;
}


/* check for commands. return 2 for exit, 1 for other command, 0 if not a command*/
static int processCmd(const std::string &command) {

	/* list of commands (static for efficiency) */
	const static std::vector<std::string>list =
	{ "EXIT", "SALIDA", "HELP", "AYUDA", "E", "S",  "F" , "N" , "X", "D",
	  "TEST", "MSIEVE", "YAFU", "V " };

	/* do 1st characters of text in command match anything in list? */
	int ix = 0;
	for (ix = 0; ix < list.size(); ix++) {
		auto len = list[ix].size();  /* get length of a string in list*/
		if (len > 1 && command.substr(0, len) == list[ix])
			break;  /* exit loop if match found */
		if (len == 1 && command == list[ix])
			break;   /* 1-letter commands with no parameters only match 1-letter input */
	}

	if (ix >= list.size())
		return 0;   /* not a command in list */

	switch (ix) {
	case 0:  /* exit*/
	case 1:  /* salida */
		return 2;  /* same command in 2 languages */
	case 2: /* Help */ {
			if (command.size() > 4)
				helpfunc(command.substr(5));
			else
				helpfunc("HELP");
			return 1;
		}
	case 3: /* ayuda */ {
			if (command.size() > 5)
				helpfunc(command.substr(6));
			else
				helpfunc("AYUDA");
			return 1;
		}
	case 4: /* E */
		{ lang = 0; return 1; }           // english
	case 5: /* S */
		{ lang = 1; return 1; }	          // spanish (Español)
	case 6: /* F */
		{ factorFlag = true; return 1; }  // do factorisation
	case 7: /* N */
		{ factorFlag = false; return 1; } // don't do factorisation
	case 8: /* X */
		{ hex = true; return 1; }         // hexadecimal output
	case 9: /* D */
		{ hex = false; return 1; }        // decimal output

	case 10: /* TEST */ {
			/* there are a number of comands that begin with TEST */
			if (command == "TEST") {
				doTests();         // do basic tests 
				return 1;
			}

			int ttype = command.data()[4] - '0';  /* for TEST2 ttype = 2, etc*/
			switch (ttype) {
			case 2: {
				if (command.size() > 5)
					doTests2(command.substr(6));         // do basic tests 
				else
					doTests2("");
				return 1;
			}
			case 3:
	#ifdef BIGNBR
				{
					doTests3();         // do basic tests 
					return 1;
				}
	#else
				return 0;
	#endif
			case 4: {
					doTests4();         // do R3 tests 
					return 1;
				}
			case 5: {
					doTests5();         // do YAFU tests 
					return 1;
				}
			case 6: {
					doTests6();         // do Msieve tests 
					return 1;
				}
			case 7: {
					/* Lucas-Lehmer test */
					if (command.size() > 5)
						doTests7(command.substr(6));
					else
						doTests7("");
					return 1;
				}
			default:
				return 0;  /* not a recognised command */
			}
		}

	case 11: /* MSIEVE */ {
			msieveParam(command);
			return 1;
		}
	case 12: /* YAFU */ {
			yafuParam(command);
			return 1;
		}
	case 13: /* V */ {
			/* will not throw an exception if input has fat finger syndrome.
			If no valid digits found, sets verbose to 0 */
			verbose = atoi(command.substr(1).data());
			std::cout << "verbose set to " << verbose << '\n';
			return 1;
		}

	default:
		return 0;   /* not a recognised command */
	}
}

int main(int argc, char *argv[]) {
	std::string expr;
	Znum Result;
	retCode rv;
	
	try {
		int version[4]; /* version info from .exe file (taken from .rc resource file) */

		hConsole = GetStdHandle(STD_OUTPUT_HANDLE);  // get handle for console window
		if (hConsole == INVALID_HANDLE_VALUE)
		{
			fprintf_s(stderr, "GetStdHandle failed with %d at line %d\n", GetLastError(), __LINE__);
			Beep(750, 1000);
			return EXIT_FAILURE;
		}

		std::string modified;  /* date & time program last modified */
		VersionInfo(argv[0], version, modified); /* get version info from .exe file */
		printf_s("%s Bigint calculator Version %d.%d \n", myTime(), version[0], version[1]);
		std::cout << "last modified on " << modified << '\n';


#ifdef __GNUC__
		printf("gcc version: %d.%d.%d\n", __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
		setlocale(LC_ALL, "en_GB.utf8");      // allows non-ascii characters to print
#endif

#ifdef _MSC_FULL_VER
/* For example, if the version number of the Visual C++ compiler is 15.00.20706.01, 
the _MSC_FULL_VER macro evaluates to 150020706 */
		long long ver = _MSC_FULL_VER;
		std::cout << "MSVC version: " << ver / 10000000 << '.';  // 1st 2 digits
		ver %= 10000000;                      // remove 1st 2 digits
		std::cout << ver / 100000 << '.' ;    // next 2 digits
		ver %= 100000;                        // remove next 2 digits
		std::cout << ver << '\n';             // last 5 digits
#endif

		auto lc = setlocale(LC_ALL, "en-EN");      // allows non-ascii characters to print
		printf("locale is now: %s\n", setlocale(LC_ALL, NULL));
		std::cout << "GMP version: " << __GNU_MP_VERSION << '.' << __GNU_MP_VERSION_MINOR
			<< '.' << __GNU_MP_VERSION_PATCHLEVEL << '\n';

#ifdef __MPIR_VERSION
		std::cout << "MPIR version: " << __MPIR_VERSION << '.' << __MPIR_VERSION_MINOR
			<< '.' << __MPIR_VERSION_PATCHLEVEL << '\n';
#endif

		processIni(argv[0]); // read .ini file if it exists

		while (true) {
			if (lang == 0) {
				printf("enter expression to be processed, or HELP, or EXIT\n");
			}
			else
				printf("ingrese la expresión para ser procesada, o AYUDA o SALIDA\n");

			PlaySound(TEXT("c:/Windows/Media/chimes.wav"), NULL, 
				SND_FILENAME | SND_NODEFAULT | SND_ASYNC | SND_NOSTOP);
			getline(std::cin, expr);    // expression may include spaces
			strToUpper(expr, expr);		// convert to UPPER CASE 

			// prevent the sleep idle time-out.
			SetThreadExecutionState(ES_CONTINUOUS | ES_SYSTEM_REQUIRED);

			int cmdCode = processCmd(expr);
			if (cmdCode == 2) break;    // EXIT command
			if (cmdCode == 1) {
				// Clear EXECUTION_STATE flags to allow the system to idle to sleep normally.
				SetThreadExecutionState(ES_CONTINUOUS);
				continue; // command has been fully processed
			}
			auto start = clock();	// used to measure execution time
			removeBlanks(expr);     // remove any spaces 

			rv = ComputeExpr(expr, Result); /* analyse expression, compute value*/

			if (rv != retCode::EXPR_OK) {
				textError(rv);   // invalid expression; print error message
			}
			else {
				std::cout << " = ";
				ShowLargeNumber(Result, 6, true, hex);   // print value of expression
				std::cout << '\n';				
				if (factorFlag) {
					doFactors(Result,false); /* factorise Result, calculate number of divisors etc */
					results.clear();  // get rid of unwanted results
				}
			}

			auto end = clock();   // measure amount of time used
			double elapsed = (double)end - start;
			PrintTimeUsed(elapsed, "time used = ");
			// Clear EXECUTION_STATE flags to allow the system to idle to sleep normally.
			SetThreadExecutionState(ES_CONTINUOUS);
		}

		// Clear EXECUTION_STATE flags to allow the system to idle to sleep normally.
		SetThreadExecutionState(ES_CONTINUOUS);

		return EXIT_SUCCESS;  // EXIT command entered
	}

	/* code below catches C++ 'throw' type exceptions */
	catch (const std::exception& e) {
		printf_s("\n*** standard exception caught, message '%s'\n", e.what());
		Beep(1200, 1000);              // sound at 1200 Hz for 1 second
		system("PAUSE");               // press any key to continue
		exit(EXIT_FAILURE);
	}

	catch (const char *str) {
		printf_s("\n*** Caught exception: <%s> \n", str);
		Beep(1200, 1000);              // sound at 1200 Hz for 1 second
		system("PAUSE");               // press any key to continue
		exit(EXIT_FAILURE);
	}

	catch (...) {
		// this executes if f() throws int or any other unrelated type
		// This catch block probably only would be executed under /EHa compiler option 
		/* most likely to be a SEH-type exception */
		printf_s("\n*** unknown exception ocurred\n");
		Beep(1200, 1000);              // sound at 1200 Hz for 1 second
		system("PAUSE");               // press any key to continue
		exit(EXIT_FAILURE);
	}
}
