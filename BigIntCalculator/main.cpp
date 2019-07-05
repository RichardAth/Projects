#include <iostream>
#include <string>
#include <vector>
#include <cctype>
#include <climits>
#include <locale>
#include <assert.h>
#include <Windows.h>
#include "factor.h"

//#define BIGNBR       // define to include bignbr tests 
#ifdef BIGNBR
#include "bignbr.h"
#include "bigint.h"
extern Znum zR, zR2, zNI, zN;
#endif

void msieveParam(std::string expupper);

#define PAREN_STACK_SIZE            100
int lang = 0;             // 0 English, 1 = Spanish
static int stackIndex=0; 
static int exprIndex;
bool hex = false;		// set true if output is in hex
bool factorFlag = true;
HANDLE hConsole;

/* list of operators, arranged in order of priority. Order is not exaxtly the
same as C or Python. */
enum class opCode {
	oper_power       =  0, 
	oper_unary_minus =  1,  // C and Python put unary minus above multiply, divide & modulus
	oper_multiply    =  2,
	oper_divide      =  3,
	oper_remainder   =  4,  // AKA modulus
	oper_comb        =  5,   // nCk, also known as binomial coefficient
	oper_plus        =  6,
	oper_minus       =  7,
	oper_shr         =  8,
	oper_shl         =  9,
	oper_not_greater = 10,
	oper_not_less    = 11,
	oper_greater     = 12,
	oper_less        = 13,
	oper_not_equal   = 14,
	oper_equal       = 15,
	oper_not         = 16,      // C and Python put bitwise not with unary minus
	oper_and         = 17,      // C and Python put AND before XOR before OR
	oper_xor         = 18,
	oper_or          = 19,
	oper_leftb       = 20,
	oper_fact,					// !   factorial 
	oper_dfact,					// !!  double factorial
	oper_prim,					// #   primorial
};

/* list of operator priority values. lower value = higher priority. Note: this 
order is not the same as C or Python */
const static int operPrio[] =
{
	0,				  // Power
	1,                // Unary minus.  (C and Python put unary minus above multiply, divide, remainder)
	2, 2, 2,          // Multiply, divide and remainder.
	3,                // combination
	4, 4,             // Plus and minus.
	5, 5,             // Shift right and left.
	6, 6, 6, 6,       // four comparison operators (equal, greater, less, etc.)
	7, 7,             // == and !=
	8,                // NOT.   (C and Python put bitwise not with unary minus)
	9,                // AND, (C and Python but AND before XOR before OR)
	10,               // XOR
	11,               // OR
};

/* error and return codes */
enum class retCode
{
	//EXPR_NUMBER_TOO_LOW,
	//EXPR_NUMBER_TOO_HIGH,
	EXPR_INTERM_TOO_HIGH = -100,
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

bool *primeFlags = NULL;
unsigned long long *primeList = NULL;
unsigned int prime_list_count = 0;
unsigned long long int primeListMax = 0;
const unsigned long long max_prime = 1000000007;  // arbitrary limit 10^9,


/* function declarations, only for functions that have forward references */
static retCode ComputeExpr(const std::string &expr, Znum &ExpressionResult);
long long MulPrToLong(const Znum &x);
static bool getBit(const unsigned long long int x, bool array[]);
void generatePrimes(unsigned long long int max_val);

void biperm(int n, mpz_t &result);   // declaration for external function
 
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

static void ShowLargeNumber(const Znum &Bi_Nbr, int digitsInGroup, bool size, bool hex) {
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

/* calculate number of divisors of n, given its list of prime factors.
NB for n=1 function returns 2, but correct value would be 1. 
From just the list of exponents we can't distinguish n=1 from n=prime */
static Znum NoOfDivisorsF(const std::vector <zFactors> exponents) {
	Znum result = 1;

	for (auto i: exponents) {
		result *= i.exponent + 1;
	}
	return result;
}

// sum of divisors is the product of (p^(e+1)-1)/(p-1) where p=prime factor and e=exponent.
static Znum SumOfDivisors( const std::vector <zFactors> primes) {
	Znum result = 1, term;

	if (primes[0].Factor == 1)
		return 1;		// special case: if original number is 1
	for (size_t i = 0; i < primes.size(); i++) {
		mpz_pow_ui(ZT(term), ZT(primes[i].Factor), primes[i].exponent + 1);  // p^(e+1)
		term = (term-1)/(primes[i].Factor-1);	                        // (p^(e+1)-1)/(p-1)
		result *= term;
	}
	return result;
}

// Find Euler's Totient as the product of p^(e-1)*(p-1) where p=prime and e=exponent.
static Znum Totient(const std::vector <zFactors> primes) {
	Znum result = 1, term;

	for (size_t i = 0; i < primes.size(); i++) {
		mpz_pow_ui(ZT(term), ZT(primes[i].Factor), primes[i].exponent - 1);  // p^(e-1)
		term = term* (primes[i].Factor- 1);	                        // (p^(e-1)-1)*(p-1)
		result *= term;
	}
	return result;
}

/* For any positive integer n, define μ(n) as the sum of the primitive nth roots of unity. 
   It has values in {−1, 0, 1} depending on the factorization of n into prime factors:
   μ(n) = 1 if n is a square-free positive integer with an even number of prime factors.
   μ(n) = −1 if n is a square-free positive integer with an odd number of prime factors.
   μ(n) = 0 if n has a squared prime factor. */
static int mobius(const std::vector <zFactors> exponents) {
	
	for (auto ex : exponents) {
		if (ex.exponent > 1)
			return 0;		// n is not square-free
	}
	auto result = exponents.size();
	if ((result & 1) == 1)
		return -1;  // odd nuber of  prime factors
	else
		return 1;    // even number of prime factors 
}

/* calculate Euler's totient for n as the product of p^(e-1)*(p-1) 
where p=prime factor and e=exponent.*/
static Znum ComputeTotient(const Znum &n) {
	std::vector <zFactors> factorlist;

	if (n == 1)
		return 1;
	auto rv = factorise(n, factorlist, nullptr);
	if (rv && !factorlist.empty()) {
		auto divisors = Totient(factorlist);
		return divisors;
	}
	else return 0;
}

/* calculate number of divisors of n, from its list of prime factors. */
static Znum ComputeNumDivs(const Znum &n) {
	std::vector <zFactors> factorlist;

	if (n == 1)
		return 1;  // 1 only has one divisor. NoOfDivisors can't handle that case
	auto rv = factorise(n, factorlist, nullptr);
	if (rv && !factorlist.empty()) {
		auto divisors = NoOfDivisorsF(factorlist);
		return divisors;
	}
	else return 0;
}

/* sum of divisors is the product of (p^(e+1)-1)/(p-1) 
 where p=prime factor and e=exponent. */
static Znum ComputeSumDivs(const Znum &n) {
	std::vector <zFactors> factorlist;

	if (n == 1)
		return 1;   // 1 only has 1 divisor. 
	auto rv = factorise(n, factorlist, nullptr);
	if (rv && !factorlist.empty()) {
		auto divisors = SumOfDivisors(factorlist);
		return divisors;
	}
	else return 0;
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
static long long ComputeNumDigits(const Znum &n, const Znum &radix)
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
		if (primeFlags == NULL || primeListMax < max_prime) {
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

/* ConcatFact(m,n): Concatenates the prime factors (base 10) of num according to the mode
mode	Order of factors	Repeated factors
0		Ascending			No
1		Descending			No
2		Ascending			Yes
3		Descending			Yes
*/
static Znum  concatFact(const Znum &mode, const Znum &num) {
	std::vector <zFactors> factorlist;
	std::string result;
	Znum rvalue = 0;
	const bool descending = ((mode & 1) == 1);
	const bool repeat = ((mode & 2) == 2);
	char *buffer = NULL;

	/* get factors of num */
	auto rv = factorise(num, factorlist, nullptr);
	if (rv && !factorlist.empty()) {
		if (descending)   /* start with largest factor */
			for (ptrdiff_t i = factorlist.size()-1; i >=0; i--) {
				buffer = mpz_get_str(NULL, 10, ZT(factorlist[i].Factor));
				if (!repeat)
					result += buffer; // concatenate factor
				else
					for (int j = 1; j <= factorlist[i].exponent; j++)
						result += buffer; // concatenate factor
				free(buffer);
				buffer = NULL;
			}
		else  /* start with smallest factor */
			for (size_t i = 0; i < factorlist.size(); i++) {
				buffer = mpz_get_str(NULL, 10, ZT(factorlist[i].Factor));
				if (!repeat)
					result += buffer;  // concatenate factor
				else
					for (int j = 1; j <= factorlist[i].exponent; j++)
						result += buffer;  // concatenate factor
				free(buffer);
				buffer = NULL;
			}
		mpz_set_str(ZT(rvalue), result.data(), 10); /* convert back from a string to a number */
		return rvalue;
	}
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

	std::vector <zFactors> factorlist;
	Znum b = 1;

	/* get factors of num */
	auto rv = factorise(num, factorlist, nullptr);
	for (size_t i = 0; i < factorlist.size(); i++) {
		if (factorlist[i].Factor <= 2)
			continue; // ignore factor 1 or 2
		if (factorlist[i].Factor % 4 == 3) { /* p = 4k+3? */
			if (factorlist[i].exponent % 2 == 1) /* exponent is odd?*/
				return 0;
		}
		else { 		/* p = 4k + 1 */
				b *= (factorlist[i].exponent + 1);
			}
	}
	return b * 4;
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
	std::vector <zFactors> factorlist;
	zFactors temp;
	auto rv = factorise(num, factorlist, nullptr);
	sqf.clear();

	for (size_t i = 0; i < factorlist.size(); i++) {
		if (factorlist[i].exponent >= 2) {
			/* copy exponent value, rounded down to multiple of 2*/
			temp.exponent = factorlist[i].exponent - (factorlist[i].exponent & 1);
			temp.Factor = factorlist[i].Factor;
			sqf.push_back(temp);
			factorlist[i].exponent -= temp.exponent;
		}
	}
	

	sq = 1;
	Znum x;
	for (size_t i = 0; i < sqf.size(); i++) {
		mpz_pow_ui(ZT(x), ZT(sqf[i].Factor), sqf[i].exponent);
		sq *= x;
			/*for (int j = 1; j <= sqf[i].exponent; j++)
				sq *= sqf[i].Factor;*/
	}
	num /= sq;
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
	std::vector <zFactors> sqf;

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
const static std::array <struct functions, 22> functionList{
	"GCD",       2,  fn_Code::fn_gcd,			// name, number of parameters, code
	"MODPOW",    3,  fn_Code::fn_modpow,
	"MODINV",    2,  fn_Code::fn_modinv,
	"TOTIENT",   1,  fn_Code::fn_totient,
	"NUMDIVS",   1,  fn_Code::fn_numdivs,
	"SUMDIVS",   1,  fn_Code::fn_sumdivs,
	"SUMDIGITS", 2,  fn_Code::fn_sumdigits,
	"NUMDIGITS", 2,  fn_Code::fn_numdigits,
	"REVDIGITS", 2,  fn_Code::fn_revdigits,
	"ISPRIME",   1,	 fn_Code::fn_isprime,
	"FactConcat",2,  fn_Code::fn_concatfact,     // FactConcat must come before F
	"F",         1,  fn_Code::fn_fib,			// fibonacci
	"LE",		 2,  fn_Code::fn_legendre,
	"L",         1,  fn_Code::fn_luc,			// Lucas Number
	"PI",		 1,  fn_Code::fn_primePi,		// prime-counting function. PI must come before P
	"P",         1,  fn_Code::fn_part,			// number of partitions
	"N",         1,  fn_Code::fn_np,				// next prime
	"B",         1,  fn_Code::fn_pp,				// previous prime
	"R2",		 1,  fn_Code::fn_r2,			// number of ways n can be expressed as sum of 2 primes
	"R3",        1,  fn_Code::fn_r3,
	"JA",		 2,  fn_Code:: fn_jacobi,
	"KR",		 2,  fn_Code::fn_kronecker
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
		if (p1 < 0 || p1 >= 60000) {
			return retCode::EXPR_INVALID_PARAM;  // note: biperm is limited to values <= 60,000
		}
		int temp = (int)MulPrToLong(p1);
		biperm(temp, ZT(result));   // calculate number of partitions
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
		result = concatFact(p1, p2);
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
	
	default:
		abort();		// if we ever get here we have a problem
	}
	return retcode;
}

/* Check whether expr matches any function name. If so ,compute value(s) of 
 parameter(s), compute function value, set Gfcode and return EXPR_OK. 
 If not, return value 1 (EXPR_FAIL).
Updates global var exprIndex */
static retCode func(const std::string &expr, const bool leftNumberFlag,
 	 int exprIndexLocal, Znum &FnValue) {
	retCode retcode;
	fn_Code fnCode = fn_Code::fn_invalid;
	long long ix;
	int NoOfArgs = 0;
	Znum args[4];     // parameter values (only up to 3 values used)
	const char *ptrExpr = expr.data() + exprIndexLocal;

	/* try to match function name */
	for (ix = 0; ix < (ptrdiff_t)functionList.size(); ix++) {
		if (_strnicmp(ptrExpr, functionList[ix].fname, strlen(functionList[ix].fname)) == 0) {
			fnCode = functionList[ix].fCode;
			NoOfArgs = functionList[ix].NoOfParams;
			break;  // found a match for the function name 
		}
	}

	if (fnCode == fn_Code::fn_invalid)
		return retCode::EXPR_FAIL;  // no function name found

	assert(NoOfArgs <= sizeof(args) / sizeof(args[0]));
		
	exprIndexLocal += (int)strlen(functionList[ix].fname); // move exprIndex past function name
	if (leftNumberFlag)
	{
		return retCode::EXPR_SYNTAX_ERROR;
	}
	if (exprIndexLocal >= expr.length() || expr[exprIndexLocal++] != '(')
	{
		return retCode::EXPR_SYNTAX_ERROR;  // no opening bracket after function name
	}

	int Bracketdepth = 0;
	for (ix = exprIndexLocal-1; ix < (long long)expr.length(); ix++) {
		if (expr[ix] == '(') Bracketdepth++;
		if (expr[ix] == ')') Bracketdepth--;
		if (Bracketdepth == 0) break;  // found matching closing bracket
	}
	if (Bracketdepth != 0)
		return retCode::EXPR_PAREN_MISMATCH;
	long long exprLength = ix - exprIndexLocal + 1;  // length of text up to closing bracket.

	/* now get value(s) of function arguments(s) */
	for (int index = 0; index < NoOfArgs; index++)
	{
		char compareChar;
		/* substring in function call below includes all characters from the parameter now
		being processed up to and including the closing bracket. */
		retcode = ComputeExpr(expr.substr(exprIndexLocal, exprLength), args[index]);
		if (retcode != retCode::EXPR_OK) {
			return retcode;   // unable to evaluate parameter
		}
		exprIndexLocal += exprIndex;  // move exprIndexLocal past expression just evaluated.
		exprLength -= exprIndex;
		/* if there's more than 1 parameter they are separated with commas. 
		The last one is followed by ')'*/
		compareChar = (index == NoOfArgs - 1 ? ')' : ',');
		if (exprIndexLocal >= expr.length() || expr[exprIndexLocal++] != compareChar)
		{
			return retCode::EXPR_SYNTAX_ERROR;  // no comma or ) when needed
		}
		exprLength--;		// subtract length of comma or closing bracket
	}

	exprIndex = exprIndexLocal;  // update global index to char after )
	retcode = ComputeFunc(fnCode, args[0], args[1], args[2], FnValue);  // evaluate function
	return retcode;
}

/* SHL: Shift left the number of bits specified on the right operand. If the
number of bits to shift is negative, this is actually a right shift. If result would
be too large an error is reported. */
static retCode ShiftLeft(const Znum first, const Znum bits, Znum &result) {
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
static bool getBit(const unsigned long long int x, bool array[])
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
	unsigned long long int numsave, count = 1, num = 3;
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
	printf("  prime %9lld is %11lld\n", count, numsave);
	primeList[count] = ULLONG_MAX;		// set end marker
	prime_list_count = (unsigned int)count;
	primeListMax = primeList[count - 1];
	return;
}

/* process one operator with 1 or 2 operands on stack. 
NOT and unary minus have 1 operand. All the others have two.
some operators can genererate an error condition e.g. EXPR_DIVIDE_BY_ZERO
otherwise return EXPR_OK. */
static retCode ComputeSubExpr(opCode stackOper, const Znum &firstArg,
	const Znum &secondArg, Znum &result)
{
	
	switch (stackOper)
	{
	case opCode::oper_comb: {  // calculate nCk AKA binomial coefficient
		if (secondArg > INT_MAX)
		return retCode::EXPR_INTERM_TOO_HIGH;
		if (secondArg < INT_MIN)
			return retCode::EXPR_INVALID_PARAM;
		long long k = MulPrToLong(secondArg);
		mpz_bin_ui(ZT(result), ZT(firstArg), k);
		return retCode::EXPR_OK;
	}
	case opCode::oper_plus: {
		result = firstArg + secondArg; 
		return retCode::EXPR_OK;
	}
	case opCode::oper_minus: {
		result = firstArg - secondArg; 
		return retCode::EXPR_OK;
	}
	case opCode::oper_unary_minus: {
		result = -secondArg; 
		return retCode::EXPR_OK;
	}
	case opCode::oper_divide: {
		if (secondArg == 0)
			return retCode::EXPR_DIVIDE_BY_ZERO;  // result would be infinity
		result = firstArg / secondArg; 
		return retCode::EXPR_OK;
	}
	case opCode::oper_multiply: {
		auto resultsize = NoOfBits(firstArg) + NoOfBits(secondArg);
		if (resultsize > 66439)  // more than 66439 bits -> more than 20,000 decimal digits
			return retCode::EXPR_INTERM_TOO_HIGH;
		result = firstArg * secondArg;
		return retCode::EXPR_OK;
	}
	case opCode::oper_remainder: {
		if (secondArg == 0)
			return retCode::EXPR_DIVIDE_BY_ZERO;  // result would be infinity
		result = firstArg % secondArg;   //return BigIntRemainder(firstArg, secondArg, result);
		return retCode::EXPR_OK;
	}
	case opCode::oper_power: {
		if (secondArg > INT_MAX)
			return retCode::EXPR_EXPONENT_TOO_LARGE;
		if (secondArg < 0)
			return retCode::EXPR_EXPONENT_NEGATIVE;
		long long exp = MulPrToLong(secondArg);
		auto resultsize = (NoOfBits(firstArg)-1)* exp;  // estimate number of bits for result
		if (resultsize > 66439)  // more than 66439 bits -> more than 20,000 decimal digits
			return retCode::EXPR_INTERM_TOO_HIGH;
		mpz_pow_ui(ZT(result), ZT(firstArg), exp);
		return retCode::EXPR_OK;
	}
	case opCode::oper_equal: {
		if (firstArg == secondArg)
			result = -1;
		else
			result = 0;
		return retCode::EXPR_OK;
	}
	case opCode::oper_not_equal: {
		if (firstArg != secondArg)
			result = -1;
		else
			result = 0;
		return retCode::EXPR_OK;
	}
	case opCode::oper_greater: {
		if (firstArg > secondArg)
			result = -1;
		else
			result = 0;
		return retCode::EXPR_OK;
	}
	case opCode::oper_not_greater: {
		if (firstArg <= secondArg)
			result = -1;
		else
			result = 0;
		return retCode::EXPR_OK;
	}
	case opCode::oper_less: {
		if (firstArg < secondArg)
			result = -1;
		else
			result = 0;
		return retCode::EXPR_OK;
	}
	case opCode::oper_not_less: {
		if (firstArg <= secondArg)
			result = -1;
		else
			result = 0;
		return retCode::EXPR_OK;
	}
	case opCode::oper_shl: {
		return ShiftLeft(firstArg, secondArg, result);
	}
	case opCode::oper_shr: {
		// invert sign of shift
		return ShiftLeft(firstArg, -secondArg, result);
	}
	case opCode::oper_not: {   // Perform binary NOT 
		//result = -1 - secondArg;  // assumes 2s complement binary numbers
		mpz_com(ZT(result), ZT(secondArg));
		return retCode::EXPR_OK;
	}
	case opCode::oper_and: {  // Perform binary AND.
		mpz_and(ZT(result), ZT(firstArg), ZT(secondArg));
		return retCode::EXPR_OK;
	}
	case opCode::oper_or: {   // Perform binary OR.
		mpz_ior(ZT(result), ZT(firstArg), ZT(secondArg));
		return retCode::EXPR_OK;
	}
	case opCode::oper_xor: {   // Perform binary XOR.
		mpz_xor(ZT(result), ZT(firstArg), ZT(secondArg));
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
/* list of operators that have format <expression> <operator> <expression> 
with corresponding codes and priority. However, operator NOT (bitwise negation) 
is in this list although it only has one operand. 
Operators ! (factorial) !! (double factorial) and # (primorial) are not in this 
list because the convention is that they follow the number or expression they 
operate on. */
const static struct oper_list operators[]  {
	{"C",	 opCode::oper_comb,	       3},
	{ "^",   opCode::oper_power,       0},
	{ "**",  opCode::oper_power,       0},     // can use ^ or ** for exponent
	{ "*",   opCode::oper_multiply,    2},
	{ "/",   opCode::oper_divide,      2},
	{ "%",   opCode::oper_remainder,   2},
	{ "+",   opCode::oper_plus,        4},
	{ "-",   opCode::oper_minus,       4},
	{ "SHL", opCode::oper_shl,         5},
	{ "<<",  opCode::oper_shl,         5},     // can use << or SHL for left shift
	{ "SHR", opCode::oper_shr,         5},
	{ ">>",  opCode::oper_shr,         5},     // can use SHR or >> for right shift
	{ "<=",  opCode::oper_not_greater, 6},
	{ ">=",  opCode::oper_not_less,    6},
	{ ">",   opCode::oper_greater,     6},	  // to avoid mismatches > and < must come after >> and <<
	{ "<",   opCode::oper_less,        6},
	{ "!=",  opCode::oper_not_equal,   7},
	{ "==",  opCode::oper_equal,       7},
	{ "NOT", opCode::oper_not,         8},      // bitwise NOT
	{ "AND", opCode::oper_and,         9},      // bitwise AND
	{ "OR",  opCode::oper_or,         11},      // bitwise OR
	{ "XOR", opCode::oper_xor,        10}   };  // bitwise exclusive or

/* search for operators e.g. '+', '*'  that have format <expression> <operator> <expression>
or <operator> <expression>. return index of operator or -1 */
static void operSearch(const std::string &expr, long long &opcode) {
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
C			binomial coefficient	0
** or ^		Exponentiate			1
*			multiply				2
/			divide					2
%			modulus (remainder)		2
-			unary minus				3
+			add						4
-			subtract				4
n SHL b		shift n left b bits		5
n SHR b		shift n right b bits	5

comparison operators
>			return zero for false 	6
>=			and -1 for true			6
<									6
<=									6
==									6
!=			not equal				6

NOT			bitwise operators		7
AND									8
OR									8
XOR									8

n!			factorial
n!!			double factorial
p#			primorial(product of all primes less or equal than p

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

/* analyse an expression supplied as a string. This function uses indirect recursion;
it calls func to detect a function name and analyse the parameters of a function, which 
in turn calls computeExpr again to evaluate each parameter to the function. 
The numerical value of the expression is returned in ExpressionResult  */
static retCode ComputeExpr(const std::string &expr, Znum &ExpressionResult) {
	static int depth=0;  // measures depth of recursion, increment on entry, decrement on return 
	bool leftNumberFlag = false;
	retCode SubExprResult;
	int startStackIndex = stackIndex;    // for nested calls this is > zero

	static Znum stackValues[PAREN_STACK_SIZE];
	static opCode stackOperators[PAREN_STACK_SIZE];
	/* stackValues is accessed by the following macros*/
#define StackTopVal stackValues[stackIndex]
#define StackP1Val stackValues[stackIndex + 1]
#define StackM1Val stackValues[stackIndex - 1]


	if (expr.empty()) {
		return retCode::EXPR_INVALID_PARAM;
	}

	exprIndex = 0;

	depth++;		// measure recursion depth
	/* exit from loop below when end of expr or a , or ) is reached 
	if an error occurs, return error code immediately to caller */
	while (exprIndex < expr.length())
	{
		int charValue;

		opCode operCode;
		long long opIndex;
		int exprLength;
		int exprIndexAux;
		Znum factorial;
		retCode retcode;

		exprLength = (int)expr.length();
		charValue = toupper(expr[exprIndex]);  // get next character of the expression

		/* find operator or evaluate number or evaluate function */

		/* first look for operators. */
		operSearch(expr.substr(exprIndex), opIndex);
		if (opIndex != -1) {
			/* found operator e.g. + that has format <expression> <operator> <expression>*/
			operCode= operators[opIndex].operCode;
			exprIndex += (int)strlen(operators[opIndex].oper);  // move index to next char after operator

			if ((operCode == opCode::oper_plus || operCode == opCode::oper_minus) && leftNumberFlag == false)
			{                    // Unary plus/minus operator
				if (operCode == opCode::oper_plus)
				{
					continue;  // process more of expr
				}
				else
				{
					if (stackIndex > startStackIndex && 
						stackOperators[stackIndex - 1] == opCode::oper_unary_minus)
					{
						stackIndex--;
						continue;   // process more of expr
					}
					if (stackIndex >= PAREN_STACK_SIZE)
					{
						depth--;				// adjust call depth
						return retCode::EXPR_TOO_MANY_PAREN;
					}
					stackOperators[stackIndex++] = opCode::oper_unary_minus; /* Unary minus */
					continue;    // process more of expr
				}
			}
			if ((leftNumberFlag == false) != (operCode == opCode::oper_not))
			{     // Missing left operator if operator is not NOT or
					// extra left operator if operator is NOT.
				depth--;				// adjust call depth
				return retCode::EXPR_SYNTAX_ERROR;
			}
			if (operCode != opCode::oper_power) {
			 // Power operator has right associativity.
				while (stackIndex > startStackIndex) {
					/* check whether we can process any stacked values and operators*/
					if (stackOperators[stackIndex - 1] == opCode::oper_leftb)
						break;
					auto p1 = operPrio[(int)stackOperators[stackIndex - 1]];
					auto p2 = operPrio[(int)operCode];  // get priorities of stacked and current operators
					if (p1 > p2)    
						break;    // stacked operator has lower priority
					if ((SubExprResult = 
						ComputeSubExpr(stackOperators[stackIndex-1],
							StackM1Val, StackTopVal, StackM1Val)) != retCode::EXPR_OK)
					{
						depth--;				// adjust call depth
						return SubExprResult;   // failed! return error code
					}
					stackIndex--;
				}
			}
			stackOperators[stackIndex++] = operCode;
			leftNumberFlag = false;
			continue;     // process more of expr
		}
		else if (charValue == '!' && expr[exprIndex + 1] == '!') { // Calculating  double factorial.
			if (leftNumberFlag == false)
			{
				depth--;				// adjust call depth
				return retCode::EXPR_SYNTAX_ERROR;
			}
			if (StackTopVal > 11081)
			{
				depth--;				// adjust call depth
				return retCode::EXPR_INTERM_TOO_HIGH;
			}
			if (StackTopVal < 0)
			{
				depth--;				// adjust call depth
				return retCode::EXPR_INVALID_PARAM;
			}

			long long temp = llabs(MulPrToLong(StackTopVal));
			mpz_2fac_ui(ZT(factorial), temp);  // get double factorial
			StackTopVal = factorial;
			exprIndex+= 2;
			continue;      // process more of expr
		}
		else if (charValue == '!') {           // Calculating factorial.
			if (leftNumberFlag == false)
			{
				depth--;				// adjust call depth
				return retCode::EXPR_SYNTAX_ERROR;
			}
			if (StackTopVal > 5984 )
			{
				depth--;				// adjust call depth
				return retCode::EXPR_INTERM_TOO_HIGH;
			}
			if (StackTopVal < 0)
			{
				depth--;				// adjust call depth
				return retCode::EXPR_INVALID_PARAM;
			}
			long long temp = llabs(MulPrToLong(StackTopVal));
			mpz_fac_ui(ZT(factorial), temp);  // get factorial
			StackTopVal = factorial;
			exprIndex++;
			continue;      // process more of expr
		}
		else if (charValue == '#') {           // Calculating primorial.
			if (leftNumberFlag == false)
			{
				depth--;				// adjust call depth
				return retCode::EXPR_SYNTAX_ERROR;
			}
			if (StackTopVal < 0 || StackTopVal > 46340)
			{
				depth--;				// adjust call depth
				return retCode::EXPR_INTERM_TOO_HIGH;
			}
			long long temp = MulPrToLong(StackTopVal);
			mpz_primorial_ui(ZT(factorial), temp);
			StackTopVal = factorial;
			exprIndex++;
			continue;      // process more of expr
		}
		else   
			/* look for function name. If found calculate the value, leave it on   
			the stack and update exprIndex past ')'. If no name found return a +ve code */
			if ((retcode = func(expr, leftNumberFlag, exprIndex, StackTopVal)) 
				<= retCode::EXPR_OK) {
				if (retcode != retCode::EXPR_OK) {
					depth--;				// adjust call depth
					return retcode;  // error processing function parameters
				}
				
				leftNumberFlag = true;
				continue; // contine procesing rest of expr
			}
		/* drop through to here only if retcode is +ve (no name found) */
		else if (charValue == '(')
		{
			if (leftNumberFlag == true)	{
				depth--;				// adjust call depth
				return retCode::EXPR_SYNTAX_ERROR;
			}
			if (stackIndex >= PAREN_STACK_SIZE)	{
				depth--;				// adjust call depth
				return retCode::EXPR_TOO_MANY_PAREN;
			}
			stackOperators[stackIndex++] = opCode::oper_leftb;
			exprIndex++;
			continue;         // process more of expr
		}
		else if (charValue == ')' || charValue == ',')
		{
			if (leftNumberFlag == false)
			{
				depth--;				// adjust call depth
				return retCode::EXPR_SYNTAX_ERROR;
			}
			while (stackIndex > startStackIndex &&
				stackOperators[stackIndex - 1] != opCode::oper_leftb)
			{  // compute value(s) of sub-expression(s) following left bracket
				if ((SubExprResult = ComputeSubExpr(stackOperators[stackIndex-1],
					StackM1Val, StackTopVal, StackM1Val)) != retCode::EXPR_OK)
				{
					depth--;				// adjust call depth
					return SubExprResult;  // return error code
				}
				stackIndex--;
			}

			if (stackIndex == startStackIndex)
			{
				break;          // jump out of while loop if done
			}
			if (charValue == ',')
			{
				depth--;				// adjust call depth
				return retCode::EXPR_PAREN_MISMATCH;
			}
			stackIndex--;    /* Discard ')' */
			StackTopVal = StackP1Val;
			leftNumberFlag = 1;
			exprIndex++;
			continue;         // process more of expr
		}
		else if (charValue >= '0' && charValue <= '9') 	{
			/* convert number from ascii to Znum */
			exprIndexAux = exprIndex;
			if (charValue == '0' && exprIndexAux < exprLength - 2 &&
				toupper(expr[exprIndexAux + 1]) == 'X')
			{  // hexadecimal
				std::vector<char> digits;
				exprIndexAux++;
				while (exprIndexAux < exprLength - 1)
				{
					charValue = expr[exprIndexAux + 1];
					if ((charValue >= '0' && charValue <= '9') ||
						(charValue >= 'A' && charValue <= 'F') ||
						(charValue >= 'a' && charValue <= 'f'))
					{
						exprIndexAux++;
						digits.push_back(charValue);
					}
					else
					{
						break;   // jump out of inner while loop
					}
				}
				digits.push_back('\0');   // null terminate string
				mpz_set_str(ZT(StackTopVal), digits.data(), 16);
				exprIndex = exprIndexAux+1;
			}
			else
			{                   // Decimal number.
				std::vector<char> digits;
				while (exprIndexAux < exprLength)
				{  // find position of last digit
					charValue = expr[exprIndexAux];
					if (charValue >= '0' && charValue <= '9')
					{
						exprIndexAux++;
						digits.push_back(charValue);
					}
					else
					{
						break;    // jump out of inner while loop
					}
				}
				// Generate big integer from decimal number
				digits.push_back('\0');   // null terminate string
				mpz_set_str(ZT(StackTopVal), digits.data(), 10);
				exprIndex = exprIndexAux;
			}
			leftNumberFlag = true;
			continue;      // process more of expr
		}

		/* If we drop through to here, we've got something we can't understand */
		depth--;				// adjust call depth
		return retCode::EXPR_SYNTAX_ERROR;
	}                              /* end while */

	if (leftNumberFlag == false)
	{
		depth--;				// adjust call depth
		return retCode::EXPR_SYNTAX_ERROR;
	}
	while (stackIndex > startStackIndex && stackOperators[stackIndex - 1] != opCode::oper_leftb)
	{
		if ((SubExprResult = ComputeSubExpr(stackOperators[stackIndex-1],
			StackM1Val, StackTopVal, StackM1Val)) != retCode::EXPR_OK)
		{
			depth--;				// adjust call depth
			return SubExprResult;   // return error code
		}
		stackIndex--;
	}
	if (stackIndex != startStackIndex) 	{
		depth--;				// adjust call depth
		return retCode::EXPR_PAREN_MISMATCH;
	}
	ExpressionResult = StackTopVal;
	depth--;				// adjust call depth
	if (depth > 0 || exprIndex >= expr.size())
		return retCode::EXPR_OK;
	else
		return retCode::EXPR_SYNTAX_ERROR;  // exit from 1st call but, still have unprocessed text.
}
#undef StackTopVal 
#undef StackP1Val 
#undef StackM1Val 

/* translate error code to text and output it*/
static void textError(retCode rc)
{
	/*
	error codes currently used:

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
	//case EXPR_NUMBER_TOO_HIGH:
	//	std::cout << (lang ? "Número muy grande \n" :
	//		"Number too high \n");
		//break;
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
		if (lang)
		{
			std::cout << ( "Error de sintaxis\n");
		}
		else
		{
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

/* convert s to UPPER CASE in d. s and d must not be the same i.e. will
not do in-place conversion */
void strToUpper(const std::string &s, std::string &d) {
	d.clear();
	d.reserve(s.size());  // may be more efficient if s is very large
	for (auto c : s) {
		d.push_back(toupper(c));
	}
}

/* remove spaces, tabs, etc  from msg */
void removeBlanks(std::string &msg) {
	for (size_t ix = 0; ix < msg.size(); ix++) {
		if ((unsigned char)msg[ix] <= 0x7f && isspace(msg[ix])) {     // look for spaces, tabs, etc
			msg.erase(ix, 1);      // remove space character
			ix--;  // adjust index to take account of removed blank
		}
	}
}

/* factorise Result, calculate number of divisors etc and print results */
static void doFactors(const Znum &Result, bool test) {
	std::vector <zFactors> factorlist;
	Znum Quad[4];
	clock_t start;
	if (test)
		start = clock();	// used to measure execution time

	/* call DA´s magic function to factorise Result */
	bool rv = factorise(Result, factorlist, Quad);
	if (rv && !factorlist.empty()) {
		if (factorlist.size() > 1 || factorlist[0].exponent > 1) {
			/* print factor list */
			std::cout << " = ";
			if (Result < 0)
				std::cout << "-";
			for (size_t i = 0; i < factorlist.size(); i++) {
				if (i > 0)
					std::cout << " * ";
				ShowLargeNumber(factorlist[i].Factor, 6, false, false);
				if (factorlist[i].exponent > 1)
					std::cout << "^" << factorlist[i].exponent;
			}
		}
		else
			std::cout << " is prime";  //number has only 1 factor

		if (abs(Result) != 1) {
			auto divisors = NoOfDivisorsF(factorlist);
			std::cout << "\nNumber of Divisors = ";
			ShowLargeNumber(divisors, 6, false, false);
		}
		else 
			std::cout << "\nNumber of Divisors = 1";   // treat n=1 as special case

		auto divisors = SumOfDivisors(factorlist);
		std::cout << "\nSum of Divisors    = ";
		ShowLargeNumber(divisors, 6, false, false);
		divisors = Totient(factorlist);
		std::cout << "\nTotient            = ";
		ShowLargeNumber(divisors, 6, false, false);
		if (Result > 0) {
			auto mob = mobius(factorlist);  // mobius only defined for +ve integers
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
			c++;
		}
		c = 'a';
		for (auto q : Quad) {
			if (q == 0)
				break;
			std::cout << "\n" << c << "= " << q;
			c++;  // change a to b, b to c, etc
		}
		std::cout << "\n";
		if (factorlist.size() > 0) {
			std::cout << "found by";
			if (counters.tdiv > 0)
				std::cout << " trial division: " << counters.tdiv;
			if (counters.prho > 0)
				std::cout << " Pollard-rho: " << counters.prho;
			if (counters.pm1 > 0)
				std::cout << " power +/- 1: " << counters.pm1;
			if (counters.ecm > 0)
				std::cout << " elliptic curve: " << counters.ecm;
			if (counters.siqs > 0)
				std::cout << " SIQS: " << counters.siqs;
			if (counters.msieve > 0)
				std::cout << " Msieve: " << counters.msieve;
			if (counters.carm > 0)
				std::cout << " Carmichael: " << counters.carm;
			if (counters.leh > 0)
				std::cout << " Lehman: " << counters.leh;
			std::cout << '\n';
		}
		if (test) {
			Znum result = 1;
			for (size_t i = 0; i < factorlist.size(); i++)
				for (int j = 1; j <= factorlist[i].exponent; j++)
					result *= factorlist[i].Factor;
			if (result != Result) {
				std::cout << "Factors expected value " << Result << " actual value " << result << '\n';
				Beep(750, 1000);
			}
			result = Quad[0] * Quad[0] + Quad[1] * Quad[1] + Quad[2] * Quad[2] + Quad[3] * Quad[3];
			if (result != Result) {
				std::cout << "Quad expected value " << Result << " actual value " << result << '\n';
				Beep(750, 1000);
			}
			auto end = clock(); double elapsed = (double)end - start;
			std::cout << "time used= " << elapsed / CLOCKS_PER_SEC << " seconds\n";
		}
	}
	else
		std::cout << " cannot be factorised\n";
}

/* perform some simple tests. Returns true if x3 is prime */
static bool factortest(const Znum x3) {
	std::vector <zFactors> factorlist;
	Znum Quad[4], result;
	long long totalFactors = 0;

	auto start = clock();	// used to measure execution time
	double end, elapsed;

	std::cout << "\nTest: factorise ";
	ShowLargeNumber(x3, 6, true, false);
	std::cout << '\n';
	factorise(x3, factorlist, Quad);  // get factors

	result = 1;
	for (size_t i = 0; i < factorlist.size(); i++)
		for (int j = 1; j <= factorlist[i].exponent; j++)
			result *= factorlist[i].Factor;  // get product of all factors
	if (result != x3) {
		std::cout << "Factors expected value " << x3 << " actual value " << result << '\n';
		Beep(750, 1000);
	}
	/* check that sum of squares is correct */
	result = Quad[0] * Quad[0] + Quad[1] * Quad[1] + Quad[2] * Quad[2] + Quad[3] * Quad[3];
	if (result != x3) {
		std::cout << "Quad expected value " << x3 << " actual value " << result << '\n';
		Beep(750, 1000);
	}
	if (factorlist.size() > 1 || factorlist[0].exponent > 1) {
		/* x3 is not prime */
		for (auto f : factorlist) {
			totalFactors += f.exponent;
		}
		std::cout << "found " << factorlist.size() << " unique factors, total "
			<< totalFactors << " factors\nfound by";
		if (counters.tdiv > 0)
			std::cout << " trial division: " << counters.tdiv;
		if (counters.prho > 0)
			std::cout << " Pollard-rho: " << counters.prho;
		if (counters.pm1 > 0)
			std::cout << " power +/- 1: " << counters.pm1;
		if (counters.ecm > 0)
			std::cout << " elliptic curve: " << counters.ecm;
		if (counters.siqs > 0)
			std::cout << " SIQS: " << counters.siqs;
		if (counters.msieve > 0)
			std::cout << " Msieve: " << counters.msieve;
		if (counters.carm > 0)
			std::cout << " Carmichael: " << counters.carm;
		if (counters.leh > 0)
			std::cout << " Lehman: " << counters.leh;
		std::cout << '\n';
		end = clock();              // measure amount of time used
		elapsed = (double)end - start;
		std::cout << "time used= " << elapsed / CLOCKS_PER_SEC << " seconds\n";
		return false;   // not prime
	}
	else
		std::cout << "is prime \n";

	end = clock();              // measure amount of time used
	elapsed = (double)end - start;
	std::cout << "time used= " << elapsed / CLOCKS_PER_SEC << " seconds\n";
	return true;    // is prime
}

static void doTests(void) {
	Znum x3, x4, result;
	int i;

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
		"gcd (12,30)",                      6,
		"modinv (17,21)",                   5,
		"20 + 32^2/4*7",                 1812,   // DAs calculator returns 56 which is wrong
		"20 + 32^2/(4*7)",                 56,
		"20 - 32^2/4*7",                -1772,
		"(20-32)^2 / 4 * 7",              252,   //DAs calculator returns 5 which is wrong
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
		"5 < 6 == 7 < 8",                  -1,   // returns true (== has lower priority)
		"5 < 6 != 7 < 8",                   0,   // returns false (!= has lower priority)
		"5 < (6 != 7) < 8",                -1,   // returns true; expr evaluated from left to right
	};

	auto start = clock();	// used to measure execution time
	for (i = 0; i < sizeof(testvalues) / sizeof(testvalues[0]); i++) {
		removeBlanks(testvalues[i].text);  // it is necessary to remove spaces
		/* but it is not necessary to convert to upper case */
		stackIndex = 0;         // ensure stack is empty 
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
		x3++;
		factortest(x3);
	}

	/* tests below show a problem with pollard-rho for certain numbers */
	factortest(99999999973789); // = 6478429 * 15435841
	factortest(183038861417);   // =  408229 *   448373
	factortest(183475587821);   // =  409477 *   448073
	factortest(181919916457);   // =  400307 *   454451
	factortest(199996999457);   // =  441361 *   453137
	factortest(204493418837);   // =  401209 *   509693

	/* exercise code specifically for power +/-1 */
	mpz_ui_pow_ui(ZT(x3), 10, 20);  // x3 = 10^20
	x3 -= 1;                        // x3 = 10^20-1
	factortest(x3);
	x3 += 2;
	factortest(x3);
	ComputeExpr("n(10^10)^2-1", x3);  // test power of large number-1
	factortest(x3);

	ComputeExpr("120#-1", x3);
	factortest(x3);
	ComputeExpr("n(10^15)^2", x3);  // test power of large number
	factortest(x3);
	ComputeExpr("n(10^6+20)^1667", x3);  // test power of large prime number
	factortest(x3);
	ComputeExpr("n(10^7)^3*n(10^8)^2", x3);  // test powers of two large prime number
	factortest(x3);
	ComputeExpr("n(3*10^5+50)*n(3*10^5+500)", x3);  // test Lehman factorisation
	factortest(x3);

	/* test using carmichael numbers. note that 1st example has no small factors  */
	long long int carmichael[] = { 90256390764228001, 7156857700403137441,  1436697831295441,
		60977817398996785 };
	for (int i = 0; i < sizeof(carmichael) / sizeof(carmichael[0]); i++) {
		factortest(carmichael[i]);
	}

	/* set x3 to large prime. see https://en.wikipedia.org/wiki/Carmichael_number */
	ComputeExpr("2967449566868551055015417464290533273077199179985304335099507"
		"5531276838753171770199594238596428121188033664754218345562493168782883", x3);
	x4 = x3 * (313 * (x3 - 1) + 1) * (353 * (x3 - 1) + 1);
	/* in general numbers > about 110 digits cannot be factorised in a reasonable time 
	but this one can, because a special algorithm just for Carmichael numbers is used. */
	factortest(x4);   // 397-digit Carmichael number
	std::cout << "factorised 397-digit Carmichael number \n";

	ComputeExpr("n(10^24)*n(10^25)*n(10^26)*n(10^27)", x3);  
	factortest(x3);

	auto end = clock();   // measure amount of time used
	double elapsed = (double)end - start;
	std::cout << "Factorisation tests completed. Time used= " << elapsed / CLOCKS_PER_SEC << " seconds\n";
}

/* generate large random number, up to 128 bits */
static void largeRand(Znum &a) {
	a = ((long long)rand() << 32) + rand();
	a <<= 64;
	a += ((long long)rand() << 32) + rand();
}

static void doTests2(void) {
	srand(756128234);
	Znum x = (Znum)rand();
	auto start = clock();	// used to measure execution time

	
	for (int i = 1; i <= 44; i++) {
		if (i == 44) {
			std::cout << "Test 44 skipped \n";
			break;  /* note: with current seed value test 44 would take about 50 minutes! 
					(about 30 minutes using Msieve) */
		}

		if (i <= 39)
			mpz_mul_2exp(ZT(x), ZT(x), 8); // shift x left 8 bits
		x += rand();
		std::cout << "\nTest # " << i << " of 44 \n";
		ShowLargeNumber(x, 6, true, false);
		std::cout << '\n';
		doFactors(x, true); /* factorise x, calculate number of divisors etc */
	}
	auto end = clock();   // measure amount of time used
	double elapsed = (double)end - start;
	std::cout << "time used= " << elapsed / CLOCKS_PER_SEC << " seconds\n";
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
		//int exp;
		//a1 = BigToLL(aBI, exp);  // a1 * 2^exp = aBI
		//mpz_div_2exp(ZT(b1), ZT(a), exp);   // b1 a/(2^exp)
		//if (b1 != a1)
		//	gmp_printf("a = %Zx a1 = %Zx, exp = %d \n", a, a1, exp);
		//assert(b1 == a1);

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
	p = 10;
	double errorv;
	Znum error, relerror;
	for (int i = 1; i < 55;  i++, p *= 1000000 ) {
		if (!ZtoBig(pBI, p))
			break;            // p too large to convert so stop
		double pdb = pBI.log();
		if (pdb > 708)
			break;
		expBigInt(amBI, pdb);
		BigtoZ(am, amBI);
		error = p - am;
		errorv = (double)(error) / (double)p ;
		relerror = (10000000000000000LL * error)/p;
		std::cout << " pdb = " << pdb/std::log(10) 
			<< " error = " << errorv << " relerror = " << relerror <<'\n';
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
this test will take about 2.5 hours, but ony about 15 mins using improved msieve! 
*/
static void doTests4(void) {
	int px = 0;
	Znum m;

	auto start = clock();	// used to measure execution time
	if (primeListMax < 1000)
		generatePrimes(393203);
	for (px = 0; px <= 69; px++) {
		std::cout << "\ntest " << px + 1 << " of 70 ";
		mpz_ui_pow_ui(ZT(m), 2, primeList[px]);
		m--;
		if (factortest(m)) /* factorise m, calculate number of divisors etc */
			std::cout << "2^" << primeList[px] << " - 1 is prime \n";
	}
	auto end = clock();              // measure amount of time used
	auto elapsed = (double)end - start;
	std::cout << "tests completed  time used= " << elapsed / CLOCKS_PER_SEC << " seconds\n";
}


int main(int argc, char *argv[]) {
	std::string expr, expupper;
	Znum Result;
	retCode rv;

	const static char helpmsg[] =
		"You can enter expressions that use the following operators, functions and parentheses:\n"
		"^ or ** : exponentiation (the exponent must be greater than or equal to zero).\n"
		"*       : multiplication\n"
		"/       : integer division\n"
		"%%       : modulus (remainder from integer division)\n"
		"+       : addition\n"
		"-       : subtraction\n"
		"SHL or <<: Shift left the number of bits specified on the right operand.\n"
		"SHR or >>: Shift right the number of bits specified on the right operand.\n"
		"<, >, <=, >= for comparisons.  The operators return zero for false and -1 for true.\n"
		"==, != for equality. The operators return zero for false and -1 for true.\n"
		"NOT, AND, XOR, OR  for bitwise operations.\n"
		"n!      : factorial (n must be greater than or equal to zero).\n"
		"n!!     : double factorial (n must be greater than or equal to zero).\n"
		"p#      : primorial (product of all primes less or equal than p).\n"
		"C       : binomial coefficient. nCk = n!/(k!*(n-k)!) but is more efficient\n"
		"B(n)    : Previous probable prime before n\n"
		"F(n)    : Fibonacci number Fn\n"
		"L(n)    : Lucas number Ln = F(n-1) + F(n+1)\n"
		"N(n)    : Next probable prime after n\n"
		"P(n)    : Unrestricted Partition Number (number of decompositions of n into sums of integers without regard to order).\n"
		"Pi(n)   : Number of primes <= n \n"
		"Gcd(m, n)       : Greatest common divisor of m and n.\n"
		"Modinv(m, n)    : inverse of m modulo n, only valid when gcd(m, n) = 1.\n"
		"Modpow(m, n, r) : finds m^n modulo r Equivalent to m^n%%r but more efficient.\n"
		"Totient(n)      : finds the number of positive integers less than n which are relatively prime to n.\n"
		"IsPrime(n)      : returns zero if n is not probable prime, -1 if it is.\n"
		"NumDivs(n)      : Number of positive divisors of n either prime or composite.\n"
		"SumDivs(n)      : Sum of positive divisors of n either prime or composite.\n"
		"NumDigits(n, r) : Number of digits of n in base r.\n"
		"SumDigits(n, r) : Sum of digits of n in base r.\n"
		"RevDigits(n, r) : finds the value obtained by writing backwards the digits of n in base r.\n"
		"ConcatFact(m,n) : Concatenates the prime factors of n according to the mode m\n"
		"R2(n)   : Number of ways n can be expressed as the sum of x^2+y^2. (order and sign of x and y are significant \n"
		"LE(a,p) : Legendre value for (a/p) \n"
		"JA(a,p) : Jacobi value for (a/p) \n"
		"KR(a,p) : Kronecker value for (a/p) \n"
		"Also the following commands: X=hexadecimal o/p, D=decimal o/p \n"
		"F = do factorisation, N = Don't factorise, S = Spanish, E=English\n"
		"HELP (this message) and EXIT\n";

	const static char ayuda[] =
		"Puedes ingresar expresiones que usen los siguientes operadores y paréntesis:\n"
		"+ para suma\n"
		"- para resta\n"
		"* para multiplicación\n"
		"/ para división entera\n"
		"%% para el resto de la división entera\n"
		"^ o ** para exponenciación(el exponente debe ser mayor o igual que cero).\n"
		"<, == , >; <= , >= , != para comparaciones.Los operadores devuelven cero si es falso y - 1 si es verdadero.\n"
		"AND, OR, XOR, NOT para lógica binaria.\n"
		"SHL  : Desplazar a la izquierda la cantidad de bits indicada en el operando derecho.\n"
		"SHR  : Desplazar a la derecha la cantidad de bits indicada en el operando derecho.\n"
		"n!   : factorial(n debe ser mayor o igual que cero).\n"
		"p#   : primorial(producto de todos los primos menores o iguales a p).\n"
		"B(n) : Número probablemente primo anterior a n\n"
		"F(n) : Número de Fibonacci Fn\n"
		"L(n) : Número de Lucas Ln = F(n-1) + F(n+1)\n"
		"N(n) : Número probablemente primo posterior a n\n"
		"P(n) : particiones irrestrictas(cantidad de descomposiciones de n en sumas de números enteros sin tener en cuenta el orden).\n"
		"Gcd(m, n)       : Máximo común divisor de estos dos números enteros.\n"
		"Modinv(m, n)    : inverso de m modulo n, sólo válido cuando gcd(m, n) = 1.\n"
		"Modpow(m, n, r) : halla m^n módulo r.\n"
		"Totient(n)      : cantidad de enteros positivos menores que n coprimos con n.\n"
		"IsPrime(n)      : returna cero si n no es un primo probable y - 1 si lo es.\n"
		"NumDivs(n)      : cantidad de divisores positivos de n primos o compuestos.\n"
		"SumDivs(n)      : suma de divisores positivos de n primos o compuestos.\n"
		"NumDigits(n, r) : cantidad de dígitos de n en base r.\n"
		"SumDigits(n, r) : suma de dígitos de n en base r.\n"
		"RevDigits(n, r) : halla el valor que se obtiene escribiendo para atrás los dígitos de n en base r.\n";

	try {
		hConsole = GetStdHandle(STD_OUTPUT_HANDLE);  // gete handle for console window
		if (hConsole == INVALID_HANDLE_VALUE)
		{
			fprintf_s(stderr, "GetStdHandle failed with %d at line %d\n", GetLastError(), __LINE__);
			Beep(750, 1000);
			return EXIT_FAILURE;
		}

		char banner[] = "Compiled on "  __DATE__ " at " __TIME__ "\n";
		printf("Version 2.0 %s", banner);
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

		while (true) {
			if (lang == 0) {
				printf("enter expression to be processed, or HELP, or EXIT\n");
			}
			else
				printf("ingrese la expresión para ser procesada, o AYUDA o SALIDA\n");

			getline(std::cin, expr);  // expression may include spaces

			strToUpper(expr, expupper);		// convert to UPPER CASE 
			if (expupper == "EXIT" || expupper == "SALIDA")
				break;

			if (expupper == "HELP" || expupper == "AYUDA") {
				if (lang == 0)
					printf(helpmsg);   // english
				else
					printf(ayuda);     // spanish
				continue;
			}
			if (expupper == "E") { lang = 0; continue; }           // english
			if (expupper == "S") { lang = 1; continue; }	       // spanish (Español)
			if (expupper == "F") { factorFlag = true; continue; }  // do factorisation
			if (expupper == "N") { factorFlag = false; continue; } // don't do factorisation
			if (expupper == "X") { hex = true; continue; }         // hexadecimal output
			if (expupper == "D") { hex = false; continue; }        // decimal output
			if (expupper == "TEST") {
				doTests();         // do basic tests 
				continue; 
			}
			if (expupper == "TEST2") {
				doTests2();         // do basic tests 
				continue;
			}
#ifdef BIGNBR
			if (expupper == "TEST3") {
				doTests3();         // do basic tests 
				continue;
			}
#endif
			if (expupper == "TEST4") {
				doTests4();         // do R3 tests 
				continue;
			}
			if (expupper.substr(0, 6) == "MSIEVE") {
				msieveParam(expupper);
				continue;
			}

			auto start = clock();	// used to measure execution time
			removeBlanks(expr);     // remove any spaces 
			stackIndex = 0;         // ensure stack is empty 
			rv = ComputeExpr(expr, Result); /* analyse expression, compute value*/

			if (rv != retCode::EXPR_OK) {
				textError(rv);   // invalid expression; print error message
				if (exprIndex < expr.size())
				{
					std::cout << "Unprocessed text:  " << expr.substr(exprIndex) << '\n';
				}
			}
			else {
				if (exprIndex < expr.size())
				{
					std::cout << "Unprocessed text:  " << expr.substr(exprIndex) << '\n';
				}
				std::cout << " = ";
				ShowLargeNumber(Result, 6, true, hex);   // print value of expression
				std::cout << '\n';				
				if (factorFlag) {
					doFactors(Result,false); /* factorise Result, calculate number of divisors etc */
				}
			}

			auto end = clock();   // measure amount of time used
			double elapsed = (double)end - start;
			std::cout << "time used= " << elapsed / CLOCKS_PER_SEC << " seconds\n";
		}

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
