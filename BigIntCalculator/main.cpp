#include <cstdio>
#include <iostream>
#include <string>
#include <vector>
#include <ctype.h>
#include <limits.h>
#include <assert.h>
#include <Windows.h>
#include "bignbr.h"
#include "factor.h"


typedef boost::multiprecision::mpz_int Znum;
/* access underlying mpz_t inside a Znum */
#define ZT(a) a.backend().data()

#define PAREN_STACK_SIZE            100
int lang = 0;             // 0 English, 1 = Spanish
static int stackIndex, exprIndex;
bool hex = false;		// set true if output is in hex
bool factorFlag = true;

/* list of operators, arranged in order of priority */
typedef enum {
	oper_comb		 =  0,   // nCk, also known as binomial coefficient
	oper_power       =  1,
	oper_multiply    =  2,
	oper_divide      =  3,
	oper_remainder   =  4,
	oper_unary_minus =  5,
	oper_plus        =  6,
	oper_minus       =  7,
	oper_shr         =  8,
	oper_shl         =  9,
	oper_not_greater = 10,
	oper_not_less    = 11,
	oper_not_equal   = 12,
	oper_equal       = 13,
	oper_greater     = 14,
	oper_less        = 15,
	oper_not         = 16,
	oper_and         = 17,
	oper_or          = 18,
	oper_xor         = 19,
	oper_leftb       = '(',
}oper_code;

enum eExprErr
{
	EXPR_NUMBER_TOO_LOW = -100,
	EXPR_NUMBER_TOO_HIGH,
	EXPR_INTERM_TOO_HIGH,
	EXPR_DIVIDE_BY_ZERO,
	EXPR_PAREN_MISMATCH,
	EXPR_SYNTAX_ERROR,
	EXPR_TOO_MANY_PAREN,
	EXPR_INVALID_PARAM,
	EXPR_ARGUMENTS_NOT_RELATIVELY_PRIME,
	EXPR_BREAK,
	EXPR_OUT_OF_MEMORY,
	EXPR_CANNOT_USE_X_IN_EXPONENT,
	EXPR_DEGREE_TOO_HIGH,
	EXPR_EXPONENT_TOO_LARGE,
	EXPR_EXPONENT_NEGATIVE,
	EXPR_LEADING_COFF_MULTIPLE_OF_PRIME,
	EXPR_CANNOT_LIFT,
	EXPR_MODULUS_MUST_BE_GREATER_THAN_ONE,
	EXPR_MODULUS_MUST_BE_PRIME_EXP,
	EXPR_BASE_MUST_BE_POSITIVE,
	EXPR_POWER_MUST_BE_POSITIVE,
	EXPR_MODULUS_MUST_BE_NONNEGATIVE,
	EXPR_VAR_OR_COUNTER_REQUIRED,
	EXPR_OK = 0,
	EXPR_FAIL = 1
};

bool *primeFlags = NULL;
unsigned long long *primeList = NULL;
unsigned int prime_list_count = 0;
unsigned long long int primeListMax = 0;
const unsigned long long max_prime = 1000000000;  // arbitrary limit 10^9,

/* function declarations, only for functions that have forward references */
static eExprErr ComputeExpr(const std::string &expr, Znum &ExpressionResult);
long long MulPrToLong(const Znum &x);
void biperm(int n, mpz_t &result);   // declaration for external function
static bool getBit(const unsigned long long int x, bool array[]);
void generatePrimes(unsigned long long int max_val);

char *hexbuffer = NULL;
/* Convert number to hexdecimal. Ensure that if number is negative the leftmost
bit of the most significant digit is set, and conversely, if the number is positive
the leftmost bit of the most significant digit is not set. This is done by 
prefixing the output with '0' or 'f' when necessary. */
char* getHex(Znum Bi_Nbr) {
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
		while (Bi_Nbr != -1) {
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
	return hexbuffer;
}

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
		rv = mpz_get_si(ZT(x)); // convert to normal integer
		return rv;
	}
	else
		throw std::range_error("big number cannot be converted to 64-bit integer");
	return 0;
}

/* calculate number of divisors of n, given its list of prime factors.
NB for n=1 function returns 2, but correct value would be 1. 
From just the list of exponents we can't distinguish n=1 from n=prime */
static Znum NoOfDivisorsF(const std::vector <int> exponents) {
	Znum result = 1;

	for (auto i: exponents) {
		result *= i + 1;
	}
	return result;
}

// sum of divisors is the product of (p^(e+1)-1)/(p-1) where p=prime factor and e=exponent.
static Znum SumOfDivisors(const std::vector <int> exponents, const std::vector <Znum> primes) {
	Znum result = 1, term;

	if (primes[0] == 1)
		return 1;		// special case: if original number is 1
	for (size_t i = 0; i < primes.size(); i++) {
		mpz_pow_ui(ZT(term), ZT(primes[i]), exponents[i] + 1);  // p^(e+1)
		term = (term-1)/(primes[i]-1);	                        // (p^(e+1)-1)/(p-1)
		result *= term;
	}
	return result;
}

// Find Euler's Totient as the product of p^(e-1)*(p-1) where p=prime and e=exponent.
static Znum Totient(const std::vector <int> exponents, const std::vector <Znum> primes) {
	Znum result = 1, term;

	for (size_t i = 0; i < primes.size(); i++) {
		mpz_pow_ui(ZT(term), ZT(primes[i]), exponents[i] - 1);  // p^(e-1)
		term = term* (primes[i] - 1);	                        // (p^(e-1)-1)*(p-1)
		result *= term;
	}
	return result;
}

/* For any positive integer n, define μ(n) as the sum of the primitive nth roots of unity. 
   It has values in {−1, 0, 1} depending on the factorization of n into prime factors:
   μ(n) = 1 if n is a square-free positive integer with an even number of prime factors.
   μ(n) = −1 if n is a square-free positive integer with an odd number of prime factors.
   μ(n) = 0 if n has a squared prime factor.
*/
static int mobius(const std::vector <int> exponents) {
	
	for (auto ex : exponents) {
		if (ex > 1)
			return 0;		// n is not square-free
	}
	auto result = exponents.size();
	if ((result & 1) == 1)
		return -1;  // odd nuber of  prime factors
	else
		return 1;    // even number of prime factors 
}

/* calculate Euler's totient for n as the product of p^(e-1)*(p-1) 
where p=prime and e=exponent.*/
static Znum ComputeTotient(const Znum n) {
	std::vector <Znum> factorlist;
	std::vector <int> exponentlist;
	Znum quads[4];  // needed, but not used
	if (n == 1)
		return 1;
	auto rv = factorise(n, factorlist, exponentlist, quads);
	if (rv && !factorlist.empty()) {
		auto divisors = Totient(exponentlist, factorlist);
		return divisors;
	}
	else return 0;
}

/* calculate number of divisors of n, from its list of prime factors. */
static Znum ComputeNumDivs(const Znum n) {
	std::vector <Znum> factorlist;
	std::vector <int> exponentlist;
	Znum quads[4];
	if (n == 1)
		return 1;  // 1 only has one divisor. NoOfDivisors can't handle that case
	auto rv = factorise(n, factorlist, exponentlist, quads);
	if (rv && !factorlist.empty()) {
		auto divisors = NoOfDivisorsF(exponentlist);
		return divisors;
	}
	else return 0;
}

// sum of divisors is the product of (p^(e+1)-1)/(p-1) where p=prime factor and e=exponent.
static Znum ComputeSumDivs(const Znum n) {
	std::vector <Znum> factorlist;
	std::vector <int> exponentlist;
	Znum quads[4];  // needed but not used
	if (n == 1)
		return 1;   // 1 only has 1 divisor. 
	auto rv = factorise(n, factorlist, exponentlist, quads);
	if (rv && !factorlist.empty()) {
		auto divisors = SumOfDivisors(exponentlist, factorlist);
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

static long long NoOfBits(const Znum &n) {
	auto result = mpz_sizeinbase(ZT(n), 2);  // calculate number of bits
	return result;
}

/* IsPrime(n): returns zero if n is not probable prime, -1 if it is. */
static int PrimalityTest(const Znum &Value) {
	static gmp_randstate_t state;
	static bool first = true;
	int rv;
	assert(Value >= 0);

	if (first) {
		gmp_randinit_default(state);
		first = false;
	}

	rv = mpz_probable_prime_p(ZT(Value), state, 16, 0);  
	/* rv is 1 if n is probably prime, or 0 if n is definitely composite.*/
	if (rv == 1)
		return -1;  // number is probably prime (less than 1 in 2^16 chance it's not prime)
	else 
		return 0;  // number is definately composite
}
/* same purpose as PrimalityTest but optimised for smaller numbers. 1st call can be very slow,
but subsequent calls are very quick. */
static int PrimalityTestSmall(const long long Value) {
	assert(Value >= 0);
	/* deal with even values first, because getBit below only handles odd values*/
	if (Value == 2) return -1;		//  2 is only even prime
	if ((Value & 1) == 0) return 0;   // even numbers > 2 are not prime

	if (Value <= max_prime) {
		if (primeFlags == NULL) {
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
static eExprErr ComputeBack(const Znum &n, Znum &p) {
	if (n < 3)
		return EXPR_INVALID_PARAM;  // 2 is the smallest prime
	if (n == 3) {
		p = 2;
		return EXPR_OK;
	}
	Znum i;		// largest odd number below n
	i = ((n & 1) == 1) ? n - 2 : n - 1;  // i >= 3

	for (; i > 0; i-=2) {
		if (PrimalityTest(i) == -1) {
			p = i;
			return EXPR_OK;
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
static Znum  concatFact(Znum mode, Znum num) {
	std::vector <Znum> factorlist;
	std::vector <int> exponentlist;
	std::string result;
	Znum rvalue;
	Znum quads[4];  // needed, but not used
	bool descending = ((mode & 1) == 1);
	bool repeat = ((mode & 2) == 2);
	char *buffer = NULL;
	/* get factors of num */
	auto rv = factorise(num, factorlist, exponentlist, quads);
	if (rv && !factorlist.empty()) {
		if (descending)   /* start with largest factor */
			for (ptrdiff_t i = factorlist.size()-1; i >=0; i--) {
				buffer = mpz_get_str(NULL, 10, ZT(factorlist[i]));
				if (!repeat)
					result += buffer; // concatenate factor
				else
					for (int j =1; j<= exponentlist[i]; j++)
						result += buffer; // concatenate factor
				free(buffer);
				buffer = NULL;
			}
		else  /* start with smallest factor */
			for (size_t i = 0; i < factorlist.size(); i++) {
				buffer = mpz_get_str(NULL, 10, ZT(factorlist[i]));
				if (!repeat)
					result += buffer;  // concatenate factor
				else
					for (int j = 1; j <= exponentlist[i]; j++)
						result += buffer;  // concatenate factor
				free(buffer);
				buffer = NULL;
			}
		mpz_set_str(ZT(rvalue), result.data(), 10); /* convert back from a string to a number */
		return rvalue;
	}
	return 0;
}

typedef enum{
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
fn_invalid = -1,
} fn_Code;

struct  functions {
	char fname[11];    // maximum name length 10 chars (allow for null terminator)
	int  NoOfParams;   // number of parameters 
	fn_Code  fCode;        // integer code for function
};

/* list of function names. No function name can begin with C because this would 
 conflict with the C operator. Longer names should come before short ones 
 to avoid mismatches */
const static std::array <struct functions, 17> functionList{
	"GCD",       2,  fn_gcd,    // name, number of parameters, code
	"MODPOW",    3,  fn_modpow,
	"MODINV",    2,  fn_modinv,
	"TOTIENT",   1,  fn_totient,
	"NUMDIVS",   1,  fn_numdivs,
	"SUMDIVS",   1,  fn_sumdivs,
	"SUMDIGITS", 2,  fn_sumdigits,
	"NUMDIGITS", 2,  fn_numdigits,
	"REVDIGITS", 2,  fn_revdigits,
	"ISPRIME",   1,	 fn_isprime,
	"FactConcat",2,  fn_concatfact,     // FactConcat must come before F
	"F",         1,  fn_fib,			// fibonacci
	"L",         1,  fn_luc,			// Lucas Number
	"PI",		 1,  fn_primePi,		// prime-counting function. PI must come before P
	"P",         1,  fn_part,			// number of partitions
	"N",         1,  fn_np,				// next prime
	"B",         1,  fn_pp,				// previous prime

};

/* Check whether expr matches any function name. If so ,compute value(s) of 
 parameter(s) leave them on the stack, set Gfcode and return EXPR_OK. 
 If not, return value 1 (EXPR_FAIL).
Updates global var exprIndex. */
static eExprErr func(const std::string &expr, const bool leftNumberFlag,
 	 int exprIndexLocal, fn_Code &Gfcode) {

	fn_Code fnCode = fn_invalid;
	long long ix;
	int NoOfArgs = 0;
	const char *ptrExpr = expr.data() + exprIndexLocal;
	/* try to match function name */
	for (ix = 0; ix < (ptrdiff_t)functionList.size(); ix++) {
		if (_strnicmp(ptrExpr, functionList[ix].fname, strlen(functionList[ix].fname)) == 0) {
			fnCode = functionList[ix].fCode;
			NoOfArgs = functionList[ix].NoOfParams;
			break;  // found a match for the function name 
		}
	}

	if (fnCode == fn_invalid)
		return EXPR_FAIL;  // no function name found
		
	exprIndexLocal += (int)strlen(functionList[ix].fname); // move exprIndex past function name
	if (leftNumberFlag)
	{
		return EXPR_SYNTAX_ERROR;
	}

	if (exprIndexLocal >= expr.length() || expr[exprIndexLocal++] != '(')
	{
		return EXPR_SYNTAX_ERROR;  // no opening bracket after function name
	}
	int Bracketdepth = 0;
	for (ix = exprIndexLocal-1; ix < (long long)expr.length(); ix++) {
		if (expr[ix] == '(') Bracketdepth++;
		if (expr[ix] == ')') Bracketdepth--;
		if (Bracketdepth == 0) break;  // found matching closing bracket
	}
	if (Bracketdepth != 0)
		return EXPR_PAREN_MISMATCH;
	long long exprLength = ix - exprIndexLocal + 1;  // length of text up to closing bracket.

	for (int index = 0; index < NoOfArgs; index++)
	{
		eExprErr retcode;
		char compareChar;
		Znum ExpressionResult;     // needed, but value is not actually used.

		if (stackIndex >= PAREN_STACK_SIZE)
		{
			return EXPR_TOO_MANY_PAREN;
		}
		/* substring in function call below includes all characters from the parameters now
		being processed up to and including the closing bracket. */
		retcode = ComputeExpr(expr.substr(exprIndexLocal, exprLength), ExpressionResult);
		if (retcode !=  EXPR_OK) { 
			return retcode; 
		}
		exprIndexLocal += exprIndex;  // exprIndexLocal now points past expression just evaluated.
		exprLength -= exprIndex;
		/* if there's more than 1 parameter they are separated with commas. 
		The last one is followed by ')'*/
		compareChar = (index == NoOfArgs - 1 ? ')' : ',');
		if (exprIndexLocal >= expr.length() || expr[exprIndexLocal++] != compareChar)
		{
			return EXPR_SYNTAX_ERROR;  // no comma or ) when needed
		}
		exprLength--;		// subtract length of comma or closing bracket
		stackIndex++;
	}

	stackIndex -= NoOfArgs;
	exprIndex = exprIndexLocal;  // update global index
	Gfcode = fnCode;  // pass function code back to caller
	return EXPR_OK;
}

/* SHL: Shift left the number of bits specified on the right operand. If the
number of bits to shift is negative, this is actually a right shift */
static enum eExprErr ShiftLeft(const Znum first, const Znum bits, Znum &result) {
	if (bits > LLONG_MAX || bits < LLONG_MIN)
		return EXPR_INVALID_PARAM;

	long long shift = MulPrToLong(bits);

	if (shift == 0) {
		result = first;
		return EXPR_OK;  // shift zero bits; nothing to do
	}

	// there is no built-in shift operator for Znums, so it is simulated
	// using multiplication or division
	if (shift > 0) {  
		mpz_mul_2exp(ZT(result), ZT(first), shift);
		return EXPR_OK;
	}
	else {   // note use of floor division
		mpz_fdiv_q_2exp(ZT(result), ZT(first), -shift);
		return EXPR_OK;
	}
}

static Znum stackValues[PAREN_STACK_SIZE];
static oper_code stackOperators[PAREN_STACK_SIZE];

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
void setBit(const unsigned long long int x, bool array[]) {
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
	/*for (int i = 0; i <= max_val / 2; i++)
		primeFlags[i] = false;*/
	//memset(primeFlags, false, max_val / 2 + 1)

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

/* list of operator priority values. lower value = higher priority */
const static int operPrio[] =
{
	0,				  // combination
	1,                // Power
	2, 2, 2,          // Multiply, divide and remainder.
	3,                // Unary minus.
	4, 4,             // Plus and minus.
	5, 5,             // Shift right and left.
	6, 6, 6, 6, 6, 6, // Six comparison operators (equal, greater, less, etc.)
	7,                // NOT.
	8, 8, 8,          // AND, OR, XOR.
};

/* process one operator with 1 or 2 operands on stack. 
NOT and unary minus have 1 operand. All the others have two.
put result back on stack.
some operators can genererate an error condition e.g. EXPR_DIVIDE_BY_ZERO
otherwise return EXPR_OK. */
static enum eExprErr ComputeSubExpr(void)
{
	auto stackOper = stackOperators[--stackIndex];
	Znum firstArg = stackValues[stackIndex];
	Znum secondArg = stackValues[stackIndex + 1];
	Znum &result = stackValues[stackIndex];

	switch (stackOper)
	{
	case oper_comb: {  // calculate nCk AKA binomial coefficient
		if (secondArg > INT_MAX)
		return EXPR_INTERM_TOO_HIGH;
		if (secondArg < INT_MIN)
			return EXPR_INVALID_PARAM;
		long long k = MulPrToLong(secondArg);
		mpz_bin_ui(ZT(result), ZT(firstArg), k);
		return EXPR_OK;
	}
	case oper_plus: {
		result = firstArg + secondArg; //BigIntAdd(firstArg, secondArg, result);
		return EXPR_OK;
	}
	case oper_minus: {
		result = firstArg - secondArg; //BigIntSubt(firstArg, secondArg, result);
		return EXPR_OK;
	}
	case oper_unary_minus: {
		result = -secondArg; //BigIntNegate(secondArg, result);
		return EXPR_OK;
	}
	case oper_divide: {
		if (secondArg == 0)
			return EXPR_DIVIDE_BY_ZERO;  // result would be infinity
		result = firstArg / secondArg; //return BigIntDivide(firstArg, secondArg, result);
		return EXPR_OK;
	}
	case oper_multiply: {
		auto resultsize = NoOfBits(firstArg) + NoOfBits(secondArg);
		if (resultsize > 66439)  // more than 66439 bits -> more than 20,000 decimal digits
			return EXPR_INTERM_TOO_HIGH;
		result = firstArg * secondArg;
		return EXPR_OK;
	}
	case oper_remainder: {
		if (secondArg == 0)
			return EXPR_DIVIDE_BY_ZERO;  // result would be infinity
		result = firstArg % secondArg;   //return BigIntRemainder(firstArg, secondArg, result);
		return EXPR_OK;
	}
	case oper_power: {
		if (secondArg > INT_MAX)
			return EXPR_EXPONENT_TOO_LARGE;
		if (secondArg < 0)
			return EXPR_EXPONENT_NEGATIVE;
		long long exp = MulPrToLong(secondArg);
		auto resultsize = (NoOfBits(firstArg)-1)* exp;  // estimate number of bits for result
		if (resultsize > 66439)  // more than 66439 bits -> more than 20,000 decimal digits
			return EXPR_INTERM_TOO_HIGH;
		mpz_pow_ui(ZT(result), ZT(firstArg), exp);
		return EXPR_OK;
	}
	case oper_equal: {
		if (firstArg == secondArg)
			result = -1;
		else
			result = 0;
		return EXPR_OK;
	}
	case oper_not_equal: {
		if (firstArg != secondArg)
			result = -1;
		else
			result = 0;
		return EXPR_OK;
	}
	case oper_greater: {
		if (firstArg > secondArg)
			result = -1;
		else
			result = 0;
		return EXPR_OK;
	}
	case oper_not_greater: {
		if (firstArg <= secondArg)
			result = -1;
		else
			result = 0;
		return EXPR_OK;
	}
	case oper_less: {
		if (firstArg < secondArg)
			result = -1;
		else
			result = 0;
		return EXPR_OK;
	}
	case oper_not_less: {
		if (firstArg <= secondArg)
			result = -1;
		else
			result = 0;
		return EXPR_OK;
	}
	case oper_shl: {
		ShiftLeft(firstArg, secondArg, result);
		return EXPR_OK;
	}
	case oper_shr: {
		secondArg = -secondArg;   // invert sign of shift
		ShiftLeft(firstArg, secondArg, result);
		return EXPR_OK;
	}
	case oper_not: {   // Perform binary NOT as result <- -1 - argument.
		result = -1 - secondArg;  // assumes 2s complement binary numbers
		return EXPR_OK;
	}
	case oper_and: {  // Perform binary AND.
		mpz_and(ZT(result), ZT(firstArg), ZT(secondArg));
		return EXPR_OK;
	}
	case oper_or: {   // Perform binary OR.
		mpz_ior(ZT(result), ZT(firstArg), ZT(secondArg));
		return EXPR_OK;
	}
	case oper_xor: {   // Perform binary XOR.
		mpz_xor(ZT(result), ZT(firstArg), ZT(secondArg));
		return EXPR_OK;
	}
	default:
		abort();	// should never get here
	}
}

struct oper_list{
	char oper[4];
	oper_code operCode;
	int operPrio;
} ;
/* list of operators that have format <expression> <operator> <expression> 
with corresponding codes and priority. However, operator NOT (bitwise negation) 
is in this list although it only has one operand. 
Operators ! (factorial) !! (double factorial) and # (primorial) are not in this 
list because the convention is that they follow the number or expression they 
operate on. */
const static struct oper_list operators[]  {
	{"C",	oper_comb,		  0},
	{ "^",  oper_power,       1},
	{ "**", oper_power,       1},     // can use ^ or ** for exponent
	{ "*", oper_multiply,     2},
	{ "/", oper_divide,       2},
	{ "%", oper_remainder,    2},
	{ "+", oper_plus,         4},
	{ "-", oper_minus,        4},
	{ "SHL", oper_shl,        5},
	{ "<<",  oper_shl,        5},     // can use << or SHL for left shift
	{ "SHR", oper_shr,        5},
	{ ">>", oper_shr,         5},     // can use SHR or >> for right shift
	{ "<=", oper_not_greater, 6},
	{ ">=", oper_not_less,    6},
	{ "!=", oper_not_equal,   6},
	{ "==", oper_equal,       6},
	{ ">", oper_greater,      6},
	{ "<", oper_less,         6},
	{ "NOT", oper_not,        7},      // bitwise NOT
	{ "AND", oper_and,        8},      // bitwise AND
	{ "OR", oper_or,          8},      // bitwise OR
	{ "XOR", oper_xor,        8}   };  // bitwise exclusive or

/* search for operators e.g. '+', '*'  that have format <expression> <operator> <expression>
or <operator> <expression>. return index of operator or -1 */
static void operSearch(const std::string &expr, long long &opcode) {
	opcode = 0;
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

Hexadecimal numbers are preceded by 0X or 0X
*/

/* analyse an expression supplied as a string. This function uses indirect recursion;
it calls func to analyse the parameters of a function, which in turn calls
computeExpr again for each parameter to the function. 
The numerical value of the expression is returned in ExpressionResult and on the stack */
static eExprErr ComputeExpr(const std::string &expr, Znum &ExpressionResult)
{
	static int depth=0;  // measures depth of recursion, increment on entry, decrement on return 

	bool leftNumberFlag = false;
	enum eExprErr SubExprResult;
	int startStackIndex = stackIndex;    // for nested calls this is > zero

	if (expr.empty()) {
		return EXPR_INVALID_PARAM;
	}

	exprIndex = 0;
	depth++;		// measure recursion depth
	/* exit from loop below when end of expr or a , or ) is reached 
	if an error occurs, return error code immediately to caller */
	while (exprIndex < expr.length())
	{
		int charValue;
		fn_Code fcode;
		oper_code opCode;
		long long opIndex;
		int exprLength;
		int exprIndexAux;
		Znum factorial;
		eExprErr retcode;
		int rv;

		exprLength = (int)expr.length();
		charValue = toupper(expr[exprIndex]);  // get next character of the expression

		/* find operator or evaluate number or evaluate function */

		/* first look for operators. */
		operSearch(expr.substr(exprIndex), opIndex);
		if (opIndex != -1) {
			/* found operator e.g. + that has format <expression> <operator> <expression>*/
			opCode= operators[opIndex].operCode;
			exprIndex += (int)strlen(operators[opIndex].oper);

			if ((opCode == oper_plus || opCode == oper_minus) && leftNumberFlag == false)
			{                    // Unary plus/minus operator
				if (opCode == oper_plus)
				{
					continue;  // process more of expr
				}
				else
				{
					if (stackIndex > startStackIndex && 
						stackOperators[stackIndex - 1] == oper_unary_minus)
					{
						stackIndex--;
						continue;   // process more of expr
					}
					if (stackIndex >= PAREN_STACK_SIZE)
					{
						depth--;				// adjust call depth
						return EXPR_TOO_MANY_PAREN;
					}
					stackOperators[stackIndex++] = oper_unary_minus; /* Unary minus */
					continue;    // process more of expr
				}
			}
			if ((leftNumberFlag == false) != (opCode == oper_not))
			{     // Missing left operator if operator is not NOT or
					// extra left operator if operator is NOT.
				depth--;				// adjust call depth
				return EXPR_SYNTAX_ERROR;
			}
			if (opCode != oper_power)
			{  // Power operator has right associativity.
				while (stackIndex > startStackIndex) {
					/* check whether we can process any stacked values and operators*/
					if (stackOperators[stackIndex - 1] == oper_leftb)
						break;
					auto p1 = operPrio[(int)stackOperators[stackIndex - 1]];
					auto p2 = operPrio[opCode];
					if (p1 > p2)    
						break;    // stacked operator has lower priority
					if ((SubExprResult = ComputeSubExpr()) != EXPR_OK)
					{
						depth--;				// adjust call depth
						return SubExprResult;   // failed! return error code
					}
				}
			}
			stackOperators[stackIndex++] = opCode;
			leftNumberFlag = false;
			continue;     // process more of expr
		}
		else if (charValue == '!' && expr[exprIndex + 1] == '!') { // Calculating  double factorial.
			if (leftNumberFlag == false)
			{
				depth--;				// adjust call depth
				return EXPR_SYNTAX_ERROR;
			}
			if (stackValues[stackIndex] > 11081)
			{
				depth--;				// adjust call depth
				return EXPR_INTERM_TOO_HIGH;
			}
			if (stackValues[stackIndex] < 0)
			{
				depth--;				// adjust call depth
				return EXPR_INVALID_PARAM;
			}

			long long temp = llabs(MulPrToLong(stackValues[stackIndex]));
			mpz_2fac_ui(ZT(factorial), temp);  // get double factorial
			stackValues[stackIndex] = factorial;
			exprIndex+= 2;
			continue;      // process more of expr
		}
		else if (charValue == '!') {           // Calculating factorial.
			if (leftNumberFlag == false)
			{
				depth--;				// adjust call depth
				return EXPR_SYNTAX_ERROR;
			}
			if (stackValues[stackIndex] > 5984 )
			{
				depth--;				// adjust call depth
				return EXPR_INTERM_TOO_HIGH;
			}
			if (stackValues[stackIndex] < 0)
			{
				depth--;				// adjust call depth
				return EXPR_INVALID_PARAM;
			}
			long long temp = llabs(MulPrToLong(stackValues[stackIndex]));
			mpz_fac_ui(ZT(factorial), temp);  // get factorial
			stackValues[stackIndex] = factorial;
			exprIndex++;
			continue;      // process more of expr
		}
		else if (charValue == '#') {           // Calculating primorial.
			if (leftNumberFlag == false)
			{
				depth--;				// adjust call depth
				return EXPR_SYNTAX_ERROR;
			}
			if (stackValues[stackIndex] < 0 || stackValues[stackIndex] > 46340)
			{
				depth--;				// adjust call depth
				return EXPR_INTERM_TOO_HIGH;
			}
			long long temp = MulPrToLong(stackValues[stackIndex]);
			mpz_primorial_ui(ZT(factorial), temp);
			stackValues[stackIndex] = factorial;
			exprIndex++;
			continue;      // process more of expr
		}
		else   /* look for function name. If found get its parameters and 
			leave them on the stack. If no name found return a +ve code */
			if ((retcode = func(expr, leftNumberFlag, exprIndex, fcode)) <= 0) {
				if (retcode != EXPR_OK) {
					depth--;				// adjust call depth
					return retcode;  // error processing function parameters
				}
				/* do any further checks needed on the parameter values, then
				evaluate the function. The result is left on the stack. 
				case numbers below must match the code numbers in functionList*/
				switch (fcode) {
				case fn_gcd: {			// GCD	
					mpz_gcd(ZT(stackValues[stackIndex]), ZT(stackValues[stackIndex + 1]), ZT(stackValues[stackIndex]));
					break;
				}
				case fn_modpow: {						// MODPOW
					mpz_powm(ZT(stackValues[stackIndex]), ZT(stackValues[stackIndex]), ZT(stackValues[stackIndex + 1]), ZT(stackValues[stackIndex + 2]));
					break;
				}
				case fn_modinv: {						// MODINV
					/* if an inverse doesn’t exist the return value is zero and rop is undefined*/
					rv = mpz_invert(ZT(stackValues[stackIndex]), ZT(stackValues[stackIndex]), ZT(stackValues[stackIndex + 1]));
					if (rv == 0) {
						depth--;				// adjust call depth
						return EXPR_ARGUMENTS_NOT_RELATIVELY_PRIME;
					}
					break;
				}

				case fn_totient: {			// totient
					if (stackValues[stackIndex] < 1) return EXPR_INVALID_PARAM;
					stackValues[stackIndex] = ComputeTotient(stackValues[stackIndex]);
					break;
				}
				case fn_numdivs: {		// NUMDIVS
					if (stackValues[stackIndex] < 1) {
						depth--;				// adjust call depth
						return EXPR_INVALID_PARAM;
					}
					stackValues[stackIndex] = ComputeNumDivs(stackValues[stackIndex]);
					break;
				}
				case fn_sumdivs: {		// SUMDIVS
					stackValues[stackIndex] = ComputeSumDivs(stackValues[stackIndex]);
					break;
				}

				case fn_sumdigits: {		// SumDigits(n, r) : Sum of digits of n in base r.
					stackValues[stackIndex] = ComputeSumDigits(stackValues[stackIndex],
						stackValues[stackIndex + 1]);
					break;
				}
				case fn_numdigits: {		// numdigits
					stackValues[stackIndex] = ComputeNumDigits(stackValues[stackIndex],
						stackValues[stackIndex + 1]);
					break;
				}
				case fn_revdigits: {	// revdigits
					stackValues[stackIndex] = ComputeRevDigits(stackValues[stackIndex],
						stackValues[stackIndex + 1]);
					break;
				}

				case fn_isprime: {  // isprime
							/* -1 indicates probably prime, 0 = composite */
					stackValues[stackIndex] = PrimalityTest(abs(stackValues[stackIndex]));
					break;
				}
				case fn_fib: {		// fibonacci
					if (stackValues[stackIndex] < 0) {
						depth--;				// adjust call depth
						return EXPR_INVALID_PARAM;
					}
					if (stackValues[stackIndex] > 95700)
					{
						depth--;				// adjust call depth
						return EXPR_INTERM_TOO_HIGH;
					}
					long long temp = MulPrToLong(ZT(stackValues[stackIndex]));
					mpz_fib_ui(ZT(stackValues[stackIndex]), temp);  // calculate fibonacci number
					break;
				}
				case fn_luc: {		// lucas number
					if (stackValues[stackIndex] < 0) {
						depth--;				// adjust call depth
						return EXPR_INVALID_PARAM;
					}
					if (stackValues[stackIndex] > 95700)
					{
						depth--;				// adjust call depth
						return EXPR_INTERM_TOO_HIGH;
					}
					long long temp = MulPrToLong(stackValues[stackIndex]);
					mpz_lucnum_ui(ZT(stackValues[stackIndex]), temp);  // calculate lucas number
					break;
				}

				case fn_part: {              // number of partitions
					if (stackValues[stackIndex] < 0 || stackValues[stackIndex] >= 60000) {
						depth--;				// adjust call depth
						return EXPR_INVALID_PARAM;  // note: biperm is limited to values <= 60,000
					}
					int temp = (int)MulPrToLong(stackValues[stackIndex]);
					biperm(temp, ZT(stackValues[stackIndex]));   // calculate number of partitions
					break;
				}
				case fn_np: {  // next prime;
					mpz_nextprime(ZT(stackValues[stackIndex]), ZT(stackValues[stackIndex]));  // get next prime
					break;
				}
				case fn_pp: {			// previous prime
					retcode = ComputeBack(stackValues[stackIndex], stackValues[stackIndex]);  // get previous prime
					if (retcode != EXPR_OK) {
						depth--;				// adjust call depth
						return retcode;   // error: number < 3
					}
					break;
				}

				case fn_primePi: {  // count primes <= n
					if (stackValues[stackIndex] > max_prime) {
						depth--;				// adjust call depth
						return EXPR_INVALID_PARAM;
					}
					stackValues[stackIndex] = primePi(stackValues[stackIndex]);
					break;
				}
				case fn_concatfact: { /*Concatenates the prime factors of n according to 
						   the mode in m */
					if (stackValues[stackIndex] < 0 || stackValues[stackIndex] > 3) {
						depth--;				// adjust call depth
						return EXPR_INVALID_PARAM;  // mode value invalid
					}
					stackValues[stackIndex] = concatFact(stackValues[stackIndex], 
						stackValues[stackIndex + 1]);
					break;
				}

				default:
					abort();		// if we ever get here we have a problem
				}

				leftNumberFlag = true;
				continue; // contine procesing rest of expr
			}

		else if (charValue == '(')
		{
			if (leftNumberFlag == true)	{
				depth--;				// adjust call depth
				return EXPR_SYNTAX_ERROR;
			}
			if (stackIndex >= PAREN_STACK_SIZE)	{
				depth--;				// adjust call depth
				return EXPR_TOO_MANY_PAREN;
			}
			stackOperators[stackIndex++] = oper_leftb;
			exprIndex++;
			continue;         // process more of expr
		}
		else if (charValue == ')' || charValue == ',')
		{
			if (leftNumberFlag == false)
			{
				depth--;				// adjust call depth
				return EXPR_SYNTAX_ERROR;
			}
			while (stackIndex > startStackIndex &&
				stackOperators[stackIndex - 1] != oper_leftb)
			{  // compute value(s) of sub-expression(s) following left bracket
				if ((SubExprResult = ComputeSubExpr()) != EXPR_OK)
				{
					depth--;				// adjust call depth
					return SubExprResult;  // return error code
				}
			}

			if (stackIndex == startStackIndex)
			{
				break;          // jump out of while loop if done
			}
			if (charValue == ',')
			{
				depth--;				// adjust call depth
				return EXPR_PAREN_MISMATCH;
			}
			stackIndex--;    /* Discard ')' */
			stackValues[stackIndex] = stackValues[stackIndex + 1];
			leftNumberFlag = 1;
			exprIndex++;
			continue;         // process more of expr
		}
		else if (charValue >= '0' && charValue <= '9') 	{
			/* convert number from ascii to BigInteger */
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
				mpz_set_str(ZT(stackValues[stackIndex]), digits.data(), 16);
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
				mpz_set_str(ZT(stackValues[stackIndex]), digits.data(), 10);
				exprIndex = exprIndexAux;
			}
			leftNumberFlag = true;
			continue;      // process more of expr
		}

		/* If we drop through to here, we've got something we can't understand */
		depth--;				// adjust call depth
		return EXPR_SYNTAX_ERROR;
	}                              /* end while */

	if (leftNumberFlag == false)
	{
		depth--;				// adjust call depth
		return EXPR_SYNTAX_ERROR;
	}
	while (stackIndex > startStackIndex && stackOperators[stackIndex - 1] != oper_leftb)
	{
		if ((SubExprResult = ComputeSubExpr()) != EXPR_OK)
		{
			depth--;				// adjust call depth
			return SubExprResult;   // return error code
		}
	}
	if (stackIndex != startStackIndex) 	{
		depth--;				// adjust call depth
		return EXPR_PAREN_MISMATCH;
	}
	ExpressionResult = stackValues[stackIndex];
	depth--;				// adjust call depth
	if (depth > 0 || exprIndex >= expr.size())
		return EXPR_OK;
	else
		return EXPR_SYNTAX_ERROR;  // on final exit but, still have unprocessed text.
}

/* translate error code to text and output it*/
static void textError(enum eExprErr rc)
{
	/*
	currently defined error codes:
	EXPR_NUMBER_TOO_LOW =		-100,
	EXPR_NUMBER_TOO_HIGH,		-99
	EXPR_INTERM_TOO_HIGH,		-98
	EXPR_DIVIDE_BY_ZERO,        -97
	EXPR_PAREN_MISMATCH,		-96
	EXPR_SYNTAX_ERROR,			-95
	EXPR_TOO_MANY_PAREN,		-94
	EXPR_INVALID_PARAM,			-93
	EXPR_ARGUMENTS_NOT_RELATIVELY_PRIME,	-92
	EXPR_BREAK,					-91
	EXPR_OUT_OF_MEMORY,
	EXPR_CANNOT_USE_X_IN_EXPONENT,
	EXPR_DEGREE_TOO_HIGH,
	EXPR_EXPONENT_TOO_LARGE,
	EXPR_EXPONENT_NEGATIVE,				-86
	EXPR_LEADING_COFF_MULTIPLE_OF_PRIME,  not used?
	EXPR_CANNOT_LIFT,					not used?
	EXPR_MODULUS_MUST_BE_GREATER_THAN_ONE,	-83
	EXPR_MODULUS_MUST_BE_PRIME_EXP,		not used?
	EXPR_BASE_MUST_BE_POSITIVE,			-81
	EXPR_POWER_MUST_BE_POSITIVE,		-80
	EXPR_MODULUS_MUST_BE_NONNEGATIVE,	-79
	EXPR_VAR_OR_COUNTER_REQUIRED,		-78
	EXPR_OK = 0
	*/
	switch (rc)
	{
	case EXPR_NUMBER_TOO_LOW:
		std::cout << (lang ? "Número muy pequeño\n" : "Number too low\n");
		break;
	case EXPR_NUMBER_TOO_HIGH:
		std::cout << (lang ? "Número muy grande \n" :
			"Number too high \n");
		break;
	case EXPR_INTERM_TOO_HIGH:
		std::cout << (lang ? "Número intermedio muy grande (más de 20000 dígitos\n" :
			"Intermediate number too high (more than 20000 digits)\n");
		break;
	case EXPR_DIVIDE_BY_ZERO:
		std::cout << (lang ? "División por cero\n" : "Division by zero\n");
		break;
	case EXPR_PAREN_MISMATCH:
		std::cout << (lang ? "Error de paréntesis\n" : "Parenthesis mismatch\n");
		break;
	case EXPR_SYNTAX_ERROR:
		if (lang)
		{
			std::cout << ( "Error de sintaxis\n");
		}
		else
		{
			std::cout << ("Syntax error\n");
		}
		break;
	case EXPR_TOO_MANY_PAREN:
		std::cout << (lang ? "Demasiados paréntesis\n" : "Too many parenthesis\n");
		break;
	case EXPR_INVALID_PARAM:
		std::cout << (lang ? "Parámetro inválido\n" : "Invalid parameter\n");
		break;
	case EXPR_ARGUMENTS_NOT_RELATIVELY_PRIME:
		std::cout << (lang ? "MCD de los argumentos no es 1\n" : "GCD of arguments is not 1\n");
		break;
	case EXPR_BREAK:
		std::cout << (lang ? "Detenido por el usuario\n" : "Stopped by use\nr");
		break;
	case EXPR_EXPONENT_NEGATIVE: {
		std::cout << "Exponent is negative\n";
		break;
	}
	case EXPR_VAR_OR_COUNTER_REQUIRED:
		if (lang)
		{
			std::cout <<  "La expresión \n";
		}
		else
		{
			std::cout << "Expression #\n";
		}
		break;
	case EXPR_BASE_MUST_BE_POSITIVE:
		std::cout << (lang ? "La base debe ser mayor que cero\n" :
			"Base must be greater than zero\n");
		break;
	case EXPR_POWER_MUST_BE_POSITIVE:
		std::cout << (lang ? "La potencia debe ser mayor que cero\n" :
			"Power must be greater than zero\n");
		break;
	case EXPR_MODULUS_MUST_BE_GREATER_THAN_ONE:
		std::cout << (lang ? "El módulo debe ser mayor que 1\n" : "Modulus must be greater than one\n");
		break;
	case EXPR_MODULUS_MUST_BE_NONNEGATIVE:
		std::cout << (lang ? "El módulo no debe ser negativo\n" :
			"Modulus must not be negative\n");
		break;
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
void doFactors(const Znum &Result) {
	std::vector <Znum> factorlist;
	std::vector<int> exponentlist;
	Znum Quad[4];

	/* call DA´s magic function to factorise Result */
	bool rv = factorise(Result, factorlist, exponentlist, Quad);
	if (rv && !factorlist.empty()) {
		if (factorlist.size() > 1 || exponentlist[0] > 1) {
			/* print factor list */
			std::cout << " = ";
			if (Result < 0)
				std::cout << "-";
			for (size_t i = 0; i < factorlist.size(); i++) {
				if (i > 0)
					std::cout << " * ";
				ShowLargeNumber(factorlist[i], 6, false, false);
				if (exponentlist[i] > 1)
					std::cout << "^" << exponentlist[i];
			}
		}
		else
			std::cout << " is prime";  //number has only 1 factor
		if (abs(Result) != 1) {
			auto divisors = NoOfDivisorsF(exponentlist);
			std::cout << "\nNumber of Divisors = ";
			ShowLargeNumber(divisors, 6, false, false);
		}
		else 
			std::cout << "\nNumber of Divisors = 1";   // treat n=1 as special case

		auto divisors = SumOfDivisors(exponentlist, factorlist);
		std::cout << "\nSum of Divisors    = ";
		ShowLargeNumber(divisors, 6, false, false);
		divisors = Totient(exponentlist, factorlist);
		std::cout << "\nTotient            = ";
		ShowLargeNumber(divisors, 6, false, false);
		if (Result > 0) {
			auto mob = mobius(exponentlist);  // mobius only defined for +ve integers
			std::cout << "\nMöbius             = " << mob;
		}

		/* show that the number is the sum of 4 or fewer squares. See
		https://www.alpertron.com.ar/4SQUARES.HTM */
		if (Result >= 0)
			std::cout << "\nn = ";
		else
			std::cout << "\n-n = ";
		char c = 'a';
		for (auto q : Quad) {
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
			c++;
		}
		std::cout << '\n';
	}
	else
		std::cout << " cannot be factorised\n";
}

int main(int argc, char *argv[]) {
	std::string expr, expupper;
	Znum Result, pval, pval2;
	eExprErr rv;

	char helpmsg[] =
		"You can enter expressions that use the following operators, functions and parentheses:\n"

		"+       : addition\n"
		"-       : subtraction\n"
		"*       : multiplication\n"
		"/       : integer division\n"
		"%%       : modulus (remainder from integer division)\n"
		"^ or ** : exponentiation (the exponent must be greater than or equal to zero).\n"
		"<, ==, >, <=, >=, != for comparisons. The operators return zero for false and -1 for true.\n"
		"AND, OR, XOR, NOT for binary logic.\n"
		"SHL or <<: Shift left the number of bits specified on the right operand.\n"
		"SHR or >>: Shift right the number of bits specified on the right operand.\n"
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
		"Also the following commands: X=hexadecimal o/p, D=decimal o/p \n"
		"F = do factorisation, N = Don't factorise, S = Spanish, E=English\n"
		"HELP (this message) and EXIT\n";

	char ayuda[] =
		"Puedes ingresar expresiones que usen los siguientes operadores y paréntesis:\n"
		"+ para suma\n"
		"- para resta\n"
		"* para multiplicación\n"
		"/ para división entera\n"
		"% para el resto de la división entera\n"
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
		setlocale(LC_ALL, "en-EN");  // allows non-ascii characters to print
		char banner[] = "Compiled on "  __DATE__ " at " __TIME__ "\n";
		printf("%s", banner);

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
					printf(helpmsg);
				else
					printf(ayuda);
				continue;
			}
			if (expupper == "E") { lang = 0; continue; }       // english
			if (expupper == "S") { lang = 1; continue; }	   // spanish (Español)
			if (expupper == "F") { factorFlag = true; continue; }  // do factorisation
			if (expupper == "N") { factorFlag = false; continue; } // don't do factorisation
			if (expupper == "X") { hex = true; continue; } // don't do factorisation
			if (expupper == "D") { hex = false; continue; } // don't do factorisation

			auto start = clock();
			removeBlanks(expr);     // remove any spaces 
			
			rv = ComputeExpr(expr, Result); /* analyse expression, compute value*/

			if (rv != EXPR_OK) {
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
				ShowLargeNumber(Result, 6, true, hex);
				std::cout << '\n';				
				if (factorFlag) {
					doFactors(Result); /* factorise Result, calculate number of divisors etc */
				}
			}

			auto end = clock();   // measure amount of time used
			double elapsed = (double)end - start;
			std::cout << "time used= " << elapsed / CLOCKS_PER_SEC << " seconds\n";
		}

		return EXIT_SUCCESS;
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
