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
#pragma fenv_access (on)

#include "pch.h"
#define WIN32_LEAN_AND_MEAN
#include <Windows.h>
#include <Mmsystem.h >   // for sound effects
#include "factor.h"

#include "diagnostic.h"

#define BIGNBR       // define to include bignbr tests 
#ifdef BIGNBR
#include "bignbr.h"
#include "bigint.h"
extern Znum zR, zR2, zNI, zN;
#endif

typedef unsigned __int32 uint32;

/* external function declaration */
void VersionInfo(const LPCSTR path, int ver[4], std::string& modified);
char* getFileName(const char* filter, HWND owner);
retCode ComputeExpr(const std::string& expr, Znum& Result, int& asgCt, bool* multiV = nullptr);

int lang = 0;             // 0 English, 1 = Spanish

bool hex = false;		// set true if output is in hex
int factorFlag = 2;     /* 0 = no factorisation, 1 = factorisation but not totient etc,
						   2 = get totient, number of divisors etc after factorisation */
/* verbose value is used to turn off or on optional messages; 
higher value = more messages */
#ifdef _DEBUG
int verbose = 1;       /* default 1 if compiled in debug mode */
#else
int verbose = 0;
#endif

HANDLE hConsole;       /* used by SetConsoleCursorPosition() function */
HWND handConsole;      /* handle to console window */

bool *primeFlags = NULL;
unsigned long long *primeList = NULL;
unsigned int prime_list_count = 0;
unsigned long long int primeListMax = 0;

std::vector <std::string> exprList;  /* expressions stored here as text once
									   they are validated */

/* name of sound file played at end of processing a command or expression
 if the elapsed time > 10 seconds */
std::string endsound = "c:/Windows/Media/Alarm09.wav";

/* name of sound file played when prompt for input is displayed */
std::string attsound = "c:/Windows/Media/chimes.wav";

/* function declarations, only for functions that have forward references */
static void textError(retCode rc);


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
	else {
		StackTrace2();
		throw std::range_error("big number cannot be converted to 64-bit integer");
	}
	return 0;
}

/* NumDigits(n,r): Number of digits of n in base r. Leading zeros are not counted  */
long long ComputeNumDigits(const Znum &n, const Znum &radix)
{
	/* more accurate than mpz_sizeinbase */
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
Znum llt(const Znum &p) {
	int rv;
	Znum n, tmp; 
	long long exp = MulPrToLong(p); /* exp = p */
	clock_t t1, t2, t3;

	/* 1st check if p is prime. */
	rv = mpz_bpsw_prp(ZT(p));  /* returns 0, 1 or 2*/
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
				printf_s("\b");    // erase previous output
			nchars = printf_s("%s llt iteration %lld %.2f%% complete", myTime(), i,
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


/* get floor(sqrt(n))*/
unsigned long long llSqrt(const unsigned long long n) {
	double root = sqrt((double)n);
	unsigned long long  iroot = llround(root);
	while (iroot*iroot > n)
		iroot--;
	return iroot;
}

// return value of bit in array corresponding to x. false (0) indicates a prime number
bool getBit(const unsigned long long int x, const bool array[])
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

	// after completing the while loop we have found all the primes < max_val
	if (verbose > 0)
		printf_s("  prime %9lld is %11lld\n", count, numsave);
	primeList[count] = ULLONG_MAX;		// set end marker
	prime_list_count = (unsigned int)count;
	primeListMax = primeList[count - 1];
	return;
}



/* translate error code to text and output it*/
static void textError(retCode rc) {
	/*
	error codes currently used:
	NUMBER_TOO_LOW,
	NUMBER_TOO_HIGH,
	INTERM_TOO_HIGH,		
	DIVIDE_BY_ZERO,        
	PAREN_MISMATCH,		
	SYNTAX_ERROR,			
	TOO_MANY_PAREN,		
	INVALID_PARAM,	
	ARGUMENTS_NOT_RELATIVELY_PRIME
	EXPONENT_TOO_LARGE,
	EXPONENT_NEGATIVE,				
	EXPR_OK = 0
	*/
	switch (rc)
	{
	case retCode::NUMBER_TOO_LOW:
		std::cout << (lang ? "Número muy pequeño\n" : "Number too low\n");
		break;
	case retCode::NUMBER_TOO_HIGH:
		std::cout << (lang ? "Número muy grande \n" :
			"Number too high \n");
		break;
	case retCode::INTERM_TOO_HIGH:
		std::cout << (lang ? "Número intermedio muy grande (más de 20000 dígitos\n" :
			"Intermediate number too high (more than 20000 digits)\n");
		break;
	case retCode::DIVIDE_BY_ZERO:
		std::cout << (lang ? "División por cero\n" : "Division by zero\n");
		break;
	case retCode::PAREN_MISMATCH:
		std::cout << (lang ? "Error de paréntesis\n" : "Parenthesis mismatch\n");
		break;
	case retCode::SYNTAX_ERROR:
		if (lang) 	{
			std::cout << ( "Error de sintaxis\n");
		}
		else {
			std::cout << ("Syntax error\n");
		}
		break;
	case retCode::TOO_MANY_PAREN:
		std::cout << (lang ? "Demasiados paréntesis\n" : "Too many parenthesis\n");
		break;
	case retCode::INVALID_PARAM:
		std::cout << (lang ? "Parámetro inválido\n" : "Invalid parameter\n");
		break;
	case retCode::ARGUMENTS_NOT_RELATIVELY_PRIME:
		std::cout << (lang ? "MCD de los argumentos no es 1\n" : "GCD of arguments is not 1\n");
		break;
	/*case EXPR_BREAK:
		std::cout << (lang ? "Detenido por el usuario\n" : "Stopped by use\nr");
		break;*/
	case retCode::EXPONENT_NEGATIVE: {
		std::cout << "Exponent is negative\n";
		break;
	}
	case retCode::EXPONENT_TOO_LARGE: {
		std::cout << "Exponent exceeds 2^31-1\n";
		break;
	}
	/*case retCode::EXPR_VAR_OR_COUNTER_REQUIRED:
		if (lang)
		{
			std::cout <<  "La expresión \n";
		}
		else
		{
			std::cout << "Expression #\n";
		}
		break;*/
	//case retCode::EXPR_BASE_MUST_BE_POSITIVE:
	//	std::cout << (lang ? "La base debe ser mayor que cero\n" :
	//		"Base must be greater than zero\n");
	//	break;
	//case retCode::EXPR_POWER_MUST_BE_POSITIVE:
	//	std::cout << (lang ? "La potencia debe ser mayor que cero\n" :
	//		"Power must be greater than zero\n");
	//	break;
	//case retCode::EXPR_MODULUS_MUST_BE_GREATER_THAN_ONE:
	//	std::cout << (lang ? "El módulo debe ser mayor que 1\n" : "Modulus must be greater than one\n");
	//	break;
	//case retCode::EXPR_MODULUS_MUST_BE_NONNEGATIVE:
	//	std::cout << (lang ? "El módulo no debe ser negativo\n" :
	//		"Modulus must not be negative\n");
	//	break;
	default:
		printf_s( "unknown error code: %d\n", (int)rc);
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
		PlaySoundA(endsound.c_str(), NULL,
			SND_FILENAME | SND_NODEFAULT | SND_ASYNC | SND_NOSTOP);

	if (elSec <= 60.0) {
		printf_s("%.4f seconds \n", elSec);  /* print time used to nearest millisecond */
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

/* remove spaces, tabs, etc  from msg (\t, \r, \n, \v   and \f)*/
static void removeBlanks(std::string &msg) {
	for (size_t ix = 0; ix < msg.size(); ix++) {
		if ((unsigned char)msg[ix] <= 0x7f && isspace(msg[ix])) {     // look for spaces, tabs, etc
			msg.erase(ix, 1);      // remove space character
			ix--;  // adjust index to take account of removed blank
		}
	}
}

/* remove initial & trailing spaces, tabs, etc from msg (\t, \r, \n, \v   and \f) */
static void removeInitTrail(std::string &msg) {
	while (!msg.empty() && isspace(msg.front())) {
		msg.erase(0, 1);      // remove 1st space character
	}
	while (!msg.empty() && isspace(msg.back())) {
		msg.resize(msg.size() - 1);  /* remove trailing spaces */
	}
}

struct summary {
	int numsize;		// number of decimal digits in number
	double time;		// time used in seconds
	int NumFacs;		// number of unique factors
	int totalFacs;		// total number of factors (=1 for a prime number)
	int testNum;		// test number (if applicable)
	int sndFac;		    // number of decimal digits in 2nd largest factor
	struct counters ctrs;
};

std::vector <summary> results;

/* print summary - 1 line per test */
static void printSummary(void) {
	long long sec, min, hour;
	double elSec;
	/* print column headings */
	printf_s("Test Num Size   time      Unique Factors Total Factors     2nd Fac");
	printf_s(" tdv prh leh crm pm1 ecm siq pwr yaf msv \n");
	for (auto res : results) {
		/* truncate elapsed time to nearest second */
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

		/* print counters showing how factors were found */
		printf_s("     %8d   %8d     %8d", res.NumFacs, res.totalFacs, res.sndFac);
		printf_s("    %3d %3d %3d %3d %3d %3d ",
			res.ctrs.tdiv, res.ctrs.prho, res.ctrs.leh, res.ctrs.carm,
			res.ctrs.pm1, res.ctrs.ecm);
		printf_s("%3d %3d %3d %3d \n", res.ctrs.siqs, res.ctrs.power, res.ctrs.yafu,
			res.ctrs.msieve);
	}
}

/* factorise Result, calculate number of divisors etc and print results */
static void doFactors(const Znum &Result, bool test) {
	fList factorlist;
	Znum Quad[4];
	clock_t start;
	summary sum;    // save summary of factorisation
	bool rv;

	if (test)
		start = clock();	// used to measure execution time

	/* call DA´s magic function to factorise Result */
	if (factorFlag > 1)
		rv = factorise(Result, factorlist, Quad);
	else
		rv = factorise(Result, factorlist, nullptr);
	if (rv && factorlist.fsize() > 0) {
		factorlist.print(Result < 0);   /* print factors */
		std::cout << '\n';
		if (factorFlag > 1) {
			auto divisors = factorlist.NoOfDivs();
			std::cout << "Number of Divisors = ";
			ShowLargeNumber(divisors, 6, false, false);

			divisors = factorlist.DivisorSum();
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
		}
		
		if (test) {
			/* recalculate result & get total number of factors */
			Znum result = factorlist.recheck(sum.totalFacs);
			if (result != Result) {
				std::cout << "Factors expected value " << Result << " actual value " << result << '\n';
				Beep(750, 1000);
			}
			if (factorFlag > 1) {
				result = Quad[0] * Quad[0] + Quad[1] * Quad[1] + Quad[2] * Quad[2] + Quad[3] * Quad[3];
				if (result != Result) {
					std::cout << "Quad expected value " << Result << " actual value " << result << '\n';
					Beep(750, 1000);
				}
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
			sum.ctrs = factorlist.getCtrs();
			results.push_back(sum);
		}
	}
	else
		std::cout << " cannot be factorised\n";
}

/* perform some simple tests. Returns true if x3 is prime 
method = 0 for standard factorisation, != 0 to use only YAFU for factorisation */
static bool factortest(const Znum &x3, const int testnum, const int method=0) {
	fList factorlist;
	Znum Quad[4], result;
	long long totalFactors = 0;
	summary sum;    // save summary of test

	auto start = clock();	// used to measure execution time
	double end, elapsed;

	sum.numsize = (int)ComputeNumDigits(x3, 10);

	std::cout << "\nTest " << testnum << ": factorise ";
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

		if (method == 0)
			factorlist.prCounts();   // print counts
		else
			sum.ctrs.yafu = sum.totalFacs;

		std::cout << "test " << testnum << " completed at ";

		end = clock();              // measure amount of time used
		elapsed = (double)end - start;
		PrintTimeUsed(elapsed, "time used = ");
		sum.time = elapsed / CLOCKS_PER_SEC;
		sum.NumFacs = (int)factorlist.fsize();
		sum.testNum = testnum;
		if (method == 0)
			sum.ctrs = factorlist.getCtrs();
		results.push_back(sum);
		return false;   // not prime
	}
	else {
		std::cout << "is prime \n";
		std::cout << "test " << testnum << " completed at ";
		sum.NumFacs = 1;
		sum.totalFacs = 1;
		sum.sndFac = 0;
		sum.testNum = testnum;
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
	int i, asgCt;
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
		"29!!!",                  72642169600,   // triple factorial
		"33!!!!",                  4996616625,
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
		"le(22, 7)",                        1,
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
		"(NOT 0x0f23456789abcde1) * -1",  0x0f23456789abcde2,
		"5 < 6 == 7 < 8",                -1,   // returns true (== has lower priority)
		"5 < 6 != 7 < 8",                 0,   // returns false (!= has lower priority)
		"5 < (6 != 7) < 8",              -1,   // returns true; expr evaluated from left to right
		"R3(49)",                        54,
		"R2(585)",                       16,
		"SQRT(1234320)",               1110,
		"NROOT(2861381721051424,5)",   1234,
		"LLT(3217)",                      1,  // 2^3217-1 is prime
		"BPSW(2^99-1)",                   0,  // not a prime number
		"BPSW(2^127-1)",                  1,  // a prime number
		"ISPRIME(2^127-1)",              -1,  // a prime number
		"aprcl(2^127-1)",                 2,  // a prime number
		"ispow(2^127-1)",                 0,  /* not a perfect power */
		"ispow(2^127)",                  -1,  /* a perfect power */
		"-not1",                          2,  // operators are processed from right to left
		"not-1",                          0,  // operators are processed from right to left
		"not5#",                        -31,  // # operator evaluated before not
		"4^5#",         1152921504606846976,  // # operator evluated before exponent
		"5!!#",                       30030,  // !! operator evaluated before #
		"5#!!",           42849873690624000,  // # operator evaluated before !!
		"$x = 99 ",                      99,  /* test user variables */
		"$y = $x+1  ",                  100,
		"modsqrt(2191, 23^3)",         1115,
		"modsqrt(4142, 29^3)",         2333,
		"modsqrt(3, 143)",               17,
		"minfact(99)",                    3,
		"maxfact(99)",                    11,
		"numfact(99)",                     2,
		"lcm(12,20)",                     60,
		"pi(500)",                        95,
		"primroot(761)",                   6,   /* primitive root */
	};

	results.clear();

	auto start = clock();	// used to measure execution time
	for (i = 0; i < sizeof(testvalues) / sizeof(testvalues[0]); i++) {

		auto  rv =ComputeExpr(testvalues[i].text, result, asgCt);
		if (rv != retCode::EXPR_OK || result != testvalues[i].expected_result) {
			std::cout << "test " << i + 1 << " failed\n" <<
				"expected " << testvalues[i].text << " = " 
				<< testvalues[i].expected_result << '\n';
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
		testcnt++;
		factortest(x3, testcnt);

		x3++;
		testcnt++;
		factortest(x3, testcnt);
	}

	/* tests below have shown a problem with pollard-rho for certain numbers */
	long long int PollRho[] = {99999999973789, 183038861417, 183475587821,
		181919916457, 199996999457, 204493418837 };
	for (int i = 0; i < sizeof(PollRho) / sizeof(PollRho[0]); i++) {
		testcnt++;
		factortest(PollRho[i], testcnt);
	}


	/* exercise code specifically for power +/-1 */
	mpz_ui_pow_ui(ZT(x3), 10, 20);  // x3 = 10^20
	x3 -= 1;                        // x3 = 10^20-1
	testcnt++;
	factortest(x3, testcnt);

	x3 += 2;
	testcnt++;
	factortest(x3, testcnt);

	testcnt++;
	ComputeExpr("n(10^10)^2-1", x3, asgCt);  // test power of large number-1
	factortest(x3, testcnt);

	testcnt++;
	ComputeExpr("120#-1", x3, asgCt);
	factortest(x3, testcnt);

	testcnt++;
	ComputeExpr("n(10^15)^2", x3, asgCt);  // test power of large number
	factortest(x3, testcnt);

	testcnt++;
	ComputeExpr("n(10^6+20)^1667", x3, asgCt);  // test power of large prime number
	factortest(x3, testcnt);

	testcnt++;
	ComputeExpr("n(10^7)^3*n(10^8)^2", x3, asgCt);  // test powers of two large prime number
	factortest(x3, testcnt);

	testcnt++;
	ComputeExpr("n(3*10^9+50)*n(3*10^10+500)", x3, asgCt);  // test Lehman factorisation
	factortest(x3, testcnt);

	testcnt++;
	ComputeExpr("n(10^15)^3*n(10^14)", x3, asgCt);  // test Lehman factorisation
	factortest(x3, testcnt);

	/* test using carmichael numbers.  */
	unsigned long long int carmichael[] = { 90256390764228001, 18118463305678126129, 
		18265521244069461529,  18349357898532971521, 18308657203978189969 };
	for (int i = 0; i < sizeof(carmichael) / sizeof(carmichael[0]); i++) {
		testcnt++;
		factortest(carmichael[i], testcnt);
	}

	ComputeExpr("16344221851913485532689", x3, asgCt);
	testcnt++;
	factortest(x3, testcnt);  /* 23 digit carmichael number */

	ComputeExpr("56897193526942024370326972321", x3, asgCt);
	testcnt++;
	factortest(x3, testcnt);  /* 29 digit pseudo-prime number */

	/* set x3 to large prime. see https://en.wikipedia.org/wiki/Carmichael_number */
	ComputeExpr("2967449566868551055015417464290533273077199179985304335099507"
		"5531276838753171770199594238596428121188033664754218345562493168782883", x3, asgCt);
	x4 = x3 * (313 * (x3 - 1) + 1) * (353 * (x3 - 1) + 1);
	/* in general numbers > about 110 digits cannot be factorised in a reasonable time 
	but this one can, because a special algorithm just for Carmichael numbers is used. */
	testcnt++;
	factortest(x4, testcnt);   // 397-digit Carmichael number
	std::cout << "factorised 397-digit Carmichael number \n";

	testcnt++;
	ComputeExpr("n(10^24)*n(10^25)*n(10^26)*n(10^27)", x3, asgCt);  
	factortest(x3, testcnt);

	testcnt++;
	ComputeExpr("n(2^97)*n(2^105)", x3, asgCt);
	factortest(x3, testcnt);

	auto end = clock();   // measure amount of time used
	double elapsed = (double)end - start;
	PrintTimeUsed(elapsed, "Factorisation tests completed. Time used= ");
	printSummary();
}


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
	while (mpz_bpsw_prp(ZT(r)) == 0)
	//while (!mpz_likely_prime_p(ZT(r), state, 0))
		r += t2;


	// 3. Compute p0 = 2(sr-2 mod r)s - 1.
	p0 = ((s*r - 2) % r)*s * 2 - 1;

	// 4. Find the first prime p in the sequence p0, p0 +2rs p0+4rs ....,
	p = p0;
	while (mpz_bpsw_prp(ZT(p)) == 0)
	//while (!mpz_likely_prime_p(ZT(p), state, 0))
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

/* generate large random number, up to 128 bits */
static void largeRand(Znum &a) {
	a = ((long long)rand() << 32) + rand();
	a <<= 64;
	a += ((long long)rand() << 32) + rand();
}

/* generate extra large number, about size*32 bits */
static void XlargeRand(Znum& a, int size) {
	a = 1;
	for (int c = 1; c <= size; c++) {
		a *= rand();
	}
}

/*  1. check basic arithmetic operators for BigIntegers
	2. test BigInteger multiplication with larger numbers
	3. BigInteger division with larger numbers
	4 & 5. Modular Multiplication using Mongomery Encoding (REDC)
*/
static void doTests3(void) {
	Znum a, a1, am, b, b1, bm, mod, p, p2, pm;
	limb aL[MAX_LEN], modL[MAX_LEN], alM[MAX_LEN], al2[MAX_LEN];
	limb bL[MAX_LEN], blM[MAX_LEN], pl[MAX_LEN], plm[MAX_LEN];
	limb one[MAX_LEN];
	int numLen = MAX_LEN-2, l2;

	BigInteger aBI, a1BI, bBI, amBI, pBI;

	auto start = clock();	// used to measure execution time

	memset(one, 0, MAX_LEN * sizeof(limb));
	one[0] = 1;                   /* set value of one to 1 */

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
		BigtoZ(am, amBI);       // convert back to Znum
		error = p - am;         // get error
		if (error != 0) {
			double e1, e2, relErrf;
			long e1l, e2l;
			e1 = mpz_get_d_2exp(&e1l, ZT(p));
			e2 = mpz_get_d_2exp(&e2l, ZT(am));
			assert(e1l == e2l);   // check power of 2 exponent is correct
			relErrf = abs(e1 - e2) / e1;   /* should be less than 10E-14 */
			relerror = (10'000'000'000'000'000LL * error) / p;  /* error ratio * 10^15 */
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
	
	modmultCallback = nullptr;      // turn off status messages from modmult
	for (int c = 1; c <= 100; c++) {
		/* set up modulus and Mongomery parameters */
		XlargeRand(mod, c);					// get large random number (up to 32 * c bits)
		mod |= 1;                       // set lowest bit (make sure mod is odd)
		GetMontgomeryParms(mod);
		numLen = MAX_LEN - 2;
		ZtoLimbs(modL, mod, MAX_LEN);    // copy value of mod to modL
		while (modL[numLen - 1] == 0)
			numLen--;                    // adjust length i.e. remove leading zeros
		memcpy(TestNbr, modL, numLen * sizeof(limb));  // set up for GetMontgomeryParms
		NumberLength = numLen;
		GetMontgomeryParms(numLen);
		XlargeRand(a, c);				     // get large random number a
		a %= mod;						 // ensure a < mod
		modmult(a, zR2, am);             // convert a to Montgomery (Znum) in am
		numLen = MAX_LEN - 2;
		ZtoLimbs(aL, a, numLen);		     // copy value of a to aL (limbs)
		//while (aL[numLen - 1] == 0)
		//	numLen--;                    // adjust length i.e. remove leading zeros
		//NumberLength = numLen;
		modmult(aL, MontgomeryMultR2, alM);  // convert a to Mongomery (limbs)
		modmult(alM, one, al2);          // convert a from Mongomery (limbs) 
		LimbstoZ(al2, a1, numLen);       // copy value to a1 (Znum)
		assert(a == a1);                 // check that all these conversions work properly
	}
	end = clock();              // measure amount of time used
	elapsed = (double)end - start;
	std::cout << "test stage 4 completed  time used= " << elapsed / CLOCKS_PER_SEC << " seconds\n";

	largeRand(b);				     // get large random number  b
	b |= 1;                         /* make sure b is odd */
	GetMontgomeryParms(b);
	modmult(b, zR2, bm);             // convert b to Montgomery (Znum)
	
	ZtoLimbs(bL, b, numLen);		 // copy value of b to bL
	numLen = MAX_LEN - 2;
	while (bL[numLen - 1] == 0)
		numLen--;                    // adjust length i.e. remove leading zeros
	NumberLength = numLen;
	memcpy(TestNbr, bL, numLen * sizeof(limb));  // set up for GetMontgomeryParms
	GetMontgomeryParms(numLen);
	modmult(bL, MontgomeryMultR2, blM);  // convert b to Mongomery (limbs)

	for (int i = 1; i < 200000000; i++) {
		modmult(alM, blM, plm);                   // p = a*b mod m (limbs)
		memcpy(alM, plm, numLen * sizeof(limb));  // a = p (limbs)
		modmult(am, bm, pm);                      // p = a*b mod m (Znum)
		am = pm;						          //a = p (Znum)
		if (i % 20000000 == 0) {
			std::cout << "test stage 5 " << i / 2000000 << "% complete \n";
		}
	}

	REDC(p, pm);					 // convert p from montgomery (Znum)
	modmult(plm, one, pl);           // convert p from Mongomery (limbs)
	LimbstoZ(pl, p2, numLen);        // convert value of p to p2 (Znum)
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
		if (factortest(m, px+1)) /* factorise m, calculate number of divisors etc */
			std::cout << "2^" << primeList[px] << " - 1 is prime \n";
	}
	auto end = clock();              // measure amount of time used
	auto elapsed = (double)end - start;
	PrintTimeUsed(elapsed, "tests completed time used = ");
	printSummary();           // print summary 1 line per test
}

/* tests using only YAFU for factorisation */
static void doTests5(void) {
	results.clear();
	int testcnt = 1, asgCt;

	Znum num = 49728103; // 7001 * 7103
	factortest(num, testcnt, 1);  // test 1

	// factorise 127 digit number   test 2
	/*  =                               280673 
	*                              2756 163353 
	*                            598990 818061 
	*                          4 527716 228491 
	*                     248158 049830 971629 
	* 33637 310674 071348 724927 955857 253537  
	*117445 937227 520353 139789 517076 610399  
	(7 factors)   */
	testcnt++;
	ComputeExpr("2056802480868100646375721251575555494408897387375737955882170045672576386016591560879707933101909539325829251496440620798637813", 
		num, asgCt);
	factortest(num, testcnt, 1);

	//factorise 57 digit number      test 3
	/* P6 = 280673
	  P12 = 598990818061
	  P10 = 2756163353
	  P13 = 4527716228491
	  P18 = 248158049830971629
	*/
	testcnt++;
	ComputeExpr("520634955263678254286743265052100815100679789130966891851", num, asgCt);
	factortest(num, testcnt, 1);

	//factorise 80 digit number (about 3 minutes)   test 4
	/* P49 = 2412329883909990626764837681505523221799246895133
       P32 = 18138544144503826077310252140817
    */
	testcnt++;
	ComputeExpr("43756152090407155008788902702412144383525640641502974083054213255054353547943661", num, asgCt);
	factortest(num, testcnt, 1);

	//factorise 85 digit number (about 7 mins)   // test 5
	/* factors are 1485325304578290487744454354798448608807999 and 
                   1263789702211268559063981919736415575710439 */
	testcnt++;
	ComputeExpr("1877138824359859508015524119652506869600959721781289179190693027302028679377371001561", num, asgCt);
	factortest(num, testcnt, 1);

	// factorise 94 digit number (about 60 mins)    test 6
	/* factors are 10910042366770069935194172108679294830357131379375349 and 
                   859735020008609871428759089831769060699941 */
	testcnt++;
	ComputeExpr("9379745492489847473195765085744210645855556204246905462578925932774371960871599319713301154409", num, asgCt);
	factortest(num, testcnt, 1);

	//factorise 100 digit number - takes many hours    test 7
	/* factor are 618162834186865969389336374155487198277265679 and
	              4660648728983566373964395375209529291596595400646068307 */
	testcnt++;
	ComputeExpr("2881039827457895971881627053137530734638790825166127496066674320241571446494762386620442953820735453", num, asgCt);
	factortest(num, testcnt, 1);

	printSummary();    // print summary - 1 line per test

	return;
}

/* tests using Msieve for factorisation. Factorise selected Mersenne numbers 
if useMsieve = false use only YAFU. Normal yafu & msieve flags are ignored. */
static void doTests6(bool useMsieve = true) {
	bool yafusave = yafu;
	bool msievesave = msieve;
	Znum m;
	msieve = useMsieve;
	yafu = !useMsieve;
	int testcnt = 1;

	results.clear();
	auto start = clock();	// used to measure execution time

	mpz_ui_pow_ui(ZT(m), 2, 277);  // get  m= 2^p
	m--;                // get 2^p -1
	factortest(m, testcnt);      // 84 digits
	/* p7  factor: 1121297
	   p38 factor: 31133636305610209482201109050392404721
	   p40 factor: 6955979459776540052280934851589652278783
	*/

	testcnt++;
	mpz_ui_pow_ui(ZT(m), 2, 293);  // get  m= 2^p
	m--;                // get 2^p -1
	factortest(m, testcnt);      // 89 digits
	/* p26 factor: 40122362455616221971122353
       p63 factor: 396645227028138890415611220710757921643910743103031701971222447 */

	testcnt++;
	mpz_ui_pow_ui(ZT(m), 2, 311);  // get  m= 2^p
	m--;                // get 2^p -1
	factortest(m, testcnt);      // 94 digits
	/* p7  factor: 5344847
	   p31 factor: 2647649373910205158468946067671
	   p57 factor: 294803681348959296477194164064643062187559537539328375831
	*/

	testcnt++;
	mpz_ui_pow_ui(ZT(m), 2, 313);  // get  m= 2^p
	m--;                // get 2^p -1
	factortest(m, testcnt);      // 95 digits
	/* p8  factor: 10960009
	   p17 factor: 14787970697180273
	   p25 factor: 3857194764289141165278097
	   p47 factor: 26693012026551688286164949958620483258358551879
	*/

	testcnt++;
	mpz_ui_pow_ui(ZT(m), 2, 349);  // get  m= 2^p
	m--;                // get 2^p -1
	factortest(m, testcnt);      // 106 digits
	/* p34 factor: 2927455476800301964116805545194017
	   p67 factor: 6725414756111955781503880188940925566051960039574573675843402666863
    */
	
	yafu = yafusave;
	msieve = msievesave;

	auto end = clock();              // measure amount of time used
	auto elapsed = (double)end - start;
	PrintTimeUsed(elapsed, "tests completed time used = ");
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
	newStr << "endsound=" << endsound << '\n';
	newStr << "attsound=" << attsound << '\n';
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
path definitions begin with yafu-path=, yafu-prog=, msieve-path=, msieve-prog=, 
helpfile=, endsound= or attsound=
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
			else if (_strnicmp("endsound=", buffer.c_str(), 9) == 0) {
				endsound = buffer.substr(9);
			}
			else if (_strnicmp("attsound=", buffer.c_str(), 9) == 0) {
				attsound = buffer.substr(9);
			}
			else inifile.push_back(buffer);  // save anything not recognised
		}
		iniStr.close();
	}
}


//search the docfile for the right entry specified by helpTopic
//just search for the heading, and print everything until
//the next heading is found
static void helpfunc(const std::string &helpTopic)
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
		newpathC = getFileName("Text\0*.TXT\0\0", handConsole);
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
		printf_s("searching for help on '%s'\n", helpTopic.data());

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
			if (strstr(helpTopic.data(), str + 1) != NULL)
				/* we get a match if the topic between [ and ] is contained 
				anywhere in helptopic */
				printtopic = true;   /* we have found the required topic*/
		}
		else {  /* not a header line */
			if (printtopic)
				printf_s("%s", str);  /* print only if within the required topic */
		}
	}

	if (feof(doc)) {
		if (printtopic)
			putchar('\n');   /* contrary to the POSIX standard, the last line of the file
							may not end with newline */
		else
			printf_s("Help for %s not found \n", helpTopic.data());
	}
	fclose(doc);
	return;
}

// machine info
double MEAS_CPU_FREQUENCY;
#ifdef _DEBUG
int VERBOSE_PROC_INFO = 1;
#else 
int VERBOSE_PROC_INFO = 0;
#endif
char CPU_ID_STR[80];

int CLSIZE;
char HAS_SSE41;
#if defined(WIN32)
char sysname[MAX_COMPUTERNAME_LENGTH + 1];
unsigned long sysname_sz;
#else
char sysname[256];
int sysname_sz;
#endif

#ifdef _MSC_VER

/* Core aware timing on Windows, courtesy of Brian Gladman */

#if defined( _WIN64 )

#define current_processor_number GetCurrentProcessorNumber

#else

unsigned long current_processor_number(void)
{
	__asm
	{
		mov     eax, 1
		cpuid
		shr     ebx, 24
		mov     eax, ebx
	}
}

#endif

static int lock_thread_to_core(void)
{
	DWORD_PTR afp, afs;

	if (GetProcessAffinityMask(GetCurrentProcess(), &afp, &afs))
	{
		afp &= (DWORD_PTR)(1LL << current_processor_number());
		if (SetThreadAffinityMask(GetCurrentThread(), afp))
			return EXIT_SUCCESS;
	}
	return EXIT_FAILURE;
}

static int unlock_thread_from_core(void)
{
	DWORD_PTR afp, afs;

	if (GetProcessAffinityMask(GetCurrentProcess(), &afp, &afs))
	{
		if (SetThreadAffinityMask(GetCurrentThread(), afp))
			return EXIT_SUCCESS;
	}
	return EXIT_FAILURE;
}

double cycles_per_second = 0.0;
double ticks_per_second = 0.0;
double cycles_per_tick = 0.0;

unsigned long long measure_processor_speed(void)
{
	unsigned long long cycles;

	lock_thread_to_core();
	cycles = __rdtsc();
	Sleep(100);
	cycles = __rdtsc() - cycles;
	unlock_thread_from_core();
	cycles_per_second = 10.0 * (double)cycles;

	if (ticks_per_second == 0.0)
	{
		LARGE_INTEGER ll;
		QueryPerformanceFrequency(&ll);
		ticks_per_second = (double)ll.QuadPart;
		cycles_per_tick = cycles_per_second / ticks_per_second;
	}
	return cycles;
}


#else

uint64 measure_processor_speed(void)
{
	uint64 cycles;
	struct timeval start, stop;
	double t_time;
	TIME_DIFF* difference;

	gettimeofday(&start, NULL);

	cycles = yafu_read_clock();
	do
	{
		gettimeofday(&stop, NULL);
		difference = my_difftime(&start, &stop);
		t_time = ((double)difference->secs +
			(double)difference->usecs / 1000000);
		free(difference);
	} while (t_time < 0.1);
	cycles = yafu_read_clock() - cycles;

	return cycles;                  /* return cycles per second  */
}

#endif

double get_tsc_time(void)
{
	if (cycles_per_second == 0.0)
		measure_processor_speed();
	return __rdtsc() / cycles_per_second;
}

double get_pfc_time(void)
{
	LARGE_INTEGER ll;

	if (ticks_per_second == 0.0)
		measure_processor_speed();
	QueryPerformanceCounter(&ll);
	return ll.QuadPart / ticks_per_second;
}

#if defined(GCC_ASM32X)
#define HAS_CPUID
#define CPUID(code, a, b, c, d) 			\
		ASM_G volatile(					\
			"movl %%ebx, %%esi   \n\t"		\
			"cpuid               \n\t"		\
			"movl %%ebx, %1      \n\t"		\
			"movl %%esi, %%ebx   \n\t"		\
			:"=a"(a), "=m"(b), "=c"(c), "=d"(d) 	\
			:"0"(code) : "%esi")
#define CPUID2(code1, code2, a, b, c, d) 			\
		ASM_G volatile(					\
			"movl %%ebx, %%esi   \n\t"		\
			"cpuid               \n\t"		\
			"movl %%ebx, %1      \n\t"		\
			"movl %%esi, %%ebx   \n\t"		\
			:"=a"(a), "=m"(b), "=c"(c), "=d"(d) 	\
			:"0"(code1), "2"(code2) : "%esi")

#elif defined(GCC_ASM64X)
#define HAS_CPUID
#define CPUID(code, a, b, c, d) 			\
		ASM_G volatile(					\
			"movq %%rbx, %%rsi   \n\t"		\
			"cpuid               \n\t"		\
			"movl %%ebx, %1      \n\t"		\
			"movq %%rsi, %%rbx   \n\t"		\
			:"=a"(a), "=m"(b), "=c"(c), "=d"(d) 	\
			:"0"(code) : "%rsi")
#define CPUID2(code1, code2, a, b, c, d)		\
		ASM_G volatile(					\
			"movq %%rbx, %%rsi   \n\t"		\
			"cpuid               \n\t"		\
			"movl %%ebx, %1      \n\t"		\
			"movq %%rsi, %%rbx   \n\t"		\
			:"=a"(a), "=m"(b), "=c"(c), "=d"(d) 	\
			:"0"(code1), "2"(code2) : "%rsi")

#elif defined(_MSC_VER)
#include <intrin.h>
#define HAS_CPUID
#define CPUID(code, a, b, c, d)	\
	{	int _z[4]; \
		__cpuid(_z, code); \
		a = _z[0]; \
		b = _z[1]; \
		c = _z[2]; \
		d = _z[3]; \
	}
#define CPUID2(code1, code2, a, b, c, d) \
	{	int _z[4]; \
		__cpuidex(_z, code1, code2); \
		a = _z[0]; \
		b = _z[1]; \
		c = _z[2]; \
		d = _z[3]; \
	}
#endif


// http://msdn.microsoft.com/en-us/library/hskdteyh.aspx
// cpuid.cpp 
// processor: x86, x64
// Use the __cpuid intrinsic to get information about a CPU
// modified for c compliers and use of CPUID macros 
//		- brb, 10/26/10

const char* szFeatures[] =
{
	"x87 FPU On Chip",
	"Virtual-8086 Mode Enhancement",
	"Debugging Extensions",
	"Page Size Extensions",
	"Time Stamp Counter",
	"RDMSR and WRMSR Support",
	"Physical Address Extensions",
	"Machine Check Exception",
	"CMPXCHG8B Instruction",
	"APIC On Chip",
	"Unknown1",
	"SYSENTER and SYSEXIT",
	"Memory Type Range Registers",
	"PTE Global Bit",
	"Machine Check Architecture",
	"Conditional Move/Compare Instruction",
	"Page Attribute Table",
	"36-bit Page Size Extension",
	"Processor Serial Number",
	"CFLUSH Extension",
	"Unknown2",
	"Debug Store",
	"Thermal Monitor and Clock Ctrl",
	"MMX Technology",
	"FXSAVE/FXRSTOR",
	"SSE Extensions",
	"SSE2 Extensions",
	"Self Snoop",
	"Multithreading Technology",
	"Thermal Monitor",
	"Unknown4",
	"Pending Break Enable"
};



/* return cpuid , cache size, etc  */
int extended_cpuid(char* CPUidstr, int* cachelinesize, char* bSSE41Extensions, int do_print)
{
	char CPUString[0x20];
	char CPUBrandString[0x40];
	int CPUInfo[4] = { -1 };
	int nSteppingID = 0;
	int nModel = 0;
	int nFamily = 0;
	int nProcessorType = 0;
	int nExtendedmodel = 0;
	int nExtendedfamily = 0;
	int nBrandIndex = 0;
	int nCLFLUSHcachelinesize = 0;
	int nLogicalProcessors = 0;
	int nAPICPhysicalID = 0;
	int nFeatureInfo = 0;
	int nCacheLineSize = 0;
	int nL2Associativity = 0;
	int nCacheSizeK = 0;
	int nPhysicalAddress = 0;
	int nVirtualAddress = 0;
	int nRet = 0;

	int nCores = 0;
	int nCacheType = 0;
	int nCacheLevel = 0;
	int nMaxThread = 0;
	int nSysLineSize = 0;
	int nPhysicalLinePartitions = 0;
	int nWaysAssociativity = 0;
	int nNumberSets = 0;

	unsigned    nIds, nExIds, i;

	char    bSSE3Instructions = 0;
	char    bMONITOR_MWAIT = 0;
	char    bCPLQualifiedDebugStore = 0;
	char    bVirtualMachineExtensions = 0;
	char    bEnhancedIntelSpeedStepTechnology = 0;
	char    bThermalMonitor2 = 0;
	char    bSupplementalSSE3 = 0;
	char    bL1ContextID = 0;
	char    bCMPXCHG16B = 0;
	char    bxTPRUpdateControl = 0;
	char    bPerfDebugCapabilityMSR = 0;
	//char    bSSE41Extensions = 0;
	char    bSSE42Extensions = 0;
	char    bPOPCNT = 0;

	char    bMultithreading = 0;

	char    bLAHF_SAHFAvailable = 0;
	char    bCmpLegacy = 0;
	char    bSVM = 0;
	char    bExtApicSpace = 0;
	char    bAltMovCr8 = 0;
	char    bLZCNT = 0;
	char    bSSE4A = 0;
	char    bMisalignedSSE = 0;
	char    bPREFETCH = 0;
	char    bSKINITandDEV = 0;
	char    bSYSCALL_SYSRETAvailable = 0;
	char    bExecuteDisableBitAvailable = 0;
	char    bMMXExtensions = 0;
	char    bFFXSR = 0;
	char    b1GBSupport = 0;
	char    bRDTSCP = 0;
	char    b64Available = 0;
	char    b3DNowExt = 0;
	char    b3DNow = 0;
	char    bNestedPaging = 0;
	char    bLBRVisualization = 0;
	char    bFP128 = 0;
	char    bMOVOptimization = 0;

	char    bSelfInit = 0;
	char    bFullyAssociative = 0;

	*bSSE41Extensions = 0;

	// __cpuid with an InfoType argument of 0 returns the number of
	// valid Ids in CPUInfo[0] and the CPU identification string in
	// the other three array elements. The CPU identification string is
	// not in linear order. The code below arranges the information 
	// in a human readable form.
	CPUID(0, CPUInfo[0], CPUInfo[1], CPUInfo[2], CPUInfo[3]);
	//__cpuid(CPUInfo, 0);
	nIds = CPUInfo[0];
	memset(CPUString, 0, sizeof(CPUString));
	*((int*)CPUString) = CPUInfo[1];
	*((int*)(CPUString + 4)) = CPUInfo[3];
	*((int*)(CPUString + 8)) = CPUInfo[2];

	// Get the information associated with each valid Id
	for (i = 0; i <= nIds; ++i)
	{
		CPUID(i, CPUInfo[0], CPUInfo[1], CPUInfo[2], CPUInfo[3]);
		//__cpuid(CPUInfo, i);

		if (do_print > 1)
		{
			printf("\nFor InfoType %d\n", i);
			printf("CPUInfo[0] = 0x%x\n", CPUInfo[0]);
			printf("CPUInfo[1] = 0x%x\n", CPUInfo[1]);
			printf("CPUInfo[2] = 0x%x\n", CPUInfo[2]);
			printf("CPUInfo[3] = 0x%x\n", CPUInfo[3]);
		}

		// Interpret CPU feature information.
		if (i == 1)
		{
			nSteppingID = CPUInfo[0] & 0xf;
			nModel = (CPUInfo[0] >> 4) & 0xf;
			nFamily = (CPUInfo[0] >> 8) & 0xf;
			nProcessorType = (CPUInfo[0] >> 12) & 0x3;
			nExtendedmodel = (CPUInfo[0] >> 16) & 0xf;
			nExtendedfamily = (CPUInfo[0] >> 20) & 0xff;
			nBrandIndex = CPUInfo[1] & 0xff;
			*cachelinesize = nCLFLUSHcachelinesize = ((CPUInfo[1] >> 8) & 0xff) * 8;
			nLogicalProcessors = ((CPUInfo[1] >> 16) & 0xff);
			nAPICPhysicalID = (CPUInfo[1] >> 24) & 0xff;
			bSSE3Instructions = (CPUInfo[2] & 0x1) || 0;
			bMONITOR_MWAIT = (CPUInfo[2] & 0x8) || 0;
			bCPLQualifiedDebugStore = (CPUInfo[2] & 0x10) || 0;
			bVirtualMachineExtensions = (CPUInfo[2] & 0x20) || 0;
			bEnhancedIntelSpeedStepTechnology = (CPUInfo[2] & 0x80) || 0;
			bThermalMonitor2 = (CPUInfo[2] & 0x100) || 0;
			bSupplementalSSE3 = (CPUInfo[2] & 0x200) || 0;
			bL1ContextID = (CPUInfo[2] & 0x300) || 0;
			bCMPXCHG16B = (CPUInfo[2] & 0x2000) || 0;
			bxTPRUpdateControl = (CPUInfo[2] & 0x4000) || 0;
			bPerfDebugCapabilityMSR = (CPUInfo[2] & 0x8000) || 0;
			*bSSE41Extensions = (CPUInfo[2] & 0x80000) || 0;
			bSSE42Extensions = (CPUInfo[2] & 0x100000) || 0;
			bPOPCNT = (CPUInfo[2] & 0x800000) || 0;
			nFeatureInfo = CPUInfo[3];
			bMultithreading = (nFeatureInfo & (1 << 28)) || 0;
		}
	}

	// Calling __cpuid with 0x80000000 as the InfoType argument
	// gets the number of valid extended IDs.
	CPUID(0x80000000, CPUInfo[0], CPUInfo[1], CPUInfo[2], CPUInfo[3]);
	//__cpuid(CPUInfo, 0x80000000);
	nExIds = CPUInfo[0];
	memset(CPUBrandString, 0, sizeof(CPUBrandString));

	// Get the information associated with each extended ID.
	for (i = 0x80000000; i <= nExIds; ++i)
	{
		CPUID(i, CPUInfo[0], CPUInfo[1], CPUInfo[2], CPUInfo[3]);
		//__cpuid(CPUInfo, i);
		if (do_print > 1)
		{
			printf("\nFor InfoType %x\n", i);
			printf("CPUInfo[0] = 0x%x\n", CPUInfo[0]);
			printf("CPUInfo[1] = 0x%x\n", CPUInfo[1]);
			printf("CPUInfo[2] = 0x%x\n", CPUInfo[2]);
			printf("CPUInfo[3] = 0x%x\n", CPUInfo[3]);
		}

		if (i == 0x80000001)
		{
			bLAHF_SAHFAvailable = (CPUInfo[2] & 0x1) || 0;
			bCmpLegacy = (CPUInfo[2] & 0x2) || 0;
			bSVM = (CPUInfo[2] & 0x4) || 0;
			bExtApicSpace = (CPUInfo[2] & 0x8) || 0;
			bAltMovCr8 = (CPUInfo[2] & 0x10) || 0;
			bLZCNT = (CPUInfo[2] & 0x20) || 0;
			bSSE4A = (CPUInfo[2] & 0x40) || 0;
			bMisalignedSSE = (CPUInfo[2] & 0x80) || 0;
			bPREFETCH = (CPUInfo[2] & 0x100) || 0;
			bSKINITandDEV = (CPUInfo[2] & 0x1000) || 0;
			bSYSCALL_SYSRETAvailable = (CPUInfo[3] & 0x800) || 0;
			bExecuteDisableBitAvailable = (CPUInfo[3] & 0x10000) || 0;
			bMMXExtensions = (CPUInfo[3] & 0x40000) || 0;
			bFFXSR = (CPUInfo[3] & 0x200000) || 0;
			b1GBSupport = (CPUInfo[3] & 0x400000) || 0;
			bRDTSCP = (CPUInfo[3] & 0x8000000) || 0;
			b64Available = (CPUInfo[3] & 0x20000000) || 0;
			b3DNowExt = (CPUInfo[3] & 0x40000000) || 0;
			b3DNow = (CPUInfo[3] & 0x80000000) || 0;
		}

		// Interpret CPU brand string and cache information.
		if (i == 0x80000002)
			memcpy(CPUBrandString, CPUInfo, sizeof(CPUInfo));
		else if (i == 0x80000003)
			memcpy(CPUBrandString + 16, CPUInfo, sizeof(CPUInfo));
		else if (i == 0x80000004)
			memcpy(CPUBrandString + 32, CPUInfo, sizeof(CPUInfo));
		else if (i == 0x80000006)
		{
			nCacheLineSize = CPUInfo[2] & 0xff;
			nL2Associativity = (CPUInfo[2] >> 12) & 0xf;
			nCacheSizeK = (CPUInfo[2] >> 16) & 0xffff;
		}
		else if (i == 0x80000008)
		{
			nPhysicalAddress = CPUInfo[0] & 0xff;
			nVirtualAddress = (CPUInfo[0] >> 8) & 0xff;
		}
		else if (i == 0x8000000A)
		{
			bNestedPaging = (CPUInfo[3] & 0x1) || 0;
			bLBRVisualization = (CPUInfo[3] & 0x2) || 0;
		}
		else if (i == 0x8000001A)
		{
			bFP128 = (CPUInfo[0] & 0x1) || 0;
			bMOVOptimization = (CPUInfo[0] & 0x2) || 0;
		}
	}
	strcpy_s(CPUidstr, 64, CPUBrandString);
	// Display all the information in user-friendly format.
	if (do_print > 0)
		printf("\n\nCPU String: %s\n", CPUString);

	if (nIds >= 1)
	{
		if (do_print > 0)
		{
			if (nSteppingID)
				printf("Stepping ID = %d\n", nSteppingID);
			if (nModel)
				printf("Model = %d\n", nModel);
			if (nFamily)
				printf("Family = %d\n", nFamily);
			if (nProcessorType)
				printf("Processor Type = %d\n", nProcessorType);
			if (nExtendedmodel)
				printf("Extended model = %d\n", nExtendedmodel);
			if (nExtendedfamily)
				printf("Extended family = %d\n", nExtendedfamily);
			if (nBrandIndex)
				printf("Brand Index = %d\n", nBrandIndex);
			if (nCLFLUSHcachelinesize)
				printf("CLFLUSH cache line size = %d\n",
					nCLFLUSHcachelinesize);
			if (bMultithreading && (nLogicalProcessors > 0))
				printf("Logical Processor Count = %d\n", nLogicalProcessors);
			if (nAPICPhysicalID)
				printf("APIC Physical ID = %d\n", nAPICPhysicalID);

			if (nFeatureInfo || bSSE3Instructions ||
				bMONITOR_MWAIT || bCPLQualifiedDebugStore ||
				bVirtualMachineExtensions || bEnhancedIntelSpeedStepTechnology ||
				bThermalMonitor2 || bSupplementalSSE3 || bL1ContextID ||
				bCMPXCHG16B || bxTPRUpdateControl || bPerfDebugCapabilityMSR ||
				*bSSE41Extensions || bSSE42Extensions || bPOPCNT ||
				bLAHF_SAHFAvailable || bCmpLegacy || bSVM ||
				bExtApicSpace || bAltMovCr8 ||
				bLZCNT || bSSE4A || bMisalignedSSE ||
				bPREFETCH || bSKINITandDEV || bSYSCALL_SYSRETAvailable ||
				bExecuteDisableBitAvailable || bMMXExtensions || bFFXSR || b1GBSupport ||
				bRDTSCP || b64Available || b3DNowExt || b3DNow || bNestedPaging ||
				bLBRVisualization || bFP128 || bMOVOptimization)
			{
				printf("\nThe following features are supported:\n");

				if (bSSE3Instructions)
					printf("\tSSE3\n");
				if (bMONITOR_MWAIT)
					printf("\tMONITOR/MWAIT\n");
				if (bCPLQualifiedDebugStore)
					printf("\tCPL Qualified Debug Store\n");
				if (bVirtualMachineExtensions)
					printf("\tVirtual Machine Extensions\n");
				if (bEnhancedIntelSpeedStepTechnology)
					printf("\tEnhanced Intel SpeedStep Technology\n");
				if (bThermalMonitor2)
					printf("\tThermal Monitor 2\n");
				if (bSupplementalSSE3)
					printf("\tSupplemental Streaming SIMD Extensions 3\n");
				if (bL1ContextID)
					printf("\tL1 Context ID\n");
				if (bCMPXCHG16B)
					printf("\tCMPXCHG16B Instruction\n");
				if (bxTPRUpdateControl)
					printf("\txTPR Update Control\n");
				if (bPerfDebugCapabilityMSR)
					printf("\tPerf\\Debug Capability MSR\n");
				if (*bSSE41Extensions)
					printf("\tSSE4.1 Extensions\n");
				if (bSSE42Extensions)
					printf("\tSSE4.2 Extensions\n");
				if (bPOPCNT)
					printf("\tPPOPCNT Instruction\n");

				i = 0;
				nIds = 1;
				while (i < (sizeof(szFeatures) / sizeof(const char*)))
				{
					if (nFeatureInfo & nIds)
					{
						printf("\t");
						printf("%s", szFeatures[i]);
						printf("\n");
					}

					nIds <<= 1;
					++i;
				}
				if (bLAHF_SAHFAvailable)
					printf("\tLAHF/SAHF in 64-bit mode\n");
				if (bCmpLegacy)
					printf("\tCore multi-processing legacy mode\n");
				if (bSVM)
					printf("\tSecure Virtual Machine\n");
				if (bExtApicSpace)
					printf("\tExtended APIC Register Space\n");
				if (bAltMovCr8)
					printf("\tAltMovCr8\n");
				if (bLZCNT)
					printf("\tLZCNT instruction\n");
				if (bSSE4A)
					printf("\tSSE4A (EXTRQ, INSERTQ, MOVNTSD, MOVNTSS)\n");
				if (bMisalignedSSE)
					printf("\tMisaligned SSE mode\n");
				if (bPREFETCH)
					printf("\tPREFETCH and PREFETCHW Instructions\n");
				if (bSKINITandDEV)
					printf("\tSKINIT and DEV support\n");
				if (bSYSCALL_SYSRETAvailable)
					printf("\tSYSCALL/SYSRET in 64-bit mode\n");
				if (bExecuteDisableBitAvailable)
					printf("\tExecute Disable Bit\n");
				if (bMMXExtensions)
					printf("\tExtensions to MMX Instructions\n");
				if (bFFXSR)
					printf("\tFFXSR\n");
				if (b1GBSupport)
					printf("\t1GB page support\n");
				if (bRDTSCP)
					printf("\tRDTSCP instruction\n");
				if (b64Available)
					printf("\t64 bit Technology\n");
				if (b3DNowExt)
					printf("\t3Dnow Ext\n");
				if (b3DNow)
					printf("\t3Dnow! instructions\n");
				if (bNestedPaging)
					printf("\tNested Paging\n");
				if (bLBRVisualization)
					printf("\tLBR Visualization\n");
				if (bFP128)
					printf("\tFP128 optimization\n");
				if (bMOVOptimization)
					printf("\tMOVU Optimization\n");
			}
		}
	}

	if (nExIds >= 0x80000004 && do_print > 0)
		printf("\nCPU Brand String: %s\n", CPUBrandString);

	if (nExIds >= 0x80000006 && do_print >0)
	{
		printf("Cache Line Size = %d\n", nCacheLineSize);
		printf("L2 Associativity = %d\n", nL2Associativity);
		printf("Cache Size = %dK\n", nCacheSizeK);
	}


	for (i = 0;; i++)
	{
		CPUID2(0x4, i, CPUInfo[0], CPUInfo[1], CPUInfo[2], CPUInfo[3]);
		//__cpuidex(CPUInfo, 0x4, i);
		if (!(CPUInfo[0] & 0xf0)) break;

		if (i == 0)
		{
			nCores = CPUInfo[0] >> 26;
			if (do_print > 0)
				printf("\n\nNumber of Cores = %d\n", nCores + 1);
		}

		nCacheType = (CPUInfo[0] & 0x1f);
		nCacheLevel = (CPUInfo[0] & 0xe0) >> 5;
		bSelfInit = (CPUInfo[0] & 0x100) >> 8;
		bFullyAssociative = (CPUInfo[0] & 0x200) >> 9;
		nMaxThread = (CPUInfo[0] & 0x03ffc000) >> 14;
		nSysLineSize = (CPUInfo[1] & 0x0fff);
		nPhysicalLinePartitions = (CPUInfo[1] & 0x03ff000) >> 12;
		nWaysAssociativity = (CPUInfo[1]) >> 22;
		nNumberSets = CPUInfo[2];

		if (do_print > 0)
		{
			printf("\n");

			printf("ECX Index %d\n", i);
			switch (nCacheType)
			{
			case 0:
				printf("   Type: Null\n");
				break;
			case 1:
				printf("   Type: Data Cache\n");
				break;
			case 2:
				printf("   Type: Instruction Cache\n");
				break;
			case 3:
				printf("   Type: Unified Cache\n");
				break;
			default:
				printf("   Type: Unknown\n");
			}

			printf("   Level = %d\n", nCacheLevel + 1);
			if (bSelfInit)
			{
				printf("   Self Initializing\n");
			}
			else
			{
				printf("   Not Self Initializing\n");
			}
			if (bFullyAssociative)
			{
				printf("   Is Fully Associatve\n");
			}
			else
			{
				printf("   Is Not Fully Associatve\n");
			}
			printf("   Max Threads = %d\n",
				nMaxThread + 1);
			printf("   System Line Size = %d\n",
				nSysLineSize + 1);
			printf("   Physical Line Partions = %d\n",
				nPhysicalLinePartitions + 1);
			printf("   Ways of Associativity = %d\n",
				nWaysAssociativity + 1);
			printf("   Number of Sets = %d\n",
				nNumberSets + 1);
		}
	}

	return  nRet;
}

// function containing system commands to get the computer name, CPU speed, etc
static void get_computer_info(char* idstr)
{
	//int ret;

	//figure out cpu freq. 0.1 seconds won't be very accurate 
	MEAS_CPU_FREQUENCY = measure_processor_speed() / 1.0e5;

#ifdef __APPLE__
	// something in extended cpuid causes a segfault on mac builds.
	// just disable it for now - this information is not critical for
	// program operation.
	strcpy(idstr, "N/A");
	CLSIZE = 0;
	L1CACHE = DEFAULT_L1_CACHE_SIZE;
	L2CACHE = DEFAULT_L2_CACHE_SIZE;
	HAS_SSE41 = 0;

#else
	
	// run an extended cpuid command to get the cache line size, and
	// optionally print a bunch of info to the screen

	extended_cpuid(idstr, &CLSIZE, &HAS_SSE41, VERBOSE_PROC_INFO);

#if defined(WIN32)

	sysname_sz = MAX_COMPUTERNAME_LENGTH + 1;
	GetComputerNameA(sysname, &sysname_sz);

#else

	ret = gethostname(sysname, sizeof(sysname) / sizeof(*sysname));
	sysname[(sizeof(sysname) - 1) / sizeof(*sysname)] = 0;	// null terminate
	if (ret != 0)
	{
		printf("error occured when getting host name\n");
		strcpy(sysname, "N/A");
	}
	sysname_sz = strlen(sysname);

#endif

#endif
	return;
}



/* evaluate 1 or more expressions, separated by commas */
static retCode ComputeMultiExpr(std::string expr, Znum result) {
	std::string subExpr;
	retCode rv = retCode::EXPR_OK;
	size_t subStart = 0, subEnd;
	int bc = 0;   /* bracket count */
	int asgCt = 0;   /* number of assignment operators */
	while (subStart < expr.size()) {
		bc = 0;
		/* find  separator or end of text*/
		for (subEnd = subStart + 1; subEnd < expr.size(); subEnd++) {
			if (expr[subEnd] == '(')
				bc++;
			if (expr[subEnd] == ')')
				bc--;
			if (bc == 0 && expr[subEnd] == ',')
				break;
		}
		if (bc != 0)
			return retCode::PAREN_MISMATCH;
		subExpr = expr.substr(subStart, subEnd - subStart);
		removeInitTrail(subExpr);  /* remove initial & trailing blanks */
		rv = ComputeExpr(subExpr, result, asgCt);
		if (rv != retCode::EXPR_OK) {
			textError(rv);   // invalid expression; print error message
			return rv;
		}
		else {
			if (asgCt == 0)
				std::cout << " = ";
			else {
				/* print names of variables assigned values */
				for (size_t ix = 0; ix < subExpr.size() && asgCt > 0; ix++) {
					putchar(subExpr[ix]);
					if (subExpr[ix] == '=')
						asgCt--;
				}
			}

			ShowLargeNumber(result, 6, true, hex);   // print value of expression
			std::cout << '\n';

		}
		subStart = subEnd + 1;  /* move past , */
	} /* end of while loop */

	return rv;
}

/* format is IF (expression) REPEAT 
          or IF (expression) STOP 
		  or IF (expression) THEN (expression) [ELSE (expression)]
          the [ELSE (expression)] part is optional
return -1 if syntax is invalid
 return 0 if expression is 0 for STOP or REPEAT
 return 1 if expression NE 0 and action is STOP
 return 2 if expression NE 0 and action is REPEAT 
 return 3 if THEN or ELSE expression evaluated succesfully */
static int ifCommand(const std::string &command) {
	int ixx, ixx2, exprLen;
	int asgCt = 0;  /* number of assignment operators*/
	std::string expr;
	Znum result;
	int bc = 1;  /* bracket count */

	/* find ( after IF */
	for (ixx = 2; ixx < command.size(); ixx++)
		if (command[ixx] == '(')
			break;
	if (ixx >= command.size()) {
		std::cout << "No ( after IF command \n";
		return -1;
	}

	/* find matching ) after IF */
	for (ixx2 = ixx + 1; ixx2 < command.size(); ixx2++) {
		if (command[ixx2] == '(')
			bc++;
		if (command[ixx2] == ')') {
			bc--;
			if (bc == 0)
				break;   /* found matching closing bracket */
		}
	}
	if (ixx2 >= command.size()) {
		std::cout << "No  matching ) after IF command \n";
		return -1;
	}

	/* evaluate expression betwen brackets */
	exprLen = ixx2 - ixx - 1;
	expr = command.substr(ixx + 1, exprLen);
	retCode rv = ComputeExpr(expr, result, asgCt);
	if (rv != retCode::EXPR_OK) {
		textError(rv);   // invalid expression; print error message
		return -1;
	}

	/* move ixx2 to next non-blank character */
	for (ixx2++; ixx2 < command.size() && isblank(command[ixx2]); ixx2++);
	if (command.substr(ixx2) == "STOP") {
		if (result != 0)
			return 1;
		else
			return 0;
	}
	if (command.substr(ixx2) == "REPEAT") {
		if (result != 0)
			return 2;
		else
			return 0;
	}
	if (command.substr(ixx2, 4) == "THEN") {
		size_t ex1Start=ixx2+4, ex1End, ex2Start, ex2End;

		/* move ex1Start to next non-blank character */
		for (ex1Start=ixx2+4; ex1Start < command.size() && isblank(command[ex1Start]); 
			ex1Start++);
	
		if (ex1Start >= command.size() || command[ex1Start] != '(') {
			std::cout << " THEN is not followed by ( \n";
			return -1;  /* no ( after THEN so invalid */
		}
		ex1End = ex1Start + 1;
		bc = 1;

		/* find matching ) */
		while (ex1End < command.size()) {
			if (command[ex1End] == '(')
				bc++;
			if (command[ex1End] == ')')
				bc--;
			if (bc == 0)
				break;
			ex1End++;
		}
		if (bc != 0) {
			std::cout << "Matching ) not found \n";
			return -1;   /* matching ) not found */
		}

		/* move ex2Start to next non-blank character if any */
		for (ex2Start = ex1End + 1; ex2Start < command.size() && isblank(command[ex2Start]);
			ex2Start++);
		if (ex2Start < command.size()) {
			/* still some unprocessed characters */

			if (command.size() < ex2Start + 4 || command.substr(ex2Start, 4) != "ELSE") {
				std::cout << "Format invalid. ELSE not found \n";
				return -1;
			}
			ex2Start += 4;   /* move past ELSE */
			for (; ex2Start < command.size() && isblank(command[ex2Start]); ex2Start++);
			if (ex2Start >= command.size() || command[ex2Start] != '(') {
				std::cout << "ELSE not followed by ( \n";
				return -1;  /* no ( after ELSE so invalid */
			}
			bc = 1;
			ex2End = ex2Start + 1;
			while (ex2End < command.size()) {
				if (command[ex2End] == '(')
					bc++;
				if (command[ex2End] == ')')
					bc--;
				if (bc == 0)
					break;
				ex2End++;
			}

			if (bc != 0) {
				std::cout << "Matching ) not found \n";
				return -1;   /* matching ) not found */
			}
		}
		if (result != 0) {
			expr = command.substr(ex1Start+1, ex1End - ex1Start-1);
		}
		else if (ex2Start < command.size())
			expr = command.substr(ex2Start+1, ex2End - ex2Start-1);
		else 
			return 3;  /* no ELSE expression to evaluate */

		retCode rv = ComputeMultiExpr(expr, result);
		return 3;  /* IF (...) THEN (...) ELSE (...) processed OK*/
	}

	std::cout << "Neither STOP, REPEAT, nor THEN found \n";
	return -1;  /* neither STOP nor REPEAT nor THEN found */
}

/* check for commands. return 2 for exit, 1 for other command, 0 if not a valid command*/
static int processCmd(const std::string &command) {

	/* list of commands (static for efficiency) */
	const static std::vector<std::string> list =
	{ "EXIT", "SALIDA", "HELP", "AYUDA", "E",     "S",  
	   "F ",   "X",     "D",    "TEST",  "MSIEVE", "YAFU", 
	   "V ",   "PRINT", "LIST", "LOOP", "REPEAT",  "IF" };

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
	case 4: /* E */ {
		 lang = 0; return 1; }           // english
	case 5: /* S */ { 
		lang = 1; return 1; }	          // spanish (Español)
	case 6: /* F */ { 
		/* will not throw an exception if input has fat finger syndrome.
			If no valid digits found, sets factorFlag to 0 */
		factorFlag = atoi(command.substr(1).data());
		std::cout << "factor set to " << factorFlag << '\n';
		return 1; }  
	case 7: /* X */ { 
		hex = true; return 1; }         // hexadecimal output
	case 8: /* D */ { 
		hex = false; return 1; }        // decimal output
	case 9: /* TEST */ {
			/* there are a number of commands that begin with TEST */
			if (command == "TEST") {
				doTests();         // do basic tests 
				return 1;
			}

			int ttype = command.data()[4] ;  /* for TEST2 ttype = 2, etc*/
			switch (ttype) {
			case '2': {
				if (command.size() > 5)
					doTests2(command.substr(6));         // do basic tests 
				else
					doTests2("");
				return 1;
			}
			case '3':
	#ifdef BIGNBR
				{
					doTests3();         // do basic tests 
					return 1;
				}
	#else
				return 0;  /* Bignbr tests omitted */
	#endif
			case '4': {
					doTests4();         // factorise Mersenne numbers 
					return 1;
				}
			case '5': {
					doTests6(false);         // do YAFU tests 
					return 1;
				}
			case '6': {
					doTests6(true);         // do Msieve tests 
					return 1;
				}
			case '7': {
					/* Lucas-Lehmer test */
					if (command.size() > 5)
						doTests7(command.substr(6));
					else
						doTests7("");
					return 1;
				}
			case '8': {
				testerrors();
				return 1;
			}
			case '9': {
				doTests9();
				return 1;
			}
			default:
				return 0;  /* not a recognised command */
			}
		}
	case 10: /* MSIEVE */ {
			msieveParam(command);
			return 1;
		}
	case 11: /* YAFU */ {
			yafuParam(command);
			return 1;
		}
	case 12: /* V */ {
			/* will not throw an exception if input has fat finger syndrome.
			If no valid digits found, sets verbose to 0 */
			verbose = atoi(command.substr(1).data());
			std::cout << "verbose set to " << verbose << '\n';
			return 1;
		}
	case 13: /* PRINT */ {
		if (command.size() > 4)
			printvars(command.substr(5));
		else
			printvars("");
		return 1;
	}
	case 14: /* LIST */ {
		if (exprList.empty())
			std::cout << "No expressions stored yet \n";
		else
			std::cout << exprList.size() << " expressions stored \n";
		for (auto exp : exprList) {
			std::cout << exp << '\n';
		}
		return 1;
	}
	case 15: /* LOOP */ {
		exprList.clear();
		return 1;
	}
	case 16: /* REPEAT */ {
		int repeat = 0;
		Znum result;
		retCode rv;
		if (command.size() > 6)
			repeat = atoi(command.substr(6).data());
		else
			repeat = 1;  /* repeat once if no count specified */

		loop1:
		for (int count = 1; count <= repeat; count++) {
			for (auto expr : exprList) {
				if (expr.size() > 2 && expr.substr(0, 2) == "IF") {
					int rv2 = ifCommand(expr);  /* analyse stored command again */
					if (rv2 == 2)
						goto loop1;  /* repeat all stored expressions again */
					else if (rv2 == 0)
						continue;  /* expr = 0; do not perform STOP or LOOP */
					else if (rv2 == 1)
						break;     /* expression is not 0 and STOP specified */
					else if (rv2 == 3)
						continue;   /* THEN or ELSE expression has been evaluated */
					else
						throw std::logic_error("unknown return code");
					abort();  /* where are we? */
				}
				/* recalculate each stored expression */
				rv = ComputeMultiExpr(expr, result);
			}
		}
		return 1;
	}
	case 17: /* IF */ {
		/*format is IF (expression) REPEAT or IF (expression) STOP 
		or IF (expression) THEN (expression[, expression ...])
		                     ELSE (expression[, expression ...]) */
		Znum result;
		int rv = ifCommand(command);  /* analyse IF command */
		if (rv == 2) {
			/* REPEAT */
			exprList.push_back(command);
			loop:
			for (auto expr : exprList) {
				if (expr.size() > 2 && expr.substr(0, 2) == "IF") {
					int rv2 = ifCommand(expr);  /* analyse stored command again */
					if (rv2 == 2)
						goto loop;  /* repeat all stored expressions again */
					else if (rv2 == 0)
						continue;  /* expr = 0; do not perform STOP or LOOP */
					else if (rv2 == 1)
						break;     /* expression is not 0 and STOP specified */
					else if (rv2 == 3)
						continue;   /* THEN or ELSE expression has been evaluated */
					else
						throw std::logic_error ("unknown return code");
						abort();  /* where are we? */
				}

				/* recalculate each stored expression */
				retCode rv = ComputeMultiExpr(expr, result);
			}
			return 1;  /* loop completed */
		}
		if (rv == 0 || rv == 3)
			/* if command processed, now save it for possible loop */
			exprList.push_back(command);  
		/* to be completed for rv =-1, rv = 1*/
		return 1;
	}
	default:
		return 0;   /* not a recognised command */
	}
}

/* initialisation code, executed once at start of program execution */
static void initialise(int argc, char *argv[]) {
	flags f = { 0,0,0, 0,0,0, 0,0,0, 0,0,0 };
#ifndef BIGNBR
	unsigned int control_word; /* save current state of fp control word here */
#endif
	int version[4]; /* version info from .exe file (taken from .rc resource file) */
	std::string modified;  /* date & time program last modified */

	f.UnEx = 1;        /* set unhandled exception 'filter' */
	f.abort = 1;       /* trap abort (as a signal) */
	f.sigterm = 1;     /* trap terminate signal */
	f.sigill = 1;      /* trap 'illegal' signal */
	f.interrupt = 1;   /* trap interrupt signal */
	//f.sigsegv = 1;     /* trap segment violation signal */
#ifdef _DEBUG
	/* only seems to work properly if compiled in debug mode */
	f.InvParam = 1;    /* trap invalid parameters on library calls */
#endif
	//f.sigfpe = 1;      /* trap floating point error signal */
	SetProcessExceptionHandlers(f);

	/* if we trap floating point errors we trap  _EM_INVALID in mpir prime test
        functions that actually work OK */
#ifndef BIGNBR
	err = _controlfp_s(&control_word, _EM_INEXACT | _EM_UNDERFLOW, MCW_EM);
	/* trap hardware FP exceptions except inexact and underflow which are
	considered to be normal, not errors. */
	if (err) {
		printf_s("could not set FP control word\n");
		exit (-1);
	}
#endif

	hConsole = GetStdHandle(STD_OUTPUT_HANDLE);  // get handle for stdout
	handConsole = GetConsoleWindow();            // get handle for console window
	if (hConsole == INVALID_HANDLE_VALUE)
	{
		fprintf_s(stderr, "GetStdHandle failed with %d at line %d\n", GetLastError(), __LINE__);
		Beep(750, 1000);
		exit (EXIT_FAILURE);
	}

	VersionInfo(argv[0], version, modified); /* get version info from .exe file */
	printf_s("%s Bigint calculator Version %d.%d.%d.%d \n", myTime(),
		version[0], version[1], version[2], version[3]);
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
	std::cout << ver / 100000 << '.';    // next 2 digits
	ver %= 100000;                        // remove next 2 digits
	std::cout << ver << '\n';             // last 5 digits

	auto lc = setlocale(LC_ALL, "en-EN");      // allows non-ascii characters to print
#endif

	printf_s("locale is now: %s\n", setlocale(LC_ALL, NULL));
	std::cout << "GMP version: " << __GNU_MP_VERSION << '.' << __GNU_MP_VERSION_MINOR
		<< '.' << __GNU_MP_VERSION_PATCHLEVEL << '\n';

#ifdef __MPIR_VERSION
	std::cout << "MPIR version: " << __MPIR_VERSION << '.' << __MPIR_VERSION_MINOR
		<< '.' << __MPIR_VERSION_PATCHLEVEL << '\n';
#endif
	std::cout << "Boost version: " << BOOST_VERSION << '\n';

	processIni(argv[0]); // read .ini file if it exists

	//get the computer name, cache sizes, etc.  store in globals
	get_computer_info(CPU_ID_STR);

	printf_s("detected %s\n", CPU_ID_STR);
	printf_s("measured cpu frequency ~= %f\n", 	MEAS_CPU_FREQUENCY);

	return;
}

/* get input from stdin. Any continuation lines are appended to 1st line.
Initial & trailing spaces are removed. Letters are converted to upper case. 
ctrl-c or ctrl-break will force the function to return, with or without input, 
but only after a 10 sec delay.*/
static void myGetline(std::string &expr) {
retry:
	getline(std::cin, expr);    // expression may include spaces
	if (breakSignal) {
		Sleep(10000);   /* wait 10 seconds */
		return;     /* Program interrupted: ctrl-c or ctrl-break */
	}
	strToUpper(expr, expr);		// convert to UPPER CASE 
	removeInitTrail(expr);       // remove initial & trailing spaces

	if (expr.empty()) {
		Beep(3000, 250);     /* audible alarm instead of error msg */
		goto retry;          /* input is zero-length; go back*/
	}

	/* check for continuation character. If found get continuation line(s) */
	while (expr.back() == '\\') {   /* ends with continuation character? */
		std::string cont;
		std::cout << "continue: ";
		getline(std::cin, cont);   /* get continuation line */
		strToUpper(cont, cont);   // convert to UPPER CASE 
		while (!cont.empty() && isspace(cont.back())) {
			cont.resize(cont.size() - 1);   /* remove trailing space */
		}
		expr.resize(expr.size() - 1); /* remove trailing \ char */
		expr += cont;    /* append continuation line to previous input */
	}

	if (breakSignal) {
		Sleep(10000);   /* wait 10 seconds */
		return;     /* Program interrupted: ctrl-c or ctrl-break */
	}

}

int main(int argc, char *argv[]) {
	std::string expr;
	Znum Result;
	retCode rv;
	int asgCt;  /* number of assignment operators */
	bool multiV = false;

	try {

		initialise(argc, argv);  /* initialisation code only executed once */

		/* start of main loop. Normal exit is via EXIT command */
		while (true) {
			if (lang == 0) {
				printf_s("enter expression to be processed, or HELP, or EXIT\n");
			}
			else
				printf_s("ingrese la expresión para ser procesada, o AYUDA o SALIDA\n");

			PlaySoundA(attsound.c_str(), NULL,
				SND_FILENAME | SND_NODEFAULT | SND_ASYNC | SND_NOSTOP);

			myGetline(expr);  /* get input from stdin */
			if (breakSignal)
				break;    /* Program interrupted: ctrl-c or ctrl-break */

			// prevent the sleep idle time-out.
			SetThreadExecutionState(ES_CONTINUOUS | ES_SYSTEM_REQUIRED);

			int cmdCode = processCmd(expr);  /* is input a command? */
			if (breakSignal)
				break;    /* Program interrupted: ctrl-c or ctrl-break */
			if (cmdCode == 2) 
				break;    // EXIT command
			if (cmdCode == 1) {
				// command has been fully processed; go back to start of loop
				// Clear EXECUTION_STATE flags to allow the system to idle to sleep normally.
				SetThreadExecutionState(ES_CONTINUOUS);
				continue;
			}
			if (cmdCode != 0) {
				fprintf_s(stderr, "Invalid return code %x from processCmd \n", cmdCode);
				continue;    // weird return code; can't process further; 
				             //  go back to start of loop
			}

			/* input is not a valid command; assume it is an expression */
			auto start = clock();	// used to measure execution time

			rv = ComputeExpr(expr, Result, asgCt, &multiV); /* analyse expression, compute value*/

			if (rv != retCode::EXPR_OK) {
				textError(rv);   // invalid expression; print error message
			}
			else {
				exprList.push_back(expr);  /* save text of expression */
				if (asgCt == 0 && !multiV) {
					std::cout << " = ";
					ShowLargeNumber(Result, 6, true, hex);   // print value of expression
					std::cout << '\n';
					if (factorFlag > 0) {
						doFactors(Result, false); /* factorise Result, calculate number of divisors etc */
						results.clear();  // get rid of unwanted results
					}
				}
				else {
					if (multiV) {
						/* expression returned multiple values */
						std::cout << " = ";
						if (roots.size() <= 31 ) 
							/* print all results if <= 31 values */
							for (auto r : roots) {
								std::cout << r << ", ";
							}
						else { /* print 1st 20 and last 10 results */
							for (int i = 0; i < 20; i++) {
								std::cout << roots[i] << ", " ;
							}
							std::cout << "\n ... \n";
							for (size_t i = roots.size()-10; i < roots.size(); i++) {
								std::cout << roots[i] << ", ";
							}
						}
						putchar('\n');
						std::cout << "found " << roots.size() << " results \n";
					}
					if (asgCt != 0)
					printvars(""); /* print variables names & values */
				}
			}

			auto end = clock();   // measure amount of time used
			double elapsed = (double)end - start;
			PrintTimeUsed(elapsed, "time used = ");
			// Clear EXECUTION_STATE flags to allow the system to idle to sleep normally.
			SetThreadExecutionState(ES_CONTINUOUS);

			if (breakSignal)
				break;     /* Program interrupted */

			/* now go back to start of loop */
		} /* end of while loop */

		/* get to here when we break out of main loop, usually by EXIT command */

		// Clear EXECUTION_STATE flags to allow the system to idle to sleep normally.
		SetThreadExecutionState(ES_CONTINUOUS);

		return EXIT_SUCCESS;  // EXIT command entered
	}

#undef max  /* remove max defined in windows.h  because of name clash */

	/* code below catches C++ 'throw' type exceptions */
	catch (const std::exception& e) {

		// Clear EXECUTION_STATE flags to allow the system to idle to sleep normally.
		SetThreadExecutionState(ES_CONTINUOUS);

		printf_s("\n*** standard exception caught, message '%s'\n", e.what());
		Beep(1200, 1000);              // sound at 1200 Hz for 1 second
		std::cout << "Press ENTER to continue...";
		std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		exit(EXIT_FAILURE);
	}

	catch (const char *str) {

		// Clear EXECUTION_STATE flags to allow the system to idle to sleep normally.
		SetThreadExecutionState(ES_CONTINUOUS);

		printf_s("\n*** Caught exception: <%s> \n", str);
		Beep(1200, 1000);              // sound at 1200 Hz for 1 second
		std::cout << "Press ENTER to continue...";
		std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		exit(EXIT_FAILURE);
	}

	catch (int e) {

		// Clear EXECUTION_STATE flags to allow the system to idle to sleep normally.
		SetThreadExecutionState(ES_CONTINUOUS);

		printf_s("\n*** Caught exception: <%d> \n", e);
		Beep(1200, 1000);              // sound at 1200 Hz for 1 second
		std::cout << "Press ENTER to continue...";
		std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		exit(EXIT_FAILURE);

	}

	catch (...) {
		// this executes if f() throws any other unrelated type
		// This catch block probably only would be executed under /EHa compiler option 
		/* most likely to be a SEH-type exception */

		// Clear EXECUTION_STATE flags to allow the system to idle to sleep normally.
		SetThreadExecutionState(ES_CONTINUOUS);

		printf_s("\n*** unknown exception ocurred\n");
		Beep(1200, 1000);              // sound at 1200 Hz for 1 second
		std::cout << "Press ENTER to continue...";
		std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		exit(EXIT_FAILURE);
	}
}
