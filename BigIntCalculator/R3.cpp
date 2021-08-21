#include "pch.h"

extern int verbose;

/* _AMD64_ 1 definition neeeded if windows.h is not included */
#ifdef _M_AMD64
#ifndef _AMD64_
#define _AMD64_ 1
#endif
#endif
#define ENABLE_INTSAFE_SIGNED_FUNCTIONS 1
#include <intsafe.h>
#include <intrin.h>

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

struct factorsS {
	int factorcount;           // number of unique prime factors
	__int64 factorlist[19][2]; // prime factor, exponent
};

extern bool *primeFlags;
extern unsigned long long int primeListMax;
extern unsigned long long *primeList;
extern unsigned int prime_list_count;
unsigned long long int gcd(unsigned long long int u, unsigned long long int v);
void generatePrimes(unsigned long long int max_val);

// calculate x^n
constexpr __int64 power(const __int64 x, unsigned int n) {
	__int64 p = x, pnew = 0;  // p is a power of x
	__int64 r = 1;  // value of x^n
	HRESULT ret = 0;

	p = x;
	r = 1;
	while (n > 0) {
		if (n % 2 == 1) {
			//r *= p;
			ret = LongLongMult(r, p, &r);
			if (ret != S_OK) {
				ThrowExs("integer overflow: ")
			}
		}
		n /= 2;
		if (n == 0) break;
		//p *= p;
		ret = LongLongMult(p, p, &pnew);
		if (ret != S_OK) {
			ThrowExs("integer overflow: ")
		};
		p = pnew;
	}
	return(r);
}

// calculate the product of all the factors in the list
unsigned __int64 FactorMult(const factorsS f) {
	unsigned __int64 prod = 1;
	HRESULT res;
	for (int i = 0; i < f.factorcount; i++) {
		/*  prod *= p^e */
		for (int k = 1; k <= f.factorlist[i][1]; k++) {
			//prod *= num.factorlist[i][0];
			res = ULongLongMult(prod, f.factorlist[i][0], &prod);
			if (res != S_OK) {
				ThrowExc("integer overflow: ")
			}
		}
	}
	return prod;
}

/* see https://en.wikipedia.org/wiki/Jacobi_symbol */
int jacobi(__int64 k, unsigned __int64 n) {
	bool neg = false;

	if (n <= 0 || (n & 1) == 0) {
		ThrowExc("invalid parameter: must be odd +ve number")
	}

	if (n == 1)
		return 1;

	while (true) {  /* loop until k = 1 or 0 */
		/* step 1*/
		k %= (long long)n;
		if (k < 0) {
			k += n;  // if k is -ve, make +ve
		}

		if (k == 0) return 0;  /* wikipedia puts this in step 3, but it's
							   needed here to avoid an infinite loop in step 2 */

							   /* step 2*/
		if (n % 8 == 1 || n % 8 == 7)
			while ((k & 1) == 0) {
				k /= 2;
			}
		else {
			while ((k & 1) == 0) {
				k /= 2;
				neg = !neg;        /* flip sign */
			}
		}

		/* step 3*/
		if (k == 0) return 0;
		if (k == 1)
			if (neg)  return -1;
			else return 1;

		/* step 4; a and b are both odd, +ve and are coprime, k < n */
		if (n % 4 == 3 && k % 4 == 3)
			neg = !neg;   /* flip sign */
		auto temp = n;
		n = k;           /* swap n & k */
		k = temp;
	}

}

// fast way to check for primes, but generatePrimes must have been called first.
// true indicates a prime number
bool isPrime2(unsigned __int64 num) {

	if (num == 1) return false;
	if (num == 2) return true;
	if ((num & 1) == 0) return false; // even numbers other than 2 are not prime
	return !primeFlags[num/2];

	//unsigned __int64 index;
	//int bit;

	////index = num / 64;
	//index = num >> 6;
	////bit = (num % 64) / 2;
	//bit = (num & 0x3f) >> 1;   // get last 6 bits, then divide by 2
	//return (((primeFlags[index] >> bit) & 1) == 0);
}


/*
* calculates (a * b) % c taking into account that a * b might overflow
*/
unsigned __int64 modMult(unsigned __int64 a, unsigned __int64 b, unsigned __int64 mod)
{
	unsigned __int64 x = 0, retval;
	HRESULT  res;

	static bool firsttime = true;
	static mpz_t al, res_t, rem;

	res = ULongLongMult(a, b, &retval);   // retval = a*b
	if (res == S_OK) {
		x = retval % mod;		// if a*b didn't overflow we do calc the short way
	}
	else {
		if (firsttime) {		// is this the 1st time modMult was called?
			mpz_inits(al, res_t, rem, NULL);		// if so, allocate storage for bigints
			firsttime = false;
		}
		mpz_set_ui(al, a);					// convert a to bigint
		mpz_mul_ui(res_t, al, b);          // res_t = a*b
		x = mpz_fdiv_r_ui(rem, res_t, mod);  //  x = res_t % mod
	}
	return x;
}

// calculate x^n%mod
static __int64 mPowerInt(const unsigned int x, unsigned int n, const unsigned int mod)
{
	__int64 p, pnew;  // p is a power of x
	__int64 r;  // value of x^n
	HRESULT ret;

	p = x;
	r = 1;

	/* loop log2(n) times.*/
	while (n > 0) {
		if (n % 2 == 1) {
			//r *= p;
			ret = LongLongMult(r, p, &r);
			if (ret != S_OK) {
				ThrowExc("integer overflow: ")
			}
			r %= mod;
		}
		n /= 2;
		if (n == 0) break;
		//p *= p   i.e. p = p^2;
		ret = LongLongMult(p, p, &pnew);
		if (ret != S_OK) {
			ThrowExc("integer overflow: ")
		}
		p = pnew % mod;
	}
	return(r);
}

// calculate a^n%mod using 'bigints'   
unsigned __int64 modPower(unsigned __int64 a, unsigned __int64 n,
	unsigned __int64 mod) {
	static mpz_t al, ml, res;
	unsigned __int64 rl;
	static bool firsttime = true;

	if (a < INT_MAX && mod < INT_MAX && n < INT_MAX)
		return mPowerInt((unsigned int)a, (unsigned int)n, (unsigned int)mod);

	if (firsttime) {		// is this the 1st time modPower was called?
		mpz_inits(al, res, ml, NULL);	// if so, allocate storage for bigints
		firsttime = false;
	}

	mpz_set_ui(al, a);    // copy a to a bigint
	mpz_set_ui(ml, mod);  // copy mod to a bigint

	mpz_powm_ui(res, al, n, ml);  // calculate res = a^n%mod
	rl = mpz_get_ui(res);		// since res <= mod, it will not overflow on conversion

	return rl;
}

// n-1 = 2^s * d with d odd by factoring powers of 2 from n-1
static bool witness(unsigned __int64 n, unsigned int s, unsigned __int64 d,
	unsigned __int64 a)
{
	unsigned __int64 x = modPower(a, d, n);   // calculate a^d%n
	unsigned __int64 y;

	while (s) {
		//y = (x * x) % n;
		y = modMult(x, x, n);		// do modular multiplication, works even if x*x would overflow
		//y %= n;
		if (y == 1 && x != 1 && x != n - 1)
			return false;
		x = y;
		--s;
	}
	if (y != 1)
		return false;
	return true;
}

/*
check whether n is prime
* this test is deterministic up to 18,446,744,073,709,551,616 = 2^64
* if n <                  1,373,653,        it is enough to test a = 2 and 3;
* if n <                  9,080,191,        it is enough to test a = 31 and 73;
* if n <              4,759,123,141,        it is enough to test a = 2, 7, and 61;
* if n <          1,122,004,669,633,        it is enough to test a = 2, 13, 23, and 1662803;
* if n <          2,152,302,898,747,        it is enough to test a = 2, 3, 5, 7, 11;
* if n <          3,474,749,660,383,        it is enough to test a = 2, 3, 5, 7, 11, 13;
* if n <        341,550,071,728,321,        it is enough to test a = 2, 3, 5, 7, 11, 13, 17.
* if n <  3,825,123,056,546,413,051,        it is enough to test a = 2, 3, 5, 7, 11, 13, 17, 19, 23.
* if n < 18,446,744,073,709,551,616 = 2^64, it is enough to test a = 2, 3, 5, 7, 11, 13, 17, 19, 23,
*                                           29, 31, and 37
*/
static bool isPrimeMR(unsigned __int64 n)
{
	if (((!(n & 1)) && n != 2)          // if n is even & > 2 return false
		|| (n < 2) 						// if n < 2 return false (i.e  if n=0 or 1)
		|| (n % 3 == 0 && n != 3))		// if n is a multiple of 3 and not equal to 3 return false
		return false;

	// if n is 2 or 3 return true
	if (n <= 3)
		return true;
	// we have now established that n is an odd number > 3
	unsigned __int64 d = n / 2;  // note: n is odd, so 2*d = n-1
	unsigned int s = 1;

	/* get the position of the least significant 1 bit in d. If d were zero,
	result would be zero, and the value in s would be undefined, otherwise
	result is non-zero and s is the bit number of the least significant 1 bit in d.*/
	auto result = _BitScanForward64((unsigned long*)&s, d);  // for gcc compiler use __builtin_ctzll instead
	d >>= s;
	s++;

	// n-1 = 2^s * d with d odd by factoring powers of 2 from n-1
	if (n < 1373653)
		return witness(n, s, d, 2) && witness(n, s, d, 3);
	if (n < 9080191)
		return witness(n, s, d, 31) && witness(n, s, d, 73);
	if (n < 4759123141)
		return witness(n, s, d, 2) && witness(n, s, d, 7) && witness(n, s, d, 61);
	if (n < 1122004669633)
		return witness(n, s, d, 2) && witness(n, s, d, 13) && witness(n, s, d, 23)
		&& witness(n, s, d, 1662803);
	if (n < 2152302898747)
		return witness(n, s, d, 2) && witness(n, s, d, 3) && witness(n, s, d, 5) && witness(n, s, d, 7)
		&& witness(n, s, d, 11);
	if (n < 3474749660383)
		return witness(n, s, d, 2) && witness(n, s, d, 3) && witness(n, s, d, 5) && witness(n, s, d, 7)
		&& witness(n, s, d, 11) && witness(n, s, d, 13);
	if (n < 341550071728321)
		return witness(n, s, d, 2) && witness(n, s, d, 3) && witness(n, s, d, 5) && witness(n, s, d, 7)
		&& witness(n, s, d, 11) && witness(n, s, d, 13) && witness(n, s, d, 17);

	if (n < 3825123056546413051)
		return witness(n, s, d, 2) && witness(n, s, d, 3) && witness(n, s, d, 5) && witness(n, s, d, 7)
		&& witness(n, s, d, 11) && witness(n, s, d, 13) && witness(n, s, d, 17) && witness(n, s, d, 19)
		&& witness(n, s, d, 23);
	// the tests below are sufficient for any 64 bit number
	return witness(n, s, d, 2) && witness(n, s, d, 3) && witness(n, s, d, 5) && witness(n, s, d, 7)
		&& witness(n, s, d, 11) && witness(n, s, d, 13) && witness(n, s, d, 17) && witness(n, s, d, 19)
		&& witness(n, s, d, 23) && witness(n, s, d, 29) && witness(n, s, d, 31) && witness(n, s, d, 37);

}

/* method to return prime divisor for n 
adapted from: 
https://www.geeksforgeeks.org/pollards-rho-algorithm-prime-factorization/ */
long long int PollardRho(long long int n)
{
	/* initialize random seed */
	std::random_device rd;   // non-deterministic generator
	std::mt19937_64 gen(rd());  // to seed mersenne twister.
	std::uniform_int_distribution<long long> dist(1, LLONG_MAX); // distribute results between 1 and MAX inclusive.

	/* no prime divisor for 1 */
	if (n == 1) return n;

	/* even number means one of the divisors is 2 */
	if (n % 2 == 0) return 2;

	/* we will pick from the range [2, N) */
	
	long long int x, y;
	x = dist(gen);  // x is in range 1 to max
	x = x % (n - 2) + 2;  // x is in range 2 to n-1
	y = x; 
	

	/* the constant in f(x).
	 * Algorithm can be re-run with a different c
	 * if it throws failure for a composite. */
	long long int c = dist(gen);  // c is in range 1 to max
		c = c% (n - 1) + 1;  // c is in range 1 to n-1

	/* Initialize candidate divisor (or result) */
	long long int d = 1;

	/* until the prime factor isn't obtained.
	   If n is prime, return n */
	while (d == 1)
	{
		/* Tortoise Move: x(i+1) = f(x(i)) */
		x = (modPower(x, 2, n) + c + n) % n;

		/* Hare Move: y(i+1) = f(f(y(i))) */
		y = (modPower(y, 2, n) + c + n) % n;
		y = (modPower(y, 2, n) + c + n) % n;

		/* check gcd of |x-y| and n */
		d = gcd(abs(x - y), n);

		/* retry if the algorithm fails to find prime factor
		 * with chosen x and c */
		if (d == n) 
			return PollardRho(n);
	}
	if (verbose >0)
		std::cout << "Pollard Rho n = " << n << " factor = " << d << '\n';

	return d;
}

static unsigned int primeFactors(unsigned __int64 tnum, factorsS &f) {
	unsigned  int count = 0, i = 0;

	while (tnum > 1) {

		/* If possible, check whether tnum is prime?*/
		if (tnum <= primeListMax) {
			if (isPrime2(tnum)) { // this check is not essential. It's to save time
				f.factorlist[count][0] = tnum;
				f.factorlist[count][1] = 1;         // set exponent for this factor to 1
				f.factorcount = count + 1;
				return count + 1;
			}
		}
		else {  // if tnum > primeListMax use another method to check if prime
			if (tnum < (primeList[i] * primeList[i])) {
				// logic here is that if tnum were not prime, it would
				// have at least 1 factor < primeList[i], which is is impossible
				// because all such factors have already been removed .
				f.factorlist[count][0] = tnum;
				f.factorlist[count][1] = 1;         // set exponent for this factor to 1
				f.factorcount = count + 1;
				return count + 1;
			}
		}

		/* perform trial division*/
		if (tnum%primeList[i] == 0) {
			f.factorlist[count][0] = primeList[i];   // store prime factor of tnum
			f.factorlist[count][1] = 1;              // set exponent for this factor to 1
			tnum /= primeList[i];
			while (tnum%primeList[i] == 0) {
				f.factorlist[count][1]++;             // increase exponent
				tnum /= primeList[i];
			}
			count++;                                // increase count of unique factors
		}

		i++;										// move index to next prime

		if (i > prime_list_count) {  // all available primes used?
			if (isPrimeMR(tnum)) {
				f.factorlist[count][0] = tnum;
				f.factorlist[count][1] = 1;
				count++;
				break;
			}

			else if (cbrt(tnum) < primeListMax) {
				/* tnum now has no factors <= primeListMax, but it is not prime. Therefore
				it must have exactly two prime factors. We can use the Pollard Rho algorithm
				to get these factors.*/
				long long factor;
				//PollardFactor(tnum, factor);
				factor = PollardRho(tnum);
#ifdef _DEBUG
				assert(tnum%factor == 0);
#endif
				if (factor > 1) {
					tnum /= factor;
					/* save smaller factor 1st in factorlist */
					if (factor <= (long long)tnum) {
						f.factorlist[count][0] = factor;
						f.factorlist[count][1] = 1;
						count++;
						f.factorlist[count][0] = tnum;
						f.factorlist[count][1] = 1;
						count++;
					}
					else {
						f.factorlist[count][0] = tnum;
						f.factorlist[count][1] = 1;
						count++;
						f.factorlist[count][0] = factor;
						f.factorlist[count][1] = 1;
						count++;
					}
					break;
				}
				else {
					ThrowExc("unable to factorise number (pollard Rho failed)");
				}
			}
			else {
				std::cout << "primeListMax = " << primeListMax << '\n';
				ThrowExc("unable to factorise number (too large)");
			}
		}

		/* drop through to here if there are still more primes to check as possible factors 
		(so go round loop again) */
	}

	/* all factors have been found */
	f.factorcount = count;
	return count;
}

/* Calculate maximum square divisor of n and divide n by this.
return adjusted n, square divisor, and factor list of divisor.*/
static void squareFree(__int64 &n, __int64 &sq, factorsS &sqf) {
	factorsS factorlist;

	primeFactors(n, factorlist);
	sqf.factorcount = 0;

	/* make n square-free, factors removed are in sqf */
	for (int i = 0; i < factorlist.factorcount; i++) {
		if (factorlist.factorlist[i][1] >= 2) {
			sqf.factorlist[sqf.factorcount][1] = factorlist.factorlist[i][1] & 0xfffffffffffffffe;
			sqf.factorlist[sqf.factorcount][0] = factorlist.factorlist[i][0];
			//factorlist.factorlist[i][1] -= sqf.factorlist[sqf.factorcount][1];
			n /= power(sqf.factorlist[sqf.factorcount][0], (int)sqf.factorlist[sqf.factorcount][1]);
			sqf.factorcount++;
		}
	}

	sq = FactorMult(sqf);  // calculate maximum square divisor of n
}

/* calculate the number of ways an integer n can be expressed as the sum of 2
squares x^2 and y^2.
see http://mathworld.wolfram.com/SumofSquaresFunction.html,
also http://oeis.org/A004018
generatePrimes must have been called first. Highest prime calculated must be
>= sqrt(n).
example: R2(4) = 4; 2^2+0, 0+2^2, (-2)^2+0, 0+(-2)^2 */
static unsigned __int64 R2(const unsigned __int64 n) {
	if (n == 0)
		return 1;
	if (n % 4 == 3)
		return 0;   // at least 1 4k+3 factor has an odd exponent

	factorsS Rfactors;
	long long b = 1;

	primeFactors(n, Rfactors);

	/* prime factors are split into 3 sets; 2, p=4k+1, p=4k+3, where k is any integer*/
	for (int i = 0; i < Rfactors.factorcount; i++) {

		if (Rfactors.factorlist[i][0] == 2)
			continue; // ignore factor 2
		if (Rfactors.factorlist[i][0] % 4 == 3) {
			/* p = 4k+3 */
			if (Rfactors.factorlist[i][1] % 2 == 1)
				return 0;   // at least 1 4k+3 factor has an odd exponent
		}
		else {
			/* p = 4k + 1 */
			b *= (Rfactors.factorlist[i][1] + 1);
		}
	}
	return 4 * b;
}

/* calculate the number of ways an integer n can be expressed as the sum of 3
squares x^2, y^2 and z^2. The order of the squares is significant. x, y and z can
be +ve, 0 or -ve See https://oeis.org/A005875 
R3(n) = 3*T(n) if n == 1,2,5,6 mod 8, 
       = 2*T(n) if n == 3 mod 8, 
	   = 0 if n == 7 mod 8 and 
	   = R(n/4) if n == 0 mod 4, 
	   where T(n) = 4 times Kronecker's function F(n). [Moreno-Wagstaff].
*/
unsigned __int64 R3(__int64 n) {
	unsigned __int64 sum = 0;
	factorsS sqf;
	__int64 sq, multiplier = 1;

	if (n < 0)
		return 0;
	if (n == 0)			 // test here necessary to avoid infinite loop
		return 1;
	if (n % 8 == 7)
		return 0;        // take short cut if possible
	while (n % 4 == 0)
		n /= 4;         // remove even factors. note that R3(4n) =R3(n)

	generatePrimes(1000003);  // this is more than enough primes to factorise any 64-bit number
	squareFree(n, sq, sqf); /* make n square-free, factors removed are in sqf */

	/* calculate R3(n). For large n this method is not very efficient, but it is simple */
	for (__int64 k = 1; k*k <= n; k++) {
		if (k*k == n)
			sum += 2;  // note: n is a perfect square
		else
			sum += 2 * R2(n - k * k);
	}
	sum += R2(n);  // note: this time (for k=0) we DON'T multiply R2 by 2
	/* we now have sum = R3(n) */

	if (sq > 1) {
		/* add back factors that were removed to make n squarefree, one prime at a time.
		see web.maths.unsw.edu.au/~mikeh/webpapers/paper63.ps specifically equation (3) */
		for (int i = 0; i < sqf.factorcount; i++) {
			auto p = sqf.factorlist[i][0];      // prime 
			auto λ = sqf.factorlist[i][1] / 2;  // exponent/2 (original exponent must be even)
			multiplier *= (power(p, (int)λ + 1) - 1) / (p - 1)
				- jacobi(-n, p)*(power(p, (int)λ) - 1) / (p - 1);
			n *= power(p, (int)λ * 2);
		}
		sum *= multiplier;
	}

	return sum;
}