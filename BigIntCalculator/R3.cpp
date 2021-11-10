#include "pch.h"

#include <map>
#include "factor.h"
#include "showtime.h"

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


extern bool *primeFlags;
extern unsigned long long int primeListMax;
extern unsigned long long *primeList;
extern unsigned int prime_list_count;
void generatePrimes(unsigned long long int max_val);
constexpr unsigned long long int gcd(unsigned long long int u, unsigned long long int v);

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
Znum power(const Znum &x, unsigned long long n) {
	Znum res;
	mpz_pow_ui(ZT(res), ZT(x), n);
	return res;
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

int jacobi(const Znum &k, const Znum &n) {
	return mpz_jacobi(ZT(k), ZT(n));
}

// fast way to check for primes, but generatePrimes must have been called first.
// true indicates a prime number
bool isPrime2(unsigned __int64 num) {

	if (num == 1) return false;
	if (num == 2) return true;
	if ((num & 1) == 0) return false; // even numbers other than 2 are not prime
	return !primeFlags[num/2];
}


/*
* calculates (a * b) % c taking into account that a * b might overflow
*/
unsigned __int64 modMult(unsigned __int64 a, unsigned __int64 b, unsigned __int64 mod)
{
	unsigned __int64 x = 0, retval=0;
	HRESULT  res =0;

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
Znum modMult(const Znum &a, const Znum &b, Znum mod) {
	Znum res;
	res = a * b;
	mpz_mod(ZT(res), ZT(res), ZT(mod));
	return res;
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

// calculate a^n%mod
unsigned __int64 modPowerBi(const Znum &a, const Znum &n, unsigned __int64 mod) {
	Znum res;
	Znum modz = mod;
	unsigned long long r1;

	mpz_powm(ZT(res), ZT(a), ZT(n), ZT(modz));
	r1 = mpz_get_ui(ZT(res));		// since res <= mod, it will not overflow on conversion
	return r1;
}
Znum modPower(const Znum &a, const Znum &n, const Znum &mod) {
	Znum res;
	mpz_powm(ZT(res), ZT(a), ZT(n), ZT(mod));
	return res;
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


/* get prime factors of tnum */
unsigned int primeFactors(unsigned __int64 tnum, factorsS &f) {
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
				factor = PollardRho(tnum);
				assert(tnum%factor == 0);

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
	while (n % 4 == 0)
		n /= 4;         // remove even factors. note that R3(4n) =R3(n)

	/* Interesting note: removing factor 4 above changes the value of n%8 
	so test for n%8 == 7 after doing that, whereas removing other factors to 
	make n square-free would not change the result of this test */
	if (n % 8 == 7)
		return 0;        // take short cut if possible

	generatePrimes(2097169);  // this is more than enough primes to factorise any 64-bit number
	             
	squareFree(n, sq, sqf); /* make n square-free, factors removed are in sqf */

	/* calculate R3(n). For large n this method is not very efficient, but it is simple */
	for (__int64 k = 1; k*k <= n; k++) {
		if (k*k == n)
			sum += 2;  // note: n is a perfect square
		else
			sum += 2 * R2(n - k * k);
		if ((k & 0xfff) == 0) {
			printf_s("%s R3: %g%% done \n", myTime(), 100.0 * double(k) / sqrt(n));
		}
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

/* calculate number of divisors of num, also returns the list of prime factors.
*/
unsigned __int64 NoOfDivisors(__int64 num, factorsS &f) {
	unsigned int count;
	unsigned __int64 result = 1;

	if (num <= 0)
		return 0;      /* num should be > 0 */

	count = primeFactors(num, f);   /* get prime factors of num */

	if (num == 1)
		return 1;  /* this is a special case */

	for (unsigned int i = 0; i < count; i++)
		result *= (f.factorlist[i][1] + 1);
	return result;
}

/* generate list of all divisors of tnum
N.B. includes non-prime factors e.g. 1, 2, 4, 8 and 16 are divisors of 16.
the value returned is the number of divisors.
*/
size_t DivisorList(unsigned __int64 tnum, std::vector <__int64> &divlist) {

	factorsS f;         /* list of factors of tnum*/
	size_t noofdivs;
	size_t ctr = 0, ct2, ccpy;

	/* get number of divisors, and factorlist */
	noofdivs = NoOfDivisors(tnum, f);
	divlist.resize(noofdivs);
	divlist[ctr++] = 1;  // 1st divisor

	for (size_t i = 0; i < f.factorcount; i++) {
		ct2 = 0;
		ccpy = ctr;
		for (int x = 1; x <= f.factorlist[i][1]; x++) {
			for (size_t j = 0; j < ctr; j++)
				divlist[j + ccpy] = divlist[j + ct2] * f.factorlist[i][0];
			ct2 += ctr;
			ccpy += ctr;
		}
		ctr = ccpy;
	}

	std::sort(divlist.begin(), divlist.end());
#ifdef _DEBUG
	/*	printf_s("divisors of %lld are: ", tnum);
		for (int i = 0; i < noofdivs; i++) {
			printf_s("%lld ", divlist[i]);
		}
		putchar('\n');*/
#endif
	return noofdivs;
}

/* return true if smallest factor of num > prime */
static bool minimumFactor(unsigned __int64 num, unsigned __int64 prime) {
	int i = 0;
	unsigned __int64 minfactor;
	for (unsigned i = 0; i < prime_list_count; i++) {
		if (primeList[i] * primeList[i] > num) {
			minfactor = num;   // num is prime
			break;
		}
		if (num%primeList[i] == 0) {
			minfactor = primeList[i];
			break;
		}
	}
	return (minfactor > prime);
}

/* store memoised value for inverse totient here */
static std::map <__int64, std::vector <unsigned __int64>> InvTot;
static long long hitcount = 0;

static void dumpCache(__int64 n, std::map <__int64, std::vector<unsigned __int64>>InvTot) {

	for (auto it : InvTot) {
		printf_s("%3lld", it.first);
		for (auto it2 : it.second) {
			printf_s(" %3lld,", it2);
		}
		putchar('\n');
	}
}

/* create a list of numbers such that the totient of each number in the list
is n. value returned is the number of items in the list, which may be zero.
If the value returned is zero the address in *result may be undefined, otherwise
result points to the list of numbers.
The numbers returned in result are always greater than n, so for large values
of n integer overflow is a possibility.
A rough limit for the maximum value is 2(n^2) and the minimum is n+1. Since there
are upper and lower limits the number of values to be returned is finite.
The algorithm used is based on https://www.insa.nic.in/writereaddata/UpLoadedFiles/IJPAM/20005a81_22.pdf
Alternative link: http://hal.inria.fr/docs/00/71/89/75/PDF/TOT.pdf

Note: this function is recursive, and uses memoisation, otherwise
stack overflow or ridiculously long execution time could be a problem.

To clear the cache afterwards call inverseTotient with n <=0
If debug = true extra information is sent to stdout.
If dump = true some contents of the cache are output.
*/
size_t inverseTotient(__int64 n, std::vector<unsigned __int64> **result, bool debug,
	int level, bool dump) {
	const unsigned __int64 sqrtn = llSqrt(n);

	std::vector <unsigned __int64>rest_solution, *r_s = NULL;
	unsigned __int64 prime_power, rest_phi;
	unsigned __int64 noofdivisors;
	std::vector <unsigned __int64> max_primes;
	unsigned __int64 prime_phi, prime, rest_sol_count;
	size_t primecount = 0;
	std::vector <__int64> divisorList;
	HRESULT rv;

	if (n <= 0) {   // n <= 0 signifies that cache memory is to be freed
		if (dump)
			dumpCache(n, InvTot);

		if (verbose > 0)
			std::cout << "Inverse Totient cache cleared (contained "
			<< InvTot.size() << " entries) \n" << "hit count = " << hitcount << '\n';
		    
		InvTot.clear();
		hitcount = 0;
		*result = nullptr;
		return 0;
	}

	if (sqrtn >= primeListMax)   // generate prime list if necessary
		generatePrimes(sqrtn + 10);

	/* If this is the 1st call: initialise cache. Each entry in the cache consists
	of a list of numbers for which the totient is n. */
	if (InvTot.size() == 0) {

		/* Note: InvTot.insert (std::make_pair(a, b)) creates an entry where
		InvTot[a] contains a vector of size b. */

		// special, set up a list with just 1 value (1)
		InvTot.insert(std::make_pair(0, 1));
		InvTot[0][0] = 1;

		// only 1 and 2 have 1 as totient value
		InvTot.insert(std::make_pair(1, 2));
		InvTot[1][0] = 1;
		InvTot[1][1] = 2;


		// 3, 4 and 6 have 2 as totient value
		InvTot.insert(std::make_pair(2, 3));
		InvTot[2][0] = 3;
		InvTot[2][1] = 4;
		InvTot[2][2] = 6;
	}

	/* is the result we need already in the cache? */
	auto cp = InvTot.find(n);

	if (cp != InvTot.end()) {
		hitcount++;
		*result = &cp->second;
		return cp->second.size();  // return number of numbers in result
	}

	if ((n & 1) == 1) {
		*result = nullptr;
		return 0;  // there are no odd totient values other than 1
	}

	if (debug || dump)
		printf_s("evaluating inv.tot(%lld) level %d\n", n, level);

	/* add result about to be calculated to cache*/
	InvTot.insert(std::make_pair(n, 0));  // initially just have key value and zero-length vector

#ifdef _DEBUG
	if (dump)
		dumpCache(n, InvTot);
#endif

	/* get all the divisors of n */
	noofdivisors = DivisorList(n, divisorList);

	max_primes.clear();              // reset size to zero
	max_primes.reserve(noofdivisors / 2);  // guess the new size
	/* make a list of all primes which are  = some divisor+1 */
	for (int i = 0; i < noofdivisors; i++) {
		/* test whether divisor+1 is prime. use fast prime check if possible */
		if (((divisorList[i] + 1 < (__int64)primeListMax) && isPrime2(divisorList[i] + 1)) ||
			((divisorList[i] + 1 >= (__int64)primeListMax) && (isPrimeMR(divisorList[i] + 1))))
			max_primes.push_back(divisorList[i] + 1);
	}
	primecount = max_primes.size();

	/* step through the list of primes, starting with the largest */
	for (ptrdiff_t ii = primecount - 1; ii >= 0; ii--) {

		prime = max_primes[ii];
		prime_phi = prime - 1; /* totient of a prime p = p-1 */
		prime_power = prime;  /* initially, prime_power = prime^1 */

		while (n%prime_phi == 0) {
			//#ifdef _DEBUG
			//		printf_s("ii = %d, prime =%lld, prime_phi=%lld, prime_power=%lld\n",
			//			ii, prime, prime_phi, prime_power);
			//
			//#endif
			rest_phi = n / prime_phi;

			/* note: this is a recursive call. If prime is 2, prime_phi is 1 and
			rest_phi = n so in this case inverseTotient will return the values
			previously calculated. This is why we start with the largest prime
			and end with 2.
			(If inverseTotient does not find the result needed in the cache, it
			will need to calculate it & add it to the cache ) */
			rest_sol_count = inverseTotient(rest_phi, &r_s, debug, level + 1, dump);
			if (rest_sol_count != 0) {
#ifdef _DEBUG
				if (debug) {
					printf_s("inv.tot of %lld are: ", rest_phi);
					for (int i = 0; i < std::min(rest_sol_count, 10ULL); i++)
						printf_s("%lld, ", (*r_s)[i]);
					printf_s("\ncount= %lld\n", rest_sol_count);
				}
#endif
				/* copy inverse totients to local variable rest_solution. This is
				a deep copy, so later changes to list referenced by r_s don't matter.*/
				rest_solution = *r_s;

				// find all numbers in rest_solution which have no prime factor <= prime
				for (int iii = 0; iii < rest_sol_count; iii++) {
					if (rest_solution[iii] > 1) {  // unless number being factored is 1
						if (minimumFactor(rest_solution[iii], prime)) {
							unsigned __int64 temp;
							rv = ULongLongMult(prime_power, rest_solution[iii],
								&temp);
							if (rv != S_OK) {
								ThrowExc("integer overflow: ")
							}
							InvTot[n].push_back(temp);  // add new inv. totient value
#ifdef _DEBUG
							if (debug) {
								printf_s("inv. totient for %lld: %lld ",
									n, temp);
								printf_s("(%lld * %lld)\n", prime_power, rest_solution[iii]);
							}
#endif
						}
					}
					else {    /* rest_solution[iii] is 1. In this case we DON'T test
							  whether it is greater than prime */
						InvTot[n].push_back(prime_power);
#ifdef _DEBUG
						if (debug)
							printf_s("inv. totient for %lld: %lld\n",
								n, prime_power);
#endif
					}
				}  // end of for loop
			} // end of "if (rest_sol_count != 0)"

			if (prime > n / prime_phi)
				break;     // try to avoid overflow for large primes

			/* calculate next power of prime */
			rv = ULongLongMult(prime_power, prime, &prime_power); //prime_power *= prime;
			if (rv != S_OK) {
				ThrowExc("integer overflow: ")
			}
			prime_phi *= prime;   // totient of p^k = p^(k-1) * (p-1)
		} // end of "while (n%prime_phi == 0)"
	}  // end of "for (ptrdiff_t ii = primecount - 1; ii >= 0; ii--)"

	/* sort values of inverse totient (not really necessary)*/
	std::sort(InvTot[n].begin(), InvTot[n].end());

	*result = &InvTot[n]; /* note that function returns a pointer into
							the cache, which therefore cannot be freed yet*/
	if (debug || dump) {
		if (InvTot.size() <= 12)
			dumpCache(n, InvTot);
		printf_s("** completed invTot(%lld) level %d\n", n, level);
	}

	return InvTot[n].size();  // return number of numbers in result
}
