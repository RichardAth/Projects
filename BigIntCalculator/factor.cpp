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
#include <intrin.h>
#include "showtime.h"
#include "factor.h"
#include "bignbr.h"
#undef min                 // use std::min

static long long Gamma[386];
static long long Delta[386];
static long long AurifQ[386];
static Znum Zfactor;
double originalTenthSecond;
int oldTimeElapsed;

static int Cos(int N) {
	switch (N % 8) 	{
	case 0:
		return 1;
	case 4:
		return -1;
	}
	return 0;   // default value 
}

static int intTotient(int N) {
	int totient, q, k;

	totient = q = N;
	if (q % 2 == 0) {
		totient /= 2;

		do {
			q /= 2;
		} while (q % 2 == 0);
	}

	if (q % 3 == 0) {
		totient = totient * 2 / 3;
		do {
			q /= 3;
		} while (q % 3 == 0);
	}

	k = 5;
	while (k * k <= q) {
		if (k % 3 != 0 && q % k == 0) {
			totient = totient * (k - 1) / k;
			do {
				q /= k;
			} while (q % k == 0);
		}
		k += 2;
	}

	if (q > 1) {
		totient = totient * (q - 1) / q;
	}
	return totient;
}

static int Moebius(int N) {
	int moebius, q, k;

	moebius = 1;
	q = N;

	if (q % 2 == 0) {
		moebius = -moebius;
		q /= 2;
		if (q % 2 == 0) {
			return 0;
		}
	}

	if (q % 3 == 0) {
		moebius = -moebius;
		q /= 3;
		if (q % 3 == 0) {
			return 0;
		}
	}
	k = 5;
	while (k * k <= q) {
		if (k % 3 != 0) {
			while (q % k == 0) {
				moebius = -moebius;
				q /= k;
				if (q % k == 0) {
					return 0;
				}
			}
		}
		k += 2;
	}

	if (q > 1) {
		moebius = -moebius;
	}
	return moebius;
}

/* power = base^exponent */
static void BigIntPowerIntExp(const Znum &Base, int exponent, Znum &Power) {
	mpz_pow_ui(ZT(Power), ZT(Base), exponent);
}

static void GetAurifeuilleFactor(fList &Factors, int L, const Znum &Base,
	const int DegreeAurif) {
	Znum x, Csal, Dsal, Nbr1;  
	int k;

	BigIntPowerIntExp(Base, L, x);   // x <- BigBase^L.
	Csal = 1;      
	Dsal = 1;       
	for (k = 1; k < DegreeAurif; k++) {
		Nbr1 = Gamma[k]; 
		Csal *= x;       
		Csal += Nbr1;     // Csal <- Csal * x + Gamma[k]
		Nbr1 = Delta[k]; 
		Dsal *= x;      
		Dsal += Nbr1;        // Dsal <- Dsal * x + Gamma[k]
	}
	Nbr1 = Gamma[k];   
	Csal = Csal*x;      
	Csal += Nbr1;        // Csal <- Csal * x + Gamma[k]
	BigIntPowerIntExp(Base, (L + 1) / 2, Nbr1);   
	Nbr1 = Nbr1 * Dsal;       // Nbr1 <- Dsal * base^((L+1)/2)
	Dsal = Csal + Nbr1; 
 
	insertBigFactor(Factors, Dsal);
	Dsal = Csal - Nbr1; 

	insertBigFactor(Factors, Dsal);
	return;
}

/* Get Aurifeuille factors.
 see https://en.wikipedia.org/wiki/Aurifeuillean_factorization */
static void InsertAurifFactors(fList &Factors, const Znum &Base,
	int Expon, int Incre)
{
	int DegreeAurif;
	if (Base >= 386) {
		return;    // Base is very big, so go out.
	}
	auto llBase = MulPrToLong(Base);

	if (Expon % 2 == 0 && Incre == -1) {
		do {
			Expon /= 2;
		} while (Expon % 2 == 0);

		Incre = llBase % 4 - 2;
	}

	if (Expon % llBase == 0
		&& Expon / llBase % 2 != 0
		&& ((llBase % 4 != 1 && Incre == 1) || (llBase % 4 == 1 && Incre == -1)))
	{
		int N1, q, L, k;
		int N = (int)llBase;
		if (N % 4 == 1) {
			N1 = N;
		}
		else {
			N1 = 2 * N;
		}
		DegreeAurif = intTotient(N1) / 2;
		for (k = 1; k <= DegreeAurif; k += 2) {
			AurifQ[k] = JacobiSymbol(N, k);
		}
		for (k = 2; k <= DegreeAurif; k += 2) {
			int t1 = k; // Calculate t2 = gcd(k, N1)
			int t2 = N1;
			while (t1 != 0) {
				int t3 = t2 % t1;
				t2 = t1;
				t1 = t3;
			}
			AurifQ[k] = (long long)Moebius(N1 / t2) * intTotient(t2) * Cos((N - 1) * k);
		}
		Gamma[0] = Delta[0] = 1;
		for (k = 1; k <= DegreeAurif / 2; k++) {
			int j;
			Gamma[k] = Delta[k] = 0;
			for (j = 0; j < k; j++) {
				Gamma[k] =
					Gamma[k]
					+ N * AurifQ[2 * k
					- 2 * j
					- 1] * Delta[j]
					- AurifQ[2 * k
					- 2 * j] * Gamma[j];
				Delta[k] =
					Delta[k]
					+ AurifQ[2 * k
					+ 1
					- 2 * j] * Gamma[j]
					- AurifQ[2 * k
					- 2 * j] * Delta[j];
			}
			Gamma[k] /= 2LL * k;
			Delta[k] = (Delta[k] + Gamma[k]) / (2LL * k + 1);
		}

		for (k = DegreeAurif / 2 + 1; k <= DegreeAurif; k++) {
			Gamma[k] = Gamma[DegreeAurif - k];
		}

		for (k = (DegreeAurif + 1) / 2; k < DegreeAurif; k++) {
			Delta[k] = Delta[DegreeAurif - k - 1];
		}

		q = Expon / (int)Base;
		L = 1;
		while (L * L <= q) {
			if (q % L == 0) {
				GetAurifeuilleFactor(Factors, L, Base, DegreeAurif);
				if (q != L * L) {
					GetAurifeuilleFactor(Factors, q / L, Base, DegreeAurif);
				}
			}
			L += 2;
		}
	}
	return;
}


static void Cunningham(fList &Factors, const Znum &Base, int Expon,
	const int increment, const Znum &Original) {
	int Expon2, k;
	Znum Nbr1, Nbr2; // Temp1;

	Expon2 = Expon;

	while (Expon2 % 2 == 0 && increment == -1) {
		Expon2 /= 2;
		BigIntPowerIntExp(Base, Expon2, Nbr1);  /* Nbr1 = Base^Expon2*/
		Nbr1 += increment; 
		insertBigFactor(Factors, Nbr1);
		InsertAurifFactors(Factors, Base, Expon2, 1);
	}

	k = 1;
	while (k * k <= Expon) {
		if (Expon % k == 0) {
			if (k % 2 != 0) { /* Only for odd exponent */
				BigIntPowerIntExp(Base, Expon / k, Nbr1); /* nbr1 = base^(expon/k) */
				Nbr1 += increment;      
				Nbr2 = gcd(Nbr1, Original);   // Nbr2 <- gcd(Base^(Expon/k)+incre, original)
				insertBigFactor(Factors, Nbr2);
				//Temp1 = Original; 
				Nbr1 = Original / Nbr2; 
				insertBigFactor(Factors, Nbr1);
				InsertAurifFactors(Factors, Base, Expon / k, increment);
			}

			if ((Expon / k) % 2 != 0) { /* Only for odd exponent */
				BigIntPowerIntExp(Base, k, Nbr1);  /* Nbr1 = Base^k*/
				Nbr1 += increment; 
				Nbr2 = gcd(Nbr1, Original);   // Nbr2 <- gcd(Base^k+incre, original)
				insertBigFactor(Factors, Nbr2);
				//Temp1 = Original; 
				Nbr1 = Original / Nbr2; 
				insertBigFactor(Factors, Nbr1);
				InsertAurifFactors(Factors, Base, k, increment);
			}
		}
		k++;
	}
	return;
}


/* check whether the number +/- 1 is a perfect power*/
static void PowerPM1Check(fList &Factors, const Znum &nbrToFactor,
	long long MaxP) {

	int Exponent = 0;
	Znum base;

	/* code below finds cases where nbr +/- 1 is a perfect power.
	If base < MaxP we may not find it here, but the factor will then
	be found anyway by trial division. */

	Exponent = (int)PowerCheck(nbrToFactor + 1, base, MaxP);
	if (Exponent != 1) {
		/* we have base^exp = nbrToFactor + 1
		i.e. nbrToFactor = base^exp -1 = (base-1) * (base^(exp-1) + .... +1)
		i.e. base-1 is a factor */
		Cunningham(Factors, base, Exponent, -1, nbrToFactor);
		return;    // number is a perfect power - 1
	}

	Exponent = (int)PowerCheck(nbrToFactor - 1, base, MaxP);
	if (Exponent != 1) {
		/* we have base^exp = nbrToFactor - 1 */
		Cunningham(Factors, base, Exponent, 1, nbrToFactor);
		return;    // number is a perfect power + 1
	}

	return;  // number is not a perfect power +/- 1
}


/* sort factors into ascending order. If two factors are equal, merge them by 
adding the 2nd exponent to the 1st then removing the second entry and moving any 
following entries up 1 position */
static void SortFactors(fList &Factors) {
	std::sort(Factors.f.begin(), Factors.f.end());
/* if two factors have the same value they are merged.  */
	for (ptrdiff_t i = 0; i < (ptrdiff_t)Factors.f.size()-1; i++) {
		if (Factors.f[i + 1] == Factors.f[i]) {
				/* factors i and i+1 are equal so merge them */
			Factors.f[i].exponent += Factors.f[i+1].exponent;  // combine the exponents
			if (Factors.f[i+1].upperBound == -1)  
				Factors.f[i].upperBound = -1; // set upperbound to show factor is prime
			else if (Factors.f[i].upperBound != -1
				&& Factors.f[i].upperBound < Factors.f[i+1].upperBound)
				/* use higher value of upperbound. */
				Factors.f[i].upperBound = Factors.f[i+1].upperBound;

			/* now remove entry i+1 & move any entries higher than i+1 down 1. */
			Factors.f.erase(Factors.f.begin() + i + 1);
			i--; /* tweak loop counter so that if there were 3 or more factors
				  with the same value they would all be merged */
		}
	}

	/* for certain numbers it is possible that a spurious factor 1 is generated. 
	In a way this is valid, but 1 is not considered a prime number so it is
	removed */
	if (Factors.f[0].Factor == 1) {
		Factors.f.erase(Factors.f.begin());
	}
	if (verbose > 1) {
		std::cout << "result after sort" << '\n';
		Factors.Xprint();
	}
}

/* Insert new factor found into factor array. */
/* assume divisor is prime. ix is index of non-prime factor which is a multiple 
of divisor. either pi is the index into the prime list of the divisor, or the  
divisor is in div */
static void insertIntFactor(fList &Factors, int pi, long long div, ptrdiff_t ix) {
	auto lastfactor = Factors.f.size();
	Znum quot, qnew;
	Znum divisor;
	if (pi >= 0)
		divisor = primeList[pi];
	else
		divisor = div;

	quot = Factors.f[ix].Factor;
	/* divide quot by divisor as many times as possible */
	auto exp = mpz_remove(ZT(qnew), ZT(quot), ZT(divisor));
	if (qnew != 1) {
		/* add new factor */
		Factors.f.resize(lastfactor + 1);  // increase size of factor list
		Factors.f[lastfactor].exponent = (int)exp * Factors.f[ix].exponent;
		Factors.f[lastfactor].Factor = divisor;
		Factors.f[lastfactor].upperBound = -1; // show that new factor is prime

		Factors.f[ix].Factor = qnew;  // adjust value of original factor
		Factors.f[ix].upperBound = pi;
	}
	else {
		/* replace residue of 1 with new factor */
		Factors.f[ix].Factor = divisor;
		Factors.f[ix].exponent *= (int)exp;
		Factors.f[ix].upperBound = -1;
	}
	SortFactors(Factors);
	if (verbose > 0) {
		/*std::cout << "result after adding factor " << divisor << '\n';
		Factors.Xprint();*/
	}
}

/* Insert new factor found into factor array. Every current factor is checked 
against the divisor. If their gcd is not 1 the existing factor is divided by 
the gcd and a new factor equal to the gcd is created. */
void insertBigFactor(fList &Factors, const Znum &divisor) {
	auto lastfactor = Factors.f.size();
	auto ipoint = lastfactor;
	zFactors temp;
	Znum g;
	if (verbose > 1) {
		std::cout << "InsertBigFactor Divisor =" << divisor << '\n';
		Factors.Xprint();
	}
	for (int i = 0; i < lastfactor; i++) {
		if (Factors.f[i].Factor == divisor)
			continue;  // factor already found, but continue checking 
		g = gcd(Factors.f[i].Factor, divisor);
		if (g != 1 && g < Factors.f[i].Factor) {
			Znum qnew;
			Factors.f.resize(lastfactor + 1);  // increase size of factor list
			/* we can replace Factor with 2 factors, Factor/g and g 
			(if Factor is a multiple of divisor, g = divisor) */
			mp_bitcnt_t fexp = mpz_remove(ZT(qnew), ZT(Factors.f[i].Factor), ZT(g));
			Factors.f[i].Factor = qnew;   //Factors[i].Factor /= g^fexp;
			Factors.f[ipoint].Factor = g;
			Factors.f[ipoint].exponent = Factors.f[i].exponent*(int)fexp;
			Factors.f[ipoint].upperBound = Factors.f[i].upperBound;
			ipoint++;
			lastfactor++;
		}
	}
	if (verbose > 1) {
		std::cout << "result before sort" << '\n';
		Factors.Xprint();
	}
	SortFactors(Factors);
}

/* called for Carmichael numbers that have no small factors. 
Return: false = No factors found, true = factors found.
 Use: Xaux for square root of -1.
      Zaux for square root of 1. */
static bool factorCarmichael(const Znum &pValue, fList &Factors)
{
	Znum randomBase = 0;  // pseudo-random number
	bool factorsFound = false;
	int countdown, ctr;
	bool sqrtOneFound = false;
	bool sqrtMinusOneFound = false;
	Znum Aux1, Aux2, Aux3, Aux4, Xaux, Zaux, Temp4;
		
	Aux1 = pValue - 1;  // Aux1 = p - 1 (p is odd, so Aux1 is even)
	DivideBigNbrByMaxPowerOf2(ctr, Aux1);  // Aux1 /= 2^ctr
	
	for (countdown = 20; countdown > 0; countdown--) {
		int i;
		randomBase = ((uint64_t)randomBase * 89547121 + 1762281733) & 0x7fffffff;
		 // Aux2 = base^Aux1 (mod p)
		mpz_powm(ZT(Aux2), ZT(randomBase), ZT(Aux1), ZT(pValue));
		// If Mult1 = 1 or Mult1 = p-1, then try next base.
		if (Aux2 ==1 || Aux2 == pValue-1) {
			continue;    // This base cannot find a factor. Try another one.
		}
		for (i = 0; i < ctr; i++)
		{              // Loop that squares number.
			// Aux3 = Aux2^2 (mod p)
			mpz_powm_ui(ZT(Aux3), ZT(Aux2), 2, ZT(pValue));
			if (Aux3 == 1)
			{            // Non-trivial square root of 1 found.
				if (!sqrtOneFound)
				{          // Save it to perform GCD later.
					Zaux = Aux2;
					sqrtOneFound = true;
				}
				else
				{          // Try to find non-trivial factor by doing GCD.
					Aux4 = Aux2 - Zaux;
					mpz_mod(ZT(Aux4), ZT(Aux4), ZT(pValue));// Aux4 = Aux2 - Zaux (mod p)
					Temp4 = gcd(pValue, Aux4);
					if ((Temp4 > 1) && (Temp4 != pValue))
					{          // Non-trivial factor found.
						insertBigFactor(Factors, Temp4);
						factorsFound = true;
					}
				}
				// Try to find non-trivial factor by doing GCD.
				Aux4 = Aux2 + 1;
				Aux4 %= pValue;             // Aux4 = Aux2+1 (mod p)
				Temp4 = gcd(pValue, Aux4);  // Temp4 = gcd(p, Aux4)
				if ((Temp4 > 1) && (Temp4 != pValue))
				{          // Non-trivial factor found.
					insertBigFactor(Factors, Temp4);
					factorsFound = true;
				}
				i = ctr;
				continue;  // Find more factors.
			}
			if (Aux3 == pValue-1)
			{            // Square root of 1 found.
				if (!sqrtMinusOneFound)
				{          // Save it to perform GCD later.
					Xaux = Aux2;
					sqrtOneFound = true;
				}
				else
				{          // Try to find non-trivial factor by doing GCD.
					Aux4 = Aux3 - Xaux;
					Aux4 %= pValue; // Aux4 = Aux3 - Xaux (mod p)
					Temp4 = gcd(pValue, Aux4);  // Temp4 = gcd(p, Aux4)
					if ((Temp4 > 1) && (Temp4 != pValue))
					{          // Non-trivial factor found.
						insertBigFactor(Factors, Temp4);
						factorsFound = true;
					}
				}
				i = ctr;
				continue;  // Find more factors.
			}
			Aux2 = Aux3;
		}
	}
	return factorsFound;
}

/*
see: http://lemire.me/blog/2013/12/26/fastest-way-to-compute-the-greatest-common-divisor/
calculate GCD of a and b. 
If using gcc compiler use __builtin_ctzll instead of _BitScanForward64
note use of long long ints to avoid overflow.

u and v are UNSIGNED, caller must use llabs or similar if necessary.
Note: mathematically gcd(a,b) = gcd(-a,b) = gcd(b,-a) = gcd (-a,-b)
Note: GCD (greatest common denominator) is also known as HCF (Highest Common Factor).
*/
constexpr unsigned long long int gcd(unsigned long long int u, unsigned long long int v)
{
	unsigned char result=0;
	unsigned long shift=0, su=0, sv=0;
	if (u == 0) return v;
	if (v == 0) return u;
	result = _BitScanForward64(&shift, u | v);  // count any common factor 2s
	result = _BitScanForward64(&su, u);
	u >>= su;             // shift u until u is odd
	do {
		result = _BitScanForward64(&sv, v);
		v >>= sv;          // shift v until v is odd
		if (u > v) {
			unsigned long long int t = v;
			v = u - v;
			u = t;
		}
		else
			v = v - u;
	} while (v != 0LL);
	return u << shift;
}

/* factorise number where we know that it only has 2 prime factors. 
see https://en.wikipedia.org/wiki/Pollard%27s_rho_algorithm 
Note: there is another (probably better) version called PollardRho. */
static void PollardFactor(const unsigned long long num, long long &factor) {
	long long x_fixed = 2, cycle_size = 2, x = 2; 
	factor = 1;
	while (factor == 1) {
		for (long long count = 1; count <= cycle_size && factor <= 1; count++) {
			/* even if x*x overflows, the function as a whole still works OK */
			x = (x*x + 1) % num;
			factor = gcd(abs(x - x_fixed), num);
		}
		if (factor == num) {
			/* there is a small possibility that PollardFactor won't work,
			even when factor is not prime */
			std::cout << "Pollard factorisation failed for num = " << num 
				<< " cycle_size = " << cycle_size << " x = " << x << " !!\n";
			factor = 1;
			return;   // factorisation failed!! 	
		}
		cycle_size *= 2;
		x_fixed = x;
	}
	if (verbose > 0) {
		std::cout << "Pollard Factor. num = " << num << " factor = " << factor
			<< " cycle_size = " << cycle_size << " x = " << x << '\n';
	}
	return;
}

/* method to return prime divisor for n
adapted from:
https://www.geeksforgeeks.org/pollards-rho-algorithm-prime-factorization/
This method generally works but for very large n it may be too slow. It uses a
truly random number generator so could give different results given the same
value of n */
long long int PollardRho(long long int n, int depth)
{
	/* initialize random seed */
	std::random_device rd;   // non-deterministic generator
	std::mt19937_64 gen(rd());  // to seed mersenne twister.
	std::uniform_int_distribution<long long> dist(1, LLONG_MAX); // distribute results between 1 and MAX inclusive.
	long long ctr = 0;     /* loop counter*/
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
	c = c % (n - 1) + 1;  // c is in range 1 to n-1

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
			return PollardRho(n, depth + 1);
		ctr++;
	}

	if (verbose > 0)
		std::cout << "Pollard Rho n = " << n << " factor = " << d
		<< " loop counter =" << ctr << " depth=" << depth << '\n';

	return d;
}

/* trial division & Pollard-Rho. Uses first 33333 primes */
static void TrialDiv(fList &Factors, const long long PollardLimit) {
	bool restart = false;  // set true if trial division has to restart
	int upperBound;
	long long testP;
	do {  /* use trial division */
		restart = false;    // change to true if any factor found
		for (int i = 0; i < Factors.f.size(); i++) {
			upperBound = Factors.f[i].upperBound;  // resume from where we left off
			if (upperBound == -1)
				continue;  // factor is prime
			/* trial division. Uses first 33333 primes */
			while (upperBound < std::min((int)prime_list_count, 33333)) {
				testP = primeList[upperBound];
				if (testP*testP > Factors.f[i].Factor) {
					Factors.f[i].upperBound = -1; // show that residue is prime
					break;
				}
				if (Factors.f[i].Factor%testP == 0) {
					insertIntFactor(Factors, upperBound, 0, i);
					restart = true;  // factor found; keep looking for more
					Factors.tdiv++;
					break;
				}
				Factors.f[i].upperBound = upperBound;
				upperBound++;
			}

			if (!restart && (Factors.f[i].upperBound != -1)
				&& Factors.f[i].Factor <= PollardLimit) {
				long long f;
				if (PrimalityTest(Factors.f[i].Factor, primeList[Factors.f[i].upperBound]) == 0)
					/* if factor is prime calling PollardFactor would waste a LOT of time*/
					Factors.f[i].upperBound = -1; // Indicate that number is prime.
				else {
					/* as factor is not prime, and it has no factors < MaxP, it must have
					just two prime factors. */
					if (verbose > 1) {
						std::cout << "factors before Pollard factorisation: ";
						Factors.Xprint();
					}
					f = PollardRho(MulPrToLong(Factors.f[i].Factor));
					/* there is a small possibility that PollardFactor won't work,
					even when factor is not prime*/
					if (f != 1) {
						insertIntFactor(Factors, -1, f, i);
						Factors.prho++;
					}
				}
			}
		}
	} while (restart);  // keep looping until no more factors found.

	if (verbose > 0) {
		std::cout << "End Trial division. " << Factors.f.size() - 1 << " factors found so far \n";
		if (Factors.f.size() > 1) {
			std::cout << "result after trial division ";
			Factors.Xprint();
		}
	}
}

/* factorise toFactor; factor list returned in Factors. */
static bool factor(const Znum &toFactor, fList &Factors) {
	long long testP;
	const long long MaxP = 393'203;  // use 1st  33333 primes
	/* larger value seems to slow down factorisation overall. */
	// MaxP must never exceed 2,097,152 to avoid overflow of PollardLimit
	const long long PollardLimit = MaxP*MaxP*MaxP;

	Factors.set(toFactor);   /* initialise factor list */
	if (toFactor <= 3)
		return true;   /* toFactor is 1 or prime so we can stop now*/

	if ((long long)primeListMax <MaxP) {  // get primes
		generatePrimes(MaxP);  // takes a while, but only needed on 1st call
	}
	
	oldTimeElapsed = 0;
	originalTenthSecond = tenths();    // record start time

	if (toFactor >= MaxP* MaxP) {
		/* may not be able to factorise entirely by trial division, so try this first */
		PowerPM1Check(Factors, toFactor, MaxP/2);  // check if toFactor is a perfect power +/- 1
		Factors.pm1 = (int)Factors.f.size() - 1;  // number of factors just found, if any
		if (Factors.f.size() > 1 && verbose > 0) {
			std::cout << "PowerPM1Check result: ";
			Factors.Xprint();
		}
	}

	// If toFactor is < PollardLimit it will be factorised completely using 
    // trial division and Pollard-Rho without ever using ECM or SIQS factorisiation. 
	TrialDiv(Factors, PollardLimit);
	/* Any small factors (up to 393,203) have now been found by trial division */

	for (ptrdiff_t i = 0; i < (ptrdiff_t)Factors.f.size(); i++) {
	/* Check whether the residue is a prime or prime power.
	Given that the residue is less than about 10^10,000 the maximum exponent is
	less than 2000.  e.g. 3000007^1826 has 10,002 digits */
		Znum Zpower;
		int expon;
		if (Factors.f[i].upperBound == -1)
			continue;         // skip if this factor is known to be prime

		Zfactor = Factors.f[i].Factor;
		testP = primeList[Factors.f[i].upperBound]; // get largest number used in trial division
		expon = (int)PowerCheck(Zfactor, Zpower,  testP - 1); // on return; Zpower^expon = Zfactor
		if (expon > 1) {    /* if factor is a perfect power*/
			Factors.f[i].Factor = Zpower;
			Factors.f[i].exponent *= expon;
			Factors.power++;
		}
		
		int result = PrimalityTest(Zpower, testP- 1);
		if (result == 0) {   // Number is prime.
			Factors.f[i].upperBound = -1; // Indicate that number is prime.
			continue;
		}
		if (result > 1) {  /* number is a pseudo-prime */
			size_t fsave = Factors.f.size();
			if (factorCarmichael(Zpower, Factors)) {
				Factors.carm += (int)(Factors.f.size() - fsave); // record any increase in number of factors;
				i = -1;			// restart loop at beginning!!
				continue;
			}
		}

		if (Zpower <= PollardLimit) {
			long long f;
			f = PollardRho(MulPrToLong(Zpower));
			if (f != 1) {
				insertIntFactor(Factors, -1, f, i);
				Factors.prho++;
				/* there is a small possibility that PollardFactor won't work,
				even when factor is not prime*/
				continue;
			}
		}

		int nooflimbs = numLimbs(Zpower);
		// get approximate size (1 limb = 64 bits)
		if (nooflimbs <=3 || (!msieve && !yafu)) {
			/* use built-in ECM & SIQS if number to factor <= 192 bits (58 digits)
			   because this is fastest for smaller numbers,
			   or if both YAFU and Msieve are turned off */
			ElipCurvNo = 1;  // start with 1st curve
			auto rv = ecm(Zpower, Factors, Zfactor);          
			// get a factor of number. result in Zfactor
			if (!rv)
				return false;  // failed to factorise number
			 //Check whether factor is not one. In this case we found a proper factor.
			if (Zfactor != 1) {
				insertBigFactor(Factors, Zfactor);
				i = -1;			// restart loop at beginning!!
			}
		}
		else {
			/* First try to factor N using Lehman algorithm. Result in Zfactor.
			This seldom achieves anything, but when it does it saves a lot of time.
			If N has 2 factors and the larger factor is < 10x smaller factor
			this should find a factor, so this complements the elliptic curve 
			method which is better for finding smaller factors. */
			for (int k = 1; k <= 10; k++) {
				LehmanZ(Zpower, k, Zfactor);
				if (Zfactor > 1) {
					Factors.leh++;     // Factor found.
					insertBigFactor(Factors, Zfactor);
					i = -1;   // success; restart loop at beginning to tidy up!	
					break;
				}
			}
			if (i == -1)
				continue;   // restart loop from beginning

			size_t fsave = Factors.f.size();
			bool rv; 
			/* msieve and yafu should never both be set. At this point one of them 
			should be set. */
			if (msieve)  
				rv = callMsieve(Zpower, Factors);
			else
				rv = callYafu(Zpower, Factors);
			if (rv) {
				i = -1;   // success; restart loop at beginning to tidy up!
				// record any increase in number of factors
				if (msieve)
					Factors.msieve += (int)(Factors.f.size() - fsave); 
				else
					Factors.yafu += (int)(Factors.f.size() - fsave); 
			}
			else {
				msieve = false;   // failed once, don't try again
				yafu = false;
				if (verbose > 0) {
					std::cout << "Msieve or YAFU failed: turn on built-in ECM/SIQS \n";
				}
				i--;
			}
		}
	}

	SortFactors(Factors);  // tidy up factor list 
	return true;
}

/* compute 4 or less values the squares of which add up to prime p,  
return values in Mult1, Mult2, Mult3 and Mult4 */
static void ComputeFourSquares(const Znum &p, Znum &Mult1, Znum &Mult2,
	Znum &Mult3, Znum &Mult4) {
	Znum a, q, K, Tmp, Tmp1, Tmp2, Tmp3, Tmp4, M1, M2, M3, M4; 
	Znum TestNbr;

	if (p == 2) {   /* Prime factor is 2 */
		Mult1 = 1;  // 2 = 1^2 + 1^2 + 0^2 + 0^2
		Mult2 = 1;
		Mult3 = 0;
		Mult4 = 0;
	}
	else {       /* Prime factor p is not 2 */
		if ((p & 3) == 1) { /* if p = 1 (mod 4) */
			q = (p-1)/4; // q = (prime-1)/4
	  
			a = 1;
			do {    // Loop that finds mult1^2 = (-1) mod p
				a++; 
				assert(a < p);
				/* Mult1 = a^q(mod p) */
				mpz_powm(ZT(Mult1), ZT(a), ZT(q), ZT(p));
				mpz_powm_ui(ZT(TestNbr), ZT(Mult1), 2, ZT(p));
			} while (TestNbr != p-1 && TestNbr != -1);

			Mult1 = abs(Mult1);
			Mult2 = 1;

			for (;;) {  
				Tmp = Mult1 * Mult1 + Mult2 * Mult2; // K <- (mult1^2 + mult2^2) / p
				K = Tmp / p;    // in this loop K is smaller each time round   
				if (K == 1) {  // are we there yet?
					Mult3 = 0;
					Mult4 = 0;
					break;     // we are finished
				}
				M1 = Mult1 % K;
				if (M1 < 0) {
					M1 += K;
				}
				M2 = Mult2 % K;
				if (M2 < 0) {
					M2 += K;
				}
				Tmp = (K+1)/2; // Tmp <- (K+1) / 2
				//subtractdivide(Tmp, -1, 2);       
				if (M1 >= Tmp) // If M1 >= K / 2 ... 
				{
					M1 -= K;
				}
				if (M2 >= Tmp) {     // If M2 >= K / 2 ...     
					M2 -= K;
				}
				Tmp = Mult1*M1 + Mult2*M2;
				Tmp2 = Tmp / K;  // Tmp2 <- (mult1*m1 + mult2*m2) / K
				Tmp = Mult1*M2 - Mult2*M1;
				Mult2 = Tmp / K;    // Mult2 <- (mult1*m2 - mult2*m1) /K
				Mult1 = Tmp2;       // Mult1 <- (mult1*m1 + mult2*m2) / K
			} /* end while */
		} /* end p = 1 (mod 4) */

		else { /* if p = 3 (mod 4) */
			q = (p-1)/2; //  q = (prime-1)/2
			Mult1 = 0;
			do {
				Mult1++;
				Tmp = Mult1*Mult1 + 1;
				Tmp = -Tmp;
				while (Tmp < 0) Tmp += p;
				mpz_powm(ZT(Tmp1), ZT(Tmp), ZT(q), ZT(p));
				       // At this moment Tmp1 = (-1 - Mult1^2)^((p-1)/2)(Mod p)
			} while (Tmp1 != 1);  // Continue loop if it is not 1.

			// After the loop finishes, Tmp1 = (-1 - Mult1^2) is a quadratic residue mod p.
			q = (p+1)/4;

			// Find Mult2 <- square root of Tmp1 = Tmp^q (mod p) 
			mpz_powm(ZT(Mult2), ZT(Tmp), ZT(q), ZT(p));

			Mult3 = 1;
			Mult4 = 0;

			for (;;) {
				// Compute K <- (Mult1^2 + Mult2^2 + Mult3^2 + Mult4^2) / p
				Tmp = Mult1*Mult1 + Mult2*Mult2;
				Tmp = Tmp + Mult3*Mult3 + Mult4*Mult4;
				assert(Tmp%p == 0);
				K = Tmp / p;    // in this loop K is smaller each time round
				assert(K > 0);
				if (K == 1) {
					break;  // we are done when K equals 1
				}


				if (ZisEven(K)) { // If K is even ...
					if (ZisEven(Mult1) != ZisEven(Mult2))
					{  // If Mult1 + Mult2 is odd...
						if (ZisEven(Mult1) == ZisEven(Mult3))
						{   // If Mult1 + Mult3 is even...
							Tmp = Mult2;  // swap mult2 and mult3
							Mult2 = Mult3;
							Mult3 = Tmp;
						}
						else {
							Tmp = Mult2; // swap mult2 and mult4
							Mult2 = Mult4;
							Mult4 = Tmp;
						}
					} // At this moment Mult1+Mult2 = even, Mult3+Mult4 = even
					Tmp1 = (Mult1 + Mult2)/2;   
					Tmp2 = (Mult1 - Mult2)/2;   
					Tmp3 = (Mult3 + Mult4)/2;   
 					Mult4 = (Mult3 - Mult4)/2 ;
					Mult3 = Tmp3;
					Mult2 = Tmp2;
					Mult1 = Tmp1;
					continue;
				} /* end if K is even */

				M1 = Mult1 % K;
				if (M1 < 0) {
					M1 += K;
				}
				M2 = Mult2 % K;
				if (M2 < 0) {
					M2 = M2 + K;
				}
				M3 = Mult3 % K;
				if (M3 < 0) {
					M3 += K;
				}
				M4 = Mult4 % K;
				if (M4 < 0) {
					M4 += K;
				}
				Tmp = (K+1)/2;  // Tmp <- (K+1) / 2
				if (M1 >= Tmp) { // If M1 >= K / 2 ... 
					M1 -= K;     // M1 = M1 - K;    
				}
				if (M2 >= Tmp) { // If M2 >= K / 2 ... 
					M2 -= K;     // M2 = M2 - K;         
				}

				if (M3 >= Tmp) {  // If M3 >= K / 2 ... 
					M3 -= K;      // M3 = M3 - K;        
				}
				if (M4 >= Tmp) { // If M4 >= K / 2 ... 
					M4 -= K;     //M4 = M4 - K;        
				}
				// Compute Tmp1 <- (Mult1*M1 + Mult2*M2 + Mult3*M3 + Mult4*M4) / K
				Tmp = Mult1*M1 + Mult2*M2 + Mult3*M3 + Mult4*M4;
				Tmp1 = Tmp / K; // BigIntDivide(Tmp, K, Tmp1);

				// Compute Tmp2 <- (Mult1*M2 - Mult2*M1 + Mult3*M4 - Mult4*M3) / K
				Tmp = Mult1*M2 - Mult2*M1 + Mult3*M4 - Mult4*M3;
				Tmp2 = Tmp / K;

				// Compute Tmp3 <- (Mult1*M3 - Mult3*M1 - Mult2*M4 + Mult4*M2) / K
				Tmp = Mult1*M3 - Mult3*M1 - Mult2*M4 + Mult4*M2;
				Tmp3 = Tmp / K;

				// Compute Mult4 <- (Mult1*M4 - Mult4*M1 + Mult2*M3 - Mult3*M2) / K
				Tmp = Mult1*M4 - Mult4*M1 + Mult2*M3 - Mult3*M2; // BigIntSubt(Tmp, Tmp4, Tmp);
				Mult4 = Tmp / K;

				Mult3 = Tmp3;
				Mult2 = Tmp2;
				Mult1 = Tmp1;
			} /* end while */
		} /* end if p = 3 (mod 4) */
	} /* end prime not 2 */
}

/* compute 3 values the squares of which add up to s *2^2r, return values in quads */
static void compute3squares(int r, const Znum &s, Znum quads[4]) {
	Znum s2, s3, r2, Tmp1, Tmp2;
	int m = 0;

	if (s == 3) {
		quads[0] = quads[1] = quads[2] = 1;
		quads[3] = 0;
		for (int ix = 0; ix <= 2; ix++)
			mpz_mul_2exp(ZT(quads[ix]), ZT(quads[ix]), r);
		return;
	}

	/* loop till we find a suitable value of s3. */
	for (Znum x = 0; ; x++) {
		assert(x*x < s);
		s2 = s - x * x;
		for (s3 = s2, m = 0; ZisEven(s3); m++) {
			s3 >>= 1;    // s3 = s2*2^m
		}
		/* we know s3 is odd, need to check whether s3 mod 4 = 1 */
		if ((s3 & 3) != 1)
			continue;
		/* In general, to establish whether or not s3 can be expressed as the sum of 2 
		squares, it would be necessary to factorise it and examine all the factors.
		However, if s3 is prime, we know it can be so expressed, and can easily find
		two squares. If s3 is not prime we keep on looking. */
#ifdef __MPIR_VERSION
		static bool first = true;
		static gmp_randstate_t rstate;
		if (first) {
			gmp_randinit_default(rstate);
			first = false;
		}

		auto rv = mpz_probable_prime_p(ZT(s3), rstate, 16, 0);
#else
		auto rv = mpz_probab_prime_p(ZT(s3), 16);
#endif
		//auto rv = mpz_bpsw_prp(ZT(s3)); /* rv = 0 for composite, 1 = probable prime, 2 = definite prime*/
		if (rv != 0) {
			/* s3 is prime of form 4k+1 */
			ComputeFourSquares(s3, quads[0], quads[1], quads[2], quads[3]);
			quads[0] = abs(quads[0]);  // ensure results are +ve
			quads[1] = abs(quads[1]);
			if (quads[0] < quads[1]) {		
				Tmp1 = quads[0];           // quads[0] < quads[1], so exchange them.
				quads[0] = quads[1];
				quads[1] = Tmp1;
			}
			assert(quads[2] == 0);
			assert(quads[3] == 0);

			/* put back factor 2 removed earlier */
			mpz_mul_2exp(ZT(quads[0]), ZT(quads[0]), m/2);
			mpz_mul_2exp(ZT(quads[1]), ZT(quads[1]), m/2);
			if ((m & 1) == 1) {
				/* if m is odd the sum needs to be multiplied by 2.
				We do this by using the formula 
				2*(q0^2 + q1^2) = (q0+q1)^2 + (q0-q1)^2 */
				Tmp1 = quads[0] + quads[1];
				Tmp2 = quads[0] - quads[1];
				quads[0] = Tmp1;
				quads[1] = Tmp2;
			}
			
			quads[2] = x;

			for (int ix = 0; ix <= 2; ix++)
				mpz_mul_2exp(ZT(quads[ix]), ZT(quads[ix]), r);

			if (quads[0] < quads[1]) {
				Tmp1 = quads[0];		// quads[0] < quads[1], so exchange them.
				quads[0] = quads[1];
				quads[1] = Tmp1;
			}

			if (quads[0] < quads[2]) {
				Tmp1 = quads[0];	// quads[0] < quads[2], so exchange them.
				quads[0] = quads[2];
				quads[2] = Tmp1;
			}

			if (quads[1] < quads[2]) {
				Tmp1 = quads[1];	// quads[1] < quads[2], so exchange them.
				quads[1] = quads[2];
				quads[2] = Tmp1;
			}

			return;
		}
	}
}

/* show that the number is the sum of 4 or fewer squares. See
https://www.alpertron.com.ar/4SQUARES.HTM */
/* uses the identity:
(a^2+b^2+c^2+d^2)*(A^2+B^2+C^2+D^2) = (aA+bB+cC+dD)^2 + (aB-bA+cD-dC)^2
                                    + (aC-bD-cA+dB)^2 + (aD-dA+bC-cB)^2 
This allows us to find the sum of squares for each factor separately then combine them */
static void ComputeFourSquares(const fList &factorlist, Znum quads[4], Znum num) {
	Znum Mult1, Mult2, Mult3, Mult4, Tmp1, Tmp2, Tmp3;
	Znum pr;
	bool twoSq = true;

	quads[0] = 1;      // 1 = 1^2 + 0^2 + 0^2 + 0^2
	quads[1] = 0;
	quads[2] = 0;
	quads[3] = 0;

	if (factorlist.f.size() == 1) {/* only 1 factor? */
		if (factorlist.f[0].Factor == 1) {   // Number to factor is 1.
			return;
		}
		if (factorlist.f[0].Factor == 0) {    // Number to factor is 0.
			quads[0] = 0;      // 0 = 0^2 + 0^2 + 0^2 + 0^2
			return;
		}
	}

	/* check whether number can be formed as sum of 1 or 2 squares */
	for (auto f: factorlist.f) {
		if ((f.Factor & 3) == 3) {  // is factor of form 4k+3?
			if ((f.exponent & 1) == 1)   // and is the exponent of that factor odd?
			twoSq = false;    // if yes, number cannot be expressed as sum of 1 or 2 squares
			break;       
		}
	}
	if (!twoSq) {  /* check whether number can be expressed as sum of 3 squares */
		int r = 0;
		while ((num & 3) == 0) {
			num /= 4;
			r++;
		}
		/* any number which is not of the form 4^r * (8k+7) can be formed as the sum of 3 squares 
		see https://en.wikipedia.org/wiki/Legendre%27s_three-square_theorem*/
		if ((numLimbs(num) < 4) && (num & 7) < 7) {
			/* use compute3squares if number is small (<= 57 digits) and can be 
			formed from 3 squares. (Each limb is up to 64 bits) */
			compute3squares(r, num, quads);  
			return;
		}
	}

	for (auto Factorx : factorlist.f) {
		if (Factorx.exponent % 2 == 0) {
			continue; /* if Prime factor appears an even number of times, no need to
					  process it in this for loop */
		}
		
		pr = Factorx.Factor;

		/* compute 4 or less values the squares of which add up to prime pr,
		return values in Mult1, Mult2, Mult3 and Mult4 */
		ComputeFourSquares(pr, Mult1, Mult2, Mult3, Mult4);
		//assert(pr == Mult1*Mult1 + Mult2*Mult2 + Mult3*Mult3 + Mult4*Mult4);

		/* use the identity:
		(a^2+b^2+c^2+d^2)*(A^2+B^2+C^2+D^2) = (aA+bB+cC+dD)^2 + (aB-bA+cD-dC)^2
		                                    + (aC-bD-cA+dB)^2 + (aD-dA+bC-cB)^2 
		note: if twoSq is true, c, d, C, & D are zero. The expression then 
		simplifies 	automatically to:  	
		       (a^2+b^2)*(A^2+B^2) = (aA+bB)^2 + (aB-bA)^2 */

		Tmp1 = Mult1*quads[0] + Mult2*quads[1] + Mult3*quads[2] + Mult4*quads[3];
		Tmp2 = Mult1*quads[1] - Mult2*quads[0] + Mult3*quads[3] - Mult4*quads[2];
		Tmp3 = Mult1*quads[2] - Mult3*quads[0] - Mult2*quads[3] + Mult4*quads[1];
		quads[3] = Mult1*quads[3] - Mult4*quads[0] + Mult2*quads[2] - Mult3*quads[1];

		quads[2] = Tmp3;
		quads[1] = Tmp2;
		quads[0] = Tmp1;
	} 

	 /* for factors that are perfect squares, multiply quads[0]-3 by sqrt(factor) */
	for (auto Factorx : factorlist.f) {
		if (Factorx.Factor >= 2) {
			mpz_pow_ui(ZT(Tmp1), ZT(Factorx.Factor), 
				Factorx.exponent / 2);
			quads[0] *= Tmp1;
			quads[1] *= Tmp1;
			quads[2] *= Tmp1;
			quads[3] *= Tmp1;
		}
	}

	quads[0] = abs(quads[0]);  // ensure results are +ve
	quads[1] = abs(quads[1]);
	quads[2] = abs(quads[2]);
	quads[3] = abs(quads[3]);

	/* Sort squares: largest in quads[0], smallest in quads[3]. This
	is equivalent to a 'bubble sort' with the loops unrolled. There are
	only 6 comparisons & exchanges for 4 items. */

	if (quads[0] < quads[1]) {		// quads[0] < quads[1], so exchange them.
		Tmp1 = quads[0];
		quads[0] = quads[1];
		quads[1] = Tmp1;
	}

	if (quads[0] < quads[2]) {	// quads[0] < quads[2], so exchange them.
		Tmp1 = quads[0];
		quads[0] = quads[2];
		quads[2] = Tmp1;
	}

	if (quads[0] < quads[3]) {	// quads[0] < quads[3], so exchange them.
		Tmp1 = quads[0];
		quads[0] = quads[3];
		quads[3] = Tmp1;
	}

	if (quads[1] < quads[2]) {	// quads[1] < quads[2], so exchange them.
		Tmp1 = quads[1];
		quads[1] = quads[2];
		quads[2] = Tmp1;
	}

	if (quads[1] < quads[3]) {	// quads[1] < quads[3], so exchange them.
		Tmp1 = quads[1];
		quads[1] = quads[3];
		quads[3] = Tmp1;
	}

	if (quads[2] < quads[3]) {	// quads[2] < quads[3], so exchange them.
		Tmp1 = quads[2];
		quads[2] = quads[3];
		quads[3] = Tmp1;
	}
	return;
}

/********************************************************************************
code below is to interface between DA's code
and the new code that uses Znums, which are really mpz_t integers from MPIR or GMP
multiprecision library, with a c++ class wrapped around them that allows them
to be used pretty much like normal integers. 
**********************************************************************************/

/* factorise number. Returns false if unable to factorise it */
bool factorise(Znum numberZ, fList &vfactors, Znum quads[]) {

	try {
		bool pos = true;
		if (numberZ == 0)
			return false;  // function factor can't factorize zero
		if (numberZ < 0) {
			pos = false;
			numberZ = -numberZ;
		}
		auto rv = factor(numberZ, vfactors);
		if (!rv)
			return false;  // failed to factorise number
		
		if (quads != nullptr) {
			ComputeFourSquares(vfactors, quads, numberZ); 
			// get a, b, c, d such that sum of their squares = numberZ
		}
		return true;
	}

	/* code below catches C++ 'throw' type exceptions */
	catch (const std::exception& e) {
		fprintf_s(stderr, "\n*** a standard exception was caught, with message\n '%s'\n", e.what());

		return false;
	}
}
