﻿#include "pch.h"
#include "factor.h"

__int64 ChineseRem(__int64 a1, __int64 n1, __int64 a2, __int64 n2);

// calculate x^n. Returns a Big Integer so should be no overflow problem.
Znum powerBi(const __int64 x, unsigned __int64 n) {
	Znum result;
	mpz_ui_pow_ui(ZT(result), x, n);
	return result;
}

// returns 2 ^ exp. exp can be any number >= 0
Znum pow2bi(unsigned __int64 exp) {
	Znum result = 1;
	mpz_mul_2exp(ZT(result), ZT(result), exp);
	return result;
}


// simulate parallel assignment. temp is used to ensure correct result
// if newa and b are the same variable. It is recommended that all variables 
// are the same type.
#define PA(a,b,newa,newb) \
temp = (newa);  \
b = (newb);     \
a = temp;

/*
the extended Euclidean algorithm is an extension to the Euclidean algorithm,
which computes, besides the greatest common divisor of integers a and b, the
coefficients of Bézout's identity, that is integers x and y such that

	  a * x + b * y = gcd ( a , b )
*/
__int64 extendedGcd(const __int64 a, const __int64 b, __int64 *x, __int64 *y) {
	__int64 s = 0, old_s = 1;;
	__int64 t = 1, old_t = 0;
	__int64 r = b, old_r = a;
	__int64 quotient, temp, qs, qt;
	HRESULT rv;

	while (r != 0) {
		quotient = old_r / r;
		PA(old_r, r, r, old_r - quotient * r)

			rv = LongLongMult(quotient, s, &qs);
		if (rv != S_OK) {
			ThrowExc("integer overflow: ")
		}
		PA(old_s, s, s, old_s - qs)

			rv = LongLongMult(quotient, t, &qt);
		if (rv != S_OK) {
			ThrowExc("integer overflow: ")
		}
		PA(old_t, t, t, old_t - qt)
	}
	//output "Bézout coefficients:", (old_s, old_t)
	//output "greatest common divisor:", old_r
	//output "quotients by the gcd:", (t, s)
	if (old_r < 0) {
		*x = -old_s;
		*y = -old_t;
		return -old_r;    // -r = gcd(a,b)
	}
	else {
		*x = old_s;
		*y = old_t;
		return old_r;     // r = gcd(a,b)
	}
}
#undef PA

/* /* find x such that x ≡ a1 (mod n1) and x ≡ a2 (mod n2)
	x will be in the range 0 to n1*n2
	This extended version handles cases where n1 and n2 are not coprime.*/
__int64 ExtChineseRem(__int64 a1, unsigned __int64 n1, __int64 a2, unsigned __int64 n2) {
/*	n1 and n2 are factorised, the 2 original congruences are expanded into a set of
	congruences, one for each prime factor and the set is then solved if possible.
	if n1 and n2 are not coprime they will have 1 or more common prime factors.
	the congruence arising from n1 must be compatible with the congruence from n2
	for the same prime factor, otherwise there is no solution.

	If there is no solution the value zero is returned. */

	factorsS n1F, n2F;
	__int64  n1fc;  /* number of prime factors of n1*/
	__int64  n2fc;  /* number of prime factors of n2*/
	unsigned __int64 cfactors[34][3] = { 0 }, minexp, x;  /* combined list based on
															factors of n1 and n2 */
	unsigned __int64 A1, A2, N1, N2;
	int n1x = 0, n2x = 0, cx = 0;     /* indices into n1factors, n2factors and cfactors */
	HRESULT ret;					/* return code from intsafe multiplication */

	n1fc = primeFactors(n1, n1F);   // get prime factors of n1
	n2fc = primeFactors(n2, n2F);   // get prime factors of n2

	/* create a combined list from the prime factors of n1 and n2 */
	while ((n1x < n1fc) || (n2x < n2fc)) {  /* exit loop only when both lists of factors
											 have been processed */

		if ((n1x < n1fc) && (n2x < n2fc) &&
			/* do n1 and n2 have a common factor? */
			(n2F.factorlist[n2x][0] == n1F.factorlist[n1x][0])) {
			/* 1st check whether the common factors can be combined */
			minexp = std::min(n2F.factorlist[n2x][1], n1F.factorlist[n1x][1]);
			minexp = power(n1F.factorlist[n1x][0], (unsigned int)minexp);
			if (a1%minexp != a2 % minexp) {
				/* factors cannot be combined */
				printf_s("ext C.R. a1=%lld, n1=%lld, a2=%lld, n2=%lld\n", a1, n1, a2, n2);
				printf_s("prime is %lld ", n2F.factorlist[n2x][0]);
				printf_s(" exponents: %lld & %lld\n",
					n1F.factorlist[n1x][1], n2F.factorlist[n2x][1]);
				printf_s("a1 mod %lld = %lld   a2 mod %lld = %lld\n", minexp, a1%minexp, minexp,
					a2%minexp);
				return 0;  // there is no solution!
			}

			/* combine the factors */
			cfactors[cx][0] = n1F.factorlist[n1x][0];  // copy prime
			cfactors[cx][1] = std::max(n1F.factorlist[n1x][1], n2F.factorlist[n2x][1]); // use max exponent
			if (n1F.factorlist[n1x][1] > n2F.factorlist[n2x][1])
				cfactors[cx][2] = a1;   // use a that corresponds to max exponent
			else {
				cfactors[cx][2] = a2;
			}
			cx++;
			n1x++;
			n2x++;
		}

		if (n1x >= n1fc) {
			while (n2x < n2fc) {  // all factor of n1 done, copy rest of n2
				cfactors[cx][0] = n2F.factorlist[n2x][0];
				cfactors[cx][1] = n2F.factorlist[n2x][1];
				cfactors[cx][2] = a2;
				cx++;
				n2x++;
			}
			break;
		}

		if (n2x >= n2fc) {
			while (n1x < n1fc) {  // all factors of n2 done, copy rest of n1
				cfactors[cx][0] = n1F.factorlist[n1x][0];
				cfactors[cx][1] = n1F.factorlist[n1x][1];
				cfactors[cx][2] = a1;
				cx++;
				n1x++;
			}
			break;
		}

		/* both n1 and n2 have unprocessed factors, copy smallest first
		so that equal prime factors will come together and can be combined */

		/* while prime factors of n1  < lowest unprocessed factor of n2,
		copy factors from n1 to combined */
		while ((n1x < n1fc) && (n1F.factorlist[n1x][0] < n2F.factorlist[n2x][0])) {
			cfactors[cx][0] = n1F.factorlist[n1x][0];
			cfactors[cx][1] = n1F.factorlist[n1x][1];
			cfactors[cx][2] = a1;
			cx++;
			n1x++;
		}

		/* while prime factors of n2  < lowest unprocessed factor of n1,
		copy factors from n2 to combined */
		while ((n2x < n2fc) && (n2F.factorlist[n2x][0] < n1F.factorlist[n1x][0])) {
			cfactors[cx][0] = n2F.factorlist[n2x][0];
			cfactors[cx][1] = n2F.factorlist[n2x][1];
			cfactors[cx][2] = a2;
			cx++;
			n2x++;
		}
	}

	/* we now have a combined list of prime factors. We now process the factors
	 by combining them until we have a solution.
	*/

	A1 = cfactors[0][2];
	N1 = power(cfactors[0][0], (unsigned int)cfactors[0][1]);
	for (int i = 1; i < cx; i++) {

		A2 = cfactors[i][2];
		N2 = power(cfactors[i][0], (unsigned int)cfactors[i][1]);
		x = ChineseRem(A1, N1, A2, N2);
		/*	printf_s("ext C.R x=%llu A1=%llu, N1=%llu, A2=%llu, N2=%llu i=%d cx=%d\n",
				x, A1, N1, A2, N2, i, cx);*/

		A1 = x;
		//N1 *= N2;
		ret = ULongLongMult(N1, N2, &N1);
		if (ret != S_OK) {
			ThrowExc("integer overflow: ")
		}
	}
	//printf_s("ext C R: x=%lld a1=%lld, n1=%lld, a2=%lld, n2=%lld\n", x, a1, n1, a2, n2);
	assert(x%n1 == a1);   // check for bugs
	assert(x%n2 == a2);
	return x;
}

/* find x such that x ≡ a1 (mod n1) and x ≡ a2 (mod n2)
	x will be in the range 0 to n1*n2 
	if there is no solution an exception will be thrown */
__int64 ChineseRem(__int64 a1, __int64 n1, __int64 a2, __int64 n2) {
	__int64 m1, m2, llgcd, llgcd2 = 0, x = 0;
	static Znum bt1, bt2, bx, bm1, bm2, t3;

/* We use the extended Euclidian algorithm to find integers m1 and m2 such that
m1*n1 + m2 * n2 = gcd(n1, n2)
A solution is given by x = a1 * m2*n2 + a2 * m1*n1, provided that n1 and n2 are co-prime

if n1 and n2 are not co-prime,
x = (a1*m2*n2 + a2 * m1*n1) / gcd(n1, n2)   ONLY if both a1 and a2 are multiples of the gcd,
otherwise there may be no solution.*/

	llgcd = extendedGcd(n1, n2, &m1, &m2);
	if (llgcd != 1) {
		if ((a1%llgcd != 0) || (a2%llgcd != 0)) {
			if (a1%llgcd != a2 % llgcd) {
				char buf[4000];  /* guess how big buffer should be; if it's too small the
								 message will be truncated */
				snprintf(buf, sizeof(buf), "Chinese Rem: no solution for %lld, %lld, %lld, %lld",
					a1, n1, a2, n2);
				ThrowExc(buf);  /* throw an exception */
			}
			else {
				x = ExtChineseRem(a1, n1, a2, n2);  // solution exists, use ext version to get it
				if (x == 0) {
					printf_s("CH: R *** x=%lld a1=%lld, n1=%lld, a2=%lld, n2=%lld\n", x, a1, n1, a2, n2);
					abort();  /* having just checked that there is a solution, logically
							   x cannot be 0 */
				}
				return x;
			}
		}
		else {   /* in this simple case, we can deal with it here although
				 n1 and n2 are not coprime. a1 and a2 can be divided exactly */
			a1 /= llgcd;
			a2 /= llgcd;
			n1 /= llgcd;        // remove common factors so n1 and n2 are coprime
			n2 /= llgcd;
			llgcd2 = extendedGcd(n1, n2, &m1, &m2);  // recalculate m1 and m2
		}
	}
	//x = a1*m2*n2 + a2*m1*n1;
	bx = (Znum)a1*m2*n2 + a2 * m1*n1;

	// use bigints to avoid overflow in intermediate products

	// ensure value of x is in correct range, 0 to n1*n2-1
	t3 = abs((Znum)(n1*n2));
	if (bx < 0) {
		mpz_fdiv_r(ZT(bx), ZT(bx), ZT(t3));
	}
	if (bx >= t3)
		mpz_fdiv_r(ZT(bx), ZT(bx), ZT(t3));

	if (llgcd != 1) {
		bx *= llgcd;
	}
	x = MulPrToLong(bx);         // convert result back to a normal integer

	return x;
}

/* find x such that x ≡ a1 (mod n1) and x ≡ a2 (mod n2)
x will be in the range 0 to n1*n2. n1 and n2 must be mutually prime
We use the extended Euclidian algorithm to find integers m1 and m2 such that
m1*n1 + m2*n2 = gcd(n1,n2)
A solution is given by x = a1*m2*n2 + a2*m1*n1, provided n1 and n2 are co-prime
N.B. if n1 and n2 are not co-prime an exception will be thrown
*/
void ChineseRem(const Znum &a1, const Znum &n1, const Znum &a2, const Znum &n2, Znum &x) {

	Znum m1, m2, gcd, t2;
	mpz_gcdext(ZT(gcd), ZT(m1), ZT(m2), ZT(n1), ZT(n2));
	if (gcd != 1) {
		char buf[4000];  /* guess how big buffer should be; if it's too small the
						 message will be truncated */
		gmp_snprintf(buf, sizeof(buf), "Chinese Rem: no solution for %Zd, %Zd, %Zd, %Zd",
			a1, n1, a2, n2);
		ThrowExc(buf);  /* throw an exception */
	}
	x = a1 * m2*n2 + a2 * m1*n1;

	t2 = abs(n1*n2);
	if (x < 0) {    // if x < 0 add a multiple of n1*n2
		//gmp_printf("x=%Zd, t2=%Zd\n", x, t2);
		mpz_fdiv_r(ZT(x), ZT(x), ZT(t2));  // ensure x is in required range
	}
	if (x >= t2) {  // if x n1*n2 subtract a multiple of n1*n2
		//gmp_printf("x=%Zd, t2=%Zd\n", x, t2);
		mpz_fdiv_r(ZT(x), ZT(x), ZT(t2));
	}
	return;
}

/* find x such that x^2 ≡ c mod p^lambda 
  the method used was partly derived from Wikipedia but the cases where c=0, 
  and where p=2 are all my own work. I have not seen descriptions that cover these
  cases let alone working code that does what this does ANYWHERE and, believe 
  me, I looked. */
std::vector <long long> ModSqrt2(long long c, const unsigned long long p, 
	const int lambda) {
	std::vector <long long> roots;
	long long r1, r2, root, x;
	//long long a, Beta;
	long long mod = power(p, lambda);
	std::vector <long long> rx;

	c %= mod;
	if (c < 0)
		c += mod;  /* ensure c is in range 0 to mod-1 */

	/* treat c=0 as special case */
	if (c == 0) {
		r1 = power(p, (lambda + 1) / 2);  /* smallest non-zero root is the
										   smallest power of p >= sqrt(mod) */
		r2 = 0;
		while (r2 < mod) {
			roots.push_back(r2);  
			r2 += r1;  /* the roots form an arithmetic progression: 0, r1, 2*r1, 3*r1 etc */
		}
		return roots;
	}

	/* must treat p=2 as a special case */
	if (p == 2 && c != 0) {
		if (c% mod == 1) {
			roots.push_back(1);  
			roots.push_back(mod - 1);  /* get 2nd root */
		}
		return roots;
	}

	// code below (derived from Wikpedia) works for p=4i+1 but the more general
	// method is better IMHO because the intermediate products are smaller so less
	// overflow possibilities to deal with
	//if (p % 4 == 1) {
	//	a = (p - 1) / 4;
	//	Beta = power(p, lambda - 1) * a;   /* β = a* p^(λ-1) */
	//	r1 = modPowerBi(powerBi(c, a)+3, Beta, mod);  /* (c^a  +3)^β (mod p^λ)*/
	//	r2 = modPower(c, (Beta + 1) / 2, mod);       /* c^((β+1)/2) (mod p^λ)      */
	//	root = modMult(r1, r2, mod);    /* get final product */
	//	roots.push_back(root);
	//	roots.push_back(mod - root);    /* get 2nd root */
	//}
	//else {

	/* this part was derived from Wikipedia 
	see https://en.wikipedia.org/wiki/Tonelli%E2%80%93Shanks_algorithm#Tonelli's_algorithm_will_work_on_mod_p^k
	The explanation there is not very clear but I got working code out of it */
		rx = primeModSqrt(c, p);   /* rx^2 ≡ c (mod p) */
		if (rx.empty()) {
			roots.clear();
			return roots;     /* there are no solutions */
		}
		if (rx.size() > 1)
			x = std::min(rx[0], rx[1]);  /* take smaller root */
		else
			x = rx[0];
		r1 = modPower(x, power(p, lambda - 1), mod);
		r2 = modPowerBi(c, (power(p, lambda) - 2 * power(p, lambda - 1) + 1) / 2, mod);
		root = modMult(r1, r2, mod);    /* get final product */
		roots.push_back(root);
		roots.push_back(mod - root);    /* get 2nd root */
	//}
	return roots;
}

std::vector <Znum> ModSqrt2(const Znum &cc, const Znum &p, const int lambda) {
	std::vector <Znum> roots;
	Znum r1, r2, root, x, c;
	//long long a, Beta;
	Znum modz = power(p, lambda);
	long long mod = MulPrToLong(modz);  /* could get overflow? */
	std::vector <Znum> rx;

	c  = cc%mod;
	if (c < 0)
		c += mod;  /* ensure c is in range 0 to mod-1 */

	/* treat c=0 as special case */
	if (c == 0) {
		r1 = power(p, (lambda + 1) / 2);  /* smallest non-zero root is the
										   smallest power of p >= sqrt(mod) */
		r2 = 0;
		while (r2 < mod) {
			roots.push_back(r2);
			r2 += r1;  /* the roots form an arithmetic progression: 0, r1, 2*r1, 3*r1 etc */
		}
		return roots;
	}

	/* must treat p=2 as a special case */
	if (p == 2 && c != 0) {
		if (c% mod == 1) {
			roots.push_back(1);
			roots.push_back(mod - 1);  /* get 2nd root */
		}
		return roots;
	}


	/* this part was derived from Wikipedia
	see https://en.wikipedia.org/wiki/Tonelli%E2%80%93Shanks_algorithm#Tonelli's_algorithm_will_work_on_mod_p^k
	The explanation there is not very clear but I got working code out of it */
	rx = primeModSqrt(c, p);   /* rx^2 ≡ c (mod p) */
	if (rx.empty()) {
		roots.clear();
		return roots;     /* there are no solutions */
	}
	if (rx.size() > 1)
		x = std::min(rx[0], rx[1]);  /* take smaller root */
	else
		x = rx[0];
	r1 = modPower(x, power(p, lambda - 1), mod);
	r2 = modPowerBi(c, (power(p, lambda) - 2 * power(p, lambda - 1) + 1) / 2, mod);
	root = modMult(r1, r2, mod);    /* get final product */
	roots.push_back(root);
	roots.push_back(mod - root);    /* get 2nd root */

	return roots;
}

/*
Square root modulo prime number using Tonelli–Shanks algorithm
Solve the equation given a and p.
	x^2 ≡ a mod p
and return list of solutions. There will be either 0, 1 or 2 solutions
see https://en.wikipedia.org/wiki/Tonelli%E2%80%93Shanks_algorithm
(renamed the solution variable n to a)
*/
std::vector <long long> primeModSqrt(long long a, const unsigned long long p) {
	std::vector <long long> result;
	unsigned long long q, s, z, m, i, e;
	unsigned long long c, t, R, b;

	a %= p;
	if (a < 0)
		a += p;   /* normalise a so it's in range 0 to p-1 */

	// Simple case
	if (a == 0) {
		result.push_back(0);
		return result;
	}

	if (p == 2) {
		result.push_back(a);  // a is 0 or 1
		return result;
	}

	/* Check solution existence on odd prime. Because p is prime the Jacobi
	symbol is the same as the Legendre symbol. */
	if (jacobi(a, p) != 1)
		return result;    // empty list; no solutions
#ifdef _DEBUG
	/* recheck existence of solution */
	assert(modPower(a, (p - 1) / 2, p) == 1);
#endif

	// Simple case
	if (p % 4 == 3) {
		R = modPower(a, (p + 1) / 4, p);
		result.push_back(R);
		result.push_back(p - R);
		return result;
	}

	// step 1: Factor p - 1 of the form q * 2 ^ s(with Q odd)
	q = p - 1;
	s = 0;
	while (q % 2 == 0) {
		s += 1;
		q /= 2;
	}

	// step 2: Select a z which is a quadratic non resudue modulo p
	z = 1;
	while (jacobi(z, p) != -1) {
		z += 1;
	}

	/* step 3 */
	m = s;
	c = modPower(z, q, p);
	t = modPower(a, q, p);
	R = modPower(a, (q + 1) / 2, p);

	// Step 4: Search for a solution
	while (t != 1) {
		if (t == 0) {
			R = 0;
			break;
		}
		// Find the lowest i such that t ^ (2 ^ i) = 1 (mod p)
		i = 0;
		e = 1;
		for (i = 0; i < m; i++) {
			if (modPower(t, e, p) == 1)
				break;
			e *= 2;  /* e = 2^i */
		}

		// Update next value to iterate
		long long temp = m - i - 1;
		assert(temp < 64 && temp >= 0);  /* pow2 exponent must be within this range */
		b = modPower(c, pow2((unsigned int)temp), p);

		/* NB multiplication below could overflow if the intermediate product
		exceeds 64 bits. modMult avoids this risk */
		R = modMult(R, b, p);		           //R = (R * b) % p;
		t = modMult(modMult(t, b, p), b, p);   // t = (t * b * b) % p; 
		 /* NB apply modulus to intermediate result to avoid risk of overflow */
		c = modMult(b, b, p);                  // c = (b * b) % p;
		m = i;
	}

#ifdef _DEBUG
	assert(modMult(R, R, p) == a);
#endif

	result.push_back(R);
	result.push_back(p - R);
	return result;
}

std::vector <Znum> primeModSqrt(const Znum &aa, const Znum &p) {
	std::vector <Znum> result;
	Znum q, z, m, i, e, a;
	Znum c, t, R, b;
	long long s;

	a = aa % p;
	if (a < 0)
		a += p;   /* normalise a so it's in range 0 to p-1 */

	// Simple case
	if (a == 0) {
		result.push_back(0);
		return result;
	}

	if (p == 2) {
		result.push_back(a);  // a is 0 or 1
		return result;
	}

	/* Check solution existence on odd prime. Because p is prime the Jacobi
	symbol is the same as the Legendre symbol. */
	if (jacobi(a, p) != 1)
		return result;    // empty list; no solutions
#ifdef _DEBUG
	/* recheck existence of solution */
	{
		Znum x1, p2;
		p2 = (p - 1) / 2;
		mpz_powm(ZT(x1), ZT(a), ZT(p2), ZT(p));
		assert(x1 == 1);
    }
#endif

	// Simple case
	if (p % 4 == 3) {
		R = modPower(a, (p + 1) / 4, p);
		result.push_back(R);
		result.push_back(p - R);
		return result;
	}

	// step 1: Factor p - 1 of the form q * 2 ^ s(with Q odd)
	q = p - 1;
	s = 0;
	while (ZisEven(q)) {
		s += 1;
		q /= 2;
	}

	// step 2: Select a z which is a quadratic non resudue modulo p
	z = 1;
	while (jacobi(z, p) != -1) {
		z += 1;
	}

	/* step 3 */
	m = s;
	c = modPower(z, q, p);
	t = modPower(a, q, p);
	R = modPower(a, (q + 1) / 2, p);

	// Step 4: Search for a solution
	while (t != 1) {
		if (t == 0) {
			R = 0;
			break;
		}
		// Find the lowest i such that t ^ (2 ^ i) = 1 (mod p)
		i = 0;
		e = 1;
		for (i = 0; i < m; i++) {
			if (modPower(t, e, p) == 1)
				break;
			e *= 2;  /* e = 2^i */
		}

		// Update next value to iterate
		long long temp = MulPrToLong (m - i - 1);  /* could get overflow? */

		b = modPower(c, pow2bi(temp), p); /* b = c^(2^(m-i-1)) mod p */

		R = (R * b) % p;
		t = (t * b * b) % p; 
		c = (b * b) % p;
		m = i;
	}

#ifdef _DEBUG
	assert((R *R%p) == a);
#endif

	result.push_back(R);
	result.push_back(p - R);
	return result;
}

/* Solve the equation given a and m.
	x^2 ≡ a mod m
to find the moduluar square root modulo m where m is not prime we need to 
find the square root modulo each prime factor p1, p2 ... of m. We can combine 
one root for each prime factor using the Chinese remainder theorem (CRT). 
If m has x unique prime factors there are generally 2^x combinations 
giving 2^x roots. The idea to use CRT came from a Stack Overflow answer */
std::vector <long long> ModSqrt(const long long aa, const unsigned long long m) {
	factorsS pFactors;
	int numFactors;
	long long cMod = 1;
	long long p, a; 
	int e;
	std::vector <long long> pRoots, cRoots, cRoots2;

	a = aa % m;
	if (a < 0)
		a += m;  /* normalise a so it's in range 0 to m-1 */

	if ( a != 0 && gcd(a, m) != 1)
		return cRoots;   /* return empty list; no solutions */

	generatePrimes(std::max(llSqrt(m), 393203ULL));
	numFactors = primeFactors(m, pFactors);

	/* get roots for each prime factor of modulus separately */
	for (int i = 0; i < numFactors; i++) {
		cRoots2.clear();
		p = pFactors.factorlist[i][0];
		e = (int)pFactors.factorlist[i][1];
		if (e != 1) {
			pRoots = ModSqrt2(a, p, e);   /*find roots such that r^2 = a mod p^e */
		}
		else 
			pRoots = primeModSqrt(a, p);   /*find roots such that r^2 = a mod p */

		if (pRoots.empty()) {
			cRoots.clear();   /* there are no solutions */
			return cRoots;
		}

		if (cRoots.empty()) {
			cMod = power(p, e);    /* 1st time round, initialise combined modulus*/
			cRoots = pRoots;      /* and combined roots */
		}
		else {
			/* generate a new set of combined roots, by taking each root of the 
			   already-combined set and combining it separately with each new root, so 
			   the combined set doubles in size each time round. */
			for (auto r1: cRoots)
				for (auto r2 : pRoots) {
					cRoots2.push_back(ChineseRem(r1, cMod, r2, power(p, e)));
				}

			cMod *= power(p, e);    /* update combined modulus & roots */
			cRoots.swap(cRoots2);
		}
	}

	/* sort roots into ascending order */
	std::sort(cRoots.begin(), cRoots.end());
#ifdef _DEBUG
	for (long long r : cRoots) {
		long long x1 = r * r%m;
		if (x1 != a)
			std::cerr << "Invalid root a = " << a << " m = " << m << " root = " << r << '\n';
	}
#endif
	return cRoots;
}

std::vector <Znum> ModSqrt(const Znum &aa, const Znum &m) {
	fList pFactors;
	Znum cMod = 1;
	Znum p, a, r;
	int e;
	std::vector <Znum> pRoots, cRoots, cRoots2;

	a = aa % m;
	if (a < 0)
		a += m;  /* normalise a so it's in range 0 to m-1 */

	if (a != 0 && gcd(a, m) != 1)
		return cRoots;   /* return empty list; no solutions */

	auto rv = factorise(m, pFactors, nullptr);
	assert(rv);    /* check factorisation worked */

	/* get roots for each prime factor of modulus separately */
	for (int i = 0; i < pFactors.fsize(); i++) {
		cRoots2.clear();
		p = pFactors.f[i].Factor;
		e = pFactors.f[i].exponent;
		if (e != 1) {
			pRoots = ModSqrt2(a, p, e);   /*find roots such that r^2 = a mod p^e */
		}
		else
			pRoots = primeModSqrt(a, p);   /*find roots such that r^2 = a mod p */

		if (pRoots.empty()) {
			cRoots.clear();   /* there are no solutions */
			return cRoots;
		}

		if (cRoots.empty()) {
			cMod = power(p, e);    /* 1st time round, initialise combined modulus*/
			cRoots = pRoots;      /* and combined roots */
		}
		else {
			/* generate a new set of combined roots, by taking each root of the
			   already-combined set and combining it separately with each new root, so
			   the combined set doubles in size each time round. */
			for (auto r1 : cRoots)
				for (auto r2 : pRoots) {
					ChineseRem(r1, cMod, r2, power(p, e), r);
					cRoots2.push_back(r);
				}

			cMod *= power(p, e);    /* update combined modulus & roots */
			cRoots.swap(cRoots2);
		}
	}

	/* sort roots into ascending order */
	std::sort(cRoots.begin(), cRoots.end());
#ifdef _DEBUG
	for (Znum r : cRoots) {
		Znum x1 = r * r%m;
		if (x1 != a)
			std::cerr << "Invalid root a = " << a << " m = " << m << " root = " << r << '\n';
	}
#endif
	return cRoots;
}