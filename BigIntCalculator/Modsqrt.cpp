﻿#include "pch.h"
#include "factor.h"

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
std::vector <Znum> ModSqrt2(const Znum &cc, const Znum &p, const int lambda) {
	std::vector <Znum> roots;
	Znum r1, r2, root, x, c;
	Znum mod = power(p, lambda);

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

	/* must treat p=2 as a special case. No roots unless c is 1 */
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
	/* the calculation is done in stages, applying modulus each time, to avoid
	 trying to generate enormous intermediate values */
	r1 = modPower(x, power(p, lambda - 1), mod);
	r2 = modPower(c, (power(p, lambda) - 2 * power(p, lambda - 1) + 1) / 2, mod);
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
std::vector <Znum> primeModSqrt(const Znum &aa, const Znum &p) {
	std::vector <Znum> result;
	Znum q, z, e, a;
	Znum c, t, R, b;
	long long i, m, s;

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
		/* check that a^((p-1)/2) ≡ 1 (mod p)*/
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

	// Tonelli-Shanks step 1: Factor p - 1 of the form q * 2 ^ s(with Q odd)
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
		long long temp = m - i - 1;  

		b = modPower(c, pow2bi(temp), p); /* b = c^(2^(m-i-1)) mod p */

		R = (R * b) % p;
		t = (t * b * b) % p; 
		c = (b * b) % p;
		m = i;           /*  i < m, so m decreases each time round */
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
std::vector <Znum> ModSqrt(const Znum &aa, const Znum &m) {
	fList pFactors;
	Znum cMod = 1;
	Znum p, a, r;
	int e;
	std::vector <Znum> pRoots, cRoots, cRoots2;

	a = aa % m;
	if (a < 0)
		a += m;  /* normalise a so it's in range 0 to m-1 */

	if (a != 0 && !isPerfectSquare(gcd(a, m)))
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