#include "pch.h"
#include "factor.h"

static std::vector <Znum> primeModSqrt(const Znum& aa, const Znum& p);

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

/* find x such that x^2 ≡ c mod prime^lambda 
  the method used was partly derived from Wikipedia but the cases where c=0, 
  and where prime=2 are all my own work. I have not seen descriptions that cover 
  these cases let alone working code that does what this does ANYWHERE and, 
  believe me, I looked. */
std::vector <Znum> ModSqrt2(const Znum &cc, const Znum &prime, const int lambda) {
	std::vector <Znum> roots;
	Znum r1, r2, root, x, c;
	Znum mod = power(prime, lambda);

	std::vector <Znum> rx;

	c  = cc%mod;
	if (c < 0)
		c += mod;  /* ensure c is in range 0 to mod-1 */

	/* treat c=0 as special case */
	if (c == 0) {
		r1 = power(prime, (lambda + 1) / 2);  /* smallest non-zero root is the
										   smallest power of prime >= sqrt(mod) */
		r2 = 0;
		while (r2 < mod) {
			roots.push_back(r2);
			r2 += r1;  /* the roots form an arithmetic progression: 0, r1, 2*r1, 3*r1 etc */
		}
		return roots;
	}

	/* must treat prime=2 as a special case. */
	if (prime == 2 && c != 0) {
		/* there are no roots unless c = 1 */
		if (c% mod == 1) {
			roots.push_back(1);
			roots.push_back(mod - 1);  /* get 2nd root */
		}
		return roots;
	}

	/* this part was derived from Wikipedia
	see https://en.wikipedia.org/wiki/Tonelli%E2%80%93Shanks_algorithm#Tonelli's_algorithm_will_work_on_mod_p^k
	The explanation there is not very clear but I got working code out of it */
	rx = primeModSqrt(c, prime);   /* rx^2 ≡ c (mod prime) */
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
	r1 = modPower(x, power(prime, lambda - 1), mod);
	r2 = modPower(c, (power(prime, lambda) - 2 * power(prime, lambda - 1) + 1) / 2, mod);
	root = modMult(r1, r2, mod);    /* get final product */
	roots.push_back(root);
	roots.push_back(mod - root);    /* get 2nd root */

	return roots;
}

/*
Square root modulo prime number using Tonelli–Shanks algorithm
Solve the equation given a and prime.
	x^2 ≡ a mod prime
and return list of solutions. There will be either 0, 1 or 2 solutions
see https://en.wikipedia.org/wiki/Tonelli%E2%80%93Shanks_algorithm
(renamed the solution variable n to a)
*/
std::vector <Znum> primeModSqrt(const Znum &aa, const Znum &prime) {
	std::vector <Znum> result;
	Znum q, z, e, a;
	Znum c, t, R, b;
	long long i, m, s;

	a = aa % prime;
	if (a < 0)
		a += prime;   /* normalise a so it's in range 0 to prime-1 */

	// Simple case
	if (a == 0) {
		result.push_back(0);
		return result;
	}

	if (prime == 2) {
		result.push_back(a);  // a is 0 or 1
		return result;
	}

	/* Check solution existence on odd prime. Because prime is prime the Jacobi
	symbol is the same as the Legendre symbol. */
	if (jacobi(a, prime) != 1)
		return result;    // empty list; no solutions
#ifdef _DEBUG
	/* recheck existence of solution */
	{
		/* check that a^((prime-1)/2) ≡ 1 (mod prime)*/
		Znum x1, p2;
		p2 = (prime - 1) / 2;
		mpz_powm(ZT(x1), ZT(a), ZT(p2), ZT(prime));
		assert(x1 == 1);
    }
#endif

	// Simple case
	if (prime % 4 == 3) {
		R = modPower(a, (prime + 1) / 4, prime);
		result.push_back(R);
		result.push_back(prime - R);
		return result;
	}

	// Tonelli-Shanks step 1: Factor prime - 1 of the form q * 2 ^ s(with Q odd)
	q = prime - 1;
	s = 0;
	while (isEven(q)) {
		s += 1;
		q /= 2;
	}

	// step 2: Select a z which is a quadratic non resudue modulo prime
	z = 1;
	while (jacobi(z, prime) != -1) {
		z += 1;
	}

	/* step 3 */
	m = s;
	c = modPower(z, q, prime);
	t = modPower(a, q, prime);
	R = modPower(a, (q + 1) / 2, prime);

	// Step 4: Search for a solution
	while (t != 1) {
		if (t == 0) {
			R = 0;
			break;
		}
		// Find the lowest i such that t ^ (2 ^ i) = 1 (mod prime)
		i = 0;
		e = 1;
		for (i = 0; i < m; i++) {
			if (modPower(t, e, prime) == 1)
				break;
			e *= 2;  /* e = 2^i */
		}

		// Update next value to iterate
		long long temp = m - i - 1;  

		b = modPower(c, pow2bi(temp), prime); /* b = c^(2^(m-i-1)) mod prime */

		R = (R * b) % prime;
		t = (t * b * b) % prime; 
		c = (b * b) % prime;
		m = i;           /*  i < m, so m decreases each time round */
	}

#ifdef _DEBUG
	assert((R *R%prime) == a);
#endif

	result.push_back(R);
	result.push_back(prime - R);
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

	/* below is a necessary condition, but satisfying this test does not guarantee 
	  that there are any roots. */
	if ((m % 4 != 2) && (mpz_kronecker(ZT(a), ZT(m)) == -1))
		return cRoots;   /* return empty list; no solutions */
	Znum m2 = m / 2;
	if ((m % 4 == 2) && (mpz_kronecker(ZT(a), ZT(m2)) == -1))
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

/* find modular square roots by brute force. Incredibly simple compared to the
  faster more sophisticated method. */
static std::vector<long long> ModSqrtBF(long long a, long long m) {
	std::vector <long long> roots;
	a = a % m;
	if (a < 0)
		a += m;  /* normalise a so it's in range 0 to m-1 */
	for (long long r = 0; r < m; r++) {
		if (r*r%m == a)
			roots.push_back(r);
	}
	return roots;
}

/* calculate modular square roots using 2 different methods & compare results.
   Return true if results match, otherwise false. */
static bool test9once(long long a, long long m) {
	std::vector <long long> r2;
	std::vector <Znum> r3;
	Znum az = a;   /* change type from 64-bit to extended precision */
	Znum mz = m;   /* change type from 64-bit to extended precision */
	bool error = false;

	r2 = ModSqrtBF(a, m);  /* brute force method */
	r3 = ModSqrt(az, mz);  /* sophisticated (faster) method */
	if (r3.size() != r2.size())
		error = true;
	else {
		for (int i = 0; i < r2.size(); i++) {
			if (r2[i] != r3[i])
				error = true;
		}
	}

	if (error) {
		printf_s("a = %lld, m = %lld, Znum results don't match! \nGot: ", a, m);
		for (auto r : r3) {
			gmp_printf("%Zd, ", r);
		}
		if (r2.empty())
			printf_s("\nActually no solutions\n");
		else {
			printf_s("\n should be: ");
			for (auto r : r2) {
				printf_s("%lld, ", r);
			}
			putchar('\n');
		}
		return false;
	}

	if (verbose > 0) {
		if (r2.empty())
			printf_s("modsqrt(%lld, %lld) has no roots \n", a, m);
		else {
			printf_s("modsqrt(%lld, %lld) = ", a, m);
			for (long long r : r2)
				printf_s("%lld, ", r);
			putchar('\n');
		}
	}
	return true;
}

/* test modular square root */
void doTests9(void) {
	bool rv = true;

	rv &= test9once(99, 107);      // roots are 45 and 62
	rv &= test9once(3, 22);        // roots are 5 & 17
	rv &= test9once(99, 100);      // no roots
	rv &= test9once(2191, 12167);  // roots are 1115, 11052
	rv &= test9once(4142, 24389);  // roots are 2333, 22056
	rv &= test9once(3, 143);       // roots are 17, 61, 82, 126
	rv &= test9once(11, 2 * 5 * 7 * 19); // roots are 121, 159 411, 639, 691, 919, 1171, 1209
	rv &= test9once(9, 44);        // roots are 3, 19, 25, 41
	rv &= test9once(0, 44);        // roots are zero, 22
	rv &= test9once(0, 4);         // roots are 0, 2
	rv &= test9once(0, 9);         // roots are 0, 3, 6
	rv &= test9once(0, 8);         // roots are 0, 4
	rv &= test9once(0, 27);        // roots are 0, 9, 18
	rv &= test9once(0, 16);        // roots are 0, 4, 8, 12
	rv &= test9once(0, 32);        // roots are 0, 8, 16, 24
	rv &= test9once(0, 176);       // roots are zero, 44, 88, 132
	rv &= test9once(1, 121550625); // roots are 1, 15491251, 51021251, 55038124,
								   // 66512501, 70529374, 106059374, 121550624,
	rv &= test9once(0, 121550625); // 11025 different roots!
	rv &= test9once(8, 28);        // roots are 6, 8, 20, 22
	rv &= test9once(6, 30);        // roots are 6, 24
	rv &= test9once(42, 66);       // roots are 30, 36

	if (rv)
		std::cout << "All modular square root tests completed successfully. \n";
	else
		std::cout << "One or more modular square root tests failed. \n";
}