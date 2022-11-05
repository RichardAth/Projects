#include "pch.h"
#include "bignbr.h"
#include "bigint.h"
#include "factor.h"
void textError(retCode rc);

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

/* find least significant 1-bit in num, equivalent to counting 0-bits, starting from 
least significant bit. If num is a power of 2, then 2^return value = num.*/
long long countZeroBits(const Znum &num) {
	mp_bitcnt_t r = mpz_scan1(ZT(num), 0);
	return r;
}

/* calculate max power of p that is a divisor of num */
long long extract(const Znum& num, const Znum p) {
	Znum residue;
	mp_bitcnt_t r = mpz_remove(ZT(residue), ZT(num), ZT(p));
	return r;
}

/* c is an even power of 2, prime = 2, c < 2^lambda */
std::vector<Znum>ModSqrtp2x(const Znum& c, const Znum& prime, const int lambda) {
	Znum increment, r1, sqrtc;
	Znum mod = power(prime, lambda);
	std::vector <Znum> roots;
	//size_t size1, size2;
	long long k = countZeroBits(c);
	assert((k & 1) == 0);  /* k must be even */
	sqrtc = power((Znum)2, k / 2);

	if ((3+k) <= lambda)
		increment = mod / (2*sqrtc);  /* increment = 2^(lambda-1-k/2)*/
	else
		increment = sqrtc*4;          /* increment = 2^(2+k/2)*/

	/* the roots are 2 arithmetic progressions */
	r1 = sqrtc;

	for (r1 = sqrtc; r1 < mod; r1 += increment) {
		roots.push_back(r1);
		roots.push_back(mod - r1);
	}
	std::sort(roots.begin(), roots.end());
	/*size1 = roots.size();
	auto last = std::unique(roots.begin(), roots.end());
	roots.erase(last, roots.end());
	size2 = roots.size();
	if (size2 != size1 && verbose > 0) {
		std::cout << "modsqrt(" << c << ", " << mod << ") discarded " << size1 - size2
			<< " duplicate roots \n" << "retained "<< size2 << " roots \n";

	}*/
	return roots;
}


/* special for prime =2, c is odd*/
std::vector<Znum>ModSqrtp2(const Znum& c, const Znum& prime, const int lambda) {
	Znum r1, r2, root, x, x2, gcdv, increment, sqrtc;
	Znum mod = power(prime, lambda);
	std::vector <Znum> roots;
	size_t size1, size2;
	assert(prime == 2);

	if (c == 1) {
		roots.push_back(1);
		roots.push_back(mod - 1);  /* get 2nd root */
		if (mod > 4) {
			roots.push_back(mod / 2 - 1);
			roots.push_back(mod / 2 + 1);
		}
	}
	else if ((lambda >= 3) && (c % 8 == 1)) {
		x2 = 1;
		for (int k = 3; k <= lambda; k++) {
			Znum i = ((x2 * x2 - c) / pow2(k)) % 2;
			if (i == -1) i = 1;
			x2 = x2 + i * pow2(k - 1);
		}
		while (x2 < 0)
			x2 += mod;
		if (x2 >= mod)
			x2 %= mod;   /* ensure x2 is in range 0 to mod-1 */

		/* x2 is 1 of 4 roots. We can easily generate the other 3 given any 1 
		of the 4 roots. */
		roots.push_back(x2);
		roots.push_back(mod - x2);
		if (x2 < mod / 2) {
			roots.push_back(mod / 2 + x2);
			roots.push_back(mod / 2 - x2);
		}
		else {
			roots.push_back(x2 - mod / 2);
			roots.push_back(mod * 3 / 2 - x2);
		}
	}

	  /* remove any duplicates */
	std::sort(roots.begin(), roots.end());
	size1 = roots.size();
	auto last = std::unique(roots.begin(), roots.end());
	roots.erase(last, roots.end());
	size2 = roots.size();
	if (size2 != size1 && verbose > 0) {
		std::cout << "modsqrt(" << c << ", " << mod << ") discarded " << size1 - size2
			<< " duplicate roots \n" << "retained " << size2 << " roots \n";

	}
	return roots;
}

/* c is a perfect square and an even power of of prime, prime > 2( */
std::vector <Znum> ModSqrt2xs(const Znum& c, const Znum& prime, const int lambda) {
	Znum mod = power(prime, lambda);
	Znum r1, r2, sqrtc, increment;
	std::vector <Znum> roots;
	//size_t size1, size2;

	long long k = extract(c, prime);   /* find k such that prime^k = c */
	assert((k & 1) == 0);   /* k must be even */

	/* the roots are 2 arithmetic progressions */
	sqrtc = power(prime, k/2);
	increment = mod / sqrtc;  /* increment = prime^(lambda - k/2) */
	for (r1 = sqrtc; r1 < mod; r1 += increment) {
		roots.push_back(r1);
		roots.push_back(mod - r1);
	}
	std::sort(roots.begin(), roots.end());
	//size1 = roots.size();
	//auto last = std::unique(roots.begin(), roots.end());
	//roots.erase(last, roots.end());
	//size2 = roots.size();
	//if (size2 != size1 && verbose > 0) {
	//	std::cout << "modsqrt(" << c << ", " << mod << ") discarded " << size1 - size2
	//		<< " duplicate roots \n";
	//}
	return roots;
}

/* c is not a multiple of prime, prime != 2 */
std::vector <Znum> ModSqrt2x(const Znum &c, const Znum &prime, const int lambda) {
	std::vector <Znum> roots;
	Znum r1, r2, root, x;
	Znum mod = power(prime, lambda);

	std::vector <Znum> rx;

	/* this part was derived from Wikipedia
	see https://en.wikipedia.org/wiki/Tonelli%E2%80%93Shanks_algorithm#Tonelli's_algorithm_will_work_on_mod_p^k
	The explanation there is not very clear but I got working code out of it. */
	rx = primeModSqrt(c, prime);   /* rx^2 ≡ c (mod prime) */
	if (rx.empty()) {
		roots.clear();
		return roots;     /* there are no solutions  */
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
	if ((root * root) % mod == c) {
		roots.push_back(root);
		roots.push_back(mod - root);    /* get 2nd root */
	}
	else {
		if (verbose > 1 || root > 0) {
			std::cerr << "discarded invalid root: c = " << c << " mod = " << mod << " root = " << root << '\n';
		}
	}
	return roots;
}

/* find x such that x^2 ≡ c mod prime^lambda
  the method used was partly derived from Wikipedia but the cases where c=0,
  and where prime=2 are all my own work. I have not seen descriptions that cover
  these cases let alone working code that does what this does ANYWHERE and,
  believe me, I looked. */
std::vector <Znum> ModSqrt2(const Znum& cc, const Znum& prime, const int lambda) {
	Znum mod = power(prime, lambda);
	Znum c, r1, r2, gcdv, sqrtgcd;
	std::vector <Znum> roots;

	c = cc % mod;
	if (c < 0)
		c += mod;  /* ensure c is in range 0 to mod-1 */

	/* treat c=0 as special case. 0 is a solution, but there may be others. */
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

	gcdv = gcd(c, mod);
	if (gcdv > 1) {
		if (isPerfectSquare(gcdv, sqrtgcd)) {
			if (prime == 2)
				roots = ModSqrtp2(c/gcdv, prime, lambda);
			else
				roots = ModSqrt2x(c / gcdv, prime, lambda);
			if (roots.empty())
				return roots;  /* no solutions*/
			r1 = roots[0];  /* we only need 1 of the roots just obtained */

			if (prime == 2)
				roots = ModSqrtp2x(gcdv, prime, lambda);
			else
				roots = ModSqrt2xs(gcdv, prime, lambda);

			for (size_t i = 0; i < roots.size(); i++)
				roots[i] = modMult(roots[i], r1, mod);
			return roots;
		}
		else return roots;  /* no solutions*/
	}
	else  /* c and mod are mutually prime*/
		if (prime == 2)
			return ModSqrtp2(c, prime, lambda);
		else
			return(ModSqrt2x(c, prime, lambda));
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
giving 2^x roots. The idea to use CRT came from a Stack Overflow answer 
 */
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
  faster more sophisticated method. N.B. very slow for larger numbers
  e.g a 30 bit modulus takes about 1 second. Each extra bit doubles the time. */
static std::vector<long long> ModSqrtBF(long long a, long long m) {
	std::vector <long long> roots;
	a = a % m;
	if (a < 0)
		a += m;  /* normalise a so it's in range 0 to m-1 */
	for (long long r = 0; r <= m/2; r++) {
		if (r * r % m == a) {
			roots.push_back(r);
			if (r > 0 &&(r*2 != m))
				roots.push_back(m - r);  /* usually, m-r is also a root */
		}
	}
	std::sort(roots.begin(), roots.end());
	return roots;
}

void test9timerx(int type, int p2d, int p3d) {
	gmp_randstate_t state;
	Znum x, m;
	std::vector <Znum> r;
	std::vector <long long> rl;
	long long rCount = 0, nrCount = 0;
	char msg[3][30] = { "Standard modsqrt: time used: ",
						"Brute force: time used: ",
						"QMES: time used: " };

	gmp_randinit_mt(state);  // use Mersenne Twister to generate pseudo-random numbers
	gmp_randseed_ui(state, 756128234);
	/* fixed seed means that the exact same tests can be repeated provided
	   that the same size of number is used each time */

	auto start = clock();	// used to measure execution time
	for (int i = 1; i <= p3d; i++) {
		mpz_urandomb(ZT(x), state, p2d);  // get random number, size=p2 bits
		mpz_urandomb(ZT(m), state, p2d);  // get random number, size=p2 bits
		switch (type) {
		case 0:
			r = ModSqrt(x, m);
			if (r.empty())
				nrCount++;
			else
				rCount += r.size();
			break;

		case 1:
			if (p2d < 64) {
				long long xl = MulPrToLong(x);
				long long ml = MulPrToLong(m);
				rl = ModSqrtBF(xl, ml);
				if (rl.empty())
					nrCount++;
				else
					rCount += rl.size();
			}
			break;
		case 2:
			r = ModSqrtQE(x, m);
			if (r.empty())
				nrCount++;
			else
				rCount += r.size();
			break;
		}
	}
	gmp_randclear(state);  // clear state - avoid memory leakage
	std::cout << "No roots = " << nrCount << " found " << rCount << " roots \n";
	auto end = clock();   // measure amount of time used
	double elapsed = (double)end - start;
	PrintTimeUsed(elapsed, msg[type]);
 }

void test9timer(const std::vector <std::string>& p) {
	int p2d = -1, p3d = -1;

	if (p.size() >= 2) {
		p2d = atoi(p[1].c_str());  /* convert to binary */
	}
	if (p2d < 10) {
		std::cout << "Use default 10 for number size in bits \n";
		p2d = 10;
	}
	if (p.size() >= 3) 
		p3d = atoi(p[2].c_str());  /* convert to binary */
	if (p3d < 5) {
		std::cout << "Use default 5 for number of tests \n";
		p3d = 5;
	}

	test9timerx(0, p2d, p3d);       /* use standard modsqrt*/
	if (p2d <= 30)
		test9timerx(1, p2d, p3d);  /* use brute force modsqrt */
	test9timerx(2, p2d, p3d);      /* use QMES modsqrt */
}

/* calculate modular square roots using 2 different methods & compare results.
   Return true if results match, otherwise false. */
static bool test9once(long long a, long long m, std::vector <long long> &r2,
	bool newb) {
	std::vector <Znum> r3;
	Znum az = a;   /* change type from 64-bit to extended precision */
	Znum mz = m;   /* change type from 64-bit to extended precision */
	bool error = false;
	r2.clear();

	r2 = ModSqrtBF(a, m);  /* brute force method */
	if (newb)
		r3 = ModSqrtQE(az, mz);
	else
		r3 = ModSqrt(az, mz);  
	if (r3.size() != r2.size())
		error = true;
	else {
		for (int i = 0; i < r2.size(); i++) {
			if (r2[i] != r3[i])
				error = true;
		}
	}

	if (error) {
		printf_s("Modsqrt(%lld, %lld), Znum results don't match! \nGot: ", a, m);
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
		if (r2.empty()) {
			if (verbose > 1)
				printf_s("modsqrt(%lld, %lld) has no roots \n", a, m);
		}
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
void doTests9(const std::string& porig) {
	std::string params = porig;
	bool rv = true;
	std::vector<long long> roots;
	std::vector<Znum> rootsZ;
	bool newb = false;
	std::vector<std::string> p;
	char seps[] = " ,\n";
	char* token = nullptr;
	char* next = nullptr;

	token = strtok_s(&params[0], seps, &next);
	while (token != nullptr) {
		p.push_back(token);
		token = strtok_s(nullptr, seps, &next);
	}

	//auto numParams = sscanf_s(params.data(), "%s,%s,%s",
		//p1, sizeof(p1), p2, sizeof(p2), p3, sizeof(p3));
	if (p.size() >= 1) {
		int timec = _strnicmp(p[0].c_str(), "time", 4);
		if (timec == 0) {
			test9timer(p);
			return;
		}

		int cmp = _strnicmp(p[0].c_str(), "new", 3);
		newb = (cmp == 0);
	}
	auto start = clock();	// used to measure execution time
	for (long long m = 2; m <= 2000; m++) {
		for (long long a = 0; a < m; a++) {
			rv &= test9once(a, m, roots, newb);
		}
		if (m % 100 == 0)
			std::cout << m + 1 << " tests completed \r";
	}

	//rv &= test9once(99, 107, roots, newb);      // roots are 45 and 62
	rv &= test9once(2191, 12167, roots, newb);  // roots are 1115, 11052
	rv &= test9once(4142, 24389, roots, newb);  // roots are 2333, 22056
	//rv &= test9once(3, 143, roots, newb);       // roots are 17, 61, 82, 126
	rv &= test9once(11, 2 * 5 * 7 * 19, roots, newb); // roots are 121, 159 411, 639, 
	                                      // 691, 919, 1171, 1209
	//rv &= test9once(0, 176, roots, newb);       // roots are zero, 44, 88, 132
	rv &= test9once(1, 121550625, roots, newb); // roots are 1, 15491251, 51021251, 
								          // 55038124, 66512501, 70529374, 
	                                      // 106059374, 121550624,
	rv &= test9once(0, 121550625, roots, newb); // 11025 different roots!
	
	if (rv)
		std::cout << "All modular square root tests completed successfully. \n";
	else
		std::cout << "One or more modular square root tests failed. \n";
	auto end = clock();   // measure amount of time used
	double elapsed = (double)end - start;
	PrintTimeUsed(elapsed, "All tests completed. Time used = ");
}
