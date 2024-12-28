#include "pch.h"

/* forward declaration */
int64_t ChineseRem(int64_t a1, int64_t n1, int64_t a2, int64_t n2);

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

	  a*x + b*y = gcd (a, b )
*/
int64_t extendedGcd(const int64_t a, const int64_t b, int64_t* x, int64_t* y) {
	int64_t s = 0, old_s = 1;;
	int64_t t = 1, old_t = 0;
	int64_t r = b, old_r = a;
	int64_t quotient, temp, qs, qt;
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

/*
The Chinese remainder theorem is a theorem of number theory, which states
that, if one knows the remainders of the division of an integer n by several
integers, then one can determine uniquely the remainder of the division of n
by the product of these integers, under the condition that the divisors are
pairwise coprime.
I.E:
Let n1, ..., nk be integers greater than 1, which are often called moduli or
divisors. Let us denote by N the product of the ni.

The Chinese remainder theorem asserts that if the ni are pairwise coprime, and
if a1, ..., ak are integers such that 0 ≤ ai < ni for every i, then there is
one and only one integer x, such that 0 ≤ x < N and the remainder of the Euclidean
division of x by ni is ai for every i.

If the moduli are NOT coprime there is a solution only if each a is a multiple of the
gcd of the moduli
If there is no solution return 0 
generatePrimes must have been called first
*/

int64_t ExtChineseRem(int64_t a1, uint64_t n1, int64_t a2, uint64_t n2) {
	/* find x such that x ≡ a1 (mod n1) and x ≡ a2 (mod n2)
	x will be in the range 0 to n1*n2
	This extended version handles cases where n1 and n2 are not coprime.
	n1 and n2 are factorised, the 2 original congruences are expanded into a set of
	congruences, one for each prime factor and the set is then solved if possible.
	if n1 and n2 are not coprime they will have 1 or more common prime factors.
	the congruence arising from n1 must be compatible with the congruence from n2
	for the same prime factor, otherwise there is no solution.

	If there is no solution the value zero is returned. */

	factorsS n1F, n2F;
	int64_t  n1fc;  /* number of prime factors of n1*/
	int64_t  n2fc;  /* number of prime factors of n2*/
	uint64_t cfactors[34][3] = { 0 }, minexp, x;  /* combined list based on
															factors of n1 and n2 */
	uint64_t A1, A2, N1, N2;
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
			minexp = min(n2F.factorlist[n2x][1], n1F.factorlist[n1x][1]);
			minexp = power(n1F.factorlist[n1x][0], (unsigned int)minexp);
			if (a1 % minexp != a2 % minexp) {
				/* factors cannot be combined */
				printf_s("ext C.R. a1=%lld, n1=%lld, a2=%lld, n2=%lld\n", a1, n1, a2, n2);
				printf_s("prime is %lld ", n2F.factorlist[n2x][0]);
				printf_s(" exponents: %lld & %lld\n",
					n1F.factorlist[n1x][1], n2F.factorlist[n2x][1]);
				printf_s("a1 mod %lld = %lld   a2 mod %lld = %lld\n", minexp, a1 % minexp, minexp,
					a2 % minexp);
				return 0;  // there is no solution!
			}

			/* combine the factors */
			cfactors[cx][0] = n1F.factorlist[n1x][0];  // copy prime
			cfactors[cx][1] = max(n1F.factorlist[n1x][1], n2F.factorlist[n2x][1]); // use max exponent
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
	assert(x % n1 == a1);   // check for bugs
	assert(x % n2 == a2);
	return x;
}


/* find x such that x ≡ a1 (mod n1) and x ≡ a2 (mod n2)
	x will be in the range 0 to n1*n2

We use the extended Euclidian algorithm to find integers m1 and m2 such that
m1*n1 + m2*n2 = gcd(n1,n2)
A solution is given by x = a1*m2*n2 + a2*m1*n1, provided that n1 and n2 are co-prime

if n1 and n2 are not co-prime,
x = (a1*m2*n2 + a2*m1*n1)/gcd(n1,n2)   ONLY if both a1 and a2 are multiples of the gcd,
otherwise there may be no solution.
If there is no solution return 0 
*/
int64_t ChineseRem(int64_t a1, int64_t n1, int64_t a2, int64_t n2) {
	int64_t m1, m2, llgcd, llgcd2 = 0, x = 0;
	static Znum bt1, bt2, bx, bm1, bm2, t3;

	llgcd = extendedGcd(n1, n2, &m1, &m2);
	if (llgcd != 1) {
		if ((a1 % llgcd != 0) || (a2 % llgcd != 0)) {
			if (a1 % llgcd != a2 % llgcd)
				return 0;   // there is no solution
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
	bx = (Znum)a1 * m2 * n2 + a2 * m1 * n1;

	// use bigints to avoid overflow in intermediate products

	// ensure value of x is in correct range, 0 to n1*n2-1
	t3 = abs((Znum)(n1 * n2));
	if (bx < 0) {
		mpz_fdiv_r(ZT(bx), ZT(bx), ZT(t3));
	}
	if (bx >= t3)
		mpz_fdiv_r(ZT(bx), ZT(bx), ZT(t3));

	if (llgcd != 1) {
		bx *= llgcd;
	}
	x = ZnumToLong(bx);         // convert result back to a normal integer

	return x;
}

/* find g such that g ≡ a (mod m) and g ≡ b (mod n)
g will be in the range 0 to m*n/gcd(m,n) -1. 
We use the extended Euclidian algorithm to find integers u and v such that
u*m + v*n = gcd(m,n)
A solution is given by g = a*v*n + b*u*m, provided m and n are co-prime
If there is no solution the value returned is -3.
*/
void ChineseRem(const Znum& a, const Znum& m, const Znum& b, const Znum& n, Znum& g) {

	Znum u, v, gcd, t2;
	/* set gcd to gcd(m, n), also m*u + n*v = gcd */
	mpz_gcdext(ZT(gcd), ZT(u), ZT(v), ZT(m), ZT(n));
	if (gcd != 1) {
		if ((a - b) % gcd != 0) {
			if (verbose > 0) {
				char buf[4000];  /* guess how big buffer should be; if it's too small the
								 message will be truncated */
				gmp_snprintf(buf, sizeof(buf), "Chinese Rem: no solution for %Zd, %Zd, %Zd, %Zd \n",
					a, m, b, n);
				std::cout << buf;
				g = -3;   /* there is no solution */
			}
			return;
		}
	}
	g = a * v * n + b * u * m;
	t2 = abs(m * n);
	if (gcd > 1) {
		g /= gcd;
		t2 /= gcd;  /* t2 = lcm(m, n) */
	}
	if (g < 0) {    // if g < 0 add a multiple of lcm(m*n)
		mpz_fdiv_r(ZT(g), ZT(g), ZT(t2));  // ensure x is in required range
	}
	if (g >= t2) {  // if g > m*n subtract a multiple of lcm(m*n)
		mpz_fdiv_r(ZT(g), ZT(g), ZT(t2));
	}
	return;
}

/* find g such that g ≡ p[0] (mod p[1], g ≡ p[2] (mod[p3]) etc.
Return -ve value for any error.
-1 if p[0] or p[2] or p[4] etc < 0
-2 if p[1] or p[3] or p[5] etc < 1 (modulus must be > 1) 
-3 if there is no solution. */
void ChineseRemV(const std::vector <Znum>& p, Znum& result) {
	Znum modulus = p[1];
	result = p[0];
	if (p[0] < 0) {
		result = -1;  /* invalid parameter value */
		return;
	}
	if (p[1] <= 1) {
		result = -2;   /* modulus must be > 1 */
		return;
	}
		
	for (unsigned int ix = 2; ix < p.size(); ix += 2) {
		if (p[ix] < 0) {
			result = -1;
			return;
		}
		if (p[ix+1] <= 1) {
			result = -2;    /* modulus must be > 1 */
			return;
		}
		ChineseRem(result, modulus, p[ix], p[ix + 1], result);
		if (result < 0)
			return;  /* there is no solution */
		modulus = lcm(modulus, p[ix + 1]);  /* get new modulus */
		if (verbose > 1)
			gmp_printf("Chinese rem. ix = %d, result %Zd, modulus = %Zd \n",
				ix, ZT(result), ZT(modulus));
	}
}
