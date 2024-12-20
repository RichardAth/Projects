﻿/*
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
#include <set>
#include <numeric>

#undef min                 // use std::min

static long long Gamma[386];
static long long Delta[386];
static long long AurifQ[386];
static int AurifCount = 0;
static Znum Zfactor;
double originalTenthSecond;
int oldTimeElapsed;

bool isEven(const Znum &a) {
    return (mpz_odd_p(ZT(a)) == 0);  /* true iff a is even (works for -ve a as well) */
}

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
 
    if (insertBigFactor(Factors, Dsal))
        AurifCount++;
    Dsal = Csal - Nbr1; 

    if (insertBigFactor(Factors, Dsal))
        AurifCount++;
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
    auto llBase = ZnumToLong(Base);

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

/* Original number to be factored = Base^Expon - increment 
increment = 1 or -1 
see https://en.wikipedia.org/wiki/Cunningham_number 
and https://en.wikipedia.org/wiki/Cunningham_Project
*/
static void Cunningham(fList &Factors, const Znum &Base, int Expon,
    const int increment, const Znum &Original) {
    int Expon2, k;
    Znum Nbr1, Nbr2; 

    Expon2 = Expon;
    AurifCount = 0;

    while (Expon2 % 2 == 0 && increment == -1) {
        Expon2 /= 2;
        BigIntPowerIntExp(Base, Expon2, Nbr1);  /* Nbr1 = Base^Expon2*/
        Nbr1 += increment;    /* Nbr1 = Base^Expon2 - 1 */
        insertBigFactor(Factors, Nbr1);     /* if Nbr1 is a factor, insert it*/
        InsertAurifFactors(Factors, Base, Expon2, 1);
    }

    k = 1;
    while (k * k <= Expon) {
        if (Expon % k == 0) {
            if (k % 2 != 0) { /* Only for odd exponent */
                BigIntPowerIntExp(Base, Expon / k, Nbr1); /* nbr1 = base^(expon/k) */
                Nbr1 += increment;      
                Nbr2 = gcd(Nbr1, Original);   // Nbr2 <- gcd(Base^(Expon/k)+incre, original)
                if (Nbr2 != Original && Nbr2 != 1) 
                    insertBigFactor(Factors, Nbr2);
                Nbr1 = Original / Nbr2;
                if (Nbr1 != Original && Nbr1 != 1)
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
    if (verbose > 0 && AurifCount > 0)
        std::cout << "AurifCount =" << AurifCount << '\n';
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
        //std::cout << "result after sort" << '\n';
        //Factors.Xprint();
    }
}

/* Insert new factor found into factor array. */
/* assume divisor is prime. ix is index of non-prime factor which is a multiple 
of divisor. either pi is the index into the prime list of the divisor, or the  
divisor value is in div. This function is used for smaller factors found by trial
 division or Pollard Rho. */
static void insertIntFactor(fList &Factors, int pi, unsigned long long div, ptrdiff_t ix) {
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
        assert(exp > 0);  /* would abort if quot were not a multiple of divisor */
        Factors.f.resize(lastfactor + 1);  // increase size of factor list
        Factors.f[lastfactor].exponent = (int)exp * Factors.f[ix].exponent;
        Factors.f[lastfactor].Factor = divisor;
        Factors.f[lastfactor].upperBound = -1; // show that new factor is prime

        Factors.f[ix].Factor = qnew;  // adjust value of original factor
        if (pi >= 0)
            Factors.f[ix].upperBound = pi;   /* show highest prime used so far in trial division */
    }
    else {
        /* replace residue of 1 with new factor */
        Factors.f[ix].Factor = divisor;
        Factors.f[ix].exponent *= (int)exp;
        Factors.f[ix].upperBound = -1;    /* new factor is prime */
    }
    SortFactors(Factors);
    if (verbose > 1) {
        std::cout << "result after adding factor " << divisor << '\n';
        Factors.Xprint();
    }
}

/* Insert new factor found into factor array. Every current factor is checked 
against the divisor. If their gcd is not 1 the existing factor is divided by 
the gcd and a new factor equal to the gcd is created. Return false if divisor
is not a factor of any existing factor */
bool insertBigFactor(fList &Factors, const Znum &divisor) {
    auto lastfactor = Factors.f.size();
    auto ipoint = lastfactor;
    bool success = false;
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
            success = true;
        }
    }
    if (success) {
        SortFactors(Factors);
        if (verbose >= 1) {
            std::cout << "Divisor = " << divisor << " result after adding factor: \n";
            Factors.Xprint();
        }
    }
    else if (verbose > 1) {
        std::cout << "insertBigFactor: could not add factor " << divisor << '\n';
    }

    return success;
}

/* n is a pseudoprime to base b. Find some divisors. This will not work if n is
   a strong pseudoprime to base b. See 
   https://en.wikipedia.org/wiki/Miller%E2%80%93Rabin_primality_test#Variants_for_finding_factors */
static bool getfactors(const Znum& n, uint32_t b, fList& Factors) {
    /* divide (n-1) by 2 until it is odd*/
    Znum two = 2;
    Znum f, kf, nm1, modpow, zb, div1, div2;
    nm1 = n - 1;
    bool factorsfound = false;
    zb = b;

    /* check that b^(n-1) = 1 (mod n) */
    mpz_powm(ZT(modpow), ZT(zb), ZT(nm1), ZT(n));
    if (modpow != 1) {
        if (verbose > 0)
            gmp_printf("%Zd is not a pseudoprime to base %d\n", n, b);
        return false;
    }

    auto shift = mpz_remove(ZT(f), ZT(nm1), ZT(two));
    /* f *(2^shift) = n-1 */
 
    for (int k = 1; k <= (1LL << shift); k <<= 1) {
        /* calculate b^kf, mod n. k doubles each time round the loop, kf also
        doubles, and modpow is squared, (mod n)*/
        kf = k * f;
        mpz_powm(ZT(modpow), ZT(zb), ZT(kf), ZT(n));
        if (modpow != 1 && modpow != nm1) {
            div1 = gcd(modpow + 1, n);
            div2 = gcd(modpow - 1, n);
            if (div1 > 1 || div2 > 1) {
                if (verbose > 0)
                    gmp_printf("pseudoprime getfactors: %Zd and %Zd are divisors of %Zd (k = %d, s = %d)\n", 
                        div1, div2, n, k, shift);
                if (div1 > 1)
                    if (insertBigFactor(Factors, div1)) {
                        factorsfound = true;
                        Factors.ct.psp++;
                    }
                if (div2 > 1)
                    if (insertBigFactor(Factors, div2)) {
                        factorsfound = true;
                        Factors.ct.psp++;
                    }
            }
        }
        else 
            break;  /* exit loop if modpow = 1 or -1, because squaring it again 
                    will give modpow = 1*/
    }
    return factorsfound;
}


/* common code to (maybe) insert factor. If factor added, factorsFound is
* set to true & the counter is incremented.
countdown, ctr, i & ref are only used if verbose > 0, to print messages. */
static void insertCarmichaelFactor(Znum &Aux4, const Znum &p, fList& Factors, 
    bool &factorsFound, const int countdown, const int ctr, const int i, int ref) {
    Znum gcdVall;
    mpz_mod(ZT(Aux4), ZT(Aux4), ZT(p));
    gcdVall = gcd(p, Aux4);
    if ((gcdVall > 1) && (gcdVall != p))
    {          // Non-trivial factor found.
        if (insertBigFactor(Factors, gcdVall)) {
            factorsFound = true;
            Factors.ct.carm++;
            if (verbose > 0)
                std::cout << "Carmichael factor found(" << ref << "), ctr = "
                << ctr << " i = " << i << " countdown = " << countdown << '\n';
        }
        else if (verbose > 1)
            std::cout << "Carmichael duplicate factor "
            << gcdVall << " found(" << ref  << "), i = " << i << " countdown = " << countdown
            << '\n';
    }
}

/* use psRand as a base. If p is a weak pseudoprime to that base, try to add 
   factors of p. p = a *2^ctr. Use: Xaux for square root of -1, mod p,
   Zaux for square root of 1 mod p.  */
static void FCarmichaelproc(const Znum& p, const Znum& a, const uint64_t PsRand, 
    Znum &Zaux, Znum &Xaux, fList &Factors, bool &sqrtOneFound, bool &sqrtMinusOneFound, 
    bool &factorsFound, const int ctr, const int countdown) {
    Znum Aux2, Aux3, Aux4, base;
    int i;

    base = PsRand;   /* change type to Znum */
    mpz_powm(ZT(Aux2), ZT(base), ZT(a), ZT(p));  /* Aux2 = base^a mod p */
    // If Aux2 = 1 or Aux2 = p-1, then try next pseudo-random number.
    if (Aux2 == 1 || Aux2 == p - 1) {
        return;    // This base cannot find a factor. Try another one.
    }

    for (i = 0; i < ctr; i++) {              // Loop that squares number.
        mpz_powm_ui(ZT(Aux3), ZT(Aux2), 2, ZT(p));   // Aux3 = Aux2^2 (mod p)
        if (Aux3 == 1)
        {            // Non-trivial square root of 1 found.
            if (!sqrtOneFound)
            {   // Save 1st non-trivial root to perform GCD later.
                Zaux = Aux2;
                sqrtOneFound = true;
            }
            else
            {          // Try to find non-trivial factor by doing GCD.
                Aux4 = Aux2 - Zaux;
                insertCarmichaelFactor(Aux4, p, Factors, factorsFound,
                    countdown, ctr, i, 1);
                Aux4 = Aux2 + Zaux;
                insertCarmichaelFactor(Aux4, p, Factors, factorsFound,
                    countdown, ctr, i, 2);
            }

            /* Aux2^2 = 1 (mod p) => (Aux2+1)*(Aux2-1) ≡ 0 (mod p)
               => Aux2+1 is a divisor of p */
            Aux4 = Aux2 + 1;
            insertCarmichaelFactor(Aux4, p, Factors, factorsFound,
                countdown, ctr, i, 3);
            return;  // Find more factors.
        }

        if (Aux3 == (p - 1))     // Aux3 = Aux2^2 (mod p)
        {   // Square root of -1 found.
            if (!sqrtMinusOneFound)
            {          // Save 1st non-trivial root to perform GCD later.
                Xaux = Aux2;
                sqrtMinusOneFound = true;
            }
            else
            {    // Try to find non-trivial factor by doing GCD.
                Aux4 = Aux2 - Xaux;   /* difference between 2 roots */
                insertCarmichaelFactor(Aux4, p, Factors, factorsFound,
                    countdown, ctr, i, -1);
            }
            return;  // Find more factors.
        }
        Aux2 = Aux3;    // Aux2 = Aux2^2 (mod p)
    }
}

/* called for pseudo-primes that have no small factors. 
Return: false = No factors found, true = factors found.
 Use: Xaux for square root of -1, mod p.
      Zaux for square root of 1 mod p. 
      We are in fact looking for bases for which p is a weak pseudo-prime. For
      Carmichael numbers, about 3/4 of numbers < p can be used. This method
      works whenever weak pseudo-prime bases can be found.

      We try  to identify a square root modulo n of −1, say R. Then, when 
      x^2 mod n ≡ n − 1 ≡ -1, we can compare the value of x against R: if x is 
      neither R nor n−R, then gcd(x − R, n) and gcd(x + R, n) are nontrivial 
      factors of n;
      R^2 ≡ -1, x^2 ≡ -1 (modulo n)
      x^2 - R^2 ≡ 0 (modulo n)
      (x-R)(x+R) ≡ 0 (modulo n)
      We also look for square roots of 1 (modulo n) in the same way.
*/
bool factorCarmichael(const Znum &p, fList &Factors, bool pseudoP)
{
    uint64_t PsRand = 0;  // pseudo-random number
    bool factorsFound = false;
    int countdown, ctr;
    bool sqrtOneFound = false;
    bool sqrtMinusOneFound = false;
    Znum a, Xaux, Zaux;
        
    a = p - 1;  // a = p - 1 (p is odd, so a is even)
    DivideZnumByMaxPowerOf2(ctr, a);  // a /= 2^ctr
    if (pseudoP) /* if p is a weak pseudoprime to base 2 */
        FCarmichaelproc(p, a, 2, Zaux, Xaux, Factors, sqrtOneFound,
            sqrtMinusOneFound, factorsFound, ctr, 21);

    /* we try up to 20 pseudo-random numbers. p is a pseudo-prime number.*/
    for (countdown = 20; countdown > 0; countdown--) {
        PsRand = ((uint64_t)PsRand * 89547121 + 1762281733) & 0x7fffffff;
        FCarmichaelproc(p, a, PsRand, Zaux, Xaux, Factors, sqrtOneFound,
            sqrtMinusOneFound, factorsFound, ctr, countdown);

         if (Factors.factComplete()) {
            if (verbose > 0)
                std::cout << "factorCarmichael has found all factors, in " << 21 - countdown
                << " cycles \n";
            break;   /* all factor are now prime, so stop. Without this test
             the outer loop would go round 20 times, probably finding the same
             factors over and over again. */
        }
    }
    if (!factorsFound && verbose > 0) {
        std::cout << "Carmichael: no factors found for " << p << '\n';
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
An alternative is to use the std::gcd function from the <numerics> library */
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
            factor = gcd( std::abs(x - x_fixed), num);
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
unsigned long long int PollardRho(unsigned long long int n, int depth)
{
    /* initialize random seed */
    std::random_device rd;   // non-deterministic generator
    std::mt19937_64 gen(rd());  // to seed mersenne twister.
    std::uniform_int_distribution<unsigned long long> dist(1, LLONG_MAX); // distribute results between 1 and MAX inclusive.
    long long ctr = 0;     /* loop counter*/
    /* no prime divisor for 1 */
    if (n == 1) return n;

    /* even number means one of the divisors is 2 */
    if (n % 2 == 0) return 2;

    /* we will pick from the range [2, N) */

    unsigned long long int x, y;
    x = dist(gen);  // x is in range 1 to max
    x = x % (n - 2) + 2;  // x is in range 2 to n-1
    y = x;


    /* the constant in f(x).
     * Algorithm can be re-run with a different c
     * if it throws failure for a composite. */
    long long int c = dist(gen);  // c is in range 1 to max
    c = c % (n - 1) + 1;  // c is in range 1 to n-1

/* Initialize candidate divisor (or result) */
    unsigned long long int d = 1;

    /* until the prime factor isn't obtained.
       If n is prime, return n */
    while (d == 1)
    {
        /* Tortoise Move: x(i+1) = f(x(i)) */
        x = (modPowerLL(x, 2, n) + c + n) % n;

        /* Hare Move: y(i+1) = f(f(y(i))) */
        y = (modPowerLL(y, 2, n) + c + n) % n;
        y = (modPowerLL(y, 2, n) + c + n) % n;

        /* check gcd of |x-y| and n */
        d = std::gcd(x>y?(x - y):(y-x), n);

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

/* trial division & Pollard-Rho. Uses first 33334 primes. Also weak
  pseudoprimes base 2 are factorised (maybe only partly) if < 190 digits. */
static void TrialDiv(fList &Factors, const unsigned long long PollardLimit) {
    bool restart = false;  // set true if trial division has to restart
    int upperBound;
    int numtype;
    unsigned long long testP;
    Znum temp;
    do {  /* use trial division */
        restart = false;    // change to true if any factor found
        for (int i = 0; i < Factors.f.size(); i++) {
            upperBound = Factors.f[i].upperBound;  // resume from where we left off
            if (upperBound == -1)
                continue;  // factor is prime
            temp = Factors.f[i].Factor;
            if (numLimbs(temp) <= 10) {  /* avoid this for numbers > 190 digits 
                  (seems faster on average not to do this for very large numbers) */
                numtype = mpz_bpsw_prp(ZT(temp));
                if (numtype == PRP_WPSP) {
                    /* number is a weak pseudoprime */
                    if (getfactors(temp, 2, Factors)) {
                        restart = true;
                        break;
                    }
                }
                else if (numtype == PRP_PRIME || numtype == PRP_PRP){
                    Factors.f[i].upperBound = -1; // show that residue is prime
                    continue;
                }
            }
            /* trial division. Uses first 33334 primes */
            //while (upperBound < std::min((int)prime_list_count,  155611)) {
            while (upperBound < std::min((int)prime_list_count, 33334)) {
                testP = primeList[upperBound];
                if (testP*testP > Factors.f[i].Factor) {
                    Factors.f[i].upperBound = -1; // show that residue is prime
                    break;
                }
                if (Factors.f[i].Factor%testP == 0) {
                    insertIntFactor(Factors, upperBound, 0, i);
                    restart = true;  // factor found; keep looking for more
                    Factors.ct.tdiv++;
                    break;
                }
                Factors.f[i].upperBound = upperBound;
                upperBound++;
            }

            if (!restart && (Factors.f[i].upperBound != -1)
                && Factors.f[i].Factor <= PollardLimit) {
                unsigned long long f;
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
                    f = PollardRho(ZnumToLong(Factors.f[i].Factor));
                    //f = SQUFOF(MulPrToLong(Factors.f[i].Factor));
                    /* there is a small possibility that PollardFactor won't work,
                    even when factor is not prime*/
                    if (f > 1) {
                        insertIntFactor(Factors, -1, f, i);
                        Factors.ct.prho++;
                    }
                }
            }
        }
    } while (restart);  // keep looping until no more factors found.

    if (verbose >= 1) {
        if (lang)
            std::cout << "fin de la división de prueba. " << Factors.f.size() - 1
            << " factores encontrados hasta ahora \n";
        else
            std::cout << "End Trial division. " << Factors.f.size() - 1 << " factors found so far \n";
        if (Factors.f.size() > 1) {
            if (lang)
                std::cout << "resultado despues de la division de prueba ";
            else
                std::cout << "result after trial division ";
            Factors.Xprint();
        }
    }
}

/* factorise toFactor; factor list returned in Factors. 
returns false only if ecm returns an error. */
static bool factor(fList &Factors) {
    Znum toFactor = Factors.n;
    unsigned long long testP;
    const unsigned long long MaxP = 393203;       // use 1st 33334 primes 
    //const unsigned long long MaxP = 2'097'143;  // use 1st 155611 primes
    /* larger value for Maxp seems to slow down factorisation overall. */
    // MaxP must never exceed 2'097'152 to avoid overflow.
    const unsigned long long PollardLimit = MaxP*MaxP*MaxP;

    if (toFactor <= 3)
        return true;   /* toFactor is 1 or prime so we can stop now*/

    if ((long long)primeListMax <MaxP) {  // get primes
        generatePrimes(MaxP);  // takes a while, but only needed on 1st call
    }
    
    oldTimeElapsed = 0;
    originalTenthSecond = tenths();    // record start time

    if (toFactor >= MaxP* MaxP) {
        /* may not be able to factorise entirely by trial division, so try this first */
        PowerPM1Check(Factors, toFactor, 2);  // check if toFactor is a perfect power +/- 1
        Factors.ct.pm1 = (int)Factors.f.size() - 1;  // number of factors just found, if any
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
            Factors.ct.powerCnt++;
        }
        
        int result = PrimalityTest(Zpower, testP- 1);
        if (result == 0) {   // Number is prime.
            Factors.f[i].upperBound = -1; // Indicate that number is prime.
            continue;
        }
        if (result > 1) {  /* number is a pseudo-prime, but is NOT prime */
            if (factorCarmichael(Zpower, Factors, result ==2)) {
                i = -1;			// restart loop at beginning!!
                continue;
            }
            //else if(result == 2)
            //    if (getfactors(Zpower, 2, Factors)) {
            //        i = -1;			// restart loop at beginning!!
            //        continue;
            //    }
        }
        if (Zpower <= PollardLimit) {
            unsigned long long f;
            f = PollardRho(ZnumToLong(Zpower));
            if (f != 1) {
                insertIntFactor(Factors, -1, f, i);
                Factors.ct.prho++;
                /* there is a small possibility that PollardFactor won't work,
                even when factor is not prime*/
                continue;
            }
        }

        int nooflimbs = numLimbs(Zpower);
        // get approximate size (1 limb = 64 bits)
        if (nooflimbs <=3 || (!msieve && !yafu && !Pari)) {
            /* use built-in ECM & SIQS if number to factor <= 192 bits (58 digits)
               because this is fastest for smaller numbers,
               or if both YAFU and Msieve are turned off */
            //ElipCurvNo = 1;  // start with 1st curve
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
                    Factors.ct.leh++;     // Factor found.
                    insertBigFactor(Factors, Zfactor);
                    i = -1;   // success; restart loop at beginning to tidy up!	
                    break;
                }
            }
            if (i == -1)
                continue;   // restart loop from beginning

            size_t fsave = Factors.f.size();
            bool rv; 
            /* one and only one of msieve, yafu and Pari should be set.*/
            if (msieve)  
                rv = callMsieve(Zpower, Factors);
            else if (yafu)
                rv = callYafu(Zpower, Factors);
            else {
                parifactor(Zpower, Factors);
                rv = true;
            }
            if (rv) {
                i = -1;   // success; restart loop at beginning to tidy up!
                // record any increase in number of factors
                if (msieve)
                    Factors.ct.msieve += (int)(Factors.f.size() - fsave); 
                else if (yafu)
                    Factors.ct.yafu += (int)(Factors.f.size() - fsave); 
                else 
                    Factors.ct.paric += (int)(Factors.f.size() - fsave);
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

/* compute 3 values the squares of which add up to prime p 
This function is intended for 8k+1 primes if we want p = x^2 + 2*y^2 
rather than p = x^2 + y^2. */
static void compute3squares(const Znum& p, Znum Mult[4]) {
    Znum q, Tmp, Tmp1, Tmp2, Tmp3, K, M1, M2, M3, M4;
    std::vector <Znum> roots;

    q = (p - 1) / 2; //  q = (prime-1)/2
    Mult[0] = 0;
    do {
        Mult[0]++;
        Tmp = Mult[0] * Mult[0] + 1;
        Tmp = -Tmp;
        while (Tmp < 0) Tmp += p;   /* Tmp = -1 - Mult[1]^2 (mod p) */
        mpz_powm(ZT(Tmp1), ZT(Tmp), ZT(q), ZT(p));
        // At this moment Tmp1 = (-1 - Mult[0]^2)^((p-1)/2)(Mod p)
    } while (Tmp1 != 1);  // Continue loop if it is not 1.

    // After the loop finishes, Tmp = (-1 - Mult[0]^2) is a quadratic residue mod p.
    Tmp1 = -1 - Mult[0] * Mult[0];
    roots = primeModSqrt(Tmp1, p);  /* use Tonelli-Shanks to get Mod sqrt */
    assert(!roots.empty());
    Mult[1] = roots[0];
    /* at this point Mult[0]^2 + Mult[1]^2 + 1 ≡ 0 (mod p)*/

    if (verbose > 1) {
        std::cout << "3 squares: prime = " << p << " initial Mult[0] = " << Mult[0]
            << " Mult[1] = " << Mult[1] << '\n';
    }

    Mult[2] = 1;
    Mult[3] = 0;
    for (;;) {
        // Compute K <- (Mult[0]^2 + Mult[1]^2 + Mult[2]^2 + Mult[3]^2) / p
        Tmp = Mult[0] * Mult[0] + Mult[1] * Mult[1];
        Tmp = Tmp + Mult[2] * Mult[2] + Mult[3] * Mult[3];
        assert(Tmp % p == 0);

        K = Tmp / p;    // in this loop K is smaller each time round
        assert(K > 0);
        if (K == 1) {
            break;  // we are done when K equals 1
        }
        if (isEven(K)) { // If K is even ...
            if (isEven(Mult[0]) != isEven(Mult[1]))
            {  // If Mult[0] + Mult[1] is odd...
                if (isEven(Mult[0]) == isEven(Mult[2])) {
                    // If Mult[0] + Mult[2] is even...
                    Mult[1].swap(Mult[2]);  // swap Mult[1] and Mult[2]
                }
                else {
                    Mult[1].swap(Mult[3]); // swap Mult[1] and Mult[3]
                }
            } // At this moment Mult[0]+Mult[1] = even, Mult[2]+Mult[3] = even
            Tmp1 = (Mult[0] + Mult[1]) / 2;
            Tmp2 = (Mult[0] - Mult[1]) / 2;
            Tmp3 = (Mult[2] + Mult[3]) / 2;
            Mult[3] = (Mult[2] - Mult[3]) / 2;
            Mult[2] = Tmp3;
            Mult[1] = Tmp2;
            Mult[0] = Tmp1;
            if (verbose > 1) {
                std::cout << "K = " << K;
                std::cout << " Mult[0]-[3] = " << Mult[0] << ", " << Mult[1] << ", " << Mult[2]
                    << ", " << Mult[3] << '\n';
            }
            continue;
        } /* end if K is even */

        M1 = Mult[0] % K;
        if (M1 < 0) {
            M1 += K;
        }
        M2 = Mult[1] % K;
        if (M2 < 0) {
            M2 = M2 + K;
        }
        M3 = Mult[2] % K;
        if (M3 < 0) {
            M3 += K;
        }
        M4 = Mult[3] % K;
        if (M4 < 0) {
            M4 += K;
        }
        Tmp = (K + 1) / 2;  // Tmp <- (K+1) / 2
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
        // Compute Tmp1 <- (Mult[0]*M1 + Mult[1]*M2 + Mult[2]*M3 + Mult[3]*M4) / K
        Tmp = Mult[0] * M1 + Mult[1] * M2 + Mult[2] * M3 + Mult[3] * M4;
        Tmp1 = Tmp / K; // BigIntDivide(Tmp, K, Tmp1);

        // Compute Tmp2 <- (Mult[0]*M2 - Mult[1]*M1 + Mult[2]*M4 - Mult[3]*M3) / K
        Tmp = Mult[0] * M2 - Mult[1] * M1 + Mult[2] * M4 - Mult[3] * M3;
        Tmp2 = Tmp / K;

        // Compute Tmp3 <- (Mult[0]*M3 - Mult[2]*M1 - Mult[1]*M4 + Mult[3]*M2) / K
        Tmp = Mult[0] * M3 - Mult[2] * M1 - Mult[1] * M4 + Mult[3] * M2;
        Tmp3 = Tmp / K;

        // Compute Mult[3] <- (Mult[0]*M4 - Mult[3]*M1 + Mult[1]*M3 - Mult[2]*M2) / K
        Tmp = Mult[0] * M4 - Mult[3] * M1 + Mult[1] * M3 - Mult[2] * M2; // BigIntSubt(Tmp, Tmp4, Tmp);
        Mult[3] = Tmp / K;

        Mult[2] = Tmp3;
        Mult[1] = Tmp2;
        Mult[0] = Tmp1;
        if (verbose > 1) {
            std::cout << "K = " << K;
            std::cout << " Mult[0]-[3] = " << Mult[0] << ", " << Mult[1] << ", " << Mult[2]
                << ", " << Mult[3] << '\n';
        }
    }
    Mult[0] = abs(Mult[0]);   // ensure results are +ve
    Mult[1] = abs(Mult[1]);
    Mult[2] = abs(Mult[2]);
    Mult[3] = abs(Mult[3]);

    /* ensure that Mult[3] is smallest */
    if (Mult[3] > Mult[2])
        Mult[3].swap(Mult[2]);
    if (Mult[3] > Mult[1])
        Mult[3].swap(Mult[1]);
    if (Mult[3] > Mult[0])
        Mult[3].swap(Mult[0]);

    /* if 2 values are equal, ensure that they are 2nd and 3rd */
    if (Mult[0] == Mult[1])
        Mult[0].swap(Mult[2]);

    assert(Mult[3] == 0);
}

/* compute 4 or less values the squares of which add up to prime p.   
return values in Mult[0] to Mult[3] 
For odd p:
if p = 1 (mod 4) p can be expressed as the sum of 2 squares. 
If p = 1 or 3 (mod 8) it can be expressed as the sum of 3 squares. Two of the 3
squares will be equal. 
If p = 7 (mod 8) it can only be expressed as the sum of 4 squares. */

static void ComputeFourSquares(const Znum &p, Znum Mult[4], const bool sqplustwosq) {
    Znum a, q, K, Tmp, Tmp1, Tmp2, Tmp3, Tmp4, M1, M2, M3, M4; 
    Znum TestNbr;
    std::vector<Znum>roots;

    if (p == 2) {   /* Prime factor is 2 */
        Mult[0] = 1;  // 2 = 1^2 + 1^2 + 0^2 + 0^2
        Mult[1] = 1;
        Mult[2] = 0;
        Mult[3] = 0;
        return;
    }
    else  /* Prime factor p is not 2 */ {
        if ((p & 3) == 1)  /* if p = 1 (mod 4) */ {
            /* if p = 1 (mod 4) p can be expressed as the sum of 2 squares. */
            if (sqplustwosq && ((p & 7) == 1)) {
                /* If p = 1 (mod 8) it can also be expressed as the sum of 3 
                squares. Two of the 3 squares will be equal. Get p as the sum 
                of 3 squares rather than 2 */
                compute3squares(p, Mult);
                Mult[3] = 0;
                return;
            }
            /* in this case p can be expressed as the sum of 2 squares */
            roots = primeModSqrt(-1, p);  /* use Legendre formula or Tonelli-Shanks to get mod sqrt */
            assert(roots.size() >= 2);
            Mult[0] = std::min(roots[0], roots[1]);   /* Mult[0]² = p-1(mod p)*/

            if (verbose > 1) {
                std::cout << "4 squares: prime = " << p << " initial Mult[0] = " << Mult[0] << '\n';
            }
            Mult[1] = 1;

            for (;;) {  
                Tmp = Mult[0] * Mult[0] + Mult[1] * Mult[1]; // K <- (Mult[0]^2 + Mult[1]^2) / p
                K = Tmp / p;    // in this loop K is smaller each time round   
                if (verbose > 1) {
                    std::cout << "K = " << K << '\n';
                }
                if (K == 1) {  // are we there yet?
                    Mult[2] = 0;
                    Mult[3] = 0;
                    break;     // we are finished
                }
                M1 = Mult[0] % K;
                if (M1 < 0) {
                    M1 += K;
                }
                M2 = Mult[1] % K;
                if (M2 < 0) {
                    M2 += K;
                }
                Tmp = (K+1)/2;   //subtractdivide(K, -1, 2);      
                if (M1 >= Tmp) // If M1 >= K / 2 ... 
                {
                    M1 -= K;
                }
                if (M2 >= Tmp) {     // If M2 >= K / 2 ...     
                    M2 -= K;
                }
                Tmp = Mult[0]*M1 + Mult[1]*M2;
                Tmp2 = Tmp / K;  // Tmp2 <- (Mult[0]*m1 + Mult[1]*m2) / K
                Tmp = Mult[0]*M2 - Mult[1]*M1;
                Mult[1] = Tmp / K;    // Mult[1] <- (Mult[0]*m2 - Mult[1]*m1) /K
                Mult[0] = Tmp2;       // Mult[0] <- (Mult[0]*m1 + Mult[1]*m2) / K
                if (verbose > 1) {
                    std::cout << "Mult[0] = " << Mult[0] << " Mult[1] = " << Mult[1] << '\n';
                }
            } /* end for */
        } /* end p = 1 (mod 4) */

        else  /* if p = 3 (mod 4) */ 	{
            q = (p-1)/2; //  q = (prime-1)/2
            Mult[0] = 0;
            do {
                Mult[0]++;
                Tmp = Mult[0]*Mult[0] + 1;
                Tmp = -Tmp;
                while (Tmp < 0) Tmp += p;   /* Tmp = -1 - Mult[0]^2 (mod p) */
                mpz_powm(ZT(Tmp1), ZT(Tmp), ZT(q), ZT(p));
                       // At this moment Tmp1 = (-1 - Mult[0]^2)^((p-1)/2)(Mod p)
            } while (Tmp1 != 1);  // Continue loop if it is not 1.

            // After the loop finishes, Tmp = (-1 - Mult[0]^2) is a quadratic residue mod p.

            q = (p+1)/4;

            // Find Mult[1] <- square root of Tmp1 = Tmp^q (mod p) (Lagrange solution)
            mpz_powm(ZT(Mult[1]), ZT(Tmp), ZT(q), ZT(p));

            /* at this point Mult[0]^2 + Mult[1]^2 + 1 ≡ 0 (mod p)*/
            if (verbose > 1) {
                std::cout << "4 squares: prime = " << p << " initial Mult[0] = " << Mult[0] 
                    << " Mult[1] = " << Mult[1] << '\n';
            }

            Mult[2] = 1;
            Mult[3] = 0;

            for (;;) {
                // Compute K <- (Mult[0]^2 + Mult[1]^2 + Mult[2]^2 + Mult[3]^2) / p
                Tmp = Mult[0]*Mult[0] + Mult[1]*Mult[1];
                Tmp = Tmp + Mult[2]*Mult[2] + Mult[3]*Mult[3];
                assert(Tmp%p == 0);
                K = Tmp / p;    // in this loop K is smaller each time round
                assert(K > 0);
                if (K == 1) {
                    break;  // we are done when K equals 1
                }


                if (isEven(K)) { // If K is even ...
                    if (isEven(Mult[0]) != isEven(Mult[1]))
                    {  // If Mult[0] + Mult[1] is odd...
                        if (isEven(Mult[0]) == isEven(Mult[2])) {  
                            // If Mult[0] + Mult[2] is even...
                            Mult[1].swap(Mult[2]);  // swap Mult[1] and Mult[2]
                        }
                        else {
                            Mult[1].swap(Mult[3]); // swap Mult[1] and Mult[3]
                        }
                    } // At this moment Mult[0]+Mult[1] = even, Mult[2]+Mult[3] = even
                    Tmp1 = (Mult[0] + Mult[1])/2;   
                    Tmp2 = (Mult[0] - Mult[1])/2;   
                    Tmp3 = (Mult[2] + Mult[3])/2;   
                    Mult[3] = (Mult[2] - Mult[3])/2 ;
                    Mult[2] = Tmp3;
                    Mult[1] = Tmp2;
                    Mult[0] = Tmp1;
                    if (verbose > 1) {
                        std::cout << "K = " << K;
                        std::cout << " Mult[0]-[3] = " << Mult[0] << ", " << Mult[1] << ", " << Mult[2]
                            << ", " << Mult[3] << '\n';
                    }
                    continue;
                } /* end if K is even */

                M1 = Mult[0] % K;
                if (M1 < 0) {
                    M1 += K;
                }
                M2 = Mult[1] % K;
                if (M2 < 0) {
                    M2 = M2 + K;
                }
                M3 = Mult[2] % K;
                if (M3 < 0) {
                    M3 += K;
                }
                M4 = Mult[3] % K;
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
                // Compute Tmp1 <- (Mult[0]*M1 + Mult[1]*M2 + Mult[2]*M3 + Mult[3]*M4) / K
                Tmp = Mult[0]*M1 + Mult[1]*M2 + Mult[2]*M3 + Mult[3]*M4;
                Tmp1 = Tmp / K; // BigIntDivide(Tmp, K, Tmp1);

                // Compute Tmp2 <- (Mult[0]*M2 - Mult[1]*M1 + Mult[2]*M4 - Mult[3]*M3) / K
                Tmp = Mult[0]*M2 - Mult[1]*M1 + Mult[2]*M4 - Mult[3]*M3;
                Tmp2 = Tmp / K;

                // Compute Tmp3 <- (Mult[0]*M3 - Mult[2]*M1 - Mult[1]*M4 + Mult[3]*M2) / K
                Tmp = Mult[0]*M3 - Mult[2]*M1 - Mult[1]*M4 + Mult[3]*M2;
                Tmp3 = Tmp / K;

                // Compute Mult[3] <- (Mult[0]*M4 - Mult[3]*M1 + Mult[1]*M3 - Mult[2]*M2) / K
                Tmp = Mult[0]*M4 - Mult[3]*M1 + Mult[1]*M3 - Mult[2]*M2; // BigIntSubt(Tmp, Tmp4, Tmp);
                Mult[3] = Tmp / K;

                Mult[2] = Tmp3;
                Mult[1] = Tmp2;
                Mult[0] = Tmp1;
                if (verbose > 1) {
                    std::cout << "K = " << K;
                    std::cout << " Mult[0]-[3] = " << Mult[0] << ", " << Mult[1] << ", " << Mult[2] 
                        << ", " << Mult[3] << '\n';
                }
            } /* end for */
        } /* end if p = 3 (mod 4) */
    } /* end prime not 2 */

    Mult[0] = abs(Mult[0]);   // ensure results are +ve
    Mult[1] = abs(Mult[1]);
    Mult[2] = abs(Mult[2]);
    Mult[3] = abs(Mult[3]);

    /* ensure that Mult[3] is smallest */
    if (Mult[3] > Mult[2])
        Mult[3].swap(Mult[2]);
    if (Mult[3] > Mult[1])
        Mult[3].swap(Mult[1]);
    if (Mult[3] > Mult[0])
        Mult[3].swap(Mult[0]);

    /* if 2 values are equal, ensure that they are 2nd and 3rd */
    if (Mult[0] == Mult[1])
        Mult[0].swap(Mult[2]);
}


/* compute 3 values the squares of which add up to s * 2^r, return values in quads */
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
    
        for (s3 = s2, m = 0; isEven(s3); m++) {
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
            ComputeFourSquares(s3, quads, false);
            if (quads[0] < quads[1]) {		
                quads[0].swap(quads[1]);  // quads[0] < quads[1], so exchange them.
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

            /* sort into ascending order */
            if (quads[0] < quads[1]) {
                quads[0].swap(quads[1]);  // quads[0] < quads[1], so exchange them.
            }

            if (quads[0] < quads[2]) {
                quads[0].swap(quads[2]);  // quads[0] < quads[2], so exchange them.
            }

            if (quads[1] < quads[2]) {
                quads[1].swap(quads[2]);  // quads[1] < quads[2], so exchange them.
            }
            if (verbose > 1) {
                std::cout << "compute3squares(" << r << ", " << s << ")\n"
                    << "quads = " << quads[0] << ", " << quads[1] 
                    << ", " << quads[2] << '\n';
            }
            return;
        }
    }
}

/* show that the number is the sum of 4 or fewer squares. See
https://www.alpertron.com.ar/4SQUARES.HTM */
/* uses the identity:
(a²+ b²+ c²+ d²)*(A²+ B²+ C²+ D²) = (aA+bB+cC+dD)² + (aB-bA+cD-dC)²
                                  + (aC-bD-cA+dB)² + (aD-dA+bC-cB)²
This allows us to find the sum of squares for each factor separately then combine them 
If the number is a perfect square just return 1 value.
If the number can be expressed as the sum of 2 squares return 2 values.
If the number can be expressed as a square + twice a square, try to find 
appropriate values.  
If the number can be expressed as the sum of 3 squares try to find appropriate values.
For large numbers this may not be practical. 
The fallback is to find 4 squares, if none of the above are feasible- */
static void ComputeFourSquares(const fList &factorlist, Znum quads[4], Znum num) {
    Znum Mult[4], Tmp1, Tmp2, Tmp3;
    Znum pr;
    bool twoSq;       /* value changed to true if num can be expressed as the sum 
                        of 2 squares */
    bool sqplustwosq; /* value changed to true if num can be expressed as the sum 
                        of a square + twice a square */

    quads[0] = 1;      /* initialise quads N.B. 1 = 1^2 + 0^2 + 0^2 + 0^2 */
    quads[1] = 0;
    quads[2] = 0;
    quads[3] = 0;

    if (factorlist.f.size() == 1) { /* only 1 factor? */
        if (factorlist.f[0].Factor == 1) {   // Number to factor is 1.
            return;   // 1 = 1^2 + 0^2 + 0^2 + 0^2
        }
        if (factorlist.f[0].Factor == 0) {    // Number to factor is 0.
            quads[0] = 0;      // 0 = 0^2 + 0^2 + 0^2 + 0^2
            return;
        }
    }

    /* check whether number can be formed as sum of 1 or 2 squares */
    twoSq = factorlist.twosq();
    if (!twoSq)
        sqplustwosq = factorlist.sqplustwosq();  /* true iff num = x^2 + 2*y^2 */
    else
        sqplustwosq = false;

    if (!twoSq) {  /* check whether number can be expressed as sum of 3 squares */
        int r = 0;
        while ((num & 3) == 0) {
            num /= 4;
            r++;
        }
        /* any number which is not of the form 4^r * (8k+7) can be formed as the sum of 3 squares 
        see https://en.wikipedia.org/wiki/Legendre%27s_three-square_theorem */
        if  (factorlist.f.size() > 1) 
            if (numLimbs(num) < 4) 
                if ((num & 7) < 7 ) {
            /* use compute3squares if number has more than 1 unique prime factor, 
            is small (<= 57 digits), and can be formed from 3 squares. 
            (Each limb is up to 64 bits). Large numbers would take too long.  */
            if (!sqplustwosq) {
                compute3squares(r, num, quads);
                return;
            }
        }
    }

    /* the method below will find 4 squares the sum of which is the required number.
    In many cases it would be possible to use just 3 squares, but the method for
    that would be too slow. */
    for (auto &Factorx : factorlist.f) {
        if (Factorx.exponent % 2 == 0) {
            continue; /* if Prime factor appears an even number of times, no need to
                      process it in this for loop */
        }
        
        pr = Factorx.Factor;

        /* compute 4 or less values the squares of which add up to prime pr,
        return values in Mult[0], Mult[1], Mult[2] and Mult[3] */
        ComputeFourSquares(pr, Mult, sqplustwosq);
        if (sqplustwosq) {
            assert(Mult[3] == 0);
            /* swap so that Mult[1] = Mult[2]. This ensures that when Mult is combined
             into quads we preserve the feature that quads[1] = quads[2] and quads[3] is zero */
            if (Mult[0] == Mult[1])
                Mult[0].swap(Mult[2]);
            else if (Mult[0] == Mult[2])
                Mult[0].swap(Mult[1]);
            else if (Mult[1] != Mult[2]) {
                std::cout << "** expected 2 of 3 values to be equal \n";
                if (verbose <= 1)
                    std::cout << "pr    = " << pr
                        << "\nMult[0] = " << Mult[0]
                        << "\nMult[1] = " << Mult[1]
                        << "\nMult[2] = " << Mult[2]
                        << "\nMult[3] = " << Mult[3] << '\n';
            }
        }
        assert(pr == Mult[0] * Mult[0] + Mult[1] * Mult[1] + Mult[2] * Mult[2] + Mult[3] * Mult[3]);
        if (verbose > 1) {
            std::cout <<   "pr    = " << pr 
                      << "\nMult[0] = " << Mult[0]
                      << "\nMult[1] = " << Mult[1]
                      << "\nMult[2] = " << Mult[2]
                      << "\nMult[3] = " << Mult[3] << '\n' ;
        }

        /* use the identity:
        (a²+ b²+ c²+ d²)*(A²+ B²+ C²+ D²) = (aA+bB+cC+dD)² + (aB-bA+cD-dC)²
                                            + (aC-bD-cA+dB)² + (aD-dA+bC-cB)² 
        note: if twoSq is true, c, d, C, & D are zero. The expression then 
        simplifies 	automatically to:  	
               (a²+b²)*(A²+B²) = (aA+bB)² + (aB-bA)² 
               
        Also, if d =0, D = 0, b = c, and B = C we effectively get
        (a²+ 2.b²)*(A²+ 2B²) = (aA+2bB)² + (aB-bA)²  + (aB-bA)² + (bB-bB)²
                                  = (aA+2bB)² + 2(aB-bA)²  
        This occurs when all odd prime factors of the form 8i+5 or 8i+7 have even exponents.
*/

        Tmp1 = Mult[0]*quads[0] + Mult[1]*quads[1] + Mult[2]*quads[2] + Mult[3]*quads[3];
        Tmp2 = Mult[0]*quads[1] - Mult[1]*quads[0] + Mult[2]*quads[3] - Mult[3]*quads[2];
        Tmp3 = Mult[0]*quads[2] - Mult[2]*quads[0] - Mult[1]*quads[3] + Mult[3]*quads[1];
        quads[3] = Mult[0]*quads[3] - Mult[3]*quads[0] + Mult[1]*quads[2] - Mult[2]*quads[1];

        quads[2] = Tmp3;
        quads[1] = Tmp2;
        quads[0] = Tmp1;
    } 

     /* for factors that are perfect squares, multiply quads[0]-[3] by sqrt(factor) */
    for (auto Factorx : factorlist.f) {
        if (Factorx.exponent >= 2) {
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

    /* firstly, put largest value in quads[0] */
    if (quads[0] < quads[1]) {		// quads[0] < quads[1], so exchange them.
        quads[0].swap(quads[1]);
    }

    if (quads[0] < quads[2]) {	// quads[0] < quads[2], so exchange them.
        quads[0].swap(quads[2]);
    }

    if (quads[0] < quads[3]) {	// quads[0] < quads[3], so exchange them.
        quads[0].swap(quads[3]);
    }

    /* put 2nd largest value in quads[1] */
    if (quads[1] < quads[2]) {	// quads[1] < quads[2], so exchange them.
        quads[1].swap(quads[2]);
    }

    if (quads[1] < quads[3]) {	// quads[1] < quads[3], so exchange them.
        quads[1].swap(quads[3]);
    }

    /* put 3nd largest value in quads[2] and smallest in quads[4] */
    if (quads[2] < quads[3]) {	// quads[2] < quads[3], so exchange them.
        quads[2].swap(quads[3]);
    }
    return;
}

/* store factors for larger numbers */
std::set<fList> savedFactors;
long long SFhitcount = 0;     /* number of 'hits' searching savedFactors */
long long SFmisscount = 0;   /* number of 'misses' searching savedFactors */

/********************************************************************************
code below is to interface between DA's code
and the new code that uses Znums, which are really mpz_t integers from MPIR or GMP
multiprecision library, with a C++ class wrapped around them that allows them
to be used pretty much like normal integers. 
**********************************************************************************/

/* factorise number. Returns false if unable to factorise it. If the number is 
more than 1 limb (64 bits) the factor list is memoised, so if the same number is
factorised again the saved value is used instead of repeating the factorisation. */
bool factorise(Znum numberZ, fList &vfactors, Znum quads[]) {

    try {
        bool pos = true;
        bool hit = false;  /* set true if numberZ factors found in cache */

        if (numberZ == 0)
            return false;  // function factor can't factorize zero
        if (numberZ < 0) {
            pos = false;
            numberZ = -numberZ;
        }
        vfactors.set(numberZ);
        if (numLimbs(numberZ) > 1) {
        /* see whether the number has already been factorised. Note that
           the find() method only looks for a match in the number to be factored */
            auto ref = savedFactors.find(vfactors);
            if (ref != savedFactors.end()) {
                vfactors = *ref;  /* use factor list obtained earlier */
                SFhitcount++;
                hit = true;
                if (verbose > 1) {
                    std::cout << "factors of " << numberZ << " found from cache \n";
                }
            }
        }
        if (!hit) {
            auto rv = factor(vfactors);
            if (!rv)
                return false;  // failed to factorise number
            if (numLimbs(numberZ) > 1) {
                savedFactors.insert(vfactors);  /* save factors just found */
                SFmisscount++;
            }
        }

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

#define nelems(x) (sizeof(x) / sizeof((x)[0]))

/* return a factor of N, using Shanks's square forms factorization method. Based on 
Wikipedia, as modified in 
https://stackoverflow.com/questions/52746812/shankss-square-form-factorization-implementation 
assume that N is not prime. 
An alternative to Pollard-Rho that seems to take about the same amount of time */
uint64_t SQUFOF(const uint64_t N) {
    uint64_t D, Po, P, Pprev, Q, Qprev, q, b, r;
    uint32_t L, B, i;
    /* all combinations of primes 3, 5, 7, 11*/
    static const int multiplier[] = { 1, 3, 5, 7, 11, 3*5, 3*7, 3*11, 5*7, 
        5*11, 7*11, 3*5*7, 3*5*11, 3*7*11, 5*7*11, 3*5*7*11 };
    /* smallest factor of each multiplier */
    static const uint64_t results[] = {
        1, 3, 5, 7, 11, 3, 3, 3, 5, 
        5, 7, 3, 3, 3, 5, 3	};
    const uint64_t s = (uint64_t)(sqrtl((double)N) + 0.5);

    if (s * s == N) 
        return s;  /* N is a perfect square */

    /* note that there is a check to prevent overflow; exit the loop if overflow 
    would occur */
    for (int k = 0; k < nelems(multiplier) && N <= UINT64_MAX / multiplier[k]; k++) {
        if (multiplier[k] == N)
            return results[k];
        D = multiplier[k] * N;  /* overflow may occur if N > 15971206990224720  */

        Po = Pprev = P = (uint64_t)sqrtl((double)D);
        Qprev = 1;
        Q = D - Po * Po;
        if (Q == 0)
            /* only happens if D is a perfect square, which implies that N is a
             multiple of multiplier[k] i.e. N has at least one small factor. */
            return results[k];

        L = 2 * (uint32_t)sqrtl(2.0 * s);
        B = 3 * L;
        for (i = 2; i < B; i++) {
            b = (uint64_t)((Po + P) / Q);
            P = b * Q - P;
            q = Q;
            Q = Qprev + b * (Pprev - P);
            r = (uint64_t)(sqrtl((double)Q) + 0.5);
            if (!(i & 1) && r * r == Q) 
                break;
            Qprev = q;
            Pprev = P;
        };
        if (i >= B) 
            continue;
        b = (uint64_t)((Po - P) / r);
        Pprev = P = b * r + P;
        Qprev = r;
        Q = (D - Pprev * Pprev) / Qprev;
        i = 0;
        do {
            b = (uint64_t)((Po + P) / Q);
            Pprev = P;
            P = b * Q - P;
            q = Q;
            Q = Qprev + b * (Pprev - P);
            Qprev = q;
            i++;
        } while (P != Pprev);
        r = gcd(N, Qprev);
        if (r != 1 && r != N) {
            if (verbose > 0)
                std::cout << "SQUFOF n = " << N << " factor = " << r << '\n';
            return r;
        }
    }

    return 1;  /* failed to find factor*/
}
