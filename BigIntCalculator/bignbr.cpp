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

/* nbr = - nbr */
void ChSignBigNbr(int nbr[], int length)
{
    int carry = 0;
    int ctr;
    for (ctr = 0; ctr < length; ctr++)
    {
        carry -= nbr[ctr];
        nbr[ctr] = carry & MAX_INT_NBR;
        carry >>= BITS_PER_INT_GROUP;
    }
}

/* nbr = - nbr */
void ChSignBigNbrB(int nbr[], int length)
{
    int carry = 0;
    int ctr;
    for (ctr = 0; ctr < length - 1; ctr++)
    {
        carry -= nbr[ctr];
        nbr[ctr] = carry & MAX_INT_NBR;
        carry >>= BITS_PER_INT_GROUP;
    }
    nbr[ctr] = carry - nbr[ctr];  /* last word does not have most significant bit masked off */
}

/* get number of limbs in Nbr */
int BigNbrLen(const long long Nbr[], int nbrLen) {
    int ix;
    for (ix = nbrLen; ix > 0; ix--) {
        if (Nbr[ix - 1] != 0)
            break;
    }
    return ix;
}

/* Prod = Factor1 *Factor2. Factor1 & Factor2 have same length, length of product
 is twice length of Factor1, Factor2*/
void MultBigNbr(const limb pFactor1[], const limb pFactor2[], limb pProd[], int nbrLen)
{
    limb* ptrProd = pProd;
    double dRangeLimb = (double)(1U << BITS_PER_GROUP);
    double dInvRangeLimb = 1.0 / dRangeLimb;
    int low = 0;
    int factor1;
    int factor2;
    double dAccumulator = 0.0;
    for (int i = 0; i < nbrLen; i++)
    {
    	for (int j = 0; j <= i; j++)
    	{
    		factor1 = *(pFactor1 + j);
    		factor2 = *(pFactor2 + i - j);
    		low += factor1 * factor2;
    		dAccumulator += (double)factor1 * (double)factor2;
    	}
    	low &= MAX_INT_NBR;    // Trim extra bits.
    	*ptrProd = low;
    	ptrProd++;
     //Subtract or add 0x20000000 so the multiplication by dVal is not nearly an integer.
     //In that case, there would be an error of +/- 1.
    	if (low < HALF_INT_RANGE)
    	{
    		dAccumulator = floor((dAccumulator + (double)FOURTH_INT_RANGE) * dInvRangeLimb);
    	}
    	else
    	{
    		dAccumulator = floor((dAccumulator - (double)FOURTH_INT_RANGE) * dInvRangeLimb);
    	}
    	low = (int)(dAccumulator - floor(dAccumulator * dInvRangeLimb) * dRangeLimb);
    }
    *ptrProd = low;
    *(ptrProd + 1) = (int)floor(dAccumulator / dRangeLimb);
}

/* Quotient = Dividend/divisor. The quotient is rounded towards zero 
(i.e. truncation division). If divisor is 0 a divide-by-zero error will occur.
the remainder is returned in rem */
void DivBigNbrByInt(const limb Dividend[], int divisor, limb Quotient[], int nbrLen, int *rem)
{
    int ctr;
    int remainder = 0;
    double dDivisor = (double)abs(divisor);  // assume divisor != 0
    double dLimb = INT_RANGE_U;
    for (ctr = nbrLen - 1; ctr >= 0; ctr--)
    {
        double dDividend, dQuotient;
        int quotient, dividend;
        dividend = (remainder << BITS_PER_INT_GROUP) + Dividend[ctr];
        dDividend = (double)remainder * dLimb + Dividend[ctr];
        dQuotient = std::floor(dDividend / dDivisor + 0.5);
        quotient = (unsigned int)dQuotient;   // quotient has correct value or 1 more.
        remainder = dividend - quotient * divisor;
        if ((unsigned int)remainder >= (unsigned int)divisor)
        {     // remainder not in range 0 <= remainder < divisor. Adjust.
            quotient--;
            remainder += divisor;
        }
        Quotient[ctr] = quotient;
    }
    *rem = remainder;
}


/* calculate NbrMod^Expon%currentPrime.
calculation uses doubles to avoid risk of overflow. */
int modPower (int NbrMod, int Expon, int currentPrime) {
    double Power = 1;
    double Square = NbrMod;
    double Modulus = currentPrime;
    while (Expon != 0) {
        if ((Expon & 1) == 1) {
            Power *= Square;
            Power -= std::floor(Power / Modulus)*Modulus;
        }
        Square *= Square;
        Square -= std::floor(Square / Modulus)*Modulus;
        Expon >>= 1;
    }
    return (int)Power;
}

/* get modular inverse of num wrt mod*/
void ModInvZnum(const Znum &num, Znum &inv, const Znum &mod) {
    auto rv = mpz_invert(ZT(inv), ZT(num), ZT(mod));
    assert(rv != 0);
}

/* returns true iff value is zero*/
bool BigNbrIsZero(const limb value[], int NumLen) {
    int ctr;
    for (ctr = 0; ctr < NumLen; ctr++) {
        if (value[ctr] != 0) {
            return false;  // Number is not zero.
        }
    }
    return true;      // Number is zero
}

/* convert value in most significant limb of Limb to a double (floating point)
note that ptrLimb points AFTER last valid value in limbs.
up to 3 most significant limbs are used. */
double getMantissa(const limb *ptrLimb, int nbrLimbs) {
    double dN = (double)*(ptrLimb - 1);
    double dInvLimb = 1 / (double)LIMB_RANGE;
    if (nbrLimbs > 1) {
        dN += (double)*(ptrLimb - 2) * dInvLimb;
    }
    if (nbrLimbs > 2) {
        dN += (double)*(ptrLimb - 3) * dInvLimb * dInvLimb;
    }
    return dN;
}

void LimbstoZ(const limb number[], Znum& numberZ, int NumLen) {
    numberZ = 0;
    for (int i = NumLen - 1; i >= 0; i--) {
        mpz_mul_2exp(ZT(numberZ), ZT(numberZ), BITS_PER_GROUP);  // shift numberZ left
        numberZ += number[i];      // add next limb
    }
}

int ZtoLimbs(limb *number, Znum numberZ, int NumLen) {
    // note: numberZ is a copy of the original. Its value is changed
    bool neg = false;
    Znum remainder;

    if (numberZ < 0) {
        neg = true;
        numberZ = -numberZ;  // make numberZ +ve
    }
    int i = 0;
    while (numberZ > 0) {
        if (i >= MAX_LEN || NumLen > MAX_LEN) {
            // number too big to convert.
            StackTrace2();
            std::string line = std::to_string(__LINE__);
            std::string mesg = "number too big : cannot convert to limbs: ";
            mesg += __func__;
            mesg += " line ";  mesg += line;
            mesg += " in file "; mesg += __FILE__;
            throw std::range_error(mesg);
        }
        //mpz_fdiv_qr_ui(ZT(quot), ZT(remainder), ZT(numberZ), LIMB_RANGE);
        /* calculating quotient and remainder separately turns
        out to be faster */
        mpz_fdiv_r_2exp(ZT(remainder), ZT(numberZ), BITS_PER_GROUP);
        number[i] = (int)MulPrToLong(remainder);
        mpz_fdiv_q_2exp(ZT(numberZ), ZT(numberZ), BITS_PER_GROUP);

        i++;
    }
    if (i < NumLen) {
        /* set any extra limbs to zero */
        std::memset(number + i, 0, (NumLen - i) * sizeof(limb));
    }
    if (neg) {
        ChSignBigNbr((int *)number, i + 1);
    }
    return i;
}

/* number = numberZ*/
int ZtoBigNbr(int number[], Znum numberZ) {
    // note: numberZ is a copy of the original. Its value is changed
    bool neg = false;
    Znum quot, remainder;

    if (numberZ < 0) {
        neg = true;
        numberZ = -numberZ;  // make numberZ +ve
    }
    int i = 0;
    while (numberZ > 0) {
        mpz_fdiv_qr_ui(ZT(quot), ZT(remainder), ZT(numberZ), LIMB_RANGE);
        //number[i] = (int)MulPrToLong(remainder);
        number[i] = (int)mpz_get_si(ZT(remainder));  // faster?? - no possibility of overflow here
        numberZ = quot;
        i++;
    }

    if (neg) {
        ChSignBigNbr(number, i);
    }
    return i;
}

/* convert integer list to Znum. */
//void ValuestoZ(Znum &numberZ, const int number[], int NumLen) {
//	numberZ = 0;
//	for (int i = NumLen - 1; i >= 0; i--) {
//		//numberZ *= LIMB_RANGE;
//		mpz_mul_2exp(ZT(numberZ), ZT(numberZ), BITS_PER_GROUP);  // shift numberZ left
//		numberZ += number[i];
//	}
//}


/* get log of Znum in base e */
double logZnum(const Znum &BigInt) {
    double BigId;
#ifdef __MPIR_VERSION
    long BiExp;   // changed for MPIR version 3.0.0
#else
    long BiExp;
#endif
    BigId = mpz_get_d_2exp(&BiExp, ZT(BigInt)); // BigId * 2^BiExp = BigInt 
    double logval = log(BigId) + BiExp * log(2);
    return logval;
}
/* get log of Znum in base b */
double logZnum(const Znum &BigInt, unsigned long long b) {
    double ln = logZnum(BigInt);   // get log in base e
    double logb = ln / std::log(b);  // get log in base b
    return logb;
}

/* returns nbrMod^Expon%currentPrime.
overflow could occur if currentPrime > 2^31
the alternative is to use modPowerLL */
static long long intModPow(long long NbrMod, long long Expon, long long currentPrime)
{
    unsigned long long power = 1;
    unsigned long long square = (unsigned long long)NbrMod;
    while (Expon != 0)
    {
        if ((Expon & 1) == 1)
        {
            power = (power * square) % (unsigned long long)currentPrime;
        }
        square = (square * square) % (unsigned long long)currentPrime;
        Expon >>= 1;
    }
    return (long long)power;
}

// This routine checks whether the number factor is a perfect power. 
// If it is not, it returns exponent 1. If it is a perfect power, it returns the  
// exponent and  the base such that base^exponent = factor.
long long PowerCheck(const Znum &factor, Znum &Base, long long upperBound) {
    /* upperbound is the largest number already tested as a factor by trial division
    i.e. factor has no factors < upperBound. This can be used to put a much
    smaller limit on maxExpon (max about 2000) */

    if (upperBound < 2)
        upperBound = 2;

    /* upperBound^maxExpon ≈ factor */
    unsigned long long maxExpon = (unsigned long long) (std::ceil(logZnum(factor) / std::log(upperBound)));

    int h;
    long long modulus, Exponent;
    unsigned long long maxPrime, j;
    int prime2310x1[] =
    { 2311, 4621, 9241, 11551, 18481, 25411, 32341, 34651, 43891, 50821 };
    // Primes of the form 2310x+1.
    bool expon2 = true, expon3 = true, expon5 = true;
    bool expon7 = true, expon11 = true;
    Znum Zmod, Root;
    std::vector<bool>ProcessExpon(maxExpon + 1);

    for (h = 0; h < sizeof(prime2310x1) / sizeof(prime2310x1[0]); h++) {
        int testprime = prime2310x1[h];
        // Zmod = mod = Bigint%testprime
        auto mod = mpz_mod_ui(ZT(Zmod), ZT(factor), testprime); // getRemainder(factor, testprime);
        if (expon2 && intModPow(mod, testprime / 2, testprime) > 1) {
            expon2 = false;
        }
        if (expon3 && intModPow(mod, testprime / 3, testprime) > 1) {
            expon3 = false;
        }
        if (expon5 && intModPow(mod, testprime / 5, testprime) > 1) {
            expon5 = false;
        }
        if (expon7 && intModPow(mod, testprime / 7, testprime) > 1) {
            expon7 = false;
        }
        if (expon11 && intModPow(mod, testprime / 11, testprime) > 1) {
            expon11 = false;
        }
    }

    maxPrime = 2 * maxExpon + 3;
    if (maxPrime > primeListMax) {
        StackTrace2();
        std::string line = std::to_string(__LINE__);
        std::string mesg = "number too big : cannot generate prime list. function : ";
        mesg += __func__;
        mesg += " line ";  mesg += line;
        mesg += " in file "; mesg += __FILE__;
        throw std::range_error(mesg);
    }

    for (h = 2; h <= maxExpon; h++) {
        ProcessExpon[h] = true;
    }

    for (size_t ix = 5, h = primeList[ix]; h < maxPrime / 2; ix++, h = primeList[ix]) {
        int processed = 0;
        for (j = 2 * h + 1; j < maxPrime; j += 2 * h) {
            if (isPrime2(j)) {
                modulus = mpz_mod_ui(ZT(Zmod), ZT(factor), j); // getRemainder(factor, j);
                if (intModPow(modulus, j / h, j) > 1) {
                    for (j = h; j <= maxExpon; j += h) {
                        ProcessExpon[j] = false;
                    }
                    break;
                }
            }
            if (++processed > 10) {
                break;
            }
        }
    }

    /* check possible exponent values. Note that largest found exponent value
    is returned although, unless this value is prime, any divisor of this
    value is also a valid exponent. */
    for (Exponent = maxExpon; Exponent >= 2; Exponent--) {
        if (Exponent % 2 == 0 && !expon2) {
            continue; // Not a square
        }
        if (Exponent % 3 == 0 && !expon3) {
            continue; // Not a cube
        }
        if (Exponent % 5 == 0 && !expon5) {
            continue; // Not a fifth power
        }
        if (Exponent % 7 == 0 && !expon7) {
            continue; // Not a 7th power
        }
        if (Exponent % 11 == 0 && !expon11) {
            continue; // Not an 11th power
        }
        if (!ProcessExpon[Exponent]) {
            continue;
        }
        if (mpz_root(ZT(Root), ZT(factor), Exponent) != 0) {
            Base = Root;   // factor is a perfect power
            return Exponent;
        }
    }

    Base = factor;   // not perfect power
    return 1;
}


// Calculate Jacobi symbol by following algorithm 2.3.5 of C&P book.
long long JacobiSymbol(long long upper, long long lower)
{
    long long tmp;
    long long a = upper % lower;
    long long m = lower;
    long long t = 1;
    while (a != 0)
    {
        while ((a & 1) == 0)
        {     // a is even.
            a >>= 1;
            if ((m & 7) == 3 || (m & 7) == 5)
            {   // m = 3 or m = 5 (mod 8)
                t = -t;
            }
        }
        tmp = a; a = m; m = tmp;   // Exchange a and m.
        if ((a & m & 3) == 3)
        {   // a = 3 and m = 3 (mod 4)
            t = -t;
        }
        a = a % m;
    }
    if (m == 1 || m == -1)
    {
        return t;
    }
    return 0;
}


// Find largest power of 2 that divides the number. Divide number
// by that power and return value of power in ShRight.
void DivideZnumByMaxPowerOf2(int &ShRight, Znum &number) {
    Znum two = 2;
    auto shift = mpz_remove(ZT(number), ZT(number), ZT(two));
    ShRight = (int)shift;
}



/* BPSW primality test:
1) If the input number is 2-SPRP composite, indicate composite and go out.
2) Perform a BPSW probable prime test on n 
   If n is not a probable prime, then n is composite.
   Otherwise, n is almost certainly prime.
   upperBound is the largest divisor already tested during trial division
Output: 0 = probable prime.
        1 = composite: not 2-Fermat pseudoprime.
        2 = composite: does not pass 2-SPRP test.
        3 = composite: does not pass BPSW test, but passes other tests.
return value 2 or 3 indicates a pseudoprime e.g. carmichael number
***********************************************************************/
int PrimalityTest(const Znum &Value, long long upperBound) {
    int i, ctr;
    Znum Mult1, Mult3, Mult4;
    const Znum two = 2;

    if (Value <= 2) {
        return 0;    // Indicate prime.
    }
    if (isEven(Value)) {
        return 1;    // Number is even and different from 2. Indicate composite.
    }

    // Perform base 2 Strong Probable Prime test (2-SPRP)
    //see https://en.wikipedia.org/wiki/Probable_prime

    Mult3 = Value - 1;
    DivideZnumByMaxPowerOf2(ctr, Mult3);  // Mult3 /= 2^ctr
    mpz_powm(ZT(Mult1), ZT(two), ZT(Mult3), ZT(Value));  // Mult1 = 2^Mult3 (mod Value)

    if (Mult1 != 1 && Mult1 != Value - 1) {
        for (i = 0; i < ctr; i++)
        {               // Loop that squares number.
            mpz_powm(ZT(Mult4), ZT(Mult1), ZT(two), ZT(Value));
            if (Mult4 == 1)
            {  // Current value is 1 but previous value is not 1 or -1: composite
                return 2;       // Composite. Not 2-strong probable prime.
            }
            if (Mult4 == Value - 1) {
                i = -1;         // Number is 2-strong pseudoprime.
                break;
            }
            Mult1 = Mult4;
        }
        if (i == ctr) {
            return 1;         // composite. Not 2-Fermat probable prime.
        }
        if (i != -1) {
            return 2;         // Composite. Not 2-strong probable prime.
        }
    }

    /* Number passes 2-SPRP test. Now perform a further test to see whether it is
    really prime or not.
    consensus is the BPSW (Baillie-Pomerance-Selfridge-Wagstaff) is better 
    than Miller-Rabin because;
    1. There are no known BPSW pseudoprimes
    2. It's faster than Miller-Rabin */

    int rv2 = mpz_bpsw_prp(ZT(Value));
    if (rv2 >= 1)
        return 0;  // 'almost' certainly prime;
    else if(rv2 == 0)
        return 3;  // composite - fails BPSW test

    /* BPSW returned error code */
    std::cout << "*** BPSW primality test returned error! *** \n";

    /* Miller-Rabin test used as backup */
#ifdef __MPIR_VERSION
    static bool first = true;
    static gmp_randstate_t rstate;
    if (first) {
        gmp_randinit_default(rstate);
        first = false;
    }

    auto rv = mpz_probable_prime_p(ZT(Value), rstate, 20, upperBound);
#else
    auto rv = mpz_probab_prime_p(ZT(Value), 16);
#endif
    if (rv == 0)
        return 3;			// composite - fails Miller-Rabin test
    else
        return 0;			// probably prime;
}

