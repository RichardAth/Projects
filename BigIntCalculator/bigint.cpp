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
#pragma fenv_access(on)

#include "pch.h"
#include <cfenv>

static BigInteger Base;     // work area
static limb approxInv[MAX_LEN];
static limb adjustedArgument[MAX_LEN];
static limb arrAux[MAX_LEN];
static int bitLengthCycle[20];

static int getNbrLimbs(const limb* bigNbr)
{
    const limb* ptrLimb = bigNbr + NumberLength;
    while (ptrLimb > bigNbr)
    {
        ptrLimb--;
        if (*ptrLimb != 0)
        {
            return (int)(ptrLimb - bigNbr + 1);
        }
    }
    return 1;
}

/*return addend1 + addend2 (used to overload + operator) */
BigInteger BigIntAdd(const BigInteger &Addend1, const BigInteger &Addend2) {
    int ctr, nbrLimbs;
    const limb *ptrBiggerAdd, *ptrSmallerAdd;
    limb *ptrSum;
    bool A1Smaller = false;
    BigInteger Sum;  // temporary variable 

    if (Addend1.nbrLimbs < Addend2.nbrLimbs) {
        A1Smaller = true;
        /* the absolute value of addend1 is less than the absolute value of addend2.*/
    }

    else if (Addend1.nbrLimbs == Addend2.nbrLimbs) {
        for (ctr = Addend1.nbrLimbs - 1; ctr >= 0; ctr--) {
            if (Addend1.limbs[ctr] != Addend2.limbs[ctr]) {
                break;
            }
        }
        if (ctr >= 0 && Addend1.limbs[ctr] < Addend2.limbs[ctr]) {
            /* the absolute value of addend1 is less than the absolute value of addend2.*/
            A1Smaller = true;
        }
    }
    if (A1Smaller) {
        nbrLimbs = Addend1.nbrLimbs;
        ptrBiggerAdd = Addend2.limbs;
        ptrSmallerAdd = Addend1.limbs;
    }
    else {
        // the absolute value of addend1 is >= the absolute value of addend2.
        nbrLimbs = Addend2.nbrLimbs;
        ptrBiggerAdd = Addend1.limbs;
        ptrSmallerAdd = Addend2.limbs;
    }
    ptrSum = Sum.limbs;

    if (Addend1.sign == Addend2.sign)
    {             // Both addends have the same sign. Sum their absolute values.
        unsigned int carry = 0;
        for (ctr = 0; ctr < nbrLimbs; ctr++) {
            carry = (carry >> BITS_PER_GROUP) + (unsigned int)ptrBiggerAdd[ctr] +
                (unsigned int)ptrSmallerAdd[ctr];
            ptrSum[ctr] = (int)(carry & MAX_INT_NBR);
        }
        if (A1Smaller)
            nbrLimbs = Addend2.nbrLimbs;
        else
            nbrLimbs = Addend1.nbrLimbs;
        for (; ctr < nbrLimbs; ctr++) {
            carry = (carry >> BITS_PER_GROUP) + (unsigned int)ptrBiggerAdd[ctr];
            ptrSum[ctr] = (int)(carry & MAX_INT_NBR);
        }
        if (carry >= LIMB_RANGE)
        {
            ptrSum[ctr] = 1;
            nbrLimbs++;
        }
    }
    else {    // addends have different signs. Subtract their absolute values.
        int borrow = 0;
        for (ctr = 0; ctr < nbrLimbs; ctr++) {
            borrow = (borrow >> BITS_PER_GROUP) + ptrBiggerAdd[ctr] - ptrSmallerAdd[ctr];
            ptrSum[ctr] = borrow & MAX_INT_NBR;
        }
        if (A1Smaller)
            nbrLimbs = Addend2.nbrLimbs;
        else
            nbrLimbs = Addend1.nbrLimbs;
        for (; ctr < nbrLimbs; ctr++) {
            borrow = (borrow >> BITS_PER_GROUP) + ptrBiggerAdd[ctr];
            ptrSum[ctr] = borrow & MAX_INT_NBR;
        }
    }

    while (nbrLimbs > 1 && Sum.limbs[nbrLimbs - 1] == 0) {
        nbrLimbs--;  // delete leading zeros.
    }

    Sum.nbrLimbs = nbrLimbs;
    if (A1Smaller)
        Sum.sign = Addend2.sign;  // use sign of addend with larger absolute value
    else
        Sum.sign = Addend1.sign;

    if (Sum.nbrLimbs == 1 && Sum.limbs[0] == 0) {
        Sum.sign = SIGN_POSITIVE; // Result is zero, so sign is +ve
    }
    return Sum;
}


/* add addend to big number */
static void addToAbsValue(limb* pLimbs, int* pNbrLimbs, int addend) {
    int ctr;
    int nbrLimbs = *pNbrLimbs;
    pLimbs[0] += addend;
    if ((unsigned int)pLimbs[0] < LIMB_RANGE)
    {     // No overflow. Go out of routine.
        return;
    }
    pLimbs[0] -= LIMB_RANGE;
    for (ctr = 1; ctr < nbrLimbs; ctr++)
    {
        //pLimbs++;        // Point to next most significant limb.
        if (pLimbs[ctr] != MAX_INT_NBR)
        {   // No overflow. Go out of routine.
            pLimbs[ctr]++;   // Add carry.
            return;
        }
        pLimbs[ctr] = 0;
    }
    (*pNbrLimbs)++;        // Result has an extra limb.
    pLimbs[ctr] = 1;   // Most significant limb must be 1.
}

/* subtract from big number */
static void subtFromAbsValue(limb* pLimbs, int* pNbrLimbs, int subt) {
    int nbrLimbs = *pNbrLimbs;
    pLimbs[0] -= subt;
    if (pLimbs[0] < 0)
    {
        int ctr;
        for (ctr = 1; ctr < nbrLimbs; ctr++)
        {
            pLimbs[ctr - 1] += LIMB_RANGE;
            if (--(pLimbs[ctr]) >= 0)
            {
                break;
            }
        }
        if (nbrLimbs > 1 && pLimbs[nbrLimbs - 1] == 0)
        {
            nbrLimbs--;
        }
    }
    *pNbrLimbs = nbrLimbs;
}

/* result += addend (used for operator overloading) */
void addbigint(BigInteger& Result, int addend) {
    int nbrLimbs = Result.nbrLimbs;
    limb* pResultLimbs = Result.limbs;
    auto sign = Result.sign;

    if (addend < 0) {
        // reverse signs of addend and result
        addend = -addend;  /* warning! if int = INT_MIN this won't work; an overflow will occur */
        if (sign == SIGN_POSITIVE) {
            sign = SIGN_NEGATIVE;  // addend and result have opposite signs
        }
        else {
            sign = SIGN_POSITIVE;   // addend and result have same sign
        }
    }

    if (sign == SIGN_POSITIVE) {   // Add abs(addend) to absolute value of pResult.
        addToAbsValue(pResultLimbs, &nbrLimbs, addend);
    }
    else {  // Subtract abs(addend) from absolute value of pResult.
        if (nbrLimbs == 1) {
            pResultLimbs[0] -= addend;
            if (pResultLimbs[0] < 0) {
                pResultLimbs[0] = -pResultLimbs[0];  // reverse sign of result
                BigIntNegate(Result);
            }
        }
        else {     // More than one limb.
            subtFromAbsValue(pResultLimbs, &nbrLimbs, addend);
        }
    }
    Result.nbrLimbs = nbrLimbs;
}


/* Dest = -Dest */
void BigIntNegate (BigInteger &pDest) {
    
    if (pDest.sign == SIGN_POSITIVE && (pDest.nbrLimbs != 1 || pDest.limbs[0] != 0))  {
        pDest.sign = SIGN_NEGATIVE;  // pDest > 0, now < 0
    }
    else {
        pDest.sign = SIGN_POSITIVE; // pDest <=0, now >= 0
    }
}


/* return minuend - subtrahend (used to overload - operator) */
BigInteger BigIntSubt(const BigInteger &Minuend, const BigInteger &Subtrahend) {
    BigInteger Difference;
    static BigInteger temp;
    temp = Subtrahend;   // copy Subtrahend to temporary variable
    BigIntNegate(temp);
    Difference = Minuend + temp;
    return Difference;
}

/* i = (i-subt)/divisor. Assume (without checking) that divisor > 0.
Does not appear to handle -ve divisor */
void subtractdivide(BigInteger& i, int subt, int divisor)
{
    int nbrLimbs = i.nbrLimbs;
    int remainder = 0;

#if 0
    char* ptrOutput = output;
    *ptrOutput++ = '2';
    *ptrOutput++ = '(';
    int2dec(&ptrOutput, i->sign);
    *ptrOutput++ = ',';
    *ptrOutput++ = ' ';
    int2dec(&ptrOutput, i->nbrLimbs);
    *ptrOutput++ = ';';
    *ptrOutput++ = ' ';
    int2dec(&ptrOutput, i->limbs[0]);
    *ptrOutput++ = ',';
    *ptrOutput++ = ' ';
    int2dec(&ptrOutput, i->limbs[1]);
    *ptrOutput++ = ',';
    *ptrOutput++ = ' ';
    *ptrOutput++ = ')';
    *ptrOutput++ = ',';
    *ptrOutput++ = ' ';
    int2dec(&ptrOutput, subt);
    *ptrOutput++ = ',';
    *ptrOutput++ = ' ';
    int2dec(&ptrOutput, divisor);
    //  databack(output);
    if ((unsigned int)i->limbs[0] >= LIMB_RANGE)
    {
        remainder = 1;
    }
#endif
    /* if subt is not zero subtract it from i*/
    if (subt > 0) {
        if (i >= 0) {       // Subtract subt from absolute value.
            subtFromAbsValue(i.limbs, &nbrLimbs, subt);
        }
        else {               // Add subt to absolute value.
            addToAbsValue(i.limbs, &nbrLimbs, subt);
        }
    }
    else if (subt < 0) {  // subt < 0
        if (i >= 0) {    // Subtract subt from absolute value.
            addToAbsValue(i.limbs, &nbrLimbs, -subt);
        }
        else {               // Add subt to absolute value.
            subtFromAbsValue(i.limbs, &nbrLimbs, -subt);
        }
    }

    // Divide number by divisor.
    i = BigIntDivideInt(i, divisor);
}

/* returns Factor1 * Factor2 (used to overload * operator)
Factor1 will be expanded to the length of Factor2 or vice versa 
Changed design from returning error code to throw an exception, because error 
codes need to be tested for after every call, passed up to the calling routine,
tested for again and so on up to the top level. In practise this wasn't done consistently. */
BigInteger BigIntMultiply(const BigInteger &Factor1, const BigInteger &Factor2)
{
    int nbrLimbsFactor1 = Factor1.nbrLimbs;
    int nbrLimbsFactor2 = Factor2.nbrLimbs;
    int nbrLimbs;
    BigInteger Product = 0;                  // temporary variable 
    limb Prodl[MAX_LEN * 2] = { 0 };    // temporary variable 

    if (Factor1 == 0 || Factor2 == 0) {    // one or both factors are zero.
        Product = 0;                       // Product is zero.
        return Product;
    }

    if (Factor1.nbrLimbs + Factor2.nbrLimbs > MAX_LEN )  // approx 23000 digits
    {
        // product out of range; throw exception
        StackTrace2();
        std::string line = std::to_string(__LINE__);
        std::string mesg = "cannot multiply: product out of range ";
        mesg += __func__;
        mesg += " line ";  mesg += line;
        mesg += " in file "; mesg += __FILE__;
        throw std::range_error(mesg);
    }

    /* Factor1 will be expanded to the length of Factor2 or vice versa. It is necessary to
    override the const-ness, but the value represented does not change */
    if (nbrLimbsFactor1 < nbrLimbsFactor2) {
        std::memset(&((BigInteger&)Factor1).limbs[nbrLimbsFactor1], 0, (nbrLimbsFactor2 - nbrLimbsFactor1) * sizeof(limb));
        nbrLimbs = nbrLimbsFactor2;
    }
    else {
        std::memset(&((BigInteger&)Factor2).limbs[nbrLimbsFactor2], 0, (nbrLimbsFactor1 - nbrLimbsFactor2) * sizeof(limb));
        nbrLimbs = nbrLimbsFactor1;
    }

    multiply(&Factor1.limbs[0], &Factor2.limbs[0], Prodl, nbrLimbs, NULL);
    nbrLimbs = nbrLimbsFactor1 + nbrLimbsFactor2;
    while (Prodl[nbrLimbs - 1] == 0)
        nbrLimbs--;     // remove leading zeros

    if (nbrLimbs > MAX_LEN)  // limit applied earlier is probably lower, this is just insurance
    {
        StackTrace2();
        std::string line = std::to_string(__LINE__);
        std::string mesg = "cannot multiply: product out of range ";
        mesg += __func__;
        mesg += " line ";  mesg += line;
        mesg += " in file "; mesg += __FILE__;
        throw std::range_error(mesg);
    }
    LimbsToBigInteger(Prodl, Product, nbrLimbs);
    if (Product.limbs[nbrLimbs - 1] == 0) {
        nbrLimbs--;  
    }
    Product.nbrLimbs = nbrLimbs;

    /* set sign of product */
    if (nbrLimbs == 1 && Product.limbs[0] == 0) {
        Product.sign = SIGN_POSITIVE;  // product is zero
    }
    else {
        if (Factor1.sign == Factor2.sign) {
            Product.sign = SIGN_POSITIVE;
        }
        else {
            Product.sign = SIGN_NEGATIVE;
        }
    }
    return Product;
}

// m *= n (used for operator overloading)
void MultBigNbrByInt(BigInteger &m, int n) {
    long long prod, carry = 0;
    int i;
    bool pos = true;

    if (n < 0) {
        pos = false;
        n = -n;   // get abs value of n
    }
    if (n == 0 || m == 0) {
        m = 0;
        return;
    }
    for (i = 0; i < m.nbrLimbs; i++) {
        prod = carry + (long long)m.limbs[i] * n;
        carry = prod >> BITS_PER_GROUP;
        m.limbs[i] = prod & MAX_VALUE_LIMB;
    }
    if (carry != 0) {
        if (i >= MAX_LEN) {
            ThrowExc("BigInteger exceeds max size");
        }
        m.limbs[i] = (int)carry;
        m.nbrLimbs++;
    }
    if (!pos) {  // if n is -ve flip sign of product 
        if (m.sign = SIGN_POSITIVE)
            m.sign = SIGN_NEGATIVE;
        else
            m.sign = SIGN_POSITIVE;
    }
}

/* calculate Dividend mod Divisor (used for operator overloading) 
uses global variable Base */
BigInteger BigIntRemainder(const BigInteger &Dividend, const BigInteger &Divisor) {
    BigInteger Remainder;  
    if (Divisor == 0) {   // If divisor = 0, then remainder is the dividend.
        Remainder = Dividend;
        return Remainder;
    }
    if (Divisor.nbrLimbs == 1) {
        /* use quick method if divisor is small */
        Remainder = getRemainder(Dividend, Divisor.limbs[0]);
        return Remainder;
    }
    Base = Dividend / Divisor;    // Get quotient of division.
    Base = Base*Divisor;   /* use long multiplication if divisor has more than 1 limb. */
    return Dividend - Base;
}

/* calculate base^exponent. Throws exception if Power would exceed max size */
void BigIntPowerIntExp(const BigInteger &base, int exponent, BigInteger &Power) {
    int mask;
    if (base == 0) {     // base = 0 -> power = 0
        Power = 0;
        return;       // base = 0, so result is zero
    }
    if (exponent == 0) {
        Power = 1LL;
        return;
    }
    Power = 1;
    for (mask = 1 << 30; mask != 0; mask >>= 1) {
        /* first look for most significant bit of exponent */
        if ((exponent & mask) != 0) {
            for (; mask != 0; mask >>= 1) {
                Power *= Power;  
                if ((exponent & mask) != 0) {
                    Power *= base; 
                }
            }
            break;  // we are finished
        }
    }
    return;
}
 
/* BigInt = e^logar.  This is the inverse function of LogBigInt. As the result
is an integer it is not exact, and also for large values only at best the 1st
15 significant digits of the result are accurate. */
void expBigInt(BigInteger &bigInt, double logar) {
    std::fenv_t envp ;
    unsigned int control_word; /* save current state of fp control word here */
    errno_t err;

    err = std::fegetenv(&envp);  /* save current floating point environment */
    err = std::feclearexcept(FE_OVERFLOW);
    /* trap hardware FP exceptions except inexact and underflow which are
    considered to be normal, not errors, and overflow which is handled by using an 
    alternative method */
    err = _controlfp_s(&control_word, _EM_INEXACT | _EM_UNDERFLOW | _EM_OVERFLOW, MCW_EM);
    if (err) {
        printf_s("could not set FP control word\n");
        std::exit(EXIT_FAILURE);
    }

    double e = std::exp(logar);  /* get e^logar directly */
    auto rv = std::fpclassify(e);  // check for overflow
    
    if (rv != FP_INFINITE && rv != FP_NAN) {
        bigInt = e;  /* note: the assignment statement uses DoubleToBigInt to
                    convert floating point to BigInteger */
    }
    else {
        /* simple approach didn't work, but more complicated method below handles
         much larger numbers. The strategy is to split the calculation into two parts.
         One part uses Big Integers only and so avoids floating point overflow. 
         The other part uses floating point, but the maximum number size is limited. 
         Mathematically this method is sound, but in practise it is less accurate. */
    
        BigInteger base=1;   
        base <<= 128;        //  new base for logs = 2^128
        double logb = logar / (std::log(2.0)*128);  // convert log to new base
        double intpf;   // integer part of log
        double fracpf = std::modf(logb, &intpf);  // split into log integer and fraction
        /* the desired result = base^(intpf+fracpf)
                 = base^intpf * base^fracpf */
        int intp = (int)std::round(intpf);         // exact conversion
        BigIntPowerIntExp(base, intp, bigInt);  /* bigInt = base^intpf*/

        /* by using a large base we make frac large, so that conversion to integer 
        doesn't introduce errors */
        double frac = std::pow(2, fracpf*128);  //  frac = base^fracpf
        BigInteger fracBI;
        fracBI = frac;    // convert to integer      
        bigInt = bigInt * fracBI;
    }
    std::feclearexcept(FE_OVERFLOW);
    std::fesetenv(&envp);  /* restore previous FP environment */
}


/* convert double dvalue to bigInt. Conversion is only accurate to about 15 
significant digits. Used for operator overloading. */
void DoubleToBigInt(BigInteger &bigInt, double dvalue) {

    if (dvalue - 0.5 > LLONG_MIN && dvalue + 0.5 < LLONG_MAX) {
        long long vv = std::llround(dvalue); // convert directly to long long if possible
        bigInt = vv;
        return;
    }
    Znum temp;
    mpz_set_d(ZT(temp), dvalue);  // convert double to Znum
    // this method has been tested to be at least as accurate as direct conversion,
    // and obviously it's easier to use standard library functions.
    ZtoBig(bigInt, temp);        // convert Znum to BigInt
    return;
}

/* estimate natural log of BigInt. */
double logBigInt (const BigInteger &pBigInt) {
    int nbrLimbs;
    double logar;
    nbrLimbs = pBigInt.nbrLimbs;
    if (nbrLimbs == 1) {
        logar = std::log((double)(pBigInt.limbs[0]));
    }
    else {
        double value = pBigInt.limbs[nbrLimbs - 2] +
            (double)pBigInt.limbs[nbrLimbs - 1] * LIMB_RANGE;
        if (nbrLimbs == 2) {
            logar = std::log(value);
        }
        else {
            logar = std::log(value + (double)pBigInt.limbs[nbrLimbs - 3] / LIMB_RANGE);
        }
        logar += (double)((nbrLimbs - 2)*BITS_PER_GROUP)*std::log(2);
    }
    return logar;
}

/* divide by 2, use right shift for speed */
void BigIntDivide2(BigInteger &pArg) {
    int nbrLimbs = pArg.nbrLimbs;
    int ctr = nbrLimbs - 1;
    unsigned int carry;
    //limb *ptrLimb = &pArg->limbs[ctr];
    limb *ptrLimb = pArg.limbs;
    carry = 0;
    for (; ctr >= 0; ctr--)
    {
        carry = (carry << BITS_PER_GROUP) + (unsigned int)ptrLimb[ctr];
        ptrLimb[ctr] = (int)(carry >> 1);
        carry &= 1;
    }
    if (nbrLimbs > 1 && pArg.limbs[nbrLimbs - 1] == 0)
    {     // Most significant limb is zero, so reduce size by one limb.
        pArg.nbrLimbs--;
    }
}

/* arg = arg*2^power. Throw exception if product is too large. 
(shiftBI does the same. May be more efficient?)*/
static void BigIntMutiplyPower2(BigInteger &pArg, int power2)
{
    int ctr;
    int nbrLimbs = pArg.nbrLimbs;
    limb *ptrLimbs = pArg.limbs;

    for (; power2 > 0; power2--) {
        /*each time round the loop multiplies arg by 2 */
        unsigned int carry = 0;
        for (ctr = 0; ctr < nbrLimbs; ctr++)
        {
            carry += (unsigned int)ptrLimbs[ctr] << 1;
            ptrLimbs[ctr] = (int)(carry & MAX_VALUE_LIMB);
            carry >>= BITS_PER_GROUP;
        }
        if (carry != 0)
        {
            ptrLimbs[ctr] = (int)carry;
            nbrLimbs++;
            if (nbrLimbs > MAX_LEN) {
                StackTrace2();
                std::string line = std::to_string(__LINE__);
                std::string mesg = "number too big : cannot do multiplication : ";
                mesg += __func__;
                mesg += " line ";  mesg += line;
                mesg += " in file "; mesg += __FILE__;
                throw std::range_error(mesg);
            }
        }
    }
    pArg.nbrLimbs = nbrLimbs;
}

/* return true if Nbr1 == Nbr2 (used for operator overloading)*/
bool TestBigNbrEqual(const BigInteger &Nbr1, const BigInteger &Nbr2) {
    int ctr;
    auto N1Limbs = Nbr1.nbrLimbs;
    auto N2Limbs = Nbr2.nbrLimbs;
    while (N1Limbs > 1)
        if (Nbr1.limbs[N1Limbs - 1] == 0)
            N1Limbs--;
        else
            break;
    while (N2Limbs > 1)
        if (Nbr2.limbs[N2Limbs - 1] == 0)
            N2Limbs--;
        else
            break;

    if (N1Limbs != N2Limbs) {        
        return false;  // Sizes of numbers are different.
    }
    if (Nbr1.sign != Nbr2.sign) { 
           // Sign of numbers are different.
        if (N1Limbs == 1 && Nbr1.limbs[0] == 0 && Nbr2.limbs[0] == 0) {              
            return true; // Both numbers are zero.
        }
        return false; // differents signs, therefore cannot be equal
    }

    // Check whether both numbers are equal.
    for (ctr = N1Limbs - 1; ctr >= 0; ctr--) {
        if (Nbr1.limbs[ctr] != Nbr2.limbs[ctr]) {
            return false;  // Numbers are different.
        }
    }        
    return true;  // Numbers are equal.
}

/* return true if Nbr1 < Nbr2 (used for operator overloading) */
bool TestBigNbrLess(const BigInteger &Nbr1, const BigInteger &Nbr2) {
    int ctr;
    auto N1Limbs = Nbr1.nbrLimbs;
    auto N2Limbs = Nbr2.nbrLimbs;
    while (N1Limbs > 1)
        if (Nbr1.limbs[N1Limbs - 1] == 0)
            N1Limbs--;
        else
            break;
    while (N2Limbs > 1)
        if (Nbr2.limbs[N2Limbs - 1] == 0)
            N2Limbs--;
        else
            break;

    if (Nbr1.sign != Nbr2.sign) {
        // Sign of numbers are different.
        if (N1Limbs == 1 && Nbr1.limbs[0] == 0 && Nbr2.limbs[0] == 0) {
            return false; // Both numbers are zero i.e Nbr1 not less than Nbr2
        }
        else return (Nbr1.sign == SIGN_NEGATIVE); /* Nbr1 < 0 & Nbr2 >= 0 */
    }

    /* numbers have same sign */
    if (N1Limbs != N2Limbs) {
        /* length of numbers is different*/
        if (Nbr1.sign == SIGN_POSITIVE)
            return N1Limbs < N2Limbs;
        else
            return N1Limbs > N2Limbs;
    }

    // numbers have same sign and length. Check whether both numbers are equal.
    for (ctr = N1Limbs - 1; ctr >= 0; ctr--) {
        if (Nbr1.limbs[ctr] < Nbr2.limbs[ctr]) {
            return true;  // Nbr1 < Nbr2.
        }
        else if (Nbr1.limbs[ctr] > Nbr2.limbs[ctr]) {
            return false;  // Nbr1 > Nbr2.
        }
    }
    return false;  // Numbers are equal.
}

/* calculate GCD of arg1 & arg2. Use Base and Power as working space */
void BigIntGcd(const BigInteger &Arg1, const BigInteger &Arg2, BigInteger &Result)
{
    int nbrLimbs1 = Arg1.nbrLimbs;
    int nbrLimbs2 = Arg2.nbrLimbs;
    int power2;
    static BigInteger Power;
    if (Arg1 == 0)
    {               // First argument is zero, so the GCD is second argument.
        Result = Arg2;    //CopyBigInt(Result, pArg2);
        return;
    }
    if (Arg2 == 0)
    {               // Second argument is zero, so the GCD is first argument.
        Result = Arg1;		//CopyBigInt(Result, pArg1);
        return;
    }
    // Reuse Base and Power temporary variables.
    Base = Arg1;     // CopyBigInt(Base, Arg1);   
    Power = Arg2;   //  CopyBigInt(Power, Arg2); 
    Base.sign = SIGN_POSITIVE;
    Power.sign = SIGN_POSITIVE;
    power2 = 0;
    while (((Base.limbs[0] | Power.limbs[0]) & 1) == 0)
    {  // Both values are even
        BigIntDivide2(Base);
        BigIntDivide2(Power);
        power2++;
    }

    while (Base != Power)	//while (TestBigNbrEqual(Base, Power) == 0)
    {    // Main GCD loop.
        if (Base.isEven()) {     // Number is even. Divide it by 2.
            BigIntDivide2(Base);
            continue;
        }
        if (Power.isEven())  {     // Number is even. Divide it by 2.
            BigIntDivide2(Power);
            continue;
        }
        Result = Base - Power; // BigIntSubt(Base, Power, Result);
        if (Result >= 0) {
            Base = Result; // CopyBigInt(Base, Result);
            BigIntDivide2(Base);
        }
        else {
            Power = Result; // CopyBigInt(Power, Result);
            Power.sign = SIGN_POSITIVE;
            BigIntDivide2(Power);
        }
    }
    Result = Base; // CopyBigInt(Result, Base);
    BigIntMutiplyPower2(Result, power2); /* Result *= 2^power     */
}


/* creates a BigInteger from a list of values.   
uses global value NumberLength for number of ints. 1st entry in list is number of values
that follow */
//void IntsToBigInteger(/*@in@*/const int *ptrValues, /*@out@*/BigInteger &bigint)
//{
//	if (NumberLength > MAX_LEN || NumberLength < 0 || ptrValues[0] > MAX_LEN) {
//		std::string line = std::to_string(__LINE__);
//		std::string mesg = "number too big : cannot convert to BigInteger: ";
//		mesg += __func__;
//		mesg += " line ";  mesg += line;
//		mesg += " in file "; mesg += __FILE__;
//		throw std::range_error(mesg);
//	}
//	limb *destLimb = bigint.limbs;
//	bigint.sign = SIGN_POSITIVE;
//	if (NumberLength == 1) {
//		destLimb[0] = ptrValues[1];
//		bigint.nbrLimbs = 1;
//	}
//	else {
//		memcpy(destLimb, ptrValues+1, ptrValues[0] * sizeof(ptrValues[0]));
//		bigint.nbrLimbs = ptrValues[0];
//		if (NumberLength > ptrValues[0])  // clear most significant limbs to zero if required
//			memset(destLimb + ptrValues[0], 0, (NumberLength - ptrValues[0]) * sizeof(ptrValues[0]));
//	}
//}


/* creates a list of values from a BigInteger, 1st entry in list is number of 
values that follow. Also uses global value NumberLength for number of ints. */
//void BigIntegerToInts(/*@out@*/int *ptrValues, /*@in@*/const BigInteger &bigint) {
//	const limb *srcLimb = bigint.limbs;
//	if (NumberLength == 1) {
//		ptrValues[0] = 1;
//		ptrValues[1] = srcLimb[0];
//	}
//	else {
//		int nbrLimbs = getNbrLimbs(bigint.limbs); //nbrLimbs <= NumberLength
//		assert(nbrLimbs == bigint.nbrLimbs);
//		ptrValues[0] = nbrLimbs;
//		std::memcpy(ptrValues + 1, srcLimb, nbrLimbs * sizeof(ptrValues[0]));
//	}
//}

/* convert limbs to BigInteger. */
void LimbsToBigInteger(/*@in@*/const limb *ptrValues, 
    /*@out@*/BigInteger &bigint, int NumLen) {
    if (NumLen > MAX_LEN || NumLen < 0 ) {
        StackTrace2();
        std::string line = std::to_string(__LINE__);
        std::string mesg = "number too big : cannot convert to BigInteger: ";
        mesg += __func__;
        mesg += " line ";  mesg += line;
        mesg += " in file "; mesg += __FILE__;
        throw std::range_error(mesg);
    }
    std::memset(bigint.limbs, 0, MAX_LEN * sizeof(limb));  /* ensure unused limbs are zero */
    if (NumLen == 1) {
        bigint.limbs[0] = ptrValues[0];
        bigint.nbrLimbs = 1;
    }
    else { 
        std::memcpy(bigint.limbs, ptrValues, NumLen * sizeof(limb));

        int nbrLimbs;   // remove any leading zeros
        for (nbrLimbs = NumLen-1; nbrLimbs > 1; nbrLimbs--) {
            if (ptrValues[nbrLimbs] != 0) {
                break;
            }
        }
        bigint.nbrLimbs = nbrLimbs+1;
    }
    bigint.sign = SIGN_POSITIVE;
}

/* Convert BigInteger to limbs. uses global value NumberLength for number of limbs. */
void BigIntegerToLimbs(/*@out@*/limb ptrValues[],
    /*@in@*/const BigInteger &bigint, int NumLen)
{
    if (NumLen == 1) {
        ptrValues[0] = bigint.limbs[0];
        ptrValues[1] = 0;
    }
    else {
        int nbrLimbs = bigint.nbrLimbs; // use lesser of bigint.nbrLimbs & NumLen
        if (nbrLimbs >= NumLen) {
            std::memcpy(ptrValues, bigint.limbs, NumLen * sizeof(limb));
            ptrValues[NumLen] = 0;
        }
        else {
            std::memcpy(ptrValues, bigint.limbs, nbrLimbs * sizeof(limb));
            /* set any extra limbs to zero */
            std::memset(ptrValues + nbrLimbs, 0, (NumLen - nbrLimbs) * sizeof(limb));
        }
    }
}


/* divide number by 2 (right shift) until it is odd. Return changed value of number
and the number of 0-bits shifted out. */
void DivideBigNbrByMaxPowerOf2(int *pShRight, limb *number, int *pNbrLimbs)
{
    int power2 = 0;
    long long mask;
    int index, index2,  shRg;
    int nbrLimbs = *pNbrLimbs;
    // Start from least significant limb (number zero).
    for (index = 0; index < nbrLimbs; index++)
    {
        if (number[index] != 0)
        {
            break;
        }
        power2 += BITS_PER_GROUP;
    }
    for (mask = 0x1; mask <= MAX_VALUE_LIMB; mask *= 2)
    {
        if ((number[index] & mask) != 0)
        {
            break;
        }
        power2++;
    }
    // Divide number by this power.
    shRg = power2 % BITS_PER_GROUP; // Shift right bit counter
    if ((number[nbrLimbs - 1] & (-(1 << shRg))) != 0)
    {   // Most significant bits set.
        *pNbrLimbs = nbrLimbs - index;
    }
    else
    {   // Most significant bits not set.
        *pNbrLimbs = nbrLimbs - index - 1;
    }
    // Move number shRg bits to the right.
    mask = (1 << shRg) - 1;
    for (index2 = index; index2 < nbrLimbs; index2++)
    {
        if (index2 < nbrLimbs - 1)
            number[index2] = ((number[index2] >> shRg) |
                (number[index2 + 1] << (BITS_PER_GROUP - shRg))) &
            MAX_VALUE_LIMB;
        else
            number[index2] = (number[index2] >> shRg);
    }
    if (index > 0)
    {   // Move limbs to final position.
         std::memmove(number, &number[index], (nbrLimbs - index) * sizeof(limb));
    }
    *pShRight = power2;
}


/* convert Znum to BigInteger. Returns false if number is too big to convert.
this function is also used to overload the assignment operator */
bool ZtoBig(BigInteger &number, Znum numberZ) {
    number.nbrLimbs = 0;
    bool neg = false;
    Znum remainder;

    if (numberZ < 0) {
        neg = true;
        numberZ = -numberZ;  // make numberZ +ve
    }
    int i = 0;
    while (numberZ > 0) {
        if (i >= MAX_LEN) 
            return false;   // number too big to convert.
        //mpz_fdiv_qr_ui(ZT(quot), ZT(remainder), ZT(numberZ), LIMB_RANGE);
        /* calculating quotient and remainder separately turns
        out to be faster */
        mpz_fdiv_r_2exp(ZT(remainder), ZT(numberZ), BITS_PER_GROUP);
        number.limbs[i] = (int)MulPrToLong(remainder);
        mpz_fdiv_q_2exp(ZT(numberZ), ZT(numberZ), BITS_PER_GROUP);
        i++;
    }
    number.nbrLimbs = i;
    if (neg) {
        number.sign = SIGN_NEGATIVE;
        numberZ = -numberZ;  // put back original value in numberZ
    }
    else
        number.sign = SIGN_POSITIVE;

    return true;
}

/* convert BigInteger to Znum */
void BigtoZ(Znum &numberZ, const BigInteger &number) {

    numberZ = 0;
    for (int i = number.nbrLimbs - 1; i >= 0; i--) {
        //numberZ *= LIMB_RANGE;
        mpz_mul_2exp(ZT(numberZ), ZT(numberZ), BITS_PER_GROUP);  // shift numberZ left
        numberZ += number.limbs[i];      // add next limb
    }
    if (number.sign == SIGN_NEGATIVE)
        numberZ = -numberZ;
}

/* convert num to long long. If num > max, truncate it (right shift)
and exp is > 0 and represents the number of discarded bits. */
long long BigToLL(const BigInteger &num, int &exp) {
    BigInteger temp = num;
    if (num.nbrLimbs == 1) {
        exp = 0;
        if (num.sign == SIGN_POSITIVE)
            return num.limbs[0];
        else
            return -num.limbs[0];
    }
    exp = num.bitLength() - 63;                // number of bits to truncate
    if (exp < 0)
        exp = 0;
    temp >>= exp;
    long long result = temp.lldata();
    return result;
}


/* inverse of BigToLL. set num = LL *2^exp */
void LLToBig(BigInteger &num, long long LL, int exp) {
    num = LL;
    if (exp == 0) { 
        return;
    }
    assert(exp > 0);
    num <<= exp;
}

/* shift first left by the number of bits specified in shiftCtr. A -ve value
in shiftCtr causes a right shift. Used for operator overloading
Right Shifts simulate 2s complement arithmetic right shift.
Mathematically, the shift result is equivalent to result = first * 2^shiftCtr,
whether ShiftCtr is +ve or -ve. */
void shiftBI(const BigInteger &first, const int shiftCtr, BigInteger &result)
{
    int delta, rem, ctr;
    long long prevLimb, curLimb;
    int ptrDest, ptrSrc;
    bool shiftleft = true;
    if (shiftCtr > 0) {
        delta = shiftCtr / BITS_PER_GROUP;
        rem = shiftCtr % BITS_PER_GROUP;
    }
    else {
        delta = (-shiftCtr) / BITS_PER_GROUP;
        rem = (-shiftCtr) % BITS_PER_GROUP;
        shiftleft = false;
    }
    int nbrLimbs = first.nbrLimbs;

    if (shiftleft) {     // Perform shift left.

        if ((first.nbrLimbs + delta) >= MAX_LEN) {
            // Shift too much to the left; would cause overflow
            StackTrace2();
            std::string line = std::to_string(__LINE__);
            std::string mesg = "cannot shift left: result out of range ";
            mesg += __func__;
            mesg += " line ";  mesg += line;
            mesg += " in file "; mesg += __FILE__;
            throw std::range_error(mesg);
        }


        result.nbrLimbs = first.nbrLimbs + delta;
        result.sign = first.sign;
        prevLimb = 0;
        ptrSrc = nbrLimbs - 1;
        ptrDest = nbrLimbs + delta;

        for (ctr = nbrLimbs - 1; ctr >= 0; ctr--)
        {  // Process starting from most significant limb.
            curLimb = first.limbs[ctr];
            result.limbs[ptrDest] = ((curLimb >> (BITS_PER_GROUP - rem))
                | (prevLimb << rem)) & MAX_INT_NBR;
            ptrDest--;
            prevLimb = curLimb;
        }

        result.limbs[ptrDest] = (prevLimb << rem) & MAX_INT_NBR;
        if (delta > 0) {
            std::memset(result.limbs, 0, delta * sizeof(limb));
        }
        //result.nbrLimbs += delta;
        if (result.limbs[result.nbrLimbs] != 0) {
            result.nbrLimbs++;
        }
    }

    else {     // Perform shift right.
        int isNegative = 0;
        if (shiftCtr > first.nbrLimbs * BITS_PER_GROUP)
        {   // Shift too much to the right. Result is zero or -1.
            if (first.sign == SIGN_POSITIVE)
                result = 0;
            else
                result = -1;

            return;
        }

        result = first;
        if (first.sign == SIGN_NEGATIVE)
        {   // If it is negative, add 1, perform shift right, and finally subtract 1 from result.
            isNegative = 1;
            result++;     
        }
        // Shift right the absolute value.
        result.limbs[nbrLimbs] = 0;
        curLimb = result.limbs[delta];
        ptrDest = 0;

        for (ctr = delta; ctr <= nbrLimbs; ctr++)
        {  // Process starting from least significant limb.
            prevLimb = curLimb;
            curLimb = result.limbs[ctr + 1];
            result.limbs[ptrDest++] = ((prevLimb >> rem)
                | (curLimb << (BITS_PER_GROUP - rem))) & MAX_INT_NBR;
        }

        result.nbrLimbs -= delta;
        if (result.nbrLimbs == 0 || result.limbs[result.nbrLimbs] != 0) {
            result.nbrLimbs++;
        }
        if (isNegative) {    // Adjust negative number.
            result--;        
        }
    }
    return;
}


// All computations are done in little-endian notation.
// Find power of 2 that divides the number.
// output: pNbrLimbs = pointer to number of limbs
//         pPower2 = reference to power of 2.
static void MultiplyBigNbrByMinPowerOf2(int &pPower2, const limb *number, int len, limb *dest)
{
    limb mostSignficLimb, oldLimb, newLimb;
    int index2, mask, shLeft;
    limb *ptrDest;

    shLeft = 0;
    mostSignficLimb = *(number + len - 1);
    for (mask = LIMB_RANGE / 2; mask > 0; mask >>= 1)
    {
        if ((mostSignficLimb & mask) != 0)
        {
            break;
        }
        shLeft++;
    }
    ptrDest = dest;
    // Multiply number by this power.
    oldLimb = 0;
    for (index2 = len; index2 > 0; index2--)
    {
        newLimb = *ptrDest;
        *(ptrDest++) = ((newLimb << shLeft) |
            (oldLimb >> (BITS_PER_GROUP - shLeft))) & MAX_VALUE_LIMB;
        oldLimb = newLimb;
    }
    *ptrDest = oldLimb >> (BITS_PER_GROUP - shLeft);
    pPower2 = shLeft;
}

/* used for operator overloading */
// After computing the number of limbs of the results, this routine finds the inverse
// of the divisor and then multiplies it by the dividend using nbrLimbs+1 limbs.
// After that, the quotient is adjusted.
BigInteger BigIntDivide(const BigInteger &Dividend, const BigInteger &Divisor) {
    BigInteger Quotient;
    double inverse;
    limb oldLimb, newLimb;
    int nbrLimbs, nbrLimbsDividend, nbrLimbsDivisor;

    // Check whether the divisor is zero.
    if (Divisor == 0) {  // Indicate error if divisor is zero.
        StackTrace2();
        std::string line = std::to_string(__LINE__);
        std::string mesg = "cannot divide by zero: ";
        mesg += __func__;
        mesg += " line ";  mesg += line;
        mesg += " in file "; mesg += __FILE__;
        throw std::range_error(mesg);
    }
    // Get number of limbs of quotient.
    nbrLimbsDividend = Dividend.nbrLimbs;
    nbrLimbsDivisor = Divisor.nbrLimbs;
    nbrLimbs = nbrLimbsDividend - nbrLimbsDivisor;
    if (nbrLimbs < 0)
    {   // Absolute value of dividend is less than absolute value of divisor.
        Quotient = 0;
        return Quotient;
    }
    if (nbrLimbs == 0)
    {   // Both divisor and dividend have the same number of limbs.
        for (nbrLimbs = nbrLimbsDividend - 1; nbrLimbs > 0; nbrLimbs--) {
            if (Dividend.limbs[nbrLimbs] != Divisor.limbs[nbrLimbs]) {
                break;
            }
        }
        if (Dividend.limbs[nbrLimbs] < Divisor.limbs[nbrLimbs])
        {   //  abs value of Dividend is less than that of divisor, so quotient is zero.
            Quotient = 0;
            return Quotient;
        }
    }
    if (nbrLimbsDividend == 1) {   // If dividend is small, perform the division directly.
        Quotient.limbs[0] = Dividend.limbs[0] / Divisor.limbs[0];
        Quotient.nbrLimbs = 1;
    }
    else if (nbrLimbsDivisor == 1)
    {   // Divisor is small: use divide by int.
        // Sign of quotient is determined later.
         Quotient = Dividend/Divisor.limbs[0];
    }
    else {
        int index;
        int bitLength;
        int bitLengthNbrCycles;
        int idx;
        int nbrLimbsQuotient;
        int power2;
        limb *ptrDest, *ptrQuotient, *ptrQuot;
        const limb *ptrDivisor, *ptrDividend;

        nbrLimbs += 3;    // Use this number of limbs for intermediate calculations.
        if (nbrLimbs > nbrLimbsDivisor) {
            std::memset(&adjustedArgument[0], 0, (nbrLimbs - nbrLimbsDivisor) * sizeof(limb));
            std::memcpy(&adjustedArgument[nbrLimbs - nbrLimbsDivisor], &Divisor.limbs[0], nbrLimbsDivisor * sizeof(limb));
        }
        else {
            std::memcpy(&adjustedArgument[0], &Divisor.limbs[nbrLimbsDivisor - nbrLimbs], nbrLimbs * sizeof(limb));
        }
        MultiplyBigNbrByMinPowerOf2(power2, adjustedArgument, nbrLimbs, adjustedArgument);
        // Initialize approximate inverse.
        inverse = MAX_VALUE_LIMB / ((double)adjustedArgument[nbrLimbs - 1] + 1);
        approxInv[nbrLimbs - 1] = 1;
        if (inverse <= 1) {
            approxInv[nbrLimbs - 2] = 0;
        }
        else {
            approxInv[nbrLimbs - 2] = (int)floor((inverse - 1)*MAX_VALUE_LIMB);
        }
        // Perform Newton approximation loop.
        // Get bit length of each cycle.
        bitLengthNbrCycles = 0;
        bitLength = nbrLimbs * BITS_PER_GROUP;
        while (bitLength >= BITS_PER_GROUP) {
            bitLengthCycle[bitLengthNbrCycles++] = bitLength;
            bitLength = (bitLength + 1) >> 1;
        }
        // Each loop increments precision.
        // Use Newton iteration: x_{n+1} = x_n(2 - x_n)
        while (--bitLengthNbrCycles >= 0) {
            limb *ptrArrAux;
            int limbLength;

            bitLength = bitLengthCycle[bitLengthNbrCycles];
            limbLength = (bitLength + 3 * (BITS_PER_GROUP)-1) / BITS_PER_GROUP;
            if (limbLength > nbrLimbs) {
                limbLength = nbrLimbs;
            }
            // Compute x(2-Nx).
            // Multiply by divisor.
            multiply(&approxInv[nbrLimbs - limbLength],
                &adjustedArgument[nbrLimbs - limbLength], arrAux, limbLength, NULL);
            // Subtract arrAux from 2.
            ptrArrAux = &arrAux[limbLength];
            for (idx = limbLength - 1; idx > 0; idx--) {
                *ptrArrAux = MAX_VALUE_LIMB - *ptrArrAux;
                ptrArrAux++;
            }
            *ptrArrAux = 1 - *ptrArrAux;
            // Multiply arrAux by approxInv.
            multiply(&arrAux[limbLength], &approxInv[nbrLimbs - limbLength], approxInv, limbLength, NULL);
             std::memmove(&approxInv[nbrLimbs - limbLength], &approxInv[limbLength - 1], limbLength * sizeof(limb));
        }
        // Multiply approxInv by argument to obtain the quotient.
        if (nbrLimbsDividend >= nbrLimbs) {
            multiply(&Dividend.limbs[nbrLimbsDividend - nbrLimbs],
                approxInv, approxInv, nbrLimbs, NULL);
        }
        else {
            std::memset(arrAux, 0, (nbrLimbs - nbrLimbsDividend) * sizeof(limb));
            std::memcpy(&arrAux[nbrLimbs - nbrLimbsDividend], Dividend.limbs, nbrLimbsDividend * sizeof(limb));
            multiply(arrAux, approxInv, approxInv, nbrLimbs, NULL);
        }             // approxInv holds the quotient.
                      // Shift left quotient power2 bits into result.
        ptrDest = &approxInv[nbrLimbs - 1];
        oldLimb = 0;
        for (index = nbrLimbs; index >= 0; index--) {
            newLimb = *ptrDest;
            *(ptrDest++) = ((newLimb << power2) |
                (oldLimb >> (BITS_PER_GROUP - power2))) & MAX_VALUE_LIMB;
            oldLimb = newLimb;
        }

        // Determine number of limbs of quotient.
        nbrLimbsQuotient = nbrLimbsDividend - nbrLimbsDivisor;
        ptrDivisor = &Divisor.limbs[nbrLimbsDivisor - 1];
        ptrDividend = &Dividend.limbs[nbrLimbsDividend - 1];
        for (idx = nbrLimbsDivisor - 1; idx > 0; idx--) {
            if (*ptrDividend != *ptrDivisor) {
                break;
            }
            ptrDividend--;
            ptrDivisor--;
        }
        if (*ptrDividend >= *ptrDivisor) {
            nbrLimbsQuotient++;
        }
        ptrQuotient = &approxInv[2 * nbrLimbs - nbrLimbsQuotient];
        if (approxInv[2 * nbrLimbs - 1] == 0)
        {  // Most significant byte is zero, so it is not part of the quotient. 
            ptrQuotient--;
        }
        ptrQuot = ptrQuotient;
        if (*(ptrQuotient - 1) > (7 << (BITS_PER_GROUP - 3))) {
            // Increment quotient.
            for (idx = 0; idx <= nbrLimbsQuotient; idx++) {
                if ((++(*(ptrQuotient + idx))) & MAX_INT_NBR) {
                    break;
                }
                *(ptrQuotient + idx) = 0;
            }
            if (idx >= nbrLimbsQuotient)
            {                // Roll back on overflow.
                for (idx = 0; idx < nbrLimbsQuotient; idx++) {
                    if (--(*(ptrQuotient + idx)) >= 0) {
                        break;
                    }
                    *(ptrQuotient + idx) = MAX_VALUE_LIMB;
                }
            }
            if (approxInv[2 * nbrLimbs - 1] != 0)
            {    // Most significant byte is not zero, so it is part of the quotient.
                ptrQuot = &approxInv[2 * nbrLimbs - nbrLimbsQuotient];
            }
            // Test whether the quotient is correct.
            // It is correct only if multiplied by the divisor, it is <= than the dividend.
            if (nbrLimbsQuotient > nbrLimbsDivisor) {
                std::memcpy(&approxInv[0], Divisor.limbs, nbrLimbsDivisor * sizeof(limb));
                std::memset(&approxInv[nbrLimbsDivisor], 0, (nbrLimbsQuotient - nbrLimbsDivisor) * sizeof(limb));
                multiply(&approxInv[0], ptrQuot, arrAux, nbrLimbsQuotient, NULL);
            }
            else {
                std::memset(&approxInv[2 * nbrLimbs], 0, (nbrLimbsDivisor - nbrLimbsQuotient) * sizeof(limb));
                multiply(Divisor.limbs, ptrQuot, arrAux, nbrLimbsDivisor, NULL);
            }
            ptrDividend = (limb *)&Dividend.limbs[Dividend.nbrLimbs - 1];
            ptrDest = &arrAux[Dividend.nbrLimbs - 1];
            for (idx = Dividend.nbrLimbs - 1; idx > 0; idx--) {
                if (*ptrDividend != *ptrDest) {
                    break;
                }
                ptrDividend--;
                ptrDest--;
            }
            if (*ptrDividend < *ptrDest) {  // Decrement quotient.
                ptrQuotient = ptrQuot;
                for (idx = 0; idx < nbrLimbsQuotient; idx++) {
                    if (--(*ptrQuotient) >= 0) {
                        break;
                    }
                    *(ptrQuotient++) = MAX_VALUE_LIMB;
                }
                if (idx == nbrLimbsQuotient) {
                    nbrLimbsQuotient--;
                }
            }
        }
        std::memcpy(&Quotient.limbs[0], ptrQuot, nbrLimbsQuotient * sizeof(limb));
        Quotient.nbrLimbs = nbrLimbsQuotient;
    }

    /* if dividend and divisor have the same sign, or if the quotient is 0, 
    quotient sign is +ve, otherwise it's -ve. */
    if (Dividend.sign == Divisor.sign || (Quotient.limbs[0] == 0 && Quotient.nbrLimbs == 1)) {
        Quotient.sign = SIGN_POSITIVE;
    }
    else {
        Quotient.sign = SIGN_NEGATIVE;
    }

    while (Quotient.nbrLimbs > 1) {
        if (Quotient.limbs[Quotient.nbrLimbs - 1] == 0)
            Quotient.nbrLimbs--;
        else break;
    }
    return Quotient;
}

/* calculate quotient and remainder of division by integer. The remainder will
have the same sign as the dividend and the quotient is rounded towards zero. */
void BigIntDivRem(const BigInteger& Dividend, const int Divisor, 
    BigInteger &Quotient, int &rem) {
    int len = Dividend.nbrLimbs;
    DivBigNbrByInt(Dividend.limbs, std::abs(Divisor), Quotient.limbs, len, &rem);
    /* now set the correct signs for quotient and remainder */
    if (Divisor >= 0)
        if (Dividend.sign == SIGN_POSITIVE) {
            Quotient.sign = SIGN_POSITIVE;
        }
        else {  /* dividend is -ve, divisor is +ve */
            Quotient.sign = SIGN_NEGATIVE;
            if (rem > 0)
                rem = -rem;
        }
    else   // Divisor is -ve
        if (Dividend.sign == SIGN_POSITIVE)
            Quotient.sign = SIGN_NEGATIVE;  /* dividend is +ve, divisor is -ve */
        else {
            Quotient.sign = SIGN_POSITIVE;  /* dividend and divisor both -ve */
            if (rem > 0)
                rem = -rem;
        }

    while (len > 1 && Quotient.limbs[len - 1] == 0)
        len--;  // remove any leading zeros
    Quotient.nbrLimbs = len;

#ifdef _DEBUG    /* (Debug only) recheck quotient and remainder using Znums */
    {  Znum divz, Qz;
    BigtoZ(divz, Dividend);
    BigtoZ(Qz, Quotient);
    if (Qz != divz / Divisor)
        std::cout << "BigIntDivideInt error: Quotient expected " << divz / Divisor
        << " got " << Qz << '\n'
        << " Dividend = " << Dividend << " divisor = " << Divisor << '\n';
    if (rem != (divz % Divisor))
        std::cout << "BigIntDivideInt error: remainder expected " << divz % Divisor
        << " got " << rem << '\n';
    }
#endif
    return;
}

/* Return Dividend/Divisor. Used for operator overloading. The remainder will
have the same sign as the dividend and the quotient is rounded towards zero */
BigInteger BigIntDivideInt(const BigInteger &Dividend, const int Divisor) {
    BigInteger Quotient;
    int rem;
    BigIntDivRem(Dividend, Divisor, Quotient, rem);
    return Quotient;
}

/* calculate BigInt modulo divisor (used for operator overloading). The remainder
will have the same sign as the dividend and the quotient is rounded towards zero
(truncation division). This function is the same as BigIntDivideInt, but returns 
the remainder rather than the quotient. */
int getRemainder(const BigInteger& Dividend, int Divisor) {
    BigInteger Quotient;
    int rem;
    BigIntDivRem(Dividend, Divisor, Quotient, rem);
    return rem;
}


bool BigIntIsZero(const BigInteger* value)
{
    if ((value->nbrLimbs == 1) && (value->limbs[0] == 0))
    {
        return true;     // Number is zero.
    }
    if (value->nbrLimbs == 0) return true;  /* no limbs (not initialised) */

    return false;      // Number is not zero.
}

/* change -ve to +ve and vice versa. Does the same job as BigIntNegate  */
void BigIntChSign(BigInteger* value)
{
    if ((value->nbrLimbs == 1) && (value->limbs[0] == 0))
    {    // Value is zero. Do not change sign.
        return;
    }
    if (value->sign == SIGN_POSITIVE)
    {
        value->sign = SIGN_NEGATIVE;
    }
    else
    {
        value->sign = SIGN_POSITIVE;
    }
}

/* change -ve to +ve */
void BigIntAbs(BigInteger* value) {
    if (value->sign == SIGN_NEGATIVE)
        value->sign = SIGN_POSITIVE;
}

void IntArray2BigInteger(const int* ptrValues, BigInteger* bigint)
{
    const int* piValues = ptrValues;
    limb* destLimb = bigint->limbs;
    int nbrLimbs = *piValues;
    piValues++;
    if (nbrLimbs > 0)
    {
        bigint->sign = SIGN_POSITIVE;
    }
    else
    {
        bigint->sign = SIGN_NEGATIVE;
        nbrLimbs = -nbrLimbs;
    }
    if (NumberLength == 1)
    {
        *destLimb = *piValues;
        bigint->nbrLimbs = 1;
    }
    else
    {
        int ctr;
        bigint->nbrLimbs = nbrLimbs;
        for (ctr = 0; ctr < nbrLimbs; ctr++)
        {
            *destLimb = *piValues;
            destLimb++;
            piValues++;
        }
        for (; ctr < NumberLength; ctr++)
        {
            *destLimb = 0;
            destLimb++;
        }
    }
}

void BigInteger2IntArray(/*@out@*/int* ptrValues, const BigInteger* bigint)
{
    int* pValues = ptrValues;
    const limb* srcLimb = bigint->limbs;
    if (NumberLength == 1)
    {
        *pValues = ((bigint->sign == SIGN_POSITIVE) ? 1 : -1);
        *(pValues + 1) = *srcLimb;
    }
    else
    {
        int nbrLimbs;
        nbrLimbs = getNbrLimbs(srcLimb);
        *pValues = ((bigint->sign == SIGN_POSITIVE) ? nbrLimbs : -nbrLimbs);
        pValues++;
        for (int ctr = 0; ctr < nbrLimbs; ctr++)
        {
            *pValues = *srcLimb;
            pValues++;
            srcLimb++;
        }
    }
}

/* set Result to 2^exponent  */
void BigIntPowerOf2(BigInteger* pResult, int exponent)
{
    unsigned int power2 = (unsigned int)exponent % (unsigned int)BITS_PER_GROUP;
    int nbrLimbs = exponent / BITS_PER_GROUP;
    if (nbrLimbs > 0)
    {
        int nbrLimbsBytes = nbrLimbs * (int)sizeof(limb);
        std::memset(pResult->limbs, 0, nbrLimbsBytes);
    }
    pResult->limbs[nbrLimbs] = (int)(1U << power2);
    pResult->nbrLimbs = nbrLimbs + 1;
    pResult->sign = SIGN_POSITIVE;
}

/* Sum = Nbr1 + Nbr2 */
void AddBigNbr(const limb* pNbr1, const limb* pNbr2, limb* pSum, int nbrLen)
{
    unsigned int carry = 0U;
    const limb* ptrNbr1 = pNbr1;
    const limb* ptrNbr2 = pNbr2;
    const limb* ptrEndSum = pSum + nbrLen;
    for (limb* ptrSum = pSum; ptrSum < ptrEndSum; ptrSum++)
    {
        unsigned int tmp;
        carry = (carry >> BITS_PER_GROUP) + (unsigned int)*ptrNbr1 +
            (unsigned int)*ptrNbr2;
        tmp = carry & MAX_INT_NBR_U;
        *ptrSum = (int)tmp;
        ptrNbr1++;
        ptrNbr2++;
    }
}

/* Diff = Nbr1 - Nbr2 */
void SubtractBigNbr(const limb* pNbr1, const limb* pNbr2, limb* pDiff, int nbrLen)
{
    unsigned int borrow = 0U;
    const limb* ptrNbr1 = pNbr1;
    const limb* ptrNbr2 = pNbr2;
    const limb* ptrEndDiff = pDiff + nbrLen;
    for (limb* ptrDiff = pDiff; ptrDiff < ptrEndDiff; ptrDiff++)
    {
        unsigned int tmp;
        borrow = (unsigned int)*ptrNbr1 - (unsigned int)*ptrNbr2 -
            (borrow >> BITS_PER_GROUP);
        tmp = borrow & MAX_INT_NBR_U;
        *ptrDiff = (int)tmp;
        ptrNbr1++;
        ptrNbr2++;
    }
}

/* nbr = nbr/2 */
void BigIntDivideBy2(BigInteger* nbr)
{
    int nbrLimbs = nbr->nbrLimbs;
    limb* ptrDest = &nbr->limbs[0];
    unsigned int curLimb = (unsigned int)*ptrDest;
    for (int ctr = 1; ctr < nbrLimbs; ctr++)
    {  // Process starting from least significant limb.
        unsigned int nextLimb = (unsigned int)*(ptrDest + 1);
        *ptrDest = (int)(((curLimb >> 1) | (nextLimb << BITS_PER_GROUP_MINUS_1)) &
            MAX_VALUE_LIMB);
        ptrDest++;
        curLimb = nextLimb;
    }
    *ptrDest = (int)((curLimb >> 1) & MAX_VALUE_LIMB);
    if ((nbrLimbs > 1) && (nbr->limbs[nbrLimbs - 1] == 0))
    {
        nbr->nbrLimbs--;
    }
}

int BigIntJacobiSymbol(const BigInteger* upper, const BigInteger* lower)
{
    int t;
    int power2;
    static BigInteger a;
    static BigInteger m;
    static BigInteger tmp;
    m = *lower; // CopyBigInt(&m, lower);               // m <- lower
    DivideBigNbrByMaxPowerOf2(&power2, m.limbs, &m.nbrLimbs);
    a = *upper % *lower; // (void)BigIntRemainder(upper, lower, &a);   // a <- upper % lower
    t = 1;
    if (upper->sign == SIGN_NEGATIVE)
    {
        a.sign = SIGN_POSITIVE;
        if ((m.limbs[0] & 3) == 3)
        {
            t = -1;
        }
    }
    while (!BigIntIsZero(&a))             // a != 0
    {
        while ((a.limbs[0] & 1) == 0)
        {     // a is even.
            BigIntDivideBy2(&a);              // a <- a / 2
            if (((m.limbs[0] & 7) == 3) || ((m.limbs[0] & 7) == 5))
            {   // m = 3 or m = 5 (mod 8)
                t = -t;
            }
        }
        tmp = a; // CopyBigInt(&tmp, &a);      // Exchange a and m.
        a = m;   // CopyBigInt(&a, &m);
        m = tmp; // CopyBigInt(&m, &tmp);
        if ((a.limbs[0] & m.limbs[0] & 3) == 3)
        {   // a = 3 and m = 3 (mod 4)
            t = -t;
        }
        tmp = a % m;  // (void)BigIntRemainder(&a, &m, &tmp);
        a = tmp;      // CopyBigInt(&a, &tmp);      // a <- a % m  
    }
    if ((m.nbrLimbs == 1) && (m.limbs[0] == 1))
    {              // Absolute value of m is 1.
        return t;
    }
    return 0;
}

enum eOper {
    OPERATION_AND = 0,
    OPERATION_OR,
    OPERATION_XOR,
};

static void ConvertToTwosComplement(BigInteger* value)
{
    int idx;
    int nbrLimbs;
    limb* ptrLimb;
    if (value->sign == SIGN_POSITIVE)
    {    // If number is positive, no conversion is needed.
        while (value->nbrLimbs > 1)
        {
            if (value->limbs[value->nbrLimbs - 1] != 0)
            {
                break;
            }
            value->nbrLimbs--;
        }
        return;
    }
    nbrLimbs = value->nbrLimbs;
    ptrLimb = &value->limbs[0];
    for (idx = 0; idx < nbrLimbs; idx++)
    {
        if (*ptrLimb != 0)
        {
            break;
        }
        ptrLimb++;
    }
    if (idx < nbrLimbs)
    {
        *ptrLimb = (int)(LIMB_RANGE - (unsigned int)*ptrLimb);
        ptrLimb++;
    }
    for (; idx < nbrLimbs; idx++)
    {
        *ptrLimb = MAX_INT_NBR - *ptrLimb;
        ptrLimb++;
    }
}

static void InternalBigIntLogical(const BigInteger* firstArgum,
    const BigInteger* secondArgum, BigInteger* result, enum eOper operation)
{
    const BigInteger* firstArg;
    const BigInteger* secondArg;
    int idx;
    int carryFirst = 0;
    int carrySecond = 0;
    int limbFirst;
    int limbSecond;
    if (firstArgum->nbrLimbs < secondArgum->nbrLimbs)
    {    // After the exchange, firstArg has not fewer limbs than secondArg.
        firstArg = secondArgum;
        secondArg = firstArgum;
    }
    else
    {
        firstArg = firstArgum;
        secondArg = secondArgum;
    }
    for (idx = 0; idx < secondArg->nbrLimbs; idx++)
    {
        limbFirst = firstArg->limbs[idx];
        limbSecond = secondArg->limbs[idx];
        if (firstArg->sign == SIGN_NEGATIVE)
        {
            carryFirst -= limbFirst;
            limbFirst = carryFirst & MAX_INT_NBR;
            carryFirst >>= 31;
        }
        if (secondArg->sign == SIGN_NEGATIVE)
        {
            carrySecond -= limbSecond;
            limbSecond = carrySecond & MAX_INT_NBR;
            carrySecond >>= 31;
        }
        if (operation == OPERATION_AND)
        {
            result->limbs[idx] = limbFirst & limbSecond;
        }
        else if (operation == OPERATION_OR)
        {
            result->limbs[idx] = limbFirst | limbSecond;
        }
        else
        {
            result->limbs[idx] = limbFirst ^ limbSecond;
        }
    }
    if (firstArg->sign == SIGN_POSITIVE)
    {
        limbFirst = 0;
    }
    else
    {
        limbFirst = -1;
    }
    for (; idx < firstArg->nbrLimbs; idx++)
    {
        limbSecond = secondArg->limbs[idx];
        if (secondArg->sign == SIGN_NEGATIVE)
        {
            carrySecond -= limbSecond;
            limbSecond = carrySecond & MAX_INT_NBR;
            carrySecond >>= 31;
        }
        if (operation == OPERATION_AND)
        {
            result->limbs[idx] = limbFirst & limbSecond;
        }
        else if (operation == OPERATION_OR)
        {
            result->limbs[idx] = limbFirst | limbSecond;
        }
        else
        {
            result->limbs[idx] = limbFirst ^ limbSecond;
        }
    }
    // Generate sign of result according to operation and
    // signs of arguments.
    if (operation == OPERATION_AND)
    {
        if ((firstArg->sign == SIGN_NEGATIVE) && (secondArg->sign == SIGN_NEGATIVE))
        {
            result->sign = SIGN_NEGATIVE;
        }
        else
        {
            result->sign = SIGN_POSITIVE;
        }
    }
    else if (operation == OPERATION_OR)
    {
        if ((firstArg->sign == SIGN_POSITIVE) && (secondArg->sign == SIGN_POSITIVE))
        {
            result->sign = SIGN_POSITIVE;
        }
        else
        {
            result->sign = SIGN_NEGATIVE;
        }
    }
    else     // XOR operation
    {
        if (firstArg->sign == secondArg->sign)
        {
            result->sign = SIGN_POSITIVE;
        }
        else
        {
            result->sign = SIGN_NEGATIVE;
        }
    }
    result->nbrLimbs = firstArg->nbrLimbs;
    ConvertToTwosComplement(result);
}

/* result = first AND second. Used for operator overloading */
void BigIntAnd(const BigInteger* firstArg, const BigInteger* secondArg, 
    BigInteger* result) {
    InternalBigIntLogical(firstArg, secondArg, result, OPERATION_AND);
}

/* result = first OR second. Used for operator overloading */
void BigIntOr(const BigInteger* firstArg, const BigInteger* secondArg, 
    BigInteger* result) {
    InternalBigIntLogical(firstArg, secondArg, result, OPERATION_OR);
}

/* result = first XOR second. Used for operator overloading */
void BigIntXor(const BigInteger* firstArg, 	const BigInteger* secondArg, 
    BigInteger* result) {
    InternalBigIntLogical(firstArg, secondArg, result, OPERATION_XOR);
}

void CompressLimbsBigInteger(/*@out@*/limb* ptrValues, const BigInteger* bigint)
{
    if (NumberLength == 1)
    {
        *ptrValues = bigint->limbs[0];
    }
    else
    {
        int numberLengthBytes = NumberLength * (int)sizeof(limb);
        int nbrLimbs = bigint->nbrLimbs;
        if (nbrLimbs > NumberLength)
        {
            std::memcpy(ptrValues, bigint->limbs, numberLengthBytes);
        }
        else
        {
            int nbrLimbsBytes = nbrLimbs * (int)sizeof(limb);
            std::memcpy(ptrValues, bigint->limbs, nbrLimbsBytes);
            nbrLimbsBytes = numberLengthBytes - nbrLimbsBytes;
            std::memset(ptrValues + nbrLimbs, 0, nbrLimbsBytes);
        }
    }
}

void UncompressLimbsBigInteger(const limb* ptrValues, /*@out@*/BigInteger* bigint)
{
    if (NumberLength == 1)
    {
        bigint->limbs[0] = *ptrValues;
        bigint->nbrLimbs = 1;
    }
    else
    {
        int nbrLimbs;
        const limb* ptrValue1;
        int numberLengthBytes = NumberLength * (int)sizeof(limb);
        std::memcpy(bigint->limbs, ptrValues, numberLengthBytes);
        ptrValue1 = ptrValues + NumberLength;
        for (nbrLimbs = NumberLength; nbrLimbs > 1; nbrLimbs--)
        {
            ptrValue1--;
            if (*ptrValue1 != 0)
            {
                break;
            }
        }
        bigint->nbrLimbs = nbrLimbs;
    }
}
