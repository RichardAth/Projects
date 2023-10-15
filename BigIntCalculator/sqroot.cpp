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


static limb approxInvSqrt[MAX_LEN];
limb approxInv[MAX_LEN];
limb adjustedArgument[MAX_LEN];
limb arrAux[MAX_LEN];
int bitLengthCycle[20];

void squareRoot(const Znum &arg, Znum & sqRoot) {
    assert(arg >= 0);  // no imaginary numbers allowed here!!
    mpz_sqrt(ZT(sqRoot), ZT(arg));
}

/* return true iff arg is a perfect square */
bool isPerfectSquare(const Znum &arg, Znum &sqRoot) {
    Znum rem;
    assert(arg >= 0);  // no imaginary numbers allowed here!!
    mpz_sqrtrem(ZT(sqRoot), ZT(rem), ZT(arg));
    return rem == 0;  // true if arg is a perfect square
}

/* faster version for small n. Do some checks before sqrt
because sqrt function is relatively slow. */
bool isPerfectSquare(__int64 x) {
    if (x < 0) return false;
    if (x == 0) return true;
    //while ((x & 0x3) == 0) 
    //	x >>= 2;   // divide by largest possible power of 4

    /* can use _BitScanForward64 instead (slightly faster) */
    unsigned long ix;
    auto result = _BitScanForward64(&ix, x);  // for gcc compiler use __builtin_ctzll instead
    if ((ix & 1) == 1)
        return false;     // if x contains an odd number of 2 factors it is 
                          // not a perfect square
    x >>= ix;

    /* at this stage, x must be odd for a perfect square,
    so sqrt(x) is odd if x is a perfect square, in which case
    (sqrt(x))(mod 8) is 1 3 5 or 7. So x (mod 8) has to be 1 */
    if ((x & 0x7) != 1)
        return false;  // not a perfect square

    long long s = std::llround(std::sqrt(x));  // slightly faster than intSqrt
    //long long s = intSqrt(x);
    return (s * s == x);
}

/* if x >= 0, returns floor(sqrt(x)) in ss */
bool isPerfectSquare(__int64 x, __int64 &ss) {
    long long s;
    if (x < 0) return false;
    if (x == 0) {
        ss = 0;
        return true;
    }
    //while ((x & 0x3) == 0) 
    //	x >>= 2;   // divide by largest possible power of 4

    /* can use _BitScanForward64 instead (slightly faster) */
    unsigned long ix;
    auto result = _BitScanForward64(&ix, x);  // for gcc compiler use __builtin_ctzll instead
    if ((ix & 1) == 1)
        return false;     // if x contains an odd number of 2 factors it is 
                          // not a perfect square
    x >>= ix;

    /* at this stage, x must be odd for a perfect square,
    so sqrt(x) is odd if x is a perfect square, in which case
    (sqrt(x))(mod 8) is 1 3 5 or 7. So x (mod 8) has to be 1 */
    if ((x & 0x7) != 1)
        return false;  // not a perfect square

    s = std::llround(std:: sqrt(x));  // slightly faster than intSqrt
    ss = s << (ix / 2);
    return (s * s == x);
}

// This routine uses Newton iteration: if x is an approximate inverse square root of N,
// a better approximation is: x(3-Nxx)/2. After the inverse square root is computed,
// the square root is found just by multiplying by N.
// The argument is multiplied by a power of 4 so the most significant limb is
// between LIMB_RANGE/4 and LIMB_RANGE-1 and there is an even number of limbs.
// At the end of the calculation, the result is divided by the power of 2.

// All computations are done in little-endian notation.
// Find power of 4 that divides the number.
// output: pNbrLimbs = pointer to number of limbs
//         pPower4 = pointer to power of 4.
static void MultiplyBigNbrByMinPowerOf4(int* pPower4, const limb* number,
    int len, limb* dest)
{
    int power4;
    unsigned int power2 = 0;
    limb mostSignficLimb;
    unsigned int shLeft;
    unsigned int shRight;
    limb prevLimb;
    limb currLimb;
    limb* ptrDest;
    const limb* ptrSrc;

    ptrSrc = number;
    if ((len & 1) != 0)
    {
        power2 = (unsigned int)BITS_PER_GROUP;
    }
    else
    {
        power2 = 0U;
    }
    shLeft = 0U;
    mostSignficLimb = *(number + len - 1);
    for (unsigned int mask = LIMB_RANGE / 2U; mask > 0U; mask >>= 1)
    {
        if (((unsigned int)mostSignficLimb & mask) != 0U)
        {
            break;
        }
        shLeft++;
    }
    power2 = (power2 + shLeft) & 0xFFFFFFFE;
    power4 = (int)power2 / 2;
    ptrDest = dest;
    if (power2 > (unsigned int)BITS_PER_GROUP)
    {
        shLeft = power2 - (unsigned int)BITS_PER_GROUP;
        *ptrDest = 0;
        ptrDest++;
    }
    else
    {
        shLeft = power2;
    }
    // Multiply number by this power.
    prevLimb = 0;
    shRight = (unsigned int)BITS_PER_GROUP - shLeft;
    for (int index2 = len; index2 >= 0; index2--)
    {
        currLimb = *ptrSrc;
        ptrSrc++;
        *ptrDest = (int)((((unsigned int)currLimb << shLeft) |
            ((unsigned int)prevLimb >> shRight)) & MAX_VALUE_LIMB);
        ptrDest++;
        prevLimb = currLimb;
    }
    *pPower4 = power4;
}

/* used for sqRoot method in BigInteger class */
void squareRoot(const limb* argument, limb* sqRoot, int len, int* pLenSqRoot)
{
    int index;
    int length = len;
    int lenInvSqrt;
    int lenInvSqrt2;
    int bitLength;
    int bitLengthNbrCycles;
    int idx;
    limb prev;
    const limb* ptrApproxInvSqrt;
    unsigned int shRight;
    unsigned int shLeft;
    limb* ptrDest;
    const limb* ptrSrc;
    double invSqrt;
    unsigned int currLimb;
    int lenBytes;
    int shRightBits;

    if ((len == 1) && (*argument == 0))
    {    // Input number is zero. Its square root is also zero.
        *pLenSqRoot = 1;
        *sqRoot = 0;
        return;
    }
    // If the number of limbs is even and the upper half of the number has all
    // limbs set to MAX_VALUE_LIMB, there could be overflow. Set the square root directly.
    if ((length % 2) == 0)
    {      // Even number of limbs.
        for (index = length / 2; index < length; index++)
        {
            if (*(argument + index) != (int)MAX_VALUE_LIMB)
            {
                break;
            }
        }
        if (index == length)
        {   // Set square root and go out.
            for (index = (length / 2) - 1; index >= 0; index--)
            {
                *(sqRoot + index) = (int)MAX_VALUE_LIMB;
            }
            *pLenSqRoot = length / 2;
            return;
        }
    }
    // Obtain logarithm of 2 of argument.
    for (index = length - 1; index > 2; index--)
    {
        if (*(argument + index) != 0)
        {
            break;
        }
    }
    lenBytes = (length + 1) / 2 * (int)sizeof(limb);
    (void)std::memset(sqRoot, 0, lenBytes);
    if (index <= 1)
    {                // Argument is small (1 limb), so compute directly its square root.
        if (index == 0)
        {
            *sqRoot = (int)std::floor(std::sqrt(*argument) + 0.000001);
        }
        else
        {   /* argument has two limbs. */
            limb square[3];   // MultBigNbr routine uses an extra limb for result.
            long long llsquare;
            double dArg = *argument + (double)*(argument + 1) * (double)LIMB_RANGE;
            dArg = std::floor(std::sqrt(dArg + 0.5)); /* get approximate square root */
            if (dArg == (double)LIMB_RANGE)
            {
                dArg = (double)MAX_VALUE_LIMB;
            }
            *sqRoot = (int)dArg;
            llsquare = (long long)*sqRoot * (*sqRoot);
            square[0] = (int)llsquare & MAX_VALUE_LIMB;
            llsquare >>= BITS_PER_GROUP;
            square[1] = llsquare & MAX_VALUE_LIMB;
            assert(llsquare <= MAX_VALUE_LIMB);
            //MultBigNbr(sqRoot, sqRoot, square, 1);
            if ((square[1] > *(argument + 1)) ||
                ((square[1] == *(argument + 1)) && (square[0] > *argument)))
            {   /*if sqRoot^2 > argument, decrease it */
                *sqRoot--;
            }
        }
        *pLenSqRoot = 1;    // Only one limb for result.
        return;
    }
    // Adjust argument.
    MultiplyBigNbrByMinPowerOf4(&shRightBits, argument, length, adjustedArgument);
    if ((length % 2) != 0)
    {
        length++;   // Make number of limbs even.
    }
    lenInvSqrt = (length + 5) / 2;
    lenBytes = lenInvSqrt * (int)sizeof(limb);
    (void)std::memset(approxInvSqrt, 0, lenBytes);
    // Initialize approximate inverse square root.
    invSqrt = (double)LIMB_RANGE /
        std::sqrt(getMantissa(adjustedArgument + length, length) * (double)LIMB_RANGE + 1.0);
    approxInvSqrt[lenInvSqrt - 1] = 1;
    invSqrt = (invSqrt - 1.0) * (double)LIMB_RANGE;
    if (invSqrt > (double)MAX_INT_NBR)
    {
        approxInvSqrt[lenInvSqrt - 2] = MAX_INT_NBR;
    }
    else
    {
        approxInvSqrt[lenInvSqrt - 2] = (int)invSqrt;
    }
    // Perform Newton approximation loop.
    // Get bit length of each cycle.
    bitLengthNbrCycles = 0;
    bitLength = lenInvSqrt * BITS_PER_GROUP;
    while (bitLength > 7)
    {
        bitLengthCycle[bitLengthNbrCycles] = bitLength;
        bitLengthNbrCycles++;
        bitLength = (bitLength + 1) / 2;
    }
    // Each loop increments precision.
    bitLengthNbrCycles--;
    while (bitLengthNbrCycles >= 0)
    {
        int limbLength;
        unsigned int prevLimb;
        limb* ptrArrAux;
        bitLength = bitLengthCycle[bitLengthNbrCycles];
        limbLength = (bitLength + (3 * BITS_PER_GROUP) - 1) / BITS_PER_GROUP;
        if (limbLength > lenInvSqrt)
        {
            limbLength = lenInvSqrt;
        }
        // Point to least significant limb for this cycle.
        ptrApproxInvSqrt = &approxInvSqrt[lenInvSqrt - limbLength];
        // Compute x(3-Nxx)/2. Start by squaring.
        multiply(ptrApproxInvSqrt, ptrApproxInvSqrt, approxInv, limbLength, NULL);
        // Multiply by argument.
        multiply(&approxInv[limbLength - 1], &adjustedArgument[length - limbLength], arrAux, limbLength, NULL);
        // Subtract arrAux from 3.
        ptrArrAux = &arrAux[limbLength];
        for (idx = limbLength - 1; idx > 0; idx--)
        {
            *ptrArrAux ^= (int)MAX_VALUE_LIMB;
            ptrArrAux++;
        }
        *ptrArrAux = 2 - *ptrArrAux;
        // Divide arrAux by 2.
        prevLimb = 0U;
        for (idx = limbLength; idx > 0; idx--)
        {
            currLimb = (unsigned int)*ptrArrAux;
            *ptrArrAux = (int)(((currLimb >> 1) | (prevLimb << BITS_PER_GROUP_MINUS_1)) &
                MAX_VALUE_LIMB);
            ptrArrAux--;
            prevLimb = currLimb;
        }
        // Multiply arrAux by approxInvSqrt.
        multiply(ptrArrAux + 1, &approxInvSqrt[lenInvSqrt - limbLength], approxInv, limbLength, NULL);
        lenBytes = limbLength * (int)sizeof(limb);
        std::memcpy(&approxInvSqrt[lenInvSqrt - limbLength], &approxInv[limbLength - 1], lenBytes);
        bitLengthNbrCycles--;
    }
    // Multiply approxInvSqrt by argument to obtain the square root.
    if (length > lenInvSqrt)
    {
        multiply(approxInvSqrt, &adjustedArgument[length - lenInvSqrt], approxInv, lenInvSqrt, NULL);
    }
    else
    {
        multiply(approxInvSqrt, adjustedArgument, approxInv, length, NULL);
    }             // approxInv holds the square root.
    lenInvSqrt2 = (length + 1) / 2;
    shLeft = (unsigned int)BITS_PER_GROUP - 3U;
    if ((unsigned int)approxInv[(2 * lenInvSqrt) - lenInvSqrt2 - 2] > (7U << shLeft))
    {                   // Increment square root.
        for (idx = (2 * lenInvSqrt) - lenInvSqrt2 - 1; idx < (2 * lenInvSqrt) - 1; idx++)
        {
            approxInv[idx] = (approxInv[idx] + 1) & MAX_INT_NBR;
            if (approxInv[idx] != 0)
            {
                break;
            }
        }
        if (idx == ((2 * lenInvSqrt) - 1))
        {                // Roll back on overflow.
            for (idx = (2 * lenInvSqrt) - lenInvSqrt2 - 1; idx < ((2 * lenInvSqrt) - 1); idx++)
            {
                approxInv[idx]--;
                if (approxInv[idx] >= 0)
                {
                    break;
                }
                approxInv[idx] = MAX_VALUE_LIMB;
            }
        }
        // Test whether square root is correct.
        // It is correct only if when squared, it is <= than the argument.
        ptrApproxInvSqrt = &approxInv[(2 * lenInvSqrt) - lenInvSqrt2 - 1];
        multiply(ptrApproxInvSqrt, ptrApproxInvSqrt, approxInvSqrt, lenInvSqrt2, NULL);
        for (idx = lenInvSqrt - 1; idx > 0; idx--)
        {
            if (adjustedArgument[idx] != approxInvSqrt[idx])
            {
                break;
            }
        }
        if (adjustedArgument[idx] < approxInvSqrt[idx])
        {                // Incorrect square root: roll back.
            for (idx = (2 * lenInvSqrt) - lenInvSqrt2 - 1; idx < (2 * lenInvSqrt) - 1; idx++)
            {
                approxInv[idx]--;
                if (approxInv[idx] >= 0)
                {
                    break;
                }
                approxInv[idx] = MAX_VALUE_LIMB;
            }
        }
    }
    length = (length + 1) / 2;
    *pLenSqRoot = length;
    // Shift right shRight bits into result.
    shRight = (unsigned int)shRightBits;
    ptrSrc = &approxInv[(2 * lenInvSqrt) - *pLenSqRoot - 1 + length - 1];
    ptrDest = sqRoot + length - 1;
    prev = 0;
    shLeft = (unsigned int)BITS_PER_GROUP - shRight;
    for (index = length; index > 0; index--)
    {
        *ptrDest = (int)((((unsigned int)prev << shLeft) |
            ((unsigned int)*ptrSrc >> shRight)) & MAX_VALUE_LIMB);
        ptrDest--;
        prev = *ptrSrc;
        ptrSrc--;
    }
}
