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

/* Perform Karatsuba multiplication in little endian order */
#include "pch.h"


#define KARATSUBA_CUTOFF 16   /* n.b. cannot increase this unless modifying 
                        ClassicalMult as well */

static limb arr[4 * MAX_LEN];      /*    3 * changed to 4 * on 16/5/2019, because
                            karatsuba exceeded arrray bound for large numbers */
static limb arrayAux[3 * MAX_LEN];
static int karatLength;
static void Karatsuba(int idxFactor1, const int length, int diffIndex);

#define PROLOG_MULTIPLICATION_DOUBLE                                    \
  factor2_i = arr[idxFactor2 + i];                                    \
  factor2_iPlus1 = arr[idxFactor2 + i + 1];                           \
  Pr = prod_iPlus0 + (uint64_t)factor2_i * factor1_0;                   \
  arrayAux[i] = (int32_t)Pr & MAX_INT_NBR;                            \
  Pr = prod_iPlus1 + (uint64_t)factor2_i * factor1_1 +                  \
       (uint64_t)factor2_iPlus1 * factor1_0 + (Pr >> BITS_PER_GROUP);   \
  arrayAux[i + 1] = (int32_t)Pr & MAX_INT_NBR;
#define MULT_MACRO_DOUBLE(m, n, p)                                      \
  Pr = prod_iPlus##p + (uint64_t)factor2_i * factor1_##p +              \
    (uint64_t)factor2_iPlus1 * factor1_##n + (Pr >> BITS_PER_GROUP);    \
  prod_iPlus##m = (int32_t)Pr & MAX_INT_NBR;
#define EPILOG_MULTIPLICATION_DOUBLE(m, n)                              \
  Pr = (uint64_t)factor2_iPlus1 * factor1_##n + (Pr >> BITS_PER_GROUP); \
  prod_iPlus##m = (uint32_t)Pr & MAX_INT_NBR;                           \
  prod_iPlus##n = (uint32_t)(Pr >> BITS_PER_GROUP);

#define PROLOG_MULTIPLICATION_SINGLE(m)                                 \
  factor2_i = arr[idxFactor2 + m];                                    \
  Pr = prod_iPlus0 + (uint64_t)factor2_i * factor1_0;                   \
  arrayAux[m] = (int32_t)Pr & MAX_INT_NBR;
#define MULT_MACRO_SINGLE(m, n)                                         \
  Pr = prod_iPlus##m + (uint64_t)factor2_i * factor1_##m +              \
       (Pr >> BITS_PER_GROUP);                                          \
  arrayAux[n] = (int32_t)Pr & MAX_INT_NBR;
#define EPILOG_MULTIPLICATION_SINGLE(m)                                 \
  arrayAux[m] = (int32_t)(Pr >> BITS_PER_GROUP);

#define M(n)                                                            \
  uint32_t factor1_##n = arr[idxFactor1 + n];                         \
  uint32_t prod_iPlus##n = 0;


/* this the entry point for external calls 
result = factor1 * factor2.  */
void multiply(const limb factor1[], const limb factor2[], limb result[], const int len, int* pResultLen)
{
    int length = len;
    // Compute length of numbers for each recursion. length is derived from len
    // such that length >= len and length = 2^x * y where y is <= KARATSUBA_CUTOFF
    if (length > KARATSUBA_CUTOFF)
    {
        int div = 1;
        while (length > KARATSUBA_CUTOFF)
        {
            div *= 2;
            length = (length + 1) / 2;
        }
        length *= div;
    }
    if (length > MAX_LEN) {
        StackTrace2();
        std::string line = std::to_string(__LINE__);
        std::string mesg = "number too big : cannot perform multiplication: ";
        mesg += __func__;
        mesg += " line ";  mesg += line;
        mesg += " in file "; mesg += __FILE__;
        throw std::range_error(mesg);
    }

    karatLength = length;
    std::memset(arr, 0, 2 * length * sizeof(limb));      /* copy both operands to arr */
    std::memcpy(&arr[0], factor1, len * sizeof(limb));
    std::memcpy(&arr[length], factor2, len * sizeof(limb));
    Karatsuba(0, length, 2 * length);           /* do the actual multiplication */
    std::memcpy(result, &arr[2 * (karatLength - length)], 2 * length * sizeof(limb));
    if (pResultLen != NULL) {
        if (karatLength > length && arr[2 * (karatLength - length) - 1] == 0) {
            *pResultLen = length * 2 - 1;
        }
        else {
            *pResultLen = length * 2;
        }
    }
}

// The return value is the sign: true: negative.
// In result the absolute value of the difference is computed.
static int absSubtract(int idxMinuend, int idxSubtrahend,
    int idxResult, int nbrLen)
{
    int sign = 0;
    limb carry;
    int i;
    for (i = nbrLen - 1; i >= 0; i--)
    {
        if (arr[idxMinuend + i] != arr[idxSubtrahend + i])
        {
            break;
        }
    }
    if (i >= 0 && arr[idxMinuend + i] < arr[idxSubtrahend + i])
    {
        sign = 1;
        i = idxMinuend;    // Exchange minuend and subtrahend.
        idxMinuend = idxSubtrahend;
        idxSubtrahend = i;
    }
    carry = 0;
    for (i = 0; i < nbrLen; i++)
    {
        carry += arr[idxMinuend + i] - arr[idxSubtrahend + i];
        arr[idxResult + i] = carry & MAX_VALUE_LIMB;
        carry >>= BITS_PER_GROUP;
    }
    return sign;
}

// Multiply two groups of nbrLen limbs. The first one starts at idxFactor1
// and the second one at idxFactor2. The 2*nbrLen limb result is stored
// starting at idxFactor1. Use arrayAux as temporary storage.
// Accumulate products by result limb.
#ifdef _USING64BITS_
static void ClassicalMult2Limbs(int idxFactor1, int idxFactor2)
{
    int i = 0;
    uint32_t factor2_i, factor2_iPlus1;
    uint64_t Pr;
    M(0); M(1);
    PROLOG_MULTIPLICATION_DOUBLE;
    EPILOG_MULTIPLICATION_DOUBLE(0, 1);
    arrayAux[2] = prod_iPlus0;
    arrayAux[3] = prod_iPlus1;
}

static void ClassicalMult3Limbs(int idxFactor1, int idxFactor2)
{
    int i = 0;
    uint32_t factor2_i, factor2_iPlus1;
    uint64_t Pr;
    M(0); M(1); M(2);
    PROLOG_MULTIPLICATION_DOUBLE;
    MULT_MACRO_DOUBLE(0, 1, 2);
    EPILOG_MULTIPLICATION_DOUBLE(1, 2);
    PROLOG_MULTIPLICATION_SINGLE(2);
    MULT_MACRO_SINGLE(1, 3);
    MULT_MACRO_SINGLE(2, 4);
    EPILOG_MULTIPLICATION_SINGLE(5);
}

static void ClassicalMult4Limbs(int idxFactor1, int idxFactor2)
{
    int i;
    uint32_t factor2_i, factor2_iPlus1;
    uint64_t Pr;
    M(0); M(1); M(2); M(3);
    for (i = 0; i < 4; i += 2)
    {
        PROLOG_MULTIPLICATION_DOUBLE;
        MULT_MACRO_DOUBLE(0, 1, 2);
        MULT_MACRO_DOUBLE(1, 2, 3);
        EPILOG_MULTIPLICATION_DOUBLE(2, 3);
    }
    arrayAux[4] = prod_iPlus0;
    arrayAux[5] = prod_iPlus1;
    arrayAux[6] = prod_iPlus2;
    arrayAux[7] = prod_iPlus3;
}

static void ClassicalMult5Limbs(int idxFactor1, int idxFactor2)
{
    int i;
    uint32_t factor2_i, factor2_iPlus1;
    uint64_t Pr;
    M(0); M(1); M(2); M(3); M(4);
    for (i = 0; i < 4; i += 2)
    {
        PROLOG_MULTIPLICATION_DOUBLE;
        MULT_MACRO_DOUBLE(0, 1, 2);
        MULT_MACRO_DOUBLE(1, 2, 3);
        MULT_MACRO_DOUBLE(2, 3, 4);
        EPILOG_MULTIPLICATION_DOUBLE(3, 4);
    }
    PROLOG_MULTIPLICATION_SINGLE(4);
    MULT_MACRO_SINGLE(1, 5);
    MULT_MACRO_SINGLE(2, 6);
    MULT_MACRO_SINGLE(3, 7);
    MULT_MACRO_SINGLE(4, 8);
    EPILOG_MULTIPLICATION_SINGLE(9);
}

static void ClassicalMult6Limbs(int idxFactor1, int idxFactor2)
{
    int i;
    uint32_t factor2_i, factor2_iPlus1;
    uint64_t Pr;
    M(0); M(1); M(2); M(3); M(4); M(5);
    for (i = 0; i < 6; i += 2)
    {
        PROLOG_MULTIPLICATION_DOUBLE;
        MULT_MACRO_DOUBLE(0, 1, 2);
        MULT_MACRO_DOUBLE(1, 2, 3);
        MULT_MACRO_DOUBLE(2, 3, 4);
        MULT_MACRO_DOUBLE(3, 4, 5);
        EPILOG_MULTIPLICATION_DOUBLE(4, 5);
    }
    arrayAux[6] = prod_iPlus0;
    arrayAux[7] = prod_iPlus1;
    arrayAux[8] = prod_iPlus2;
    arrayAux[9] = prod_iPlus3;
    arrayAux[10] = prod_iPlus4;
    arrayAux[11] = prod_iPlus5;
}

static void ClassicalMult7Limbs(int idxFactor1, int idxFactor2)
{
    int i;
    uint32_t factor2_i, factor2_iPlus1;
    uint64_t Pr;
    M(0); M(1); M(2); M(3); M(4); M(5); M(6);
    for (i = 0; i < 6; i += 2)
    {
        PROLOG_MULTIPLICATION_DOUBLE;
        MULT_MACRO_DOUBLE(0, 1, 2);
        MULT_MACRO_DOUBLE(1, 2, 3);
        MULT_MACRO_DOUBLE(2, 3, 4);
        MULT_MACRO_DOUBLE(3, 4, 5);
        MULT_MACRO_DOUBLE(4, 5, 6);
        EPILOG_MULTIPLICATION_DOUBLE(5, 6);
    }
    PROLOG_MULTIPLICATION_SINGLE(6);
    MULT_MACRO_SINGLE(1, 7);
    MULT_MACRO_SINGLE(2, 8);
    MULT_MACRO_SINGLE(3, 9);
    MULT_MACRO_SINGLE(4, 10);
    MULT_MACRO_SINGLE(5, 11);
    MULT_MACRO_SINGLE(6, 12);
    EPILOG_MULTIPLICATION_SINGLE(13);
}

static void ClassicalMult8Limbs(int idxFactor1, int idxFactor2)
{
    int i;
    uint32_t factor2_i, factor2_iPlus1;
    uint64_t Pr;
    M(0); M(1); M(2); M(3); M(4); M(5); M(6); M(7);
    for (i = 0; i < 8; i += 2)
    {
        PROLOG_MULTIPLICATION_DOUBLE;
        MULT_MACRO_DOUBLE(0, 1, 2);
        MULT_MACRO_DOUBLE(1, 2, 3);
        MULT_MACRO_DOUBLE(2, 3, 4);
        MULT_MACRO_DOUBLE(3, 4, 5);
        MULT_MACRO_DOUBLE(4, 5, 6);
        MULT_MACRO_DOUBLE(5, 6, 7);
        EPILOG_MULTIPLICATION_DOUBLE(6, 7);
    }
    arrayAux[8] = prod_iPlus0;
    arrayAux[9] = prod_iPlus1;
    arrayAux[10] = prod_iPlus2;
    arrayAux[11] = prod_iPlus3;
    arrayAux[12] = prod_iPlus4;
    arrayAux[13] = prod_iPlus5;
    arrayAux[14] = prod_iPlus6;
    arrayAux[15] = prod_iPlus7;
}

static void ClassicalMult9Limbs(int idxFactor1, int idxFactor2)
{
    int i;
    uint32_t factor2_i, factor2_iPlus1;
    uint64_t Pr;
    M(0); M(1); M(2); M(3); M(4); M(5); M(6); M(7); M(8);
    for (i = 0; i < 8; i += 2)
    {
        PROLOG_MULTIPLICATION_DOUBLE;
        MULT_MACRO_DOUBLE(0, 1, 2);
        MULT_MACRO_DOUBLE(1, 2, 3);
        MULT_MACRO_DOUBLE(2, 3, 4);
        MULT_MACRO_DOUBLE(3, 4, 5);
        MULT_MACRO_DOUBLE(4, 5, 6);
        MULT_MACRO_DOUBLE(5, 6, 7);
        MULT_MACRO_DOUBLE(6, 7, 8);
        EPILOG_MULTIPLICATION_DOUBLE(7, 8);
    }
    PROLOG_MULTIPLICATION_SINGLE(8);
    MULT_MACRO_SINGLE(1, 9);
    MULT_MACRO_SINGLE(2, 10);
    MULT_MACRO_SINGLE(3, 11);
    MULT_MACRO_SINGLE(4, 12);
    MULT_MACRO_SINGLE(5, 13);
    MULT_MACRO_SINGLE(6, 14);
    MULT_MACRO_SINGLE(7, 15);
    MULT_MACRO_SINGLE(8, 16);
    EPILOG_MULTIPLICATION_SINGLE(17);
}

static void ClassicalMult10Limbs(int idxFactor1, int idxFactor2)
{
    int i;
    uint32_t factor2_i, factor2_iPlus1;
    uint64_t Pr;
    M(0); M(1); M(2); M(3); M(4); M(5); M(6); M(7); M(8); M(9);
    for (i = 0; i < 10; i += 2)
    {
        PROLOG_MULTIPLICATION_DOUBLE;
        MULT_MACRO_DOUBLE(0, 1, 2);
        MULT_MACRO_DOUBLE(1, 2, 3);
        MULT_MACRO_DOUBLE(2, 3, 4);
        MULT_MACRO_DOUBLE(3, 4, 5);
        MULT_MACRO_DOUBLE(4, 5, 6);
        MULT_MACRO_DOUBLE(5, 6, 7);
        MULT_MACRO_DOUBLE(6, 7, 8);
        MULT_MACRO_DOUBLE(7, 8, 9);
        EPILOG_MULTIPLICATION_DOUBLE(8, 9);
    }
    arrayAux[10] = prod_iPlus0;
    arrayAux[11] = prod_iPlus1;
    arrayAux[12] = prod_iPlus2;
    arrayAux[13] = prod_iPlus3;
    arrayAux[14] = prod_iPlus4;
    arrayAux[15] = prod_iPlus5;
    arrayAux[16] = prod_iPlus6;
    arrayAux[17] = prod_iPlus7;
    arrayAux[18] = prod_iPlus8;
    arrayAux[19] = prod_iPlus9;
}

static void ClassicalMult11Limbs(int idxFactor1, int idxFactor2)
{
    int i;
    uint32_t factor2_i, factor2_iPlus1;
    uint64_t Pr;
    M(0); M(1); M(2); M(3); M(4); M(5); M(6); M(7); M(8); M(9); M(10);
    for (i = 0; i < 10; i += 2)
    {
        PROLOG_MULTIPLICATION_DOUBLE;
        MULT_MACRO_DOUBLE(0, 1, 2);
        MULT_MACRO_DOUBLE(1, 2, 3);
        MULT_MACRO_DOUBLE(2, 3, 4);
        MULT_MACRO_DOUBLE(3, 4, 5);
        MULT_MACRO_DOUBLE(4, 5, 6);
        MULT_MACRO_DOUBLE(5, 6, 7);
        MULT_MACRO_DOUBLE(6, 7, 8);
        MULT_MACRO_DOUBLE(7, 8, 9);
        MULT_MACRO_DOUBLE(8, 9, 10);
        EPILOG_MULTIPLICATION_DOUBLE(9, 10);
    }
    PROLOG_MULTIPLICATION_SINGLE(10);
    MULT_MACRO_SINGLE(1, 11);
    MULT_MACRO_SINGLE(2, 12);
    MULT_MACRO_SINGLE(3, 13);
    MULT_MACRO_SINGLE(4, 14);
    MULT_MACRO_SINGLE(5, 15);
    MULT_MACRO_SINGLE(6, 16);
    MULT_MACRO_SINGLE(7, 17);
    MULT_MACRO_SINGLE(8, 18);
    MULT_MACRO_SINGLE(9, 19);
    MULT_MACRO_SINGLE(10, 20);
    EPILOG_MULTIPLICATION_SINGLE(21);
}

static void ClassicalMult12Limbs(int idxFactor1, int idxFactor2)
{
    int i;
    uint32_t factor2_i, factor2_iPlus1;
    uint64_t Pr;
    M(0); M(1); M(2); M(3); M(4); M(5); M(6); M(7); M(8); M(9);
    M(10); M(11);
    for (i = 0; i < 12; i += 2)
    {
        PROLOG_MULTIPLICATION_DOUBLE;
        MULT_MACRO_DOUBLE(0, 1, 2);
        MULT_MACRO_DOUBLE(1, 2, 3);
        MULT_MACRO_DOUBLE(2, 3, 4);
        MULT_MACRO_DOUBLE(3, 4, 5);
        MULT_MACRO_DOUBLE(4, 5, 6);
        MULT_MACRO_DOUBLE(5, 6, 7);
        MULT_MACRO_DOUBLE(6, 7, 8);
        MULT_MACRO_DOUBLE(7, 8, 9);
        MULT_MACRO_DOUBLE(8, 9, 10);
        MULT_MACRO_DOUBLE(9, 10, 11);
        EPILOG_MULTIPLICATION_DOUBLE(10, 11);
    }
    arrayAux[12] = prod_iPlus0;
    arrayAux[13] = prod_iPlus1;
    arrayAux[14] = prod_iPlus2;
    arrayAux[15] = prod_iPlus3;
    arrayAux[16] = prod_iPlus4;
    arrayAux[17] = prod_iPlus5;
    arrayAux[18] = prod_iPlus6;
    arrayAux[19] = prod_iPlus7;
    arrayAux[20] = prod_iPlus8;
    arrayAux[21] = prod_iPlus9;
    arrayAux[22] = prod_iPlus10;
    arrayAux[23] = prod_iPlus11;
}

static void ClassicalMult13Limbs(int idxFactor1, int idxFactor2)
{
    int i;
    uint32_t factor2_i, factor2_iPlus1;
    uint64_t Pr;
    M(0); M(1); M(2); M(3); M(4); M(5); M(6); M(7); M(8); M(9);
    M(10); M(11); M(12);
    for (i = 0; i < 12; i += 2)
    {
        PROLOG_MULTIPLICATION_DOUBLE;
        MULT_MACRO_DOUBLE(0, 1, 2);
        MULT_MACRO_DOUBLE(1, 2, 3);
        MULT_MACRO_DOUBLE(2, 3, 4);
        MULT_MACRO_DOUBLE(3, 4, 5);
        MULT_MACRO_DOUBLE(4, 5, 6);
        MULT_MACRO_DOUBLE(5, 6, 7);
        MULT_MACRO_DOUBLE(6, 7, 8);
        MULT_MACRO_DOUBLE(7, 8, 9);
        MULT_MACRO_DOUBLE(8, 9, 10);
        MULT_MACRO_DOUBLE(9, 10, 11);
        MULT_MACRO_DOUBLE(10, 11, 12);
        EPILOG_MULTIPLICATION_DOUBLE(11, 12);
    }
    PROLOG_MULTIPLICATION_SINGLE(12);
    MULT_MACRO_SINGLE(1, 13);
    MULT_MACRO_SINGLE(2, 14);
    MULT_MACRO_SINGLE(3, 15);
    MULT_MACRO_SINGLE(4, 16);
    MULT_MACRO_SINGLE(5, 17);
    MULT_MACRO_SINGLE(6, 18);
    MULT_MACRO_SINGLE(7, 19);
    MULT_MACRO_SINGLE(8, 20);
    MULT_MACRO_SINGLE(9, 21);
    MULT_MACRO_SINGLE(10, 22);
    MULT_MACRO_SINGLE(11, 23);
    MULT_MACRO_SINGLE(12, 24);
    EPILOG_MULTIPLICATION_SINGLE(25);
}

static void ClassicalMult14Limbs(int idxFactor1, int idxFactor2)
{
    int i;
    uint32_t factor2_i, factor2_iPlus1;
    uint64_t Pr;
    M(0); M(1); M(2); M(3); M(4); M(5); M(6); M(7); M(8); M(9);
    M(10); M(11); M(12); M(13);
    for (i = 0; i < 14; i += 2)
    {
        PROLOG_MULTIPLICATION_DOUBLE;
        MULT_MACRO_DOUBLE(0, 1, 2);
        MULT_MACRO_DOUBLE(1, 2, 3);
        MULT_MACRO_DOUBLE(2, 3, 4);
        MULT_MACRO_DOUBLE(3, 4, 5);
        MULT_MACRO_DOUBLE(4, 5, 6);
        MULT_MACRO_DOUBLE(5, 6, 7);
        MULT_MACRO_DOUBLE(6, 7, 8);
        MULT_MACRO_DOUBLE(7, 8, 9);
        MULT_MACRO_DOUBLE(8, 9, 10);
        MULT_MACRO_DOUBLE(9, 10, 11);
        MULT_MACRO_DOUBLE(10, 11, 12);
        MULT_MACRO_DOUBLE(11, 12, 13);
        EPILOG_MULTIPLICATION_DOUBLE(12, 13);
    }
    arrayAux[14] = prod_iPlus0;
    arrayAux[15] = prod_iPlus1;
    arrayAux[16] = prod_iPlus2;
    arrayAux[17] = prod_iPlus3;
    arrayAux[18] = prod_iPlus4;
    arrayAux[19] = prod_iPlus5;
    arrayAux[20] = prod_iPlus6;
    arrayAux[21] = prod_iPlus7;
    arrayAux[22] = prod_iPlus8;
    arrayAux[23] = prod_iPlus9;
    arrayAux[24] = prod_iPlus10;
    arrayAux[25] = prod_iPlus11;
    arrayAux[26] = prod_iPlus12;
    arrayAux[27] = prod_iPlus13;
}

static void ClassicalMult15Limbs(int idxFactor1, int idxFactor2)
{
    int i;
    uint32_t factor2_i, factor2_iPlus1;
    uint64_t Pr;
    M(0); M(1); M(2); M(3); M(4); M(5); M(6); M(7); M(8); M(9);
    M(10); M(11); M(12); M(13); M(14);
    for (i = 0; i < 14; i += 2)
    {
        PROLOG_MULTIPLICATION_DOUBLE;
        MULT_MACRO_DOUBLE(0, 1, 2);
        MULT_MACRO_DOUBLE(1, 2, 3);
        MULT_MACRO_DOUBLE(2, 3, 4);
        MULT_MACRO_DOUBLE(3, 4, 5);
        MULT_MACRO_DOUBLE(4, 5, 6);
        MULT_MACRO_DOUBLE(5, 6, 7);
        MULT_MACRO_DOUBLE(6, 7, 8);
        MULT_MACRO_DOUBLE(7, 8, 9);
        MULT_MACRO_DOUBLE(8, 9, 10);
        MULT_MACRO_DOUBLE(9, 10, 11);
        MULT_MACRO_DOUBLE(10, 11, 12);
        MULT_MACRO_DOUBLE(11, 12, 13);
        MULT_MACRO_DOUBLE(12, 13, 14);
        EPILOG_MULTIPLICATION_DOUBLE(13, 14);
    }
    PROLOG_MULTIPLICATION_SINGLE(14);
    MULT_MACRO_SINGLE(1, 15);
    MULT_MACRO_SINGLE(2, 16);
    MULT_MACRO_SINGLE(3, 17);
    MULT_MACRO_SINGLE(4, 18);
    MULT_MACRO_SINGLE(5, 19);
    MULT_MACRO_SINGLE(6, 20);
    MULT_MACRO_SINGLE(7, 21);
    MULT_MACRO_SINGLE(8, 22);
    MULT_MACRO_SINGLE(9, 23);
    MULT_MACRO_SINGLE(10, 24);
    MULT_MACRO_SINGLE(11, 25);
    MULT_MACRO_SINGLE(12, 26);
    MULT_MACRO_SINGLE(13, 27);
    MULT_MACRO_SINGLE(14, 28);
    EPILOG_MULTIPLICATION_SINGLE(29);
}

static void ClassicalMult16Limbs(int idxFactor1, int idxFactor2)
{
    int i;
    uint32_t factor2_i, factor2_iPlus1;
    uint64_t Pr;
    M(0); M(1); M(2); M(3); M(4); M(5); M(6); M(7); M(8); M(9);
    M(10); M(11); M(12); M(13); M(14); M(15);
    for (i = 0; i < 16; i += 2)
    {
        PROLOG_MULTIPLICATION_DOUBLE;
        MULT_MACRO_DOUBLE(0, 1, 2);
        MULT_MACRO_DOUBLE(1, 2, 3);
        MULT_MACRO_DOUBLE(2, 3, 4);
        MULT_MACRO_DOUBLE(3, 4, 5);
        MULT_MACRO_DOUBLE(4, 5, 6);
        MULT_MACRO_DOUBLE(5, 6, 7);
        MULT_MACRO_DOUBLE(6, 7, 8);
        MULT_MACRO_DOUBLE(7, 8, 9);
        MULT_MACRO_DOUBLE(8, 9, 10);
        MULT_MACRO_DOUBLE(9, 10, 11);
        MULT_MACRO_DOUBLE(10, 11, 12);
        MULT_MACRO_DOUBLE(11, 12, 13);
        MULT_MACRO_DOUBLE(12, 13, 14);
        MULT_MACRO_DOUBLE(13, 14, 15);
        EPILOG_MULTIPLICATION_DOUBLE(14, 15);
    }
    arrayAux[16] = prod_iPlus0;
    arrayAux[17] = prod_iPlus1;
    arrayAux[18] = prod_iPlus2;
    arrayAux[19] = prod_iPlus3;
    arrayAux[20] = prod_iPlus4;
    arrayAux[21] = prod_iPlus5;
    arrayAux[22] = prod_iPlus6;
    arrayAux[23] = prod_iPlus7;
    arrayAux[24] = prod_iPlus8;
    arrayAux[25] = prod_iPlus9;
    arrayAux[26] = prod_iPlus10;
    arrayAux[27] = prod_iPlus11;
    arrayAux[28] = prod_iPlus12;
    arrayAux[29] = prod_iPlus13;
    arrayAux[30] = prod_iPlus14;
    arrayAux[31] = prod_iPlus15;
}

#endif

static void ClassicalMult(int idxFactor1, int idxFactor2, int nbrLen) {
#ifdef _USING64BITS_
    uint64_t product;

    switch (nbrLen)     /* currently must have nbrLen <= 16 */
    {
    case 1:
        product = (int64_t)arr[idxFactor1] * arr[idxFactor2];
        // Store least significant limb.
        arr[idxFactor1] = (int32_t)product & MAX_INT_NBR;
        // Store most significant limb.
        arr[idxFactor2] = (int32_t)(product >> BITS_PER_GROUP) & MAX_INT_NBR;
        return;
    case 2:
        ClassicalMult2Limbs(idxFactor1, idxFactor2);  /* arrayAux = Factor1 * Factor2 */
        break;
    case 3:
        ClassicalMult3Limbs(idxFactor1, idxFactor2);  /* arrayAux = Factor1 * Factor2 */
        break;
    case 4:
        ClassicalMult4Limbs(idxFactor1, idxFactor2);  /* arrayAux = Factor1 * Factor2 */
        break;
    case 5:
        ClassicalMult5Limbs(idxFactor1, idxFactor2);  /* arrayAux = Factor1 * Factor2 */
        break;
    case 6:
        ClassicalMult6Limbs(idxFactor1, idxFactor2);  /* arrayAux = Factor1 * Factor2 */
        break;
    case 7:
        ClassicalMult7Limbs(idxFactor1, idxFactor2);  /* arrayAux = Factor1 * Factor2 */
        break;
    case 8:
        ClassicalMult8Limbs(idxFactor1, idxFactor2);  /* arrayAux = Factor1 * Factor2 */
        break;
    case 9:
        ClassicalMult9Limbs(idxFactor1, idxFactor2);   /* arrayAux = Factor1 * Factor2 */
        break;
    case 10:
        ClassicalMult10Limbs(idxFactor1, idxFactor2);  /* arrayAux = Factor1 * Factor2 */
        break;
    case 11:
        ClassicalMult11Limbs(idxFactor1, idxFactor2);  /* arrayAux = Factor1 * Factor2 */
        break;
    case 12:
        ClassicalMult12Limbs(idxFactor1, idxFactor2);  /* arrayAux = Factor1 * Factor2 */
        break;
    case 13:
        ClassicalMult13Limbs(idxFactor1, idxFactor2);  /* arrayAux = Factor1 * Factor2 */
        break;
    case 14:
        ClassicalMult14Limbs(idxFactor1, idxFactor2);  /* arrayAux = Factor1 * Factor2 */
        break;
    case 15:
        ClassicalMult15Limbs(idxFactor1, idxFactor2);  /* arrayAux = Factor1 * Factor2 */
        break;
    case 16:
        ClassicalMult16Limbs(idxFactor1, idxFactor2);  /* arrayAux = Factor1 * Factor2 */
        break;

    default: std::abort();   // should never get here
    }
#else
    limb *ptrFactor1, *ptrFactor2;
    int prodCol, fact1Col;
    double dRangeLimb = (double)LIMB_RANGE;
    double dInvRangeLimb = 1 / dRangeLimb;
    int low = 0;              // Low limb of sums of multiplications.
    double dAccumulator = 0;  // Approximation to the sum of multiplications.
    int factor1, factor2;
    for (prodCol = 0; prodCol < 2 * nbrLen - 1; prodCol++)
    {    // Process each limb of product (least to most significant limb).
        if (prodCol < nbrLen)
        {   // Processing first half (least significant) of product.
            ptrFactor2 = &arr[idxFactor2 + prodCol];
            ptrFactor1 = &arr[idxFactor1];
            fact1Col = prodCol;
        }
        else
        {  // Processing second half (most significant) of product.
            ptrFactor2 = &arr[idxFactor2 + nbrLen - 1];
            ptrFactor1 = &arr[idxFactor1 + prodCol - nbrLen + 1];
            fact1Col = 2 * (nbrLen - 1) - prodCol;
        }
        for (; fact1Col >= 0; fact1Col--)
        {
            factor1 = (ptrFactor1++)->x;
            factor2 = (ptrFactor2--)->x;
            low += factor1 * factor2;
            dAccumulator += (double)factor1 * (double)factor2;
        }
        low &= MAX_VALUE_LIMB;    // Trim extra bits.
        arrayAux[prodCol] = low;
        // Subtract or add 0x20000000 so the multiplication by dVal is not nearly an integer.
        // In that case, there would be an error of +/- 1.
        if (low < HALF_INT_RANGE)
        {
            dAccumulator = floor(dAccumulator * dInvRangeLimb + 0.25);
        }
        else
        {
            dAccumulator = floor(dAccumulator * dInvRangeLimb - 0.25);
        }
        low = (unsigned int)(dAccumulator - floor(dAccumulator * dInvRangeLimb) * dRangeLimb);
    }
    arrayAux[prodCol] = low;
#endif
    /* copy product back into factor1 and Factor2 */
    std::memcpy(&arr[idxFactor1], &arrayAux[0], 2 * nbrLen * sizeof(limb));
    return;
}

// Recursive Karatsuba function.
static void Karatsuba(int idxFactor1, const int nbrLen, int diffIndex)
{
    /* check that array index is within bounds */
    assert(idxFactor1 + nbrLen * 2 <= sizeof(arr) / sizeof(arr[0]));

    int idxFactor2 = idxFactor1 + nbrLen;
    int i;
    unsigned int carry1First, carry1Second;
    unsigned int carry2Second;
    limb *ptrResult, *ptrHigh, tmp;
    int middle;
    int sign;
    int halfLength;
    
    // Check if one of the factors is equal to zero.
    if (BigNbrIsZero (&arr[idxFactor1], nbrLen) || BigNbrIsZero(&arr[idxFactor2], nbrLen))
    {    // One of the factors is equal to zero.
        for (i = nbrLen - 1; i >= 0; i--) 	{
            arr[idxFactor1 + i] = arr[idxFactor2 + i] = 0;
        }
        return;
    }

    if (nbrLen <= KARATSUBA_CUTOFF) {
        // Below cutoff: perform standard classical multiplcation.
        ClassicalMult(idxFactor1, idxFactor2, nbrLen);
        return;
    }

    // Length > KARATSUBA_CUTOFF: Use Karatsuba multiplication.
    // It uses three half-length multiplications instead of four.
    //  x*y = (xH*b + xL)*(yH*b + yL)
    //  x*y = (b + 1)*(xH*yH*b + xL*yL) + (xH - xL)*(yL - yH)*b
    // The length of b is stored in variable halfLength.
    // Since the absolute values of (xH - xL) and (yL - yH) fit in
    // a single limb, there will be no overflow.

    // At this moment the order is: xL, xH, yL, yH.
    // Exchange high part of first factor with low part of 2nd factor.
    halfLength = nbrLen >> 1;
    for (i = idxFactor1 + halfLength; i<idxFactor2; i++)
    {
        tmp = arr[i];
        arr[i] = arr[i + halfLength];
        arr[i + halfLength] = tmp;
    }
    // At this moment the order is: xL, yL, xH, yH.
    // Get absolute values of (xH-xL) and (yL-yH) and the signs.
    sign = absSubtract(idxFactor1, idxFactor2, diffIndex, halfLength);
    sign ^= absSubtract(idxFactor2 + halfLength, idxFactor1 + halfLength,
        diffIndex + halfLength, halfLength);
    middle = diffIndex;
    diffIndex += nbrLen;
    Karatsuba(idxFactor1, halfLength, diffIndex); // Multiply both low parts.
    Karatsuba(idxFactor2, halfLength, diffIndex); // Multiply both high parts.
    Karatsuba(middle, halfLength, diffIndex);     // Multiply the differences.
                                                  // Process all carries at the end.
                                                  // Obtain (b+1)(xH*yH*b + xL*yL) = xH*yH*b^2 + (xL*yL+xH*yH)*b + xL*yL
                                                  // The first and last terms are already in correct locations.
    ptrResult = &arr[idxFactor1 + halfLength];
    carry1First = carry1Second = carry2Second = 0;
    for (i = halfLength; i > 0; i--)
    {
        // The sum of three ints overflows an unsigned int variable,
        // so two adds are required. Also carries must be separated in
        // order to avoid overflow:
        // 00000001 + 7FFFFFFF + 7FFFFFFF = FFFFFFFF
        unsigned int accum1Lo = carry1First + *ptrResult + *(ptrResult + halfLength);
        unsigned int accum2Lo;
        carry1First = accum1Lo >> BITS_PER_GROUP;
        accum2Lo = carry2Second + (accum1Lo & MAX_VALUE_LIMB) +
            *(ptrResult - halfLength);
        carry2Second = accum2Lo >> BITS_PER_GROUP;
        accum1Lo = carry1Second + (accum1Lo & MAX_VALUE_LIMB) +
            *(ptrResult + nbrLen);
        carry1Second = accum1Lo >> BITS_PER_GROUP;
        *(ptrResult + halfLength) = accum1Lo & MAX_VALUE_LIMB;
        *ptrResult = accum2Lo & MAX_VALUE_LIMB;
        ptrResult++;
    }
    *(ptrResult + halfLength) += carry1First + carry1Second;
    *ptrResult += carry1First + carry2Second;
    // Process carries.
    ptrResult = &arr[idxFactor1];
    carry1First = 0;
    for (i = 2 * nbrLen; i > 0; i--)
    {
        carry1First += *ptrResult;
        *(ptrResult++) = carry1First & MAX_VALUE_LIMB;
        carry1First >>= BITS_PER_GROUP;
    }
    // Compute final product.
    ptrHigh = &arr[middle];
    ptrResult = &arr[idxFactor1 + halfLength];
    if (sign != 0)
    {            // (xH-xL) * (yL-yH) is negative.
        int borrow = 0;
        for (i = nbrLen; i > 0; i--)
        {
            borrow += *ptrResult - *(ptrHigh++);
            *(ptrResult++) = borrow & MAX_VALUE_LIMB;
            borrow >>= BITS_PER_GROUP;
        }
        for (i = halfLength; i > 0; i--)
        {
            borrow += *ptrResult;
            *(ptrResult++) = borrow & MAX_VALUE_LIMB;
            borrow >>= BITS_PER_GROUP;
        }
    }
    else
    {            // (xH-xL) * (yL-yH) is positive or zero.
        unsigned int carry = 0;
        for (i = nbrLen; i > 0; i--)
        {
            carry += (unsigned int)*ptrResult + (unsigned int)*(ptrHigh++);
            *(ptrResult++) = (int)(carry & MAX_VALUE_LIMB);
            carry >>= BITS_PER_GROUP;
        }
        for (i = halfLength; i > 0; i--)
        {
            carry += (unsigned int)*ptrResult;
            *(ptrResult++) = (int)(carry & MAX_VALUE_LIMB);
            carry >>= BITS_PER_GROUP;
        }
    }
}