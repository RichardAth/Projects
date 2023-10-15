 //
// This file is part of Alpertron Calculators.
//
// Copyright 2017-2021 Dario Alejandro Alpern
//
// Alpertron Calculators is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Alpertron Calculators is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Alpertron Calculators.  If not, see <http://www.gnu.org/licenses/>.
//
#include <cassert>
#include "pch.h"

#include "quadmodLL.h"
//#include "commonstruc.h"

extern BigInteger MontgomeryMultR1BI, MontgomeryMultR2BI;

struct stQuad
{
    BigInteger Solution1[400];
    BigInteger Solution2[400];
    BigInteger Increment[400];
};

static const int maxFactors = 400;  /* modulus cannot have > 400 unique factors */
static BigInteger Quadr;
static BigInteger Linear;
static BigInteger Const;
static BigInteger sqrRoot;
static BigInteger ValAOdd;
static BigInteger ValBOdd;
static BigInteger ValCOdd;
static BigInteger Mult;
static BigInteger currentSolution;
static BigInteger discriminant;
static BigInteger SqrtDisc;
BigInteger prime;
static BigInteger bigBase;
static BigInteger K1;
static BigInteger L;
static BigInteger Q;
static BigInteger V;
static int Exponents[maxFactors];
static BigInteger Aux[maxFactors];
static BigInteger tmp1;
static BigInteger tmp2;
static BigInteger ValNn;
static int nbrFactors;
static bool sol1Invalid;
static bool sol2Invalid;
BigInteger LastModulus;
static BigInteger Aux0;
static BigInteger Aux1;
static BigInteger Aux2;
static BigInteger* pGcdAll;
const static BigInteger* pValNn;
fList Nfactors;

/* function pointers */
pSolution SolutionP = nullptr;
static pShowSolutionsModPrime ShowSolutionsModPrime = nullptr;;
static pShowNoSolsModPrime ShowNoSolsModPrime = nullptr;

struct sFactors astFactorsMod[maxFactors+1];
int factorsMod[20000];

extern char* ptrOutput;
extern int SolNbr;

struct stQuad quad;

#undef min  /* remove min defined in windows.h */

/* mainly for debugging */
void printFactors(sFactors flist[]) {
    int nbrFactors = std::min(flist[0].multiplicity, (int)maxFactors);
    printf_s("factor list: \n");
    for (int i = 1; i <= nbrFactors; i++) {
        printf_s("% 4d ", i);
        IntArray2BigInteger(flist[i].ptrFactor, &prime);
        PrintBigInteger(&prime, 0);
        if (flist[i].multiplicity != 1) {
            printf_s("^%d", flist[i].multiplicity);
        }
        /* values below only set for larger factors*/
        printf_s(" upperBound = %d, Type = %d \n", flist[i].upperBound, flist[i].Type);
    }
    std::putchar('\n');
}

/* the same factor list is returned in both factorsMod and astFactorsMod,
in different formats 
Uses global variable Aux0 */
static void factor(const BigInteger* pValN, int factorsMod[], sFactors astFactorsMod[]) {
    Znum N;
    int numFactors;
    int* pFactorsMod = factorsMod;
    BigtoZ(N, *pValN);   /* convert N to Znum */
    bool rv = factorise(N, Nfactors, nullptr);
    assert(rv);
    numFactors = (int)Nfactors.fsize();
    int numLimbs;
    astFactorsMod[0].multiplicity = numFactors; 
    /* if necessary, reduce numFactors so that array bound is not exceeded 
    This would be a problem as array of factors is truncated */
    numFactors = std::min((int)maxFactors, numFactors);

    for (int i = 1; i <= numFactors; i++) {
        ZtoBig(Aux0, Nfactors.f[i - 1].Factor);
        numLimbs = Aux0.nbrLimbs;
        *pFactorsMod = numLimbs;
        pFactorsMod++;
        BigIntegerToLimbs(pFactorsMod, Aux0, numLimbs);
        astFactorsMod[i].ptrFactor = pFactorsMod - 1;
        pFactorsMod += numLimbs;
        //BigIntegerToLimbs(astFactorsMod[i].ptrFactor, Aux0, numLimbs);
        astFactorsMod[i].multiplicity = Nfactors.f[i - 1].exponent;
    }
}

// Use Chinese remainder theorem to obtain the solutions.
/* for each solution <solution()> is called.
The input parameters are in quad, astFactorsMod, Aux, exponents,  GcdAll, pValNn, 
... */
static void PerformChineseRemainderTheorem(const BigInteger *pValA, const BigInteger *pValB, 
    const BigInteger *pValC, const BigInteger *pValN)
{
    int T1;
    int expon;
    do {
        const struct sFactors* pstFactor;
        Aux[0] = quad.Increment[0] * (Exponents[0] / 2); // multint(&Aux[0], &quad.Increment[0], Exponents[0] / 2);
        if ((Exponents[0] & 1) != 0)
        {
            Aux[0] += quad.Solution2[0]; // BigIntAdd(&Aux[0],  &quad.Solution2[0], &Aux[0]);
        }
        else
        {
            Aux[0] += quad.Solution1[0]; // BigIntAdd(&Aux[0], &quad.Solution1[0], &Aux[0]);
        }
        currentSolution = Aux[0];  // CopyBigInt(&currentSolution, &Aux[0]);
        nbrFactors = astFactorsMod[0].multiplicity;
        pstFactor = &astFactorsMod[1];
        IntArray2BigInteger(pstFactor->ptrFactor, &prime);
        (void)BigIntPowerIntExp(prime, pstFactor->multiplicity, Mult);  /* Mult = prime^exp */
        for (T1 = 1; T1 < nbrFactors; T1++) {
            pstFactor++;
            if (pstFactor->multiplicity == 0) {
                Aux[T1] = 0;   // intToBigInteger(&Aux[T1], 0);
                continue;
            }
            expon = Exponents[T1];
            Aux[T1] = quad.Increment[T1] * (expon / 2); // multint(&Aux[T1], &quad.Increment[T1], expon / 2);
            if ((expon & 1) != 0) {
                Aux[T1] += quad.Solution2[T1]; // BigIntAdd(&Aux[T1], &quad.Solution2[T1], &Aux[T1]);
            }
            else {
                Aux[T1] += quad.Solution1[T1]; //BigIntAdd(&Aux[T1], &quad.Solution1[T1], &Aux[T1]);
            }
            NumberLength = *pstFactor->ptrFactor;
            IntArray2BigInteger(pstFactor->ptrFactor, &prime);
            (void)BigIntPowerIntExp(prime, pstFactor->multiplicity, K1);  /* K1 = prime^exp */
            prime = K1;    // CopyBigInt(&prime, &K1);
            for (int E = 0; E < T1; E++)
            {
                int NumberLengthBytes;
                Q = Aux[T1] - Aux[E];  // BigIntSubt(&Aux[T1], &Aux[E], &Q);
                IntArray2BigInteger(astFactorsMod[E + 1].ptrFactor, &bigBase);
                (void)BigIntPowerIntExp(bigBase, astFactorsMod[E + 1].multiplicity, L);
                NumberLength = prime.nbrLimbs;
                NumberLengthBytes = NumberLength * (int)sizeof(limb);
                std::memcpy(TestNbr, prime.limbs, NumberLengthBytes);
                TestNbr[NumberLength] = 0;
                GetMontgomeryParms(NumberLength);
                BigIntModularDivision(&Q, &L, &prime, &Aux[T1]);
            }
            L = Aux[T1] % prime;   //(void)BigIntRemainder(&Aux[T1], &prime, &L);
            Aux[T1] = L;           // CopyBigInt(&Aux[T1], &L);
            // currentSolution <- Aux[T1] * Mult + currentSolution
            L = Aux[T1] * Mult;   //(void)BigIntMultiply(&Aux[T1], &Mult, &L);
            currentSolution += L; //BigIntAdd(&currentSolution, &L, &currentSolution);
            Mult *= K1;           // (void)BigIntMultiply(&K1, &Mult, &Mult);
        }   /* end for */
        V = 0;  // intToBigInteger(&V, 0);
        K1 = V - *pGcdAll; // BigIntSubt(&V, pGcdAll, &K1);
        // Perform loop while V < GcdAll.
        while (K1.sign == SIGN_NEGATIVE)
        {
            // The solution is V*ValNn + currentSolution
            K1 = V * *pValNn;      // (void)BigIntMultiply(&V, pValNn, &K1);
            K1 += currentSolution; // BigIntAdd(&K1, &currentSolution, &K1);
            SolutionP(&K1, pValA, pValB, pValC, pValN);
            V++;                   // addbigint(&V, 1);      // V <- V + 1
            K1 = V - *pGcdAll;     // BigIntSubt(&V, pGcdAll, &K1);
        }
        for (T1 = nbrFactors - 1; T1 >= 0; T1--)
        {
            IntArray2BigInteger(astFactorsMod[T1 + 1].ptrFactor, &bigBase);
            BigIntPowerIntExp(bigBase, astFactorsMod[T1 + 1].multiplicity, prime);
            K1 = quad.Solution1[T1] - quad.Solution2[T1];
            // BigIntSubt(&common.quad.Solution1[T1], &common.quad.Solution2[T1], &K1);
            if (K1 == 0)    //((K1.nbrLimbs == 1) && (K1.limbs[0] == 0))
            {     // common.quad.Solution1[T1] == common.quad.Solution2[T1]
                Exponents[T1] += 2;
            }
            else
            {     // common.quad.Solution1[T1] != common.quad.Solution2[T1]
                Exponents[T1]++;
            }
            L = quad.Increment[T1] * Exponents[T1]; // multadd(&L, Exponents[T1], &quad.Increment[T1], 0);   // 
            K1 = prime * 2;  // multadd(&K1, 2, &prime, 0);               K1 <- 2 * prime
            L -= K1;        // BigIntSubt(&L, &K1, &L);
            if (L.sign == SIGN_NEGATIVE)  {
                break;
            }
            Exponents[T1] = 0;
        }   /* end for */
    } while (T1 >= 0);
}

/* adjust number of limbs (remove leading zero limbs) */
static void setNbrLimbs(BigInteger* pBigNbr) {
    pBigNbr->nbrLimbs = NumberLength;
    pBigNbr->sign = SIGN_POSITIVE;
    while (pBigNbr->nbrLimbs > 1)   {
        if (pBigNbr->limbs[pBigNbr->nbrLimbs - 1] != 0)  {
            break;
        }
        pBigNbr->nbrLimbs--;
    }
}

// Solve Bx + C = 0 (mod N).
// uses global variables Aux[0], NumberLength, Common.Quad.Increment, Aux2
// quad.solution1, quad.solution2, factorsMod, astFactorsMod, etc
static void SolveModularLinearEquation(const BigInteger * pValA, const BigInteger* pValB,
    const BigInteger* pValC, const BigInteger* pValN)
{
    int NumberLengthBytes;
    int powerOf2;
    int solutionNbr = 0;
    int* ptrFactorsMod = factorsMod;
    struct sFactors* pstFactor = &astFactorsMod[1];
    BigInteger* ptrSolution1 = quad.Solution1;
    BigInteger* ptrSolution2 = quad.Solution2;
    ValNn = *pValN;

    if (verbose > 1) {
        printf_s("SolveModularLinearEquation: b = ");
        PrintBigInteger(pValB, 0);
        printf_s(" c = ");
        PrintBigInteger(pValC, 0);
        printf_s(" n = ");
        PrintBigInteger(pValN, 0);
        std::putchar('\n');
    }
    BigIntGcd(*pValB, *pValN, Aux[0]);
    if (Aux[0] != 1)  //((Aux[0].nbrLimbs != 1) || (Aux[0].limbs[0] != 1))
    {         // ValB and ValN are not coprime. Go out.
        return;
    }
    // Calculate z <- -ValC / ValB (mod ValN)
    // Modular division routines used work for power of 2 or odd numbers.
    // This requires to compute the quotient in two steps.
    // N = r*2^k (r = odd)
    DivideBigNbrByMaxPowerOf2(&powerOf2, ValNn.limbs, &ValNn.nbrLimbs);
   
    NumberLength = ValNn.nbrLimbs;
    NumberLengthBytes = NumberLength * (int)sizeof(limb);
    if (ValNn != 1)   //((pValN->nbrLimbs != 1) || (pValN->limbs[0] != 1))
    {       // ValN is not 1.
        quad.Increment[solutionNbr] = ValNn;  // CopyBigInt(&quad.Increment[solutionNbr], pValN);
        Exponents[solutionNbr] = 1;
        std::memcpy(TestNbr, ValNn.limbs, NumberLengthBytes);
        TestNbr[NumberLength] = 0;
        // Perform division using odd modulus r.
        GetMontgomeryParms(NumberLength);
        // ptrSolution1 <- ValC / |ValB|
        BigIntModularDivision(pValC, pValB, &ValNn, ptrSolution1);
        // ptrSolution1 <- -ValC / ValB (mod N)
        if (*ptrSolution1 != 0)     //(!BigIntIsZero(ptrSolution1))
        {
            *ptrSolution1 = *pValN - *ptrSolution1; // BigIntSubt(pValN, ptrSolution1, ptrSolution1);
        }
        *ptrSolution2 = *ptrSolution1; //BigIntCopy(ptrSolution2, ptrSolution1);
        BigInteger2IntArray(ptrFactorsMod, pValN);
        pstFactor->ptrFactor = ptrFactorsMod;
        pstFactor->multiplicity = 1;
        pstFactor++;    /* advance to next factor */
        ptrFactorsMod += *ptrFactorsMod;  /* advance pointer to next factor */
        ptrFactorsMod++;
        ptrSolution1++;
        ptrSolution2++;
        solutionNbr++;
    }
    // Perform division using power of 2.
    if (powerOf2 > 0)
    {
        BigIntPowerOf2(ptrSolution1, powerOf2);  /* multiply Solution1 by 2^powerOf2 */
        quad.Increment[solutionNbr] = *ptrSolution1; // CopyBigInt(&quad.Increment[solutionNbr], ptrSolution1);
        Exponents[solutionNbr] = 1;
        BigInteger2IntArray(ptrFactorsMod, ptrSolution1); /* FactorsMod = Solution1 */
        pstFactor->ptrFactor = ptrFactorsMod;
        pstFactor->multiplicity = 1;
        GetMontgomeryParmsPowerOf2(powerOf2);
        // Use Aux2 as a temporary area.
        // ptrSolution1 <- 1 / |ValB|
        ComputeInversePower2(pValB->limbs, ptrSolution1->limbs, Aux2.limbs);
        // ptrSolution1 <- |ValC| / |ValB|
        modmult(ptrSolution1->limbs, pValC->limbs, ptrSolution1->limbs);
        NumberLengthBytes = NumberLength * (int)sizeof(int);
        // ptrSolution1 <- -ValC / ValB
        if (pValB->sign == pValC->sign)
        {  /* change  sign of Solution1 */
            std::memset(Aux2.limbs, 0, NumberLengthBytes);
            SubtractBigNbr(Aux2.limbs, ptrSolution1->limbs, ptrSolution1->limbs, NumberLength);
        }
        // Discard bits outside number in most significant limb.
        ptrSolution1->limbs[NumberLength - 1] &= (1 << (powerOf2 % BITS_PER_GROUP)) - 1;
        ptrSolution1->nbrLimbs = NumberLength;
        ptrSolution1->sign = SIGN_POSITIVE;
        *ptrSolution2 = *ptrSolution1; // CopyBigInt(ptrSolution2, ptrSolution1);
        solutionNbr++;
    }
    astFactorsMod[0].multiplicity = solutionNbr;
    if (verbose > 1) {
        printf_s("solution = ");
        PrintBigInteger(ptrSolution1, 0);
        printf_s("\nquad.Increment[0] = ");
        PrintBigInteger(&quad.Increment[0], 0);
        printf_s("\nquad.Increment[1] = ");
        PrintBigInteger(&quad.Increment[1], 0);
        printf_s("\nsolutionNbr = %d \n", solutionNbr);
    } 
    PerformChineseRemainderTheorem(pValA, pValB, pValC, pValN);
}

// Compute sqrRoot <- sqrt(ValCOdd) mod 2^expon.
// To compute the square root, compute the inverse of sqrt,
// so only multiplications are used.
// f(x) = invsqrt(x), f_{n+1}(x) = f_n * (3 - x*f_n^2)/2
// uses external variables sqrRoot, ValCOdd, tmp1, tmp2
static void ComputeSquareRootModPowerOf2(int expon, int bitsCZero)
{
    int lenBytes;
    int correctBits;
    int nbrLimbs;
    int rem;
    // First approximation to inverse of square root.
    // If value is ...0001b, the inverse of square root is ...01b.
    // If value is ...1001b, the inverse of square root is ...11b.
    sqrRoot.limbs[0] = (((ValCOdd.limbs[0] & 15) == 1) ? 1 : 3);
    correctBits = 2;
    nbrLimbs = 1;
    while (correctBits < expon)
    {   // Compute f(x) = invsqrt(x), f_{n+1}(x) = f_n * (3 - x*f_n^2)/2
        correctBits *= 2;
        nbrLimbs = (correctBits / BITS_PER_GROUP) + 1;
        multiply(sqrRoot.limbs, sqrRoot.limbs, tmp2.limbs, nbrLimbs, nullptr); // tmp2 = sqrRoot*sqrRoot
        multiply(tmp2.limbs, ValCOdd.limbs, tmp1.limbs, nbrLimbs, nullptr);   //tmp1 = tmp2 * ValCOdd
        ChSignBigNbr(tmp1.limbs, nbrLimbs);    /* tmp1 = -tmp1 */
        lenBytes = nbrLimbs * (int)sizeof(limb);
        std::memset(tmp2.limbs, 0, lenBytes);  /* tmp2 = 0 */
        tmp2.limbs[0] = 3;                      /* tmp2 = 3 */
        AddBigNbr(tmp2.limbs, tmp1.limbs, tmp2.limbs, nbrLimbs);  /* tmp2 += tmp1 */
        multiply(tmp2.limbs, sqrRoot.limbs, tmp1.limbs, nbrLimbs, nullptr);  /* tmp1 = tmp2*sqrRoot */
        std::memcpy(sqrRoot.limbs, tmp1.limbs, lenBytes);     /* sqrRoot = tmp1 */
        DivBigNbrByInt(tmp1.limbs, 2, sqrRoot.limbs, nbrLimbs, &rem); /* sqrRoot = Tmp1/2 */
        correctBits--;
    }
    // Get square root of ValCOdd from its inverse by multiplying by ValCOdd.
    multiply(ValCOdd.limbs, sqrRoot.limbs, tmp1.limbs, nbrLimbs, nullptr);
    lenBytes = nbrLimbs * (int)sizeof(limb);
    std::memcpy(sqrRoot.limbs, tmp1.limbs, lenBytes);  /* sqrRoot = tmp1*/
    setNbrLimbs(&sqrRoot);
    if ((bitsCZero / 2) > 0)
        sqrRoot <<= (bitsCZero / 2);
    //for (int ctr = 0; ctr < (bitsCZero / 2); ctr++)
    //{  /* sqrRoot *= 2^(bitsCZero / 2) */
    //    sqrRoot *= 2; //BigIntMultiplyBy2(&sqrRoot);
    //}
}

// Find quadratic solution of Quadr*x^2 + Linear*x + Const = 0 (mod 2^exponent)
// when Quadr is even and Linear is odd. In this case there is a unique solution.
// uses global variables Quadr, Linear, Const, Aux0, Aux1, Aux2
// see Lemma 3, and Algorithm 1
static void findQuadraticSolution(BigInteger* pSolution, int exponent)
{
    int bytesLen;
    int expon = exponent;
    int bitMask = 1;
    limb* ptrSolution = pSolution->limbs;
    BigIntPowerOf2(&Aux0, expon);    /* Aux0 = 2^exponent */
    bytesLen = Aux0.nbrLimbs * (int)sizeof(limb);
    std::memset(pSolution->limbs, 0, bytesLen);
    while (expon > 0)
    {
        expon--;
        BigIntPowerOf2(&Aux2, expon);   /* Aux2 *= 2^expon */
        Aux2--; //addbigint(&Aux2, -1);   // Aux2 <- 2^expon -1 (used as a bitmask)
        if ((Const.limbs[0] & 1) != 0)
        {        // Const is odd.
            *ptrSolution |= bitMask;
            // Const <- Quadr/2 + floor(Linear/2) + floor(Const/2) + 1
            if (Const.sign == SIGN_NEGATIVE)
            {
                Const--;   // addbigint(&Const, -1);
            }
            BigIntDivideBy2(&Const);          // floor(Const/2)
            Const++;   // addbigint(&Const, 1);             // floor(Const/2) + 1
            Aux1 = Linear;     // CopyBigInt(&Aux1, &Linear);
            if (Aux1.sign == SIGN_NEGATIVE)
            {
                Aux1--;  // addbigint(&Aux1, -1);
            }
            BigIntDivideBy2(&Aux1);           // floor(Linear/2)
            Const += Aux1;   // BigIntAdd(&Const, &Aux1, &Const);
            Aux1 = Quadr;     // CopyBigInt(&Aux1, &Quadr);
            BigIntDivideBy2(&Aux1);            // Quadr/2
            Const += Aux1;   // BigIntAdd(&Const, &Aux1, &Const);

            // Linear <- 2*Quadr + Linear and Quadr <- 2*Quadr.
            Quadr *= 2; //BigIntMultiplyBy2(&Quadr);         // Quadr*2
            Linear += Quadr;  // BigIntAdd(&Linear, &Quadr, &Linear);
            Linear &= Aux2; // BigIntAnd(&Linear, &Aux2, &Linear);   // Reduce mod 2^expon
        }
        else
        {        // Const is even.
            BigIntDivideBy2(&Const);           // Const/2
            Quadr *= 2; // BigIntMultiplyBy2(&Quadr);         // Quadr*2
        }
        Const &= Aux2; // (&Const, &Aux2, &Const);    // Reduce mod 2^expon
        Quadr &= Aux2; // BigIntAnd(&Quadr, &Aux2, &Quadr);    // Reduce mod 2^expon
        bitMask *= 2;
        if (bitMask < 0)  {  /* bitmask overflow?*/
            bitMask = 1;
            ptrSolution++;   /* if so, move to next limb */
        }
    }
    NumberLength = Aux0.nbrLimbs;
    setNbrLimbs(pSolution);
}

// Solve Ax^2 + Bx + C = 0 (mod 2^expon).
static bool SolveQuadraticEqModPowerOf2(int exponent, int factorIndex,
    const BigInteger* pValA, const BigInteger* pValB, const BigInteger* pValC, 
    BigInteger& prime) {
    int expon = exponent;
    int bitsAZero;
    int bitsBZero;
    int bitsCZero;

    if (verbose > 1) {
        printf_s("SolveQuadraticEqModPowerOf2: expon = %d, factorIndex = %d \n a = ", exponent, factorIndex);
        PrintBigInteger(pValA, 0);
        printf_s("\n b = ");
        PrintBigInteger(pValB, 0);
        printf_s("\n c = ");
        PrintBigInteger(pValC, 0);
        printf_s("\n prime = ");
        PrintBigInteger(&prime, 0);
        std::putchar('\n');
    }

    // ax^2 + bx + c = 0 (mod 2^expon)
    // This follows the paper Complete solving the quadratic equation mod 2^n
    // of Dehnavi, Shamsabad and Rishakani.
    // https://arxiv.org/pdf/1711.03621.pdf
    // Get odd part of A, B and C and number of bits to zero.
    ValAOdd = *pValA;   // CopyBigInt(&ValAOdd, pValA);
    DivideBigNbrByMaxPowerOf2(&bitsAZero, ValAOdd.limbs, &ValAOdd.nbrLimbs);
    ValBOdd = *pValB;   // CopyBigInt(&ValBOdd, pValB);
    DivideBigNbrByMaxPowerOf2(&bitsBZero, ValBOdd.limbs, &ValBOdd.nbrLimbs);
    ValCOdd = *pValC;   // CopyBigInt(&ValCOdd, pValC);
    DivideBigNbrByMaxPowerOf2(&bitsCZero, ValCOdd.limbs, &ValCOdd.nbrLimbs);
    if ((bitsAZero > 0) && (bitsBZero > 0) && (bitsCZero > 0))
    {
        int minExpon = bitsAZero;
        if (minExpon < bitsBZero)   {
            minExpon = bitsBZero;
        }
        if (minExpon < bitsCZero)  {
            minExpon = bitsCZero;
        }
        bitsAZero -= minExpon;
        bitsBZero -= minExpon;
        bitsCZero -= minExpon;
        expon -= minExpon;
    }
    if (((bitsAZero == 0) && (bitsBZero == 0) && (bitsCZero == 0)) ||
        ((bitsAZero > 0) && (bitsBZero > 0) && (bitsCZero == 0)))
    {
        return false;   // No solutions, so go out. (Lemma 2)
    }
    if ((bitsAZero == 0) && (bitsBZero > 0))
    {  /* see corollary 1 */
        // The solution in this case requires square root.
        // compute s = ((b/2)^2 - a*c)/a^2, q = odd part of s,
        // r = maximum exponent of power of 2 that divides s.
        tmp1 = *pValB;  // CopyBigInt(&tmp1, pValB);
        BigIntDivideBy2(&tmp1);
        tmp1 *= tmp1;   // (void)BigIntMultiply(&tmp1, &tmp1, &tmp1);  // (b/2)^2
        tmp2 = *pValA * *pValC;     // (void)BigIntMultiply(pValA, pValC, &tmp2);  // a*c
        tmp1 -= tmp2; // BigIntSubt(&tmp1, &tmp2, &tmp1);      // (b/2)^2 - a*c
        BigIntPowerOf2(&K1, expon);  /* K1 = 2^expon */
        K1--;         //addbigint(&K1, -1);   (K1 is used as a bitmask)
        ValCOdd = tmp1 & K1; // BigIntAnd(&tmp1, &K1, &ValCOdd);   // r = (b/2) - a*c mod 2^n
        NumberLength = K1.nbrLimbs;
        if (NumberLength > ValAOdd.nbrLimbs)
        {   /* extra leading limbs are set to zero */
            int lenBytes = (NumberLength - ValAOdd.nbrLimbs) * (int)sizeof(int);
            std::memset(&ValAOdd.limbs[ValAOdd.nbrLimbs], 0, lenBytes);
        }
        if (NumberLength > tmp2.nbrLimbs)
        {     /* extra leading limbs are set to zero */
            int lenBytes = (NumberLength - tmp2.nbrLimbs) * (int)sizeof(int);
            std::memset(&tmp2.limbs[tmp2.nbrLimbs], 0, lenBytes);
        }
        ComputeInversePower2(ValAOdd.limbs, tmp2.limbs, tmp1.limbs);
        ValCOdd *= tmp2;     // (void)BigIntMultiply(&ValCOdd, &tmp2, &ValCOdd);
        ValCOdd &= K1; // BigIntAnd(&ValCOdd, &K1, &ValCOdd);      // ((b/2) - a*c)/a mod 2^n
        ValCOdd *= tmp2;        // (void)BigIntMultiply(&ValCOdd, &tmp2, &ValCOdd);
        ValCOdd &= K1; // BigIntAnd(&ValCOdd, &K1, &ValCOdd);    // s = ((b/2) - a*c)/a^2 mod 2^n
        if (ValCOdd == 0)        //(BigIntIsZero(&ValCOdd))
        {         // s = 0, so its square root is also zero.
            sqrRoot = 0; // intToBigInteger(&sqrRoot, 0);
            expon -= expon / 2;
        }
        else  {
            DivideBigNbrByMaxPowerOf2(&bitsCZero, ValCOdd.limbs, &ValCOdd.nbrLimbs);
            // At this moment, bitsCZero = r and ValCOdd = q.
            if (((ValCOdd.limbs[0] & 7) != 1) || (bitsCZero & 1))  {
                return false;          // q != 1 or p2(r) == 0, so go out.
            }
            if (expon < 2) {         // Modulus is 2.
                sqrRoot = (bitsCZero > 0) ? 0 : 1; // intToBigInteger(&sqrRoot, (bitsCZero > 0) ? 0 : 1);
            }
            else { // Compute sqrRoot as the square root of ValCOdd.
                expon -= bitsCZero / 2;
                ComputeSquareRootModPowerOf2(expon, bitsCZero);
                expon--;
                if (expon == (bitsCZero / 2)) {
                    expon++;
                }
            }
        }

        // x = sqrRoot - b/2a.
        BigIntPowerOf2(&K1, expon);  /* K1 *= 2^expon  */
        K1--;          //addbigint(&K1, -1);
        NumberLength = K1.nbrLimbs;
        ComputeInversePower2(ValAOdd.limbs, tmp2.limbs, tmp1.limbs);
        setNbrLimbs(&tmp2);
        tmp1 = *pValB;  // CopyBigInt(&tmp1, pValB);
        BigIntDivideBy2(&tmp1);               // b/2
        tmp1 *= tmp2;       // (void)BigIntMultiply(&tmp1, &tmp2, &tmp1);  // b/2a
        BigIntNegate(tmp1);                  // -b/2a
        tmp1 &= K1; // BigIntAnd(&tmp1, &K1, &tmp1);         // -b/2a mod 2^expon
        tmp2 = tmp1 + sqrRoot;  // BigIntAdd(&tmp1, &sqrRoot, &tmp2);
        quad.Solution1[factorIndex] = tmp2 & K1; // BigIntAnd(&tmp2, &K1, &quad.Solution1[factorIndex]);
        tmp2 = tmp1 - sqrRoot;   // BigIntSubt(&tmp1, &sqrRoot, &tmp2);
        quad.Solution2[factorIndex] = tmp2 & K1; // BigIntAnd(&tmp2, &K1, &quad.Solution2[factorIndex]);
    }
    else if ((bitsAZero == 0) && (bitsBZero == 0))  {
        Quadr = *pValA *2;   // CopyBigInt(&Quadr, pValA);
        //BigIntMultiplyBy2(&Quadr);         // 2a
        Linear = *pValB;  // CopyBigInt(&Linear, pValB);        // b
        Const = *pValC;   // CopyBigInt(&Const, pValC);
        BigIntDivideBy2(&Const);           // c/2
        findQuadraticSolution(&quad.Solution1[factorIndex], expon - 1);
        quad.Solution1[factorIndex] *= 2; // BigIntMultiplyBy2(&quad.Solution1[factorIndex]);

        Quadr = *pValA *2;   // CopyBigInt(&Quadr, pValA);
        // BigIntMultiplyBy2(&Quadr);         // 2a
        Linear = Quadr + *pValB;  // BigIntAdd(&Quadr, pValB, &Linear); // 2a+b
        Const = *pValA;    // CopyBigInt(&Const, pValA);
        Const += *pValB;   // BigIntAdd(&Const, pValB, &Const);
        Const += *pValC;    //BigIntAdd(&Const, pValC, &Const);
        BigIntDivideBy2(&Const);           // (a+b+c)/2
        findQuadraticSolution(&quad.Solution2[factorIndex], expon - 1);
        quad.Solution2[factorIndex] *= 2; // BigIntMultiplyBy2(&quad.Solution2[factorIndex]);
        quad.Solution2[factorIndex]++; // addbigint(&quad.Solution2[factorIndex], 1);
    }
    else  {
        Quadr = *pValA;    // CopyBigInt(&Quadr, pValA);
        Linear = *pValB;   // CopyBigInt(&Linear, pValB);
        Const = *pValC;    // CopyBigInt(&Const, pValC);
        findQuadraticSolution(&quad.Solution1[factorIndex], expon);
        sol2Invalid = true;
    }
    BigIntPowerOf2(&Q, expon);         // Store increment. Q *= 2^expon
    return true;
}

/* uses global variable prime, Aux[3] to Aux[9], TestNbr, discriminant, 
result returned in sqrRoot */
static void ComputeSquareRootModPowerOfP(int nbrBitsSquareRoot)  {
    int correctBits;

    if (verbose > 1) {
        printf_s("ComputeSquareRootModPowerOfP: Aux[3] = ");
        PrintBigInteger(&Aux[3], 0);
        printf_s("\n prime = ");
        PrintBigInteger(&prime, 0);
        std::putchar('\n');
    }

    NumberLength = prime.nbrLimbs;
    int NumberLengthBytes = NumberLength * (int)sizeof(limb);
    std::memcpy(TestNbr, prime.limbs, NumberLengthBytes);  /* TestNbr = prime*/
    TestNbr[NumberLength] = 0;
    GetMontgomeryParms(NumberLength);
    if (verbose > 1) {
        printf_s("MontgomeryMultR1 = ");
        PrintBigInteger(&MontgomeryMultR1BI, 0);
        printf_s("\nMontgomeryMultR2 = ");
        PrintBigInteger(&MontgomeryMultR2BI, 0);
        std::putchar('\n');
    }
    Q = prime;    // CopyBigInt(&Q, &prime);
    if (verbose > 1) {
        printf_s("Q = ");
        PrintBigInteger(&Q, 0);
        printf_s("\n");
    }
    if ((prime.limbs[0] & 3) == 3)
    {                 // prime mod 4 = 3
        subtractdivide(Q, -1, 4);   // Q <- (prime+1)/4.
        BigIntModularPower(&Aux[3], &Q, &SqrtDisc);
    }
    else {
        limb* toConvert;
        // Convert discriminant to Montgomery notation.
        CompressLimbsBigInteger(Aux[5].limbs, &Aux[3]);
        modmult(Aux[5].limbs, MontgomeryMultR2, Aux[6].limbs);  // u
        if ((prime.limbs[0] & 7) == 5)
        {  // prime mod 8 = 5: use Atkin's method for modular square roots.
           // Step 1. v <- (2u)^((p-5)/8) mod p
           // Step 2. i <- (2uv^2) mod p
           // Step 3. square root of u <- uv (i-1)
          // Step 1.
            subtractdivide(Q, 5, 8);   // Q <- (prime-5)/8.
            AddBigNbrMod(Aux[6].limbs, Aux[6].limbs, Aux[7].limbs); // 2u
            modPow(Aux[7].limbs, Q.limbs, Q.nbrLimbs, Aux[8].limbs);
            // At this moment Aux[7].limbs is v in Montgomery notation.

            // Step 2.
            modmult(Aux[8].limbs, Aux[8].limbs, Aux[9].limbs);  // v^2
            modmult(Aux[7].limbs, Aux[9].limbs, Aux[9].limbs);  // i

            // Step 3.
            SubtBigNbrMod(Aux[9].limbs, MontgomeryMultR1, Aux[9].limbs); // i-1
            modmult(Aux[8].limbs, Aux[9].limbs, Aux[9].limbs);           // v*(i-1)
            modmult(Aux[6].limbs, Aux[9].limbs, Aux[9].limbs);           // u*v*(i-1)
            toConvert = Aux[9].limbs;
        }
        else
        { // prime = 1 (mod 8). Use Shanks' method for modular square roots.
          // Step 1. Select e >= 3, q odd such that p = 2^e * q + 1.
          // Step 2. Choose x at random in the range 1 < x < p such that jacobi (x, p) = -1.
          // Step 3. Set z <- x^q (mod p).
          // Step 4. Set y <- z, r <- e, x <- a^((q-1)/2) mod p, v <- ax mod p, w <- vx mod p.
          // Step 5. If w = 1, the computation ends with +/-v as the square root.
          // Step 6. Find the smallest value of k such that w^(2^k) = 1 (mod p)
          // Step 7. Set d <- y^(2^(r-k-1)) mod p, y <- d^2 mod p, r <- k, v <- dv mod p, w <- wy mod p.
          // Step 8. Go to step 5.
            int x;
            int e;
            int r;
            // Step 1.
            Q--;  // subtractdivide(Q, 1, 1);   // Q <- (prime-1).
            if (verbose > 1) {
                printf_s("Q = ");
                PrintBigInteger(&Q, 0);
                printf_s("\n");
            }
            DivideBigNbrByMaxPowerOf2(&e, Q.limbs, &Q.nbrLimbs);
            if (verbose > 1) {
                printf_s("step 1: Q = ");
                PrintBigInteger(&Q, 0);
                printf_s(" e = %d \n", e);
            }
            // Step 2.
            x = 1;
            do
            {
                x++;
                Aux[3] = x;  // intToBigInteger(&Aux[3], x);
            } while (BigIntJacobiSymbol(&Aux[3], &prime) >= 0);
            if (verbose > 1) {
                printf_s("step2: Aux[3] = %d \n", x);
            }

            // Step 3.
            // Get z <- x^q (mod p) in Montgomery notation.
            modPowBaseInt(x, Q.limbs, Q.nbrLimbs, Aux[4].limbs);  // z
            if (verbose > 1) {
                Aux[4].nbrLimbs = NumberLength;
                printf_s("step 3: Aux[4] = ");
                PrintBigInteger(&Aux[4], 0);
                std::putchar('\n');
            }
            // Step 4.
            NumberLengthBytes = NumberLength * (int)sizeof(limb);
            std::memcpy(Aux[5].limbs, Aux[4].limbs, NumberLengthBytes); // y
            r = e;
            K1 = Q;    // CopyBigInt(&K1, &Q);
            subtractdivide(K1, 1, 2);  /* K1 = (K-1)/2*/
            modPow(Aux[6].limbs, K1.limbs, K1.nbrLimbs, Aux[7].limbs); // x
            modmult(Aux[6].limbs, Aux[7].limbs, Aux[8].limbs);         // v
            modmult(Aux[8].limbs, Aux[7].limbs, Aux[9].limbs);         // w
            if (verbose > 1) {
                Aux[6].nbrLimbs = NumberLength;
                Aux[7].nbrLimbs = NumberLength;
                Aux[8].nbrLimbs = NumberLength;
                Aux[9].nbrLimbs = NumberLength;
                printf_s("step 4: Aux[5] = ");
                PrintBigInteger(&Aux[5], 0);
                printf_s("\nAux[6] = ");
                PrintBigInteger(&Aux[6], 0);
                printf_s("\nAux[7] = ");
                PrintBigInteger(&Aux[7], 0);
                printf_s("\nAux[8] = ");
                PrintBigInteger(&Aux[8], 0);
                printf_s("\nAux[9] = ");
                PrintBigInteger(&Aux[9], 0);
                std::putchar('\n');
            }
            // Step 5
            while (memcmp(Aux[9].limbs, MontgomeryMultR1, NumberLengthBytes) != 0)
            {
                // Step 6
                int k = 0;
                std::memcpy(Aux[10].limbs, Aux[9].limbs, NumberLengthBytes);
                do
                {
                    k++;
                    modmult(Aux[10].limbs, Aux[10].limbs, Aux[10].limbs);
                } while (memcmp(Aux[10].limbs, MontgomeryMultR1, NumberLengthBytes) != 0);
                // Step 7
                std::memcpy(Aux[11].limbs, Aux[5].limbs, NumberLengthBytes); // d
                for (int ctr = 0; ctr < (r - k - 1); ctr++)
                {
                    modmult(Aux[11].limbs, Aux[11].limbs, Aux[11].limbs);
                }
                modmult(Aux[11].limbs, Aux[11].limbs, Aux[5].limbs);   // y
                r = k;
                modmult(Aux[8].limbs, Aux[11].limbs, Aux[8].limbs);    // v
                modmult(Aux[9].limbs, Aux[5].limbs, Aux[9].limbs);     // w
            }
            toConvert = Aux[8].limbs;
        }
        // Convert from Montgomery to standard notation.
        NumberLengthBytes = NumberLength * (int)sizeof(limb);
        std::memset(Aux[4].limbs, 0, NumberLengthBytes); // Convert power to standard notation.
        Aux[4].limbs[0] = 1;
        modmult(Aux[4].limbs, toConvert, toConvert);
        UncompressLimbsBigInteger(toConvert, &SqrtDisc);
    }
    if (verbose > 1) {
        printf_s("SqrtDisc = ");
        PrintBigInteger(&SqrtDisc, 0);
        printf_s("\n");
    }
    // Obtain inverse of square root stored in SqrtDisc (mod prime).
    tmp2 = 1;  // intToBigInteger(&tmp2, 1);
    BigIntModularDivision(&tmp2, &SqrtDisc, &prime, &sqrRoot);
    correctBits = 1;
    Q = prime; // CopyBigInt(&Q, &prime);
    // Obtain nbrBitsSquareRoot correct digits of inverse square root.
    while (correctBits < nbrBitsSquareRoot)
    {   // Compute f(x) = invsqrt(x), f_{n+1}(x) = f_n * (3 - x*f_n^2)/2
        correctBits *= 2;
        Q = Q * Q;                  //(void)BigIntMultiply(&Q, &Q, &Q);           // Square Q.
        tmp1 = sqrRoot * sqrRoot;   //(void)BigIntMultiply(&sqrRoot, &sqrRoot, &tmp1);
        tmp2 = tmp1 % Q;            //(void)BigIntRemainder(&tmp1, &Q, &tmp2);
        tmp1 = tmp2 * discriminant; // (void)BigIntMultiply(&tmp2, &discriminant, &tmp1);
        tmp2 = tmp1 % Q;            // (void)BigIntRemainder(&tmp1, &Q, &tmp2);
        tmp1 = 3;                   // intToBigInteger(&tmp1, 3);
        tmp2 = tmp1 - tmp2;         // BigIntSubt(&tmp1, &tmp2, &tmp2);
        tmp1 = tmp2 * sqrRoot;      // (void)BigIntMultiply(&tmp2, &sqrRoot, &tmp1);
        if ((tmp1.limbs[0] & 1) != 0)  // if tmp1 is odd? 
        {
            tmp1 += Q;             // BigIntAdd(&tmp1, &Q, &tmp1);
        }
        BigIntDivide2(tmp1);
        sqrRoot = tmp1 % Q;        // (void)BigIntRemainder(&tmp1, &Q, &sqrRoot);
    }
    // Get square root of discriminant from its inverse by multiplying by discriminant.
    if (sqrRoot.sign == SIGN_NEGATIVE) {
        sqrRoot += Q;         // BigIntAdd(&sqrRoot, &Q, &sqrRoot);
    }
    sqrRoot *= discriminant;  // (void)BigIntMultiply(&sqrRoot, &discriminant, &sqrRoot);
    sqrRoot %= Q;             // (void)BigIntRemainder(&sqrRoot, &Q, &sqrRoot);
    if (verbose > 1) {
        printf_s("sqrRoot = ");
        PrintBigInteger(&sqrRoot, 0);
        printf_s("\n");
    }
}

// Solve Ax^2 + Bx + C = 0 (mod p^expon).
/* uses global variables prime, discriminant, V */
static bool SolveQuadraticEqModPowerOfP(int expon, int factorIndex,
    const BigInteger* pValA, const BigInteger* pValB)
{
    int correctBits;
    int nbrLimbs;
    int bitsAZero;
    int ctr;
    int deltaZeros;
    int NumberLengthBytes;

    if (verbose > 1) {
        printf_s("SolveQuadraticEqModPowerOfP: expon = %d \na = ", expon);
        PrintBigInteger(pValA, 0);
        printf_s("\nb = ");
        PrintBigInteger(pValB, 0);
        printf_s("\nprime = ");
        PrintBigInteger(&prime, 0);
        printf_s("\ndiscriminant = ");
        PrintBigInteger(&discriminant, 0);
        printf_s("\nV = ");
        PrintBigInteger(&V, 0);
        std::putchar ('\n');
    }

    // Number of bits of square root of discriminant to compute: expon + bits_a + 1,
    // where bits_a is the number of least significant bits of a set to zero.
    // To compute the square root, compute the inverse of sqrt, so only multiplications are used.
    // f(x) = invsqrt(x), f_{n+1}(x) = f_n * (3 - x*f_n^2)/2
    // Get maximum power of prime which divide ValA.
    ValAOdd = *pValA;  // CopyBigInt(&ValAOdd, pValA);
    bitsAZero = 0;
    for (;;)  {
        tmp1 = ValAOdd % prime;   // (void)BigIntRemainder(&ValAOdd, &prime, &tmp1);
        if (ValAOdd.sign == SIGN_NEGATIVE)   {
            ValAOdd += prime; // BigIntAdd(&ValAOdd, &prime, &ValAOdd);
        }
        if (tmp1 != 0)    // (!BigIntIsZero(&tmp1))
        {
            break;
        }
        ValAOdd /= prime;   // (void)BigIntDivide(&ValAOdd, &prime, &ValAOdd);
        bitsAZero++;
    }
    discriminant %= V;      // (void)BigIntRemainder(&discriminant, &V, &discriminant);
    // Get maximum power of prime which divides discriminant.
    if (discriminant == 0)  // (BigIntIsZero(&discriminant))
    {      // Discriminant is zero.
        deltaZeros = expon;
    }
    else {      // Discriminant is not zero.
        deltaZeros = 0;
        for (;;) {
            tmp1 = discriminant % prime;  // (void)BigIntRemainder(&discriminant, &prime, &tmp1);
            if (tmp1 != 0)   /* (!BigIntIsZero(&tmp1)) */ {
                break;
            }
            discriminant /= prime;  //(void)BigIntDivide(&discriminant, &prime, &discriminant);
            deltaZeros++;
        }
    }
    if (((deltaZeros & 1) != 0) && (deltaZeros < expon))
    {          // If delta is of type m*prime^n where m is not multiple of prime
               // and n is odd, there is no solution, so go out.
        return false;
    }
    deltaZeros >>= 1;
    // Compute inverse of -2*A (mod prime^(expon - deltaZeros)).
    ValAOdd = *pValA + *pValA; // BigIntAdd(pValA, pValA, &ValAOdd);
    (void)BigIntPowerIntExp(prime, expon - deltaZeros, tmp1); /* tmp1 = prime^(expon-deltaZeros) */
    if (verbose > 1) {
        printf_s("prime^(expon-deltaZeros) = ");
        PrintBigInteger(&tmp1, 0);
        std::putchar('\n');
    }
    ValAOdd %= tmp1;  //(void)BigIntRemainder(&ValAOdd, &tmp1, &ValAOdd);
    nbrLimbs = tmp1.nbrLimbs;
    if (ValAOdd.sign == SIGN_NEGATIVE) {
        ValAOdd.sign = SIGN_POSITIVE;           // Negate 2*A
    }
    else if (ValAOdd != 0)  // (!BigIntIsZero(&ValAOdd))
    {
        ValAOdd = tmp1 - ValAOdd;   // BigIntSubt(&tmp1, &ValAOdd, &ValAOdd);  // Negate 2*A
    }
    else {     // Nothing to do.
    }

    tmp2 = 1;        // intToBigInteger(&tmp2, 1);
    NumberLength = tmp1.nbrLimbs;
    NumberLengthBytes = NumberLength * (int)sizeof(limb);
    std::memcpy(TestNbr, tmp1.limbs, NumberLengthBytes);
    TestNbr[NumberLength] = 0;
    GetMontgomeryParms(NumberLength);
    BigIntModularDivision(&tmp2, &ValAOdd, &tmp1, &Aux[0]);
    ValAOdd = Aux[0];       // CopyBigInt(&ValAOdd, &Aux[0]);
    if (verbose > 1) {
        printf_s("ValAOdd = ");
        PrintBigInteger(&ValAOdd, 0);
        std::putchar('\n');
    }
    if (discriminant == 0)  // (BigIntIsZero(&discriminant))
    {     // Discriminant is zero.
        int lenBytes = nbrLimbs * (int)sizeof(limb);
        std::memset(sqrRoot.limbs, 0, lenBytes);
        sqrRoot.nbrLimbs = 1;            /* added 31/10/22 */
        sqrRoot.sign = SIGN_POSITIVE;    /* added 31/10/22 */
    }
    else
    {      // Discriminant is not zero.
      // Find number of digits of square root to compute.
        int nbrBitsSquareRoot = expon + bitsAZero - deltaZeros;
        (void)BigIntPowerIntExp(prime, nbrBitsSquareRoot, tmp1);
        nbrLimbs = tmp1.nbrLimbs;
        discriminant %= tmp1;    // (void)BigIntRemainder(&discriminant, &tmp1, &discriminant);
        if (discriminant.sign == SIGN_NEGATIVE) {
            discriminant += tmp1;  // BigIntAdd(&discriminant, &tmp1, &discriminant);
        }
        if (nbrLimbs > discriminant.nbrLimbs)
        {
            int lenBytes = (nbrLimbs - discriminant.nbrLimbs) * (int)sizeof(limb);
            std::memset(&discriminant.limbs[nbrLimbs], 0, lenBytes);
        }
        Aux[3] = discriminant % prime;     // (void)BigIntRemainder(&discriminant, &prime, &Aux[3]);
        if (Aux[3].sign == SIGN_NEGATIVE)  {
            Aux[3] += prime;  // BigIntAdd(&Aux[3], &prime, &Aux[3]);
        }
        if (BigIntJacobiSymbol(&Aux[3], &prime) != 1) {
            return false;         // Not a quadratic residue, so go out.
        }
        // Compute square root (mod prime) of discriminant in Aux[3]. result in sqrRoot
        ComputeSquareRootModPowerOfP(nbrBitsSquareRoot);
        //if (verbose > 1) {
        //    printf_s("sqrRoot = ");
        //    PrintBigInteger(&sqrRoot, 0);
        //    std::putchar('\n');
        //}
        // Multiply by square root of discriminant by prime^deltaZeros.
        for (ctr = 0; ctr < deltaZeros; ctr++)  {
            sqrRoot *= prime;    // (void)BigIntMultiply(&sqrRoot, &prime, &sqrRoot);
        }
    }
    correctBits = expon - deltaZeros;
    (void)BigIntPowerIntExp(prime, correctBits, Q);      // Store increment.
    // Compute x = (b + sqrt(discriminant)) / (-2a) and x = (b - sqrt(discriminant)) / (-2a)
    tmp1 = *pValB + sqrRoot;  // BigIntAdd(pValB, &sqrRoot, &tmp1);
    for (ctr = 0; ctr < bitsAZero; ctr++) {
        tmp2 = tmp1 % prime;   // (void)BigIntRemainder(&tmp1, &prime, &tmp2);
        if (tmp2 != 0)    // (!BigIntIsZero(&tmp2))
        {   // Cannot divide by prime, so go out.
            sol1Invalid = true;
            break;
        }
        tmp1 /= prime;     // (void)BigIntDivide(&tmp1, &prime, &tmp1);
    }
    tmp1 *= ValAOdd;    // (void)BigIntMultiply(&tmp1, &ValAOdd, &tmp1);
    quad.Solution1[factorIndex] = tmp1 % Q; // (void)BigIntRemainder(&tmp1, &Q, &quad.Solution1[factorIndex]);
    if (quad.Solution1[factorIndex].sign == SIGN_NEGATIVE)  {
        quad.Solution1[factorIndex] += Q; // BigIntAdd(&quad.Solution1[factorIndex], &Q, &quad.Solution1[factorIndex]);
    }

    tmp1 = *pValB - sqrRoot; // BigIntSubt(pValB, &sqrRoot, &tmp1);
    for (ctr = 0; ctr < bitsAZero; ctr++)  {
        tmp2 = tmp1 % prime;  // (void)BigIntRemainder(&tmp1, &prime, &tmp2);
        if (tmp2 != 0)    // (!BigIntIsZero(&tmp2))
        {   // Cannot divide by prime, so go out.
            sol2Invalid = true;
            break;
        }
        tmp1 /= prime;   // (void)BigIntDivide(&tmp1, &prime, &tmp1);
    }
    tmp1 *= ValAOdd;       // (void)BigIntMultiply(&tmp1, &ValAOdd, &tmp1);
    quad.Solution2[factorIndex] = tmp1 % Q;  // (void)BigIntRemainder(&tmp1, &Q, &quad.Solution2[factorIndex]);
    if (quad.Solution2[factorIndex].sign == SIGN_NEGATIVE) {
        quad.Solution2[factorIndex] += Q;  // BigIntAdd(&quad.Solution2[factorIndex], &Q, &quad.Solution2[factorIndex]);
    }
    return true;
}

/* Solve aX^2 + bX + C = 0 (mod p^e)
uses global variables prime, Q, V, quad.solution1, quad.solution2 ... */
static void QuadraticTermMultipleOfP(int expon, int factorIndex,
    const BigInteger* pValA, const BigInteger* pValB, const BigInteger* pValC)
{
    // Perform Newton approximation:
    // x_{n+1} = x_n - (a*x_n^2 + b*x_n + c) / (2*a_x + b)      
    BigInteger* ptrSolution = &quad.Solution1[factorIndex];
    int NumberLengthBytes;

    if (verbose > 1) {
        printf_s("QuadraticTermMultipleOfP: expon = %d \nA = ", expon);
        PrintBigInteger(pValA, 0);
        printf_s("\nB = ");
        PrintBigInteger(pValB, 0);
        printf_s("\nC = ");
        PrintBigInteger(pValC, 0);
        printf_s("\nprime = ");
        PrintBigInteger(&prime, 0);
        std::putchar('\n');
    }

    NumberLength = prime.nbrLimbs;
    NumberLengthBytes = NumberLength * (int)sizeof(limb);
    std::memcpy(TestNbr, prime.limbs, NumberLengthBytes);  /* TestNbr = prime */
    TestNbr[NumberLength] = 0;
    GetMontgomeryParms(NumberLength);
    BigIntModularDivision(pValC, pValB, &prime, ptrSolution); /* Solution = C/B (mod prime) */
    BigIntNegate(*ptrSolution);  
    if (ptrSolution->sign == SIGN_NEGATIVE)
    {
        *ptrSolution += prime; // BigIntAdd(ptrSolution, &prime, ptrSolution);
    }
    for (int currentExpon = 2; currentExpon < (2 * expon); currentExpon *= 2)
    {
        (void)BigIntPowerIntExp(prime, currentExpon, V);   /* V = prime^currentExpon */
        Q = *pValA * *ptrSolution;   // (void)BigIntMultiply(pValA, ptrSolution, &Q);// a*x_n
        L = Q;                       // CopyBigInt(&L, &Q);  
        Q += *pValB;                 // BigIntAdd(&Q, pValB, &Q);    // a*x_n + b
        Q %= V;                      //(void)BigIntRemainder(&Q, &V, &Q);
        Q *= *ptrSolution;           // (void)BigIntMultiply(&Q, ptrSolution, &Q);   // a*x_n^2 + b*x_n
        Q += *pValC;                 // BigIntAdd(&Q, pValC, &Q);   // a*x_n^2 + b*x_n + c
        Q %= V;                      //(void)BigIntRemainder(&Q, &V, &Q);    // Numerator. 
        L *= 2;                     // multint(&L, &L, 2);                   // 2*a*x_n
        L += *pValB;                // BigIntAdd(&L, pValB, &L);      // 2*a*x_n + b
        L %= V;                     // (void)BigIntRemainder(&L, &V, &L);    // Denominator
        NumberLength = V.nbrLimbs;
        NumberLengthBytes = NumberLength * (int)sizeof(limb);
        std::memcpy(TestNbr, V.limbs, NumberLengthBytes);
        TestNbr[NumberLength] = 0;
        GetMontgomeryParms(NumberLength);
        BigIntModularDivision(&Q, &L, &V, &Aux1);    /* Aux1 = Q/L (mod V) */
        *ptrSolution -= Aux1;   // BigIntSubt(ptrSolution, &Aux1, ptrSolution);
        *ptrSolution %= V;      // (void)BigIntRemainder(ptrSolution, &V, ptrSolution);
        if (ptrSolution->sign == SIGN_NEGATIVE)
        {
            *ptrSolution += V;  //BigIntAdd(ptrSolution, &V, ptrSolution);
        }
    }

    (void)BigIntPowerIntExp(prime, expon, Q);  /* Q = prime^expon */
    *ptrSolution %= Q;           //(void)BigIntRemainder(ptrSolution, &Q, ptrSolution);
    /* make Solution 2 = Solution 1*/
    quad.Solution2[factorIndex] = quad.Solution1[factorIndex];  // CopyBigInt(&quad.Solution2[factorIndex], &quad.Solution1[factorIndex]);
    if (verbose > 1) {
        printf_s("Solution = ");
        PrintBigInteger(ptrSolution, 0);
        std::putchar('\n');
    }
}

/* uses global variables sol1Invalid, sol2Invalid, Aux[0], discriminant, prime,
quad.solution1, quad.solution2, */
static bool QuadraticTermNotMultipleOfP(int expon, int factorIndex,
    const BigInteger* pValA, const BigInteger* pValB, const BigInteger* pValC)
{
    bool solutions;

    if (verbose > 1) {
        printf_s("QuadraticTermNotMultipleOfP: expon = %d, factorIndex = %d \n a = ", expon, factorIndex);
        PrintBigInteger(pValA, 0);
        printf_s("\n b = ");
        PrintBigInteger(pValB, 0);
        printf_s("\n c = ");
        PrintBigInteger(pValC, 0);
    }

    sol1Invalid = false;
    sol2Invalid = false;
    // Compute discriminant = ValB^2 - 4*ValA*ValC.
    Aux[0] = *pValB * *pValB;        // (void)BigIntMultiply(pValB, pValB, &Aux[0]);          /* Aux[0] = b^2  */
    discriminant = *pValA * *pValC;  // (void)BigIntMultiply(pValA, pValC, &discriminant);    /* discriminant =a.c */
    discriminant *= 4;               // multint(&discriminant, &discriminant, 4);             /* discriminant = 4 * a.c */
    discriminant = Aux[0] - discriminant;  // BigIntSubt(&Aux[0], &discriminant, &discriminant);  /*  discriminant = b^2 - 4*a.c */

    if (verbose > 1) {
        printf_s("\n discriminant = ");
        PrintBigInteger(&discriminant, 0);
        std::putchar('\n');
    }

    if (prime == 2)     //((prime.nbrLimbs == 1) && (prime.limbs[0] == 2))
    {         /* Prime p is 2 */
        solutions = SolveQuadraticEqModPowerOf2(expon, factorIndex,
            pValA, pValB, pValC, prime);
    }
    else  {                // Prime is not 2
        solutions = SolveQuadraticEqModPowerOfP(expon, factorIndex, pValA, pValB);
    }
    if (!solutions || (sol1Invalid && sol2Invalid))
    {     // Both solutions are invalid. Go out.
        if (ShowNoSolsModPrime != NULL)
        {
            ShowNoSolsModPrime(expon);
        }
        return false;
    }
    if (sol1Invalid)
    {     // quad.Solution1 is invalid. Overwrite it with quad.Solution2.
        quad.Solution1[factorIndex] = quad.Solution2[factorIndex]; // (&quad.Solution1[factorIndex], &quad.Solution2[factorIndex]);
    }
    else if (sol2Invalid)
    {     // quad.Solution2 is invalid. Overwrite it with quad.Solution1.
        quad.Solution2[factorIndex] = quad.Solution1[factorIndex]; // CopyBigInt(&quad.Solution2[factorIndex], &quad.Solution1[factorIndex]);
    }
    else {              // Nothing to do.
    }

    if (verbose > 1) {
        std::cout << "solutions = " << solutions << " sol1Invalid = " << sol1Invalid
            << " sol2Invalid = " << sol2Invalid << '\n';
        if (!sol1Invalid) {
            printf_s("solution 1 = ");
            PrintBigInteger(&quad.Solution1[factorIndex], 0);
            std::putchar('\n');
        }
        if (!sol2Invalid) {
            printf_s("solution 2 = ");
            PrintBigInteger(&quad.Solution2[factorIndex], 0);
            std::putchar('\n');
        }
    }

    // BigIntSubt(&quad.Solution2[factorIndex], &quad.Solution1[factorIndex], &Aux[0]);
    // if (Aux[0].sign == SIGN_NEGATIVE)
    if (quad.Solution2[factorIndex] < quad.Solution1[factorIndex])
    {     // quad.Solution2 is less than quad.Solution1, so exchange them.
        Aux[0] = quad.Solution1[factorIndex]; // CopyBigInt(&Aux[0], &common.quad.Solution1[factorIndex]);
        quad.Solution1[factorIndex] = quad.Solution2[factorIndex]; // CopyBigInt(&quad.Solution1[factorIndex], &quad.Solution2[factorIndex]);
        quad.Solution2[factorIndex] = Aux[0];  // CopyBigInt(&quad.Solution2[factorIndex], &Aux[0]);
    }
    if (ShowSolutionsModPrime != NULL) {
        ShowSolutionsModPrime(factorIndex, expon, &Q);
    }
    return true;
}

void SetCallbacksForSolveEquation(pSolution solutionCback
    //, pShowSolutionsModPrime showSolutionsModPrime, pShowNoSolsModPrime showNoSolsModPrime
)
{
    SolutionP = solutionCback;
    //ShowSolutionsModPrime = showSolutionsModPrime;
    //ShowNoSolsModPrime = showNoSolsModPrime;
}

// Solve Ax² + Bx + C = 0 (mod N).
void SolveEquation(const BigInteger* pValA, const BigInteger* pValB,
    const BigInteger* pValC, const BigInteger* pValN,
    BigInteger* pGcdAllParm, const BigInteger* pValNnParm)
{
    const struct sFactors* pstFactor;
    pGcdAll = pGcdAllParm;
    pValNn = pValNnParm;  /* pValNn is used in PerformChineseRemainderTheorem */

    if (verbose > 1) {
        printf_s("SolveEquation: a = ");
        PrintBigInteger(pValA, 0);
        printf_s(" b = ");
        PrintBigInteger(pValB, 0);
        printf_s(" c = ");
        PrintBigInteger(pValC, 0);
        printf_s(" n = ");
        PrintBigInteger(pValN, 0);
        printf_s("\nGcdAllParm = ");
        PrintBigInteger(pGcdAllParm, 0);
        printf_s(" ValNnParm = ");
        PrintBigInteger(pValNnParm, 0);
        std::putchar('\n');
    }

    Aux[0] = *pValA % *pValN;   //(void)BigIntRemainder(pValA, pValN, &Aux[0]);
    if (Aux[0] == 0)            //(BigIntIsZero(&Aux[0]))
    {           // Linear equation.
        SolveModularLinearEquation(pValA, pValB, pValC, pValN);
        return;
    }
    if (LastModulus != *pValN)  //(!BigIntEqual(&LastModulus, pValN))
    {     // Last modulus is different from ValN; need to factorise N

        LastModulus = *pValN;  //CopyBigInt(&LastModulus, pValN);
        NumberLength = pValN->nbrLimbs;
        //BigInteger2IntArray(nbrToFactor, pValN);  /* copy N to nbrToFactor */
        //factor(pValN, nbrToFactor, factorsMod, astFactorsMod);  /* get factors of N */
        factor(pValN, factorsMod, astFactorsMod);
        /* the same factor list is returned in both factorsMod and astFactorsMod,
           in different formats */

        if (verbose > 1)
            printFactors(astFactorsMod);  /* print factors */

    }
    Q = 0;  // intToBigInteger(&Q, 0);  /* set Q to 0 */
    nbrFactors = astFactorsMod[0].multiplicity;
    if (nbrFactors > 400) {
        /* array index would exceed maximum */
        copyStr(&ptrOutput, "Modulus too large (>400 unique factors) \n");
        SolNbr = -1;  /* small frig so that we get the ouput text we want */
        return;
    }
    pstFactor = &astFactorsMod[1];
    for (int factorIndex = 0; factorIndex < nbrFactors; factorIndex++) {
        int expon = pstFactor->multiplicity;
        if (expon == 0) {
            quad.Solution1[factorIndex] = 0; // intToBigInteger(&quad.Solution1[factorIndex], 0);
            quad.Solution2[factorIndex] = 0; // intToBigInteger(&quad.Solution2[factorIndex], 0);
            quad.Increment[factorIndex] = 1; // intToBigInteger(&quad.Increment[factorIndex], 1);
            pstFactor++;    /* move pointer to next factor*/
            continue;
        }
        NumberLength = *pstFactor->ptrFactor;
        IntArray2BigInteger(pstFactor->ptrFactor, &prime);  /* convert factor to BigInteger */
        (void)BigIntPowerIntExp(prime, expon, V);  /* V = prime Factor^exponent */
        L = *pValA % prime;    // (void)BigIntRemainder(pValA, &prime, &L);
        //if (BigIntIsZero(&L) && !((prime.nbrLimbs == 1) && (prime.limbs[0] == 2)))
        if (L == 0 && prime != 2)
        {     // ValA multiple of prime, means linear equation mod prime. Also prime is not 2.
            //if ((BigIntIsZero(pValB)) && !(BigIntIsZero(pValC)))
            if (*pValB == 0 && *pValC != 0) {
                return;  // There are no solutions: B=0 and C!=0
            }
            QuadraticTermMultipleOfP(expon, factorIndex, pValA, pValB, pValC);
        }
        else {           /* If quadratic equation mod p */
            if (!QuadraticTermNotMultipleOfP(expon, factorIndex, pValA, pValB, pValC))
            {
                return; /* no solutions */
            }
        }

        if (verbose > 1) {
            printf_s("Q = ");
            PrintBigInteger(&Q, 0);
            std::putchar('\n');
        }

        quad.Increment[factorIndex] = Q; //CopyBigInt(&quad.Increment[factorIndex], &Q);  /* Increment = Q */
        Exponents[factorIndex] = 0;
        pstFactor++;    /* move pointer to next factor */
    }   // end for
    PerformChineseRemainderTheorem(pValA, pValB, pValC, pValN);
}
