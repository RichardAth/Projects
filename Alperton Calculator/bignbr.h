#pragma once
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
#ifndef _BIGNBR_H
#define _BIGNBR_H
#define MAX_LEN 2500        // 20000 digits
#define BITS_PER_GROUP 31
#define BITS_PER_INT_GROUP 31
#define HALF_INT_RANGE (1 << (BITS_PER_INT_GROUP - 1))
#define MAX_INT_NBR ((int)((1U << BITS_PER_INT_GROUP)-1))
#define LIMB_RANGE (1U<<BITS_PER_GROUP)
#define SMALL_NUMBER_BOUND 32768
#define _USING64BITS_ 1
struct mylimb
{
	int x;      // maximum value in each limb is MAX_VALUE_LIMB
};
#define MAX_VALUE_LIMB (LIMB_RANGE-1)
typedef struct mylimb limb;
typedef enum eBoolean
{
	FALSE = 0,
	TRUE,
} boolean;

enum eSign
{
	SIGN_POSITIVE = 0,
	SIGN_NEGATIVE,
};

typedef struct BigInteger
{
	limb limbs[MAX_LEN];
	int nbrLimbs;
	enum eSign sign;
} BigInteger;

extern limb TestNbr[MAX_LEN];
extern limb MontgomeryMultR2[MAX_LEN];
extern limb MontgomeryMultR1[MAX_LEN];
extern int NumberLength, NumberLengthR1;
extern int groupLen;

//#include "expression.h"
void multiply(const limb *factor1, const limb *factor2, limb *result, int len, int *pResultLen);
void squareRoot(limb *argument, limb *sqRoot, int len, int *pLenSqRoot);
void Dec2Bin(const char *decimal, limb *binary, int digits, int *bitGroups);
void Bin2Dec(const limb *binary, char *decimal, int nbrLimbs, int groupLen);
void Bin2Hex(const limb *binary, char *decimal, int nbrLimbs, int groupLen);
void int2dec(char **pOutput, long long nbr);
void int2hex(char **pOutput, int nbr);
void GetMontgomeryParms(int len);
void AddBigNbrModN(limb *Nbr1, limb *Nbr2, limb *Sum, limb *TestNbr, int NumberLength);
void SubtBigNbrModN(limb *Nbr1, limb *Nbr2, limb *Sum, limb *TestNbr, int NumberLength);
void SubtBigNbrMod(limb *Nbr1, limb *Nbr2, limb *Sum);
void modmult(const limb *factor1, const limb *factor2, limb *product);
void modmultInt(limb *factorBig, int factorInt, limb *result);
void modmultIntExtended(limb *factorBig, int factorInt, limb *result, const limb *pTestNbr, int nbrLen);
void AddBigNbrMod(limb *Nbr1, limb *Nbr2, limb *Sum);
void modPowBaseInt(int base, const limb *exp, int nbrGroupsExp, limb *power);
void modPow(const limb *base, const limb *exp, int nbrGroupsExp, limb *power);
void modPowLimb(const limb *base, const limb *exp, limb *power);
void AdjustModN(limb *Nbr, limb *TestNbr, int NumberLength);
int fsquares(void);
void AddBigInt(limb *pAddend1, limb *pAddend2, limb *pSum, int nbrLimbs);
void SubtractBigInt(limb *pAddend1, limb *pAddend2, limb *pSum, int nbrLimbs);
void BigIntAdd (const BigInteger *pAddend1, const BigInteger *pAddend2, BigInteger *pSum);
void BigIntSubt(const BigInteger *pAddend1, const BigInteger *pAddend2, BigInteger *pSum);
void BigIntNegate(const BigInteger *pSrc, BigInteger *pDest);
enum eExprErr BigIntDivide(const BigInteger *pDividend, const BigInteger *pDivisor, 
	BigInteger *pQuotient);
enum eExprErr BigIntMultiply(const BigInteger *pFactor1, const BigInteger *pFactor2, BigInteger *pProduct);
enum eExprErr BigIntRemainder(const BigInteger *pDividend, const BigInteger *pDivisor, 
	BigInteger *pRemainder);
enum eExprErr BigIntPower(BigInteger *pBase, BigInteger *pExponent, BigInteger *pPower);
enum eExprErr BigIntPowerIntExp(const BigInteger *pBase, int expon, BigInteger *pPower);
void BigInteger2Dec(const BigInteger *pBigInt, char *decimal, int groupLen);
void BigInteger2Hex(const BigInteger *pBigInt, char *decimal, int groupLen);
void BigIntGcd(const BigInteger *pArg1, const BigInteger *pArg2, BigInteger *pResult);
void BigIntGeneralModularDivision(BigInteger *Num, BigInteger *Den, BigInteger *mod, BigInteger *quotient);
void BigIntModularDivision(BigInteger *Num, BigInteger *Den, BigInteger *mod, BigInteger *quotient);
void BigIntModularDivisionPower2(BigInteger *Num, BigInteger *Den, BigInteger *mod, BigInteger *quotient);
void BigIntModularDivisionSaveTestNbr(BigInteger *Num, BigInteger *Den, BigInteger *mod, BigInteger *quotient);
void multint(BigInteger *pResult, const BigInteger *pMult, int iMult);
void multadd(BigInteger *pResult, int iMult, const BigInteger *pMult, int addend);
void addmult(BigInteger *pResult, const BigInteger *pMult1, int iMult1, const BigInteger *pMult2, int iMult2);
void BigIntPowerOf2(BigInteger *pResult, int expon);
int getRemainder(const BigInteger *pBigInt, int divisor);
void subtractdivide(BigInteger *pBigInt, int subt, int divisor);
void addbigint(BigInteger *pResult, int addend);
boolean TestBigNbrEqual(const BigInteger *pNbr1, const BigInteger *pNbr2);
void CopyBigInt(BigInteger *pDest, const BigInteger *pSrc);
void ModInvBigNbr(limb *num, limb *inv, limb *mod, int NumberLength);
int modInv(int NbrMod, int currentPrime);
int getNbrLimbs(const limb *bigNbr);
void BigIntDivide2(BigInteger *pArg);
int PowerCheck(const BigInteger *pBigNbr, BigInteger *pBase);
int BpswPrimalityTest(const BigInteger *pBigNbr);
void UncompressBigInteger(/*@in@*/const int *ptrValues, /*@out@*/BigInteger *bigint);
void CompressBigInteger(/*@out@*/int *ptrValues, /*@in@*/const BigInteger *bigint);
void UncompressLimbsBigInteger(/*@in@*/const limb *ptrValues, /*@out@*/BigInteger *bigint);
void CompressLimbsBigInteger(/*@out@*/limb *ptrValues, /*@in@*/const BigInteger *bigint);
void NbrToLimbs(int nbr, /*@out@*/limb *limbs, int len);
void ComputeInversePower2(/*@in@*/const limb *value, /*@out@*/limb *result, /*@out@*/limb *aux);
int BigNbrIsZero(const limb *value);
void intToBigInteger(BigInteger *bigint, int value);
void longToBigInteger(BigInteger *bigint, long long value);
void expBigNbr(BigInteger *pBigNbr, double logar);
double logBigNbr(const BigInteger *pBigNbr);
double logLimbs(const limb *pBigNbr, int nbrLimbs);
double getMantissa(const limb *ptrLimb, int nbrLimbs);
void UncompressIntLimbs(/*@in@*/const int *ptrValues, /*@out@*/limb *bigint, int nbrLen);
void CompressIntLimbs(/*@out@*/int *ptrValues, /*@in@*/const limb *bigint, int nbrLen);
int checkOne(const limb *value, int nbrLimbs);
int checkMinusOne(const limb *value, int nbrLimbs);
void DivideBigNbrByMaxPowerOf2(int *pShRight, limb *number, int *pNbrLimbs);
void BigIntModularPower(const BigInteger *base, BigInteger *exponent, BigInteger *power);
enum eExprErr BigIntGeneralModularPower(BigInteger *base, BigInteger *exponent, BigInteger *mod, BigInteger *power);

void ChSignBigNbr      (int *nbr, int length);
void ChSignBigNbrB     (int *nbr, int length);
void AddBigNbr         (const int *pNbr1, const int *pNbr2, int *pSum, int nbrLen);
void SubtractBigNbr    (const int *pNbr1, const int *pNbr2, int *pDiff, int nbrLen);
void AddBigNbrB        (const int *pNbr1, const int *pNbr2, int *pSum, int nbrLen);
void SubtractBigNbrB   (const int *pNbr1, const int *pNbr2, int *pDiff, int nbrLen);
void AddBigIntModN     (const int *pNbr1, const int *pNbr2, int *pSum, const int *pMod, int nbrLen);
void SubtractBigNbrModN(const int *pNbr1, const int *pNbr2, int *pDiff, const int *pMod, int nbrLen);
void MultBigNbrByInt   (const int *bigFactor, int factor, int *bigProduct, int nbrLen);
void MultBigNbrByIntB  (const int *bigFactor, int factor, int *bigProduct, int nbrLen);
void DivBigNbrByInt    (const int *pDividend, int divisor, int *pQuotient, int nbrLen);
int RemDivBigNbrByInt  (const int *pDividend, int divisor, int nbrLen);
void MultBigNbr        (const int *pFactor1, const int *pFactor2, int *pProd, int nbrLen);
void IntToBigNbr       (int value, int *bigNbr, int nbrLength);
int BigNbrToBigInt     (const BigInteger *pBigNbr, int *pBigInt);
void BigIntToBigNbr    (BigInteger *pBigNbr, const int *pBigInt, int nbrLenBigInt);
//int BigNbrIsZero(const limb *pNbr);
void GcdBigNbr         (const int *pNbr1, const int *pNbr2, int *pGcd, int nbrLen);
void AdjustBigIntModN  (int *Nbr, int *Mod, int nbrLen);
void MultBigNbrModN    (int *Nbr1, int *Nbr2, int *Prod, int *Mod, int nbrLen);
void MultBigNbrByIntModN(int *Nbr1, int Nbr2, int *Prod, int *Mod, int nbrLen);
int intDoubleModPow(int NbrMod, int Expon, int currentPrime);
void ModInvBigInt(int *num, int *inv, int *mod, int NumberLength);
void IntToBigNbr(int value, int *bigNbr, int nbrLength);
int JacobiSymbol(int upper, int lower);
int BigIntJacobiSymbol(BigInteger *upper, BigInteger *lower);
void DivideBigNbrByMaxPowerOf4(int *pPower4, limb *value, int *pNbrLimbs);


typedef void(*mmCback)(void);
#endif
