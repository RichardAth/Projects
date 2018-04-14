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
#include <vector>
#define MAX_LEN 2500        // approximately 20000 digits
#define BITS_PER_GROUP 31
#define BITS_PER_INT_GROUP 31
#define HALF_INT_RANGE (1 << (BITS_PER_INT_GROUP - 1))
#define MAX_INT_NBR ((int)((1U << BITS_PER_INT_GROUP)-1))
#define LIMB_RANGE (1U<<BITS_PER_GROUP)
#define SMALL_NUMBER_BOUND 32768
#define _USING64BITS_ 1
#define Max_Factors 2500  // number of unique factors
struct limb
{
	int x=0;      // maximum value in each limb is MAX_VALUE_LIMB
};
#define MAX_VALUE_LIMB (LIMB_RANGE-1)

enum eSign
{
	SIGN_POSITIVE = 0,
	SIGN_NEGATIVE,
};


class BigInteger {
public:
	limb limbs[MAX_LEN];
	int nbrLimbs = 1;
	enum eSign sign = SIGN_POSITIVE;
	/* overload assignment operator here */
	BigInteger & operator = (const BigInteger &other) {
		if (&other == this)
			return *this;		// if lhs == rhs do nothing
		nbrLimbs = other.nbrLimbs;
		sign = other.sign;
		memcpy(limbs, other.limbs, nbrLimbs * sizeof(int));
		while (nbrLimbs > 1 && limbs[nbrLimbs - 1].x == 0) {
			nbrLimbs--;  // remove any leading zeros
		}
		return *this;
	}
	BigInteger & operator = (const int value) {
		if (value >= 0) {
			limbs[0].x = value;
			sign = SIGN_POSITIVE;
		}
		else {
			limbs[0].x = -value;
			sign = SIGN_NEGATIVE;
		}
		nbrLimbs = 1;
		return *this;
	}
	BigInteger & operator = (long long value) {
		int noOfLimbs = 0;
		sign = SIGN_POSITIVE;
		if (value < 0) {
			sign = SIGN_NEGATIVE;
			value = -value;
		}

		do {
			limbs[noOfLimbs++].x = (int)value & MAX_VALUE_LIMB;
			value >>= BITS_PER_GROUP;
		} while (value != 0);

		nbrLimbs = noOfLimbs;
		return *this;
	}
	/* other operators are overloaded later */
} ;

extern limb TestNbr[MAX_LEN];
extern limb MontgomeryMultR2[MAX_LEN];
extern limb MontgomeryMultR1[MAX_LEN];
extern int NumberLength, NumberLengthR1;
extern int groupLen;
bool BigNbrIsZero(const limb *value, int nbrLength);
void multiply(const limb *factor1, const limb *factor2, limb *result, int len, int *pResultLen);
void squareRoot(limb *argument, limb *sqRoot, int len, int *pLenSqRoot);
void int2dec(char **pOutput, long long nbr);
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
void AdjustModN(limb *Nbr, const limb *TestNbr, int NumberLength);
//void BigIntAdd(const BigInteger &pAddend1, const BigInteger &pAddend2, BigInteger &pSum);
BigInteger &BigIntAdd(const BigInteger &Addend1, const BigInteger &Addend2);
static BigInteger &operator +(const BigInteger &a, const BigInteger &b) {
	return BigIntAdd(a, b);
}
static BigInteger &operator +=(BigInteger &a, const BigInteger &b) {
	a = BigIntAdd(a, b);
	return a;
}
//void BigIntSubt(const BigInteger &pAddend1, const BigInteger &pAddend2, BigInteger &pSum);
BigInteger &BigIntSubt(const BigInteger &Minuend, const BigInteger &Subtrahend);
static BigInteger &operator -(const BigInteger &a, const BigInteger &b) {
	return BigIntSubt(a, b);
}
static BigInteger &operator -=(BigInteger &a, const BigInteger &b) {
	a =  BigIntSubt(a, b);
	return a;
}

void BigIntNegate(BigInteger &pDest);
//void BigIntDivide(const BigInteger &pDividend, const BigInteger &pDivisor,
//	BigInteger &pQuotient);
BigInteger BigIntDivide(const BigInteger &Dividend, const BigInteger &Divisor);
static BigInteger operator / (const BigInteger &a, const BigInteger &b) {
	return BigIntDivide(a, b);
}

//void BigIntMultiply(const BigInteger &pFactor1, const BigInteger &pFactor2, BigInteger &pProduct);
BigInteger BigIntMultiply(const BigInteger &Factor1, const BigInteger &Factor2);
static BigInteger operator * (const BigInteger &a, const BigInteger &b) {
	return BigIntMultiply(a, b);
}
static BigInteger operator *= (BigInteger &a, const BigInteger &b) {
	return a = BigIntMultiply(a, b);
}

//void BigIntRemainder(const BigInteger &Dividend, const BigInteger &Divisor,
//	BigInteger &Remainder);
BigInteger BigIntRemainder(const BigInteger &Dividend, const BigInteger &Divisor);
static BigInteger operator %(const BigInteger &Dividend, const BigInteger &Divisor) {
	return BigIntRemainder(Dividend, Divisor);
}

void BigIntPowerIntExp(const BigInteger &pBase, int expon, BigInteger &pPower);
void BigIntGcd(const BigInteger &pArg1, const BigInteger &pArg2, BigInteger &pResult);
void BigIntModularDivision(BigInteger &Num, BigInteger &Den, BigInteger &mod, BigInteger &quotient);

int getRemainder(const BigInteger &pBigInt, int divisor);
static int operator %(const BigInteger &a, const int divisor) {
	return getRemainder(a, divisor);
}
void subtractdivide(BigInteger &pBigInt, int subt, int divisor);
void addbigint(BigInteger &pResult, int addend);
static BigInteger operator +=(BigInteger &a, const int b) {
	addbigint(a, b);
	return a;
}
static BigInteger operator ++(BigInteger &a, int x) {
	addbigint(a, 1);
	return a;
}
static BigInteger operator --(BigInteger &a, int x) {
	addbigint(a, -1);
	return a;
}
static BigInteger operator -=(BigInteger &a, const int b) {
	addbigint(a, -b);
	return a;
}
bool TestBigNbrEqual(const BigInteger &Nbr1, const BigInteger &Nbr2);
bool TestBigNbrLess(const BigInteger &Nbr1, const BigInteger &Nbr2);
/* overload comparison operators here */
static bool operator ==(const BigInteger &a, const BigInteger &b) {
	return TestBigNbrEqual(a, b);
}
static bool operator ==(const BigInteger &a, const long long b) {
	if (b == 0)
		return BigNbrIsZero(a.limbs, a.nbrLimbs);
	BigInteger Btemp;
	Btemp = b;
	return TestBigNbrEqual(a, Btemp);
}

static bool operator !=(const BigInteger &a, const BigInteger &b) {
	return !TestBigNbrEqual(a, b);
}
static bool operator !=(const BigInteger &a, const long long b) {
	BigInteger Btemp;
	Btemp = b;
	return !TestBigNbrEqual(a, Btemp);
}

static bool operator <(const BigInteger &a, const BigInteger &b) {
	return TestBigNbrLess(a, b);
}
static bool operator <(const BigInteger &a, const long long b) {
	if (b == 0) {
		return (a.sign == SIGN_NEGATIVE);
	}
	BigInteger Btemp;
	Btemp = b;
	return TestBigNbrLess(a, Btemp);
}

static bool operator >(const BigInteger &a, const BigInteger &b) {
	return TestBigNbrLess(b, a);
}
static bool operator >(const BigInteger &a, const long long b) {
	BigInteger Btemp;
	Btemp = b;
	return TestBigNbrLess(Btemp, a);
}

static bool operator <=(const BigInteger &a, const BigInteger &b) {
	return !TestBigNbrLess(b, a);
}
static bool operator <=(const BigInteger &a, const long long b) {
	BigInteger Btemp;
	Btemp = b;
	return !TestBigNbrLess(Btemp, a);
}
static bool operator >=(const BigInteger &a, const BigInteger &b) {
	return !TestBigNbrLess(a, b);
}
static bool operator >=(const BigInteger &a, const long long b) {
	if (b == 0)
		return (a.sign == SIGN_POSITIVE);
	BigInteger Btemp;
	Btemp = b;
	return !TestBigNbrLess(a, Btemp);
}

//BigInteger &CopyBigInt(BigInteger &pDest, const BigInteger &pSrc);
void ModInvBigNbr(limb *num, limb *inv, limb *mod, int NumberLength);
//int modInv(int NbrMod, int currentPrime);
void BigIntDivide2(BigInteger &pArg);
int PowerCheck(const BigInteger &pBigNbr, BigInteger &pBase);
int BpswPrimalityTest(const BigInteger *pBigNbr);
void UncompressBigInteger(/*@in@*/const int *ptrValues, /*@out@*/BigInteger &bigint);
void CompressBigInteger(/*@out@*/int *ptrValues, /*@in@*/const BigInteger &bigint);
void UncompressLimbsBigInteger(/*@in@*/const limb *ptrValues, /*@out@*/BigInteger &bigint);
void CompressLimbsBigInteger(/*@out@*/limb *ptrValues, /*@in@*/const BigInteger &bigint);
void ComputeInversePower2(/*@in@*/const limb *value, /*@out@*/limb *result, /*@out@*/limb *aux);

//void intToBigInteger(BigInteger &bigint, const int value);
//void longToBigInteger(BigInteger &bigint, long long value);
void expBigNbr(BigInteger &pBigNbr, double logar);
double logBigNbr(const BigInteger &pBigNbr);
double logLimbs(const limb *pBigNbr, int nbrLimbs);
double getMantissa(const limb *ptrLimb, int nbrLimbs);
void UncompressIntLimbs(/*@in@*/const int *ptrValues, /*@out@*/limb *bigint, int nbrLen);
void CompressIntLimbs(/*@out@*/int *ptrValues, /*@in@*/const limb *bigint, int nbrLen);
bool checkOne(const limb *value, int nbrLimbs);
bool checkMinusOne(const limb *value, int nbrLimbs);
void DivideBigNbrByMaxPowerOf2(int *pShRight, limb *number, int *pNbrLimbs);

void ChSignBigNbr(int nbr[], int length);
void ChSignBigNbrB(int nbr[], int length);
void AddBigNbr(const int Nbr1[], const int Nbr2[], int Sum[], int nbrLen);
void SubtractBigNbr(const int Nbr1[], const int Nbr2[], int Diff[], int nbrLen);
void AddBigNbrB(const int Nbr1[], const int Nbr2[], int Sum[], int nbrLen);
void SubtractBigNbrB(const int Nbr1[], const int Nbr2[], int Diff[], int nbrLen);
void AddBigIntModN(const int Nbr1[], const int Nbr2[], int Sum[], const int Mod[], int nbrLen);
void SubtractBigNbrModN(const int Nbr1[], const int Nbr2[], int Diff[], const int Mod[], int nbrLen);
void MultBigNbrByInt(const int bigFactor[], int factor, int bigProduct[], int nbrLen);
void MultBigNbrByIntB(const int bigFactor[], int factor, int bigProduct[], int nbrLen);
void DivBigNbrByInt(const int Dividend[], int divisor, int Quotient[], int nbrLen);
int RemDivBigNbrByInt(const int Dividend[], int divisor, int nbrLen);
void MultBigNbr(const int Factor1[], const int Factor2[], int Prod[], int nbrLen);
void IntToBigNbr(int value, int bigNbr[], int nbrLength);
int BigIntToBigNbr(const BigInteger &pBigInt, int BigNbr[]);
void BigNbrToBigInt(BigInteger &pBigInt, const int BigNbr[], int nbrLenBigInt);
void GcdBigNbr(const int *pNbr1, const int *pNbr2, int *pGcd, int nbrLen);
void MultBigNbrModN(int Nbr1[], int Nbr2[], int Prod[], const int Mod[], int nbrLen);
void MultBigNbrByIntModN(int Nbr1[], int Nbr2, int Prod[], const int Mod[], int nbrLen);
int intDoubleModPow(int NbrMod, int Expon, int currentPrime);
void ModInvBigInt(int *num, int *inv, int *mod, int NumberLength);
//void IntToBigNbr(int value, int *bigNbr, int nbrLength);
int JacobiSymbol(int upper, int lower);

typedef void(*mmCback)(void);
