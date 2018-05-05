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
#include "mpir.h"
#include "boost/multiprecision/gmp.hpp" 
#define ZT(a) a.backend().data()
#define ZisEven(a) (mpz_even_p(ZT(a)) != 0)  /* true iff a is even (works for -ve a as well) */
typedef boost::multiprecision::mpz_int Znum;
long long MulPrToLong(const Znum &x);

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


bool BigNbrIsZero(const limb *value, int nbrLength);
void squareRoot(const limb *argument, limb *sqRoot, int len, int *pLenSqRoot);
double getMantissa(const limb *ptrLimb, int nbrLimbs);

class BigInteger {
public:                  // ideally, values should be private
	limb limbs[MAX_LEN];
	int nbrLimbs = 1;
	enum eSign sign = SIGN_POSITIVE;

public:
	int bitLength() const {
		unsigned int mask;
		int bitCount;
		const int lastLimb = nbrLimbs - 1;
		int bitLen = lastLimb*BITS_PER_GROUP;
		unsigned int limb = (unsigned int)(limbs[lastLimb].x);
		mask = 1;
		for (bitCount = 0; bitCount < BITS_PER_GROUP; bitCount++) {
			if (limb < mask) {
				break;
			}
			mask *= 2;
		}
		return bitLen + bitCount;
	}
	bool isEven() const {
		return ((limbs[0].x & 1) == 0);
	}
	BigInteger sqRoot() const {
		BigInteger sqrRoot;
		squareRoot((*this).limbs, sqrRoot.limbs, (*this).nbrLimbs, &sqrRoot.nbrLimbs);
		sqrRoot.sign = SIGN_POSITIVE;
		return sqrRoot;
	}
	long long lldata() const { /* convert to long long */
		int noOfLimbs; 
		long long rv = 0;
		for (noOfLimbs = nbrLimbs - 1; noOfLimbs >= 0; noOfLimbs--) {
			rv *= LIMB_RANGE;
			rv += limbs[noOfLimbs].x;
		}
		if (sign == SIGN_NEGATIVE)
			rv = -rv;
		return rv;  // no check for overflow; just return last 64 bits
	}
	double fdata() const { /* convert to double */
		double rv = getMantissa(limbs+nbrLimbs, nbrLimbs);
		if (sign == SIGN_NEGATIVE)
			rv = -rv;
		return rv;  
	}
	/* overload assignment operator here 
	there are 5 overloads, for assignments from BigIntegers, Integers, long long,
	double and Znum */
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
	BigInteger & operator = (double value) {
		DoubleToBigInt (*this, value);
		return *this;
	}
	BigInteger & operator = (Znum x) {
		ZtoBig(*this, x);
		return *this;
	}
	/* other operators are overloaded later */
	BigInteger(long long i=0) {   // constructor
		*this = i;
	}
	friend BigInteger BigIntAdd      (const BigInteger &Addend1, const BigInteger &Addend2);
	friend BigInteger BigIntSubt     (const BigInteger &Minuend, const BigInteger &Subtrahend);
	friend static void BigIntNegate  (BigInteger &pDest);
	friend BigInteger BigIntDivide   (const BigInteger &Dividend, const BigInteger &Divisor);
	friend BigInteger BigIntDivideInt(const BigInteger &Dividend, const int Divisor);
	friend BigInteger BigIntMultiply (const BigInteger &Factor1, const BigInteger &Factor2);
	friend BigInteger BigIntRemainder(const BigInteger &Dividend, const BigInteger &Divisor);
	/* calculate base^expon. Throw exception if result is out of range */
	friend void BigIntPowerIntExp    (const BigInteger &Base, int expon, BigInteger &Power);
	friend void BigIntGcd (const BigInteger &pArg1, const BigInteger &pArg2, BigInteger &Result);
	friend void BigIntModularDivision(const BigInteger &Num, const BigInteger &Den, 
		const BigInteger &mod, BigInteger &quotient);
	friend void subtractdivide(BigInteger &BigInt, int subt, int divisor);
	friend void addbigint     (BigInteger &Result, int addend); // Result += addend

	friend int getRemainder    (const BigInteger &pBigInt, int divisor);  // BigInt%divisor
	friend bool TestBigNbrEqual(const BigInteger &Nbr1, const BigInteger &Nbr2);
	friend bool TestBigNbrLess (const BigInteger &Nbr1, const BigInteger &Nbr2);
	friend void BigIntDivide2  (BigInteger &Arg);   // arg /=2;
	friend void expBigInt      (BigInteger &BigInt, double logar); /* BigInt = e^logar */
	friend void DoubleToBigInt (BigInteger &bigInt, double dvalue);
	friend double logBigNbr(const BigInteger &BigInt); /* natural log of BigInt */
	friend static void BigIntMutiplyPower2(BigInteger &pArg, int power2);
	friend void IntsToBigInteger(/*@in@*/const int *ptrValues, /*@out@*/BigInteger &bigint);
	//friend void BigIntegerToInts(/*@out@*/int *ptrValues, /*@in@*/const BigInteger &bigint);
	friend void LimbsToBigInteger(/*@in@*/const limb *ptrValues, 
		/*@out@*/BigInteger &bigint, int NumberLength);
	friend void BigIntegerToLimbs(/*@out@*/limb *ptrValues, 
		/*@in@*/const BigInteger &bigint, int NumberLength);
	//friend int PowerCheck(const BigInteger &pBigNbr, BigInteger &pBase);
	friend void DoubleToBigInt(BigInteger &bigInt, double dvalue);
	friend bool ZtoBig(BigInteger &number, Znum numberZ);
	friend void BigtoZ(Znum &numberZ, const BigInteger &number);

	BigInteger operator +  (const BigInteger &b) const {
		return BigIntAdd(*this, b);
	}
	BigInteger &operator += (const BigInteger &b) {
		*this = BigIntAdd(*this, b);
		return *this;
	}
	BigInteger operator -  (const BigInteger &b) const {
		return BigIntSubt(*this, b);
	}
	BigInteger &operator -= (const BigInteger &b) {
		*this = BigIntSubt(*this, b);
		return *this;
	}
	BigInteger operator /  (const BigInteger &b) const {
		return BigIntDivide(*this, b);
	}
	BigInteger &operator /= (const BigInteger &b) {
		return *this = BigIntDivide(*this, b);
	}
	BigInteger  operator *  (const BigInteger &b) const {
		return BigIntMultiply(*this, b);
	}
	BigInteger &operator *= (const BigInteger &b) {
		return *this = BigIntMultiply(*this, b);
	}
	BigInteger  operator %  (const BigInteger &Divisor) const {
		return BigIntRemainder(*this, Divisor);
	}

	/* note that only a few operators are oveloaded for BigInt <op> int */
	int         operator %  (int divisor) const {
		return getRemainder(*this, divisor);
	}
	BigInteger operator / (int divisor) const {
		return BigIntDivideInt(*this, divisor);
	}
	BigInteger &operator /= (int divisor) {
		subtractdivide(*this, 0, divisor);
		return *this;
	}
	BigInteger &operator += (const int b) {
		addbigint(*this, b);
		return *this;
	}
	BigInteger &operator ++ (int x) {
		addbigint(*this, 1);
		return *this;
	}
	BigInteger &operator -- (int x) {
		addbigint(*this, -1);
		return *this;
	}
	BigInteger &operator -= (const int b) {
		addbigint(*this, -b);
		return *this;
	}

	/* overload comparison operators here */
	bool operator ==(const BigInteger &b) const {
		return TestBigNbrEqual(*this, b);
	}
	bool operator ==(const long long b) const {
		if (b == 0)
			return BigNbrIsZero((*this).limbs, (*this).nbrLimbs);
		BigInteger Btemp;
		Btemp = b;
		return TestBigNbrEqual(*this, Btemp);
	}
	bool operator !=(const BigInteger &b) const {
		return !TestBigNbrEqual(*this, b);
	}
	bool operator !=(const long long b) const {
		BigInteger Btemp;
		Btemp = b;
		return !TestBigNbrEqual(*this, Btemp);
	}
	bool operator < (const BigInteger &b) const {
		return TestBigNbrLess(*this, b);
	}
	bool operator < ( const long long b) const {
		if (b == 0) {
			return ((*this).sign == SIGN_NEGATIVE);
		}
		BigInteger Btemp = b;
		return TestBigNbrLess(*this, Btemp);
	}
	bool operator > (const BigInteger &b) const {
		return TestBigNbrLess(b, *this);
	}
	bool operator > (const long long b) const {
		BigInteger Btemp = b;
		return TestBigNbrLess(Btemp, *this);
	}
	bool operator <=(const BigInteger &b) const {
		return !TestBigNbrLess(b, *this);
	}
	bool operator <=(const long long b) const {
		BigInteger Btemp = b;
		return !TestBigNbrLess(Btemp, *this);
	}
	bool operator >=(const BigInteger &b) const {
		return !TestBigNbrLess(*this, b);
	}
	bool operator >=(const long long b) const {
		if (b == 0)
			return ((*this).sign == SIGN_POSITIVE);
		BigInteger Btemp = b;
		return !TestBigNbrLess(*this, Btemp);
	}
} ;

extern limb TestNbr[MAX_LEN];
extern limb MontgomeryMultR2[MAX_LEN];
extern limb MontgomeryMultR1[MAX_LEN];
extern int NumberLength, NumberLengthR1;
extern int groupLen;
long long PowerCheck(const Znum &BigInt, Znum &Base, long long upperBound);
double logBigNbr(const Znum &BigInt);
void BigIntPowerIntExp(const Znum &Base, int exponent, Znum &Power);
void multiply(const limb *factor1, const limb *factor2, limb *result, int len, int *ResultLen);
void int2dec(char **pOutput, long long nbr);
void GetMontgomeryParms(int len);
void AddBigNbrModNB (const limb *Nbr1, const limb *Nbr2, limb *Sum, const limb *TestNbr, int NumberLength);
void SubtBigNbrModN(const limb *Nbr1, const limb *Nbr2, limb *Sum, const limb *TestNbr, int NumberLength);
void SubtBigNbrMod (const limb *Nbr1, const limb *Nbr2, limb *Sum);
void modmult(const limb *factor1, const limb *factor2, limb *product);
void modmultInt(limb *factorBig, int factorInt, limb *result);
void modmultIntExtended(limb *factorBig, int factorInt, limb *result, const limb *pTestNbr, int nbrLen);
void AddBigNbrMod(limb *Nbr1, limb *Nbr2, limb *Sum);
void modPowBaseInt(int base, const limb *exp, int nbrGroupsExp, limb *power);
void modPow(const limb *base, const limb *exp, int nbrGroupsExp, limb *power);
void AdjustModN(limb *Nbr, const limb *TestNbr, int NumberLength);

void ModInvBigNbr(limb *num, limb *inv, limb *mod, int NumberLength);

//int BpswPrimalityTest(const BigInteger &pBigNbr);
int PrimalityTest(const Znum &Value, long long upperBound);

void ComputeInversePower2(/*@in@*/const limb *value, /*@out@*/limb *result, /*@out@*/limb *aux);
double logLimbs(const limb *pBigNbr, int nbrLimbs);
//void UncompressIntLimbs(/*@in@*/const int *ptrValues, /*@out@*/limb *bigint, int nbrLen);
//void CompressIntLimbs(/*@out@*/int *ptrValues, /*@in@*/const limb *bigint, int nbrLen);
bool checkOne(const limb *value, int nbrLimbs);
bool checkMinusOne(const limb *value, int nbrLimbs);
//void DivideBigNbrByMaxPowerOf2(int *pShRight, limb *number, int *pNbrLimbs);
void DivideBigNbrByMaxPowerOf2(int &ShRight, Znum &number);
void ChSignBigNbr(int nbr[], int length);
void ChSignBigNbrB(int nbr[], int length);
void AddBigNbr(const int Nbr1[], const int Nbr2[], int Sum[], int nbrLen);
void SubtractBigNbr(const int Nbr1[], const int Nbr2[], int Diff[], int nbrLen);
void AddBigNbrB(const int Nbr1[], const int Nbr2[], int Sum[], int nbrLen);
void SubtractBigNbrB(const int Nbr1[], const int Nbr2[], int Diff[], int nbrLen);
void AddBigNbrModN(const int Nbr1[], const int Nbr2[], int Sum[], const int Mod[], int nbrLen);
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
void MultBigNbrModN(Znum Nbr1, Znum Nbr2, Znum Prod, const Znum Mod);
void MultBigNbrByIntModN(int Nbr1[], int Nbr2, int Prod[], const int Mod[], int nbrLen);
void MultBigNbrByIntModN(Znum Nbr1, int Nbr2, Znum Prod, const Znum Mod);
int intDoubleModPow(int NbrMod, int Expon, int currentPrime);
void ModInvBigInt(int *num, int *inv, int *mod, int NumberLength);

long long JacobiSymbol(long long upper, long long lower);

typedef void(*mmCback)(void);
