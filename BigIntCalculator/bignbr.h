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
#ifdef __GNUC__
#include <gmp.h>
typedef unsigned long mpir_ui;  // for compatibility with MPIR
#else
#include <mpir.h>
#endif
#include "boost/multiprecision/gmp.hpp"
#define ZT(a) a.backend().data()
#define ZisEven(a) (mpz_even_p(ZT(a)) != 0)  /* true iff a is even (works for -ve a as well) */
typedef boost::multiprecision::mpz_int Znum;

#define MAX_LEN 2500        // approximately 20000 digits
#define BITS_PER_GROUP 31
#define BITS_PER_INT_GROUP 31
#define HALF_INT_RANGE (1 << (BITS_PER_INT_GROUP - 1))
#define MAX_INT_NBR ((int)((1U << BITS_PER_INT_GROUP)-1))
#define LIMB_RANGE (1U<<BITS_PER_GROUP)
#define SMALL_NUMBER_BOUND 32768
#define _USING64BITS_ 1

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

extern long long lModularMult;
bool BigNbrIsZero(const limb *value, int nbrLength);
//void squareRoot(const limb *argument, limb *sqRoot, int len, int *pLenSqRoot);
void squareRoot(const Znum &arg, Znum & sqRoot);
bool isPerfectSquare(const Znum &arg, Znum &sqRoot);
bool isPerfectSquare(__int64 x);
double getMantissa(const limb *ptrLimb, int nbrLimbs);

class BigInteger {

public:				// should be private members
	limb limbs[MAX_LEN];
	int nbrLimbs = 1;
	enum eSign sign = SIGN_POSITIVE;

public:
	int bitLength() const {
		const int lastLimb = nbrLimbs - 1;
		int bitLen = lastLimb*BITS_PER_GROUP;
		unsigned int limb = (unsigned int)(limbs[lastLimb].x);
		if (limb != 0) {
			/* use _BitScanReverse instead of loop because it's faster */
			unsigned long bitCount;
			_BitScanReverse(&bitCount, limb);
			bitLen += (bitCount + 1);
		}
		return bitLen;
	}
	bool isEven() const {
		return ((limbs[0].x & 1) == 0);
	}
	/*BigInteger sqRoot() const {
		BigInteger sqrRoot;
		squareRoot((*this).limbs, sqrRoot.limbs, (*this).nbrLimbs, &sqrRoot.nbrLimbs);
		sqrRoot.sign = SIGN_POSITIVE;
		return sqrRoot;
	}*/
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
	double fpoint() const { /* convert to double */
		double rv = getMantissa(limbs+nbrLimbs, nbrLimbs);
		if (sign == SIGN_NEGATIVE)
			rv = -rv;
		if (nbrLimbs > 1) {
			rv *= pow(LIMB_RANGE, nbrLimbs - 1);
		}
		auto c = fpclassify(rv);
		if (c == FP_INFINITE) {
			std::string line = std::to_string(__LINE__);
			std::string mesg = "floating point overflow ";
			mesg += __func__;
			mesg += " line ";  mesg += line;
			mesg += " in file "; mesg += __FILE__;
			throw std::range_error(mesg);
		}
		return rv;  
	}

	/* overload assignment operator here. There are 5 overloads, 
	for assignments from BigIntegers, Integers, long long,
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
		auto rv = ZtoBig(*this, x);
		if (!rv) {
			std::string line = std::to_string(__LINE__);
			std::string mesg = "Number too big to convert ";
			mesg += __func__;
			mesg += " line ";  mesg += line;
			mesg += " in file "; mesg += __FILE__;
			throw std::range_error(mesg);
		}
		return *this;
	}
	/* other operators are overloaded later */

	 // constructor
	BigInteger(long long i=0) {  
		*this = i;
	}
	friend BigInteger BigIntAdd      (const BigInteger &Addend1, const BigInteger &Addend2);
	friend BigInteger BigIntSubt     (const BigInteger &Minuend, const BigInteger &Subtrahend);
	friend static void BigIntNegate  (BigInteger &pDest);
	friend BigInteger BigIntDivide   (const BigInteger &Dividend, const BigInteger &Divisor);
	friend BigInteger BigIntDivideInt(const BigInteger &Dividend, const int Divisor);
	friend BigInteger BigIntMultiply (const BigInteger &Factor1, const BigInteger &Factor2);
	friend void MultBigNbrByInt(BigInteger &m, int n); // m *= n
	friend BigInteger BigIntRemainder(const BigInteger &Dividend, const BigInteger &Divisor);
	/* calculate base^expon. Throw exception if result is out of range */
	//friend void BigIntPowerIntExp    (const BigInteger &Base, int expon, BigInteger &Power);
	//friend void BigIntGcd (const BigInteger &pArg1, const BigInteger &pArg2, BigInteger &Result);
	//friend void BigIntModularDivision(const BigInteger &Num, const BigInteger &Den, 
		//const BigInteger &mod, BigInteger &quotient);
	friend void subtractdivide(BigInteger &BigInt, int subt, int divisor);
	friend void addbigint     (BigInteger &Result, int addend); // Result += addend

	friend int getRemainder    (const BigInteger &pBigInt, int divisor);  // BigInt%divisor
	friend bool TestBigNbrEqual(const BigInteger &Nbr1, const BigInteger &Nbr2);
	friend bool TestBigNbrLess (const BigInteger &Nbr1, const BigInteger &Nbr2);
	//friend void BigIntDivide2  (BigInteger &Arg);   // arg /=2;
	//friend void expBigInt      (BigInteger &BigInt, double logar); /* BigInt = e^logar */
	friend void DoubleToBigInt (BigInteger &bigInt, double dvalue);
	//friend double logBigNbr(const BigInteger &BigInt); /* natural log of BigInt */
	//friend static void BigIntMutiplyPower2(BigInteger &pArg, int power2);
	//friend void IntsToBigInteger(/*@in@*/const int *ptrValues, /*@out@*/BigInteger &bigint);
	//friend void BigIntegerToInts(/*@out@*/int *ptrValues, /*@in@*/const BigInteger &bigint);
	friend void LimbsToBigInteger(/*@in@*/const limb *ptrValues, 
		/*@out@*/BigInteger &bigint, int NumLen);
	friend void BigIntegerToLimbs(/*@out@*/limb *ptrValues, 
		/*@in@*/const BigInteger &bigint, int NumLen);
	//friend int PowerCheck(const BigInteger &pBigNbr, BigInteger &pBase);
	friend void DoubleToBigInt(BigInteger &bigInt, double dvalue);
	friend bool ZtoBig(BigInteger &number, Znum numberZ);
	friend void BigtoZ(Znum &numberZ, const BigInteger &number);
	//friend int BigIntToBigNbr(const BigInteger &pBigInt, int BigNbr[]);
	//friend void BigNbrToBigInt(BigInteger &pBigInt, const int BigNbr[], int nbrLenBigInt);
	//friend void ModInvBigNbr(int *num, int *inv, int *mod, int NumLen);
	friend void shiftBI(const BigInteger &first, const int shiftCtr, BigInteger &result);
	friend double logBigNbr(const BigInteger &BigInt);
	friend void expBigInt(BigInteger &bigInt, double logar);

	BigInteger  operator +  (const BigInteger &b) const {
		return BigIntAdd(*this, b);
	}
	BigInteger &operator += (const BigInteger &b) {
		*this = BigIntAdd(*this, b);
		return *this;
	}
	BigInteger  operator -  (const BigInteger &b) const {
		return BigIntSubt(*this, b);
	}
	BigInteger &operator -= (const BigInteger &b) {
		*this = BigIntSubt(*this, b);
		return *this;
	}
	BigInteger  operator /  (const BigInteger &b) const {
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
	BigInteger &operator %= (const BigInteger &b) {
		return *this = BigIntRemainder(*this, b);
	}

	int         operator %  (int divisor) const {
		return getRemainder(*this, divisor);
	}
	BigInteger operator / (int divisor) const {
		return BigIntDivideInt(*this, divisor);
	}
	BigInteger &operator /= (int divisor) {
		*this = BigIntDivideInt(*this, divisor);
		return *this;
	}
	BigInteger operator * (int m) const {
		BigInteger prod =*this;
		MultBigNbrByInt(prod, m);
		return prod;
	}
	BigInteger &operator *= (int b) {
		MultBigNbrByInt(*this, b);
		return *this;
	}

	BigInteger &operator += (const int b) {
		addbigint(*this, b);
		return *this;
	}
	BigInteger operator + (const int b) const {
		BigInteger sum = *this;
		sum += b;
		return sum;
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
	BigInteger operator - (const int b) const {
		BigInteger sum = *this;
		sum -= b;
		return sum;
	}

	BigInteger operator << (const int b) const {
		BigInteger temp;
		shiftBI(*this, b, temp);
		return temp;
	}
	BigInteger &operator <<= (const int b) {
		BigInteger temp;
		shiftBI(*this, b, temp);
		*this = temp;
	}
	BigInteger operator >>= (const int b) {
		BigInteger temp;
		shiftBI(*this, -b, temp);
		*this = temp;
	}
	BigInteger &operator >> (const int b) const {
		BigInteger temp;
		shiftBI(*this, -b, temp);
		return temp;
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

void Bin2Dec(const BigInteger &BigInt, char *decimal, int groupLength);
/* output operators declaration */
std::ostream  &operator<<(std::ostream  &s, const BigInteger &n);
std::wostream &operator<<(std::wostream &s, const BigInteger &n);

extern BigInteger TestNbrBI;
extern limb * const TestNbr;
extern limb *const MontgomeryMultR2;
extern limb *const MontgomeryMultR1;
#define NumberLength TestNbrBI.nbrLimbs
//extern int groupLen;

void ModInvBigNbr(Znum &num, Znum &inv, Znum &mod);
//void FactoringSIQSx(const Znum &NbrToFactor, Znum &Factor);
void FactoringSIQS(const Znum &NbrToFactor, Znum &Factor);
void multiply(const limb *factor1, const limb *factor2, limb *result, int len, int *ResultLen);
void int2dec(char **pOutput, long long nbr);
void Bin2Dec(const limb *binary, char *decimal, int nbrLimbs, int groupLength);
void GetMontgomeryParms(int len);
void GetMontgomeryParms(const Znum &Nval);
void AddBigNbrModNB (const limb Nbr1[], const limb Nbr2[], limb Sum[], const limb TestNbr[], int NumLen);
void SubtBigNbrModN(const limb Nbr1[], const limb Nbr2[], limb Sum[], const limb TestNbr[], int NumLen);
void SubtBigNbrModN(const BigInteger &Nbr1, const BigInteger &Nbr2, BigInteger &Diff);
//void SubtBigNbrMod (const limb *Nbr1, const limb *Nbr2, limb *Sum);
void modmult(const limb *factor1, const limb *factor2, limb *product);
void modmult(const BigInteger &factor1, const BigInteger &factor2, BigInteger &product);
void modmult(const Znum &a, const Znum &b, Znum &result);
void REDC(Znum &t, const Znum &T);
void modmultInt(const limb *factorBig, int factorInt, limb *result);
void modmultInt(const BigInteger &factorBig, int factorInt, BigInteger &result);
//void modmultIntExtended(const limb *factorBig, int factorInt, limb *result, const limb *pTestNbr, int nbrLen);
//void AddBigNbrMod(const limb *Nbr1, const limb *Nbr2, limb *Sum);
//void modPowBaseInt(int base, const limb *exp, int nbrGroupsExp, limb *power);
//void modPow(const limb *base, const limb *exp, int nbrGroupsExp, limb *power);
//void AdjustModN(limb *Nbr, const limb *TestNbr, int NumLen);
void AdjustModN(Znum &Nbr, const Znum &Modulus);
void ModInvBigNbr(const limb *num, limb *inv, const limb *mod, int NumLen);
void ModInvBigNbr(const BigInteger &num, BigInteger &inv, const BigInteger &mod);
void ValuestoZ(Znum &numberZ, const int number[], int NumLen);
//static void ComputeInversePower2(/*@in@*/const limb *value, /*@out@*/limb *result, /*@out@*/limb *aux);
//double logLimbs(const limb *pBigNbr, int nbrLimbs);
//void UncompressIntLimbs(/*@in@*/const int *ptrValues, /*@out@*/limb *bigint, int nbrLen);
//void CompressIntLimbs(/*@out@*/int *ptrValues, /*@in@*/const limb *bigint, int nbrLen);
//bool checkOne(const limb *value, int nbrLimbs);
//bool checkMinusOne(const limb *value, int nbrLimbs);
//void DivideBigNbrByMaxPowerOf2(int *pShRight, limb *number, int *pNbrLimbs);

void ChSignBigNbr(int nbr[], int length);
void ChSignBigNbrB(int nbr[], int length);
int BigNbrLen(const long long Nbr[], int nbrLen);
//void AddBigNbr(const int Nbr1[], const int Nbr2[], int Sum[], int nbrLen);
//void AddBigNbr(const Znum &Nbr1, const Znum &Nbr2, Znum &Sum);
//void SubtractBigNbr(const int Nbr1[], const int Nbr2[], int Diff[], int nbrLen);
//void SubtractBigNbr(const Znum &Nbr1, const Znum &Nbr2, Znum &Sum);
//void AddBigNbrB(const int Nbr1[], const int Nbr2[], int Sum[], int nbrLen);
//void AddBigNbrB(const Znum &Nbr1, const Znum &Nbr2, Znum &Sum);
//void SubtractBigNbrB(const int Nbr1[], const int Nbr2[], int Diff[], int nbrLen);
//void SubtractBigNbrB(const Znum &Nbr1, const Znum &Nbr2, Znum &Sum);
//void AddBigNbrModN(const int Nbr1[], const int Nbr2[], int Sum[], const int Mod[], int nbrLen);
void AddBigNbrModN(const Znum &Nbr1, const Znum &Nbr2, Znum &Diff, const Znum &Mod);
//void SubtractBigNbrModN(const int Nbr1[], const int Nbr2[], int Diff[], const int Mod[], int nbrLen);
void SubtractBigNbrModN(const Znum &Nbr1, const Znum &Nbr2, Znum &Diff, const Znum &Mod);
//void MultBigNbrByInt(const int bigFactor[], int factor, int bigProduct[], int nbrLen);
//void MultBigNbrByInt(const Znum &bigFactor, int factor, Znum &bigProd);
//void MultBigNbrByIntB(const int bigFactor[], int factor, int bigProduct[], int nbrLen);
//void MultBigNbrByIntB(const Znum &bigFactor, int factor, Znum &bigProd);
void DivBigNbrByInt(const int Dividend[], int divisor, int Quotient[], int nbrLen);
//void DivBigNbrByInt(const Znum &Dividend, int divisor, Znum &Quotient);
//int RemDivBigNbrByInt(const int Dividend[], int divisor, int nbrLen);
mpir_ui RemDivBigNbrByInt(const Znum &Dividend, mpir_ui divisor);
//void MultBigNbr(const int Factor1[], const int Factor2[], int Prod[], int nbrLen);
//void MultBigNbr(const Znum &Fact1, const Znum &Fact2, Znum &Prod);
//void IntToBigNbr(int value, int bigNbr[], int nbrLength);
//void GcdBigNbr(const int *pNbr1, const int *pNbr2, int *pGcd, int nbrLen);
//void MultBigNbrModN(const int Nbr1[], const int Nbr2[], int Prod[], const int Mod[], int nbrLen);
void MultBigNbrModN(const Znum &Nbr1, const Znum &Nbr2, Znum &Prod, const Znum &Mod);
//void MultBigNbrByIntModN(const int Nbr1[], int Nbr2, int Prod[], const int Mod[], int nbrLen);
void MultBigNbrByIntModN(const Znum &Nbr1, int Nbr2, Znum &Prod, const Znum &Mod);
int modPower(int NbrMod, int Expon, int currentPrime);
// calculate a^n%mod using 'bigints'   
unsigned __int64 modPower(unsigned __int64 a, unsigned __int64 n,
	unsigned __int64 mod);
int ZtoLimbs(limb *number, Znum numberZ, int NumLen);
int ZtoBigNbr(int number[], Znum numberZ);
void LimbstoZ(const limb *number, Znum &numberZ, int NumLen);

typedef void(*mmCback)(void);
extern mmCback modmultCallback;
