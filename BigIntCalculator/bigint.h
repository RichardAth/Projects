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

/* BigInteger has been made into a class, with operator overload so that normal
operators do what you would expect. However, these Bigintegers are now only used
by modmult functions which in turn are only used by the built-in ECM function. 
Oct 2022 update: Solving Modular Quadratic Equations makes extensive use of
BIgIntegers. As it uses modular arithmetic in obscure ways, and  is is fast
enough, switching to Znums does not seem worthwhile. */

enum eSign
{
	SIGN_POSITIVE = 0,
	SIGN_NEGATIVE,
};
void squareRoot(const limb* argument, limb* sqRoot, int len, int* pLenSqRoot);


class BigInteger {

public:				// should be private members
	limb limbs[MAX_LEN];
	int nbrLimbs = 1;
	enum eSign sign = SIGN_POSITIVE;

public:
	int bitLength() const {
		const int lastLimb = nbrLimbs - 1;
		int bitLen = lastLimb * BITS_PER_GROUP;
		unsigned int limb = (unsigned int)(limbs[lastLimb]);
		if (limb != 0) {
			/* use _BitScanReverse instead of loop because it's faster */
			unsigned long bitCount;
			_BitScanReverse(&bitCount, limb);
			bitLen += (bitCount + 1);
		}
		return bitLen;
	}
	bool isEven() const {
		return ((limbs[0] & 1) == 0);
	}
	BigInteger sqRoot() const {
		BigInteger sqrRoot;
		squareRoot((*this).limbs, sqrRoot.limbs, (*this).nbrLimbs, &sqrRoot.nbrLimbs);
		sqrRoot.sign = SIGN_POSITIVE;
		return sqrRoot;
	}
	BigInteger abs() const {
		BigInteger result = *this;
		if (result.sign == SIGN_NEGATIVE)
			result.sign = SIGN_POSITIVE;
		return result;
	}
	long long lldata() /* convert to long long */ const {
		int noOfLimbs;
		long long rv = 0;
		for (noOfLimbs = nbrLimbs - 1; noOfLimbs >= 0; noOfLimbs--) {
			rv *= LIMB_RANGE;
			rv += limbs[noOfLimbs];
		}
		if (sign == SIGN_NEGATIVE)
			rv = -rv;
		return rv;  // no check for overflow; just return last 64 bits
	}
	double fpoint() /* convert to double */ const {
		double rv = getMantissa(limbs + nbrLimbs, nbrLimbs);
		if (sign == SIGN_NEGATIVE)
			rv = -rv;
		if (nbrLimbs > 1) {
			rv *= std::pow(LIMB_RANGE, nbrLimbs - 1);
		}
		auto c = fpclassify(rv);
		if (c == FP_INFINITE) {
			StackTrace2();
			std::string line = std::to_string(__LINE__);
			std::string mesg = "floating point overflow ";
			mesg += __func__;
			mesg += " line ";  mesg += line;
			mesg += " in file "; mesg += __FILE__;
			throw std::range_error(mesg);
		}
		return rv;
	}
	double log() const {
		return logBigInt(*this);
	}

	/* overload assignment operator here. There are 5 overloads,
	for assignments from BigIntegers, Integers, long long,
	double and Znum */
	BigInteger & operator = (const BigInteger &other) {
		if (&other == this)
			return *this;		// if lhs == rhs do nothing
		nbrLimbs = other.nbrLimbs;
		sign = other.sign;
		std::memcpy(limbs, other.limbs, nbrLimbs * sizeof(int));
		while (nbrLimbs > 1 && limbs[nbrLimbs - 1] == 0) {
			nbrLimbs--;  // remove any leading zeros
		}
		return *this;
	}
	BigInteger & operator = (const int value) {
		if (value >= 0) {
			limbs[0] = value;
			sign = SIGN_POSITIVE;
		}
		else {
			limbs[0] = -value;
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
			limbs[noOfLimbs++] = (int)value & MAX_VALUE_LIMB;
			value >>= BITS_PER_GROUP;
		} while (value != 0);

		nbrLimbs = noOfLimbs;
		return *this;
	}
	BigInteger & operator = (double value) {
		DoubleToBigInt(*this, value);
		return *this;
	}
	BigInteger & operator = (Znum x) {
		auto rv = ZtoBig(*this, x);
		if (!rv) {
			StackTrace2();
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
	BigInteger(long long i = 0) {
		*this = i;
	}

	friend BigInteger BigIntAdd(const BigInteger &Addend1, const BigInteger &Addend2);
	friend BigInteger BigIntSubt(const BigInteger &Minuend, const BigInteger &Subtrahend);
	friend void BigIntNegate(BigInteger &pDest);  /* Dest = -Dest */
	friend BigInteger BigIntDivide(const BigInteger &Dividend, const BigInteger &Divisor);
	friend BigInteger BigIntDivideInt(const BigInteger &Dividend, const int Divisor);
	friend BigInteger BigIntMultiply(const BigInteger &Factor1, const BigInteger &Factor2);
	friend void MultBigIntByInt(BigInteger &m, int n); // m *= n
	friend BigInteger BigIntRemainder(const BigInteger &Dividend, const BigInteger &Divisor);
	/* calculate base^expon. Throw exception if result is out of range */
	friend void BigIntPowerIntExp    (const BigInteger &Base, int expon, BigInteger &Power);
	friend void BigIntGcd (const BigInteger &pArg1, const BigInteger &pArg2, BigInteger &Result);
	//friend void BigIntModularDivision(const BigInteger &Num, const BigInteger &Den, 
		//const BigInteger &mod, BigInteger &quotient);
	friend void subtractdivide(BigInteger &BigInt, int subt, int divisor);
	friend void addbigint(BigInteger &Result, int addend); // Result += addend

	friend int getRemainder(const BigInteger &pBigInt, int divisor);  // BigInt%divisor
	friend bool TestBigIntEqual(const BigInteger &Nbr1, const BigInteger &Nbr2);
	friend bool TestBigIntLess(const BigInteger &Nbr1, const BigInteger &Nbr2);
	friend void BigIntDivide2  (BigInteger &Arg);   // arg /=2;
	friend void expBigInt      (BigInteger &BigInt, double logar); /* BigInt = e^logar */
	friend void DoubleToBigInt(BigInteger &bigInt, double dvalue);
	friend double logBigInt(const BigInteger &BigInt); /* natural log of BigInt */
	friend static void BigIntMutiplyPower2(BigInteger &pArg, int power2);
	//friend void IntsToBigInteger(/*@in@*/const int *ptrValues, /*@out@*/BigInteger &bigint);
	//friend void BigIntegerToInts(/*@out@*/int *ptrValues, /*@in@*/const BigInteger &bigint);
	friend void LimbsToBigInteger(/*@in@*/const limb *ptrValues,
		/*@out@*/BigInteger &bigint, int NumLen);
	friend void BigIntegerToLimbs(/*@out@*/limb ptrValues[],
		/*@in@*/const BigInteger &bigint, int NumLen);
	//friend int PowerCheck(const BigInteger &pBigNbr, BigInteger &pBase);
	friend bool ZtoBig(BigInteger &number, Znum numberZ);
	friend void BigtoZ(Znum &numberZ, const BigInteger &number);
	friend long long BigToLL(const BigInteger &num, int &exp);
	friend void LLToBig(BigInteger &num, long long LL, int exp);
	//friend int BigIntToBigNbr(const BigInteger &pBigInt, int BigNbr[]);
	//friend void BigNbrToBigInt(BigInteger &pBigInt, const int BigNbr[], int nbrLenBigInt);
	
	friend void shiftBI(const BigInteger &first, const int shiftCtr, BigInteger &result);
	friend void BigIntAnd(const BigInteger* firstArg, const BigInteger* secondArg, 
		BigInteger* result);
	friend void BigIntOr(const BigInteger* firstArg, const BigInteger* secondArg,
		BigInteger* result);
	friend void BigIntXor(const BigInteger* firstArg, const BigInteger* secondArg,
		BigInteger* result);

	/* overloads where 2nd operand is a bigInteger */
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

	/* overloads where 2nd operand is an int */
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
		BigInteger prod = *this;
		MultBigIntByInt(prod, m);
		return prod;
	}
	BigInteger &operator *= (int b) {
		MultBigIntByInt(*this, b);
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
		if (b <= 0x3FFFFFFF && b >= -0x3FFFFFFF) {
			addbigint(*this, -b);
			return *this;
		}
		else {
			BigInteger bb = b;
			*this = BigIntSubt(*this, bb);
			return *this;
		}
	}
	BigInteger operator - (const int b) const {
		BigInteger sum = *this;
		sum -= b;
		return sum;
	}

	/* unary - */
	BigInteger operator - () {
		BigInteger result = *this;
		BigIntNegate(result);
		return result;
	}

	/* bitwise and */
	BigInteger operator & (const BigInteger& b) const {
		BigInteger result;
		BigIntAnd(this, &b, &result);
		return result;
	}
	BigInteger operator &= (const BigInteger& b) {
		BigIntAnd(this, &b, this);
		return *this;
	}

	/* bitwise or */
	BigInteger operator | (const BigInteger& b) const {
		BigInteger result;
		BigIntOr(this, &b, &result);
		return result;
	}
	BigInteger operator |= (const BigInteger& b) {
		BigIntOr(this, &b, this);
		return *this;
	}

	/* bitwise xor */
	BigInteger operator ^ (const BigInteger& b) const {
		BigInteger result;
		BigIntXor(this, &b, &result);
		return result;
	}
	BigInteger operator ^= (const BigInteger& b) {
		BigIntXor(this, &b, this);
		return *this;
	}

	/* shift operators */
	BigInteger operator << (const int b) const {
		BigInteger temp;
		shiftBI(*this, b, temp);
		return temp;
	}
	BigInteger &operator <<= (const int b) {
		BigInteger temp;
		shiftBI(*this, b, temp);
		*this = temp;
		return *this;
	}
	BigInteger operator >>= (const int b) {
		BigInteger temp;
		shiftBI(*this, -b, temp);
		*this = temp;
		return *this;
	}
	BigInteger &operator >> (const int b) const {
		BigInteger temp;
		shiftBI(*this, -b, temp);
		return temp;
	}

	/* overload comparison operators here */
	bool operator ==(const BigInteger &b) const {
		return TestBigIntEqual(*this, b);
	}
	bool operator ==(const long long b) const {
		if (b == 0)
			return BigNbrIsZero((*this).limbs, (*this).nbrLimbs);
		BigInteger Btemp;
		Btemp = b;
		return TestBigIntEqual(*this, Btemp);
	}
	bool operator !=(const BigInteger &b) const {
		return !TestBigIntEqual(*this, b);
	}
	bool operator !=(const long long b) const {
		BigInteger Btemp;
		Btemp = b;
		return !TestBigIntEqual(*this, Btemp);
	}
	bool operator < (const BigInteger &b) const {
		return TestBigIntLess(*this, b);
	}
	bool operator < (const long long b) const {
		if (b == 0) {
			return ((*this).sign == SIGN_NEGATIVE);
		}
		BigInteger Btemp = b;
		return TestBigIntLess(*this, Btemp);
	}
	bool operator > (const BigInteger &b) const {
		return TestBigIntLess(b, *this);
	}
	bool operator > (const long long b) const {
		BigInteger Btemp = b;
		return TestBigIntLess(Btemp, *this);
	}
	bool operator <=(const BigInteger &b) const {
		return !TestBigIntLess(b, *this);
	}
	bool operator <=(const long long b) const {
		BigInteger Btemp = b;
		return !TestBigIntLess(Btemp, *this);
	}
	bool operator >=(const BigInteger &b) const {
		return !TestBigIntLess(*this, b);
	}
	bool operator >=(const long long b) const {
		if (b == 0)
			return ((*this).sign == SIGN_POSITIVE);
		BigInteger Btemp = b;
		return !TestBigIntLess(*this, Btemp);
	}
};  /* end of BigInteger class definition */

bool BigIntIsZero(const BigInteger* value);
//void BigIntChSign(BigInteger* value);   /* does same job as BigIntNegate */
void BigIntAbs(BigInteger* value);  /* if value is -ve, change its sign */
void Bin2Dec(const BigInteger &BigInt, char *decimal, int groupLength);
void PrintBigInteger(const BigInteger* pBigInt, int groupLength);
void copyStr(char** pptrString, const char* stringToCopy);
void BigInteger2Dec(char** ppDecimal, const BigInteger* pBigInt, int groupLength);
void IntArray2BigInteger(const int* ptrValues, /*@out@*/BigInteger* bigint);
void DivideBigNbrByMaxPowerOf2(int* pShRight, limb* number, int* pNbrLimbs);
void BigInteger2IntArray(/*@out@*/int* ptrValues, const BigInteger* bigint);
void BigIntPowerOf2(BigInteger* pResult, int exponent);
void AddBigNbr     (const limb* pNbr1, const limb* pNbr2, limb* pSum, int nbrLen);
void SubtractBigNbr(const limb* pNbr1, const limb* pNbr2, limb* pDiff, int nbrLen);
void BigIntDivideBy2(BigInteger* nbr);

void CompressLimbsBigInteger(/*@out@*/limb* ptrValues, const BigInteger* bigint);
void UncompressLimbsBigInteger(const limb* ptrValues, /*@out@*/BigInteger* bigint);
int BigIntJacobiSymbol(const BigInteger* upper, const BigInteger* lower);
void BigIntModularPower(const BigInteger* base, const BigInteger* exponent,
	BigInteger* power);
void BigIntModularDivision(const BigInteger* Num, const BigInteger* Den,
	const BigInteger* mod, BigInteger* quotient);

/* output operators declaration */
std::ostream  &operator<<(std::ostream  &s, const BigInteger &n);
std::wostream &operator<<(std::wostream &s, const BigInteger &n);

extern BigInteger TestNbrBI;
#define NumberLength TestNbrBI.nbrLimbs
