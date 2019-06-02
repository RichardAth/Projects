#pragma once
#undef _CRT_SECURE_NO_WARNINGS   // prevent warning message
#define _CRT_SECURE_NO_WARNINGS 1
#include <cstdio>
#include <iostream>
#include <strstream>
#include <cmath>
#include <cfloat>
#include <assert.h>
#include <intrin.h>


typedef struct {
	uint64_t i[2];
} _S2;

typedef struct {
	unsigned int i[4];
} _S4;

typedef struct {
	unsigned short s[8];
} _S8;

typedef struct {
	unsigned char b[16];
} _S16;


#define HI64(n) (n).s2.i[1]
#define LO64(n) (n).s2.i[0]

/* declarations for asm functions
For X64 architecture use Int128x64.asm */
extern "C" {
	void int128add(void *dst, const void *x);       // dst += x
	void int128sub(void *dst, const void *x);       // dst -= x
	void int128mul(void *dst, const void *x);       // dst *= x
	void int128div(void *dst, void *x);             // dst /= x
	void int128rem(void *dst, void *x);             // dst %=x 
	void int128neg(void *x);                        // x = -x
	void int128inc(void *x);                        // x++
	void int128dec(void *x);                        // x--
	void int128shr(void *x, int shft);              // x >>= shft
	void int128shl(void *x, int shft);              // x <<= shft
	int  int128cmp(const void *x1, const void *x2); /* 3 way compare. 
													return  -ve, 0, or  +ve */
	/* separate unsigned versions are needed for the following: */
	void uint128div(void *dst, const void *x);      // dst /= x
	void uint128rem(void *dst, const void *x);      // dst %=x 
	void uint128shr(void *x, int shft);             // x >>= shft
	int  uint128cmp(const void *x1, const void *x2);/* 3 way compare. 
													return  -ve, 0, or  +ve */
};


class _int128 {
	/* note that the there are no private fields */
public:
	union {
		_S2  s2;
		_S4  s4;
		_S8  s8;
		_S16 s16;
	};


	// constructors
	inline _int128() {}
	inline _int128(const unsigned __int64 &n) {
		HI64(*this) = 0;
		LO64(*this) = n;
	}
	inline _int128(const __int64 &n) {
		HI64(*this) = n < 0 ? -1 : 0;
		LO64(*this) = n;
	}
	inline _int128(unsigned long n) {
		HI64(*this) = 0;
		LO64(*this) = n;
	}
	inline _int128(long n) {
		HI64(*this) = n < 0 ? -1 : 0;
		LO64(*this) = n;
	}
	inline _int128(unsigned int n) {
		HI64(*this) = 0;
		LO64(*this) = n;
	}
	inline _int128(int n) {
		HI64(*this) = n < 0 ? -1 : 0;
		LO64(*this) = n;
	}
	inline _int128(unsigned short n) {
		HI64(*this) = 0;
		LO64(*this) = n;
	}
	inline _int128(short n) {
		HI64(*this) = n < 0 ? -1 : 0;
		LO64(*this) = n;
	}
	explicit inline _int128(const unsigned __int64 &hi, const unsigned __int64 &lo) {
		HI64(*this) = hi;
		LO64(*this) = lo;
	}
	
	/* get approximate value in floating point */
	double ToDouble() const {
		bool neg = false;
		_int128 absThis;
		if ((*this).isNegative()) {
			absThis = *this;
			int128neg(&absThis);   // split into 2 statements to avoid error C2593
			neg = true;
		}
		else
			absThis = *this;

		/* absThis = abs(*this). Sign is in neg */
	
		double l = (double)LO64(absThis);
		if (HI64(absThis) != 0) {
			double h = (double)HI64(absThis);
			l += h * pow(2.0, 64);
		}
		if (!neg)
			return l;
		else
			return -l;
	}

	// type operators
	operator unsigned __int64() const {
		return LO64(*this);
	}
	operator __int64() const {
		return LO64(*this);
	}
	operator unsigned long() const {
		return (unsigned long)LO64(*this);
	}
	operator long() const {
		return (long)LO64(*this);
	}
	operator unsigned int() const {
		return (unsigned int)LO64(*this);
	}
	operator int() const {
		return (int)LO64(*this);
	}
	inline operator bool() const {
		return LO64(*this) || HI64(*this);
	}
	operator double() const {
		return this->ToDouble();
	}

	// assign operators
	inline _int128 &operator++() {   // prefix-form
		int128inc(this);
		return *this;
	}
	inline _int128 &operator--() {   // prefix-form
		int128dec(this);
		return *this;
	}

	inline _int128 operator++(int) { // postfix-form
		const _int128 result(*this);
		int128inc(this);
		return result;
	}
	inline _int128 operator--(int) { // postfix-form
		const _int128 result(*this);
		int128dec(this);
		return result;
	}

	inline bool isNegative() const {
		return ((int)s4.i[3] < 0);
	}
	inline bool isZero() const {
		return LO64(*this) == 0 && HI64(*this) == 0;
	}

	/* count number of bits set */
	int popcnt() const {
		auto u = (int)__popcnt64(HI64(*this));  // count 1-bits in top half
		auto l = (int)__popcnt64(LO64(*this));  // count 1-bits in bottom half 
		return (int)u + l;    /* return total 1 bits */
	}

	/* get bit number of most significant 1 bit
	bits are numbered from 0 to 127 */
	int msb() const {
		DWORD index;
		if (HI64(*this) != 0) {
			auto rv = _BitScanReverse64(&index, HI64(*this));
			return index + 64;
		}
		else if (LO64(*this) != 0) {
			auto rv = _BitScanReverse64(&index, LO64(*this));
			return index;
		}
		else return -1;
	}

	/* get bit number of least significant 1 bit
	bits are numbered from 0 to 127 */
	int lsb() const {
		DWORD index;
		if (LO64(*this) != 0) {
			auto rv = _BitScanForward64(&index, LO64(*this));
			return index;
		}
		else if (HI64(*this) != 0) {
			auto rv = _BitScanForward64(&index, HI64(*this));
			return index + 64;
		}
		else return -1;
	}

};

class _uint128 {
public:
	union {
		_S2  s2;
		_S4  s4;
		_S8  s8;
		_S16 s16;
	};


	// constructors
	inline _uint128() {}

	inline _uint128(const _int128 &n) {
		HI64(*this) = HI64(n);
		LO64(*this) = LO64(n);
	}
	inline _uint128(const unsigned __int64 &n) {
		HI64(*this) = 0;
		LO64(*this) = n;
	}
	inline _uint128(const __int64 &n) {
		HI64(*this) = n < 0 ? -1 : 0;
		LO64(*this) = n;
	}
	inline _uint128(unsigned long n) {
		HI64(*this) = 0;
		LO64(*this) = n;
	}
	inline _uint128(long n) {
		HI64(*this) = n < 0 ? -1 : 0;
		LO64(*this) = n;
	}
	inline _uint128(unsigned int n) {
		HI64(*this) = 0;
		LO64(*this) = n;
	}
	inline _uint128(int n) {
		HI64(*this) = n < 0 ? -1 : 0;
		LO64(*this) = n;
	}
	inline _uint128(unsigned short n) {
		HI64(*this) = n < 0 ? -1 : 0;
		LO64(*this) = n;
	}
	inline _uint128(short n) {
		HI64(*this) = n < 0 ? -1 : 0;
		LO64(*this) = n;
	}
	inline _uint128(const unsigned __int64 &hi, const unsigned __int64 &lo) {
		HI64(*this) = hi;
		LO64(*this) = lo;
	}
	
	/* get approximate value in floating point */
	double toDouble() const {
		double l = (double)LO64(*this);
		if (HI64(*this) != 0) {
			double h = (double)HI64(*this);
			return h * pow(2.0, 64) + l;
		}
		else return l;
	}

	// type operators
	inline operator _int128() const {
		return *(_int128*)(void*)this;
	}
	inline operator unsigned __int64() const {
		return LO64(*this);
	}
	inline operator __int64() const {
		return LO64(*this);
	}
	inline operator unsigned long() const {
		return (unsigned long)LO64(*this);
	}
	inline operator long() const {
		return (long)LO64(*this);
	}
	inline operator unsigned int() const {
		return (unsigned int)LO64(*this);
	}
	inline operator int() const {
		return (int)LO64(*this);
	}
	inline operator bool() const {
		return LO64(*this) || HI64(*this);
	}
	inline operator double() const {
		return this->toDouble();
	}

	inline _uint128 &operator++() {   // prefix-form
		int128inc(this);
		return *this;
	}
	inline _uint128 &operator--() {   // prefix-form
		int128dec(this);
		return *this;
	}

	inline _uint128 operator++(int) { // postfix-form
		const _uint128 result(*this);
		int128inc(this);
		return result;
	}
	inline _uint128 operator--(int) { // postfix-form
		const _uint128 result(*this);
		int128dec(this);
		return result;
	}

	inline bool isNegative() const {
		return false;
	}
	inline bool isZero() const {
		return LO64(*this) == 0 && HI64(*this) == 0;
	}

	/* count number of bits set */
	int popcnt() const {
		auto u = (int)__popcnt64(HI64(*this));  // count 1-bits in top half
		auto l = (int)__popcnt64(LO64(*this));  // count 1-bits in bottom half 
		return (int)u + l;    /* return total 1 bits */
	}

	/* get bit number of most significant 1 bit
	bits are numbered from 0 to 127 */
	int msb() const {
		DWORD index;
		if (HI64(*this) != 0) {
			auto rv = _BitScanReverse64(&index, HI64(*this));
			return index + 64;
		}
		else if (LO64(*this) != 0) {
			auto rv = _BitScanReverse64(&index, LO64(*this));
			return index;
		}
		else return -1;  // indicate no bits set
	}

	/* get bit number of least significant 1 bit
	bits are numbered from 0 to 127 */
	int lsb() const {
		DWORD index;
		if (LO64(*this) != 0) {
			auto rv = _BitScanForward64(&index, LO64(*this));
			return index;
		}
		else if (HI64(*this) != 0) {
			auto rv = _BitScanForward64(&index, HI64(*this));
			return index + 64;
		}
		else return -1;  // indicate no bits set
	}
};




// 4 version of all 5 binary arithmetic operators,
// 3 binary logical operators and 6 compare-operators
//    signed   op signed
//    signed   op unsigned
//    unsigned op signed
//    unsigned op unsigned
//  For +, -, *, &, |, ^, ==, != the called function is the same
//  regardless of signed/unsigned combinations.
//  For /, %, <, >, <=, >= however the signed function is used
//  only for the "signed op signed" combination.
//  For left shift (<<) there is no difference for
//  signed and unsigned function, but for right shift (>>)
//  the leftmost bit (bit 127) indicates the sign, and will
//  be copied to all new bits comming in from left for _int128
//  and 0-bits will be shifted in for _uint128 (because there
//  is no sign).
//  For assign-operators (+=,-=...) the same rules apply.
//  Vesions for built in integral types are then defined
//  on top of these

// 4 basic combination of operator+ (128-bit integers - dont care about signed/unsigned)
inline _int128 operator+(const _int128 &lft, const _int128 &rhs) {
	_int128 result(lft);
	int128add(&result, &rhs);
	return result;
}
inline _int128 operator+(const _int128 &lft, const _uint128 &rhs) {
	_int128 result(lft);
	int128add(&result, &rhs);
	return result;
}
inline _uint128 operator+(const _uint128 &lft, const _int128 &rhs) {
	_uint128 result(lft);
	int128add(&result, &rhs);
	return result;
}
inline _uint128 operator+(const _uint128 &lft, const _uint128 &rhs) {
	_uint128 result(lft);
	int128add(&result, &rhs);
	return result;
}

// 4 basic combination of operator- (128-bit integers - dont care about signed/unsigned)
inline _int128 operator-(const _int128 &lft, const _int128 &rhs) {
	_int128 result(lft);
	int128sub(&result, &rhs);
	return result;
}
inline _int128 operator-(const _int128 &lft, const _uint128 &rhs) {
	_int128 result(lft);
	int128sub(&result, &rhs);
	return result;
}
inline _uint128 operator-(const _uint128 &lft, const _int128 &rhs) {
	_uint128 result(lft);
	int128sub(&result, &rhs);
	return result;
}
inline _uint128 operator-(const _uint128 &lft, const _uint128 &rhs) {
	_uint128 result(lft);
	int128sub(&result, &rhs);
	return result;
}

// 4 basic combination of operator* (128-bit integers - dont care about signed/unsigned)
inline _int128 operator*(const _int128 &lft, const _int128 &rhs) {
	_int128 result(lft);
	int128mul(&result, &rhs);
	return result;
}
inline _int128 operator*(const _int128 &lft, const _uint128 &rhs) {
	_int128 result(lft);
	int128mul(&result, &rhs);
	return result;
}
inline _uint128 operator*(const _uint128 &lft, const _int128 &rhs) {
	_uint128 result(lft);
	int128mul(&result, &rhs);
	return result;
}
inline _uint128 operator*(const _uint128 &lft, const _uint128 &rhs) {
	_uint128 result(lft);
	int128mul(&result, &rhs);
	return result;
}

// 4 basic combination of operator/ - signed division only if both are signed
inline _int128 operator/(const _int128 &lft, const _int128 &rhs) {
	_int128 result(lft), tmp(rhs);
	if (tmp.isZero())
		throw std::exception("divide by zero");
	int128div(&result, &tmp);
	return result;
}
inline _int128 operator/(const _int128 &lft, const _uint128 &rhs) {
	_int128 result(lft);
	if (rhs.isZero())
		throw std::exception("divide by zero");
	uint128div(&result, &rhs);
	return result;
}
inline _uint128 operator/(const _uint128 &lft, const _int128 &rhs) {
	_uint128 result(lft);
	if (rhs.isZero())
		throw std::exception("divide by zero");
	uint128div(&result, &rhs);
	return result;
}
inline _uint128 operator/(const _uint128 &lft, const _uint128 &rhs) {
	_uint128 result(lft);
	if (rhs.isZero())
		throw std::exception("divide by zero");
	uint128div(&result, &rhs);
	return result;
}

// 4 basic combination of operator% - signed % only if both are signed
inline _int128 operator%(const _int128 &lft, const _int128 &rhs) {
	_int128 result(lft), tmp(rhs);
	if (rhs.isZero())
		throw std::exception("divide by zero");
	int128rem(&result, &tmp);
	return result;
}
inline _int128 operator%(const _int128 &lft, const _uint128 &rhs) {
	_int128 result(lft);
	if (rhs.isZero())
		throw std::exception("divide by zero");
	uint128rem(&result, &rhs);
	return result;
}
inline _uint128 operator%(const _uint128 &lft, const _int128 &rhs) {
	_uint128 result(lft);
	if (rhs.isZero())
		throw std::exception("divide by zero");
	uint128rem(&result, &rhs);
	return result;
}
inline _uint128 operator%(const _uint128 &lft, const _uint128 &rhs) {
	_uint128 result(lft);
	if (rhs.isZero())
		throw std::exception("divide by zero");
	uint128rem(&result, &rhs);
	return result;
}

// 2 version of unary - (dont care about signed/unsigned)
inline _int128 operator-(const _int128 &x) { // unary minus
	_int128 result(x);
	int128neg(&result);
	return result;
}
inline _uint128 operator-(const _uint128 &x) {
	_uint128 result(x);
	int128neg(&result);
	return result;
}

// Basic bit operators
// 4 basic combinations of operator& (bitwise AND)
inline _int128 operator&(const _int128 &lft, const _int128 &rhs) {
	return _int128(HI64(lft) & HI64(rhs), LO64(lft) & LO64(rhs));
}
inline _int128 operator&(const _int128 &lft, const _uint128 &rhs) {
	return _int128(HI64(lft) & HI64(rhs), LO64(lft) & LO64(rhs));
}
inline _uint128 operator&(const _uint128 &lft, const _int128 &rhs) {
	return _uint128(HI64(lft) & HI64(rhs), LO64(lft) & LO64(rhs));
}
inline _uint128 operator&(const _uint128 &lft, const _uint128 &rhs) {
	return _int128(HI64(lft) & HI64(rhs), LO64(lft) & LO64(rhs));
}

// 4 basic combinations of operator| (inclusive OR)
inline _int128 operator|(const _int128 &lft, const _int128 &rhs) {
	return _int128(HI64(lft) | HI64(rhs), LO64(lft) | LO64(rhs));
}
inline _int128 operator|(const _int128 &lft, const _uint128 &rhs) {
	return _int128(HI64(lft) | HI64(rhs), LO64(lft) | LO64(rhs));
}
inline _uint128 operator|(const _uint128 &lft, const _int128 &rhs) {
	return _uint128(HI64(lft) | HI64(rhs), LO64(lft) | LO64(rhs));
}
inline _uint128 operator|(const _uint128 &lft, const _uint128 &rhs) {
	return _uint128(HI64(lft) | HI64(rhs), LO64(lft) | LO64(rhs));
}

// 4 basic combinations of operator^ (bitwise XOR)
inline _int128 operator^(const _int128 &lft, const _int128 &rhs) {
	return _int128(HI64(lft) ^ HI64(rhs), LO64(lft) ^ LO64(rhs));
}
inline _int128 operator^(const _int128 &lft, const _uint128 &rhs) {
	return _int128(HI64(lft) ^ HI64(rhs), LO64(lft) ^ LO64(rhs));
}
inline _uint128 operator^(const _uint128 &lft, const _int128 &rhs) {
	return _uint128(HI64(lft) ^ HI64(rhs), LO64(lft) ^ LO64(rhs));
}
inline _uint128 operator^(const _uint128 &lft, const _uint128 &rhs) {
	return _uint128(HI64(lft) ^ HI64(rhs), LO64(lft) ^ LO64(rhs));
}

// 2 versions of operator~
inline _int128 operator~(const _int128 &n) {
	return _int128(~HI64(n), ~LO64(n));
}
inline _uint128 operator~(const _uint128 &n) {
	return _uint128(~HI64(n), ~LO64(n));
}

// 2 version of operator>> (arithmetic shift for signed, logical shift for unsigned)
inline _int128 operator>>(const _int128 &lft, const int shft) {
	_int128 copy(lft);
	if (shft != 0)
		int128shr(&copy, shft);
	return copy;
}
inline _uint128 operator>>(const _uint128 &lft, const int shft) {
	_uint128 copy(lft);
	if (shft != 0)
		uint128shr(&copy, shft);
	return copy;
}

// 2 version of operator<< (dont care about signed/unsigned)
inline _int128 operator<<(const _int128 &lft, const int shft) {
	_int128 copy(lft);
	if (shft != 0)
		int128shl(&copy, shft);
	return copy;
}
inline _int128 operator<<(const _uint128 &lft, const int shft) {
	_uint128 copy(lft);
	if (shft != 0)
		int128shl(&copy, shft);
	return copy;
}


// 4 basic combinations of operator==. (dont care about signed/unsigned)
inline bool operator==(const _int128 &lft, const _int128 &rhs) {
	return (LO64(lft) == LO64(rhs)) && (HI64(lft) == HI64(rhs));
}
inline bool operator==(const _int128 &lft, const _uint128 &rhs) {
	return (LO64(lft) == LO64(rhs)) && (HI64(lft) == HI64(rhs));
}
inline bool operator==(const _uint128 &lft, const _int128 &rhs) {
	return (LO64(lft) == LO64(rhs)) && (HI64(lft) == HI64(rhs));
}
inline bool operator==(const _uint128 &lft, const _uint128 &rhs) {
	return (LO64(lft) == LO64(rhs)) && (HI64(lft) == HI64(rhs));
}

// 4 basic combinations of operator!= (dont care about signed/unsigned)
inline bool operator!=(const _int128 &lft, const _int128 &rhs) {
	return (LO64(lft) != LO64(rhs)) || (HI64(lft) != HI64(rhs));
}
inline bool operator!=(const _int128 &lft, const _uint128 &rhs) {
	return (LO64(lft) != LO64(rhs)) || (HI64(lft) != HI64(rhs));
}
inline bool operator!=(const _uint128 &lft, const _int128 &rhs) {
	return (LO64(lft) != LO64(rhs)) || (HI64(lft) != HI64(rhs));
}
inline bool operator!=(const _uint128 &lft, const _uint128 &rhs) {
	return (LO64(lft) != LO64(rhs)) || (HI64(lft) != HI64(rhs));
}

// 4 basic combinations of operator> (signed compare only if both are signed)
inline bool operator>(const _int128 &lft, const _int128 &rhs) {
	return int128cmp(&lft, &rhs) > 0;
}
inline bool operator>(const _int128 &lft, const _uint128 &rhs) {
	return uint128cmp(&lft, &rhs) > 0;
}
inline bool operator>(const _uint128 &lft, const _int128 &rhs) {
	return uint128cmp(&lft, &rhs) > 0;
}
inline bool operator>(const _uint128 &lft, const _uint128 &rhs) {
	return uint128cmp(&lft, &rhs) > 0;
}

// 4 basic combinations of operator>= (signed compare only if both are signed)
inline bool operator>=(const _int128 &lft, const _int128 &rhs) {
	return int128cmp(&lft, &rhs) >= 0;
}
inline bool operator>=(const _int128 &lft, const _uint128 &rhs) {
	return uint128cmp(&lft, &rhs) >= 0;
}
inline bool operator>=(const _uint128 &lft, const _int128 &rhs) {
	return uint128cmp(&lft, &rhs) >= 0;
}
inline bool operator>=(const _uint128 &lft, const _uint128 &rhs) {
	return uint128cmp(&lft, &rhs) >= 0;
}

// 4 basic combinations of operator< (signed compare only if both are signed)
inline bool operator<(const _int128 &lft, const _int128 &rhs) {
	return int128cmp(&lft, &rhs) < 0;
}
inline bool operator<(const _int128 &lft, const _uint128 &rhs) {
	return uint128cmp(&lft, &rhs) < 0;
}
inline bool operator<(const _uint128 &lft, const _int128 &rhs) {
	return uint128cmp(&lft, &rhs) < 0;
}
inline bool operator<(const _uint128 &lft, const _uint128 &rhs) {
	return uint128cmp(&lft, &rhs) < 0;
}

// 4 basic combinations of operator<= (signed compare only if both are signed)
inline bool operator<=(const _int128 &lft, const _int128 &rhs) {
	return int128cmp(&lft, &rhs) <= 0;
}
inline bool operator<=(const _int128 &lft, const _uint128 &rhs) {
	return uint128cmp(&lft, &rhs) <= 0;
}
inline bool operator<=(const _uint128 &lft, const _int128 &rhs) {
	return uint128cmp(&lft, &rhs) <= 0;
}
inline bool operator<=(const _uint128 &lft, const _uint128 &rhs) {
	return uint128cmp(&lft, &rhs) <= 0;
}

// Assign operators
// operator+= (dont care about sign)
inline _int128 &operator+=(_int128 &lft, const _int128 &rhs) {
	int128add(&lft, &rhs);
	return lft;
}
inline _int128 &operator+=(_int128 &lft, const _uint128 &rhs) {
	int128add(&lft, &rhs);
	return lft;
}
inline _uint128 &operator+=(_uint128 &lft, const _int128 &rhs) {
	int128add(&lft, &rhs);
	return lft;
}
inline _uint128 &operator+=(_uint128 &x, const _uint128 &rhs) {
	int128add(&x, &rhs);
	return x;
}

// operator-= (dont care about sign)
inline _int128 &operator-=(_int128 &lft, const _int128 &rhs) {
	int128sub(&lft, &rhs);
	return lft;
}
inline _int128 &operator-=(_int128 &lft, const _uint128 &rhs) {
	int128sub(&lft, &rhs);
	return lft;
}
inline _uint128 &operator-=(_uint128 &lft, const _int128 &rhs) {
	int128sub(&lft, &rhs);
	return lft;
}
inline _uint128 &operator-=(_uint128 &lft, const _uint128 &rhs) {
	int128sub(&lft, &rhs);
	return lft;
}

// operator*= (dont care about sign)
inline _int128 &operator*=(_int128 &lft, const _int128 &rhs) {
	int128mul(&lft, &rhs);
	return lft;
}
inline _int128 &operator*=(_int128 &lft, const _uint128 &rhs) {
	int128mul(&lft, &rhs);
	return lft;
}
inline _uint128 &operator*=(_uint128 &lft, const _int128 &rhs) {
	int128mul(&lft, &rhs);
	return lft;
}
inline _uint128 &operator*=(_uint128 &lft, const _uint128 &rhs) {
	int128mul(&lft, &rhs);
	return lft;
}

// operator/= (use signed div only if both are signed)
inline _int128 &operator/=(_int128 &lft, const _int128 &rhs) {
	_int128 tmp(rhs);
	if (rhs.isZero())
		throw std::exception("divide by zero");
	int128div(&lft, &tmp);
	return lft;
}
inline _int128 &operator/=(_int128 &lft, const _uint128 &rhs) {
	if (rhs.isZero())
		throw std::exception("divide by zero");
	uint128div(&lft, &rhs);
	return lft;
}
inline _uint128 &operator/=(_uint128 &lft, const _int128 &rhs) {
	if (rhs.isZero())
		throw std::exception("divide by zero");
	uint128div(&lft, &rhs);
	return lft;
}
inline _uint128 &operator/=(_uint128 &lft, const _uint128 &rhs) {
	if (rhs.isZero())
		throw std::exception("divide by zero");
	uint128div(&lft, &rhs);
	return lft;
}

// operator%= (use signed % only if both are signed)
inline _int128 &operator%=(_int128 &lft, const _int128 &rhs) {
	_int128 tmp(rhs);
	if (rhs.isZero())
		throw std::exception("divide by zero");
	int128rem(&lft, &tmp);
	return lft;
}
inline _int128 &operator%=(_int128 &lft, const _uint128 &rhs) {
	if (rhs.isZero())
		throw std::exception("divide by zero");
	uint128rem(&lft, &rhs);
	return lft;
}
inline _uint128 &operator%=(_uint128 &lft, const _int128 &rhs) {
	if (rhs.isZero())
		throw std::exception("divide by zero");
	uint128rem(&lft, &rhs);
	return lft;
}
inline _uint128 &operator%=(_uint128 &lft, const _uint128 &rhs) {
	if (rhs.isZero())
		throw std::exception("divide by zero");
	uint128rem(&lft, &rhs);
	return lft;
}

// operator&= (dont care about sign)
inline _int128 &operator&=(_int128 &lft, const _int128 &rhs) {
	LO64(lft) &= LO64(rhs); HI64(lft) &= HI64(rhs);
	return lft;
}
inline _int128 &operator&=(_int128 &lft, const _uint128 &rhs) {
	LO64(lft) &= LO64(rhs); HI64(lft) &= HI64(rhs);
	return lft;
}
inline _uint128 &operator&=(_uint128 &lft, const _int128 &rhs) {
	LO64(lft) &= LO64(rhs); HI64(lft) &= HI64(rhs);
	return lft;
}
inline _uint128 &operator&=(_uint128 &lft, const _uint128 &rhs) {
	LO64(lft) &= LO64(rhs); HI64(lft) &= HI64(rhs);
	return lft;
}

// operator|= (dont care about sign)
inline _int128 &operator|=(_int128 &lft, const _int128 &rhs) {
	LO64(lft) |= LO64(rhs); HI64(lft) |= HI64(rhs);
	return lft;
}
inline _int128 &operator|=(_int128 &lft, const _uint128 &rhs) {
	LO64(lft) |= LO64(rhs); HI64(lft) |= HI64(rhs);
	return lft;
}
inline _uint128 &operator|=(_uint128 &lft, const _int128 &rhs) {
	LO64(lft) |= LO64(rhs); HI64(lft) |= HI64(rhs);
	return lft;
}
inline _uint128 &operator|=(_uint128 &lft, const _uint128 &rhs) {
	LO64(lft) |= LO64(rhs); HI64(lft) |= HI64(rhs);
	return lft;
}

// operator^= (dont care about sign)
inline _int128 &operator^=(_int128 &lft, const _int128 &rhs) {
	LO64(lft) ^= LO64(rhs); HI64(lft) ^= HI64(rhs);
	return lft;
}
inline _int128 &operator^=(_int128 &lft, const _uint128 &rhs) {
	LO64(lft) ^= LO64(rhs); HI64(lft) ^= HI64(rhs);
	return lft;
}
inline _uint128 &operator^=(_uint128 &lft, const _int128 &rhs) {
	LO64(lft) ^= LO64(rhs); HI64(lft) ^= HI64(rhs);
	return lft;
}
inline _uint128 &operator^=(_uint128 &lft, const _uint128 &rhs) {
	LO64(lft) ^= LO64(rhs); HI64(lft) ^= HI64(rhs);
	return lft;
}

inline _int128 &operator>>=(_int128 &lft, int shft) {
	if (shft != 0)
		int128shr(&lft, shft);
	return lft;
}
inline _uint128 &operator>>=(_uint128 &lft, const int shft) {
	if (shft != 0) {
		uint128shr(&lft, shft);
	}
	return lft;
}
inline _int128 &operator<<=(_int128 &lft, int shft) {
	if (shft != 0)
		int128shl(&lft, shft);
	return lft;
}
inline _uint128 &operator<<=(_uint128 &lft, int shft) {
	if (shft != 0)
		int128shl(&lft, shft);
	return lft;
}

// Now all combinations of binary operators for lft = 128-bit and rhs is  
// integral type
// operator + for built in integral types as second argument
inline _int128  operator+(const _int128  &lft, __int64 rhs) {
	return lft + (_int128)rhs;
}
inline _int128  operator+(const _int128  &lft, unsigned __int64 rhs) {
	return lft + (_uint128)rhs;
}
inline _int128  operator+(const _int128  &lft, long rhs) {
	return lft + (_int128)rhs;
}
inline _int128  operator+(const _int128  &lft, unsigned long rhs) {
	return lft + (_uint128)rhs;
}
inline _int128  operator+(const _int128  &lft, int rhs) {
	return lft + (_int128)rhs;
}
inline _int128  operator+(const _int128  &lft, unsigned int rhs) {
	return lft + (_uint128)rhs;
}
inline _int128  operator+(const _int128  &lft, short rhs) {
	return lft + (_int128)rhs;
}
inline _int128  operator+(const _int128  &lft, unsigned short rhs) {
	return lft + (_uint128)rhs;
}

inline _uint128 operator+(const _uint128 &lft, __int64 rhs) {
	return lft + (_int128)rhs;
}
inline _uint128 operator+(const _uint128 &lft, unsigned __int64 rhs) {
	//return lft + (_uint128)rhs;
	_uint128 sum;
	auto c = _addcarry_u64(0, LO64(lft), rhs, &LO64(sum));
	HI64(sum) = HI64(lft) + c;
	return sum;
}
inline _uint128 operator+(const _uint128 &lft, long rhs) {
	return lft + (_int128)rhs;
}
inline _uint128 operator+(const _uint128 &lft, unsigned long rhs) {
	//return lft + (_uint128)rhs;
	_uint128 sum;
	auto c = _addcarry_u64(0, LO64(lft), rhs, &LO64(sum));
	HI64(sum) = HI64(lft) + c;
	return sum;
}
inline _uint128 operator+(const _uint128 &lft, int rhs) {
	return lft + (_int128)rhs;
}
inline _uint128 operator+(const _uint128 &lft, unsigned int rhs) {
	//return lft + (_uint128)rhs;
	_uint128 sum;
	auto c = _addcarry_u64(0, LO64(lft), rhs, &LO64(sum));
	HI64(sum) = HI64(lft) + c;
	return sum;
}
inline _uint128 operator+(const _uint128 &lft, short rhs) {
	return lft + (_int128)rhs;
}
inline _uint128 operator+(const _uint128 &lft, unsigned short rhs) {
	//return lft + (_uint128)rhs;
	_uint128 sum;
	auto c = _addcarry_u64(0, LO64(lft), rhs, &LO64(sum));
	HI64(sum) = HI64(lft) + c;
	return sum;
}


// operator - for built in integral types as second argument
inline _int128  operator-(const _int128  &lft, __int64 rhs) {
	return lft - (_int128)rhs;
}
inline _int128  operator-(const _int128  &lft, unsigned __int64 rhs) {
	return lft - (_uint128)rhs;
}
inline _int128  operator-(const _int128  &lft, long rhs) {
	return lft - (_int128)rhs;
}
inline _int128  operator-(const _int128  &lft, unsigned long rhs) {
	return lft - (_uint128)rhs;
}
inline _int128  operator-(const _int128  &lft, int rhs) {
	return lft - (_int128)rhs;
}
inline _int128  operator-(const _int128  &lft, unsigned int rhs) {
	return lft - (_uint128)rhs;
}
inline _int128  operator-(const _int128  &lft, short rhs) {
	return lft - (_int128)rhs;
}
inline _int128  operator-(const _int128  &lft, unsigned short rhs) {
	return lft - (_uint128)rhs;
}

inline _uint128 operator-(const _uint128 &lft, __int64 rhs) {
	return lft - (_int128)rhs;
}
inline _uint128 operator-(const _uint128 &lft, unsigned __int64 rhs) {
	return lft - (_uint128)rhs;
}
inline _uint128 operator-(const _uint128 &lft, long rhs) {
	return lft - (_int128)rhs;
}
inline _uint128 operator-(const _uint128 &lft, unsigned long rhs) {
	return lft - (_uint128)rhs;
}
inline _uint128 operator-(const _uint128 &lft, int rhs) {
	return lft - (_int128)rhs;
}
inline _uint128 operator-(const _uint128 &lft, unsigned int rhs) {
	return lft - (_uint128)rhs;
}
inline _uint128 operator-(const _uint128 &lft, short rhs) {
	return lft - (_int128)rhs;
}
inline _uint128 operator-(const _uint128 &lft, unsigned short rhs) {
	return lft - (_uint128)rhs;
}


// operator * for built in integral types as second argument
inline _int128  operator*(const _int128  &lft, __int64 rhs) {
	return lft * (_int128)rhs;
}
inline _int128  operator*(const _int128  &lft, unsigned __int64 rhs) {
	return lft * (_uint128)rhs;
}
inline _int128  operator*(const _int128  &lft, long rhs) {
	return lft * (_int128)rhs;
}
inline _int128  operator*(const _int128  &lft, unsigned long rhs) {
	return lft * (_uint128)rhs;
}
inline _int128  operator*(const _int128  &lft, int rhs) {
	return lft * (_int128)rhs;
}
inline _int128  operator*(const _int128  &lft, unsigned int rhs) {
	return lft * (_uint128)rhs;
}
inline _int128  operator*(const _int128  &lft, short rhs) {
	return lft * (_int128)rhs;
}
inline _int128  operator*(const _int128  &lft, unsigned short rhs) {
	return lft * (_uint128)rhs;
}

inline _uint128 operator*(const _uint128 &lft, __int64 rhs) {
	return lft * (_int128)rhs;
}
inline _uint128 operator*(const _uint128 &lft, unsigned __int64 rhs) {
	return lft * (_uint128)rhs;
}
inline _uint128 operator*(const _uint128 &lft, long rhs) {
	return lft * (_int128)rhs;
}
inline _uint128 operator*(const _uint128 &lft, unsigned long rhs) {
	return lft * (_uint128)rhs;
}
inline _uint128 operator*(const _uint128 &lft, int rhs) {
	return lft * (_int128)rhs;
}
inline _uint128 operator*(const _uint128 &lft, unsigned int rhs) {
	return lft * (_uint128)rhs;
}
inline _uint128 operator*(const _uint128 &lft, short rhs) {
	return lft * (_int128)rhs;
}
inline _uint128 operator*(const _uint128 &lft, unsigned short rhs) {
	return lft * (_uint128)rhs;
}


// operator / for built in integral types as second argument
inline _int128  operator/(const _int128  &lft, __int64 rhs) {
	return lft / (_int128)rhs;
}
inline _int128  operator/(const _int128  &lft, unsigned __int64 rhs) {
	return lft / (_int128)rhs;
}
inline _int128  operator/(const _int128  &lft, long rhs) {
	return lft / (_int128)rhs;
}
inline _int128  operator/(const _int128  &lft, unsigned long rhs) {
	return lft / (_int128)rhs;
}
inline _int128  operator/(const _int128  &lft, int rhs) {
	return lft / (_int128)rhs;
}
inline _int128  operator/(const _int128  &lft, unsigned int rhs) {
	return lft / (_int128)rhs;
}
inline _int128  operator/(const _int128  &lft, short rhs) {
	return lft / (_int128)rhs;
}
inline _int128  operator/(const _int128  &lft, unsigned short rhs) {
	return lft / (_int128)rhs;
}

inline _uint128 operator/(const _uint128 &lft, __int64 rhs) {
	return lft / (_int128)rhs;
}
inline _uint128 operator/(const _uint128 &lft, unsigned __int64 rhs) {
	return lft / (_uint128)rhs;
}
inline _uint128 operator/(const _uint128 &lft, long rhs) {
	return lft / (_int128)rhs;
}
inline _uint128 operator/(const _uint128 &lft, unsigned long rhs) {
	return lft / (_uint128)rhs;
}
inline _uint128 operator/(const _uint128 &lft, int rhs) {
	return lft / (_int128)rhs;
}
inline _uint128 operator/(const _uint128 &lft, unsigned int rhs) {
	return lft / (_uint128)rhs;
}
inline _uint128 operator/(const _uint128 &lft, short rhs) {
	return lft / (_int128)rhs;
}
inline _uint128 operator/(const _uint128 &lft, unsigned short rhs) {
	return lft / (_uint128)rhs;
}


// operator % for built in integral types as second argument
inline _int128  operator%(const _int128  &lft, __int64 rhs) {
	return lft % (_int128)rhs;
}
inline _int128  operator%(const _int128  &lft, unsigned __int64 rhs) {
	return lft % (_int128)rhs;
}
inline _int128  operator%(const _int128  &lft, long rhs) {
	return lft % (_int128)rhs;
}
inline _int128  operator%(const _int128  &lft, unsigned long rhs) {
	return lft % (_int128)rhs;
}
inline _int128  operator%(const _int128  &lft, int rhs) {
	return lft % (_int128)rhs;
}
inline _int128  operator%(const _int128  &lft, unsigned int rhs) {
	return lft % (_int128)rhs;
}
inline _int128  operator%(const _int128  &lft, short rhs) {
	return lft % (_int128)rhs;
}
inline _int128  operator%(const _int128  &lft, unsigned short rhs) {
	return lft % (_int128)rhs;
}

inline _uint128 operator%(const _uint128 &lft, __int64 rhs) {
	return lft % (_int128)rhs;
}
inline _uint128 operator%(const _uint128 &lft, unsigned __int64 rhs) {
	return lft % (_uint128)rhs;
}
inline _uint128 operator%(const _uint128 &lft, long rhs) {
	return lft % (_int128)rhs;
}
inline _uint128 operator%(const _uint128 &lft, unsigned long rhs) {
	return lft % (_uint128)rhs;
}
inline _uint128 operator%(const _uint128 &lft, int rhs) {
	return lft % (_int128)rhs;
}
inline _uint128 operator%(const _uint128 &lft, unsigned int rhs) {
	return lft % (_uint128)rhs;
}
inline _uint128 operator%(const _uint128 &lft, short rhs) {
	return lft % (_int128)rhs;
}
inline _uint128 operator%(const _uint128 &lft, unsigned short rhs) {
	return lft % (_uint128)rhs;
}


// operator & for built in integral types as second argument
inline _int128  operator&(const _int128  &lft, __int64 rhs) {
	return lft & (_int128)rhs;
}
inline _int128  operator&(const _int128  &lft, unsigned __int64 rhs) {
	return lft & (_int128)rhs;
}
inline _int128  operator&(const _int128  &lft, long rhs) {
	return lft & (_int128)rhs;
}
inline _int128  operator&(const _int128  &lft, unsigned long rhs) {
	return lft & (_int128)rhs;
}
inline _int128  operator&(const _int128  &lft, int rhs) {
	return lft & (_int128)rhs;
}
inline _int128  operator&(const _int128  &lft, unsigned int rhs) {
	return lft & (_int128)rhs;
}
inline _int128  operator&(const _int128  &lft, short rhs) {
	return lft & (_int128)rhs;
}
inline _int128  operator&(const _int128  &lft, unsigned short rhs) {
	return lft & (_int128)rhs;
}

inline _uint128 operator&(const _uint128 &lft, __int64 rhs) {
	return lft & (_int128)rhs;
}
inline _uint128 operator&(const _uint128 &lft, unsigned __int64 rhs) {
	return lft & (_uint128)rhs;
}
inline _uint128 operator&(const _uint128 &lft, long rhs) {
	return lft & (_int128)rhs;
}
inline _uint128 operator&(const _uint128 &lft, unsigned long rhs) {
	return lft & (_uint128)rhs;
}
inline _uint128 operator&(const _uint128 &lft, int rhs) {
	return lft & (_int128)rhs;
}
inline _uint128 operator&(const _uint128 &lft, unsigned int rhs) {
	return lft & (_uint128)rhs;
}
inline _uint128 operator&(const _uint128 &lft, short rhs) {
	return lft & (_int128)rhs;
}
inline _uint128 operator&(const _uint128 &lft, unsigned short rhs) {
	return lft & (_uint128)rhs;
}


// operator | for built in integral types as second argument
inline _int128  operator|(const _int128  &lft, __int64 rhs) {
	return lft | (_int128)rhs;
}
inline _int128  operator|(const _int128  &lft, unsigned __int64 rhs) {
	return lft | (_int128)rhs;
}
inline _int128  operator|(const _int128  &lft, long rhs) {
	return lft | (_int128)rhs;
}
inline _int128  operator|(const _int128  &lft, unsigned long rhs) {
	return lft | (_int128)rhs;
}
inline _int128  operator|(const _int128  &lft, int rhs) {
	return lft | (_int128)rhs;
}
inline _int128  operator|(const _int128  &lft, unsigned int rhs) {
	return lft | (_int128)rhs;
}
inline _int128  operator|(const _int128  &lft, short rhs) {
	return lft | (_int128)rhs;
}
inline _int128  operator|(const _int128  &lft, unsigned short rhs) {
	return lft | (_int128)rhs;
}

inline _uint128 operator|(const _uint128 &lft, __int64 rhs) {
	return lft | (_int128)rhs;
}
inline _uint128 operator|(const _uint128 &lft, unsigned __int64 rhs) {
	return lft | (_uint128)rhs;
}
inline _uint128 operator|(const _uint128 &lft, long rhs) {
	return lft | (_int128)rhs;
}
inline _uint128 operator|(const _uint128 &lft, unsigned long rhs) {
	return lft | (_uint128)rhs;
}
inline _uint128 operator|(const _uint128 &lft, int rhs) {
	return lft | (_int128)rhs;
}
inline _uint128 operator|(const _uint128 &lft, unsigned int rhs) {
	return lft | (_uint128)rhs;
}
inline _uint128 operator|(const _uint128 &lft, short rhs) {
	return lft | (_int128)rhs;
}
inline _uint128 operator|(const _uint128 &lft, unsigned short rhs) {
	return lft | (_uint128)rhs;
}


// operator ^ for built in integral types as second argument
inline _int128  operator^(const _int128  &lft, __int64 rhs) {
	return lft ^ (_int128)rhs;
}
inline _int128  operator^(const _int128  &lft, unsigned __int64 rhs) {
	return lft ^ (_int128)rhs;
}
inline _int128  operator^(const _int128  &lft, long rhs) {
	return lft ^ (_int128)rhs;
}
inline _int128  operator^(const _int128  &lft, unsigned long rhs) {
	return lft ^ (_int128)rhs;
}
inline _int128  operator^(const _int128  &lft, int rhs) {
	return lft ^ (_int128)rhs;
}
inline _int128  operator^(const _int128  &lft, unsigned int rhs) {
	return lft ^ (_int128)rhs;
}
inline _int128  operator^(const _int128  &lft, short rhs) {
	return lft ^ (_int128)rhs;
}
inline _int128  operator^(const _int128  &lft, unsigned short rhs) {
	return lft ^ (_int128)rhs;
}

inline _uint128 operator^(const _uint128 &lft, __int64 rhs) {
	return lft ^ (_int128)rhs;
}
inline _uint128 operator^(const _uint128 &lft, unsigned __int64 rhs) {
	return lft ^ (_uint128)rhs;
}
inline _uint128 operator^(const _uint128 &lft, long rhs) {
	return lft ^ (_int128)rhs;
}
inline _uint128 operator^(const _uint128 &lft, unsigned long rhs) {
	return lft ^ (_uint128)rhs;
}
inline _uint128 operator^(const _uint128 &lft, int rhs) {
	return lft ^ (_int128)rhs;
}
inline _uint128 operator^(const _uint128 &lft, unsigned int rhs) {
	return lft ^ (_uint128)rhs;
}
inline _uint128 operator^(const _uint128 &lft, short rhs) {
	return lft ^ (_int128)rhs;
}
inline _uint128 operator^(const _uint128 &lft, unsigned short rhs) {
	return lft ^ (_uint128)rhs;
}

// Compare operators where second argument is built in integral type
// operator == for built in integral types as second argument
inline bool operator==(const _int128 &lft, __int64 rhs) {
	return lft == _int128(rhs);
}
inline bool operator==(const _int128 &lft, unsigned __int64 rhs) {
	return lft == _int128(rhs);
}
inline bool operator==(const _int128 &lft, long rhs) {
	return lft == _int128(rhs);
}
inline bool operator==(const _int128 &lft, unsigned long rhs) {
	return lft == _int128(rhs);
}
inline bool operator==(const _int128 &lft, int rhs) {
	return lft == _int128(rhs);
}
inline bool operator==(const _int128 &lft, unsigned int rhs) {
	return lft == _int128(rhs);
}
inline bool operator==(const _int128 &lft, short rhs) {
	return lft == _int128(rhs);
}
inline bool operator==(const _int128 &lft, unsigned short rhs) {
	return lft == _int128(rhs);
}

inline bool operator==(const _uint128 &lft, __int64 rhs) {
	return lft == _int128(rhs);
}
inline bool operator==(const _uint128 &lft, unsigned __int64 rhs) {
	return lft == _uint128(rhs);
}
inline bool operator==(const _uint128 &lft, long rhs) {
	return lft == _int128(rhs);
}
inline bool operator==(const _uint128 &lft, unsigned long rhs) {
	return lft == _uint128(rhs);
}
inline bool operator==(const _uint128 &lft, int rhs) {
	return lft == _int128(rhs);
}
inline bool operator==(const _uint128 &lft, unsigned int rhs) {
	return lft == _uint128(rhs);
}
inline bool operator==(const _uint128 &lft, short rhs) {
	return lft == _int128(rhs);
}
inline bool operator==(const _uint128 &lft, unsigned short rhs) {
	return lft == _uint128(rhs);
}


// operator != for built in integral types as second argument
inline bool operator!=(const _int128 &lft, __int64 rhs) {
	return lft != _int128(rhs);
}
inline bool operator!=(const _int128 &lft, unsigned __int64 rhs) {
	return lft != _int128(rhs);
}
inline bool operator!=(const _int128 &lft, long rhs) {
	return lft != _int128(rhs);
}
inline bool operator!=(const _int128 &lft, unsigned long rhs) {
	return lft != _int128(rhs);
}
inline bool operator!=(const _int128 &lft, int rhs) {
	return lft != _int128(rhs);
}
inline bool operator!=(const _int128 &lft, unsigned int rhs) {
	return lft != _int128(rhs);
}
inline bool operator!=(const _int128 &lft, short rhs) {
	return lft != _int128(rhs);
}
inline bool operator!=(const _int128 &lft, unsigned short rhs) {
	return lft != _int128(rhs);
}

inline bool operator!=(const _uint128 &lft, __int64 rhs) {
	return lft != _int128(rhs);
}
inline bool operator!=(const _uint128 &lft, unsigned __int64 rhs) {
	return lft != _uint128(rhs);
}
inline bool operator!=(const _uint128 &lft, long rhs) {
	return lft != _int128(rhs);
}
inline bool operator!=(const _uint128 &lft, unsigned long rhs) {
	return lft != _uint128(rhs);
}
inline bool operator!=(const _uint128 &lft, int rhs) {
	return lft != _int128(rhs);
}
inline bool operator!=(const _uint128 &lft, unsigned int rhs) {
	return lft != _uint128(rhs);
}
inline bool operator!=(const _uint128 &lft, short rhs) {
	return lft != _int128(rhs);
}
inline bool operator!=(const _uint128 &lft, unsigned short rhs) {
	return lft != _uint128(rhs);
}


// operator > for built in integral types as second argument
inline bool operator>(const _int128 &lft, __int64 rhs) {
	return lft > _int128(rhs);
}
inline bool operator>(const _int128 &lft, unsigned __int64 rhs) {
	return lft > _uint128(rhs);
}
inline bool operator>(const _int128 &lft, long rhs) {
	return lft > _int128(rhs);
}
inline bool operator>(const _int128 &lft, unsigned long rhs) {
	return lft > _uint128(rhs);
}
inline bool operator>(const _int128 &lft, int rhs) {
	return lft > _int128(rhs);
}
inline bool operator>(const _int128 &lft, unsigned int rhs) {
	return lft > _uint128(rhs);
}
inline bool operator>(const _int128 &lft, short rhs) {
	return lft > _int128(rhs);
}
inline bool operator>(const _int128 &lft, unsigned short rhs) {
	return lft > _uint128(rhs);
}

inline bool operator>(const _uint128 &lft, __int64 rhs) {
	return lft > _int128(rhs);
}
inline bool operator>(const _uint128 &lft, unsigned __int64 rhs) {
	return lft > _uint128(rhs);
}
inline bool operator>(const _uint128 &lft, long rhs) {
	return lft > _int128(rhs);
}
inline bool operator>(const _uint128 &lft, unsigned long rhs) {
	return lft > _uint128(rhs);
}
inline bool operator>(const _uint128 &lft, int rhs) {
	return lft > _int128(rhs);
}
inline bool operator>(const _uint128 &lft, unsigned int rhs) {
	return lft > _uint128(rhs);
}
inline bool operator>(const _uint128 &lft, short rhs) {
	return lft > _int128(rhs);
}
inline bool operator>(const _uint128 &lft, unsigned short rhs) {
	return lft > _uint128(rhs);
}


// operator >= for built in integral types as second argument
inline bool operator>=(const _int128 &lft, __int64 rhs) {
	return lft >= _int128(rhs);
}
inline bool operator>=(const _int128 &lft, unsigned __int64 rhs) {
	return lft >= _uint128(rhs);
}
inline bool operator>=(const _int128 &lft, long rhs) {
	return lft >= _int128(rhs);
}
inline bool operator>=(const _int128 &lft, unsigned long rhs) {
	return lft >= _uint128(rhs);
}
inline bool operator>=(const _int128 &lft, int rhs) {
	return lft >= _int128(rhs);
}
inline bool operator>=(const _int128 &lft, unsigned int rhs) {
	return lft >= _uint128(rhs);
}
inline bool operator>=(const _int128 &lft, short rhs) {
	return lft >= _int128(rhs);
}
inline bool operator>=(const _int128 &lft, unsigned short rhs) {
	return lft >= _uint128(rhs);
}

inline bool operator>=(const _uint128 &lft, __int64 rhs) {
	return lft >= _int128(rhs);
}
inline bool operator>=(const _uint128 &lft, unsigned __int64 rhs) {
	return lft >= _uint128(rhs);
}
inline bool operator>=(const _uint128 &lft, long rhs) {
	return lft >= _int128(rhs);
}
inline bool operator>=(const _uint128 &lft, unsigned long rhs) {
	return lft >= _uint128(rhs);
}
inline bool operator>=(const _uint128 &lft, int rhs) {
	return lft >= _int128(rhs);
}
inline bool operator>=(const _uint128 &lft, unsigned int rhs) {
	return lft >= _uint128(rhs);
}
inline bool operator>=(const _uint128 &lft, short rhs) {
	return lft >= _int128(rhs);
}
inline bool operator>=(const _uint128 &lft, unsigned short rhs) {
	return lft >= _uint128(rhs);
}


// operator < for built in integral types as second argument
inline bool operator<(const _int128 &lft, __int64 rhs) {
	return lft < _int128(rhs);
}
inline bool operator<(const _int128 &lft, unsigned __int64 rhs) {
	return lft < _uint128(rhs);
}
inline bool operator<(const _int128 &lft, long rhs) {
	return lft < _int128(rhs);
}
inline bool operator<(const _int128 &lft, unsigned long rhs) {
	return lft < _uint128(rhs);
}
inline bool operator<(const _int128 &lft, int rhs) {
	return lft < _int128(rhs);
}
inline bool operator<(const _int128 &lft, unsigned int rhs) {
	return lft < _uint128(rhs);
}
inline bool operator<(const _int128 &lft, short rhs) {
	return lft < _int128(rhs);
}
inline bool operator<(const _int128 &lft, unsigned short rhs) {
	return lft < _uint128(rhs);
}

inline bool operator<(const _uint128 &lft, __int64 rhs) {
	return lft < _int128(rhs);
}
inline bool operator<(const _uint128 &lft, unsigned __int64 rhs) {
	return lft < _uint128(rhs);
}
inline bool operator<(const _uint128 &lft, long rhs) {
	return lft < _int128(rhs);
}
inline bool operator<(const _uint128 &lft, unsigned long rhs) {
	return lft < _uint128(rhs);
}
inline bool operator<(const _uint128 &lft, int rhs) {
	return lft < _int128(rhs);
}
inline bool operator<(const _uint128 &lft, unsigned int rhs) {
	return lft < _uint128(rhs);
}
inline bool operator<(const _uint128 &lft, short rhs) {
	return lft < _int128(rhs);
}
inline bool operator<(const _uint128 &lft, unsigned short rhs) {
	return lft < _uint128(rhs);
}


// operator <= for built in integral types as second argument
inline bool operator<=(const _int128 &lft, __int64 rhs) {
	return lft <= _int128(rhs);
}
inline bool operator<=(const _int128 &lft, unsigned __int64 rhs) {
	return lft <= _uint128(rhs);
}
inline bool operator<=(const _int128 &lft, long rhs) {
	return lft <= _int128(rhs);
}
inline bool operator<=(const _int128 &lft, unsigned long rhs) {
	return lft <= _uint128(rhs);
}
inline bool operator<=(const _int128 &lft, int rhs) {
	return lft <= _int128(rhs);
}
inline bool operator<=(const _int128 &lft, unsigned int rhs) {
	return lft <= _uint128(rhs);
}
inline bool operator<=(const _int128 &lft, short rhs) {
	return lft <= _int128(rhs);
}
inline bool operator<=(const _int128 &lft, unsigned short rhs) {
	return lft <= _uint128(rhs);
}

inline bool operator<=(const _uint128 &lft, __int64 rhs) {
	return lft <= _int128(rhs);
}
inline bool operator<=(const _uint128 &lft, unsigned __int64 rhs) {
	return lft <= _uint128(rhs);
}
inline bool operator<=(const _uint128 &lft, long rhs) {
	return lft <= _int128(rhs);
}
inline bool operator<=(const _uint128 &lft, unsigned long rhs) {
	return lft <= _uint128(rhs);
}
inline bool operator<=(const _uint128 &lft, int rhs) {
	return lft <= _int128(rhs);
}
inline bool operator<=(const _uint128 &lft, unsigned int rhs) {
	return lft <= _uint128(rhs);
}
inline bool operator<=(const _uint128 &lft, short rhs) {
	return lft <= _int128(rhs);
}
inline bool operator<=(const _uint128 &lft, unsigned short rhs) {
	return lft <= _uint128(rhs);
}

// Assign operators where second argument is built in integral type
// operator += for built in integral types as second argument
inline _int128  &operator+=(_int128  &lft, __int64 rhs) {
	return lft += (_int128)rhs;
}
inline _int128  &operator+=(_int128  &lft, unsigned __int64 rhs) {
	return lft += (_uint128)rhs;
}
inline _int128  &operator+=(_int128  &lft, long rhs) {
	return lft += (_int128)rhs;
}
inline _int128  &operator+=(_int128  &lft, unsigned long rhs) {
	return lft += (_uint128)rhs;
}
inline _int128  &operator+=(_int128  &lft, int rhs) {
	return lft += (_int128)rhs;
}
inline _int128  &operator+=(_int128  &lft, unsigned int rhs) {
	return lft += (_uint128)rhs;
}
inline _int128  &operator+=(_int128  &lft, short rhs) {
	return lft += (_int128)rhs;
}
inline _int128  &operator+=(_int128  &lft, unsigned short rhs) {
	return lft += (_uint128)rhs;
}

inline _uint128 &operator+=(_uint128 &lft, __int64 rhs) {
	return lft += (_int128)rhs;
}
inline _uint128 &operator+=(_uint128 &lft, unsigned __int64 rhs) {
	//return lft += (_uint128)rhs;
	auto c = _addcarry_u64(0, LO64(lft), rhs, &LO64(lft));
	HI64(lft) += c;
	return lft;
}
inline _uint128 &operator+=(_uint128 &lft, long rhs) {
	return lft += (_int128)rhs;
}
inline _uint128 &operator+=(_uint128 &lft, unsigned long rhs) {
	return lft += (_uint128)rhs;
}
inline _uint128 &operator+=(_uint128 &lft, int rhs) {
	return lft += (_int128)rhs;
}
inline _uint128 &operator+=(_uint128 &lft, unsigned int rhs) {
	return lft += (_uint128)rhs;
}
inline _uint128 &operator+=(_uint128 &lft, short rhs) {
	return lft += (_int128)rhs;
}
inline _uint128 &operator+=(_uint128 &lft, unsigned short rhs) {
	return lft += (_uint128)rhs;
}


// operator -= for built in integral types as second argument
inline _int128  &operator-=(_int128  &lft, __int64 rhs) {
	return lft -= (_int128)rhs;
}
inline _int128  &operator-=(_int128  &lft, unsigned __int64 rhs) {
	return lft -= (_uint128)rhs;
}
inline _int128  &operator-=(_int128  &lft, long rhs) {
	return lft -= (_int128)rhs;
}
inline _int128  &operator-=(_int128  &lft, unsigned long rhs) {
	return lft -= (_uint128)rhs;
}
inline _int128  &operator-=(_int128  &lft, int rhs) {
	return lft -= (_int128)rhs;
}
inline _int128  &operator-=(_int128  &lft, unsigned int rhs) {
	return lft -= (_uint128)rhs;
}
inline _int128  &operator-=(_int128  &lft, short rhs) {
	return lft -= (_int128)rhs;
}
inline _int128  &operator-=(_int128  &lft, unsigned short rhs) {
	return lft -= (_uint128)rhs;
}

inline _uint128 &operator-=(_uint128 &lft, __int64 rhs) {
	return lft -= (_int128)rhs;
}
inline _uint128 &operator-=(_uint128 &lft, unsigned __int64 rhs) {
	return lft -= (_uint128)rhs;
}
inline _uint128 &operator-=(_uint128 &lft, long rhs) {
	return lft -= (_int128)rhs;
}
inline _uint128 &operator-=(_uint128 &lft, unsigned long rhs) {
	return lft -= (_uint128)rhs;
}
inline _uint128 &operator-=(_uint128 &lft, int rhs) {
	return lft -= (_int128)rhs;
}
inline _uint128 &operator-=(_uint128 &lft, unsigned int rhs) {
	return lft -= (_uint128)rhs;
}
inline _uint128 &operator-=(_uint128 &lft, short rhs) {
	return lft -= (_int128)rhs;
}
inline _uint128 &operator-=(_uint128 &lft, unsigned short rhs) {
	return lft -= (_uint128)rhs;
}


// operator*= for built in integral types as second argument
inline _int128  &operator*=(_int128  &lft, __int64 rhs) {
	return lft *= (_int128)rhs;
}
inline _int128  &operator*=(_int128  &lft, unsigned __int64 rhs) {
	return lft *= (_uint128)rhs;
}
inline _int128  &operator*=(_int128  &lft, long rhs) {
	return lft *= (_int128)rhs;
}
inline _int128  &operator*=(_int128  &lft, unsigned long rhs) {
	return lft *= (_uint128)rhs;
}
inline _int128  &operator*=(_int128  &lft, int rhs) {
	return lft *= (_int128)rhs;
}
inline _int128  &operator*=(_int128  &lft, unsigned int rhs) {
	return lft *= (_uint128)rhs;
}
inline _int128  &operator*=(_int128  &lft, short rhs) {
	return lft *= (_int128)rhs;
}
inline _int128  &operator*=(_int128  &lft, unsigned short rhs) {
	return lft *= (_uint128)rhs;
}

inline _uint128 &operator*=(_uint128 &lft, __int64 rhs) {
	return lft *= (_int128)rhs;
}
inline _uint128 &operator*=(_uint128 &lft, unsigned __int64 rhs) {
	return lft *= (_uint128)rhs;
}
inline _uint128 &operator*=(_uint128 &lft, long rhs) {
	return lft *= (_int128)rhs;
}
inline _uint128 &operator*=(_uint128 &lft, unsigned long rhs) {
	return lft *= (_uint128)rhs;
}
inline _uint128 &operator*=(_uint128 &lft, int rhs) {
	return lft *= (_int128)rhs;
}
inline _uint128 &operator*=(_uint128 &lft, unsigned int rhs) {
	return lft *= (_uint128)rhs;
}
inline _uint128 &operator*=(_uint128 &lft, short rhs) {
	return lft *= (_int128)rhs;
}
inline _uint128 &operator*=(_uint128 &lft, unsigned short rhs) {
	return lft *= (_uint128)rhs;
}


// operator /= for built in integral types as second argument
inline _int128  &operator/=(_int128  &lft, __int64 rhs) {
	return lft /= (_int128)rhs;
}
inline _int128  &operator/=(_int128  &lft, unsigned __int64 rhs) {
	return lft /= (_uint128)rhs;
}
inline _int128  &operator/=(_int128  &lft, long rhs) {
	return lft /= (_int128)rhs;
}
inline _int128  &operator/=(_int128  &lft, unsigned long rhs) {
	return lft /= (_uint128)rhs;
}
inline _int128  &operator/=(_int128  &lft, int rhs) {
	return lft /= (_int128)rhs;
}
inline _int128  &operator/=(_int128  &lft, unsigned int rhs) {
	return lft /= (_uint128)rhs;
}
inline _int128  &operator/=(_int128  &lft, short rhs) {
	return lft /= (_int128)rhs;
}
inline _int128  &operator/=(_int128  &lft, unsigned short rhs) {
	return lft /= (_uint128)rhs;
}

inline _uint128 &operator/=(_uint128 &lft, __int64 rhs) {
	return lft /= (_int128)rhs;
}
inline _uint128 &operator/=(_uint128 &lft, unsigned __int64 rhs) {
	return lft /= (_uint128)rhs;
}
inline _uint128 &operator/=(_uint128 &lft, long rhs) {
	return lft /= (_int128)rhs;
}
inline _uint128 &operator/=(_uint128 &lft, unsigned long rhs) {
	return lft /= (_uint128)rhs;
}
inline _uint128 &operator/=(_uint128 &lft, int rhs) {
	return lft /= (_int128)rhs;
}
inline _uint128 &operator/=(_uint128 &lft, unsigned int rhs) {
	return lft /= (_uint128)rhs;
}
inline _uint128 &operator/=(_uint128 &lft, short rhs) {
	return lft /= (_int128)rhs;
}
inline _uint128 &operator/=(_uint128 &lft, unsigned short rhs) {
	return lft /= (_uint128)rhs;
}


// operator %= for built in integral types as second argument
inline _int128  &operator%=(_int128  &lft, __int64 rhs) {
	return lft %= (_int128)rhs;
}
inline _int128  &operator%=(_int128  &lft, unsigned __int64 rhs) {
	return lft %= (_uint128)rhs;
}
inline _int128  &operator%=(_int128  &lft, long rhs) {
	return lft %= (_int128)rhs;
}
inline _int128  &operator%=(_int128  &lft, unsigned long rhs) {
	return lft %= (_uint128)rhs;
}
inline _int128  &operator%=(_int128  &lft, int rhs) {
	return lft %= (_int128)rhs;
}
inline _int128  &operator%=(_int128  &lft, unsigned int rhs) {
	return lft %= (_uint128)rhs;
}
inline _int128  &operator%=(_int128  &lft, short rhs) {
	return lft %= (_int128)rhs;
}
inline _int128  &operator%=(_int128  &lft, unsigned short rhs) {
	return lft %= (_uint128)rhs;
}

inline _uint128 &operator%=(_uint128 &lft, __int64 rhs) {
	return lft %= (_int128)rhs;
}
inline _uint128 &operator%=(_uint128 &lft, unsigned __int64 rhs) {
	return lft %= (_uint128)rhs;
}
inline _uint128 &operator%=(_uint128 &lft, long rhs) {
	return lft %= (_int128)rhs;
}
inline _uint128 &operator%=(_uint128 &lft, unsigned long rhs) {
	return lft %= (_uint128)rhs;
}
inline _uint128 &operator%=(_uint128 &lft, int rhs) {
	return lft %= (_int128)rhs;
}
inline _uint128 &operator%=(_uint128 &lft, unsigned int rhs) {
	return lft %= (_uint128)rhs;
}
inline _uint128 &operator%=(_uint128 &lft, short rhs) {
	return lft %= (_int128)rhs;
}
inline _uint128 &operator%=(_uint128 &lft, unsigned short rhs) {
	return lft %= (_uint128)rhs;
}


// operator &= for built in integral types as second argument
inline _int128  &operator&=(_int128  &lft, __int64 rhs) {
	return lft &= (_int128)rhs;
}
inline _int128  &operator&=(_int128  &lft, unsigned __int64 rhs) {
	return lft &= (_uint128)rhs;
}
inline _int128  &operator&=(_int128  &lft, long rhs) {
	return lft &= (_int128)rhs;
}
inline _int128  &operator&=(_int128  &lft, unsigned long rhs) {
	return lft &= (_uint128)rhs;
}
inline _int128  &operator&=(_int128  &lft, int rhs) {
	return lft &= (_int128)rhs;
}
inline _int128  &operator&=(_int128  &lft, unsigned int rhs) {
	return lft &= (_uint128)rhs;
}
inline _int128  &operator&=(_int128  &lft, short rhs) {
	return lft &= (_int128)rhs;
}
inline _int128  &operator&=(_int128  &lft, unsigned short rhs) {
	return lft &= (_uint128)rhs;
}

inline _uint128 &operator&=(_uint128 &lft, __int64 rhs) {
	return lft &= (_int128)rhs;
}
inline _uint128 &operator&=(_uint128 &lft, unsigned __int64 rhs) {
	return lft &= (_uint128)rhs;
}
inline _uint128 &operator&=(_uint128 &lft, long rhs) {
	return lft &= (_int128)rhs;
}
inline _uint128 &operator&=(_uint128 &lft, unsigned long rhs) {
	return lft &= (_uint128)rhs;
}
inline _uint128 &operator&=(_uint128 &lft, int rhs) {
	return lft &= (_int128)rhs;
}
inline _uint128 &operator&=(_uint128 &lft, unsigned int rhs) {
	return lft &= (_uint128)rhs;
}
inline _uint128 &operator&=(_uint128 &lft, short rhs) {
	return lft &= (_int128)rhs;
}
inline _uint128 &operator&=(_uint128 &lft, unsigned short rhs) {
	return lft &= (_uint128)rhs;
}


// operator |= for built in integral types as second argument
inline _int128  &operator|=(_int128  &lft, __int64 rhs) {
	return lft |= (_int128)rhs;
}
inline _int128  &operator|=(_int128  &lft, unsigned __int64 rhs) {
	return lft |= (_uint128)rhs;
}
inline _int128  &operator|=(_int128  &lft, long rhs) {
	return lft |= (_int128)rhs;
}
inline _int128  &operator|=(_int128  &lft, unsigned long rhs) {
	return lft |= (_uint128)rhs;
}
inline _int128  &operator|=(_int128  &lft, int rhs) {
	return lft |= (_int128)rhs;
}
inline _int128  &operator|=(_int128  &lft, unsigned int rhs) {
	return lft |= (_uint128)rhs;
}
inline _int128  &operator|=(_int128  &lft, short rhs) {
	return lft |= (_int128)rhs;
}
inline _int128  &operator|=(_int128  &lft, unsigned short rhs) {
	return lft |= (_uint128)rhs;
}

inline _uint128 &operator|=(_uint128 &lft, __int64 rhs) {
	return lft |= (_int128)rhs;
}
inline _uint128 &operator|=(_uint128 &lft, unsigned __int64 rhs) {
	return lft |= (_uint128)rhs;
}
inline _uint128 &operator|=(_uint128 &lft, long rhs) {
	return lft |= (_int128)rhs;
}
inline _uint128 &operator|=(_uint128 &lft, unsigned long rhs) {
	return lft |= (_uint128)rhs;
}
inline _uint128 &operator|=(_uint128 &lft, int rhs) {
	return lft |= (_int128)rhs;
}
inline _uint128 &operator|=(_uint128 &lft, unsigned int rhs) {
	return lft |= (_uint128)rhs;
}
inline _uint128 &operator|=(_uint128 &lft, short rhs) {
	return lft |= (_int128)rhs;
}
inline _uint128 &operator|=(_uint128 &lft, unsigned short rhs) {
	return lft |= (_uint128)rhs;
}


// operator ^= for built in integral types as second argument
inline _int128  &operator^=(_int128  &lft, __int64 rhs) {
	return lft ^= (_int128)rhs;
}
inline _int128  &operator^=(_int128  &lft, unsigned __int64 rhs) {
	return lft ^= (_uint128)rhs;
}
inline _int128  &operator^=(_int128  &lft, long rhs) {
	return lft ^= (_int128)rhs;
}
inline _int128  &operator^=(_int128  &lft, unsigned long rhs) {
	return lft ^= (_uint128)rhs;
}
inline _int128  &operator^=(_int128  &lft, int rhs) {
	return lft ^= (_int128)rhs;
}
inline _int128  &operator^=(_int128  &lft, unsigned int rhs) {
	return lft ^= (_uint128)rhs;
}
inline _int128  &operator^=(_int128  &lft, short rhs) {
	return lft ^= (_int128)rhs;
}
inline _int128  &operator^=(_int128  &lft, unsigned short rhs) {
	return lft ^= (_uint128)rhs;
}

inline _uint128 &operator^=(_uint128 &lft, __int64 rhs) {
	return lft ^= (_int128)rhs;
}
inline _uint128 &operator^=(_uint128 &lft, unsigned __int64 rhs) {
	return lft ^= (_uint128)rhs;
}
inline _uint128 &operator^=(_uint128 &lft, long rhs) {
	return lft ^= (_int128)rhs;
}
inline _uint128 &operator^=(_uint128 &lft, unsigned long rhs) {
	return lft ^= (_uint128)rhs;
}
inline _uint128 &operator^=(_uint128 &lft, int rhs) {
	return lft ^= (_int128)rhs;
}
inline _uint128 &operator^=(_uint128 &lft, unsigned int rhs) {
	return lft ^= (_uint128)rhs;
}
inline _uint128 &operator^=(_uint128 &lft, short rhs) {
	return lft ^= (_int128)rhs;
}
inline _uint128 &operator^=(_uint128 &lft, unsigned short rhs) {
	return lft ^= (_uint128)rhs;
}

_int128   _strtoi128(const char    *str, char    **end, int radix);
_uint128  _strtoui128(const char    *str, char    **end, int radix);
_int128   _wcstoi128(const wchar_t *str, wchar_t **end, int radix);
_uint128  _wcstoui128(const wchar_t *str, wchar_t **end, int radix);

char     *_i128toa(_int128   value, char    *str, int radix);
char     *_ui128toa(_uint128  value, char    *str, int radix);
wchar_t  *_i128tow(_int128   value, wchar_t *str, int radix);
wchar_t  *_ui128tow(_uint128  value, wchar_t *str, int radix);
struct fc {
	bool minus = false;
	bool plus = false;
	bool zero = false;
	bool blank = false;
	bool hash = false;
};

int sPrintf128(char *buffer, const int buflen, const char* FormatStr, const _int128  &n);
int sPrintf128(char *buffer, const int buflen, const char* FormatStr, const _uint128  &n);


inline char radixLetter(unsigned int c) {
	return (c < 10) ? ('0' + c) : ('a' + (c - 10));
}

/* return true if ch is a valid octal digit */
inline bool iswodigit(wchar_t ch) {
	return ('0' <= ch) && (ch < '8');
}

extern const _int128  _I128_MIN, _I128_MAX;
extern const _uint128 _UI128_MAX;

std::istream  &operator>>(std::istream  &s, _int128  &n);
std::ostream  &operator<<(std::ostream  &s, const _int128  &n);
std::istream  &operator>>(std::istream  &s, _uint128 &n);
std::ostream  &operator<<(std::ostream  &s, const _uint128 &n);

std::wistream &operator>>(std::wistream &s, _int128  &n);
std::wostream &operator<<(std::wostream &s, const _int128  &n);
std::wistream &operator>>(std::wistream &s, _uint128 &n);
std::wostream &operator<<(std::wostream &s, const _uint128 &n);

typedef struct {
	_int128 quot;
	_int128 rem;
} _int128div_t;

_int128div_t divrem(const _int128 &a, const _int128 &b);

typedef struct {
	_uint128 quot;
	_uint128 rem;
} _uint128div_t;

_uint128div_t divrem(const _uint128 &a, const _uint128 &b);


/* calculate x^n.   Beware of overflow!
e.g. 2^128, 3^81, 4^64 etc will cause an overflow */
_int128 power(const _int128 &x, int n);

/* find nth root of a  i.e. return r such that abs(r^n) <= abs(a)
if n < 1 there is no solution.
If n is even and a is negative there is no solution.
if n = 1 or a=0 or a = 1, r = a */
_int128 nthroot(const _int128 &aa, int n);

/* calculate GCD of u and v. */
_uint128 gcd(_uint128 u, _uint128 v);
_int128 extendedGcd(const _int128 &u, const _int128 &v, _int128 &x, _int128 &y);

/* get modular inverse of a wrt m */
_int128 modMultInv(_int128 a, _int128 m);

/* get absolute value of x */
static inline _int128 abs(const _int128 &x) {
	if (x.isNegative())
		return -x;
	else
		return x;
}

/* convert double x to _int128. Number is truncated i.e. fractional part
is discarded, like casting to int. */
_int128 doubleTo128(const double x);

/* multiply 64 bit x 64 bit to get 128 bits. This is a more efficient alternative
to casting each of the 64 bit values to 128 bits before multiplication, while
still ensuring that overflow will not occur. */

/* this version is for unsigned integers. */
#define ui128mult(p128, a64, b64)     \
	LO64(p128) =_umul128(a64, b64, &HI64(p128))
/* this version is for signed integers. */
#define i128mult(p128, a64, b64)     \
	LO64(p128) =_mul128(a64, b64, (long long*)&HI64(p128))