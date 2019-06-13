/***************************************************************
File: Int128Str.cpp

defines functions _strtoi128, _strtoui128, _i128toa, _ui128toa,
				  _i128tow, _ui128tow, _wcstoi128, _wcstoui128

also: abs, doubleTo128, gcd, extendedGcd, nthroot, power, modMultInv, divrem
****************************************************************/

#include "pch.h"

//#define _CRT_SECURE_NO_WARNINGS
#include <assert.h>

#pragma warning(disable : 4073)

const _int128  _I128_MIN(0x8000000000000000, 0x0000000000000000);
const _int128  _I128_MAX(0x7fffffffffffffff, 0xffffffffffffffff);
const _uint128 _UI128_MAX(0xffffffffffffffff, 0xffffffffffffffff);

// Conversion from _int128/_uint128 to string

// Number of digits that should be appended to the string for each loop.
// Also used as maximal number of digits to scan, that fits in 32 bit UINT
// digitCount[R] = floor(logR(_UI32_MAX+1))
// (index = radix = [2..36])
static const unsigned short digitCount[37] = {
   0, 0,32,20,16,13,12,11
 ,10,10, 9, 9, 8, 8, 8, 8
 , 8, 7, 7, 7, 7, 7, 7, 7
 , 6, 6, 6, 6, 6, 6, 6, 6
 , 6, 6, 6, 6, 6
};

// Highest power of radix that fits in 32 bit (index = radix)
// For radix 2,4,8,16,32 special code is used which doesn't use this table.
// For all other values the element is used to get as many digits as possible
// by modulus and division. (powRadix[r] == pow(r,digitCount[r])
static const _uint128 powRadix[] = {
  0          //  not used
 ,0          //  not used
 ,0          //  not used
 ,0xcfd41b91 //  3486784401 =  3^20
 ,0          //  not used
 ,0x48c27395 //  1220703125 =  5^13
 ,0x81bf1000 //  2176782336 =  6^12
 ,0x75db9c97 //  1977326743 =  7^11
 ,0          //  not used
 ,0xcfd41b91 //  3486784401 =  9^10
 ,0x3b9aca00 //  1000000000 = 10^9
 ,0x8c8b6d2b //  2357947691 = 11^9
 ,0x19a10000 //   429981696 = 12^8
 ,0x309f1021 //   815730721 = 13^8
 ,0x57f6c100 //  1475789056 = 14^8
 ,0x98c29b81 //  2562890625 = 15^8
 ,0          //  not used
 ,0x18754571 //   410338673 = 17^7
 ,0x247dbc80 //   612220032 = 18^7
 ,0x3547667b //   893871739 = 19^7
 ,0x4c4b4000 //  1280000000 = 20^7
 ,0x6b5a6e1d //  1801088541 = 21^7
 ,0x94ace180 //  2494357888 = 22^7
 ,0xcaf18367 //  3404825447 = 23^7
 , 0xb640000 //   191102976 = 24^6
 , 0xe8d4a51 //   244140625 = 25^6
 ,0x1269ae40 //   308915776 = 26^6
 ,0x17179149 //   387420489 = 27^6
 ,0x1cb91000 //   481890304 = 28^6
 ,0x23744899 //   594823321 = 29^6
 ,0x2b73a840 //   729000000 = 30^6
 ,0x34e63b41 //   887503681 = 31^6
 ,0          //  not used
 ,0x4cfa3cc1 //  1291467969 = 33^6
 ,0x5c13d840 //  1544804416 = 34^6
 ,0x6d91b519 //  1838265625 = 35^6
 ,0x81bf1000 //  2176782336 = 36^6
};

/* convert integer v to a text string using radix.
there would be a buffer overflow problem here if str were too small */
#define ULTOSTR(v, str, radix)        \
{ if(sizeof(Ctype) == sizeof(char))   \
    _ultoa(v, (char*)str, radix);     \
  else                                \
    _ultow(v, (wchar_t*)str, radix);  \
}

#define STRLEN(str)      ((sizeof(Ctype)==sizeof(char))?strlen((char*)str):wcslen((wchar_t*)str))
#define STRCPY(dst, src) ((sizeof(Ctype)==sizeof(char))?(Ctype*)strcpy((char*)dst, (char*)src):(Ctype*)wcscpy((wchar_t*)dst, (wchar_t*)src))
#define STRREV(str)      ((sizeof(Ctype)==sizeof(char))?(Ctype*)_strrev((char*)str):(Ctype*)_wcsrev((wchar_t*)str))

/* Function template for conversion of 128 bit integer to an ascii string
These functions can write past the end of a buffer that is too small. To prevent
buffer overruns, ensure that buffer is large enough to hold the converted digits
plus the trailing null-character and a sign character. */
template<class Int128Type, class Ctype>
Ctype *int128toStr(Int128Type value, Ctype *str, unsigned int radix) {
	if ((radix < 2) || (radix > 36)) {
		errno = EINVAL;
		str[0] = 0;
		return str;
	}
	errno = 0;
	bool setSign = false;
	if (value.isZero()) {
		str[0] = '0';
		str[1] = 0;
		return str;
	}

	Ctype *s = str;
	switch (radix) {
	case 2:
	case 4:
	case 16:
	{ const unsigned int dc = digitCount[radix];
	for (int i = 3; i >= 0; i--) {
		if (value.s4.i[i]) {
			Ctype tmpStr[40];
#pragma warning( disable : 4996)
			/* there would be a buffer overflow problem here if tmpStr was too small */
			ULTOSTR(value.s4.i[i], tmpStr, radix);  // convert 32 bits to text in tmpStr
			const size_t l = STRLEN(tmpStr);
			if (s != str) {
				for (size_t i = dc - l; i--;) *(s++) = '0'; // fill up with zeroes, if not leading digits
			}
			//STRCPY(s, tmpStr);         // changed to avoid internal compiler error!!!
			/* there would be a buffer overflow problem here if s was too small */
			if (sizeof(Ctype) == sizeof(char))
				strcpy((char*)s, (char *)tmpStr);
			else
				wcscpy((wchar_t*)s, (wchar_t*)tmpStr);
			s += l;
		}
		else if (s != str) {
			for (size_t i = dc; i--;) *(s++) = '0'; // fill up with zeroes, if not leading digits
		}
	}
	*s = 0;
	}
	return str;

	case 8: // Get 3 bits/digit giving 30 bits/loop, ie 10 digits/loop
	case 32: // Get 5 bits/digit giving 30 bits/loop too! which is 6 digits/loop
	{	const unsigned int shft = (radix == 32) ? 5 : 3;
		const unsigned int mask = (1 << shft) - 1;
		const unsigned int dpl = 30 / shft;
		_uint128   v = value;
		for (;;) {
			unsigned int v30 = v.s4.i[0] & ((1 << 30) - 1);
			v >>= 30;
			unsigned int count;
			for (count = 0; v30; count++, v30 >>= shft) {
				*(s++) = radixLetter(v30 & mask);
			}
			if (v.isZero()) break;
			while (count++ < dpl) *(s++) = '0';
		}
	}
	break;

	case 10:
		if (value.isNegative()) {
			value = -value;
			setSign = true;
		}
		// NB continue case
	default:
		const unsigned int dc = digitCount[radix];
		const _uint128 &divisor = powRadix[radix];
		_uint128        v = value;
		for (;;) {
			const unsigned int c = v % divisor;
			Ctype tmpStr[40];
			/* there would be a buffer overflow problem here if tmpStr were too small */
			ULTOSTR(c, tmpStr, radix);
			STRREV(tmpStr);
			//STRCPY(s, tmpStr);  // changed to avoid internal compiler error!!!
			if (sizeof(Ctype) == sizeof(char))
				strcpy((char*)s, (char *)tmpStr);
			else
				wcscpy((wchar_t*)s, (wchar_t*)tmpStr);
			size_t l = STRLEN(tmpStr);
			s += l;
			v /= divisor;
			if (v) {
				while (l++ < dc) *(s++) = '0'; // append zeroes
			}
			else {
				break;
			}
		}
		if (setSign) 
			*(s++) = '-';
		break;
	}

	*s = 0;
	return STRREV(str);
}

char *_i128toa(_int128 value, char *str, int radix) {
	return int128toStr<_int128, char>(value, str, radix);
}

wchar_t *_i128tow(_int128 value, wchar_t *str, int radix) {
	return int128toStr<_int128, wchar_t>(value, str, radix);
}

char*_ui128toa(_uint128 value, char *str, int radix) {
	return int128toStr<_uint128, char>(value, str, radix);
}

wchar_t *_ui128tow(_uint128 value, wchar_t *str, int radix) {
	return int128toStr<_uint128, wchar_t>(value, str, radix);
}

// Conversion from string to _int128/_uint128

static inline bool isRadixDigit(wchar_t ch, unsigned int radix, unsigned int &value) {
	if (!iswalnum(ch)) return false;
	const unsigned int v = iswdigit(ch) ?
		(ch - '0') : (ch - (iswupper(ch) ? 'A' : 'a') + 10);
	if (v >= radix) return false;
	value = v;
	return true;
}

/* template function to convert ascii to 128 byte integer.
max signed integer is 340 282366 920938 463463 374607 431768 211455
On return errno may be set to EINVAL or ERANGE.

The function is called by "strtoint128 <itype, Ctype, withsign>(str, end, radix)
where: itype = _int128 or _uint128    (output either signed or unsigned 128 bit integer)
       ctype = char or wchar_t        (input is 8 bit or wide characters)
	   withsign = true or false       (true for signed, false for unsigned)
	   str   = pointer to char or wchar_t string to be converted
	   end   = pointer to pointer to char or wchar_t. 
	           On return is set to character after last digit processed.
	   radix = zero or positive number between 2 and 36.
	   if radix = 0 input can be octal, decimal or hex 
	                (0X or 0x prefix for hex, 
                     0 prefix for octal, 
		             no prefix for decimal) */
template<class Int128Type, class Ctype, bool withSign>
Int128Type strtoint128(const Ctype *s, Ctype **end, unsigned int radix) {

	if ((s == NULL) || ((radix != 0) && ((radix < 2) || (radix > 36)))) {
		errno = EINVAL;
		return 0;
	}
	errno = 0;
	bool negative = false;
	bool gotDigit = false;

	while (iswspace(*s)) 
		s++; // skip whitespace

	if (*s == '-') {   // read optional sign
		s++;
		negative = true;
	}
	else if (*s == '+') {
		s++;
	}

	unsigned int digit;
	if (radix == 0) { // first determine if radix is 8, 10 or 16
		if (*s == '0') {
			gotDigit = true;
			s++;
			if (end)
				*end = (Ctype*)s;
			if ((*s == 'x') || (*s == 'X')) {
				radix = 16; s++;
			}
			else {
				radix = 8;
			}
			if (!isRadixDigit(*s, radix, digit)) { // we've scanned "0[x]",
				return 0;                           // if no more digits, then 0
			}
		}
		else if (iswdigit(*s)) {
			radix = 10;
		}
		else {
			return 0; // nothing recognized
		}
	}

	Int128Type result128;
	bool       overflow = false;
	const unsigned int maxDigitCount32 = digitCount[radix];

	if (isRadixDigit(*(s++), radix, digit)) {
		gotDigit = true;
		while ((digit == 0) && isRadixDigit(*s, radix, digit))
			s++; // skip leading zeroes

		bool firstChunk = true;
		unsigned int result32 = digit;

		if ((radix & -(int)radix) == radix) { // is radix 2,4,8,16 or 32
			const unsigned int maxBitCount = withSign ? 127 : 128;
			unsigned int       totalBitCount;
			const unsigned int bitsPerDigit = 32 / maxDigitCount32;
			const unsigned int maxBitCount32 = bitsPerDigit * maxDigitCount32;
			for (unsigned int bitCount32 = bitsPerDigit;; result32 = bitCount32 = 0) {
				while (isRadixDigit(*(s++), radix, digit)) {
					result32 <<= bitsPerDigit;
					result32 |= digit;
					if ((bitCount32 += bitsPerDigit) == maxBitCount32) break;
				}
				if (firstChunk) {
					result128 = result32;
					firstChunk = false;
					totalBitCount = bitCount32;
				}
				else if (bitCount32) {
					if ((totalBitCount += bitCount32) > maxBitCount) {
						const unsigned char mask = ~(((withSign && !negative) ? 0x7f : 0xff) >> (bitCount32 & 7));
						if (result128.s16.b[15 - (bitCount32 >> 3)] & mask) {
							overflow = true;
							if (bitCount32 == maxBitCount32) {
								while (isRadixDigit(*(s++), radix, digit));
							}
							break;
						}
					}
					result128 <<= bitCount32;
					if (result32) 
						result128 |= result32;
				}
				else {
					break;
				}
				if (bitCount32 < maxBitCount32) 
					break;
			}
		}
		else {
			/* radix is not 2, 4, 8, 16 or 32 */
			for (unsigned int digitCount32 = 1, p32 = radix;; p32 = 1, result32 = digitCount32 = 0) {
				while (isRadixDigit(*(s++), radix, digit)) {
					result32 *= radix;
					result32 += digit;
					p32 *= radix;
					if (++digitCount32 == maxDigitCount32)
						break;
				}

				if (firstChunk) {
					result128 = result32;
					firstChunk = false;
				}
				else if (digitCount32) {
					const unsigned int lastH = result128.s4.i[3];
					result128 *= p32;
					if (result32)
						result128 += result32;
					if (lastH && ((result128.s4.i[3] == 0) || (withSign && !negative&&result128.isNegative()))) {
						overflow = true;
						if (digitCount32 == maxDigitCount32) {
							while (isRadixDigit(*(s++), radix, digit));
						}
						break;
					}
				}
				else {
					break;
				}
				if (digitCount32 < maxDigitCount32)
					break;
			}
		}
	}

	if (!gotDigit)
		return 0;

	if (end)
		*end = (Ctype*)s - 1;

	if (overflow) {
		errno = ERANGE;
		//return withSign ? (negative ? _I128_MIN : _I128_MAX) : _UI128_MAX;
		if (withSign)
			return negative ? _I128_MIN : _I128_MAX;
		else
			return _UI128_MAX;
	}
	//return negative ? -result128 : result128;
	if (negative)
		return -result128;
	else
		return result128;
}

/*str   = pointer to char or wchar_t string
  end   = pointer to pointer to char or wchar_t. On return is st to character after last digit processed.
  radix = zero or positive number between 2 and 36. 
  if radix = 0 input can be octal (0X or 0x prefix for hex, 
                                   0 prefix for octal, 
								   no prefix for decimal) */
_int128 _strtoi128(const char *str, char **end, int radix) {
	return strtoint128<_int128, char, true >(str, end, radix);
}

_uint128 _strtoui128(const char *str, char **end, int radix) {
	return strtoint128<_uint128, char, false>(str, end, radix);
}

_int128 _wcstoi128(const wchar_t *str, wchar_t **end, int radix) {
	return strtoint128<_int128, wchar_t, true >(str, end, radix);
}

_uint128 _wcstoui128(const wchar_t *str, wchar_t **end, int radix) {
	return strtoint128<_uint128, wchar_t, false>(str, end, radix);
}


/* << definitions for output streams (std::cout etc) 
definitions below are crude but they work. Only decimal output
is supported for direct output via cout & wcout*/
std::ostream  &operator<<(std::ostream  &s, const _int128  &n) {
	static char buf[50];  // a 128 bit integer is at most 39 decimal digits, so 50 is plenty
	return s << _i128toa(n, buf, 10);
}

std::ostream  &operator<<(std::ostream  &s, const _uint128 &n) {
	static char buf[50];  // a 128 bit integer is at most 39 decimal digits, so 50 is plenty
	return s << _ui128toa(n, buf, 10);
}

std::wostream &operator<<(std::wostream &s, const _int128  &n) {
	static WCHAR buf[50];  // a 128 bit integer is at most 39 decimal digits, so 50 is plenty
	return s << _i128tow(n, buf, 10);
}

std::wostream &operator<<(std::wostream &s, const _uint128 &n) {
	static WCHAR buf[50];  // a 128 bit integer is at most 39 decimal digits, so 50 is plenty
	return s << _ui128tow(n, buf, 10);
}

/*
Convert a 128-bit integer to an ascii string
The number can be output in octal, decimal or hexadecimal
Process format string of form "%[flags][width][.precision]type"
Width specification

The width argument is a non-negative integer that controls the minimum number of 
characters that are output. If the number of characters in the output value is 
less than the specified width, blanks are added to the left or the right of the 
values, depending on whether or not the left-alignment flag (-) is specified, 
until the minimum width is reached.

The precision specifies the minimum number of digits to be printed. If the 
number of digits in the argument is less than precision, the output value is 
padded on the left with zeros. The value is not truncated when the number 
of digits exceeds precision.

Flag characters
Flag 	Meaning 									Default
- 		Left align the result within the given 		Right align.
		field width. 	
+ 		Use a sign (+ or -) to prefix the output 	Sign appears only for negative  (-).
		value if it is of a signed type. 			signed values
0 		Leading zeros are added until the minimum 	No padding.
        width is reached. If both 0 and - appear,
		the 0 is ignored. If 0 is specified for an 
		integer format and a precision specification 
		is also present the  0 is ignored.
 
blank (' ') 	Use a blank to prefix the output    No blank appears.
       value if it is signed and positive. The 
	   blank is ignored if both the blank and 
	   + flags appear. 	
# 	   When it's used with the o, x, or X format,   No prefix appears.
       the # flag uses 0, 0x, or 0X, respectively, 
	   to prefix any nonzero output value. 
*/
/* to be completed. Invalid format string may cause abort
Process format string of form "%[flags][width][.precision]type" 
 */
static void getFormat(const char FormatStr[], fc &flags, int &precision, int &width, int &radix) {
	int ix = 0;
	auto len = strnlen_s(FormatStr, 50);
	assert(len < 50);  // assert fails if string is not null-terminated

	/* set default values */
	width = 0;
	precision = 0;
	radix = 10;
	flags.minus = flags.plus = flags.blank = flags.hash = flags.zero = false;

	if (FormatStr[0] == '%')
		ix++;

	do {
		if (FormatStr[ix] == '-') { flags.minus = true; ix++;  continue; }
		else if (FormatStr[ix] == '+') { flags.plus = true; ix++;  continue; }
		else if (FormatStr[ix] == ' ') { flags.blank = true; ix++;  continue; }
		else if (FormatStr[ix] == '#') { flags.hash = true; ix++;  continue; }
		else if (FormatStr[ix] == '0') { flags.zero = true; ix++;  break; }
		else break;
	} while (ix < len);

	assert(ix < len);
	while (ix < len && isdigit(FormatStr[ix])) {
		width *= 10;
		width += FormatStr[ix] - '0';
		ix++;
	}
	assert(ix < len);

	if (FormatStr[ix] == '.') {
		/* precision specified */
		while (ix < len && isdigit(FormatStr[ix])) {
			precision *= 10;
			precision += FormatStr[ix] - '0';
			ix++;
		}
	}
	assert(ix < len);
	 /* process type */
	char type = FormatStr[ix];

	if (type == 'd' || type == 'i' || type == 'u') {
		radix = 10;
		return;
	}

	if (type == 'o') {
		radix = 8;
		return;
	}
	if (type == 'x' || type == 'X') {
		radix = 16;
		return;
	}

	std::cout << "format string error in: " << FormatStr << '\n';
}

/* Convert integer n to an ascii string. modelled on sprintf_s. 
Types d, i and u are all treated the same. Appropriate code for signed
or unsigned integers is executed based on the parameter type. 
The value returned is the length of the string placed in buffer,
or -1 if an error occurred.
*/
int sPrintf128(char *buffer, const int buflen, const char* FormatStr, const _int128  &n){
	std::string sbuffer;
	fc flags = { false, false, false, false, false };
	int precision=0, width=0, radix=10;

	getFormat(FormatStr, flags, precision, width, radix);

	/* ensure parameters are valid */
	if (radix != 8 && radix != 10 && radix != 16)
		radix = 10;  // ensure radix is valid
	if (width > buflen - 1)
		width = buflen - 1;
	if (width < 0)
		width = 0;
	if (precision > buflen - 1)
		precision = buflen - 1;
	if (precision < 0)
		precision = 0;

	/* potential buffer overflow problem here if buffer is too small,
	but radix can only be 8, 10 or 16, so 50 is enough */
	sbuffer.resize(50);
	if (radix == 10)
		_i128toa(n, &sbuffer[0], radix);  // convert int to chars.
	else
		/* for hex or octal print as if unsigned */
		_ui128toa(n, &sbuffer[0], radix);  // convert int to chars

	sbuffer.shrink_to_fit();

	int ip = 0;
	if (sbuffer[0] == '-')
		ip = 1;      // If ve 1st char is -move insertion point after sign
	else if (flags.plus) {
		sbuffer.insert(0, 1, '+');  // insert + as 1st character
		ip=1;    // move insertion point after sign
		}
		else if (flags.blank) {
			sbuffer.insert(0, 1, ' ');  // insert space in place of sign
			ip=1;
		}

	/* insert octal or hex prefix after sign if required  */
	if (flags.hash && !n.isZero()) {
		if (radix == 8) {
			sbuffer.insert(ip, 1, '0');
			ip++;
		}
		if (radix == 16) {
			sbuffer.insert(ip, "0x");
			ip+=2;
		}
	}
	
	if (precision > 0) {
		while (sbuffer.length() < precision) {
			sbuffer.insert(ip, 1, '0');  // insert zeros
		}
	}

	if (flags.zero && precision == 0) {
		while (sbuffer.length() < width) {
			sbuffer.insert(ip, 1, '0');  // insert zeros
		}
	}

	if (width > 0) {
		if (flags.minus) {
			/* number is left justified. If < minimum width pad on right with spaces */
			ptrdiff_t blanks = width - sbuffer.length();
			if (blanks > 0)
				sbuffer.append(blanks, ' ');
		}
		else {
			/* number is right justified.If < minimum width pad on left with spaces*/
			ptrdiff_t blanks = width - sbuffer.length();
			if (blanks > 0)
				sbuffer.insert(0, blanks, ' ');
		}
	}

	auto rv = strcpy_s(buffer, buflen, sbuffer.c_str());
	if (rv == 0)   // did copy work OK?
		return  (int)sbuffer.length();  // yes, return length
	else return -1;                     // no, return error
}

int sPrintf128(char *buffer, const int buflen, const char* FormatStr, const _uint128  &n) {

	std::string sbuffer;
	fc flags;
	int precision = 0, width = 0, radix = 10;

	getFormat(FormatStr, flags, precision, width, radix);

	/* ensure parameters are valid */
	if (radix != 8 && radix != 10 && radix != 16)
		radix = 10;  // ensure radix is valid
	if (width > buflen - 1)
		width = buflen - 1;
	
	if (precision > buflen - 1)
		precision = buflen - 1;
	
	/* /* potential buffer overflow problem here if buffer is too small,
	but radix can only be 8, 10 or 16, so 50 is enough */
	sbuffer.resize(50);
	_ui128toa(n, &sbuffer[0], radix);  // convert unsigned int to chars. 

	sbuffer.shrink_to_fit();
	
	int ip = 0;
	if (sbuffer[0] == '-')
		ip = 1;      // move insertion point after sign
	else if (flags.plus) {
		sbuffer.insert(0, 1, '+');
		ip = 1;    // move insertion point after sign
	}
	else if (flags.blank) {
		sbuffer.insert(0, 1, ' ');  // insert space in place of sign
		ip = 1;
	}

	/* insert octal or hex prefix after sign if required  */
	if (flags.hash && !n.isZero()) {
		if (radix == 8) {
			sbuffer.insert(ip, 1, '0');
			ip++;
		}
		if (radix == 16) {
			sbuffer.insert(ip, "0x");
			ip += 2;
		}
	}

	if (precision > 0) {
		while (sbuffer.length() < precision) {
			sbuffer.insert(ip, 1, '0');  // insert zeros
		}
	}

	if (flags.zero && precision == 0) {
		while (sbuffer.length() < width) {
			sbuffer.insert(ip, 1, '0');  // insert zeros
		}
	}

	if (width > 0) {
		if (flags.minus) {
			/* number is left justified. If < minimum width pad on right with spaces */
			ptrdiff_t blanks = width - sbuffer.length();
			if (blanks > 0)
				sbuffer.append(blanks, ' ');
		}
		else {
			/* number is right justified.If < minimum width pad on left with spaces*/
			ptrdiff_t blanks = width - sbuffer.length();
			if (blanks > 0)
				sbuffer.insert(0, blanks, ' ');
		}
	}

	auto rv = strcpy_s(buffer, buflen, sbuffer.c_str());
	if (rv == 0)
		return  (int)sbuffer.length();
	else return -1;
}


/* calculate x^n.   Beware of overflow!
e.g. 2^128, 3^81, 4^64 etc will cause an overflow */
_int128 power(const _int128 &x, int n) {
	_int128 p = x;
	_int128 r = 1;  // to hold result

	while (n > 0) {
		if (n % 2 == 1) {
			r *= p;
		}
		n /= 2;
		if (n == 0) break;
		p *= p;
	}
	return r;
}

/* convert double x to _int128. Number is truncated i.e. fractional part
is discarded, like casting to int.
e.g -1.9 is converted to -1
-0.9 is converted to 0
0.9 is converted to 0
1.9 is converted to 1, etc
NaN is converted to zero */
_int128 doubleTo128(const double x) {
	int sign;           // bit  63
	int exponent;       // bits 62 to 52. Bias = 1023
	long long mantissa; // bits 51 to 0

	/* struct below allows components of double to be separated */
	union {
		double dx;         // used to set up values
		long long llx;     // for debugging
		struct IEEE_754 {
			unsigned long long mantissa : 52;
			unsigned long long exponent : 11;
			unsigned long long     sign : 1;
		} dd;
	} floatd;

	_int128 rv;

	auto cat = std::fpclassify(x);
	if (cat == FP_ZERO || cat == FP_NAN || cat == FP_SUBNORMAL) {
		rv = 0;  // special case 
		return rv;
	}
	if (cat == FP_INFINITE || x >= pow(2, 127)) {
		return _I128_MAX;  // out of range return max possible value
	}
	if (x < -pow(2, 127)) {
		return _I128_MIN;  // out of range return min possible value
	}

	floatd.dx = x;
	sign = floatd.dd.sign;
	exponent = (int)floatd.dd.exponent - 1023;
	mantissa = floatd.dd.mantissa + (1LL << 52);

	rv = mantissa;
	/* shift rv so that it contains just the integer part */
	if (exponent > 52)
		rv <<= exponent - 52;
	else
		rv >>= 52 - exponent;
	if (sign)
		rv = -rv;

#ifdef _DEBUG
	//std::cout << "x = " << x;
	//printf_s(" llx = %llx ", floatd.llx);
	//std::cout << "sign = " << sign;
	//std::cout << " exponent = " << exponent;
	//printf_s(" mantissa = %llx \n", mantissa);
#endif

	return rv;
}


/* find nth root of a  i.e. return r such that abs(r^n) <= abs(a)
This is the inverse of the power function
if n < 1 there is no solution.
If n is even and a is negative there is no solution (unless using complex numbers).
if n = 1 or a=0 or a = 1, r = a */
_int128 nthroot(const _int128 &aa, int n) {
	bool neg = false;

	if ((n < 1) || (aa.isNegative() && ((n & 1) == 0)))
		/* either n is not +ve, or asking for even root of -ve number */
		throw std::exception("nthroot invalid parameters");

	if ((n == 1) || (aa == 1) || (aa == 0))
		return aa;

	_int128 root;
	if (aa.isNegative()) {
		root = -aa;    /* get abs(aa) */
		neg = true;
	}
	else
		root = aa;

	if (n >= 127)
		if (!neg)
			return 1;
		else
			return -1;

	/* for large n, calculate this way to avoid overflow */
	if (n >= 81) {
		// 3^81 is > max _int128, so result  < 3.
		if (root >= power(2, n))
			if (!neg)
				return 2;
			else
				return -2;
		else
			if (!neg)
				return 1;
			else return -1;
	}

	/* using floating point should give a much more accurate starting value */
	double adouble = root;  // convert aa to double 
	double rootd = pow(adouble, 1.0 / n);  // get nth root as floating point
	root =doubleTo128(rootd);    // convert root back to _int128

#ifdef _DEBUG
	//auto msb = a.msb();  // get most significant bit 
	//std::cout << "aa = " << aa << " msb = " << msb << '\n';
	//std::cout << "a  = " << a << '\n';
#endif
	if (pow(floor(rootd) + 1, n) > pow(2, 127))
		// can't do loop below because overflow would occur, so just return
		// best approximation.
		if (!neg)
			return root;
		else
			return -root;

	_int128 rnext = root;

	/* use newton-raphson method to find root. Generally the value obtained
	above is very accurate, but for square roots of integers > 2^126
	a 64 bit floating point cannot hold the exact value. */
	do {
		/* this check is needed because sometimes root does not converge to one
		integer but cycles through several different integers. This checks
		whether root satisfies the fundamental definition of a solution */
		if (power(root, n) <= aa && (power(root + 1, n) > aa))
			break;

		_int128 pan1 = power(root, n - 1);
		assert(pan1 > 0);         // check for overflow
		rnext = (aa / pan1 + root * (n - 1)) / n;
		if (rnext == root)
			break;
#ifdef _DEBUG
		//std::cout << "rnext = " << rnext << '\n';
#endif
		root = rnext;
	} while (true);

	if (!neg)
		return root;
	else
		return -root;
}

/*
see: http://lemire.me/blog/2013/12/26/fastest-way-to-compute-the-greatest-common-divisor/
calculate GCD of u and v. Use of division is avoided by using right shifts
and subtraction instead. Although more iterations are needed, it is faster
overall because shift is much faster than divide.
u and v are UNSIGNED, caller must use abs or similar if necessary.
Note: mathematically gcd(a,b) = gcd(-a,b) = gcd(b,-a) = gcd (-a,-b)
Note: GCD (greatest common denominator) is also known as HCF (Highest Common Factor).
*/
_uint128 gcd(_uint128 u, _uint128 v) {
	int shift, su, sv;
	if (u == 0) return v;
	if (v == 0) return u;

	/* shift represents number of factor 2s common to u and v */
	shift = (u | v).lsb();

	/* su represents factor 2s in u */
	su = u.lsb();
	u >>= su;          // remove factor 2s so that u is now odd 

	do {
		sv = v.lsb();  // remove factor 2s so that v is now odd
		v >>= sv;

		if (u > v) {
			auto t = v;
			v = u - v;
			u = t;
		}
		else
			v = v - u;
	} while (!v.isZero());

	return u << shift;  // put back factor 2s common to u and v
}

// simulate parallel assignment. temp is used to ensure correct result
// if newa and b are the same variable. It is recommended that all variables 
// are the same type.
#define PA(a,b,newa,newb) \
temp = (newa);  \
b = (newb);     \
a = temp;

/*
the extended Euclidean algorithm is an extension to the Euclidean algorithm,
which computes, besides the greatest common divisor of integers a and b, the
coefficients of Bézout's identity, that is integers x and y such that

	  a * x + b * y = gcd ( a , b )
*/
_int128 extendedGcd(const _int128 &a, const _int128 &b, _int128 &x, _int128 &y) {
	_int128 s = 0, old_s = 1;;
	_int128 t = 1, old_t = 0;
	_int128 r = b, old_r = a;
	_int128 quotient, temp, qs, qt;


	while (r != 0) {
		quotient = old_r / r;
		PA(old_r, r, r, old_r - quotient * r)
			qs = quotient * s;

		PA(old_s, s, s, old_s - qs)
			qt = quotient * t;
		PA(old_t, t, t, old_t - qt)
	}

	if (old_r < 0) {
		x = -old_s;
		y = -old_t;
		return -old_r;    // -r = gcd(a,b)
	}
	else {
		x = old_s;
		y = old_t;
		return old_r;     // r = gcd(a,b)
	}
}
#undef PA

/* the modular multiplicative inverse of an integer a modulo m is an integer x such that
(a* x) ≡ 1 (mod m). Given a and m this function calculates x.
N.B. a and m must be mutually prime i.e.no common factors i.e. gcd(a, m) is 1.
If they are, the function returns the inverse, otherwise it throws an exception.
*/
_int128 modMultInv(_int128 a, _int128 m) {
	_int128 llgcd, x, y;

	if (a > 0 && m > 0)
		a %= m;    // if a > m reduce it;
	llgcd = extendedGcd(a, m, x, y);  // use extended euclidian algorithm.
	if (llgcd != 1) {
		// there is no inverse
		throw std::range_error("Mod Inv error: Input parameters are not mutually prime");
	}
	if (x < 0) x += m; // ensure x in range 0 to m-1

#ifdef _DEBUG
	_int128 t;
	t = (a*x) % m;       // OK check provided a*x doesn't overflow
	assert(t == 1);    // check for bugs
#endif
	return x;
}

/* based on standard C div() function */
_int128div_t divrem(const _int128 &a, const _int128 &b) {
	static _int128div_t result;

	result.quot = a / b;               // quotient
	result.rem = a - result.quot*b;    // remainder (multiplication is faster than % operator)
	return result;
}

_uint128div_t divrem(const _uint128 &a, const _uint128 &b) {
	static _uint128div_t result;

	result.quot = a / b;               // quotient
	result.rem = a - result.quot*b;    // remainder (multiplication is faster than % operator)
	return result;
}