Adapted from https://github.com/JesperMikkelsen/Big-Numbers

This is an attempt to integrate 128 bit integers as fully as possible into
Visual Studio C++. The objective that you just declare _int128 and _uint128
variables and use them as normal. However library functions such as printf
can't be enhanced to handle 128 byte integers directly.

This VS project consists of:

An assembly module  Int128x64.asm that implements 128 bit arithmetic.
The C/C++ declarations for them are:
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

	The overhead for a function call is significant here because the compiler cannot
	inline the functions or do any optimisation. Although these are actually pure 
	functions there is no way to indicate this to the compiler.

A header file Int128.h that creates classes for signed and unsigned 128 bit integers
These classes overload all the arithmetic, comparison, assignment and bit manipulation
operators. Also constructors etc so that 128 bit integers and standard integers can be 
used together and conversions are performed automatically.

For addition, the _addcarry_u64 intrinsic function is used because it is faster 
than calling the assembly routine. 

The & (bitwise and), | (bitwise or), and ^ (exclusive or) operators are also 
overloaded for 128 bit variables but don't requre assembly subroutines to
implement them.

std::cout does simple conversion of 128 bit integers to decimal, with - sign if negative.
No formatting control is available with std::cout.

An _int128 or _uint128 number has the following member functions:
bool   x.isNegative()
bool   x.isZero()
int    x.msb()        // bit number of most significant 1 bit of x
int    x.lsb()        // bit number of least significant 1 bit of x
int    x.popcnt()     // number of bits in x set to 1
double x.toDouble()   // convert x to double (can also be invoked by casting to double:
						double x = i.toDouble(); is equivalent to
						double x = (double)i;  )

File Int128Str.cpp contains the following functions:

/* convert integer to text (ascii or wide character) 
(design based on _itoa() & _itow() functions.)
str   = pointer to char or wchar_t string to that, on return, contains the digits
radix = positive number between 2 and 36. 
These functions can write past the end of a buffer that is too small. To prevent 
buffer overruns, ensure that buffer is large enough to hold the converted digits 
plus the trailing null-character and a sign character. */
char     *_i128toa( _int128 value,    char *str, int radix)
wchar_t  *_i128tow( _int128 value, wchar_t *str, int radix)
char    *_ui128toa(_uint128 value,    char *str, int radix)
wchar_t *_ui128tow(_uint128 value, wchar_t *str, int radix)

/* convert a 128 bit integer to text. modelled on sprintf_s */
int sPrintf128(char *buffer, const int buflen, const char* FormatStr, const _int128  &n)
int sPrintf128(char *buffer, const int buflen, const char* FormatStr, const _uint128  &n)

/* convert text (ascii or wide character) to integer 
str   = pointer to char or wchar_t string to be converted
end   = pointer to pointer to char or wchar_t. 
        On return is set to character after last digit processed.
radix = zero or positive number between 2 and 36.
        if radix = 0 input can be octal, decimal or hex 
	                (0X or 0x prefix for hex, 
                     0 prefix for octal, 
		             no prefix for decimal) */
_int128 _strtoi128(const char *str, char **end, int radix)
_uint128 _strtoui128(const char *str, char **end, int radix)
_int128 _wcstoi128(const wchar_t *str, wchar_t **end, int radix)
_uint128 _wcstoui128(const wchar_t *str, wchar_t **end, int radix)

/* miscellaneous functions */
_int128 abs(const _int128 &x)
 _int128div_t divrem(const  _int128 &a, const  _int128 &b)
_uint128div_t divrem(const _uint128 &a, const _uint128 &b)
_int128 doubleTo128(const double x)    //convert double to signed 128 bit integer

_int128 extendedGcd(const _int128 &a, const _int128 &b, _int128 &x, _int128 &y)
_uint128 gcd(_uint128 u, _uint128 v)
_int128 modMultInv(_int128 a, _int128 m)  // modular inverse of a wrt m
_int128 nthroot(const _int128 &x, int n)  // calculate nth root of x
_int128 power(const _int128 &x, int n)     // calculate x^n

macros to multiply 64 bit x 64 bit to get 128 bits. These are a more efficient 
alternative to casting each of the 64 bit values to 128 bits before 
multiplication, while still ensuring that overflow will not occur. The code
generated uses an intrinsic function instead of assembler. It has just one 
multiply instruction and will usually be inlined, which makes it faster.

	for unsigned integers:
ui128mult(p128, a64, b64)    // (p128 = a64 * b64)

	for signed integers:
i128mult(p128, a64, b64)