/* trn.c                 Thomas R. Nicely          2016.04.01.0400
*
* in this version a lot of unneeded code has been commented out or simply deleted.
 *
 * Freeware copyright (C) 2016 Thomas R. Nicely <http://www.trnicely.net>.
 * see https://web.archive.org/web/20161030053956/http://trnicely.net/misc/trn.zip
 * Released into the public domain by the author, who disclaims any legal
 * liability arising from its use. The routine iEvalExprMPZ is included
 * under the terms of the GNU GPL; see the header of that routine (at
 * the end of the code) for details.
 *
 * Common custom routines. Callable from co-compiled or linked codes.
 *
 * Header file is trn.h. It is assumed the compiler and
 * libraries have __non-trivial__ support for the long double
 * data type, and that the processor is Intel 386 compatible
 * (the macro __i386__ is defined).
 *
 * The macro __DJGPP__ is used to detect DJGPP.
 * The macro __DMC__ is used to detect 32-bit Digital Mars.
 * The macro __LINUX__ is used to detect LINUX GNU/Linux 10.0.
 * The macro __CYGWIN__ is used to detect Cygwin.
 * The macro __MINGW__ is used to detect MinGW.
 * The macro __BORLANDC__ is used to detect Borland/Inprise C.
 *   Version 5.00 or later is presumed.
 * The macro __MSVC__ is used to detect Microsoft Visual C++.
 *   Support for Microsoft Visual C++ is at alpha level.
 * The macro __WIN32__ is used to detect the Win32 API, which
 *   includes Cygwin, MinGW, Digital Mars in its default mode,
 *   and Borland C in its default mode.
 *
 * The macro __GMP__ indicates support for GMP, the GNU Multiple
 * Precision library. It is assumed to be present on all platforms
 * except Digital Mars, Borland C, and MSVC. To compile without GMP
 * on other platforms, add the command-line parameter "-D__NOGMP__"
 * (or the equivalent). If you _do_ have GMP with Digital Mars,
 * Borland C, or MSVC, activate the directive below defining
 *  __GMP__---and send me a copy of your libgmp files!
 *
 * The macro __MPFR__ indicates support for MPFR, the GNU Multiple
 * Precision Floating-point library with reliable rounding, version
 * 2.2.1 or later. It is assumed to be present on all platforms
 * except DJGPP, Digital Mars, Borland C, and MSVC. To compile
 * without MPFR on other platforms, add the command-line parameter
 * "-D__NOMPFR__" (or the equivalent). If you _do_ have MPFR with
 * DJGPP, Digital Mars, Borland C, or MSVC, activate the directive
 * below defining  __MPFR__---and send me a copy of your libmpfr
 * files!
 *
 * If support for DJGPP/Borland C style conio console functions
 * is required for compilers and platforms other than DJGPP and
 * Borland C, include conio3.h and compile and link conio3.c (q.v.).
 *
 * MinGW and Cygwin codes are targeted to run in an ordinary Windows
 * DOS box, _not_ within the MSYS/Cygwin UNIX emulation environments,
 * where the executables exhibit different behavior. Note that Cygwin
 * code compiled with the -mno-cygwin option exhibits its own
 * eccentricities, different from those of standalone Cygwin or
 * standalone MinGW.
 *
 * NOTE: The mathematical routines herein do not, in general,
 * attempt to treat all the exceptions and conditions spelled out
 * in documents such as IEEE754, IEEE854, and the C and C++
 * standards. These include detection, handling, and signalling
 * of overflow, underflow, denormals, infinities, and NaNs.
 *
 */
#define _CRT_SECURE_NO_WARNINGS

#if !defined(_TRN_H_)
#include "trn.h"
#endif

 /* M_EPSILON1 is the convergence tolerance in several calculations. */

#define M_EPSILON1 LDBL_EPSILON


/* The following external variables may be accessed from linked codes
   by means of global "extern" declarations. For the arrays, use,
   e.g., the global declaration "extern unsigned long ulPrime16[]". */

unsigned long ulDmax = 0;  /* tracks global max of Lucas-Selfridge |D| */
unsigned long ulPrime16[6545];  /* array of 16-bit primes < 65538 */
int iPrime16Initialized = 0;  /* Is ulPrime16 initialized? */
//int iNoDivWarnings = 0;  /* declare extern in the calling code and set to 1
//			 to suppress warnings of non-exact
//			 integer divisions */
//long double ldZ[66];  /* Zeta function values zeta(2..65) */

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
/* 80-bit long double floating point arithmetic routines for DJGPP,   */
/* MinGW, and Cygwin, using inline assembly (GNU-style AT&T syntax)   */
/* to compute the results in the x87 FPU.                             */
/*                                                                    */
/* NOTES:                                                             */
/*                                                                    */
/* (1) These functions are already available in Borland C, Digital    */
/*     Mars, and SUSE 10.x Linux.                                     */
/* (2) Borland C, Digital Mars, and MSVC cannot parse the GNU-style   */
/*     AT&T assembler syntax employed.                                */
/* (3) MSVC does not support 80-bit long doubles, period, so it is    */
/*     incompatible with such functions.                              */
/* (4) An x387 or later FPU is assumed to be present.                 */
/* (5) Function identifiers have been changed, more than once, to     */
/*     resolve clashes with various compiler library functions.       */
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
//long double __FABSL(long double ldArg)
//{
//	if (ldArg >= 0)return(ldArg);
//	return(-ldArg);
//}
/**********************************************************************/
/* The following compilers are syntactically incompatible with AT&T
   assembler syntax */
   /**********************************************************************/
#if !defined(__BORLANDC__) && !defined(__DMC__) && !defined(__MSVC__)
/**********************************************************************/
long double __CEILL(long double ldArg)
{
	unsigned short uhCWSave, uhCWTemp;
	long double ldResult;

	asm("fstcw %0" : "=m" (uhCWSave) : );
	uhCWTemp = uhCWSave & 0xFBFF;
	uhCWTemp = uhCWTemp | 0x0800;
	asm("fldcw %0" : : "m" (uhCWTemp));
	asm("frndint" : "=t" (ldResult) : "0" (ldArg));
	asm("fldcw %0" : : "m" (uhCWSave));
	return(ldResult);
}
/**********************************************************************/
long double __EXPL(long double ldExp)
{
	long double ldExp2, ldInt, ldFrac, ldMant, ldResult;

	ldExp2 = M_LOG_E_BASE2 * ldExp;
	ldInt = __FLOORL(ldExp2);
	ldFrac = ldExp2 - ldInt;
	asm("f2xm1" : "=t" (ldMant) : "0" (ldFrac));
	ldMant += 1.0L;
	asm("fscale" : "=t" (ldResult) : "0" (ldMant), "u" (ldInt));
	return(ldResult);
}
/**********************************************************************/
long double __FLOORL(long double ldArg)
{
	unsigned short uhCWSave, uhCWTemp;
	long double ldResult;

	asm("fnstcw %0" : "=m" (uhCWSave) : );
	uhCWTemp = uhCWSave & 0xF7FF;
	uhCWTemp = uhCWTemp | 0x0400;
	asm("fldcw %0" : : "m" (uhCWTemp));
	asm("frndint" : "=t" (ldResult) : "0" (ldArg));
	asm("fldcw %0" : : "m" (uhCWSave));
	return(ldResult);
}
/**********************************************************************/
long double __FMODL(long double ldTop, long double ldBottom)
{
	long double ldRem, ldNumerator;

	if (ldBottom == 0)
	{
		fprintf(stderr,
			"\n ERROR: Zero modulus passed to __FMODL.\n");
		signal(SIGFPE, SIG_DFL);
		raise(SIGFPE);
		exit(EXIT_FAILURE);
	}

	ldNumerator = ldTop;
	while (1)
	{
		asm("fprem" : "=t" (ldRem) : "0" (ldNumerator), "u" (ldBottom));
		if (__FABSL(ldRem) <= __FABSL(ldBottom))break;
		ldNumerator = ldRem;
	}
	return(ldRem);
}
/**********************************************************************/
long double __LOGL(long double ldArg)
{
	long double ldResult, ldLn_2 = M_LN2;

	if (ldArg <= 0)
	{
		fprintf(stderr,
			"\n ERROR: Non-positive argument passed to __LOGL.\n");
		signal(SIGFPE, SIG_DFL);
		raise(SIGFPE);
		exit(EXIT_FAILURE);
	}

	asm("fyl2x" : "=t" (ldResult) : "0" (ldArg), "u" (ldLn_2) : "st(1)");
	return(ldResult);
}
/**********************************************************************/
long double __LOG10L(long double ldArg)
{
	long double ldResult, ldLog10_2 = M_LOG_2_BASE10;

	if (ldArg <= 0)
	{
		fprintf(stderr,
			"\n ERROR: Non-positive argument passed to __LOG10L.\n");
		signal(SIGFPE, SIG_DFL);
		raise(SIGFPE);
		exit(EXIT_FAILURE);
	}

	asm("fyl2x" : "=t" (ldResult) : "0" (ldArg), "u" (ldLog10_2) : "st(1)");
	return(ldResult);
}
/**********************************************************************/
long double __LOG2L(long double ldArg)
{
	long double ldResult, ldOne = 1.0L;

	if (ldArg <= 0)
	{
		fprintf(stderr,
			"\n ERROR: Non-positive argument passed to __LOG2L.\n");
		signal(SIGFPE, SIG_DFL);
		raise(SIGFPE);
		exit(EXIT_FAILURE);
	}

	asm("fyl2x" : "=t" (ldResult) : "0" (ldArg), "u" (ldOne) : "st(1)");
	return(ldResult);
}
/**********************************************************************/
long double __POWL(long double ldBase, long double ldExp)
{
	long double ld2Exp, ldInt, ldFrac, ldMant, ldResult;

	if (ldBase <= 0)
	{
		fprintf(stderr,
			"\n ERROR: Non-positive base passed to __POWL.\n");
		signal(SIGFPE, SIG_DFL);
		raise(SIGFPE);
		exit(EXIT_FAILURE);
	}

	/* Evaluate as 2^(ldExp*log(ldBase,2)); do exponent expression first */

	asm("fyl2x" : "=t" (ld2Exp) : "0" (ldBase), "u" (ldExp) : "st(1)");

	/* Separate exponent result into integer and fractional parts */

	ldInt = __FLOORL(ld2Exp);
	ldFrac = ld2Exp - ldInt;
	asm("f2xm1" : "=t" (ldMant) : "0" (ldFrac));  /* 2^(fr part) - 1 */
	ldMant += 1.0L;  /* 2^(fr part) */

	/* Now multiply by 2^(integer part */

	asm("fscale" : "=t" (ldResult) : "0" (ldMant), "u" (ldInt));

	return(ldResult);
}
/**********************************************************************/
long double __SQRTL(long double ldArg)
{
	long double ldResult;

	if (ldArg < 0)
	{
		fprintf(stderr,
			"\n ERROR: Negative argument passed to __SQRTL.\n");
		signal(SIGFPE, SIG_DFL);
		raise(SIGFPE);
		exit(EXIT_FAILURE);
	}
	asm("fsqrt" : "=t" (ldResult) : "0" (ldArg));
	return(ldResult);
}
/**********************************************************************/
#endif  /* not Borland C and not DMC and not MSVC */
/**********************************************************************/
#ifndef __MSVC__
/**********************************************************************/
long double __NEARBYINTL(long double ldArg)
{
	/* The nearbyintl function of C99 is simulated for 80-bit x87 long
	   doubles. This is done for two reasons: (1) DJGPP 2.03 does not
	   support nearbyintl; (2) GCC 4.14, at least as implemented in SUSE
	   Linux 10.0, has a bug in nearbyintl which unpredictably interferes
	   with the cast of its value to an unsigned long long.

	   Thus, -4.5 and -3.5 must both be rounded to -4. No error trapping is
	   included. Assuming a mantissa of 64 bits (as in the x87 hardware),
	   the return will be unreliable for values of |ld| > 2^63. Be aware
	   also that some implementations of strtold (also _atold, printf, etc.)
	   have problems correctly converting the 19th and succeeding significant
	   decimal digits (if present). */

#if defined(__BORLANDC__) || defined(__DMC__)

	long double ld2, ldFrac, ldInt;

	ld2 = floorl(ldArg);
	ldFrac = modfl(ldArg, &ldInt);
	if (fabsl(ldFrac) == 0.5L)  /* Round to nearest even integer */
	{
		if (fabsl(fmodl(ld2, 2)) < 0.5L)  /* Is ld2 even? */
			return(ld2);
		else
			return(ld2 + 1);
	}
	return(floorl(ldArg + 0.5L));

#else

	unsigned short uhCWSave, uhCWTemp;
	long double ldResult;

	asm("fstcw %0" : "=m" (uhCWSave) : );  /* store FPU control word */
	uhCWTemp = uhCWSave & 0xF3FF;  /* clear rounding bits ==> nearest/even */
	asm("fldcw %0" : : "m" (uhCWTemp));  /* load new control word */
	asm("frndint" : "=t" (ldResult) : "0" (ldArg));  /* bankers' rounding */
	asm("fldcw %0" : : "m" (uhCWSave));  /* restore previous control word */
	return(ldResult);

#endif
}
/**********************************************************************/
long double __NEXTTOWARDL(long double ld, long double ldDirection)
{
	/* Returns the rational number ld1 with the following three properties:

	   (1) ld1 is representable exactly as a long double;
	   (2) ld1 lies between ld and ldDirection;
	   (3) No other such number lies between ld and ld1.

	   If no such number can be found, ld is returned.

	   Do _not_ use HUGE_VALL or similar manifest constants as ldDirection
	   (they might be defined as infinities or NaNs); however, LDBL_MAX and
	   -LDBL_MAX are OK.

	   Under MSVC, this function will behave like nexttoward (the double
	   precision version).
	*/

	int i, iSign;
	long double eps, ld1;

	if (ld == 0)
	{
		if (ldDirection > 0)return(LDBL_MIN);
		if (ldDirection < 0)return(-LDBL_MIN);
	}

	if (ldDirection > ld)iSign = +1;
	else if (ldDirection < ld)iSign = -1;
	else return(ld);

#ifdef __MSVC__
	eps = fabs(ld)*1e-17;
#else
	eps = fabsl(ld)*1e-20L;
#endif

	ld1 = ld;
	for (i = 1; i < 100; i++)
	{
		ld1 = ld + iSign * i*eps;
		if (ld1 != ld)return(ld1);
	}

	return(ld);
}
/**********************************************************************/
char *szLDtoHex(char *sz, long double ld)
{
	/*
	 * Converts a long double to a string containing the ten-byte (Intel
	 * x86/x87) hex representation. ASSUMPTIONS: short integers are 16 bits;
	 * long longs are 64 bits; long doubles are 80-bit format (IEEE
	 * 754 double extended), or 96 bits including 16 bits of padding.
	 *
	 * Returning the result as both a function argument and function value
	 * may appear redundant. However, returning the result from a local
	 * automatic variable is not allowed. Returning the result from a local
	 * static variable is allowed, but would produce unexpected results in
	 * a statement such as
	 *
	 * printf("\n 10==>%s   -10==>%s\n", szLDtoHex2(10), szLDtoHex2(-10));
	 *
	 * A single value (the result for 10) is printed twice! This will also
	 * occur if the present function is used and the same target string
	 * is referenced in both calls. Looks like a bug in printf to me---but
	 * that's the way it works in both DJGPP 2.03 and SUSE Linux 10.0.
	 *
	 */

	short h;
	int64_t *pll, ll;

	pll = (int64_t *)(&ld);
	ll = *(pll);
	h = *((short *)(++pll));
	sprintf(sz, "0x%04hX%016" PRIX64, h, ll);
	return(sz);
}
/**********************************************************************/
#endif  /* not MSVC */
/**********************************************************************/
//long double ldRand64(void)
//{
//	/* Generates a pseudorandom long double, between zero and one (exclusive),
//	   with 64 bits of precision. The initializer seed64 will be invoked on the
//	   first call. Note that this code is not yet suitable for a multithreaded
//	   or parallelized environment.
//
//	   If the compiler fails to support long doubles with a 64-bit mantissa
//	   (e.g., Microsoft Visual C++), the return value will be only double
//	   precision (53 bits).
//
//	   This algorithm has not been tested for rigorous conformance to the
//	   standards for pseudorandom number generation. For practical
//	   purposes it appears to perform acceptably.
//
//	   This 64-bit linear congruential generator is based on those
//	   developed by the Florida State SPRNG group, and by Beck and
//	   Brooks at UCRL/LLNL as part of an effort to improve Monte Carlo
//	   simulations in radiation and neutron transport studies. See
//	   the following references:
//
//	o <http://nuclear.llnl.gov/CNP/rng/rngman/node4.html> (accessed
//	  27 April 2010).
//	o "The RNG Random Number Library," Bret R. Beck and Eugene D. Brooks III
//	  (8 December 2000), University of California, Lawrence Livermore National
//	  Laboratory, Livermore, California 94550, UCRL-MA-141673. See
//	  <http://nuclear.llnl.gov/CNP/rng/rngman/node3.html>.
//	o M. Mascagni, S. A. Cuccaro, D. V. Pryor and M. L. Robinson (1995),
//	  "A Fast, High Quality, and Reproducible Parallel Lagged-Fibonacci
//	  Pseudorandom Number Generator," Journal of Computational Physics
//	  119:211-219.
//	o M. Mascagni (1999), "SPRNG: A Scalable Library for Pseudorandom Number
//	  Generation," to appear in Proceedings of the Third International
//	  Conference on Monte Carlo and Quasi Monte Carlo Methods in Scientific
//	  Computing, J. Spanier et al., editors, Springer Verlag, Berlin. See
//	  <http://sprng.cs.fsu.edu/RNG/>.
//	*/
//
//#define MULTIPLIER __ULL(2862933555777941757)
//#define ADDEND     __ULL(3037000493)
//
//	static uint64_t ull = 0;
//	//char sz[256];
//
//	while (ull == 0)ull = seed64();
//	ull = MULTIPLIER * ull + ADDEND;
//	return(ull / (1 + (long double)UINT64_MAX));
//
//#undef MULTIPLIER
//#undef ADDEND
//}
/**********************************************************************/
//static uint64_t seed64(void)
//{
//	/* Generates a pseudorandom 64-bit unsigned integer to be used as the
//	   seed for the 64-bit pseudorandom number generator ldRand64. Not
//	   intended for use by any other routine except ldRand64. */
//
//	long NB, BA;
//	uint64_t ull, RM;
//
//	/* If available, use the time stamp counter maintained on the Pentium
//	   and compatible chips. */
//
//#ifdef __DJ204__
//	return(_rdtsc());
//#elif defined(__WIN32__) && defined(WINVER) && (WINVER >= 0x0400)
//	   /* Thanks to Huang Yuanbing <bailuzhou(at)163(dot)com> for this
//		  code fragment, adapted from his PrimeNumber.cpp code (12 Feb 2009). */
//	static LARGE_INTEGER s_freq;
//	LARGE_INTEGER performanceCount;
//	if (s_freq.QuadPart && QueryPerformanceFrequency(&s_freq))
//	{
//		QueryPerformanceCounter(&performanceCount);
//		return(performanceCount.QuadPart);
//	}
//#elif defined(__WIN32__)  /* use GetTickCount() */
//	srand(time(NULL));
//	return((GetTickCount() + (uint64_t)rand())*((uint64_t)time(NULL) +
//		rand()) + rand());
//#endif
//
//	/* At this point neither Windows nor the Pentium time stamp is available. */
//
//#ifdef __LINUX__  /* BSD gettimeofday is present in most Linux versions */
//	static struct timeval tv;
//	gettimeofday(&tv, NULL);
//#endif
//
//	RM = (uint64_t)RAND_MAX;
//	NB = 1;
//
//	while (1)  /* Compute the number of bits NB in RAND_MAX */
//	{
//		RM = RM / 2;
//		if (!RM)break;
//		NB++;
//	}
//	srand((unsigned int)time(NULL));
//	ull = rand();
//	BA = NB;  /* Number of bits already assigned random values */
//	while (BA < 64 - NB)
//	{
//		ull = (ull << NB) + rand();
//		BA += NB;
//	}
//	ull = (ull << (64 - BA)) + (rand() >> (NB - 64 + BA)) + (uint64_t)time(NULL);
//#ifdef __LINUX__
//	ull += tv.tv_usec;  /* enhances the period of ldRand64 */
//#endif
//
//	return(ull);
//}
/**********************************************************************/
long double ___strtold(const char *sz, char **ppch)
{
	/* Convert a string to the equivalent long double value; specifically,
	   to the long double whose value is nearer to the string value than
	   any other long double. Intended to give greater and more consistent
	   precision than strtold, _strold, and _atold.

	   NOTES:

	   (1) It is assumed that the long doubles are stored in IEEE 754
		   extended double format (Intel x86/x87 ten-byte extended
		   or temporary reals). Otherwise, use strtold and hope for
		   the best. MSVC has no 80-bit long doubles, so simply return
		   the result of strtod.
	   (2) Some of the mpf's are redundant, but have been retained for
		   readability.
	   (3) Input strings containing syntax errors are handled by the
		   default (possibly less precise) strtold or strtod functions;
		   furthermore, if ppch is not NULL, its value is also set by the
		   default strtold (or strtod).
	   (4) The original function name __strtold was changed to ___strtold
		   to avoid a linker clash in Cygwin code compiled with the
		   -mno-cygwin option (to produce MinGW-like code). The linker
		   complains that an intrinsic routine also named __strtold exists,
		   although it appears to be dysfunctional, and I cannot find a
		   prototype for it; see the "hack" comment in trn.h, where a
		   #define should make the use of the old __strtold transparent---
		   except, of course, when -mno-cygwin is used.
	   (5) If GMP is not available, the highest precision native string-to-fp
		   function available is called. Note that the DEFINEs in trn.h
		   always make strtold an alias for this function, except in
		   Cygwin and MSVC, which have no such function for long doubles.
		   Furthermore, those DEFINEs, together with the IFNDEF __GMP__
		   block below, also allow ___strtold to be called unconditionally.
	   (6) MPFR is used if available.
	*/

	char ch; // *pch;
	unsigned short uh, *puh;
	int i, iSign, iExp, iRet;
	unsigned long ulMSL, ulLSL, *pMSL, *pLSL;
	//int64_t ll;
	uint64_t *pull; // ull;
	long double ld;

	static char szLocal[512];
	static int iFirst = 1;
#ifdef __GMP__
	static mpf_t mpf, mpf2, mpf3, mpf4, mpf5, mpfHalf;
#endif

#ifdef __MPFR__
	mpfr_t mpfr;
	mpfr_init2(mpfr, 128);
	mpfr_strtofr(mpfr, sz, ppch, 10, GMP_RNDN);
	ld = mpfr_get_ld(mpfr, GMP_RNDN);
	mpfr_clear(mpfr);
	return(ld);
#endif

	if (!sz || !(*sz))return(0.0L);

	/* Use the intrinsic functions to establish the default return values. */

#if defined(__CYGWIN__) || defined(__MSVC__)
	ld = (long double)strtod(sz, ppch);  /* Cygwin and MSVC have _no_ strtold */
#else
	ld = strtold(sz, ppch);  /* normal default for ppch and ld */
#endif

#if !defined(__GMP__) || defined(__MSVC__)
	return(ld);
#endif

#ifdef __GMP__

	if (iFirst)
	{
		mpf_init2(mpf, 512);
		mpf_init2(mpf2, 512);
		mpf_init2(mpf3, 512);
		mpf_init2(mpf4, 512);
		mpf_init2(mpf5, 512);
		mpf_init2(mpfHalf, 32);
		mpf_set_ui(mpfHalf, 1);
		mpf_div_2exp(mpfHalf, mpfHalf, 1);
		iFirst = 0;
	}

	strncpy(szLocal, sz, 511);  /* sz is const */
	szLocal[511] = 0;

	/* Isolate the token to be converted. We can't use strtok as it
	   may modify static buffers employed by the calling code; and
	   since the first argument sz is "const", we must make a local
	   copy szLocal for modification. Stop at the first space or
	   non-printing character. */

	szTrimMWS(szLocal);
	i = 0;
	for (i = 0; i < strlen(szLocal); i++)
	{
		ch = szLocal[i];
		if (!isgraph(ch))
		{
			szLocal[i] = 0;
			break;
		}
	}

	iRet = __mpf_set_str(mpf, szLocal, 10);
	iSign = mpf_sgn(mpf);
	if (iRet == -1 || iSign == 0)return(ld);  /* Handles input syntax errors */

	mpf_set(mpf2, mpf);
	mpf_abs(mpf2, mpf2);
	iExp = 0;
	if (mpf_cmp_ui(mpf2, 1) > 0)
	{
		while (mpf_cmp_ui(mpf2, 2) >= 0)
		{
			mpf_div_2exp(mpf2, mpf2, 1);
			iExp++;
		}
	}
	else
	{
		while (mpf_cmp_ui(mpf2, 1) < 0)
		{
			mpf_mul_2exp(mpf2, mpf2, 1);
			iExp--;
		}
	}
	mpf_mul_2exp(mpf2, mpf2, 63);
	mpf_add(mpf2, mpf2, mpfHalf);
	mpf_floor(mpf2, mpf2);

	/* mpf2 now contains the integer value of the 64-bit mantissa in
	   the long double representation of pchToken. */

	mpf_div_2exp(mpf3, mpf2, 32);
	mpf_floor(mpf3, mpf3);  /* mpf3 = high 64-32 bits of mpf2 */
	ulMSL = (unsigned long)mpf_get_ui(mpf3);  /* least significant 32 bits */
	mpf_mul_2exp(mpf4, mpf3, 32);  /* mpf2 with low 32 bits zeroed */
	mpf_sub(mpf5, mpf2, mpf4);  /* mpf5 = low 64-32 bits of mpf2 */
	ulLSL = (unsigned long)mpf_get_ui(mpf5);

	/* ulMSL and ulLSL are now the most significant and least significant
	   32-bit unsigned integers of the 64-bit integer mantissa. */

	pMSL = (unsigned long *)(&ld);
	pLSL = pMSL++;
	*pMSL = ulMSL;
	*pLSL = ulLSL;

	uh = iExp + 0x3FFF;  /* Exponent bias */
	if (iSign < 0)uh += 0x8000;  /* Sign bit incorporated */

	/* uh now contains the value of the 16-bit sign+biased exponent
	   field of the long double representation of pchToken. */

	pull = (uint64_t *)(&ld);
	puh = ((unsigned short *)(++pull));
	*puh = uh;

	return(ld);
#endif  /* GMP */
}
/**********************************************************************/
//char *__szLD(char *sz, long double ld, char *szFMT)
//{
//	/* Converts ld to an ASCII string, using the format specifier szFMT.
//	   If szFMT is NULL, empty, or lacks the signature 'L' for long doubles,
//	   a default format of "%.19Le" is used (or "%.15le" for MinGW, which
//	   is unable to process "L", and MSVC, which has no 80-bit long doubles).
//
//	   The ASCII string is stored in sz, which must have sufficient
//	   allocation. A pointer to the result is also returned.
//
//	   GMP, an mpf, and gmp_sprintf are used to process the long double
//	   argument. If GMP is not available, the compiler's native sprintf is
//	   called. MPFR is used if available.
//
//	   This routine should exhibit precision >= that of most compilers'
//	   sprintf or printf routines, and should be particularly useful in MinGW
//	   with GMP. It is of little value for MSVC.
//
//	   WARNING: Since the printing value is cast to a double (53-bit mantissa)
//	   for MSVC and MinGW, if |ld| is < DBL_MIN or > DBL_MAX, incorrect
//	   values may be returned on those platforms.
//
//	   NOTE: If the first 'L' or 'l' following the '%' character in an
//	   (invalid) format specifier is _not_ the size specifier, this routine
//	   may crash or return garbage; e.g., szFMT=="%.12e is the length".
//	*/
//
//	char szFMT2[128], *pch, *pch2;
//#ifdef __GMP__
//	static int iFirst = 1;
//	static mpf_t mpf;
//#endif
//
//	/* Edit the format string. Replace the 'L' following '%' with an 'F'
//	   for gmp_printf, or with an 'l' for MinGW without GMP or MSVC, or
//	   leave it alone for other platforms without GMP. If szFMT is NULL
//	   or empty, use the default format. */
//
//	if ((szFMT == NULL) || (*szFMT == 0))
//	{
//#ifdef __MSVC__
//		strcpy(szFMT2, "%.15le");
//#elif defined(__GMP__)
//		strcpy(szFMT2, "%.19Fe");
//#elif defined(__MINGW__)
//		strcpy(szFMT2, "%.15le");
//#else
//		strcpy(szFMT2, "%.19Le");
//#endif
//	}
//	else
//	{
//		strcpy(szFMT2, szFMT);
//		pch = strchr(szFMT2, '%');
//		if (pch)
//		{
//			pch2 = strchr(pch, 'L');
//			if (!pch2)pch2 = strchr(pch, 'l');
//			if (pch2)
//			{
//#ifdef __MSVC__
//				*pch2 = 'l';
//#elif defined(__GMP__)
//				*pch2 = 'F';
//#elif defined(__MINGW__)
//				*pch2 = 'l';
//#endif
//			}
//			else
//			{
//#ifdef __MSVC__
//				strcpy(szFMT2, "%.15le");
//#elif defined(__GMP__)
//				strcpy(szFMT2, "%.19Fe");
//#elif defined(__MINGW__)
//				strcpy(szFMT2, "%.15le");
//#else
//				strcpy(szFMT2, "%.19Le");
//#endif
//			}
//		}
//		else
//		{
//#ifdef __MSVC__
//			strcpy(szFMT2, "%.15le");
//#elif defined __GMP__
//			strcpy(szFMT2, "%.19Fe");
//#elif defined(__MINGW__)
//			strcpy(szFMT2, "%.15le");
//#else
//			strcpy(szFMT2, "%.19Le");
//#endif
//		}
//	}
//
//#ifdef __MPFR__
//	mpfr_t mpfr;
//	mpf_t mpfTemp;
//	mpfr_init2(mpfr, 128);
//	mpf_init2(mpf, 128);
//	mpfr_set_ld(mpfr, ld, GMP_RNDN);
//	mpfr_get_f(mpf, mpfr, GMP_RNDN);
//	gmp_sprintf(sz, szFMT2, mpf);
//	mpf_clear(mpf);
//	mpfr_clear(mpfr);
//	return(sz);
//#endif
//
//	/* Now print ld to sz using gmp_sprintf or sprintf and the edited format. */
//
//#if (defined(__MINGW__) && !defined(__GMP__)) \
//    || defined (__MSVC__)
//
///* MSVC and MinGW cannot print 80-bit long doubles correctly; edit for
//   overflow and cast to double. */
//
//	if (ld > DBL_MAX)
//		ld = DBL_MAX;
//	else if (ld < -DBL_MAX)
//		ld = -DBL_MAX;
//
//	sprintf(sz, szFMT2, (double)ld);  /* MSVC has no 80-bit ld's */
//
//#elif defined(__GMP__)
//
//	if (iFirst)  /* Initialize static mpf on the first invocation */
//	{
//		mpf_init(mpf);
//		iFirst = 0;
//	}
//	__mpf_set_ld(mpf, ld);
//	gmp_sprintf(sz, szFMT2, mpf);
//
//#else
//
//	sprintf(sz, szFMT, ld);
//
//#endif
//
//	return(sz);
//}
/**********************************************************************/
//char *__szLL(char *sz, int64_t ll)
//{
//	/* Converts ll to an ASCII string of minimum length. It is assumed
//	 * that ll is a 64-bit signed integer.
//	 *
//	 * The ASCII string is stored in sz, which must have sufficient
//	 * allocation. A pointer to the result is also returned.
//	 *
//	 * Targeted at printf-challenged environments, such as MinGW.
//	 *
//	 */
//
//	static char sz2[32];
//	unsigned long ul1, ul2, ul3;
//	uint64_t ull;
//
//	sz[0] = 0;
//	if (ll == 0)
//	{
//		strcpy(sz, "0");
//		return(sz);
//	}
//	if (ll < 0)
//	{
//		strcpy(sz, "-");
//		ll = -ll;
//	}
//
//	if (ll < 1e9)
//		sprintf(sz2, "%lu", (unsigned long)ll);
//	else if (ll < 1e18)
//	{
//		ul1 = (unsigned long)(ll / __ULL(1000000000));
//		ul2 = (unsigned long)(ll % __ULL(1000000000));
//		sprintf(sz2, "%lu%09lu", ul1, ul2);
//	}
//	else
//	{
//		ul1 = (unsigned long)(ll / __ULL(1000000000000000000));
//		ull = ll % __ULL(1000000000000000000);
//		ul2 = (unsigned long)(ull / __ULL(1000000000));
//		ul3 = (unsigned long)(ull % __ULL(1000000000));
//		sprintf(sz2, "%lu%09lu%09lu", ul1, ul2, ul3);
//	}
//	strcat(sz, sz2);
//
//	return(sz);
//}
/**********************************************************************/
//char *__szULL(char *sz, uint64_t ull)
//{
//	/* Converts ull to an ASCII string of minimum length. It is assumed
//	 * that ull is a 64-bit unsigned integer.
//	 *
//	 * The ASCII string is stored in sz, which must have sufficient
//	 * allocation. A pointer to the result is also returned.
//	 *
//	 * Targeted at printf-challenged environments, such as MinGW.
//	 *
//	 */
//
//	static char sz2[32];
//	unsigned long ul1, ul2, ul3;
//	uint64_t ull2;
//
//	sz[0] = 0;
//	if (ull == 0)
//	{
//		strcpy(sz, "0");
//		return(sz);
//	}
//
//	if (ull < 1e9)
//		sprintf(sz2, "%lu", (unsigned long)ull);
//	else if (ull < 1e18)
//	{
//		ul1 = (unsigned long)(ull / __ULL(1000000000));
//		ul2 = (unsigned long)(ull % __ULL(1000000000));
//		sprintf(sz2, "%lu%09lu", ul1, ul2);
//	}
//	else
//	{
//		ul1 = (unsigned long)(ull / __ULL(1000000000000000000));
//		ull2 = ull % __ULL(1000000000000000000);
//		ul2 = (unsigned long)(ull2 / __ULL(1000000000));
//		ul3 = (unsigned long)(ull2 % __ULL(1000000000));
//		sprintf(sz2, "%lu%09lu%09lu", ul1, ul2, ul3);
//	}
//	strcat(sz, sz2);
//
//	return(sz);
//}
/**********************************************************************/
/**********************************************************************/
#ifdef __GMP__
/**********************************************************************/
/**********************************************************************/
/*                    GMP mpz BIGINT functions                        */

/**********************************************************************/
/*                         GMP mpf functions                          */
/**********************************************************************/
/**********************************************************************/
int __mpf_set_str(mpf_t mpf, char *sz, int iBase)
{
	/* Fixes a bug (feature?) in the GNU GMP floating-point library (noted
	   in both version 4.01 and version 4.14). A leading plus sign in sz is
	   not properly recognized, causing an erroneous value of zero to be
	   assigned to mpf. The fix used is to check the first visible character
	   in sz; if it is a '+', replace it with a space. */

	char ch;
	int i, j = -1, iRet;

	for (i = 0; i < strlen(sz); i++)
	{
		ch = sz[i];
		if (ch == '+')
		{
			sz[i] = ' ';
			j = i;  /* save the change */
			break;
		}
		if (isgraph(ch))break; /* printing characters except space */
	}
	iRet = mpf_set_str(mpf, sz, iBase);
	if (j != -1)sz[i] = '+';  /* restore sz */
	return(iRet);
}
/***********************************************************************/

#endif  /* GMP available */
#ifdef __MPFR__
/**********************************************************************/
/**********************************************************************/
/*                         GMP + MPFR functions                       */
/**********************************************************************/
/**********************************************************************/
#undef LIT_LI2
#undef LIT_C2
#undef LIT_C3
#undef LIT_C4
#undef LIT_R2
#undef DEF_PREC_BITS
#undef MAX_PREC_BITS
#define LIT_LI2 "1.0451637801174927848445888891946131365226155781512015758329091440750132052103595301727174056263833563"
#define LIT_C2  "0.6601618158468695739278121100145557784326233602847334133194484233354056423044952771437600314138398679"
#define LIT_C3  "0.6351663546042712072066965912725224173420656873323724508997344604867846131161391882080291386764046176"
#define LIT_C4  "0.3074948787583270931233544860710768530221785199506639282983083962608887672966929948138402646817149384"
#define LIT_R2  "1.5410090161871318832885037866275465435308992182709284694298743638692436905783372048127122002248898594"
#define DEF_PREC_BITS  128  /* 128 bits or 38 SDD */
#define MAX_PREC_BITS  320  /* 320 bits or 96 SDD */
/**********************************************************************/
/**********************************************************************/
void mpfrLIRZ(mpfr_t mpfrx, mpfr_t mpfrLi, mpfr_t mpfrHL2,
	mpfr_t mpfrHL3, mpfr_t mpfrHL4, mpfr_t mpfrR, int iLocalPrec)
{
	/* The logarithmic integral expressions Li(x), HL2(x), HL3(x), and HL4(x),
	   approximating pi(x), pi_2(x), pi_3(x), and pi_4(x) respectively---the
	   counts from 0 to x of the primes, twin-prime pairs, basic prime triplets,
	   and basic prime quadruplets---are approximated using a series obtained as
	   explained below.

	   In addition, Riemann's prime counting function R(x) is approximated,
	   using a truncated Gram series; see "Prime numbers and computer methods
	   for factorization," Hans Reisel (Birkhauser, Boston, 1994), pp. 50-51,
	   especially eqn 2.27.

	   NOTE: The domain of the algorithm is x > 1; the accuracy degrades
	   near x=1, a singular point of these functions; and these functions
	   are rarely called with arguments x < 2. Consequently, for x < 2,
	   this routine returns artificial values of zero.

	   The default precision is 128 bits (iLocalPrec=128) or 36 significant
	   decimal digits (after typical rounding errors). The routine can be
	   called with values of iLocalPrec as high as 320, for 96 significant
	   decimal digits, or with iLocalPrec as low as 32 (9 SDD). Higher
	   precisions could be obtained by a slight modification of the
	   code---consisting mainly of providing more precise values for the
	   literal constants LIT_LI2, LIT_C2, LIT_C3, LIT_C4, and LIT_R2, and
	   changing the #define for MAX_PREC_BITS.

	   In general, one can expect cumulative rounding errors of two or
	   three SDD. Thus 128 bits will actually yield 35 or 36 correct
	   decimal digits.

	   The fact that Li(x) is asymptotic to pi(x) is one statement of the
	   Prime Number Theorem. The use of L2(x), L3(x), and L4(x) as
	   approximations to pi_2(x), pi_3(x), and pi_4(x) is a consequence
	   of the prime k-tuples conjecture of Hardy and Littlewood (ca. 1922).

	   The technique is as follows. In I4=int(1/(ln t)^4, t) substitute
	   u=ln(t), yielding the integral int(exp(u)/u^4, u). Substitute the
	   Maclaurin series for exp(u). Integrate the first five terms separately
	   to yield

	   I4=-1/(3u^3) - 1/(2u^2) - 1/(2u) + (ln u)/6 + u/24
											  + sum(u^k/(k*(k+3)!), k, 2, inf).

	   Replace u by ln(t) and apply the limits t=2 to t=x to produce the
	   result for the integral in L4.  Note that the terms in the resulting
	   series are related by

	   T(k+1)=T(k)*(ln t)*k/((k+1)(k+4).

	   Iterate the series until the ratio of successive terms is < M_EPSILON1
	   (LDBL_EPSILON/4, approximately 2.71e-20 on an x386 system).

	   Once I4 is evaluated, I3, I2, and I1 can be obtained using (reverse)
	   integration by parts on I1(t)=int(1/(ln t), t):

	   I1(t)=t/(ln t) + I2(t)=t/(ln t) + t/(ln t)^2 + 2*I3(t)
			=t/(ln t) + t/(ln t)^2 + 2t/(ln t)^3 + 6*I4(t) ,
	   or

	   I3(t)=t/(ln t)^3 + 3*I4(t)

	   I2(t)=t/(ln t)^2 + 2*I3(t)

	   I1(t)=t/(ln t) + I2(t).

	   Now apply the limits 2 to x. Add ldLi2 to L1 to account for
	   the lower limit being 0 rather than 2. Multiply I4, I3, and I2
	   by the appropriate Hardy-Littlewood constants to obtain the
	   estimates HL4(x), HL3(x), and HL2(x)for pi_4(x), pi_3(x),
	   and pi_2(x).

	   R(x), the Riemann prime-counting function, is calculated using
	   a different algorithm (truncated Gram series); see the discussion
	   in the routine ldRPCF.
	*/

	static int iFirst = 1;
	static mpf_t mpf;
	static mpfr_t mpfrLogx, mpfrLogx2, mpfrLogx3, mpfrLI2, mpfrLN2, mpfrLN2_2,
		mpfrLN2_3, mpfr1_LN2, mpfr1_LN2_2, mpfr1_LN2_3, mpfrLN2_24, mpfrC2,
		mpfrC3, mpfrC4, mpfrI2, mpfrI3, mpfrI4, mpfr1, mpfr2, mpfrTerm1,
		mpfrTerm2, mpfrDelta, mpfrEpsilon, mpfr3, mpfr4, mpfr5, mpfrt,
		mpfrFactor, mpfrZeta;
	int iGlobalPrec;
	unsigned long ul;

	iGlobalPrec = mpfr_get_default_prec();
	if (iLocalPrec < 32)iLocalPrec = DEF_PREC_BITS;
	if (iLocalPrec > MAX_PREC_BITS)iLocalPrec = MAX_PREC_BITS;
	mpfr_set_default_prec(iLocalPrec);

	if (iFirst)
	{
		mpfr_init_set_str(mpfrLI2, LIT_LI2, 10, GMP_RNDN);
		mpfr_init(mpfrLN2);
		mpfr_const_log2(mpfrLN2, GMP_RNDN);
		mpfr_init(mpfrLN2_24);
		mpfr_div_ui(mpfrLN2_24, mpfrLN2, 24, GMP_RNDN);
		mpfr_init(mpfrLN2_2);
		mpfr_sqr(mpfrLN2_2, mpfrLN2, GMP_RNDN);
		mpfr_init(mpfrLN2_3);
		mpfr_mul(mpfrLN2_3, mpfrLN2_2, mpfrLN2, GMP_RNDN);
		mpfr_init(mpfr1_LN2);
		mpfr_ui_div(mpfr1_LN2, 1, mpfrLN2, GMP_RNDN);
		mpfr_init(mpfr1_LN2_2);
		mpfr_sqr(mpfr1_LN2_2, mpfr1_LN2, GMP_RNDN);
		mpfr_init(mpfr1_LN2_3);
		mpfr_mul(mpfr1_LN2_3, mpfr1_LN2_2, mpfr1_LN2, GMP_RNDN);
		mpfr_init_set_str(mpfrC2, LIT_C2, 10, GMP_RNDN);
		mpfr_init_set_str(mpfrC3, LIT_C3, 10, GMP_RNDN);
		mpfr_init_set_str(mpfrC4, LIT_C4, 10, GMP_RNDN);
		mpfr_init_set_ui(mpfrI2, 0, GMP_RNDN);
		mpfr_init_set_ui(mpfrI3, 0, GMP_RNDN);
		mpfr_init_set_ui(mpfrI4, 0, GMP_RNDN);
		mpfr_init_set_ui(mpfrLogx, 0, GMP_RNDN);
		mpfr_init_set_ui(mpfrLogx2, 0, GMP_RNDN);
		mpfr_init_set_ui(mpfrLogx3, 0, GMP_RNDN);
		mpfr_init_set_ui(mpfr1, 0, GMP_RNDN);
		mpfr_init_set_ui(mpfr2, 0, GMP_RNDN);
		mpfr_init_set_ui(mpfr3, 0, GMP_RNDN);
		mpfr_init_set_ui(mpfr4, 0, GMP_RNDN);
		mpfr_init_set_ui(mpfr5, 0, GMP_RNDN);
		mpfr_init_set_ui(mpfrTerm1, 0, GMP_RNDN);
		mpfr_init_set_ui(mpfrTerm2, 0, GMP_RNDN);
		mpfr_init_set_ui(mpfrDelta, 0, GMP_RNDN);
		mpfr_init_set_ui(mpfrt, 0, GMP_RNDN);
		mpfr_init_set_ui(mpfrFactor, 0, GMP_RNDN);
		mpfr_init_set_ui(mpfrZeta, 0, GMP_RNDN);
		iFirst = 0;
	}

	mpfr_init_set_ld(mpfrEpsilon, powl(2, -iLocalPrec - 3), GMP_RNDN);
	mpf_init2(mpf, iLocalPrec);
	mpfr_set_str(mpfrLi, LIT_LI2, 10, GMP_RNDN);
	mpfr_set_ui(mpfrHL2, 0, GMP_RNDN);
	mpfr_set_ui(mpfrHL3, 0, GMP_RNDN);
	mpfr_set_ui(mpfrHL4, 0, GMP_RNDN);
	mpfr_set_str(mpfrR, LIT_R2, 10, GMP_RNDN);

	if (mpfr_cmp_ui(mpfrx, 2) == 0)return;
	if (mpfr_cmp_ui(mpfrx, 2) < 0)
	{
		mpfr_set_ui(mpfrLi, 0, GMP_RNDN);
		mpfr_set_ui(mpfrR, 0, GMP_RNDN);
		return;
	}

	mpfr_log(mpfrLogx, mpfrx, GMP_RNDN);
	mpfr_sqr(mpfrLogx2, mpfrLogx, GMP_RNDN);
	mpfr_mul(mpfrLogx3, mpfrLogx2, mpfrLogx, GMP_RNDN);

	/* Initialize I4 */

	mpfr_mul_ui(mpfrI4, mpfrLogx3, 3, GMP_RNDN);
	mpfr_si_div(mpfrI4, -1, mpfrI4, GMP_RNDN);       //  -1/(3*ldLx3)

	mpfr_div_ui(mpfr1, mpfr1_LN2_3, 3, GMP_RNDN);
	mpfr_add(mpfrI4, mpfrI4, mpfr1, GMP_RNDN);      //  +1/(3*M_LN2_CUBED)

	mpfr_mul_2ui(mpfr1, mpfrLogx2, 1, GMP_RNDN);
	mpfr_ui_div(mpfr1, 1, mpfr1, GMP_RNDN);
	mpfr_sub(mpfrI4, mpfrI4, mpfr1, GMP_RNDN);      //  -1/(2*ldLx2)

	mpfr_mul_2ui(mpfr1, mpfrLN2_2, 1, GMP_RNDN);
	mpfr_ui_div(mpfr1, 1, mpfr1, GMP_RNDN);
	mpfr_add(mpfrI4, mpfrI4, mpfr1, GMP_RNDN);      //  +1/(2*M_LN2_SQUARED)

	mpfr_mul_2ui(mpfr1, mpfrLogx, 1, GMP_RNDN);
	mpfr_ui_div(mpfr1, 1, mpfr1, GMP_RNDN);
	mpfr_sub(mpfrI4, mpfrI4, mpfr1, GMP_RNDN);      //  -1/(2*ldLx)

	mpfr_mul_2ui(mpfr1, mpfrLN2, 1, GMP_RNDN);
	mpfr_ui_div(mpfr1, 1, mpfr1, GMP_RNDN);
	mpfr_add(mpfrI4, mpfrI4, mpfr1, GMP_RNDN);      //   +1/(2*M_LN2)

	mpfr_div(mpfr1, mpfrLogx, mpfrLN2, GMP_RNDN);
	mpfr_log(mpfr1, mpfr1, GMP_RNDN);
	mpfr_div_ui(mpfr1, mpfr1, 6, GMP_RNDN);         //
	mpfr_add(mpfrI4, mpfrI4, mpfr1, GMP_RNDN);      //   +logl(ldLx/M_LN2)/6

	/* Initialize term1 and term2 */

	mpfr_div_ui(mpfrTerm1, mpfrLogx, 24, GMP_RNDN);
	mpfr_set(mpfrTerm2, mpfrLN2_24, GMP_RNDN);

	for (ul = 1; ; ul++)
	{
		mpfr_add(mpfrI4, mpfrI4, mpfrTerm1, GMP_RNDN);
		mpfr_sub(mpfrI4, mpfrI4, mpfrTerm2, GMP_RNDN);
		mpfr_mul(mpfrDelta, mpfrEpsilon, mpfrI4, GMP_RNDN);
		mpfr_abs(mpfrDelta, mpfrDelta, GMP_RNDN);
		mpfr_set_ui(mpfr1, ul + 1, GMP_RNDN);
		mpfr_mul_ui(mpfr1, mpfr1, ul + 4, GMP_RNDN);
		mpfr_ui_div(mpfr1, ul, mpfr1, GMP_RNDN);
		mpfr_mul(mpfrTerm1, mpfrTerm1, mpfr1, GMP_RNDN);
		mpfr_mul(mpfrTerm1, mpfrTerm1, mpfrLogx, GMP_RNDN);
		if (mpfr_cmp(mpfrTerm1, mpfrDelta) < 0)break;
		if (mpfr_cmp(mpfrTerm2, mpfrDelta) > 0)
		{
			mpfr_mul(mpfrTerm2, mpfrTerm2, mpfrLN2, GMP_RNDN);
			mpfr_mul(mpfrTerm2, mpfrTerm2, mpfr1, GMP_RNDN);
		}
		else
			mpfr_set_ui(mpfrTerm2, 0, GMP_RNDN);
	}

	mpfr_mul(mpfrHL4, mpfrI4, mpfrC4, GMP_RNDN);
	mpfr_mul_ui(mpfrHL4, mpfrHL4, 27, GMP_RNDN);
	mpfr_div_2ui(mpfrHL4, mpfrHL4, 1, GMP_RNDN);  // HL4 = 27/2*c_4*I4

	mpfr_mul_ui(mpfrI3, mpfrI4, 3, GMP_RNDN);  // I3 = 3*I4 + x/(ln x)^3 - 2/(ln 2)^3
	mpfr_div(mpfr1, mpfrx, mpfrLogx3, GMP_RNDN);
	mpfr_add(mpfrI3, mpfrI3, mpfr1, GMP_RNDN);
	mpfr_mul_2ui(mpfr1, mpfr1_LN2_3, 1, GMP_RNDN);
	mpfr_sub(mpfrI3, mpfrI3, mpfr1, GMP_RNDN);
	mpfr_mul(mpfrHL3, mpfrI3, mpfrC3, GMP_RNDN);
	mpfr_mul_ui(mpfrHL3, mpfrHL3, 9, GMP_RNDN);
	mpfr_div_2ui(mpfrHL3, mpfrHL3, 1, GMP_RNDN);  // HL3 = 9/2*c_3*I3

	mpfr_mul_2ui(mpfrI2, mpfrI3, 1, GMP_RNDN);  // I2 = 2*I3 + x/(ln x)^2 - 2/(ln 2)^2
	mpfr_div(mpfr1, mpfrx, mpfrLogx2, GMP_RNDN);
	mpfr_add(mpfrI2, mpfrI2, mpfr1, GMP_RNDN);
	mpfr_mul_2ui(mpfr1, mpfr1_LN2_2, 1, GMP_RNDN);
	mpfr_sub(mpfrI2, mpfrI2, mpfr1, GMP_RNDN);
	mpfr_mul(mpfrHL2, mpfrI2, mpfrC2, GMP_RNDN);
	mpfr_mul_2ui(mpfrHL2, mpfrHL2, 1, GMP_RNDN);  // HL2 = 2*c_2*I2

	mpfr_div(mpfrLi, mpfrx, mpfrLogx, GMP_RNDN);  //  Li = I2 + x/(ln x) - 2/(ln 2) +
	mpfr_add(mpfrLi, mpfrLi, mpfrI2, GMP_RNDN);   //          + Li(2)
	mpfr_mul_2ui(mpfr1, mpfr1_LN2, 1, GMP_RNDN);
	mpfr_sub(mpfrLi, mpfrLi, mpfr1, GMP_RNDN);
	mpfr_add(mpfrLi, mpfrLi, mpfrLI2, GMP_RNDN);

	/* Compute the Riemann prime-counting function R(x), using a truncated
	   Gram series. See "Prime numbers and computer methods for
	   factorization," Hans Reisel (Birkhauser, Boston, 1994), pp. 50-51,
	   especially eqn 2.27. */

	mpfr_set(mpfrt, mpfrLogx, GMP_RNDN);
	mpfr_set(mpfrFactor, mpfrLogx, GMP_RNDN);
	mpfr_set_ui(mpfr1, 2, GMP_RNDN);
	mpfr_zeta(mpfrZeta, mpfr1, GMP_RNDN);
	mpfr_div(mpfrR, mpfrFactor, mpfrZeta, GMP_RNDN);
	for (ul = 2; ; ul++)
	{
		mpfr_mul(mpfrFactor, mpfrFactor, mpfrt, GMP_RNDN);
		mpfr_mul_ui(mpfrFactor, mpfrFactor, ul - 1, GMP_RNDN);
		mpfr_set_ui(mpfr1, ul, GMP_RNDN);
		mpfr_sqr(mpfr1, mpfr1, GMP_RNDN);
		mpfr_div(mpfrFactor, mpfrFactor, mpfr1, GMP_RNDN);
		mpfr_set_ui(mpfr1, ul + 1, GMP_RNDN);
		mpfr_zeta(mpfrZeta, mpfr1, GMP_RNDN);
		mpfr_div(mpfrTerm1, mpfrFactor, mpfrZeta, GMP_RNDN);
		mpfr_add(mpfrR, mpfrR, mpfrTerm1, GMP_RNDN);
		mpfr_div(mpfr1, mpfrTerm1, mpfrR, GMP_RNDN);
		if (mpfr_cmp(mpfr1, mpfrEpsilon) < 0)break;
	}
	mpfr_add_ui(mpfrR, mpfrR, 1, GMP_RNDN);

	mpfr_set_default_prec(iGlobalPrec);
	return;
}
/**********************************************************************/
void mpfrRRF(mpfr_t mpfrOut, mpfr_t mpfrIn, int iLocalPrec)
{
	/* Compute the Riemann prime-counting function R(x), also known as
	   Riemann's R-function, using a truncated Gram series. See "Prime
	   numbers and computer methods for factorization," Hans Reisel
	   (Birkhauser, Boston, 1994), pp. 50-51, especially eqn 2.27. */

	static int iFirst = 1;
	static mpfr_t mpfrLogx, mpfr1, mpfrTerm1, mpfrEpsilon, mpfrt,
		mpfrFactor, mpfrZeta;
	int iGlobalPrec;
	unsigned long ul;

	iGlobalPrec = mpfr_get_default_prec();
	if (iLocalPrec < 32)iLocalPrec = DEF_PREC_BITS;
	if (iLocalPrec > MAX_PREC_BITS)iLocalPrec = MAX_PREC_BITS;
	mpfr_set_default_prec(iLocalPrec);

	if (iFirst)
	{
		mpfr_init_set_ui(mpfrLogx, 0, GMP_RNDN);
		mpfr_init_set_ui(mpfr1, 0, GMP_RNDN);
		mpfr_init_set_ui(mpfrTerm1, 0, GMP_RNDN);
		mpfr_init_set_ui(mpfrt, 0, GMP_RNDN);
		mpfr_init_set_ui(mpfrFactor, 0, GMP_RNDN);
		mpfr_init_set_ui(mpfrZeta, 0, GMP_RNDN);
		mpfr_init_set_ld(mpfrEpsilon, powl(2, -iLocalPrec - 3), GMP_RNDN);
		iFirst = 0;
	}
	mpfr_set_str(mpfrOut, LIT_R2, 10, GMP_RNDN);

	if (mpfr_cmp_ui(mpfrIn, 2) == 0)return;
	if (mpfr_cmp_ui(mpfrIn, 2) < 0)
	{
		mpfr_set_ui(mpfrOut, 0, GMP_RNDN);
		mpfr_set_default_prec(iGlobalPrec);
		return;
	}

	mpfr_log(mpfrLogx, mpfrIn, GMP_RNDN);
	mpfr_set(mpfrt, mpfrLogx, GMP_RNDN);
	mpfr_set(mpfrFactor, mpfrLogx, GMP_RNDN);
	mpfr_set_ui(mpfr1, 2, GMP_RNDN);
	mpfr_zeta(mpfrZeta, mpfr1, GMP_RNDN);
	mpfr_div(mpfrOut, mpfrFactor, mpfrZeta, GMP_RNDN);
	for (ul = 2; ; ul++)
	{
		mpfr_mul(mpfrFactor, mpfrFactor, mpfrt, GMP_RNDN);
		mpfr_mul_ui(mpfrFactor, mpfrFactor, ul - 1, GMP_RNDN);
		mpfr_set_ui(mpfr1, ul, GMP_RNDN);
		mpfr_sqr(mpfr1, mpfr1, GMP_RNDN);
		mpfr_div(mpfrFactor, mpfrFactor, mpfr1, GMP_RNDN);
		mpfr_set_ui(mpfr1, ul + 1, GMP_RNDN);
		mpfr_zeta(mpfrZeta, mpfr1, GMP_RNDN);
		mpfr_div(mpfrTerm1, mpfrFactor, mpfrZeta, GMP_RNDN);
		mpfr_add(mpfrOut, mpfrOut, mpfrTerm1, GMP_RNDN);
		mpfr_div(mpfr1, mpfrTerm1, mpfrOut, GMP_RNDN);
		if (mpfr_cmp(mpfr1, mpfrEpsilon) < 0)break;
	}
	mpfr_add_ui(mpfrOut, mpfrOut, 1, GMP_RNDN);

	mpfr_set_default_prec(iGlobalPrec);
	return;
}
/**********************************************************************/
/**********************************************************************/
#endif  /* MPFR available */
/**********************************************************************/
/**********************************************************************/
/*              Prime number generation and testing                   */
/**********************************************************************/
/**********************************************************************/
//uint64_t ullGCD(uint64_t N1, uint64_t N2)
//{
//	/* Return the GCD of two unsigned 64-bit integers (Euclidean algorithm).
//	   The GMP documentation states that the binary algorithm is faster,
//	   but timing comparisons contradict this, at least for 64-bit integers. */
//
//	uint64_t R;
//
//	if (N1*N2 == 0)return(0);
//	while (1)
//	{
//		R = N1 % N2;
//		if (R == 0)break;
//		N1 = N2;  /* a <-- b      */
//		N2 = R;   /* b <-- a - qb */
//	}
//	return(N2);
//}
/**********************************************************************/
void vGenPrimesDiv(unsigned long *ulPrime, unsigned long *nPrimes,
	unsigned long *ulUB)         /* revised version of vGenPrimes */
{
	/* Using only trial prime divisors to the square root of N, generate all
	   the primes <= *ulUB. If this would be < *nPrimes, continue to the
	   (*nPrimes)th prime. The values of the primes are stored in the
	   array ulPrime, with markers ulPrime[0]=ulPrime[*nPrimes + 1]=0.

	   Before this function is called, the array ulPrime must have been
	   previously declared and allocated with sufficient storage for this
	   number of primes, plus two additional elements used to mark the
	   beginning and end of the array. Thus at least 4*(pi(*ulUB) + 2) or
	   4*(2 + *nPrimes) bytes must have been allocated in the calling code;
	   this can be accomplished by malloc(4*(2 + *nPrimes)) and/or
	   malloc(4*(2 + ceil(ldLi(*ulUB)))) in the calling code.

	   The return values are *nPrimes=(the actual number of primes generated)
	   and *ulUB=(the actual value of the largest, or (*nPrime)th, prime
	   generated).

	   *ulUB may not exceed the largest 32-bit prime, UINT32_MAX - 4,
	   and *nPrimes may not exceed the number of 32-bit primes (203280221).
	   As a minimum, the primes 2, 3, and 5 are generated, requiring five
	   elements in ulPrime (20 bytes).

	   Note that this function differs from the old vGenPrimes and the
	   new vGenPrimesSieve in that the lower bounds of *nPrimes=6543 and
	   *ulUB=65537 are no longer present.
	*/

	unsigned long ulN, ul, ul2, ulSqrt, ulDivisor;

	ulPrime[0] = 0;  /* 0th prime---just a marker */
	ulPrime[1] = 2;  /* 1st prime */
	ulPrime[2] = 3;  /* 2nd prime */
	ulPrime[3] = 5;  /* 3rd prime */
	if (*ulUB < 7)*ulUB = 5;
	if (*nPrimes < 3)*nPrimes = 3;
	if (*ulUB > UINT32_MAX - 4)*ulUB = UINT32_MAX - 4;  /* Largest 32-bit prime */
	if (*nPrimes > NUM_32BIT_PRIMES)*nPrimes = NUM_32BIT_PRIMES;
	ul = 4;
	ulN = 7;
	while (1)
	{
		ul2 = 2;
		ulSqrt = ceil(sqrt(ulN + 0.5));  /* faster than integer sqrt */
		while (1)
		{
			ulDivisor = ulPrime[ul2++];
			if (ulDivisor >= ulSqrt)
			{
				ulPrime[ul++] = ulN;
				break;
			}
			if (ulN % ulDivisor == 0)break;
		}
		ulN += 2;
		if ((ulN > *ulUB) && (ul > *nPrimes))break;  /* post-increments */
	}
	ulPrime[ul] = 0;  /* end marker */
	*nPrimes = ul - 1;
	*ulUB = ulPrime[ul - 1];

	return;
}
/**********************************************************************/
void vGenPrimes16(void)
{
	/* Generate from scratch all the primes < 2^16 + 2. There are 6543 such
	   primes, having a sum of 202353624. The array of primes ulPrime16 is
	   global (static) to this code, and must be accessed from the calling
	   code by means of the global declaration "extern unsigned long
	   ulPrime16[]". Do not allocate it or otherwise declare ulPrime16
	   in the calling code. */

	unsigned long nPrimes, ulUB;

	if (iPrime16Initialized)return;
	nPrimes = 6543;
	ulUB = 65537UL;
	vGenPrimesDiv(ulPrime16, &nPrimes, &ulUB);
	if (ulPrime16[6543] == 65537UL)
		iPrime16Initialized = 1;
	else
	{
		fprintf(stderr,
			"\n ***FATAL ERROR: ulPrime16 initialization failed in "
			"vGenPrimes16().\n");
		exit(EXIT_FAILURE);
	}

	return;
}
/**********************************************************************/
//void vGenPrimesSieve(unsigned long *ulPrime, unsigned long *nPrimes,
//	unsigned long *ulUB)
//{
//	/* Using the Sieve of Eratosthenes, generate all the primes <= *ulUB.
//	   If this would be < *nPrimes, continue to the (*nPrimes)th prime.
//	   Before this function is called, the array ulPrime must have been
//	   previously declared and allocated with sufficient storage for this
//	   number of primes, plus two additional elements used to mark the
//	   beginning and end of the array. Thus at least 4*(pi(*ulUB) + 2) or
//	   4*(2 + *nPrimes) bytes must have been allocated in the calling code;
//	   this can be accomplished by malloc(4*(2 + *nPrimes)) and/or
//	   malloc(4*(2 + ceil(ldLi(*ulUB)))) in the calling code.
//
//	   The return values are *nPrimes=(the actual number of primes generated)
//	   and *ulUB=(the actual value of the largest, or (*nPrime)th, prime
//	   generated).
//
//	   *ulUB may not exceed the largest 32-bit prime, UINT32_MAX - 4,
//	   and *nPrimes may not exceed the number of 32-bit primes (203280221).
//	   Minimum values are also enforced: *ulUB >= 65537 and *nPrimes >= 6543.
//
//	   This algorithm bootstraps some base primes by first calling
//	   vGenPrimesDiv, which does not use sieving, to establish a base set.
//	   An array allocated for use as a sieve is then used to generate
//	   succeeding blocks of primes by assigning successive odd integers to
//	   the elements of the sieving array. The base primes, bootstrapped by
//	   use of trial divisors, are required for sieving the first block. Each
//	   sieving then generates more than enough new primes to sieve the
//	   succeeding block.
//	*/
//
//	unsigned char *uchSieve;
//	unsigned long nPrimes0, ulLB0, ulUB0, ulPMA, ulSpan, ulOffset,
//		ulOffsetMax, ulN, ulSieveBytes, ulLB0In, ulUB0In, nP;
//	long double ldSB, ldSB1, ldSB2, lognP;
//
//	/* Adjust the bounds. */
//
//	if (*ulUB > UINT32_MAX - 4)*ulUB = UINT32_MAX - 4;  /* Largest 32-bit prime */
//	if (*nPrimes > NUM_32BIT_PRIMES)*nPrimes = NUM_32BIT_PRIMES;
//	if (*ulUB < 65537UL)*ulUB = 65537UL;
//	if (*nPrimes < 6543)*nPrimes = 6543;
//
//	/* Allocate memory for the sieving array. The size in bytes needs to
//	   be half the size of the interval to be sieved; a smaller size
//	   will cause iterations of the sieving process. */
//
//	ulPMA = __ulPhysicalMemoryAvailable() & 0xFFFFFFFE;  /* ulPMA is even */
//	ulSieveBytes = ulPMA / 2;  /* Use no more than half the available memory */
//	/* Use even less if this exceeds half the projected upper bound */
//	ldSB1 = (*ulUB + 1) / 2 + 1;
//	/* The following formulas, due to J. B. Rosser (1983) and G. Robin (1983),
//	   provide upper bounds for the (*nPrimes)th prime. See "The new book
//	   of prime number records," Paulo Ribenboim (Springer, 1995). */
//	nP = *nPrimes;
//	if (nP > 6543)
//	{
//		lognP = logl(nP);
//		if (nP < 7022)
//			ldSB2 = nP * (lognP + logl(lognP) + 8) / 2 + 2;
//		else
//			ldSB2 = nP * (lognP + logl(lognP) - 0.9385) / 2 + 2;
//	}
//	if (ldSB2 > ldSB1)ldSB = ldSB2; else ldSB = ldSB1;
//	if (ldSB < ulSieveBytes)ulSieveBytes = ceill(ldSB);
//	uchSieve = (unsigned char *)malloc(ulSieveBytes);
//	while (!uchSieve)
//	{
//		ulSieveBytes = ulSieveBytes / 2;
//		if (ulSieveBytes < 32000)
//		{
//			fprintf(stderr,
//				"\n ERROR: uchSieve malloc failed in vGenPrimesSieve.\n");
//			exit(EXIT_FAILURE);
//		}
//		uchSieve = (unsigned char *)malloc(ulSieveBytes);
//	}
//	if (!uchSieve)
//	{
//		fprintf(stderr, "\n ERROR: uchSieve malloc failed in vGenPrimesSieve.\n");
//		exit(EXIT_FAILURE);
//	}
//	ulSpan = 2 * ulSieveBytes - 2;  /* ulSpan is even */
//
//	/* Bootstrap the primes necessary for the initial sieve (to at least
//	   the next prime beyond the square root of the first sieve upper bound). */
//
//	nPrimes0 = 0;
//	ulUB0 = 1442 + ceil(sqrt(3 + 2 * ulSieveBytes));  /* No gap > 1442 below 1e18 */
//	ulUB0 = ulUB0 + (ulUB0 & 1);
//	vGenPrimesDiv(ulPrime, &nPrimes0, &ulUB0);
//
//	if ((nPrimes0 >= *nPrimes) && (ulUB0 >= *ulUB))
//	{
//		*nPrimes = nPrimes0;
//		*ulUB = ulUB0;
//		goto DONE_vGP;
//	}
//
//	ulLB0 = 3;  /* Bootstrap primes will be regenerated---avoids complications */
//	nPrimes0 = 1;  /* Even prime 2 is an exception */
//	while (1)
//	{
//		if (ulLB0 + (uint64_t)ulSpan > (uint64_t)UINT32_MAX - 4)
//		{                      /* overflow trap */
//			ulUB0 = UINT32_MAX - 4;
//			if (ulUB0 == ulLB0)
//			{
//				nPrimes0++;
//				*nPrimes = nPrimes0;
//				*ulUB = ulUB0;
//				break;
//			}
//			if (ulUB0 < ulLB0)
//			{
//				*nPrimes = nPrimes0;
//				*ulUB = ulUB0;
//				break;
//			}
//			ulSpan = ulUB0 - ulLB0;
//		}
//		else  /* normal case */
//		{
//			ulUB0 = ulLB0 + ulSpan;
//			ulOffsetMax = ulSpan / 2;
//		}
//		ulLB0In = ulLB0;
//		ulUB0In = ulUB0;
//		vSieve(uchSieve, &ulLB0, &ulUB0, ulPrime);
//		if ((ulLB0In != ulLB0) || (ulUB0In != ulUB0))
//		{
//			fprintf(stderr, "\n ERROR, Bounds anomaly in call to vSieve:");
//			fprintf(stderr, "\n LB==> %lu to %lu    UB==> %lu to %lu\n",
//				ulLB0In, ulLB0, ulUB0In, ulUB0);
//			exit(EXIT_FAILURE);
//		}
//		ulOffset = 0;
//		ulN = ulLB0;
//		while (1)
//		{
//			if (uchSieve[ulOffset])ulPrime[++nPrimes0] = ulN;
//			ulN += 2;  /* ulN is odd */
//			if ((ulN > *ulUB) && (nPrimes0 >= *nPrimes))
//			{
//				*nPrimes = nPrimes0;
//				*ulUB = ulPrime[nPrimes0];
//				goto DONE_vGP;
//			}
//			ulOffset++;
//			if (ulOffset > ulOffsetMax)break;
//		}
//		ulLB0 = ulN;
//	}
//
//DONE_vGP:
//
//	free(uchSieve);
//	return;
//}
/**********************************************************************/
//void vSieve(unsigned char *uchSieve, unsigned long *ulLB,
//	unsigned long *ulUB, unsigned long *ulPrime)
//{
//	/* A deprecated routine which now reverts to vSieveForDivisors. The
//	   reference array ulPrime is no longer needed, and may be replaced
//	   by a dummy variable in the call. The original routine is still
//	   coded below (conditionally excluded by an #if 0) for reference,
//	   but appears to offer little if any advantage. */
//
//	vSieveForDivisors(uchSieve, *ulLB, *ulUB, ceill(sqrtl(*ulUB + 0.5L)));
//	return;
//
//#if 0
//
//	/* Applies the Sieve of Eratosthenes to the interval between the odd
//	   positive integers *ulLB and *ulUB (*ulLB >= 3; *ulUB <= 2^32-1;
//	   *ulLB <= *ulUB). Elements of the sieving array uchSieve represent
//	   successive odd integers from *ulLB to *ulUB inclusive; uchSieve is
//	   assumed pre-allocated with at least 1 + (*ulUB - *ulLB) elements
//	   (unsigned chars). Sieving is carried out by elements of ulPrime up
//	   to ceil(sqrt(*ulUB)); ulPrime is assumed previously populated with
//	   primes up to the next prime beyond ceil(sqrt(*ulUB)), with
//	   ulPrime[0]=0 (or 1) and ulPrime[i]=the ith prime.
//
//	   Upon return, the prime elements in uchSieve will be represented
//	   by 1's and the composite elements by 0's.
//	*/
//
//	unsigned long ul, ulSqrt, ulUBLocal, ulDiv, ulRem, ulOffset, ulOffsetMax;
//
//	/* Must have *ulLB <= 3 <= *ulUB <= 2^32-1 with *ulLB and *ulUB odd. */
//
//	if (*ulUB < *ulLB)return;
//	if (!(*ulLB & 1))(*ulLB)++;
//	if (!(*ulUB & 1))(*ulUB)--;
//	if (*ulUB < *ulLB)return;
//
//	ulOffsetMax = (*ulUB - *ulLB) / 2;
//	for (ul = 0; ul <= ulOffsetMax; ul++)uchSieve[ul] = 1;
//	if (*ulLB == 1)uchSieve[0] = 0;
//	ulSqrt = floorl(sqrtl(*ulUB + 0.5L));
//	ul = 2;     /* evens already sieved out; start with divisor 3 */
//	while (1)  /* Sieve by ulPrime[ul++] in each loop */
//	{
//		ulDiv = ulPrime[ul++];
//		if (ulDiv == 0)
//		{
//			fprintf(stderr, "\n ERROR: ulDiv=0 in vSieve for ul=%lu", ul - 1);
//			fprintf(stderr, "\n        *ulLB=%lu  *ulUB=%lu\n", *ulLB, *ulUB);
//			fprintf(stderr, "\n        ulPrime[%lu]=%lu\n", ul - 2, ulPrime[ul - 2]);
//			exit(EXIT_FAILURE);
//		}
//		if (ulDiv > ulSqrt)break;
//		ulRem = *ulLB % ulDiv;
//		if (ulRem == 0)
//			ulOffset = 0;
//		else if (ulRem == 2 * (ulRem / 2))           /* ulRem even */
//			ulOffset = (2 * ulDiv - ulRem) / 2;     /* Assumes ulDiv < INT32_MAX */
//		else
//			ulOffset = (ulDiv - ulRem) / 2;       /* ulRem odd */
//		if (ulDiv >= *ulLB)ulOffset += ulDiv;  /* do not mark the prime itself */
//		while (ulOffset <= ulOffsetMax)
//		{
//			uchSieve[ulOffset] = 0;
//			ulOffset += ulDiv;
//		}
//	}
//
//	return;
//
//#endif
//}
/**********************************************************************/
//void vSieveForDivisors(unsigned char *uchPrime, unsigned long ulStart,
//	unsigned long ulEnd, unsigned long ulMaxDiv)
//{
//	/* Secondary sieve, for intervals from as little as ulStart=3 to
//	 * as much as ulEnd=4294967295. This is intended primarily for
//	 * generating trial divisors > 2^16, for use in primary sieving of
//	 * intervals exceeding the unsigned long range (ullEnd > 2^32).
//	 * It can function as a "partial" sieve, with trial divisors cut
//	 * off at ulMaxDiv for speed. This would leave a few composites
//	 * marked as primes; if they are used as trial divisors in a
//	 * primary sieve, this does not affect the final result, but may
//	 * result in improved performance overall.
//	 *
//	 * The interval [ulStart, ulEnd], with odd positive integer bounds,
//	 * is sieved for primes, using the algorithm of Eratosthenes and trial
//	 * divisors up to sqrt(ulEnd); i.e., for each prime q <= sqrt(ulEnd),
//	 * all integer multiples mq, ulStart <= mq <= ulEnd and m > 1, are
//	 * marked as composite. The sieved interval is represented by the byte
//	 * array uchPrime2, each byte representing an odd number---
//	 * uchPrime2[0]=ulStart and uchPrime[ulMaxIndex2]=ulEnd. For a given
//	 * prime divisor q, the offset ulIndex of the first odd mq in the
//	 * interval is computed, then the positions of successive odd multiples
//	 * are found simply by repeated addition of q.
//	 *
//	 * ulMaxDiv specifies an optional maximum value for the trial divisors
//	 * used within the secondary sieve. If ulMaxDiv is 0 or 1, all prime trial
//	 * divisors < 2^16 are available; otherwise, the sieving process is
//	 * terminated with the largest prime not exceeding ulMaxDiv. Such
//	 * termination has the effect of leaving some composites (those with
//	 * all prime factors exceeding ulMaxDiv) marked as primes. Experience
//	 * indicates that a proper choice of ulMaxDiv will save more time in
//	 * the secondary sieve than is lost by the presence of the spurious
//	 * composites in the primary sieve (their presence has no effect on
//	 * the final result, other than the time loss).
//	 *
//	 * uchPrime is the (unsigned character) array encoding the primes;
//	 *   byte 0 represents ulStart; byte (ulEnd - ulStart)/2 represents
//	 *   ulEnd; the values are either 0 (composite) or 1 (prime or almost
//	 *   prime). The proper space must have been allocated for uchPrime
//	 *   prior to calling this routine.
//	 * ulStart and ulEnd are odd positive integers, 2 < ulStart <= ulEnd < 2^32.
//	 * ulMaxDiv is the maximum trial divisor to be used in sieving. A
//	 *   value of 0 or 1 is treated as the largest 32-bit prime, 4294967291.
//	 *   If ulMaxDiv < sqrt(ulEnd), some of the "primes" returned may be
//	 *   composites with all prime factors > ulMaxDiv ("almost" primes).
//	 * ulPrime16 is an array of unsigned longs representing the primes
//	 *   to 65537 inclusive (all 16-bit primes plus one overflow value).
//	 *   Its use is internal, and the end user need not be aware of its
//	 *   existence. It provides the trial divisors for sieving.
//	 */
//
//	unsigned long ul = 2, ulMaxIndex, r, ulIndex, ulDivisor, ulSqrt; // ulQuot;
//
//	if (ulStart > ulEnd)return;
//	if ((ulStart & 1) == 0)ulStart++;
//	if ((ulEnd & 1) == 0)ulEnd--;
//	ulMaxIndex = (ulEnd - ulStart) / 2;
//	if (ulMaxDiv < 2)ulMaxDiv = MAX_32BIT_PRIME;
//
//	vGenPrimes16();  /* Will actually be done only once */
//	ulSqrt = floor(sqrt(ulEnd + 0.5));
//
//	/* Initially, all the odd integers in [ulStart, ulEnd] are assumed
//	   prime. */
//
//	for (ulIndex = 0; ulIndex <= ulMaxIndex; ulIndex++)uchPrime[ulIndex] = 1;
//	if (ulStart == 1)uchPrime[0] = 0;
//	while (1)
//	{
//		ulDivisor = ulPrime16[ul++];
//		if (ulDivisor == MIN_32BIT_PRIME)break;  /* 65537 = trial divisor overflow */
//		if (ulDivisor > ulSqrt)return;  /* surpassed sqrt of upper bound */
//		if (ulDivisor > ulMaxDiv)break;
//
//		/* Calculate the offset in uchPrime of the first odd multiple of
//		   ulDivisor in the sieve interval. */
//
//		r = ulStart % ulDivisor;
//		ulIndex = r ? ((r & 1) ? ((ulDivisor - r) >> 1) : (((ulDivisor << 1) - r) >> 1)) : 0;
//#if 0
//		/* The preceding statement is optimized for speed; it assumes
//		 * ulDivisor < INT32_MAX. Its less cryptic equivalent follows. */
//
//		if (r == 0)
//			ulIndex = 0;
//		else if (r & 1 == 1)               /* if r is odd */
//			ulIndex = (ulDivisor - r) / 2;
//		else                              /* if r is even and not zero */
//			ulIndex = (2 * ulDivisor - r) / 2;
//#endif
//		/* Now avoid striking out the prime divisor itself. */
//
//		if (ulDivisor == ulStart + 2 * ulIndex)ulIndex += ulDivisor;
//
//		while (1)
//		{
//			if (ulIndex > ulMaxIndex)break;
//			if (uchPrime[ulIndex])uchPrime[ulIndex] = 0;
//			ulIndex += ulDivisor;
//		}
//	}
//
//	return;
//}
/**********************************************************************/
//void vSieveULL(unsigned char *uchPrime, uint64_t ullStart, uint64_t ullEnd)
//{
//	/* Primary sieve, for intervals from as little as ullStart=1 to
//	 * as much as ullEnd=2^64-1 (approximately 1.8e19).
//	 *
//	 * uchPrime is the array representing the (successive odd integers of)
//	 *   the sieved interval (integer 2 is not represented or detected).
//	 *   It must be pre-allocated by the calling routine with at least
//	 *   (ullEnd - ullStart)/2 + 1 bytes of storage.
//	 * ullStart and ullEnd are odd positive integers,
//	 *   0 <= ullStart <= ullEnd < UINT64_MAX (approximately 1.8*10^19);
//	 *   however, ullEnd - ullStart must be < 2*UINT32_MAX - 1.
//	 *
//	 * The interval [ullStart, ullEnd], with odd positive integer bounds,
//	 * is sieved for primes, using the algorithm of Eratosthenes and trial
//	 * divisors up to sqrt(ullEnd); i.e., for each prime q <= sqrt(ullEnd),
//	 * all integer multiples mq, ullStart <= mq <= ullEnd and m > 1, are
//	 * marked as composite. The sieved interval is represented by the byte
//	 * array uchPrime, each byte representing an odd number---
//	 * uchPrime[0]=ullStart and uchPrime[ulMaxIndex]=ullEnd. For a given
//	 * prime divisor q, the offset ulIndex of the first odd mq in the
//	 * interval is computed, then the positions of successive odd multiples
//	 * are found simply by repeated addition of q.
//	 *
//	 * If trial divisors exceeding 2^16 are required (ullEnd > 2^32), they
//	 * are generated using a secondary sieve, which is tuned to return a
//	 * few composite trial divisors (this may improve efficiency), and a
//	 * call to vSieveForDivisors.
//	 *
//	 */
//
//	unsigned char *uchPrime2;
//	long slSS;
//	unsigned long ul, ulMaxIndex, r, ulIndex, ulDivisor, ulMaxIndex2, ulSqrt2;
//	uint64_t ullOffset, ull1, ull2, ullQuot;
//
//	if (ullStart > ullEnd)return;
//	if ((ullStart & 1) == 0)ullStart++;
//	if ((ullEnd & 1) == 0)ullEnd--;
//	if (ullEnd - ullStart > __ULL(2)*UINT32_MAX - __ULL(2))
//	{
//		fprintf(stderr, "\n ERROR: Sieving interval too large in vSieveULL.\n");
//		exit(EXIT_FAILURE);
//	}
//	ulMaxIndex = (ullEnd - ullStart) / 2;
//
//	/* Initially, all the odd integers (except 1) in [ullStart, ullEnd]
//	   are assumed prime. */
//
//	for (ulIndex = 0; ulIndex <= ulMaxIndex; ulIndex++)uchPrime[ulIndex] = 1;
//	if (ullStart == 1)uchPrime[0] = 0;
//
//	/* Main loop, for divisors 3 through 65521. */
//
//	vGenPrimes16();
//	ul = 2;
//	ulSqrt2 = ulSqrt(ullEnd) + 1;
//	while (1)
//	{
//		ulDivisor = ulPrime16[ul++];
//		if (ulDivisor == MIN_32BIT_PRIME)break;    /* divisors exceed ulPrime16 */
//		if (ulDivisor > ulSqrt2)return;  /* reached sqrt of upper bound */
//
//		/* Calculate the offset in uchPrime of the first odd multiple of
//		   ulDivisor in the sieve interval. */
//
//		r = ullStart % ulDivisor;
//
//		if (r == 0)
//			ulIndex = 0;
//		else if ((r & 1) == 1)  /* if r is odd */
//			ulIndex = (ulDivisor - r) / 2;
//		else  /* if r is even and not zero */
//			ulIndex = (2 * ulDivisor - r) / 2;
//
//		/*
//		 * Faster but cryptic version of the above, for speed; both assume
//		 * ulDivisor < UINT32_MAX/2.
//		 *
//		 * ulIndex=r?((r&1)?((ulDivisor-r)>>1):(((ulDivisor<<1)-r)>>1)):0;
//		 *
//		 */
//
//		 /* Avoid striking out the prime divisor itself */
//
//		if (ulDivisor == ullStart + 2 * ulIndex)ulIndex += ulDivisor;
//
//		while (1)
//		{
//			if (ulIndex > ulMaxIndex)break;
//			if (uchPrime[ulIndex])uchPrime[ulIndex] = 0;
//			ulIndex += ulDivisor;
//		}
//	}
//
//	ull1 = MIN_32BIT_PRIME;  /* this is 65537 */
//
//	/* If we arrive here, the prime divisors have exceeded 2^16, due to
//	 * ullEnd exceeding 2^32. A secondary sieve is called to generate
//	 * the necessary trial divisors. For efficiency reasons, these may
//	 * include some composite integers as well. Allocate memory for
//	 * the secondary sieve.
//	 */
//
//	slSS = __ulMem() / 2;
//	slSS = 1000000UL * (slSS / 1000000UL);
//	if (slSS > 100000000L)slSS = 100000000UL;
//	while (1)
//	{
//		if (slSS < 100000L)
//		{
//			fprintf(stderr, "\n ERROR: Insufficient free memory for uchPrime2.");
//			exit(EXIT_FAILURE);
//		}
//		uchPrime2 = (unsigned char *)malloc(slSS);
//		if (uchPrime2)break;
//		slSS = slSS - 1000000UL;
//	}
//
//	while (1)
//	{
//		if (ull1 > MAX_32BIT_PRIME)break;  /* this is 4294967291 */
//		if (ull1 > ulSqrt2)break;  /* reached sqrt of upper bound */
//		ull2 = ull1 + __ULL(2)*slSS - 2;
//		if (ull2 > MAX_32BIT_PRIME)ull2 = MAX_32BIT_PRIME;
//
//		/* Initialize the secondary sieve array. Call the secondary sieve
//		 * with an upper limit of 1000 on the secondary trial divisors; this
//		 * will cause some composites (those with no prime factors < 1000)
//		 * to be returned as trial divisors by the secondary sieve, but the
//		 * hope is that the few spurious composites will be offset by the
//		 * time saved in skipping more than 97 % of the secondary trial
//		 * divisors. The presence of composite trial divisors has no effect
//		 * on the eventual value of the uchPrime array.
//		 */
//
//		ulMaxIndex2 = (ull2 - ull1) / 2;
//		for (ulIndex = 0; ulIndex <= ulMaxIndex2; ulIndex++)
//			uchPrime2[ulIndex] = 1;
//		if (ull1 == 1)uchPrime2[0] = 0;
//		vSieveForDivisors(uchPrime2, ull1, ull2, 1000);
//
//		for (ul = 0; ul <= ulMaxIndex2; ul++)
//		{
//			if (uchPrime2[ul] == 0)continue;  /* not a prime divisor */
//			ulDivisor = ull1 + 2 * ul;
//			if (ulDivisor > ulSqrt2)break;  /* reached sqrt of upper bound */
//
//			/* Calculate the offset in uchPrime of the first odd multiple of
//			   ulDivisor in the sieve interval. */
//
//			r = ullStart % ulDivisor;
//
//			if (r == 0)
//				ullOffset = 0;
//			else if ((r & 1) == 1)  /* if r is odd */
//				ullOffset = (ulDivisor - r) / 2;
//			else  /* if r is even and not zero */
//				ullOffset = (__ULL(2)*ulDivisor - r) / 2;  /* Avoids UL overflow */
//
//			  /* Avoid striking out the prime divisor itself */
//
//			if (ulDivisor == ullStart + 2 * ullOffset)ullOffset += ulDivisor;
//
//			while (1)
//			{
//				if (ullOffset > ulMaxIndex)break;
//				if (uchPrime[(long)ullOffset])uchPrime[(long)ullOffset] = 0;
//				ullOffset += ulDivisor;
//			}
//		}
//		ull1 = ull2 + 2;  /* Proceed to next block */
//	}
//
//	free(uchPrime2);
//	return;
//}
/**********************************************************************/

//int iIsPrime32(unsigned long ulN)
//{
//	/* Returns 1 if ulN is prime, zero otherwise. No sieving is used. The
//	   routine simply checks for prime divisors up to the sqrt of ulN. */
//
//	unsigned long ulSqrtN, ul = 2, ulDiv;
//
//	if ((ulN < 3) || ((ulN & 1) == 0))return(ulN == 2 ? 1 : 0);
//
//	if (!iPrime16Initialized)vGenPrimes16();
//	ulSqrtN = ulSqrt(ulN);
//
//	while (1)
//	{
//		ulDiv = ulPrime16[ul++];
//		if (ulDiv > ulSqrtN)return(1);
//		if (ulN%ulDiv == 0)return(0);
//	}
//}
/**********************************************************************/
#ifdef __GMP__
/**********************************************************************/
//int iIsPrime64(uint64_t ullN, unsigned long ulMaxDivisor)
//{
//	/* Returns 1 if ullN is prime, zero otherwise. No sieving is used.
//	   The routine checks for prime divisors up to the smaller of the
//	   sqrt of ullN or ulMaxDivisor. If no prime divisor is found, and
//	   N > ulMaxDivisor^2 exceeds 2^32, the strong BPSW primality test
//	   is invoked. If 0 or 1 is specified for ulMaxDivisor, a default
//	   value of 1000 is used. */
//
//	static long        iFirst = 1;
//	unsigned long	   ulSqrtN, ul = 2, ulDiv;
//	mpz_t              mpzN;
//
//	if ((ullN < 3) || ((ullN & 1) == 0))return(ullN == 2 ? 1 : 0);
//	if (!iPrime16Initialized)vGenPrimes16();
//
//	if (ulMaxDivisor < 2)ulMaxDivisor = 1000UL;
//	if (ulMaxDivisor > 65536UL)ulMaxDivisor = 65536UL;
//
//	ulSqrtN = ulSqrt(ullN);
//	while (1)
//	{
//		ulDiv = ulPrime16[ul++];
//		if (ulDiv > ulSqrtN)return(1);
//		if (ulDiv > ulMaxDivisor)break;
//		if (ullN%ulDiv == 0)return(0);
//	}
//
//	/* If there are no small prime divisors, we use the strong BPSW test
//	   for primality. */
//
//	if (iFirst)
//	{
//		mpz_init2(mpzN, 512);
//		iFirst = 0;
//	}
//	__mpz_set_ull(mpzN, ullN);
//	return(iBPSW2(mpzN));
//}
/**********************************************************************/
int iPrP(mpz_t mpzN, unsigned long ulNMR, unsigned long ulMaxDivisor)
{
	/* Returns 1 if mpzN is a probable prime according to the strong
	 * modified Baillie-PSW test.
	 *
	 * Returns 0 if mpzN is definitely composite (and mpzN > 1).
	 * Returns 0 if mpzN is not prime (and mpzN < 2).
	 *
	 * ulNMR indicates the total number of Miller's tests to be used;
	 * the default is one (with B=2). If ulNMR > 1, ulNMR - 1 extra
	 * Miller's tests will be carried out (ulNMR <= 6543), along with
	 * an additional floor(ulNMR/5) extra strong Lucas tests.
	 *
	 * ulMaxDivisor specifies the upper bound for small prime trial divisors
	 * If ulMaxDivisor < (# of binary digits in N) is specified, the
	 * small prime divisors checked in Miller's test are used (up to
	 * qMax = # of binary digits in N).
	 *
	 * This test consists of the strong Baillie-PSW test enhanced as
	 * follows: (1) The domain of trial divisors may be altered;
	 * (2) The number of Miller's tests may be increased by specifying
	 * ulNMR > 1; (3) if the total number of tests ulNMR > 4, an
	 * additional extra strong Lucas test is performed after each five
	 * Miller's tests.
	 *
	 */

	int iComp2;
	unsigned long qMax;
	unsigned long ulDiv, ul;

	/* First eliminate all N < 3 and all even N. */

	iComp2 = mpz_cmp_si(mpzN, 2);
	if (iComp2 < 0)return(0);
	if (iComp2 == 0)return(1);
	if (mpz_even_p(mpzN))return(0);

	/* Any small prime divisors will be found in iMiller. */

	qMax = mpz_sizeinbase(mpzN, 2);
	if (ulMaxDivisor > qMax)
	{
		ulDiv = ulPrmDiv(mpzN, ulMaxDivisor);
		if (ulDiv == 1)return(1);
		if (ulDiv > 1)return(0);
	}

	if (iMiller(mpzN, 2) == 0)return(0);  /* Carry out Miller's test with base 2 */

	/* Now N is a prime, or a base-2 strong pseudoprime with no small
	   prime divisors. Apply the strong Lucas-Selfridge primality test. */

	if (iStrongLucasSelfridge(mpzN) == 0)return(0);

	/* The following is in addition to the strong Baillie-PSW test.
	   Additional Miller's tests (numbering ulNMR - 1) can be
	   mandated, a strategy rumored to be in use by Mathematica
	   In addition, after each five Miller's tests, we perform
	   an extra strong Lucas test. */

	if (ulNMR < 2)return(1);
	if (ulNMR > 6543)ulNMR = 6543;
	if (!iPrime16Initialized)vGenPrimes16();
	for (ul = 2; ul <= ulNMR; ul++)
	{
		if (iMiller(mpzN, ulPrime16[ul]) == 0)return(0);
		if (ul % 5 == 0)
			if (iExtraStrongLucas(mpzN, ulPrime16[ul / 5 + 1]) == 0)return(0);
	}

	return(1);
}
/**********************************************************************/
unsigned long ulPrmDiv(mpz_t mpzN, unsigned long ulMaxDivisor)
{
	/* Returns the smallest proper prime divisor (p <= ulMaxDivisor) of N.

	   If N < 2, return 0.
	   If N is prime and "small", return 1. "Small" means N < approximately
		 ulMaxDivisor^2.
	   If N is prime and "large", return 0. "Large" means N > approximately
		 ulMaxDivisor^2.
	   If N is composite and its smallest prime divisor p <= ulMaxDivisor,
		 return p.
	   If N is composite but its smallest prime divisor p > ulMaxDivisor,
		 return 0. In this case N will be "large", as above.

	   A return of 0 indicates "no conclusion"; N might be < 2
	   (questionable input), or N might be "large" (either prime
	   or composite) and have no prime divisor p <= ulMaxDivisor.

	   A return of 1 indicates that N is a "small" prime.

	   A return > 1 indicates that N is composite, and the
	   returned value is the smallest proper prime divisor of N.

	   If ulMaxDivisor is zero or one, a default value of 1000 is used.
	*/

	int iComp2, i;
	unsigned long ul, ulDiv, ulBase;
	static mpz_t mpzSqrt;
	static int d[8] = { 1,7,11,13,17,19,23,29 }, iFirst = 1;

	/* First eliminate all N < 3 and all even N. */

	iComp2 = mpz_cmp_si(mpzN, 2);
	if (iComp2 < 0)return(0);  /* No conclusion (N < 2) */
	if (iComp2 == 0)return(1);   /* Prime (N=2) */
	if (mpz_even_p(mpzN))return(2);  /* Composite (even) */

	if (!iPrime16Initialized)vGenPrimes16();
	if (ulMaxDivisor < 2)ulMaxDivisor = 1000UL;
	if (ulMaxDivisor > MAX_32BIT_PRIME)ulMaxDivisor = MAX_32BIT_PRIME;

	if (iFirst)
	{
		mpz_init(mpzSqrt);
		iFirst = 0;
	}
	mpz_sqrt(mpzSqrt, mpzN);

	ul = 2;  /* first trial divisor will be 3, the 2nd prime */
	while (1)
	{
		ulDiv = ulPrime16[ul++];
		if (ulDiv > ulMaxDivisor)return(0);  /* No conclusion */
		if (ulDiv > 65536UL)break;
		if (mpz_cmp_ui(mpzSqrt, ulDiv) < 0)return(1);  /* Prime */
		if (mpz_divisible_ui_p(mpzN, ulDiv))return(ulDiv);  /* Composite */
	}

	/* Once the 16-bit divisors have been exhausted, use trial divisors of
	   the form 30n + d, where d=1,7,11,13,17,19,23,29. */

	ulBase = 30 * (ulDiv / 30);
	while (1)
	{
		for (i = 0; i < 8; i++)
		{
			ulDiv = ulBase + d[i];
			if (ulDiv > ulMaxDivisor)return(0);  /* No conclusion */
			if (mpz_cmp_ui(mpzSqrt, ulDiv) < 0)return(1);  /* Prime */
			if (mpz_divisible_ui_p(mpzN, ulDiv))return(ulDiv);  /* Composite */
		}
		ulBase += 30;
	}

	return(0);  /* No conclusion */
}
/**********************************************************************/
int iMillerRabin(mpz_t mpzN, const long iB)
{
	return(iMiller(mpzN, iB));
}
/**********************************************************************/
int iMiller(mpz_t mpzN, long iB)
{
	/* Test N for primality using the Miller's strong probable prime
	   test with base B. See Gary Miller's famous paper ("Riemann's
	   hypothesis and tests for primality," Journal of Computation and
	   System Science, 1976, Volume 13, pp 300-317).

	   Returns 1 if N is a prime or a base-B strong probable prime.
	   Returns 0 if N is definitely not a prime (composite or < 2).

	   NOTE 1: Some will not regard this as a "pure" Miller's test with
	   base B, since certain adjustments are made, prior to applying the
	   algorithm, in order to work around invalid input values and
	   improve performance:

	   1) N < 3 and even N are screened out first.
	   2) Multiples of the small primes (to qMax=# of binary digits in N)
		  are returned as composite.
	   3) Values of B < 2 are replaced by B=2.
	   4) If N divides B exactly, B is replaced by B+1.

	   If such adjustments are not made, a third return value (e.g., -1)
	   must be allowed, indicating invalid input and an indeterminate result,
	   and complicating the calling source code.

	   NOTE 2: Not all authorities agree exactly on the conditions for
	   Miller's test. Some allow B < 2 (Rosen, "Elementary number theory and
	   its applications," 3rd ed., 1993, pp. 195-200), although there are good
	   reasons to exclude such values. On the other hand, some require
	   1 < B < N (Ribenboim, "The new book of prime number records,"
	   3rd ed., 1996, pp. 113-116, 143-146). As far as I know, no one
	   allows N to divide B, as this produces "false negatives"; e.g.,
	   N=3 and B=6 fails Miller's test, thus indicating N=3 as composite.
	   In practice, there appears to be no good reason to use any B < 2,
	   and in fact its value is almost always chosen to be a small
	   (positive) prime number. Furthermore, in my opinion, refusing to
	   first screen out even values of N and N < 3 gratuitously adds
	   unnecessary complication to the test.
	*/

	static long iFirst = 1;
	static mpz_t mpzB, mpzNm1, mpzd, mpzRem, mpzSqrt, mpzOne;
	long iComp2, iBits, s, j, q; // digits;
	unsigned long qMax;

	/* Allocate the static variables used in Miller's test. */

	if (iFirst)
	{
		mpz_init(mpzB);  /* Never needs more than one limb */
		iBits = mp_bits_per_limb * (1 + mpz_size(mpzN));
		if (iBits < 512)iBits = 512;
		mpz_init2(mpzNm1, iBits);
		mpz_init2(mpzOne, iBits);
		mpz_set_si(mpzOne, 1);
		mpz_init2(mpzd, iBits);
		mpz_init2(mpzRem, 2 * iBits);  /* must contain products */
		if (!iPrime16Initialized)vGenPrimes16();
		iFirst = 0;
	}

	/* First take care of all N < 3 and all even N. */

	iComp2 = mpz_cmp_si(mpzN, 2);
	if (iComp2 < 0)return(0);        /* N < 2 is by convention not prime */
	if (iComp2 == 0)return(1);         /* N=2 is prime */
	if (mpz_even_p(mpzN))return(0);  /* Even N > 2 is composite */

	/* Try small prime divisors from 3 to an UB qMax determined by the size
	   of N (qMax >= 31). */

	mpz_sqrt(mpzSqrt, mpzN);
	qMax = mpz_sizeinbase(mpzN, 2);  /* Number of binary digits in N */
	if (qMax < 36)qMax = 36;
	j = 2;  /* First trial divisor is 3, the second prime */
	while (1)
	{
		q = ulPrime16[j++];
		if (q > qMax)break;
		if (mpz_cmp_si(mpzN, q) == 0)return(1);
		if (mpz_cmp_si(mpzSqrt, q) < 0)return(1);
		if (mpz_divisible_ui_p(mpzN, q))return(0);
	}

	/* Check for valid input. Miller's test requires B > 1, and N must not
	   divide B exactly. Choose B=2 and B<--B+1 if these problems arise.
	   This is technically a deviation from the pure Miller's test, but
	   avoids the necessity of handling an error return of -1. */

	if (iB < 2)iB = 2;
	mpz_set_si(mpzB, iB);
	if (mpz_divisible_p(mpzB, mpzN))mpz_add_ui(mpzB, mpzB, 1);

	/* Now compute d and s, where d is odd and N - 1 = (2^s)*d. */

	mpz_sub_ui(mpzNm1, mpzN, 1);
	s = mpz_scan1(mpzNm1, 0);
	mpz_tdiv_q_2exp(mpzd, mpzNm1, s);

	/* Now proceed with the Miller's algorithm. First, if B^d is
	   congruent to 1 mod N, N is a strong probable prime to base B. */

	mpz_powm(mpzRem, mpzB, mpzd, mpzN);
	if (mpz_cmp_si(mpzRem, 1) == 0)return(1);

	/* Now calculate B^((2^j)*d), for j=0,1,...,s-1 by successive
	   squaring. If any of these is congruent to -1 mod N, N is a
	   sprp base B. Start with j=0 and B^d, already computed.
	   Miller's test uses repeated modular squaring in place of repeated
	   modular exponentiation for speed (squaring is an order of
	   magnitude faster). */

	if (mpz_cmp(mpzRem, mpzNm1) == 0)return(1);  /* j=0 case */
	for (j = 1; j < s; j++)
	{
		mpz_mul(mpzRem, mpzRem, mpzRem);
		mpz_mod(mpzRem, mpzRem, mpzN);
		if (mpz_cmp(mpzRem, mpzNm1) == 0)return(1);
		if (mpz_cmp(mpzRem, mpzOne) == 0)return(0);
	}

	return(0);
}
/**********************************************************************/

int iLucasSelfridge(mpz_t mpzN)
{
	/* Test mpzN for primality using Lucas's test with Selfridge's parameters.
	   Returns 1 if mpzN is prime or a Lucas-Selfridge pseudoprime. Returns
	   0 if mpzN is definitely composite. Note that a Lucas-Selfridge test
	   typically requires three to seven times as many bit operations as a
	   single Miller's test. The frequency of Lucas-Selfridge pseudoprimes
	   appears to be roughly four times that of base-2 strong pseudoprimes;
	   the Baillie-PSW test is based on the hope (verified by the author,
	   May, 2005, for all N < 10^13; and by Martin Fuller, January, 2007,
	   for all N < 10^15) that the two tests have no common pseudoprimes. */

	int iComp2, iP, iJ, iSign;
	long lDabs, lD, lQ;
	unsigned long ulMaxBits, ulNbits, ul, ulGCD;
	mpz_t mpzU, mpzV, mpzNplus1, mpzU2m, mpzV2m, mpzQm, mpz2Qm,
		mpzT1, mpzT2, mpzT3, mpzT4, mpzD;

#undef RETURN
#define RETURN(n)           \
  {                         \
  mpz_clear(mpzU);          \
  mpz_clear(mpzV);          \
  mpz_clear(mpzNplus1);     \
  mpz_clear(mpzU2m);        \
  mpz_clear(mpzV2m);        \
  mpz_clear(mpzQm);         \
  mpz_clear(mpz2Qm);        \
  mpz_clear(mpzT1);         \
  mpz_clear(mpzT2);         \
  mpz_clear(mpzT3);         \
  mpz_clear(mpzT4);         \
  mpz_clear(mpzD);          \
  return(n);                \
  }

	/* This implementation of the algorithm assumes N is an odd integer > 2,
	   so we first eliminate all N < 3 and all even N. As a practical matter,
	   we also need to filter out all perfect square values of N, such as
	   1093^2 (a base-2 strong pseudoprime); this is because we will later
	   require an integer D for which Jacobi(D,N) = -1, and no such integer
	   exists if N is a perfect square. The algorithm as written would
	   still eventually return zero in this case, but would require
	   nearly sqrt(N)/2 iterations. */

	iComp2 = mpz_cmp_si(mpzN, 2);
	if (iComp2 < 0)return(0);
	if (iComp2 == 0)return(1);
	if (mpz_even_p(mpzN))return(0);
	if (mpz_perfect_square_p(mpzN))return(0);

	/* Allocate storage for the mpz_t variables. Most require twice
	   the storage of mpzN, since multiplications of order O(mpzN)*O(mpzN)
	   will be performed. */

	ulMaxBits = 2 * mpz_sizeinbase(mpzN, 2) + mp_bits_per_limb;
	mpz_init2(mpzU, ulMaxBits);
	mpz_init2(mpzV, ulMaxBits);
	mpz_init2(mpzNplus1, ulMaxBits);
	mpz_init2(mpzU2m, ulMaxBits);
	mpz_init2(mpzV2m, ulMaxBits);
	mpz_init2(mpzQm, ulMaxBits);
	mpz_init2(mpz2Qm, ulMaxBits);
	mpz_init2(mpzT1, ulMaxBits);
	mpz_init2(mpzT2, ulMaxBits);
	mpz_init2(mpzT3, ulMaxBits);
	mpz_init2(mpzT4, ulMaxBits);
	mpz_init(mpzD);

	/* Find the first element D in the sequence {5, -7, 9, -11, 13, ...}
	   such that Jacobi(D,N) = -1 (Selfridge's algorithm). Although
	   D will nearly always be "small" (perfect square N's having
	   been eliminated), an overflow trap for D is present. */

	lDabs = 5;
	iSign = 1;
	while (1)
	{
		lD = iSign * lDabs;
		iSign = -iSign;
		ulGCD = mpz_gcd_ui(NULL, mpzN, lDabs);
		/* if 1 < GCD < N then N is composite with factor lDabs, and
		   Jacobi(D,N) is technically undefined (but often returned
		   as zero). */
		if ((ulGCD > 1) && mpz_cmp_ui(mpzN, ulGCD) > 0)RETURN(0);
		mpz_set_si(mpzD, lD);
		iJ = mpz_jacobi(mpzD, mpzN);
		if (iJ == -1)break;
		lDabs += 2;
		if (lDabs > ulDmax)ulDmax = lDabs;  /* tracks global max of |D| */
		if (lDabs > INT32_MAX - 2)
		{
			fprintf(stderr,
				"\n ERROR: D overflows signed long in Lucas-Selfridge test.");
			fprintf(stderr, "\n N=");
			mpz_out_str(stderr, 10, mpzN);
			fprintf(stderr, "\n |D|=%ld\n\n", lDabs);
			exit(EXIT_FAILURE);
		}
	}

	iP = 1;         /* Selfridge's choice */
	lQ = (1 - lD) / 4;  /* Required so D = P*P - 4*Q */

	/* NOTE: The conditions (a) N does not divide Q, and
	   (b) D is square-free or not a perfect square, are included by
	   some authors; e.g., "Prime numbers and computer methods for
	   factorization," Hans Riesel (2nd ed., 1994, Birkhauser, Boston),
	   p. 130. For this particular application of Lucas sequences,
	   these conditions were found to be immaterial. */

	mpz_add_ui(mpzNplus1, mpzN, 1); /* must compute U_(N - Jacobi(D,N)) */

	/* mpzNplus1 is always even, so the accumulated values U and V
	   are initialized to U_0 and V_0 (if the target index were odd,
	   U and V would be initialized to U_1=1 and V_1=P). In either case,
	   the values of U_2m and V_2m are initialized to U_1 and V_1;
	   the FOR loop calculates in succession U_2 and V_2, U_4 and
	   V_4, U_8 and V_8, etc. If the corresponding bits of N+1 are
	   on, these values are then combined with the previous totals
	   for U and V, using the composition formulas for addition
	   of indices. */

	mpz_set_ui(mpzU, 0);           /* U=U_0 */
	mpz_set_ui(mpzV, 2);           /* V=V_0 */
	mpz_set_ui(mpzU2m, 1);         /* U_1 */
	mpz_set_si(mpzV2m, iP);        /* V_1 */
	mpz_set_si(mpzQm, lQ);
	mpz_set_si(mpz2Qm, 2 * lQ);

	ulNbits = mpz_sizeinbase(mpzNplus1, 2);
	for (ul = 1; ul < ulNbits; ul++)  /* zero bit off, already accounted for */
	{
		/* Formulas for doubling of indices (carried out mod N). Note that
		 * the indices denoted as "2m" are actually powers of 2, specifically
		 * 2^(ul-1) beginning each loop and 2^ul ending each loop.
		 *
		 * U_2m = U_m*V_m
		 * V_2m = V_m*V_m - 2*Q^m
		 */
		mpz_mul(mpzU2m, mpzU2m, mpzV2m);
		mpz_mod(mpzU2m, mpzU2m, mpzN);
		mpz_mul(mpzV2m, mpzV2m, mpzV2m);
		mpz_sub(mpzV2m, mpzV2m, mpz2Qm);
		mpz_mod(mpzV2m, mpzV2m, mpzN);
		if (mpz_tstbit(mpzNplus1, ul))
		{
			/* Formulas for addition of indices (carried out mod N);
			 *
			 * U_(m+n) = (U_m*V_n + U_n*V_m)/2
			 * V_(m+n) = (V_m*V_n + D*U_m*U_n)/2
			 *
			 * Be careful with division by 2 (mod N)!
			 */
			mpz_mul(mpzT1, mpzU2m, mpzV);
			mpz_mul(mpzT2, mpzU, mpzV2m);
			mpz_mul(mpzT3, mpzV2m, mpzV);
			mpz_mul(mpzT4, mpzU2m, mpzU);
			mpz_mul_si(mpzT4, mpzT4, lD);
			mpz_add(mpzU, mpzT1, mpzT2);
			if (mpz_odd_p(mpzU))mpz_add(mpzU, mpzU, mpzN);
			mpz_fdiv_q_2exp(mpzU, mpzU, 1);
			mpz_add(mpzV, mpzT3, mpzT4);
			if (mpz_odd_p(mpzV))mpz_add(mpzV, mpzV, mpzN);
			mpz_fdiv_q_2exp(mpzV, mpzV, 1);
			mpz_mod(mpzU, mpzU, mpzN);
			mpz_mod(mpzV, mpzV, mpzN);
		}
		/* Calculate Q^m for next bit position, doubling the exponent.
		   The irrelevant final iteration is omitted. */
		if (ul < ulNbits - 1)  /* Q^m not needed for MSB. */
		{

			mpz_mul(mpzQm, mpzQm, mpzQm);
			mpz_mod(mpzQm, mpzQm, mpzN);  /* prevents overflow */
			mpz_add(mpz2Qm, mpzQm, mpzQm);
		}
	}

	/* If U_(N - Jacobi(D,N)) is congruent to 0 mod N, then N is
	   a prime or a Lucas pseudoprime; otherwise it is definitely
	   composite. */

	if (mpz_sgn(mpzU) == 0)RETURN(1);
	RETURN(0);
}
/**********************************************************************/
int iStrongLucasSelfridge(mpz_t mpzN)
{
	/* Test N for primality using the strong Lucas test with Selfridge's
	   parameters. Returns 1 if N is prime or a strong Lucas-Selfridge
	   pseudoprime (in which case N is also a pseudoprime to the standard
	   Lucas-Selfridge test). Returns 0 if N is definitely composite.

	   The running time of the strong Lucas-Selfridge test is, on average,
	   roughly 10 % greater than the running time for the standard
	   Lucas-Selfridge test (3 to 7 times that of a single Miller's test).
	   However, the frequency of strong Lucas pseudoprimes appears to be
	   only (roughly) 30 % that of (standard) Lucas pseudoprimes, and only
	   slightly greater than the frequency of base-2 strong pseudoprimes,
	   indicating that the strong Lucas-Selfridge test is more computationally
	   effective than the standard version. */

	int iComp2, iP, iJ, iSign;
	long lDabs, lD, lQ;
	unsigned long ulMaxBits, uldbits, ul, ulGCD, r, s;
	mpz_t mpzU, mpzV, mpzNplus1, mpzU2m, mpzV2m, mpzQm, mpz2Qm,
		mpzT1, mpzT2, mpzT3, mpzT4, mpzD, mpzd, mpzQkd, mpz2Qkd;

#undef RETURN
#define RETURN(n)           \
  {                         \
  mpz_clear(mpzU);          \
  mpz_clear(mpzV);          \
  mpz_clear(mpzNplus1);     \
  mpz_clear(mpzU2m);        \
  mpz_clear(mpzV2m);        \
  mpz_clear(mpzQm);         \
  mpz_clear(mpz2Qm);        \
  mpz_clear(mpzT1);         \
  mpz_clear(mpzT2);         \
  mpz_clear(mpzT3);         \
  mpz_clear(mpzT4);         \
  mpz_clear(mpzD);          \
  mpz_clear(mpzd);          \
  mpz_clear(mpzQkd);        \
  mpz_clear(mpz2Qkd);       \
  return(n);                \
  }

	/* This implementation of the algorithm assumes N is an odd integer > 2,
	   so we first eliminate all N < 3 and all even N. As a practical matter,
	   we also need to filter out all perfect square values of N, such as
	   1093^2 (a base-2 strong pseudoprime); this is because we will later
	   require an integer D for which Jacobi(D,N) = -1, and no such integer
	   exists if N is a perfect square. The algorithm as written would
	   still eventually return zero in this case, but would require
	   nearly sqrt(N)/2 iterations. */

	iComp2 = mpz_cmp_si(mpzN, 2);
	if (iComp2 < 0)return(0);
	if (iComp2 == 0)return(1);
	if (mpz_even_p(mpzN))return(0);
	if (mpz_perfect_square_p(mpzN))return(0);

	/* Allocate storage for the mpz_t variables. Most require twice
	   the storage of mpzN, since multiplications of order O(mpzN)*O(mpzN)
	   will be performed. */

	ulMaxBits = 2 * mpz_sizeinbase(mpzN, 2) + mp_bits_per_limb;
	mpz_init2(mpzU, ulMaxBits);
	mpz_init2(mpzV, ulMaxBits);
	mpz_init2(mpzNplus1, ulMaxBits);
	mpz_init2(mpzU2m, ulMaxBits);
	mpz_init2(mpzV2m, ulMaxBits);
	mpz_init2(mpzQm, ulMaxBits);
	mpz_init2(mpz2Qm, ulMaxBits);
	mpz_init2(mpzT1, ulMaxBits);
	mpz_init2(mpzT2, ulMaxBits);
	mpz_init2(mpzT3, ulMaxBits);
	mpz_init2(mpzT4, ulMaxBits);
	mpz_init(mpzD);
	mpz_init2(mpzd, ulMaxBits);
	mpz_init2(mpzQkd, ulMaxBits);
	mpz_init2(mpz2Qkd, ulMaxBits);

	/* Find the first element D in the sequence {5, -7, 9, -11, 13, ...}
	   such that Jacobi(D,N) = -1 (Selfridge's algorithm). Theory
	   indicates that, if N is not a perfect square, D will "nearly
	   always" be "small." Just in case, an overflow trap for D is
	   included. */

	lDabs = 5;
	iSign = 1;
	while (1)
	{
		lD = iSign * lDabs;
		iSign = -iSign;
		ulGCD = mpz_gcd_ui(NULL, mpzN, lDabs);
		/* if 1 < GCD < N then N is composite with factor lDabs, and
		   Jacobi(D,N) is technically undefined (but often returned
		   as zero). */
		if ((ulGCD > 1) && mpz_cmp_ui(mpzN, ulGCD) > 0)RETURN(0);
		mpz_set_si(mpzD, lD);
		iJ = mpz_jacobi(mpzD, mpzN);
		if (iJ == -1)break;
		lDabs += 2;
		if (lDabs > ulDmax)ulDmax = lDabs;  /* tracks global max of |D| */
		if (lDabs > INT32_MAX - 2)
		{
			fprintf(stderr,
				"\n ERROR: D overflows signed long in Lucas-Selfridge test.");
			fprintf(stderr, "\n N=");
			mpz_out_str(stderr, 10, mpzN);
			fprintf(stderr, "\n |D|=%ld\n\n", lDabs);
			exit(EXIT_FAILURE);
		}
	}

	iP = 1;         /* Selfridge's choice */
	lQ = (1 - lD) / 4;  /* Required so D = P*P - 4*Q */

	/* NOTE: The conditions (a) N does not divide Q, and
	   (b) D is square-free or not a perfect square, are included by
	   some authors; e.g., "Prime numbers and computer methods for
	   factorization," Hans Riesel (2nd ed., 1994, Birkhauser, Boston),
	   p. 130. For this particular application of Lucas sequences,
	   these conditions were found to be immaterial. */

	   /* Now calculate N - Jacobi(D,N) = N + 1 (even), and calculate the
		  odd positive integer d and positive integer s for which
		  N + 1 = 2^s*d (similar to the step for N - 1 in Miller's test).
		  The strong Lucas-Selfridge test then returns N as a strong
		  Lucas probable prime (slprp) if any of the following
		  conditions is met: U_d=0, V_d=0, V_2d=0, V_4d=0, V_8d=0,
		  V_16d=0, ..., etc., ending with V_{2^(s-1)*d}=V_{(N+1)/2}=0
		  (all equalities mod N). Thus d is the highest index of U that
		  must be computed (since V_2m is independent of U), compared
		  to U_{N+1} for the standard Lucas-Selfridge test; and no
		  index of V beyond (N+1)/2 is required, just as in the
		  standard Lucas-Selfridge test. However, the quantity Q^d must
		  be computed for use (if necessary) in the latter stages of
		  the test. The result is that the strong Lucas-Selfridge test
		  has a running time only slightly greater (order of 10 %) than
		  that of the standard Lucas-Selfridge test, while producing
		  only (roughly) 30 % as many pseudoprimes (and every strong
		  Lucas pseudoprime is also a standard Lucas pseudoprime). Thus
		  the evidence indicates that the strong Lucas-Selfridge test is
		  more effective than the standard Lucas-Selfridge test, and a
		  Baillie-PSW test based on the strong Lucas-Selfridge test
		  should be more reliable. */


	mpz_add_ui(mpzNplus1, mpzN, 1);
	s = mpz_scan1(mpzNplus1, 0);
	mpz_tdiv_q_2exp(mpzd, mpzNplus1, s);

	/* We must now compute U_d and V_d. Since d is odd, the accumulated
	   values U and V are initialized to U_1 and V_1 (if the target
	   index were even, U and V would be initialized instead to U_0=0
	   and V_0=2). The values of U_2m and V_2m are also initialized to
	   U_1 and V_1; the FOR loop calculates in succession U_2 and V_2,
	   U_4 and V_4, U_8 and V_8, etc. If the corresponding bits
	   (1, 2, 3, ...) of t are on (the zero bit having been accounted
	   for in the initialization of U and V), these values are then
	   combined with the previous totals for U and V, using the
	   composition formulas for addition of indices. */

	mpz_set_ui(mpzU, 1);                      /* U=U_1 */
	mpz_set_ui(mpzV, iP);                     /* V=V_1 */
	mpz_set_ui(mpzU2m, 1);                    /* U_1 */
	mpz_set_si(mpzV2m, iP);                   /* V_1 */
	mpz_set_si(mpzQm, lQ);
	mpz_set_si(mpz2Qm, 2 * lQ);
	mpz_set_si(mpzQkd, lQ);  /* Initializes calculation of Q^d */

	uldbits = mpz_sizeinbase(mpzd, 2);
	for (ul = 1; ul < uldbits; ul++)  /* zero bit on, already accounted for */
	{
		/* Formulas for doubling of indices (carried out mod N). Note that
		 * the indices denoted as "2m" are actually powers of 2, specifically
		 * 2^(ul-1) beginning each loop and 2^ul ending each loop.
		 *
		 * U_2m = U_m*V_m
		 * V_2m = V_m*V_m - 2*Q^m
		 */
		mpz_mul(mpzU2m, mpzU2m, mpzV2m);
		mpz_mod(mpzU2m, mpzU2m, mpzN);
		mpz_mul(mpzV2m, mpzV2m, mpzV2m);
		mpz_sub(mpzV2m, mpzV2m, mpz2Qm);
		mpz_mod(mpzV2m, mpzV2m, mpzN);
		/* Must calculate powers of Q for use in V_2m, also for Q^d later */
		mpz_mul(mpzQm, mpzQm, mpzQm);
		mpz_mod(mpzQm, mpzQm, mpzN);  /* prevents overflow */
		mpz_mul_2exp(mpz2Qm, mpzQm, 1);
		if (mpz_tstbit(mpzd, ul))
		{
			/* Formulas for addition of indices (carried out mod N);
			 *
			 * U_(m+n) = (U_m*V_n + U_n*V_m)/2
			 * V_(m+n) = (V_m*V_n + D*U_m*U_n)/2
			 *
			 * Be careful with division by 2 (mod N)!
			 */
			mpz_mul(mpzT1, mpzU2m, mpzV);
			mpz_mul(mpzT2, mpzU, mpzV2m);
			mpz_mul(mpzT3, mpzV2m, mpzV);
			mpz_mul(mpzT4, mpzU2m, mpzU);
			mpz_mul_si(mpzT4, mpzT4, lD);
			mpz_add(mpzU, mpzT1, mpzT2);
			if (mpz_odd_p(mpzU))mpz_add(mpzU, mpzU, mpzN);
			mpz_fdiv_q_2exp(mpzU, mpzU, 1);
			mpz_add(mpzV, mpzT3, mpzT4);
			if (mpz_odd_p(mpzV))mpz_add(mpzV, mpzV, mpzN);
			mpz_fdiv_q_2exp(mpzV, mpzV, 1);
			mpz_mod(mpzU, mpzU, mpzN);
			mpz_mod(mpzV, mpzV, mpzN);
			mpz_mul(mpzQkd, mpzQkd, mpzQm);  /* Calculating Q^d for later use */
			mpz_mod(mpzQkd, mpzQkd, mpzN);
		}
	}

	/* If U_d or V_d is congruent to 0 mod N, then N is a prime or a
	   strong Lucas pseudoprime. */

	if (mpz_sgn(mpzU) == 0)RETURN(1);
	if (mpz_sgn(mpzV) == 0)RETURN(1);

	/* NOTE: Ribenboim ("The new book of prime number records," 3rd ed.,
	   1995/6) omits the condition Vр0 on p.142, but includes it on
	   p. 130. The condition is NECESSARY; otherwise the test will
	   return false negatives---e.g., the primes 29 and 2000029 will be
	   returned as composite. */

	   /* Otherwise, we must compute V_2d, V_4d, V_8d, ..., V_{2^(s-1)*d}
		  by repeated use of the formula V_2m = V_m*V_m - 2*Q^m. If any of
		  these are congruent to 0 mod N, then N is a prime or a strong
		  Lucas pseudoprime. */

	mpz_mul_2exp(mpz2Qkd, mpzQkd, 1);  /* Initialize 2*Q^(d*2^r) for V_2m */
	for (r = 1; r < s; r++)
	{
		mpz_mul(mpzV, mpzV, mpzV);
		mpz_sub(mpzV, mpzV, mpz2Qkd);
		mpz_mod(mpzV, mpzV, mpzN);
		if (mpz_sgn(mpzV) == 0)RETURN(1);
		/* Calculate Q^{d*2^r} for next r (final iteration irrelevant). */
		if (r < s - 1)
		{
			mpz_mul(mpzQkd, mpzQkd, mpzQkd);
			mpz_mod(mpzQkd, mpzQkd, mpzN);
			mpz_mul_2exp(mpz2Qkd, mpzQkd, 1);
		}
	}

	/* Otherwise N is definitely composite. */

	RETURN(0);
}
/**********************************************************************/
int iExtraStrongLucas(mpz_t mpzN, long lB)
{
	/* Test N for primality using the extra strong Lucas test with base B,
	   as formulated by Zhaiyu Mo and James P. Jones ("A new primality test
	   using Lucas sequences," preprint, circa 1997), and described by Jon
	   Grantham in "Frobenius pseudoprimes," (preprint, 16 July 1998),
	   available at <http://www.pseudoprime.com/pseudo1.ps>.

	   Returns 1 if N is prime or an extra strong Lucas pseudoprime (base B).
	   Returns 0 if N is definitely composite.

	   Even N and N < 3 are eliminated before applying the Lucas test.

	   In this implementation of the algorithm, Q=1, and B is an integer
	   in 2 < B < INT32_MAX (2147483647 on 32-bit machines); the default value
	   is B=3. B is incremented as necessary if the values of B and N are
	   inconsistent with the hypotheses of Jones and Mo: P=B, Q=1,
	   D=P*P - 4*Q, GCD(N,2D)=1, Jacobi(D,N) <> 0.

	   Since the base B is used solely to calculate the discriminant
	   D=B*B - 4, negative values of B are redundant. The bases B=0 and
	   B=1 are excluded because they produce huge numbers of pseudoprimes,
	   and B=2 is excluded because the resulting D=0 fails the Jones-Mo
	   hypotheses.

	   Note that the choice Q=1 eliminates the computation of powers of Q
	   which appears in the weak and strong Lucas tests.

	   The running time of the extra strong Lucas-Selfridge test is, on
	   average, roughly 80 % that of the standard Lucas-Selfridge test
	   or 2 to 6 times that of a single Miller's test. This is superior
	   in speed to both the standard and strong Lucas-Selfridge tests. The
	   frequency of extra strong Lucas pseudoprimes also appears to be
	   about 80 % that of the strong Lucas-Selfridge test and 30 % that of
	   the standard Lucas-Selfridge test, comparable to the frequency of
	   spsp(2).

	   Unfortunately, the apparent superior peformance of the extra strong
	   Lucas test is offset by the fact that it is not "backwards compatible"
	   with the Lucas-Selfridge tests, due to the differing choice of
	   parameters: P=B and Q=1 in the extra strong test, while P=1 and
	   Q=(1 - D)/4 in the standard and strong Lucas-Selfridge tests (with D
	   chosen from the sequence 5, -7, 9, ...). Thus, although every extra
	   strong Lucas pseudoprime to base B is also both a strong and standard
	   Lucas pseudoprime with parameters P=B and Q=1, the extra strong
	   pseudoprimes do *NOT* constitute a proper subset of the Lucas-Selfridge
	   standard and strong pseudoprimes. As a specific example, 4181 is an
	   extra strong Lucas pseudoprime to base 3, but is neither a standard
	   nor strong Lucas-Selfridge pseudoprime.

	   As a result, the corresponding Baillie-PSW test is fatally flawed.
	   Regardless of the base chosen for the extra strong Lucas test, it
	   appears that there exist numerous N for which the corresponding
	   extra strong Lucas pseudoprimes (xslpsp) will also be strong
	   pseudoprimes to base 2 (or any other particular Miller's base).
	   For example, 6368689 is both spsp(2) and xslpsp(3); 8725753
	   is both spsp(2) and xslpsp(11); 80579735209 is spsp(2) and
	   simultaneously xslpsp for the bases 3, 5, and 7; 105919633 is
	   spsp(3) and xslpsp(11); 1121176981 is spsp(19) and xslpsp(31);
	   and so on. Perhaps some combination of the extra strong test
	   and multiple Miller's tests could match the performance of the
	   Lucas-Selfridge BPSW tests, but the prospects do not look bright.
	*/

	int iComp2, iJ;
	long lD, lP, lQ;
	unsigned long ulMaxBits, uldbits, ul, ulGCD, r, s;
	mpz_t mpzU, mpzV, mpzM, mpzU2m, mpzV2m, mpzT1, mpzT2, mpzT3, mpzT4,
		mpzD, mpzd, mpzTwo, mpzMinusTwo;

#undef RETURN
#define RETURN(n)           \
  {                         \
  mpz_clear(mpzU);          \
  mpz_clear(mpzV);          \
  mpz_clear(mpzM);          \
  mpz_clear(mpzU2m);        \
  mpz_clear(mpzV2m);        \
  mpz_clear(mpzT1);         \
  mpz_clear(mpzT2);         \
  mpz_clear(mpzT3);         \
  mpz_clear(mpzT4);         \
  mpz_clear(mpzD);          \
  mpz_clear(mpzd);          \
  mpz_clear(mpzTwo);        \
  mpz_clear(mpzMinusTwo);   \
  return(n);                \
  }

	/* This implementation of the algorithm assumes N is an odd integer > 2,
	   so we first eliminate all N < 3 and all even N. */

	iComp2 = mpz_cmp_si(mpzN, 2);
	if (iComp2 < 0)return(0);
	if (iComp2 == 0)return(1);
	if (mpz_even_p(mpzN))return(0);

	/* Allocate storage for the mpz_t variables. Most require twice
	   the storage of mpzN, since multiplications of order O(mpzN)*O(mpzN)
	   will be performed. */

	ulMaxBits = 2 * mpz_sizeinbase(mpzN, 2) + mp_bits_per_limb;
	mpz_init2(mpzU, ulMaxBits);
	mpz_init2(mpzV, ulMaxBits);
	mpz_init2(mpzM, ulMaxBits);
	mpz_init2(mpzU2m, ulMaxBits);
	mpz_init2(mpzV2m, ulMaxBits);
	mpz_init2(mpzT1, ulMaxBits);
	mpz_init2(mpzT2, ulMaxBits);
	mpz_init2(mpzT3, ulMaxBits);
	mpz_init2(mpzT4, ulMaxBits);
	mpz_init(mpzD);
	mpz_init2(mpzd, ulMaxBits);
	mpz_init_set_si(mpzTwo, 2);
	mpz_init_set_si(mpzMinusTwo, -2);

	/* The parameters specified by Zhaiyu Mo and James P. Jones,
	   as set forth in Grantham's paper, are P=B, Q=1, D=P*P - 4*Q,
	   with (N,2D)=1 so that Jacobi(D,N) <> 0. As explained above,
	   bases B < 3 are excluded. */

	if (lB < 3)
		lP = 3;
	else
		lP = lB;
	lQ = 1;

	/* We check to make sure that N and D are relatively prime. If not,
	   then either 1 < (D,N) < N, in which case N is composite with
	   divisor (D,N); or N = (D,N), in which case N divides D and may be
	   either prime or composite, so we increment the base B=P and
	   try again. */

	while (1)
	{
		lD = lP * lP - 4 * lQ;
		ulGCD = mpz_gcd_ui(NULL, mpzN, labs(lD));
		if (ulGCD == 1)break;
		if (mpz_cmp_ui(mpzN, ulGCD) > 0)RETURN(0);
		lP++;
	}

	/* Now calculate M = N - Jacobi(D,N) (M even), and calculate the
	   odd positive integer d and positive integer s for which
	   M = 2^s*d (similar to the step for N - 1 in Miller's
	   test). The extra strong Lucas-Selfridge test then returns N as
	   an extra strong Lucas probable prime (eslprp) if any of the
	   following conditions is met: U_d=0 and V_dрс2; or V_d=0; or
	   V_2d=0, V_4d=0, V_8d=0, V_16d=0, ..., etc., ending with
	   V_{2^(s-2)*d}=V_{M/4}р0 (all equalities mod N). Thus d is the
	   highest index of U that must be computed (since V_2m is
	   independent of U), compared to U_M for the standard Lucas
	   test; and no index of V beyond M/4 is required, compared to
	   M/2 for the standard and strong Lucas tests. Furthermore,
	   since Q=1, the powers of Q required in the standard and
	   strong Lucas tests can be dispensed with. The result is that
	   the extra strong Lucas test has a running time shorter than
	   that of either the standard or strong Lucas-Selfridge tests
	   (roughly two to six times that of a single Miller's test).
	   The extra strong test also produces fewer pseudoprimes.
	   Unfortunately, the pseudoprimes produced are *NOT* a subset
	   of the standard or strong Lucas-Selfridge pseudoprimes (due
	   to the incompatible parameters P and Q), and consequently the
	   extra strong test does not combine with a single Miller's test
	   to produce a Baillie-PSW test of the reliability level of the
	   BPSW tests based on the standard or strong Lucas-Selfridge tests. */

	mpz_set_si(mpzD, lD);
	iJ = mpz_jacobi(mpzD, mpzN);
	assert(iJ != 0);
	if (iJ == 1)
		mpz_sub_ui(mpzM, mpzN, 1);
	else
		mpz_add_ui(mpzM, mpzN, 1);

	s = mpz_scan1(mpzM, 0);
	mpz_tdiv_q_2exp(mpzd, mpzM, s);

	/* We must now compute U_d and V_d. Since d is odd, the accumulated
	   values U and V are initialized to U_1 and V_1 (if the target
	   index were even, U and V would be initialized instead to U_0=0
	   and V_0=2). The values of U_2m and V_2m are also initialized to
	   U_1 and V_1; the FOR loop calculates in succession U_2 and V_2,
	   U_4 and V_4, U_8 and V_8, etc. If the corresponding bits
	   (1, 2, 3, ...) of t are on (the zero bit having been accounted
	   for in the initialization of U and V), these values are then
	   combined with the previous totals for U and V, using the
	   composition formulas for addition of indices. */

	mpz_set_ui(mpzU, 1);                       /* U=U_1 */
	mpz_set_si(mpzV, lP);                      /* V=V_1 */
	mpz_set_ui(mpzU2m, 1);                     /* U_1 */
	mpz_set_si(mpzV2m, lP);                    /* V_1 */

	uldbits = mpz_sizeinbase(mpzd, 2);
	for (ul = 1; ul < uldbits; ul++)  /* zero bit on, already accounted for */
	{
		/* Formulas for doubling of indices (carried out mod N). Note that
		 * the indices denoted as "2m" are actually powers of 2, specifically
		 * 2^(ul-1) beginning each loop and 2^ul ending each loop.
		 *
		 * U_2m = U_m*V_m
		 * V_2m = V_m*V_m - 2*Q^m
		 */
		mpz_mul(mpzU2m, mpzU2m, mpzV2m);
		mpz_mod(mpzU2m, mpzU2m, mpzN);
		mpz_mul(mpzV2m, mpzV2m, mpzV2m);
		mpz_sub_ui(mpzV2m, mpzV2m, 2);
		mpz_mod(mpzV2m, mpzV2m, mpzN);
		if (mpz_tstbit(mpzd, ul))
		{
			/* Formulas for addition of indices (carried out mod N);
			 *
			 * U_(m+n) = (U_m*V_n + U_n*V_m)/2
			 * V_(m+n) = (V_m*V_n + D*U_m*U_n)/2
			 *
			 * Be careful with division by 2 (mod N)!
			 */
			mpz_mul(mpzT1, mpzU2m, mpzV);
			mpz_mul(mpzT2, mpzU, mpzV2m);
			mpz_mul(mpzT3, mpzV2m, mpzV);
			mpz_mul(mpzT4, mpzU2m, mpzU);
			mpz_mul_si(mpzT4, mpzT4, lD);
			mpz_add(mpzU, mpzT1, mpzT2);
			if (mpz_odd_p(mpzU))mpz_add(mpzU, mpzU, mpzN);
			mpz_fdiv_q_2exp(mpzU, mpzU, 1);
			mpz_add(mpzV, mpzT3, mpzT4);
			if (mpz_odd_p(mpzV))mpz_add(mpzV, mpzV, mpzN);
			mpz_fdiv_q_2exp(mpzV, mpzV, 1);
			mpz_mod(mpzU, mpzU, mpzN);
			mpz_mod(mpzV, mpzV, mpzN);
		}
	}

	/* N first passes the extra strong Lucas test if V_dр0, or if V_dрс2
	   and U_dр0.  U and V are tested for divisibility by N, rather than
	   zero, in case the previous FOR is a zero-iteration loop.*/

	if (mpz_divisible_p(mpzV, mpzN))RETURN(1);
	if (mpz_divisible_p(mpzU, mpzN))
	{
		if (mpz_congruent_p(mpzV, mpzTwo, mpzN))RETURN(1);
		if (mpz_congruent_p(mpzV, mpzMinusTwo, mpzN))RETURN(1);
	}

	/* Otherwise, we must compute V_2d, V_4d, V_8d, ..., V_{2^(s-2)*d}
	   by repeated use of the formula V_2m = V_m*V_m - 2*Q^m. If any of
	   these are congruent to 0 mod N, then N is a prime or an extra
	   strong Lucas pseudoprime. */

	for (r = 1; r < s - 1; r++)
	{
		mpz_mul(mpzV, mpzV, mpzV);
		mpz_sub_ui(mpzV, mpzV, 2);
		mpz_mod(mpzV, mpzV, mpzN);
		if (mpz_sgn(mpzV) == 0)RETURN(1);
	}

	/* Otherwise N is definitely composite. */

	RETURN(0);
}
/**********************************************************************/
#endif  /* GMP available */
/**********************************************************************/

/**********************************************************************/
/*                       string editing                               */
/**********************************************************************/
/**********************************************************************/
char *szTrimMWS(char *pch)
{
	return(szTrimLWS(szTrimTWS(pch)));
}
/**********************************************************************/
char *szTrimTWS(char *pch)
{
	char            *pchStart;
	//long            sl;
	unsigned long   ulLen;

	if (*pch == 0)return(pch);
	pchStart = pch;
	ulLen = strlen(pch);
	pch = pchStart + ulLen - 1;
	while (1)
	{
		if (isgraph(*pch))
			return(pchStart);
		else
			*pch = 0;
		if (pch == pchStart)
		{
			*pch = 0;
			return(pch);
		}
		pch--;
	}
}
/**********************************************************************/
char *szTrimLWS(char *pch)
{
	char            *pchStart;

	pchStart = pch;
	while (1)
	{
		if (*pch == 0)
		{
			*pchStart = 0;
			return(pch);
		}
		else if (isgraph(*pch))
		{
			if (pch == pchStart)
				return(pch);
			else
			{
				memmove(pchStart, pch, strlen(pch) + 1);
				return(pch);
			}
		}
		else
			pch++;
	}
}
/********************************************************/

/* Returns the number of seconds elapsed since some fixed event, dependent 
  upon the function call and platform. */
double lfSeconds2(void)
{
	/* NOTE: The peculiar design of this routine is due to MinGW's
	   inablity to correctly parse the preprocessor directive
	   #if (CLOCKS_PER_SEC <= 1000) and similar directives.

	   Returns the number of seconds elapsed since some fixed event,
	   dependent upon the function call and platform. The clock()
	   based routine normally returns the number of seconds since
	   either the beginning of program execution or the first call
	   to the function. The gettimeofday and time(NULL) based routines
	   return the number of seconds since the beginning of the UNIX epoch
	   (00:00:00 GMT 1 Jan 1970). The Win32 GetTickCount returns the
	   number of milliseconds since system boot (in reality, resolution
	   may be as poor as +/- 55 ms); QueryPerformanceCounter can
	   theoretically achieve nanosecond resolution. The granularity of the
	   clock() routine is generally either 0.01 sec or 0.055 sec. The
	   granularity of gettimeofday is nominally 1 microsecond, but
	   in reality 0.01 second is more common. The granularity of
	   time(NULL) is 1 second.

	   PORTABILITY AND BUGS: The clock() and time(NULL) functions are
	   part of standard C. The gettimeofday function is not part of
	   standard C, but is available on most platforms. GetTickCount()
	   and QueryPerformanceCounter(*LARGE_INTEGER) are part of Win32.
	   GetTickCount has maximum compatibility, supported at least
	   since Windows 3.0 (1991). QueryPerformanceCounter is a more
	   recent addition to the Windows API, available since at least
	   December, 2002; it is designed to take advantage of the
	   Time Stamp Counter (TSC) on Pentium compatible processors,
	   and has a potential resolution level in the nanosecond range.

	   Known bugs in clock() include the rollover problem, which
	   will usually cause INT32_MAX to rollover to INT32_MIN after
	   2^31 ticks. This is a major problem on systems (including many
	   GNU/Linux systems) which comply with the P*SIX standard
	   CLOCKS_PER_SECOND=1000000; then the first rollover occurs after
	   less than 36 minutes. Rollover results in clock() failing to be
	   monotonic increasing, so that simply differencing clock() values
	   may not reflect the true time difference (and may even generate
	   a ridiculous negative time difference). Rollover can generally
	   be ignored on systems where CLOCKS_PER_SECOND <= 1000, as
	   rollover will take at least 24.85 days. Otherwise, it must be
	   trapped in the routine, and this can become quite problematical
	   because of the possibility of multiple rollovers and masked
	   rollovers. Cygwin's clock() is avoided, as it yields erroneous
	   values (typically twice the correct value) on some systems.

	   Bugs in gettimeofday have been reported by several users; these
	   are either "backward jumps" in value in rare instances, or
	   anomalous values returned at local midnight and then quickly
	   self-correcting. More recent versions of gettimeofday (starting
	   with GNU/Linux 2.4) appear to be more reliable, but I have
	   observed the midnight anomaly on my own Windows systems, using
	   the gettimeofday in DJGPP 2.03. It appears to have no rollover
	   problem, although one may be coming in 2038, when the UNIX
	   epoch attains 2^31 seconds.

	   The only bugs of which I am aware in time(NULL) are a midnight
	   rollover anomaly, similar to the one exhibited by gettimeofday,
	   on some Windows systems; and the same Y2K type problem
	   looming in 2038. Its huge disadvantage is the poor granularity.

	   As a last resort, time(NULL) is returned, with a granularity
	   of only one second.

	   The clock() routine has a correction factor to compensate for
	   DJGPP's use of 91/5 PC clock ticks per second (the correct
	   value is 1193180/65536 = 18.2046819336). */

	static double lftPrevious = -1e308, lft, dt;

	/* lftPrevious stores the value from the previous call. If the
	   next computed value is slightly *less* than lftPrevious, it
	   is almost certainly the result of an interfering OS layer
	   blocking direct access to the hardware clock, and lftPrevious
	   is returned instead, in an effort to keep the function
	   monotonic non-decreasing (and time intervals non-negative).
	   If lftPrevious is significantly less than the computed lft,
	   it is probably due to a system clock rollover, for which
	   there is no simple remedy; in that case, the inconsistent
	   present value is returned (with apologies). NOTE: -DBL_MAX
	   is avoided as the value of the initializer above due to a
	   bug in the DJGPP compilers. */

#if defined(__WIN32__)  /* MinGW, Cygwin, DM, Borland */
#if defined(WINVER) && (WINVER >= 0x0400)
	   /* Thanks to Huang Yuanbing <bailuzhou(at)163(dot)com> for this
		  code fragment, adapted from his PrimeNumber.cpp code
		  (12 Feb 2009). The preceding conditional is a fudged test
		  for the availability of the QueryPerformanceXXX routines, to
		  be replaced by a more discriminating one if and when found. */
	static LARGE_INTEGER s_freq;
	LARGE_INTEGER performanceCount;
	if (s_freq.QuadPart == 0 && !QueryPerformanceFrequency(&s_freq))
	{
		/* Performance counter has failed; fall back on GetTickCount */
		lft = GetTickCount() / 1000.0;
		dt = lft - lftPrevious;
		if ((dt >= 0) || (dt < -1))
		{
			lftPrevious = lft;
			return(lft);
		}
		else
			return(lftPrevious);  /* guard against non-monotonic clock */
	}
	QueryPerformanceCounter(&performanceCount);
	return(performanceCount.QuadPart / (double)s_freq.QuadPart);
#else
	lft = GetTickCount() / 1000.0;
	dt = lft - lftPrevious;
	if ((dt >= 0) || (dt < -1))
	{
		lftPrevious = lft;
		return(lft);
	}
	else
		return(lftPrevious);  /* guard against non-monotonic clock */
#endif
#elif defined(CLOCKS_PER_SEC)
	static unsigned long ulCPS = CLOCKS_PER_SEC;
	if (ulCPS <= 1000)
	{
		lft = clock() / ((double)CLOCKS_PER_SEC);
#ifdef __DJGPP__
		if (ulCPS == 91)lft = 0.9996439766*lft;
#endif
		dt = lft - lftPrevious;
		if ((dt >= 0) || (dt < -1))
		{
			lftPrevious = lft;
			return(lft);
		}
		else
			return(lftPrevious);
	}
	else
	{
#if defined(__LINUX__)
		static struct timeval tv;
		gettimeofday(&tv, NULL);
		lft = tv.tv_sec + tv.tv_usec / 1000000.0;
		dt = lft - lftPrevious;
		if ((dt >= 0) || (dt < -1))
		{
			lftPrevious = lft;
			return(lft);
		}
		else
			return(lftPrevious);
#else
		return(time(NULL));  /* last resort */
#endif
	}
#endif
}
/**********************************************************************/
void vAtExit(void)
{
	__OBSL__;
	return;
}
/**********************************************************************/


