INTRODUCTION

This is a calculator and factorisation program based on Dario Alpert's program.
See https://www.alpertron.com.ar/ECM.HTM

It factorises numbers or numeric expressions using fast algorithms ECM and SIQS.
ECM = Lenstra elliptic-curve factorization
see https://en.wikipedia.org/wiki/Lenstra_elliptic-curve_factorization
SIQS = Self-initializing quadratic sieve
see http://www.mersennewiki.org/index.php/Self-Initializing_Quadratic_Sieve
or https://en.wikipedia.org/wiki/Quadratic_sieve

HISTORY

the program was converted from browser-based to a console window. The advantage
is that it is somewhat faster, and the calculator part has been carefully tested.
The disadvantage is that you have to download and build it. In order to do that
you need GMP or MPIR multi-precision library and the Boost multi-precision library.

the following changes were made:

C++ rather than C

The calculator part was extensively rewritten to use GMP/MPIR extended precision
and the Boost multiprecision library. Wherever possible GMP functions were used
instead of DAs version.

Hexadecimal output is turned on with the X command and turned off with D.
If numbers are displayed in hexadecimal they are preceded by 0x. If the number is
negative it is displayed in 2s complement form, and the 1st digit will be in the 
range 8 to f. If necessary an extra f is added at the beginning e.g. -129 is 
displayed as 0xf7f. If the number is positive the 1st digit will be in the range 
0 to 7. If necessary an extra 0 is added at the beginning e.g 128 is displayed as 
0x080

GMP functions used include:
mpz_probable_prime_p        test if prime
mpz_fac_ui                  factorial
mpz_2fac_ui                 double factorial
mpz_fib_ui                  fibonacci
mpz_gcd                     Greatest Common Denominator
mpz_invert                  modular inversion
mpz_lucnum_ui               Lucas number
mpz_nextprime               get next prime
mpz_powm                    modular exponent
mpz_primorial_ui            primorial
mpz_pow_ui                  exponent
mpz_mul_2exp                left shift
mpz_fdiv_q_2exp             right shift
mpz_xor                     bitwise exclusive or
mpz_ior                     bitwise inclusive or
mpz_and                     bitwise and

Standard arithmetic operators for extended precision values use the Boost library,  
which implements them using GMP/MPIR functions.

The PARTITION function was rewritten as DAs version initially had a bug and an 
alternative that worked was already available.

FACTORIZATION

As mentioned above it uses algorithms ECM and SIQS.

The factoriser is essentially DAs program, with an interface function that converts 
GMP/MPIR extended precision numbers to DAs BigIntegers and vice versa. The progress 
messages it produces have been modified to work with a console window instead of a 
Web Browser. 

A couple of checks for array index overflows were added. If an error is detected 
an exception is thrown, an error message is output, and the program continues 
without factorising the number. I suspect there are other places checks should be 
added, but realistically any number small enough to be factorised in a reasonable 
time will not exceed the array bounds.

Profiling using Visual Studio suggests that changing from DAs BigIntegers to 
GMP/MPIR would not improve performance much. This would be a huge task.

the time taken to factorise a large number is unpredictable e.g.

1000 000000 000000 000000 008970 000000 000000 000000 014661 000000 000000 000000 006097 (76 digits)
= 1 000000 000000 000000 000007 * 10 000000 000000 000000 000013 * 100 000000 000000 000000 000067
took 291 seconds

40526 919504 877216 755680 601905 432322 134980 384796 226602 145184 481280 000000 000000 (77 digits)
 = 2^53 * 3^27 * 5^13 * 7^9 * 11^5 * 13^4 * 17^3 * 19^3 * 23^2 * 29 * 31 * 37 * 41 * 43 * 47 * 53
 took 0.054 seconds. The larger number was facorised in about 1/5000 of the time!

 The time required depends mainly on the size of the 2nd largest factor, but you don't
 know what that is in advance.

 ALTERNATIVES

 Use DA's web page.
 Advantage:     ready to use, no installation, possibility to split calculation over
                several processors, possible to use any previously known factors.
 Disadvantage:  some bugs in calculator. No way to interface to other programs

 Use Python interpreter
 Advantages:    very easy to install. 
                more user-friendly,
                line-by line interpreter available, so can be used as a (programmable)
                calculator 
                very easy to write small programs
                supports big integers as standard.
disadvantages:  no functions such as gcd, fibonacci, factorial, primorial etc
                no factorisation

Operators and Functions
calculator      C/C++   Python      notes
 +                +      +
 -                -      -
 *                *      *
 /                /      //         calculator and C/C++ uses truncation division
                                    python uses floor division (/ operator in python generates
                                    a floating point number)
 %                %      %          python uses floor division
^ or **          N/A     **         has right-to-left associativity. This is mathematically
                                    correct but seems strange to programmers. e.g.
                                    2**3**4 = 2**(3**4) = 2417851639229258349412352 = 2^81
                                    (2**3)**4 = 2**(3*4) = 4096 = 2^12
AND               &      &          bitwise and
OR                |      |          bitwise or
XOR               ^      ^          bitwise exclusive or. Do not confuse XOR with exponentiation.
NOT               ~      ~          bitwise not
SHL or <<         <<     <<         bitwise left shift
SHR or >>         >>     >>         bitwise arithmetic right shift
 C                N/A    N/A        binomial coeffiecient. nCk = n!/(k!*(n-k)!) but is more efficient

comparision operators <, <=, ==, !=, >, and >= are the same in all three. The calculator
returns -1 for true, 0 for false. This allows this allows AND, OR, XOR and NOT to operate
on the returned values as if they were boolean variables.

#(primorial) !(factorial) and !!(double factorial) operators have the highest priority, 
above ^ (exponentiation)
The normal rules for operator precedence and use of brackets to over-ride the default order
of evaluation apply.

the following functions and operators are only in the calculator
n!                                  factorial
n!!                                 double factorial, not to be confused with (n!)!
n#                                  primorial. (DA's calculator requires n to be prime)

B(n)                                Previous probable prime before n
F(n)                                Fibonacci number Fn
L(n)                                Lucas number Ln = Fn-1 + Fn+1
N(n)                                Next probable prime after n
PI(n)                               the number of prime numbers less than or equal to n
P(n)                                Unrestricted Partition Number (number of 
                                    decompositions of n into sums of integers without 
                                    regard to order).
Gcd(m,n)                            Greatest common divisor of m and n.
Modinv(m,n)                         inverse of m modulo n, only valid when gcd(m,n)=1.
Modpow(m,n,r)                       finds m^n modulo r. more efficient than (m^n)%r
Totient(n)                          finds the number of positive integers less than n 
                                    which are relatively prime to n.
IsPrime(n)                          returns zero if n is not probable prime, -1 if it is.
NumDivs(n)                          Number of positive divisors of n either prime or composite.
SumDivs(n)                          Sum of all positive divisors of n both prime and composite.
NumDigits(n,r)                      Number of digits of n in base r.
SumDigits(n,r)                      Sum of digits of n in base r.
RevDigits(n,r)                      finds the value obtained by writing backwards 
                                    the digits of n in base r. 

Some functions and operators limit the range of their parameters:
^ or ** (exponent)                  0 <= exponent <= 2^31-1, also result is estimated 
                                    before calculation and if it appears to be > 20,000
                                    digits an error will be reported.
*                                   result is estimated before calculation and if it 
                                    appears to be > 20,000 digits an error will be reported.
/ (division)                        divisor must not be zero
% (modulus)                         modulus must not be zero
! (factorial)                       0 < number <= 5984 (limits result to 20,000 digits)	
!! (double factorial)               0 < number <= 11081	(limits result to 20,000 digits)
# (primorial)                       0 < number <= 46340 (limits result to 20,000 digits)
nCk (binomial coefficient)          (-2^31 <= k <= 2^31-1)	
<< and >> (shift operators)         -2^63 <= bitshift value <= 2^63-1

totient(n)                          n >= 1 (mathematically undefined for n < 1)
PI(n)                               n <= 10^9 (this is because this calculation can be 
                                    very slow)
F(n)                                n <= 95700 (limits result to 20,000 digits)
L(n)                                n <= 95700 (limits result to 20,000 digits)
P(n)                                n < 60000
