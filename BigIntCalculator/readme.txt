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
instead of DAs version of bigIntegers.

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

As mentioned above ECM and SIQS are used for factorisation.

At all stages a factor list is maintained that contains all known factors, whether
prime or not, that have not yet been split into smaller factors. Initially the
factor list just contains the number to be factorised. Whenever a new factor is found,
existing factors that are multiples of that factor are split into two factors. If
this process creates equal factors they are merged. 

Step 1: if number to factor is > than the the square of the largest prime to be 
used in trial division, check if the number is a perfect power +/- 1 and if so 
attempt to factorise it. This step may return factors that are not prime factors.

Step 2: Try to factorise each factor (there may already be more than 1) by trial division 
using a list of primes. If a factor is found that is: 
    less than the cube of the largest prime in the list, 
    and has no factors in the list of primes, 
    and the factor is not prime 
then the factor can only have two prime factors. In this case factorise it using 
Pollard's Rho algorithm. 

Step 3: For each factor that is not already known to be prime:
A. If the factor is a perfect power replace it with the number which is the root 
   and adjust the the factor's exponent accordingly. Although such numbers are  
   rare this is necessary because SIQS and ECM do not work for numbers which are
   perfect powers.
B. Test whether the number is prime: 
   If it is a Carmichael number factorise it using a specific algorithm. 
   If it is prime mark it as such. 
   Otherwise either:
      factorise it using built-in ECM, SIQS and Lehman algorithms. SIQS is only
      used for numbers between 30 and 95 digits. This method is alway used for 
      numbers less than 2^192 (58 digits)
  OR
      factorise it using Msieve. Using modified Msieve that makes more use of
   elliptic curve factorisation than the standard version (even with the e option)
   does, this is nearly always the fastest for larger numbers that would use SIQS 
   for factorisation.
   OR
      factorise it using YAFU


The built-in factoriser is essentially DAs program, with an interface function 
that converts GMP/MPIR extended precision numbers to DAs BigIntegers and vice 
versa. The progress messages it produces have been modified to work with a 
console window instead of a Web Browser. 

A couple of checks for overflows were added. If an error is detected an exception 
is thrown, an error message is output, and the program continues without 
factorising the number. Realistically any number small enough to be factorised in
a reasonable time will not cause an overflow.

Profiling using Visual Studio suggests that changing from DAs BigIntegers to 
GMP/MPIR would not improve performance. This would be a huge task.

the time taken to factorise a large number is unpredictable e.g.

1000 000000 000000 000000 008970 000000 000000 000000 014661 000000 000000 000000 006097 (76 digits)
= 1 000000 000000 000000 000007 * 10 000000 000000 000000 000013 * 100 000000 000000 000000 000067
took 291 seconds (227 seconds using Msieve, 90 seconds using YAFU)

40526 919504 877216 755680 601905 432322 134980 384796 226602 145184 481280 000000 000000 (77 digits)
= 2^53 * 3^27 * 5^13 * 7^9 * 11^5 * 13^4 * 17^3 * 19^3 * 23^2 * 29 * 31 * 37 * 41 * 43 * 47 * 53
took 0.054 seconds. The larger number was factorised in about 1/5000 of the time!

The time required depends mainly on the size of the 2nd largest factor, but you don't
know what that is in advance.

A feature was added to allow use of Msieve as an alternative to the original ECM & SIQS
factorisation. My conclusion was that for larger numbers (> about 50 digits) Msieve
is usually significantly faster, provided that the -e option is used with Msieve.
The -e option causes greater use of ECM within Msieve. The default in mseive is 
to use ECM only for factors < 15 digits. 

   The MSIEVE command controls factorisation using Msieve:

   MSIEVE ON    Turns on factorisation using Msieve (and turns YAFU off)
   MSIEVE OFF   Turns off factorisation using Msieve (revert to built-in ECM and SIQS)
   MSIEVE PATH  Displays path used when starting Msieve
   MSIEVE PATH SET Change pathe used when starting Msieve
   MSIEVE LOG   Displays path & file name for Msieve output file
   MSIEVE E ON  Turns on -e option in Msieve; perform 'deep' ECM
   MSIEVE E OFF Turns off -e option in Msieve
   MSIEVE N ON  Turns on -n option in Msieve; use NFS instead of SIQS
   MSIEVE N OFF Turns off -n option in Msieve

The source code for Msieve was downloaded from https://sourceforge.net/projects/msieve/files/
(use green buton)

For a precompiled version try https://download.mersenne.ca/msieve

To build Msieve from source requires GMP-ECM library functions, which in turn 
require Pthreads. 
download GMP-ECM from https://gforge.inria.fr/projects/ecm/. Note that only the
libecm part is requred for Msieve or YAFU.
The dll and lib files for pthreads were downloaded from SourceForge.

In addition I modified the function choose_max_digits within the source file gmp_ecm.c 
to increase the use of ECM still further. The code change is:

	if (obj->flags & MSIEVE_FLAG_DEEP_ECM) {
		if (bits > 200) {           // 200 bits = about 60 digits
			if (bits < 240)         // 240 bits = about 72 digits
				max_digits = 20;    // increased from 15
			else if (bits < 280)    // 280 bits = about 84 digits 
				max_digits = 25;    // increased from 20
			else if (bits < 320)    // 320 bits = about 96 digits
				max_digits = 30;    // increased from 25
			else if (bits < 360)    // 360 bits = about 108 digits
				max_digits = 30;
			else if (bits < 400)    // 400 bits = about 120 digits
				max_digits = 35;
			else
				max_digits = 40;
		}

Msieve is off by default but can be turned on by the "MSIEVE ON" command (which 
also turns YAFU off) and turned back off by the MSIEVE OFF command. The path to 
access msieve should be changed to whatever is appropriate.

The command MSIEVE PATH will display the current path and verify that the msieve.exe 
file exists.
The command MSIEVE PATH SET will display a windows explorer style window where 
you can navigate to the correct folder then click on the Msieve exe file. The new path 
to the Msieve file will then be saved in the BigIntCalculator.ini file.


Additionally, the option to use YAFU (Yet Another Factorisation Utility) instead of 
Msieve or the built-in ECM and SIQS has been added. The conclusion is that for 
numbers over about 60 digits YAFU is generally significantly faster than Msieve 
but for larger numbers > about 95 digits it relies on ggnfs.

The source for YAFU was obtained from https://sourceforge.net/projects/yafu/files/1.24/
(use green button to download version 1.34)
For a precompiled version try https://download.mersenne.ca/YAFU

Also need  GGNFS got precompiled from 
    https://mersenneforum.org/attachment.php?attachmentid=18244&d=1525946072
    can also try https://download.mersenne.ca/GGNFS

To build YAFU from source requires Msieve and GMP-ECM library functions, which in
turn require Pthreads. 

   The YAFU command controls factorisation using YAFU:

   YAFU ON   Turns on factorisation using YAFU (and turns Msieve off)
   YAFU OFF  Turns off factorisation using YAFU (revert to built-in ECM and SIQS)
   YAFU PATH Displays or changes path used when starting YAFU
   YAFU LOG  Displays path & file name for YAFU output file
   YAFU PLAN <name> 
             where <name> is NONE, NOECM, LIGHT, NORMAL, or DEEP.
             These correspond to the -plan options documented in YAFU's
             docfile.txt file. The default is NORMAL. This is a convenient way
             to control the amount of testing using ECM before switching to SIQS 
             or NFS.


YAFU is on by default but can be turned off by the "YAFU OFF" command and turned 
back on by the YAFU ON command (which also turns Msieve off). The default path to 
access YAFU is hard coded in yafu.cpp and should be changed to whatever is appropriate.

The command YAFU PATH will display the current path and verify that the yafu-X64.exe 
file exists.
The command YAFU PATH SET will display a windows explorer style window where 
you can navigate to the correct folder then click on the yafu exe file. The new path 
to the YAFU file will then be saved in the BigIntCalculator.ini file.

The command YAFU INI will check whether the YAFU.ini file exist. If it does it is
read and any ggnfs_dir parameter is found. If ggnfs_dir is found the path that it 
contains is checked by checking for the existence of the gnfs-lasieve4iXXe.exe files,
where XX is between 11 and 16. Note that YAFU INI requires the YAFU path to already 
be set correctly (by YAFU PATH or otherwise).

The command YAFU INI I allows the ggnfs_dir parameter to be altered or created. It
will display a windows explorer style window where you can navigate to the correct 
folder then click on any file. The new path to the ggnfs files will then be saved
in the ggnfs_dir parameter of the YAFU.ini file.


 ALTERNATIVES

 Use DA's web page.
 Advantage:     ready to use, no installation, possibility to split calculation over
                several processors, possible to use any previously known factors.
 Disadvantage:  some bugs in calculator. No way to interface to other programs

 Use Msieve directly
 Advantage:	     It's faster. For a 94 digit Mersenne number 2^311-1 it took
                 about 5 mins vs 1 hour using DA's code. 
 Disadvantage:   Build is tricky & Calculator is basic, and result isn't sent to
                 console window. default is that input & output are from/to files, 
				 not console window.

Use YAFU         Generally faster than Msieve. For a 94 digit Mersenne number 
                 2^311-1 it took about 80 sec.
                 needs GGNFS for large numbers > about 95 digits. 
                 try https://mersenneforum.org/attachment.php?attachmentid=18244&d=1525946072
                 https://www.mersenneforum.org/showthread.php?t=22215&page=5
                 https://www.mersenneforum.org/attachment.php?attachmentid=22535&d=1591570404
                 https://mersenneforum.org/showthread.php?t=25304
                 or https://download.mersenne.ca/

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
 +                +      +          unary + and - are supported and have higher priority than 
 -                -      -          other operators.
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
 C                N/A    N/A        binomial coeffiecient. nCk = n!/(k!*(n-k)!) but is 
                                    more efficient

comparision operators <, <=, ==, !=, >, and >= are similar in all three. The calculator
returns -1 for true, 0 for false. This allows this allows AND, OR, XOR and NOT to operate
on the returned values as if they were boolean variables.

#(primorial) !(factorial) and !..!(multi-factorial) operators have the highest priority, 
above NOT and unary -, which are above ^ (exponentiation)
The normal rules for operator precedence and use of brackets to over-ride the default order
of evaluation apply.

the following functions and operators are only in the calculator
n!                                  factorial
n!..!                               multi-factorial, not to be confused with (n!)!
n#                                  primorial. 

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
                                    NB. If n is -ve, the value is only defined if the 
                                    inverse of m with respect to r exitst, i.e. only 
                                    if gcd(m, r) is 1
Totient(n)                          finds the number of positive integers less than n 
                                    which are relatively prime to n.
IsPrime(n)                          returns zero if n is not probable prime, -1 if it is.
NumDivs(n)                          Number of positive divisors of n either prime or composite.
SumDivs(n)                          Sum of all positive divisors of n both prime and composite.
NumDigits(n,r)                      Number of digits of n in base r.
SumDigits(n,r)                      Sum of digits of n in base r.
RevDigits(n,r)                      finds the value obtained by writing backwards 
                                    the digits of n in base r. 
R2(n)                               The number of ways n can be formed as the sum of x^2 + y^2
                                    where x and y are negative, zero, or positive. The order is 
									significant e.g. R2(1) = 4. (0 + 1^2, 0 + (-1)^2, 1^2 +0,
									(-1)^2 + 0)
R3(n)                               The number of ways n can be formed as the sum of x^2 + y^2 + z^2
                                    where x, y and z are negative, zero, or positive.
									WARNING: for large n this is very slow
llt(n)                              Do Lucas-Lehmer primality test on 2^n-1.
                                    Return 0 if 2^n-1 is composite, 1 if prime
sqrt(n)                             Calculate floor(sqrt(n))
nroot(x, n)                         Calculate nth root of x
numfact(n)                          returns the number of uniqe factors in n
minfact(n)                          returns the value of the smallest factor of n
maxfact(n)                          returns the value of the largest factor of n

Some functions and operators limit the range of their parameters:
^ or ** (exponent)                  0 <= exponent <= 2^31-1, also result is estimated 
                                    before calculation and if it appears to be > 20,000
                                    digits an error will be reported.
*                                   result is estimated before calculation and if it 
                                    appears to be > 20,000 digits an error will be reported.
/ (division)                        divisor must not be zero
% (modulus)                         modulus must not be zero
! (factorial & multi-factorial)     0 < number <= 11081 (limits result to 20,000 digits)	

# (primorial)                       0 < number <= 46340 (limits result to 20,000 digits)
nCk (binomial coefficient)          (-2^31 <= k <= 2^31-1)	
<< and >> (shift operators)         -2^63 <= bitshift value <= 2^63-1
									result is estimated before calculation and if it 
                                    appears to be > 20,000 digits an error will be reported.
totient(n)                          n >= 1 (mathematically undefined for n < 1)
PI(n)                               n <= 10^9 (this is because this calculation can be 
                                    very slow)
F(n)                                n <= 95700 (limits result to 20,000 digits)
L(n)                                n <= 95700 (limits result to 20,000 digits)
P(n)                                n <= 1000000 (implementation limit, P(100000) takes 5
                                    minutes to calculate!)
sqrt(n)                             n >= 0
nroot(x, n)                         if n is even, x >= 0. (If n is odd -ve x is OK)

There are a number of built-in test commands. In all tests the results are
checked for correctness and a summary of the factorisation results is printed
after the tests are completed:

TEST    1st tests most of the calculator functions, then factorises a series of
        numbers. The numbers are chosen to test for various special cases as well
        as normal factorisation. This test is good for regression testing.

TEST2   test factorisation using pseudo-random numbers of a specified size.
        Command format is TEST2 [num1[,num[,num3]]] where num1 is the number of 
        tests,  num2 is the size of the numbers to be factored in bits. If num3  
        NE 0 the number to be factored consists of 2 approximately same-sized 
        'strong' prime factors, otherwise it is a random number that can contain 
        any number of factors. If the values for num2 and num3 are the same in 
        two TEST2 commands the same sequence of numbers to be factored is 
        generated. This command is useful for benchmarking factorisation.

TEST3    Tests for the built-in bigintegers. These are based on DA's biginteger 
         functions but made into a C++ class with arithmetic , shift, compare 
         operators, etc implemented. However the built-in bigintegers are now 
         only used in the the built-in ECM & SIQS, which themselves are normally 
         only used for numbers up to 2^192, so this test is not normally included.

TEST4     test factorisation of mersenne numbers. See https://en.wikipedia.org/wiki/Mersenne_prime
          this test will take over an hour.

TEST5     tests using only YAFU for factorisation. Note that these tests bypass the
          trial division etc normally used and rely on YAFU for all the factorisation.

TEST6     tests using only Msieve for factorisation. Factorise selected Mersenne numbers.

TEST7     tests the Lucas-Lehmer function. Format is TEST7 [num]. All primes <= num are
          tested to see whether 2^p-1 is prime or not. The default value for num is
          12000.   Sample output:

            test7 87000
            18:47:15 2^2 -1 is prime ***
            18:47:15 2^3 -1 is prime ***
            18:47:15 2^5 -1 is prime ***
            18:47:15 2^7 -1 is prime ***
            18:47:15 2^13 -1 is prime ***
            18:47:15 2^17 -1 is prime ***
            18:47:15 2^19 -1 is prime ***
            18:47:15 2^31 -1 is prime ***
            18:47:15 2^61 -1 is prime ***
            18:47:15 2^89 -1 is prime ***
            18:47:15 2^107 -1 is prime ***
            18:47:15 2^127 -1 is prime ***
            18:47:15 2^521 -1 is prime ***
            18:47:15 2^607 -1 is prime ***
            18:47:16 2^1279 -1 is prime ***
            18:47:16 2^2203 -1 is prime ***
            18:47:16 2^2281 -1 is prime ***
            18:47:17 2^3217 -1 is prime ***
            18:47:18 2^4253 -1 is prime ***
            18:47:18 2^4423 -1 is prime ***
            18:47:53 2^9689 -1 is prime ***
            18:47:56 2^9941 -1 is prime ***
            18:48:15 2^11213 -1 is prime ***
            18:55:42 2^19937 -1 is prime ***
            18:57:54 2^21701 -1 is prime ***
            19:00:14 2^23209 -1 is prime ***
            20:33:19 2^44497 -1 is prime ***
            08:18:34 2^86243 -1 is prime ***
            Found 28 Mersenne primes  out of 8450 numbers tested
            2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279, 2203, 2281, 
            3217, 4253, 4423, 9689, 9941, 11213, 19937, 21701, 23209, 44497, 86243,
            0.33% primes found by llt
            99.67% composites found by llt
            08:40:38 test 7 completed time used =  13h 53 min 23.314sec

          This has found the 1st 28 Mersenne primes in less than 14 hours. The 
          last of these was originally found in 1982 using a mult-million $ Cray
          supercomputer. Talk about standing on the shoulders of giants!

Test example:

enter expression to be processed, or HELP, or EXIT
TEST2 5
Use default 48 for number size in bits
 
 <******** output from indivdial tests omitted *****************>

10:48:34 All tests completed. Time used = 0.172 seconds
Test Num Size   time      Unique Factors Total Factors     2nd Fac
   1       14  0:00:00.01             4          4            4
   2       15  0:00:00.04             4          4            3
   3       15  0:00:00.03             3          3            5
   4       15  0:00:00.03             2          2            5
   5       14  0:00:00.02             5          6            4


   note: The 'Size' column gives the size in decimal digits of the number to be
   factored. The '2nd Fac' column gives the size of the 2nd largest factor, 
   which correlates much more closely with time required than the 'Size' does.
   The 'time' column gives the time required for each factorisation to one-
   hundredth of a second.

   Added 5/6/2021

   The 'engine'at the heart of the calculator was largely rewritten; 
   It was divided into 3 parts:
   1.   'Tokenise' all terms in the expression i.e. each number, operator, 
        Function name, bracket & comma is turned into a token. Also check that
        the opening and closing brackets pair up correctly.
    2.  Convert to Reverse Polish. This uses the well-known 'shunting' algorithm,
        but a recursive call to the reverse polish function is made for each
        function parameter (this also takes care of nested function calls).
        Also there is a tweak for the factorial, double factorial and primorial
        functions because the operator follows the number rather than precedes it.
        Some syntax checks are made but there is no guarantee that all syntax
        errors will be detected.
    3.  Calculate the value of the reverse polish sequence. If there is more than
        one number on the stack at the end, or at any time there are not enough
        numbers on the stack to perform an operation an error is reported.