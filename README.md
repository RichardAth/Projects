This is a calculator and factorisation program based on Dario Alpern's program.
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

Changes include:

C++ rather than C

The calculator part was extensively rewritten to use GMP/MPIR extended precision
and the Boost multiprecision library. Wherever possible GMP functions were used
instead of DAs version.

FACTORIZATION

As mentioned above it uses algorithms ECM and SIQS.

The factoriser is essentially DAs program, with an interface function that converts 
GMP/MPIR extended precision numbers to DAs BigIntegers and vice versa. The progress 
messages it produces have been modified to work with a console window instead of a 
Web Browser. 

My attempt to speed up factorisation failed completely, so I took an entirely different approach.

A feature was added to allow use of Msieve as an alternative to the original ECM & SIQS
factorisation. My conclusion was that for larger numbers (> about 50 digits) Msieve
is usually significantly faster, provided that the -i option is used with Msieve.
The -e option causes greater use of ECM within Msieve. The default is to use ECM
only for factors < 15 digits. In addition I modified the function choose_max_digits
within the source file gmp_ecm.c to increase the use of ECM still further.

Use of Msieve is controlled by Msieve on/off commands. Note: msive must be installed
separately and the path to it must be correct.

Also an interface to YAFU, which is generally faster than Msieve. Again YAFU must be installed separately
and the path to it must be correct. For numbers > 95 digits YAFU neeeeds GGNFS to be installed as well.
Use of YAFU is controlled by YAFU ON/OFF commands.

Turning YAFU on turns Msieve off and vice versa. If both are off the original ECM/SIQS is used.

To build Msieve requires some components of GMP-ECM, as well as Pthreads

To build YAFU requires components of Msieve, GMP-ECM and Pthreads.

Some other unrelated stuff is also included in the solution:

Int128 - allows use of 128-bit integers with visual studio. I had intended to use this in the 
factorisation but never actually did so.

Baillie-PSW - a copy of a program I found that demonstrates various primality-testing functions


For more information see readme.txt file
