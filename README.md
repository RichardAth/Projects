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
