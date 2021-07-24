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

#include "pch.h"
#include "factor.h"

void squareRoot(const Znum &arg, Znum & sqRoot) {
	assert(arg >= 0);  // no imaginary numbers allowed here!!
	mpz_sqrt(ZT(sqRoot), ZT(arg));
}

/* return true iff arg is a perfect square */
bool isPerfectSquare(const Znum &arg, Znum &sqRoot) {
	Znum rem;
	assert(arg >= 0);  // no imaginary numbers allowed here!!
	mpz_sqrtrem(ZT(sqRoot), ZT(rem), ZT(arg));
	return rem == 0;  // true if arg is a perfect square
}

/* faster version for small n. Do some checks before sqrt
because sqrt function is relatively slow. */
bool isPerfectSquare(__int64 x) {
	if (x < 0) return false;
	if (x == 0) return true;
	//while ((x & 0x3) == 0) 
	//	x >>= 2;   // divide by largest possible power of 4

	/* can use _BitScanForward64 instead (slightly faster) */
	unsigned long ix;
	auto result = _BitScanForward64(&ix, x);  // for gcc compiler use __builtin_ctzll instead
	if ((ix & 1) == 1)
		return false;     // if x contains an odd number of 2 factors it is 
						  // not a perfect square
	x >>= ix;

	/* at this stage, x must be odd for a perfect square,
	so sqrt(x) is odd if x is a perfect square, in which case
	(sqrt(x))(mod 8) is 1 3 5 or 7. So x (mod 8) has to be 1 */
	if ((x & 0x7) != 1)
		return false;  // not a perfect square

	long long s = llround(sqrt(x));  // slightly faster than intSqrt
	//long long s = intSqrt(x);
	return (s * s == x);
}