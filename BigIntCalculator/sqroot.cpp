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

#include <stdlib.h>
#include "bignbr.h"

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