#pragma once

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

#include "mpir.h"
#include "boost/multiprecision/gmp.hpp" 
#include "bignbr.h"
#include "showtime.h"

//#ifdef __EMSCRIPTEN__
//void getCunn(char *url, char *factorsFromServer);
//#endif
struct sFactors
{
	int *ptrFactor;
	int multiplicity;
	int upperBound;
};
typedef boost::multiprecision::mpz_int Znum;


void FactoringSIQSx(const limb *pNbrToFactor, limb *pFactor);
extern int lang;


/* access underlying mpz_t inside an bigint */
#define ZT(a) a.backend().data()

bool factorise(const Znum numberZ, std::vector <Znum> &factorlist,
	std::vector <int> &exponentlist, Znum Quad[]);


