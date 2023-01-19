#pragma once
//
// This file is part of Alpertron Calculators.
//
// Copyright 2017-2021 Dario Alejandro Alpern
//
// Alpertron Calculators is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Alpertron Calculators is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Alpertron Calculators.  If not, see <http://www.gnu.org/licenses/>.
//
#include "pch.h"

typedef void (*pShowSolutionsModPrime)(int factorIndex, int expon,
	const Znum* pIncrement);
typedef void (*pSolution)(Znum* value);
typedef void (*pShowNoSolsModPrime)(int expon);

void SolveEquation(Znum* pValA, const Znum* pValB,
	const Znum* pValC, Znum* pValN,
	Znum* GcdAll, Znum* pValNn);
void SetCallbacksForSolveEquation(pSolution solutionCback,
	pShowSolutionsModPrime showSolutionsModPrime, pShowNoSolsModPrime showNoSolsModPrime);

