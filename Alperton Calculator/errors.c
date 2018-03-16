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

#define _CRT_SECURE_NO_DEPRECATE

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "highlevel.h"
#include "bignbr.h"
#include "expression.h"

int expressionNbr;
/* convert error code to text */
void textError(char *ptrOutput, enum eExprErr rc)
{
	/*
	currently defined error codes:
	EXPR_NUMBER_TOO_LOW =		-100,
	EXPR_NUMBER_TOO_HIGH,		-99
	EXPR_INTERM_TOO_HIGH,		-98
	EXPR_DIVIDE_BY_ZERO,        -97
	EXPR_PAREN_MISMATCH,		-96
	EXPR_SYNTAX_ERROR,			-95
	EXPR_TOO_MANY_PAREN,		-94
	EXPR_INVALID_PARAM,			-93
	EXPR_ARGUMENTS_NOT_RELATIVELY_PRIME,	-92
	EXPR_BREAK,					-91
	EXPR_OUT_OF_MEMORY,
	EXPR_CANNOT_USE_X_IN_EXPONENT,
	EXPR_DEGREE_TOO_HIGH,
	EXPR_EXPONENT_TOO_LARGE,
	EXPR_EXPONENT_NEGATIVE,				not used?
	EXPR_LEADING_COFF_MULTIPLE_OF_PRIME,  not used?
	EXPR_CANNOT_LIFT,					not used?
	EXPR_MODULUS_MUST_BE_GREATER_THAN_ONE,	-83
	EXPR_MODULUS_MUST_BE_PRIME_EXP,		not used?
	EXPR_BASE_MUST_BE_POSITIVE,			-81
	EXPR_POWER_MUST_BE_POSITIVE,		-80
	EXPR_MODULUS_MUST_BE_NONNEGATIVE,	-79
	EXPR_VAR_OR_COUNTER_REQUIRED,		-78
	EXPR_OK = 0
	*/
	switch (rc)
	{
	case EXPR_NUMBER_TOO_LOW:
		strcpy(ptrOutput, lang ? "Número muy pequeño\n" : "Number too low\n");
		break;
	case EXPR_NUMBER_TOO_HIGH:
		strcpy(ptrOutput, lang ? "Número muy grande (más de 10000 dígitos)\n" :
			"Number too high (more than 10000 digits)\n");
		break;
	case EXPR_INTERM_TOO_HIGH:
		strcpy(ptrOutput, lang ? "Número intermedio muy grande (más de 20000 dígitos\n" :
			"Intermediate number too high (more than 20000 digits)\n");
		break;
	case EXPR_DIVIDE_BY_ZERO:
		strcpy(ptrOutput, lang ? "División por cero\n" : "Division by zero\n");
		break;
	case EXPR_PAREN_MISMATCH:
		strcpy(ptrOutput, lang ? "Error de paréntesis\n" : "Parenthesis mismatch\n");
		break;
	case EXPR_SYNTAX_ERROR:
		if (lang)
		{
			strcpy(ptrOutput, "Error de sintaxis");
			if (expressionNbr > 0)
			{
				ptrOutput += strlen(ptrOutput);
				strcpy(ptrOutput, " en la expresión ");
				ptrOutput += strlen(ptrOutput);
				*ptrOutput++ = (char)(expressionNbr + '0');
				*ptrOutput++ = 0;
			}
		}
		else
		{
			strcpy(ptrOutput, "Syntax error");
			if (expressionNbr > 0)
			{
				ptrOutput += strlen(ptrOutput);
				strcpy(ptrOutput, " in expression #");
				ptrOutput += strlen(ptrOutput);
				*ptrOutput++ = (char)(expressionNbr + '0');
				*ptrOutput++ = 0;
			}
		}
		break;
	case EXPR_TOO_MANY_PAREN:
		strcpy(ptrOutput, lang ? "Demasiados paréntesis\n" : "Too many parenthesis\n");
		break;
	case EXPR_INVALID_PARAM:
		strcpy(ptrOutput, lang ? "Parámetro inválido\n" : "Invalid parameter\n");
		break;
	case EXPR_ARGUMENTS_NOT_RELATIVELY_PRIME:
		strcpy(ptrOutput, lang ? "MCD de los argumentos no es 1\n" : "GCD of arguments is not 1\n");
		break;
	case EXPR_BREAK:
		strcpy(ptrOutput, lang ? "Detenido por el usuario\n" : "Stopped by use\nr");
		break;
	case EXPR_VAR_OR_COUNTER_REQUIRED:
		if (lang)
		{
			strcpy(ptrOutput, "La expresión ");
			ptrOutput += strlen(ptrOutput);
			*ptrOutput++ = (char)(expressionNbr + '0');
			strcpy(ptrOutput, " debe incluir la variable <var>x</var> y/o el contador <var>c</var>\n");
		}
		else
		{
			strcpy(ptrOutput, "Expression #");
			ptrOutput += strlen(ptrOutput);
			*ptrOutput++ = (char)(expressionNbr + '0');
			strcpy(ptrOutput, " must include the variable <var>x</var> and/or the counter <var>c</var>\n");
		}
		break;
	case EXPR_BASE_MUST_BE_POSITIVE:
		strcpy(ptrOutput, lang ? "La base debe ser mayor que cero\n" :
			"Base must be greater than zero\n");
		break;
	case EXPR_POWER_MUST_BE_POSITIVE:
		strcpy(ptrOutput, lang ? "La potencia debe ser mayor que cero\n" :
			"Power must be greater than zero\n");
		break;
	case EXPR_MODULUS_MUST_BE_GREATER_THAN_ONE:
		strcpy(ptrOutput, lang ? "El módulo debe ser mayor que 1\n" : "Modulus must be greater than one\n");
		break;
	case EXPR_MODULUS_MUST_BE_NONNEGATIVE:
		strcpy(ptrOutput, lang ? "El módulo no debe ser negativo\n" :
			"Modulus must not be negative\n");
		break;
	default:
		sprintf(ptrOutput, "unknown error code: %d\n", (int)rc);
		break;
	}
}