#define _CRT_SECURE_NO_DEPRECATE
#define __EMSCRIPTEN__
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "highlevel.h"
#include "bignbr.h"
#include "expression.h"
#include "batch.h"
#include "showtime.h"


int valuesProcessed;         // kludge *****************
static char input[160] = { 0 };
static char helpmsg[] =
"You can enter expressions that use the following operators, functions and parentheses:\n"

"+ for addition\n"
"- for subtraction\n"
"* for multiplication\n"
"/ for integer division\n"
"%% for modulus(remainder of the integer division)\n"
"^ or ** for exponentiation(the exponent must be greater than or equal to zero).\n"
"<, == , >; <= , >= , != for comparisons.The operators return zero for false and -1 for true.\n"
"AND, OR, XOR, NOT for binary logic.\n"
"SHL: Shift left the number of bits specified on the right operand.\n"
"SHR : Shift right the number of bits specified on the right operand.\n"
"n!: factorial(n must be greater than or equal to zero).\n"
"p# : primorial(product of all primes less or equal than p).\n"
"B(n) : Previous probable prime before n\n"
"F(n) : Fibonacci number Fn\n"
"L(n) : Lucas number Ln = F(n-1) + F(n+1)\n"
"N(n) : Next probable prime after n\n"
"P(n) : Unrestricted Partition Number(number of decompositions of n into sums of integers without regard to order).\n"
"Gcd(m, n)       : Greatest common divisor of these two integers.\n"
"Modinv(m, n)    : inverse of m modulo n, only valid when gcd(m, n) = 1.\n"
"Modpow(m, n, r) : finds m^n modulo r.\n"
"Totient(n)      : finds the number of positive integers less than n which are relatively prime to n.\n"
"IsPrime(n)      : returns zero if n is not probable prime, -1 if it is.\n"
"NumDivs(n)      : Number of positive divisors of n either prime or composite.\n"
"SumDivs(n)      : Sum of positive divisors of n either prime or composite.\n"
"NumDigits(n, r) : Number of digits of n in base r.\n"
"SumDigits(n, r) : Sum of digits of n in base r.\n"
"RevDigits(n, r) : finds the value obtained by writing backwards the digits of n in base r.\n";

static char ayuda[] =
"Puedes ingresar expresiones que usen los siguientes operadores y paréntesis:\n"

"+ para suma\n"
"- para resta\n"
"* para multiplicación\n"
"/ para división entera\n"
"% para el resto de la división entera\n"
"^ o ** para exponenciación(el exponente debe ser mayor o igual que cero).\n"
"<, == , >; <= , >= , != para comparaciones.Los operadores devuelven cero si es falso y - 1 si es verdadero.\n"
"AND, OR, XOR, NOT para lógica binaria.\n"
"SHL  : Desplazar a la izquierda la cantidad de bits indicada en el operando derecho.\n"
"SHR  : Desplazar a la derecha la cantidad de bits indicada en el operando derecho.\n"
"n!   : factorial(n debe ser mayor o igual que cero).\n"
"p#   : primorial(producto de todos los primos menores o iguales a p).\n"
"B(n) : Número probablemente primo anterior a n\n"
"F(n) : Número de Fibonacci Fn\n"
"L(n) : Número de Lucas Ln = F(n-1) + F(n+1)\n"
"N(n) : Número probablemente primo posterior a n\n"
"P(n) : particiones irrestrictas(cantidad de descomposiciones de n en sumas de números enteros sin tener en cuenta el orden).\n"
"Gcd(m, n)       : Máximo común divisor de estos dos números enteros.\n"
"Modinv(m, n)    : inverso de m modulo n, sólo válido cuando gcd(m, n) = 1.\n"
"Modpow(m, n, r) : halla m^n módulo r.\n"
"Totient(n)      : cantidad de enteros positivos menores que n coprimos con n.\n"
"IsPrime(n)      : returna cero si n no es un primo probable y - 1 si lo es.\n"
"NumDivs(n)      : cantidad de divisores positivos de n primos o compuestos.\n"
"SumDivs(n)      : suma de divisores positivos de n primos o compuestos.\n"
"NumDigits(n, r) : cantidad de dígitos de n en base r.\n"
"SumDigits(n, r) : suma de dígitos de n en base r.\n"
"RevDigits(n, r) : halla el valor que se obtiene escribiendo para atrás los dígitos de n en base r.\n";

/* tofactortext is the 'expression to evaluate' e.g. the number to factor
the number to factor is converted to a biginteger, 
ptrOutput is a pointer to a pointer into the output buffer*/
enum eExprErr BatchProcessing(char *tofactorText, BigInteger *tofactor, char **ptrOutput) {
	enum eExprErr rv = EXPR_OK;

	if (tofactorText != NULL) {
		printf("Expression is: %s = ", tofactorText);
		rv = ComputeExpression(tofactorText, 1, tofactor);
	}
	else {
		if (lang == 0) {
			printf("enter expression to be processed, or HELP\n");
		}
		else
			printf("ingrese la expresión para ser procesada, o AYUDA\n");
		gets_s(input, sizeof(input));
		if (strcmp(input, "help") == 0 || strcmp(input, "HELP") == 0
			|| strcmp(input, "ayuda") == 0 || strcmp(input, "AYUDA") == 0) {
			if (lang == 0) {
				printf(helpmsg);
				printf("enter expression to be processed\n");
			}
			else {
				printf(ayuda);
				printf("ingrese la expresión para ser procesada\n");
			}
			gets_s(input, sizeof(input));
		}
		clock_t beginTime = clock();
		originalTenthSecond = tenths();
		rv = ComputeExpression(input, 1, tofactor);
	}
	if (rv != EXPR_OK) {
		textError(*ptrOutput, rv);
	}
	else {     // output value of expression
		BigInteger2Dec(tofactor, *ptrOutput, 6);
	}
	strcat(*ptrOutput, "\n");
	*ptrOutput += strlen(*ptrOutput);
	valuesProcessed = 1;
	return rv;
}