#pragma once

#define COPYRIGHT_SPANISH "Hecho por Dar�o Alpern. Actualizado el 18 de febrero de 2018."
#define COPYRIGHT_ENGLISH "Written by Dario Alpern. Last updated on 18 February 2018."

#ifdef __EMSCRIPTEN__
int stamp(void);
#endif
void databack(char *data);

enum eExprErr
{
	EXPR_NUMBER_TOO_LOW = -100,
	EXPR_NUMBER_TOO_HIGH,
	EXPR_INTERM_TOO_HIGH,
	EXPR_DIVIDE_BY_ZERO,
	EXPR_PAREN_MISMATCH,
	EXPR_SYNTAX_ERROR,
	EXPR_TOO_MANY_PAREN,
	EXPR_INVALID_PARAM,
	EXPR_ARGUMENTS_NOT_RELATIVELY_PRIME,
	EXPR_BREAK,
	EXPR_OUT_OF_MEMORY,
	EXPR_CANNOT_USE_X_IN_EXPONENT,
	EXPR_DEGREE_TOO_HIGH,
	EXPR_EXPONENT_TOO_LARGE,
	EXPR_EXPONENT_NEGATIVE,
	EXPR_LEADING_COFF_MULTIPLE_OF_PRIME,
	EXPR_CANNOT_LIFT,
	EXPR_MODULUS_MUST_BE_GREATER_THAN_ONE,
	EXPR_MODULUS_MUST_BE_PRIME_EXP,
	EXPR_BASE_MUST_BE_POSITIVE,
	EXPR_POWER_MUST_BE_POSITIVE,
	EXPR_MODULUS_MUST_BE_NONNEGATIVE,
	EXPR_VAR_OR_COUNTER_REQUIRED,
	EXPR_OK = 0
};

extern char inputString[1000000];
extern char output[3000000];
extern BigInteger valueX;
extern int counterC;
extern int expressionNbr;
enum eExprErr ComputeGaussianExpression(char *expr, BigInteger *ExpressionResult);
enum eExprErr ComputeExpression(char *expr, int type, BigInteger *ExpressionResult);
enum eExprErr ComputeExpr(char *expr, BigInteger *ExpressionResult);
void partition(int val, BigInteger *pResult);
void textError(char *output, enum eExprErr rc);
void PerformFactorization(const BigInteger *tofactor);

#define FACTORIZATION_FUNCTIONS 1