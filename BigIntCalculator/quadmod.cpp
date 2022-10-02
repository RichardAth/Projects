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
#include "bignbr.h"
#include "bigint.h"


#include "main.h"
#include "expression.h"
#include "factor.h"
#include "quadmodLL.h"

int groupLen;

char output[300000];

/* external functions */
void textError(retCode rc);
retCode ComputeExpr(const std::string& expr, Znum& Result, int& asgCt, 
    bool* multiV = nullptr);


retCode ComputeExpression(const char* exprA, BigInteger *result, bool dummy) {
    std::string exp = exprA;
    Znum value;
    int asgct;
    bool multiv;
    retCode rc = ComputeExpr(exp, value, asgct, &multiv);
    *result = value;
    return rc;
}


static BigInteger discriminant;
static BigInteger sqrtDiscriminant;
static BigInteger Aux0;
static BigInteger Aux1;
static BigInteger Aux2;
static BigInteger GcdAll;
static BigInteger ValNn;

BigInteger ValA;
BigInteger ValB;
BigInteger ValC;
BigInteger ValN;
int SolNbr;
char* ptrOutput;
extern int factorsMod[20000];

static int Show(const BigInteger* num, const char* str, int t)
{
    if (!BigIntIsZero(num))
    {     // num is not zero.
        if (((t & 1) != 0) && (num->sign == SIGN_POSITIVE))
        {
            copyStr(&ptrOutput, " +");
        }
        if (num->sign == SIGN_NEGATIVE)
        {
            copyStr(&ptrOutput, " -");
        }
        if ((num->nbrLimbs != 1) || (num->limbs[0] != 1))
        {    // num is not 1 or -1.
            *ptrOutput = ' ';
            ptrOutput++;
            Bin2Dec(&ptrOutput, num->limbs, num->nbrLimbs, groupLen);
        }
        copyStr(&ptrOutput, str);
        return t | 1;
    }
    return t;
}

void Show1(const BigInteger* num, int t)
{
    int u = Show(num, "", t);
    if (((u & 1) == 0) || ((num->nbrLimbs == 1) && (num->limbs[0] == 1)))
    {
        *ptrOutput = ' ';
        ptrOutput++;
        Aux1 = *num; // CopyBigInt(&Aux1, num);
        Aux1.sign = SIGN_POSITIVE;
        BigInteger2Dec(&ptrOutput, &Aux1, groupLen);
    }
}

/* store solution */
void Solution(BigInteger* value)
{
    SolNbr++;
    int2dec(&ptrOutput, SolNbr);
    copyStr(&ptrOutput, " x = ");
    BigInteger2Dec(&ptrOutput, value, groupLen);
    copyStr(&ptrOutput, "\n");

    /* temp */
    printf("solution nbr %4d ", SolNbr);
    PrintBigInteger(value, 0);
    putchar('\n');
}

// Solve Ax^2 + Bx + C = 0 in integers.
// uses global variables ValA, ValB, ValC
static void SolveIntegerEquation(const BigInteger& ValA, const BigInteger& ValB, 
    const BigInteger& ValC) {
    if (ValA == 0)           //(BigIntIsZero(&ValA))
    {      // Not Quadratic Equation
        if (ValB == 0)       // (BigIntIsZero(&ValB))
        {    // Constant Equation
            if (ValC == 0) //(BigIntIsZero(&ValC))
            {  // 0 = 0
                copyStr(&ptrOutput, "The equation is satisfied by any integer x \n");
            }
            else {
                return;  /* no solution */
            }
        }
        else {        // Linear Equation: ValB * x + ValC = 0.
            Aux0 = ValC % ValB; //(void)BigIntRemainder(&ValC, &ValB, &Aux0);
            if (Aux0 == 0)      //(BigIntIsZero(&Aux0))
            {      // ValC is multiple of ValB: solution is -ValC / ValB.
                Aux0 = ValC / ValB; //(void)BigIntDivide(&ValC, &ValB, &Aux0);
                Aux0 = -Aux0;       //BigIntNegate(&Aux0, &Aux0);
                Solution(&Aux0);
            }
            else {      // No solutions on integers.
                return;
            }
        }
    }
    else
    {                       // Quadratic Equation
      // The discriminant is b^2 - 4.a.c  = (ValB * ValB - 4 * ValA * ValC).
        // the soulution formula is (-b +/- sqrt(b^2 - 4.a.c))/2a
        Aux0 = ValB * ValB;         // (void)BigIntMultiply(&ValB, &ValB, &Aux0);
        Aux1 = ValA * ValC;         // (void)BigIntMultiply(&ValA, &ValC, &Aux1);
        Aux1 *= 4;                  // multint(&Aux1, &Aux1, 4);
        discriminant = Aux0 - Aux1; // BigIntSubt(&Aux0, &Aux1, &discriminant);
        if (discriminant.sign == SIGN_NEGATIVE)
        {        // No integer solutions.
            return;
        }
        // Set sqrtDiscriminant to square root of discriminant.
        //squareRoot(discriminant.limbs, sqrtDiscriminant.limbs, discriminant.nbrLimbs, &sqrtDiscriminant.nbrLimbs);
        sqrtDiscriminant = discriminant.sqRoot();
        sqrtDiscriminant.sign = SIGN_POSITIVE;
        Aux0 = sqrtDiscriminant * sqrtDiscriminant; // (void)BigIntMultiply(&sqrtDiscriminant, &sqrtDiscriminant, &Aux0);
        if (Aux0 != discriminant)
            //BigIntSubt(&Aux0, &discriminant, &Aux0);
            //if (!BigIntIsZero(&Aux0))
        {  // discriminant has no integer square root.
            return;
        }
        Aux0 = ValA * 2;                // multint(&Aux0, &ValA, 2);    // Denominator
        Aux1 = sqrtDiscriminant - ValB; // BigIntSubt(&sqrtDiscriminant, &ValB, &Aux1);
        Aux2 = Aux1 % Aux0;             //(void)BigIntRemainder(&Aux1, &Aux0, &Aux2);
        if (Aux2 == 0)                  //(BigIntIsZero(&Aux2))
        {      // (sqrtDiscriminant-ValB)/(2*ValA) is integer: it is a solution.
            Aux2 = Aux1 / Aux0;        //(void)BigIntDivide(&Aux1, &Aux0, &Aux2);
            Solution(&Aux2);
        }
        //BigIntNegate(&sqrtDiscriminant, &sqrtDiscriminant);
        BigIntChSign(&sqrtDiscriminant);
        Aux1 = sqrtDiscriminant - ValB; //BigIntSubt(&sqrtDiscriminant, &ValB, &Aux1);
        Aux2 = Aux1 % Aux0;             //(void)BigIntRemainder(&Aux1, &Aux0, &Aux2);
        if (Aux2 == 0)                  //(BigIntIsZero(&Aux2))
        {      // (-sqrtDiscriminant-ValB)/(2*ValA) is integer: it is a solution.
            Aux2 = Aux1 / Aux0;         //(void)BigIntDivide(&Aux1, &Aux0, &Aux2);
            Solution(&Aux2);
        }
    }
    SolNbr = 1;
}

void textErrorQuadMod(char** pptrOutput, retCode rc)
{
    if (rc == retCode::EXPR_MODULUS_MUST_BE_NONNEGATIVE)
    {
        copyStr(pptrOutput, lang ? "No debe ser negativo" :
            "Must not be negative");
    }
    else
    {
        textError(rc);
    }
}

/* solve Quadratic modular equations of the form a⁢x² + b⁢x + c ≡ 0 (mod n) where
the integer unknown x is in the range 0 ≤ x < n. In particular, it can find
modular square roots by setting a = -1, b = 0, c = number whose root we want to find
and n = modulus.  The modulus has to be factored, so if it is too large
the factorisation can take an excessive length of time
Coefficients are in ValA, ValB, ValC and ValN=modulus
*/
static void ModulusIsNotZero(BigInteger& ValA, BigInteger& ValB, BigInteger& ValC, BigInteger& ValN) {
    ValN.sign = SIGN_POSITIVE;
    ValA %= ValN; // (void)BigIntRemainder(&ValA, &ValN, &ValA);
    ValB %= ValN; // (void)BigIntRemainder(&ValB, &ValN, &ValB);
    ValC %= ValN; //  (void)BigIntRemainder(&ValC, &ValN, &ValC);
    BigIntGcd(ValA, ValB, Aux0);
    Aux0.sign = SIGN_POSITIVE;
    BigIntGcd(ValC, Aux0, GcdAll); /* get GCD of (a, b, c)*/
    GcdAll.sign = SIGN_POSITIVE;
    Aux0 = ValC % GcdAll;   // (void)BigIntRemainder(&ValC, &GcdAll, &Aux0);
    if (Aux0 != 0)          //(!BigIntIsZero(&Aux0))
    {  // ValC must be multiple of gcd(ValA, ValB).
       // Otherwise go out because there are no solutions.
        return;
    }
    BigIntGcd(ValN, GcdAll, Aux0);
    GcdAll = Aux0;      //CopyBigInt(&GcdAll, &Aux0);
    // Divide all coefficients by gcd(ValA, ValB, ValC, ValN).
    ValA /= GcdAll;     //(void)BigIntDivide(&ValA, &GcdAll, &ValA);
    ValB /= GcdAll;     //(void)BigIntDivide(&ValB, &GcdAll, &ValB);
    ValC /= GcdAll;     //(void)BigIntDivide(&ValC, &GcdAll, &ValC);
    ValN /= GcdAll;     //(void)BigIntDivide(&ValN, &GcdAll, &ValN);
    ValNn = ValN;       //CopyBigInt(&ValNn, &ValN);
    if (ValN == 1)      //((ValNn.nbrLimbs == 1) && (ValNn.limbs[0] == 1))
    {     // All values from 0 to GcdAll - 1 are solutions.
        if (GcdAll > 5)   //((GcdAll.nbrLimbs > 1) || (GcdAll.limbs[0] > 5))
        {
            copyStr(&ptrOutput, "All values of x between 0 and ");
            GcdAll--;    //addbigint(&GcdAll, -1);
            BigInteger2Dec(&ptrOutput, &GcdAll, groupLen);
            copyStr(&ptrOutput, " are solutions.\n");
        }
        else
        {
            long long Gcd;
            Gcd = GcdAll.lldata();   /* change type */
            for (int ctr = 0; ctr < Gcd; ctr++)
            {
                Aux0 = ctr; //intToBigInteger(&Aux0, ctr);
                Solution(&Aux0);
            }
        }
        return;
    }

    /* N != 1 */
    SetCallbacksForSolveEquation(Solution, NULL, NULL);
    SolveEquation(&ValA, &ValB, &ValC, &ValN, &GcdAll, &ValNn);
}

/* solve Quadratic modular equations of the form a⁢x² + b⁢x + c ≡ 0 (mod n) where
the integer unknown x is in the range 0 ≤ x < n. In particular, it can find
modular square roots by setting a = -1, b = 0, c = number whose root we want to find
and n = modulus. The parameters are supplied as text, which can be numbers or
numerical expressions. The modulus has to be factored, so if it is too large
the factorisation can take an excessive length of time */
void quadmodText(const char* aText, const char* bText, const char* cText,
    const char* modText, int groupLength)
{
    groupLen = groupLength; /* used for formatting the output */
    char* ptrBeginSol;
    retCode rc;     /* return code*/
    ptrOutput = output;   /* buffer for output */

    rc = ComputeExpression(aText, &ValA, false);  /* convert a from text to BigInteger */
    if (rc != retCode::EXPR_OK)
    {
        copyStr(&ptrOutput, lang ? "Coeficiente cuadrático: " : "Quadratic coefficient: ");
        textErrorQuadMod(&ptrOutput, rc);
        //copyStr(&ptrOutput, "/n");
    }
    rc = ComputeExpression(bText, &ValB, false);  /* convert b from text to BigInteger */
    if (rc != retCode::EXPR_OK)
    {
        copyStr(&ptrOutput, lang ? "Coeficiente lineal: " : "Linear coefficient: ");
        textErrorQuadMod(&ptrOutput, rc);
        copyStr(&ptrOutput, "/n");
    }
    rc = ComputeExpression(cText, &ValC, false);   /* convert c from text to BigInteger */
    if (rc != retCode::EXPR_OK)
    {
        copyStr(&ptrOutput, lang ? "Término independiente: " : "Constant coefficient: ");
        textErrorQuadMod(&ptrOutput, rc);
        copyStr(&ptrOutput, "/n");
    }
    rc = ComputeExpression(modText, &ValN, false);   /* convert mod from text to BigInteger */
    if ((rc == retCode::EXPR_OK) && (ValN.sign == SIGN_NEGATIVE))
    {
        rc = retCode::EXPR_MODULUS_MUST_BE_NONNEGATIVE;
    }
    if (rc != retCode::EXPR_OK)
    {
        copyStr(&ptrOutput, lang ? "Módulo: " : "Modulus: ");
        textErrorQuadMod(&ptrOutput, rc);
        copyStr(&ptrOutput, "/n");
    }
    if (ptrOutput == output)
    {    // No errors found. send ax² + bx +c = 0 (mod n) to output buffer.
        //  values of a, b , c and n are displayed in decimal with appropriate signs 

        int u = Show(&ValA, "x²", 2);  /* send a as text to output, then "x²" */
        u = Show(&ValB, " x", u);      /* send +b or -b as text, then "x" */
        Show1(&ValC, u);               /* send +c or -c as text */
        copyStr(&ptrOutput, " = 0 (mod ");
        BigInteger2Dec(&ptrOutput, &ValN, groupLen);
        copyStr(&ptrOutput, ")\n");

        SolNbr = 0;
        ptrBeginSol = ptrOutput;  /* save output pointer */
        copyStr(&ptrOutput, "\n");
        if (ValN == 0) {
            // Mod zero => Equation in integer numbers
            SolveIntegerEquation(ValA, ValB, ValC);
        }
        else {
            ModulusIsNotZero(ValA, ValB, ValC, ValN);
        }
        if (SolNbr == 0)
        {
            ptrOutput = ptrBeginSol;
            copyStr(&ptrOutput, lang ? "\nNo hay soluciones.\n" : "\nThere are no solutions.\n");
        }
        else
        {
            copyStr(&ptrOutput, "\n");
        }
    }
    //copyStr(&ptrOutput, "<p>");
    copyStr(&ptrOutput, lang ? COPYRIGHT_SPANISH : COPYRIGHT_ENGLISH);
    copyStr(&ptrOutput, "\n");
}


int quadModEqn(const std::string &command) {
    char A[50], B[50], C[50], N[50];
    printf("Quadratic Modular Equation Solver\n");
    printf("solve equations of the form Ax² + Bx + C = 0 (mod N) where"
        " the unknown integer x is in the range 0 <= x < N. \n");
    printf("In particular, it can find modular square roots by setting "
        "A = -1, B = 0, C = number whose root we want to find and N = modulus.\n");

    printf("Enter value for A: ");
    fgets(A, sizeof(A), stdin);
    if (A[strlen(A) - 1] == '\n') A[strlen(A) - 1] = '\0'; /* remove trailing \n */

    printf("Enter value for B: ");
    fgets(B, sizeof(B), stdin);
    if (B[strlen(B) - 1] == '\n') B[strlen(B) - 1] = '\0'; /* remove trailing \n */

    printf("Enter value for C: ");
    fgets(C, sizeof(C), stdin);
    if (C[strlen(C) - 1] == '\n') C[strlen(C) - 1] = '\0'; /* remove trailing \n */

    printf("Enter value for N: ");
    fgets(N, sizeof(N), stdin);
    if (N[strlen(N) - 1] == '\n') N[strlen(N) - 1] = '\0'; /* remove trailing \n */

    quadmodText(A, B, C, N, 6);
    printf("%s\n", output);
    //system("PAUSE");   /* press any key to continue */
    return 0;
}



