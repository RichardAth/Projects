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

#include "main.h"
#include "expression.h"
#include "quadmodLL.h"

int groupLen;

char output[300000];

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
int SolNbr = 0;
char* ptrOutput;
static std::vector <Znum> roots;

/* interfac matcher between DA's code and expression evaluation function */
retCode ComputeExpression(const char* exprA, BigInteger* result, bool dummy) {
    std::string exp = exprA;
    Znum value;
    int asgct;
    bool multiv;
    retCode rc = ComputeExpr(exp, value, asgct, &multiv);
    *result = value;  /* copy value of expression */
    return rc;  /* pass back return code too */
}

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

/* uses global variable Aux1 */
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

#ifdef _DEBUG
/* verify that solution is valid */
bool checkSolution(const BigInteger* sol, const BigInteger* pValA,
    const BigInteger* pValB, const BigInteger* pValC, const BigInteger* pValN) {
    ValNn = *pValA;
    ValNn *= *sol;
    ValNn *= *sol; 
    ValNn += *sol * *pValB; 
    ValNn += *pValC;
    ValNn %= *pValN;   /*ValNn = (a*sol^2 +b*sol +c) modulo n */
    if (ValNn != 0) {   
        std::cout << "invalid solution: a= " << *pValA
            << "\n    b = " << *pValB
            << "\n    c = " << *pValC
            << "\n    n = " << *pValN
            << "\n  sol = " << *sol << '\n';

        return false;   /* not a valid solution */
    }
    else return true;
}
#endif

/* store solution in output buffer*/
void Solution(const BigInteger* value, const BigInteger * pValA, 
    const BigInteger * pValB, const BigInteger * pValC, const BigInteger * pValN)
{
    SolNbr++;
    int2dec(&ptrOutput, SolNbr);  /* copy solution to output buffer */
    copyStr(&ptrOutput, " x = ");
    BigInteger2Dec(&ptrOutput, value, groupLen);
    copyStr(&ptrOutput, "\n");

    if (verbose > 1) {
        printf("solution nbr %4d ", SolNbr);
        PrintBigInteger(value, 0);
        printf("\na = ");
        PrintBigInteger(pValA, 0);
        printf(", b= ");
        PrintBigInteger(pValB, 0);
        printf(", c= ");
        PrintBigInteger(pValC, 0); 
        printf(", n= ");
        PrintBigInteger(pValN, 0);
        putchar('\n');
    }
#ifdef _DEBUG
    checkSolution(value, pValA, pValB, pValC, pValN);
#endif
}

/* store solution in roots vector */
void solms(const BigInteger* value, const BigInteger* pValA,
    const BigInteger* pValB, const BigInteger* pValC, const BigInteger* pValN) {
    Znum vZ;

    BigtoZ(vZ, *value);
    roots.push_back(vZ);
    if (verbose > 1) {
        SolNbr++;
        gmp_printf("solution nbr %4d %Zd ", SolNbr, vZ);
        printf("\na = ");
        PrintBigInteger(pValA, 0);
        printf(", b= ");
        PrintBigInteger(pValB, 0);
        printf(", c= ");
        PrintBigInteger(pValC, 0);
        printf(", n= ");
        PrintBigInteger(pValN, 0);
        putchar('\n');
    }
#ifdef _DEBUG
    checkSolution(value, pValA, pValB, pValC, pValN);
#endif
}

// Solve Ax^2 + Bx + C = 0 in integers.
// uses global variables Aux0, Aux1, Aux2, discriminant, sqrtDiscriminant
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
                Solution(&Aux0, &ValA, &ValB, &ValC, &ValN);
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
            Solution(&Aux2, &ValA, &ValB, &ValC, &ValN);
        }
        //BigIntNegate(&sqrtDiscriminant, &sqrtDiscriminant);
        BigIntChSign(&sqrtDiscriminant);
        Aux1 = sqrtDiscriminant - ValB; //BigIntSubt(&sqrtDiscriminant, &ValB, &Aux1);
        Aux2 = Aux1 % Aux0;             //(void)BigIntRemainder(&Aux1, &Aux0, &Aux2);
        if (Aux2 == 0)                  //(BigIntIsZero(&Aux2))
        {      // (-sqrtDiscriminant-ValB)/(2*ValA) is integer: it is a solution.
            Aux2 = Aux1 / Aux0;         //(void)BigIntDivide(&Aux1, &Aux0, &Aux2);
            Solution(&Aux2, &ValA, &ValB, &ValC, &ValN);
        }
    }
    SolNbr = 1;
}

void textErrorQuadMod(retCode rc){
    textError(rc);
}

/* solve Quadratic modular equations of the form a⁢x² + b⁢x + c ≡ 0 (mod n) where
the integer unknown x is in the range 0 ≤ x < n. In particular, it can find
modular square roots by setting a = -1, b = 0, c = number whose root we want to 
find and n = modulus.  The modulus has to be factored, so if it is too large the 
factorisation can take an excessive length of time
Coefficients are in ValA, ValB, ValC and ValN=modulus
Uses global variable Aux0
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
    Aux0.sign = SIGN_POSITIVE;
    GcdAll = Aux0;      //CopyBigInt(&GcdAll, &Aux0);
    if (verbose > 1 && GcdAll > 1)
        std::cout << "ModulusIsNotZero: GcdAll = " << GcdAll << '\n';
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
            for (int ctr = 0; ctr < Gcd; ctr++) {
                Aux0 = ctr; //intToBigInteger(&Aux0, ctr);
                SolutionP(&Aux0, &ValA, &ValB, &ValC, &ValN);
            }
        }
        return;
    }

    /* N != 1 */

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
        std::cout << (lang ? "Coeficiente cuadrático: " : "Quadratic coefficient: ") ;
        textError(rc);
        return;
    }
    rc = ComputeExpression(bText, &ValB, false);  /* convert b from text to BigInteger */
    if (rc != retCode::EXPR_OK)
    {
        std::cout << (lang ? "Coeficiente lineal: " : "Linear coefficient: ");
        textError(rc);
        return;
    }
    rc = ComputeExpression(cText, &ValC, false);   /* convert c from text to BigInteger */
    if (rc != retCode::EXPR_OK)
    {
        std::cout << (lang ? "Término independiente: " : "Constant coefficient: ");
        textError(rc);
        return;
    }
    rc = ComputeExpression(modText, &ValN, false);   /* convert mod from text to BigInteger */
    if ((rc == retCode::EXPR_OK) && (ValN.sign == SIGN_NEGATIVE))
    {
        rc = retCode::EXPR_MODULUS_MUST_BE_NONNEGATIVE;
    }
    if (rc != retCode::EXPR_OK)
    {
        std::cout << (lang ? "Módulo: " : "Modulus: ");
        textError(rc);
        return;
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
            SetCallbacksForSolveEquation(Solution);
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

    copyStr(&ptrOutput, lang ? COPYRIGHT_SPANISH : COPYRIGHT_ENGLISH);
    copyStr(&ptrOutput, "\n");
}

/* solve quadratic modular equation a⁢x² + b⁢x + c ≡ 0 (mod n). 
expressions for a, b, c, and n are entered as text and converted to numbers,
then the equation is solved. */
int quadModEqn(const std::string &command) {
    char A[50], B[50], C[50], N[50];
    if (lang) {
        printf("Resolución de ecuaciones cuadráticas modulares \n");
        printf("Esta aplicación resuelve ecuaciones de la forma ax² + bx + c ≡ 0"
            " (mod n) donde la incógnita entera x se encuentra en el rango"
            " 0 ≤ x < n.\n");
         printf("En particular, puede hallar raíces cuadradasmodulares si ingresa "
            " a = -1, b = 0, c = número cuya raíz se desea hallar y n = módulo.\n");

         printf("ingrese el valor para A: ");
    }
    else { /* english */
        printf("Quadratic Modular Equation Solver\n");
         printf("solve equations of the form Ax² + Bx + C ≡ 0 (mod N) where"
            " the unknown integer x is in the range 0 <= x < N. \n");
        printf("In particular, it can find modular square roots by setting "
            "A = -1, B = 0, C = number whose root we want to find and N = modulus.\n");

        printf("Enter value for A: ");
    }
    fgets(A, sizeof(A), stdin);
    if (A[strlen(A) - 1] == '\n') A[strlen(A) - 1] = '\0'; /* remove trailing \n */

    printf((lang)?  "ingrese el valor para B: " : "Enter value for B: ");
    fgets(B, sizeof(B), stdin);
    if (B[strlen(B) - 1] == '\n') B[strlen(B) - 1] = '\0'; /* remove trailing \n */

    printf((lang) ? "ingrese el valor para C: " : "Enter value for C: ");
    fgets(C, sizeof(C), stdin);
    if (C[strlen(C) - 1] == '\n') C[strlen(C) - 1] = '\0'; /* remove trailing \n */

    printf((lang) ? "ingrese el valor para N: " : "Enter value for N: ");
    fgets(N, sizeof(N), stdin);
    if (N[strlen(N) - 1] == '\n') N[strlen(N) - 1] = '\0'; /* remove trailing \n */

    quadmodText(A, B, C, N, 6);
    printf("%s\n", output);
    //system("PAUSE");   /* press any key to continue */
    return 0;
}

/* find modular square root. Solve the equation given aa and m.
	x^2 ≡ aa mod m */
std::vector <Znum> ModSqrtQE(const Znum& aa, const Znum& m) {

    roots.clear();
    Znum am = aa%m;
    if (am < 0)
        am += m;  /* normalise a so it's in range 0 to m-1 */
    if (verbose > 1 && am != aa) {
        std::cout << "modsqrt(" << aa << ", " << m << ") \n changed to: \nmodsqrt("
            << am << ", " << m << ")\n";
    }

    BigInteger a, b, c, n;

    if (m == 1) {
        /* every number ≡ 0 (mod 1)*/
        roots.push_back(0);
        return roots;
    }

    a = -1;
    b = 0;
    c = am;
    n = m;
    /* solve -x² +c = 0 (mod n) */
    SetCallbacksForSolveEquation(solms);
    SolNbr = 0;
    ModulusIsNotZero(a, b, c, n);  /* get the roots */
    std::sort(roots.begin(), roots.end());
    return roots;
}

/* test quadratic modular equation solver */

/* verify that solution is valid */
bool checkSolution(const Znum &sol, const Znum &A, const Znum &B, const Znum &C, 
    const Znum &N) {
    Znum result;
    result = A *sol * sol;
    result += sol * B;
    result += C;
    result %= N;   /*ValNn = (a*sol^2 +b*sol +c) modulo n */
    if (result != 0) {
        std::cout << "invalid solution: a= " << A
            << "\n    b = " << B
            << "\n    c = " << C
            << "\n    n = " << N
            << "\n  sol = " << sol << '\n';
        std::cout << "result = " << result << " should be 0 \n";
        return false;   /* not a valid solution */
    }
    else return true;
}

void doTestsA(const std::string& params) {
    long long p1;  // number of tests; must be greater than 0
    long long p2;  // size of numbers to be factored in bits (>= 48)
    gmp_randstate_t state;
    Znum a, b, c, n;

    auto numParams = sscanf_s(params.data(), "%lld,%lld", &p1, &p2);
    if (p1 <= 0 || numParams < 1) {
        std::cout << "Use default 20 for number of tests \n";
        p1 = 20;
    }
    if (p2 < 20 || numParams < 2) {
        std::cout << "Use default 20 for number size in bits \n";
        p2 = 20;
    }
    gmp_randinit_mt(state);  // use Mersenne Twister to generate pseudo-random numbers
    gmp_randseed_ui(state, 756128234);
    /* fixed seed means that the exact same tests can be repeated provided
       that the same size of number is used each time */

    auto start = clock();	// used to measure execution time
    SetCallbacksForSolveEquation(solms);

    for (int i = 1; i <= p1; i++) {
        mpz_urandomb(ZT(a), state, p2);  // get random number, size=p2 bits
        mpz_urandomb(ZT(b), state, p2);  // get random number, size=p2 bits
        mpz_urandomb(ZT(c), state, p2);  // get random number, size=p2 bits
        mpz_urandomb(ZT(n), state, p2);  // get random number, size=p2 bits
        SolNbr = 0;

        if (n == 0) {
            ValA = a;
            ValB = b;
            ValC = c;
            SolveIntegerEquation(ValA, ValB, ValC);
            if (roots.size() > 0) {
                std::cout << "roots are: ";
                for (int j = 0; j < roots.size(); j++) {
                    std::cout << roots[j] << ", ";

                }
                putchar('\n');
                for (int j = 0; j < roots.size(); j++)
                    checkSolution(roots[j], a, b, c, n);
            }
        }
        else {
            while (true) {
                roots.clear();
                std::cout << "solve: " << a << "x² + "
                    << b << "x + " << c << " = 0 (mod " << n << ") \n";
                ValA = a;
                ValB = b;
                ValC = c;
                ValN = n;
                ModulusIsNotZero(ValA, ValB, ValC, ValN);
                if (roots.size() > 0) {
                    std::cout << "roots are: ";
                    for (int j = 0; j < roots.size(); j++) {
                        std::cout << roots[j] << ", ";
                        
                    }
                    putchar('\n'); 
                    for (int j = 0; j < roots.size(); j++)
                        checkSolution(roots[j], a, b, c, n);
                }
                if (roots.size() > 0 || i > p1)
                    break;
                i++;
                c++;       /* try again, wth only c changed */
                c %= n;
            }
        }
    }

    auto end = clock();   // measure amount of time used
    double elapsed = (double)end - start;
    PrintTimeUsed(elapsed, "All tests completed. Time used = ");
    gmp_randclear(state);  // clear state - avoid memory leakage
}