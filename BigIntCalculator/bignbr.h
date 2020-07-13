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

#ifdef __GNUC__
#include <gmp.h>
typedef unsigned long mpir_ui;  // for compatibility with MPIR
#else
#include <mpir.h>
#endif
#include "boost/multiprecision/gmp.hpp"
#define ZT(a) a.backend().data()
#define ZisEven(a) (mpz_even_p(ZT(a)) != 0)  /* true iff a is even (works for -ve a as well) */
typedef boost::multiprecision::mpz_int Znum;

long long MulPrToLong(const Znum &x);

extern bool *primeFlags;
extern unsigned long long *primeList;
extern unsigned int prime_list_count;
extern unsigned long long int primeListMax;
bool isPrime2(unsigned __int64 num);

#define MAX_LEN 2500        // approximately 20000 digits
#define BITS_PER_GROUP 31
#define BITS_PER_INT_GROUP 31
#define HALF_INT_RANGE (1 << (BITS_PER_INT_GROUP - 1))
#define MAX_INT_NBR ((int)((1U << BITS_PER_INT_GROUP)-1))
#define LIMB_RANGE (1U<<BITS_PER_GROUP)
#define SMALL_NUMBER_BOUND 32768
#define _USING64BITS_ 1

struct limb
{
	int x=0;      // maximum value in each limb is MAX_VALUE_LIMB
};
#define MAX_VALUE_LIMB (LIMB_RANGE-1)

extern long long lModularMult;  // count of number of modular multiplications
bool BigNbrIsZero(const limb *value, int nbrLength);
//void squareRoot(const limb *argument, limb *sqRoot, int len, int *pLenSqRoot);
void squareRoot(const Znum &arg, Znum & sqRoot);
bool isPerfectSquare(const Znum &arg, Znum &sqRoot);
bool isPerfectSquare(__int64 x);
double getMantissa(const limb *ptrLimb, int nbrLimbs);


extern limb * const TestNbr;
extern limb *const MontgomeryMultR2;
extern limb *const MontgomeryMultR1;

//extern int groupLen;

void ModInvBigNbr(Znum &num, Znum &inv, Znum &mod);
//void FactoringSIQSx(const Znum &NbrToFactor, Znum &Factor);
void FactoringSIQS(const Znum &NbrToFactor, Znum &Factor);
void multiply(const limb *factor1, const limb *factor2, limb *result, int len, int *ResultLen);
void int2dec(char **pOutput, long long nbr);
void Bin2Dec(const limb *binary, char *decimal, int nbrLimbs, int groupLength);
void GetMontgomeryParms(int len);
void GetMontgomeryParms(const Znum &Nval);
void AddBigNbrModNB (const limb Nbr1[], const limb Nbr2[], limb Sum[], const limb TestNbr[], int NumLen);
void SubtBigNbrModN(const limb Nbr1[], const limb Nbr2[], limb Sum[], const limb TestNbr[], int NumLen);

void modmult(const limb *factor1, const limb *factor2, limb *product);

void modmult(const Znum &a, const Znum &b, Znum &result);
void REDC(Znum &t, const Znum &T);
void modmultInt(const limb *factorBig, int factorInt, limb *result);

void AdjustModN(Znum &Nbr, const Znum &Modulus);
void ModInvBigNbr(const limb *num, limb *inv, const limb *mod, int NumLen);

void ValuestoZ(Znum &numberZ, const int number[], int NumLen);

void ChSignBigNbr(int nbr[], int length);
void ChSignBigNbrB(int nbr[], int length);
int BigNbrLen(const long long Nbr[], int nbrLen);
void AddBigNbrModN(const Znum &Nbr1, const Znum &Nbr2, Znum &Diff, const Znum &Mod);
void SubtractBigNbrModN(const Znum &Nbr1, const Znum &Nbr2, Znum &Diff, const Znum &Mod);
void DivBigNbrByInt(const int Dividend[], int divisor, int Quotient[], int nbrLen);
mpir_ui RemDivBigNbrByInt(const Znum &Dividend, mpir_ui divisor);
void MultBigNbrModN(const Znum &Nbr1, const Znum &Nbr2, Znum &Prod, const Znum &Mod);
void MultBigNbrByIntModN(const Znum &Nbr1, int Nbr2, Znum &Prod, const Znum &Mod);
int modPower(int NbrMod, int Expon, int currentPrime);
// calculate a^n%mod using 'bigints'   
unsigned __int64 modPower(unsigned __int64 a, unsigned __int64 n,
	unsigned __int64 mod);
int ZtoLimbs(limb *number, Znum numberZ, int NumLen);
int ZtoBigNbr(int number[], Znum numberZ);
void LimbstoZ(const limb *number, Znum &numberZ, int NumLen);

void DivideBigNbrByMaxPowerOf2(int &ShRight, Znum &number);
int PrimalityTest(const Znum &Value, long long upperBound);
long long JacobiSymbol(long long upper, long long lower);
long long PowerCheck(const Znum &BigInt, Znum &Base, long long upperBound);
double logBigNbr(const Znum &BigInt);
double logBigNbr(const Znum &BigInt, unsigned long long b);

typedef void(*mmCback)(void);
extern mmCback modmultCallback;
