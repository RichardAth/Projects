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


//#include "boost/multiprecision/gmp.hpp"
//typedef boost::multiprecision::mpz_int Znum;
//#define ZT(a) a.backend().data()

//#define ZisEven(a) (mpz_even_p(ZT(a)) != 0)  /* true iff a is even (works for -ve a as well) */

bool isPrime2(unsigned __int64 num);

#define MAX_LEN 2500        // approximately 20000 digits
#define BITS_PER_GROUP 31
#define BITS_PER_GROUP_MINUS_1 30
#define BITS_PER_INT_GROUP 31
#define HALF_INT_RANGE (1 << (BITS_PER_INT_GROUP - 1))
#define FOURTH_INT_RANGE      0x20000000
#define HALF_INT_RANGE_U      0x40000000U
#define MAX_INT_NBR ((int)((1U << BITS_PER_INT_GROUP)-1))
#define MAX_INT_NBR_U         0x7FFFFFFFU
#define LIMB_RANGE (1U<<BITS_PER_GROUP)
#define SMALL_NUMBER_BOUND 32768
#define _USING64BITS_ 1

//struct limb
//{
//	int x=0;      // maximum value in each limb is MAX_VALUE_LIMB
//};
#define MAX_VALUE_LIMB (LIMB_RANGE-1)

typedef int limb;

extern long long lModularMult;  // count of number of modular multiplications used to control status display
bool BigNbrIsZero(const limb value[], int nbrLength);

void squareRoot(const Znum &arg, Znum & sqRoot);
bool isPerfectSquare(const Znum &arg, Znum &sqRoot);
bool isPerfectSquare(__int64 x);
bool isPerfectSquare(__int64 x, __int64& s);
double getMantissa(const limb *ptrLimb, int nbrLimbs);


extern limb * const TestNbr;
extern limb *const MontgomeryMultR2;
extern limb *const MontgomeryMultR1;

//extern int groupLen;

void ModInvZnum(const Znum &num, Znum &inv, const Znum &mod);

void multiply(const limb factor1[], const limb factor2[], limb result[], int len, int* ResultLen);
void int2dec(char **pOutput, long long nbr);
void Bin2Dec(const limb binary[], char* decimal, int nbrLimbs, int groupLength);
void Bin2Dec(char** ppDecimal, const limb* binary, int nbrLimbs, int groupLength);
void ComputeInversePower2(const limb* value, limb* result, limb* tmp);
void GetMontgomeryParmsPowerOf2(int powerOf2);
void GetMontgomeryParms(int len);
void GetMontgomeryParms(const Znum &Nval);
void AddBigNbrModNB (const limb Nbr1[], const limb Nbr2[], limb Sum[], 
	const limb TestNbr[], int NumLen);
void AddBigNbrModN(const limb* Nbr1, const limb* Nbr2, limb* Sum,
	const limb* mod, int nbrLen);
#define AddBigNbrMod(Nbr1, Nbr2, Sum) AddBigNbrModN(Nbr1, Nbr2, Sum, TestNbr, NumberLength) 
void SubtBigNbrModN(const limb Nbr1[], const limb Nbr2[], limb Sum[], const limb TestNbr[], int NumLen);
#define SubtBigNbrMod(Nbr1, Nbr2, Sum) SubtBigNbrModN(Nbr1, Nbr2, Sum, TestNbr, NumberLength) 

void modmult(const limb factor1[], const limb factor2[], limb product[]);

void modmult(const Znum &a, const Znum &b, Znum &result);
void REDC(Znum &t, const Znum &T);
void modmultInt(const limb factorBig[], int factorInt, limb result[]);

//void AdjustModN(Znum &Nbr, const Znum &Modulus);
void ModInvBigNbr(const limb num[], limb inv[], const limb mod[], int NumLen);
void modPow(const limb* base, const limb* exp, int nbrGroupsExp, limb* power);
void modPowBaseInt(int base, const limb* exp, int nbrGroupsExp, limb* power);


//void ValuestoZ(Znum &numberZ, const int number[], int NumLen);

void ChSignBigNbr(int nbr[], int length);
void ChSignBigNbrB(int nbr[], int length);
int BigNbrLen(const long long Nbr[], int nbrLen);
void MultBigNbr(const limb* pFactor1, const limb* pFactor2, limb* pProd, int nbrLen);
void DivBigNbrByInt(const limb Dividend[], int divisor, limb Quotient[], int nbrLen);
int modPower(int NbrMod, int Expon, int currentPrime);

int ZtoLimbs(limb *number, Znum numberZ, int NumLen);
int ZtoBigNbr(int number[], Znum numberZ);
void LimbstoZ(const limb number[], Znum& numberZ, int NumLen);

void DivideBigNbrByMaxPowerOf2(int &ShRight, Znum &number);
int PrimalityTest(const Znum &Value, long long upperBound);
long long JacobiSymbol(long long upper, long long lower);
long long PowerCheck(const Znum &BigInt, Znum &Base, long long upperBound);
double logZnum(const Znum &BigInt);
double logZnum(const Znum &BigInt, unsigned long long b);

typedef void(*mmCback)(void);
extern mmCback modmultCallback;
