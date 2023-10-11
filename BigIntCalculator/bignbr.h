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


/* MAX_LEN decreased from 2500 to 2304 so that Karatsuba multiplication
can handle maximum length numbers */
#define MAX_LEN 2304        // approximately 21500 digits
#define BITS_PER_GROUP 31
#define BITS_PER_GROUP_MINUS_1 30
#define BITS_PER_INT_GROUP 31
#define HALF_INT_RANGE (1 << (BITS_PER_INT_GROUP - 1))
//#define FOURTH_INT_RANGE      0x20000000
#define HALF_INT_RANGE_U      0x40000000U
#define INT_RANGE_U 0x80000000U
#define MAX_INT_NBR ((int)((1U << BITS_PER_INT_GROUP)-1))
#define MAX_INT_NBR_U         0x7FFFFFFFU
#define LIMB_RANGE (1U<<BITS_PER_GROUP)
#define SMALL_NUMBER_BOUND 32768
#define _USING64BITS_ 1

#define MAX_VALUE_LIMB (LIMB_RANGE-1)

typedef int limb;

bool BigNbrIsZero(const limb value[], int nbrLength);
int BigNbrLen(const long long Nbr[], int nbrLen);
void ChSignBigNbr(int nbr[], int length);
void ChSignBigNbrB(int nbr[], int length);
void DivBigNbrByInt(const limb Dividend[], int divisor, limb Quotient[], int nbrLen, int *remainder);
void DivideZnumByMaxPowerOf2(int& ShRight, Znum& number);
double getMantissa(const limb *ptrLimb, int nbrLimbs);
long long JacobiSymbol(long long upper, long long lower);
void LimbstoZ(const limb number[], Znum& numberZ, int NumLen);
/* get log of Znum in base e */
double logZnum(const Znum& BigInt);
/* get log of Znum in base b */
double logZnum(const Znum& BigInt, unsigned long long b);
void ModInvZnum(const Znum &num, Znum &inv, const Znum &mod);
//void AdjustModN(Znum &Nbr, const Znum &Modulus);
//void ValuestoZ(Znum &numberZ, const int number[], int NumLen);
int modPower(int NbrMod, int Expon, int currentPrime);
void MultBigNbr(const limb Factor1[], const limb Factor2[], limb Prod[], int nbrLen);
long long PowerCheck(const Znum& BigInt, Znum& Base, long long upperBound);
int PrimalityTest(const Znum& Value, long long upperBound);
int ZtoBigNbr(int number[], Znum numberZ);
int ZtoLimbs(limb *number, Znum numberZ, int NumLen);


/* defined in sqroot.cpp */
void squareRoot(const Znum& arg, Znum& sqRoot);
bool isPerfectSquare(const Znum& arg, Znum& sqRoot);
bool isPerfectSquare(__int64 x);
bool isPerfectSquare(__int64 x, __int64& s);

/* defined in karatsuba.cpp 
result = factor1 * factor2. */
void multiply(const limb factor1[], const limb factor2[], limb result[], int len, int* ResultLen);

/* defined in baseconv.cpp */
void int2dec(char** pOutput, long long nbr);
//void Bin2Dec(const limb binary[], char* decimal, int nbrLimbs, int groupLength);
void Bin2DecV2(char** ppDecimal, const limb* binary, int nbrLimbs, int groupLength);

/* defined in modmult.cpp */
void ComputeInversePower2(const limb* value, limb* result, limb* tmp);
void GetMontgomeryParmsPowerOf2(int powerOf2);
void GetMontgomeryParms(int len);
void modmult(const limb factor1[], const limb factor2[], limb product[]);
void modmultInt(const limb factorBig[], int factorInt, limb result[]);
void AddBigNbrModNB(const limb Nbr1[], const limb Nbr2[], limb Sum[],
    const limb TestNbr[], int NumLen);
void AddBigNbrModN(const limb* Nbr1, const limb* Nbr2, limb* Sum,
    const limb* mod, int nbrLen);
#define AddBigNbrMod(Nbr1, Nbr2, Sum) AddBigNbrModN(Nbr1, Nbr2, Sum, TestNbr, NumberLength) 
void SubtBigNbrModN(const limb Nbr1[], const limb Nbr2[], limb Sum[], const limb TestNbr[], int NumLen);
#define SubtBigNbrMod(Nbr1, Nbr2, Sum) SubtBigNbrModN(Nbr1, Nbr2, Sum, TestNbr, NumberLength) 
void ModInvBigNbr(const limb num[], limb inv[], const limb mod[], int NumLen);
void modPow(const limb* base, const limb* exp, int nbrGroupsExp, limb* power);
void modPowBaseInt(int base, const limb* exp, int nbrGroupsExp, limb* power);
typedef void(*mmCback)(void);
extern mmCback modmultCallback;
extern limb* const TestNbr;
extern limb* const MontgomeryMultR2;
extern limb* const MontgomeryMultR1;
extern long long lModularMult;  // count of number of modular multiplications.
                                // used to control status display

void modmult(const Znum& a, const Znum& b, Znum& result);
void REDC(Znum& t, const Znum& T);
void GetMontgomeryParms(const Znum& Nval);
