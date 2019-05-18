#include <cstdint>
#include <cstdlib>
#include <iostream>

#include "bignbr.h"

// Compute Nbr <- Nbr mod Modulus.
void AdjustModN(Znum &Nbr, const Znum &Modulus) {
	mpz_mod(ZT(Nbr), ZT(Nbr), ZT(Modulus));
}

Znum zN,     // (copy of Nval)
zR,      // R is a power of 2 such that R > N
zR2,     // = R ^ 2(mod N)
zNI,     // = -N^(-1) mod R i.e. NI*N ≡ −1 mod R;
zNI2;

long long zRexp;
Znum zR1;
Znum zR_1;  // = zR-1

// given the modulus N (which must be odd)
// Let R be a power of 2 such that R > N
// Compute R1  and NI using the formulas:
// R1 = R mod N
// NI = -N^(-1) mod R i.e. NI*N ≡ −1 mod R;
// NB. if N < 8 limbs only calculate NI2 as - N ^ (-1) mod 2 ^ 63
// also compute R2 = R^2(mod N), Rexp (2^Rexp = R)
// REDC Uses external global variables zR_1, zNI, zN (copy of Nval)
void GetMontgomeryParms(const Znum &Nval) {
	zN = Nval;        // save value for use in REDC and ModMultInt
	assert((zN & 1) != 0);		// N must be odd

	auto index = mpz_sizeinbase(ZT(zN), 2);
	if (index > 1) {
		// (2^Rexp = R).
		zRexp = (index +BITS_PER_GROUP-1)/BITS_PER_GROUP;	 
		zRexp *= BITS_PER_GROUP;   // exp is rounded up to a multiple of 63
		               // to match what is used with limbs
		zRexp--;
		zR = 1;
		mpz_mul_2exp(ZT(zR), ZT(zR), zRexp);   // zR is a power of 2 such that zR > zN
		zR_1 = zR - 1;              // used as a bit mask, assumes that R is a power of 2
		
		// R1 = R mod N
		zR1 = zR % Nval;
				
		mpz_invert(ZT(zNI), ZT(Nval), ZT(zR));  // zNI*Nval ≡ 1 mod zR;
		zNI = zR - zNI;				// zNI*Nval ≡ -1 mod zR ≡ (zR-1) mod zR;
		
		if (zRexp >= 8 * BITS_PER_GROUP) {
			zNI2 = zNI;
		}
		else {
			Znum temp = 1LL << 31;
			mpz_invert(ZT(zNI2), ZT(Nval), ZT(temp));
			zNI2 = temp - zNI2;
		}

		mpz_powm_ui(ZT(zR2), ZT(zR), 2, ZT(zN)); // zR2 = zR^2(mod zN)
#ifdef _DEBUG
		std::cout << "N     = " << zN
			    << "\nNI    = " << zNI
			    << "\nNI2   = " << zNI2
			    << "\nR1    = " << zR1
			    << "\nR2    = " << zR2 
			    << "\nzRexp = " << zRexp << '\n';
#endif
	}
	else {
		std::cout << "Modulus must not be zero\n";
		exit(EXIT_FAILURE);
	}
}


// The REDC algorithm :
// see https://en.wikipedia.org/wiki/Montgomery_modular_multiplication
/*
input:  Integers R and N with gcd(R, N) = 1,
(R is a power of 2 such that R > N)
Integer NI in [0, R − 1] such that N.NI ≡ −1 mod R,
Integer T in the range [0, R.N − 1]
output: Integer S in the range [0, N − 1] such that S ≡ TR^(−1) mod N
For large numbers the REDC algorithm is:
m <- ((T mod R)NI) mod R     // mod R is performed by masking out upper bits
t <- (T + mN) / R            // division by R is performed by right shift
// the bits shifted out are zeros
if t >= N then
return t - N
else
return t
end if
Conversion of a out of Montgomery form is done by computing REDC(aR mod N).*/
void REDC(Znum &t, const Znum &T) {
	static Znum m;

	mpz_and(ZT(m), ZT(T), ZT(zR_1)); // m = T (mod R)  
	m *= zNI;						 // m = (T mod R)NI
									 // mpz_mul(m, m, ZT(zNI)); 

	mpz_and(ZT(m), ZT(m), ZT(zR_1)); // m = ((T mod R)NI) mod R

	m = (T + m *zN);  // minor change from pseudocode above, ensures correct result
					  // even when t and T are the same variable
	mpz_tdiv_q_2exp(ZT(t), ZT(m), zRexp); // t = (T + m*N) /R

	if (t >= zN)
		t -= zN;	   // (ensure that t is in correct range)

					   /* we should have t ≡ T.zR^(-1)   (mod zN)
					   i.e.              t.zR ≡ T        (mod zN) */

	return;
}

/* result= a*b mod zN.
a and b must be in Montgomery form. result is also in Montgomery form */
void modmult(const Znum &a, const Znum &b, Znum &result) {
	static Znum T;
	T = a * b;  // T is the ordinary product of the two mumbers
	REDC(result, T);
}

/* Sum = Nbr1 + Nbr2 (mod m) */
void AddBigNbrModNB(const Znum &Nbr1, const Znum &Nbr2, Znum &Sum, const Znum &m) {
	Sum = Nbr1 + Nbr2;
	if (Sum > m)
		Sum -= m;
}

/* Diff = Nbr1-Nbr2 (mod m)*/
void SubtBigNbrModN(const Znum &Nbr1, const Znum &Nbr2, Znum &Diff, const Znum &m) {
	Diff = Nbr1 - Nbr2;
	if (Diff < 0)
		Diff += m;
}

/* result = FactorBig* factorInt (mod zN) */
void modmultInt(const Znum &factorBig, int factorInt, Znum &result) {
	result = (factorBig * factorInt) % zN;
}

/***********************************************************************/
/* NAME: ModInvBigNbr                                                  */
/* PURPOSE: Find the inverse of num modulo M.  num.inv ≡  1 (mod m)    */
/***********************************************************************/
void ModInvBigNbr(const Znum &num, Znum &inv, const Znum &mod) {
	mpz_invert(ZT(inv), ZT(num), ZT(mod));
}