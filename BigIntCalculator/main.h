#pragma once

#include "boost/multiprecision/gmp.hpp" 
typedef boost::multiprecision::mpz_int Znum;

/* defined in main.cpp */
long long MulPrToLong(const Znum& x);

/* ComputeNumDigits(n,r): Number of digits of n in base r. */
long long ComputeNumDigits(const Znum& n, const Znum& radix);
void ShowLargeNumber(const Znum& Bi_Nbr, int digitsInGroup, bool size, bool hex);
void writeIni(void);
void generatePrimes(unsigned long long int max_val);
unsigned long long llSqrt(const unsigned long long n);
bool getBit(const unsigned long long int x, const bool array[]);


extern bool* primeFlags;
extern unsigned long long* primeList;
extern unsigned int prime_list_count;
extern unsigned long long int primeListMax;
extern int verbose;
extern int lang;    // 0 English, 1 = Spanish

/* defined in various places: not in main.cpp */
extern std::string YafuPath;
extern std::string yafuprog;
extern std::string MsievePath;
extern std::string MsieveProg;
extern bool breakSignal;
extern std::vector <Znum> roots;   /* used by functions that return multiple values */
extern bool msieve;
extern bool yafu;

void delfile(const std::string& path, const char* FileName);
bool isPerfectSquare(const Znum &num);
void msieveParam(const std::string& expupper);   /*process Msieve commands */
void yafuParam(const std::string& command);      /*process YAFU commands */
void printvars(std::string name);
void doTests9(void);  /* modular square root test */

// calculate a^n%mod   
unsigned __int64 modPowerLL(unsigned __int64 a, unsigned __int64 n,
	unsigned __int64 mod);

// calculate a^n%mod using 'bigints'   
Znum             modPower(const Znum& a, const Znum& n, const Znum& mod);
unsigned __int64 modPowerBi(const Znum& a, const Znum& n, unsigned __int64 mod);
constexpr __int64 power(const __int64 x, unsigned int n);
Znum              power(const Znum& x, unsigned long long n);
int jacobi(__int64 k, unsigned __int64 n);
int jacobi(const Znum& k, const Znum& n); 
unsigned __int64 modMult(unsigned __int64 a, unsigned __int64 b, unsigned __int64 mod);
Znum             modMult(const Znum& a, const Znum& b, Znum mod);

unsigned __int64 R3h(Znum n);
Znum Hclassno12(const Znum& n);
Znum classno(const Znum& n, int flag);
Znum tau(const Znum& n);
Znum stirling(const Znum& n, const Znum& m, const Znum& flag);