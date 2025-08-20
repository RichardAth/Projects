#pragma once
/* R3.cpp */
struct uint128_t {
	uint64_t hi;
	uint64_t lo;
};
bool isPrime2(unsigned __int64 num);
size_t inverseTotient(__int64 n, std::vector<unsigned __int64>** result, bool debug,
	int level, bool dump);
long long R4alt(int num);

// calculate a^n%mod using 'bigints'   
Znum             modPower(const Znum& a, const Znum& n, const Znum& mod);
unsigned __int64 modPowerBi(const Znum& a, const Znum& n, unsigned __int64 mod);
constexpr __int64 power(const __int64 x, unsigned int n);
Znum              power(const Znum& x, unsigned long long n);
int jacobi(__int64 k, unsigned __int64 n);
int jacobi(const Znum& k, const Znum& n);
uint128_t divide_uint128_by_uint64(uint128_t dividend, uint64_t divisor,
	uint64_t* remainder);

Znum             modMult(const Znum& a, const Znum& b, const Znum& mod);
extern unsigned __int64 R2(const unsigned __int64 n);
extern unsigned __int64 R3(__int64 n);
/* get prime factors of tnum, using trial division */
unsigned int primeFactors(unsigned __int64 tnum, factorsS& f);