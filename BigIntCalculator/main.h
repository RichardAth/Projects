#pragma once


/* defined in main.cpp */
/* error and return codes, errors are -ve, OK is 0, FAIL is +1 */
enum class retCode
{
	NUMBER_TOO_LOW = -100,
	NUMBER_TOO_HIGH,
	INTERIM_TOO_HIGH,
	DIVIDE_BY_ZERO,
	PAREN_MISMATCH,
	SYNTAX_ERROR,
	TOO_MANY_PAREN,
	INVALID_PARAM,
	ARGUMENTS_NOT_RELATIVELY_PRIME,
	//EXPR_BREAK,
	//EXPR_OUT_OF_MEMORY,
	//EXPR_CANNOT_USE_X_IN_EXPONENT,
	//EXPR_DEGREE_TOO_HIGH,
	EXPONENT_TOO_LARGE,
	EXPONENT_NEGATIVE,
	//EXPR_LEADING_COFF_MULTIPLE_OF_PRIME,
	//EXPR_CANNOT_LIFT,
	//EXPR_MODULUS_MUST_BE_GREATER_THAN_ONE,
	//EXPR_MODULUS_MUST_BE_PRIME_EXP,
	EXPR_BASE_MUST_BE_POSITIVE,
	//EXPR_POWER_MUST_BE_POSITIVE,
	EXPR_MODULUS_MUST_BE_NONNEGATIVE,
	//EXPR_VAR_OR_COUNTER_REQUIRED,
	EXPR_OK = 0,
	EXPR_FAIL = 1
};

long long MulPrToLong(const Znum& x);
void writeIni(void);
void generatePrimes(unsigned long long int max_val);
unsigned long long llSqrt(const unsigned long long n);
bool getBit(const unsigned long long int x, const bool array[]);
void textError(retCode rc);
/* remove initial & trailing spaces, tabs, etc from msg (\t, \r, \n, \v   and \f) */
void removeInitTrail(std::string& msg);
/* remove any spaces between 2 digits, also multiple consecutive spaces reduced to 1 space */
void removeIntSpace(std::string& msg);
void strToUpper(const std::string& s, std::string& d);
void PrintTimeUsed(double elapsed, const std::string& msg = "");
void ErrorDisp(const char* lpszFunction);
void ShowLargeNumber(const Znum& Bi_Nbr, int digitsInGroup, bool size, bool hexPrFlag);
/* ComputeNumDigits(n,r): Number of digits of n in base r. */
long long ComputeNumDigits(const Znum& n, const Znum& radix);
bool factortest(const Znum& x3, const int testnum, const int method = 0);
void printSummary(void);


extern bool* primeFlags;
extern unsigned long long* primeList;
extern unsigned int prime_list_count;
extern unsigned long long int primeListMax;
extern int verbose;
extern long long lltTdivCnt;  /* count no of Mersenne numbers (partly) factored by
                           trial division */
extern long long lltCmpCnt;   /* count no of Mersenne numbers determined to be composite
                           using Lucas-Lehmer test */
extern long long lltPriCnt;   /* count no of Mersenne numbers determined to be prime
                           using Lucas-Lehmer test */
extern int lang;    // 0 English, 1 = Spanish
extern int groupSize;
extern bool hexPrFlag;
extern HWND handConsole;      /* handle to console window */
extern std::string helpFilePath;

/* defined in various places: not in main.cpp */

void doTests14(const std::vector<std::string>& p);

/* yafu.cpp */
bool callYafu(const Znum& num, fList& Factors);
void yafuParam(const std::vector<std::string>& p);      /*process YAFU commands */
extern std::string YafuPath;
extern std::string yafuprog;
extern std::string outPath;
char* getFileName(const char* filter, HWND owner, bool MustExist = true);
bool fileStatus(const std::string& fileName);
bool changepathPP(std::string& path, std::string& prog);
bool changepath2(std::string& path);
extern bool yafu;

/* msieve.cpp */
bool callMsieve(const Znum& num, fList& Factors);
void msieveParam(const std::vector<std::string>& p);   /*process Msieve commands */
extern std::string MsievePathS;
extern std::string MsieveProg;
extern std::string MsieveLogPath;
extern bool msieve;

/* evalexpr.cpp */
retCode ComputeExpr(const std::string& expr, Znum& Result, int& asgCt, bool* multiV = nullptr);
/* evaluate 1 or more expressions, separated by commas */
retCode ComputeMultiExpr(std::string expr, Znum result);
void printvars(std::string name);
extern std::vector <Znum> roots;   /* used by functions that return multiple values */
Znum llt(const Znum& p);
Znum R4(Znum num);

/* libpariinterface.cpp */
void pariParam(const std::vector<std::string>& p);
void parifactor(const Znum& n, fList& factors);
extern bool Pari;
Znum R3h(Znum n);
Znum Hclassno12(const Znum& n);
Znum classno(const Znum& n, int flag);
Znum tau(const Znum& n);
Znum stirling(const Znum& n, const Znum& m, const Znum& flag);
Znum quaddisc(const Znum& n);
Znum eulerfrac(const Znum& n);

/* quadmod.cpp */
void doTestsA(const std::vector<std::string>& params);   /* quadratic modular equation solver */
int quadModEqn(const std::vector<std::string>& p);  /* Quadratic Modular Equation Solver */
std::vector <Znum> ModSqrtQE(const Znum& aa, const Znum& m); /* modular square root*/

/* Modsqrt.cpp */
std::vector <Znum> primeModSqrt(const Znum& aa, const Znum& p);
std::vector <Znum> ModSqrt(const Znum& aa, const Znum& m);
void doTests9(const std::vector<std::string> & p);  /* modular square root test */

/* fileversioninfo.cpp */
void VersionInfo(const LPCSTR path, int ver[4], std::string& modified);
DWORD getComCtlVer(void);

/* mpz_prp.cpp */
/*************************************************************/
/*************************************************************/
/* These are the definitions for the probable prime routines */
/*************************************************************/
/*************************************************************/
#define PRP_ERROR -1
#define PRP_COMPOSITE 0 /* composite, not a pseudoprime */
#define PRP_SPSP 1     /* strong pseudoprime */
#define PRP_WPSP 2      /* weak pseudo-prime*/
#define PRP_PRIME 3     /* definate prime */
#define PRP_PRP 4       /* probable prime */

#define APRTCLE_ERROR -1
#define APRTCLE_COMPOSITE 0
#define APRTCLE_PRIME 3       /* definate prime */

int mpz_bpsw_prp(const mpz_t n); /* Baillie-Pomerance-Selfridge-Wagstaff probablistic primality test*/
int mpz_aprtcle(const mpz_t N, const int verbose);  /* APR-CL prime testing */

/* R3.cpp */
bool isPrime2(unsigned __int64 num);
size_t inverseTotient(__int64 n, std::vector<unsigned __int64>** result, bool debug,
	int level, bool dump);
long long R4alt(int num);
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
Znum             modMult(const Znum& a, const Znum& b, const Znum& mod);
extern unsigned __int64 R2(const unsigned __int64 n);
extern unsigned __int64 R3(__int64 n);
/* get prime factors of tnum, using trial division */
unsigned int primeFactors(unsigned __int64 tnum, factorsS& f);

/* help.cpp */
void helpfunc(const std::vector<std::string>& command);

/* partition.cpp */
void biperm(int n, Znum& result);

/* computerinfo.cpp */
void get_computer_info(char* CPUidstr, double &MEAS_CPU_FREQUENCY);
