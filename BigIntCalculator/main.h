#pragma once


/* defined in main.cpp */
/* error and return codes, errors are -ve, OK is 0, FAIL is +1 */
enum class retCode
{
	NUMBER_TOO_LOW = -100,
	NUMBER_TOO_HIGH,
	INTERM_TOO_HIGH,
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
void PrintTimeUsed(double elapsed, const std::string& msg = "");
void ErrorDisp(const char* lpszFunction);
void ShowLargeNumber(const Znum& Bi_Nbr, int digitsInGroup, bool size, bool hexPrFlag);
/* ComputeNumDigits(n,r): Number of digits of n in base r. */
long long ComputeNumDigits(const Znum& n, const Znum& radix);


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

extern std::string YafuPath;
extern std::string yafuprog;
extern std::string outPath;
extern std::string MsievePathS;
extern std::string MsieveProg;
extern std::string MsieveLogPath;
extern bool breakSignal;
extern std::vector <Znum> roots;   /* used by functions that return multiple values */
extern bool msieve;
extern bool yafu;
extern bool Pari;

const char* myTime(void);
//void delfile(const std::string& path, const char* FileName);
void msieveParam(const std::vector<std::string>& p);   /*process Msieve commands */
void yafuParam(const std::vector<std::string>& p);      /*process YAFU commands */
void pariParam(const std::vector<std::string>& p);      /*process YAFU commands */
int quadModEqn(const std::vector<std::string>& p);  /* Quadratic Modular Equation Solver */
std::vector <Znum> primeModSqrt(const Znum& aa, const Znum& p);
std::vector <Znum> ModSqrtQE(const Znum& aa, const Znum& m); /* modular square root*/
std::vector <Znum> ModSqrt(const Znum& aa, const Znum& m);
void printvars(std::string name);
void doTests9(const std::vector<std::string> & p);  /* modular square root test */
void doTestsA(const std::vector<std::string> & params);   /* quadratic modular equation solver */
void VersionInfo(const LPCSTR path, int ver[4], std::string& modified);
char* getFileName(const char* filter, HWND owner, bool MustExist = true);
/* check file status. Print date & time modified, return false if file not found */
bool fileStatus(const std::string& fileName);
bool changepathPP(std::string& path, std::string& prog);
bool changepath2(std::string& path);
DWORD getComCtlVer(void);
retCode ComputeExpr(const std::string& expr, Znum& Result, int& asgCt, bool* multiV = nullptr);
/* evaluate 1 or more expressions, separated by commas */
retCode ComputeMultiExpr(std::string expr, Znum result);
int mpz_bpsw_prp(const mpz_t n); /* Baillie-Pomerance-Selfridge-Wagstaff probablistic primality test*/
int mpz_aprtcle(const mpz_t N, const int verbose);  /* APR-CL prime testing */
bool isPrime2(unsigned __int64 num);
void helpfunc(const std::vector<std::string>& command);
void biperm(int n, Znum& result);
size_t inverseTotient(__int64 n, std::vector<unsigned __int64>** result, bool debug,
	int level, bool dump);

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
Znum             modMult(const Znum& a, const Znum& b, const Znum &mod);

Znum R3h(Znum n);
Znum Hclassno12(const Znum& n);
Znum classno(const Znum& n, int flag);
Znum tau(const Znum& n);
Znum stirling(const Znum& n, const Znum& m, const Znum& flag);
void get_computer_info(char* CPUidstr, double &MEAS_CPU_FREQUENCY);
Znum llt(const Znum& p);