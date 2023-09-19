#pragma once


/* defined in main.cpp */
long long MulPrToLong(const Znum& x);


void writeIni(void);
void generatePrimes(unsigned long long int max_val);
unsigned long long llSqrt(const unsigned long long n);
bool getBit(const unsigned long long int x, const bool array[]);
void textError(retCode rc);



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

/* defined in various places: not in main.cpp */
extern long long SFhitcount, SFmisscount;
extern std::string YafuPath;
extern std::string yafuprog;
extern std::string outPath;
extern std::string MsievePath;
extern std::string MsieveProg;
extern std::string MsieveLogPath;
extern bool breakSignal;
extern std::vector <Znum> roots;   /* used by functions that return multiple values */
extern bool msieve;
extern bool yafu;
extern bool Pari;
extern bool hexPrFlag;

const char* myTime(void);
void delfile(const std::string& path, const char* FileName);
bool isPerfectSquare(const Znum &num);
void msieveParam(const std::vector<std::string>& p);   /*process Msieve commands */
void yafuParam(const std::vector<std::string>& p);      /*process YAFU commands */
void pariParam(const std::vector<std::string>& p);      /*process YAFU commands */
int quadModEqn(const std::vector<std::string>& p);  /* Quadratic Modular Equation Solver */
void printvars(std::string name);
void doTests9(const std::vector<std::string> & p);  /* modular square root test */
void doTestsA(const std::vector<std::string> & params);   /* quadratic modular equation solver */
void PrintTimeUsed(double elapsed, const std::string& msg = "");
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

/* remove initial & trailing spaces, tabs, etc from msg (\t, \r, \n, \v   and \f) */
void removeInitTrail(std::string& msg);
/* remove any spaces between 2 digits, also multiple consecutive spaces reduced to 1 space */
void removeIntSpace(std::string& msg);

void helpfunc(const std::vector<std::string>& command);


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
void get_computer_info(char* CPUidstr, double MEAS_CPU_FREQUENCY);
Znum llt(const Znum& p);