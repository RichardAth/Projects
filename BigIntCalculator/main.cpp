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
#pragma fenv_access (on)  /* allow access to floating point environment */

#include "pch.h"
#include <strsafe.h>

#include <Mmsystem.h >   // for sound effects
#include "diagnostic.h"
#include "resource.h"
#include <commctrl.h>  /* requires Comctl32.lib. In visual studio:
 project -> properties -> linker -> input -> additional dependencies*/
/* use manifest to ensure that latest Comctl32.dll is used, otherwise the old 
dll for Windows 7 would be used */
#pragma comment(linker,"\"/manifestdependency:type='win32' \
name='Microsoft.Windows.Common-Controls' version='6.0.0.0' \
processorArchitecture='*' publicKeyToken='6595b64144ccf1df' language='*'\"")

#define BIGNBR       // define to include bignbr tests 
#ifdef BIGNBR
extern Znum zR, zR2, zNI, zN;
#endif

EXTERN_C IMAGE_DOS_HEADER __ImageBase;
#define HINST_THISCOMPONENT ((HINSTANCE)&__ImageBase)
HINSTANCE calcHandle = HINST_THISCOMPONENT;  /* instance handle for this program */
HINSTANCE cH2 = GetModuleHandle(NULL);  /* should contain the same value */

int lang = 0;             // 0 English, 1 = Spanish
bool VSrun = false;       /* set to true if program started from Visual Studio */

bool hexPrFlag = false;		// set true if output is in hex
int factorFlag = 2;     /* 0 = no factorisation, 1 = factorisation but not totient etc,
                           2 = get totient, number of divisors etc after factorisation */
int groupSize = 6;      /* for large numbers, insert space between 6 digit groups */
/* verbose value is used to turn off or on optional messages; 
higher value = more messages */
#ifdef _DEBUG
int verbose = 1;       /* default 1 if compiled in debug mode */
#else
int verbose = 0;
#endif

HANDLE hConsole;       /* used by SetConsoleCursorPosition() function */
HWND handConsole;      /* handle to console window */

bool *primeFlags = NULL;
unsigned long long *primeList = NULL;
unsigned int prime_list_count = 0;
unsigned long long int primeListMax = 0;
// machine info
double MEAS_CPU_FREQUENCY;
char CPU_ID_STR[80] = { '\0' };

std::vector <std::string> exprList;  /* expressions stored here as text once
                                       they are validated */

/* name of sound file played at end of processing a command or expression
 if the elapsed time > 10 seconds */
std::string endsound = "c:/Windows/Media/Alarm09.wav";

/* name of sound file played when prompt for input is displayed */
std::string attsound = "c:/Windows/Media/chimes.wav";


/* display last error in a message box. lpszFunction should be a pointer to text 
containing the function name. 
See https://learn.microsoft.com/en-us/windows/win32/debug/retrieving-the-last-error-code */
void ErrorDisp(const char *lpszFunction)
{
    // Retrieve the system error message for the last-error code

    LPVOID lpMsgBuf;
    LPVOID lpDisplayBuf;
    DWORD dw = GetLastError(); /* GetLastError gets the last error that was set 
                                  by a Windows API function */
    if (dw == 0)
        return;   /* no error so just return */

    FormatMessage(
        FORMAT_MESSAGE_ALLOCATE_BUFFER |
        FORMAT_MESSAGE_FROM_SYSTEM |
        FORMAT_MESSAGE_IGNORE_INSERTS,
        NULL,
        dw,
        MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
        (LPTSTR)&lpMsgBuf,
        0, NULL);

    // Display the error message

    lpDisplayBuf = (LPVOID)LocalAlloc(LMEM_ZEROINIT,
        (lstrlen((LPCTSTR)lpMsgBuf) + strlen(lpszFunction) + 40) * sizeof(TCHAR));
    StringCchPrintf((LPTSTR)lpDisplayBuf,
        LocalSize(lpDisplayBuf) / sizeof(TCHAR),
        TEXT("%S failed with error %d: %s"),
        lpszFunction, dw, lpMsgBuf);
    MessageBox(NULL, (LPCTSTR)lpDisplayBuf, TEXT("Error"), MB_OK);

    LocalFree(lpMsgBuf);
    LocalFree(lpDisplayBuf);
    // ExitProcess(dw);
}


/* get time in format hh:mm:ss */
const char * myTime(void) {
    static char timestamp[10];   // time in format hh:mm:ss
    struct tm newtime;

    const time_t current = time(NULL);  // time as seconds elapsed since midnight, January 1, 1970
    localtime_s(&newtime, &current);    // convert time to tm structure
    /* convert time to hh:mm:ss */
    strftime(timestamp, sizeof(timestamp), "%H:%M:%S", &newtime);
    return timestamp;
}


/* Convert number to hexdecimal. Ensure that if number is negative the leftmost
bit of the most significant digit is set, and conversely, if the number is positive
the leftmost bit of the most significant digit is not set. This is done by 
prefixing the output with '0' or 'f' when necessary. */
static char* getHex(Znum Bi_Nbr) {
    static char *hexbuffer = NULL;  // IMPORTANT. This must be a static variable!
    std::string obuff;
    if (Bi_Nbr >= 0) {
        hexbuffer = mpz_get_str(NULL, 16, ZT(Bi_Nbr));
        if (!(hexbuffer[0] >= '0' && hexbuffer[0] <= '7')) {
            /* we need to insert a zero. move everything down 1 byte*/
            hexbuffer = (char *)realloc(hexbuffer, strlen(hexbuffer) + 2);
            assert(hexbuffer != NULL);
            memmove(&hexbuffer[1], &hexbuffer[0], strlen(hexbuffer) + 1); 
            hexbuffer[0] = '0';
        }
    }
    else {
        if (Bi_Nbr == -1)
            obuff = 'f';		// -1 has to be a special case		
        else while (Bi_Nbr != -1) {
            Znum q, r;
            long long rll;
            mpz_fdiv_q_2exp(ZT(q), ZT(Bi_Nbr), 4);	// q = Bi_Nbr/16 rounded towards -inf
            mpz_fdiv_r_2exp(ZT(r), ZT(Bi_Nbr), 4);	// r = Bi_Nbr - q*16 (0 <= r <= 15)
            rll = MulPrToLong(r);
            obuff += (rll <= 9) ? '0' + (char)rll : 'a' + (char)rll -10;  // get char '0' to 'f'
            Bi_Nbr = q;
        }

        std::reverse(std::begin(obuff), std::end(obuff));  // get digits into correct order
        obuff += '\0';				// null-terminate string
        if (obuff[0] >= '0' && obuff[0] <= '7') {
            /* we need to insert an f */
            obuff.insert(obuff.begin(), 'f');
        }
        /* copy text from local STL string to global C-style string */
        hexbuffer = (char *)malloc(obuff.size() + 1);
        assert(hexbuffer != NULL);
        strncpy_s(hexbuffer, obuff.size(), obuff.c_str(), _TRUNCATE);
    }
    return hexbuffer;  // hexbuffer must be a static variable
}

/* output value of Nbr as ascii text to stdout. DigitsInGroup controls the spacing
between groups of digits. If size = true and the number is large show the number 
of digits. If hexPrFlag is true output in hex instead of decimal. */
void ShowLargeNumber(const Znum &Bi_Nbr, int digitsInGroup, bool size, bool hexPrFlag) {
    std::string nbrOutput = "";
    char* buffer = NULL;
    size_t msglen, index = 0;

    if (!hexPrFlag) {
        // convert to null-terminated ascii string, base 10, with sign if -ve
        // corrrectly-sized buffer is allocated automatically.
        buffer = mpz_get_str(NULL, 10, ZT(Bi_Nbr));
    }
    else {
        buffer = getHex(Bi_Nbr);
    }
    msglen = strnlen(buffer, 500000);  // arbitrary limit of 500000
    if (buffer[0] == '-') {      // if number is -ve
        nbrOutput = "-";
        index = 1;
    }
    if (hexPrFlag) 
        nbrOutput += "0x";
    for (; index < msglen; index++) {
        if (digitsInGroup > 0 && (msglen - index) % digitsInGroup == 0) {
            // divide digits into groups, if there are more than 6 digits
            if ((index > 0 && buffer[0] != '-') || (index > 1)) {
                nbrOutput += " ";  // put space after group of digits 
            }
        }
        nbrOutput += buffer[index];
    }
    if (buffer[0] == '-')   // if number is -ve
        msglen--;
    free(buffer);		// avoid memory leakage
    std::cout << nbrOutput;
    if (msglen > 24 && size)
        if (lang) std::cout << " (" << msglen << " dígitos)";
        else std::cout << " (" << msglen << " digits)";
}

/* convert biginteger to normal. Checks for overflow */
long long MulPrToLong(const Znum &x) {
    long long rv;
    // note: do not use mpz_fits_slong_p because it checks whether x fits a 32 bit integer, rather than a 64 bit integer.
    if (x >= LLONG_MIN && x <= LLONG_MAX) { // is x value OK for normal integer?
        //rv = mpz_get_si(ZT(x)); // convert to normal integer
        rv = ZT(x)->_mp_d[0];     // accessing the limb directly seems to be a lot faster
        if (ZT(x)->_mp_size < 0)  // than calling mpz_get_si
            rv = -rv;
        if (ZT(x)->_mp_size == 0)
            rv = 0;
        return rv;
    }
    else {
        StackTrace2();
        throw std::range_error("big number cannot be converted to 64-bit integer");
    }
    return 0;
}

/* NumDigits(n,r): Number of digits of n in base r. Leading zeros are not counted  */
long long ComputeNumDigits(const Znum &n, const Znum &radix)
{
    /* more accurate than mpz_sizeinbase */
    Znum result = n;
    long long digits = 0;
    
    while (result > 0) {
        result /= radix;
        digits++;
    }
    return digits;
}

/* get 'width' of n in bits */
static long long NoOfBits(const Znum &n) {
    auto result = mpz_sizeinbase(ZT(n), 2);  // calculate number of bits
    return result;
}


long long lltTdivCnt = 0;  /* count no of Mersenne numbers (partly) factored by
                           trial division */
long long lltCmpCnt = 0;   /* count no of Mersenne numbers determined to be composite
                           using Lucas-Lehmer test */
long long lltPriCnt = 0;   /* count no of Mersenne numbers determined to be prime
                           using Lucas-Lehmer test */


/* get floor(sqrt(n))*/
unsigned long long llSqrt(const unsigned long long n) {
    double root = sqrt((double)n);
    unsigned long long  iroot = llround(root);
    while (iroot*iroot > n)
        iroot--;
    return iroot;
}

// return value of bit in array corresponding to x. false (0) indicates a prime number
bool getBit(const unsigned long long int x, const bool array[])
{
    return array[x/2];
}

/* sets bit in 'array' corresponding to 'x' */
/* assume 'num' is odd. no bits are stored for even numbers */
static void setBit(const unsigned long long int x, bool array[]) {
    array[x / 2] = true;
    return;
}

void generatePrimes(unsigned long long int max_val) {
    unsigned long long int numsave=0, count = 1, num = 3;
    size_t plist_size = 0;  // size of prime list
    unsigned long long int sqrt_max_val;

    if (primeListMax >= max_val)
        return;  // nothing to do if prime list already generated

    max_val += 63 - (max_val + 63) % 64;  // round up to next multiple of 64
    sqrt_max_val = llSqrt(max_val);

    /* initialise flags */
    
    if (primeFlags != NULL) delete []primeFlags;	// if generatePrimes has been called before
                                                // clear out old prime list
    //primeFlags = (uint32_t*)calloc((max_val / 16) + 1, 1);
    primeFlags = new bool[max_val / 2 + 1]{ 0 };
    assert(primeFlags != NULL);
    
    // allocate storage for primeList if required
    {
        if (verbose > 0)
            if (lang)
                printf_s("número esperado de primos es %.0f \n", 
                    (double)max_val / (log((double)max_val) - 1));
            else
                printf_s("Expected no of primes is %.0f\n",
                (double)max_val / (log((double)max_val) - 1));
        if (primeList != NULL) free(primeList);
        plist_size = (size_t)((double)max_val / (log((double)max_val) - 1)) * 102 / 100;
        // add 2% for safety
        primeList = (unsigned long long *)malloc(plist_size * sizeof(long long));
        assert(primeList != NULL);
        prime_list_count = 1;
        primeList[0] = 2;  // put 1st prime into list
    }

    while (num <= max_val)
    {
        if (!getBit(num, primeFlags))
        {   /* we have found a prime number */
            primeList[count] = num;
            count++;
            numsave = num;

            if (num <= sqrt_max_val)  /* check whether starting value for i is already
                                      more than max_val, while avoiding an integer overflow*/
                for (unsigned long long int i = num*num; i <= max_val; i += (num * 2))
                { /* mark all odd multiples of num in range as not prime */
                  /* starting value is num*num rather than num*3 because the flag bits for
                  multiples of num less than num*num are already set and using num*num
                  makes the program run significantly faster */
                    setBit(i, primeFlags);
                }
        }
        num += 2;	// advance to next possible prime
    }

    // after completing the while loop we have found all the primes < max_val
    if (verbose > 0) {
        if (lang)
            printf_s("  el primo %9lld es %11lld\n", count, numsave);
        else 
            printf_s("  prime %9lld is %11lld\n", count, numsave);
    }
    primeList[count] = ULLONG_MAX;		// set end marker
    prime_list_count = (unsigned int)count;
    primeListMax = primeList[count - 1];
    return;
}


/* translate error code to text and output it*/
void textError(retCode rc) {
    /*
    error codes currently used:
    NUMBER_TOO_LOW,
    NUMBER_TOO_HIGH,
    INTERM_TOO_HIGH,		
    DIVIDE_BY_ZERO,        
    PAREN_MISMATCH,		
    SYNTAX_ERROR,			
    TOO_MANY_PAREN,		
    INVALID_PARAM,	
    ARGUMENTS_NOT_RELATIVELY_PRIME
    EXPONENT_TOO_LARGE,
    EXPONENT_NEGATIVE,				
    EXPR_OK = 0
    */
    switch (rc)
    {
    case retCode::NUMBER_TOO_LOW:
        std::cout << (lang ? "Número muy pequeño\n" : "Number too low\n");
        break;
    case retCode::NUMBER_TOO_HIGH:
        std::cout << (lang ? "Número muy grande \n" :
            "Number too high \n");
        break;
    case retCode::INTERM_TOO_HIGH:
        std::cout << (lang ? "Número intermedio muy grande (más de 20000 dígitos\n" :
            "Intermediate number too high (more than 20000 digits)\n");
        break;
    case retCode::DIVIDE_BY_ZERO:
        std::cout << (lang ? "División por cero\n" : "Division by zero\n");
        break;
    case retCode::PAREN_MISMATCH:
        std::cout << (lang ? "Error de paréntesis\n" : "Parenthesis mismatch\n");
        break;
    case retCode::SYNTAX_ERROR:
        if (lang) 	{
            std::cout << ( "Error de sintaxis\n");
        }
        else {
            std::cout << ("Syntax error\n");
        }
        break;
    case retCode::TOO_MANY_PAREN:
        std::cout << (lang ? "Demasiados paréntesis\n" : "Too many parenthesis\n");
        break;
    case retCode::INVALID_PARAM:
        std::cout << (lang ? "Parámetro inválido\n" : "Invalid parameter\n");
        break;
    case retCode::ARGUMENTS_NOT_RELATIVELY_PRIME:
        std::cout << (lang ? "MCD de los argumentos no es 1\n" : "GCD of arguments is not 1\n");
        break;
    /*case EXPR_BREAK:
        std::cout << (lang ? "Detenido por el usuario\n" : "Stopped by use\nr");
        break;*/
    case retCode::EXPONENT_NEGATIVE: {
        std::cout << (lang? "Exponente no debe ser negativo\n" : "Exponent must not be negative\n");
        break;
    }
    case retCode::EXPONENT_TOO_LARGE: {
        std::cout << (lang? "El exponente es mayor que 2^31-1\n": "Exponent exceeds 2^31-1\n");
        break;
    }
    /*case retCode::EXPR_VAR_OR_COUNTER_REQUIRED:
        if (lang)
        {
            std::cout <<  "La expresión \n";
        }
        else
        {
            std::cout << "Expression #\n";
        }
        break;*/
    case retCode::EXPR_BASE_MUST_BE_POSITIVE:
        std::cout << (lang ? "La base debe ser mayor que un\n" :
            "Base must be greater than one\n");
        break;
    //case retCode::EXPR_POWER_MUST_BE_POSITIVE:
    //	std::cout << (lang ? "La potencia debe ser mayor que cero\n" :
    //		"Power must be greater than zero\n");
    //	break;
    //case retCode::EXPR_MODULUS_MUST_BE_GREATER_THAN_ONE:
    //	std::cout << (lang ? "El módulo debe ser mayor que 1\n" : "Modulus must be greater than one\n");
    //	break;
    case retCode::EXPR_MODULUS_MUST_BE_NONNEGATIVE:
        std::cout << (lang ? "El módulo no debe ser negativo\n" :
            "Modulus must not be negative\n");
        break;
    default:
        printf_s( "unknown error code: %d\n", (int)rc);
        break;
    }
}

/* convert s to UPPER CASE in d.  d could be the same string as s 
similar to _strupr in c */
static void strToUpper(const std::string &s, std::string &d) {
    if (&d != &s)
        d.resize(s.size());  // resize d unless it's the same string as s
    for (size_t ix = 0; ix < s.size(); ix++)
        d[ix] = toupper(s[ix]);
}

/* print elapsed time. If > 60 seconds print in hour min sec format */
void PrintTimeUsed(double elapsed, const std::string &msg) {

    if (msg.size() > 1)
        std::cout << myTime() << ' ' << msg;
    auto elSec = elapsed / CLOCKS_PER_SEC; // convert ticks to seconds

    if (elSec > 10.0)
        PlaySoundA(endsound.c_str(), NULL,
            SND_FILENAME | SND_NODEFAULT | SND_ASYNC | SND_NOSTOP);

    if (elSec <= 60.0) {
        if (lang)
            printf_s("%.4f segundos \n", elSec);
        else
            printf_s("%.4f seconds \n", elSec);  /* print time used to nearest millisecond */
    }
    else {
        /* round to nearest second */
        long long sec = (long long)std::floor(elSec); // convert to an integer
        long long min = sec / 60;  // min >= 1
        elSec -= min * 60;         // get seconds
        long long hour = min / 60; // hour may be zero
        min %= 60;                 // get minutes
        if (hour > 0)
            std::cout << ' ' << hour << "h " << min << " min " << elSec << "sec \n";
        else 
            std::cout << ' ' << min << " min " << elSec << " sec \n";
    }
}

/* remove spaces, tabs, etc  from msg (\t, \r, \n, \v   and \f)*/
static void removeBlanks(std::string &msg) {
    for (size_t ix = 0; ix < msg.size(); ix++) {
        if ((unsigned char)msg[ix] <= 0x7f && isspace(msg[ix])) {     // look for spaces, tabs, etc
            msg.erase(ix, 1);      // remove space character
            ix--;  // adjust index to take account of removed blank
        }
    }
}

/* remove initial & trailing spaces, tabs, etc from msg (\t, \r, \n, \v   and \f) */
void removeInitTrail(std::string &msg) {
    while (!msg.empty() && isspace(msg.front())) {
        msg.erase(0, 1);      // remove 1st space character
    }
    while (!msg.empty() && isspace(msg.back())) {
        msg.resize(msg.size() - 1);  /* remove trailing spaces */
    }
}

/* remove any spaces between 2 digits, also multiple consecutive spaces reduced to 1 space */
void removeIntSpace(std::string& msg) {
    if (msg.empty())
        return;
    for (int i = 1; i < (msg.size() - 1);) {
        if (isdigit(msg[i - 1]) && isspace(msg[i]) && 
            (isdigit(msg[i + 1]) || isspace(msg[i+1]))) {
            msg.erase(i, 1);  /* remove space between 2 digits*/
        }
        else
            i++;  /* advance to next char if no space */
    }
}

struct summary {
    int numsize;		// number of decimal digits in number
    double time;		// time used in seconds
    int NumFacs;		// number of unique factors
    int totalFacs;		// total number of factors (=1 for a prime number)
    int testNum;		// test number (if applicable)
    int sndFac;		    // number of decimal digits in 2nd largest factor
    struct counters ctrs;
};

std::vector <summary> results;

/* print summary - 1 line per test */
static void printSummary(void) {
    long long sec, min, hour;
    double elSec;
    /* print column headings */
    printf_s("Test Num Size   time      Unique Factors Total Factors     2nd Fac");
    printf_s(" tdv prh leh crm pm1 ecm siq pwr yaf msv \n");
    for (auto res : results) {
        /* truncate elapsed time to nearest second */
        sec = (long long)std::floor(res.time); // convert to an integer
        min = sec / 60;  // min may be 0
        elSec = res.time - min * 60;         // get seconds including fractional part
        hour = min / 60; // hour may be zero
        min %= 60;       // get minutes
        printf_s("%4d %8d %2lld:%02lld:", res.testNum, res.numsize,
            hour, min);
        if (elSec < 10.0)
            printf_s("0%.2f ", elSec);  // insert leading zero
        else
            printf_s("%.2f ", elSec);

        /* print counters showing how factors were found */
        printf_s("     %8d   %8d     %8d", res.NumFacs, res.totalFacs, res.sndFac);
        printf_s("    %3d %3d %3d %3d %3d %3d ",
            res.ctrs.tdiv, res.ctrs.prho, res.ctrs.leh, res.ctrs.carm,
            res.ctrs.pm1, res.ctrs.ecm);
        printf_s("%3d %3d %3d %3d \n", res.ctrs.siqs, res.ctrs.power, res.ctrs.yafu,
            res.ctrs.msieve);
    }

    if (verbose > 0) {
        std::cout << "saved Factors list hit  count = " << SFhitcount << '\n'
                  << "saved Factors list miss count = " << SFmisscount << '\n';;
    }
}

/* factorise Result, calculate number of divisors etc and print results */
static void doFactors(const Znum &Result, bool test) {
    fList factorlist;
    Znum Quad[4];
    clock_t start;
    summary sum;    // save summary of factorisation
    bool rv;

    if (test)
        start = clock();	// used to measure execution time

    /* call DA´s magic function to factorise Result */
    if (factorFlag > 1)
        rv = factorise(Result, factorlist, Quad);
    else
        rv = factorise(Result, factorlist, nullptr);
    if (rv && factorlist.fsize() > 0) {
        factorlist.print(Result < 0);   /* print factors */
        std::cout << '\n';
        if (factorFlag > 1) {
            auto divisors = factorlist.NoOfDivs();
            std::cout << (lang? "Cantidad de Divisores = " : "Number of Divisors = ");
            ShowLargeNumber(divisors, groupSize, false, false);

            divisors = factorlist.DivisorSum();
            std::cout << (lang? "\nSuma de divisores     = " : "\nSum of Divisors    = ");
            ShowLargeNumber(divisors, groupSize, false, false);
            divisors = factorlist.totient();
            std::cout << (lang ? "\nPhi de Euler          = " : "\nTotient            = ");
            ShowLargeNumber(divisors, groupSize, false, false);
            if (Result > 0) {
                auto mob = factorlist.mob();  // mobius only defined for +ve integers
                std::cout << "\nMöbius             = " << mob;
    
            }

            /* show that the number is the sum of 4 or fewer squares. See
            https://www.alpertron.com.ar/4SQUARES.HTM */
            if (Result >= 0)
                std::cout << "\nn = ";
            else
                std::cout << "\n-n = ";
            char c = 'a';
            for (auto q : Quad) {  // print "n = a² + b² + c² + d²"
                if (q == 0)
                    break;
                if (c > 'a') std::cout << "+ ";  // precede number with + unless its the 1st number
                std::cout << c << "² ";
                c++;    // change a to b, b to c, etc
            }
            c = 'a';
            for (auto q : Quad) {
                if (q == 0)
                    break;
                std::cout << "\n" << c << "= " << q;
                c++;  // change a to b, b to c, etc
            }
            std::cout << "\n";
            factorlist.prCounts();  // print counts
        }
        
        if (test) {
            /* recalculate result & get total number of factors */
            Znum result = factorlist.recheck(sum.totalFacs);
            if (result != Result) {
                std::cout << "Factors expected value " << Result << " actual value " << result << '\n';
                Beep(750, 1000);
            }
            if (factorFlag > 1) {
                result = Quad[0] * Quad[0] + Quad[1] * Quad[1] + Quad[2] * Quad[2] + Quad[3] * Quad[3];
                if (result != Result) {
                    std::cout << "Quad expected value " << Result << " actual value " << result << '\n';
                    Beep(750, 1000);
                }
                else {
                    if (factorlist.twosq()) {
                        if (Quad[2] != 0 || Quad[3] != 0)
                            std::cout << "expected c = d= 0; got: \n"
                            << "number = " << Result << '\n'
                            << "a= " << Quad[0] << '\n'
                            << "b= " << Quad[1] << '\n'
                            << "c= " << Quad[2] << '\n'
                            << "d= " << Quad[3] << '\n';
                    }
                    //else if (factorlist.sqplustwosq()) {
                    //	if (Quad[1] != Quad[2] || Quad[3] != 0)
                    //		std::cout << "expected b = c and d= 0; got: \n"
                    //		<< "number = " << Result << '\n'
                    //		<< "a= " << Quad[0] << '\n'
                    //		<< "b= " << Quad[1] << '\n'
                    //		<< "c= " << Quad[2] << '\n'
                    //		<< "d= " << Quad[3] << '\n';
                    //}
                }
            }
            auto end = clock(); 
            double elapsed = (double)end - start;
            PrintTimeUsed(elapsed, lang? "Tiempo transcurrido = " :"time used = ");

            /* store info for summary */
            sum.time = elapsed / CLOCKS_PER_SEC;
            sum.numsize = (int)ComputeNumDigits(result, 10);
            sum.NumFacs = (int)factorlist.fsize();
            /* get number of digits in 2nd largest factor */
            sum.sndFac = factorlist.sndFac();
            sum.ctrs = factorlist.getCtrs();
            results.push_back(sum);
        }
    }
    else
        std::cout << (lang? "no se puede factorizar\n" : " cannot be factorised\n");
}

/* perform some simple tests. Returns true if x3 is prime 
method = 0 for standard factorisation, != 0 to use only YAFU for factorisation */
static bool factortest(const Znum &x3, const int testnum, const int method=0) {
    fList factorlist;
    Znum Quad[4], result;
    long long totalFactors = 0;
    summary sum;    // save summary of test

    auto start = clock();	// used to measure execution time
    double end, elapsed;

    sum.numsize = (int)ComputeNumDigits(x3, 10);

    if (lang)
        std::cout << '\n' << myTime() << "Prueba " << testnum << ": factoriza ";
    else
        std::cout << '\n' << myTime() << " Test " << testnum << ": factorise ";
    ShowLargeNumber(x3, groupSize, true, false);
    std::cout << '\n';
    if (method == 0) {
        factorise(x3, factorlist, Quad);  // get factors
    }
    if (method != 0) {
        factorlist.set(x3);
        callYafu(x3, factorlist);
    }

    /* recalculate result & get total number of factors */
    result = factorlist.recheck(sum.totalFacs);
    if (result != x3) {
        std::cout << "Factors expected value " << x3 << " actual value " << result << '\n';
        Beep(750, 1000);
    }
    
    /* get number of digits in 2nd largest factor */
    sum.sndFac = factorlist.sndFac();

    if (method == 0) {
        /* check that sum of squares is correct */
        result = Quad[0] * Quad[0] + Quad[1] * Quad[1] + Quad[2] * Quad[2] + Quad[3] * Quad[3];
        if (result != x3) {
            std::cout << "Quad expected value " << x3 << " actual value " << result << '\n'
                << "a= " << Quad[0] << '\n'
                << "b= " << Quad[1] << '\n'
                << "c= " << Quad[2] << '\n'
                << "d= " << Quad[3] << '\n';
            Beep(750, 1000);
        }
        else {
            if (factorlist.twosq()) {
                if (Quad[2] != 0 || Quad[3] != 0)
                    std::cout << "expected c = d= 0; got: \n"
                    << "number = " << x3 << '\n'
                    << "a= " << Quad[0] << '\n'
                    << "b= " << Quad[1] << '\n'
                    << "c= " << Quad[2] << '\n'
                    << "d= " << Quad[3] << '\n';
            }
            else if (factorlist.sqplustwosq()) {
                if ((Quad[0] != Quad[1] && Quad[1] != Quad[2]) || Quad[3] != 0)
                    std::cout << "expected a= b or b = c and d= 0; got: \n"
                    << "number = " << x3 << '\n'
                    << "a= " << Quad[0] << '\n'
                    << "b= " << Quad[1] << '\n'
                    << "c= " << Quad[2] << '\n'
                    << "d= " << Quad[3] << '\n';
            }
            else {
                result = x3;
                while (isEven(result)) result >>= 1;
                if ((numLimbs(result) < 4) && ((result & 7) < 7)) {
                    if (Quad[3] != 0)
                        std::cout << "expected d= 0; got: \n"
                        << "number = " << x3 << '\n'
                        << "a= " << Quad[0] << '\n'
                        << "b= " << Quad[1] << '\n'
                        << "c= " << Quad[2] << '\n'
                        << "d= " << Quad[3] << '\n';
                }
            }
        }
    }

    if (!factorlist.isPrime() ) {
        /* x3 is not prime */

        if (lang)
            std::cout  << factorlist.fsize() << " factores únicos encontrados, total "
            << sum.totalFacs << " factores\n";
        else
        std::cout << "found " << factorlist.fsize() << " unique factors, total "
            << sum.totalFacs << " factors\n";

        if (method == 0)
            factorlist.prCounts();   // print counts
        else
            sum.ctrs.yafu = sum.totalFacs;

        if (lang)
            std::cout << "prueba " << testnum << " terminada as las ";
        else
            std::cout << "test " << testnum << " completed at ";

        end = clock();              // measure amount of time used
        elapsed = (double)end - start;
        PrintTimeUsed(elapsed, lang ? "Tiempo transcurrido = " : "time used = ");
        sum.time = elapsed / CLOCKS_PER_SEC;
        sum.NumFacs = (int)factorlist.fsize();
        sum.testNum = testnum;
        if (method == 0)
            sum.ctrs = factorlist.getCtrs();
        results.push_back(sum);
        return false;   // not prime
    }
    else {
        std::cout << "is prime \n";
        std::cout << "test " << testnum << " completed at ";
        sum.NumFacs = 1;
        sum.totalFacs = 1;
        sum.sndFac = 0;
        sum.testNum = testnum;
        end = clock();              // measure amount of time used
        elapsed = (double)end - start;
        PrintTimeUsed(elapsed, lang ? "Tiempo transcurrido = " : "time used = ");
        sum.time = elapsed / CLOCKS_PER_SEC;
        results.push_back(sum);
        return true;    // is prime
    }
}

static void doTests(void) {
    Znum x3, x4, result, result1;
    int i, asgCt;
    int testcnt = 0;
    struct test {
        std::string text;        // text of expression to be evaluated
        long long expected_result;   // can only do tests that return a value <2^63
    };

    static test testvalues [] 
    {
        "2 - 3 + 4",                        3,
        "2 - (3+4)",                       -5,  // + and - have same priority, left-to-right evaluation
        "-5/2",                            -2,  // uses truncation division
        "-5%2",                            -1,
        "2-5/2",                            0,  // division has higher priority than -
        "(2-5)/2",                         -1,
        "7 * 11 / 2",                      38,  // * and / have same priority, left-to-right evaluation
        "7 * (11/2)",                      35,
        "-2^6",                            64,   // unary - has higher priority than ^
        "0-2^6",                          -64,   // ^ has higher priority than -
        "gcd (12,30)",                      6,
        "modinv (17,21)",                   5,
        "20 + 32^2/4*7",                 1812,    
        "20 + 32^2/(4*7)",                 56,
        "20 - 32^2/4*7",                -1772,
        "(20-32)^2 / 4 * 7",              252,   
        "(20-32)^2 / (4*7)",                5,
        "19!",             121645100408832000,   // factorial
        "19!!",  	                654729075,   // double factorial
        "29!!!",                  72642169600,   // triple factorial
        "33!!!!",                  4996616625,
        "37#",                  7420738134810,   // primorial
        "b(123456789)",	            123456761,   // prime before n
        "f(60)",	            1548008755920,   // fibonacci
        "l(60)",                3461452808002,   // lucas number
        "n(123456789)",             123456791,   // prime after n
        "p(150)",	              40853235313,   // number of partitions
        "modinv(7, 19)",                   11,   // modular inverse
        "7* 11 % 19",                       1,   // verfy modinv in previous test
        "modinv(11,19)",                    7,
        "modpow(8, 7, 6)",                  2,   // modular power
        "modpow(8, -2, 5)",                 4,   // negative exponent OK in some cases
        "totient(201)",		              132,
        "carmichael(201)",                 66,
        "numdivs(7!)",                     60,
        "sumdivs(7!)",                  19344,
        "numdigits(123456789, 6)",         11,
        "sumdigits(123456789, 6)",         19,
        "revdigits(1234567890, 10)",  987654321, 
        "factconcat(2, 11!)", 22222222333355711,
        "le(22, 7)",                        1,    /* legendre */
        "17c7",                         19448,    // binomial coefficient
        "(17!) / ((17-7)!*7!)",         19448,
        "4 ^ 3 ^ 2",                   262144,    // NB expoentiation is right to left evaluation
        "4 ^ (3^2)",                   262144,
        "(4^3) ^ 2",                     4096,
        "4 ^ (3*2)",                     4096,
        "gcd(12, gcd(30, 40))",             2,   // nested function calls
        "7 * (11%19)",                     77,   // override normal left to right evaluation
        "(0x12345678) AND 0x018",        0x18,   // brackets round 1st number are needed so that A of AND is not considered part of number
        "0x12345678 OR 0x1",       0x12345679,
        "0x12345678 XOR 0x10",     0x12345668,
        "(NOT 0x0f23456789abcde1) * -1",  0x0f23456789abcde2,
        "5 < 6 == 7 < 8",                -1,   // returns true (== has lower priority)
        "5 < 6 != 7 < 8",                 0,   // returns false (!= has lower priority)
        "5 < (6 != 7) < 8",              -1,   // returns true; expr evaluated from left to right
        "R3(49)",                        54,
        "R3H(49)",                       54,
        "R2(585)",                       16,
        "SQRT(1234320)",               1110,
        "NROOT(2861381721051424,5)",   1234,
        "LLT(3217)",                      1,  // 2^3217-1 is prime
        "BPSW(2^99-1)",                   0,  // not a prime number
        "BPSW(2^127-1)",                  1,  // a prime number
        "ISPRIME(2^127-1)",              -1,  // a prime number
        "aprcl(2^127-1)",                 2,  // a prime number
        "ispow(2^127-1)",                 0,  /* not a perfect power */
        "ispow(2^127)",                  -1,  /* a perfect power */
        "-not1",                          2,  // operators are processed from right to left
        "not-1",                          0,  // operators are processed from right to left
        "not5#",                        -31,  // # operator evaluated before not
        "4^5#",         1152921504606846976,  // # operator evluated before exponent
        "5!!#",                       30030,  // !! operator evaluated before #
        "5#!!",           42849873690624000,  // # operator evaluated before !!
        "$x = 99 ",                      99,  /* test user variables */
        "$y = $x+1  ",                  100,
        "modsqrt(2191, 23^3)",         1115,
        "modsqrt(4142, 29^3)",         2333,
        "modsqrt(3, 143)",               17,
        "modsqrt(9, 27)",                 3,
        "modsqrt(17, 32)",                7,
        "minfact(99)",                    3,
        "maxfact(99)",                    11,
        "numfact(99)",                     2,
        "lcm(12,20)",                     60,
        "pi(500)",                        95,
        "primroot(761)",                   6,   /* primitive root */
        "hclass(999)",                   384,   /* hurwitz class number*/
        "classno(1000001)",               94,   /* class number */
        "gf(21)",                47297536000,   /* gauss factorial */
        "carmichael(497)",             210,     /* carmichael function */
        "numdivs(116)",                  6,     /* number of divisors*/
        "sumdivs(116)",                210,
        "invtot(132)",                 161,
        "popcnt(123456789)",            16,     /* number of 1-bits */
        "tau(9)",                  -113643,  /* Ramanujan's tau function */
    };

    results.clear();

    auto start = clock();	// used to measure execution time
    for (i = 0; i < sizeof(testvalues) / sizeof(testvalues[0]); i++) {

        auto  rv =ComputeExpr(testvalues[i].text, result, asgCt);
        if (rv != retCode::EXPR_OK || result != testvalues[i].expected_result) {
            std::cout << "test " << i + 1 << " failed\n" <<
                "expected " << testvalues[i].text << " = " 
                << testvalues[i].expected_result << '\n';
            std::cout << "got " << result << '\n';
            Beep(750, 1000);
        }
    }
    std::cout << i << (lang ? "  pruebas completadas\n" : " tests completed\n");

    for (Znum i = 1000; i <= 100000000000000000; ) {
        Znum x1, x2;
        mpz_nextprime(ZT(x1), ZT(i));  // get next prime
        i *= 10;
        mpz_nextprime(ZT(x2), ZT(i));  // get next prime
        x3 = x1*x2;
        testcnt++;
        factortest(x3, testcnt);

        x3++;
        testcnt++;
        factortest(x3, testcnt);
    }

    /* tests below have shown a problem with pollard-rho for certain numbers */
    long long int PollRho[] = {99999999973789, 183038861417, 183475587821,
        181919916457, 199996999457, 204493418837 };
    for (int i = 0; i < sizeof(PollRho) / sizeof(PollRho[0]); i++) {
        testcnt++;
        factortest(PollRho[i], testcnt);
    }


    /* exercise code specifically for power +/-1 */
    mpz_ui_pow_ui(ZT(x3), 10, 20);  // x3 = 10^20
    x3 -= 1;                        // x3 = 10^20-1
    testcnt++;
    factortest(x3, testcnt);

    x3 += 2;  // x3 = 10^20+1
    testcnt++;
    factortest(x3, testcnt);

    testcnt++;
    ComputeExpr("n(10^10)^2-1", x3, asgCt);  // test power of large number-1
    /* x3 = 10000000019^2 -1 */
    factortest(x3, testcnt);

    testcnt++;
    ComputeExpr("120#-1", x3, asgCt);
    factortest(x3, testcnt);

    testcnt++;
    ComputeExpr("n(10^15)^2", x3, asgCt);  // test power of large number
    factortest(x3, testcnt);

    testcnt++;
    ComputeExpr("n(10^6+20)^1667", x3, asgCt);  // test power of large prime number
    factortest(x3, testcnt);

    testcnt++;
    ComputeExpr("n(10^7)^3*n(10^8)^2", x3, asgCt);  // test powers of two large prime number
    factortest(x3, testcnt);

    testcnt++;
    ComputeExpr("n(3*10^9+50)*n(3*10^10+500)", x3, asgCt);  // test Lehman factorisation
    factortest(x3, testcnt);

    testcnt++;
    ComputeExpr("n(10^15)^3*n(10^14)", x3, asgCt);  // test Lehman factorisation
    factortest(x3, testcnt);

    /* test using carmichael numbers.  */
    unsigned long long int carmichael[] = { 90256390764228001, 18118463305678126129, 
        18265521244069461529,  18349357898532971521, 18308657203978189969 };
    for (int i = 0; i < sizeof(carmichael) / sizeof(carmichael[0]); i++) {
        testcnt++;
        factortest(carmichael[i], testcnt);
    }

    ComputeExpr("16344221851913485532689", x3, asgCt);
    testcnt++;
    factortest(x3, testcnt);  /* 23 digit carmichael number */

    ComputeExpr("56897193526942024370326972321", x3, asgCt);
    testcnt++;
    factortest(x3, testcnt);  /* 29 digit pseudo-prime number */

    /* set x3 to large prime. see https://en.wikipedia.org/wiki/Carmichael_number */
    ComputeExpr("2967449566868551055015417464290533273077199179985304335099507"
        "5531276838753171770199594238596428121188033664754218345562493168782883", x3, asgCt);
    x4 = x3 * (313 * (x3 - 1) + 1) * (353 * (x3 - 1) + 1);
    /* in general numbers > about 110 digits cannot be factorised in a reasonable time 
    but this one can, because a special algorithm just for Carmichael numbers is used. */
    testcnt++;
    factortest(x4, testcnt);   // 397-digit Carmichael number
    std::cout << "factorised 397-digit Carmichael number \n";

    testcnt++;
    ComputeExpr("n(10^24)*n(10^25)*n(10^26)*n(10^27)", x3, asgCt);  
    factortest(x3, testcnt);

    testcnt++;
    ComputeExpr("n(2^97)*n(2^105)", x3, asgCt);
    factortest(x3, testcnt);

    testcnt++;
    /* test reduction of number x3 to squares a^2 + 2*b^2 */
    ComputeExpr("n(10^7)^2 + 2*n(3*10^6)^2", x3, asgCt);
    factortest(x3, testcnt);

    testcnt++;
    /* test reduction of mumber x3 (> 2^64) to squares a^2 + 2*b^2 */
    ComputeExpr("n(2^34)^2 + 2*n(2^34+200)^2", x3, asgCt);
    factortest(x3, testcnt);

    testcnt++;
    /* test reduction of mumber x3 (> 2^64) to squares a^2 + 2*b^2 */
    ComputeExpr("n(10^27)^2 + 2*n(10^32)^2", x3, asgCt);
    factortest(x3, testcnt);

    testcnt++;
    /* test reduction of mumber x3 (> 2^192) to squares  */
    ComputeExpr("n(10^20+477)*n(10^24)*n(10^22)", x3, asgCt);
    factortest(x3, testcnt);

    auto end = clock();   // measure amount of time used
    double elapsed = (double)end - start;
    PrintTimeUsed(elapsed, "Factorisation tests completed. Time used= ");
    printSummary();
}


//find a random 'strong' prime of size 'bits'
//follows the Handbook of applied cryptography
static void gordon(Znum &p, gmp_randstate_t &state, const long long bits) {
    /*
    SUMMARY: a strong prime p is generated.
    1. Generate two large random primes s and t of roughly equal bitlength (see Note 4.54).
    2. Select an integer i0. Find the first prime in the sequence 2it + 1,
    for i = i0; i0 + 1; i0 + 2; : : : (see Note 4.54). Denote this prime by r = 2it+ 1.
    3. Compute p0 = 2(sr-2 mod r)s - 1.
    4. Select an integer j0. Find the first prime in the sequence p0 +2jrs,
    for j = j0; j0 + 1; j0 + 2; : : : (see Note 4.54). Denote this prime by
    p = p0 + 2jrs. Note: if p0 is has a common factor with r this sequence will 
    never find a prime. 
    5. Return(p).

  4.54 Note (implementing Gordon’s algorithm)
    (i) The primes s and t required in step 1 can be probable primes generated by
    Algorithm 4.44. The Miller-Rabin test (Algorithm 4.24) can be used to test
    each candidate for primality in steps 2 and 4, after ruling out candidates
    that are divisible by a small prime less than some bound B. See Note 4.45 for
    guidance on selecting B. Since the Miller-Rabin test is a probabilistic
    primality test, the output of this implementation of Gordon’s algorithm is
    a probable prime.
    (ii) By carefully choosing the sizes of primes s, t and parameters i0, j0,
    one can control the exact bitlength of the resulting prime p. Note that the
    bitlengths of r and s will be about half that of p, while the bitlength of
    t will be slightly less than that of r.
    */
    Znum r, s, t, t2, p0;

    //1. s and t should be about half the bitlength of p
    for (;;) {
        mpz_urandomb(ZT(s), state, bits / 2 - 2);
        mpz_nextprime(ZT(s), ZT(s));        // make s prime
        if (NoOfBits(s) >= (bits / 2 - 5))
            break;  /* try again if s is too small */
    }
    for (;;) {
    mpz_urandomb(ZT(t), state, bits/2 -2);
    mpz_nextprime(ZT(t), ZT(t));        // make t prime
    if (NoOfBits(t) >= (bits / 2 - 5))
        break;  /* try again if t is too small */
    }


    // 2 Find the first prime r in the sequence 2t + 1, 4t+1 6t+1 ...
    t2 = t * 2;
    r = t2 + 1;
    while (mpz_bpsw_prp(ZT(r)) == 0)
    //while (!mpz_likely_prime_p(ZT(r), state, 0))
        r += t2;


    // 3. Compute p0 = 2(sr-2 mod r)s - 1.
    p0 = ((s*r - 2) % r)*s * 2 - 1;

    // 4. Find the first prime p in the sequence p0, p0 +2rs p0+4rs ....,
    p = p0;
    while (gcd(p, r*s) != 1)
        p+=2;  /* if p has any common factors with r or s we need to make an
               adjustment. Otherwise we would never find a prime. */
    while (mpz_bpsw_prp(ZT(p)) == 0)
    //while (!mpz_likely_prime_p(ZT(p), state, 0))
        p += 2 * r*s;
    return;
}

/* generate RSA-style difficult to factor number, size = bits +/- 2 */
static void get_RSA(Znum &x, gmp_randstate_t &state, const long long bits) {
    Znum p, q;
    x = 0;

    /* keep generating random strong primes till we get a pair that produces
    a product of about the right size. */
    while (abs(NoOfBits(x)-bits) >2) {
        gordon(p, state, bits / 2);
        gordon(q, state, bits / 2);
        x = p * q;
    } 
    return;
}

/* test factorisation using pseudo-random numbers of a specified size.
Command format is TEST 2 [p1[,p2[,p3]]] where 
p1 is the number of tests, 
p2 is the size of the numbers to be factored in bits,
if p3  NE 0 the number(s) to be factored consist(s) of 2 approximately same-sized 
prime factors, if p3 = 0 it(they) will be (a) random number(s) that can contain 
any number of factors. 
if p3 <= 1 use fixed random seed value
if p3 = 2 use truly random seed value
if p3 > 2  use p3 as the seed value */
static void doTests2(const std::vector<std::string> &p) {
    long long p1=0;  // number of tests; must be greater than 0
    long long p2=0;  // size of numbers to be factored in bits (>= 48)
    long long p3=0;  // optional, if non-zero generate RSA-style difficult to factor number
    gmp_randstate_t state;
    Znum x;

 
    if (p.size() >= 3)
        p1 = atoi(p[2].data());
    if (p.size() >= 4)
        p2 = atoi(p[3].data());
    if (p.size() >= 5)
        p3 = atoi(p[4].data());

    if (p1 <= 0 ) {
        std::cout << "Use default 2 for number of tests \n";
        p1 = 2;
    }
    if (p2 < 48) {
        std::cout << "Use default 48 for number size in bits \n";
        p2 = 48;
    }
    if (p2 >= 500)
        std::cout << "**warning: factoring such large numbers could take weeks! \n";
    gmp_randinit_mt(state);  // use Mersenne Twister to generate pseudo-random numbers
    if (p3 <= 1)
        gmp_randseed_ui(state, 756128234);
    /* fixed seed means that the exact same tests can be repeated provided
       that the same size of number is used each time */
    if (p3 == 2) {
        std::random_device rd;   // non-deterministic generator
        unsigned long long seedval = rd();
        std::cout << "random generator seed value = " << seedval << '\n';
        gmp_randseed_ui(state, seedval);
    }
    else 
        gmp_randseed_ui(state, p3); /* use supplied value as random seed value*/
    

    auto start = clock();	// used to measure execution time

    results.clear();

    for (int i = 1; i <= p1; i++) {
        std::cout << '\n' << myTime() << "  Test " << i << " of " << p1 << '\n';
        if (p3 == 0)
            mpz_urandomb(ZT(x), state, p2);  // get random number, size=p2 bits
        else
            get_RSA(x, state, p2);  // get RSA type number size p2 bits
        ShowLargeNumber(x, groupSize, true, false);
        std::cout << '\n';
        doFactors(x, true); /* factorise x, calculate number of divisors etc */
        results.back().testNum = i;
    }

    auto end = clock();   // measure amount of time used
    double elapsed = (double)end - start;
    PrintTimeUsed(elapsed, "All tests completed. Time used = ");
    printSummary();
    gmp_randclear(state);  // clear state - avoid memory leakage
}

#ifdef BIGNBR

/* generate large random number, up to 128 bits */
static void largeRand(Znum &a) {
    a = ((long long)rand() << 32) + rand();
    a <<= 64;
    a += ((long long)rand() << 32) + rand();
}

/* generate extra large number, about size*32 bits */
static void XlargeRand(Znum& a, int size) {
    a = 1;
    for (int c = 1; c <= size; c++) {
        a *= rand();
    }
}

/*  1. check basic arithmetic operators for BigIntegers
    2. test BigInteger multiplication with larger numbers
    3. BigInteger division with larger numbers
    4 & 5. Modular Multiplication using Mongomery Encoding (REDC)
*/
static void doTests3(void) {
    Znum a, a1, am, b, b1, bm, mod, p, p2, pm;
    limb aL[MAX_LEN], modL[MAX_LEN], alM[MAX_LEN], al2[MAX_LEN];
    limb one[MAX_LEN];
    int numLen = MAX_LEN-2, l2;

    BigInteger aBI, a1BI, bBI, amBI, pBI;

    auto start = clock();	// used to measure execution time

    std::memset(one, 0, MAX_LEN * sizeof(limb));
    one[0] = 1;                   /* set value of one to 1 */

    srand(421040034);               // seed random number generator 
    
    /* check basic arithmetic operators for BigIntegers */
    a = 1291;        // set starting values (increased each time round loop)
    b = 131;
    for (int ctr = 1; ctr <= 1700; ctr++) {
        /* with the values currently used for a and b, ctr cannot go much beyond 1700
        or overflow would occur. */
        aBI = a;                      // copy value of a to aBI (BigInteger)
        bBI = b;
        BigtoZ(a1, aBI);              // copy value to a1 (Znum)
        assert(a1 == a);              // verify conversion to & from Biginteger

        pBI = aBI + bBI;
        BigtoZ(p, pBI);
        assert(p == (a + b));         // verify addition        

        pBI = aBI;
        pBI += bBI;
        assert(p == (a + b));         // verify addition

        pBI = aBI - bBI;
        BigtoZ(p, pBI);
        assert(p == (a - b));         // verify subtraction

        pBI = aBI;
        pBI -= bBI;
        assert(p == (a - b));         // verify subtraction

        pBI = -aBI;           /* unary - */
        BigtoZ(p, pBI);
        assert(p = -a);

        pBI = aBI * bBI;
        BigtoZ(p, pBI);
        Znum prod;
        Znum diff;
        prod = a * b;
        diff = p ^ prod;   // p XOR prod (0 if they are equal)
        if (diff != 0) {              // verify multiplication
            std::cout << "a = " << a << " b = " << b << '\n';
            gmp_printf("p     = %Zx \na * b = %Zx \n", p, prod);
            gmp_printf("diffference = %Zx \n", diff);
        }

        pBI = aBI;
        pBI *= bBI;
        BigtoZ(p, pBI);               // p = pBI = aBI * bBI
        assert(p == (a * b));         // verify multiplication

        pBI = aBI;
        pBI *= INT_MAX;
        BigtoZ(p, pBI);               // p = pBI = aBI * INT_MAX
        assert(p == a * INT_MAX);	  // verify multiplication by int

        pBI = aBI / bBI;
        BigtoZ(p, pBI);                // p = pBI = aBI / bBI
        if (p != (a / b)) {              // verify division
            std::cout << "a = " << a << "\nb = " << b << '\n';
            std::cout << "p = " << p << "\na/b = " << a / b << '\n';
        }

        pBI = aBI;
        pBI /= bBI;            
        BigtoZ(p, pBI);                  // p = pBI = aBI / bBI
        if (p != (a / b)) {              // verify division
            std::cout << "a = " << a << " b = " << b << '\n';
            std::cout << "p = " << p << "\na/b = " << a / b << '\n';
        }

        pBI = aBI % bBI;
        BigtoZ(p, pBI);               // p = pBI = aBI % bBI
        if (p != (a%b))               // verify modulus
            std::cout << "p = " << p << "\na%b = " << a % b << '\n';

        pBI = aBI;
        pBI %= bBI;
        BigtoZ(p, pBI);                 // p = pBI = aBI % bBI
        if (p != (a%b))               // verify modulus
            std::cout << "p = " << p << "\na%b = " << a % b << '\n';

        /* check BigtoLL function */
        int exp;
        a1 = BigToLL(aBI, exp);  // a1 * 2^exp = aBI
        mpz_div_2exp(ZT(b1), ZT(a), exp);   // b1 = a/(2^exp)
        if (b1 != a1)
            gmp_printf("a = %Zx a1 = %Zx, exp = %d \n", a, a1, exp);
        assert(b1 == a1);

        pBI = aBI << 4;      // check left shift
        BigtoZ(p, pBI);      // p = pBI = aBI << 4
        p2 = a << 4;
        if (p != p2) {
            gmp_printf("a = %#Zx p = %#Zx expected %#Zx \n", a, p, p2);
        }


        pBI = aBI << 56;    // check left shift
        BigtoZ(p, pBI);     // p = pBI = aBI << 56
        p2 = a << 56;
        if (p != p2) {
            gmp_printf("a = %#Zx p = %#Zx expected %#Zx \n", a, p, p2);
        }

        pBI = aBI << 567;      // check left shift
        BigtoZ(p, pBI);        // p = pBI = aBI << 567
        p2 = a << 567;
        if (p != p2) {
            gmp_printf("a = %#Zx p = %#Zx expected %#Zx \n", a, p, p2);
        }

        pBI = aBI + INT_MAX;
        BigtoZ(p, pBI);              // p = pBI = aBI + max;
        assert(p == a + INT_MAX);    // veryify addition of int

        pBI = aBI + INT_MIN;
        BigtoZ(p, pBI);              // p = pBI = aBI + (-max);
        b1 = a + INT_MIN;
        if (p != b1) {
            gmp_printf("a =%Zd p = %Zd expected %Zd  \n", a, p, b1);
        }
        assert(p == b1);       // veryify addition of -ve int

        pBI = aBI - INT_MAX;
        BigtoZ(p, pBI);               // p = pBI = aBI - max;
        assert(p == a - INT_MAX);;    // veryify subtraction of int

        pBI = aBI - (long long)INT_MIN;
        BigtoZ(p, pBI);               // p = pBI = aBI - max;
        b1 = a - (long long)INT_MIN;
        if (p != b1) {
            gmp_printf("a =%Zd p = %Zd expected %Zd  \n", a, p, b1);
        }

        pBI = aBI & bBI.abs();
        BigtoZ(p, pBI);
        if (p != (a & abs(b))) {
            gmp_printf("logical and of %Zx & %Zx failed; \ngot %Zx, expected %Zx \n",
                a, abs(b), p, (a & abs(b)));
        }
        assert(p == (a & abs(b)));

        pBI = aBI;
        pBI &= bBI.abs();
        assert(p == (a & abs(b)));

        a *= 1237953;                 // increase a & b, then repeat
        b *= -129218;
    }

    auto end = clock();              // measure amount of time used
    double elapsed = (double)end - start;
    std::cout << "test stage 1 completed  time used= " << elapsed / CLOCKS_PER_SEC << " seconds\n";
    
    /* test multiplication with larger numbers */
    aBI = 0x7fffffffffffffff;
    a1BI = aBI;
    a = 0x7fffffffffffffff;
    for (bBI = 1; bBI <= 256; bBI++) {
        aBI = aBI * a1BI;
        a = a * 0x7fffffffffffffffLL;
        BigtoZ(a1, aBI);

        if (a1 != a) {
            gmp_printf("aBI = %Zx \n", ZT(a1));
            gmp_printf("  a = %Zx \n", ZT(a));
        }
    }

    /* check multiplication of large numbers */
    aBI = 17341;
    a = 17341;
    for (int i = 1; i < 13; i++) {
        if (aBI.bitLength() > MAX_LEN*BITS_PER_GROUP / 2)
            break;  // exit loop if a is too large to square it
        aBI = aBI * aBI;
        a = a * a;
        BigtoZ(a1, aBI);

        /*gmp_printf("aBI = %Zx \n", ZT(a1));
        gmp_printf("  a = %Zx \n", ZT(a));*/
        assert(a1 == a);
    }

    end = clock();              // measure amount of time used
    elapsed = (double)end - start;
    std::cout << "test stage 2 completed  time used= " << elapsed / CLOCKS_PER_SEC << " seconds\n";

    
    /* check conversion to & from floating point */
    /*for (int i = 1; i <= 30; i++) {
        expBigInt(amBI, i);
        std::cout << "e^=" << amBI << " expected " << std::exp(i) << '\n';
    }*/

    p = 10;

    Znum error, relerror;
    for (int i = 1; i < 99;  i++, p *= 1000000 ) {
        if (!ZtoBig(pBI, p))
            break;            // p too large to convert so stop
        double pdb = pBI.log();  // get natural log of p
        /*if (pdb > 708)
            break;*/
        expBigInt(amBI, pdb);   // convert back from log
        BigtoZ(am, amBI);       // convert back to Znum
        error = p - am;         // get error
        if (error != 0) {
            double e1, e2, relErrf;
            long e1l, e2l;
            e1 = mpz_get_d_2exp(&e1l, ZT(p));
            e2 = mpz_get_d_2exp(&e2l, ZT(am));
            assert(e1l == e2l);   // check power of 2 exponent is correct
            relErrf = abs(e1 - e2) / e1;   /* should be less than 10E-14 */
            if (relErrf > 10E-13) {
                relerror = (10'000'000'000'000'000LL * error) / p;  /* error * 10^15 */
                /* print log base 10 of p, then error ratio */
                std::cout << " pdb = " << pdb / std::log(10)
                    << " error ratio = " << relErrf << " error*10^15 = " << relerror << '\n';
            }
        }
    }

    /* check division with large numbers */
    p = 12345678901;
    p2 = 11;
    for (l2 = 1; l2 < 4560; l2++) {
        if (!ZtoBig(pBI, p))
            break;                 // exit loop if p is too big to fit BigInteger (approximately 20000 digits)
        if (!ZtoBig(amBI, p2))
            break;
        bBI = pBI / amBI;         // calculçate p/p2 using Bigintegers
        BigtoZ(b1, bBI);          // convert quotient back to Znum in b1
        if (b1 != p / p2) {
            std::cout << "p      = " << p
                << "\np2     = " << p2
                << "\nexpected " << p / p2
                << "\ngot      " << b1;
            break;
        }
        bBI = pBI % amBI;
        BigtoZ(b1, bBI);
        if (b1 != p % p2) {
            std::cout << "p      = " << p
                << "\np2     = " << p2
                << "\nexpected " << p % p2
                << "\ngot      " << b1;
            break;
        }
        p *= (long long)rand() * 2 + 3;
        p2 *= rand();
    }
    std::cout << "division test " << l2 << " cycles\n";
    end = clock();              // measure amount of time used
    elapsed = (double)end - start;
    std::cout << "test stage 3 completed  time used= " << elapsed / CLOCKS_PER_SEC << " seconds\n";

    /* test to compare modular multiplication using DA's code against GMP/MPIR
    BigIntegers. Both use Mongomery notation for the integers to avoid slow
    division operations. The conclusion is that GMP takes about twice as long
    as DA's code. */
    
    modmultCallback = nullptr;      // turn off status messages from modmult
    for (int c = 1; c <= 100; c++) {
        /* set up modulus and Mongomery parameters */
        XlargeRand(mod, c);					// get large random number (up to 32 * c bits)
        mod |= 1;                       // set lowest bit (make sure mod is odd)
        GetMontgomeryParms(mod);   /* set up zR, zR1, zR2, zNI */
        numLen = MAX_LEN - 2;
        ZtoLimbs(modL, mod, MAX_LEN);    // copy value of mod to modL
        while (modL[numLen - 1] == 0)
            numLen--;                    // adjust length i.e. remove leading zeros
        std::memcpy(TestNbr, modL, numLen * sizeof(limb));  // set up for GetMontgomeryParms
        NumberLength = numLen;
        GetMontgomeryParms(numLen);
        XlargeRand(a, c);				     // get large random number a
        a %= mod;						 // ensure a < mod
        modmult(a, zR2, am);             // convert a to Montgomery (Znum) in am
        numLen = MAX_LEN - 2;
        ZtoLimbs(aL, a, numLen);		     // copy value of a to aL (limbs)
        while (aL[numLen - 1] == 0)
            numLen--;                    // adjust length i.e. remove leading zeros
        //NumberLength = numLen;
        modmult(aL, MontgomeryMultR2, alM);  // convert a to Mongomery (limbs)
        modmult(alM, one, al2);          // convert a from Mongomery (limbs) 
        LimbstoZ(al2, a1, numLen);       // copy value to a1 (Znum)
        assert(a == a1);                 // check that all these conversions work properly
    }
    end = clock();              // measure amount of time used
    elapsed = (double)end - start;
    std::cout << "test stage 4 completed  time used= " << elapsed / CLOCKS_PER_SEC << " seconds\n";

    /* test square root function */
    a = 1;
    for (int i = 1; i <= 20; i++) {
        a *= rand();
        b = sqrt(a);         /* get square root of Znum, use Boost Library function */
        aBI = a;             /* convert Znum to BigInteger*/
        bBI = aBI.sqRoot();  /* get square root of big integer, use DA's function */
        BigtoZ(bm, bBI);     /* convert BigInteger to Znum */
        std::cout << "sqrt(" << a << ") = " << b << '\n';
        assert(bm == b);      /* check that both square roots have the same value */
    }

    //largeRand(mod);				     // get large random number  b
    //mod |= 1;                         /* make sure b is odd */
    //std::cout << "Mongomery modular multiplication. modulus = " << mod << '\n';
    //GetMontgomeryParms(mod);      /* set up zR, zR1, zR2, zNI */
    //numLen = MAX_LEN - 2;
    //ZtoLimbs(modL, mod, MAX_LEN);    // copy value of mod to modL
    //while (modL[numLen - 1] == 0)
    //    numLen--;                    // adjust length i.e. remove leading zeros
    //memcpy(TestNbr, modL, numLen * sizeof(limb));  // set up for GetMontgomeryParms
    //NumberLength = numLen;
    //GetMontgomeryParms(numLen);

    //a %= mod;
    //modmult(a, zR2, am);      // convert a to montgomery
    //ZtoLimbs(aL, a, numLen);

    //largeRand(b);
    //b %= mod;
    //modmult(b, zR2, bm);             // convert b to Montgomery (Znum)
    //numLen = MAX_LEN - 2;
    //ZtoLimbs(bL, b, numLen);		 // copy value of b to bL
    //
    //while (bL[numLen - 1] == 0)
    //    numLen--;                    // adjust length i.e. remove leading zeros
    ////NumberLength = numLen;

    ////memcpy(TestNbr, bL, numLen * sizeof(limb));  // set up for GetMontgomeryParms
    ////GetMontgomeryParms(numLen);
    //modmult(bL, MontgomeryMultR2, blM);  // convert b to Mongomery (limbs)

    //for (int i = 1; i < 200000000; i++) {
    //    modmult(alM, blM, plm);                   // p = a*b modulo mod (limbs)
    //    memcpy(alM, plm, numLen * sizeof(limb));  // a = p (limbs)
    //    modmult(am, bm, pm);                      // p = a*b modulo mod (Znum)
    //    am = pm;						          //a = p (Znum)
    //    REDC(p, pm);					 // convert p from montgomery (Znum)
    //    modmult(plm, one, pl);           // convert p from Mongomery (limbs)
    //    LimbstoZ(pl, p2, numLen);        // convert value of p to p2 (Znum)
    //    assert(p2 == p);
    //    if (i % 20000000 == 0) {
    //        std::cout << "test stage 5 " << i / 2000000 << "% complete \n";
    //    }
    //    if (p == 0)
    //        break;
    //}

    ////REDC(p, pm);					 // convert p from montgomery (Znum)
    ////modmult(plm, one, pl);           // convert p from Mongomery (limbs)
    ////LimbstoZ(pl, p2, numLen);        // convert value of p to p2 (Znum)
    ////assert(p2 == p);

    end = clock();              // measure amount of time used
    elapsed = (double)end - start;
    std::cout << "tests completed  time used= " << elapsed / CLOCKS_PER_SEC << " seconds\n";
}
#endif

/* tests for r3 function. This test exploits the fact that R3 and R3h calculate
the same value by completely different methods. R3 is too slow to be practical
for numbers greater than about 11 digits, but R3h requires that pari/GP has 
been installed. */
/* see https://oeis.org/A005875 
Command format is TEST 11 [p1[,p2[,p3]]] where
p1 is the number of tests,
p2 is the size of the numbers to be processed in bits,
if p3 <= 1 use fixed random seed value (default)
if p3 = 2 use truly random seed value
if p3 > 2  use p3 as the seed value*/
static void doTestsB(const std::vector<std::string> &p) {
    int i;
    auto start = clock();	// used to measure execution time
    long long p1 = 0;  // number of tests; must be greater than 0, default is 20
    long long p2 = 0;  // size of numbers to be tested, in bits (default is 32, maximum is 48)
    long long p3 = 0;
    long long xl;
    unsigned long long rv;
    gmp_randstate_t state;  /* use gmp/mpir random number generator */
    Znum x, rv3;
    /* convert parameters to binary integers */
    if (p.size() >= 3)
        p1 = atoi(p[2].data());
    if (p.size() >= 4)
        p2 = atoi(p[3].data());
    if (p.size() >= 5)
        p3 = atoi(p[4].data());
    /* check whether parameter values are within acceptable ranges */
    if (p1 <= 0) {
        std::cout << "Use default 20 for number of tests \n";
        p1 = 20;
    }
    if (p2 <= 7 || p2 >41) {
        std::cout << "Use default 32 for number size in bits \n";
        p2 = 32;
    }
    /* use value of p3 to seed the the 'random' number generator */
    gmp_randinit_mt(state);  // use Mersenne Twister to generate pseudo-random numbers
    if (p3 <= 1)
        gmp_randseed_ui(state, 756128234);
    /* fixed seed means that the exact same tests can be repeated provided
       that the same size of number is used each time */
    else if (p3 == 2) {
        std::random_device rd;   // non-deterministic generator
        unsigned long long seedval = rd();
        std::cout << "random generator seed value = " << seedval << '\n';
        gmp_randseed_ui(state, seedval);
        }
        else
            gmp_randseed_ui(state, p3); /* use supplied value as random seed value*/

    for (i = 1; i <= p1; i++) {
        mpz_urandomb(ZT(x), state, p2);  // get random number, size=p2 bits
        xl = MulPrToLong(x);
        rv = R3(xl);
        rv3 = R3h(x);
        if (rv3 != rv || verbose > 1 || p2 >= 37) {
            std::cout << myTime() <<  " R3(" << xl << ") =" << rv 
                << "; R3h(" << x << ") = " << rv3 << '\n';
        }
    }
    std::cout << "R3 " << p1 << " tests completed \n";
    auto end = clock();              // measure amount of time used
    auto elapsed = (double)end - start;
    PrintTimeUsed(elapsed, " time used = ");
    return;
}

/* test factorisation of mersenne numbers see https://en.wikipedia.org/wiki/Mersenne_prime
this test will take about 2.5 hours 
*/
static void doTests4(void) {
    int px = 0;
    const int pmax = 75;
    Znum m;

    results.clear();
    auto start = clock();	// used to measure execution time
    if (primeListMax < 1000)
        generatePrimes(393203);

    for (px = 0; px <= pmax; px++) {
        if (px >= 70 && !yafu) {
            std::cout << "remaining tests skipped \n";
            break;
        }
        std::cout << '\n' << myTime() << " test " << px + 1 << " of " << pmax + 1;
        mpz_ui_pow_ui(ZT(m), 2, primeList[px]);  // get  m= 2^p
        m--;                // get 2^p -1
        if (factortest(m, px+1)) /* factorise m, calculate number of divisors etc */
            std::cout << "2^" << primeList[px] << " - 1 is prime \n";
    }
    auto end = clock();              // measure amount of time used
    auto elapsed = (double)end - start;
    PrintTimeUsed(elapsed, "tests completed time used = ");
    printSummary();           // print summary 1 line per test
}

/* tests using only YAFU for factorisation */
static void doTests5(void) {
    results.clear();
    int testcnt = 1, asgCt;

    Znum num = 49728103; // 7001 * 7103
    factortest(num, testcnt, 1);  // test 1

    // factorise 127 digit number   test 2
    /*  =                               280673 
    *                              2756 163353 
    *                            598990 818061 
    *                          4 527716 228491 
    *                     248158 049830 971629 
    * 33637 310674 071348 724927 955857 253537  
    *117445 937227 520353 139789 517076 610399  
    (7 factors)   */
    testcnt++;
    ComputeExpr("2056802480868100646375721251575555494408897387375737955882170045672576386016591560879707933101909539325829251496440620798637813", 
        num, asgCt);
    factortest(num, testcnt, 1);

    //factorise 57 digit number      test 3
    /* P6 = 280673
      P12 = 598990818061
      P10 = 2756163353
      P13 = 4527716228491
      P18 = 248158049830971629
    */
    testcnt++;
    ComputeExpr("520634955263678254286743265052100815100679789130966891851", num, asgCt);
    factortest(num, testcnt, 1);

    //factorise 80 digit number (about 3 minutes)   test 4
    /* P49 = 2412329883909990626764837681505523221799246895133
       P32 = 18138544144503826077310252140817
    */
    testcnt++;
    ComputeExpr("43756152090407155008788902702412144383525640641502974083054213255054353547943661", num, asgCt);
    factortest(num, testcnt, 1);

    //factorise 85 digit number (about 7 mins)   // test 5
    /* factors are 1485325304578290487744454354798448608807999 and 
                   1263789702211268559063981919736415575710439 */
    testcnt++;
    ComputeExpr("1877138824359859508015524119652506869600959721781289179190693027302028679377371001561", num, asgCt);
    factortest(num, testcnt, 1);

    // factorise 94 digit number (about 60 mins)    test 6
    /* factors are 10910042366770069935194172108679294830357131379375349 and 
                   859735020008609871428759089831769060699941 */
    testcnt++;
    ComputeExpr("9379745492489847473195765085744210645855556204246905462578925932774371960871599319713301154409", num, asgCt);
    factortest(num, testcnt, 1);

    //factorise 100 digit number - takes many hours    test 7
    /* factor are 618162834186865969389336374155487198277265679 and
                  4660648728983566373964395375209529291596595400646068307 */
    testcnt++;
    ComputeExpr("2881039827457895971881627053137530734638790825166127496066674320241571446494762386620442953820735453", num, asgCt);
    factortest(num, testcnt, 1);

    printSummary();    // print summary - 1 line per test

    return;
}

/* tests using Msieve for factorisation. Factorise selected Mersenne numbers 
if useMsieve = FALSE use only YAFU. Normal yafu & msieve flags are ignored. */
static void doTests6(bool usesMsieve = true) {
    bool yafusave = yafu;
    bool msievesave = msieve;
    bool Parisave = Pari;
    Znum m;
    msieve = usesMsieve;
    yafu = !usesMsieve;
    Pari = false;
    int testcnt = 1;

    results.clear();
    auto start = clock();	// used to measure execution time

    mpz_ui_pow_ui(ZT(m), 2, 277);  // get  m= 2^p
    m--;                // get 2^p -1
    factortest(m, testcnt);      // 84 digits
    /* p7  factor: 1121297
       p38 factor: 31133636305610209482201109050392404721
       p40 factor: 6955979459776540052280934851589652278783
    */

    testcnt++;
    mpz_ui_pow_ui(ZT(m), 2, 293);  // get  m= 2^p
    m--;                // get 2^p -1
    factortest(m, testcnt);      // 89 digits
    /* p26 factor: 40122362455616221971122353
       p63 factor: 396645227028138890415611220710757921643910743103031701971222447 */

    testcnt++;
    mpz_ui_pow_ui(ZT(m), 2, 311);  // get  m= 2^p
    m--;                // get 2^p -1
    factortest(m, testcnt);      // 94 digits
    /* p7  factor: 5344847
       p31 factor: 2647649373910205158468946067671
       p57 factor: 294803681348959296477194164064643062187559537539328375831
    */

    testcnt++;
    mpz_ui_pow_ui(ZT(m), 2, 313);  // get  m= 2^p
    m--;                // get 2^p -1
    factortest(m, testcnt);      // 95 digits
    /* p8  factor: 10960009
       p17 factor: 14787970697180273
       p25 factor: 3857194764289141165278097
       p47 factor: 26693012026551688286164949958620483258358551879
    */

    testcnt++;
    mpz_ui_pow_ui(ZT(m), 2, 349);  // get  m= 2^p
    m--;                // get 2^p -1
    factortest(m, testcnt);      // 106 digits
    /* p34 factor: 2927455476800301964116805545194017
       p67 factor: 6725414756111955781503880188940925566051960039574573675843402666863
    */
    
    yafu = yafusave;
    msieve = msievesave;
    Pari = Parisave;

    auto end = clock();              // measure amount of time used
    auto elapsed = (double)end - start;
    PrintTimeUsed(elapsed, "tests completed time used = ");
    printSummary();           // print 1 line per test
}

/* Lucas-Lehmer test*/
static void doTests7(const std::vector<std::string> &p) {
    std::vector <long long> mPrimes;
    int i=0;
    lltPriCnt = lltCmpCnt = lltTdivCnt = 0;  /* reset counters */
#ifdef _DEBUG
    int limit = 2000;  /* find 1st 15 Mersenne  primes, test 303 primes */
#else
    int limit = 12000; /* find 1st 23 Mersenne primes, test 1438 primes */
#endif
    auto start = clock();	// used to measure execution time

    if (p.size() >= 3)
        limit = atoi(p[2].data());
    if (limit < 0 || limit > 120000) {
        std::cout << "limit out of range; use 12000 \n";
        limit = 12000;
    }
    generatePrimes(limit);    /* make prime list if not already done */

    for (i = 0; primeList[i] < limit; i++) {
        Znum p = primeList[i];
        Znum rv = llt(p); /* Return 0 if 2^p-1 is composite, 1 if prime  */
        if (rv == 1) {
            std::cout << myTime() << " 2^" << primeList[i] << " -1 is prime *** \n";
            mPrimes.push_back(primeList[i]);
        }
        else if (verbose > 0 || (i & 0x3f) == 0)
            /* \r instead of usual \n means that each messsage overwrites the 
            previous one */
            std::cout << myTime() << " 2^" << primeList[i] << " -1 is NOT PRIME \r";
    }

    /* print the results */
    /* see https://oeis.org/A000043 */
    std::cout << "Found " << mPrimes.size() << " Mersenne primes  out of " << i 
        << " numbers tested \n" ;
    for (auto p : mPrimes) {
        std::cout << p << ", ";
    }
    putchar('\n');
    long long other = (long long)i - (lltPriCnt+ lltTdivCnt+ lltCmpCnt);
    printf_s("%5.2f%% primes found by llt \n", 100.0 *lltPriCnt / i);
    if (lltTdivCnt != 0)
        printf_s("%5.2f%% composites found by trial division \n", 100.*lltTdivCnt / i);
    if (lltCmpCnt != 0)
        printf_s("%5.2f%% composites found by llt \n", 100.0*lltCmpCnt / i);
    if (other != 0)
        printf_s("%5.2f%% other \n", 100.0*other / i);

    auto end = clock();              // measure amount of time used
    auto elapsed = (double)end - start;
    PrintTimeUsed(elapsed, "test 7 completed time used = ");
}


std::vector<std::string> inifile;  // copy some contents of .ini here
std::string iniPath;          // path to .ini file

  /* can be overwritten from the .ini file */
std::string helpFilePath = "docfile.txt";
std::string PariPath = "C:/Program Files (x86)/Pari64-2-13-2/libpari.dll";  

/* (re)write the BigIntCalculator.ini file
initially a .new file is created, then any .old file is deleted
then the current .ini if any is renamed as .old,
then the .new file is renamed as .ini */
void writeIni(void) {
    std::string newFname = iniPath + "BigIntCalculator.new";
    std::string iniFname = iniPath + "BigIntCalculator.ini";
    std::string oldFname = iniPath + "BigIntCalculator.old";

    std::ofstream newStr(newFname, std::ios::out);  // open .new file for output
    if (!newStr.is_open()) {
        std::cerr << "cannot open BigIntCalculator.new \n";
        return;
    }
    /* copy comments etc from current ini file to new one */
    for (auto text : inifile) {
        newStr << text << '\n';
    }
    newStr << "yafu-path=" << YafuPath << '\n';
    newStr << "yafu-prog=" << yafuprog << '\n';
    newStr << "yafu-out=" << outPath << '\n';
    newStr << "msieve-path=" << MsievePathS << '\n';
    newStr << "msieve-prog=" << MsieveProg << '\n';
    newStr << "msieve-log=" << MsieveLogPath << '\n';
    newStr << "helpfile=" << helpFilePath << '\n';
    newStr << "paripath=" << PariPath << '\n';
    newStr << "endsound=" << endsound << '\n';
    newStr << "attsound=" << attsound << '\n';
    newStr.close();

      // delete any previous .old
    int rc = remove(oldFname.c_str());
    if (rc != 0 && errno != ENOENT) {
        perror("could not remove BigIntCalculator.old file ");
    }
    // rename .ini as .old
    int rv = rename(iniFname.c_str(), oldFname.c_str());   
    if (rv == 0 || errno == ENOENT) {
        int rv2 = rename(newFname.c_str(), iniFname.c_str());   // .new -> .ini
        if (rv2 != 0)
            perror("unable to rename BigIntCalculator.new as BigIntCalculator.ini");
        else if (verbose > 0)
            std::cout << myTime() << " BigIntCalculator.ini written to disk \n";
    }
    else
        perror("unable to rename BigIntCalculator.ini as BigIntCalculator.old");
}

/* read the .ini file and update paths. 
path definitions begin with yafu-path=, yafu-prog=, msieve-path=, msieve-prog=, 
helpfile=, endsound=, attsound= or paripath=
Paths are not case-sensitive. 
Anything else is saved and is copied if the .ini file is updated 

arg is arg[0] of the call to main, which conveniently contains the full path
for the .exe file. We use the same path for the .ini file. 
If the .ini file can't be opened a new one is created */
static void processIni(const char * arg) {
    std::string iniFname;
    std::string buffer;

    iniPath = arg;  // save path in a global variable
    auto b = iniPath.find_last_of("\\/"); // find last / character
    if (b != std::string::npos) {
        iniPath.erase(b + 1);  // remove everything after /
        iniFname = iniPath + "BigIntCalculator.ini";
    }
    else 
        iniFname = "BigIntCalculator.ini";

    std::ifstream iniStr(iniFname, std::ios::in);  // open .ini file
    if (!iniStr.is_open()) {
        /* can't open .ini file - make one from scratch */
        const time_t currtime = time(NULL);  // time as seconds elapsed since midnight, January 1, 1970
        struct tm mytime;
        char timestamp[23];   // date & time in format "dd/mm/yyyy at hh:mm:ss"

        localtime_s(&mytime, &currtime);  // convert to tm format
        /* convert to dd/mm/yyyy hh:mm:ss */
        strftime(timestamp, sizeof(timestamp), "%d/%m/%C%y at %H:%M:%S", &mytime);
        buffer = "%file originally created on ";
        buffer += timestamp;     //  "dd/mm/yyyy at hh:mm:ss"
        inifile.push_back(buffer);
        writeIni();  // create a new .ini file using hard coded initial values

        std::cout << " cannot open BigIntCalculator.ini - create new file \n";
    }
    else {    /* read in .ini file */
        while (std::getline(iniStr, buffer)) {

            if (_strnicmp("yafu-path=", buffer.c_str(), 10) == 0) {
                YafuPath = buffer.substr(10); // copy path following '=' character
            }
            else if (_strnicmp("yafu-prog=", buffer.c_str(), 10) == 0) {
                yafuprog = buffer.substr(10); // copy path following '=' character
            }
            else if (_strnicmp("yafu-out=", buffer.c_str(), 9) == 0) {
                outPath = buffer.substr(9); // copy path following '=' character
            }
            else if (_strnicmp("msieve-path=", buffer.c_str(), 12) == 0) {
                MsievePathS = buffer.substr(12); // copy path following '=' character
            }
            else if (_strnicmp("msieve-prog=", buffer.c_str(), 12) == 0) {
                MsieveProg = buffer.substr(12); // copy path following '=' character
            }
            else if (_strnicmp("msieve-log=", buffer.c_str(), 11) == 0) {
                MsieveLogPath = buffer.substr(11); // copy path following '=' character
            }
            else if (_strnicmp("helpfile=", buffer.c_str(), 9) == 0) {
                helpFilePath = buffer.substr(9);
            }
            else if (_strnicmp("endsound=", buffer.c_str(), 9) == 0) {
                endsound = buffer.substr(9);
            }
            else if (_strnicmp("attsound=", buffer.c_str(), 9) == 0) {
                attsound = buffer.substr(9);
            }
            else if (_strnicmp("paripath=", buffer.c_str(), 9) == 0) {
                PariPath = buffer.substr(9);
            }
            else inifile.push_back(buffer);  // save anything not recognised
        }
        iniStr.close();
    }
}

/* format is IF (expression) REPEAT 
          or IF (expression) STOP 
          or IF (expression) THEN (expression) [ELSE (expression)]
          the [ELSE (expression)] part is optional
return -1 if syntax is invalid
 return 0 if expression is 0 for STOP or REPEAT
 return 1 if expression NE 0 and action is STOP
 return 2 if expression NE 0 and action is REPEAT 
 return 3 if THEN or ELSE expression evaluated succesfully */
static int ifCommand(const std::string &command) {
    int ixx, ixx2, exprLen;
    int asgCt = 0;  /* number of assignment operators*/
    std::string expr;
    Znum result;
    int bc = 1;  /* bracket count */

    /* find ( after IF */
    for (ixx = 2; ixx < command.size(); ixx++)
        if (command[ixx] == '(')
            break;
    if (ixx >= command.size()) {
        std::cout << "No ( after IF command \n";
        return -1;
    }

    /* find matching ) after IF */
    for (ixx2 = ixx + 1; ixx2 < command.size(); ixx2++) {
        if (command[ixx2] == '(')
            bc++;
        if (command[ixx2] == ')') {
            bc--;
            if (bc == 0)
                break;   /* found matching closing bracket */
        }
    }
    if (ixx2 >= command.size()) {
        std::cout << "No  matching ) after IF command \n";
        return -1;
    }

    /* evaluate expression betwen brackets */
    exprLen = ixx2 - ixx - 1;
    expr = command.substr(ixx + 1, exprLen);
    retCode rv = ComputeExpr(expr, result, asgCt);
    if (rv != retCode::EXPR_OK) {
        textError(rv);   // invalid expression; print error message
        return -1;
    }

    /* move ixx2 to next non-blank character */
    for (ixx2++; ixx2 < command.size() && isblank(command[ixx2]); ixx2++);
    if (command.substr(ixx2) == "STOP") {
        if (result != 0)
            return 1;
        else
            return 0;
    }
    if (command.substr(ixx2) == "REPEAT") {
        if (result != 0)
            return 2;
        else
            return 0;
    }
    if (command.substr(ixx2, 4) == "THEN") {
        size_t ex1Start=ixx2+4, ex1End, ex2Start, ex2End;

        /* move ex1Start to next non-blank character */
        for (ex1Start=ixx2+4; ex1Start < command.size() && isblank(command[ex1Start]); 
            ex1Start++);
    
        if (ex1Start >= command.size() || command[ex1Start] != '(') {
            std::cout << " THEN is not followed by ( \n";
            return -1;  /* no ( after THEN so invalid */
        }
        ex1End = ex1Start + 1;
        bc = 1;

        /* find matching ) */
        while (ex1End < command.size()) {
            if (command[ex1End] == '(')
                bc++;
            if (command[ex1End] == ')')
                bc--;
            if (bc == 0)
                break;
            ex1End++;
        }
        if (bc != 0) {
            std::cout << "Matching ) not found \n";
            return -1;   /* matching ) not found */
        }

        /* move ex2Start to next non-blank character if any */
        for (ex2Start = ex1End + 1; ex2Start < command.size() && isblank(command[ex2Start]);
            ex2Start++);
        if (ex2Start < command.size()) {
            /* still some unprocessed characters */

            if (command.size() < ex2Start + 4 || command.substr(ex2Start, 4) != "ELSE") {
                std::cout << "Format invalid. ELSE not found \n";
                return -1;
            }
            ex2Start += 4;   /* move past ELSE */
            for (; ex2Start < command.size() && isblank(command[ex2Start]); ex2Start++);
            if (ex2Start >= command.size() || command[ex2Start] != '(') {
                std::cout << "ELSE not followed by ( \n";
                return -1;  /* no ( after ELSE so invalid */
            }
            bc = 1;
            ex2End = ex2Start + 1;
            while (ex2End < command.size()) {
                if (command[ex2End] == '(')
                    bc++;
                if (command[ex2End] == ')')
                    bc--;
                if (bc == 0)
                    break;
                ex2End++;
            }

            if (bc != 0) {
                std::cout << "Matching ) not found \n";
                return -1;   /* matching ) not found */
            }
        }
        if (result != 0) {
            expr = command.substr(ex1Start+1, ex1End - ex1Start-1);
        }
        else if (ex2Start < command.size())
            expr = command.substr(ex2Start+1, ex2End - ex2Start-1);
        else 
            return 3;  /* no ELSE expression to evaluate */

        retCode rv = ComputeMultiExpr(expr, result);
        return 3;  /* IF (...) THEN (...) ELSE (...) processed OK*/
    }

    std::cout << "Neither STOP, REPEAT, nor THEN found \n";
    return -1;  /* neither STOP nor REPEAT nor THEN found */
}

// Description:
//   Creates a tooltip for an item in a dialog box. 
// Parameters:
//   idTool - identifier of an dialog box item.
//   hDlg - window handle of the dialog box.
//   pszText - string to use as the tooltip text.
//   hParnt = parent instance handle
// Returns:
//   The handle to the tooltip.
//
static HWND CreateToolTip(int toolID, HWND hDlg, LPWSTR pszText, HINSTANCE hParnt)
{
    if (!toolID || !hDlg || !pszText)
    {
        return FALSE;
    }
    // Get the window of the tool.
    HWND hwndTool = GetDlgItem(hDlg, toolID);

    // Create the tooltip. 
    HWND hwndTip = CreateWindowEx(NULL, TOOLTIPS_CLASS, NULL,
        WS_POPUP | TTS_ALWAYSTIP /* | TTS_BALLOON  */ ,
        CW_USEDEFAULT, CW_USEDEFAULT,
        CW_USEDEFAULT, CW_USEDEFAULT,
        hDlg, NULL,
        hParnt, NULL);

    if (!hwndTool || !hwndTip)
    {
        return (HWND)NULL;
    }

    auto rv = SendMessage(hwndTip, TTM_ACTIVATE, TRUE, 0);
    if (rv != TRUE)
        ErrorDisp(__FUNCTION__ " TTM ACTIVATE");

    /* typedef struct tagTOOLINFOA {
  UINT      cbSize;  // size in bytes 
  UINT      uFlags;
  HWND      hwnd;   // handle to the window that contains the tool
  UINT_PTR  uId;    // tool identifier or handle 
  RECT      rect;   // ignored if TTF_IDISHWND is set
  HINSTANCE hinst;  // Handle to the instance that contains the string 
                    // resource for the tool.
  LPSTR     lpszText;  // text for the  tooltip
  LPARAM    lParam;    // application defined
  void      *lpReserved;
} TTTOOLINFOA, *PTOOLINFOA, *LPTTTOOLINFOA;  */

    // Associate the tooltip with the tool.
    TOOLINFO toolInfo = { 0 };
    toolInfo.cbSize = sizeof(toolInfo);
    toolInfo.hwnd = hDlg;  // handle to the window that contains the tool
    /* indicate that uId is the tool handle */
    toolInfo.uFlags = TTF_IDISHWND | TTF_SUBCLASS;
    toolInfo.uId = (UINT_PTR)hwndTool;
    toolInfo.lpszText = pszText;  /* text for the  tooltip */
    rv = SendMessage(hwndTip, TTM_ADDTOOL, 0, (LPARAM)&toolInfo);
    if (rv != TRUE)
        ErrorDisp(__FUNCTION__ " TTM_ADDTOOL");
    return hwndTip;
}

/* process message from dialog box e.g. WM_COMMAND, WM_INITDIALOG,
* additionalInfo1 = Menu identifier, or control notification code + control identifier
additionalInfo2 = 0 or handle to control window.
returns TRUE or FALSE. 
Note that TRUE and FALSE are macros  used with C-stye BOOL variables, but true and false
are C++ language words used with C++ bool variables. */
static INT_PTR SetDialogAct(HWND DiBoxHandle,
    UINT message,
    WPARAM additionalInfo1,
    LPARAM additionalInfo2) {

    static HWND hh;  /* tooltip handle */
    int wpHi = HIWORD(additionalInfo1);   /* 0 or control notification code */
    int wpLo = LOWORD(additionalInfo1);   /* menu or control identifier */
    // plan name can be NONE, NOECM, LIGHT, NORMAL, DEEP
    char planText[][7] = { "none", "noECM", "light", "normal", "deep" };
    int planTextSize = sizeof(planText) / sizeof(planText[0]);
    extern int pvalue;  // YAFU PLAN parameter value 4 = PLAN NORMAL (default)
    extern bool eopt;   // set -e option in Msieve: perform 'deep' ECM, seek factors > 15 digits
    int good;
    int temp;

    switch (message) {
    case WM_DESTROY:       /* 0x02 */
    case WM_MOVE:          /* 0x03 */
    case WM_ACTIVATE:      /* 0x06 */
    case WM_SETFOCUS:      /* 0x07 */
    case WM_KILLFOCUS:     /* 0x08 */
    case WM_PAINT:         /* 0x0f */
    case WM_CLOSE:         /* 0x10 */
    case WM_ERASEBKGND:    /* 0x14 */
    case WM_SHOWWINDOW:    /* 0x18 */
    case WM_ACTIVATEAPP:   /* 0x1c */
    case WM_CANCELMODE:    /* 0x1f */
    case WM_SETCURSOR:     /* 0x20 */
    case WM_MOUSEACTIVATE:  /* 0x21 */
    case WM_GETMINMAXINFO:  /* 0x24 */
    case WM_SETFONT:       /* 0x30 */
    case WM_WINDOWPOSCHANGING:  /* 0x46 */
    case WM_WINDOWPOSCHANGED:  /* 0x47 */
    case WM_NOTIFY:        /* 0x4e */
    case WM_HELP:          /* 0x53  (F1 key pressed) */
    case WM_NOTIFYFORMAT:  /* 0x55 */
    case WM_GETICON:       /* 0x75 */
    case WM_NCDESTROY:     /* 0x82 */
    case WM_NCHITTEST:     /* 0x84 */
    case WM_NCPAINT:       /* 0x85 */
    case WM_NCACTIVATE:    /* 0x86 */
    case 0x90:             /* can't find any documentation for this */
    case WM_NCMOUSEMOVE:   /* 0xa0 */
    case WM_NCLBUTTONDOWN: /* 0xa1 */
    case WM_NCLBUTTONUP:   /* 0xa2 */
    case WM_NCLBUTTONDBLCLK: /* 0xa3 */
    case 0xae:             /* can't find any documentation for this */
    case WM_SYSCOMMAND:    /* 0x112 */
    case WM_INITMENU:      /* 0x116 */
    case WM_MENUSELECT:    /* 0x11f */
    case WM_ENTERIDLE:     /* 0x121 */
    case WM_CHANGEUISTATE: /* 0x127 */
    case WM_UPDATEUISTATE: /* 0x128 */
    case WM_QUERYUISTATE:  /* 0x129*/
    case WM_CTLCOLOREDIT:  /* 0x133 */
    case WM_CTLCOLORLISTBOX: /* 0x134 */
    case WM_CTLCOLORBTN:   /* 0x135 */
    case WM_CTLCOLORDLG:   /* ox136 */
    case WM_CTLCOLORSTATIC: /* 0x138 */
    case WM_MOUSEFIRST:     /* 0x200 */
    case WM_LBUTTONDOWN:    /* 0x201 */
    case WM_LBUTTONUP:      /* 0x202 */
    case WM_LBUTTONDBLCLK:  /* 0x203 */
    case WM_ENTERMENULOOP:  /* 0x211 */
    case WM_EXITMENULOOP:   /* 0x212 */
    case WM_CAPTURECHANGED: /* 0x215 */
    case WM_MOVING:         /* 0x216 */
    case WM_ENTERSIZEMOVE:  /* 0x231 */
    case WM_EXITSIZEMOVE:   /* 0x232 */
    case WM_IME_SETCONTEXT:  /* 0x281 */
    case WM_IME_NOTIFY:      /* 0x282 */
    case WM_NCMOUSELEAVE:    /* 0x2a2 */
    case WM_PRINTCLIENT:     /* 0x318 */
    case WM_DWMNCRENDERINGCHANGED:  /* 0x31f */
    case WM_USER:            /* 0x400 */
        return FALSE;

    default:
        printf_s("SetDialogAct: unknown message %x wpHi = %x, wpLo = %x"
            " additionalinfo2 = %llx \n", 
            message, wpHi, wpLo, additionalInfo2);
        return FALSE;

    case WM_INITDIALOG:    /* 0x110 */ {
        BOOL rv = TRUE;
        LRESULT rr;

        /* tooltip help message */
        static wchar_t grouphelp[] = L"help message \n";

        /* initialise combi box */
        HWND hwYafuPlan = GetDlgItem(DiBoxHandle, YafuPlan);
        for (int i = 0; i < planTextSize; i++)
            SendMessageA(hwYafuPlan, (UINT)CB_ADDSTRING, (WPARAM)0, 
                (LPARAM)planText[i]);
        /* highlight value currently selected */
        SendMessageA(hwYafuPlan, (UINT)CB_SETCURSEL, (WPARAM)pvalue-1, (LPARAM)0);

         /* set other controls to current values */       
        SetDlgItemInt(DiBoxHandle, verboseValue, verbose, FALSE);
        SetDlgItemInt(DiBoxHandle, post_eval_process, factorFlag, FALSE);
        SetDlgItemInt(DiBoxHandle, group_size_int, groupSize, FALSE);

        /* set tick boxes to current values */
        CheckDlgButton(DiBoxHandle, Msieve_E_option, eopt);
        CheckDlgButton(DiBoxHandle, hexPrint, hexPrFlag);
        CheckDlgButton(DiBoxHandle, sel_language, lang);

        if (yafu)
            rv = CheckRadioButton(DiBoxHandle, useYAFU, usepari, useYAFU);
        if (msieve)
            rv = CheckRadioButton(DiBoxHandle, useYAFU, usepari, useMsieve);
        if (Pari)
            rv = CheckRadioButton(DiBoxHandle, useYAFU, usepari, usepari);
        if (!yafu && !msieve && !Pari)
            rv = CheckRadioButton(DiBoxHandle, useYAFU, usepari, Builtin);
        if (rv == FALSE) {
            ErrorDisp(__FUNCTION__);
        }

        hh = CreateToolTip(group_size_int, DiBoxHandle, grouphelp, calcHandle);
        /* the amount of time the tooltip window remains visible */
        auto poptime = SendMessage(hh, TTM_GETDELAYTIME, TTDT_AUTOPOP, 0);
        /*  the amount of time the pointer must remain stationary within a tool's 
        bounding rectangle before the tooltip window appears.*/
        auto initTime = SendMessage(hh, TTM_GETDELAYTIME, TTDT_INITIAL, 0);
        /* the amount of time it takes for subsequent tooltip windows to appear 
        as the pointer moves from one tool to another*/
        auto rstime = SendMessage(hh, TTM_GETDELAYTIME, TTDT_RESHOW, 0);
        rr = SendMessage(hh, TTM_ACTIVATE, TRUE, 0);
        if (rr != TRUE)
            ErrorDisp(__FUNCTION__ " TTM ACTIVATE");

        /* return true so system sets focus to 1st control */
        return TRUE;
        }

    case WM_COMMAND:      /* 0x111 control selected by user */
        std::vector<std::string> p;
  
        switch (wpLo) {  /* switch according to control selected */
            int b_ck;    /* button status: checked, indeterminate or unchecked */
            LRESULT rr;

        case IDOK:       /* OK     */
        case IDCANCEL:   /* cancel */
            EndDialog(DiBoxHandle, wpLo);    /* close dialog box */
            return TRUE;

        case Builtin:  /* 'builtin' radio button. wpHi contains code (BN_CLICKED, etc)*/
            /* Dont't use YAFU, Msieve or Pari for factoring */
            yafu = false;
            msieve = false;
            Pari = false;
            return FALSE;

        case useYAFU:     /* yafu radio button. wpHi contains code (BN_CLICKED, etc)*/
            yafu = true;
            msieve = false;
            Pari = false;
            return FALSE;

        case useMsieve:   /* Msieve radio button. wpHi contains code (BN_CLICKED, etc)*/
            yafu = false;
            msieve = true;
            Pari = false;
            return FALSE;

        case usepari:  /* pari radio button. wpHi contains code (BN_CLICKED, etc)*/
            yafu = false;
            msieve = false;
            Pari = true;
            return FALSE;

        case IDC_BUTTON1:   /* check YAFU path */
            p.resize(2);
            p[0] = "YAFU";
            p[1] = "PATH";
            yafuParam(p);
            return FALSE;

        case setYAFUpath:
            p.resize(3);
            p[0] = "YAFU";
            p[1] = "PATH";
            p[2] = "SET";
            yafuParam(p);
            return FALSE;

        case YafuPlan:  /* change YAFU plan*/
            switch (wpHi) {
            case CBN_SELCHANGE: {
                int itemIndex = (int)SendMessageA(GetDlgItem(DiBoxHandle, YafuPlan),
                    (UINT)CB_GETCURSEL, (WPARAM)0, (LPARAM)0);
                pvalue = itemIndex + 1;
                break;
            }
            case CBN_DBLCLK:
            case CBN_SETFOCUS:
            case CBN_KILLFOCUS:
            case CBN_EDITCHANGE:
            case CBN_EDITUPDATE:
            case CBN_DROPDOWN:
            case CBN_CLOSEUP:
            case CBN_SELENDOK:
            case CBN_SELENDCANCEL:
                break;
            }
            return FALSE;

        case Yafulog:
            p.resize(2);
            p[0] = "YAFU";
            p[1] = "OUT";
            yafuParam(p);
            return FALSE;

        case YAFU_out_path:
            p.resize(3);
            p[0] = "YAFU";
            p[1] = "OUT";
            p[2] = "SET";
            yafuParam(p);
            return FALSE;

        case Yafu_GGNFS:
            p.resize(2);
            p[0] = "YAFU";
            p[1] = "INI";
            yafuParam(p);
            return FALSE;

        case YAFU_GGNFS_change:
            p.resize(3);
            p[0] = "YAFU";
            p[1] = "INI";
            p[2] = "I";
            yafuParam(p);
            return FALSE;

        case MSievePath:
            p.resize(2);
            p[0] = "MSIEVE";
            p[1] = "PATH";
            msieveParam(p);
            return FALSE;

        case Msieve_path_set:
            p.resize(3);
            p[0] = "MSIEVE";
            p[1] = "PATH";
            p[2] = "SET";
            msieveParam(p);
            return FALSE;

        case Msieve_E_option:   /* set/unset eopt*/
            b_ck = IsDlgButtonChecked(DiBoxHandle, Msieve_E_option);
            if (b_ck == BST_CHECKED)
                eopt = true;
            if (b_ck == BST_UNCHECKED)
                eopt = false;
            return FALSE;

        case Check_Msieve_Log:
            p.resize(2);
            p[0] = "MSIEVE";
            p[1] = "LOG";
            msieveParam(p);
            return FALSE;

        case Set_Msieve_Log:
            p.resize(3);
            p[0] = "MSIEVE";
            p[1] = "LOG";
            p[2] = "SET";
            msieveParam(p);
            return FALSE;

        case Pari_path:
            p.resize(2);
            p[0] = "PARI";
            p[1] = "PATH";
            pariParam(p);
            return FALSE;

        case Set_Pari_path:
            p.resize(3);
            p[0] = "PARI";
            p[1] = "PATH";
            p[2] = "SET";
            pariParam(p);
            return FALSE;

        case verboseValue:   /* set verbose */
            switch (wpHi) {
            case EN_SETFOCUS:
            case EN_KILLFOCUS:
            default:
                return FALSE;

            case EN_CHANGE:
            case EN_UPDATE:
                temp = GetDlgItemInt(DiBoxHandle, verboseValue, &good, FALSE);
                if (good == TRUE)
                    verbose = temp;
                return FALSE;
            }
            return FALSE;

        case post_eval_process:  /* set factorflag */
            switch (wpHi) {
            case EN_SETFOCUS:
            case EN_KILLFOCUS:
            default:
                return FALSE;

            case EN_CHANGE:
            case EN_UPDATE:
                temp = GetDlgItemInt(DiBoxHandle, post_eval_process, &good, FALSE);
                if (good == TRUE)
                    factorFlag = temp;
                return FALSE;
            }
            return FALSE;

        case hexPrint:       /* set/unset hex*/
            b_ck = IsDlgButtonChecked(DiBoxHandle, hexPrint);
            if (b_ck == BST_CHECKED)
                hexPrFlag = true;
            if (b_ck == BST_UNCHECKED)
                hexPrFlag = false;
            return FALSE;

        case sel_language:
            b_ck = IsDlgButtonChecked(DiBoxHandle, sel_language);
            if (b_ck == BST_CHECKED)
                lang = TRUE;    /* select Spanish */
            if (b_ck == BST_UNCHECKED)
                lang = FALSE;   /* select English */
            return FALSE;

        case group_size_int:
            switch (wpHi) {
            case EN_SETFOCUS:
                rr = SendMessage(hh, TTM_ACTIVATE, TRUE, 0);
                if (rr != TRUE)
                    ErrorDisp(__FUNCTION__ " TTM_ACTIVATE");
                /* TTM_POPUP requires Comclt32.dll version 6.0. */
                rr = SendMessage(hh, TTM_POPUP, 0, 0);
                if (rr != TRUE)
                    ErrorDisp(__FUNCTION__ " TTM_POPUP");
                return false;

            case EN_KILLFOCUS:
            default:
                return FALSE;

            case EN_CHANGE:
            case EN_UPDATE:
                temp = GetDlgItemInt(DiBoxHandle, group_size_int, &good, FALSE);
                if (good == TRUE)
                    groupSize = temp;
                return FALSE;
            }
            return FALSE;

        case help_file_path_check:
            std::cout << "path = " << helpFilePath << '\n';
            fileStatus(helpFilePath);
            return FALSE;

        case help_file_path_set:
            char* newpathC;
            newpathC = getFileName("Text\0*.TXT\0\0", handConsole);
            if (newpathC == NULL) {
                std::cout << "command cancelled \n";
            }
            else {
                helpFilePath = newpathC; /* copy new path for doc file to permanent storage */
                writeIni();       /* update the .ini file*/
                std::cout << "new path = " << helpFilePath << '\n';
                fileStatus(helpFilePath);
            }
            return FALSE;

        default:  /* unknown control*/
            std::cout << "SetDialog WM_COMMAND  wpHi = " << wpHi << " wpLo = " << wpLo
                << " info2 = " << additionalInfo2 << '\n';
        }
        return FALSE;
    }
    return FALSE;
}

/* set up dialog box for SET command. When the user select an action it is processed
in the SetDialogAct function. */
static long long setdiag(void) {
    auto rv = DialogBoxParamW(GetModuleHandle(nullptr), MAKEINTRESOURCE(Change_settings),
        handConsole, SetDialogAct, (LPARAM)99L);
    /* control is returned here when the dialog box is closed. */
    if (rv != IDOK && rv != IDCANCEL) {
        std::cout << "rv = " << rv << '\n';
        ErrorDisp(__FUNCTION__);
        return rv;
    }

    return rv;
}

/* check for commands. return 2 for exit, 1 for other command, 0 if not a valid command.
Note: command syntax checking is not rigorous; typing errors may produce unexpected 
results. command is a copy of the input buffer, not a reference. This is intentional
because strtok_s overwrites part of the buffer. */
static int processCmd(std::string command) {

    /* list of commands (static for efficiency) */
    const static std::vector<std::string> list =
    { "EXIT",  "SALIDA", "HELP", "AYUDA", "E",      "S",  
       "F",    "X",      "D",    "TEST",  "MSIEVE", "YAFU", 
       "V",    "PRINT",  "LIST", "LOOP",  "REPEAT", "IF",
       "PARI", "QMES",   "SET"};

    std::vector<std::string> p;   /* each parameter is stored separately in an element of p */
    int p1 = INT_MIN;             /* if 1st parameter is numeric, store its value here */
    const char seps[] = ", \n";   /* separators between parameters; either , or space */
    char* token = nullptr;
    char* next = nullptr;         /* used by strtok_s */

    if (command.empty()) return 0;  /* buffer is empty; not a command */

    /* separate params text into an array of tokens, by finding the separator characters */
    token = strtok_s(&command[0], seps, &next);
    while (token != nullptr) {
        p.push_back(token);
        token = strtok_s(nullptr, seps, &next);
    }
    if (p.size() >= 2 && isdigit(p[1][0]))
        p1 = std::atoi(p[1].data());  /* if p[1] is a decimal number set p1 to value */

    /* do 1st characters of text in command match anything in list? */
    int ix = 0;
    for (ix = 0; ix < list.size(); ix++) {
        if (p[0] == list[ix])
            break;
    }

    if (ix >= list.size())
        return 0;   /* not a command in list */

    switch (ix) {
    case 0:  /* exit*/
    case 1:  /* salida */
        return 2;  /* same command in 2 languages */
    case 2: /* Help */ {
            helpfunc(p);
            return 1;
        }
    case 3: /* ayuda */ {
            helpfunc(p);
            return 1;
        }
    case 4: /* E */ {
         lang = 0; return 1; }           // english
    case 5: /* S */ { 
        lang = 1; return 1; }	          // spanish (Español)
    case 6: /* F */ { 
         if (p.size() >= 2)
            factorFlag = p1;
         std::cout << (lang ? "factor establecido como " : "factor set to ") << factorFlag << '\n';
        return 1; }  
    case 7: /* X */ { 
        hexPrFlag = true; return 1; }         // hexadecimal output
    case 8: /* D */ { 
        hexPrFlag = false; return 1; }        // decimal output
    case 9: /* TEST */ {
            /* there are a number of commands that begin with TEST */
         if (p.size() == 1) {
             doTests();         // do basic tests 
             return 1;
         }

         char ttype = toupper(p[1][0]);  /* print help message? */
         if (ttype == 'H') {
             std::cout << "Test command format: \n"
                 << "      Test    (with no parameters) do basic tests of calculator and factorisation \n"
                 << "      Test 2 [p1[,p2[,p3]]]        test factorisation, where p1 = no of tests, \n"
                 << "                                   p2 = number size in bits \n"
                 << "      Test 3                       tests for builtin BigIntegers \n"
                 << "      Test 4                       test factorisation of mersenne numbers. \n"
                 << "      Test 5                       test using YAFU for factorisation \n"
                 << "      Test 6                       test using Msieve for factorisation \n"
                 << "      Test 7  [p1]                 tests the Lucas-Lehmer function for Mersenne \n"
                 << "                                   numbers up to 2^p1 -1 \n"
                 << "      Test 8                       test error handling \n"
                 << "      Test 9                       test modular square root (use \"TEST 9 H\" for more info.\n"
                 << "      Test 10 [p1[,p2]]            test quadratic modular equation solver \n"
                 << "                   where p1 is the number of tests and p2 is the number size in bits \n"
                 << "      Test 11 [p1[,p2[,p3]]]       test R3 & R3h functions, where p1 = no of tests, \n"
                 << "                                   p2 = number size in bits \n";

             return 1;
         }

         switch (p1)  {
         case 2: /* test factorisation */ {
             doTests2(p);
             return 1;
         }
         case 3:
    #ifdef BIGNBR
             {
                doTests3();         // do basic tests 
                return 1;
             }
    #else
             return 0;  /* Bignbr tests omitted */
    #endif
         case 4: /* factorise Mersenne numbers*/ {
             doTests4(); 
             return 1;
         }
         case 5: /* do YAFU tests*/ {
             doTests6(false);  
             return 1;
         }
         case 6: /* do Msieve tests*/ {
             doTests6(true);    
             return 1;
         }
         case 7: /* Lucas-Lehmer test */ {
             doTests7(p);
             return 1;
         }
         case 8: /* test error handling */ {
             testerrors();  
             return 1;
         }
         case 9: /* test modular square root*/ {  /* command format is:
                    TEST 9
                    or
                    TEST 9 new
                    or
                    TEST 9 time x y
                    where x is the size in bits of the numbers to test (default =10)
                          y is the number of tests (default = 5) */
            doTests9(p);   
            return 1;
         }
         case 10: /* test quadratic modular equation solver */ {
       
            doTestsA(p);   /* test quadratic modular equation solver  */
            return 1;
         }

         case 11:   /* R3 tests*/ {
             doTestsB(p);
             return 1;
         }

         default:
             return 0;  /* not a recognised command */
         }
        }
    case 10: /* MSIEVE */ {
            msieveParam(p);
            return 1;
        }
    case 11: /* YAFU */ {
            yafuParam(p);
            return 1;
        }
    case 12: /* V */ {
            if (p.size() >= 2)
                verbose = p1;
            std::cout << "verbose set to " << verbose << '\n';
            return 1;
        }
    case 13: /* PRINT */ {
            if (p.size() >= 2)
                printvars(p[1]);
            else
                printvars("");
            return 1;
        }
    case 14: /* LIST */ {
            if (exprList.empty())
                std::cout << "No expressions stored yet \n";
            else
                std::cout << exprList.size() << " expressions stored \n";
            for (auto exp : exprList) {
                std::cout << exp << '\n';
            }
            return 1;
        }
    case 15: /* LOOP */ {
            exprList.clear();
            return 1;
        }
    case 16: /* REPEAT */ {
            int repeat = 0;
            Znum result;
            retCode rv;
            if (p.size() >=2)
                repeat = p1;
            else
                repeat = 1;  /* repeat once if no count specified */

            loop1:
            for (int count = 1; count <= repeat; count++) {
                for (auto expr : exprList) {
                    if (expr.size() > 2 && expr.substr(0, 2) == "IF") {
                        int rv2 = ifCommand(expr);  /* analyse stored command again */
                        if (rv2 == 2)
                            goto loop1;  /* repeat all stored expressions again */
                        else if (rv2 == 0)
                            continue;  /* expr = 0; do not perform STOP or LOOP */
                        else if (rv2 == 1)
                            break;     /* expression is not 0 and STOP specified */
                        else if (rv2 == 3)
                            continue;   /* THEN or ELSE expression has been evaluated */
                        else
                            throw std::logic_error("unknown return code");
                        abort();  /* where are we? */
                    }
                    /* recalculate each stored expression */
                    rv = ComputeMultiExpr(expr, result);
                }
            }
            return 1;
        }
    case 17: /* IF */ {
            /*format is IF (expression) REPEAT or IF (expression) STOP 
            or IF (expression) THEN (expression[, expression ...])
                                 ELSE (expression[, expression ...]) */
            Znum result;
            int rv = ifCommand(command);  /* analyse IF command */
            if (rv == 2) {
                /* REPEAT */
                exprList.push_back(command);
                loop:
                for (auto expr : exprList) {
                    if (expr.size() > 2 && expr.substr(0, 2) == "IF") {
                        int rv2 = ifCommand(expr);  /* analyse stored command again */
                        if (rv2 == 2)
                            goto loop;  /* repeat all stored expressions again */
                        else if (rv2 == 0)
                            continue;  /* expr = 0; do not perform STOP or LOOP */
                        else if (rv2 == 1)
                            break;     /* expression is not 0 and STOP specified */
                        else if (rv2 == 3)
                            continue;   /* THEN or ELSE expression has been evaluated */
                        else
                            throw std::logic_error ("unknown return code");
                            abort();  /* where are we? */
                    }

                    /* recalculate each stored expression */
                    retCode rv = ComputeMultiExpr(expr, result);
                }
                return 1;  /* loop completed */
            }
            if (rv == 0 || rv == 3)
                /* if command processed, now save it for possible loop */
                exprList.push_back(command);  
            /* to be completed for rv =-1, rv = 1*/
            return 1;
        }
    case 18: /* PARI */ {
            pariParam(p);
            return 1;
        }
    case 19: /* QMES */ {
            quadModEqn(p);  /* Quadratic Modular Equation Solver */
            return 1;
        }
    case 20:  /* SET */
        {
            setdiag();   /* change settings using a dialog box */
            return 1;
        }
    default:
        return 0;   /* not a recognised command */
    }
}

/* initialisation code, executed once at start of program execution */
static void initialise(int argc, char *argv[]) {
    char VSversion[100] = { 0 };  /* visual studio version*/
    size_t vslen = 0;             /* no of chars in visual studio version */
    flags f = { 0,0,0, 0,0,0, 0,0,0, 0,0,0 };
#ifndef BIGNBR
    unsigned int control_word; /* save current state of fp control word here */
#endif
    int version[4]; /* version info from .exe file (taken from .rc resource file) */
    std::string modified;  /* date & time program last modified */

    f.UnEx = 1;        /* set unhandled exception 'filter' */
    f.abort = 1;       /* trap abort (as a signal) */
    f.sigterm = 1;     /* trap terminate signal */
    f.sigill = 1;      /* trap 'illegal' signal */
    f.interrupt = 1;   /* trap interrupt signal */
    //f.sigsegv = 1;     /* trap segment violation signal */
#ifdef _DEBUG
    /* only seems to work properly if compiled in debug mode */
    f.InvParam = 1;    /* trap invalid parameters on library calls */
#endif
    //f.sigfpe = 1;      /* trap floating point error signal */
    SetProcessExceptionHandlers(f);

    /* if we trap floating point errors we trap  _EM_INVALID in mpir prime test
        functions that actually work OK */
#ifndef BIGNBR
    err = _controlfp_s(&control_word, _EM_INEXACT | _EM_UNDERFLOW, MCW_EM);
    /* trap hardware FP exceptions except inexact and underflow which are
    considered to be normal, not errors. */
    if (err) {
        printf_s("could not set FP control word\n");
        exit (-1);
    }
#endif

    hConsole = GetStdHandle(STD_OUTPUT_HANDLE);  // get handle for stdout
    handConsole = GetConsoleWindow();            // get handle for console window
    if (hConsole == INVALID_HANDLE_VALUE)
    {
        fprintf_s(stderr, "GetStdHandle failed with %d at line %d\n", GetLastError(), __LINE__);
        Beep(750, 1000);
        exit (EXIT_FAILURE);
    }

    VersionInfo(argv[0], version, modified); /* get version info from .exe file */
    printf_s(lang? "Bigint calculadora versão %d.%d.%d.%d \n" : 
           "%s Bigint calculator Version %d.%d.%d.%d \n", 
        myTime(), version[0], version[1], version[2], version[3]);
    std::cout << (lang? "última modificação em " : "last modified on ") << modified << '\n';

#ifdef __GNUC__
    printf("gcc version: %d.%d.%d\n", __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
    setlocale(LC_ALL, "en_GB.utf8");      // allows non-ascii characters to print
#endif

#ifdef _MSC_FULL_VER
/* For example, if the version number of the Visual C++ compiler is 15.00.20706.01,
the _MSC_FULL_VER macro evaluates to 150020706 */
    long long ver = _MSC_FULL_VER;
    std::cout << "MSVC version: " << ver / 10000000 << '.';  // 1st 2 digits
    ver %= 10000000;                      // remove 1st 2 digits
    std::cout << ver / 100000 << '.';    // next 2 digits
    ver %= 100000;                        // remove next 2 digits
    std::cout << ver << '\n';             // last 5 digits

    auto lc = setlocale(LC_ALL, "en-EN.utf8");      // allows non-ascii characters to print
#endif

    printf_s("locale is now: %s\n", setlocale(LC_ALL, NULL));
    std::cout << "GMP version: " << __GNU_MP_VERSION << '.' << __GNU_MP_VERSION_MINOR
        << '.' << __GNU_MP_VERSION_PATCHLEVEL << '\n';

#ifdef __MPIR_VERSION
    std::cout << "MPIR version: " << __MPIR_VERSION << '.' << __MPIR_VERSION_MINOR
        << '.' << __MPIR_VERSION_PATCHLEVEL << '\n';
#endif
    std::cout << "Boost version: " << BOOST_VERSION << '\n';

    unsigned long long comCtlVer = getComCtlVer();  /* get version of ComCtl32.dll */

    getenv_s(&vslen, VSversion, "VisualStudioEdition");
    if (vslen != 0) {
        VSrun = true;      /* program started from visual studio */
        if (verbose > 0) {
            std::cout << "Visual Studio version: " << VSversion << '\n';
        }
    }
    processIni(argv[0]); // read .ini file if it exists

    //get the computer name, cache sizes, etc.  store in globals
    get_computer_info(CPU_ID_STR, MEAS_CPU_FREQUENCY);

    printf_s("detected %s\n", CPU_ID_STR);
    printf_s("measured cpu frequency ~= %.0f\n", 	MEAS_CPU_FREQUENCY);

    INITCOMMONCONTROLSEX ccset = { 0, ICC_WIN95_CLASSES };
    ccset.dwSize = sizeof(ccset);

    BOOL iccRV = InitCommonControlsEx(&ccset);
    if (iccRV == FALSE)
        ErrorDisp(__FUNCTION__);
    return;
}

/* get input from stdin. Any continuation lines are appended to 1st line.
Initial & trailing spaces are removed. Letters are converted to upper case. 
ctrl-c or ctrl-break will force the function to return, with or without input, 
but only after a 5 sec delay.*/
static void myGetline(std::string &expr) {

    std::getline(std::cin, expr);    // expression may include spaces
    if (std::cin.fail() || std::cin.bad()) {
        expr.erase();
        return;   /* error reading from stdin */
    }

    if (breakSignal) {
        Sleep(5000);   /* wait 5 seconds */
        return;     /* Program interrupted: ctrl-c or ctrl-break */
    }
    strToUpper(expr, expr);		// convert to UPPER CASE 
    removeInitTrail(expr);       // remove initial & trailing spaces

    if (expr.empty()) {
        return;
    }

    /* check for continuation character. If found get continuation line(s) */
    while (expr.back() == '\\') {   /* ends with continuation character? */
        std::string cont;
        std::cout << (lang ? "continuar: " : "continue: ");
        std::getline(std::cin, cont);   /* get continuation line */
        if (std::cin.fail() || std::cin.bad()) {
            expr.erase();
            return;
        }
        strToUpper(cont, cont);   // convert to UPPER CASE 
        while (!cont.empty() && isspace(cont.back())) {
            cont.resize(cont.size() - 1);   /* remove trailing space */
        }
        expr.resize(expr.size() - 1); /* remove trailing \ char */
        expr += cont;    /* append continuation line to previous input */
    }

    if (breakSignal) {
        Sleep(5000);   /* wait 5 seconds */
        return;     /* Program interrupted: ctrl-c or ctrl-break */
    }

}

int main(int argc, char *argv[]) {
    std::string expr;
    Znum Result;
    retCode rv;
    int asgCt;  /* number of assignment operators */
    bool multiV = false;

    try {

        initialise(argc, argv);  /* initialisation code only executed once */

        /* start of main loop. Normal exit is via EXIT command */
        while (true) {
            if (lang == 0) {
                printf_s("enter expression to be processed, or HELP, EXIT, SET or QMES\n");
            }
            else
                printf_s("ingrese la expresión para ser procesada, o AYUDA o SALIDA o SET o QMES\n");

            PlaySoundA(attsound.c_str(), NULL,
                SND_FILENAME | SND_NODEFAULT | SND_ASYNC | SND_NOSTOP);

            myGetline(expr);  /* get input from stdin */
            if (breakSignal)
                break;    /* Program interrupted: ctrl-c or ctrl-break */

            if (expr.empty()) {
                Sleep(1000);       
                continue;            /* no input */
            }

            // prevent the sleep idle time-out.
            SetThreadExecutionState(ES_CONTINUOUS | ES_SYSTEM_REQUIRED);

            int cmdCode = processCmd(expr);  /* is input a command? */
            if (breakSignal)
                break;    /* Program interrupted: ctrl-c or ctrl-break */
            if (cmdCode == 2) 
                break;    // EXIT command
            if (cmdCode == 1) {
                // command has been fully processed; go back to start of loop
                // Clear EXECUTION_STATE flags to allow the system to idle to sleep normally.
                SetThreadExecutionState(ES_CONTINUOUS);
                continue;
            }
            if (cmdCode != 0) {
                fprintf_s(stderr, "Invalid return code %x from processCmd \n", cmdCode);
                continue;    // weird return code; can't process further; 
                             //  go back to start of loop
            }

            /* input is not a valid command; assume it is an expression */
            auto start = clock();	// used to measure execution time
            removeIntSpace(expr);   /* remove spaces between digits */
            rv = ComputeExpr(expr, Result, asgCt, &multiV); /* analyse expression, compute value*/

            if (rv != retCode::EXPR_OK) {
                textError(rv);   // invalid expression; print error message
            }
            else {
                exprList.push_back(expr);  /* save text of expression */
                if (!multiV) {
                    std::cout << "result = ";
                    ShowLargeNumber(Result, groupSize, true, hexPrFlag);   // print value of expression
                    std::cout << '\n';
                    if (factorFlag > 0) {
                        doFactors(Result, false); /* factorise Result, calculate number of divisors etc */
                        results.clear();  // get rid of unwanted results
                    }
                }
                else {
                    /* expression returned multiple values */
                    std::cout << " = ";
                    if (roots.size() <= 31)
                        /* print all results if <= 31 values */
                        for (auto r : roots) {
                            std::cout << r << ", ";
                        }
                    else { /* print 1st 20 and last 10 results */
                        for (int i = 0; i < 20; i++) {
                            std::cout << roots[i] << ", ";
                        }
                        std::cout << "\n ... \n";
                        for (size_t i = roots.size() - 10; i < roots.size(); i++) {
                            std::cout << roots[i] << ", ";
                        }
                    }
                    putchar('\n');
                    std::cout << "found " << roots.size() << " results \n";
                }
                if (asgCt != 0)
                    printvars(""); /* print variables names & values */
            }

            auto end = clock();   // measure amount of time used
            double elapsed = (double)end - start;
            PrintTimeUsed(elapsed, lang? "Tiempo transcurrido = " : "time used = ");
            // Clear EXECUTION_STATE flags to allow the system to idle to sleep normally.
            SetThreadExecutionState(ES_CONTINUOUS);

            if (breakSignal)
                break;     /* Program interrupted */

            /* now go back to start of loop */
        } /* end of while loop */

        /* get to here when we break out of main loop, usually by EXIT command */

        // Clear EXECUTION_STATE flags to allow the system to idle to sleep normally.
        SetThreadExecutionState(ES_CONTINUOUS);
        if (lang)
            std::cout << "¡Adiós!";
        else
            std::cout << "Goodbye \n";
        if (!VSrun)
            system("PAUSE");  /* "press any key to continue ..." (unless program 
                             started from Visual Studio) */
        return EXIT_SUCCESS;  // EXIT command entered
    }

#undef max  /* remove max defined in windows.h  because of name clash */

    /* code below catches C++ 'throw' type exceptions */
    catch (const std::exception& e) {

        // Clear EXECUTION_STATE flags to allow the system to idle to sleep normally.
        SetThreadExecutionState(ES_CONTINUOUS);

        printf_s("\n*** standard exception caught, message '%s'\n", e.what());
        Beep(1200, 1000);              // sound at 1200 Hz for 1 second
        std::cout << "Press ENTER to continue...";
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        exit(EXIT_FAILURE);
    }

    catch (const char *str) {

        // Clear EXECUTION_STATE flags to allow the system to idle to sleep normally.
        SetThreadExecutionState(ES_CONTINUOUS);

        printf_s("\n*** Caught exception: <%s> \n", str);
        Beep(1200, 1000);              // sound at 1200 Hz for 1 second
        std::cout << "Press ENTER to continue...";
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        exit(EXIT_FAILURE);
    }

    catch (int e) {

        // Clear EXECUTION_STATE flags to allow the system to idle to sleep normally.
        SetThreadExecutionState(ES_CONTINUOUS);

        printf_s("\n*** Caught exception: <%d> \n", e);
        Beep(1200, 1000);              // sound at 1200 Hz for 1 second
        std::cout << "Press ENTER to continue...";
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        exit(EXIT_FAILURE);

    }

    catch (...) {
        // this executes if f() throws any other unrelated type
        // This catch block probably only would be executed under /EHa compiler option 
        /* most likely to be a SEH-type exception */

        // Clear EXECUTION_STATE flags to allow the system to idle to sleep normally.
        SetThreadExecutionState(ES_CONTINUOUS);

        printf_s("\n*** unknown exception ocurred\n");
        Beep(1200, 1000);              // sound at 1200 Hz for 1 second
        std::cout << "Press ENTER to continue...";
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        exit(EXIT_FAILURE);
    }
}
