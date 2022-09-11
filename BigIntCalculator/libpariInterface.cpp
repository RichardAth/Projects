#include "pch.h"
#include <windows.h>
#include "factor.h"
//#include "pari.h"
extern std::string PariPath;

/* stuff below copied from pari headers, to avoid including humungous pari.h etc*/
typedef long long* GEN;
typedef unsigned long long ulong;

#define BITS_IN_LONG 64
#define TYPnumBITS   7
#define SIGNnumBITS  2
#define LGnumBITS (BITS_IN_LONG - 1 - TYPnumBITS)
#define TYPSHIFT (BITS_IN_LONG - TYPnumBITS)
#define SIGNSHIFT (BITS_IN_LONG - SIGNnumBITS)
#define EXPOnumBITS (BITS_IN_LONG - SIGNnumBITS)
#define HIGHEXPOBIT (1ULL<<(EXPOnumBITS-1))
#define EXPOBITS    ((1ULL<<EXPOnumBITS)-1)
#define LGBITS      ((1ULL<<LGnumBITS)-1)

#define lg(x)         ((int64_t)(((ulong)((x)[0])) & LGBITS))
#define typ(x)        ((int64_t)(((ulong)((x)[0])) >> TYPSHIFT))
#define signe(x)      (((int64_t)((x)[1])) >> SIGNSHIFT)
#define lgefint(x)      ((int64_t)(((ulong)((x)[1])) & LGBITS))
#define gmael1(m,x1)             (((GEN*)    (m))[x1])
#define gel(m,x)     gmael1(m,x)
#define expo(x)       ((int64_t) ((((ulong)((x)[1])) & EXPOBITS) - HIGHEXPOBIT))
#define realprec(x)   ((int64_t)(((ulong)((x)[0])) & LGBITS))
#define int_MSW(x) ((x)+lgefint((x))-1)
#define int_precW(x) ((x)-1)

enum {
    t_INT = 1,
    t_REAL = 2,
    t_INTMOD = 3,
    t_FRAC = 4,
    t_FFELT = 5,
    t_COMPLEX = 6,
    t_PADIC = 7,
    t_QUAD = 8,
    t_POLMOD = 9,
    t_POL = 10,
    t_SER = 11,
    t_RFRAC = 13,
    t_QFR = 15,
    t_QFI = 16,
    t_VEC = 17,
    t_COL = 18,
    t_MAT = 19,
    t_LIST = 20,
    t_STR = 21,
    t_VECSMALL = 22,
    t_CLOSURE = 23,
    t_ERROR = 24,
    t_INFINITY = 25
};
/* end of stuff from pari headers */

#define ZT(a) a.backend().data()  /* access mpz_t within a Znum (Boost mpz_int)*/

bool fRunTimeLinkSuccess = false;    /* true when links to parilib set up */
bool stackinit = false;              /* true when pari stack is initialised */
bool Pari = false;                   /* true if parilib used for factorisation */

/* typedefs for access to libpari functions */
typedef      GEN(__cdecl* G_GG)     (GEN x, GEN y);  /* used for all funtions of the form GEN f(GEN x, GEN y) */
typedef      GEN(__cdecl* G_G)      (GEN x);         /* used for all funtions of the form GEN f(GEN x) */
typedef     void(__cdecl* v_v) (void);               /* used for functions of the form void f(void) */
typedef     void(__cdecl* pari_initX)(size_t, ulong);    /*void    pari_init(size_t parisize, ulong maxprime);*/
//typedef     void(__cdecl* pari_stack_initX)(pari_stack* s, size_t size, void** data);
//typedef  int64_t(__cdecl* pari_stack_newX)(pari_stack* s);
typedef     void(__cdecl* paristack_setsizeX)(size_t, size_t); /* void paristack_setsize(size_t rsize, size_t vsize) */
typedef     void(__cdecl* pari_init_optsX)(size_t parisize, ulong maxprime, ulong init_opts);
typedef     void(__cdecl* pari_init_primesX)(ulong maxprime);
typedef    ulong(__cdecl* hclassno6uX)(ulong D);   /* ulong hclassno6u(ulong D)*/
typedef      GEN(__cdecl* dvmdiiX)  (GEN x, GEN y, GEN* z);
typedef      int(__cdecl* invmodX)  (GEN a, GEN b, GEN* res);  /* int  invmod(GEN a, GEN b, GEN *res); */
typedef      GEN(__cdecl* itorX)    (GEN x, int64_t prec);     /* GEN    itor(GEN x, int64_t prec); */
typedef   double(__cdecl* rtodblX)  (GEN x);      /* double  rtodbl(GEN x); */
typedef      GEN(__cdecl* dbltorX)  (double x);   /* GEN     dbltor(double x); */
typedef      GEN(__cdecl* divruX)   (GEN x, ulong y);

typedef      GEN(__cdecl* expX)     (GEN x, int64_t prec);  /* GEN gexp(GEN x, int64_t prec)*/
typedef      GEN(__cdecl* logX)     (GEN x, int64_t prec);  /* GEN glog(GEN x, int64_t prec)*/
typedef   char* (__cdecl* GENtostrX)(GEN x);        /* char*   GENtostr(GEN x);*/
typedef      GEN(__cdecl* utoiposX) (ulong x);     /* GEN utoipos(ulong x)*/
typedef      GEN(__cdecl* utoinegX) (ulong x);
typedef      GEN(__cdecl* cgetiposX)(int64_t x); /* GEN    cgetipos(int64_t x);*/
typedef      GEN(__cdecl* cgetinegX)(int64_t x);
typedef      GEN(__cdecl* qfbclassnox)(GEN x, int64_t flag);  /* GEN     qfbclassno0(GEN x,int64_t flag);*/
typedef ulong** avmaX;
typedef     void(__cdecl* set_avmaX)(ulong av);  /* void   set_avma(ulong av);*/
typedef      GEN(__cdecl* stirlingX)(int64_t n, int64_t m, int64_t flag); /* GEN stirling(int64_t n, int64_t m, int64_t flag)*/

HINSTANCE hinstLib;

/* function pointers to access libpari functions */
static pari_initX        pari_init_ref;
//static pari_stack_initX  pari_stack_init_ref;
//static pari_stack_newX   pari_stack_new_ref;
static pari_init_optsX   pari_init_opts_ref;
static pari_init_primesX pari_init_primes_ref;
static G_GG              addii_ref;     /* a + b */
static G_GG              subii_ref;     /* a - b */
static G_GG              mulii_ref;     /* a * b */
static dvmdiiX           dvmdii_ref;
static hclassno6uX       hclassno6u_ref;
static G_G               hclassno_ref;
static invmodX           invmod_ref;
static itorX             itor_ref;
static rtodblX           rtodbl_ref;    /* floating point to double */
static dbltorX           dbltor_ref;
static G_GG              addrr_ref;     /* a + b */
static G_GG              subrr_ref;     /* a - b */
static G_GG              divrr_ref;     /* a / b  */
static G_GG              mulrr_ref;     /* a * b */
static divruX            divru_ref;    
static v_v               pari_print_version_ref;
//parimainstackX    parimainstack_ref;
static expX              exp_ref;       /*  exp(a) */
static logX              log_ref;       /* log(a)  */
static GENtostrX         GENtostr_ref;
static utoiposX          utoipos_ref;
static utoinegX          utoineg_ref;
static G_G               floorr_ref;     /* floor(x) */
static cgetiposX         cgetipos_ref;
static cgetinegX         cgetineg_ref;
static qfbclassnox       qfbclassno_ref;
static avmaX             avma_ref;
static set_avmaX         set_avma_ref;
static G_G               tau_ref;
static stirlingX         stirling_ref;
static paristack_setsizeX paristack_setsize_ref;
static G_G               factor_ref;
static v_v               pari_close_ref;

/* set up pointers to parilib functions, initialise pari stack.
pari library functions are accessed in this way because linking statically to
libpari.dll doesn't work, presumably because it was compiled with Msys2/gcc and
we are using windows + Visual Studio.
for every function we want to use, we have to set up a function pointer.
For some reason dynamic linking works OK. */
static void specinit() {

    if (!fRunTimeLinkSuccess) {
        // Get a handle to the pari DLL module.
        //hinstLib = LoadLibraryA("C:/Program Files (x86)/Pari64-2-13-2/libpari.dll");
        hinstLib = LoadLibraryA(PariPath.c_str());

        // If the handle is valid, try to get the function addresses.
        if (hinstLib == nullptr) {
            std::cerr << "could not access dll for PARI \n";
            system("PAUSE");
            abort();
        }

        else {
            /* set up function pointers. There over 6000 accessible libpari functions. This is
            a selection of functions that might be useful. */
            pari_init_ref = (pari_initX)GetProcAddress(hinstLib, "pari_init");
            paristack_setsize_ref = (paristack_setsizeX)GetProcAddress(hinstLib, "paristack_setsize");
            /* pari_stack_init_ref  = (pari_stack_initX)GetProcAddress(hinstLib, "pari_stack_init");
             pari_stack_new_ref   = (pari_stack_newX)GetProcAddress(hinstLib, "pari_stack_new");*/
            pari_init_opts_ref = (pari_init_optsX)GetProcAddress(hinstLib, "pari_init_opts");
            pari_init_primes_ref = (pari_init_primesX)GetProcAddress(hinstLib, "pari_init_primes");
            addii_ref = (G_GG)GetProcAddress(hinstLib, "addii");
            subii_ref = (G_GG)GetProcAddress(hinstLib, "subii");
            mulii_ref = (G_GG)GetProcAddress(hinstLib, "mulii");
            dvmdii_ref = (dvmdiiX)GetProcAddress(hinstLib, "dvmdii");
            hclassno6u_ref = (hclassno6uX)GetProcAddress(hinstLib, "hclassno6u");
            hclassno_ref = (G_G)GetProcAddress(hinstLib, "hclassno");
            invmod_ref = (invmodX)GetProcAddress(hinstLib, "invmod");
            itor_ref = (itorX)GetProcAddress(hinstLib, "itor");
            rtodbl_ref = (rtodblX)GetProcAddress(hinstLib, "rtodbl");
            dbltor_ref = (dbltorX)GetProcAddress(hinstLib, "dbltor");
            addrr_ref = (G_GG)GetProcAddress(hinstLib, "addrr");
            subrr_ref = (G_GG)GetProcAddress(hinstLib, "subrr");
            divrr_ref = (G_GG)GetProcAddress(hinstLib, "divrr");
            mulrr_ref = (G_GG)GetProcAddress(hinstLib, "mulrr");
            divru_ref = (divruX)GetProcAddress(hinstLib, "divru");
            pari_print_version_ref = (v_v)GetProcAddress(hinstLib, "pari_print_version");
            exp_ref = (expX)GetProcAddress(hinstLib, "gexp");
            log_ref = (logX)GetProcAddress(hinstLib, "glog");
            utoipos_ref = (utoiposX)GetProcAddress(hinstLib, "utoipos");
            utoineg_ref = (utoinegX)GetProcAddress(hinstLib, "utoineg");
            GENtostr_ref = (GENtostrX)GetProcAddress(hinstLib, "GENtostr");
            floorr_ref = (G_G)GetProcAddress(hinstLib, "floorr");
            cgetipos_ref = (cgetiposX)GetProcAddress(hinstLib, "cgetipos");
            cgetineg_ref = (cgetinegX)GetProcAddress(hinstLib, "cgetineg");
            qfbclassno_ref = (qfbclassnox)GetProcAddress(hinstLib, "qfbclassno0");
            avma_ref = (avmaX)GetProcAddress(hinstLib, "avma");
            set_avma_ref = (set_avmaX)GetProcAddress(hinstLib, "set_avma");
            tau_ref = (G_G)GetProcAddress(hinstLib, "ramanujantau");
            stirling_ref = (stirlingX)GetProcAddress(hinstLib, "stirling");
            factor_ref = (G_G)GetProcAddress(hinstLib, "factor");
            pari_close_ref = (v_v)GetProcAddress(hinstLib, "pari_close");

            /* check that all function pointers were set up successfully */
            if (nullptr == pari_init_ref || nullptr == paristack_setsize_ref ||
                /*           nullptr == pari_stack_init_ref ||
                           nullptr == pari_stack_new_ref ||*/
                nullptr == pari_init_opts_ref || nullptr == pari_init_primes_ref ||
                nullptr == addii_ref || nullptr == subii_ref ||
                nullptr == mulii_ref || nullptr == dvmdii_ref ||
                nullptr == invmod_ref || nullptr == itor_ref ||
                nullptr == rtodbl_ref || nullptr == dbltor_ref ||
                nullptr == addrr_ref || nullptr == subrr_ref ||
                nullptr == divrr_ref || nullptr == divru_ref ||
                nullptr == pari_print_version_ref ||
                nullptr == hclassno_ref || nullptr == exp_ref ||
                nullptr == log_ref || nullptr == utoipos_ref ||
                nullptr == utoineg_ref || nullptr == GENtostr_ref ||
                nullptr == floorr_ref || nullptr == cgetipos_ref ||
                nullptr == cgetineg_ref || nullptr == hclassno6u_ref ||
                nullptr == qfbclassno_ref || nullptr == avma_ref ||
                nullptr == set_avma_ref || nullptr == tau_ref ||
                nullptr == stirling_ref || nullptr == factor_ref ||
                nullptr == pari_close_ref) {
                fRunTimeLinkSuccess = false;
                std::cerr << "PARI dynamic linking failed \n";
                system("PAUSE");
                abort();
            }

            fRunTimeLinkSuccess = true;
        }
    }

    if (!stackinit) {
        pari_init_ref(8000000, 500000);  /* stack size 8 MB, maxprime */
        paristack_setsize_ref(8'000'000, 80'000'000);  /* change parisizemax to 80 MB */
        stackinit = true;
        if (verbose > 0) {
            pari_print_version_ref();  /* print some info from libpari dll */
            printf("Pari stack initialised. avma = %p\n", *avma_ref);
        }
    }
}

/* convert t_INT in x to mpz_t in value */
static void InttoMP(const GEN x, mpz_t value) {
#ifdef _DEBUG
    assert(typ(x) == t_INT);
#endif
    int64_t s = signe(x);
    int64_t i, lx = lgefint(x);   /* lx = effective length of x */

    if (s == 0) {
        mpz_set_ui(value, 0);   /* x = 0 */
        return;
    }

    mpz_set_ui(value, 0);   /* value = 0 */
    GEN y = int_MSW(x);  /* y points to most significant word of x */
    for (i = 2; i < lx; i++, y = int_precW(y)) {
        mpz_mul_2exp(value, value, 64); /* move more significant bits left 64 bits*/
        mpz_add_ui(value, value, *y);   /* add next 64-bit limb */
    }
    if (s < 0)
        mpz_neg(value, value);  /* if x  is -ve, flip sign */
}

/* convert GEN to MPIR/GMP multi-precision
for type t_INT (integer) value is set to the value, denom is set to 1,
    and val_d is set to the value as a double.
for type t_FRAC (rational number) value is set to the numerator, denom is set to
    the denominator and val_d is value/denom as a floating point.
for type t_REAL val_d is set to the value as a floating point, value is set to
    the integer part of val_d, and denom is set to 1.
 Note: 1. overflow of val_d is possible.
       2. for any other GEN type abort() is called.
       3. value and denom must be initialised before GENtoMP is called.*/
static void GENtoMP(const GEN x, mpz_t value, mpz_t denom, double& val_d) {
    int64_t typx = typ(x); /* get type of object in x */

    switch (typx) {
    case t_INT: {
        int64_t s = signe(x);
        int64_t lx = lgefint(x);   /* lx = effective length of x */

        mpz_set_ui(denom, 1);
        if (s == 0) {
            mpz_set_ui(value, 0);   /* x = 0 */
            val_d = 0.0;
            return;
        }
        InttoMP(x, value);
        val_d = mpz_get_d(value);  /* convert to floating point */
        return;
    }
    case t_FRAC: {
        GEN num = gel(x, 1);      /* numerator */
        GEN den = gel(x, 2);     /* denominator */
        int64_t s = signe(num);

        if (s == 0) {
            mpz_set_ui(value, 0);  /* x = 0 */
            mpz_set_ui(denom, 1);
            val_d = 0;
            return;
        }
        InttoMP(num, value);

        s = signe(den);
        if (s == 0) {
            mpz_set_ui(denom, 0);  /* this would be weird! Mathematically the value is infinite!*/
            return;
        }
        InttoMP(den, denom);
        val_d = mpz_get_d(value) / mpz_get_d(denom);
        return;
    }
    case t_REAL: {
        GEN xi = floorr_ref(x);   /* xi = integer part of x*/
        mpz_set_ui(denom, 1);
        InttoMP(xi, value);/* convert (GEN)xi to (mpz_t)value */
        val_d = rtodbl_ref(x);  /* also convert x to normal floating point */
        return;
    }

    default:
        abort();  /* unsupported GEN type */
    }
}

/* convert mpz_t to GEN */
GEN MPtoGEN(const mpz_t num) {
    ptrdiff_t numlimbs = mpz_size(num);
    GEN rv;
    if (numlimbs == 0) {
        rv = utoipos_ref(0);  /* num = 0 */
        return rv;
    }
    if (mpz_sgn(num) < 0)
        rv = cgetineg_ref(numlimbs + 2);
    else
        rv = cgetipos_ref(numlimbs + 2);

    for (int i = 0; i < numlimbs; i++)
        rv[i + 2] = num->_mp_d[i];  /* copy limbs from num to rv */
    return rv;
}

/* calculate R3 using Hurwitz class number. Let libpari do the heavy lifting. 
Note: n is intentionally a copy of the original value. This copy is modified. */
Znum R3h(Znum n) {
    mpz_t num, denom, result;
    double hxd;
    GEN x, hx;
    Znum r;

    if (n == 0)
        return 1;   /* easiest to make n = 0 a special case*/

    specinit();   /* initialise as required*/
    ulong* av = *avma_ref;  /* save address of available memory */
    mpz_inits(num, denom, result, nullptr);

    /* note: r3(4n) = r3(n) */
    while ((n & 3) == 0)
        n >>= 2;  /* if n is a multiple of 4, divide n by 4*/

    int mod = (int)MulPrToLong(n & 7);

    /* note that in general h(n) may not be an integer. hclassno returns num and
    denom and the ratio is the class number */
    switch (mod) {
    case 1:  /* n%4 = 1 or 2 */
    case 2:
    case 5:
    case 6: {
        /* get 12 * h(4n) */
        Znum n4 = n * 4;
        x = MPtoGEN(ZT(n4));
        hx = hclassno_ref(x);
        GENtoMP(hx, num, denom, hxd);    /* convert hx from GEN to mpz_t */
        mpz_mul_ui(result, num, 12);
        mpz_div(result, result, denom);  /* get 12*h(n) */
        r = result;          /* convert to Znum */
        break;
    }

    case 3:  /* n%8 = 3 */
        /* get 24 * h(n) */ {
        x = MPtoGEN(ZT(n));
        hx = hclassno_ref(x);
        GENtoMP(hx, num, denom, hxd);
        mpz_mul_ui(result, num, 24);     /* result = num * 24 */
        mpz_div(result, result, denom);  /* result =(num * 24)/denom  = 24*h(n)  */
        r = result;         /* convert to Znum */
        break;
    }
    case 7:  /* n%8 = 7 */
        r = 0;
        break;

    default:   /* n%8 = 0 or 4 */
        abort();   /* should never happen */
    }

    mpz_clears(num, denom, result, nullptr);  /* avoid memory leakage */
    ptrdiff_t diff = av - *avma_ref;
    if (verbose > 1)
        printf("used %lld bytes on pari stack \n", (long long)diff);
    set_avma_ref((ulong)av);      /* recover memory used */
    return r;            /* return result */
}

/* convert a T_Real to an mpf_t. The precison of the mpf_t is set to match the T_Real.
value must be initialised before calling TrealToMP. */
void TrealToMP(const GEN x, mpf_t value) {
    int64_t typx = typ(x);
    int64_t sign, exp, prec, i;
    if (typx != t_REAL) {
        std::cerr << "invalid conversion attempted; T_Real expected \n";
        return;
    }
    sign = signe(x);
    exp = expo(x);
    prec = realprec(x);  /* this is really the total length of x in 64-bit words */
    assert(prec >= 3);
    mpf_set_prec(value, (prec - 2) * 64);  /* set precision (convert length in
                                           words to length in bits) */
    mpf_set_ui(value, 0);
    for (i = 2; i < prec; i++) {    /* copy mantissa from x to value */
        mpf_mul_2exp(value, value, 64);
        mpf_add_ui(value, value, x[i]);
    }

    /* set exponent */
    int64_t shift = (prec - 2) * 64 - exp - 1;
    if (shift >= 0)
        mpf_div_2exp(value, value, shift);
    else
        mpf_mul_2exp(value, value, -shift);

    if (sign < 0)
        mpf_neg(value, value);  /* if x is -ve make value -ve */
    return;
}

/* get Hurwitz class number * 12. The Hurwitz class number is not always an integer
but h(n) * 12 always is. */
Znum Hclassno12(const Znum &n) {

    double rvd;
    Znum num, denom;

    specinit();  /* initialise as required*/
    ulong* av = *avma_ref;

    GEN ng = MPtoGEN(ZT(n));
    GEN retval = hclassno_ref(ng);
    GENtoMP(retval, ZT(num), ZT(denom), rvd);

    ptrdiff_t diff = av - *avma_ref;
    if (verbose > 1)
        printf("used %lld bytes on pari stack \n", (long long)diff);
    set_avma_ref((ulong)av);      /* recover memory used */

    return (num * 12) / denom;  /* division is always exact; remainder = 0 */
}

/* get class number. flag = 0 for Shanks method, 1 to use Euler products */
Znum classno(const Znum& n, int flag) {
    Znum num;

    specinit();    /* initialise as required*/
    ulong* av = *avma_ref;

    GEN ng = MPtoGEN(ZT(n));
    /* flag = 0 for Shanks method, 1 to use Euler products */
    GEN retval = qfbclassno_ref(ng, flag);
    InttoMP(retval, ZT(num));
    ptrdiff_t diff = av - *avma_ref;
    if (verbose > 1)
        printf("used %lld bytes on pari stack \n", (long long)diff);
    set_avma_ref((ulong)av);      /* recover memory used */
    return num;
}

/* ramanujantau(n): compute the value of Ramanujan's tau function at n, assuming the GRH.
Algorithm in O(n^{1/2+eps}). */
Znum tau(const Znum& n) {
    Znum num;

    specinit();    /* initialise as required*/
    ulong* av = *avma_ref;

    GEN ng = MPtoGEN(ZT(n));
    GEN retval = tau_ref(ng);
    InttoMP(retval, ZT(num));
    ptrdiff_t diff = av - *avma_ref;
    if (verbose > 1)
        printf("used %lld bytes on pari stack \n", (long long)diff);
    set_avma_ref((ulong)av);      /* recover memory used */
    return num;
}

/* if flag=1 return the Stirling number of the first kind s(n, k), 
if flag=2, return the Stirling number of the second kind S(n, k). */
Znum stirling(const Znum& n, const Znum& m, const Znum& flag) {
    Znum num;

    specinit();     /* initialise as required*/
    ulong* av = *avma_ref;

    /* convert m, n, and flag to 64-bit integers */
    uint64_t ln = MulPrToLong(n);
    uint64_t lm = MulPrToLong(m);
    long long f = MulPrToLong(flag);
    GEN retval = stirling_ref(ln, lm, f);  /* get stirling number */
    InttoMP(retval, ZT(num));  /* convert result to Znums */
    ptrdiff_t diff = av - *avma_ref;
    if (verbose > 1)
        printf("used %lld bytes on pari stack \n", (long long)diff);
    set_avma_ref((ulong)av);      /* recover memory used */
    return num;
}

/* process "PARI ... " command */
void pariParam(const std::string& command) {
    std::string param = command.substr(4);  /* remove "PARI */

    while (param[0] == ' ')
        param.erase(0, 1);              /* remove leading space(s) */

    if (param == "ON" || param == "on") {
        Pari = true;  /* switch to using parilib for factorisation */
        yafu = false;
        msieve = false;
    }
    else if (param == "OFF" || param == "off")
        Pari = false;
    else if (param == "close" || param == "CLOSE") {
        if (stackinit) {
            pari_close_ref();    /* free pari stack */
            stackinit = false;
            if (verbose > 0)
                std::cout << "Pari stack closed \n";
        }
    }
    else
        std::cout << "Invalid Pari command \n";
    return;
}

/* factorise n using parilib factor() function.*/
void parifactor(const Znum& n, fList &factors) {
    GEN nG, flM, fG, prime, exp, expv;
    long long numFactors, typfl, typf, typexp;
    Znum pZ, expZ;

    specinit();    /* initialise as required*/
    ulong* av = *avma_ref;  /* save address of available memory */

    if (verbose > 1) {
        std::cout << "factorise " << n << " using parilib \n";
    }

    nG = MPtoGEN(ZT(n));   /* convert n to t_INT format */
    flM = factor_ref(nG);  /* get factors of n using parilib */
    
    typfl = typ(flM);
    assert(typfl == t_MAT);  /* factors are returned in a matrix*/
  
    fG = gel(flM, 1);     /* the matrix contains 2 columns*/
    typf = typ(fG);       /* 1st column is the (prime?) factors */
    assert(typf == t_COL);
    numFactors = lg(fG)-1;
    assert(numFactors >= 1);

    exp = gel(flM, 2);     /* 2nd column is the exponents */
    typexp = typ(exp);
    assert(typexp == t_COL);

    for (int i = 1; i <= numFactors; i++) {
        prime = gel(fG, i);     /* get factor (should be prime) */
        long long tp = typ(prime);
        assert(tp == t_INT);
        expv = gel(exp, i);    /* and corresponding exponent */
        long long texp = typ(expv);
        assert(texp == t_INT);
        InttoMP(prime, ZT(pZ));   /* convert factor to Znum format */
        insertBigFactor(factors, pZ);     /* add factor to factor list */
        if (verbose > 1) {
            InttoMP(expv, ZT(expZ));
            std::cout << "factor " << pZ << "^" << expZ << "\n";
        }
    }
    set_avma_ref((ulong)av);      /* recover memory used */
}