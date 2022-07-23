#include "pch.h"
#include <windows.h>
#include "pari.h"

#define ZT(a) a.backend().data()  /* access mpz_t within a Znum (Boost mpz_int)*/

bool fRunTimeLinkSuccess = false;

/* typedefs for access to libpari functions */
typedef    void(__cdecl* pari_initX)(size_t parisize, ulong maxprime);    /*void    pari_init(size_t parisize, ulong maxprime);*/
typedef    void(__cdecl* pari_stack_initX)(pari_stack* s, size_t size, void** data);
typedef int64_t(__cdecl* pari_stack_newX)(pari_stack* s);
typedef    void(__cdecl* pari_init_optsX)(size_t parisize, ulong maxprime, ulong init_opts);
typedef    void(__cdecl* pari_init_primesX)(ulong maxprime);
typedef   ulong(__cdecl* hclassno6uX)(ulong D);   /* ulong hclassno6u(ulong D)*/
typedef     GEN(__cdecl* hclassnoX)(GEN x);        /* GEN     hclassno(GEN x);*/

typedef GEN(__cdecl* addiiX)(GEN x, GEN y);   /* GEN addii(GEN x, GEN y) etc*/
typedef GEN(__cdecl* subiiX)(GEN x, GEN y);
typedef GEN(__cdecl* muliiX)(GEN x, GEN y);
typedef GEN(__cdecl* dvmdiiX)(GEN x, GEN y, GEN* z);
typedef int(__cdecl* invmodX)(GEN a, GEN b, GEN* res);  /* int  invmod(GEN a, GEN b, GEN *res); */
/* GEN    itor(GEN x, int64_t prec); */
typedef GEN(__cdecl* itorX)(GEN x, int64_t prec);
/* double  rtodbl(GEN x); */
typedef double(__cdecl* rtodblX)(GEN x);
/* GEN     dbltor(double x); */
typedef GEN(__cdecl* dbltorX)(double x);
/* GEN    addrr(GEN x, GEN y); etc */
typedef GEN(__cdecl* addrrX) (GEN x, GEN y);
typedef GEN(__cdecl* subrrX) (GEN x, GEN y);
typedef GEN(__cdecl* mulrrX) (GEN x, GEN y);
typedef GEN(__cdecl* divrrX) (GEN x, GEN y);
typedef GEN(__cdecl* divruX) (GEN x, ulong y);
typedef void(__cdecl* pari_print_versionX)(void);   /* void pari_print_version(void)*/
//typedef pari_mainstack* (__cdecl* parimainstackX)();
typedef GEN(__cdecl* expX)(GEN x, int64_t prec);  /* GEN gexp(GEN x, int64_t prec)*/
typedef GEN(__cdecl* logX)(GEN x, int64_t prec);  /* GEN glog(GEN x, int64_t prec)*/
typedef char* (__cdecl* GENtostrX)(GEN x);        /* char*   GENtostr(GEN x);*/
typedef GEN(__cdecl* utoiposX)(ulong x);     /* GEN utoipos(ulong x)*/
typedef GEN(__cdecl* utoinegX)(ulong x);
typedef GEN(__cdecl* floorrX)(GEN x);    /* GEN floorr(GEN x)*/
typedef GEN(__cdecl* cgetiposX)(int64_t x); /* GEN    cgetipos(int64_t x);*/
typedef GEN(__cdecl* cgetinegX)(int64_t x);
typedef GEN(__cdecl* qfbclassnox)(GEN x, int64_t flag);  /* GEN     qfbclassno0(GEN x,int64_t flag);*/
typedef ulong** avmaX;
typedef void(__cdecl* set_avmaX)(ulong av);  /* void   set_avma(ulong av);*/
typedef GEN(__cdecl* tauX)(GEN n);     /* GEN ramanujantau(GEN n)*/
typedef GEN(__cdecl* stirlingX)(int64_t n, int64_t m, int64_t flag); /* GEN stirling(int64_t n, int64_t m, int64_t flag)*/

HINSTANCE hinstLib;

/* function pointers to access libpari functions */
static pari_initX        pari_init_ref;
static pari_stack_initX  pari_stack_init_ref;
static pari_stack_newX   pari_stack_new_ref;
static pari_init_optsX   pari_init_opts_ref;
static pari_init_primesX pari_init_primes_ref;
static addiiX            addii_ref;
static subiiX            subii_ref;
static muliiX            mulii_ref;
static dvmdiiX           dvmdii_ref;
static hclassno6uX       hclassno6u_ref;
static hclassnoX         hclassno_ref;
static invmodX           invmod_ref;
static itorX             itor_ref;
static rtodblX           rtodbl_ref;
static dbltorX           dbltor_ref;
static addrrX            addrr_ref;
static subrrX            subrr_ref;
static divrrX            divrr_ref;
static mulrrX            mulrr_ref;
static divruX            divru_ref;
static pari_print_versionX pari_print_version_ref;
//parimainstackX    parimainstack_ref;
static expX              exp_ref;
static logX              log_ref;
static GENtostrX         GENtostr_ref;
static utoiposX          utoipos_ref;
static utoinegX          utoineg_ref;
static floorrX           floorr_ref;
static cgetiposX         cgetipos_ref;
static cgetinegX         cgetineg_ref;
static qfbclassnox       qfbclassno_ref;
static avmaX             avma_ref;
static set_avmaX         set_avma_ref;
static tauX              tau_ref;
static stirlingX         stirling_ref;

/* pari library functions are accessed in this way because linking statically to
libpari.dll doesn't work, presumably because it was compiled with Msys2/gcc and
we are using windows + Visual Studio.
for every function we want to use, we have to set up a function pointer.
For some reason dynamic linking works OK. */
static void specinit()
{
    // Get a handle to the pari DLL module.
    hinstLib = LoadLibraryA("C:/Program Files (x86)/Pari64-2-13-2/libpari.dll");

    // If the handle is valid, try to get the function addresses.
    if (hinstLib == nullptr) {
        std::cerr << "could not access dll for PARI \n";
        system("PAUSE");
        abort();
    }

    else {
        /* set up function pointers. There over 6000 accessible libpari functions. This is
        a selection of functions that might be useful. */
        pari_init_ref        = (pari_initX)GetProcAddress(hinstLib, "pari_init");
        pari_stack_init_ref  = (pari_stack_initX)GetProcAddress(hinstLib, "pari_stack_init");
        pari_stack_new_ref   = (pari_stack_newX)GetProcAddress(hinstLib, "pari_stack_new");
        pari_init_opts_ref   = (pari_init_optsX)GetProcAddress(hinstLib, "pari_init_opts");
        pari_init_primes_ref = (pari_init_primesX)GetProcAddress(hinstLib, "pari_init_primes");
        addii_ref            = (addiiX)GetProcAddress(hinstLib, "addii");
        subii_ref            = (subiiX)GetProcAddress(hinstLib, "subii");
        mulii_ref            = (muliiX)GetProcAddress(hinstLib, "mulii");
        dvmdii_ref           = (dvmdiiX)GetProcAddress(hinstLib, "dvmdii");
        hclassno6u_ref       = (hclassno6uX)GetProcAddress(hinstLib, "hclassno6u");
        hclassno_ref         = (hclassnoX)GetProcAddress(hinstLib, "hclassno");
        invmod_ref           = (invmodX)GetProcAddress(hinstLib, "invmod");
        itor_ref             = (itorX)GetProcAddress(hinstLib, "itor");
        rtodbl_ref           = (rtodblX)GetProcAddress(hinstLib, "rtodbl");
        dbltor_ref           = (dbltorX)GetProcAddress(hinstLib, "dbltor");
        addrr_ref            = (addrrX)GetProcAddress(hinstLib, "addrr");
        subrr_ref            = (subrrX)GetProcAddress(hinstLib, "subrr");
        divrr_ref            = (divrrX)GetProcAddress(hinstLib, "divrr");
        mulrr_ref            = (mulrrX)GetProcAddress(hinstLib, "mulrr");
        divru_ref            = (divruX)GetProcAddress(hinstLib, "divru");
        pari_print_version_ref = (pari_print_versionX)GetProcAddress(hinstLib, "pari_print_version");
        exp_ref              = (expX)GetProcAddress(hinstLib, "gexp");
        log_ref              = (logX)GetProcAddress(hinstLib, "glog");
        utoipos_ref          = (utoiposX)GetProcAddress(hinstLib, "utoipos");
        utoineg_ref          = (utoinegX)GetProcAddress(hinstLib, "utoineg");
        GENtostr_ref         = (GENtostrX)GetProcAddress(hinstLib, "GENtostr");
        floorr_ref           = (floorrX)GetProcAddress(hinstLib, "floorr");
        cgetipos_ref         = (cgetiposX)GetProcAddress(hinstLib, "cgetipos");
        cgetineg_ref         = (cgetinegX)GetProcAddress(hinstLib, "cgetineg");
        qfbclassno_ref       =(qfbclassnox)GetProcAddress(hinstLib, "qfbclassno0");
        avma_ref             = (avmaX)GetProcAddress(hinstLib, "avma");
        set_avma_ref         = (set_avmaX)GetProcAddress(hinstLib, "set_avma");
        tau_ref              = (tauX)GetProcAddress(hinstLib, "ramanujantau");
        stirling_ref         =(stirlingX)GetProcAddress(hinstLib, "stirling");


        /* check that all function pointers were set up successfully */
        if (nullptr == pari_init_ref ||
            nullptr == pari_stack_init_ref ||
            nullptr == pari_stack_new_ref ||
            nullptr == pari_init_opts_ref ||
            nullptr == pari_init_primes_ref ||
            nullptr == addii_ref ||
            nullptr == subii_ref ||
            nullptr == mulii_ref ||
            nullptr == dvmdii_ref ||
            nullptr == invmod_ref ||
            nullptr == itor_ref ||
            nullptr == rtodbl_ref ||
            nullptr == dbltor_ref ||
            nullptr == addrr_ref ||
            nullptr == subrr_ref ||
            nullptr == divrr_ref ||
            nullptr == divru_ref ||
            nullptr == pari_print_version_ref ||
            nullptr == hclassno_ref ||
            nullptr == exp_ref ||
            nullptr == log_ref ||
            nullptr == utoipos_ref ||
            nullptr == utoineg_ref ||
            nullptr == GENtostr_ref ||
            nullptr == floorr_ref ||
            nullptr == cgetipos_ref ||
            nullptr == cgetineg_ref ||
            nullptr == hclassno6u_ref ||
            nullptr == qfbclassno_ref ||
            nullptr == avma_ref ||
            nullptr == set_avma_ref ||
            nullptr == tau_ref ||
            nullptr == stirling_ref) {
            fRunTimeLinkSuccess = false;
            std::cerr << "PARI dynamic linking failed \n";
            system("PAUSE");
            abort();
        }

        fRunTimeLinkSuccess = true;

        pari_init_ref(8000000, 500000);  /* stack size, maxprime */

        if (verbose > 0) {
            pari_print_version_ref();  /* print some info from libpari dll */
            printf("avma = %p\n", *avma_ref);
        }
    }
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
    int64_t typx = typ(x);

    switch (typx) {
    case t_INT: {
        int64_t s = signe(x);
        int64_t i, lx = lgefint(x);

        mpz_set_ui(denom, 1);
        if (s == 0) {
            mpz_set_ui(value, 0);
            return;
        }

        mpz_set_ui(value, 0);   /* value = 0 */
        int64_t* y = x + lx - 1;  /* y points to most significant word of x */
        for (i = 2; i < lx; i++, y--) {
            mpz_mul_2exp(value, value, 64); /* move more significant bits left 64 bits*/
            mpz_add_ui(value, value, *y);   /* add next 64-bit limb */
        }
        if (s < 0)
            mpz_neg(value, value);  /* if x  is -ve, flip sign */

        val_d = mpz_get_d(value);  /* convert to floating point */
        return;
    }
    case t_FRAC: {
        GEN num = gel(x, 1);      /* numerator */
        GEN den = gel(x, 2);     /* denominator */
        int64_t s = signe(num);
        int64_t i, lnum = lgefint(num), lden = lgefint(den);
        if (s == 0) {
            mpz_set_ui(value, 0);
            mpz_set_ui(denom, 1);
            return;
        }

        mpz_set_ui(value, 0);   /* value = 0 */
        int64_t* y = num + lnum - 1;  /* y points to most significant word of num */
        for (i = 2; i < lnum; i++, y--) {
            mpz_mul_2exp(value, value, 64);
            mpz_add_ui(value, value, *y);
        }
        if (s < 0)
            mpz_neg(value, value);  /* if num  is -ve, flip sign */

        s = signe(den);
        if (s == 0) {
            mpz_set_ui(denom, 0);  /* this would be weird! Mathematically the value is infinite!*/
            return;
        }

        mpz_set_ui(denom, 0);   /* denom = 0 */
        y = den + lden - 1;  /* y points to most significant word of denom */
        for (i = 2; i < lden; i++, y--) {
            mpz_mul_2exp(denom, denom, 64);
            mpz_add_ui(denom, denom, *y);
        }
        if (s < 0)
            mpz_neg(denom, denom);  /* if denom  is -ve, flip sign */

        val_d = mpz_get_d(value) / mpz_get_d(denom);
        return;
    }
    case t_REAL: {
        GEN xi = floorr_ref(x);   /* xi = integer part of x*/
        int64_t s = signe(xi);
        int64_t i, lx = lgefint(xi);
        mpz_set_ui(denom, 1);

        /* convert (GEN)xi to (mpz_t)value */
        if (s == 0) {
            mpz_set_ui(value, 0);
        }
        else {
            mpz_set_ui(value, 0);   /* value = 0 */
            int64_t* y = xi + lx - 1;  /* y points to most significant word of xi */
            for (i = 2; i < lx; i++, y--) {
                mpz_mul_2exp(value, value, 64); /* move more significant bits left 64 bits*/
                mpz_add_ui(value, value, *y);   /* add next 64-bit limb */
            }
            if (s < 0)
                mpz_neg(value, value);  /* if xi  is -ve, flip sign */
        }

        val_d = rtodbl_ref(x);  /* also convert x to normal floating point */

        return;
    }
    default:
        abort();
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
assume that R3(n) < 2^64. Note: n is a copy of the original value. This copy
is modified. */
uint64_t R3h(Znum n) {
    mpz_t num, denom, result;
    mpz_inits(num, denom, result, nullptr);
    double hxd;
    GEN x, hx;
    uint64_t r;

    if (fRunTimeLinkSuccess == false)
        specinit();
    ulong* av = *avma_ref;  /* save address of available memory */

    if (n == 0)
        return 1;   /* easiest to make n = 0 a special case*/
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
        r = mpz_get_ui(result);
        break;
 
    }

    case 3:  /* n%8 = 3 */
        /* get 24 * h(n) */
        x = MPtoGEN(ZT(n));
        hx = hclassno_ref(x);
        GENtoMP(hx, num, denom, hxd);
        mpz_mul_ui(result, num, 24);     /* result = num * 24 */
        mpz_div(result, result, denom);  /* result =(num * 24)/denom  = 24*h(n)  */
        r = mpz_get_ui(result);         /* convert to 64-bit integer */
        break;

    case 7:  /* n%8 = 7 */
        return 0;

    default:   /* n%8 = 0 or 4 */
        abort();   /* should never happen */
    }

    mpz_clears(num, denom, result, nullptr);
    ptrdiff_t diff = av - *avma_ref;
    if (verbose > 1)
        printf("used %lld bytes on pari stack \n", (long long)diff);
    set_avma_ref((ulong)av);      /* recover memory used */
    return r;            /* return result */
}

/* get Hurwitz class number * 12. The Hurwitz class number is not always an integer
but h(n) * 12 always is. */
Znum Hclassno12(const Znum &n) {

    double rvd;
    Znum num, denom;

    if (fRunTimeLinkSuccess == false)
        specinit();
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

/* get class number flag = 0 for Shanks method, 1 to use Euler products */
Znum classno(const Znum& n, int flag) {
    double rvd;
    Znum num, denom;

    if (fRunTimeLinkSuccess == false)
        specinit();
    ulong* av = *avma_ref;

    GEN ng = MPtoGEN(ZT(n));
    /* flag = 0 for Shanks method, 1 to use Euler products */
    GEN retval = qfbclassno_ref(ng, flag);
    GENtoMP(retval, ZT(num), ZT(denom), rvd);
    ptrdiff_t diff = av - *avma_ref;
    if (verbose > 1)
        printf("used %lld bytes on pari stack \n", (long long)diff);
    set_avma_ref((ulong)av);      /* recover memory used */
    return num / denom;
}

/* ramanujantau(n): compute the value of Ramanujan's tau function at n, assuming the GRH.
Algorithm in O(n^{1/2+eps}). */
Znum tau(const Znum& n) {
    double rvd;
    Znum num, denom;

    if (fRunTimeLinkSuccess == false)
        specinit();
    ulong* av = *avma_ref;

    GEN ng = MPtoGEN(ZT(n));
    GEN retval = tau_ref(ng);
    GENtoMP(retval, ZT(num), ZT(denom), rvd);
    ptrdiff_t diff = av - *avma_ref;
    if (verbose > 1)
        printf("used %lld bytes on pari stack \n", (long long)diff);
    set_avma_ref((ulong)av);      /* recover memory used */
    return num / denom;
}

/* if flag=1 return the Stirling number of the first kind s(n, k), 
if flag=2, return the Stirling number of the second kind S(n, k). */
Znum stirling(const Znum& n, const Znum& m, const Znum& flag) {
    double rvd;
    Znum num, denom;

    if (fRunTimeLinkSuccess == false)
        specinit();
    ulong* av = *avma_ref;

    /* convert m, n, and flag to 64-bit integers */
    uint64_t ln = MulPrToLong(n);
    uint64_t lm = MulPrToLong(m);
    long long f = MulPrToLong(flag);
    GEN retval = stirling_ref(ln, lm, f);  /* get stirling number */
    GENtoMP(retval, ZT(num), ZT(denom), rvd);  /* convert result to Znums */
    ptrdiff_t diff = av - *avma_ref;
    if (verbose > 1)
        printf("used %lld bytes on pari stack \n", (long long)diff);
    set_avma_ref((ulong)av);      /* recover memory used */
    return num / denom;
}