/* generated from ecm-h.in by gen_ecm_h.bat */
/* ecm.h - public interface for libecm.
Copyright 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011
Paul Zimmermann, Alexander Kruppa, David Cleaver, Cyril Bouvier.
This file is part of the ECM Library.
The ECM Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 3 of the License, or (at your
option) any later version.
The ECM Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.
You should have received a copy of the GNU Lesser General Public License
along with the ECM Library; see the file COPYING.LIB.  If not, see
http://www.gnu.org/licenses/ or write to the Free Software Foundation, Inc.,
51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA. */
#ifndef _ECM_H 
#define _ECM_H 1
#include <stdio.h> /* for FILE */
#include <gmp.h> 
#define ECM_VERSION "7.0.4"
#ifdef __cplusplus 
extern "C" {
#endif  
#define EC_W_NBUFS 9 /* for Hessian form */
/* More ec forms */
#define ECM_EC_TYPE_MONTGOMERY 1
#define ECM_EC_TYPE_WEIERSTRASS 2
#define ECM_EC_TYPE_HESSIAN 3
#define ECM_EC_TYPE_WEIERSTRASS_COMPLETE 4
/* which type of law used */
#define ECM_LAW_AFFINE 1
#define ECM_LAW_HOMOGENEOUS 2
typedef struct 
{  
int type; 
int law; 
mpz_t a4; /* for MONTGOMERY: b*y^2=x^3+A*x^2+x 
for WEIERSTRASS: y^2=x^3+A*x+B
for HESSIAN: U^3+V^3+W^3=3*A*U*V*W */
mpz_t a1, a3, a2, a6;  /* for complete WEIERSTRASS */
mpz_t buf[EC_W_NBUFS]; /* used in the addition laws */
int disc; /* in case E is known to have CM by Q(sqrt(disc)) */
mpz_t sq[10]; /* for CM curves, we might have squareroots */
} __ell_curve_struct; 
typedef __ell_curve_struct ell_curve_t[1];
typedef struct 
{  
mpz_t x; 
mpz_t y; 
mpz_t z; 
} __ell_point_struct; 
typedef __ell_point_struct ell_point_t[1];
typedef struct 
{  
int method; /* factorization method, default is ecm */
mpz_t x, y;        /* starting point (if non zero) */
int param; /* (ECM only) What parametrization do we used */ 
mpz_t sigma; /* (ECM only) The parameter for the parametrization */ 
/* May contains A */
int sigma_is_A; /* if  1, 'parameter' contains A (Montgomery form),
if 0, 'parameter' contains sigma (Montgomery form),
if -1, 'parameter' contains A, and the input curve is in
Weierstrass form y^2 = x^3 + A*x + B, with y in 'go'. */
__ell_curve_struct *E; /* the curve, particularly useful for CM ones */
mpz_t go; /* initial group order to preload (if NULL: do nothing),
or y for Weierstrass form if sigma_is_A = -1. */
double B1done; /* step 1 was already done up to B1done */
mpz_t B2min; /* lower bound for stage 2 (default is B1) */
mpz_t B2; /* step 2 bound (chosen automatically if < 0.0) */
unsigned __int64 k;/* number of blocks in stage 2 */
int S; /* degree of the Brent-Suyama's extension for stage 2 */
int repr; /* representation for modular arithmetic: ECM_MOD_MPZ=mpz,         
ECM_MOD_MODMULN=modmuln (Montgomery's quadratic multiplication),
ECM_MOD_REDC=redc (Montgomery's subquadratic multiplication),
ECM_MOD_GWNUM=Woltman's gwnum routines (tbd),
> 16 : special base-2 representation        
MOD_DEFAULT: automatic choice */
int nobase2step2; /* disable special base-2 code in ecm stage 2 only */
int verbose; /* verbosity level: 0 no output, 1 normal output,   
2 diagnostic output */
FILE *os; /* output stream (for verbose messages) */
FILE *es; /* error  stream (for error   messages) */
char *chkfilename; /* Filename to write stage 1 checkpoints to */
char *TreeFilename; /* Base filename for storing product tree of F */
double maxmem; /* Maximal amount of memory to use in stage 2, in bytes.
0. means no limit (optimise only for speed) */
double stage1time; /* Time to add for estimating expected time to find fac.*/
gmp_randstate_t rng; /* State of random number generator */
int use_ntt; /* set to 1 to use ntt poly code in stage 2 */
int (*stop_asap) (void); /* Pointer to function, if it returns 0, contine 
normally, otherwise exit asap. May be NULL */
/* The batch mode is used for stage 1 when param=1 or param=2)*/
mpz_t batch_s; /* s is the product of primes up to B1 for batch mode */
double batch_last_B1_used; /* Last B1 used in batch mode. Used to avoid */
/* computing s when B1 = batch_last_B1_used */
int gpu; /* do we use the GPU for stage 1. */
/* If different from 0, the GPU is used */
/* Else, the parameters beginning by gpu_* have no meaning */
int gpu_device; /* Which device do we use */
int gpu_device_init; /* Is the device initialized?*/
unsigned int gpu_number_of_curves; 
double gw_k; /* use for gwnum stage 1 if input has form k*b^n+c */
unsigned __int64 gw_b;  /* use for gwnum stage 1 if input has form k*b^n+c */
unsigned __int64 gw_n;  /* use for gwnum stage 1 if input has form k*b^n+c */
signed __int64 gw_c;    /* use for gwnum stage 1 if input has form k*b^n+c */
} __ecm_param_struct; 
typedef __ecm_param_struct ecm_params[1];
typedef __ecm_param_struct *ecm_params_ptr;
#define ECM_MOD_NOBASE2 -1
#define ECM_MOD_DEFAULT 0
#define ECM_MOD_MPZ 1
#define ECM_MOD_BASE2 2
#define ECM_MOD_MODMULN 3
#define ECM_MOD_REDC 4
/* values <= -16 or >= 16 have a special meaning */
const char *ecm_version();
int ecm_factor (mpz_t, mpz_t, double, ecm_params);
void ecm_init (ecm_params);
void ecm_clear (ecm_params);
/* the following interface is not supported */
int ecm (mpz_t, mpz_t, mpz_t, int*, mpz_t, mpz_t, mpz_t, double *, double, mpz_t, mpz_t,
unsigned long long, const int, int, int, int, int, int, 
ell_curve_t, FILE* os, FILE* es,
char*, char *, double, double, gmp_randstate_t, int (*)(void), mpz_t, 
double *, double, unsigned __int64, unsigned __int64, signed __int64);
int pp1 (mpz_t, mpz_t, mpz_t, mpz_t, double *, double, mpz_t, mpz_t, 
signed __int64, int, int, int, FILE*, FILE*, char*,
char *, double, gmp_randstate_t, int (*)(void));
int pm1 (mpz_t, mpz_t, mpz_t, mpz_t, double *, double, mpz_t, 
mpz_t, unsigned long long, int, int, int, FILE*, 
FILE*, char *, char*, double, gmp_randstate_t, int (*)(void));
/* different methods implemented */
#define ECM_ECM 0
#define ECM_PM1 1
#define ECM_PP1 2
/* return value of ecm, pm1, pp1 */
#define ECM_FACTOR_FOUND_STEP1 1 /* should be positive */
#define ECM_FACTOR_FOUND_STEP2 2 /* should be positive */
#define ECM_NO_FACTOR_FOUND 0 /* should be zero */
#define ECM_ERROR -1 /* should be non-zero */
#define ECM_FACTOR_FOUND_P(x) ((x) > 0)
#define ECM_ERROR_P(x) ((x) < 0)
#define ECM_DEFAULT_B1_DONE 1.0
#define ECM_IS_DEFAULT_B1_DONE(x) (x <= 1.0)
/* Different parametrizations used in stage 1 of ECM */
#define ECM_PARAM_DEFAULT -1
#define ECM_PARAM_SUYAMA 0
#define ECM_PARAM_BATCH_SQUARE 1
#define ECM_PARAM_BATCH_2 2
#define ECM_PARAM_BATCH_32BITS_D 3
/* we keep 4 as spare */
#define ECM_PARAM_WEIERSTRASS 5
#define ECM_PARAM_HESSIAN 6
#define ECM_PARAM_TORSION 7
/* stage 2 bound */
#define ECM_DEFAULT_B2 -1
#define ECM_IS_DEFAULT_B2(x) (mpz_cmp_si (x, ECM_DEFAULT_B2) == 0)
#define ECM_DEFAULT_K 0 /* default number of blocks in stage 2. 0 = automatic
choice */ 
#define ECM_DEFAULT_S 0 /* polynomial is chosen automatically */
/* Apple uses '\r' for newlines */
#define IS_NEWLINE(c) (((c) == '\n') || ((c) == '\r'))
#ifdef __cplusplus 
}  
#endif  
#endif /* _ECM_H */
