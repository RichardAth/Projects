/* init.c */
#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_NONSTDC_NO_DEPRECATE
#include <stdio.h>
#include "stack.h"
#include "integer.h"
#include "calc.h"
#ifdef _WIN32
#include "ytab.h"
#else
#include "y.tab.h"
#endif
#include "fun.h"
#include "wrappers.h"


void QSORTMPIX();

const static struct { /* Built-ins which return POLYI */
	char *name;
	POLYI(*func)();
	int argTypes[3];
} builtinsp[] = {
	{ "primitive", PRIMITIVEPI_W ,{ 1, POLY,0 } } ,
	{ "deriv", DERIVPI_W ,{ 1, POLY,0 } } ,
	{ "gcdpi", GCDPI_W ,{ 2, POLY , POLY } } ,
	{ "cyclotomic", CYCLOTOMIC_W ,{ 1, NUM,0 } } ,
	{ "fg", FGPI_W ,{ 2, POLY , POLY } } ,
	{ 0, 0 }
};

const static struct { /* Built-ins which return *MPI */
	char *name;
	MPI *(*func)();
	int argTypes[10];
	/* Type of function is well defined in builtins.txt */
} builtins[] = {
	{ "gcd", GCD_W,{ 2, NUM, NUM } },
	{ "lcm", LCM_W,{ 2, NUM, NUM } },
	{ "lcma", LCM_ARRAY_W,{ 1, ARR } },
	{ "gcdv", EUCLIDI_W,{ 4, NUM, NUM, VARADR, VARADR } },
	{ "jacobi", JACOB_W,{ 2, NUM, NUM } },
	{ "peralta", SQRTMX_W,{ 2, NUM, NUM } },
	{ "gcda", GCD_ARRAY_W,{ 1, ARR } },
	{ "gcdav", GCD_ARRAYVX_W,{ 2, ARR, ARRADR } },
	{ "congr", CONGRX_W,{ 4, NUM, NUM, NUM, VARADR } },
	{ "chinese", CHINESEX_W,{ 5, NUM, NUM, NUM, NUM, VARADR } },
	{ "chinesea", CHINESE_ARRAYX_W,{ 3, ARR, ARR, VARADR } },
	{ "mthroot", BIG_MTHROOTX_W,{ 2, NUM, NUM } },
	{ "fund", FUND_UNITX_W,{ 3, NUM, VARADR, VARADR } },
	{ "pell", PELX_W,{ 4, NUM, NUM, VARADR, VARADR } },
	{ "surd", SURDX_W,{ 9, NUM, NUM, NUM, NUM, ARRADR, ARRADR, ARRADR, ARRADR, ARRADR } },
	{ "mpower", MPOWERX_W,{ 3, NUM, NUM, NUM } },
	{ "nprime", Nextprime_W,{ 1, NUM } },
	{ "inv", INVERSEMX_W,{ 2, NUM, NUM } },
	{ "factor", FACTORX_W,{ 1, NUM } },
	{ "tau", DIVISORX_W,{ 1, NUM } },
	{ "mobius", MOBIUSX_W,{ 1, NUM } },
	{ "euler", EULERX_W,{ 1, NUM } },
	{ "sigma", SIGMAX_W,{ 1, NUM } },
	{ "lprimroot", LPRIMROOTX_W,{ 1, NUM } },
	/*{ "orderp", ORDERPX_W,  {2, NUM, NUM} }, */
	/*	{ "orderq", ORDERQX_W,  {} }, */
	{ "orderm", ORDERMX_W,{ 2, NUM, NUM } },
	/*	{ "lucasu", LUCASUX_W,  {} },
	{ "lucasb", LUCASBX_W,  {} }, */
	{ "lucas", LUCASX_W,{ 1, NUM } },
	{ "length", LENGTHX_W,{ 1, NUM } },
	/* 	{ "mult32", MULT32X_W,  {} }, */
	{ "rsae", RSAE_W,{ 2, NUM, NUM } },
	{ "nprimeap", NEXTPRIMEAPX_W,{ 3, NUM, NUM, NUM } },
	{ "pollard", POLLARD_W,{ 1, NUM } },
	{ "elliptic", EFACTORX_W,{ 3, NUM, NUM, NUM } },
	{ "perfectpower", PERFECT_POWER_W,{ 1, NUM } },
	/*	{ "random", RANDOMI_W,  {} }, */
	{ "leastqnr", LEASTQNRX_W,{ 1, NUM } },
	{ "content", CONTENTPI2_W,{ 1, POLY } },
	{ "sqroot", SQROOT_W,{ 5, NUM, NUM, ARRADR, VARADR, VARADR } },
	{ "congq", QUADRATIC_W,{ 5, NUM, NUM, NUM, NUM, ARRADR } },
	{ "absmod", HALFMOD_W,{ 2, NUM, NUM } },
	{ "ceil", CEILINGI_W,{ 2, NUM, NUM } },
	{ "resultant", SUBRESULTANT_W,{ 2, POLY, POLY } },
	{ "discriminant", DISCRIMINANTPI_W,{ 1, POLY } },
	{ "primes", PRIME_GENERATOR_W,{ 2, NUM, NUM } },
	{ "sturmsequence", STURM_SEQUENCE_W,{ 3, POLY, NUM, NUM } },
	{ "classnop", POS_W,{ 1, NUM } },
	{ "classnon", NEG_W,{ 2, NUM, NUM } },
	{ "nearint", NEARINT_W,{ 2, NUM, NUM } },
	{ "reduceneg", REDUCE_NEG_W,{ 3, NUM, NUM, NUM } },
	{ "reducepos", REDUCE_POS_W,{ 3, NUM, NUM, NUM } },
	{ "classnop0", POS0_W,{ 1, NUM } },
	{ "tableneg", TABLENEG_W,{ 2, NUM, NUM } },
	{ "tablepos", TABLEPOS_W,{ 2, NUM, NUM } },
	{ "davison", DAVISON_W,{ 3, NUM, NUM, NUM } },
	{ "raney", RANEY1_W,{ 4, NUM, NUM, NUM, NUM } },
	{ "unimodular", UNIMODULAR_W,{ 4, NUM, NUM, NUM, NUM } },
	{ "sigmak", SIGMAK_W,{ 2, NUM, NUM } },
	{ "tauprimepower", TAU_PRIMEPOWER_W,{ 2, NUM, NUM } },
	{ "ramanujan", TAU_COMPOSITE_W,{ 1, NUM } },
	{ "repdefinite", REP_DEFINITE_W,{ 5, NUM, NUM, NUM, NUM, NUM } },
	{ "euclid1", EUCLIDI1_W,{ 2, NUM, NUM } },
	{ "rcfperiod", CFRAC_PERIOD_W,{ 1, NUM } },
	{ "nscfperiod", NSCF_PERIOD_W,{ 1, NUM } },
	{ "cfracn", CFRACN_W,{ 1, NUM } },
	{ "spiralinverse", SPIRAL_INVERSE_W,{ 2, NUM, NUM } },
	{ "nicfperiod", NICF_PERIOD_W,{ 1, NUM } },
	{ "nscfperiod0", NSCF_PERIOD0_W,{ 4, NUM, NUM, VARADR, VARADR } },
	{ "rcfperiod0", RCF_PERIOD0_W,{ 4, NUM, NUM, VARADR, VARADR } },
	{ "nicfperiod0", NICF_PERIOD0_W,{ 4, NUM, NUM, VARADR, VARADR } },
	{ "rcfperiod1", CFRAC_PERIOD1_W,{ 1, NUM } },
	{ "nscfperiod1", NSCF_PERIOD1_W,{ 1, NUM } },
	{ "nicfperiod1", NICF_PERIOD1_W,{ 1, NUM } },
	{ "divisorpminus1", DIVISOR_PMINUS1_W,{ 2, NUM, ARRADR } },
	{ "tangent", TANGENT_W,{ 1, NUM } },
	{ "ternary", TERNARY_W,{ 2, NUM, ARRADR } },
	{ "pell4", PEL4X_W,{ 4, NUM, NUM, VARADR, VARADR } },
	{ "stolt", STOLT_W,{ 4, NUM, NUM, NUM, NUM } },
	{ "minarray", MIN_MPI_ARRAY_W,{ 1, ARR } },
	{ "lprimefactor", LPRIME_FACTORX_W,{ 1, NUM } },
	{ "conwaycycles", CONWAY_CYCLESX_W,{ 2, NUM, NUM } },
	{ "conwaycycles0", CONWAY_CYCLES0X_W,{ 2, NUM, NUM } },
	{ 0, 0 }
};

const static struct { /* Built-ins which return void */
	char *name;
	void(*func)();
	int argTypes[10];
} builtinsv[] = {
	{ "serret", SERRET_W,{ 3, NUM, VARADR, VARADR } },
	{ "collatz", COLLATZ_W,{ 2, NUM, NUM } },
	{ "mthrootr", MTHROOTX_W,{ 4, NUM, NUM, NUM, NUM } },
	{ "miller", MILLERX_W,{ 2, NUM, NUM } },
	{ "juggler", JUGGLER_W,{ 2, NUM, NUM } },
	{ "hermite", HERMITE,{ 0 } },
	{ "mlll", MLLL,{ 0 } },
	{ "smith", SMITH,{ 0 } },
	{ "encode", ENCODE_W,{ 2, NUM, NUM } },
	{ "decode", DECODEX_W,{ 3, NUM, NUM, NUM } },
	{ "egcd", EXTGCDX,{ 0 } },
	{ "fp", FINCKE_POHSTX,{ 0 } },
	{ "improvep", IMPROVEPX,{ 0 } },
	/*{ "qsort", QSORTMPIX_W,  {} },
	{ "qsortmat", QSORTMATIX_W,  {} }, */
	{ "sgcd", SCHNORRGCD_W,{ 1, NUM } },
	{ "inhomfp", SHORTESTTTX,{ 0 } },
	/*{ "fibmin", FIB_MIN_W,  {} },
	{ "printww", PRINTWW_W,  {} },
	{ "printdefect", PRINT_DEFECT_W,  {} },*/
	{ "addcubicr", ADD_CUBICRX,{ 0 } },
	{ "powercubicr", POWER_CUBICRX,{ 0 } },
	{ "shallit", SHALLIT,{ 0 } },
	{ "lucasmin", LUCAS_MIN,{ 0 } },
	{ "lllgcd", LLLGCDX,{ 0 } },
	{ "jacobigcd", JACOBIGCDX,{ 0 } },
	{ "shermite", SCHNORRHERMITE_W,{ 1, NUM } },
	{ "lllhermite", LLLHERMITE1X,{ 0 } },
	/*{ "gcd33", GCD33_W,  {} },
	{ "gcd4", GCD4,  {} },
	{ "gcd5", GCD5,  {} },
	{ "gcd6", GCD6_W,  {} },
	{ "gcd10", GCD10_W,  {} },
	{ "gcd11", GCD11_W,  {} },*/
	{ "euclid", EUCLID_W,{ 7, NUM, NUM, ARRADR, ARRADR, ARRADR, ARRADR, VARADR } },
	{ "convergents", CONVERGENTS_W,{ 3, ARR, ARRADR, ARRADR } },
	{ "lagrange", LAGRANGE_W,{ 3, POLY, ARRADR, NUM } },
	{ "axb", AXB1,{ 0 } },
	/*{ "axb1", AXB_W,  {} },
	{ "testaxb", TESTAXB_W,  {} }, */
	/*{ "changel",CHANGELX_W,  {} },
	{ "fermatq",FERMAT_QUOTIENT_W,  {} }, */
	{ "lllgcd0", LLLGCD0X,{ 0 } },
	{ "slv", SLVECTORX,{ 0 } },
	/*{ "gcdconj", GCD_CONJ_W,  {} },
	{ "gcd3", GCD3_W,  {} },
	{ "conj4", GCDCONJECTURE4_W,  {} },
	{ "conj5", GCDCONJECTURE5_W,  {} },
	{ "conj6", GCDCONJECTURE6_W,  {} },
	{ "conj7", GCDCONJECTURE7_W,  {} },
	{ "conjm", GCDCONJECTUREM_W,  {} }, */
	{ "absnearint", ABS_NEAREST_INTRX,{ 0 } },
	{ "cycle", CYCLEX,{ 0 } },
	{ "addcubicm", ADD_CUBICMX,{ 0 } },
	{ "powercubicm", POWER_CUBICMX,{ 0 } },
	{ "ordercubicm", ORDER_CUBICMX,{ 0 } },
	{ "powercubicr", POWER_CUBICRX,{ 0 } },
	{ "ordercubicr", ORDER_CUBICRX,{ 0 } },
	{ "sturm", STURM_W,{ 1, POLY } },
	{ "rootexp", ROOTEXPANSION,{ 2, POLY, NUM } },
	{ "intlog", INTLOG,{ 4, NUM, NUM, NUM, NUM } },
	{ "intlog1", INTLOG1,{ 4, NUM, NUM, NUM, NUM } },
	{ "log1", LOG1,{ 5, NUM, NUM, NUM, NUM, NUM } },
	{ "intlog2", INTLOG2,{ 4, NUM, NUM, NUM, NUM } },
	{ "logg", LOGG,{ 5, NUM, NUM, NUM, NUM, NUM } },
	{ "log2", LOGGG,{ 6, NUM, NUM, NUM, NUM, NUM, NUM } },
	{ "sqroot1", SQROOT1_W,{ 3, NUM, NUM, NUM } },
	{ "sqroot2", SQROOT2_W,{ 2, NUM, NUM } },
	{ "sqroot3", SQROOT3_W,{ 3, NUM, NUM, NUM } },
	{ "cornacchia", CORNACCHIA_W,{ 3, NUM, NUM, NUM } },
	{ "patz", PATZ_W,{ 2, NUM, NUM } },
	{ "shankslog", SHANKSLOG,{ 4, NUM, NUM, NUM, NUM } },
	/*{ "binaryform", BINARYFORM_W,  {4, NUM, NUM, NUM, NUM} },*/
	{ "gauss", GAUSS_W,{ 7, NUM, NUM, NUM, NUM, VARADR, VARADR, VARADR } },
	{ "binform", BINFORM_W,{ 5, NUM, NUM, NUM, NUM, NUM } },
	{ "log", LOG_W,{ 6, NUM, NUM, NUM, NUM, ARRADR, VARADR } },
	{ "testlog1", TESTLOG1_W,{ 4, NUM, NUM, NUM, NUM } },
	{ "testlog", TESTLOG_W,{ 5, NUM, NUM, NUM, NUM, NUM } },
	{ "twoadicsqrt", TWOADICSQRT_W,{ 3, NUM, NUM, ARRADR } },
	{ "padicsqrt", PADICSQRT_W,{ 4, NUM, NUM, NUM, ARRADR } },
	{ "powerd", POWERD_W,{ 6, NUM, NUM, NUM, NUM, VARADR, VARADR } },
	{ "midptcount", MIDPT_COUNT_W,{ 2, NUM, NUM } },
	{ "cfraccount", CFRAC_COUNT_W,{ 2, NUM, NUM } },
	{ "algcycle", ALG_CYCLEX,{ 0 } },
	{ "spiral", SPIRAL_W,{ 3, NUM, VARADR, VARADR } },
	{ "rcfcount", TIME_COUNT_RCF_W,{ 2, NUM, NUM } },
	{ "nicfcount", TIME_COUNT_NICF_W,{ 2, NUM, NUM } },
	{ "nscfcount", TIME_COUNT_NSCF_W,{ 2, NUM, NUM } },
	{ "rcfperiodcount", PERIOD_TIME_COUNT_RCF_W,{ 2, NUM, NUM } },
	{ "nscfperiodcount", PERIOD_TIME_COUNT_NSCF_W,{ 2, NUM, NUM } },
	{ "nicfperiodcount", PERIOD_TIME_COUNT_NICF_W,{ 2, NUM, NUM } },
	{ "rcfcount00", TIME_COUNT_RCF00_W,{ 2, NUM, NUM } },
	{ "nicfcount00", TIME_COUNT_NICF00_W,{ 2, NUM, NUM } },
	{ "nscfcount00", TIME_COUNT_NSCF00_W,{ 2, NUM, NUM } },
	{ "rcfcount000", TIME_COUNT_RCF000_W,{ 2, NUM, NUM } },
	{ "nscfcount000", TIME_COUNT_NSCF000_W,{ 2, NUM, NUM } },
	{ "nicfcount000", TIME_COUNT_NICF000_W,{ 2, NUM, NUM } },
	{ "rcfcountd", TIME_COUNT_RCF00_D_W,{ 2, NUM, NUM } },
	{ "nicfcountd", TIME_COUNT_NICF00_D_W,{ 2, NUM, NUM } },
	{ "nscfcountd", TIME_COUNT_NSCF00_D_W,{ 2, NUM, NUM } },
	{ "rcfperiodcount1", PERIOD_TIME_COUNT_RCF1_W,{ 2, NUM, NUM } },
	{ "nscfperiodcount1", PERIOD_TIME_COUNT_NSCF1_W,{ 2, NUM, NUM } },
	{ "nicfperiodcount1", PERIOD_TIME_COUNT_NICF1_W,{ 2, NUM, NUM } },
	{ "rcfcountqimproved", TIME_COUNT_RCF_QIMPROVED_W,{ 2, NUM, NUM } },
	{ "nicfcountqimproved", TIME_COUNT_NICF_QIMPROVED_W,{ 2, NUM, NUM } },
	{ "sortarray", SORT_ARRAY_MPI_W,{ 3, ARR, ARRADR, NUM } },
	{ "carmichael", CARMICHAEL_W,{ 1, NUM } },
	{ "carnielli", CARNIELLI_CYCLEX,{ 0 } },
	{ "lupei", LUPEI_CYCLEX,{ 0 } },
	{ "bernoulli", BERNOULLI_W,{ 3, NUM, VARADR, VARADR } },
	{ "twocycle", TWO_CYCLEX,{ 0 } },
	{ "mcycle", M_CYCLEX,{ 0 } },
	{ "pqcycle", PQ_CYCLEX,{ 0 } },
	{ "mmcycle", MM_CYCLEX,{ 0 } },
	{ "exceptionals", EXCEPTIONALS_W,{ 1, NUM } },
	{ "exceptionals0", EXCEPTIONALS0_W,{ 1, NUM } },
	{ "exceptionals10", EXCEPTIONALS10_W,{ 1, NUM } },
	{ "gplus", GPLUS_W,{ 3, POLY, POLY, POLY } },
	{ "gzero", GZERO_W,{ 3, POLY, POLY, POLY } },
	{ "gminus", GMINUS_W,{ 3, POLY, POLY, POLY } },
	{ "nagell", NAGELL_W,{ 2, NUM, NUM } },
	{ "frattini", FRATTINI_W,{ 2, NUM, NUM } },
	{ "kashihara", KASHIHARA_W,{ 2, NUM, NUM } },
	{ "branch", BRANCH_W,{ 3, NUM, NUM, NUM } },
	{ "conwayrangetest", CONWAY_RANGE_TEST_W,{ 2, NUM, NUM } },
	{ "conwayrangetest1", CONWAY_RANGE_TEST1_W,{ 4, NUM, NUM, NUM, NUM } },
	{ 0, 0 }
};

void init()
/* install built-ins in table */
{
	int i;
	Symbol *s;

	for (i = 0; builtins[i].name; i++)
	{
		s = installFunc(builtins[i].name, BLTIN, builtins[i].argTypes);
		s->u.ptr = builtins[i].func;
	}
	for (i = 0; builtinsp[i].name; i++)
	{
		s = installFunc(builtinsp[i].name, BLTINP, builtinsp[i].argTypes);
		s->u.ptrp = builtinsp[i].func;
	}
	for (i = 0; builtinsv[i].name; i++)
	{
		s = installFunc(builtinsv[i].name, BLTINV, builtinsv[i].argTypes);
		s->u.ptrv = builtinsv[i].func;
	}
}
