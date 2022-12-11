/* calc.h */
#include "stack.h"

/* Constants for parse.y */
#define RET_NULL 0
#define RET_MPI 1
#define RET_MPIA 2
#define RET_POLY 3
#define FUNC_FAIL 4

extern int rettype;
/* end constants for parse.y */
typedef struct _Symbol { /* symbol table entry */
						 /* Ref: The Unix Programming Environment, Kernighan and Pike p. 247. */
	char *name;
	int type; /* ARRAY, VAR, BLTIN, BLTINV, UNDEF */
	union {
		MPIA symarr;            /* if ARRAY  */
		MPI *symval;             /* if VAR    */
		POLYI sympval;             /* if POLY */
		MPI *(*ptr)();    	 /* if BLTIN  */
		POLYI(*ptrp)();       /* if BLTINP */
		void(*ptrv)();    	 /* if BLTINV */
	} u;
	int *argTypes;  /* only for BLTIN, BLTINP, BLTINV */
	struct _Symbol *next;   /* to link to another */
} Symbol;

typedef struct _Argument { /* symbol table entry */
						   /* Ref: The Unix Programming Environment, Kernighan and Pike p. 247. */
	int type; /* ARRAY, VARADR, NUMBER, POLY, ARRAYADR */
	int defined; /* only used for ARRADR and VARADR. 0 for undefined.
				 * 1 for defined */
	union {
		MPIA *arrayAdr;        /* if ARRADR */
		MPIA array;            /* if ARR  */
		MPI **varAdr;         /* if VARADR */
		MPI *num;             /* if NUM    */
		POLYI poly;             /* if POLY */
	} u;
} *Argument;

Symbol *install(char *s, int typ), *lookup(char *s, int typ);
/* add function to symlist. s = name, t=type */
Symbol *installFunc(const char *s, int typ, const int *argTypes);
Argument createArg(void *arg, int type, int defined);
void freeArg(Argument Arg);

/* constants to be known  to checkArgs, freeDatum, freeArg and createArg */
#define POLY 512
#define ARRADR 513
#define VARADR  514
#define NUM 515
#define ARR 516 /* note there is a different between this and ARRAY */


/* A note on the new addition 'int *argTypes'.  This a pointer to a
* dynamically allocated array of integers.  The first entry argTypes[0]
* contains the number of integers _left_ in the array.
* ie. If argTypes[0] == 3 then argTypes[3] contains the last argType.
* These argTypes can be compared against the symbols declared in parse.y
* eg. ARRAY, POLY, NUMBER etc. There is a new one however. ARRAYADR. This
* is for functions that have the format f(&a[]).
*/
