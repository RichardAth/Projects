/* integer.h */
#ifndef __INTEGER_H__
#define __INTEGER_H__

#include <stddef.h>
/*
* If DEBUG is defined in the Makefile, the nett change in global variable
* nettbytes is  calculated at the end of the program.  nettbytes should be
* zero at program end if proper allocation and deallocation have been done.
* If investigating a program crash due to suspected misfreeing, turn on the
* Michael Forbes debugging device in int/i3.c by defining FORBES in Makefile.
* (Should only be done on mainframes as it requires a large pointer table.)
* Banking procedures are used to save on the mallocing, reallocing and freeing.
* If DEBUG is not switched on, to save a little time, malloc, realloc and free
* are used on the mainframe version, as space considerations are unnecessary.
* However mmalloc and realloc1 are used in the PC version, as space is at a
* premium and we would like an exit(1) to occur when no space can be allocated.
*/

#define MIN_ARRAY_SIZE 11 /* enough for slots 0 to 10 */
#define MAX_ARRAY_SIZE 65536
#ifdef DEBUG
#define ffree(x,y) sfree(x,y,__FILE__,__LINE__)
#define SAVEPI(x) ffree((char *)x, sizeof(struct TERMI))
#define SAVEPR(x) ffree((char *)x, sizeof(struct TERMR))
#define SAVEPCI(x) ffree((char *)x, sizeof(struct TERMCI))
#define SAVEPCR(x) ffree((char *)x, sizeof(struct TERMCR))
#define SAVEPm(x) ffree((char *)x, sizeof(struct TERMm))
#define rrealloc(x,y,z) realloc1(x,y,z)
#define mmalloc(x) malloc1(x)
#define ccalloc(x,y) calloc1(x,y)
#else
#define SAVEPI SAVEI
#define SAVEPR SAVER
#define SAVEPCI SAVECI
#define SAVEPCR SAVECR
#define SAVEPm SAVEm
#define ffree(x,y) free(x)
#ifdef _WIN32
#define rrealloc(x,y,z) realloc1(x,y,z)
#define mmalloc(x) malloc1(x)
#define ccalloc(x,y) calloc1(x,y)
#else
#define rrealloc(x,y,z) realloc(x,y)
#define mmalloc(x) malloc(x)
#define ccalloc(x,y) calloc(x,y)
#endif
#endif

#define R0 65536 /* base = 2^16 = 65536 */
#define T0 16 /* used in int/i1.c and int/i2.c */
#define Z0 98 /* number of rows in array reciprocal in primes.h */
#define Y0 2048 /* number of primes in the array PRIME in primes.h */
#define O5 5
#define O32 32
#define O31 31
#define USL unsigned long
#define USI unsigned int

#define  get_element(M, row, col) (M[row][col>>5] & ((USL)1<<(31-(col&31)))?1:0)
#define elt(m, r, c) ((m)->V[(r) * (m)->C + (c)])

typedef struct _MPI
{
	unsigned long *V;
	int S; /* S = -1,0,-1, with S=0 corresponding to D=0 and V[0]=0.*/
	unsigned int D; /* length of array V - 1. */
	struct _MPI *NEXT;/* only used for BANKING, when it's = NULL. */
}
MPI;

typedef struct
{
	MPI *N;
	MPI *D;
}
MPR;
/* Here gcd(*N, *D) = 1 and *D > 0. */

typedef struct
{
	MPI *R;
	MPI *I;
}
MPCI;

typedef struct
{
	MPR *R;
	MPR *I;
}
MPCR;

typedef struct
{
	MPI ***V;
	unsigned int R;
	unsigned int C;
}
MPMATI;

typedef struct
{
	MPR **V; /* the vector way of presenting a matrix - used in CMAT */
	unsigned int R;
	unsigned int C;
}
MPMATR;

struct TERMI
{
	int DEG;
	MPI *COEF;
	struct TERMI *NEXT;
};

/* THE MPI array.  Added by Sean Seefried */
typedef struct _MPIA {
	MPI **A;
	unsigned size; /* A defacto standard is 0 <= size <= 2^16-1 */
	unsigned slots; /* The number of slots allocated.  Sometimes different to
					* size. eg. when array built slots = MIN_ARRAY_SIZE
					* while size = 0 */
} *MPIA;

typedef struct TERMI *POLYI;
/* A variable of type POLYI is thus the head pointer to a linear linked-list
* which represents a polynomial with MPI coefficients.
* each structure in the linked list corresponds to a component monomial.
* the zero polynomial is represented by the null pointer NULL.
* Our constructions are based on the account in Scientific Pascal by
* H. Flanders, pp. 175-189.
*/

#endif

