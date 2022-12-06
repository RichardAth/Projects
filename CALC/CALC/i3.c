/* i3.c */
#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_NONSTDC_NO_DEPRECATE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "integer.h"
#include "fun.h"


#ifdef DEBUG
long int nettbytes = 0;
#endif

/*
For comments on the memory banking algorithm, please see the appropriate
routines later in this file.
Peter Adams,  1/4/1993
*/

void *realloc1(void *ptr, size_t nbytes, long mbytes)
/*
* Reallocates nbytes of memory using realloc().
* mbytes = nett increase in memory. Exits if call to realloc fails.
*/
{
	register char *s;

	if ((s = realloc(ptr, nbytes)) == NULL)
	{
		fprintf(stderr, "out of memory!\n\
		realloc request for %u bytes failed.\n", (USI)nbytes);
		exit(1);

	}

#ifdef DEBUG
	nettbytes += mbytes;
#endif

	return(s); /* returns only if s points to a valid block of memory.*/
}

void *malloc1(size_t nbytes)
/*
* mmalloc -- allocates nbytes of memory using malloc().
* Exits if call to malloc fails.
*/
{
	register char *s;

	if ((s = malloc(nbytes)) == NULL)
	{
		fprintf(stderr, "out of memory!\n\
		malloc request for %lu bytes failed.\n", (USL)nbytes);
		exit(1);

	}
#ifdef DEBUG
	nettbytes += nbytes;
#endif
	return(s); /*returns only if s points to a valid block of memory.*/
}

void *calloc1(size_t m, size_t n)
/* allocates m * n of contiguous memory using calloc().
* Exits if call to calloc fails.
*/
{
	register char *s;

	if ((s = calloc(m, n)) == NULL)
	{
		fprintf(stderr, "out of memory!\n\
		calloc request for %ld bytes failed.\n", (long)(m * n));
		exit(1);

	}
#ifdef DEBUG
	nettbytes += m * n;
#endif

	return(s); /*returns only if s points to a valid block of memory.*/
}

void sfree(void *ptr, size_t nbytes, char *errorfile, int errorline)
/*
* frees memory using free, keeping track of nettbytes.
*/
{
	free(ptr);

#ifdef DEBUG
	nettbytes -= nbytes;
#endif

	return;
}

unsigned int EQONEI(MPI *Mptr)
/*
* returns 1 if *Mptr == 1, 0 otherwise.
*/
{
	if (Mptr->S == 1 && Mptr->D == 0 && Mptr->V[0] == 1)
		return (1);
	else
		return (0);
}

unsigned int EQZEROI(MPI *Mptr)
/*
* returns 1 if *Mptr == 0, 0 otherwise.
*/
{
	if (Mptr->S == 0)
		return (1);
	else
		return (0);
}

unsigned int EQMINUSONEI(MPI *Mptr)
/*
* returns 1 if *Mptr == -1, 0 otherwise.
*/
{
	if (Mptr->S == -1 && Mptr->D == 0 && Mptr->V[0] == 1)
		return (1);
	else
		return (0);
}

unsigned int EQONER(MPR *Mptr)
/*
* returns 1 if *Mptr == 1, 0 otherwise.
*/
{
	if (EQUALI(Mptr->N, Mptr->D))
		return (1);
	else
		return (0);
}

unsigned int EQMINUSONER(MPR *Mptr)
/*
* returns 1 if *Mptr == -1, 0 otherwise.
*/
{
	if (EQMINUSONEI(Mptr->N) * EQONEI(Mptr->D))
		return (1);
	else
		return (0);
}

unsigned int EQZEROR(MPR *Mptr)
/*
* returns 1 if *Mptr == 0, 0 otherwise.
*/
{
	if (Mptr->N->S == 0)
		return (1);
	else
		return (0);
}

unsigned int EQMINUSONECR(MPCR *Mptr)
/*
* returns 1 if *Mptr == -1, 0 otherwise.
*/
{
	if (EQMINUSONER(Mptr->R) && Mptr->I->N->S == 0)
		return (1);
	else
		return (0);
}

unsigned int EQONECR(MPCR *Mptr)
/*
* returns 1 if *Mptr == 1, 0 otherwise.
*/
{
	if (EQONER(Mptr->R) && Mptr->I->N->S == 0)
		return (1);
	else
		return (0);
}

unsigned int EQZEROCR(MPCR *Mptr)
/*
* returns 1 if *Mptr == 0, 0 otherwise.
*/
{
	if (Mptr->R->N->S == 0 && Mptr->I->N->S == 0)
		return (1);
	else
		return (0);
}

unsigned int EQMINUSONECI(MPCI *Mptr)
/*
* returns 1 if *Mptr == -1, 0 otherwise.
*/
{
	if (EQMINUSONEI(Mptr->R) && Mptr->I->S == 0)
		return (1);
	else
		return (0);
}

unsigned int EQONECI(MPCI *Mptr)
/*
* returns 1 if *Mptr == 1, 0 otherwise.
*/
{
	if (EQONEI(Mptr->R) && Mptr->I->S == 0)
		return (1);
	else
		return (0);
}

unsigned int EQZEROCI(MPCI *Mptr)
/*
* returns 1 if *Mptr == 0, 0 otherwise.
*/
{
	if (Mptr->R->S == 0 && Mptr->I->S == 0)
		return (1);
	else
		return (0);
}
