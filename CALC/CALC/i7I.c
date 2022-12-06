/* i7I.c */
#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_NONSTDC_NO_DEPRECATE
#include <stdio.h>
#include "integer.h"
#include "fun.h"

MPI *DOTRI(MPMATI *Mptr, USI i, USI j)
/*
* *Aptr is the dot product of rows i and j in *Mptr.
*/
{
	unsigned int k;
	MPI *H, *Aptr, *TempI;

	Aptr = ZEROI();
	for (k = 0; k <= Mptr->C - 1; k++)
	{
		H = MULTI(Mptr->V[i][k], Mptr->V[j][k]);
		TempI = Aptr;
		Aptr = ADDI(Aptr, H);
		FREEMPI(H);
		FREEMPI(TempI);
	}
	return (Aptr);
}
