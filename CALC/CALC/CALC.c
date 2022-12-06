/* binary.c */

#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_NONSTDC_NO_DEPRECATE
#include <stdio.h>
#include "integer.h"
#include "fun.h"

/*
main()
{
unsigned long i, j;
int s;
printf("Enter an unsigned int: ");
s=scanf("%lu", &i);
printf("%u\n", binary(i));
*/


/*

for (i = 0; i < 16; i++)
{
for (j = 0; j < 16; j++)
{
if (i != 0 || j != 0)
printf("%u,", BINARY(16 * i + j));
}
printf("\n");
}
}

*/
unsigned int binary(USL i)
/*
* returns the number of binary digits of i.
*/
{
	static unsigned int bin[256] = {
		0,1,2,2,3,3,3,3,4,4,4,4,4,4,4,4,
		5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
		6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
		6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
		7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,
		7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,
		7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,
		7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,
		8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,
		8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,
		8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,
		8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,
		8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,
		8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,
		8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,
		8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8 };

	unsigned int j;

	if ((j = (i >> 24)))
		return (24 + bin[j]);
	else if ((j = (i >> 16)))
		return (16 + bin[j]);
	else if ((j = (i >> 8)))
		return (8 + bin[j]);
	else /* 0 <= i < 255 */
		return (bin[i]);
}

unsigned int BINARY(USL i)
/*
* returns the location of the rightmost bitset of i (non-zero).
*/
{
	static unsigned int BIN[255] = {
		0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
		4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
		5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
		4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
		6,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
		4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
		5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
		4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
		7,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
		4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
		5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
		4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
		6,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
		4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
		5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
		4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0 };

	unsigned int j;

	if ((j = (unsigned int)(i & (USL)255)))
		return (BIN[j - 1]);
	else if ((j = (unsigned int)((i >> 8) & (USL)255)))
		return (8 + BIN[j - 1]);
	else if ((j = (unsigned int)((i >> 16) & (USL)255)))
		return (16 + BIN[j - 1]);
	else
	{
		j = (unsigned int)((i >> 24) & (USL)255);
		return (24 + BIN[j - 1]);
	}
}

unsigned int BINARYB(MPI *N)
/*
* returns the number of binary digits of N.
*/
{
	return (T0 * (N->D) + binary(N->V[N->D]));
}
