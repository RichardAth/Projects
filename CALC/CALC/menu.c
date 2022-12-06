/* menu.c */
#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_NONSTDC_NO_DEPRECATE
#include <stdio.h>
#include <stdlib.h>
#include "integer.h"
#include "fun.h"

#define MEMINCR 256
#define M0 15

char *Gets()
/*
* As for Fgets() (see below) but from standard input.
*/
{
	return Fgets(stdin);
}

char *Fgets(FILE *f)
/*
* Reads one line from file into a string allocated by malloc.  Return NULL on
* error, otherwise a pointer to the string.
* Byte, June 1988, page 316.
*/
{
	char *r;
	unsigned int n, m;
	int c;

	n = 0;
	m = MEMINCR;
	r = (char *)malloc((m + 1) * sizeof(unsigned int));
	while ((c = fgetc(f)) != '\n') {
		if (c == EOF)
			return NULL;
		if (--m == 0) {
			if ((r = (char*)realloc(r, n + MEMINCR + 1)) == NULL)
				return NULL;
			m = MEMINCR;
		}
		r[n++] = c;
	}
	r[n] = '\0';
	if ((r = (char*)realloc(r, n + 1)) == NULL)
		return NULL;
	return r;
}

void SelOpt()
/*
* This function simply prints the "SELECT OPTION: " prompt.  It is here only
* for consistency, ie whenever the words SELECT OPTION are printed it can be
* guaranteed that they are coming from this function, and the way they are
* printed will be exactly the same every time.
*/
{
	printf("\nSELECT OPTION: ");
}

unsigned int GetYN()
/*
* Gets a character in from the keyboard, making sure it's a y or an n (either
* case).  If at first the user doesn't succeed, he tries, tries again.
* 0 is returned if n, 1 if y.
*/
{
	int c;

	printf("Enter y or n : ");
	c = GetCharFlush();
	while (c != 'y' && c != 'Y' && c != 'n' && c != 'N') {
		printf("Try again. Enter y or n: ");
		c = GetCharFlush();
	}
	return c == 'y' || c == 'Y';
}

int whitespace(FILE *f)
/* removes spaces or newlines from a stream. */
{
	int ch;

	while ((ch = fgetc(f)) != EOF && (ch == ' ' || ch == '\n'))
		;
	return ch;
}

void FFlush(FILE *f)
{
	int ch;
	while ((ch = fgetc(f)) != EOF && ch != '\n')
		;
}

void Flush()
{
	int ch;
	while ((ch = getchar()) != EOF && ch != '\n')
		;
}

int GetCharFlush()
{
	char ch;
	int s;
	s = scanf(" %c", &ch);
	Flush();
	return ch;
}

void GetReturn()
{
	printf("HIT RETURN: ");
	Flush();
}

void TryAgain()
{
	printf("Hit return and Try again: ");
	Flush();
}


void clearscreen()
/* Calls Unix to clear the screen. */
{
	system("/usr/ucb/clear");
}

int scanf01(USI *jptr)
/*
* returns 1 if unsigned int *jptr < M0 is successfully entered, 0 otherwise.
*/
{
	int n;

	n = scanf("%u", jptr);
	Flush();
	if (n == 1)
	{
		if (*jptr < M0)
			return 1;
		else
			printf("integer %u exceeded array bound %u, try again:\n", *jptr, M0 - 1);
	}
	else
		printf("try again:\n");
	GetReturn();
	return 0;
}

int scanf011(USI *jptr, USI *kptr)
/*
* returns 1 if unsigned int *jptr < M0 and unsigned int *kptr are successfully
* entered, 0 otherwise.
*/
{
	int n;

	n = scanf("%u%u", jptr, kptr);
	Flush();
	if (n == 2)
	{
		if (*jptr < M0)
			return 1;
		else
			printf("first integer %u exceeded array bound %u, try again:\n", *jptr, M0 - 1);
	}
	else
		printf("try again:\n");
	GetReturn();
	return 0;
}

int scanf02(USI *jptr, USI *kptr)
/*
* returns 1 if unsigned ints *jptr, *kptr < M0 are successfully entered,
* 0 otherwise.
*/
{
	int n;

	n = scanf("%u%u", jptr, kptr);
	Flush();
	if (n == 2)
	{
		if (*jptr < M0 && *kptr < M0)
			return 1;
		else
			printf("first integer %u or second integer %u exceeded array bound %u, try again:\n", *jptr, *kptr, M0 - 1);
	}
	else
		printf("try again:\n");
	GetReturn();
	return 0;
}

int scanf021(USI *jptr, USI *kptr, USI *mptr)
/*
* returns 1 if unsigned ints *jptr, *kptr, each less than M0
* and unsigned int *mptr are successfully entered, 0 otherwise.
*/
{
	int n;

	n = scanf("%u%u%u", jptr, kptr, mptr);
	Flush();
	if (n == 3)
	{
		if (*jptr < M0 && *kptr < M0)
			return 1;
		else
			printf("first integer %u or second integer %u exceeded array bound %u, try again:\n", *jptr, *kptr, M0 - 1);
	}
	else
		printf("try again:\n");
	GetReturn();
	return 0;
}

int scanf21(USI *jptr, USI *kptr, USI *mptr)
/*
* returns 1 if unsigned ints *jptr, *kptr, with *mptr < M0 are successfully
* entered, 0 otherwise.
*/
{
	int n;

	n = scanf("%u%u%u", jptr, kptr, mptr);
	Flush();
	if (n == 3)
	{
		if (*mptr < M0)
			return 1;
		else
			printf("third integer %u exceeded array bound %u, try again:\n", *mptr, M0 - 1);
	}
	else
		printf("try again:\n");
	GetReturn();
	return 0;
}

int scanf121(USI *hptr, USI *jptr, USI *kptr, USI *mptr)
/*
* returns 1 if arbitrary unsigned ints *hptr, *mptr and
* unsigned ints *jptr, *kptr, each less than M0 are successfully entered,
* 0 otherwise.
*/
{
	int n;

	n = scanf("%u%u%u%u", hptr, jptr, kptr, mptr);
	Flush();
	if (n == 4)
	{
		if (*jptr < M0 && *kptr < M0)
			return 1;
		else
			printf("second integer %u or third integer %u exceeded array bound %u, try again:\n", *jptr, *kptr, M0 - 1);
	}
	else
		printf("try again:\n");
	GetReturn();
	return 0;
}

int scanf221(USI *hptr, USI *iptr, USI *jptr, USI *kptr, USI *mptr)
/*
* returns 1 if unsigned ints *hptr, *iptr, *mptr and
* unsigned ints *jptr, *kptr, each less than M0 are successfully entered,
* 0 otherwise.
*/
{
	int n;

	n = scanf("%u%u%u%u%u", hptr, iptr, jptr, kptr, mptr);
	Flush();
	if (n == 5)
	{
		if (*jptr < M0 && *kptr < M0)
			return 1;
		else
			printf("third integer %u or fourth integer %u exceeded array bound %u, try again:\n", *jptr, *kptr, M0 - 1);
	}
	else
		printf("try again:\n");
	GetReturn();
	return 0;
}

int scanf321(USI *gptr, USI *hptr, USI *iptr, USI *jptr, USI *kptr, USI *mptr)
/*
* returns 1 if unsigned ints *gptr, *hptr, *iptr, *mptr and
* unsigned ints *jptr, *kptr, each less than M0 are successfully entered,
* 0 otherwise.
*/
{
	int n;

	n = scanf("%u%u%u%u%u%u", gptr, hptr, iptr, jptr, kptr, mptr);
	Flush();
	if (n == 6)
	{
		if (*jptr < M0 && *kptr < M0)
			return 1;
		else
			printf("fourth integer %u or fifth integer %u exceeded array bound %u, try again:\n", *jptr, *kptr, M0 - 1);
	}
	else
		printf("try again:\n");
	GetReturn();
	return 0;
}

int scanf03(USI *jptr, USI *kptr, USI *lptr)
/*
* returns 1 if unsigned ints *jptr, *kptr, *lptr, each less than M0 are
* successfully entered, 0 otherwise.
*/
{
	int n;

	n = scanf("%u%u%u", jptr, kptr, lptr);
	Flush();
	if (n == 3)
	{
		if (*jptr < M0 && *kptr < M0 && *lptr < M0)
			return 1;
		else
			printf("one of %u, %u and %u exceeded array bound %u, try again:\n", *jptr, *kptr, *lptr, M0 - 1);
	}
	else
		printf("try again:\n");
	GetReturn();
	return 0;
}

int scanf031(USI *jptr, USI *kptr, USI *lptr, USI *mptr)
/*
* returns 1 if unsigned ints *jptr, *kptr, *lptr, each less than M0 and
* unsigned int *mptr are successfully entered, 0 otherwise.
*/
{
	int n;

	n = scanf("%u%u%u%u", jptr, kptr, lptr, mptr);
	Flush();
	if (n == 4)
	{
		if (*jptr < M0 && *kptr < M0 && *lptr < M0)
			return 1;
		else
			printf("one of %u, %u and %u exceeded array bound %u, try again:\n", *jptr, *kptr, *lptr, M0 - 1);
	}
	else
		printf("try again:\n");
	GetReturn();
	return 0;
}

int scanf131(USI *iptr, USI *jptr, USI *kptr, USI *lptr, USI *mptr)
/*
* returns 1 if unsigned ints *jptr, *kptr, *lptr, each less than M0 and
* unsigned int *mptr, *iptr are successfully entered, 0 otherwise.
*/
{
	int n;

	n = scanf("%u%u%u%u%u", iptr, jptr, kptr, lptr, mptr);
	Flush();
	if (n == 5)
	{
		if (*jptr < M0 && *kptr < M0 && *lptr < M0)
			return 1;
		else
			printf("one of %u, %u and %u exceeded array bound %u, try again:\n", *jptr, *kptr, *lptr, M0 - 1);
	}
	else
		printf("try again:\n");
	GetReturn();
	return 0;
}

int scanf05(USI *iptr, USI *jptr, USI *kptr, USI *lptr, USI *mptr)
/*
* returns 1 if unsigned ints *iptr, *jptr, *kptr, *lptr, *mptr with
* *iptr, *jptr, *kptr, *lptr and *mptr < M0 are successfully entered,
* 0 otherwise.
*/
{
	int n;

	n = scanf("%u%u%u%u%u", iptr, jptr, kptr, lptr, mptr);
	Flush();
	if (n == 5)
	{
		if (*iptr < M0 && *jptr < M0 && *kptr < M0 && *lptr < M0 && *mptr < M0)
			return 1;
		else
			printf("one of %u, %u, %u, %u and %u exceeded array bound %u, try again:\n", *iptr, *jptr, *kptr, *lptr, *mptr, M0 - 1);
	}
	else
		printf("try again:\n");
	GetReturn();
	return 0;
}

int scanf051(USI *iptr, USI *jptr, USI *kptr, USI *lptr, USI *mptr, USI *nptr)
/*
* returns 1 if unsigned ints *iptr, *jptr, *kptr, *lptr, *mptr, *nptr with,
* *iptr, *jptr, *kptr, *lptr and *mptr < M0 are successfully entered,
* 0 otherwise.
*/
{
	int n;

	n = scanf("%u%u%u%u%u%u", iptr, jptr, kptr, lptr, mptr, nptr);
	Flush();
	if (n == 6)
	{
		if (*iptr < M0 && *jptr < M0 && *kptr < M0 && *lptr < M0 && *mptr < M0)
			return 1;
		else
			printf("one of %u, %u, %u, %u and %u exceeded array bound %u, try again:\n", *iptr, *jptr, *kptr, *lptr, *mptr, M0 - 1);
	}
	else
		printf("try again:\n");
	GetReturn();
	return 0;
}

int scanf10(USI *jptr)
/*
* returns 1 if unsigned int *jptr successfully, 0 otherwise.
*/
{
	int n;

	n = scanf("%u", jptr);
	Flush();
	if (n == 1)
		return 1;
	else
		printf("try again:\n");
	GetReturn();
	return 0;
}

int scanfd10(int *jptr)
/*
* returns 1 if int *jptr successfully entered, 0 otherwise.
*/
{
	int n;

	n = scanf("%d", jptr);
	Flush();
	if (n == 1)
		return 1;
	else
		printf("try again:\n");
	GetReturn();
	return 0;
}

int scanf11(USI *jptr, USI *kptr)
/*
* returns 1 if  unsigned ints *jptr, *kptr with *kptr < M0 are successfully
* entered, 0 otherwise.
*/
{
	int n;

	n = scanf("%u%u", jptr, kptr);
	Flush();
	if (n == 2)
	{
		if (*kptr < M0)
			return 1;
		else
			printf("second integer %u exceeded array bound %u, try again:\n", *kptr, M0 - 1);
	}
	else
		printf("try again:\n");
	GetReturn();
	return 0;
}

int scanf12(USI *jptr, USI *kptr, USI *lptr)
/*
* returns 1 if unsigned ints *jptr, *kptr, *lptr with *kptr and *lptr < M0
* are successfully entered, 0 otherwise.
*/
{
	int n;

	n = scanf("%u%u%u", jptr, kptr, lptr);
	Flush();
	if (n == 3)
	{
		if (*kptr < M0 && *lptr < M0)
			return 1;
		else
			printf("second integer %u or third integer %u exceeded array bound %u, try again:\n", *kptr, *lptr, M0 - 1);
	}
	else
		printf("try again:\n");
	GetReturn();
	return 0;
}

int scanf012(USI *jptr, USI *kptr, USI *lptr)
/*
* returns 1 if unsigned ints *jptr, *kptr, *lptr with *jptr < M0 are
* successfully entered, 0 otherwise.
*/
{
	int n;

	n = scanf("%u%u%u", jptr, kptr, lptr);
	Flush();
	if (n == 3)
	{
		if (*jptr < M0)
			return 1;
		else
			printf("first integer %u exceeded array bound %u, try again:\n", *jptr, M0 - 1);
	}
	else
		printf("try again:\n");
	GetReturn();
	return 0;
}

int scanf0112(USI *jptr, USI *kptr, USI *lptr, USI *mptr)
/*
* returns 1 if unsigned ints *jptr, *kptr, *lptr, *mptr with
* jptr < M0, *lptr < M0, *mptr < M0 are successfully entered, 0 otherwise.
*/
{
	int n;

	n = scanf("%u%u%u%u", jptr, kptr, lptr, mptr);
	Flush();
	if (n == 4)
	{
		if (*jptr < M0 && *lptr < M0 && *mptr < M0)
			return 1;
		else
			printf("one of %u, %u and %u exceeded array bound %u, try again:\n", *jptr, *lptr, *mptr, M0 - 1);
	}
	else
		printf("try again:\n");
	GetReturn();
	return 0;
}

int scanf0122(USI *jptr, USI *kptr, USI *lptr, USI *mptr, USI *nptr)
/*
* returns 1 if unsigned ints *jptr, *kptr, *lptr, *mptr, *nptr with
* *jptr < M0, *mptr < M0, *nptr < M0 are successfully entered, 0 otherwise.
*/
{
	int n;

	n = scanf("%u%u%u%u%u", jptr, kptr, lptr, mptr, nptr);
	Flush();
	if (n == 5)
	{
		if (*jptr < M0 && *mptr < M0 && *nptr < M0)
			return 1;
		else
			printf("one of %u, %u and %u exceeded array bound %u, try again:\n", *jptr, *mptr, *nptr, M0 - 1);
	}
	else
		printf("try again:\n");
	GetReturn();
	return 0;
}

int scanf01221(USI *jptr, USI *kptr, USI *lptr, USI *mptr, USI *nptr, USI *optr)
/*
* returns 1 if unsigned ints *jptr, *kptr, *lptr, *mptr, *nptr, *optr
* with *jptr < M0, *mptr < M0, *nptr < M0 successfully entered, 0 otherwise.
*/
{
	int n;

	n = scanf("%u%u%u%u%u%u", jptr, kptr, lptr, mptr, nptr, optr);
	Flush();
	if (n == 6)
	{
		if (*jptr < M0 && *mptr < M0 && *nptr < M0)
			return 1;
		else
			printf("one of %u, %u and %u exceeded array bound %u, try again:\n", *jptr, *mptr, *nptr, M0 - 1);
	}
	else
		printf("try again:\n");
	GetReturn();
	return 0;
}

int scanf22(USI *iptr, USI *jptr, USI *kptr, USI *lptr)
/*
* returns 1 if unsigned ints *iptr, *jptr, *kptr, *lptr with *kptr and
* *lptr < M0 are successfully entered, 0 otherwise.
*/
{
	int n;

	n = scanf("%u%u%u%u", iptr, jptr, kptr, lptr);
	Flush();
	if (n == 4)
	{
		if (*kptr < M0 && *lptr < M0)
			return 1;
		else
			printf("third integer %u or fourth integer %u exceeded array bound %u, try again:\n", *kptr, *lptr, M0 - 1);
	}
	else
		printf("try again:\n");
	GetReturn();
	return 0;
}

int scanf20(USI *iptr, USI *jptr)
/*
* returns 1 if unsigned ints *iptr, *jptr successfully entered, 0 otherwise.
*/
{
	int n;

	n = scanf("%u%u", iptr, jptr);
	Flush();
	if (n == 2)
		return 1;
	else
		printf("try again:\n");
	GetReturn();
	return 0;
}

int scanf30(USI *iptr, USI *jptr, USI *kptr)
/*
* returns 1 if unsigned ints *iptr, *jptr, *kptr are succesfully entered,
* 0 otherwise.
*/
{
	int n;

	n = scanf("%u%u%u", iptr, jptr, kptr);
	Flush();
	if (n == 3)
		return 1;
	else
		printf("try again:\n");
	GetReturn();
	return 0;
}

int scanf40(USI *iptr, USI *jptr, USI *kptr, USI *lptr)
/*
* returns 1 if unsigned ints *iptr, *jptr, *kptr, *lptr are successfully
* entered, 0 otherwise.
*/
{
	int n;

	n = scanf("%u%u%u%u", iptr, jptr, kptr, lptr);
	Flush();
	if (n == 4)
		return 1;
	else
		printf("try again:\n");
	GetReturn();
	return 0;
}

int scanf23(USI *iptr, USI *jptr, USI *kptr, USI *lptr, USI *mptr)
/*
* returns 1 if unsigned ints *iptr, *jptr, *kptr, *lptr, *mptr with
* *kptr, *lptr and *mptr < M0 are successfully entered, 0 otherwise.
*/
{
	int n;

	n = scanf("%u%u%u%u%u", iptr, jptr, kptr, lptr, mptr);
	Flush();
	if (n == 5)
	{
		if (*kptr < M0 && *lptr < M0 && *mptr < M0)
			return 1;
		else
			printf("one of %u, %u and %u exceeded array bound %u, try again:\n", *kptr, *lptr, *mptr, M0 - 1);
	}
	else
		printf("try again:\n");
	GetReturn();
	return 0;
}

int scanf231(USI *iptr, USI *jptr, USI *kptr, USI *lptr, USI *mptr, USI *nptr)
/*
* returns 1 if unsigned ints *iptr, *jptr, *kptr, *lptr, *mptr, *nptr with
* *kptr, *lptr and *mptr < M0 are successfully entered, 0 otherwise.
*/
{
	int n;

	n = scanf("%u%u%u%u%u%u", iptr, jptr, kptr, lptr, mptr, nptr);
	Flush();
	if (n == 6)
	{
		if (*kptr < M0 && *lptr < M0 && *mptr < M0)
			return 1;
		else
			printf("one of %u, %u and %u exceeded array bound %u, try again:\n", *kptr, *lptr, *mptr, M0 - 1);
	}
	else
		printf("try again:\n");
	GetReturn();
	return 0;
}

int scanf13(USI *iptr, USI *jptr, USI *kptr, USI *lptr)
/*
* returns 1 if unsigned ints *iptr, *jptr, *kptr, *lptr with
* *jptr, *kptr and *lptr < M0 are successfully entered, 0 otherwise.
*/
{
	int n;

	n = scanf("%u%u%u%u", iptr, jptr, kptr, lptr);
	Flush();
	if (n == 4)
	{
		if (*jptr < M0 && *kptr < M0 && *lptr < M0)
			return 1;
		else
			printf("one of %u, %u and %u exceeded array bound %u, try again:\n", *jptr, *kptr, *lptr, M0 - 1);
	}
	else
		printf("try again:\n");
	GetReturn();
	return 0;
}

FILE *Fopen_r(char *filename)
/*
* An attempt is made at opening the file specified for reading. If successful,
* a pointer to the file is returned, otherwise a null pointer is returned.
*/
{
	FILE *f;
	if ((f = fopen(filename, "r")) == NULL)
	{
		fprintf(stderr, "cmat: ");
		perror(filename);
		exit(1);
	}
	return f;
}
