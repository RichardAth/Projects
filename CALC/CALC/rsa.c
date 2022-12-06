/* program: rsa.c */

#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_NONSTDC_NO_DEPRECATE
#include <stdio.h>
#include <string.h>
#include "integer.h"
#include "fun.h"

void doit(FILE *f, char *valstr, MPI *Eptr, MPI *Rptr)
/*
* The string consisting of up to 4 juxtaposed ascii numbers is sent
* to file f.
*/
{
	char *s;
	MPI *Temp, *Mptr;
	unsigned int n;

	printf("juxtaposition of 4 ascii numbers :%s\n", valstr);
	s = valstr;
	Mptr = ZEROI();
	while (*s != '\0')
	{
		Temp = Mptr;
		Mptr = MULT_I(Temp, 10);
		FREEMPI(Temp);
		n = (unsigned int)(*s - '0');
		Temp = Mptr;
		Mptr = ADD0_I(Temp, n);
		FREEMPI(Temp);
		s++;
	}
	Temp = MPOWER(Mptr, Eptr, Rptr);
	printf("encoded number : "); PRINTI(Temp); printf("\n");
	FPRINTI(f, Temp);
	fprintf(f, "\n");
	FREEMPI(Temp);
	FREEMPI(Mptr);
}

void undoit(FILE *f, char *unencstr, MPI *Dptr, MPI *Rptr)
/*
* The numbers in file "encoded.out" are inputted as strings.
* Each string is decomposed into up to 4 ascii numbers.
* The corresponding characters are concatenated to reconstruct the
* original message.
*/
{
	char *s, *s1;
	MPI *Temp, *Mptr;
	unsigned int u;
	char tmpstr[100], valstr[500];
	int newval;

	while (1)
	{
		Mptr = FINPUTI(f, &u);
		if (Mptr->S == -1 && Mptr->V[0] == 1)
		{
			FREEMPI(Mptr);
			break;
		}
		if (Mptr->S == -1 && Mptr->V[0] == 2)
		{
			FREEMPI(Mptr);
			strcat(unencstr, "\n");
			continue;
		}
		printf("Encoded number =  :"); PRINTI(Mptr); printf("\n");
		Temp = MPOWER(Mptr, Dptr, Rptr);
		FREEMPI(Mptr);
		SPRINTI(valstr, Temp);
		FREEMPI(Temp);
		printf("decoded juxtaposition of 4 ascii numbers = %s\n", valstr);
		newval = 0;
		s = valstr;
		s1 = tmpstr;
		while (*s != '\0')
		{
			newval = newval * 10 + (*s) - '0';
			if (newval>31)
			{
				printf("ascii number %d\n", newval);
				*s1++ = newval;
				newval = 0;
			}
			s++;
		}
		if (newval)
		{
			printf("ascii number %d\n", newval);
			*s1++ = newval;
		}
		*s1 = '\0';
		strcat(unencstr, tmpstr);
	}
	return;
}

void ENCODE(MPI *Eptr, MPI *Rptr)
/*
* *Rptr is the encryption modulus, *Eptr is the enciphering key.
* The user is prompted to input a string of <= 500 characters (not control
* characters).  Each character is converted to its ASCII equivalent.
* The resulting numbers lie in the range 32-126. They are concatenated in
* blocks of 4. Each block is then raised to the exponent *Eptr (mod *Rptr)
* and sent to the file "encoded.out", with one number per line.
* 355143^2 > 126126126126 > 355142^2. So *Rptr is assumed to be composed of
* two primes, each greater than 355142.
*/
{
	char str[500], valstr[50], tmp[50], em[20];
	char *cp;
	int ctr, n;
	FILE * outfile;
	MPI *Temp, *Temp1;
	int resp;
	char junk[200];
	char fname[100];
	int ok, s;
	FILE *ifile;

	ifile = 0;
	strcpy(em, "encoded.out");
	n = 0;
	do
	{
		printf("Please type 1 for keyboard input, 2 for file....");
		s = scanf("%d", &resp);
		(void)fgets(junk, 200, stdin);	 /* read in the newline and any other junk */
	} while (resp != 1 && resp != 2);
	if (resp == 2)
	{
		do
		{
			printf("What is the name of the file to be encoded? ");
			s = scanf("%s", fname);
			(void)fgets(junk, 200, stdin);	 /* read in the newline and any other junk */
			if (!(ifile = fopen(fname, "r")))
				ok = 0;
			else
				ok = 1;
			if (!ok)
				printf("Missing file %s!  Please try again.\n", fname);
		} while (!ok);
	}
	outfile = fopen(em, "w");
	while (1)
	{
		if (resp == 1)
		{
			printf("Please type some text, followed by return. \n");
			fgets(str, 500, stdin);
		}
		else
		{
			if (!fgets(str, 500, ifile))
				break;
			cp = str;
			while (*cp && *cp != '\n')
				cp++;
			*cp = 0;
		}
		cp = str;
		ctr = 0;
		valstr[0] = '\0';
		while (*cp)
		{
			sprintf(tmp, "%d", *cp);
			strcat(valstr, tmp);
			ctr++;
			if (ctr == 4)
			{
				doit(outfile, valstr, Eptr, Rptr);
				n++;
				ctr = 0;
				valstr[0] = '\0';
			}
			cp++;
		}
		if (ctr)
		{
			doit(outfile, valstr, Eptr, Rptr);
			n++;
		}
		Temp = MINUS_ONEI();
		Temp1 = ADDI(Temp, Temp);
		FPRINTI(outfile, Temp1);
		fprintf(outfile, "\n");
		FREEMPI(Temp);
		FREEMPI(Temp1);
		if (resp == 1)
			break;
	}
	Temp = MINUS_ONEI();
	FPRINTI(outfile, Temp);
	fprintf(outfile, "\n");
	FREEMPI(Temp);
	fclose(outfile);
	printf("counter = %u\n", n);
	printf("encoded message sent to file 'encoded.out'\n");
	return;
}

void DECODE(MPI *Dptr, MPI *Rptr)
/*
* Inputs each number from the file "encoded.out", raises it to the
* deciphering modulus *Dptr (mod *Rptr) and splits up the resulting
* number into 4 ASCII numbers. The decrypted message is printed on the
* screen, as well as being sent to a file called "decoded.out".
*/
{
	FILE * outfile, *infile;
	char unencstr[500], dm[20], em[20];

	strcpy(em, "encoded.out");
	infile = fopen(em, "r");
	unencstr[0] = 0;
	printf("about to undoit\n");
	undoit(infile, unencstr, Dptr, Rptr);
	fclose(infile);
	strcpy(dm, "decoded.out");
	printf("Unencoded string (sent to file 'decoded.out'):\n%s\n", unencstr);
	outfile = fopen(dm, "w");
	fprintf(outfile, "%s", unencstr);
	fclose(outfile);
	return;
}

MPI *RSAEX(MPI *Pptr, MPI *Qptr)
/* Checks that p and q are primes > 355142, before calling RSAE below. */
{
	MPI *X, *Y;

	X = CHANGE(355142UL);
	if (RSV(Pptr, X) <= 0 || RSV(Qptr, X) <= 0) {
		FREEMPI(X);
		return(NULL);
	}
	FREEMPI(X);
	X = LUCASX(Pptr);
	Y = LUCASX(Qptr);
	if (X == NULL || Y == NULL) {
		return(NULL);
	}
	else {
		return(RSAE(Pptr, Qptr));
	}
}

MPI *RSAE(MPI *Pptr, MPI *Qptr)
/*
* The least integer e which is relatively prime to (p-1)*(q-1) and
* which satisfies 32^e >= 16*((pq)->D+1) is returned as an MPI*.
*/
{
	MPI *Rptr, *Tmp1, *Tmp2, *Tmp3, *Eptr, *G;
	USI e = 3;

	Rptr = MULTI(Pptr, Qptr);
	Tmp1 = SUB0_I(Pptr, 1);
	Tmp2 = SUB0_I(Qptr, 1);
	Tmp3 = MULTI(Tmp1, Tmp2);
	FREEMPI(Tmp1);
	FREEMPI(Tmp2);
	while (1)
	{
		Eptr = CHANGE(e);
		G = GCD(Eptr, Tmp3);
		if (EQONEI(G) && (5 * e >= 16 * (1 + Rptr->D)))
		{
			FREEMPI(G);
			break;
		}
		FREEMPI(G);
		FREEMPI(Eptr);
		e = e + 2;
	}
	FREEMPI(Tmp3);
	FREEMPI(Rptr);
	return (Eptr);
}
