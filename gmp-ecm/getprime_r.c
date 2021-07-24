/* Dynamic Eratosthenes sieve.
 
  Copyright 2001-2016 Paul Zimmermann and Alexander Kruppa.
  Imported from CADO-NFS, which imported it from GMP-ECM.
  (use 'unsigned long' instead of 'double'; thread-safe)

  This file is part of the ECM Library.

  The ECM Library is free software; you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation; either version 2.1 of the License, or (at your
  option) any later version.

  The ECM Library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
  License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with the ECM Library; see the file COPYING.LIB.  If not, write to
  the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,
  MA 02110-1301, USA.
*/

/* compile with -DMAIN to use as a standalone program */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "getprime_r.h"

/* provided for in cado.h, but we want getprime.c to be standalone */
#ifndef ASSERT
#define ASSERT(x)
#endif

/* This function returns successive odd primes, starting with 3.
   To perform a loop over all primes <= B1, do the following
   (compile this file with -DMAIN to count primes):

      prime_info_t pi;
      prime_info_init (pi);
      for (p = 2; p <= B1; p = getprime_mt (pi))
         {
            ...
         }

      prime_info_clear (pi);
*/

void
prime_info_init (prime_info_t i)
{
  i->offset = 0;
  i->current = -1;
  i->primes = NULL;
  i->nprimes = 0;
  i->sieve = NULL;
  i->len = 0;
  i->moduli = NULL;
}

void
prime_info_clear (prime_info_t i)
{
  free (i->primes);
  free (i->sieve);
  free (i->moduli);
}

/* this function is thread-safe */
ecm_uint
getprime_mt (prime_info_t i)
{
  if (i->len)
    {
      unsigned char *ptr = i->sieve + i->current;
      while (!*++ptr);
      i->current = ptr - i->sieve;
    }
  else
    i->current = 0;

  if (i->current < i->len) /* most calls will end here */
    return i->offset + 2 * i->current;

  /* otherwise we have to sieve */
  i->offset += 2 * i->len;

  /* first enlarge sieving table if too small */
  if ((ecm_uint) i->len * i->len < i->offset)
    {
      free (i->sieve);
      i->len *= 2;
      i->sieve = (unsigned char *) malloc ((i->len + 1 )
                                           * sizeof (unsigned char));
      /* assume this "small" malloc will not fail in normal usage */
      ASSERT(i->sieve != NULL);
      i->sieve[i->len] = 1; /* End mark */
    }

  /* now enlarge small prime table if too small */
  if ((i->nprimes == 0) ||
      ((ecm_uint) i->primes[i->nprimes - 1] * (ecm_uint)
       i->primes[i->nprimes - 1] < i->offset + i->len))
      {
	if (i->nprimes == 0) /* initialization */
	  {
	    i->nprimes = 1;
	    i->primes = (ecm_uint*) malloc (i->nprimes
                                                * sizeof(ecm_uint));
	    /* assume this "small" malloc will not fail in normal usage */
	    ASSERT(i->primes != NULL);
	    i->moduli = (ecm_uint*) malloc (i->nprimes
                                                * sizeof(ecm_uint));
	    /* assume this "small" malloc will not fail in normal usage */
	    ASSERT(i->moduli != NULL);
	    i->len = 1;
	    i->sieve = (unsigned char *) malloc((i->len + 1) *
                                       sizeof(unsigned char)); /* len=1 here */
	    /* assume this "small" malloc will not fail in normal usage */
	    ASSERT(i->sieve != NULL);
	    i->sieve[i->len] = 1; /* End mark */
	    i->offset = 5;
	    i->sieve[0] = 1; /* corresponding to 5 */
	    i->primes[0] = 3;
	    i->moduli[0] = 1; /* next odd multiple of 3 is 7, i.e. next to 5 */
	    i->current = -1;
	    return 3;
	  }
	else
	  {
	    ecm_uint k, p, j, ok;

	    k = i->nprimes;
	    i->nprimes *= 2;
	    i->primes = (ecm_uint*) realloc (i->primes, i->nprimes *
                                                 sizeof(ecm_uint));
	    i->moduli = (ecm_uint*) realloc (i->moduli, i->nprimes *
                                                 sizeof(ecm_uint));
	    /* assume those "small" realloc's will not fail in normal usage */
	    ASSERT(i->primes != NULL && i->moduli != NULL);
	    for (p = i->primes[k-1]; k < i->nprimes; k++)
	      {
		/* find next (odd) prime > p */
		do
		  {
		    for (p += 2, ok = 1, j = 0; (ok != 0) && (j < k); j++)
		      ok = p % i->primes[j];
		  }
		while (ok == 0);
		i->primes[k] = p;
		/* moduli[k] is the smallest m such that offset + 2*m = k*p */
		j = i->offset % p;
		j = (j == 0) ? j : p - j; /* -offset mod p */
		if ((j % 2) != 0)
		  j += p; /* ensure j is even */
		i->moduli[k] = j / 2;
	      }
	  }
      }

  /* now sieve for new primes */
  {
    ecm_int k;
    ecm_uint j, p;
    
    memset (i->sieve, 1, sizeof(unsigned char) * (i->len + 1));
    for (j = 0; j < i->nprimes; j++)
      {
	p = i->primes[j];
	for (k = i->moduli[j]; k < i->len; k += p)
	  i->sieve[k] = 0;
	i->moduli[j] = k - i->len; /* for next sieving array */
      }
  }

  {
    unsigned char *ptr = i->sieve - 1;
    while (!*++ptr);
    i->current = ptr - i->sieve;
  }

  ASSERT(i->current < i->len); /* otherwise we found a prime gap >= sqrt(x)
                                  around x */
  return i->offset + 2 * i->current;
}

#ifdef MAIN
int
main (int argc, char *argv[])
{
  unsigned long p, B;
  unsigned long pi = 0;
  prime_info_t i;

  if (argc != 2)
    {
      fprintf (stderr, "Usage: getprime <bound>\n");
      exit (EXIT_FAILURE);
    }

  B = strtoul (argv[1], NULL, 0);
  
  prime_info_init (i);

  for (pi = 0, p = 2; p <= B; p = getprime_mt (i), pi++);
  printf ("pi(%lu)=%lu\n", B, pi);

  prime_info_clear (i); /* free the tables */

  return 0;
}
#endif
