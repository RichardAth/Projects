/* Implements algorithm polyeval and remainder tree using middle product.

Copyright 2003, 2004, 2005, 2006, 2007, 2008, 2009 Laurent Fousse,
Alexander Kruppa, Paul Zimmermann.

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

#define _CRT_SECURE_NO_DEPRECATE

#include <stdlib.h>
#include <string.h> /* for strlen */
#include "ecm-impl.h"

#ifdef HAVE_UNISTD_H
# include <unistd.h> /* for unlink */
#endif


#ifndef MAX
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#endif

/* #define DEBUG_TREEDATA */

extern unsigned int Fermat;

#if defined(DEBUG) || defined(DEBUG_TREEDATA)
void
print_vect (listz_t t, unsigned int l)
{
    unsigned int i;

    fprintf (ECM_STDOUT, "[");
    for (i = 0; i < l; i++)
    {
        mpz_out_str (ECM_STDOUT, 10, t[i]);
        if (i != l - 1)
            fprintf (ECM_STDOUT, ", ");
        else
            fprintf (ECM_STDOUT, "]");
    }
}
#endif

/* Computes TUpTree as described in ref[1]. k is the degree of the
 * polynomial at the root of the tree. sh is the shift we need to
 * apply to find the actual coefficients of the polynomial at the root
 * of the tree.
 */

void
TUpTree (listz_t b, listz_t *Tree, unsigned long long k, listz_t tmp, int dolvl,
         unsigned long long sh, mpz_t n, FILE *TreeFile)
{
    unsigned long long m, l;

    m = k / 2;
    l = k - m;
    
    if (k == 1)
      return;
   
#ifdef DEBUG
    fprintf (ECM_STDOUT, "In TupTree, k = %d.\n", k);

    fprintf (ECM_STDOUT, "b = ");
    print_vect (b, k);
    fprintf (ECM_STDOUT, "\nThe polynomials at that level are: ");
    print_vect (Tree[0] + sh, k);
    fprintf (ECM_STDOUT, "\n");
#endif

    if (dolvl == 0 || dolvl == -1)
      {
        if (TreeFile != NULL)
          {
            list_inp_raw (tmp + k, TreeFile, l);
#ifdef DEBUG_TREEDATA
            printf ("Read from file: ");
            print_vect (tmp + k, l);
#endif
            TMulGen (tmp + l, m - 1, tmp + k, l - 1, b, k - 1, tmp + k + l, n);
            list_inp_raw (tmp + k, TreeFile, m);
#ifdef DEBUG_TREEDATA
            print_vect (tmp + k, m);
            printf ("\n");
#endif
            TMulGen (tmp, l - 1, tmp + k, m - 1, b, k - 1, tmp + k + m, n);
          }
        else
          {
#ifdef DEBUG_TREEDATA
            printf ("Got from Tree: ");
            print_vect (Tree[0] + sh, l);
            print_vect (Tree[0] + sh + l, m);
            printf ("\n");
#endif
            TMulGen (tmp + l, m - 1, Tree[0] + sh, l - 1, b, k - 1, tmp + k, n);
            TMulGen (tmp, l - 1, Tree[0] + sh + l, m - 1, b, k - 1, tmp + k, n);
          }

#if defined(DEBUG) || defined (DEBUG_TREEDATA)
        fprintf (ECM_STDOUT, "And the result at that level (before correction) is:");
        print_vect (tmp, k);
        fprintf (ECM_STDOUT, "\n");
#endif

        /* GMP-ECM specific: leading coefficients in the product tree
        * are implicit ones, so we need some extra work here.
        */

        list_add (tmp, tmp, b + m, l);
        list_add (tmp + l, tmp + l, b + l, m);

        list_mod (b, tmp, k, n); /* reduce both parts simultaneously */

#ifdef DEBUG
        fprintf (ECM_STDOUT, "And the result at this level is:");
        print_vect (b, k);
        fprintf (ECM_STDOUT, "\n");
#endif
      }
    
    if (dolvl > 0 || dolvl == -1)
      {
        if (dolvl > 0)
          dolvl--;
        TUpTree (b,     Tree + 1, l, tmp, dolvl, sh,     n, TreeFile);
        TUpTree (b + l, Tree + 1, m, tmp, dolvl, sh + l, n, TreeFile);
      }
}

static unsigned long long
TUpTree_space (unsigned long long k)
{

    unsigned long long m, l;
    unsigned long long r1, r2;

    m = k / 2;
    l = k - m;
    
    if (k == 1)
      return 0;
   
    r1 = TMulGen_space (l - 1, m - 1, k - 1) + l;
    if (m != l)
      {
        r2 = TMulGen_space (m - 1, l - 1, k - 1) + k;
        r1 = MAX (r1, r2);
      }

    r2 = TUpTree_space (l);
    r1 = MAX (r1, r2);
    
    if (m != l)
      {
        r2 = TUpTree_space (m);
        r1 = MAX (r1, r2);
      }

    return r1;
}

/* This is the documentation of the (now removed) polyeval() function.
   Algorithm polyeval from section 3.7 of Peter Montgomery's dissertation.
Input: 
   G - an array of k elements of R, G[i], 0 <= i < k
       representing the coefficients of a polynomial G(x) of degree < k
   Tree - the product tree produced by PolyFromRoots
   Tree[0][0..k-1] (degree k/2)
   Tree[1][0..k-1] (degree k/4), ...,
   Tree[lgk-1][0..k-1] (degree 1)
Output: the sequence of values of G(a[i]) are stored in G[i] for 0 <= i < k
Remark: we need an auxiliary (k+1)-th cell G[k] in G.
The memory used is M(k) = max(3*floor(k/2)+list_mul_mem(floor(k/2)),
                              k+list_mul_mem(ceil(k/2)),
                              floor(k/2) + M(ceil(k/2))).
Since list_mul_mem(k) >= 2*k, the maximum is the 1st.
*/

/* Same as polyeval. Needs invF as extra argument.
   Return non-zero iff an error occurred.
*/
int
polyeval_tellegen (listz_t b, unsigned long long k, listz_t *Tree, listz_t tmp,
                   unsigned long long sizeT, listz_t invF, mpz_t n, 
                   char *TreeFilename)
{
    unsigned long long tupspace;
    unsigned long long tkspace;
    int allocated = 0, 
        r = 0; /* return value, 0 = no error */
    listz_t T;

    ASSERT(Tree != NULL || TreeFilename != NULL);
    
    tupspace = TUpTree_space (k) + k;
    tkspace = 2 * k - 1 + list_mul_mem (k);

    tupspace = MAX (tupspace, tkspace);
    
    if (TreeFilename != NULL)
      tupspace += (k + 1) / 2;

    if (sizeT >= tupspace)
        T = tmp;
    else
      {
        outputf (OUTPUT_DEVVERBOSE, "polyeval_tellegen: allocating extra temp"
                 " space, want %d but T has only %d\n", tupspace, sizeT);
        T = init_list (tupspace);
	if (T == NULL)
	  return ECM_ERROR;
        allocated = 1;
      }
    
#ifdef TELLEGEN_DEBUG
    fprintf (ECM_STDOUT, "In polyeval_tellegen, k = %d.\n", k);
    fprintf (ECM_STDOUT, "Required memory: %d.\n", 
	     TMulGen_space (k - 1, k - 1, k - 1));
#endif

    if (Fermat)
      {
        /* Schoenhage-Strassen can't do a half product faster than a full */
        F_mul (T, invF, b, k, DEFAULT, Fermat, T + 2 * k);
        list_mod (T, T + k - 1, k, n);
      }
    else
      {
        /* need space 2k-1+list_mul_mem(k) in T */
        list_mul_high (T, invF, b, k);
        list_mod (T, T + k - 1, k, n);
      }
    list_revert (T, k);
    if (TreeFilename != NULL)
      {
        unsigned int lgk, i;
        FILE *TreeFile;
	char *fullname = (char *) malloc (strlen (TreeFilename) + 1 + 2 + 1);
        ASSERT_ALWAYS(fullname != NULL);

	lgk = ceil_log2 (k);
        for (i = 0; i < lgk; i++)
          {
            sprintf (fullname, "%s.%d", TreeFilename, i);
            
	    TreeFile = fopen (fullname, "rb");
            if (TreeFile == NULL)
              {
                outputf (OUTPUT_ERROR, 
                         "Error opening file %s for product tree of F\n",
                         fullname);
                r = ECM_ERROR;
                goto clear_T;
              }
            TUpTree (T, NULL, k, T + k, i, 0, n, TreeFile);
            fclose (TreeFile);
            auto err = _unlink (fullname);
			assert(err == 0);
          }
        free (fullname);
      }
    else
      TUpTree (T, Tree, k, T + k, -1, 0, n, NULL);
    list_swap (b, T, k); /* more efficient than list_set, since T is not
                            needed anymore */

clear_T:
    if (allocated)
      clear_list (T, tupspace);

    return r;
}
