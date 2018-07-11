#include <stdlib.h>
#include <assert.h>

#ifdef __GNUC__
#include "gmp.h"
#else
#include "mpir.h"
#endif

mpz_t bistoredp[60000];
mpz_t bi1;

static void bipermwork(int n, mpz_t &result) {
	int gk;

	if (n <= 1) {   //  p(1) =1, by convention p(0) is 1
		mpz_set_ui(result, 1);
		mpz_init(bistoredp[0]);
		mpz_set_ui(bistoredp[0], 1);  // set p(0) to 1
		mpz_init(bistoredp[1]);
		mpz_set_ui(bistoredp[1], 1);  // set p(1) to 1
		return;
	}
	mpz_inits(bi1, NULL);
	mpz_set_ui(bi1, 0);   // set total to 0;
	for (int k = 0;; ) {
		if (k <= 0) k = 1 - k;
		else k = -k;			// generate next generalised pentagonal number.
		gk = (k*(3 * k - 1)) / 2;
		if (gk > n) 
			break;		// if true, the rest of the infinite number of terms are 0
								// because p(n) is 0 if n <0
		if ((k % 2) == 0)
			mpz_sub(bi1, bi1, bistoredp[n - gk]);
		else
			mpz_add(bi1, bi1, bistoredp[n - gk]);
	}
	mpz_init(bistoredp[n]);
	mpz_set(bistoredp[n], bi1);
	mpz_set(result, bi1);
	return;
}

/* Unrestricted Partition Number (number of decompositions of n into sums 
of integers without regard to order) */
void biperm(int n, mpz_t &result) {
	static int maxn = 0;

	assert(n < sizeof(bistoredp) / sizeof(bistoredp[0]));
	/* make sure we don't exceed maximum array index*/
	/* because the calculation of p(n) uses smaller values of p(n) it is necessary 
	to calculate all of p(1) to p(n) in sequence. Values are stored so they are only 
	calculated once.*/
	if (n > maxn) {
		for (int i = maxn+1; i <= n; i++) {
			bipermwork(i, result);
		}
		maxn = n;
		return;
	}
	else {
		mpz_set(result, bistoredp[n]);
	}
}