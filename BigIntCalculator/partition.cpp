#include "pch.h"


std::vector <Znum>bistoredp;   /* calculated values of p(n) */
Znum bi1;

/* calculate partition number of n. Assumes that p(n) has already been calculated 
for all smaller values of n*/
static void bipermwork(int n, Znum &result) {
	int gk;

	if (n <= 1) {   //  p(1) =1, by convention p(0) is 1
		result = 1;
		bistoredp[0] = 1;   // set p(0) to 1
		bistoredp[1] = 1;   // set p(0) to 1
		return;
	}

	bi1 = 0;                // set total to 0
	for (int k = 0;; ) {
		/* k = 0, 1, -1, 2, -2, 3, -3 etc... */
		if (k <= 0) k = 1 - k;
		else k = -k;			// generate next generalised pentagonal number.
		gk = (k*(3 * k - 1)) / 2;
		if (gk > n) 
			break;		// if true, the rest of the infinite number of terms are 0
								// because p(n) is 0 if n <0
		if ((k % 2) == 0)
			bi1 -= bistoredp[n - gk];
		else
			bi1 += bistoredp[n - gk];
	}
	bistoredp[n] = bi1;
	result = bi1;
	return;
}

/* Unrestricted Partition Number (number of decompositions of n into sums 
of integers without regard to order) 
see https://oeis.org/A000041
also https://en.wikipedia.org/wiki/Partition_(number_theory) 
also https://en.wikipedia.org/wiki/Partition_function_(number_theory)
also https://en.wikipedia.org/wiki/Pentagonal_number */
void biperm(int n, Znum &result) {
	static int maxn = 0;

	if (n >= bistoredp.size())
		bistoredp.resize(n + 1);
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
		result = bistoredp[n];
	}
}