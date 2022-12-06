/* i5m.c */
/* modular arithmetic */
#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_NONSTDC_NO_DEPRECATE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "integer.h"
#include "fun.h"

extern unsigned int reciprocal[][Z0];
extern long int nettbytes;
void DUMPMPI(MPI *Mptr, char *name)
{
	unsigned int i;
	if (Mptr == NULL) {
		printf("NULL\n");
		return;
	}

	printf("%s:", name);
	printf("Mptr->S = %d, Mptr->D = %u, ", Mptr->S, Mptr->D);
	for (i = 0; i <= Mptr->D; i++)
		printf("Mptr->V[%u] = %lu, ", i, Mptr->V[i]);
	printf("\n");
	return;
}

void PRINT_MPI(MPI *Mptr, char *name)
{
	if (Mptr == NULL) {
		printf("NULL\n");
		return;
	}
	printf("%s = ", name);
	PRINTI(Mptr);
	printf("\n");
	return;
}

MPI *ADDM(MPI *Aptr, MPI *Bptr, MPI *Mptr)
/*
* *Cptr = (*Aptr + *Bptr) mod (*Mptr).
* Here 0 <= *Aptr, *Bptr, *Cptr < *Mptr.
*/
{
	MPI *Cptr, *TempI;

	Cptr = ADD0I(Aptr, Bptr);
	if (RSV(Cptr, Mptr) >= 0)
	{
		TempI = Cptr;
		Cptr = SUB0I(Cptr, Mptr);
		FREEMPI(TempI);
	}
	return (Cptr);
}

unsigned long ADDm(USL a, USL b, USL m)
/*
* returns a + b mod m, if m > 0, where 0 <= a,b < m < 2^32.
* a + b in GF(4) if m = 0.
*/
{
	unsigned long c;

	if (a == 0)
		return b;
	if (b == 0)
		return a;
	if (m)
	{
		c = (a < m - b) ? a + b : a - (m - b);
		return c;
	}
	else
	{
		if (a == 1 && b == 1)
			return (USL)0;
		else if (a == 1 && b == 2)
			return (USL)3;
		else if (a == 1 && b == 3)
			return (USL)2;
		else if (a == 2 && b == 1)
			return (USL)3;
		else if (a == 2 && b == 2)
			return (USL)0;
		else if (a == 2 && b == 3)
			return (USL)1;
		else if (a == 3 && b == 1)
			return (USL)3;
		else if (a == 3 && b == 2)
			return (USL)1;
		else  /*(a == 3 && b == 3)*/
			return (USL)0;
	}
}

MPI *MINUSM(MPI *Aptr, MPI *Mptr)
/*
* *Bptr = -(*Aptr) mod (*Mptr).
* Here 0 <= *Aptr, *Bptr < *Mptr.
*/
{

	if (Aptr->S == 0)
		return ZEROI();
	else
		return SUB0I(Mptr, Aptr);
}

unsigned long MINUSm(USL a, USL m)
/*
* returns -a mod m if m > 0; here 0 <= a < m < 2^32.
* -a  = a in GF(4) if m = 0.
*/
{
	if (m)
	{
		if (a == 0)
			return (USL)0;
		else
			return (m - a);
	}
	else
		return a;
}

MPI *SUBM(MPI *Aptr, MPI *Bptr, MPI *Mptr)
/*
* *Cptr = (*Aptr - *Bptr) mod (*Mptr).
* Here 0 <= *Aptr, *Bptr, *Cptr < *Mptr.
*/
{
	MPI *Cptr, *TempI;

	Cptr = MINUSM(Bptr, Mptr);
	TempI = Cptr;
	Cptr = ADDM(Aptr, Cptr, Mptr);
	FREEMPI(TempI);
	return (Cptr);
}

unsigned long SUBm(USL a, USL b, USL m)
/*
* returns a - b mod m if m > 0; here 0 <= a, b < m < 2^32.
* a - b = a + b in GF(4) if m = 0.
*/
{
	if (m)
	{
		if (a >= b)
			return (a - b);
		else
		{
			return (MINUSm(b - a, m));
		}
	}
	else
		return ADDm(a, b, (USL)0);
}

MPI *MULTM(MPI *Aptr, MPI *Bptr, MPI *Mptr)
/*
* *Cptr = (*Aptr * *Bptr) mod (*Mptr).
* Here 0 <= *Aptr, *Bptr, *Cptr < *Mptr.
*/
{
	MPI *TempI, *Cptr;

	Cptr = MULTI(Aptr, Bptr);
	TempI = Cptr;
	Cptr = MOD0(Cptr, Mptr);
	FREEMPI(TempI);
	return (Cptr);
}

unsigned long MULTm(USL a, USL b, USL m)
/*
* returns a * b mod m  if m > 0; here 0 <= a, b < m < 2^32.
* a * b in GF(4) if m = 0.
*/
{
	if (a == 0 || b == 0)
		return (USL)0;
	if (m)
		return (RUSSm(0, a, b, m));
	else
	{
		if (a == 1)
			return b;
		else if (b == 1)
			return a;
		else if (a == 2 && b == 2)
			return (USL)3;
		else if (a == 2 && b == 3)
			return (USL)1;
		else if (a == 3 && b == 2)
			return (USL)1;
		else  /*(a == 3 && b == 3)*/
			return (USL)2;
	}
}

MPI *INVERSEM(MPI *Nptr, MPI *Mptr)
/*
* here gcd(*Nptr, *Mptr) = 1, 1 <= *Nptr < *Mptr.
* *Kptr satisfies 1 <= *Kptr < *Mptr with (*Kptr) * (*Nptr) = 1  mod *Mptr.
* See Knuth, p. 325.
*/
{
	MPI *A, *B, *C, *H, *L, *Q, *Kptr;
	MPI *TempI;
	unsigned long i;

	if (EQONEI(Nptr))
		return ONEI();
	A = COPYI(Mptr);
	B = COPYI(Nptr);
	C = MOD(A, B);
	Kptr = ONEI();
	L = ZEROI();
	i = 1;
	while (C->S > 0)
	{
		Q = INT0(A, B);
		FREEMPI(A);
		A = B;
		B = C;
		C = MOD0(A, B);
		TempI = Q;
		Q = MULTI(Q, Kptr);
		FREEMPI(TempI);
		H = ADD0I(L, Q);
		FREEMPI(Q);
		FREEMPI(L);
		L = Kptr;
		Kptr = H;
		i = 1 - i;
	}
	FREEMPI(L);
	FREEMPI(A);
	FREEMPI(B);
	FREEMPI(C);
	if (i == 0)
	{
		TempI = Kptr;
		Kptr = SUB0I(Mptr, Kptr);
		FREEMPI(TempI);
	}
	return (Kptr);
}

unsigned long GCDm(USL m, USL n)
/*
* returns the gcd of m and n.
*/
{
	unsigned long c;

	if (n == 0)
		return m;
	c = m % n;
	while (c)
	{
		m = n;
		n = c;
		c = m % n;
	}
	return n;
}

unsigned long INVERSEm(USL n, USL m)
/*
*  If 1 <= n < m < 2^32, gcd(n, m) = 1; returns k, 1 <= k < m with
* n * k congruent to 1 mod m.
* returns n^{-1} in GF(4) if m = 0.
*/
{
	unsigned long a, b, c, h, i, k, l, q;
	unsigned int t;

	if (n == 1)
		return ((USL)1);
	if (m)
	{
		if (n == m - 1)
			return (m - 1);
		if (m <= Z0 + 1)
		{
			t = reciprocal[m - 2][n - 1];
			return (USL)t;
		}
		a = m;
		b = n;
		c = a % b;
		l = 0;
		k = 1;
		i = 1;
		while (c)
		{
			q = a / b;
			a = b;
			b = c;
			c = a % b;
			h = l + q * k; /* no overflow - see Knuth p. 595 */
			l = k;
			k = h;
			i = 1 - i;
		}
		if (i)
			return (k);
		else
			return (m - k);
		return (k);
	}
	else
	{
		if (n == 2)
			return (USL)3;
		else /* (n == 3)*/
			return (USL)2;
	}
}

MPI *DIVM(MPI *Aptr, MPI *Bptr, MPI *Mptr)
/*
* *Cptr = (*Aptr / *Bptr) mod (*Mptr).
* Here 0 <= *Aptr, *Bptr, *Cptr < *Mptr, gcd(*Bptr, *Mptr) = 1.
*/
{
	MPI *K, *Cptr;

	K = INVERSEM(Bptr, Mptr);
	Cptr = MULTM(Aptr, K, Mptr);
	FREEMPI(K);
	return (Cptr);
}

unsigned long DIVm(USL a, USL b, USL m)
/*
* returns a / b mod m if m > 0, a / b in GF(4) if m = 0.
* here 0 <= a, b < m < R0, gcd(b, m) = 1 if m > 0, b != 0 if m = 0.
*/
{
	unsigned long k;

	k = INVERSEm(b, m);
	return (MULTm(a, k, m));
}

unsigned long POWERm(USL a, MPI *Bptr, USL m)
/*
* returns a^ *Bptr mod m.
* here 0 <= a < m < R0 and 0 <= *Bptr.
*/
{
	unsigned long x, z;
	MPI *Y, *TempI;

	x = a;
	Y = COPYI(Bptr);
	z = 1;
	while (Y->S)
	{
		while (!(Y->V[0] & 1))
		{
			TempI = Y;
			Y = INT0_(Y, (USL)2);
			FREEMPI(TempI);
			x = MULTm(x, x, m);
		}
		TempI = Y;
		Y = SUB0_I(Y, (USL)1);
		FREEMPI(TempI);
		z = MULTm(z, x, m);
	}
	FREEMPI(Y);
	return (z);
}

unsigned long POWER_m(USL a, USL y, USL m)
/*
* returns a^ y mod m.
* Here 0 <= a < m < R0 and 0 <= y is an unsigned int.
*/
{
	unsigned long x, z;

	x = a;
	z = 1;
	while (y)
	{
		while (y % 2 == 0)
		{
			y = y / 2;
			x = MULTm(x, x, m);
		}
		y = y - 1;
		z = MULTm(z, x, m);
	}
	return (z);
}

MPI *MPOWER(MPI *Aptr, MPI *Bptr, MPI *Cptr)
/*
* *Aptr, *Bptr, *Cptr, *Eptr  are MPI'S, *Cptr > 0, *Bptr >= 0,
* Eptr = (*Aptr)^ *Bptr (mod *Cptr).
*/
{
	MPI *X, *Y, *Z, *TempI;

	X = MOD(Aptr, Cptr);
	Y = COPYI(Bptr);
	Z = ONEI();
	while (Y->S)
	{
		while (!(Y->V[0] & (USL)1))
		{
			TempI = Y;
			Y = INT0_(Y, (USL)2);
			FREEMPI(TempI);
			TempI = X;
			X = MULTI(X, X);
			FREEMPI(TempI);
			TempI = X;
			X = MOD0(X, Cptr);
			FREEMPI(TempI);
		}
		TempI = Y;
		Y = SUB0_I(Y, (USL)1);
		FREEMPI(TempI);
		TempI = Z;
		Z = MULTM(Z, X, Cptr);
		FREEMPI(TempI);
	}
	FREEMPI(X);
	FREEMPI(Y);
	return (Z);
}

MPI *MPOWER_M(MPI *Aptr, USL b, MPI *Cptr)
/*
* *Aptr, *Cptr, *Eptr  are MPI'S, *Cptr > 0, b an unsigned int.
* *Eptr = (*Aptr)^ b (mod *Cptr).
*/
{
	MPI *X, *Z, *TempI;
	unsigned long y;

	X = MOD(Aptr, Cptr);
	y = b;
	Z = ONEI();
	while (y)
	{
		while (y % 2 == 0)
		{
			y = y / 2;
			TempI = X;
			X = MULTI(X, X);
			FREEMPI(TempI);
			TempI = X;
			X = MOD0(X, Cptr);
			FREEMPI(TempI);
		}
		y = y - 1;
		TempI = Z;
		Z = MULTM(Z, X, Cptr);
		FREEMPI(TempI);
	}
	FREEMPI(X);
	return (Z);
}

MPI *MPOWER_(long a, MPI *Bptr, MPI *Cptr)
/*
* 0 < |a| < R0.
*  *Cptr, *Eptr  are MPI'S, *Cptr > 0, *Bptr >= 0,
* *Eptr = a^ *Bptr (mod *Cptr).
*/
{
	MPI *X, *Y;

	X = CHANGEI(a);
	Y = MPOWER(X, Bptr, Cptr);
	FREEMPI(X);
	return (Y);
}

MPI *BINOMIAL(USI n, USI m)
/*
* returns n choose m, where n >= m are unsigned integers.
*/
{
	unsigned int j;
	MPI *G, *TempI, *Aptr;

	Aptr = ONEI();
	for (j = 1; j <= m; j++)
	{
		G = CHANGE((USL)n - j + 1);
		TempI = Aptr;
		Aptr = MULTI(G, Aptr);
		FREEMPI(G);
		FREEMPI(TempI);
		TempI = Aptr;
		Aptr = INT0_(Aptr, (USL)j);
		FREEMPI(TempI);
	}
	return (Aptr);
}

unsigned int LENGTHm(USL n)
/*
* returns the number of decimal digits in the unsigned int n.
*/
{
	unsigned int i = 0;

	do
	{
		n = n / 10;
		i++;
	} while (n);
	return (i);
}

unsigned long RUSSm(USL a, USL b, USL c, USL p)
/*
*  input: unsigned long ints a, b, c, p, with p > 0.
* output: a + b * c (mod p).
* Russian Peasant multiplication algorithm. Uses the identities
* RUSSm(a, b, 2 * c, p) = RUSSm(a, 2 * b, c, p),
* RUSSm(a, b, c + 1, p) = RUSSm(a + b, b, c, p).
* If a, b, c and p are less than 2^32, * so is RUSSm(a, b, c, p).
* From H. Luneburg, "On the Rational Normal Form of * Endomorphisms",
* B.I. WissenSchaftsverlag, Mannheim/Wien/Zurich, 1987, pp 18-19.
* Luneburg's restriction to 2*p<2^32 removed by krm on 18/4/94.
*/
{
	USL t;

	a = a % p;
	b = b % p;
	c = c % p;
	while (c)
	{
		while (!(c & 1))
		{
			c = c >> 1;
			t = p - b;
			b = (b < t) ? (b << 1) % p : (b - t) % p;
		}
		c--;
		t = p - b;
		a = (a < t) ? a + b : a - t;
	}
	return a;
}

unsigned long RANDOMm(USL x)
/*
* input: unsigned int x, output:a "random number" a * x + c (mod m).
* a = 1001, m = R0 = 65536, c = 65;
* From H. Luneburg, "On the Rational Normal Form of Endomorphisms",
* B.I. WissenSchaftsverlag, Mannheim/Wien/Zurich, 1987.
* See Knuth Vol 2, Theorem A, p. 16.
*/
{
	unsigned long a, c, m;

	m = R0;
	a = 1001;
	c = 65;
	return RUSSm(c, a, x, m);
}

unsigned long *FFT(USL N, USL *a, USL w, USL p)
/*
* Fast Fourier Transform.
* Here N is a power of 2, a is an array of N elements from Z_p and
* w is a primitive N-th root of unity mod p and N | p - 1.
* See "Elements of Algebra and Algebraic Computing", J.D. Lipson, p.298.
* Outputs a(w^k), 0 <= k < N, where a(x)=\sum_{k=0}^{N-1} a[k]x^k.
*/
{
	USL *A, *B, *C, *b, *c, i, k, n, s, t, u;
	A = (USL *)mmalloc(N * sizeof(USL));
	if (N == 1)
	{
		A[0] = a[0];
		return (A);
	}
	else
	{
		n = N >> 1;
		/*printf("level N=%lu\n", N);*/
		b = (USL *)mmalloc(n * sizeof(USL));
		for (i = 0; i < n; i++)
			b[i] = a[2 * i];
		/*printf("creating c of length %lu\n", n);*/
		c = (USL *)mmalloc(n * sizeof(USL));
		for (i = 0; i < n; i++)
			c[i] = a[2 * i + 1];
		t = RUSSm(0, w, w, p);
		B = FFT(n, b, t, p);
		C = FFT(n, c, t, p);
		for (k = 0; k < n; k++)
		{
			u = POWER_m(w, k, p);
			s = RUSSm(0, u, C[k], p);
			A[k] = ADDm(B[k], s, p);
			A[k + n] = SUBm(B[k], s, p);
		}
		ffree((USL *)b, n * sizeof(USL));
		ffree((USL *)c, n * sizeof(USL));
		return (A);
	}
}

unsigned long *FFI(USL N, USL *b, USL w, USL p)
/*
* Here N is a power of 2, b is an array of N elements from Z_p and
* w is a primitive N-th root of unity mod p and N | p - 1.
* Fast Fourier Interpolation.
* Outputs a(x)=\sum_{k=0}^{N-1} a[k]x^k, where a(w^k)=b[k], 0 <= k < N.
* See "Elements of Algebra and Algebraic Computing", J.D. Lipson, p.303.
*/
{
	USL *c, i, t;

	c = FFT(N, b, INVERSEm(w, p), p);
	t = INVERSEm(N, p);
	for (i = 0; i < N; i++)
		c[i] = RUSSm(0, t, c[i], p);
	return (c);
}

unsigned long *FFP(USL N, USL *a, USL *b, USL m, USL n, USL w, USL p)
/*
* Input: arrays a and b of USL's mod p, representing polynomials of
* degrees m, n, respectively; N = 2^e > m + n, N | p - 1.
* w is a primitive N-th root of unity mod p.
* Output: array c mod p, representing a(x)b(x).
*/
{
	USL *aa, *bb, *c, *A, *B, *C, i;

	/* pad a and b with zeros. */
	aa = (USL *)ccalloc(N, sizeof(USL));
	bb = (USL *)ccalloc(N, sizeof(USL));
	for (i = 0; i <= m; i++)
		aa[i] = a[i];
	for (i = 0; i <= n; i++)
		bb[i] = b[i];
	A = FFT(N, aa, w, p);
	B = FFT(N, bb, w, p);
	ffree((USL *)aa, N * sizeof(USL));
	ffree((USL *)bb, N * sizeof(USL));
	C = (USL *)mmalloc(N * sizeof(USL));
	for (i = 0; i < N; i++)
		C[i] = RUSSm(0, A[i], B[i], p);
	ffree((USL *)A, N * sizeof(USL));
	ffree((USL *)B, N * sizeof(USL));
	c = FFI(N, C, w, p);
	ffree((USL *)C, N * sizeof(USL));
	c = (USL *)rrealloc(c, (m + n + 1) * sizeof(USL), -((long)(N - (m + n + 1)) * sizeof(USL)));
	return (c);
}

MPI *CRA(USL n, USL *a, USL *m)
/*
* Garner's Chinese Remainder Theorem. Page 180, "Algorithms for computer
* algebra", K.O. Geddes, S.R. Czapor, G. Labahn".
* Here gcd(m[i],m[j])=1 if i != j.
* Returns the least remainder mod(m[0].....m[n]) of u=a[i]mod(m[i]),0<=i<=n.
*/
{
	USL *g, *nu, i, temp, product;
	long int k, j;
	MPI *U, *UU, *TEMP;

	/* Compute g[k]=(m[0]...m[k-1])^{-1}mod(m[k]),1<=k<=n. */
	g = (USL *)mmalloc(n * sizeof(USL));
	nu = (USL *)mmalloc((n + 1) * sizeof(USL));
	for (k = 1; k <= n; k++)
	{
		product = m[0] % m[k];
		for (i = 1; i < k; i++)
			product = RUSSm(0, product, m[i], m[k]);
		g[k] = INVERSEm(product, m[k]);
	}
	/* compute the mixed radix coefficients nu[k], 0<=k<=n. */
	nu[0] = a[0];
	for (k = 1; k <= n; k++)
	{
		temp = nu[k - 1];
		for (j = k - 2; j >= 0; j--)
			temp = RUSSm(nu[j], temp, m[j], m[k]);
		nu[k] = RUSSm(0, SUBm(a[k], temp, m[k]), g[k], m[k]);
	}

	/* convert to base R0 */
	U = CHANGE(nu[n]);
	for (k = n - 1; k >= 0; k--)
	{
		TEMP = U;
		U = MULT_II(U, m[k]);
		FREEMPI(TEMP);
		TEMP = U;
		UU = CHANGE(nu[k]);
		U = ADD0I(U, UU);
		FREEMPI(UU);
		FREEMPI(TEMP);
	}
	ffree((USL *)g, n * sizeof(USL));
	ffree((USL *)nu, (n + 1) * sizeof(USL));
	return (U);
}

USL fp[4] = { 2013265921UL, 2281701377UL, 3221225473UL, 3489660929UL };
USL lprimroot[4] = { 31, 3, 5, 3 };

MPI *FFM(MPI *a, MPI *b, USL K)
/*
* Returns the product of a=(a[0],...,a[m])_B and b=(b[0],...,b[n])_B, B=2^16,
* m = a->D, n = b->D, using the Discrete Fast Fourier Transform.
* Let M = min(m,n). Then a(x)*b(x)=c(x), where 0<= c[k] < (M+1)B^2.
* Using the CRA mod fp[i] for 0 <= i <= K-1, enables us to reconstruct c(B),
* provided that fp[0]...fp[K-1] >= (M+1)B^2. We also need m+n < N = 2^e,
* where N | fp[i] - 1.
* If m<B and n<B, then M<B, then K=3 primes suffice as fp[i]>=B; take e>=17.
* If m<2^26 and n<2^26, then M<2^26, then K=4 primes suffice; take e>=27.
*      N = 2^27
* fp[0] = 2013265921, lprimroot = 31, primitive Nth root =  440564289;
* fp[1] = 2281701377, lprimroot =  3, primitive Nth root =  129140163;
* fp[2] = 3221225473, lprimroot =  5, primitive Nth root =  229807484;
* fp[3] = 3489660929, lprimroot =  3, primitive Nth root = 1392672017.
* See "Elements of Algebra and Algebraic Computing", J.D. Lipson, p.310.
*/
{
	USL N, i, t, **A, *w, *temp, r;
	long int j;
	MPI **c, *U, *TEMP;
	unsigned int m, n;

	m = a->D;
	n = b->D;
	if (m >= 67108864 || n >= 67108864)
	{
		fprintf(stderr, "overflow in FFM\n");
		exit(1);
	}
	N = 1;
	t = m + n;
	while (N <= t)
		N = 2 * N;
	A = (USL **)mmalloc(K * sizeof(USL *));
	w = (USL *)mmalloc(K * sizeof(USL));
	for (i = 0; i < K; i++)
	{
		w[i] = POWER_m(lprimroot[i], (fp[i] - 1) / N, fp[i]);
		A[i] = FFP(N, a->V, b->V, m, n, w[i], fp[i]);
	}
	/* Each A[i] has length m + n + 1 and represents a(x)*b(x) mod fp[i]. */
	c = (MPI **)mmalloc((t + 1) * sizeof(MPI *));
	temp = (USL *)mmalloc(K * sizeof(USL));
	r = K - 1;
	for (j = 0; j <= t; j++)
	{
		for (i = 0; i < K; i++)
			temp[i] = A[i][j];
		c[j] = CRA(r, temp, fp);
	}
	U = COPYI(c[t]);
	/* bug site: c[i] may be >= R0 */
	for (j = (long)t - 1; j >= 0; j--)
	{
		TEMP = U;
		U = MULT_II(U, (long)65536);
		FREEMPI(TEMP);
		TEMP = U;
		U = ADD0I(U, c[j]);
		FREEMPI(TEMP);
	}
	for (j = 0; j <= t; j++)
		FREEMPI(c[j]);
	ffree((MPI **)c, (t + 1) * sizeof(MPI *));
	ffree((USL *)w, K * sizeof(USL));
	for (i = 0; i < K; i++)
		ffree((USL *)A[i], (t + 1) * sizeof(USL));
	ffree((USL **)A, K * sizeof(USL *));
	ffree((USL *)temp, K * sizeof(USL));
	return (U);
}

MPI *RANDOMI(MPI *X, MPI *P, MPI *T) {
	/* If Z=MOD(P*X,T), returns Z if Z < T/2,  Z-T otherwise. */
	/* Usually take P=96069 or 1+4n, n odd and T a power of 2. */
	MPI *Z, *W, *U, *V; int t;

	V = MOD(X, T);
	Z = MULTM(P, V, T);
	FREEMPI(V);
	W = MULT_II(Z, 2);
	t = RSV(W, T);
	FREEMPI(W);
	if (t == 1) {
		U = SUBI(Z, T);
	}
	else
		U = COPYI(Z);
	FREEMPI(Z);
	return (U);
}


void RANDOM_MATRIX(MPI *M, MPI *N, MPI *X, MPI *P, MPI *T) {
	/* returns a random m x n matrix of integers, |A[i][j]| < T/2. */
	USI i, j, m, n;
	MPMATI *A;
	MPI *Temp;

	m = (USI)CONVERTI(M);
	n = (USI)CONVERTI(N);
	A = BUILDMATI(m, n);
	Temp = X;
	for (i = 0; i < m; i++)
		for (j = 0; j < n; j++)
		{
			Temp = RANDOMI(Temp, P, T);
			A->V[i][j] = Temp;
		}
	PRINTMATI(0, A->R - 1, 0, A->C - 1, A);
	FREEMATI(A);
	return;
}

MPMATI *JCOLSMATI(MPMATI *Mptr, MPMATI *Nptr)
/*
* returns [*Mptr | *Nptr]
*/
{
	unsigned int i, j;
	MPMATI *Lptr;

	Lptr = BUILDMATI(Mptr->R, Mptr->C + Nptr->C);
	for (i = 0; i < Lptr->R; i++)
	{
		for (j = 0; j < Lptr->C; j++)
		{
			if (j < Mptr->C)
				Lptr->V[i][j] = COPYI(Mptr->V[i][j]);
			else
				Lptr->V[i][j] = COPYI(Nptr->V[i][j - Mptr->C]);
		}
	}
	return (Lptr);
}

MPMATI *RANDOM_MATRIXA(USI m, USI n, MPI *X, MPI *P, MPI *T)
{
	/* returns a random m x (n+1) augmented matrix of integers, |A[i][j]| < T/2. */
	/* which is consistent. */
	USI i, j;
	MPMATI *A, *B, *C, *D;
	MPI *Temp;

	A = BUILDMATI(m, n);
	Temp = X;
	for (i = 0; i < m; i++)
		for (j = 0; j < n; j++)
		{
			Temp = RANDOMI(Temp, P, T);
			A->V[i][j] = Temp;
		}
	C = BUILDMATI(n, 1);
	for (i = 0; i < n; i++)
	{
		Temp = RANDOMI(Temp, P, T);
		C->V[i][0] = Temp;
	}
	B = MULTMATI(A, C);
	FREEMATI(C);
	D = JCOLSMATI(A, B);
	FREEMATI(A);
	FREEMATI(B);
	return (D);

}

MPMATI *RANDOM_MATRIXA3(USI m, USI n, MPI *X, MPI *P, MPI *T)
{
	/*
	* returns a random m x (n+1) augmented matrix of integers,
	* |A[i][j]| < T/2 which is consistent, with X components
	* 0,-1,1.
	*/
	USI i, j;
	MPMATI *A, *B, *C, *D;
	MPI *Temp, *Temp1, *Temp2;

	Temp1 = THREEI();
	A = BUILDMATI(m, n);
	Temp = X;
	for (i = 0; i < m; i++)
		for (j = 0; j < n; j++)
		{
			Temp = RANDOMI(Temp, P, T);
			A->V[i][j] = Temp;
		}
	C = BUILDMATI(n, 1);
	Temp2 = COPYI(Temp);
	for (i = 0; i < n; i++)
	{
		Temp = Temp2;
		Temp2 = RANDOMI(Temp2, P, T);
		C->V[i][0] = HALFMOD(Temp2, Temp1);
		FREEMPI(Temp);
	}
	B = MULTMATI(A, C);
	FREEMATI(C);
	D = JCOLSMATI(A, B);
	FREEMATI(A);
	FREEMATI(B);
	FREEMPI(Temp1);
	return (D);
}
