/* qres.c */
#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_NONSTDC_NO_DEPRECATE
#include <stdio.h>
#include <stdlib.h>
#include "integer.h"
#include "fun.h"

int JACOBIB(MPI *N, MPI *M)
/*
* Evaluates the Jacobi symbol (N/M); M odd, M > 0.
*/
{
	int sign = 1;
	unsigned long parity, t;
	MPI *Tmp, *NN, *MM;

	NN = MOD(N, M);
	if (NN->S == 0)
	{
		FREEMPI(NN);
		return (0);
	}
	MM = COPYI(M);
	while (1)
	{
		parity = 0;
		while ((NN->V[0] & 1) == 0)
		{
			Tmp = NN;
			NN = INT0_(NN, (USL)2);
			FREEMPI(Tmp);
			parity = 1 - parity;
		}
		if (parity)
		{
			t = MM->V[0] & 7;
			if (t == 3 || t == 5)
				sign = -sign;
		}
		if (EQONEI(NN))
			break;
		if ((MM->V[0] & 3) == 3 && (NN->V[0] & 3) == 3)
			sign = -sign;
		Tmp = NN;
		NN = MOD0(MM, NN);
		FREEMPI(MM);
		if (NN->S == 0)
		{
			FREEMPI(Tmp);
			FREEMPI(NN);
			return (0);
		}
		MM = Tmp;
	}
	FREEMPI(NN);
	FREEMPI(MM);
	return (sign);
}

int JACOBI(USL n, USL m)
/*
* Evaluates the Jacobi symbol (n/m); m odd, 0 < n < m.
*/
{
	int sign = 1;
	unsigned int parity, temp;
	unsigned long t;

	while (1)
	{
		parity = 0;
		while (n % 2 == 0)
		{
			n = n / 2;
			parity = 1 - parity;
		}
		if (parity)
		{
			t = m % 8;
			if (t == 3 || t == 5)
				sign = -sign;
		}
		if (n == 1)
			break;
		if (m % 4 == 3 && n % 4 == 3)
			sign = -sign;
		temp = n;
		n = m %	n;
		if (n == 0)
			return (0);
		m = temp;
	}
	return (sign);
}

void MULTXm(USL a, USL b, USL c, USL d, USL *eptr, USL *fptr, USL x, USL p)
/*
* *eptr+(*fptr)*sqrt(x)=(a+b*sqrt(x))(c+dsqrt(x)) (mod p). Used in R.C. Peralta's
* probabilistic method of finding sqrt(x) (mod p).
*/
{
	unsigned long u, v;

	u = (b * d) % p;
	v = (a * c) % p;
	u = (u * x) % p;
	*eptr = ADDm(v, u, p);
	u = (b * c) % p;
	v = (a * d) % p;
	*fptr = ADDm(v, u, p);
	return;
}

void POWERXm(USL a, USL b, USL n, USL *eptr, USL *fptr, USL x, USL p)
/*
* *eptr+(*fptr)*sqrt(x)=(a+b*sqrt(x))^n (mod p). Used in R.C. Peralta's
* probabilistic method of finding sqrt(x) (mod p).
*/
{
	unsigned long x1, x2, z1, z2;

	x1 = a;
	x2 = b;
	z1 = 1;
	z2 = 0;
	while (n)
	{
		while (n % 2 == 0)
		{
			n = n / 2;
			MULTXm(x1, x2, x1, x2, &x1, &x2, x, p);
		}
		n = n - 1;
		MULTXm(z1, z2, x1, x2, &z1, &z2, x, p);
	}
	*eptr = z1;
	*fptr = z2;
	return;
}

unsigned long SQRTm(USL x, USL p)
/*
* Calculates sqrt(x) (mod p) using "A simple and fast probabilistic algorithm
* for computing square roots modulo a prime number", I.E.E.E. Trans. Inform.
* Theory, IT-32, 1986, 846-847, R. Peralta.
* Here x is a quadratic residue mod p, 0 < x < p.
*/
{
	unsigned long r, z, u, v;

	z = (p - 1) / 2;
	if (p % 4 == 3)
		return (POWER_m(x, (z + 1) / 2, p));
	for (r = 1; r <= z; r++)
	{
		if (x == (r * r) % p)
			return (r);
		POWERXm(r, (USL)1, z, &u, &v, x, p);
		if (u == 0)
			break;
	}
	return (INVERSEm(v, p));
}

void MULTXM(MPI *a, MPI *b, MPI *c, MPI *d, MPI **eptr, MPI **fptr, MPI *x, MPI *p)
/*
* *eptr+(*fptr)*sqrt(x)=(a+b*sqrt(x))(c+dsqrt(x)) (mod p). Used in R.C. Peralta's
* probabilistic method of finding sqrt(x) (mod p).
*/
{
	MPI *u, *v, *tmp1, *tmp2, *tmp3;

	tmp1 = MULTI(b, d);
	tmp2 = MOD(tmp1, p);
	FREEMPI(tmp1);
	tmp3 = MULTI(tmp2, x);
	FREEMPI(tmp2);
	u = MOD(tmp3, p);
	FREEMPI(tmp3);
	tmp1 = MULTI(a, c);
	v = MOD(tmp1, p);
	FREEMPI(tmp1);
	*eptr = ADDM(v, u, p);
	tmp2 = u;
	tmp1 = MULTI(b, c);
	u = MOD(tmp1, p);
	FREEMPI(tmp1);
	FREEMPI(tmp2);
	tmp2 = v;
	tmp1 = MULTI(a, d);
	v = MOD(tmp1, p);
	FREEMPI(tmp1);
	FREEMPI(tmp2);
	*fptr = ADDM(v, u, p);
	FREEMPI(u);
	FREEMPI(v);
	return;
}

void POWERXM(MPI *a, MPI *b, MPI *n, MPI **eptr, MPI **fptr, MPI *x, MPI *p)
/*
* *eptr+(*fptr)*sqrt(x)=(a+b*sqrt(x))^n (mod p). Used in R.C. Peralta's
* probabilistic method of finding sqrt(x) (mod p).
*/
{
	MPI *x1, *x2, *x3, *x4, *x5, *x6, *z1, *z2, *tmp, *tmp1, *tmp2, *N;

	x1 = COPYI(a);
	x2 = COPYI(b);
	z1 = ONEI();
	z2 = ZEROI();
	N = COPYI(n);
	while (N->S)
	{
		while ((N->V[0]) % 2 == 0)
		{
			tmp = N;
			N = INT0_(N, (USL)2);
			FREEMPI(tmp);
			tmp1 = x1;
			tmp2 = x2;
			MULTXM(x1, x2, x1, x2, &x3, &x4, x, p);
			x1 = x3;
			x2 = x4;
			FREEMPI(tmp1);
			FREEMPI(tmp2);
		}
		tmp = N;
		N = SUB0_I(N, (USL)1);
		FREEMPI(tmp);
		tmp1 = z1;
		tmp2 = z2;
		MULTXM(z1, z2, x1, x2, &x5, &x6, x, p);
		z1 = x5;
		z2 = x6;
		FREEMPI(tmp1);
		FREEMPI(tmp2);
	}
	*eptr = z1;
	*fptr = z2;
	FREEMPI(x1);
	FREEMPI(x2);
	FREEMPI(N);
	return;
}

MPI *SQRTM(MPI *x, MPI *p)
/*
* Calculates sqrt(x) (mod p) using "A simple and fast probabilistic algorithm
* for computing square roots modulo a prime number", I.E.E.E. Trans. Inform.
* Theory, IT-32, 1986, 846-847, R. Peralta.
* Here x is a quadratic residue mod p. x can be negative.
*/
{
	MPI *z, *u, *v, *tmp, *tmp1, *tmp2, *tmp3;
	unsigned int t;
	unsigned long r;

	z = INT0_(p, (USL)2);
	if ((p->V[0]) % 4 == 3)
	{
		tmp = ADD0_I(z, (USL)1);
		FREEMPI(z);
		tmp1 = INT0_(tmp, (USL)2);
		FREEMPI(tmp);
		tmp2 = MPOWER(x, tmp1, p);
		FREEMPI(tmp1);
		return (tmp2);
	}
	tmp = ONEI();
	for (r = 1; r < R0; r++)
	{
		tmp1 = CHANGE(r);
		tmp2 = MULTI(tmp1, tmp1);
		tmp3 = MOD0(tmp2, p);
		FREEMPI(tmp2);
		t = EQUALI(x, tmp3);
		FREEMPI(tmp3);
		if (t)
		{
			FREEMPI(tmp);
			FREEMPI(z);
			/*		FREEMPI(x); deleted on 10th Dec 2000 */
			return (tmp1);
		}
		POWERXM(tmp1, tmp, z, &u, &v, x, p);
		FREEMPI(tmp1);
		if (u->S == 0)
		{
			tmp2 = INVERSEM(v, p);
			FREEMPI(u);
			FREEMPI(v);
			FREEMPI(tmp);
			FREEMPI(z);
			break;
		}
		FREEMPI(u);
		FREEMPI(v);
	}
	return (tmp2);
}

MPI *LUCASU(MPI *N, MPI *Q, MPI *M)
/*
* returns U_n (mod m), where U_{k+1}=U_k-qU_{k-1}, U_0=0,U_1=1.
* D=1-4q != 0.
* We use the recurrence relations:
* U_{k+1}=U_k-qU_{k-1},
* U_{2k}=-2qU_{k-1}U_k+U_kU_k,
* U_{2k-1}=U_kU_k-qU_{k-1}U_{k-1}
* See Niven, Zuckermann and Montgomery, p.202.
* Also D. Bressoud, Factorization and Primality Testing, p.179-191 and
* C. Pomerance, Kent State Lecture Notes, 19-22.
*/
{
	unsigned int *b, h;
	int i;
	MPI *Y, *R, *S, *T, *U, *V, *Temp0, *Temp, *Temp1, *Temp2;

	if (N->S == 0)
		return (ZEROI());
	if (EQONEI(N))
		return(ONEI());
	h = 0;
	b = (USI *)mmalloc(10000 * sizeof(USI));
	Y = COPYI(N);
	while (Y->S)
	{
		b[h] = (Y->V[0]) % 2;
		Temp = Y;
		Y = INT0_(Y, 2);
		FREEMPI(Temp);
		h++;
	}
	R = ZEROI();
	S = ONEI();
	FREEMPI(Y);
	Temp = MINUSI(Q);
	T = MOD(Temp, M);
	FREEMPI(Temp);
	V = ADDM(T, T, M);
	i = h - 1;
	while (1)
	{
		U = S;
		Temp1 = MULTM(V, R, M);
		Temp2 = ADDM(Temp1, S, M);
		S = MULTM(Temp2, S, M);
		FREEMPI(Temp1);
		FREEMPI(Temp2);
		Temp0 = MULTM(R, R, M);
		Temp1 = MULTM(Temp0, T, M);
		FREEMPI(Temp0);
		Temp2 = MULTM(U, U, M);
		FREEMPI(U);
		Temp = R;
		R = ADDM(Temp2, Temp1, M);
		FREEMPI(Temp);
		FREEMPI(Temp1);
		FREEMPI(Temp2);
		i--;
		if (i == -1)
			break;
		if (b[i])
		{
			Temp1 = MULTM(T, R, M);
			FREEMPI(R);
			R = S;
			S = ADDM(S, Temp1, M);
			FREEMPI(Temp1);
		}
	}
	FREEMPI(R);
	FREEMPI(T);
	FREEMPI(V);
	ffree(b, 10000 * sizeof(USI));
	return (S);
}

MPI *LUCASB(MPI *N)
/*
* Let N > 1 be odd. Then LUCASB(N) determines an integer D(!=-3) of the form
* 4m+1,  such that the Jacobi symbol j(D,N)= -1. Then with Q=(1-d)/4, the Lucas
* pseudoprime test examines L=U_{(N+1)/2}mod N. If L=0, N is a Lucas probable
* prime, otherwise N is composite. See "The Pseudoprimes to 25.10^9",
* Mathematics of computation, 35 (1980) 1003-1026. At the end of this paper
* it is conjectured that if N is a strong base 2 pseudoprime and a Lucas
* probable prime, then N is in fact a prime. A $30 prize is offered for a
* counterexample.
* if LUCASB(N) = 1, then N is a Lucas probable prime, else N is composite.
*/
{
	unsigned int h;
	MPI *S, *T, *X, *Q, *L, *Temp, *Temp1;
	int d, t, s;

	X = BIG_MTHROOT(N, 2);
	Temp = MULTI(X, X);
	FREEMPI(X);
	if (EQUALI(Temp, N))
	{
		/* N is a perfect square */;
		FREEMPI(Temp);
		return (ZEROI());
	}
	FREEMPI(Temp);
	s = 5;
	t = -7;
	for (h = 0; 1; h++)
	{
		S = CHANGEI(s);
		if (JACOBIB(S, N) == -1)
		{
			FREEMPI(S);
			d = s;
			break;
		}
		FREEMPI(S);
		s = s + 4;
		T = CHANGEI(t);
		if (JACOBIB(T, N) == -1)
		{
			FREEMPI(T);
			d = t;
			break;
		}
		FREEMPI(T);
		t = t - 4;
	}
	Q = CHANGEI((1 - d) / 4);
	Temp = ADD0_I(N, 1);
	Temp1 = INT0_(Temp, 2);
	FREEMPI(Temp);
	L = LUCASU(Temp1, Q, N);
	FREEMPI(Q);
	FREEMPI(Temp1);
	t = L->S;
	FREEMPI(L);
	if (t == 0)/* N is a Lucas probable prime */
		return (ONEI());
	else /* N is composite */
		return (ZEROI());
}

MPI *LUCAS(MPI *N)
/* Here N is odd and > 1.
* If LUCAS(N) returns 1, then N is a strong base 2 pseudoprime and a Lucas
* probable prime; if LUCAS(N) returns 0, then N is composite.
* See "The Pseudoprimes to 25.10^9", Mathematics of computation, 35 (1980)
* 1003-1026. At the end of this paper it is conjectured that if N is a strong
* base 2 pseudoprime and a Lucas probable prime, then N is in fact a prime.
* A $30 prize is offered for a counterexample.
*/
{
	unsigned int l;
	int t;
	MPI *L;

	if (N->D == 0 && N->S == 1 && N->V[0] == 2) {
		return (ONEI());
	}
	l = MILLER(N, 2);
	if (l == 0)
		return (ZEROI());
	L = LUCASB(N);
	t = L->S;
	FREEMPI(L);
	if (t)
		return (ONEI());
	else
		return (ZEROI());
}

MPI *LEASTQNR(MPI *P)
/*
* Returns NP, the least quadratic non-residue (mod P) if NP < R0.
* Else returns 0.
*/
{
	MPI *Q, *S, *T;
	unsigned int flag = 0;
	unsigned long i;

	Q = INT0_(P, (USL)2);
	T = MPOWER_((long)2, Q, P);
	S = SUBI(T, P);
	if (EQMINUSONEI(S))
	{
		FREEMPI(Q);
		FREEMPI(S);
		FREEMPI(T);
		return (CHANGE((USL)2));
	}
	FREEMPI(S);
	FREEMPI(T);
	for (i = 3; i < R0; i = i + 2)
	{
		T = MPOWER_((long)i, Q, P);
		S = SUBI(T, P);
		if (EQMINUSONEI(S))
		{
			FREEMPI(S);
			FREEMPI(T);
			flag = 1;
			break;
		}
		FREEMPI(S);
		FREEMPI(T);
	}
	FREEMPI(Q);
	if (flag)
		return (CHANGE(i));
	else
		return(ZEROI());
}
