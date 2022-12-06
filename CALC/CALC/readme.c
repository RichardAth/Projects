/* program "readme.c" */
/*
*  List and description of CALC functions.
*/

#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_NONSTDC_NO_DEPRECATE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "integer.h"
#include "fun.h"

void readme()
{
	char *p, name[20];
	char quit_str[] = "quit";
	char absmod_str[] = "absmod";
	char euclid_str[] = "euclid";
	char gcd_str[] = "gcd";
	char gcdv_str[] = "gcdv";
	char gcda_str[] = "gcda";
	char gcdav_str[] = "gcdav";
	char egcd_str[] = "egcd";
	char sgcd_str[] = "sgcd";
	char lllgcd_str[] = "lllgcd";
	char lllgcd0_str[] = "lllgcd0";
	char lllhermite_str[] = "lllhermite";
	char axb_str[] = "axb";
	char fp_str[] = "fp";
	char inhomfp_str[] = "inhomfp";
	char slv_str[] = "slv";
	char lcm_str[] = "lcm";
	char lcma_str[] = "lcma";
	char length_str[] = "length";
	char pollard_str[] = "pollard";
	char nprime_str[] = "nprime";
	char nprimeap_str[] = "nprimeap";
	char jacobi_str[] = "jacobi";
	char peralta_str[] = "peralta";
	char cong_str[] = "congr";
	char chinese_str[] = "chinese";
	char chinesea_str[] = "chinesea";
	char mthroot_str[] = "mthroot";
	char mthrootr_str[] = "mthrootr";
	char fund_str[] = "fund";
	char pell_str[] = "pell";
	char surd_str[] = "surd";
	char mpower_str[] = "mpower";
	char inv_str[] = "inv";
	char elliptic_str[] = "elliptic";
	char factor_str[] = "factor";
	char tau_str[] = "tau";
	char sigma_str[] = "sigma";
	char mobius_str[] = "mobius";
	char euler_str[] = "euler";
	char lprimroot_str[] = "lprimroot";
	char orderm_str[] = "orderm";
	char lucas_str[] = "lucas";
	char serret_str[] = "serret";
	char collatz_str[] = "collatz";
	char cycle_str[] = "cycle";
	char miller_str[] = "miller";
	char hermite_str[] = "hermite";
	char improvep_str[] = "improvep";
	char mlll_str[] = "mlll";
	char smith_str[] = "smith";
	char rsae_str[] = "rsae";
	char encode_str[] = "encode";
	char decode_str[] = "decode";
	char addcubicr_str[] = "addcubicr";
	char powercubicr_str[] = "powercubicr";
	char ordercubicr_str[] = "ordercubicr";
	char addcubicm_str[] = "addcubicm";
	char powercubicm_str[] = "powercubicm";
	char ordercubicm_str[] = "ordercubicm";
	char convergents_str[] = "convergents";
	char lagrange_str[] = "lagrange";
	char perfectpower_str[] = "perfectpower";
	char leastqnr_str[] = "leastqnr";
	char examples_str[] = "examples";
	char content_str[] = "content";
	char primitive_str[] = "primitive";
	char sturm_str[] = "sturm";
	char rootexp_str[] = "rootexp";
	char log_str[] = "log";
	char log1_str[] = "log1";
	char log2_str[] = "log2";
	char sqroot_str[] = "sqroot";
	char cornacchia_str[] = "cornacchia";
	char patz_str[] = "patz";
	char congq_str[] = "congq";
	char binform_str[] = "binform";
	char ceil_str[] = "ceil";
	char testlog_str[] = "testlog";
	char resultant_str[] = "resultant";
	char discriminant_str[] = "discriminant";
	char deriv_str[] = "deriv";
	char primes_str[] = "primes";
	char sturmsequence_str[] = "sturmsequence";
	char cyclotomic_str[] = "cyclotomic";
	char classnop_str[] = "classnop";
	char classnon_str[] = "classnon";
	char nearint_str[] = "nearint";
	char reduceneg_str[] = "reduceneg";
	char reducepos_str[] = "reducepos";
	char classnop0_str[] = "classnop0";
	char tableneg_str[] = "tableneg";
	char tablepos_str[] = "tablepos";
	char davison_str[] = "davison";
	char raney_str[] = "raney";
	char unimodular_str[] = "unimodular";
	char twoadicsqrt_str[] = "twoadicsqrt";
	char padicsqrt_str[] = "padicsqrt";
	char sigmak_str[] = "sigmak";
	char ramanujan_str[] = "ramanujan";
	char repdefinite_str[] = "repdefinite";
	char powerd_str[] = "powerd";
	char euclid1_str[] = "euclid1";
	char rcf_period_str[] = "rcfperiod";
	char nscf_period_str[] = "nscfperiod";
	char alg_cycle_str[] = "algcycle";
	char spiral_str[] = "spiral";
	char spiralinverse_str[] = "spiralinverse";
	char nicf_period_str[] = "nicfperiod";
	char rcf_period0_str[] = "rcfperiod0";
	char nscf_period0_str[] = "nscfperiod0";
	char nicf_period0_str[] = "nicfperiod0";
	char carmichael_str[] = "carmichael";
	char carnielli_str[] = "carnielli";
	char lupei_str[] = "lupei";
	char tangent_str[] = "tangent";
	char bernoulli_str[] = "bernoulli";
	char partition_str[] = "partition";
	char twocycle_str[] = "twocycle";
	char mcycle_str[] = "mcycle";
	char pqcycle_str[] = "pqcycle";
	char mmcycle_str[] = "mmcycle";
	char exceptionals[] = "exceptionals";
	char exceptionals0[] = "exceptionals0";
	char fg[] = "fg";
	char exceptionals10[] = "exceptionals10";
	char gplus[] = "gplus";
	char gzero[] = "gzero";
	char gminus[] = "gminus";
	char nagell[] = "nagell";
	char frattini[] = "frattini";
	char kashihara[] = "kashihara";
	char branch[] = "branch";
	char ternary[] = "ternary";
	char pell4[] = "pell4";
	char lprimefactor[] = "lprimefactor";
	char minarray[] = "minarray";
	char conwaycycles[] = "conwaycycles";
	char conwayrangetest1[] = "conwayrangetest1";
	while (1)
	{
		printf("          **************************\n");
		printf("          * List of CALC functions *\n");
		printf("          **************************\n");
		printf("\n");
		printf("absmod, euclid, gcd, gcdv, gcda, gcdav, egcd, sgcd, lllgcd, lllgcd0, \n");
		printf("lllhermite, axb, fp, slv, inhomfp, lcm, lcma, length, \n");
		printf("pollard, nprime, nprimeap, jacobi, peralta, congr, chinese, chinesea, \n");
		printf("mthroot, mthrootr, fund, pell, surd, mpower, inv, elliptic, \n");
		printf("factor, tau, sigma, mobius, euler, lprimroot, orderm, lucas, \n");
		printf("serret, collatz, cycle, miller, hermite, improvep, mlll, smith, \n");
		printf("rsae, encode, decode, addcubicr, powercubicr, ordercubicr, \n");
		printf("addcubicm, powercubicm, ordercubicm, convergents, lagrange, \n");
		printf("perfectpower, leastqnr, content, primitive, sturm, rootexp, log, \n");
		printf("sqroot, cornacchia, patz, congq, binform, ceil, testlog,\n");
		printf("resultant, discriminant, deriv, primes, sturmsequence, cyclotomic,\n");
		printf("classnop, classnon, nearint, reduceneg, reducepos, classnop0,\n");
		printf("tableneg, tablepos, davison, raney, unimodular, twoadicsqrt,\n");
		printf("padicsqrt, sigmak, ramanujan, repdefinite, powerd, euclid1,\n");
		printf("rcfperiod, nscfperiod, nicfperiod, algcycle, spiral, spiralinverse,\n");
		printf("rcfperiod0, nscfperiod0, nicfperiod0, carmichael, carnielli, lupei,\n");
		printf("tangent, bernoulli, partition, twocycle, mcycle, pqcycle, mmcycle,\n");
		printf("exceptionals, exceptionals0, fg, exceptionals10, gplus, gzero, gminus,\n");
		printf("nagell, frattini, kashihara, branch, ternary, pell4, lprimefactor, minarray,\n");
		printf("conwaycycles,conwayrangetest1\n");
		printf("          **************************\n");
		printf("Type a function name or type qu to return to command line:");
		p = Fgets(stdin);
		strcpy(name, p);
		if (strncmp(name, quit_str, 2) == 0)
			break;
		else if (strcmp(name, euclid_str) == 0)
		{
			printf("          **************************\n");
			printf("Usage: euclid(a,b,&q[],&r[],&s[],&t[],&m)\n");
			printf("Description: m=n+1, where the arrays q[0]=NULL,...,q[n],q[n+1]=NULL,\n");
			printf("r[0]=a,r[1]=b,...,r[n+1],\n");
			printf("s[0]=1,s[1]=0,...,s[n+1],\n");
			printf("t[0]=0,t[1]=1,...,t[n+1],\n");
			printf("from Euclid's algorithm are returned and printed in euclid.out:\n");
			printf("r[k]=r[k+1]*q[k+1]+r[k+2], 0 < r[k+2] < r[k+1],\n");
			printf("s[k]=-q[k-1]*s[k-1]+s[k-2],\n");
			printf("t[k]=-q[k-1]*t[k-1]+t[k-2]\n");
			printf("r[n]=gcd(a,b)=s[n]*a+t[n]*b.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, absmod_str) == 0)
		{
			printf("          **************************\n");
			printf("Usage: z=absmod(a,b), b>= 1;\n");
			printf("z=r=a(mod b) if r<=b/2,\n");
			printf("z=r-b) if r>b/2\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, gcd_str) == 0)
		{
			printf("          **************************\n");
			printf("Usage: z=gcd(m,n)\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, gcdv_str) == 0)
		{
			printf("          **************************\n");
			printf("Usage: z=gcdv(m,n,&s,&t)\n");
			printf("As well as returning z=gcd(m,n),\n");
			printf("this gives numbers s and t satisfying z=s*m+t*n,\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, gcdav_str) == 0) {
			printf("\n********************************\n");
			printf("Usage: z=gcdav(a[],&b[])\n");
			printf("z = gcd(a[0],...,a[n-1]). Also gives integers\n");
			printf("b[0],...,b[n-1] satisfying z = b[0]a[0]+...+b[n-1]a[n-1].\n");
			printf("********************************\n");
			GetReturn();
		}
		else if (strcmp(name, gcda_str) == 0)
		{
			printf("          **************************\n");
			printf("Usage: z=gcda(a[])\n");
			printf("z=gcd(a[0],...,a[n-1]),\n");
			printf("where values for a[0],...,a[n-1]\n");
			printf("already exist.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, egcd_str) == 0)
		{
			printf("          **************************\n");
			printf("Usage: z=egcda()\n");

			printf("This does the same as gcda, except that\n");
			printf("one can input the numbers from a file as well.\n");
			printf("The file should have as its first line n, where n\n");
			printf("is the number of integers. The integers should be listed\n");
			printf("on separate lines.  A small multiplier vector P is found\n");
			printf("by LLL based Algorithm 1 of Havas, Majewski and Matthews.\n");
			printf("An m x m$ matrix  whose rows are X_1,...,X_{m-1},P \n");
			printf("is sent to an output file egcdmat.out.\n");
			printf("Here X_1,...,X_{m-1} form a LLL reduced basis\n");
			printf("for the lattice L defined by the equation\n");
			printf("x_1d_1+...+x_md_m=0.\n");
			printf("The inhomogeneous version of the Fincke--Pohst algorithm \n");
			printf("can then be used as an option to find a shortest \n");
			printf("multiplier vector by solving the inequality\n");
			printf("||P - x_1X_1-...-x_{m-1}X_{m-1}||^2 <= ||P||^2\n");
			printf("in integers x_1,...,x_{m-1}.\n");
			printf("Each time a shorter multiplier vector\n");
			printf("$Q=P-x_1X_1-...-x_{m-1}X_{m-1}$ is found,\n");
			printf("$P$ is replaced by $Q$, until the shortest $Q$ is found.\n");
			printf("The multipliers are sent to an output file called egcdmult.out.\n");
			printf("There is an option for finding all the shortest multipliers.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, sgcd_str) == 0)
		{
			printf("          **************************\n");
			printf("sgcd(): usage: sgcd(N)\n");
			printf("This performs the LLL algorithm on [I_n|NA],\n");
			printf("where A is a column vector of positive integers.\n");
			printf("If N is sufficiently large, the last column will be\n");
			printf("reduced to +-NdE_n, where d=gcd(a[1],...,a[n]).\n");
			printf("Output is sent to sgcdbas.out.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, lllgcd_str) == 0)
		{
			printf("          **************************\n");
			printf("Usage: lllgcd(). This performs a modification of LLL\n");
			printf("which is really a limiting form of sgcd(N) for large N.\n");
			printf("It is superior to egcd() in that it avoids inputting a\n");
			printf("large initial unimodular matrix and instead builds one\n");
			printf("from the identity matrix at the outset. \n");
			printf("The underlying algorithm is explained in a paper at\n");
			printf("http://www.maths.uq.edu.au/~krm/lll.html.\n");
			printf("The multipliers are sent to an output file called lllgcdmult.out.\n");
			printf("The matrix B of that paper is sent to lllgcdmat.out. \n");
			printf("Note: if the shortest vector option is chosen, the last\n");
			printf("row of B has been replaced by this vector.\n");
			printf("There is an option for finding all the shortest multipliers.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, lllgcd0_str) == 0)
		{
			printf("          **************************\n");
			printf("Usage: lllgcd0().  This sits on top of lllgcd()\n");
			printf("and finds a multiplier that in general is better than\n");
			printf("that delivered by lllgcd().\n");
			printf("See the paper at http://www.maths.uq.edu.au/~krm/lll.html.\n");
			printf("The outputs go to lllgcd0mat.out and lllgcd0mult.out.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, lllhermite_str) == 0)
		{
			printf("          **************************\n");
			printf("Usage: lllhermite().  This is a new LLL based algorithm\n");
			printf("for the Hermite normal form HNF(A) of A.\n");
			printf("HNF(A) is sent to lllhermitebas.out, \n");
			printf("while the unimodular transformation matrix P\n");
			printf("is sent to lllhermitetrans.out.\n");
			printf("See the paper at http://www.maths.uq.edu.au/~krm/lll.html.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, axb_str) == 0)
		{
			printf("          **************************\n");
			printf("Usage: axb(). This solves AX=B, where A is mxn,\n");
			printf("X is nx1, B is mx1. The method is MLLL based,\n");
			printf("as in an example of M. Pohst. A short basis is found\n");
			printf("for the nullspace of A and this is used to size-reduce\n");
			printf("an initial large solution.  The option of finding all \n");
			printf("the shortest X is also given. The short solution\n");
			printf("together with all shortest, if requested, are sent to\n");
			printf("axb.out, while the short basis for the nullspace of A\n");
			printf("is sent to axbbas.out.\n");
			printf("See the paper at http://www.maths.uq.edu.au/~krm/lll.html.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, fp_str) == 0)
		{
			printf("          **************************\n");
			printf("Usage: fp(), This is the simplest form of the Fincke-Pohst algorithm.\n");
			printf("This takes an integer matrix with LI rows and a positive integer C\n");
			printf("as input and finds all lattice vectors X with ||X||^2 <= C.\n");
			printf("The vectors found are sent to a file called fp.out. \n");
			printf("It is a good idea to start with a matrix whose rows are LLL reduced.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, inhomfp_str) == 0)
		{
			printf("          **************************\n");
			printf("Usage: inhomfp(). The inhomogeneous version of the Fincke-Pohst algorithm.\n");
			printf("Input: An m x M integer matrix A whose first m-1 rows are LI\n");
			printf("and which are preferably in LLL reduced form.\n");
			printf("A positive integer C is also entered.\n");
			printf("L is the lattice spanned by the first m-1 rows of A\n");
			printf("and P is the last row of A.\n");
			printf("Output: All lattice vectors X of L such that ||X-P||^2 <= C.\n");
			printf("The vectors found are sent to a file called inhomfp.out.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, slv_str) == 0)
		{
			printf("          **************************\n");
			printf("Usage: slv(). This finds the shortest vectors in the lattice\n");
			printf("spanned by the LI family b[1],...,b[m]. \n");
			printf("It applies Fincke-Pohst to examine all nonzero X in L\n");
			printf("satisfying ||X||^2<=C=||b[1]||^2. If it finds an X\n");
			printf("shorter than b[1], the new bound C=||X||^2 is chosen.\n");
			printf("At the end, Fincke-Pohst has only the shortest vectors to enumerate.\n");
			printf("Only the vectors with highest positive coefficient of b[i]  are given.\n");
			printf("The output is sent to slv.out.  \n");
			printf("It is suggested that the input matrix is LLL reduced.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, lcm_str) == 0)
		{
			printf("          **************************\n");
			printf("lcm: usage: z=lcm(x,y) \n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, lcma_str) == 0)
		{
			printf("          **************************\n");
			printf("lcma: usage: z=lcma(a[]) \n");
			printf("z=lcm(a[0],...,a[n-1]) (n is size of array.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, length_str) == 0)
		{
			printf("          **************************\n");
			printf("length: usage: z=length(n)\n");
			printf("z is the number of decimal digits of n.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, pollard_str) == 0)
		{
			printf("          **************************\n");
			printf("pollard: usage: z=pollard(x) \n");
			printf("This attempts to return a factor of a composite x\n");
			printf("using Pollard's p-1 method.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, nprime_str) == 0)
		{
			printf("          **************************\n");
			printf("nprime: usage: z=nprime(x)\n");
			printf("This finds the first integer after x which passes the\n");
			printf("strong base 2 pseudoprime test and the Lucas pseudoprime test.\n");
			printf("This integer is likely to be prime.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, nprimeap_str) == 0)
		{
			printf("          **************************\n");
			printf("nprimeap: usage: z=nprimeap(a,b,m)\n");
			printf("This finds the first p, p=b(mod a), m <= p,\n");
			printf("which passes the strong base 2 pseudoprime test\n");
			printf("and the Lucas pseudoprime test. \n");
			printf("Here a must be even, b odd, 1 <= b < a, gcd(a,b)=1, b <= m. \n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, jacobi_str) == 0)
		{
			printf("          **************************\n");
			printf("jacobi: usage: z=jacobi(x,y)\n");
			printf("z is the value of the Jacobi symbol (x/y).\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, peralta_str) == 0)
		{
			printf("          **************************\n");
			printf("peralta: usage: z=peralta(a,p)\n");
			printf("Peralta's algorithm is used to return a square root z of a (mod p).\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, cong_str) == 0)
		{
			printf("          **************************\n");
			printf("congr: usage: x=congr(a,b,m,&n)\n");
			printf("Returns the solution x of the congruence ax=b(mod m).\n");
			printf("Also n=m/gcd(a,m) is returned.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, chinese_str) == 0)
		{
			printf("          **************************\n");
			printf("chinese: usage: x=chinese(a,b,m,n,&l)\n");
			printf("Returns the solution x(mod l) of the system of congruences\n");
			printf("x=a(mod m) and x=b(mod n).\n");
			printf("Also l=lcm(m,n) is returned.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, chinesea_str) == 0)
		{
			printf("          **************************\n");
			printf("chinesea: usage: x=chinesea(a[],m[],&l)\n");
			printf("Returns the solution x(mod l) of the system of congruences\n");
			printf("x=a[i](mod m[i]), 0<= i<n. (n = smaller of two array sizes)\n");
			printf("Also l=lcm(m[0],...,m[n-1]) is returned.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, mthroot_str) == 0)
		{
			printf("          **************************\n");
			printf("mthroot: usage: z=mthroot(x,m)\n");
			printf("The integer part of the m-th root of x is returned.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, mthrootr_str) == 0)
		{
			printf("          **************************\n");
			printf("mthrootr: usage: mthrootr(x,y,m,r)\n");
			printf("The m-th root of x/y is computed to r decimal places.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, fund_str) == 0)
		{
			printf("          **************************\n");
			printf("fund: usage: z=fund(d,&x,&y)\n");
			printf("x+y*omega is the fundamental unit of Q(sqrt{d}).\n");
			printf("z=Norm(x+y*omega) is also returned.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, pell_str) == 0)
		{
			printf("          **************************\n");
			printf("pell: usage: z=pell(d,e,&x,&y)\n");
			printf("The continued fraction expansion of sqrt{d}\n");
			printf("is periodic after the first term:\n");
			printf("a[0],a[1],...,a[n-1],2a[0],a[1],...,a[n-1],2a[0],....\n");
			printf("Also the section a[1],...,a[n-1] is palindromic.\n");
			printf("We print a[0] and half the palindrome iff e is nonzero,\n");
			printf("sending the output to a file called pell.out.\n");
			printf("The least solution x and y of Pell's equation\n");
			printf("x^2-dy^2=epsilon is returned, where epsilon=1 or -1.\n");
			printf("z=epsilon is also returned.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, surd_str) == 0)
		{
			printf("          **************************\n");
			printf("surd: usage: z=surd(d,t,u,v,&a[],&u[],&v[],&p[],&q[])\n");
			printf("The continued fraction of (u+t*sqrt{d})/v\n");
			printf("is determined up to the end of the period.\n");
			printf("a[i] is the ith partial quotient.\n");
			printf("u[i]+sqrt(d))/v[i] is the ith complete convergent.\n");
			printf("p[i]/q[i] is the ith convergent.\n");
			printf("The period length z is returned\n");
			printf("Output is sent to file surd.out.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, mpower_str) == 0)
		{
			printf("          **************************\n");
			printf("mpower: usage: z=mpower(a,b,c)\n");
			printf("z=a^b(mod c) is returned.\n");
			printf("a,b,c integers, b>=0,c>0.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, inv_str) == 0)
		{
			printf("          **************************\n");
			printf("inv: usage: z=inv(a,m)\n");
			printf("z=a^{-1}(mod m) is returned.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, elliptic_str) == 0)
		{
			printf("          **************************\n");
			printf("elliptic: usage: z=elliptic(n,m,p)\n");
			printf("The elliptic curve method is used to try to find\n");
			printf("a factor z of a composite number n.\n");
			printf("Here 1 <= p < 2^32 and 10 < m < 1279.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, factor_str) == 0)
		{
			printf("          **************************\n");
			printf("factor: usage: z=factor(n)\n");
			printf("The factorization of n is performed using the\n");
			printf("multiple polynomial quadratic sieve if length(n) <= 55;\n");
			printf("otherwise the elliptic curve method is used.\n");
			printf("z=omega(n), the number of distinct prime factors of n.\n");
			printf("One can nominate the number of elliptic curves to be used.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, tau_str) == 0)
		{
			printf("          **************************\n");
			printf("tau: usage: z=tau(n)\n");
			printf("The divisor function z=tau(n) is returned.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, sigma_str) == 0)
		{
			printf("          **************************\n");
			printf("sigma: usage: z=sigma(n)\n");
			printf("z=sigma(n), the sum of the divisors of n, is returned.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, mobius_str) == 0)
		{
			printf("          **************************\n");
			printf("mobius: usage: z=mobius(n)\n");
			printf("The Moebius function z=mu(n) is returned.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, euler_str) == 0)
		{
			printf("          **************************\n");
			printf("euler: usage: z=euler(n)\n");
			printf("Euler's function z=phi(n) is returned.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, lprimroot_str) == 0)
		{
			printf("          **************************\n");
			printf("lprimroot: usage: z=lprimroot(n)\n");
			printf("The least primitive root mod p, an odd prime, is returned.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, orderm_str) == 0)
		{
			printf("          **************************\n");
			printf("orderm: usage: z=orderm(a,n)\n");
			printf("The order of a(mod n) is returned.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, lucas_str) == 0)
		{
			printf("          **************************\n");
			printf("lucas: usage: z=lucas(n)\n");
			printf("n is subjected to a strong pseudoprime test to base 2,\n");
			printf("together with a Lucas pseudoprime test.\n");
			printf("If z=0 is returned, n is composite;\n");
			printf("if z=1 is returned, n is a Lucas probable prime,\n");
			printf("as well as a base 2 strong pseudoprime.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, serret_str) == 0)
		{
			printf("          **************************\n");
			printf("serret: usage: serret(p,&x,&y)\n");
			printf("Here p is a prime of the form 4n+1.\n");
			printf("Serret's algorithm finds integers x and y such that p=x^2+y^2. \n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, collatz_str) == 0)
		{
			printf("          **************************\n");
			printf("collatz: usage: collatz(x, n)\n");
			printf("Collatz' 3x+1 conjecture is tested.\n");
			printf("The iterates x,T(x),.. are printed iff n is nonzero.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, cycle_str) == 0)
		{
			printf("          **************************\n");
			printf("cycle: usage cycle()\n");
			printf("Tries to find all cycles for the generalised Collatz mapping T,\n");
			printf("which arise from starting numbers p, |p|<= RANGE/2.\n");
			printf("Trajectories containing an iterate whose magnitude\n");
			printf("exceeds a prescribed value INFINITY are deemed non-cycling.\n");
			printf("T(x)=int(m[i]*x/d)+X[i] if x=i (mod d).\n");
			printf("The d nonzero moduli m[i] and d shifts X[i] are entered from the keyboard.\n");
			printf("(See http://www.numbertheory.org/3x+1/)\n");
			printf("Output is sent to cycle.out\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, miller_str) == 0)
		{
			printf("          **************************\n");
			printf("miller: usage: miller(m,b)\n");
			printf("Here m>1,b>1 and m does not divide b.\n");
			printf("Miller's test to base b is applied to m.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, hermite_str) == 0)
		{
			printf("          **************************\n");
			printf("hermite: usage: hermite()\n");
			printf("The Hermite normal form HNF(A) of an integer matrix A\n");
			printf("is found using the Kannan-Bachem algorithm.\n");
			printf("A can be entered either from the keyboard or from a file\n");
			printf("such as \n");
			printf("        2 3\n");
			printf("        1 -4 5\n");
			printf("        3  2 1\n");
			printf("in the case of a 2 x 3 matrix.\n");
			printf("The hermite normal form is sent to a file called hermite.out.\n");
			printf("A unimodular integer matrix P such that PA = HNF(A) is sent to hermitep.out, if desired.\n");
			printf("The last line of this file contains the value of rank(A).\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, improvep_str) == 0)
		{
			printf("          **************************\n");
			printf("improvep: usage: improvep()\n");
			printf("The unimodular matrix contained in the file hermitep.out is improved\n");
			printf("using the LLL based method, followed by Gauss lattice reduction.\n");
			printf("The output is sent to a file improvep.out.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, mlll_str) == 0)
		{
			printf("          **************************\n");
			printf("mlll: usage: mlll()\n");
			printf("The MLLL algorithm of M. Pohst is applied to an integer matrix\n");
			printf("whose first row is non-zero.\n");
			printf("The reduced matrix is sent to mlllbas.out\n");
			printf("the transforming matrix is sent to mllltran.out.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, smith_str) == 0)
		{
			printf("          **************************\n");
			printf("smith: usage: smith()\n");
			printf("The Smith normal form SNF(A) of an integer matrix A is found\n");
			printf("using a pivoting strategy due to G. Havas.\n");
			printf("A cutoff value is requested. When coefficients grow\n");
			printf("above this value in size, the MLLL algorithm is used to\n");
			printf("reduce coefficient explosion. Unimodular matrices P and Q\n");
			printf("are found such that PAQ=SNF(A).\n");
			printf("The invariant factors d[1],...,d[r],\n");
			printf("together with P and Q, are sent to files\n");
			printf("smith.out, smithp.out, smithq.out, respectively, if desired.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, rsae_str) == 0)
		{
			printf("          **************************\n");
			printf("rsae: usage: e=rsae(p,q)\n");
			printf("Here p and q are distinct odd primes\n");
			printf("and e is the least integer such that\n");
			printf("gcd(e,(p-1)(q-1))=1 and 32^e>pq.\n");
			printf("e is for use as an RSA encryption modulus.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, encode_str) == 0)
		{
			printf("          **************************\n");
			printf("encode: usage: encode(e,n)\n");
			printf("n=p*q, a product of primes each greater than 355142.\n");
			printf("(p*q>126126126126.)\n");
			printf("e is the RSA encryption modulus, found using rsae() above.\n");
			printf("A string of non-control characters when entered from the keyboard\n");
			printf("or from a file consisting of lines each containing less than 500 characters.\n");
			printf("These characters have ascii values in the range 32-126.\n");
			printf("The message string is encoded using the RSA algorithm:\n");
			printf("every 4 characters are converted to ascii, joined as strings\n");
			printf("and the resulting large number m is encoded as n=m^e(mod p*q).\n");
			printf("The encoded numbers are sent to encoded.out, which is terminated by an entry -1.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, decode_str) == 0)
		{
			printf("          **************************\n");
			printf("decode: usage: decode(e,p,q)\n");
			printf("This calculates the decryption modulus d,\n");
			printf("then decodes each number n in encoded.out.\n");
			printf("m=n^d(mod pq). The ascii characters are split off\n");
			printf("and the original string of message characters is recreated.\n");
			printf("The decoded message is sent to decoded.out, as well as to the screen.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, addcubicr_str) == 0)
		{
			printf("          **************************\n");
			printf("addcubicr: usage: addcubicr()\n");
			printf("The sum of two points on the elliptic curve y^2+a_1xy+a_3y=x^3+a_2x^2+a_4x+a_6 over Q\n");
			printf("is calculated. The discriminant is also calculated.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, powercubicr_str) == 0)
		{
			printf("          **************************\n");
			printf("powercubicr: usage: powercubicr()\n");
			printf("The point nP, where P is on the elliptic curve over Q:\n");
			printf("y^2+a_1xy+a_3y=x^3+a_2x^2+a_4x+a_6, is calculated.\n");
			printf("The discriminant is also calculated.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, ordercubicr_str) == 0)
		{
			printf("          **************************\n");
			printf("ordercubicr: usage: ordercubicr()\n");
			printf("Finds the order of P on the elliptic curve over Q:\n");
			printf("y^2+a_1xy+a_3y=x^3+a_2x^2+a_4x+a_6.\n");
			printf("The discriminant is also calculated.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, addcubicm_str) == 0)
		{
			printf("          **************************\n");
			printf("addcubicm: usage: addcubicm()\n");
			printf("The sum of two points on the elliptic curve\n");
			printf("y^2+a_1xy+a_3y=x^3+a_2x^2+a_4x+a_6 over Z_p\n");
			printf("is calculated. The discriminant is also calculated.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, powercubicm_str) == 0)
		{
			printf("          **************************\n");
			printf("powercubicm: usage: powercubicm()\n");
			printf("The point nP, where P is on the elliptic curve over Z_p:\n");
			printf("y^2+a_1xy+a_3y=x^3+a_2x^2+a_4x+a_6, is calculated.\n");
			printf("The discriminant is also calculated.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, ordercubicm_str) == 0)
		{
			printf("          **************************\n");
			printf("ordercubicm: usage: ordercubicm()\n");
			printf("Finds the order of P on the elliptic curve over Z_p:\n");
			printf("y^2+a_1xy+a_3y=x^3+a_2x^2+a_4x+a_6.\n");
			printf("The discriminant is also calculated.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, convergents_str) == 0)
		{
			printf("          **************************\n");
			printf("convergents: usage: convergents(a[],&p[],&q[])\n");
			printf("The convergents p[0]/q[0],...,p[n]/q[n] of the continued fraction\n");
			printf("[a[0];a[1],...,a[n]] are returned as arrays p[] and q[]\n");
			printf("and are sent to a file convergents.out.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, lagrange_str) == 0)
		{
			printf("          **************************\n");
			printf("Usage: lagrange(\"poly\",&a[],m)\n");
			printf("""poly"" is a polynomial\n");
			printf("with integer coefficients, having no rational roots\n");
			printf("and having exactly one real positive root x > 1.\n");
			printf("The method of Lagrange (1797) is used to find the\n");
			printf("the first m+1 partial quotients a[0],...,a[m] of x.\n");
			printf("(See Knuth, Art of computer programming, volume2, problem 13, 4.5.3.\n");
			printf("Also S. Lang and H. Trotter, 'Continued fractions for some algebraic numbers',\n");
			printf("J. fur Math. 255 (1972) 112-134; Addendum 267 (1974) ibid. 219-220.\n");
			printf("Also D.G. Cantor, P.H. Galyean, H.G. Zimmer,\n");
			printf("'A continued fraction algorithm for real algebraic numbers',\n");
			printf("Math. Comp. 26, July 1972, 785-791.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, perfectpower_str) == 0)
		{
			printf("          **************************\n");
			printf("perfectpower: usage: z=perfectpower(n)\n");
			printf("Here n > 1.  z = x if n = x^k, for some x, k > 1,\n");
			printf("NULL otherwise.\n");
			printf("See E. Bach and J. Sorenson,\n");
			printf("'Sieve algorithms for perfect power testing'\n");
			printf("Algorithmica 9 (1993) 313-328.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, leastqnr_str) == 0)
		{
			printf("          **************************\n");
			printf("leastqnr: usage: z=leastqnr(p):\n");
			printf("Returns n_p, the least quadratic non-residue (mod p),\n");
			printf("if n_p< 65536, otherwise returns 0.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, examples_str) == 0)
		{
			printf("          **************************\n");
			printf("SOME EXAMPLES\n");
			printf("To find z=gcd(4,6), type z=gcd(4,6)\n");
			printf("Then z is stored and its value is printed on typing z,\n");
			printf("as follows:\n");
			printf("z\n");
			printf("2\n");
			GetReturn();

			printf("To find z=gcd(4,6) together with integers u,v\n");
			printf("satisfying z=4u+6v, type z=gcdv(4,6,&u,&v)\n");
			printf("The values of u,v are stored and their values printed,\n");
			printf("by typing u,v, as in \n");
			printf("u\n");
			printf("-1\n");
			printf("v\n");
			printf("1\n");
			GetReturn();

			printf("To find z=gcd(4,6,9), type\n");
			printf("x[0]=4;x[1]=6;x[2]=9\n");
			printf("z=gcda(x[])\n");
			GetReturn();

			printf("To find z=gcd(4,6,9), together with integers\n");
			printf("b[0],b[1],b[2] satisfying z=4*b[0]+6*b[1]+9*b[2],\n");
			printf("type\n");
			printf("x[0]=4;x[1]=6;x[2]=9\n");
			printf("z=gcdav(x[],&b[])\n");
			printf("The values of b[0],b[1],b[2] are stored\n");
			printf("and their values printed by typing printa(b,2) \n");
			GetReturn();
			printf("\nCalc can also parse polynomials.\n");
			printf("Typing\n> (X+2)^2\ngives\nX^2 + 4X + 4\n");
			printf("If you assign the polynomial to a variable with\n");
			printf("> z=X^3+3X^2+2X+1\nyou can evaluate the polynomial\n");
			printf("at various integer points by typing (for example)\n");
			printf("> z(2)\nwhich gives\n25\n");
			printf("One can also enter arrays in the following manner\n");
			printf("> a[] = {1, 2, 3, 4, 5, 6}\n");
			printf("Note: Any previous definition of the array will be");
			printf("erased\n");
			printf("Or one can enter arrays like so:\n");
			printf("> a[1]=2\n> a[3]=2\n> a[7]=1000\n");
			printf("One can view the contents of the array by typing\n");
			printf("> a[]\nwhich gives\n[0]:0\n[1]:2\n[2]:0\n[3]:2\n");
			printf("[4]:0\n[5]:0\n[6]:0\n[7]:1000\n");
			printf("Note that the undefined subscripts have been\n");
			printf("initialised to zero.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, content_str) == 0) {
			printf("\n***************************************\n");
			printf("Usage: content(f(X))\n");
			printf("Returns the content of a polynomial f(X).\n");
			printf("eg. content(3X^2+6) returns 3\n");
			printf("*****************************************\n");
			GetReturn();
		}
		else if (strcmp(name, primitive_str) == 0) {
			printf("\n***************************************\n");
			printf("Usage: primitive(f(X))\n");
			printf("Returns the primitive part of a polynomial f(X).\n");
			printf("eg. primitive(3X^2+6) returns X^2+2\n");
			printf("*****************************************\n");
			GetReturn();
		}
		else if (strcmp(name, sturm_str) == 0) {
			printf("\n***************************************\n");
			printf("Usage: sturm(f(X))\n");
			printf("Prints rational open intervals that are guaranteed to\n");
			printf("contain only one real root of f(X).\n");
			printf("It is not guaranteed to produce correct results if the\n");
			printf("polynomial entered has rational roots.  The output is \n");
			printf("sent to the file sturm.out and to the screen.\n");
			printf("***************************************\n");
			GetReturn();
		}
		else if (strcmp(name, rootexp_str) == 0) {
			printf("\n***************************************\n");
			printf("Usage: rootexp(f(X), m)\n");
			printf("Finds m partial quotients of the continued fraction \n");
			printf("expansion of all real roots of a squarefree polynomial f(X)\n");
			printf("(having no rational roots)\n");
			printf("using Lagrange's method and methods presented in a\n");
			printf("paper by Cantor, Galyean and Zimmer called A Continued \n");
			printf("Fraction Algorithm for Real Algebraic Numbers.  The \n");
			printf("output is sent to a file rootexp.out and to the screen.\n");
			printf("***************************************\n");
			GetReturn();
		}
		else if (strcmp(name, log2_str) == 0) {
			printf("\n***************************************\n");
			printf("Usage: log2(a,b,d,r,s,e)\n");
			printf("This is a discrete version of Shank's log_b(a) algorithm.\n");
			printf("d > 1, 1< =r < s are integers. The larger is s-r, the more\n");
			printf("accurate the output. eg s=2r, but this is slow when r=500.\n");
			printf("But s=r+10 should be adequate.\n");
			printf("A certain number (usually at least r if d=10) of \n");
			printf("number of partial quotients are printed.\n");
			printf("The correctness is not 100 percent guaranteed.\n");
			printf("e=0 prints only the partial quotients,\n");
			printf("while e nonzero prints convergents and decimal expansion.\n");
			printf("Output is sent to log2.out.\n");
			printf("***************************************\n");
			GetReturn();
		}
		else if (strcmp(name, log_str) == 0) {
			printf("\n***************************************\n");
			printf("Usage: log(a,b,d,r,&a[],&l)\n");
			printf("This is a discrete version of Shank's log_b(a) algorithm.\n");
			printf("a>b>1, d>1, r>0 are integers. \n");
			printf("l  is the number of partial quotients a[s].\n");
			printf("The correctness is not 100 percent guaranteed.\n");
			printf("Output is sent to log.out.\n");
			printf("See http://www.numbertheory.org/pdfs/log.pdf.\n");
			printf("***************************************\n");
			GetReturn();
		}
		else if (strcmp(name, testlog_str) == 0) {
			printf("\n***************************************\n");
			printf("Usage: testlog(a,b,d,m,n)\n");
			printf("This runs the modified Shank's log_b(a) algorithm\n");
			printf("for r=m,...,n. \n");
			printf("Here a>b>1, d>1 are integers. \n");
			printf("Output is sent to testlog.out.\n");
			printf("With d=10 and m=n, we expect to get at least m partial quotients\n");
			printf("The idea is to run it for a range (m-t,m+t) with d=10\n");
			printf("to get a good idea of the correct partial quotients.\n");
			printf("See http://www.numbertheory.org/pdfs/log.pdf.\n");
			printf("***************************************\n");
			GetReturn();
		}
		else if (strcmp(name, log1_str) == 0) {
			printf("\n***************************************\n");
			printf("Usage: log1(a,b,d,r,e), where d>=2, r>=1.\n");
			printf("d=2 is fast, but d = say 10 gives more output.\n");
			printf("This is Algorithm 1 of http://www.maths.uq.edu.au/~krm/log.pdf.\n");
			printf("A certain number (approximately r) of \n");
			printf("number of partial quotients are printed.\n");
			printf("The correctness is not 100 percent guaranteed.\n");
			printf("So far a counter-example has not emerged.\n");
			printf("e=0 prints only the partial quotients,\n");
			printf("while e nonzero prints convergents, the integers A[i] and  the decimal expansion of log_b(a).\n");
			printf("Output is sent to log1.out.\n");
			printf("***************************************\n");
			GetReturn();
		}
		else if (strcmp(name, sqroot_str) == 0) {
			printf("\n***************************************\n");
			printf("Usage: z=sqroot(a,n,&s[],&m,&l)\n");
			printf("Returns the solutions of x^2=a (mod n)\n");
			printf("as x=s[i] or -s[i] (mod m), 0 <= s[i] <= m/2.\n");
			printf("z is the number of solutions mod n.\n");
			printf("If there is no solution, z=0 is returned\n");
			printf("together with s[0]=NULL and m=NULL.\n");
			printf("l is the number of solutions s[i] returned.\n");
			printf("***************************************\n");
			GetReturn();
		}
		else if (strcmp(name, cornacchia_str) == 0) {
			printf("\n***************************************\n");
			printf("Usage: z=cornacchia(a,b,m)\n");
			printf("Returns the positive primitive solutions (x,y) of a*x^2+b*y^2=m, where if a=b=1, x>=y.\n");
			printf("Here a>0,b>0,m>a+b, gcd(a,b)=1=gcd(a,m).\n");
			printf("***************************************\n");
			GetReturn();
		}
		else if (strcmp(name, patz_str) == 0) {
			printf("\n***************************************\n");
			printf("Usage: patz(d,n)\n");
			printf("Returns the positive, primitive fundamental solutions (x,y) of x^2-d*y^2=n and -n,\n");
			printf("in the case of solubility, where d>0 and not a perfect square.\n");
			printf("Also output is sent to patz.out.\n");
			printf("***************************************\n");
			GetReturn();
		}
		else if (strcmp(name, congq_str) == 0) {
			printf("\n***************************************\n");
			printf("Usage: z=congq(a,b,c,n,&s[])\n");
			printf("Solves the congruence ax^2+bx+c=0(mod n), a nonzero, n>0.\n");
			printf("z= number of solutions mod n.\n");
			printf("***************************************\n");
			GetReturn();
		}
		else if (strcmp(name, binform_str) == 0) {
			printf("\n***************************************\n");
			printf("Usage: binform(a,b,c,n,e)\n");
			printf("Solves the diophantine equation ax^2+bxy+cy^2=n, n non-zero,\n");
			printf("where D=b^2-4ac>0 and is not a perfect square.\n");
			printf("One solution from each class is printed,\n");
			printf("together with the corresponding solution n\n");
			printf("of the congruence n^2=D (mod 4|N|), -|N|<n<=|N|.\n");
			printf("e=1 is verbose, e=0 is terse\n");
			printf("***************************************\n");
			GetReturn();
		}
		else if (strcmp(name, ceil_str) == 0)
		{
			printf("          **************************\n");
			printf("ceil: usage: ceil(a,b)\n");
			printf("Returns the least integer not less than a/b\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, resultant_str) == 0)
		{
			printf("          **************************\n");
			printf("resultant: usage: r=resultant(p,q)\n");
			printf("Returns the resultant of non-constant polynomials p and q\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, discriminant_str) == 0)
		{
			printf("          **************************\n");
			printf("discriminant: usage: r=discriminant(p)\n");
			printf("Returns the discriminant of polynomial p\n");
			printf("not of the from aX+b\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, deriv_str) == 0)
		{
			printf("          **************************\n");
			printf("deriv: usage: q=deriv(p)\n");
			printf("Returns the derivative of polynomial p\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, primes_str) == 0)
		{
			printf("          **************************\n");
			printf("primes: usage: c=primes(m,n)\n");
			printf("where 1<=m<10^10, 1<=n<10^10, \n");
			printf("prints the primes if any, in the interval [m,n]\n");
			printf("and returns their number\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, sturmsequence_str) == 0)
		{
			printf("          **************************\n");
			printf("primes: usage: c=sturmdequence(f,b,e)\n");
			printf("where f is a polynomial, squarefree, deg(f)>1\n");
			printf("Also f(b) is non-zero.\n");
			printf("The Sturm polynomials are listed, along with\n");
			printf("their values at x=b. c is the no. of sign-changes.\n");
			printf("e=0 suppresses printing.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, cyclotomic_str) == 0)
		{
			printf("          **************************\n");
			printf("cyclotomic: usage: p=cyclotomic(n)\n");
			printf("where 1 <= n < 65536\n");
			printf("Returns the nth cyclotomic polynomial\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, classnop_str) == 0)
		{
			printf("          **************************\n");
			printf("classnop: usage: h=classnop(d)\n");
			printf("where 1 < d < 10^6 is squarefree.\n");
			printf("Returns the class-number of the real quadratic field\n");
			printf("Q(sqrt(d) and the sign of the fundamental unit.\n");
			printf("A complete set of reduced binary forms is given\n");
			printf("corresponding to the classes of ideals.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, classnon_str) == 0)
		{
			printf("          **************************\n");
			printf("classnon: usage: h=classnon(d,e)\n");
			printf("where d < 0 and 1 < |d| < 10^6 is squarefree\n");
			printf("and d=0 or 1(mod 4).\n");
			printf("    This is Henri Cohen's Algorithm 5.3.5, p. 228,\n");
			printf("for finding the class number h(d) of binary quadratic forms\n");
			printf("of discriminant d, when d<0.\n");
			printf("If e=1, we print only the primitive forms.\n");
			printf("h(d) is returned in each case.\n\n");
			printf("If d is the discriminant of an imaginary quadratic field K,\n");
			printf("then the primitive forms class-number h(d) is also\n");
			printf("the class number of K.\n\n");
			printf("    Davenport's Higher Arithmetic has a table of forms,\n");		printf("which lists the imprimitive ones with an asterisk.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, nearint_str) == 0)
		{
			printf("          **************************\n");
			printf("Usage: z=nearint(a,b), b > 0;\n");
			printf("z the nearest integer t to a/b,\n");
			printf("where z = t if a/b = 1/2 + t, t an integer.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, reduceneg_str) == 0)
		{
			printf("          **************************\n");
			printf("Usage: n=reduceneg(a,b,c)), a > 0, c > 0 and b^2-4ac < 0.\n");
			printf("This is Gauss's algorithm for reducing a positive\n");
			printf("definite binary quadratic form. See L.E. Dickson, \n");
			printf("Introduction to the theory of numbers, page 69. \n");
			printf("The reduced form (A,B,C) satisfies -A<B<=A, C>=A, \n");
			printf("with B>=0 if C=A.\n");
			printf("The number of steps taken in the reduction is returned.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, reducepos_str) == 0)
		{
			printf("          **************************\n");
			printf("Usage: n=reducepos(a,b,c), d=b^2-4ac > 0, d not a square.\n");
			printf("We also assume that d < 10^6.\n");
			printf("We use the PQa continued fraction algorithm to find\n");
			printf("an equivalent reduced form and thence a cycle of reduced forms.\n");
			printf("A unimodular tranforming matrix transforming (a,b,c) to a reduced form is constructed.\n");
			printf("The cycle-length is returned. A file reducepos.out is also created.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, classnop0_str) == 0)
		{
			printf("          **************************\n");
			printf("classnop0: usage: h=classnop0(d)\n");
			printf("where 1 < d < 10^6 is not a perfect square.\n");
			printf("d = 0 or 1 (mod 4).\n");
			printf("h is the number of classes of binary quadratic forms of discriminant d\n");
			printf("A complete set of reduced binary forms is given.\n");
			printf("We determine if the Pell equation x^2-d*y^2=-4 has\n");
			printf("has a solution, by using the fact that the equation\n");
			printf("is soluble iff at least one of the above cycles is odd.\n");
			printf("If there is no solution, the reduced forms (-a,b,-c)\n");
			printf("have to be counted as well. \n");
			printf("(See G.B. Mathews, Theory of Numbers, 80-81.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, tableneg_str) == 0)
		{
			printf("          **************************\n");
			printf("Usage: h=tableneg(m,n),1<=m<=n<10^6.\n");
			printf("Calculates h(-d) for all squarefree d with m<=d<=n.\n");
			printf("The number of squarefree d  in the range is returned.\n");
			printf("The output is also sent to tableneg.out.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, tablepos_str) == 0)
		{
			printf("          **************************\n");
			printf("Usage: h=tablepos(m,n),2<=m<=n<10^6.\n");
			printf("Calculates h(d) for all squarefree d with m<=d<=n.\n");
			printf("The number of squarefree d  in the range is returned.\n");
			printf("The output is also sent to tablepos.out.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, davison_str) == 0)
		{
			printf("          **************************\n");
			printf("Usage: h=davison(l,m,n),l,m>=1, 10^5>=n>=0\n");
			printf("h partial quotients a[i] of e^{l/m} are found.\n");
			printf("We cannot predict the value of h.\n");
			printf("The a[i] are also sent to davison.out.\n");
			printf("The program stops if 10^6 partial quotients a[i] are found.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, raney_str) == 0)
		{
			printf("          **************************\n");
			printf("Usage: t=raney(p,q,r,s),p,q,r,s>=0, p*s-q*r nonzero.\n");
			printf("Let A=[p,q;r,s]. Then we assume A!=I_2, A!=[0,1;1,0].\n");
			printf("With L=[1,0;1,1] and R=[1,1;0,1], we express A uniquely as\n");
			printf("a product of non-negative powers of L and R, (at least one is positive).\n");
			printf("followed by a row-balanced B.\n");
			printf("B=[a,b;c,d] is row-balanced if (a<c & b>d) or (c<a & d>b) and a,b,c>=0.\n");
			printf("The number k of powers of L and R is returned.\n");
			printf("The output is also sent to raney.out.\n");
			printf("See G.N. Raney, Math, Annalen 206 (1973) 265-283.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, unimodular_str) == 0)
		{
			printf("Usage: t=unimodular(p,q,r,s),p,q,r,s>=0, p*s-q*r = 1 or -1.\n");
			printf("This algorithm expresses a unimodular matrix \n");
			printf("A !=I_2 or U=[0,1;1,0] with non-negative coefficients\n");
			printf("as a product of one of the following forms:\n"); printf("P, UP, PU, or UPU, where P is a product of matrices\n");
			printf("of the form U[a]=[a,1;1,0], a > 0.\n");
			printf("The representation is unique. \n");
			printf("See Kjell Kolden, 'Continued fractions and linear substitutions'\n");
			printf("Arch. Math. Naturvid. 50 (1949), 141-196.\n");
			printf("The number t of matrices in the product is returned.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, twoadicsqrt_str) == 0)
		{
			printf("Usage: twoadicsqrt(b,n), b>0, b=8k+1, n > 0.\n");
			printf("Finds n terms a[0]...a[n-1] of the 2-adic square root x of a, x=1 (mod 4).\n");
			printf("Output also sent to 2-adic.out.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, padicsqrt_str) == 0)
		{
			printf("Usage: padicsqrt(b,p,n.&a[]), b>0, b a quadratic residue (mod p), n > 0.\n");
			printf("Finds a square root u of b (mod p), 0 < u < p.\n");
			printf("Then finds n terms a[0]...a[n-1] of the p-adic square root x of b, x=u (mod p).\n");
			printf("Output also sent to p-adic.out.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, sigmak_str) == 0)
		{
			printf("Usage: u=sigmak(k,n), 2^16>k>0, n > 0.\n");
			printf("Returns the sum of the kth powers of the divisors of n.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, ramanujan_str) == 0)
		{
			printf("Usage: u=ramanujan(n), 2^16>n>0.\n");
			printf("Returns Ramanujan's tau function.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, repdefinite_str) == 0)
		{
			printf("          **************************\n");
			printf("Usage: repdefinite(a,b,c,m,print_flag), a > 0, c > 0 and b^2-4ac < 0.\n");
			printf("This algorithm of Gauss solves the diophantine equation\n");
			printf("ax^2+bxy+cy^2=m, where d=b^2-4ac<0, a>0, c>0, m>0.\n");
			printf("See  L.E. Dickson, Introduction to the Theory of Numbers, 74-75 \n");
			printf("output is sent to repdefinite.out\n");
			printf("print_flag=0 lists only the solutions,\n");
			printf("print_flag=1 lists unimodular transformations\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, powerd_str) == 0)
		{
			printf("          **************************\n");
			printf("powerd: usage: powerd(a,b,d,n,&aa,&bb)\n");
			printf("(a+b*sqrt(d))^n=aa+bb*sqrt(d) is returned.\n");
			printf("a,b,c integers, d>0,n>=0.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, euclid1_str) == 0)
		{
			printf("          **************************\n");
			printf("euclid1: usage: z=euclid1(a,b)\n");
			printf("z is the length of Euclid's algorithm for a /b.\n");
			printf("a and b are positive integers.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, nscf_period_str) == 0)
		{
			printf("          **************************\n");
			printf("nscfperiod: usage: z=nscfperiod(d)\n");
			printf("The period-length z of the nearest square\n");
			printf("continued fraction expansion of sqrt{d}\n");
			printf("is returned");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, rcf_period_str) == 0)
		{
			printf("          **************************\n");
			printf("rcfperiod: usage: z=rcfperiod(d)\n");
			printf("The period-length z of the regular\n");
			printf("continued fraction expansion of sqrt{d}\n");
			printf("is returned");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, nicf_period_str) == 0)
		{
			printf("          **************************\n");
			printf("nicfperiod: usage: z=nicfperiod(d)\n");
			printf("The period-length z of the nearest integer\n");
			printf("continued fraction expansion of sqrt{d}\n");
			printf("is returned");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, alg_cycle_str) == 0)
		{
			printf("          **************************\n");
			printf("algcycle: usage: algcycle()\n");
			printf("We search for cycles of the d-branched function\n");
			printf("obtained by dividing by sqrt(d), d > 0.\n");
			printf("x=u+v*sqrt(d) -> ");
			printf("((M[i]+N[i]sqrt(d))x - (X[i]+Y[i]sqrt(d)))/sqrt(d)\n");
			printf("if sqrt(d) divides x-i, 0<= i < d, i.e. if d divides u-i.  \n");
			printf("Trajectories are declared divergent if |u| or |v| ");
			printf(" becomes greater than 2^{16xINFINITY}\n");
			printf("and where the trajectories start from spiral(n), ");
			printf("0 <= n <= RANGE.\n");
			printf("Output is sent to alg_cycle.out\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, spiral_str) == 0)
		{
			printf("          **************************\n");
			printf("spiral: usage: spiral(n,&x,&y)\n");
			printf("If n>=0, (x,y) is the point on the spiral\n");
			printf("defined on page 99 of Concrete Mathematics by\n");
			printf("Graham, Knuth and Pataschnik.");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, spiralinverse_str) == 0)
		{
			printf("          **************************\n");
			printf("spiral: usage: n=spiralinverse(x,y)\n");
			printf("n gives rise to the point (x,y) on the spiral\n");
			printf("defined on page 99 of Concrete Mathematics by\n");
			printf("Graham, Knuth and Pataschnik.");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, nscf_period0_str) == 0)
		{
			printf("          **************************\n");
			printf("nscfperiod0: usage: z=nscfperiod0(d,e,&x,&y)\n");
			printf("The period-length z of the nearest square\n");
			printf("continued fraction expansion of sqrt{d}\n");
			printf("is returned, together with the fundamental solution of the Pell equation");
			printf("x^2-dy^2=+/-1.\n");
			printf("z=nscfperiod0(d,1,&x,&y) prints the type of midpoint criteria met.\n");
			printf("z=nscfperiod0(d,0,&x,&y) prints only the period.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, rcf_period0_str) == 0)
		{
			printf("          **************************\n");
			printf("rcfperiod0: usage: z=rcfperiod0(d,e,&x,&y)\n");
			printf("The period-length z of the regular\n");
			printf("continued fraction expansion of sqrt{d}\n");
			printf("is returned, together with the fundamental solution of the Pell equation");
			printf("x^2-dy^2=+/-1.\n");
			printf("z=rcfperiod0(d,1,&x,&y) prints the type of midpoint criteria met.\n");
			printf("z=rcfperiod0(d,0,&x,&y) prints only the period.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, nicf_period0_str) == 0)
		{
			printf("          **************************\n");
			printf("nicfperiod0: usage: z=nicfperiod0(d,e,&x,&y)\n");
			printf("The period-length z of the nearest integer\n");
			printf("continued fraction expansion of sqrt{d}\n");
			printf("is returned, together with the fundamental solution of the Pell equation");
			printf("x^2-dy^2=+/-1.\n");
			printf("z=nicfperiod0(d,1,&x,&y) prints the type of midpoint criteria met.\n");
			printf("z=nicfperiod0(d,0,&x,&y) prints only the period.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, carmichael_str) == 0)
		{
			printf("          **************************\n");
			printf("carmichael: usage: carmichael(n)\n");
			printf("This finds the solutions, if any of phi(x)=n.\n");
			printf("If the number of primes p such that p-1 divides n is >100,\n");
			printf("we stop when two solutions are produced.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, carnielli_str) == 0)
		{
			printf("          **************************\n");
			printf("carnielli: usage: carnielli()\n");
			printf("This looks for cycles for Walter Carnielli's generalization of the 3x+1 mapping.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, lupei_str) == 0)
		{
			printf("          **************************\n");
			printf("lupei: usage: lupei()\n");
			printf("This looks for cycles for Lu Pei's generalization of the 3x-1 mapping.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, tangent_str) == 0)
		{
			printf("          **************************\n");
			printf("tangent: usage: tangent(n)\n");
			printf("This calculates the tangent function t(n) for n <= 2000.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, bernoulli_str) == 0)
		{
			printf("          **************************\n");
			printf("bernoulli: usage: bernoulli(n,&x,&y)\n");
			printf("This calculates the n-th Bernoulli number x/y for n <= 4000.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, partition_str) == 0)
		{
			printf("          **************************\n");
			printf("partition: usage: partition(n)\n");
			printf("This calculates the n-th partition number p(n) for n <= 65535.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, twocycle_str) == 0)
		{
			printf("          **************************\n");
			printf("twocycle: usage twocycle()\n");
			printf("Tries to find all cycles for the 2-branched generalised Collatz mapping T,\n");
			printf("which arise from starting numbers p, |p|<= RANGE/2.\n");
			printf("Trajectories containing an iterate whose magnitude\n");
			printf("exceeds a prescribed value INFINITY are deemed non-cycling.\n");
			printf("T(x)=int(ai*x/mi)+xi if x=i (mod 2).\n");
			printf("The 2 nonzero divisors mi, multipliers ai and shifts xi are entered from the keyboard.\n");
			printf("(See http://www.numbertheory.org/3x+1/)\n");
			printf("Output is sent to two_cycle.out\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, mcycle_str) == 0)
		{
			printf("          **************************\n");
			printf("mcycle: usage mcycle( )\n");
			printf("Tries to find all cycles for the 2-branched m-Collatz mapping of Benoit Cloitre,\n");
			printf("which arise from starting numbers p, |p|<= RANGE/2.\n");
			printf("Trajectories containing an iterate whose magnitude\n");
			printf("exceeds a prescribed value INFINITY are deemed non-cycling.\n");
			printf("f(x)=int((m+1)*x/m).\n");
			printf("T(x)=f(x)/2 if f(x) is even,\n");
			printf("T(x)=(3*f(x)+1)/2 if f(x) is odd.\n");
			printf("Output is sent to m_cycle.out\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, pqcycle_str) == 0)
		{
			printf("          **************************\n");
			printf("pqcycle: usage pqcycle( )\n");
			printf("Tries to find all cycles for the pq-Collatz mapping of Benoit Cloitre,\n");
			printf("which arise from starting numbers pp, |pp|<= RANGE/2.\n");
			printf("Trajectories containing an iterate whose magnitude\n");
			printf("exceeds a prescribed value INFINITY are deemed non-cycling.\n");
			printf("f(x)=int((p*x/q), p > q > 1, p not a multiple of q.\n");
			printf("T(x)=f(x)/2 if f(x) is even,\n");
			printf("T(x)=(3*f(x)+1)/2 if f(x) is odd.\n");
			printf("Output is sent to m_cycle.out\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, mmcycle_str) == 0)
		{
			printf("          **************************\n");
			printf("mcycle: usage mmcycle( )\n");
			printf("Tries to find all cycles for the 2-branched m-Collatz mapping of Benoit Cloitre,\n");
			printf("m >= 7, which arise from starting numbers p, |p|<= RANGE/2.\n");
			printf("f(x)=int((m+1)*x/m).\n");
			printf("T(x)=f(x)/2 if f(x) is even,\n");
			printf("T(x)=(3*f(x)+1)/2 if f(x) is odd.\n");
			printf("Output is sent to mm_cycle.out and is for eventual latex tabular output\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, exceptionals) == 0)
		{
			printf("          **************************\n");
			printf("exceptionals: usage exceptionals(c)\n");
			printf("Finds all exceptional solutions (k,x,y) with k<= c\n");
			printf("Output is sent to exceptionals.out\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, exceptionals0) == 0)
		{
			printf("          **************************\n");
			printf("exceptionals0: usage exceptionals0(c)\n");
			printf("Finds all exceptional solutions (k,x,y) with k<= c\n");
			printf("Output is sent only to exceptionals0.out and not to the screen\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, fg) == 0)
		{
			printf("          **************************\n");
			printf("fg: usage h=fg(p,q), p and q polynomials\n");
			printf("computes h=p(q(X))\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, exceptionals10) == 0)
		{
			printf("          **************************\n");
			printf("exceptionals10: usage exceptionals10(c)\n");
			printf("Finds all Type 1 exceptional solutions (k,x,y) with k<= c\n");
			printf("Output is sent only to exceptionals10.out\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, gplus) == 0)
		{
			printf("          **************************\n");
			printf("gplus: usage gplus(k,x,y)\n");
			printf("This is the function g+ of paper dujella.pdf\n");
			printf(" (k1,x1,y1)=g+(k,x,y).\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, gzero) == 0)
		{
			printf("          **************************\n");
			printf("gzero: usage gzero(k,x,y)\n");
			printf("This is the function g0 of paper dujella.pdf\n");
			printf(" (k1,x1,y1)=g0(k,x,y).\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, gminus) == 0)
		{
			printf("          **************************\n");
			printf("gminus: usage gminus(k,x,y)\n");
			printf("This is the function g- of paper dujella.pdf\n");
			printf(" (k1,x1,y1)=g-(k,x,y).\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, nagell) == 0)
		{
			printf("          **************************\n");
			printf("nagell: usage nagell(d,n)\n");
			printf("This solves x^2-dy^2=n using nagells' method\n");
			printf("and finds the fundamental solutions\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, frattini) == 0)
		{
			printf("          **************************\n");
			printf("frattini: usage frattini(d,n)\n");
			printf("This solves x^2-dy^2=n using frattini's method\n");
			printf("and finds a sequence of non-negative solutions u+v\\sqrt(D) which generate the non-negative solutions x+y\\sqrt(D)=(u+v\\sqrt(D))epsilon^n, n>=0\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, kashihara) == 0)
		{
			printf("          **************************\n");
			printf("kashihara: usage kashihara(m,n), where m<=n\n");
			printf("This solves x^2-(z^2-1)y^2=2-z^2 using nagell's method\n");
			printf(" and looks for z in the range m<=z<=n for which there are\n");
			printf(" >= 6 solution classes.  So far for m=2, n=10^7, only one \n");
			printf(" found is z=33539, which was known to Kashihara\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, branch) == 0)
		{
			printf("          **************************\n");
			printf("branch: usage branch(n,TYPE, FLAG).\n");
			printf("This constructs the polynomial triple (k(t),x(t),y(t))\n");
			printf("corresponding to a path given by the reverse digits\n");
			printf("of the negative 3-adic expansion of N, where the root node\n");
			printf("is (t,t,0) if TYPE=0, (t,t^2-t+1,t-1) if TYPE=-1,\n");
			printf("(t,t^2+t+1,t+1) if TYPE=1. If FLAG is nonzero, we print\n");
			printf("(k(t),x(t),y(t)), otherwise only k(t).\n");
			printf("Dujella's unicity conjecture is equivalent to the statement\n");
			printf("that the values of k(t) as t ranges over t >= 2 for TYPE =0,\n");
			printf("t >= 1 for TYPE 1 and t>=2 for TYPE -1, are distinct.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, ternary) == 0)
		{
			printf("          **************************\n");
			printf("ternary: usage z=ternary(n,&a[]).\n");
			printf("This contructs the ternary expansion of a positive integer,\n");
			printf("with digits a[i] from 0, 1 or -1.\n");
			printf("Returns the number of digits.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, pell4) == 0)
		{
			printf("          **************************\n");
			printf("pell4: usage pell4(d,e,&x1,&y1).\n");
			printf("This finds the smallest positive integer solution(x1,y1) of,\n");
			printf("the equation x^2-dy^2=4, where d > 1 is not a square.\n");
			printf("We use the fact if d > 16, a positive solution with gcd(x,y)=1 is a\n");
			printf("convergent of the continued fraction expansion of sqrt(d).\n");
			printf("If there is no such convergent, we use the midpoint method\n");
			printf("for finding the least positive solution (X1,Y1) of X^2-dY^2=1 \n");
			printf("and then (x1,y1)=(2X1,2Y1).\n");
			printf("e=1 displays some partial quotients.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, lprimefactor) == 0)
		{
			printf("          **************************\n");
			printf("lprimefactor: usage z=lprimefactor(n).\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, minarray) == 0)
		{
			printf("          **************************\n");
			printf("minarray: usage z=minarray(a[],n).\n");
			printf("This finds the minimum integer in a list a[0],...,a[n-1].\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, conwaycycles) == 0)
		{
			printf("          **************************\n");
			printf("minarray: usage z=conwaycycles(a,b).\n");
			printf("Here a and b are positive integers. \n");
			printf("We perform Conway's conjecture for his least prime factor sequence\n");
			printf("and believe that all sequences will eventually a cycle of length z.\n");
			printf(" I believe that apart from cycles of length 1\n");
			printf("there are just 6 cycles.\n");
			printf("          **************************\n");
			GetReturn();
		}
		else if (strcmp(name, conwayrangetest1) == 0)
		{
			printf("          **************************\n");
			printf("conwayrangetest1: usage z=conwayrangetest1(m1,m2,n1,n2).\n");
			printf("Here m1<=m2, n1<=n2 are positive integers. \n");
			printf("We test Conway's conjecture for his least prime factor sequence\n");
			printf("in the range [m1,m2], [n1,n2] and believe that all sequences will eventually cycle.\n");
			printf(" I believe that apart from cycles of length 1\n");
			printf("there are just 6 cycles (31st May 2016).\n");
			printf("          **************************\n");
			GetReturn();
		}
		else {
			continue;
		}
	}
	return;
}

