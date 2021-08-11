Freeware copyright (c) 2009 Thomas R. Nicely <http://www.trnicely.net>.
 see https://web.archive.org/web/20131208185421/http://trnicely.net/misc/bpsw1.zip
 Released into the public domain by the author, who disclaims any legal
 liability arising from its use.

 Dr Nicely died Wednesday, September 11, 2019, 
 see https://www.lynchburg.edu/news/2019/09/remembering-dr-thomas-nicely-legendary-teacher-who-discovered-pentium-bug/
 
 NOTE: Most of the functions called in this code are defined and
 implemented in the separately downloadable support files trn.h and
 trn.c, available at <http://www.trnicely.net/misc/trn.zip>.
 N.B. Now need to use 
       https://web.archive.org/web/20161030053956/http://trnicely.net/misc/trn.zip
 Command-line compilation of the code is carried out by a command
 such as
 gcc bpsw1.c trn.c -std=gnu99 -D__NOMPFR__ -lm -lgmp  -obpsw1.exe
 
 with trn.c and trn.h present in the current directory or on
 the search path. This will be dependent on the environment and
 configuration of your system. Note that recompilation does require
 GMP, but does _not_ require MPFR.
 
 SYNTAX: bpsw1 LB UB|dN [UF]
 
 The bounds may be expressed in floating point notation, e.g.,
 bpsw1 1e50 1e5. If arg2 is less than arg1, it is interpreted as
 an increment, and UB = arg1 + arg2. If arg2 is negative, the
 bounds are (arg1 + arg2) and arg1. The optional third argument
 UF sets the screen update frequency (default is every 10000
 integers).
 
 The purpose of this code is not simply to compute pi(x),
 which can be done much more efficiently, but to test and
 illustrate the various primality testing routines called.
 
 This code, and its support routines trn.c and trn.h, implement
 the standard and strong versions of the Lucas-Selfridge and
 Baillie-PSW primality tests, as well as the extra strong Lucas
 test; see <http://www.trnicely.net/misc/bpsw.html> for details.
 The GMP mpz_probab_prime_p function, employing Miller-Rabin
 tests with 13 different bases, is called to determine the
 "true" primality of each odd number between the specified
 bounds, and this result is then compared with those from the
 standard BPSW test, strong BPSW test, and extra strong Lucas
 test (base 3), as well as the standard and strong Lucas tests
 and the Miller-Rabin test with base 2. Discrepancies are
 reported as counts of psueudoprimes for each category.
 As of this date, there is no known integer N for which this
 implementation of either the standard or strong BPSW test
 will return a false primality result---there is no known
 strong or standard BPSW pseudoprime. The author has
 directly verified that none exists for N < 10^13; Martin Fuller
 has verified that none exists for N < 10^15. If a BPSW
 pseudoprime is detected (a significant and highly unexpected
 event), it will be reported directly to the screen.

 Note: according to 
 https://web.archive.org/web/20191121062007/http://www.trnicely.net/misc/bpsw.html
 there is no pseudoprime below 2^64 (approximately 1.845*10^19)
 
 See the documentation in the individual modules for additional
 details, including bibliographies. These modules are part of
 the accompanying support files trn.h and trn.c.