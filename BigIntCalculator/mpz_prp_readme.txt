Copyright 2011-2015 David Cleaver

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.


v1.0 Posted to SourceForge on 2013/07/04

v1.1 Posted to SourceForge on 2013/12/27
[The following fix was recommended by Dana Jacobsen and verified by Jon Grantham]
     - Bug fix: Removed unnecessary vl==0 check in mpz_extrastronglucas_prp
[The following improvements/fixes were recommended by Laurent Desnogues in 2013/08]
     - Speed improvement 1: Removed extraneous NormalizeJS calls in ARPCL
     - Speed improvement 2: Removed/consolidated calls to mpz_mod in APRCL
       (these improvements make the APRCL code about 1.5-2.2x faster)
     - Bug fix: Final test in APRCL routine is now correct

v1.2 Posted to SourceForge on 2015/03/07
  - Minor change to code to remove "warning: array subscript is above array bounds"
    encountered while compiling with the options ( -O3 -Wall )


The PRP functions presented here are based on the paper:
Grantham, Jon. Frobenius Pseudoprimes. Math. Comp. 70 (2001), 873-891.


**********************************************************************************
APR-CL (also known as APRT-CLE) is a prime proving algorithm developed by:
L. Adleman, C. Pomerance, R. Rumely, H. Cohen, and H. W. Lenstra
APRT-CLE = Adleman-Pomerance-Rumely Test Cohen-Lenstra Extended version
You can find all the details of this implementation in the Cohen & Lenstra paper:
   H. Cohen and A. K. Lenstra, "Implementation of a new primality test",
   Math. Comp., 48 (1987) 103--121

----------------------------------------------------------------------------------

This C/GMP version is a conversion of Dario Alpern's Java based APRT-CLE code
His code was based on Yuji Kida's UBASIC code

Based on APRT-CLE Written by Dario Alejandro Alpern (Buenos Aires - Argentina)
From: Last updated September 10th, 2011. See http://www.alpertron.com.ar/ECM.HTM

On 2012/11/12 Dario Alpern has approved the conversion, from Java to C/GMP, of 
his implementation of the APR-CL algorithm, and that it be licensed under the LGPL.

----------------------------------------------------------------------------------

With improvements based on Jason Moxham's APRCL v1.15 code, from 2003/01/01

On 2013/04/14 Toby Moxham has approved the APR-CL code and data tables, 
originally written by his brother Jason Moxham on 2003/01/01, to be released 
under the LGPL.
*********************************************************************************

There are two versions of the APR-CL algorithm that are available.

The first is called mpz_aprcl on SourceForge and should work fine on any system that
can compile C, GMP and use a 64-bit data type.  I used the 'long long' data type in
the code, but you can change that through the available typedef in the three files:
mpz_aprcl.c
mpz_aprcl.h
jacobi_sum.h

In those three files, change the following typdef's:
#ifndef HAVE_U64_T
#define HAVE_U64_T
typedef long long s64_t;
typedef unsigned long long u64_t;
#endif
to whatever 64-bit type is available on your system.
The 64-bit version can prove numbers up to 6021 decimal digits long.  Note, however, it
took me two weeks to verify a 6000 digit prime number.  For composite numbers it is quite
a bit faster.  At the end of this file I am including a list of times to prove various
numbers prime with APR-CL and probably-prime by BPSW.

The APR-CL method does not produce a certificate of primality.  It runs through a series
of computations and at the end will tell you if the number is prime or composite.

The second version of this project is called mpz_aprcl32 on SourceForge.  It should
work on systems that can compile C, GMP, and 32-bit code but that do not have a 64-bit type.
The limitation of the 32-bit version is that it can "only" prove the primality of numbers
up to 3827 decimal digits.

Above those limits, the code will fall back to a BPSW test which will tell you if the number
is probably prime, or if it is composite.  Each APR-CL project has the full set of mpz_prp
functions included with it so that you can use probable prime tests or the full APR-CL prime
test by including just this one project into your own.


Time, in seconds, to run the APR-CL test on each given number
on one core of a dual Xeon 5335 2.0GHz, 12GB FB-DDR2 PC2-5300, computer.
(Not necessarily best observed time, just an indication of relative performance)

            :        Number           :   BPSW   :  APRCL v1.0 :  APRCL v1.1 :
  10 digits : 2^32+15                 :  0.0000s :     0.0000s :     0.0000s :
  20 digits : 2^64+13                 :  0.0000s :     0.0000s :     0.0000s :
  30 digits : (2^101+1)/3             :  0.0000s :     0.0156s :     0.0000s :
  40 digits : 4*10^39+7               :  0.0000s :     0.0156s :     0.0156s :
  50 digits : (2^167+1)/3             :  0.0000s :     0.0469s :     0.0313s :
  60 digits : (2^199+1)/3             :  0.0000s :     0.0625s :     0.0313s :
  70 digits : 10^69+9                 :  0.0000s :     0.0938s :     0.0469s :
  80 digits : 11^76*6+1               :  0.0000s :     0.1563s :     0.0781s :
  90 digits : 13^80-2                 :  0.0000s :     0.2031s :     0.0938s :
 100 digits : 5^143+2                 :  0.0000s :     0.3125s :     0.1563s :
 110 digits : 7^129+4                 :  0.0000s :     0.3750s :     0.2031s :
 120 digits : 13^107*3+2              :  0.0000s :     0.5469s :     0.2969s :
 130 digits : (10^129*89-17)/9        :  0.0000s :     0.7031s :     0.3438s :
 140 digits : 7^164*8+5               :  0.0000s :     0.9375s :     0.4531s :
 150 digits : 3^313*4+5               :  0.0000s :     1.1719s :     0.5469s :
 160 digits : (81^86-86^81)/20435     :  0.0000s :     1.4688s :     0.7031s :
 170 digits : 11^162*16+1             :  0.0156s :     1.7188s :     0.8438s :
 180 digits : 5^256*3+2               :  0.0000s :     2.1250s :     1.0000s :
 190 digits : 11^181*14+3             :  0.0156s :     2.4219s :     1.1563s :
 200 digits : 7^235*9+2               :  0.0000s :     3.0469s :     1.4844s :
 210 digits : 2^688*687-1             :  0.0000s :     3.6250s :     1.7031s :
 220 digits : (68^123+123^68)/335203  :  0.0156s :     4.1094s :     2.5000s :
 230 digits : (5^333-3)/2074          :  0.0156s :     4.6719s :     2.7656s :
 240 digits : 12^196*196^12+1         :  0.0000s :     5.2500s :     3.0938s :
 250 digits : 5^356*8-1               :  0.0156s :     6.4844s :     3.6563s :
 260 digits : 44^104*104^44+1         :  0.0000s :     6.9063s :     3.8750s :
 270 digits : (5^385*11+3)/2          :  0.0000s :     9.0781s :     5.2031s :
 280 digits : 3^583*17+2              :  0.0000s :    10.2813s :     5.6563s :
 290 digits : 5^413*18-1              :  0.0000s :    12.1563s :     6.3594s :
 300 digits : 424^114+3               :  0.0156s :    12.7344s :     6.7344s :
 310 digits : 3^648-2                 :  0.0156s :    16.2031s :     8.1719s :
 320 digits : 395^123+2               :  0.0156s :    17.8125s :     8.8750s :
 350 digits : 1160^114+7              :  0.0000s :    23.6563s :    12.6563s :
 400 digits : (291^163-1)/290         :  0.0156s :    39.3281s :    20.0938s :
 450 digits : 232^190+7               :  0.0313s :    58.0313s :    29.4375s :
 500 digits : 1014^166+7              :  0.0156s :    88.1406s :    50.7969s :
 550 digits : 10^549*9-7              :  0.0313s :   147.4688s :    83.1406s :
 600 digits : 1432^190+7              :  0.0625s :   210.3281s :   109.1094s :
 650 digits : 2^2159+375              :  0.0625s :   265.9844s :   137.3594s :
 700 digits : (157^319+319^157)/28    :  0.0781s :   403.0781s :   209.6406s :
 750 digits : 10^749*2+89             :  0.0938s :   519.9531s :   260.3125s :
 800 digits : (10^799*61-7)/9         :  0.1094s :   709.7500s :   358.6406s :
 850 digits : 2^2821-183              :  0.1406s :   974.4688s :   483.2500s :
 900 digits : (24^653-1)/23           :  0.1406s :  1063.8125s :   524.5625s :
 950 digits : 10^949*4-9              :  0.1719s :  1422.8125s :  1091.0156s :
1000 digits : 10^999+7                :  0.1719s :  1588.3125s :  1195.0313s :
1100 digits : 2^3653+41               :  0.2031s :  2237.8438s :  1574.7031s :
1200 digits : 10^1199*5+9             :  0.3281s :  3347.6094s :  2110.3281s :
1300 digits : 2^4318+165              :  0.3906s :  4324.2344s :  2610.4844s :
1400 digits : (187^617-1)/186         :  0.5000s :  6247.0469s :  3650.9844s :
1500 digits : 2^4972*1779-1           :  0.3281s :  8380.2500s :  4621.9844s :
1600 digits : (12^1483+1)/13          :  0.6406s : 10285.5469s :  5607.5469s :
1700 digits : 2^5644-227              :  0.8281s : 13942.1250s :  7221.0625s :
1800 digits : 10^1800-87              :  0.8750s : 16868.7188s :  8556.3281s :
1900 digits : (10^1900+3)/7           :  1.0938s : 20384.2813s :  9946.6719s :
2000 digits : (2^6643*113+1)/115      :  1.1094s : 24359.6719s : 18041.1094s :
