This a test program based on Dario Alpert's code. It is capable of computing integer
expressions incorporating numerous additional operators and functions such as
#(primorial) !(factorial) and !!(double factorial)

B(n)                                Previous probable prime before n
F(n)                                Fibonacci number Fn
L(n)                                Lucas number Ln = Fn-1 + Fn+1
N(n)                                Next probable prime after n
P(n)                                Unrestricted Partition Number (number of 
                                    decompositions of n into sums of integers without 
                                    regard to order).
Gcd(m,n)                            Greatest common divisor of m and n.
Modinv(m,n)                         inverse of m modulo n, only valid when gcd(m,n)=1.
Modpow(m,n,r)                       finds m^n modulo r. more efficient than (m^n)%r
Totient(n)                          finds the number of positive integers less than n 
                                    which are relatively prime to n.
IsPrime(n)                          returns zero if n is not probable prime, -1 if it is.
NumDivs(n)                          Number of positive divisors of n either prime or composite.
SumDivs(n)                          Sum of all positive divisors of n both prime and composite.
NumDigits(n,r)                      Number of digits of n in base r.
SumDigits(n,r)                      Sum of digits of n in base r.
RevDigits(n,r)                      finds the value obtained by writing backwards 
                                    the digits of n in base r. 
                                    
Large numbers (100 digits or more) can be factorised
