#pragma once

#include "targetver.h"

#include <iostream>
#define __STDC_WANT_LIB_EXT1__ 1  /* ask for printf_s etc */
#include <cstdio>
#include <cstring>
#include <sys/types.h>
#include <sys/stat.h>
#include <vector>
#include <cctype>
#include <ctime>
#include <stdexcept>

#include <sstream>
#include <limits>
#include <locale>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <random>
#include <algorithm>
#ifdef __GNUC__
#include "gmp.h"
#else
#include "mpir.h"
#endif

#include <Windows.h>

#include "boost/multiprecision/gmp.hpp" 
typedef boost::multiprecision::mpz_int Znum;
#define ZT(a) a.backend().data()  /* access mpz_t within a Znum (Boost mpz_int)*/

/* following 5 includes must be in this order; */
#include "showtime.h"
#include "diagnostic.h"
#include "bignbr.h"
#include "bigint.h"
#include "factor.h"
#include "main.h"

#define ENABLE_INTSAFE_SIGNED_FUNCTIONS 1
#include <intsafe.h>
#include <intrin.h>
