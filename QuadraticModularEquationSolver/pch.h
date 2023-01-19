// pch.h: This is a precompiled header file.
// Files listed below are compiled only once, improving build performance for future builds.
// This also affects IntelliSense performance, including code completion and many code browsing features.
// However, files listed here are ALL re-compiled if any one of them is updated between builds.
// Do not add files here that you will be updating frequently as this negates the performance advantage.

#ifndef PCH_H
#define PCH_H

// add headers that you want to pre-compile here
#include "framework.h"
#pragma once
#include <iostream>
#define __STDC_WANT_LIB_EXT1__ 1  /* ask for printf_s etc */
#include <cstdio>
#include <string>
#include <vector>
#include <cctype>
#include <ctime>
#include <stdexcept>
#include <iostream>
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

#include "boost/multiprecision/gmp.hpp" 
typedef boost::multiprecision::mpz_int Znum;

#define ENABLE_INTSAFE_SIGNED_FUNCTIONS 1
#include <intsafe.h>
#include <intrin.h>

#endif //PCH_H
