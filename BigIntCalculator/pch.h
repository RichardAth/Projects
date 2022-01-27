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
#include "main.h"

#define ENABLE_INTSAFE_SIGNED_FUNCTIONS 1
#include <intsafe.h>
#include <intrin.h>