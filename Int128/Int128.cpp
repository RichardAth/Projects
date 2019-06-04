// Int128.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"
#include <tchar.h>
#include <limits.h>
#include <assert.h>
#include <time.h>

/* convert unsigned/signed number to ascii string */
std::string ToString(const _uint128 &n, int radix = 10) {
	static char buffer[50];
	if (radix == 16)
		sPrintf128(buffer, sizeof(buffer), "%#32x", n);	
	else
		sPrintf128(buffer, sizeof(buffer), "%20d", n);
	
	return std::string(buffer);
}
std::string ToString(const _int128 &n, int radix = 10) {

	static char buffer[50];
	if (radix == 16)
		sPrintf128(buffer, sizeof(buffer), "%#32x", n);
	else
		sPrintf128(buffer, sizeof(buffer), "%20d", n);

	return std::string(buffer);
}

/* convert bool to ascii string */
char * boolToStr(bool b) {
	static char t[] = "true",
		f[] = "false";

	if (b)
		return t;
	else
		return f;
}

#define RADIX 0
//#define RADIX 16

/* set value in signed int128 from keyboard input */
static _int128 inputInt128(TCHAR const * const msg) {
	char buffer[200];
	char * endptr;
	wprintf_s(L"%s", msg);
	std::cin >> buffer;
	_int128 result = _strtoi128(buffer, &endptr, RADIX);
	if (errno == ERANGE) {
		throw std::runtime_error("out of range");
	}

	return result;
}

/* set value in unsigned int128 from keyboard input */
static _uint128 inputUint128(TCHAR const * const msg) {
	char buffer[200];
	char * endptr;
	wprintf_s(L"%s", msg);
	std::cin >> buffer;
	_uint128 result = _strtoui128(buffer, &endptr, RADIX);
	if (errno == ERANGE) {
		throw std::runtime_error("out of range");;
	}
	if (strlen(endptr) != 0) {
		std::cout << "unprocessed characters: " << endptr << '\n';
	}

	return result;
}

template<class _itype> void loopItype(_itype(*inputItype)(const TCHAR *)) {
	for (;;) {
		const _itype x = inputItype(L"Enter x:");  // inputItype = inputUint128 or inputInt128
		const _itype y = inputItype(L"Enter y:");

		_itype s = x + y;
		_itype d = x - y;
		_itype p = x * y;
		_itype q = x / y;
		_itype r = x % y;

		/*std::cout << "x = " << x << '\n';
		std::cout << "y = " << y << '\n';*/
		_tprintf(L"x       = %hs\n", ToString(x).c_str());
		_tprintf(L"y       = %S\n", ToString(y).c_str());

		_tprintf(L"x+y = s = %hs\n", ToString(s).c_str());
		_tprintf(L"x-y = d = %hs\n", ToString(d).c_str());
		_tprintf(L"x*y = p = %hs\n", ToString(p).c_str());
		_tprintf(L"x/y = q = %hs\n", ToString(q).c_str());
		_tprintf(L"x%%y = r = %hs\n", ToString(r).c_str());

		_itype c1 = (s + d) / _int128(2);
		_itype c2 = (s - d) / _int128(2);
		_itype c3 = q * y + r;

		_tprintf(L"(s+d)/2 = %hs (check = x)\n", ToString(c1).c_str());  // (x+y) + (x-y)=2x; 2x/2 = x
		_tprintf(L"(s-d)/2 = %hs (check = y)\n", ToString(c2).c_str());  // (x+y) -(x-y) = 2y; 2y/2 = y;
		_tprintf(L"(q*y+r) = %hs (check = x)\n", ToString(c3).c_str());             // (x/y)*y + r = x
		_tprintf(L"______________________________________________\n");

		const bool lt = x < y;
		const bool le = x <= y;
		const bool gt = x > y;
		const bool ge = x >= y;

		_tprintf(L"x<y:%hs, x<=y:%hs, x>y:%hs, x>=y:%hs\n\n",
			boolToStr(lt), boolToStr(le), boolToStr(gt), boolToStr(ge));
		s = x + 10ULL;
		d = s-10;
		assert(d == x);
	}
}

long long dobenchmark(void) {
	std::string oper;
	int operval = -1;
#ifdef _DEBUG
	const int tries = 9'000'000;   // compensate partly for debug being slower
#else
	const int tries = 90'000'000;
#endif
	_uint128 x, y, a, b;
	unsigned long long xll, yll, zll, all, bll;
	double ratio;
	clock_t starttime, midtime, endtime;
	long long ctr = 0;

	srand(474644141);

	do {
		std::cout << "select operator (e to exit) \n";
		std::cin >> oper;
		if (oper == "e" || oper == "E")
			break;

		operval = -1;
		if (oper == "*")
			operval = 0;
		else if (oper == "/")
			operval = 1;
		else if (oper == "%")
			operval = 2;
		else if (oper == "+")
			operval = 3;
		else if (oper == "-")
			operval = 4;
		else if (oper == "<<")
			operval = 5;
		else if (oper == ">>")
			operval = 6;
		else if (oper == "<")
			operval = 7;
		else if (oper == ">")
			operval = 8;

#ifdef _DEBUG
		std::cout << "operval = " << operval << '\n';
#endif

		switch (operval) {
		case 0:    // multiply
			xll = rand();
			yll = rand();
			x = xll;
			y = yll;
			starttime = clock();
			for (ctr = 1; ctr <=tries; ctr++) {
				a = x * y;         // use _uint128
				b = a * x;
				a = b * y;
			}
			midtime = clock();
			std::cout << "_uint128 part completed \n";
			for (ctr = 1; ctr <=tries*10; ctr++) {
				all = xll * yll;   // use long long do 10x number of operations
				bll = all * xll;
				all = bll * yll;
				zll += all;
			}
			endtime = clock();
			/* multiply time for _uint128 by 10 to compensate for fewer cycles*/
			ratio = double(midtime - starttime)*10 / double(endtime - midtime);
			std::cout << "multiply ratio = " << ratio;
			std::cout << " elapsed time = " << (double)(endtime - starttime) / CLOCKS_PER_SEC << 
				" seconds \n";
			if (a != all)
				std::cout << "error in calc x = " << x << '\n';;
			break;

		case 1:    // divide
			xll = rand()*rand()*rand();
			yll = rand();
			x = xll;
			y = yll;
			starttime = clock();
			for (ctr = 1; ctr <=tries; ctr++) {
				a = x / y;         // use _uint128
				b = a / y;
				a = b / y;
			}
			midtime = clock();
			std::cout << "_uint128 part completed \n";
			for (ctr = 1; ctr <=tries*10; ctr++) {
				all = xll / yll;   // use long long do 10x number of operations
				bll = all / xll;
				all = bll / yll;
				zll += all;
			} 
			endtime = clock();
			/* multiply time for _uint128 by 10 to compensate for fewer cycles*/
			ratio = double(midtime - starttime) * 10 / double(endtime - midtime);
			std::cout << "divide ratio = " << ratio;
			std::cout << " elapsed time = " << (double)(endtime - starttime) / CLOCKS_PER_SEC <<
				" seconds \n";
			if (a != all)
				std::cout << "error in calc x = " << x << '\n';
			break;

		case 2:		// modulus
			xll = rand()*rand()*rand();
			yll = rand();
			x = xll;
			y = yll;
			starttime = clock();
			for (ctr = 1; ctr <=tries; ctr++) {
				a = x % y;         // use _uint128
				b = a % y;
				a = b % y;
			}
			midtime = clock();
			std::cout << "_uint128 part completed \n";
			for (ctr = 1; ctr <=tries*10; ctr++) {
				all = xll % yll;   // use long long do 10x number of operations
				bll = all % xll;
				all = bll % yll;
				zll += all;
			}
			endtime = clock();
			/* multiply time for _uint128 by 10 to compensate for fewer cycles*/
			ratio = double(midtime - starttime) * 10 / double(endtime - midtime);
			std::cout << "modulus ratio = " << ratio;
			std::cout << " elapsed time = " << (double)(endtime - starttime) / CLOCKS_PER_SEC <<
				" seconds \n";
			if (a != all)
				std::cout << "error in calc x = " << x << '\n';
			break;

		case 3:	   // add
			xll = rand();
			yll = rand();
			x = xll;
			y = yll;
			starttime = clock();
			for (ctr = 1; ctr <=tries; ctr++) {
				a = x + y;         // use _uint128
				b = a + x;
				a = b + y;
			}
			midtime = clock();
			std::cout << "_uint128 part completed \n";
			for (ctr = 1; ctr <=tries*10; ctr++) {
				all = xll + yll;   // use long long do 10x number of operations
				bll = all + xll;
				all = bll + yll;
				zll += all;
			}
			endtime = clock();
			/* multiply time for _uint128 by 10 to compensate for fewer cycles*/
			ratio = double(midtime - starttime) * 10 / double(endtime - midtime);
			std::cout << "add ratio = " << ratio;
			std::cout << " elapsed time = " << (double)(endtime - starttime) / CLOCKS_PER_SEC <<
				" seconds \n";
			if (a != all)
				std::cout << "error in calc x = " << x << '\n';;
			break;

		case 4:    // subtract
			xll = rand();
			yll = rand();
			x = xll;
			y = yll;
			starttime = clock();
			for (ctr = 1; ctr <=tries; ctr++) {
				a = x - y;         // use _uint128
				b = a - x;
				a = b - y;
			}
			midtime = clock();
			std::cout << "_uint128 part completed \n";
			for (ctr = 1; ctr <=tries*10; ctr++) {
				all = xll - yll;   // use long long do 10x number of operations
				bll = all - xll;
				all = bll - yll;
				zll += all;
			}
			endtime = clock();
			/* multiply time for _uint128 by 10 to compensate for fewer cycles*/
			ratio = double(midtime - starttime) * 10 / double(endtime - midtime);
			std::cout << "subtract ratio = " << ratio;
			std::cout << " elapsed time = " << (double)(endtime - starttime) / CLOCKS_PER_SEC <<
				" seconds \n";
			if (a != all)
				std::cout << "error in calc x = " << x << '\n';;
			break;

		case 5:    // left shift
			xll = rand();
			x = xll;
			starttime = clock();
			for (ctr = 1; ctr <=tries; ctr++) {
				a = x << 4;         // use _uint128
				b = a << 5;
				a = b << 6;
			}
			midtime = clock();
			std::cout << "_uint128 part completed \n";
			for (ctr = 1; ctr <=tries*10; ctr++) {
				all = xll << 4;   // use long long do 10x number of operations
				bll = all << 5;
				all = bll << 6;
				zll += all;
			}
			endtime = clock();
			/* multiply time for _uint128 by 10 to compensate for fewer cycles*/
			ratio = double(midtime - starttime) * 10 / double(endtime - midtime);
			std::cout << "left shift ratio = " << ratio;
			std::cout << " elapsed time = " << (double)(endtime - starttime) / CLOCKS_PER_SEC <<
				" seconds \n";
			if (a != all)
				std::cout << "error in calc x = " << x << '\n';;
			break;

		case 6:    // right shift 
			xll = rand();
			x = xll;
			starttime = clock();
			for (ctr = 1; ctr <= tries; ctr++) {
				a = x >> 4;         // use _uint128
				b = a >> 5;
				a = b >> 6;
			}
			midtime = clock();
			std::cout << "_uint128 part completed \n";
			for (ctr = 1; ctr <= tries * 10; ctr++) {
				all = xll >> 4;   // use long long do 10x number of operations
				bll = all >> 5;
				all = bll >> 6;
				zll += all;
			}
			endtime = clock();
			/* multiply time for _uint128 by 10 to compensate for fewer cycles*/
			ratio = double(midtime - starttime) * 10 / double(endtime - midtime);
			std::cout << "right shift ratio = " << ratio;
			std::cout << " elapsed time = " << (double)(endtime - starttime) / CLOCKS_PER_SEC <<
				" seconds \n";
			if (a != all)
				std::cout << "error in calc x = " << x << '\n';;
			break;
		case 7:    // less than
			xll = rand();
			yll = rand();
			x = xll;
			y = yll;
			starttime = clock();
			for (ctr = 1; ctr <= tries; ctr++) {
				a = x < y;         // use _uint128
				b = a < x;
				a = b < y;
			}
			midtime = clock();
			std::cout << "_uint128 part completed \n";
			for (ctr = 1; ctr <= tries * 10; ctr++) {
				all = xll < yll;   // use long long do 10x number of operations
				bll = all < xll;
				all = bll < yll;
				zll += all;
			}
			endtime = clock();
			/* multiply time for _uint128 by 10 to compensate for fewer cycles*/
			ratio = double(midtime - starttime) * 10 / double(endtime - midtime);
			std::cout << "'less than' ratio = " << ratio;
			std::cout << " elapsed time = " << (double)(endtime - starttime) / CLOCKS_PER_SEC <<
				" seconds \n";
			if (a != all)
				std::cout << "error in calc x = " << x << '\n';;
			break;

		case 8:    // greater than
			xll = rand();
			yll = rand();
			x = xll;
			y = yll;
			starttime = clock();
			for (ctr = 1; ctr <= tries; ctr++) {
				a = x > y;         // use _uint128
				b = a > x;
				a = b > y;
			}
			midtime = clock();
			std::cout << "_uint128 part completed \n";
			for (ctr = 1; ctr <= tries * 10; ctr++) {
				all = xll > yll;   // use long long do 10x number of operations
				bll = all > xll;
				all = bll > yll;
				zll += all;
			}
			endtime = clock();
			/* multiply time for _uint128 by 10 to compensate for fewer cycles*/
			ratio = double(midtime - starttime) * 10 / double(endtime - midtime);
			std::cout << "'greater than' ratio = " << ratio;
			std::cout << " elapsed time = " << (double)(endtime - starttime) / CLOCKS_PER_SEC <<
				" seconds \n";
			if (a != all)
				std::cout << "error in calc x = " << x << '\n';;
			break;

		default:
			std::cout << "invalid operator \n";
		}
	} while (oper != "e" && oper != "E");

	return zll;
}

int main(int argc, TCHAR **argv) {
	try {
		char yn;
		do {
			std::cout << "Do benchmark tests? answer Y or N \n";
			std::cin >> yn;
			yn = toupper(yn);
		} while (yn != 'Y' && yn != 'N');

		if (yn == 'Y') {
			dobenchmark();
			return 0;
		}

		for (double x = -10.0; x <= 10.0; x++) {
			_int128 x128;
			x128 = doubleTo128(x);
			std::cout << "x = " << x << ", x128 = " << x128 << '\n';
		}
		for (int i = -10; i <= 100; i += 2) {
			double x = pow(2, i);
			_int128 x128;
			x128 = doubleTo128(x);
			std::cout << "x = " << x << ", x128 = " << x128 << '\n';
		}

		_int128 x128 = -1;
		char buffer[50];
		/* note that % character is optional in format string. It's allowed
		for compatibility with standard printf */
		sPrintf128(buffer, sizeof(buffer), "#x", x128);
		printf_s("%s \n", buffer);
		sPrintf128(buffer, sizeof(buffer), "+d", x128);
		printf_s("%s \n", buffer);

		long long a = 12345678901, b = 45678901234;
		ui128mult(x128, a, b);   // multiply 64 bit x 64 bit to get 128 bits.
		//assert(x128 == 5639370470903950470ULL);
		assert(LO64(x128) == 10534724974170115354ULL);
		assert(HI64(x128) == 30);

		x128 = 0;
		i128mult(x128, a, b);   // multiply 64 bit x 64 bit to get 128 bits.
		assert(LO64(x128) == 10534724974170115354ULL);
		assert(HI64(x128) == 30);

		x128 = 0;
		b = -1;
		i128mult(x128, a, b);   // multiply 64 bit x 64 bit to get 128 bits.
		assert(x128 == -12345678901);


		loopItype<_uint128>(inputUint128);  /* test unsigned 128 bit integers */
	}
	catch (const std::exception& e) {
		std::cout << "exception occurred: " << e.what() << '\n';
	}

	try {
		loopItype<_int128 >(inputInt128);   /* test signed 128 bit integers */
	}
	catch (const std::exception& e) {
		std::cout << "exception occurred: " << e.what() << '\n';
	}

	try {
		_uint128 x6, y, z;
		_int128 z2, a, b, r;

		/* test gcd function */
		do {
			x6 = inputUint128(L"enter value for x:");
			y = inputUint128(L"enter value for y:");
			z = gcd(x6, y);
			z2 = extendedGcd(x6, y, a, b);
			/* should have x6*a + y*b = gcd */
			r = x6 * a + y * b;

			if (r != z2)
				std::cout << "**error. Expected " << z2 << " got " << r << '\n';

			std::cout << "gcd(" << x6 << ',' << y << ") = " << z << '\n';
			std::cout << "extended gcd. z2 = " << z2 << " a = " << a
				<< " b = " << b << '\n';
		} while (x6 != 0);

		_uint128  x1(0xaaaaaaaabbbbbbbb, 0xfabcdef12345678);
		_uint128  x2 = 0xccccdddddddd;

		_uint128 x3 = x1 / x2;
		_uint128 x4 = x1 % x2;
		_uint128 x5 = x3 * x2 + x4;  // x1/x2*x2+x4 = x1

		_tprintf(L"x3.%016I64x:%016I64x (%016I64x%016I64x)\n",
			HI64(x3), LO64(x3), HI64(x3), LO64(x3));
		_tprintf(L"x4.%016I64x:%016I64x (%016I64x%016I64x)\n",
			HI64(x4), LO64(x4), HI64(x4), LO64(x4));
		_tprintf(L"x5.%016I64x:%016I64x (%016I64x%016I64x)\n",
			HI64(x5), LO64(x5), HI64(x5), LO64(x5));


		_int128 x(0x1234567890, 0x1234567890123456);
		y = 0x43215456543;

		_tprintf(L"maxInt:%hs (= %hs)\n", ToString(x1, 10).c_str(), ToString(x1, 16).c_str());
		_tprintf(L"minInt:%hs\n", ToString(x2, 10).c_str());
		_tprintf(L"maxUInt:%hs (= %hs)\n", ToString(x3, 10).c_str(), ToString(x3, 16).c_str());
		//return 0;

		for (int i = 0; i < 129; i++) {
			_uint128 z = x >> i;
			_tprintf(L"i:%3d, x>>%3d:%32hs ", i, i, ToString(z, 16).c_str());
			_tprintf(L" (1 bits =  %d) \n", z.popcnt());
		}
		for (int i = 0; i < 129; i++) {
			_uint128 z = x << i;
			_tprintf(L"i:%3d, x<<%3d:%32hs  ", i, i, ToString(z, 16).c_str());
			_tprintf(L" (msb =  %d) ", z.msb());
			_tprintf(L" (lsb =  %d) \n", z.lsb());
		}
		//return 0;


		_int128 s = x + y;
		_int128 d = x - y;
		_int128 p = x * y;
		_int128 q = x / y;
		r = x % y;

		_tprintf(L"%-20hs + %-20hs = %20hs\n", ToString(x, RADIX).c_str(),
			ToString(y, RADIX).c_str(), ToString(s).c_str());
		_tprintf(L"%-20hs - %-20hs = %20hs\n", ToString(x, RADIX).c_str(),
			ToString(y, RADIX).c_str(), ToString(d).c_str());
		_tprintf(L"%-20hs * %-20hs = %20hs\n", ToString(x, RADIX).c_str(),
			ToString(y, RADIX).c_str(), ToString(p).c_str());
		_tprintf(L"%-20hs / %-20hs = %20hs\n", ToString(x, RADIX).c_str(),
			ToString(y, RADIX).c_str(), ToString(q, RADIX).c_str());
		_tprintf(L"%-20hs %% %-20hs = %20hs\n", ToString(x, RADIX).c_str(),
			ToString(y, RADIX).c_str(), ToString(r, RADIX).c_str());

		x <<= 40;
		x += 123456789012;
		printf("x = %s \n", ToString(x, RADIX).c_str());
		for (int i = 2; i <= 126; i++) {
			auto root = nthroot(x, i);
			printf("%dth root = %s \n", i, ToString(root, RADIX).c_str());
		}
		r = modMultInv(x, y);
		std::cout << "modular inverse of " << x << " wrt " << y << " is " << r << '\n';

		a = -987654;
		b = 12345;
		auto dr = divrem(a, b);
		q = dr.quot;
		r = dr.rem;
		std::cout << "a, b, q, r = " << a << ", " << b << ", " << q << ", " << r << '\n';
		return 0;
	}
	catch (const std::exception& e) {
		std::cout << "exception occurred: " << e.what() << '\n';
		return -1;
	}
}



// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
