#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdio>
#include "factor.h"
const char * myTime(void);  // get time as hh:mm:ss

#ifndef _DEBUG
static std::string Path = "C:\\Users\\admin99\\Source\\Repos\\RichardAth\\Projects\\"
"bin\\x64\\Release";
#else
static std::string Path = "C:\\Users\\admin99\\Source\\Repos\\RichardAth\\Projects\\"
"bin\\x64\\Debug";
#endif
//static std::string Path = "C:\\Users\\admin99\\Documents\\Downloads_long_term_storage"
//"\\yafu-1.34-src\\yafu-1.34.3\\bin\\x64\\release";
static std::string yafuprog = "yafu-x64.exe ";
static std::string logPath = "C:/users/admin99/factors.txt";
//static std::string logPath = "factors.txt";
std::string batfilename = "YafurunTemp.bat";
bool yafu = true;  // true: use YAFU. false: use built-in ECM and SIQS or Msieve

/*process YAFU commands */
void yafuParam(std::string command) {
	std::string param = command.substr(4);  /* remove "YAFU" */
	while (param[0] == ' ')
		param.erase(0, 1);              /* remove leading space(s) */

	if (param == "ON") {
		yafu = true;
		msieve = false;
	}
	else if (param == "OFF") {
		yafu = false;
	}
	else if (param == "PATH") {
		std::cout << "path = " << Path << '\n';
		/* todo; allow command to change path */
	}
	else if (param == "LOG") {
		std::cout << "log file = " << logPath << '\n';
		/* todo; allow command to change log file name */
	}
	
	else {
		std::cout << "invalid YAFU command (use ON, OFF, PATH or LOG \n";
	}

	/* to be completed */
}

/* generate short batch file to run YAFU */
void genfile(const std::string &numStr) {
	FILE * batfile;
	std::string buffer;

	int rv = fopen_s(&batfile, batfilename.data(), "w");
	if (rv != 0) {
		if (errno != 0) {
			perror("error opening batch file");
			abort();
		}
	}

	buffer = "pushd " + Path + '\n';  // change directory
	rv = fputs(buffer.data(), batfile);
	assert(rv >= 0);

	buffer = yafuprog;
	buffer += " factor(" + numStr + ")";
	buffer += " -p";                // set batch priority
	for (int i = 1; i <= verbose; i++)
		buffer += " -v";                    // set verbose mode for YAFU

	buffer += " | wintee " + logPath + '\n';
	if (verbose > 0) {
		//static char time[10];
		//_strtime_s(time, sizeof(time));  // get current time as "hh:mm:ss"
		std::cout << myTime() << " command is: \n" << buffer ;  
	}
	rv = fputs(buffer.data(), batfile);
	assert(rv >= 0);

	rv = fputs("popd \n", batfile);    // change directory back
	assert(rv >= 0);

	fclose(batfile);
}

bool callYafu(Znum num, std::vector<zFactors>&Factors) {
	int rv;
	std::string command = Path + yafuprog;
	std::string numStr;
	std::string buffer;
	std::string logfile = logPath;
	//FILE *log;
	int fcount = 0;


	size_t numdigits = mpz_sizeinbase(ZT(num), 10);  // get number of decimal digits in num
	numStr.resize(numdigits + 5);             // resize buffer
	mpz_get_str(&numStr[0], 10, ZT(num));     // convert num to decimal (ascii) digits
	numStr.resize(strlen(&numStr[0]));        // get exact size of string in bufer
	genfile(numStr);

	//command += " factor(" + numStr + ")";      // append number to be factored to command line

	int rc = remove(logfile.data());
	if (rc != 0 && errno != ENOENT) {
		perror("could not remove old YAFU log file ");
	}

	rv = system(batfilename.data());             // start YAFU;

	/* get control back when YAFU has finished */
	if (rv == -1) {
		std::cout << "cannot start YAFU errno = " << errno << '\n';
		return false;
	}
	else if (rv != 0) {
		std::cout << "cannot start YAFU return code = " << rv << '\n';
		return false;
	}

	std::ifstream logStr(logfile, std::ios::in);  // open log file for input
	if (!logStr.is_open()) {
		std::cout << "cannot open YAFU log file \n";
		return false;
	}

	while (true) {            // read log file, find prime factors, ignore everything else
		Znum f;
		int i;
		if (!std::getline(logStr, buffer))   // read 1 line into buffer
			break;                   // exit loop when end of file reached 

		if (buffer[0] == 'P') {      // ignore log entry unless it's a prime factor
						 // note: some other log entries also have a 'p' in this position
			int  len = 0;
			for (i = 1; isdigit(buffer[i]); i++)  // get number of ascii digits in factor
				len = len * 10 + (buffer[i] - '0');
			if (len == 0)
				continue;                    // p not followed by a number
			i += 3;
			auto rv = mpz_set_str(ZT(f), buffer.substr(i, len).c_str(), 10);
			insertBigFactor(Factors, f);     // put f into factor list
			fcount++;                        // increase count of factors found
		}
	}

		logStr.close();
		if (fcount > 0)
			return true;    // success
		else
			return false;   // failure
}