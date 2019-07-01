#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <Windows.h>
#include "factor.h"

//
bool msieve = true;
std::string callPath = "C:/Users/admin99/Documents/Downloads_long_term_storage"
"/msieve-code-r1030-trunk/bin/x64/Release/msieve.exe";
std::string logPath = "C:\\users\\admin99\\msieve.log ";
std::string options = " -e ";   // perform 'deep' ECM, seek factors > 15 digits
std::string redirectOP = " >\\.\\pipe\\StdOutPipe 2>\\.\\pipe\\StdErrPipe ";

/* process MSIEVE commands */
void msieveParam(std::string expupper) {
	std::string param = expupper.substr(6);  /* remove "MSIEVE */
	while (param[0] == ' ')
		param.erase(0, 1);              /* remove leading space */

	if (param == "ON")
		msieve = true;
	else if (param == "OFF")
		msieve = false;
	else if (param == "PATH")
		std::cout << "path = " << callPath << '\n';
	/* todo; allow command to change path */
	else if (param == "LOG")
		std::cout << "log file = " << logPath << '\n';
	/* todo; allow command to change log file name */

	/* to be completed */
}


int openPipe(const char PipeName[], HANDLE *hPipe) {

	while (1) {
		*hPipe = CreateNamedPipeA(PipeName,                     // name
			PIPE_ACCESS_DUPLEX,                                // open mode
			PIPE_TYPE_BYTE | PIPE_READMODE_BYTE | PIPE_WAIT,   // pipe mode 
			// FILE_FLAG_FIRST_PIPE_INSTANCE is not needed but forces CreateNamedPipe(..) 
			// to fail if the pipe already exists...
			1,                                             // maximum number of instances             
			1024 * 16,                                     // out buffer size
			1024 * 16,                                     // in buffer size
			NMPWAIT_USE_DEFAULT_WAIT,                      // default timeout
			NULL);                                         // security attributes

	    // Break if the pipe handle is valid. 
		if (hPipe != INVALID_HANDLE_VALUE)
			break;

		// Exit if an error other than ERROR_PIPE_BUSY occurs. 
		if (GetLastError() != ERROR_PIPE_BUSY) 	{
			printf("Could not open pipe. GLE=%d\n", GetLastError());
			return -1;
		}

		// All pipe instances are busy, so wait for 20 seconds. 
		if (!WaitNamedPipeA(PipeName, 20000)) {
			printf("Could not open pipe: 20 second wait timed out.");
			return -1;
		}
	}

	return 0;
}

/* use msieve to factorise num. msieve places its results in a log file.
All entries in the log file start with a time stamp such as "Sun Jun 30 16:07:40 2019". 
For prime factors this is followed by text such as "p9 factor: 193707721" 
where p designates a prime factor, the 9 is the number of decimal digits in the factor
and 193707721 is the factor itself in decimal.

The whole interface to msieve is a bit of a kludge but it works. Integrating them
properly would be so tedious it's not worth it. I had to make several kludges
to build msieve, and the prebuilt msieve wouldn't work */
bool callMsieve(Znum num, std::vector<zFactors>&Factors) {
	//HANDLE logP;
	std::string command = callPath + options + " -l " + logPath;
	std::string numStr;
	size_t numdigits = mpz_sizeinbase(ZT(num), 10);  // get number of decimal digits in num
	int rv, fcount = 0;
	FILE *log;
	std::string buffer;

	rv = fopen_s(&log, logPath.data(), "w");      
	if (rv != 0) {
		if (errno != 0)
			perror(NULL);
	}
	else
		fclose(log);     // if log exists, erase its contents

	numStr.resize(numdigits + 5);             // resize buffer
	mpz_get_str(&numStr[0], 10, ZT(num));     // convert num to decimal (ascii) digits
	numStr.resize(strlen(&numStr[0]));        // get exact size of string in bufer
	command += " ";
	command += numStr;                        // append number to be factored to command line
	//std::cout << "command is: \n" << command << '\n';  // temp
	rv = system(command.data());             // start msieve;

	/* get control back when msieve has finished */
	if (rv ==  -1) {
		std::cout << "cannot start msieve errno = " << errno << '\n';
		return false;
	}
	else if (rv != 0) {
		std::cout << "cannot start msieve return code = "<<  rv << '\n';
		return false;
	}

	std::ifstream logStr(logPath, std::ios::in);  // open log file for input
	if (!logStr.is_open()) {
		std::cout << "cannot open msieve log file \n";
		return false;
	}

/* msieve has placed its results in a log file.
All entries in the log file start with a time stamp such as "Sun Jun 30 16:07:40 2019". 
For prime factors this is followed by text such as "p9 factor: 193707721" 
where p designates a prime factor, the 9 is the number of decimal digits in the factor
and 193707721 is the factor itself in decimal.*/

	while (true){           // read log file, find prime factors
		Znum f;
		int i;
		if (!std::getline(logStr, buffer))
			break;          // exit loop when end of file reached 
		if (buffer.size() < 38)
			continue;
		if (buffer[26] == 'p') {  // ignore log entry unless it's a prime factor
			//std::cout << buffer << '\n';   // temp
			int  len = 0;
			for (i = 27; isdigit(buffer[i]); i++)
				len = len * 10 + (buffer[i] - '0');
			i += 9;        // move past text " factor: "
			//std::cout << "factor = " << buffer.substr(i, len).c_str() << '\n';
			/* convert factor from ascii to Znum */
			auto rv = mpz_set_str(ZT(f), buffer.substr(i, len).c_str(), 10);
			//std::cout << "factor = " << f << '\n';
			insertBigFactor(Factors, f);     // put f into factor list
			fcount++;
		}
	}

	fclose(log);
	if (fcount > 0)
		return true;    // success
	else
		return false;   // failure
}