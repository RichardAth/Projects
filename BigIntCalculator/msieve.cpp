#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdio>
#define WIN32_LEAN_AND_MEAN
#include <Windows.h>
#include "factor.h"

/* output from Msieve -h option:

Msieve v. 1.54 (SVN 945)

usage: C:/Users/admin99/Documents/Downloads_long_term_storage/msieve-code-r1030-trunk/bin/x64/Release/msieve.exe [options] [one_number]

numbers starting with '0' are treated as octal,
numbers starting with '0x' are treated as hexadecimal

options:
   -s <name> save intermediate results to <name>
			 instead of the default msieve.dat
   -l <name> append log information to <name>
			 instead of the default msieve.log
   -i <name> read one or more integers to factor from
			 <name> (default worktodo.ini) instead of
			 from the command line
   -m        manual mode: enter numbers via standard input
   -q        quiet: do not generate any log information,
			 only print any factors found
   -d <min>  deadline: if still sieving after <min>
			 minutes, shut down gracefully (default off)
   -r <num>  stop sieving after finding <num> relations
   -p        run at idle priority
   -v        verbose: write log information to screen
			 as well as to logfile
   -z        you are Paul Zimmermann
   -t <num>  use at most <num> threads

 elliptic curve options:
   -e        perform 'deep' ECM, seek factors > 15 digits

 quadratic sieve options:
   -c        client: only perform sieving

 number field sieve options:

		   [nfs_phase] "arguments"

 where the first part is one or more of:
   -n        use the number field sieve (80+ digits only;
			 performs all NFS tasks in order)
   -nf <name> read from / write to NFS factor base file
			 <name> instead of the default msieve.fb
   -np       perform only NFS polynomial selection
   -np1      perform stage 1 of NFS polynomial selection
   -nps      perform NFS polynomial size optimization
   -npr      perform NFS polynomial root optimization
   -ns       perform only NFS sieving
   -nc       perform only NFS combining (all phases)
   -nc1      perform only NFS filtering
   -nc2      perform only NFS linear algebra
   -ncr      perform only NFS linear algebra, restarting
			 from a previous checkpoint
   -nc3      perform only NFS square root

 the arguments are a space-delimited list of:
 polynomial selection options:
   polydegree=X    select polynomials with degree X
   min_coeff=X     minimum leading coefficient to search
				   in stage 1
   max_coeff=X     maximum leading coefficient to search
				   in stage 1
   stage1_norm=X   the maximum norm value for stage 1
   stage2_norm=X   the maximum norm value for stage 2
   min_evalue=X    the minimum score of saved polyomials
   poly_deadline=X stop searching after X seconds (0 means
				   search forever)
   X,Y             same as 'min_coeff=X max_coeff=Y'
 line sieving options:
   X,Y             handle sieve lines X to Y inclusive
 filtering options:
   filter_mem_mb=X  try to limit filtering memory use to
					X megabytes
   filter_maxrels=X limit the filtering to using the first
					X relations in the data file
   filter_lpbound=X have filtering start by only looking
					at ideals of size X or larger
   target_density=X attempt to produce a matrix with X
					entries per column
   max_weight=X     have filtering start by looking at ideals
					of max weight >= X
   X,Y              same as 'filter_lpbound=X filter_maxrels=Y'
 linear algebra options:
   skip_matbuild=1  start the linear algebra but skip building
					the matrix (assumes it is built already)
   la_block=X       use a block size of X (512<=X<=65536)
   la_superblock=X  use a superblock size of X
   cado_filter=1    assume filtering used the CADO-NFS suite
 square root options:
   dep_first=X start with dependency X, 1<=X<=64
   dep_last=Y  end with dependency Y, 1<=Y<=64
   X,Y         same as 'dep_first=X dep_last=Y'

*/

bool msieve = true;  // true: use Msieve. false: use built-in ECM and SIQS
bool eopt = true;    // set -e option in Msieve: perform 'deep' ECM, seek factors > 15 digits
bool nopt = true;    // set -n option in Msieve: use the number field sieve (80+ digits only;
				     //        performs all NFS tasks in order)

std::string callPath = "C:/Users/admin99/Documents/Downloads_long_term_storage"
"/msieve-code-r1030-trunk/bin/x64/Release/msieve.exe ";
std::string logPath = "C:\\users\\admin99\\msieve.log ";
std::string options = " -e ";   // perform 'deep' ECM, seek factors > 15 digits
std::string redirectOP = " >\\.\\pipe\\StdOutPipe 2>\\.\\pipe\\StdErrPipe ";

/* process MSIEVE commands */
void msieveParam(std::string command) {
	std::string param = command.substr(6);  /* remove "MSIEVE */
	while (param[0] == ' ')
		param.erase(0, 1);              /* remove leading space(s) */

	if (param == "ON")
		msieve = true;
	else if (param == "OFF")
		msieve = false;
	else if (param == "PATH") {
		std::cout << "path = " << callPath << '\n';
		/* todo; allow command to change path */
	}
	else if (param == "LOG") {
		std::cout << "log file = " << logPath << '\n';
		/* todo; allow command to change log file name */
	}
	else if (param == "E ON")
		eopt = true;
	else if (param == "E OFF")
		eopt = false;
	else if (param == "N ON")
		nopt = true;
	else if (param == "N OFF")
		nopt = false;

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

/* use Msieve to factorise num. Msieve places its results in a log file.
All entries in the log file start with a time stamp such as "Sun Jun 30 16:07:40 2019". 
For prime factors this is followed by text such as "p9 factor: 193707721" 
where p designates a prime factor, the 9 is the number of decimal digits in the factor
and 193707721 is the factor itself in decimal.

The whole interface to Msieve is a bit of a kludge but it works. Integrating them
properly would be so tedious it's not worth it. I had to make several kludges
to build msieve, and the prebuilt msieve wouldn't work */
bool callMsieve(Znum num, std::vector<zFactors>&Factors) {
	/* set up command to invoke Msieve */
	std::string command = callPath;
	std::string numStr;
	int rv, fcount = 0;
	FILE *log;
	std::string buffer;

	if (eopt)
		command += options;       // add -e option
	if (nopt)
		command += " -n ";
	command += " -l " + logPath + " ";

	size_t numdigits = mpz_sizeinbase(ZT(num), 10);  // get number of decimal digits in num
	numStr.resize(numdigits + 5);             // resize buffer
	mpz_get_str(&numStr[0], 10, ZT(num));     // convert num to decimal (ascii) digits
	numStr.resize(strlen(&numStr[0]));        // get exact size of string in bufer
	command += numStr;                        // append number to be factored to command line
	//std::cout << "command is: \n" << command << '\n';  // temp

	/* open earlier Msieve log file, if any exists. Replace it with an empty file */
	rv = fopen_s(&log, logPath.data(), "w");
	if (rv != 0) {
		if (errno != 0)
			perror(NULL);
	}
	else
		fclose(log);     // if log exists, erase its contents 

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

	while (true){            // read log file, find prime factors, ignore everything else
		Znum f;
		int i;
		if (!std::getline(logStr, buffer))   // read 1 line into buffer
			break;                   // exit loop when end of file reached 
		if (buffer.size() < 38)
			continue;
		if (buffer[26] == 'p') {      // ignore log entry unless it's a prime factor
			             // note: some other log entries also have a 'p' in this position
			int  len = 0;  
			for (i = 27; isdigit(buffer[i]); i++)  // get number of ascii digits in factor
				len = len * 10 + (buffer[i] - '0');
			if (len == 0)
				continue;                    // p not followed by a number
			if (buffer.substr(i, 9) != " factor: ")
				continue;                    // pnn not followed by " factor: "
			i += 9;                          // move past text " factor: "
			/* convert factor from ascii to Znum */
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