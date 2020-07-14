#include "pch.h"
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
int pvalue = 4;    // 4 = PLAN NORMAL (default)

/*process YAFU commands */
void yafuParam(const std::string &command) {
	std::string param = command.substr(4);  /* remove "YAFU" */
	const std::vector<std::string> paramList = {
		"ON", "OFF", "PATH", "LOG", "PLAN"
	};

	ptrdiff_t ix = 0;

	while (param[0] == ' ')
		param.erase(0, 1);              /* remove leading space(s) */

	/* match parameter value to a list entry */
	for (ix = 0; ix < paramList.size(); ix++) {
		if (param.size() < paramList[ix].size()) continue;
		if (paramList[ix].compare(param.substr(0, paramList[ix].size())) == 0)
			break;
	}

	switch (ix) {
	case 0:{    // YAFU ON
			yafu = true;
			msieve = false;
			break;
		}
	case 1: {    // YAFU OFF
			yafu = false;
			break;
		}
	case 2: {     // YAFU PATH
			std::cout << "path = " << Path << '\n';
			/* todo; allow command to change path */
			break;
		}
	case 3: {     // YAFU LOG
			std::cout << "log file = " << logPath << '\n';
			/* todo; allow command to change log file name */
			break;
		}
	case 4: {   // YAFU PLAN
		// plan name can be NONE, NOECM, LIGHT, NORMAL, DEEP
		// CUSTOM is not supported, NORMAL is default
		param.erase(0, 4);  // get rid of "PLAN"
		while (param[0] == ' ')
			param.erase(0, 1);              /* remove leading space(s) */
		if (param == "NONE")         pvalue = 1;
		else if (param == "NOECM")   pvalue = 2;
		else if (param == "LIGHT")   pvalue = 3;
		else if (param == "NORMAL")  pvalue = 4;
		else if (param == "DEEP")    pvalue = 4;
		else std::cout
			<< "YAFU PLAN value invalid; use NONE, NOECM, LIGHT, NORMAL, or DEEP\n";
		break;
	}

	default: {
			std::cout << "invalid YAFU command (use ON, OFF, PATH, LOG or PLAN \n";
			break;
		}
	}
	/* to be completed */
}

/* generate short batch file to run YAFU */
static void genfile(const std::string &numStr) {
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

	if (pvalue != 4) {
		switch (pvalue) {
		case 1: buffer += " -plan none"; break;
		case 2: buffer += " -plan noecm"; break;
		case 3: buffer += " -plan light"; break;
		case 4: buffer += " -plan normal"; break;
		case 5: buffer += " -plan deep"; break;
		default: break;  // ignore invalid pvalue
		}
	}
	buffer += " -of " + logPath + '\n';   

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

bool callYafu(const Znum &num, fList &Factors) {
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

	/* the format of the text file containing the factors is as follows:
	The number to be factored is in brackets followed by a / 
	The factors follow separated by / characters. If a factor 
	is repeated this is indicated by a ^ for exponentiation followed by the 
	exponent.
	There is no / after the last factor.
	e.g. (531069)/3/7/11^3/19   
	The file consists of one line. Before calling YAFU the old file is deleted 
	because otherwise YAFU would append the new result at the end of the file. */

	std::ifstream factorStr(logfile, std::ios::in);  // open factor file
	if (!factorStr.is_open()) {
		std::cout << "cannot open YAFU factor file \n";
		return false;
	}
	if (!std::getline(factorStr, buffer)) {  // read 1 line into buffer
		std::cout << "cannot read YAFU factor file \n";
		factorStr.close();
		return false;
	}

	auto pos1 = buffer.find("/");
	size_t pos2;
	ptrdiff_t len;
	if (pos1 == std::string::npos) {
		std::cout << "cannot find any factors in YAFU factor file \n";
		factorStr.close();
		return false;
	}
	auto factors = buffer.substr(pos1 + 1);

	while (true) {
		Znum f;
		pos2 = factors.find_first_of("/^");  // find end of factor text string
		if (pos2 != std::string::npos)
			len = pos2;
		else
			len = factors.size();
		auto rv = mpz_set_str(ZT(f), factors.substr(0, len).c_str(), 10);
		assert(rv == 0);
		insertBigFactor(Factors, f);     // put f into factor list
		fcount++;                        // increase count of factors found
		if (pos2 == std::string::npos)
			break;       // pos2 at end of buffer; no more factors
		if (factors[pos2] != '/')
			pos2 = factors.find("/");  // move pos2 past exponent following factor
		if (pos2 != std::string::npos)
			factors = factors.substr(pos2 + 1);  // remove factor just processed
		else break;    // pos2 at end of buffer; no more factors
	}
	factorStr.close();
	if (fcount > 0)
		return true;    // success
	else
		return false;   // failure
}