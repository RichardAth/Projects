#include "pch.h"

#include "showtime.h"

extern HWND handConsole;      /* handle to console window */

/* Note: the path specified here can be overwritten by a path specified in the .ini file */
#ifndef _DEBUG
std::string YafuPath = "C:\\Users\\admin99\\Source\\Repos\\RichardAth\\Projects\\"
"bin\\x64\\Release";
#else
std::string YafuPath = "C:\\Users\\admin99\\Source\\Repos\\RichardAth\\Projects\\"
"bin\\x64\\Debug";
#endif

std::string yafuprog = "yafu-x64.exe";
static std::string logPath = "C:/users/admin99/factors.txt";

static std::string batfilename = "YafurunTemp.bat";
bool yafu = true;  // true: use YAFU. false: use built-in ECM and SIQS or Msieve
int pvalue = 4;    // 4 = PLAN NORMAL (default)

/* delete file specified by path + FileName */
void delfile(const std::string &path,  const char * FileName) {
	std::string fname;
	struct __stat64 fileStat;

	if (!path.empty())
		fname = path + "\\" + FileName;
	else
		fname = FileName;

	int err = _stat64(fname.data(), &fileStat);
	if (0 != err) 
		return;

	auto  fsize = fileStat.st_size/1024 ;
	std::cout << FileName << " size is " << fsize << " KB \n";

	int rc = remove(fname.data());
	if (rc != 0 && errno != ENOENT) {
		perror("could not remove file ");
	}
	else std::cout << "removed file: " << FileName << '\n';
}

/* use windows explorer-type dialogue box to select file */
char * getFileName(const char *filter, HWND owner) {

	OPENFILENAMEA ofn;
	static char * fileNameStr = NULL;
	static char afilename[MAX_PATH] = "";
	static char fileName[MAX_PATH] = "";

	fileNameStr = NULL;
	ZeroMemory(&ofn, sizeof(ofn));
	ofn.lStructSize = sizeof(OPENFILENAMEA);	// length of the structure
	ofn.hwndOwner = owner;						// handle of owner window
	ofn.lpstrFilter = filter;					// filters to select specific file types
	ofn.lpstrFile = fileName;					// on return, fileName contains the name of the file
	ofn.nMaxFile = MAX_PATH;					// length of buffer to contain the file name
	ofn.nMaxFileTitle = 0;						// length of lpstrFileTitle buffer
	ofn.lpstrFileTitle = NULL;					// file name and extension (without path) of selected file
	ofn.lpstrInitialDir = NULL;					// initial directory
	ofn.nFilterIndex = 1;						// do not use the custom filter
	ofn.Flags = OFN_PATHMUSTEXIST | OFN_FILEMUSTEXIST | OFN_HIDEREADONLY | OFN_EXPLORER;
	ofn.lpstrDefExt = "";						// default file extension


	if (GetOpenFileNameA(&ofn)) {
		fileNameStr = afilename;
		strcpy_s(afilename, sizeof(afilename), fileName);
	}

	return fileNameStr;
}

/* replace path and/or program name, 
return true if change actually made. 
otherwise return false */
bool changepath(std::string &path, std::string &prog) {
	bool rewrite = false;
	char * newpathC;
	std::string newpath;
	std::string newprog;
	std::string::size_type n;

	newpathC = getFileName("Program Files\0*.exe\0\0", handConsole);
	if (newpathC == NULL) {
		std::cout << "command cancelled \n";
		return false;
	}
	newpath = newpathC;  // copy C-style string to std::string
	n = newpath.rfind('\\');   // look for '\'
	if (n == std::string::npos) {
		std::cout << "command cancelled \n";
		return false;
	}
	else {
		newprog = newpath.substr(n + 1);  // copy characters after last '\' 
		newpath.erase(n);  // erase last \ and any following characters
		if (newpath != path) {
			std::cout << "old path: " << path << '\n';
			std::cout << "new path: " << newpath << '\n';
			path = newpath;
			rewrite = true;
		}
		if (!prog.empty() && newprog != prog) {
			std::cout << "old prog: " << prog << '\n';
			std::cout << "new prog: " << newprog << '\n';
			prog = newprog;
			rewrite = true;
		}
		//if (rewrite)
		//	writeIni();  // write any changes to .ini file
	}
	return rewrite;
}

/* check file status. Print date & time modified */
void fileStatus(const std::string &progname) {
	struct __stat64 fileStat;
	struct tm ftimetm;
	char timestamp[23];   // date & time in format "dd/mm/yyyy at hh:mm:ss"
	int err = _stat64(progname.data(), &fileStat);
	if (err == 0) {
		auto ftime = fileStat.st_mtime;  // time last modified in time_t format
		localtime_s(&ftimetm, &ftime);   // convert to tm format
		/* convert to dd/mm/yyyy hh:mm:ss */
		strftime(timestamp, sizeof(timestamp), "%d/%m/%C%y at %H:%M:%S", &ftimetm);
		std::cout << progname << "  modified on ";
		std::cout << timestamp << '\n';
	}
	else
		std::cout << progname << " not found \n";
}

/* check or replace yafu.ini file */
static void inifile(std::string &param) {
	std::string yafuini = "yafu.ini";
	std::string iniFname = YafuPath + "\\" + yafuini;
	std::string buffer;
	std::string ggnfsPath;
	int lineNo = 1;
	int ggnfsLine = -1;
	std::vector<std::string> inifile;  // copy contents of yafu.ini here

	std::ifstream iniStr(iniFname, std::ios::in);  // open yafu.ini file
	if (!iniStr.is_open()) {
		std::cout << " cannot open yafu.ini \n";
	}
	else {    /* read in yafu.ini file */
		while (std::getline(iniStr, buffer)) {
			if (_strnicmp("ggnfs_dir=", buffer.c_str(), 10) == 0) {
				ggnfsLine = lineNo;
				ggnfsPath = buffer.substr(10); // copy path following '=' character
			}

			if (_strnicmp("%print", buffer.c_str(), 6) == 0)
				std::cout << lineNo << ": " << buffer << '\n';

			inifile.push_back(buffer);
			lineNo++;
		}
		iniStr.close();
		std::cout << "yafu.ini contains " << lineNo - 1 << " lines \n";

		if (ggnfsLine > 0) {
			/* ggnfs_dir was found */
			std::cout << ggnfsLine << ": " << inifile[ggnfsLine - 1] << '\n';
			for (int i = 11; i <= 16; i++) {
				char ia[3];  /* 11 to 16 as ascii text */

				_itoa_s(i, ia, sizeof(ia), 10);   /* convert i from binary to ascii */
				std::string progname = ggnfsPath;
				progname += "gnfs-lasieve4i";
				progname += +ia;
				progname += "e.exe";
				fileStatus(progname);
			}
		}
		else
			std::cout << "ggnfs_dir parameter not found \n";
	}

	/* should we change/add the ggnfs_dir parameter?*/
	if (!param.empty() && toupper(param[0]) == 'I') {
		std::string newFname = YafuPath + "\\" + "yafu.new";
		std::string oldFname = YafuPath + "\\" + "yafu.old";

		/* replace the ggnfs_dir parameter */
		std::string newprog;
		if (!changepath(ggnfsPath, newprog))
			return;  // bail out if change path cancelled

		buffer = "ggnfs_dir=" + ggnfsPath + '/' + '\n';

		std::ofstream newStr(newFname, std::ios::out);  // open yafu.new file for output
		if (!newStr.is_open()) {
			std::cout << "cannot open yafu.new \n";
			return;
		}
		/* copy everything from the yafu.ini to the yafu.new file except 
		   ggnfs_dir= parameter */
		for (size_t ix = 0; ix < (lineNo - 1); ix++) {
			if (ix != (ggnfsLine - 1))
				newStr << inifile[ix] << '\n';
			else
				newStr << buffer;   // replace ggnfs_dir parameter with new value
		}

		if (ggnfsLine <= 0)
			newStr << buffer;    // append ggnfs_dir parameter at end of file

		newStr.close();

		delfile(YafuPath, oldFname.c_str());  // delete any previous yafu.old
		int rv = rename(iniFname.c_str(), oldFname.c_str());   // yafu.ini -> yafu.old
		if (rv == 0)
			rename(newFname.c_str(), iniFname.c_str());   // yafu.new -> yafu.ini
		else
			perror("unable to rename yafu.ini as yafu.old");
	}
}

/*process YAFU commands */
void yafuParam(const std::string &command) {
	std::string param = command.substr(4);  /* remove "YAFU" */
	const std::vector<std::string> paramList = {
		"ON", "OFF", "PATH", "LOG", "PLAN", "TIDY", "INI"
	};

	size_t ix = 0;

	while (param[0] == ' ')
		param.erase(0, 1);              /* remove leading space(s) */

	/* match parameter value to a list entry */
	for (ix = 0; ix < paramList.size(); ix++) {
		if (param.size() < paramList[ix].size()) continue;
		if (paramList[ix].compare(param.substr(0, paramList[ix].size())) == 0)
			break;
	}

	switch (ix)  /*  switch according to parameter value */ {
	case 0: /* YAFU ON   */ {   
			yafu = true;
			msieve = false;
			Pari = false;
			break;
		}
	case 1: /* YAFU OFF  */ {    
			yafu = false;
			break;
		}
	case 2: /* YAFU PATH */ { 
		param.erase(0, 4);  // get rid of "PATH"
		while (param[0] == ' ')
			param.erase(0, 1);              /* remove leading space(s) */

		if (param   == "SET") {
			if (changepath(YafuPath, yafuprog))
				writeIni();  // rewrite .ini file
		}
		else {
			std::cout << "path = " << YafuPath << '\n';
		}

		fileStatus(YafuPath + '\\' + yafuprog);
		break;

	}
	case 3: /* YAFU LOG  */ {  
			std::cout << "log file = " << logPath << '\n';
			/* todo; allow command to change log file name */
			break;
		}
	case 4: /* YAFU PLAN */ { 
		// plan name can be NONE, NOECM, LIGHT, NORMAL, DEEP
		// CUSTOM is not supported, NORMAL is default
		param.erase(0, 4);  // get rid of "PLAN"
		while (param[0] == ' ')
			param.erase(0, 1);              /* remove leading space(s) */
		if (param == "NONE")         pvalue = 1;
		else if (param == "NOECM")   pvalue = 2;
		else if (param == "LIGHT")   pvalue = 3;
		else if (param == "NORMAL")  pvalue = 4;
		else if (param == "DEEP")    pvalue = 5;
		else std::cout
			<< "YAFU PLAN value invalid; use NONE, NOECM, LIGHT, NORMAL, or DEEP\n";
		break;
	}
	case 5: /* YAFU TIDY */ {
		delfile(YafuPath, "factor.log");
		delfile(YafuPath, "session.log");
		delfile(YafuPath, "siqs.dat");
		delfile(YafuPath, "nfs.log");
		delfile(YafuPath, "nfs.job");
		delfile(YafuPath, "ggnfs.log");
		delfile(YafuPath, "YAFU_get_poly_score.out");
		break;
	}
	case 6: /* YAFU INI  */ {
		param.erase(0, 3);  // get rid of "INI"
		while (param[0] == ' ')
			param.erase(0, 1);              /* remove leading space(s) */

		inifile(param);  // process "YAFU INI .... " command
		break;
	}

	default: {
			std::cout << "invalid YAFU command (use ON, OFF, PATH, LOG, PLAN, TIDY or INI \n";
			break;
		}
	}
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

	buffer = "pushd " + YafuPath + '\n';  // change directory
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
	std::string command = YafuPath + yafuprog;
	std::string numStr;
	std::string buffer;
	std::string logfile = logPath;
	//FILE *log;
	int fcount = 0;


	size_t numdigits = mpz_sizeinbase(ZT(num), 10);  // get number of decimal digits in num
	numStr.resize(numdigits + 5);             // resize buffer
	mpz_get_str(&numStr[0], 10, ZT(num));     // convert num to decimal (ascii) digits
	numStr.resize(strlen(&numStr[0]));        // get exact size of string in buffer
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