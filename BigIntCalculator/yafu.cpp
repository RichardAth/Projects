#include "pch.h"

#include "json.h"

/* forward declarations */
static void process_value_s(json_value* value, int depth, const char* name,
	const int index);

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
std::string outPath = "C:/users/admin99/factors.txt";

static std::string batfilename = "YafurunTemp.bat";
bool yafu = true;  // true: use YAFU. false: use built-in ECM and SIQS or Msieve
bool useOldYafu = false;
int pvalue = 4;    // 4 = PLAN NORMAL (default)

/* delete file specified by path + FileName */
void delfile(const std::string& path, const char* FileName) {
	std::string fname;
	struct __stat64 fileStat;

	if (!path.empty()) {
		fname = path;
		fname.append("\\");
		fname.append(FileName);
	}
	else
		fname = FileName;

	int err = _stat64(fname.data(), &fileStat);
	if (0 != err) {
		if (errno != ENOENT) {
			std::cout << "could not remove file " << fname.data() << ' '
#pragma warning(suppress : 4996)
				<< strerror(errno) << '\n';
		}
		return;
	}

	auto  fsize = fileStat.st_size/1024 ;
	std::cout << FileName << " size is " << fsize << " KB \n";

	int rc = remove(fname.data());
	if (rc != 0) {
		if (errno != ENOENT) {
			std::cout << "could not remove file " << fname.data() << ' '
#pragma warning(suppress : 4996)
				<< strerror(errno) << '\n';
		}
	}
	else std::cout << "removed file: " << FileName << '\n';
}

/* use windows explorer-type dialogue box to select file */
char * getFileName(const char *filter, HWND owner, bool MustExist) {

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
	ofn.Flags = OFN_PATHMUSTEXIST | OFN_HIDEREADONLY | OFN_EXPLORER;
	if (MustExist)
		ofn.Flags |= OFN_FILEMUSTEXIST;
	ofn.lpstrDefExt = "";						// default file extension


	if (GetOpenFileNameA(&ofn)) {
		fileNameStr = afilename;
		strcpy_s(afilename, sizeof(afilename), fileName);
	}

	return fileNameStr;
}

/* replace path and/or program name, 
return true if change actually made. otherwise return false */
bool changepathPP(std::string &path, std::string &prog) {
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

/* replace path, return true if change actually made, otherwise return false */
bool changepath2(std::string& path) {
	bool rewrite = false;
	char* newpathC;
	std::string newpath;
	std::string newprog;
	std::string::size_type n;

	newpathC = getFileName("Text Files\0*.TXT\0log files\0*.LOG\0All\0*.*\0", handConsole, false);
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
		if (newpath != path) {
			std::cout << "old path: " << path << '\n';
			std::cout << "new path: " << newpath << '\n';
			path = newpath;
			rewrite = true;
		}
	}
	return rewrite;
}

/* check file status. Print date & time modified, return false if file not found */
bool fileStatus(const std::string &fileName) {
	struct __stat64 fileStat;
	struct tm ftimetm;
	char timestamp[23];   // date & time in format "dd/mm/yyyy at hh:mm:ss"
	int err = _stat64(fileName.data(), &fileStat);
	if (err == 0) {
		auto ftime = fileStat.st_mtime;  // time last modified in time_t format
		localtime_s(&ftimetm, &ftime);   // convert to tm format
		/* convert to dd/mm/yyyy hh:mm:ss */
		strftime(timestamp, sizeof(timestamp), "%d/%m/%C%y at %H:%M:%S", &ftimetm);
		std::cout << fileName << "  modified on ";
		std::cout << timestamp << '\n';
		return true;
	}
	else {
		std::cout << "could not access " << fileName << " "
#pragma warning(suppress : 4996)
			<< strerror(errno)<<  "\n";
		return false;
	}
}

/* check or replace yafu.ini file */
static void inifile(const std::string &param) {
	std::string yafuini = "yafu.ini";
	std::string iniFname;

	std::string buffer;
	std::string ggnfsPath;
	int lineNo = 1;
	int ggnfsLine = -1;
	std::vector<std::string> inifile;  // copy contents of yafu.ini here
	std::string curPath(MAX_PATH, 0);
	std::string newFname;
	std::string oldFname;

	if (useOldYafu)
		iniFname = YafuPath + "\\" + yafuini;
	else {
		DWORD rv2 = GetCurrentDirectoryA(MAX_PATH, (LPSTR)curPath.data());
		assert(rv2 != 0);
		while (curPath.back() == 0)
			curPath.pop_back();  /* remove trailing null(s) */
		iniFname = curPath + "\\" + yafuini;
	}

	std::ifstream iniStr(iniFname, std::ios::in);  // open yafu.ini file
	if (!iniStr.is_open()) {
		std::cout << " cannot open " << iniFname <<  '\n';
	}

	else {    /* read in yafu.ini file */
		if (verbose > 0) {
			std::cout << iniFname << " opened \n";
		}
		while (std::getline(iniStr, buffer)) {
			if (_strnicmp("ggnfs_dir=", buffer.c_str(), 10) == 0) {
				ggnfsLine = lineNo;
				ggnfsPath = buffer.substr(10); // copy path following '=' character
			}

			if ((_strnicmp("%print", buffer.c_str(), 6) == 0) || (verbose > 0))
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


		/* should we change/add the ggnfs_dir parameter?*/
		if (!param.empty() && toupper(param[0]) == 'I') {
			if (useOldYafu) {
				newFname = YafuPath + "\\" + "yafu.new";
				oldFname = YafuPath + "\\" + "yafu.old";
			}
			else {
				newFname = curPath + "\\" + "yafu.new";
				oldFname = curPath + "\\" + "yafu.old";
			}

			/* replace the ggnfs_dir parameter */
			std::string newprog;
			if (!changepathPP(ggnfsPath, newprog))
				return;  // bail out if change path cancelled

			buffer = "ggnfs_dir=" + ggnfsPath + '/' + '\n';

			std::ofstream newStr(newFname, std::ios::out);  // open yafu.new file for output
			if (!newStr.is_open()) {
				std::cout << "cannot open " << newFname << '\n';
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

			delfile("", oldFname.c_str());  // delete any previous yafu.old
			int rv = rename(iniFname.c_str(), oldFname.c_str());   // yafu.ini -> yafu.old
			if (rv == 0)
				rename(newFname.c_str(), iniFname.c_str());   // yafu.new -> yafu.ini
			else
				perror("unable to rename yafu.ini as yafu.old");
		}
	}
}

/*process YAFU commands */ 
void yafuParam(const std::vector<std::string>& p) {
	const std::vector<std::string> paramList = {
		"ON", "OFF", "PATH", "OUT", "PLAN", "TIDY", "INI"
	};
	size_t ix = 0;
	std::string curPath(MAX_PATH, 0);  /* path to current directory*/

	if(p.size() < 2) {
		std::cout << "invalid YAFU command (use ON, OFF, PATH, OUT, PLAN, TIDY or INI \n";
		return;
	}

	/* match parameter value to a list entry */
	for (ix = 0; ix < paramList.size(); ix++) {
		if (paramList[ix] == p[1]) 
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
		std::cout << "path = " << YafuPath << '\n';
		if (p.size() >= 3 && p[2] == "SET") {
			if (changepathPP(YafuPath, yafuprog))
				writeIni();  // rewrite .ini file
		}
		fileStatus(YafuPath + '\\' + yafuprog);
		break;

	}
	case 3: /* YAFU OUT  */ {  
		std::cout << "YAFU output file = " << outPath << '\n';
		if (p.size() >=3 && p[2] == "SET") {
			if (changepath2(outPath))
				writeIni();  // rewrite .ini file
		}
		fileStatus(outPath);
		break;
	}
	case 4: /* YAFU PLAN */ { 
		// plan name can be NONE, NOECM, LIGHT, NORMAL, DEEP
		// CUSTOM is not supported, NORMAL is default
		if (p.size() >= 3) {
			if (p[2] == "NONE")         pvalue = 1;
			else if (p[2] == "NOECM")   pvalue = 2;
			else if (p[2] == "LIGHT")   pvalue = 3;
			else if (p[2] == "NORMAL")  pvalue = 4;
			else if (p[2] == "DEEP")    pvalue = 5;
			else std::cout
				<< "YAFU PLAN value invalid; use NONE, NOECM, LIGHT, NORMAL, or DEEP\n";
			break;
		}
	}
	case 5: /* YAFU TIDY */ {
		if (useOldYafu) {
			delfile(YafuPath, "factor.log");
			delfile(YafuPath, "session.log");
			delfile(YafuPath, "factor.json");
			delfile(YafuPath, "siqs.dat");
			delfile(YafuPath, "nfs.job");
			delfile(YafuPath, "ggnfs.log");
			delfile(YafuPath, "nfs.dat");
			delfile(YafuPath, "nfs.fb");
			delfile(YafuPath, "nfs.log");
			delfile(YafuPath, "nfs.dat.mat");
			delfile(YafuPath, "YAFU_get_poly_score.out");
		}
		else {
			DWORD rv2 = GetCurrentDirectoryA(MAX_PATH, (LPSTR)curPath.data());
			assert(rv2 != 0);
			while (curPath.back() == 0)
				curPath.pop_back();  /* remove trailing null(s) */
			delfile(curPath, "factor.log");
			delfile(curPath, "session.log");
			delfile(curPath, "factor.json");
			delfile(curPath, "siqs.dat");
			delfile(curPath, "nfs.job");
			delfile(curPath, "ggnfs.log");
			delfile(curPath, "nfs.dat");
			delfile(curPath, "nfs.fb");
			delfile(curPath, "nfs.log");
			delfile(curPath, "nfs.dat.mat");
			delfile(curPath, "YAFU_get_poly_score.out");
		}
		break;
	}
	case 6: /* YAFU INI  */ {
		if (p.size() >= 3)
			inifile(p[2]);  // process "YAFU INI .... " command
		else
			inifile("");
		break;
	}

	default: {
			std::cout << "invalid YAFU command (use ON, OFF, PATH, OUT, PLAN, TIDY or INI \n";
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
	buffer += " -of " + outPath + '\n';   

	if (verbose > 0) {
		std::cout << myTime() << " command is: \n" << buffer ;  
	}

	rv = fputs(buffer.data(), batfile);
	assert(rv >= 0);

	rv = fputs("popd \n", batfile);    // change directory back
	assert(rv >= 0);

	fclose(batfile);
}

bool callYafuOld(const Znum &num, fList &Factors) {
	int rv;
	std::string command = YafuPath + yafuprog;
	std::string numStr;
	std::string buffer;
	std::string logfile = outPath;
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
	delfile(YafuPath, "nfs.dat");  /* if earlier run leaves this file undeleted
	                               it would cause problems */
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

Znum ToBeFactored;
std::vector <Znum> factors;

/* check whether the product of all the factors is equal to the number to be
factored */
static bool sanityCheck(const Znum& ToBeFactored, const std::vector<Znum>& factors, 
	const Znum &num) {
	Znum residue = ToBeFactored;
	Znum remainder;

	if (ToBeFactored != num) {
		gmp_printf("number to be factored from json record = %Zd\n"
			"original number to be factored        = %Zd\n"
			"We appear to be processing the wrong json record \n",
			ToBeFactored, num);
	}

	for (auto f : factors) {
		/* divide by each of the factors in turn, checking the remainder each time */
		mpz_fdiv_qr(ZT(residue), ZT(remainder), ZT(residue), ZT(f));
		assert(remainder == 0);
	}
	if (residue != 1) {
		gmp_printf("YAFU: factors OK but residue after removing factors  = %Zd \n",
			ZT(residue));
		return false;
	}
	else {
		if (verbose > 0)
			printf("YAFU: all factors found \n");
		return true;
	}
}

/* save specified values. At the moment, we are only looking for the number
to be factored and the factors found. These are decimal mumbers saved as text
strings */
static void save_value(json_value* value, const int index, int depth) {
	json_type t = value->type;
	char* sp;
	Znum numValue;
	int rv;

	assert(t == json_string); /* at the moment, only strings can be processed. */
	sp = value->u.string.ptr;
	rv = mpz_set_str(ZT(numValue), sp, 10); /* convert ascii decimal to binary */
	if (rv != 0) {
		fprintf(stderr, "YAFU: json: invalid value: not a decimal number: %s \n", sp);
		return;
	}
	switch (index) {
	case 0:
		if (verbose > 0)
			gmp_printf("YAFU: json: number to be factored = %Zd \n", ZT(numValue));
		ToBeFactored = numValue;
		break;
	case 1:
		if (verbose > 0)
			gmp_printf("YAFU: json: factor = %Zd \n", ZT(numValue));
		factors.push_back(numValue);
		break;

	default:
		abort();   /* WTF? */
	}
}

static void process_object_s(json_value* value, int depth, const char* name,
	const int index)
{
	int length, x;
	if (value == NULL) {
		return;
	}
	length = value->u.object.length;
	for (x = 0; x < length; x++) {
		if (strcmp(value->u.object.values[x].name, name) != 0)
			continue;  /* ignore object unless name matches */
		//print_depth_shift(depth);
		//printf("object[%d].name = %s\n", x, value->u.object.values[x].name);
		process_value_s(value->u.object.values[x].value, depth + 1, name, index);
	}
}

static void process_array_s(json_value* value, int depth, const char* name,
	const int index)
{
	int length, x;
	if (value == NULL) {
		return;
	}
	length = value->u.array.length;
	//printf("array\n");
	for (x = 0; x < length; x++) {
		process_value_s(value->u.array.values[x], depth, name, index);
	}
}

static void process_value_s(json_value* value, int depth, const char* name,
	const int index) {
	if (value == NULL) {
		return;
	}
	switch (value->type) {
	case json_none:
	case json_null:
		break;
	case json_integer:
	case json_double:
	case json_string:
	case json_boolean:
		save_value(value, index, depth);
		break;

	case json_object:
		process_object_s(value, depth + 1, name, index);
		break;
	case json_array:
		process_array_s(value, depth + 1, name, index);
		break;

	default:
		printf("YAFU: invalid json type\n");
		break;
	}
}

/* returns -1 for read error, +1 for parsing error, otherwise
treats each line as a separate json record; parses it and stores
selected fields of last record in a decoded form. */
static int process_file_s(FILE* fp, int* counter, const char* name_list[],
	const int name_list_size, const Znum& num) {
	char* string_p = NULL;
	char buffer[4096] = "";
	json_char* json = NULL;
	json_value* value = NULL;

	*counter = 0;  /* counts number of records in file */
	for (;;) {
		string_p = fgets(buffer, sizeof(buffer), fp); /* get 1 record */
		if (string_p != NULL) {
			//printf("%s \n", buffer);  /* print raw data */
			json = (json_char*)buffer;
			json_value_free(value);  /* avoid memory leakage */
			value = json_parse(json, strlen(buffer));
			if (value == NULL) {
				fprintf(stderr, "Unable to parse data\n");
				return(1);
			}
			(*counter)++;  /* count number of records */
		}
		else {   /* error or end-of-file */
			if (feof(fp)) {
				/* process last record read */
				factors.clear();
				for (int index = 0; index < name_list_size; index++)
					process_value_s(value, 0, name_list[index], index);
				json_value_free(value);  /* avoid memory leakage */
				if (!factors.empty())
					sanityCheck(ToBeFactored, factors, num);
				break;    /* normal exit */
			}
			else {
				perror("YAFU: fgets error");
				return(-1);
			}
		}
	}
	return(0);  /* normal exit */
}

/* returns -1 for read error, +1 for parsing error, 0 for succcess*/
static int process_json_file_main(const Znum& num) {
	FILE* fp;
	char filename[MAX_PATH] = "";
	struct stat filestatus;
	int file_size;
	int counter = 0;
	int rv;
	const char* name_list[] = { "input-decimal", "factors-prime" };
	DWORD rv2 = GetCurrentDirectoryA(MAX_PATH, filename);
	assert(rv2 != 0);

	rv = strcat_s(filename, "\\factor.json");
	if (stat(filename, &filestatus) != 0) {
		fprintf(stderr, "File %s not found\n", filename);
		return -1;
	}
	file_size = filestatus.st_size;
	errno_t errNo = fopen_s(&fp, filename, "rt");
	if (fp == NULL) {
		perror(NULL);
		fprintf(stderr, "Unable to open %s\n", filename);
		fclose(fp);
		return -1;
	}
	rv = process_file_s(fp, &counter, name_list, 2, num); /* read & process file */
	if (counter > 20 || verbose > 0)
		printf("YAFU: factor.json file contains %d records, size = %.1f Kb\n", counter,
		(double)file_size / 1024.0);
	fclose(fp);
	return rv;
}

bool callYafu(const Znum& num, fList& Factors) {
	int rv;
	std::string numStr;
	std::string buffer;
	int fcount = 0;

	if (useOldYafu) {
		return callYafuOld(num, Factors);
	}

	size_t numdigits = mpz_sizeinbase(ZT(num), 10);  // get number of decimal digits in num
	numStr.resize(numdigits + 5);             // resize buffer
	mpz_get_str(&numStr[0], 10, ZT(num));     // convert num to decimal (ascii) digits
	numStr.resize(strlen(&numStr[0]));        // get exact size of string in buffer

	buffer = YafuPath + "\\" + yafuprog;
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
	buffer += '\n';
	if (verbose > 0) {
		std::cout << myTime() << " command is: \n" << buffer;
	}
	delfile("", "nfs.dat"); /* if earlier run leaves this file undeleted
	                               it would cause problems */
	rv = system(buffer.data());             // start YAFU;

	/* get control back when YAFU has finished */
	if (rv == -1) {
#pragma warning(suppress : 4996)
		std::cout << "cannot start YAFU: " << strerror(errno) << '\n';
		return false;
	}
	else if (rv != 0) {
		std::cout << "cannot start YAFU return code = " << rv << '\n';
		return false;
	}

	rv = process_json_file_main(num);
	if (rv != 0)
		return false;
	for (auto f : factors) {
		insertBigFactor(Factors, f);
		fcount++;
	}
	return (fcount > 0);
}