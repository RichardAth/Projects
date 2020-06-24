/*----------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Ben Buhrow. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

Some parts of the code (and also this header), included in this 
distribution have been reused from other sources. In particular I 
have benefitted greatly from the work of Jason Papadopoulos's msieve @ 
www.boo.net/~jasonp, Scott Contini's mpqs implementation, and Tom St. 
Denis Tom's Fast Math library.  Many thanks to their kind donation of 
code to the public domain.
       				   --bbuhrow@gmail.com 7/28/10
----------------------------------------------------------------------*/
#define _CRT_SECURE_NO_WARNINGS
#include "yafu.h"
#include "soe.h"
#include "calc.h"
#include "yafu_string.h"
#include "util.h"
#include "factor.h"
#include "gmp.h"

/* stuff to handle floating point exceptions */
#include <signal.h>

#pragma float_control(precise, on)
#pragma fenv_access(on)
#pragma float_control(except, on)

void __cdecl handle_fpe(int sig, int num);   // Prototype
DWORD DumpStackTrace(EXCEPTION_POINTERS *ep);
void DumpStackTrace2(void);

//#if defined(_MSC_VER)
//	#include <gmp-ecm\config.h>
//#else
//	#include <config.h>
//#endif
/* Use the above if linking against GMP-ECM revision 2344 or earlier, and the
following if linking against 2345+ (all 6.* versions are the old way) */

// note that visual studio builds require a recent SVN, in order to capture
// changes implemented by Brian Gladman to put ECM_VERSION in ecm.h.
#if defined(_MSC_VER)
	#include <ecm.h>
#else
	#include <ecm.h>
#endif

// the number of recognized command line options
//#define NUMOPTIONS 67
// maximum length of command line option strings
#define MAXOPTIONLEN 20

// command line options visible to driver.c
//static const char OptionArray[][MAXOPTIONLEN] = {
//	"B1pm1", "B1pp1", "B1ecm",  "rhomax", "B2pm1",        //   0 - 4
//	"B2pp1", "B2ecm", "qssave", "siqsB",  "siqsTF",       //   5 - 9 
//	"siqsR", "siqsT", "siqsNB", "siqsM",  "logfile",      //  10 -14 
//	"batchfile", "seed", "sigma", "session", "threads",   //  15 -19 
//	"v", "silent", "pfile", "pscreen", "forceDLP",        //  20 -24 
//	"fmtmax", "noopt", "vproc", "noecm", "ggnfs_dir",     //  25 -29
//	"tune_info", "pretest_ratio", "xover", "one", "op",   //  30 -34
//	"of", "ou", "plan", "pretest", "no_expr",             //  35 -39
//	"o", "a", "r", "ggnfsT", "job",                       //  40 -44
//	"ns", "np", "nc", "psearch", "R",                     //  45 -49
//	"pbatch", "ecm_path", "siever", "ncr", "lathreads",   //  50 -54
//	"nc2", "nc3", "p", "work", "nprp",                    //  55 -59
//	"ext_ecm", "testsieve", "nt", "aprcl_p", "aprcl_d",   //  60 -64
//	"filt_bump", "nc1"};                                  //  65 -66


// indication of whether or not an option needs a corresponding argument
// 0 = no argument
// 1 = argument required
// 2 = argument optional
//static const int needsArg[] = {
//	1,1,1,1,1,   //   0 - 4
//	1,1,1,1,1,   //   5 - 9 
//	1,1,1,1,1,   //  10 -14
//	1,1,1,1,1,   //  15 -19
//	0,0,0,0,0,   //  20 -24 
//	1,0,0,0,1,   //  25 -29
//	1,1,1,0,1,   //  30 -34
//	1,1,1,2,0,   //  35 -39
//	1,0,0,1,1,   //  40 -44
//	2,2,0,1,0,   //  45 -49
//	1,1,1,0,1,   //  50 -54
//	0,0,0,1,1,   //  55 -59
//	1,1,1,1,1,   //  60 -64
//	1,0 };       //  65 -66

/* the order of the entries in option must match the case statement in 
applyOpt2 */
static struct {
	char o[MAXOPTIONLEN];  // option value
	int a; // indication of whether or not an option needs a corresponding 
			// argument
			// 0 = no argument
			// 1 = argument required
			// 2 = argument optional
	int type;  // 1 for one integer argument
			   // 2 for 2 integer arguments
			   // 3 for 1 floating point argument
} option[] = {
	"B1pm1",         1, 1,
	"B1pp1",         1, 1, 
	"B1ecm",         1, 1,
	"rhomax",        1, 1,
	"B2pm1",         1, 1,     //   0 - 4
	"B2pp1",         1, 1,
	"B2ecm",         1, 1,
	"qssave",        1, 0,
	"siqsB",         1, 1,
	"siqsTF",        1, 1,     //   5 - 9 
	"siqsR",         1, 1,
	"siqsT",         1, 1,
	"siqsNB",        1, 1,
	"siqsM",         1, 1,
	"logfile",       1, 0,      //  10 -14 
	"batchfile",     1, 0, 
	"seed",          1, 2,
	"sigma",         1, 1,
	"session",       1, 0,
	"threads",       1, 1,  //  15 -19 
	"v",             0, 0,
	"silent",        0, 0,
	"pfile",         0, 0,
	"pscreen",       0, 0,
	"forceDLP",      0, 0,      //  20 -24 
	"fmtmax",        1, 1,
	"noopt",         0, 0,
	"vproc",         0, 0,
	"noecm",         0, 0,
	"ggnfs_dir",     1, 0,     //  25 -29
	"tune_info",     1, 0,
	"pretest_ratio", 1, 3,
	"xover",         1, 3,
	"one",           0, 0,
	"op",            1, 0,  //  30 -34
	"of",            1, 0,
	"ou",            1, 0,
	"plan",          1, 0,
	"pretest",       2, 1,
	"no_expr",       0, 0,    //  35 -39
	"o",             1, 0,
	"a",             0, 0,
	"r",             0, 0,
	"ggnfsT",        1, 1,
	"job",           1, 0,     //  40 -44
	"ns",            2, 2,
	"np",            2, 2,
	"nc",            0, 0,
	"psearch",       1, 0,
	"R",             0, 0,        //  45 -49
	"pbatch",        1, 1,
	"ecm_path",      1, 0,
	"siever",        1, 1,
	"ncr",           0, 0,
	"lathreads",     1, 1,   //  50 -54
	"nc2",           0, 0,
	"nc3",           0, 0,
	"p",             0, 0,
	"work",          1, 3,
	"nprp",          1, 1,         //  55 -59
	"ext_ecm",       1, 1,
	"testsieve",     1, 1,
	"nt",            1, 0,
	"aprcl_p",       1, 1,
	"aprcl_d",       1, 1,    //  60 -64
	"filt_bump",     1, 3,
	"nc1",           0, 0,  //  65 -66
};
static const int NumOptions = sizeof(option) / sizeof(option[0]);

// function to read the .ini file and populate options
static void readINI(fact_obj_t *fobj);
static void apply_tuneinfo(fact_obj_t *fobj, const char *arg);

// functions to populate the global options with default values, and to free
// those which allocate memory
static void set_default_globals(void);
static void free_globals(void);

// function containing system commands to get the computer name
static void get_computer_info(char *idstr);

// function to print the splash screen to file/screen
static void print_splash(int is_cmdline_run, FILE *logfile, char *idstr, fact_obj_t *fobj);

// functions to make a batchfile ready to execute, and to process batchfile lines
static void prepare_batchfile(char *input_exp);
static char * process_batchline(char *input_exp, const char *indup, int *code);
static void finalize_batchline();

// functions to process all incoming arguments
static int process_arguments(int argc, char **argv, char *input_exp, fact_obj_t *fobj);
//static void applyOpt(const char *opt, const char *arg, fact_obj_t *fobj, int j);
static void applyOpt2(const char *opt, const char *arg, fact_obj_t *fobj, int j);
static unsigned process_flags(int argc, char **argv, fact_obj_t *fobj);

/* get time in format hh:mm:ss */
char * myTime(void) {
	static char timestamp[10];   // time in format hh:mm:ss
	const time_t current = time(NULL);
	strftime(timestamp, sizeof(timestamp), "%H:%M:%S", localtime(&current));
	return timestamp;
}

/* the first parameter should have value SIG_FPE
the 2nd parameter can have any of the following values:
FPE_INVALID			0x81   This exception occurs when the result of an operation is 
                           ill-defined, such as (0.0 / 0.0). If trapping is enabled, 
						   the floating-point-invalid condition is signalled. 
						   Otherwise, a quiet NaN is returned. 
FPE_DENORMAL		0x82
FPE_ZERODIVIDE		0x83   This exception occurs when a float is divided by zero. 
                           If trapping is enabled, the  divide-by-zero condition 
						   is signalled. Otherwise, the appropriate infinity is 
						   returned. 
FPE_OVERFLOW		0x84   This exception occurs when the result of an operation 
                           is too large to be represented as a float in its 
						   format. If trapping is enabled, the  floating-point-
						   overflow exception is signalled. Otherwise, the 
						   operation results in the appropriate infinity. 
FPE_UNDERFLOW		0x85   This exception occurs when the result of an operation 
                           is too small to be represented as a normalized float in 
						   its format. If trapping is enabled, the floating-point-
						   underflow condition is signalled. Otherwise, the 
						   operation results in a denormalized float or zero. 
FPE_INEXACT			0x86   This exception occurs when the result of an operation 
                           is not exact, i.e. the result was rounded. If trapping 
						   is enabled, the floating-point-inexact condition is 
						   signalled. Otherwise, the rounded result is returned. 
FPE_UNEMULATED		0x87
FPE_SQRTNEG			0x88
FPE_STACKOVERFLOW	0x8a
FPE_STACKUNDERFLOW	0x8b
FPE_EXPLICITGEN     0x8c // raise(SIGFPE);
*/
void handle_fpe(int sig, int num) {
	unsigned int control_word;
	errno_t err = _controlfp_s(&control_word, _EM_INEXACT | _EM_UNDERFLOW,
		MCW_EM);  /* allow hardware FP exceptions exept inexact & underflow*/
	if (err) {
		printf_s("could not set FP control word\n");
		exit(1);
	}
	/* reset SIGFPE handler  */
	if (signal(SIGFPE, (void(__cdecl *)(int))handle_fpe) == SIG_ERR) {
		printf_s("could not install handler on SIGFPE\n");
		exit(1);
	}
	_fpreset();
	if (num == _FPE_INEXACT)
		return;
	/* strictly speaking, should not call I/O routines here, but it seems to work OK */
	fprintf_s(stderr, "\nreceived FPE signal %d %#x; continue \n", sig, num);
	DumpStackTrace2();
	if (num == _FPE_UNDERFLOW) {
		return;  // for underflow, continue after recording it
	}

#ifdef _DEBUG
	__debugbreak();     // try to enter debuggger (to look at call stack)
	Beep(1200, 1000);   // beep at 1200 Hz for 1 second
	system("PAUSE");    // press any key to continue
#endif

	return;
}

/* if an error occurs exit with code 1 or -1 */
int main(int argc, char *argv[])
{
	uint32 insize = GSTR_MAXSIZE;
	char *input_exp, *ptr;
	str_t str;
	mpz_t tmp;
	int nooutput, slog, is_cmdline_run=0;
	ptrdiff_t offset;
	FILE *logfile;
	fact_obj_t *fobj;
	char indup[GSTR_MAXSIZE];

	//the input expression
	input_exp = (char *)malloc(GSTR_MAXSIZE*sizeof(char));
	mallocCheck(input_exp)
	/*indup = (char *)malloc(GSTR_MAXSIZE*sizeof(char));
	mallocCheck(indup)*/
	//strcpy_s(input_exp, sizeof(input_exp), "");
	input_exp[0] = '\0';

	sInit(&str);

	//set defaults for various things
	set_default_globals();

	/* detect floating point errors */
	if (signal(SIGFPE, (void(__cdecl *)(int))handle_fpe) == SIG_ERR) {
		printf_s("could not install handler on SIGFPE\n");
		return -1;
	}
	unsigned int control_word;
	errno_t err = _controlfp_s(&control_word, _EM_INEXACT | _EM_UNDERFLOW, MCW_EM);
	/* unblock hardware FP exceptions exept inexact */
	if (err) {
		printf_s("could not set FP control word\n");
		return -1;
	}

	// a factorization object that gets passed around to any factorization routine
	// called out in the input expression.  if no factorization routine is specified,
	// this is not used.  the factor object must be initialized prior to parsing
	// input options in order to make those options stick.
	fobj = (fact_obj_t *)malloc(sizeof(fact_obj_t));
	mallocCheck(fobj)
	init_factobj(fobj);

	//get the computer name, cache sizes, etc.  store in globals
	get_computer_info(CPU_ID_STR);	

	//now check for an .ini file, which will override these defaults
	//command line arguments will override the .ini file
	readINI(fobj);	

	//check/process input arguments
	is_cmdline_run = process_arguments(argc, argv, input_exp, fobj);
	
	// now that we've processed arguments, spit out vproc info if requested
#ifndef __APPLE__
	if (VERBOSE_PROC_INFO)
		extended_cpuid(CPU_ID_STR, &CLSIZE, &HAS_SSE41, VERBOSE_PROC_INFO);
#endif

	// get the batchfile ready, if requested
	if (USEBATCHFILE)
	{
		prepare_batchfile(input_exp);		
		
		//batchfile jobs are command line in nature
		is_cmdline_run = 1;

		//remember the input expression
		strcpy_s(indup, sizeof(indup), input_exp);
	}

	//never run silently when run interactively, else the results of
	//calculations will never be displayed.
	if (!is_cmdline_run && Vflag < 0)
		Vflag = 0;

	//session log
	logfile = fopen(sessionname,"a");
	if (logfile == NULL)
	{
		printf("fopen error: %s\n", strerror(errno));
		printf("couldn't open %s for appending\n",sessionname);
		slog = 0;
	}
	else
		slog = 1;	
		
	// print the splash screen, to the logfile and depending on options, to the screen
	print_splash(is_cmdline_run, logfile, CPU_ID_STR, fobj);

	// bigint used in this routine
	mpz_init(tmp);
	
	//start the calculator
	//right now this just allocates room for user variables
	calc_init();				
		
	if (USERSEED)
	{
		logprint(logfile,"User random seed:  %u\n\n",g_rand.low);
	}
	else
	{		
		logprint(logfile,"New random seeds: %u, %u\n\n",g_rand.hi,g_rand.low);
	}
	fflush(logfile);

	//printf("WARNING: constant seed is set\n");
	//g_rand.hi = 123;
	//g_rand.low = 123;
	srand(g_rand.low);
	gmp_randinit_default(gmp_randstate);
	gmp_randseed_ui(gmp_randstate, g_rand.low);

#if BITS_PER_DIGIT == 64
	LCGSTATE = (uint64)g_rand.hi << 32 | (uint64)g_rand.low;
#else
	LCGSTATE = g_rand.low;
#endif	


	//command line
	while (1)
	{		
		reset_factobj(fobj);		

		//handle a batch file, if passed in.
		if (USEBATCHFILE)
		{
			int code;
			input_exp = process_batchline(input_exp, indup, &code);
			if (code == 1)
			{
				finalize_batchline();
				break;
			}
			else if (code == 2)
				continue;
		}
		else if (!is_cmdline_run)
		{
			int c = fgetc(stdin);
			if (c == EOF)
				break; // ^D quits yafu (but, for reasons I've not investigated, doesn't print the proper newline)
			ungetc(c, stdin);

			// get command from user
			fgets(input_exp,GSTR_MAXSIZE,stdin);
			while (1)
			{
				if (input_exp[strlen(input_exp) - 1] == '\r' || input_exp[strlen(input_exp) - 1] == '\n')
				{
					//replace with a null char and continue
					printf("\n");
					fflush(stdout);
					input_exp[strlen(input_exp) - 1] = '\0';
					break;
				}
				else
				{
					//last char is not a carriage return means
					//the input is longer than allocated.
					//reallocate and get another chunk
					insize += GSTR_MAXSIZE;
					input_exp = (char *)realloc(input_exp,insize*sizeof(char));
					if (input_exp == NULL)
					{
						printf("couldn't reallocate string when parsing\n");
						exit(-1);
					}
					fgets(input_exp+strlen(input_exp),GSTR_MAXSIZE,stdin);
				}
			}	
		}
		else
		{
			// input expression already read in.  nothing to do here.

		}
		
		//search for substring help in input
		ptr = strstr(input_exp,"help");
		if (ptr != NULL)
			helpfunc(input_exp);
		else if ((strcmp(input_exp,"quit") == 0) || (strcmp(input_exp,"exit") == 0))
			break;
		else
		{
			logprint(logfile,"Processing expression: %s\n\n",input_exp);
			toStr(input_exp,&str);

			preprocess(&str);
			strcpy_s(input_exp, GSTR_MAXSIZE * sizeof(char), str.s);

			//detect an '=' operation, and verify the destination of the = is valid
			//pass on everything to the right of the = to the calculator
			if ((ptr = strchr(str.s,'=')) != NULL)
			{
				*ptr = '\0';
				if (invalid_dest(str.s))
				{
					printf("invalid destination %s\n",str.s);
					offset = ptr-str.s+1;
					sFree(&str);
					sInit(&str);
					memcpy(str.s,input_exp+offset,(GSTR_MAXSIZE-offset) * sizeof(char));
					strcpy_s(input_exp, GSTR_MAXSIZE * sizeof(char), "ans");
					str.nchars = (int)strlen(str.s) + 1;
				}
				else
				{
					offset = ptr-str.s+1;
					sFree(&str);
					sInit(&str);
					memcpy(str.s,input_exp+offset,(GSTR_MAXSIZE-offset) * sizeof(char));

					input_exp[offset-1] = '\0';
					str.nchars = (int)strlen(str.s) + 1;
				}
			}
			else
				strcpy_s(input_exp, GSTR_MAXSIZE * sizeof(char), "ans");

			//look for a trailing semicolon
			if (str.s[str.nchars-2] == ';')
			{
				nooutput = 1;
				str.s[str.nchars-2] = '\0';
				str.nchars--;
			}
			else
				nooutput = 0;

			// new factorization
			fobj->refactor_depth = 0;

			if (!calc(&str, fobj))
			{
				if (strcmp(str.s,"") != 0)
				{
					clock_t start, stop;
					double t;

					start = clock();
					mpz_set_str(tmp, str.s, 0);
					stop = clock();
					t = (double)(stop - start)/(double)CLOCKS_PER_SEC;
					//printf("str2hexz in = %6.4f seconds.\n",t);

					if (set_uvar(input_exp,tmp))
						new_uvar(input_exp,tmp);
					if (nooutput == 0)
					{
						if (OBASE == DEC)
						{
							int sz = (int)mpz_sizeinbase(tmp, 10) + 10;
							if (gstr1.alloc < sz)
							{
								gstr1.s = (char *)realloc(gstr1.s, sz * sizeof(char));
								gstr1.alloc = sz;
							}
							if (Vflag >= 0)
								printf("\n%s = %s\n\n",input_exp, mpz_get_str(gstr1.s, 10, tmp));
						}
						else if (OBASE == HEX)
						{
							ptrdiff_t sz = mpz_sizeinbase(tmp, 16) + 10;
							if (gstr1.alloc < sz)
							{
								gstr1.s = (char *)realloc(gstr1.s, sz * sizeof(char));
								gstr1.alloc = (int)sz;
							}
							if (Vflag >= 0)
								printf("\n%s = %s\n\n",input_exp, mpz_get_str(gstr1.s, 16, tmp));
						}
					}
				}
			}						
		}

#if defined(WIN32)
		fflush(stdin);	//important!  otherwise scanf will process printf's output
		
#else
		if (!is_cmdline_run)
		{
			fflush(stdin);	//important!  otherwise scanf will process printf's output
			fflush(stdout);
		}
#endif

		input_exp = (char *)realloc(input_exp,GSTR_MAXSIZE*sizeof(char));
		if (input_exp == NULL)
		{
			printf("couldn't reallocate string during cleanup\n");
			exit(-1);
		}
		input_exp[0] = '\0';

		if (is_cmdline_run)
		{
			if (USEBATCHFILE)
			{
				// the line from the batchfile finished.  make the temporary file
				// created in processs_batchline the new batchfile, with the line
				// we just finished removed.
				finalize_batchline();
			}
			else
				break;
		}
		else
			printf(">> ");

	}

	if (slog)
		fclose(logfile);

	calc_finalize();
	free_globals();
	sFree(&str);
	free(input_exp);
	//free(indup);	
	mpz_clear(tmp);
	free_factobj(fobj);
	free(fobj);

	return 0;
}

//check if optbuf contains a valid option. If so set j to index for option
// if not valid return -1;
static int checkopt(char optbuf[]) {
	if (optbuf == NULL)
		return -1;
	for (int j = 0; j < NumOptions; j++)
	{
		if (strcmp(option[j].o, optbuf) == 0)
			return j;
	}
	return - 1;   // invalid option
}

// function to read the .ini file and populate options
static void readINI(fact_obj_t *fobj)
{
	FILE *doc;
	char str[1024];
	char *key;
	char *value;
	ptrdiff_t len;

	doc = fopen("yafu.ini","r");

	if (doc == NULL) {
		perror("could not open yafu.ini");
		return;
	}

	while (fgets(str, 1024, doc) != NULL)
	{
#ifdef _DEBUG
		printf("yafu.ini; %.160s", str);
#else
		//if first character is a % sign, skip this line
		if (strncmp(str, "%print", 6) == 0) {
			printf("yafu.ini; %.160s", str+1);
			continue;
		}
#endif
		if (str[0] == '%')
			continue;
		if (str[0] == '\0')
			continue;  // ignore 0-length string

		//if last character of line is newline, remove it
		do {
			len = strlen(str);
			if (str[len - 1] == '\n')         // newline
				str[len - 1] = '\0';
			else if (str[len - 1] == '\r')    // carriage return
				str[len - 1] = '\0';
			else
				break;
		} while (len > 0);

		//read keyword by looking for an equal sign
		// = is replaced by a null character, so key is null-terminated
		key = strtok(str,"=");
		int ix = checkopt(key);
		if (key == NULL || ix < 0)
		{
			printf("Invalid line in yafu.ini, use Keyword=Value pairs \n"
				"%s\n"
				"See docfile.txt for valid keywords", str);
			continue;  // don't process this invalid line.
		}

		//read value
		value = strtok((char *)0, "=");

		//apply the option... same routine command line options use
		applyOpt2(key, value, fobj, ix);
	}

	fclose(doc);

#ifdef _DEBUG
	/* cant' use Vflag because it's not set up yet*/
		printf("yafu.ini processed\n");
#endif
	return;
}

//search the docfile for the right entry
//just search for the heading, and print everything until
//the next heading is found
static void helpfunc(const char *s)
{
	FILE *doc;
	const char *func;
	char str[1024];
	int printtopic = 0;
	int j;

	j=4;
	if (s[j] == '\0')
		j = 0;
	else
		while (isspace(s[j])) j++;		//skip white space
	func = s + j;

	//func now points to a string with the desired help topic
	//open the doc file and search matching topics
	doc = fopen("docfile.txt","r");
	if (doc == NULL)
	{
		printf("fopen error: %s\n", strerror(errno));
		printf("documentation file not found\n");
		return;
	}

	printf("searching for help on '%s'\n",func);

	while (!feof(doc))
	{
		//read a line
		fgets(str, 1024, doc);
		//is this a header?
		//printf("(%d) %s",strlen(str),str);
		if ((str[0] == '[') && (str[strlen(str) - 2] == ']'))
		{
			//printf("in printtopic if\n");
			if (printtopic == 1)
				break;
			printtopic = 0;
			//does it match our topic?
			str[strlen(str) - 2] = '\0';
			if (strstr(func,str+1) != NULL)
				printtopic = 1;
		}
		else
		{
			if (printtopic)
				printf("%s", str);
		}
	}
	fclose(doc);
	return;
}

//return 1 if invalid, 0 otherwise
static int invalid_dest(char *dest)
{
	int i;

	if (getFunc(dest,&i) >= 0)
		return 1;	//is a function name

	//global vars are ok
	if (strcmp(dest,"POLLARD_STG1_MAX") == 0) {
		return 0;}
	else if (strcmp(dest,"POLLARD_STG2_MAX") == 0) {
		return 0;}
	else if (strcmp(dest,"WILL_STG1_MAX") == 0) {
		return 0;}
	else if (strcmp(dest,"WILL_STG2_MAX") == 0) {
		return 0;}
	else if (strcmp(dest,"ECM_STG1_MAX") == 0) {
		return 0;}
	else if (strcmp(dest,"ECM_STG2_MAX") == 0) {
		return 0;}
	else if (strcmp(dest,"BRENT_MAX_IT") == 0) {
		return 0;}
	else if (strcmp(dest,"IBASE") == 0) {
		return 0;}
	else if (strcmp(dest,"OBASE") == 0) {
		return 0;}
	else if (strcmp(dest,"QS_DUMP_CUTOFF") == 0) {
		return 0;}
	else if (strcmp(dest,"NUM_WITNESSES") == 0) {
		return 0;}
	else if (strcmp(dest,"LOGFLAG") == 0) {
		return 0;}
	else if (strcmp(dest,"VFLAG") == 0) {
		return 0;}
	else if (strcmp(dest,"PRIMES_TO_FILE") == 0) {
		return 0;}
	else if (strcmp(dest,"PRIMES_TO_SCREEN") == 0) {
		return 0;}

	//check starting char not lower case letter or _ or `
	if (dest[0] < '_' || dest[0] > 'z') 
		return 1; // not lower case letter or _ or `

	return 0;
}

//check that num consists of only numeric or alphanumeric characters. 
// return 0 if valid, 1 if invalid
//int invalid_num(char *num) {
//	int i=0;
//	int nchars = strlen(num);
//	
//	if (nchars == 0) return 1;
//
//	if (num[0] == '-')
//		i++;
//	
//	//check for 0x, 0d, 0b, or 0o.  nchars must be > 3-i in this case
//	if (num[i] == '0' && num[i+1] == 'x' && ((nchars-i) > 2))
//	{
//		//num is hex, and can have lower or upper case alpha characters
//		i += 2;
//		for (;i<nchars;i++)
//		{ 
//			if (num[i] > 102)	//102 == f
//				return 1;
//			else if (num[i] < 48) 
//				return 1;
//			else if (num[i] > 57 && num[i] < 65)
//				return 1;
//			else if (num[i] > 70 && num[i] < 97)	//97 == a
//				return 1;
//		}
//	}
//	else if (num[i] == '0' && num[i+1] == 'd' && ((nchars-i) > 2))
//	{
//		//num is dec, and can have only digits 0 - 9
//		i += 2;
//		for (;i<nchars;i++)
//		{ 
//			if (num[i] < '0' || num[i] > '9') 
//				return 1;
//		}
//	}
//	else if (num[i] == '0' && num[i+1] == 'b' && ((nchars-i) > 2))
//	{
//		//num is bin, and can have only digits 0 - 1
//		i += 2;
//		for (;i<nchars;i++)
//		{ 
//			if (num[i] < '0' || num[i] > '1') 
//				return 1;
//		}
//	}
//	else if (num[i] == '0' && num[i+1] == 'o' && ((nchars-i) > 2))
//	{
//		//num is oct, and can have only digits 0 - 7
//		i += 2;
//		for (;i<nchars;i++)
//		{ 
//			if (num[i] < '0' || num[i] > '7')
//				return 1;
//		}
//	}
//	else
//	{
//		//no base designator, go by IBASE
//		if (IBASE == HEX)
//		{
//			//num is hex, and can have only upper case alpha characters
//			for (;i<nchars;i++)
//			{ 
//				if (num[i] < 48) 
//					return 1;
//				else if (num[i] > 57 && num[i] < 65)
//					return 1;
//				else if (num[i] > 70)	//70 == F
//					return 1;
//			}
//		}
//		else if (IBASE == DEC)
//		{
//			//num is dec, and can have only digits 0 - 9
//			for (;i<nchars;i++)
//			{ 
//				if (num[i] < 48 || num[i] > 57) 
//					return 1;
//			}
//		}
//		else if (IBASE == BIN)
//		{
//			//num is bin, and can have only digits 0 - 1
//			for (;i<nchars;i++)
//			{ 
//				if (num[i] < 48 || num[i] > 49) 
//					return 1;
//			}
//		}
//		else if (IBASE == OCT)
//		{
//			//num is oct, and can have only digits 0 - 7
//			for (;i<nchars;i++)
//			{ 
//				if (num[i] < 48 || num[i] > 55) 
//					return 1;
//			}
//		}
//	}
//
//	return 0;
//}

// function to make a batchfile ready to execute 
static void prepare_batchfile(char *input_exp)
{
	char *ptr;
	
	//look for @ symbol in input expression
	ptr = strchr(input_exp,'@');
	if (ptr == NULL)
	{
		printf("missing variable indicator (@) in input expression\n");
		printf("ignoring any input expression: interpreting batchfile lines as input expressions\n");
		sprintf_s(input_exp, GSTR_MAXSIZE * sizeof(char), "@");
	}

	return;
}

// functions to process all incoming arguments
static int process_arguments(int argc, char **argv, char *input_exp, fact_obj_t *fobj)
{
	int is_cmdline_run=0;
	FILE *in = stdin;

	//now check for and handle any incoming arguments, whatever
	//their source.  
	if (argc > 1)
	{
		//process arguments

		if (argv[1][0] == '-')
		{
			//then there are just flags, no expression.  start up normally
			//after processing the flags
			process_flags(argc-1, &argv[1], fobj);
		}
		else
		{
			//assume its a command.  execute once, after processing any
			//flags.  an invalid command will be handled by calc.
			is_cmdline_run = 1;
			strcpy_s(input_exp, GSTR_MAXSIZE * sizeof(char), argv[1]);
			if (argc > 2)
				process_flags(argc-2, &argv[2], fobj);
		}		
	}
	else
	{
		//else, need to check for input from a redirect or pipe.
		//this is done differently on unix like systems vs. windows
#if defined(WIN32)	//not complete, but ok for now
		fseek(in,-1,SEEK_END);
		if (ftell(in) >= 0)
		{
			rewind(in);
			fgets(input_exp,1024,in);
			is_cmdline_run = 1;
		}

#else
		if (isatty(fileno(in)) == 0)
		{			
			fgets(input_exp,1024,in);
			is_cmdline_run = 1;
		}
#endif
	}

	return is_cmdline_run;

}

// function to print the splash screen to file/screen
static void print_splash(int is_cmdline_run, FILE *logfile, char *idstr, 
	fact_obj_t *fobj)
{
	if (Vflag >= 0)
		printf("\n\n");

	if (Vflag > 0 || !is_cmdline_run)
	{	
		logprint(NULL,"System/Build Info: \n");
	}
	logprint(logfile,"System/Build Info: \n");
	fflush(stdout);

	if (Vflag > 0 || !is_cmdline_run)
#ifdef _MSC_MPIR_VERSION
		printf("Using GMP-ECM %s, Powered by MPIR %s\n", ECM_VERSION,
_MSC_MPIR_VERSION);
		fprintf(logfile,"Using GMP-ECM %s, Powered by MPIR %s\n", ECM_VERSION,
_MSC_MPIR_VERSION);
#else
	#ifdef ECM_VERSION
		printf("Using GMP-ECM %s, Powered by GMP %d.%d.%d\n", ECM_VERSION, 
			__GNU_MP_VERSION,__GNU_MP_VERSION_MINOR,__GNU_MP_VERSION_PATCHLEVEL);
		fprintf(logfile,"Using GMP-ECM %s, Powered by GMP %d.%d.%d\n", ECM_VERSION,
		__GNU_MP_VERSION,__GNU_MP_VERSION_MINOR,__GNU_MP_VERSION_PATCHLEVEL);
	#else
		printf("Using GMP-ECM, Powered by GMP\n");
		fprintf(logfile,"Using GMP-ECM, Powered by GMP\n");
	#endif

#endif

	fflush(stdout);

	fprintf(logfile,"cached %u primes. pmax = %u\n",szSOEp,spSOEprimes[szSOEp-1]);
	fprintf(logfile,"detected %s\ndetected L1 = %d bytes, L2 = %d bytes, CL = %d bytes\n",
		idstr,L1CACHE,L2CACHE,CLSIZE);
	fprintf(logfile,"measured cpu frequency ~= %f\n",
		MEAS_CPU_FREQUENCY);
	fprintf(logfile,"using %u random witnesses for Rabin-Miller PRP checks\n\n",
			NUM_WITNESSES);

	fflush(logfile);

	if (Vflag > 0 || !is_cmdline_run)
	{		
		printf("detected %s\ndetected L1 = %d bytes, L2 = %d bytes, CL = %d bytes\n",
			idstr,L1CACHE,L2CACHE,CLSIZE);
		printf("measured cpu frequency ~= %f\n",
			MEAS_CPU_FREQUENCY);
		printf("using %u random witnesses for Rabin-Miller PRP checks\n\n",
			NUM_WITNESSES);

		printf("===============================================================\n");
		printf("======= Welcome to YAFU (Yet Another Factoring Utility) =======\n");
		printf("======= version:          " VERSION_STRING "                        =======\n");
		printf("=======             bbuhrow@gmail.com                   =======\n");
		printf("=======     Type help at any time, or quit to quit      =======\n");
		printf("===============================================================\n");
		printf("cached %u primes. pmax = %u\n\n",szSOEp,spSOEprimes[szSOEp-1]);
		char * path = _getcwd(NULL, 0);   // get current working directory
		printf("Current directory is: %s\n", path);
		free(path);
		printf("QS to NFS crossover at %1.0f digits\n",
			fobj->autofact_obj.qs_gnfs_xover);
		printf("%s compiled on %s \n", __FILE__, __DATE__);
		printf("\n>> ");

	}

	return;
}

// function containing system commands to get the computer name
static void get_computer_info(char *idstr)
{
	//int ret;

	//figure out cpu freq in order to scale qs time estimations
	//0.1 seconds won't be very accurate, but hopefully close
	//enough for the rough scaling we'll be doing anyway.
    MEAS_CPU_FREQUENCY = measure_processor_speed() / 1.0e5;
	
#ifdef __APPLE__
	// something in extended cpuid causes a segfault on mac builds.
	// just disable it for now - this information is not critical for
	// program operation.
	strcpy(idstr, "N/A");
	CLSIZE = 0;
	L1CACHE = DEFAULT_L1_CACHE_SIZE;
	L2CACHE = DEFAULT_L2_CACHE_SIZE;
	HAS_SSE41 = 0;

#else
	//read cache sizes
	yafu_get_cache_sizes(&L1CACHE, &L2CACHE);

	// run an extended cpuid command to get the cache line size, and
	// optionally print a bunch of info to the screen
	extended_cpuid(idstr, &CLSIZE, &HAS_SSE41, VERBOSE_PROC_INFO);

	#if defined(WIN32)

		sysname_sz = MAX_COMPUTERNAME_LENGTH + 1;
		GetComputerName(sysname,&sysname_sz);
	
	#else

		ret = gethostname(sysname,sizeof(sysname) / sizeof(*sysname));
		sysname[(sizeof(sysname)-1)/sizeof(*sysname)] = 0;	// null terminate
		if (ret != 0)
		{
			printf("error occured when getting host name\n");
			strcpy(sysname, "N/A");
		}
		sysname_sz = strlen(sysname);
	
	#endif

#endif
	return;
}

static void yafu_set_idle_priority(void) {

#if defined(WIN32) || defined(_WIN64)
	SetPriorityClass(GetCurrentProcess(),
			IDLE_PRIORITY_CLASS);
#else
	nice(100);
#endif
}

// function to populate the global options with default values
static void set_default_globals(void)
{
	uint64 limit, i;
	uint32 seed_p[6542], num_sp;  // there are 6542 primes < 65537
	
	Vflag = 0;
	VERBOSE_PROC_INFO = 0;
	LOGFLAG = 1;

	NUM_WITNESSES = 20;
	
	PRIMES_TO_FILE = 0;
	PRIMES_TO_SCREEN = 0;
	GLOBAL_OFFSET = 0;
	
	USEBATCHFILE = 0;
	USERSEED = 0;
	THREADS = 1;
	LATHREADS = 0;

	strcpy_s(sessionname, sizeof(sessionname), "session.log");	

	// initial limit of cache of primes.
	szSOEp = 1000000;	

	//set some useful globals
	zInit(&zZero);
	zInit(&zOne);
	zInit(&zTwo);
	zInit(&zThree);
	zInit(&zFive);
	zOne.val[0] = 1;
	zTwo.val[0] = 2;
	zThree.val[0] = 3;
	zFive.val[0] = 5;

	//global strings, used mostly for logprint stuff
	sInit(&gstr1);
	sInit(&gstr2);
	sInit(&gstr3);

	//global i/o base
	IBASE = DEC;
	OBASE = DEC;

	//find, and hold globally, primes less than some N
	//bootstrap the process by finding some initial sieve primes.
	//if the requested offset+range is large we may need to find more - 
	//we can use these primes to accomplish that.
	num_sp = tiny_soe(65537, seed_p);  // there are 6542 primes < 65537
	PRIMES = GetPRIMESRange(seed_p, num_sp, NULL, 0, szSOEp, &limit);

	//save a batch of sieve primes too.
	spSOEprimes = (uint32 *)malloc((size_t) (limit * sizeof(uint32)));
	mallocCheck(spSOEprimes)
	for (i=0; i<limit; i++)
		spSOEprimes[i] = (uint32)PRIMES[i];

	szSOEp = (uint32)limit;
	NUM_P = limit;
	P_MIN = 0; 
	P_MAX = PRIMES[NUM_P-1];

	// random seeds
	get_random_seeds(&g_rand);	

	return;
}

// function to free the global options which allocate memory
static void free_globals(void)
{
	zFree(&zZero);
	zFree(&zOne);
	zFree(&zTwo);
	zFree(&zThree);
	zFree(&zFive);
	free(spSOEprimes);
	free(PRIMES);
	sFree(&gstr1);
	sFree(&gstr2);
	sFree(&gstr3);

	return;
}

static void finalize_batchline()
{
	rename(batchfilename,"_bkup");
	rename("__tmpbatchfile",batchfilename);
	remove("_bkup");
	remove("__tmpbatchfile");

	return;
}

// function to process batchfile lines
static char * process_batchline(char *input_exp, const char *indup, int *code)
{
	int nChars, j, i;
	char *line, tmpline[GSTR_MAXSIZE], *ptr, *ptr2;
	FILE *batchfile, *tmpfile;

	//try to open the file
	batchfile = fopen(batchfilename,"r");	

	if (batchfile == NULL)
	{
		printf("fopen error: %s\n", strerror(errno));
		printf("couldn't open %s for reading\n",batchfilename);
		exit(-1);
	}	

	//load the next line of the batch file and get the expression
	//ready for processing
	line = (char *)malloc(GSTR_MAXSIZE * sizeof(char));
	mallocCheck(line)
	strcpy_s(line, GSTR_MAXSIZE * sizeof(char), "");
	strcpy_s(input_exp, GSTR_MAXSIZE * sizeof(char), "");

	// read a line - skipping blank lines
	do
	{
		while (1)
		{
			ptr = fgets(tmpline,GSTR_MAXSIZE,batchfile);
			strcpy_s(line + strlen(line), 
				GSTR_MAXSIZE * sizeof(char)- strlen(line),
				tmpline);
			
			// stop if we didn't read anything
			if (feof(batchfile))
			{
				printf("eof; done processing batchfile\n");
				fclose(batchfile);
				*code = 1;
				free(line);
				return input_exp;
			}

			if (ptr == NULL)
			{
				printf("fgets returned null; done processing batchfile\n");		
				fclose(batchfile);
				*code = 1;
				free(line);
				return input_exp;
			}

			// if we got the end of the line, stop reading
			if ((line[strlen(line)-1] == 0xa) ||
				(line[strlen(line)-1] == 0xd))
				break;

			// else reallocate the buffer and get some more
			line = (char *)realloc(line, (strlen(line) + GSTR_MAXSIZE) * sizeof(char));
		} 

		// remove LF an CRs from line
		nChars = 0;
		for (j=0; j<strlen(line); j++)
		{
			switch (line[j])
			{
			case 13:
			case 10:
				break;  // don't copy cr or lf
			default:
				line[nChars++] = line[j];
				break;
			}
		}
		line[nChars++] = '\0';

	} while (strlen(line) == 0);	

	// copy everything in the file after the line we just read to
	// a temporary file.  if the expression we just read finishes, 
	// the temporary file will become the batch file (effectively 
	// eliminating the expression from the batch file).
	tmpfile = fopen("__tmpbatchfile", "w");
	
	if (tmpfile == NULL)
	{
		printf("fopen error: %s\n", strerror(errno));
		printf("couldn't open __tmpbatchfile for reading\n");
		exit(-1);
	}	

	while (!feof(batchfile))
	{
		ptr2 = fgets(tmpline,GSTR_MAXSIZE,batchfile);
		if (ptr2 == NULL)
			break;

		if (strlen(tmpline) == 0)
			continue;

		fputs(tmpline,tmpfile);
	}
	fclose(tmpfile);

	// close the batchfile
	fclose(batchfile);		

	//ignore blank lines
	if (strlen(line) == 0)
	{
		*code = 2;
		free(line);
		return input_exp;
	}

	//ignore comment lines
	if (((line[0] == '/') && (line[1] == '/')) || (line[0] == '%'))
	{
		*code = 2;
		free(line);
		return input_exp;
	}

	//substitute the batchfile line into the '@' symbol in the input expression
	nChars = 0;
	if ((strlen(indup) + strlen(line)) >= GSTR_MAXSIZE)
		input_exp = (char *)realloc(input_exp, strlen(indup) + strlen(line) + 2);

	for (i=0; i<strlen(indup); i++)
	{
		if (indup[i] == '@')
		{
			for (j=0; j<strlen(line); j++)
				input_exp[nChars++] = line[j];
		}
		else				
			input_exp[nChars++] = indup[i];
	}
	input_exp[nChars++] = '\0';

	printf("=== Starting work on batchfile expression ===\n");
	printf("%s\n",input_exp);
	printf("=============================================\n");
	fflush(stdout);

	free(line);
	*code = 0;
	return input_exp;;
}

// functions to process all incoming arguments
static unsigned process_flags(int argc, char **argv, fact_obj_t *fobj)
{
    int ch = 0, ix, j;
	char optbuf[MAXOPTIONLEN];
	char argbuf[80];

    //argument loop
	ix = 0;
	while (ix < argc)
	{
		//read in the option
		ch = argv[ix][0];
		if (ch != '-')
		{
			printf("no switch detected\n");
			exit(1);
		}
		auto ec =strcpy_s(optbuf, sizeof(optbuf), argv[ix]);
		if (ec != 0) {
			printf("invalid option %s\n", argv[ix]);
			exit(1);  // argv too long??
		}

		//check if it's valid. If so set j to index for option
		j = checkopt(optbuf + 1);

		if (j < 0) {
			printf("invalid option %s\n", optbuf);
			exit(1);
		}

		//check to see if this option requires an argument
		if (option[j].a == 1) {
			ix++;
			if (ix == argc) 	{
				printf("argument expected for %s\n",optbuf);
				exit(1);
			}
			ec = strcpy_s(argbuf, sizeof(argbuf), argv[ix]);
			if (ec != 0) {
				printf("invalid option argument %s\n", argv[ix]);
				exit(1);  // argv too long??
			}

			//now apply -option argument
			//printf("applying option %s with argument %s\n",optbuf+1,argbuf);
			applyOpt2(optbuf+1, argbuf, fobj, j);
		}
		else if (option[j].a == 2)
		{	// Argument is optional. Check to see if an argument was supplied
			if (((ix+1) == argc) || argv[ix+1][0] == '-')
			{
				// no option supplied.  use default option
				applyOpt2(optbuf+1, NULL, fobj, j);
			}
			else {  // an option was supplied, pass it on
				ix++;
				ec = strcpy_s(argbuf, sizeof(argbuf), argv[ix]);
				if (ec != 0) {
					printf("invalid option argument %s\n", argv[ix]);
					exit(1);  // argv too long??
				}
				//now apply -option argument
				applyOpt2(optbuf+1, argbuf, fobj, j);
			}

		}
		else {			//apply -option without argument
			applyOpt2(optbuf+1, NULL, fobj, j);
		}
		ix++;
	}

    return 1;
}

/* check whether arg is numeric */
static int isnumeric(const char arg[]) {
	for (int i = 0; i < (int)strlen(arg); i++)
	{
		if (!isdigit(arg[i]))
			return FALSE;
	}
	return TRUE;
}

/* functions to process all incoming arguments 
arguments can come from the command line or the yafu.ini file 
opt is the option as a Null-terminated ascii string, 
arg is the argument as a Null-terminated ascii string, or NULL if no argument
If an error occurs exit with code 1
*/
//static void applyOpt(const char *opt, const char *arg, fact_obj_t *fobj,
//	int ix)
//{
//	char **ptr;
//	//int i;            // not used any more
//	//z tmp;            // not used
//
//	//zInit(&tmp);
//
//	ptr = NULL;
//	if (strcmp(opt,OptionArray[0]) == 0) // B1pm1
//	{
//		if (!isnumeric(arg)) {  //argument should be all numeric
//			printf("expected numeric input for option %s\n",opt);
//			exit(1);
//		}
//
//		fobj->pm1_obj.B1 = strtoul(arg,ptr,10);
//		if (fobj->pm1_obj.stg2_is_default)
//		{
//			//stg2 hasn't been specified yet, so set it to the default value
//			fobj->pm1_obj.B2 = 100 * (uint64)fobj->pm1_obj.B1;
//			//else, we have already specified a B2, so don't overwrite it with
//			//the default
//		}
//	}
//	else if (strcmp(opt,OptionArray[1]) == 0)  // B1pp1
//	{
//		if (!isnumeric(arg)) {  //argument should be all numeric
//			printf("expected numeric input for option %s\n", opt);
//			exit(1);
//		}
//
//		fobj->pp1_obj.B1 = strtoul(arg,ptr,10);
//		if (fobj->pp1_obj.stg2_is_default)
//		{
//			//stg2 hasn't been specified yet, so set it to the default value
//			fobj->pp1_obj.B2 = 50 * (uint64)fobj->pp1_obj.B1;
//			//else, we have already specified a B2, so don't overwrite it with
//			//the default
//		}
//	}
//	else if (strcmp(opt,OptionArray[2]) == 0)  // B1ecm
//	{
//		if (!isnumeric(arg)) {  //argument should be all numeric
//			printf("expected numeric input for option %s\n", opt);
//			exit(1);
//		}
//
//		fobj->ecm_obj.B1 = strtoul(arg,ptr,10);
//		if (fobj->pp1_obj.stg2_is_default)
//		{
//			//stg2 hasn't been specified yet, so set it to the default value
//			fobj->ecm_obj.B2 = 25 * (uint64)fobj->ecm_obj.B1;
//			//else, we have already specified a B2, so don't overwrite it with
//			//the default
//		}
//	}
//	else if (strcmp(opt,OptionArray[3]) == 0) // rhomax
//	{
//		if (!isnumeric(arg)) {  //argument should be all numeric
//			printf("expected numeric input for option %s\n", opt);
//			exit(1);
//		}
//
//		fobj->rho_obj.iterations = strtoul(arg,ptr,10);
//	}
//	else if (strcmp(opt,OptionArray[4]) == 0)  // B2pm1
//	{
//		if (!isnumeric(arg)) {  //argument should be all numeric
//			printf("expected numeric input for option %s\n", opt);
//			exit(1);
//		}
//
//		fobj->pm1_obj.B2 = strto_uint64(arg,ptr,10);
//		fobj->pm1_obj.stg2_is_default = 0;
//	}
//	else if (strcmp(opt,OptionArray[5]) == 0)   // B2pp1
//	{
//		if (!isnumeric(arg)) {  //argument should be all numeric
//			printf("expected numeric input for option %s\n", opt);
//			exit(1);
//		}
//
//		fobj->pp1_obj.B2 = strto_uint64(arg,ptr,10);
//		fobj->pp1_obj.stg2_is_default = 0;
//	}
//	else if (strcmp(opt,OptionArray[6]) == 0)  // B2ecm
//	{
//		if (!isnumeric(arg)) {  //argument should be all numeric
//			printf("expected numeric input for option %s\n", opt);
//			exit(1);
//		}
//
//		fobj->ecm_obj.B2 = strto_uint64(arg,ptr,10);
//		fobj->ecm_obj.stg2_is_default = 0;
//	}
//	else if (strcmp(opt,OptionArray[7]) == 0)  // qssave
//	{
//		//argument is a string
//	
//		if (strlen(arg) < sizeof(fobj->qs_obj.siqs_savefile))
//			strcpy_s(fobj->qs_obj.siqs_savefile, sizeof(fobj->qs_obj.siqs_savefile),
//				arg);
//		else
//			printf("*** argument to savefile too long, ignoring ***\n");
//	}
//	else if (strcmp(opt,OptionArray[8]) == 0)   // siqsB
//	{
//		if (!isnumeric(arg)) {  //argument should be all numeric
//			printf("expected numeric input for option %s\n", opt);
//			exit(1);
//		}
//
//		fobj->qs_obj.gbl_override_B = strtoul(arg,ptr,10);
//		fobj->qs_obj.gbl_override_B_flag = 1;
//	}
//	else if (strcmp(opt,OptionArray[9]) == 0)   // siqsTF
//	{
//		if (!isnumeric(arg)) {  //argument should be all numeric
//			printf("expected numeric input for option %s\n", opt);
//			exit(1);
//		}
//
//		fobj->qs_obj.gbl_override_tf = strtoul(arg,ptr,10);
//		fobj->qs_obj.gbl_override_tf_flag = 1;
//	}
//	else if (strcmp(opt,OptionArray[10]) == 0)  // siqsR
//	{
//		if (!isnumeric(arg)) {  //argument should be all numeric
//			printf("expected numeric input for option %s\n", opt);
//			exit(1);
//		}
//
//		fobj->qs_obj.gbl_override_rel = strtoul(arg,ptr,10);
//		fobj->qs_obj.gbl_override_rel_flag = 1;
//	}
//	else if (strcmp(opt,OptionArray[11]) == 0)    // siqsT
//	{
//		if (!isnumeric(arg)) {  //argument should be all numeric
//			printf("expected numeric input for option %s\n", opt);
//			exit(1);
//		}
//
//		fobj->qs_obj.gbl_override_time = strtoul(arg,ptr,10);
//		fobj->qs_obj.gbl_override_time_flag = 1;
//	}
//	else if (strcmp(opt,OptionArray[12]) == 0)    // siqsNB
//	{
//		if (!isnumeric(arg)) {  //argument should be all numeric
//			printf("expected numeric input for option %s\n", opt);
//			exit(1);
//		}
//
//		fobj->qs_obj.gbl_override_blocks = strtoul(arg,ptr,10);
//		fobj->qs_obj.gbl_override_blocks_flag = 1;
//	}
//	else if (strcmp(opt,OptionArray[13]) == 0)   // siqsM
//	{
//		if (!isnumeric(arg)) {  //argument should be all numeric
//			printf("expected numeric input for option %s\n", opt);
//			exit(1);
//		}
//
//		fobj->qs_obj.gbl_override_lpmult = strtoul(arg,ptr,10);
//		fobj->qs_obj.gbl_override_lpmult_flag = 1;
//	}
//	else if (strcmp(opt,OptionArray[14]) == 0)    // logfile
//	{
//		//argument is a string
//		if (strlen(arg) < sizeof(fobj->flogname))
//			strcpy_s(fobj->flogname, sizeof(fobj->flogname), arg);
//		else
//			printf("*** argument to logfile too long, ignoring ***\n");
//	}
//	else if (strcmp(opt,OptionArray[15]) == 0)  // batchfile
//	{
//		//argument is a string
//		if (strlen(arg) < 80)
//		{
//			strcpy_s(batchfilename, 80, arg);
//			USEBATCHFILE = 1;
//		}
//		else
//			printf("*** argument to batchfile too long, ignoring ***\n");
//	}
//	else if (strcmp(opt,OptionArray[16]) == 0)    // seed
//	{
//		USERSEED = 1;
//		sscanf(arg,"%u,%u",&g_rand.hi,&g_rand.low);
//	}
//	else if (strcmp(opt,OptionArray[17]) == 0)   // sigma
//	{
//	if (!isnumeric(arg)) {  //argument should be all numeric
//		printf("expected numeric input for option %s\n", opt);
//		exit(1);
//	}
//
//		fobj->ecm_obj.sigma = strtoul(arg,NULL,10);
//	}
//	else if (strcmp(opt,OptionArray[18]) == 0)    // session
//	{
//		//argument is a string
//		if (strlen(arg) < sizeof(sessionname))
//		{
//			strcpy_s(sessionname, sizeof(sessionname), arg);
//		}
//		else
//			printf("*** argument to sessionname too long, ignoring ***\n");
//	}
//	else if (strcmp(opt,OptionArray[19]) == 0)    // threads
//	{
//		if (!isnumeric(arg)) {  //argument should be all numeric
//			printf("expected numeric input for option %s\n", opt);
//			exit(1);
//		}
//
//		THREADS = strtoul(arg,NULL,10);
//	}
//	else if (strcmp(opt,OptionArray[20]) == 0) // v option
//	{
//		Vflag++;
//	}
//	else if (strcmp(opt,OptionArray[21]) == 0)  // silent
//	{
//		Vflag = -1;
//	}
//	else if (strcmp(opt,OptionArray[22]) == 0)  // pfile
//	{
//		PRIMES_TO_FILE = 1;
//	}
//	else if (strcmp(opt,OptionArray[23]) == 0)   // pscreen
//	{
//		PRIMES_TO_SCREEN = 1;
//	}
//	else if (strcmp(opt,OptionArray[24]) == 0)  // forceDLP
//	{
//		fobj->qs_obj.gbl_force_DLP = 1;
//	}
//	else if (strcmp(opt,OptionArray[25]) == 0)  // fmtmax
//	{
//		if (!isnumeric(arg)) {  //argument should be all numeric
//			printf("expected numeric input for option %s\n", opt);
//			exit(1);
//		}
//
//		fobj->div_obj.fmtlimit = strtoul(arg,NULL,10);
//	}
//	else if (strcmp(opt,OptionArray[26]) == 0)     // noopt
//	{
//		fobj->qs_obj.no_small_cutoff_opt = 1;
//	}
//	else if (strcmp(opt,OptionArray[27]) == 0)    // vproc
//	{
////#ifdef __APPLE__
////		printf("extended cpuid not supported\n");
////#else
//		VERBOSE_PROC_INFO++;
////#endif
//	}
//	else if (strcmp(opt,OptionArray[28]) == 0)   // noecm
//	{
//		fobj->autofact_obj.yafu_pretest_plan = PRETEST_NOECM;
//	}
//	else if (strcmp(opt,OptionArray[29]) == 0)  // ggnfs-dir
//	{
//		//argument is a string
//		if (strlen(arg) < sizeof(fobj->nfs_obj.ggnfs_dir))
//			strcpy_s(fobj->nfs_obj.ggnfs_dir, sizeof(fobj->nfs_obj.ggnfs_dir), 
//				arg);
//		else
//			printf("*** argument to ggnfs_dir too long, ignoring ***\n");
//	}
//	else if (strcmp(opt,OptionArray[30]) == 0)   // tune-info
//	{
//		//parse the tune_info string and if it matches the current OS and CPU, 
//		//set the appropriate globals
//		apply_tuneinfo(fobj, arg);
//	}
//	else if (strcmp(opt,OptionArray[31]) == 0) //argument "pretest_ratio"
//	{
//		/* convert text to double. Would strtod be better? */
//		sscanf(arg, "%lf", &fobj->autofact_obj.target_pretest_ratio);
//	}
//	else if (strcmp(opt,OptionArray[32]) == 0)  //argument "xover"
//	{  /* convert text to double. Would strtod be better? */
//		sscanf(arg, "%lf", &fobj->autofact_obj.qs_gnfs_xover);
//		fobj->autofact_obj.prefer_xover = 1;
//	}
//	else if (strcmp(opt,OptionArray[33]) == 0) //argument "one"
//	{
//		fobj->autofact_obj.want_only_1_factor = 1;
//	}
//	else if (strcmp(opt,OptionArray[34]) == 0)
//	{	//argument "op".  argument is a string
//		if (strlen(arg) < sizeof(fobj->autofact_obj.op_str))
//		{
//			strcpy_s(fobj->autofact_obj.op_str, sizeof(fobj->autofact_obj.op_str),
//				arg);
//			fobj->autofact_obj.want_output_primes = 1;
//		}
//		else
//			printf("*** argument to -op too long, ignoring ***\n");
//	}
//	else if (strcmp(opt,OptionArray[35]) == 0)
//	{	//option "of".  argument is a string
//		if (strlen(arg) < sizeof(fobj->autofact_obj.of_str))
//		{
//			strcpy_s(fobj->autofact_obj.of_str, sizeof(fobj->autofact_obj.of_str), 
//				arg);
//			fobj->autofact_obj.want_output_factors = 1;
//		}
//		else
//			printf("*** argument to -of too long, ignoring ***\n");
//	}
//	else if (strcmp(opt,OptionArray[36]) == 0)
//	{	//argument "out".  argument is a string
//		if (strlen(arg) < sizeof(fobj->autofact_obj.ou_str))
//		{
//			strcpy_s(fobj->autofact_obj.ou_str, sizeof(fobj->autofact_obj.ou_str),
//				arg);
//			fobj->autofact_obj.want_output_unfactored = 1;
//		}
//		else
//			printf("*** argument to -ou too long, ignoring ***\n");
//	}
//	else if (strcmp(opt, OptionArray[37]) == 0)
//	{	//argument "plan".  argument is a string
//		if (strlen(arg) < sizeof(fobj->autofact_obj.plan_str))
//		{
//			strcpy_s(fobj->autofact_obj.plan_str, 
//				sizeof(fobj->autofact_obj.plan_str), arg);
//
//			// test for recognized options.  
//			if (strcmp(fobj->autofact_obj.plan_str, "none") == 0)
//				fobj->autofact_obj.yafu_pretest_plan = PRETEST_NONE;
//			else if (strcmp(fobj->autofact_obj.plan_str, "noecm") == 0)
//				fobj->autofact_obj.yafu_pretest_plan = PRETEST_NOECM;
//			else if (strcmp(fobj->autofact_obj.plan_str, "light") == 0)
//				fobj->autofact_obj.yafu_pretest_plan = PRETEST_LIGHT;
//			else if (strcmp(fobj->autofact_obj.plan_str, "deep") == 0)
//				fobj->autofact_obj.yafu_pretest_plan = PRETEST_DEEP;
//			else if (strcmp(fobj->autofact_obj.plan_str, "normal") == 0)
//				fobj->autofact_obj.yafu_pretest_plan = PRETEST_NORMAL;
//			else if (strcmp(fobj->autofact_obj.plan_str, "custom") == 0)
//				fobj->autofact_obj.yafu_pretest_plan = PRETEST_CUSTOM;
//			else			
//			{
//				printf("*** unknown plan option, ignoring ***\n");
//				strcpy_s(fobj->autofact_obj.plan_str,
//					sizeof(fobj->autofact_obj.plan_str), "normal");
//			}
//		}
//		else
//			printf("*** argument to -plan too long, ignoring ***\n");
//	}
//	else if (strcmp(opt,OptionArray[38]) == 0)
//	{
//		//argument "pretest"
//		if (arg == NULL)
//		{
//			// no argument, use the value "1" to signify doing 
//			// pretesting to the bounds computed by factor()
//			fobj->autofact_obj.only_pretest = 1;
//		}
//		else
//		{
//			if (!isnumeric(arg)) {  //argument should be all numeric
//				printf("expected numeric input for option %s\n", opt);
//				exit(1);
//			}
//			// an optional argument to pretest is interpreted as
//			// a maximum t-level to pretest to
//			fobj->autofact_obj.only_pretest = strtoul(arg, NULL, 10);
//			
//			// default behavior
//			if (fobj->autofact_obj.only_pretest == 0)
//				fobj->autofact_obj.only_pretest = 1;
//		}
//	}
//	else if (strcmp(opt,OptionArray[39]) == 0)
//	{
//		//argument "no_expr"
//		fobj->autofact_obj.want_output_expressions = 0;
//	}	
//	else if (strcmp(opt,OptionArray[40]) == 0)
//	{
//		//argument "o".  Indicates output filename ggnfs sieving.
//		char *cptr;
//
//		if (strlen(arg) < sizeof(fobj->nfs_obj.outputfile))
//		{
//			char tmp[GSTR_MAXSIZE];
//			strcpy_s(fobj->nfs_obj.outputfile,
//				sizeof(fobj->nfs_obj.outputfile), arg);
//			strcpy_s(tmp, GSTR_MAXSIZE, fobj->nfs_obj.outputfile);
//			cptr = strchr(tmp, '.');   // look for '.'
//			if (cptr == NULL)
//			{
//				//no . in provided filename
//				sprintf_s(fobj->nfs_obj.logfile, 
//					sizeof(fobj->nfs_obj.logfile), "%s.log",fobj->nfs_obj.outputfile);
//				sprintf_s(fobj->nfs_obj.fbfile, 
//					sizeof(fobj->nfs_obj.fbfile), "%s.fb",fobj->nfs_obj.outputfile);
//			}
//			else
//			{				
//				cptr[0] = '\0';  // remove '.'
//				sprintf_s(fobj->nfs_obj.logfile, 
//					sizeof(fobj->nfs_obj.logfile), "%s.log", tmp);
//				sprintf_s(fobj->nfs_obj.fbfile, 
//					sizeof(fobj->nfs_obj.fbfile), "%s.fb", tmp);
//			}
//		}
//		else
//			printf("*** argument to -o too long, ignoring ***\n");
//	}
//	else if (strcmp(opt,OptionArray[41]) == 0)
//	{
//		//argument "a".  Indicates algebraic side special Q.
//		fobj->nfs_obj.sq_side = 1;
//	}
//	else if (strcmp(opt,OptionArray[42]) == 0)
//	{
//		//argument "r".  Indicates rational side special Q.
//		//fobj->nfs_obj.sq_side = 0;
//		fobj->nfs_obj.sq_side = -1;
//	}
//	else if (strcmp(opt,OptionArray[43]) == 0)
//	{
//		//argument "ggnfsT".  Indicates timeout (in seconds) for NFS job.
//		fobj->nfs_obj.timeout = strtoul(arg,NULL,10);
//	}
//	else if (strcmp(opt,OptionArray[44]) == 0)
//	{
//		//argument "job".  Indicates input .job file automated NFS.
//		if (strlen(arg) < sizeof(fobj->nfs_obj.job_infile))
//		{
//			strcpy_s(fobj->nfs_obj.job_infile, 
//				sizeof(fobj->nfs_obj.job_infile), arg);
//		}
//		else
//			printf("*** argument to -job too long, ignoring ***\n");
//	}
//	else if (strcmp(opt,OptionArray[45]) == 0)
//	{
//		char **nextptr = &(char *)arg;
//
//		//argument "ns".  do nfs sieving
//		fobj->nfs_obj.nfs_phases |= NFS_PHASE_SIEVE;
//		
//		if (arg != NULL)
//		{
//			// if an argument was supplied, parse the start and range of 
//			//special Q in the format X,Y
//			fobj->nfs_obj.startq = strtoul(arg, nextptr,10);
//
//			if (*nextptr[0] != ',')
//			{
//				printf("format of sieving argument is START,STOP\n");
//				exit(1);
//			}
//			fobj->nfs_obj.rangeq = strtoul(*nextptr + 1,NULL,10);
//
//			if (fobj->nfs_obj.startq >= fobj->nfs_obj.rangeq)
//			{
//				printf("format of sieving argument is START,STOP; STOP must be > START\n");
//				exit(1);
//			}
//			fobj->nfs_obj.rangeq = fobj->nfs_obj.rangeq - fobj->nfs_obj.startq;
//		}
//		else
//		{
//			fobj->nfs_obj.startq = 0;
//			fobj->nfs_obj.rangeq = 0;
//		}
//
//	}
//	else if (strcmp(opt,OptionArray[46]) == 0)
//	{		//argument "np".  do poly finding.
//		char **nextptr = &(char *)arg;
//
//		fobj->nfs_obj.nfs_phases |= NFS_PHASE_POLY;
//
//		if (arg != NULL)
//		{
//			// if an argument was supplied, parse the start and stop coefficient range in the
//			// format X,Y
//			fobj->nfs_obj.polystart = strtoul(arg, nextptr, 10);
//
//			if (*nextptr[0] != ',')
//			{
//				printf("format of poly select argument is START,STOP\n");
//				exit(1);
//			}
//			fobj->nfs_obj.polyrange = strtoul(*nextptr + 1,NULL,10);
//
//			if (fobj->nfs_obj.polystart >= fobj->nfs_obj.polyrange)
//			{
//				printf("format of poly select argument is START,STOP; STOP must be > START\n");
//				exit(1);
//			}
//			fobj->nfs_obj.polyrange = fobj->nfs_obj.polyrange - fobj->nfs_obj.polystart;
//		}
//		else
//		{
//			fobj->nfs_obj.polystart = 0;
//			fobj->nfs_obj.polyrange = 0;
//		}
//	}
//	else if (strcmp(opt,OptionArray[47]) == 0)
//	{	//argument "nc".  Do post processing, starting with filtering
//		fobj->nfs_obj.nfs_phases |= NFS_PHASE_FILTER;
//		fobj->nfs_obj.nfs_phases |= NFS_PHASE_LA;
//		fobj->nfs_obj.nfs_phases |= NFS_PHASE_SQRT;
//	}
//	else if (strcmp(opt,OptionArray[48]) == 0)
//	{	//argument "psearch".  modify poly search methodology
//		if (strlen(arg) < 1024)
//		{
//			if (strcmp(arg, "wide") == 0)
//				fobj->nfs_obj.poly_option = 1;
//			else if (strcmp(arg, "deep") == 0)
//				fobj->nfs_obj.poly_option = 2;
//			else if (strcmp(arg, "fast") == 0)
//				fobj->nfs_obj.poly_option = 0;
//			else
//			{
//				printf("option -psearch recognizes arguments 'deep', 'wide', or 'fast'.\n  see docfile.txt for details\n"); 
//				exit(1);
//			}
//
//		}
//		else
//			printf("*** argument to -psearch too long, ignoring ***\n");
//
//	}
//	else if (strcmp(opt,OptionArray[49]) == 0)
//	{
//		//argument "R".  nfs restart flag
//		fobj->nfs_obj.restart_flag = 1;
//	}
//	else if (strcmp(opt,OptionArray[50]) == 0)
//	{
//		//argument "pbatch".  Indicates size of blocks of leading coefficients to
//		//distribute to each thread in threaded NFS poly selection.
//		fobj->nfs_obj.polybatch = strtoul(arg,NULL,10);
//		if (fobj->nfs_obj.polybatch == 0)
//			fobj->nfs_obj.polybatch = 250;
//	}
//	else if (strcmp(opt,OptionArray[51]) == 0)
//	{
//		// argument "ecm_path"
//		//argument is a string
//		if (strlen(arg) < sizeof(fobj->ecm_obj.ecm_path))
//			strcpy_s(fobj->ecm_obj.ecm_path, sizeof(fobj->ecm_obj.ecm_path), 
//				arg);
//		else
//			printf("*** argument to ecm_path too long, ignoring ***\n");
//	}
//	else if (strcmp(opt,OptionArray[52]) == 0)
//	{
//		// argument "siever"
//		if (!isnumeric(arg)) {  //argument should be all numeric
//			printf("expected numeric input for option %s\n", opt);
//			exit(1);
//		}
//
//		fobj->nfs_obj.siever = strtoul(arg,NULL,10);
//	}
//	else if (strcmp(opt,OptionArray[53]) == 0)
//	{	//parameter  "ncr".  linear algebra restart flag
//		fobj->nfs_obj.nfs_phases |= NFS_PHASE_LA_RESUME;
//	}
//	else if (strcmp(opt,OptionArray[54]) == 0)
//	{	// parameter "lathreads"
//		if (!isnumeric(arg)) {  //argument should be all numeric
//			printf("expected numeric input for option %s\n", opt);
//			exit(1);
//		}
//
//		LATHREADS = strtoul(arg,NULL,10);
//	}
//	else if (strcmp(opt,OptionArray[55]) == 0)
//	{
//		//argument "nc2".  do linear algebra.
//		fobj->nfs_obj.nfs_phases |= NFS_PHASE_LA;
//	}
//	else if (strcmp(opt,OptionArray[56]) == 0)
//	{
//		//argument "nc3".  do nfs sqrt
//		fobj->nfs_obj.nfs_phases |= NFS_PHASE_SQRT;
//	}
//	else if (strcmp(opt, OptionArray[57]) == 0)
//	{
//		//argument "p".  set to idle priority.
//		//TODO: check to see if ggnfs and ecm binaries called through system
//		//retain idle priority
//		yafu_set_idle_priority();
//	}
//	else if (strcmp(opt,OptionArray[58]) == 0)
//	{
//		//argument "work"
//		sscanf(arg, "%lf", &fobj->autofact_obj.initial_work);
//	}
//	else if (strcmp(opt,OptionArray[59]) == 0)
//	{
//		//argument "nprp"
//		NUM_WITNESSES = strtoul(arg,NULL,10);
//	}
//	else if (strcmp(opt,OptionArray[60]) == 0)
//	{
//		// parameter "ecm_ext"
//		if (!isnumeric(arg)) {  //argument should be all numeric
//			printf("expected numeric input for option %s\n", opt);
//			exit(1);
//		}
//
//		fobj->ecm_obj.ecm_ext_xover = strtoul(arg,NULL,10);
//	}
//	else if (strcmp(opt,OptionArray[61]) == 0)
//	{
//		//argument "testsieve"
//		fobj->nfs_obj.snfs_testsieve_threshold = strtoul(arg,NULL,10);
//	}
//	else if (strcmp(opt,OptionArray[62]) == 0)
//	{
//		// argument "nt"
//		if (arg == NULL)
//		{
//			printf("expected argument for option %s\n", opt);
//			exit(1);
//		}
//		else
//			strcpy_s(fobj->nfs_obj.filearg, sizeof(fobj->nfs_obj.filearg), 
//				arg);
//	}
//	else if (strcmp(opt,OptionArray[63]) == 0)
//	{
//		// argument "aprcl_p", setting the threshold below which numbers
//		// are proved prime using APR-CL
//		if (!isnumeric(arg)) {  //argument should be all numeric
//			printf("expected numeric input for option %s\n", opt);
//			exit(1);
//		}
//		fobj->aprcl_prove_cutoff = strtoul(arg,NULL,10);
//		if (fobj->aprcl_prove_cutoff > 6021)
//		{
//			printf("APR-CL primality proving is possible only for numbers less"
//				" than 6021 digits... setting limit to 6021 digits\n");
//			fobj->aprcl_prove_cutoff = 6021;
//		}
//	}
//	else if (strcmp(opt,OptionArray[64]) == 0)
//	{
//		// argument "aprcl_d", setting the threshold above which numbers
//		// that are proved prime using APR-CL have additional verbosity enabled
//		if (!isnumeric(arg)) {  //argument should be all numeric
//			printf("expected numeric input for option %s\n", opt);
//			exit(1);
//		}
//		fobj->aprcl_display_cutoff = strtoul(arg,NULL,10);
//	}
//	else if (strcmp(opt,OptionArray[65]) == 0)
//	{
//		//argument "filt_bump"
//		sscanf(arg, "%lf", &fobj->nfs_obj.filter_min_rels_nudge);
//		fobj->nfs_obj.filter_min_rels_nudge = 1 + fobj->nfs_obj.filter_min_rels_nudge / 100;
//	}
//	else if (strcmp(opt,OptionArray[66]) == 0)
//	{
//		//argument "nc1".  do msieve filtering.
//		fobj->nfs_obj.nfs_phases |= NFS_PHASE_FILTER;
//	}
//	else
//	{
//		printf("invalid option %s\n",opt);
//		exit(1);
//	}
//
//	//zFree(&tmp);
//	return;
//}

/* functions to process all incoming arguments
arguments can come from the command line or the yafu.ini file
opt is the option as a Null-terminated ascii string,
arg is the argument as a Null-terminated ascii string, or NULL if no argument
If an error occurs exit with code 1
*/
static void applyOpt2(const char *opt, const char *arg, fact_obj_t *fobj,
	int ix)
{
	char **ptr = NULL;
	unsigned long long arg1Int=0, arg2Int=0;
	double arg1F=0.0;
	int numCount = 0;   // number of numbers extracted from arg

	if (option[ix].a == 1 && arg == NULL) {
		printf("argument expected for %s\n", opt);
		exit(1);
	}

	if (option[ix].a != 0 && arg != NULL) {
		int type = option[ix].type;  // 1 for one integer argument
			   // 2 for 2 integer arguments
			   // 3 for 1 floating point argument
		switch (type) {
		default:
			break;  // do nothing if type is not 1, 2 or 3
		case 1:  // type = 1 for one integer argument
			if (!isnumeric(arg)) {  //argument should be all numeric
				printf("expected numeric input for option %s\n", opt);
				exit(1);
			}
			arg1Int = strtoul(arg, ptr, 10);
			numCount = 1;
			break;

		case 2: // type = 2 for two integer arguments
			 numCount = sscanf(arg, "%llu,%llu", &arg1Int, &arg2Int);
			 if (numCount < 2) {
				 printf("expected number,number for option %s\n", opt);
			 }
			 if (numCount <= 1)
				 arg2Int = 0;
			 if (numCount <= 0)
				 arg1Int = 0;
			 break;

		case 3: // type = 3 for one floating point argument
			numCount = sscanf(arg, "%lf", &arg1F);
			break;
		}
	}


	switch (ix) {
		case 0:		 // B1pm1 - 1 integer argument
		{
			fobj->pm1_obj.B1 = (uint32)arg1Int;
			if (fobj->pm1_obj.stg2_is_default)
			{
				//stg2 hasn't been specified yet, so set it to the default value
				fobj->pm1_obj.B2 = 100 * (uint64)fobj->pm1_obj.B1;
				//else, we have already specified a B2, so don't overwrite it with
				//the default
			}
			break; 
		}

		case 1:  // B1pp1 - 1 integer argument
		{
			fobj->pp1_obj.B1 = (uint32)arg1Int;
			if (fobj->pp1_obj.stg2_is_default)
			{
				//stg2 hasn't been specified yet, so set it to the default value
				fobj->pp1_obj.B2 = 50 * (uint64)fobj->pp1_obj.B1;
				//else, we have already specified a B2, so don't overwrite it with
				//the default
			}
			break;
		}

		case 2:  // B1ecm  - 1 integer argument
		{
			fobj->ecm_obj.B1 = (uint32)arg1Int;
			if (fobj->pp1_obj.stg2_is_default)
			{
				//stg2 hasn't been specified yet, so set it to the default value
				fobj->ecm_obj.B2 = 25 * (uint64)fobj->ecm_obj.B1;
				//else, we have already specified a B2, so don't overwrite it with
				//the default
			}
			break;
		}

		case 3:		 // rhomax - 1 integer argument
		{
			fobj->rho_obj.iterations = (uint32)arg1Int;
			break;
		}

		case 4:  // B2pm1  - 1 integer argument
		{
			fobj->pm1_obj.B2 = arg1Int;
			fobj->pm1_obj.stg2_is_default = 0;
			break;
		}

		case 5:	   // B2pp1 - 1 integer argument
		{
			fobj->pp1_obj.B2 = arg1Int;
			fobj->pp1_obj.stg2_is_default = 0;
			break;
		}

		case 6:		 // B2ecm - 1 integer argument
		{
			fobj->ecm_obj.B2 = arg1Int;
			fobj->ecm_obj.stg2_is_default = 0;
			break;
		}

		case 7:  // qssave  argument is a string
		{	
			if (strlen(arg) < sizeof(fobj->qs_obj.siqs_savefile))
				strcpy_s(fobj->qs_obj.siqs_savefile, sizeof(fobj->qs_obj.siqs_savefile),
					arg);
			else
				printf("*** argument to savefile too long, ignoring ***\n");
			break;
		}

		case 8:	   // siqsB - 1 integer argument
		{
			fobj->qs_obj.gbl_override_B = (uint32)arg1Int;
			fobj->qs_obj.gbl_override_B_flag = 1;
			break;
		}

		case 9:  // siqsTF - 1 integer argument
 		{
			fobj->qs_obj.gbl_override_tf = (uint32)arg1Int;
			fobj->qs_obj.gbl_override_tf_flag = 1;
			break;
		}

		case 10:  // siqsR - 1 integer argument
		{
			fobj->qs_obj.gbl_override_rel = (uint32)arg1Int;
			fobj->qs_obj.gbl_override_rel_flag = 1;
			break;
		}

		case 11:    // siqsT - 1 integer argument
		{
			fobj->qs_obj.gbl_override_time = (uint32)arg1Int;
			fobj->qs_obj.gbl_override_time_flag = 1;
			break;
		}

		case 12:   // siqsNB - 1 integer argument
		{
			fobj->qs_obj.gbl_override_blocks = (uint32)arg1Int;
			fobj->qs_obj.gbl_override_blocks_flag = 1;
			break;
		}

		case 13:	   // siqsM - 1 integer argument
		{
			fobj->qs_obj.gbl_override_lpmult = (uint32)arg1Int;
			fobj->qs_obj.gbl_override_lpmult_flag = 1;
			break;
		}

		case 14:  // logfile   argument is a string
		{	
			if (strlen(arg) < sizeof(fobj->flogname))
				strcpy_s(fobj->flogname, sizeof(fobj->flogname), arg);
			else
				printf("*** argument to logfile too long, ignoring ***\n");
			break;
		}

		case 15:  // batchfile   argument is a string
		{	
			if (strlen(arg) < 80) {
				strcpy_s(batchfilename, 80, arg);
				USEBATCHFILE = 1;
			}
			else
				printf("*** argument to batchfile too long, ignoring ***\n");
			break;
		}

		case 16:	    // seed - needs 2 integer parameters
		{
			USERSEED = 1;
			//sscanf(arg, "%u,%u", &g_rand.hi, &g_rand.low);
			if (numCount < 2) {
				printf("expected number,number for option %s\n", opt);
				exit(1);
			}
			g_rand.hi = (uint32)arg1Int;
			g_rand.low = (uint32)arg2Int;
			break;
		}

		case 17:   // sigma - 1 integer argument
		{
			fobj->ecm_obj.sigma = (uint32)arg1Int;
			break;
		}

		case 18:    // session. argument is a string
		{	
			if (strlen(arg) < sizeof(sessionname))
			{
				strcpy_s(sessionname, sizeof(sessionname), arg);
			}
			else
				printf("*** argument to sessionname too long, ignoring ***\n");
			break;
		}

		case 19:    // threads - 1 integer argument
		{
			THREADS = (uint32)arg1Int;
			break;
		}

		case 20: // v option
		{
			Vflag++;
			break;
		}

		case 21: // silent
		{
			Vflag = -1;
			break;
		}

		case 22:	  // pfile
		{
			PRIMES_TO_FILE = 1;
			break;
		}

		case 23:	   // pscreen
		{
			PRIMES_TO_SCREEN = 1;
			break;
		}

		case 24:	  // forceDLP
		{
			fobj->qs_obj.gbl_force_DLP = 1;
			break;
		}

		case 25:	  // fmtmax - 1 integer argument
		{
			fobj->div_obj.fmtlimit = (uint32)arg1Int;
			break;
		}

		case 26: // noopt
   		{
			fobj->qs_obj.no_small_cutoff_opt = 1;
			break;
		}

		case 27:  // vproc
		{
			//#ifdef __APPLE__
			//		printf("extended cpuid not supported\n");
			//#else
			VERBOSE_PROC_INFO++;
			//#endif
			break;
		}

		case 28:   // noecm
		{
			fobj->autofact_obj.yafu_pretest_plan = PRETEST_NOECM;
			break;
		}

		case 29: // ggnfs-dir argument is a string
		{
			if (strlen(arg) < sizeof(fobj->nfs_obj.ggnfs_dir))
				strcpy_s(fobj->nfs_obj.ggnfs_dir, sizeof(fobj->nfs_obj.ggnfs_dir),
					arg);
			else
				printf("*** argument to ggnfs_dir too long, ignoring ***\n");
			break;
		}

		case 30:   // tune-info
		{	//parse the tune_info string and if it matches the current OS and CPU, 
			//set the appropriate globals
			apply_tuneinfo(fobj, arg);
			break;
		}

		case 31: //option "pretest_ratio" - 1 floating point argument
		{
			/* convert text to double. Would strtod be better? */
			//sscanf(arg, "%lf", &fobj->autofact_obj.target_pretest_ratio);
			fobj->autofact_obj.target_pretest_ratio = arg1F;
			break;
		}

		case 32:  //option "xover"- 1 floating point argument
		{  
			fobj->autofact_obj.qs_gnfs_xover = arg1F;
			fobj->autofact_obj.prefer_xover = 1;
			break;
		}

		case 33:	 //argument "one"
		{
			fobj->autofact_obj.want_only_1_factor = 1;
			break;
		}

		case 34:      //option  "op".  argument is a string
		{	
			if (strlen(arg) < sizeof(fobj->autofact_obj.op_str))
			{
				strcpy_s(fobj->autofact_obj.op_str, sizeof(fobj->autofact_obj.op_str),
					arg);
				fobj->autofact_obj.want_output_primes = 1;
			}
			else
				printf("*** argument to -op too long, ignoring ***\n");
			break;
		}

		case 35:  //option "of".  argument is a string
		{	
			if (strlen(arg) < sizeof(fobj->autofact_obj.of_str))
			{
				strcpy_s(fobj->autofact_obj.of_str, sizeof(fobj->autofact_obj.of_str),
					arg);
				fobj->autofact_obj.want_output_factors = 1;
			}
			else
				printf("*** argument to -of too long, ignoring ***\n");
			break;
		}

		case 36:   //option "ou".  argument is a string
		{	
			if (strlen(arg) < sizeof(fobj->autofact_obj.ou_str))
			{
				strcpy_s(fobj->autofact_obj.ou_str, sizeof(fobj->autofact_obj.ou_str),
					arg);
				fobj->autofact_obj.want_output_unfactored = 1;
			}
			else
				printf("*** argument to -ou too long, ignoring ***\n");
			break;
		}

		case 37:     //option "plan".  argument is a string
		{	
			if (strlen(arg) < sizeof(fobj->autofact_obj.plan_str))
			{
				strcpy_s(fobj->autofact_obj.plan_str,
					sizeof(fobj->autofact_obj.plan_str), arg);

				// test for recognized options.  
				if (strcmp(fobj->autofact_obj.plan_str, "none") == 0)
					fobj->autofact_obj.yafu_pretest_plan = PRETEST_NONE;
				else if (strcmp(fobj->autofact_obj.plan_str, "noecm") == 0)
					fobj->autofact_obj.yafu_pretest_plan = PRETEST_NOECM;
				else if (strcmp(fobj->autofact_obj.plan_str, "light") == 0)
					fobj->autofact_obj.yafu_pretest_plan = PRETEST_LIGHT;
				else if (strcmp(fobj->autofact_obj.plan_str, "deep") == 0)
					fobj->autofact_obj.yafu_pretest_plan = PRETEST_DEEP;
				else if (strcmp(fobj->autofact_obj.plan_str, "normal") == 0)
					fobj->autofact_obj.yafu_pretest_plan = PRETEST_NORMAL;
				else if (strcmp(fobj->autofact_obj.plan_str, "custom") == 0)
					fobj->autofact_obj.yafu_pretest_plan = PRETEST_CUSTOM;
				else
				{
					printf("*** unknown plan option, set to normal ***\n");
					strcpy_s(fobj->autofact_obj.plan_str,
						sizeof(fobj->autofact_obj.plan_str), "normal");
					fobj->autofact_obj.yafu_pretest_plan = PRETEST_NORMAL;
				}
			}
			else
				printf("*** argument to -plan too long, ignoring ***\n");
			break;
		}

		case 38:    //option "pretest" - 1 optional integer argument
		{
			if (arg == NULL) {
				// no argument, use the value "1" to signify doing 
				// pretesting to the bounds computed by factor()
				fobj->autofact_obj.only_pretest = 1;
			}
			else {
				fobj->autofact_obj.only_pretest = (uint32)arg1Int;

				// default behavior
				if (fobj->autofact_obj.only_pretest == 0)
					fobj->autofact_obj.only_pretest = 1;
			}
			break;
		}

		case 39:	//argument "no_expr"
		{
			fobj->autofact_obj.want_output_expressions = 0;
			break;
		}

		case 40:   //argument "o".  Indicates output filename ggnfs sieving.
		{
			char *cptr;

			if (strlen(arg) < sizeof(fobj->nfs_obj.outputfile))
			{
				char tmp[GSTR_MAXSIZE];
				strcpy_s(fobj->nfs_obj.outputfile,
					sizeof(fobj->nfs_obj.outputfile), arg);
				strcpy_s(tmp, GSTR_MAXSIZE, fobj->nfs_obj.outputfile);
				cptr = strchr(tmp, '.');   // look for '.'
				if (cptr == NULL)
				{
					//no . in provided filename
					sprintf_s(fobj->nfs_obj.logfile,
						sizeof(fobj->nfs_obj.logfile), "%s.log", fobj->nfs_obj.outputfile);
					sprintf_s(fobj->nfs_obj.fbfile,
						sizeof(fobj->nfs_obj.fbfile), "%s.fb", fobj->nfs_obj.outputfile);
				}
				else
				{
					cptr[0] = '\0';  // remove '.'
					sprintf_s(fobj->nfs_obj.logfile,
						sizeof(fobj->nfs_obj.logfile), "%s.log", tmp);
					sprintf_s(fobj->nfs_obj.fbfile,
						sizeof(fobj->nfs_obj.fbfile), "%s.fb", tmp);
				}
			}
			else
				printf("*** argument to -o too long, ignoring ***\n");
			break;
		}

		case 41:  //argument "a".  Indicates algebraic side special Q.
		{
			fobj->nfs_obj.sq_side = 1;
			break;
		}

		case 42:   //argument "r".  Indicates rational side special Q.
		{
			//fobj->nfs_obj.sq_side = 0;
			fobj->nfs_obj.sq_side = -1;
			break;
		}

		case 43:  /*argument "ggnfsT".  Indicates timeout (in seconds) for NFS job.
				 - 1 integer argument */
		{
			fobj->nfs_obj.timeout = (uint32)arg1Int;
			break;
		}

		case 44:	//argument "job".  Indicates input .job file automated NFS.
		{
			if (strlen(arg) < sizeof(fobj->nfs_obj.job_infile))
			{
				strcpy_s(fobj->nfs_obj.job_infile,
					sizeof(fobj->nfs_obj.job_infile), arg);
			}
			else
				printf("*** argument to -job too long, ignoring ***\n");
			break;
		}

		case 45:	//argument "ns".  do nfs sieving
		{
			fobj->nfs_obj.nfs_phases |= NFS_PHASE_SIEVE;

			if (arg != NULL) 	{
				// if an argument was supplied, parse the start and range of 
				//special Q in the format X,Y
				if (numCount < 2) {
					printf("format of sieving argument is START,STOP\n");
					exit(1);
				}
				fobj->nfs_obj.startq = (uint32)arg1Int;
				fobj->nfs_obj.rangeq = (uint32)arg2Int;

				if (fobj->nfs_obj.startq >= fobj->nfs_obj.rangeq)
				{
					printf("format of sieving argument is START,STOP; STOP must be > START\n");
					exit(1);
				}
				fobj->nfs_obj.rangeq = fobj->nfs_obj.rangeq - fobj->nfs_obj.startq;
			}
			else
			{
				fobj->nfs_obj.startq = 0;
				fobj->nfs_obj.rangeq = 0;
			}
			break;
		}

		case 46://argument "np".  do poly finding. 2 optional integer parameters
		{		
			fobj->nfs_obj.nfs_phases |= NFS_PHASE_POLY;

			if (arg != NULL) {
				// if an argument was supplied, parse the start and stop coefficient range in the
				// format X,Y
				fobj->nfs_obj.polystart = (uint32)arg1Int;

				if (numCount < 2)
				{
					printf("format of poly select argument is START,STOP\n");
					exit(1);
				}
				fobj->nfs_obj.polyrange = (uint32)arg2Int;

				if (fobj->nfs_obj.polystart >= fobj->nfs_obj.polyrange)
				{
					printf("format of poly select argument is START,STOP; STOP must be > START\n");
					exit(1);
				}
				fobj->nfs_obj.polyrange = fobj->nfs_obj.polyrange - fobj->nfs_obj.polystart;
			}
			else
			{
				fobj->nfs_obj.polystart = 0;
				fobj->nfs_obj.polyrange = 0;
			}
			break;
		}

		case 47:	//argument "nc".  Do post processing, starting with filtering
		{
			fobj->nfs_obj.nfs_phases |= NFS_PHASE_FILTER;
			fobj->nfs_obj.nfs_phases |= NFS_PHASE_LA;
			fobj->nfs_obj.nfs_phases |= NFS_PHASE_SQRT;
			break;
		}

		case 48:	//argument "psearch".  modify poly search methodology
		{
			if (strlen(arg) < 1024)
			{
				if (strcmp(arg, "wide") == 0)
					fobj->nfs_obj.poly_option = 1;
				else if (strcmp(arg, "deep") == 0)
					fobj->nfs_obj.poly_option = 2;
				else if (strcmp(arg, "fast") == 0)
					fobj->nfs_obj.poly_option = 0;
				else
				{
					printf("option -psearch recognizes arguments 'deep', 'wide', or 'fast'.\n  see docfile.txt for details\n");
					exit(1);
				}

			}
			else
				printf("*** argument to -psearch too long, ignoring ***\n");
			break;
		}

		case 49:	//argument "R".  nfs restart flag
		{
			fobj->nfs_obj.restart_flag = 1;		
			break;
		}

		case 50: // option "pbatch".  Indicates size of blocks of leading 
			// coefficients to distribute to each thread in threaded NFS 
			// poly selection. - 1 integer argument 
		{
			fobj->nfs_obj.polybatch = (uint32)arg1Int;
			if (fobj->nfs_obj.polybatch == 0)
				fobj->nfs_obj.polybatch = 250;
			break;
		}

		case 51:   // option "ecm_path" argument is a string
		{
			if (strlen(arg) < sizeof(fobj->ecm_obj.ecm_path))
				strcpy_s(fobj->ecm_obj.ecm_path, sizeof(fobj->ecm_obj.ecm_path),
					arg);
			else
				printf("*** argument to ecm_path too long, ignoring ***\n");
			break;
		}

		case 52:		//  "siever" - 1 integer argument 
		{
			fobj->nfs_obj.siever = (uint32)arg1Int;
			break;
		}

		case 53:   //parameter  "ncr".  linear algebra restart flag
		{	
			fobj->nfs_obj.nfs_phases |= NFS_PHASE_LA_RESUME;
			break;
		}

		case 54:  // parameter "lathreads" - 1 integer argument 
		{	
			LATHREADS = (int)arg1Int;
			break;
		}

		case 55:  //argument "nc2".  do linear algebra.
		{
			fobj->nfs_obj.nfs_phases |= NFS_PHASE_LA;
			break;
		}

		case 56:	//argument "nc3".  do nfs sqrt
		{
			fobj->nfs_obj.nfs_phases |= NFS_PHASE_SQRT;
			break;
		}

		case 57:   //argument "p".  set to idle priority.
			//TODO: check to see if ggnfs and ecm binaries called through system
			//retain idle priority
		{
			yafu_set_idle_priority();
			break;
		}

		case 58:	//"work" - 1 floating point argument
		{
			//sscanf(arg, "%lf", &fobj->autofact_obj.initial_work);
			fobj->autofact_obj.initial_work = arg1F;
			break;
		}

		case 59:	//argument "nprp" - 1 integer argument 
		{
			NUM_WITNESSES = (uint32)arg1Int;
			break;
		}

		case 60:	// parameter "ecm_ext" - 1 integer argument 
		{
			fobj->ecm_obj.ecm_ext_xover = (uint32)arg1Int;
			break;
		}

		case 61:	//"testsieve" - 1 integer argument 
		{
			fobj->nfs_obj.snfs_testsieve_threshold = (uint32)arg1Int;
			break;
		}

		case 62:	// argument "nt"
		{
			if (arg == NULL) {
				printf("expected argument for option %s\n", opt);
				exit(1);
			}
			else
				strcpy_s(fobj->nfs_obj.filearg, sizeof(fobj->nfs_obj.filearg),
					arg);
			break;
		}

		case 63:// argument "aprcl_p", setting the threshold below which numbers
				// are proved prime using APR-CL
			// - 1 integer argument 
		{
			fobj->aprcl_prove_cutoff = (uint32)arg1Int;
			if (fobj->aprcl_prove_cutoff > 6021)
			{
				printf("APR-CL primality proving is possible only for numbers less"
					" than 6021 digits... setting limit to 6021 digits\n");
				fobj->aprcl_prove_cutoff = 6021;
			}
			break;
		}

		case 64:// argument "aprcl_d", setting the threshold above which numbers
				// that are proved prime using APR-CL have additional verbosity enabled
				// 1 integer parameter
		{
			fobj->aprcl_display_cutoff = (int)arg1Int;
			break;
		}

		case 65:	// "filt_bump" 1 floating point parameter
		{
			//sscanf(arg, "%lf", &fobj->nfs_obj.filter_min_rels_nudge);
			fobj->nfs_obj.filter_min_rels_nudge = arg1F;
			fobj->nfs_obj.filter_min_rels_nudge = 1 + fobj->nfs_obj.filter_min_rels_nudge / 100;
			break;
		}

		case 66:	//"nc1".  do msieve filtering.
		{
			fobj->nfs_obj.nfs_phases |= NFS_PHASE_FILTER;
			break;		
		}

		default:
		{	/* invalid value of ix; indicates a logic error */
			printf("invalid option %s\n", opt);
			exit(1);
		}
	}
	return;
}

static void apply_tuneinfo(fact_obj_t *fobj, const char *arg)
{
	int i,j;
	char cpustr[80], osstr[80];

	//read up to the first comma - this is the cpu id string
	j=0;
	for (i=0; i<strlen(arg); i++)
	{
		if (arg[i] == '\n') break;   // LF?
		if (arg[i] == '\r') break;  // CR?
		if (arg[i] == ',') break;
		cpustr[j++] = arg[i];
	}
	cpustr[j] = '\0';
	i++;

	//read up to the next comma - this is the OS string
	j=0;
	for ( ; i<strlen(arg); i++)
	{
		if (arg[i] == '\n') break; // NL
		if (arg[i] == '\r') break; // CR
		if (arg[i] == ',') break;
		osstr[j++] = arg[i];
	}
	osstr[j] = '\0';

	//printf("found OS = %s and CPU = %s in tune_info field\n",osstr, cpustr);


#if defined(_WIN64)
//#ifdef _DEBUG
//	printf("CPU_ID_STR = \"%s\"\n", CPU_ID_STR);
//	printf("cpustr     = \"%s\"\n", cpustr);
//#endif
	if (strcmp(cpustr,CPU_ID_STR) == 0) 
		if (strcmp(osstr, "WIN64") == 0)
		{
			//printf("Applying tune_info entry for %s - %s\n",osstr,cpustr);
			double xoversave = fobj->autofact_obj.qs_gnfs_xover;
			sscanf(arg + i + 1, "%lg, %lg, %lg, %lg, %lg, %lg",
				&fobj->qs_obj.qs_multiplier, &fobj->qs_obj.qs_exponent,
				&fobj->nfs_obj.gnfs_multiplier, &fobj->nfs_obj.gnfs_exponent, 
				&fobj->autofact_obj.qs_gnfs_xover, &fobj->nfs_obj.gnfs_tune_freq);
			fobj->qs_obj.qs_tune_freq = fobj->nfs_obj.gnfs_tune_freq;
			if (fobj->autofact_obj.prefer_xover == 1)
				fobj->autofact_obj.qs_gnfs_xover = xoversave;
		}
#elif defined(WIN32)
	if ((strcmp(cpustr,CPU_ID_STR) == 0) && (strcmp(osstr, "WIN32") == 0))
	{
		//printf("Applying tune_info entry for %s - %s\n",osstr,cpustr);
		sscanf(arg + i + 1, "%lg, %lg, %lg, %lg, %lg, %lg",
			&fobj->qs_obj.qs_multiplier, &fobj->qs_obj.qs_exponent,
			&fobj->nfs_obj.gnfs_multiplier, &fobj->nfs_obj.gnfs_exponent, 
			&fobj->autofact_obj.qs_gnfs_xover, &fobj->nfs_obj.gnfs_tune_freq);
		fobj->qs_obj.qs_tune_freq = fobj->nfs_obj.gnfs_tune_freq;
	}
#elif BITS_PER_DIGIT == 64
	if ((strcmp(cpustr,CPU_ID_STR) == 0) && (strcmp(osstr, "LINUX64") == 0))
	{
		//printf("Applying tune_info entry for %s - %s\n",osstr,cpustr);
		
		sscanf(arg + i + 1, "%lg, %lg, %lg, %lg, %lg, %lg",
			&fobj->qs_obj.qs_multiplier, &fobj->qs_obj.qs_exponent,
			&fobj->nfs_obj.gnfs_multiplier, &fobj->nfs_obj.gnfs_exponent, 
			&fobj->autofact_obj.qs_gnfs_xover, &fobj->nfs_obj.gnfs_tune_freq);
		fobj->qs_obj.qs_tune_freq = fobj->nfs_obj.gnfs_tune_freq;
	}
#else 
	if ((strcmp(cpustr,CPU_ID_STR) == 0) && (strcmp(osstr, "LINUX32") == 0))
	{
		//printf("Applying tune_info entry for %s - %s\n",osstr,cpustr);
		
		sscanf(arg + i + 1, "%lg, %lg, %lg, %lg, %lg, %lg",
			&fobj->qs_obj.qs_multiplier, &fobj->qs_obj.qs_exponent,
			&fobj->nfs_obj.gnfs_multiplier, &fobj->nfs_obj.gnfs_exponent, 
			&fobj->autofact_obj.qs_gnfs_xover, &fobj->nfs_obj.gnfs_tune_freq);
		fobj->qs_obj.qs_tune_freq = fobj->nfs_obj.gnfs_tune_freq;
	}
#endif	

	//printf("QS_MULTIPLIER = %lg, QS_EXPONENT = %lg\nNFS_MULTIPLIER = %lg, NFS_EXPONENT = %lg\nXOVER = %lg, TUNE_FREQ = %lg\n",
	//	fobj->qs_obj.qs_multiplier, fobj->qs_obj.qs_exponent,
	//	fobj->nfs_obj.gnfs_multiplier, fobj->nfs_obj.gnfs_exponent, 
	//	fobj->autofact_obj.qs_gnfs_xover, fobj->qs_obj.qs_tune_freq);

	return;
}

//function get_random_seeds courtesy of Jason Papadopoulos
static void get_random_seeds(rand_t *r) {

	uint32 tmp_seed1, tmp_seed2;

	/* In a multithreaded program, every msieve object
	   should have two unique, non-correlated seeds
	   chosen for it */

	//in YAFU, make them available everywhere, by putting them in
	//a global structure that holds them.

#ifndef WIN32

	FILE *rand_device = fopen("/dev/urandom", "r");

	if (rand_device != NULL) {

		/* Yay! Cryptographic-quality nondeterministic randomness! */

		fread(&tmp_seed1, sizeof(uint32), (size_t)1, rand_device);
		fread(&tmp_seed2, sizeof(uint32), (size_t)1, rand_device);
		fclose(rand_device);
	}
	else

#endif
	{
		/* <Shrug> For everyone else, sample the current time,
		   the high-res timer (hopefully not correlated to the
		   current time), and the process ID. Multithreaded
		   applications should fold in the thread ID too */

		uint64 high_res_time = yafu_read_clock();
		tmp_seed1 = ((uint32)(high_res_time >> 32) ^
			     (uint32)time(NULL)) * 
			    (uint32)getpid();
		tmp_seed2 = (uint32)high_res_time;
	}

	/* The final seeds are the result of a multiplicative
	   hash of the initial seeds */

	r->low = tmp_seed1 * ((uint32)40499 * 65543);
	r->hi = tmp_seed2 * ((uint32)40499 * 65543);
}

