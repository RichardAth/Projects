#pragma once

extern bool breakSignal;
extern bool fpeReRegister;

/* flags used to control action of SetProcessExceptionHandlers() function */
struct flags {
	unsigned char UnEx : 1;        /* set unhandled exception 'filter' */
	unsigned char PureCall : 1;    /* Catch pure virtual function calls. */
	unsigned char InvParam : 1;    /* set invalid parameter handler */
	unsigned char New : 1;         /* Catch new operator memory allocation exceptions*/
	unsigned char abort : 1;       /* catch abort signal */
	unsigned char interrupt : 1;   /* catch ctrl-c & ctrl-break signals */
	unsigned char sigterm : 1;     /* Catch a termination request signal */
	unsigned char term : 1;        /* set up terminate handler */
	unsigned char unexpected : 1;  /* set up 'unexpected' handler */
	unsigned char sigfpe : 1;      /* catch floating point error signal */
	unsigned char sigill : 1;      /* catch illegal instruction signal */
	unsigned char sigsegv : 1;     /* Catch illegal storage access error signal */
};

static void __cdecl UnexpectedHandler();   /* 'unexpected' handler */
static void __cdecl TerminateHandler();    /* terminate handler */
static void __cdecl PureCallHandler();     /* pure virtual function calls handler*/

static int __cdecl NewHandler(size_t);     /*  new operator memory allocation exceptions */

/* signal handlers */
static void SigabrtHandler(int);     /* abort signal handler */
static void SigintHandler(int);      /* interrupt signal handler */
static void SigillHandler(int);      /* illegal instruction handler */
static void SigsegvHandler(int);     /* illegal storage access handler */
static void SigtermHandler(int);
void SigfpeHandler(int sig, int num);  /* floating point error handler */

unsigned long StackTrace(const EXCEPTION_POINTERS *ep);
unsigned long StackTrace2(void);
void InvalidParameterHandler(const wchar_t* expression,
	const wchar_t* function,
	const wchar_t* file,
	unsigned int line,
	size_t pReserved);


int handle_program_memory_depletion(size_t memsize);
void SetProcessExceptionHandlers(flags f);

void createMiniDump(const EXCEPTION_POINTERS* pExcPtrs);

/* Structured Exception Handler. called if the __except section of a function is
entered, . This can happen for a variety of causes, e.g. signals, divide errors etc.
*/
int filter(unsigned int code, const struct _EXCEPTION_POINTERS *ep);

/* does same as filter, but doesn't need code value */
long filter2(struct _EXCEPTION_POINTERS *ep);

/* call this to generate one of a variety of errors; test error handling */
void testerrors(void);

/* re-register the FP error signal handler */
void ReRegister(void);