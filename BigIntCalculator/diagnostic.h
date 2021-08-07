#pragma once

struct flags {
	unsigned char UnEx : 1;
	unsigned char PureCall : 1;
	unsigned char InvParam : 1;
	unsigned char New : 1;
	unsigned char abort : 1;
	unsigned char interrupt : 1;
	unsigned char sigterm : 1;
	unsigned char term : 1;
	unsigned char unexpected : 1;
	unsigned char sigfpe : 1;
	unsigned char sigill : 1;
	unsigned char sigsegv : 1;
};

static void __cdecl UnexpectedHandler();
static void __cdecl TerminateHandler();
static void __cdecl PureCallHandler();

static int __cdecl NewHandler(size_t);

static void SigabrtHandler(int);
static void SigintHandler(int);
static void SigillHandler(int);
static void SigsegvHandler(int);
static void SigtermHandler(int);
void SigfpeHandler(int sig, int num);

int filter(unsigned int code, struct _EXCEPTION_POINTERS *ep);
unsigned long StackTrace(EXCEPTION_POINTERS *ep);
unsigned long StackTrace2(void);
void InvalidParameterHandler(const wchar_t* expression,
	const wchar_t* function,
	const wchar_t* file,
	unsigned int line,
	size_t pReserved);


int handle_program_memory_depletion(size_t memsize);
void SetProcessExceptionHandlers(flags f);

void createMiniDump(EXCEPTION_POINTERS* pExcPtrs);

/* Structured Exception Handler. called if the __except section of a function is
entered, . This can happen for a variety of causes, e.g. signals, divide errors etc.
*/
int filter(unsigned int code, struct _EXCEPTION_POINTERS *ep);

/* does same as filter, but doesn't need code value */
long filter2(struct _EXCEPTION_POINTERS *ep);

/* call this to generate one of a variety of errors; test error handling */
void testerrors(void);