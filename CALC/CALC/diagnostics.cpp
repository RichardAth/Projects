#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <cstdint>
#include <cfloat>
#include <array>
#include <vector>
#include <string>
#include <iostream>
#include <Windows.h>
#include <setjmp.h>
#include <signal.h>

#include <excpt.h>
#include <new.h>
#include <eh.h>
#include <sstream>

#include <Psapi.h>
/* To ensure correct resolution of symbols, add Psapi.lib and dbghelp.lib
 In visual studio:
 project -> properties -> linker -> input -> additional dependencies */

#include <dbghelp.h>
#include <algorithm>
#include <iterator>

/* Kludge */
//#ifdef __cplusplus 
//EXTERN_C{
//#endif
#include <minidumpapiset.h>


#ifdef __cplusplus 
EXTERN_C{
#endif


static void __cdecl UnexpectedHandler();
static void __cdecl TerminateHandler();
static void __cdecl PureCallHandler();
long filter2(struct _EXCEPTION_POINTERS *ep);

static void __cdecl InvalidParameterHandler(const wchar_t* expression,
		const wchar_t* function, const wchar_t* file,
		unsigned int line, uintptr_t pReserved);

static int __cdecl NewHandler(size_t);

static void SigabrtHandler(int);
static void SigfpeHandler(int /*code*/, int subcode);
static void SigintHandler(int);
static void SigillHandler(int);
static void SigsegvHandler(int);
static void SigtermHandler(int);
int filter(unsigned int code, struct _EXCEPTION_POINTERS *ep);


struct module_data{
	std::string image_name;
	std::string module_name;
	void *base_address;
	DWORD load_size;
};

class symbol {
	typedef IMAGEHLP_SYMBOL64 sym_type;
	sym_type *sym;
	static const int max_name_len = 1024;

public:
	symbol(HANDLE process, DWORD64 address) : sym((sym_type *)::operator new(sizeof(*sym) + max_name_len)) {
		memset(sym, '\0', sizeof(*sym) + max_name_len);
		sym->SizeOfStruct = sizeof(*sym);
		sym->MaxNameLength = max_name_len;
		DWORD64 displacement;

		SymGetSymFromAddr64(process, address, &displacement, sym);
	}

	std::string name() { return std::string(sym->Name); }
	std::string undecorated_name() {
		if (*sym->Name == '\0')
			return "<couldn't map PC to fn name>";
		std::vector<char> und_name(max_name_len);
		UnDecorateSymbolName(sym->Name, &und_name[0], max_name_len, UNDNAME_COMPLETE);
		return std::string(&und_name[0], strlen(&und_name[0]));
	}
};


class get_mod_info {
	HANDLE process;
	static const int buffer_length = 4096;
public:
	get_mod_info(HANDLE h) : process(h) {}

	module_data operator()(HMODULE module) {
		module_data ret;
		char temp[buffer_length];
		MODULEINFO mi;

		GetModuleInformation(process, module, &mi, sizeof(mi));
		ret.base_address = mi.lpBaseOfDll;
		ret.load_size = mi.SizeOfImage;

		GetModuleFileNameExA(process, module, temp, sizeof(temp));
		ret.image_name = temp;
		GetModuleBaseNameA(process, module, temp, sizeof(temp));
		ret.module_name = temp;
		std::vector<char> img(ret.image_name.begin(), ret.image_name.end());
		std::vector<char> mod(ret.module_name.begin(), ret.module_name.end());
		SymLoadModule64(process, 0, &img[0], &mod[0], (DWORD64)ret.base_address, ret.load_size);
		return ret;
	}
};


// if you use C++ exception handling: install a translator function
// with set_se_translator(). In the context of that function (but *not*
// afterwards), you can either do your stack dump, or save the CONTEXT
// record as a local copy. Note that you must do the stack dump at the
// earliest opportunity, to avoid the interesting stack-frames being gone
// by the time you do the dump.
DWORD DumpStackTrace(EXCEPTION_POINTERS *ep)
{
	HANDLE process = GetCurrentProcess();
	HANDLE hThread = GetCurrentThread();
	int frame_number = 0;
	DWORD offset_from_symbol = 0;
	IMAGEHLP_LINE64 line = { 0 };
	std::vector<module_data> modules;
	DWORD cbNeeded;
	std::vector<HMODULE> module_handles(1);

	// Load the symbols:
	// WARNING: You'll need to replace <pdb-search-path> with either NULL
	// or some folder where your clients will be able to find the .pdb file.
	if (!SymInitialize(process, nullptr, false))
		//throw(std::logic_error("Unable to initialize symbol handler"));
	{
		std::cout << "dumpstack failure; Unable to initialize symbol handler\n";
		abort();
	}
	DWORD symOptions = SymGetOptions();
	symOptions |= SYMOPT_LOAD_LINES | SYMOPT_UNDNAME;
	SymSetOptions(symOptions);
	EnumProcessModules(process, &module_handles[0],
		(DWORD)module_handles.size() * sizeof(HMODULE), &cbNeeded);
	module_handles.resize(cbNeeded / sizeof(HMODULE));
	EnumProcessModules(process, &module_handles[0],
		(DWORD)module_handles.size() * sizeof(HMODULE), &cbNeeded);
	std::transform(module_handles.begin(), module_handles.end(), std::back_inserter(modules), get_mod_info(process));
	void *base = modules[0].base_address;

	// Setup stuff:
	CONTEXT* context = ep->ContextRecord;
#ifdef _M_X64
	STACKFRAME64 frame;
	frame.AddrPC.Offset = context->Rip;
	frame.AddrPC.Mode = AddrModeFlat;
	frame.AddrStack.Offset = context->Rsp;
	frame.AddrStack.Mode = AddrModeFlat;
	frame.AddrFrame.Offset = context->Rbp;
	frame.AddrFrame.Mode = AddrModeFlat;
#else
	STACKFRAME64 frame;
	frame.AddrPC.Offset = context->Eip;
	frame.AddrPC.Mode = AddrModeFlat;
	frame.AddrStack.Offset = context->Esp;
	frame.AddrStack.Mode = AddrModeFlat;
	frame.AddrFrame.Offset = context->Ebp;
	frame.AddrFrame.Mode = AddrModeFlat;
#endif
	line.SizeOfStruct = sizeof line;
	IMAGE_NT_HEADERS *h = ImageNtHeader(base);
	DWORD image_type = h->FileHeader.Machine;
	int n = 0;

	// Build the string:
	std::ostringstream builder;
	do {
		if (frame.AddrPC.Offset != 0) {
			std::string fnName = symbol(process, frame.AddrPC.Offset).undecorated_name();
			builder << fnName;
			if (SymGetLineFromAddr64(process, frame.AddrPC.Offset, &offset_from_symbol, &line))
				builder << "  " /*<< line.FileName*/ << "(" << line.LineNumber << ")\n";
			else builder << "\n";
			if (fnName == "main")
				break;  /* stop when we get to top level of all code */
			//if (fnName == "RaiseException") {
			//	// This is what we get when compiled in Release mode:
			//	printf_s("%s", "the program has crashed.\n\n");
			//	break;
			//}
		}
		else
			builder << "(No Symbols: PC == 0)";
		if (!StackWalk64(image_type, process, hThread, &frame, context, NULL,
			SymFunctionTableAccess64, SymGetModuleBase64, NULL))
			break;
		if (++n > 10)
			break;
	} while (frame.AddrReturn.Offset != 0);


	SymCleanup(process);

	// Display the string:
	printf_s("%s %s", "the program has crashed. Stack Trace; \n",
		builder.str().c_str());
	return EXCEPTION_EXECUTE_HANDLER;
}

DWORD DumpStackTrace2(void)
{
	HANDLE process = GetCurrentProcess();
	HANDLE hThread = GetCurrentThread();
	int frame_number = 0;
	DWORD offset_from_symbol = 0;
	IMAGEHLP_LINE64 line = { 0 };
	std::vector<module_data> modules;
	DWORD cbNeeded;
	int NumModules;
	std::vector<HMODULE> module_handles(5);

	// Load the symbols:
	// WARNING: You may need to replace nullptr with some folder where 
	// your clients will be able to find the .pdb file.
	if (!SymInitialize(process, nullptr, false))
		//throw(std::logic_error("Unable to initialize symbol handler"));
	{
		std::cout << "dumpstack failure; Unable to initialize symbol handler\n";
		abort();
	}

	DWORD symOptions = SymGetOptions();
	symOptions |= SYMOPT_LOAD_LINES | SYMOPT_UNDNAME;
	SymSetOptions(symOptions);

	EnumProcessModules(process, &module_handles[0],
		sizeof(module_handles), &cbNeeded);
	module_handles.resize(cbNeeded / sizeof(HMODULE));
	EnumProcessModules(process, &module_handles[0],
		(DWORD)module_handles.size() * sizeof(HMODULE), &cbNeeded);

	NumModules = cbNeeded / sizeof(HMODULE);

	std::transform(module_handles.begin(), module_handles.end(),
		std::back_inserter(modules), get_mod_info(process));


	void *base = modules[0].base_address;
	// Setup stuff:

	CONTEXT ct;
	RtlCaptureContext(&ct);

#ifdef _M_X64
	STACKFRAME64 frame;
	frame.AddrPC.Offset = ct.Rip;
	frame.AddrPC.Mode = AddrModeFlat;
	frame.AddrStack.Offset = ct.Rsp;
	frame.AddrStack.Mode = AddrModeFlat;
	frame.AddrFrame.Offset = ct.Rbp;
	frame.AddrFrame.Mode = AddrModeFlat;
#else
	STACKFRAME64 frame;
	frame.AddrPC.Offset = ct.Eip;
	frame.AddrPC.Mode = AddrModeFlat;
	frame.AddrStack.Offset = ct.Esp;
	frame.AddrStack.Mode = AddrModeFlat;
	frame.AddrFrame.Offset = ct.Ebp;
	frame.AddrFrame.Mode = AddrModeFlat;
#endif
	line.SizeOfStruct = sizeof line;
	IMAGE_NT_HEADERS *h = ImageNtHeader(base);
	DWORD image_type = h->FileHeader.Machine;
	int n = 0;   // count of stack entries processed

	// Build the string:
	char * builder;
	size_t bIx = 0;
	const size_t bIxMax = 10000;
	builder = (char *)calloc(bIxMax, 1);
	do {
		if (frame.AddrPC.Offset != 0) {
			std::string fnName = symbol(process, frame.AddrPC.Offset).undecorated_name();
			strcat_s(builder, bIxMax - bIx, fnName.c_str());
			if (SymGetLineFromAddr64(process, frame.AddrPC.Offset, &offset_from_symbol, &line)) {
				bIx = strlen(builder);
				strcat_s(builder, bIxMax - bIx, "  (line ");
				bIx = strlen(builder);
				sprintf_s(builder + bIx, bIxMax - bIx, "%4d) \n", line.LineNumber);
			}
			else strcat_s(builder, bIxMax - bIx, "\n");

			if (fnName == "main")
				break;    // exit when top level in stack is reached

			//if (fnName == "RaiseException") {
			//	// This is what we get when compiled in Release mode:
			//	printf_s("%s", "the program has crashed.\n\n");
			//	break;
			//}
		}
		else
			strcat_s(builder, bIxMax - bIx, "(No Symbols: PC == 0)\n");
		if (!StackWalk64(image_type, process, hThread, &frame, &ct, NULL,
			SymFunctionTableAccess64, SymGetModuleBase64, NULL))
			break;
		if (++n > 20)
			break;   // exit if more than 20 stack entries
	} while (frame.AddrReturn.Offset != 0);

	SymCleanup(process);

	// Display the string:
	printf_s("%s%s", "the program has crashed; Stack Trace: \n",
		builder);
	free(builder);
	return EXCEPTION_EXECUTE_HANDLER;
}

void myInvalidParameterHandler(const wchar_t* expression,
	const wchar_t* function,
	const wchar_t* file,
	unsigned int line,
	uintptr_t pReserved)
{
	wprintf(L"Invalid parameter detected in function %s."
		L" File: %s Line: %d\n", function, file, line);
	wprintf(L"Expression: %s\n", expression);
}

/* convert Floating Point error code to text */
const char* FPEcodeToTxt(int num) {
	static char buf[80];

	switch (num) {
	case FPE_INVALID:
		return "invalid";
	case FPE_DENORMAL:
		return "denormal";
	case FPE_ZERODIVIDE:
		return "Divide by zero";
	case FPE_OVERFLOW:
		return "Overflow";
	case FPE_UNDERFLOW:
		return "underflow";
	case FPE_INEXACT:
		return "inexact";
	case FPE_UNEMULATED:
		return "unemulated";
	case FPE_SQRTNEG:
		return "square root of -ve number";
	case FPE_STACKOVERFLOW:
		return "stack overflow";
	case FPE_STACKUNDERFLOW:
		return "stack underflow";
	case FPE_EXPLICITGEN:
		return "raise(SIGFPE) called";
	case _FPE_MULTIPLE_TRAPS:
		return "multiple traps";
	case _FPE_MULTIPLE_FAULTS:
		return "multiple faults";
	default:
		sprintf_s(buf, sizeof(buf), "unknown floating point exception: code %#x ",
			num);
		return buf;
	}

	return NULL;  /* error; should never get here */
}

/* the first parameter should have value SIG_FPE
the 2nd parameter can have any of the following values:
FPE_INVALID			0x81   This exception occurs when the result of an operation is
						   ill-defined, such as (0.0 / 0.0). If trapping is enabled,
						   the floating-point-invalid condition is signalled.
						   Otherwise, a quiet NaN is returned.
FPE_DENORMAL		0x82   operation results in a denormalized float?
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
static void SigfpeHandler(int sig, int num) {
	const char *msg = NULL;
	if (num == FPE_INEXACT)
		return;

	msg = FPEcodeToTxt(num);
	/* strictly speaking, should not call I/O routines here, but it seems to work OK */
	fprintf_s(stderr, "\nreceived FPE signal %d %s; continue \n", sig, msg);
	DumpStackTrace2();
	//__debugbreak();     // try to enter debuggger (to look at call stack)
	Beep(1200, 1000);   // beep at 1200 Hz for 1 second
	system("PAUSE");    // press any key to continue

	return;
}

static void SigabrtHandler(int sig) {
	fprintf(stderr, "Abort occurred\n");
	DumpStackTrace2();
	//__debugbreak();     // try to enter debuggger (to look at call stack)
	Beep(1200, 1000);   // beep at 1200 Hz for 1 second
	system("PAUSE");    // press any key to continue
	exit(EXIT_FAILURE);
}

static void SigillHandler(int sig) {
	fprintf(stderr, "Illegal Instruction \n");
	DumpStackTrace2();
	//__debugbreak();     // try to enter debuggger (to look at call stack)
	Beep(1200, 1000);   // beep at 1200 Hz for 1 second
	system("PAUSE");    // press any key to continue
	exit(EXIT_FAILURE);
}

static void SigsegvHandler(int sig) {
	fprintf(stderr, "Segment Violation \n");
	DumpStackTrace2();
	//__debugbreak();     // try to enter debuggger (to look at call stack)
	Beep(1200, 1000);   // beep at 1200 Hz for 1 second
	system("PAUSE");    // press any key to continue
	exit(EXIT_FAILURE);
}

static void SigintHandler(int sig) {
	fprintf(stderr, "Interrupt \n");
	DumpStackTrace2();
	//__debugbreak();     // try to enter debuggger (to look at call stack)
	Beep(1200, 1000);   // beep at 1200 Hz for 1 second
	system("PAUSE");    // press any key to continue
	exit(EXIT_FAILURE);
}

static void SigtermHandler(int sig) {
	fprintf(stderr, "termination \n");
	DumpStackTrace2();
	//__debugbreak();     // try to enter debuggger (to look at call stack)
	Beep(1200, 1000);   // beep at 1200 Hz for 1 second
	system("PAUSE");    // press any key to continue
	exit(EXIT_FAILURE);
}

int handle_program_memory_depletion(size_t memsize)
{
	// Your code
	system("PAUSE");   /* temporary */
	return 0;
}

void SetProcessExceptionHandlers()
{
	SetUnhandledExceptionFilter(filter2);
	// Catch pure virtual function calls.
	// Because there is one _purecall_handler for the whole process, 
	// calling this function immediately impacts all threads. The last 
	// caller on any thread sets the handler. 
	// http://msdn.microsoft.com/en-us/library/t296ys27.aspx
	//_set_purecall_handler(PureCallHandler);

	_invalid_parameter_handler oldHandler, newHandler;
	newHandler = myInvalidParameterHandler;
	oldHandler = _set_invalid_parameter_handler(newHandler);

	// Catch new operator memory allocation exceptions
	_set_new_handler(handle_program_memory_depletion);

	// Set up C++ signal handlers

	_set_abort_behavior(_CALL_REPORTFAULT, _CALL_REPORTFAULT);

	// Catch an abnormal program termination
	signal(SIGABRT, SigabrtHandler);

	// Catch interrupt handler
	signal(SIGINT, SigintHandler);

	// Catch a termination request
	signal(SIGTERM, SigtermHandler);
}

void SetThreadExceptionHandlers()
{

	// Catch terminate() calls. 
	// In a multithreaded environment, terminate functions are maintained 
	// separately for each thread. Each new thread needs to install its own 
	// terminate function. Thus, each thread is in charge of its own termination handling.
	// http://msdn.microsoft.com/en-us/library/t6fk7h29.aspx
	//set_terminate(TerminateHandler);

	// Catch unexpected() calls.
	// In a multithreaded environment, unexpected functions are maintained 
	// separately for each thread. Each new thread needs to install its own 
	// unexpected function. Thus, each thread is in charge of its own unexpected handling.
	// http://msdn.microsoft.com/en-us/library/h46t5b69.aspx  
	//set_unexpected(UnexpectedHandler);

	/* Catch a floating point error. 
	Note: some errors e.g. floating point errors will be caught by the filter
	function if a signal handler is not registered */
	typedef void(*sigh)(int); /* force compiler to accept SigfpeHandler although
							  it actually has 2 parameters, not 1 */
	//signal(SIGFPE, (sigh)SigfpeHandler);

	// Catch an illegal instruction
	signal(SIGILL, SigillHandler);

	// Catch illegal storage access errors
	signal(SIGSEGV, SigsegvHandler);

}

static void createMiniDump(EXCEPTION_POINTERS* pExcPtrs)
{
	HMODULE hDbgHelp = NULL;
	HANDLE hFile = NULL;
	MINIDUMP_EXCEPTION_INFORMATION mei;
	MINIDUMP_CALLBACK_INFORMATION mci;

	// Load dbghelp.dll
	hDbgHelp = LoadLibrary(TEXT("dbghelp.dll"));
	if (hDbgHelp == NULL)
	{
		// Error - couldn't load dbghelp.dll
		return;
	}

	// Create the minidump file
	hFile = CreateFile(
		TEXT("crashdump.dmp"),
		GENERIC_WRITE,
		0,
		NULL,
		CREATE_ALWAYS,
		FILE_ATTRIBUTE_NORMAL,
		NULL);

	if (hFile == INVALID_HANDLE_VALUE)
	{
		// Couldn't create file
		return;
	}

	// Write minidump to the file
	mei.ThreadId = GetCurrentThreadId();
	mei.ExceptionPointers = pExcPtrs;
	mei.ClientPointers = FALSE;
	mci.CallbackRoutine = NULL;
	mci.CallbackParam = NULL;

	typedef BOOL(WINAPI *LPMINIDUMPWRITEDUMP)(
		HANDLE hProcess,
		DWORD ProcessId,
		HANDLE hFile,
		MINIDUMP_TYPE DumpType,
		CONST PMINIDUMP_EXCEPTION_INFORMATION ExceptionParam,
		CONST PMINIDUMP_USER_STREAM_INFORMATION UserEncoderParam,
		CONST PMINIDUMP_CALLBACK_INFORMATION CallbackParam);

	LPMINIDUMPWRITEDUMP pfnMiniDumpWriteDump =
		(LPMINIDUMPWRITEDUMP)GetProcAddress(hDbgHelp, "MiniDumpWriteDump");
	if (!pfnMiniDumpWriteDump)
	{
		// Bad MiniDumpWriteDump function
		return;
	}

	HANDLE hProcess = GetCurrentProcess();
	DWORD dwProcessId = GetCurrentProcessId();

	BOOL bWriteDump = pfnMiniDumpWriteDump(
		hProcess,
		dwProcessId,
		hFile,
		MiniDumpNormal,
		&mei,
		NULL,
		&mci);

	if (!bWriteDump)
	{
		// Error writing dump.
		return;
	}

	// Close file
	CloseHandle(hFile);

	// Unload dbghelp.dll
	FreeLibrary(hDbgHelp);
	fprintf_s(stderr, "minidump file created; crashdump.dmp \n");
}

/* convert an error code to text */
const char *getText(const int errorcode) {
	static char msg[80];

	switch (errorcode) {
	case EXCEPTION_ACCESS_VIOLATION: {
			return "Access Violation";
		}
	case EXCEPTION_DATATYPE_MISALIGNMENT: {
			return "Datatype Misalignment";
		}
	case EXCEPTION_BREAKPOINT: {
			return "Breakpoint";
		}
	case EXCEPTION_SINGLE_STEP: {
			return "Single Step";
		}
	case EXCEPTION_ARRAY_BOUNDS_EXCEEDED: {
			return "Array Bounds Exceeded";
		}
	case EXCEPTION_FLT_DENORMAL_OPERAND: {
			return "Float Denormal Operand";
		}
	case EXCEPTION_FLT_DIVIDE_BY_ZERO: {
			return "Float divide by zero";
		}
	case EXCEPTION_FLT_INEXACT_RESULT: {
			return "Float inexact result";
		}
	case EXCEPTION_FLT_INVALID_OPERATION: {
				return "Float invalid operation";
			}
	case EXCEPTION_FLT_OVERFLOW: {
			return "Float overflow";
		}
	case EXCEPTION_FLT_STACK_CHECK: {
			return "Float stack check";
		}
	case EXCEPTION_FLT_UNDERFLOW: {
			return "Float underflow";
		}
	case EXCEPTION_INT_DIVIDE_BY_ZERO: {
			return "Integer divide by zero";
		}
	case EXCEPTION_INT_OVERFLOW: {
			return "Integer Overflow";
		}
	case EXCEPTION_PRIV_INSTRUCTION: {
			return "Priviledged Instruction";
		}
	case EXCEPTION_IN_PAGE_ERROR: {
			return "In-page error";
		}
	case  EXCEPTION_ILLEGAL_INSTRUCTION: {
			return "Illegal instruction";
		}
	case EXCEPTION_NONCONTINUABLE_EXCEPTION: {
			return "Non-continuable exception";
		}
	case EXCEPTION_STACK_OVERFLOW: {
			return "Stack overflow";
		}
	case EXCEPTION_INVALID_DISPOSITION: {
			return "Invalid disposition";
		}
	case EXCEPTION_GUARD_PAGE: {
			return "Guard page";
		}
	case EXCEPTION_INVALID_HANDLE: {
			return "Invalid handle";
		}
	case CONTROL_C_EXIT: {
			return "Control-C";
		}
	case 123:          /* C++ throw exception, not using debugger */
	case 0xe06d7363:  /* C++ throw exception debug mode */ {
			return "C++ throw ";
		}

	default: {
			sprintf_s(msg, sizeof(msg), "unknown exception code = %X ", errorcode);
			break;
		}
	}

	return msg;
}
/* Structured Exception Handler. called if the __except section of a function is 
entered. This can happen for a variety of causes, e.g. signals, divide errors etc. */
int filter(unsigned int code, struct _EXCEPTION_POINTERS *ep)
{
	const char *msg;
	unsigned long flags;

	int code2 = ep->ExceptionRecord->ExceptionCode;
	if (code2 != code) {
		fprintf(stderr, "???\n");
	}
	flags = ep->ExceptionRecord->ExceptionFlags;

	/* convert code to a text message in msg */
	msg = getText(code);
	fprintf_s(stderr, "Entered SEH exception handler: %s ", msg);

	if ((flags && EXCEPTION_NONCONTINUABLE) != 0)
		fprintf_s(stderr, "(non-continuable) \n)");
	else
		fprintf_s(stderr, "(continuable) \n");

	createMiniDump(ep);
	DumpStackTrace(ep);  /* print call stack */
	return EXCEPTION_EXECUTE_HANDLER;  /* continue execution of the __except block */
}

/* interface matcher; separates exception code */
long filter2(struct _EXCEPTION_POINTERS *ep) {
	int code = ep->ExceptionRecord->ExceptionCode;
	return filter(code, ep);
}

/* call this to generate one of a variety of errors; test error handling */
void testerrors(void) {
	printf("Choose an exception type:\n");
	printf("0 - SEH exception (access violation)\n");
	printf("1 - terminate\n");
	printf("2 - unexpected\n");
	printf("3 - pure virtual method call\n");
	printf("4 - invalid parameter\n");
	printf("5 - new operator fault\n");
	printf("6 - SIGABRT\n");
	printf("7 - SIGFPE\n");
	printf("8 - SIGILL\n");
	printf("9 - SIGINT\n");
	printf("10 - SIGSEGV\n");
	printf("11 - SIGTERM\n");
	printf("12 - RaiseException\n");
	printf("13 - throw C++ typed exception\n");
	printf("14 - SEH exception (divide error) \n");
	printf("Your choice >  ");

	int ExceptionType = 0;
	scanf_s("%d", &ExceptionType);

	switch (ExceptionType) 	{
	case 0: /* SEH */	{
		// Access violation
		int *p = 0;
#pragma warning(disable : 6011)
		// warning C6011: Dereferencing NULL pointer 'p'
		*p = 0;
#pragma warning(default : 6011)   
		break;
	}
	case 1: /* terminate */	{
		// Call terminate
		terminate();
		break;
	}
	case 2: /* unexpected */ {
		// Call unexpected
		unexpected();
		break;
	}
	case 3: /* pure virtual method call */	{
		// pure virtual method call
		//CDerived derived;
		break;
	}
	case 4: /* invalid parameter */ {
		char* formatString;
		// Call printf_s with invalid parameters.
		formatString = NULL;
#pragma warning(disable : 6387)
		// warning C6387: 'argument 1' might be '0': this does
		// not adhere to the specification for the function 'printf'
		printf(formatString);
#pragma warning(default : 6387)   

		break;
	}
	case 5: /* new operator fault */ {
		// Cause memory allocation error
		const int BIG_NUMBER = 0x7fffffff;  /* try to allocate 1 terabyte*/
		size_t *pi = new size_t[BIG_NUMBER];
		delete[]pi;
		break;
	}
	case 6: /* SIGABRT */ {
		// Call abort
		abort();
		break;
	}
	case 7: /* SIGFPE */ {
		// floating point exception ( /fp:except compiler option)
		double x = 0, y = 1, z;
		z = y / x;  /* generate divide error */
		printf("z= %g \n", z); /* stop optimiser from removing test code */
		break;
	}
	case 8: /* SIGILL */ {
		raise(SIGILL);
		break;
	}
	case 9: /* SIGINT */ {
		raise(SIGINT);
		break;
	}
	case 10: /* SIGSEGV */ {
		raise(SIGSEGV);
		break;
	}
	case 11: /* SIGTERM */ 	{
		raise(SIGTERM);
		break;
	}
	case 12: /* RaiseException */ 	{
		// Raise noncontinuable software exception
		RaiseException(124, EXCEPTION_NONCONTINUABLE, 0, NULL);
		break;
	}

	case 13: /* throw */	{
		// Throw typed C++ exception.
		//throw std::range_error("test for throw");
		throw 13;
		break;
	}
	case 14: /* SEH exception */ {
		int x = 1, y = 0, z;
		z = x / (2 - x - x);  /* divide by zero */
		printf("z=%d", z);
		break;
	}
	default: 	{
		printf("Unknown exception type %d specified. \n", ExceptionType);
	
		getchar();
		break;
	}

		break;
	}
}

#ifdef __cplusplus
}
#endif
