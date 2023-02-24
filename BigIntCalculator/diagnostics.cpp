#include "pch.h"

#include <setjmp.h>
#include <signal.h>

#include <excpt.h>
#include <new.h>
#include <eh.h>
#include <sstream>
#include <string.h>

#include <Psapi.h>
/* To ensure correct resolution of symbols, add Psapi.lib and dbghelp.lib
 In visual studio:
 project -> properties -> linker -> input -> additional dependencies */
#include <dbghelp.h>
#include <algorithm>
#include <iterator>

 
#include <minidumpapiset.h>
#include "diagnostic.h"

bool breakSignal = false;
extern HWND handConsole;      /* handle to console window */

/* used for minidump */
struct module_data {
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

/* walkbaack through call stack */
DWORD StackTrace(const EXCEPTION_POINTERS *ep)
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
		std::cerr << "dumpstack failure; Unable to initialize symbol handler\n";
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
				builder << "  " << line.FileName << " (" << line.LineNumber << ")\n";
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
			builder << "(No Symbols: PC.offset == 0)\n";
		if (!StackWalk64(image_type, process, hThread, &frame, context, NULL,
			SymFunctionTableAccess64, SymGetModuleBase64, NULL))
			break;
		if (++n > 20)
			break;
	} while (frame.AddrReturn.Offset != 0);


	SymCleanup(process);

	// Display the string:
	fprintf_s(stderr, " Stack Trace; \n%s", builder.str().c_str());
	return EXCEPTION_EXECUTE_HANDLER;
}

/* does the same as StackTrace but doesn't need the exception pointer 
parameter*/
DWORD StackTrace2(void)
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
		std::cerr << "dumpstack failure; Unable to initialize symbol handler\n";
		exit(EXIT_FAILURE);
	}

	DWORD symOptions = SymGetOptions();
	symOptions |= SYMOPT_LOAD_LINES | SYMOPT_UNDNAME;
	SymSetOptions(symOptions); /* set load lines & undecorated names */

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
			sprintf_s(builder+bIx, bIxMax-bIx, "%-16s ", fnName.c_str());
			bIx = strlen(builder);
			if (SymGetLineFromAddr64(process, frame.AddrPC.Offset, &offset_from_symbol, &line)) {
				strcat_s(builder, bIxMax, line.FileName);
				bIx = strlen(builder);
				sprintf_s(builder+bIx, bIxMax-bIx, " (%4d) \n", line.LineNumber);
				bIx = strlen(builder);
			}
			else strcat_s(builder, bIxMax, "\n");

			if (fnName == "main")
				break;    // exit when top level in stack is reached

			//if (fnName == "RaiseException") {
			//	// This is what we get when compiled in Release mode:
			//	printf_s("%s", "the program has crashed.\n\n");
			//	break;
			//}
			bIx = strlen(builder);
		}
		else {
			strcat_s(builder, bIxMax - bIx, "(No Symbols: PC.offset == 0)\n");
			bIx = strlen(builder);
		}
		if (!StackWalk64(image_type, process, hThread, &frame, &ct, NULL,
			SymFunctionTableAccess64, SymGetModuleBase64, NULL))
			break;
		if (++n > 20)
			break;   // exit if more than 20 stack entries
	} while (frame.AddrReturn.Offset != 0);

	SymCleanup(process);

	// Display the string:
	fprintf_s(stderr, "%s%s", " Stack Trace: \n",
		builder);
	free(builder);
	return EXCEPTION_EXECUTE_HANDLER;
}

/* to handle calls to library functions that have invalid parameters 
(only works when compiled in debug mode???) */
void InvalidParameterHandler(const wchar_t* expression,
	const wchar_t* function,
	const wchar_t* file,
	unsigned int line,
	size_t pReserved)
{
	wprintf(L"Invalid parameter detected in: \n"
		    L"function    %s.\n"
		    L" File:      %s \n"
		    L"Line:       %d\n", function, file, line);
	wprintf(L"Expression: %s\n", expression);
	StackTrace2();
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
						   Also occurs if, when attempting to convert to an integer,
						   it cannot be converted because it is a NaN, infinite or 
						   outside the integer range.
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
void SigfpeHandler(int sig, int num) {
	const char *msg = NULL;
	if (num == FPE_INEXACT || num == FPE_UNDERFLOW)
		return;  /* ignore these signals */

	msg = FPEcodeToTxt(num);
	/* strictly speaking, should not call I/O routines here, but it seems to work OK */
	fprintf_s(stderr, "\nreceived FPE signal %d %s \n", sig, msg);
	StackTrace2();
	EXCEPTION_POINTERS *ep = (EXCEPTION_POINTERS *)_pxcptinfoptrs;
	createMiniDump(ep);
	Beep(1200, 1000);   // beep at 1200 Hz for 1 second
	std::cout << "Press ENTER to continue...";
#undef max   /* remove max defined in windows.h because of name clash */
	std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	return;
}

void SigabrtHandler(int sig) {
	fprintf(stderr, "Abort signal occurred\n");
	StackTrace2();
	EXCEPTION_POINTERS *ep = (EXCEPTION_POINTERS *)_pxcptinfoptrs;
	createMiniDump(ep);
	Beep(1200, 1000);   // beep at 1200 Hz for 1 second
	Beep(2400, 1000);   // beep at 2400 Hz for 1 second
	Sleep(2000);   /* wait 2 seconds*/
	std::cout << "Press ENTER to continue...";
	std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	// Clear EXECUTION_STATE flags to allow the system to idle to sleep normally.
	SetThreadExecutionState(ES_CONTINUOUS);
	exit(EXIT_FAILURE);
}

void SigillHandler(int sig) {
	fprintf(stderr, "Illegal Instruction \n");
	StackTrace2();
	//__debugbreak();     // try to enter debuggger (to look at call stack)
	EXCEPTION_POINTERS *ep = (EXCEPTION_POINTERS *)_pxcptinfoptrs;
	createMiniDump(ep);
	Beep(1200, 1000);   // beep at 1200 Hz for 1 second
	std::cout << "Press ENTER to continue...";
	std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	// Clear EXECUTION_STATE flags to allow the system to idle to sleep normally.
	SetThreadExecutionState(ES_CONTINUOUS);
	exit(EXIT_FAILURE);
}

void SigsegvHandler(int sig) {
	fprintf(stderr, "Segment Violation \n");
	StackTrace2();
	//__debugbreak();     // try to enter debuggger (to look at call stack)
	EXCEPTION_POINTERS *ep = (EXCEPTION_POINTERS *)_pxcptinfoptrs;
	createMiniDump(ep);
	Beep(1200, 1000);   // beep at 1200 Hz for 1 second
	std::cout << "Press ENTER to continue...";
	std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	// Clear EXECUTION_STATE flags to allow the system to idle to sleep normally.
	SetThreadExecutionState(ES_CONTINUOUS);
	exit(EXIT_FAILURE);
}

/* handles ctrl-C (SIGINT) & ctrl-Break (SIGBREAK) 
Also activated when the console window is closed 
When a CTRL+C interrupt occurs, Win32 operating systems generate a new thread to 
specifically handle that interrupt. This can cause a single-thread application 
to become multithreaded and cause unexpected behavior. */
void SigintHandler(int sig) {
	breakSignal = true;     /* Main program can check this. Allows it to minimise
							strange behaviour before it terminates */
	/* use a message box so as not to mess up any output from main program which 
	is still running. */
	int r = MessageBoxA(handConsole, "Press YES to terminate program, NO to continue", 
		"Interrupt", MB_YESNO);
	switch (r){
	case IDYES: {  /* YES selected */
		// Clear EXECUTION_STATE flags to allow the system to idle to sleep normally.
		SetThreadExecutionState(ES_CONTINUOUS);
		exit (EXIT_FAILURE);
	}
	case IDNO: { /* NO selected */
		/* re-register signal handler */
		signal(SIGINT, SigintHandler);       // interrupt (Ctrl - C)
		signal(SIGBREAK, SigintHandler);     // Ctrl - Break sequence
		return;  /* main program continues, but cannot get any more input from stdin */
	}

	default:
		abort();    /* should never get to here!! */
	}
}

void SigtermHandler(int sig) {
	fprintf(stderr, "termination \n");
	StackTrace2();
	//__debugbreak();     // try to enter debuggger (to look at call stack)
	EXCEPTION_POINTERS *ep = (EXCEPTION_POINTERS *)_pxcptinfoptrs;
	createMiniDump(ep);
	Beep(1200, 1000);   // beep at 1200 Hz for 1 second
	std::cout << "Press ENTER to continue...";
	std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	// Clear EXECUTION_STATE flags to allow the system to idle to sleep normally.
	SetThreadExecutionState(ES_CONTINUOUS);
	exit(EXIT_FAILURE);
}

int handle_program_memory_depletion(size_t memsize)
{
	// Your code
	std::cout << "Press ENTER to continue...";
	std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	return 0;
}

/* set up exception handlers; caller specifies which are required using flags */
void SetProcessExceptionHandlers(flags f)
{
//#ifndef _DEBUG
	/* send errors to calling process */
	SetErrorMode(SEM_FAILCRITICALERRORS | SEM_NOGPFAULTERRORBOX);
//#endif

/* if you use C++ exception handling: install a translator function
with set_se_translator() or use SetUnhandledExceptionFilter(). In the 
context of that function (but *not* afterwards), you can either do 
your stack dump, or save the CONTEXT record as a local copy. Note 
that you must do the stack dump at the earliest opportunity, to avoid 
the interesting stack-frames being gone by the time you do the dump.  */
	if (f.UnEx) {
		/* The system does display the critical-error-handler message box.*/
		SetErrorMode(SEM_FAILCRITICALERRORS); 
		// Suppress the abort message (debug only?, has no effect in release)
		_set_abort_behavior(0, _WRITE_ABORT_MSG);
		SetUnhandledExceptionFilter(filter2);
	}

	/* Catch pure virtual function calls.
	Because there is one _purecall_handler for the whole process, 
	calling this function immediately impacts all threads. The last 
	caller on any thread sets the handler. 
	http://msdn.microsoft.com/en-us/library/t296ys27. */
	//if (f.PureCall)
		//_set_purecall_handler(PureCallHandler);

	_invalid_parameter_handler oldHandler, newHandler;
	if (f.InvParam) {
		newHandler = InvalidParameterHandler;
		oldHandler = _set_invalid_parameter_handler(newHandler);
	}

	// Catch new operator memory allocation exceptions
	if (f.New)
		_set_new_handler(handle_program_memory_depletion);

	// Set up C++ signal handlers

	if (f.abort) {
		// Suppress the abort message (debug only?, has no effect in release)
		_set_abort_behavior(0, _WRITE_ABORT_MSG);
		// Catch an abnormal program termination
		signal(SIGABRT, SigabrtHandler);
	}

	/* Catch interrupt handler. When a CTRL+C interrupt occurs, Win32 operating 
	systems generate a new thread to specifically handle that interrupt. This can 
	cause a single-thread application, to become multithreaded and cause 
	unexpected behavior.
	*/
	if (f.interrupt) {
		signal(SIGINT,   SigintHandler);     // interrupt (Ctrl - C)
		signal(SIGBREAK, SigintHandler);     // Ctrl - Break sequence
	}

	// Catch a termination request
	if (f.sigterm)
		signal(SIGTERM, SigtermHandler);


	// Catch terminate() calls. 
	// In a multithreaded environment, terminate functions are maintained 
	// separately for each thread. Each new thread needs to install its own 
	// terminate function. Thus, each thread is in charge of its own termination handling.
	// http://msdn.microsoft.com/en-us/library/t6fk7h29.aspx
	/*if (f.term)
		set_terminate(TerminateHandler);*/

	// Catch unexpected() calls.
	// In a multithreaded environment, unexpected functions are maintained 
	// separately for each thread. Each new thread needs to install its own 
	// unexpected function. Thus, each thread is in charge of its own unexpected handling.
	// http://msdn.microsoft.com/en-us/library/h46t5b69.aspx  
	/*if (f.unexpected)
		set_unexpected(UnexpectedHandler);*/

	/* Catch a floating point error.
	Note: some errors e.g. floating point errors will be caught by the filter
	function if a signal handler is not registered. If both a SIGFPE signal
	function and a filter function are used, both seem to get called for each 
	floating point error. */
	typedef void(*sigh)(int); /* force compiler to accept SigfpeHandler although
							  it actually has 2 parameters, not 1 */
	if (f.sigfpe)
		signal(SIGFPE, (sigh)SigfpeHandler);

	// Catch an illegal instruction
	if (f.sigill)
		signal(SIGILL, SigillHandler);

	// Catch illegal storage access errors
	if (f.sigsegv)
		signal(SIGSEGV, SigsegvHandler);

}

/* pExcPtrs is a pointer to a EXCEPTION_POINTERS structure 
describing the client exception that caused the minidump to be generated. If 
the value of this parameter is NULL, no exception information is included in 
the minidump file 
the dump file will be placed in the TEMP directory.
The file name is crashdumpHHMMSS.dmp so that each dump has a unique name and will
not overwrite earlier dumps. */
void createMiniDump(const EXCEPTION_POINTERS* pExcPtrs)
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
	char buf[1024];
	DWORD rv = GetTempPathA(sizeof(buf), buf);
	if (rv > MAX_PATH || rv == 0)
		/* error - couldn't get path*/
		return;

	strcat_s(buf, "crashdump");
	char timestamp[10];   // time in format hhmmss
	struct tm newtime;
	const time_t current = time(NULL);  // time as seconds elapsed since midnight, January 1, 1970
	localtime_s(&newtime, &current);    // convert time to tm structure
	/* convert time to hhmmss */
	strftime(timestamp, sizeof(timestamp), "%H%M%S", &newtime);
	strcat_s(buf, timestamp);
	strcat_s(buf, ".dmp");
	// Create the minidump file
	hFile = CreateFileA(
		buf,
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
	mei.ExceptionPointers = (EXCEPTION_POINTERS *)pExcPtrs;
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
		hProcess,      /* A handle to the process for which the information is 
					   to be generated.*/
		dwProcessId,   /* The identifier of the process for which the information 
					   is to be generated*/
		hFile,         /* A handle to the file in which the information is to be written*/
		MiniDumpWithDataSegs,  /* includes all of the data sections from loaded modules 
							    in order to capture global variable contents. */
		&mei,          /* A pointer to a MINIDUMP_EXCEPTION_INFORMATION structure 
					   describing the client exception that caused the minidump 
					   to be generated. If the value of this parameter is NULL, 
					   no exception information is included in the minidump file*/
		NULL,          /* no user-defined information is included in the minidump file */
		&mci);         /* no callbacks are performed */

	if (!bWriteDump)
	{
		// Error writing dump.
		return;
	}

	// Close file
	CloseHandle(hFile);

	// Unload dbghelp.dll
	FreeLibrary(hDbgHelp);
	fprintf_s(stderr, "minidump file created; %s \n", buf);
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
			return "Privileged Instruction";
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

/* Structured Exception Handler Filter. called if the __except section of a 
function is entered. This can happen for a variety of causes, e.g. signals, 
divide errors etc. 

return Value 	                   Meaning
EXCEPTION_EXECUTE_HANDLER 0x1      Return from UnhandledExceptionFilter and
								   execute the associated exception handler.
								   This usually results in process termination.

EXCEPTION_CONTINUE_EXECUTION       Return from UnhandledExceptionFilter and
0xffffffff                         continue execution from the point of the
								   exception. Note that the filter function is
								   free to modify the continuation state by
								   modifying the exception information supplied
								   through its LPEXCEPTION_POINTERS parameter.
								   BUT in practise control jumps straight to main.

EXCEPTION_CONTINUE_SEARCH 0x0      Proceed with normal execution of
								   UnhandledExceptionFilter. That means obeying
								   the SetErrorMode flags, or invoking the
								   Application Error pop-up message box.
*/
int filter(unsigned int code, const struct _EXCEPTION_POINTERS *ep)
{
	const char *msg;
	unsigned long flags;

	// Clear EXECUTION_STATE flags to allow the system to idle to sleep normally.
	SetThreadExecutionState(ES_CONTINUOUS);

	flags = ep->ExceptionRecord->ExceptionFlags;

	/* convert code to a text message in msg */
	msg = getText(code);
	fprintf_s(stderr, "Entered SEH exception handler: %s ", msg);

	if ((flags && EXCEPTION_NONCONTINUABLE) != 0)
		fprintf_s(stderr, "(non-continuable) \n)");
	else
		fprintf_s(stderr, "(continuable) \n");

	createMiniDump(ep);
	StackTrace(ep);  /* print call stack */

	/* floating point error handling is well-defined even without any exception
	handling, so it should be relatively safe to continue execution 
	BUT in practise control jumps straight to main. */
	if (code == EXCEPTION_FLT_DENORMAL_OPERAND ||
		code == EXCEPTION_FLT_DIVIDE_BY_ZERO ||
		code == EXCEPTION_FLT_INEXACT_RESULT ||
		code == EXCEPTION_FLT_OVERFLOW ||
		code == EXCEPTION_FLT_UNDERFLOW ||
		code == EXCEPTION_FLT_INVALID_OPERATION)
		return EXCEPTION_CONTINUE_EXECUTION;

	return EXCEPTION_EXECUTE_HANDLER;  /* continue execution of the __except block */
}

/* interface matcher; gets the exception code from the exception record. 
This version is useful for SetUnhandledExceptionFilter() 
*/
long filter2(struct _EXCEPTION_POINTERS *ep) {
	int rcode;
	int code = ep->ExceptionRecord->ExceptionCode;

	//if (code == EXCEPTION_FLT_DENORMAL_OPERAND ||
	//	code == EXCEPTION_FLT_DIVIDE_BY_ZERO ||
	//	code == EXCEPTION_FLT_INEXACT_RESULT ||
	//	code == EXCEPTION_FLT_OVERFLOW ||
	//	code == EXCEPTION_FLT_UNDERFLOW ||
	//	code == EXCEPTION_FLT_INVALID_OPERATION)
	//	return EXCEPTION_CONTINUE_EXECUTION;

	rcode = filter(code, ep);
	return EXCEPTION_EXECUTE_HANDLER;  /* override return code from filter */
}

/* call this to generate one of a variety of errors; test error handling */
void testerrors(void) {

	// Clear EXECUTION_STATE flags to allow the system to idle to sleep normally.
	SetThreadExecutionState(ES_CONTINUOUS);

	printf_s("Choose an exception type:\a\n");
	printf_s("0 - SEH exception (access violation)\n");
	printf_s("1 - terminate\n");
	printf_s("2 - unexpected\n");
	printf_s("3 - pure virtual method call\n");
	printf_s("4 - invalid parameter\n");
	printf_s("5 - new operator fault\n");
	printf_s("6 - SIGABRT\n");
	printf_s("7 - SIGFPE (divide error) \n");
	printf_s("8 - SIGILL\n");
	printf_s("9 - SIGINT\n");
	printf_s("10 - SIGSEGV\n");
	printf_s("11 - SIGTERM\n");
	printf_s("12 - RaiseException\n");
	printf_s("13 - throw C++ typed exception\n");
	printf_s("14 - SEH exception (divide error) \n");
	printf_s("15 - SIGFPE (overflow) \n");
	printf_s("16 - SIGFPE (invalid) \n");
	printf_s("17 - SIGFPE (raise) \n");
	printf_s("18 - SIGBREAK \n");
	printf_s("Your choice >  ");

	int ExceptionType = 9999;
	char buffer[256];
	fgets(buffer, sizeof(buffer), stdin);
	int rv = sscanf_s(buffer, "%d", &ExceptionType);


	switch (ExceptionType) {
	case 0: /* SEH */ {
			// Access violation
			int *p = 0;
#pragma warning(disable : 6011)
			// warning C6011: Dereferencing NULL pointer 'p'
			*p = 0;
#pragma warning(default : 6011)   
			break;
		}
	case 1: /* terminate */ {
			terminate();  // Call terminate
			break;
		}
	case 2: /* unexpected */ {
			unexpected();    // Call unexpected
			break;
		}
	case 3: /* pure virtual method call */ {
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
			int rc = printf_s(formatString);
#pragma warning(default : 6387)   
			printf_s("return code from printf_s is %d \n", rc);
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
			abort();    // Call abort
			break;
		}
	case 7: /* SIGFPE */ {
			// floating point exception ( /fp:except compiler option)
			double x = 0, y = 1, z;
			z = y / x;  /* generate divide error */
			printf_s("z= %g \n", z); /* stop optimiser from removing test code */
			y++;
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
	case 11: /* SIGTERM */ {
			raise(SIGTERM);
			break;
		}
	case 12: /* RaiseException */ {
			// Raise noncontinuable software exception
			RaiseException(124, EXCEPTION_NONCONTINUABLE, 0, NULL);
			break;
		}

	case 13: /* throw */ {
			// Throw typed C++ exception.
			throw std::range_error("test for throw");
			//throw 13;
			break;
		}
	case 14: /* SEH exception */ {
			int x = 1, y = 0, z;
			z = x / (2 - x - x);  /* divide by zero */
			printf_s("z=%d\n", z);
			break;
		}
	case 15: /* floating point overflow */ {
			double z = std::pow(DBL_MAX, 2);    /* generate overflow */
			printf_s("z=%g\n", z);                /* if overflow not trapped z = inf. */
			z = DBL_MAX;
			printf_s("z=%g\n", z);
			z *= 100.0;                           /* if overflow not trapped z = inf. */
			z -= 1.0;
			printf_s("z=%g\n", z);
			break;
		}
	case 16: /* floating point invalid 
			 This exception is raised if the given operands are invalid for the 
			 operation to be performed. Examples are (see IEEE 754, section 7):

             Addition or subtraction: &infin; - &infin;. (But &infin; + &infin; = &infin;).
             Multiplication: 0  x infin .
             Division: 0/0 or  infin/infin .
             Remainder: x REM y, where y is zero or x is infinite.
             Square root if the operand is less then zero. More generally, any 
			 mathematical function evaluated outside its domain produces this exception.
             Conversion of a floating-point number to an integer or decimal string, 
			 when the number cannot be represented in the target format (due to 
			 overflow, infinity, or NaN).
             Conversion of an unrecognizable input string.
             Comparison via predicates involving < or >, when one or other of the
			 operands is NaN. */ {
			double z = std::acos(2);
			printf_s("z=%g\n", z);    /* if error not trapped, Z = NaN */
			break;
		}
	case 17: /* raise FPE */ {
			raise(SIGFPE);  
			break;
		}
	case 18: /* raise Break signal */ {
		raise(SIGBREAK);
		break;
	}

	default: {
			printf_s("Unknown exception type %d specified. \n", ExceptionType);
			break;
		}

		break;
	}
}