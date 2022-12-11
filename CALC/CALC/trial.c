/* trial.c */

#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_NONSTDC_NO_DEPRECATE

#include <stdio.h>
#include <io.h>
//#include <stdlib.h>
#include <setjmp.h>
//#include <string.h>
#include <signal.h>
#include <windows.h> // for EXCEPTION_ACCESS_VIOLATION
#include <excpt.h>
#include <float.h>

#include "integer.h"
//#include "fun.h"
/* cant include fun.h due to name clash; copy required declarations below */
void FREEMPIBANK();
void Flush();
void readme();
void clean_symtab();
void INITMPIBANK();
void init();
void Parse(char *s);
void testerrors(void);
void SetProcessExceptionHandlers();
void SetThreadExceptionHandlers();


#include "primes.h"
#define MAXLEN 1000 /* Maximum length of line */
char *progname;
jmp_buf begin;

#ifdef DEBUG
extern long int nettbytes;
#endif

/* Structured Exception Handler */
int filter(unsigned int code, struct _EXCEPTION_POINTERS *ep);

int main(int argc, char ** argv)
{
	/* __try and __except are the Microsoft way of catching errors in C programs*/
	__try
	{
		char *ns = NULL;
		char *command = NULL;
		static char quit_str[] = "exit";
		static char exit_str[] = "quit";
		static char help_str[] = "help";
		int quit;
		char *cps, *cpe;
		int tty;

		INITMPIBANK();
#ifdef UNIX
		setlinebuf(stdout);
#endif
#ifdef DEBUG
		printf("before testing: nettbytes = %ld\n", nettbytes);
#endif
		//tty = isatty(0);                /* is stdin a terminal ?*/
		tty = _isatty(_fileno(stdin));    /* is stdin a terminal ?*/
		progname = argv[0];
		init();             /*install functions in symbol table */

		/* detect floating point errors */
		unsigned int control_word;
		int err = _controlfp_s(&control_word, _EM_INEXACT | _EM_UNDERFLOW, MCW_EM);
		/* allow hardware FP exceptions except inexact and underflow which are
		considered to be normal, not errors. */
		if (err) {
			printf_s("could not set FP control word\n");
			return -1;
		}

		if (ns == NULL)
			if (!(ns = malloc(sizeof(char) * 1000))) {
				printf("Not enough memory\n");
				exit(-1); /* return failure */
			}
		if (command == NULL)
			if (!(command = malloc(sizeof(char) * 1000))) {
				printf("Not enough memory\n");
				exit(-1); /* return failure */
			}
		if (argc == 2) {
			int i;
			char dummy[4096];
			for (i = 0; i < 4095; i++)dummy[i] = 0;
			strcpy(dummy, argv[1]);
			for (i = 0; i < strlen(dummy); i++)if (dummy[i] == ';')dummy[i] = '\n';
			dummy[strlen(argv[1])] = '\n';
			Parse(dummy);
			exit(0);
		}
		printf("\n");
		printf("\t\t\tCALC\n");
		printf("\n");
		printf("\t\tA NUMBER THEORY CALCULATOR\n");
		printf("\t\tK.R.MATTHEWS, 20th September 2016\n\n");
		printf("\n\t\t see http://www.numbertheory.org/calc/krm_calc.html for more info \n");
		;
		printf("Type exit to quit, help for information:\n");

		//SetProcessExceptionHandlers();
		//SetThreadExceptionHandlers();
		//testerrors();    /* temporary */

		while (1)
		{
			setjmp(begin);
			quit = 0;

			ns[0] = 0;
			fflush(stdout);
			if (tty)printf("> ");
			if (scanf("%999[^\n]", ns) == EOF)break;	/* exit on end of file */
			strcpy(command, ns);
			Flush();
			cps = ns;
			while (1) {
				cpe = cps;
				while (*cpe && *cpe != ';')
					cpe++;
				*cpe = '\n';
				*(cpe + 1) = '\0';
				if (strncmp(cps, quit_str, 4) == 0 || strncmp(cps, exit_str, 4) == 0) {
					quit = 1;
					break;
				}
				if (strncmp(cps, help_str, 4) == 0) {
					quit = 2;
					break;
				}
				Parse(cps);
				strcpy(ns, command);
				if (!(*cpe)) break;
				cps = cpe + 1;
				while (*cps && (*cps == ' ' || *cps == '\n'))
					cps++;
				if (!(*cps))
					break;
			}
			if (quit == 1)
				break;
			if (quit == 2) {
				readme();   // help message
				continue;
			}
		}
		clean_symtab();
		FREEMPIBANK();
#ifdef DEBUG
		/*  FREEBANKI(); */
		printf("after testing: nettbytes = %ld\n", nettbytes);
#endif
		exit(0); /* return success */
	}

	/* the filter function does general diagnostic collecting.  */
	__except (filter(GetExceptionCode(), GetExceptionInformation())) {
		/* the filter function may or may not return control to allow the code below 
		to execute. Application-specific cleaning up could go here */
		system("PAUSE");   /* press any key to continue */
		return EXIT_FAILURE;
	}
}
void warning(char *s, char *t) {
	if (t)
		fprintf(stderr, "%s ", t);
	fprintf(stderr, "%s\n", s);
}

int execerror(char *s, char *t) {
	warning(s, t);
	longjmp(begin, 0);
}


